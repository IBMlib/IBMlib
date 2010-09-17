ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Implements internal state update interface (module particle state)
c     and active motion components of tracer dynamics for sandeel eggs/larvae   
c
c     init_particle_state:       parameters from input file
c     update_tracer_state(...):  arguments specific to each provider
c     close_particle_state:      civilized module close-down
c
c     write_particle_state(...): write state to file(s)   
c
c
c     TODO: * more accurate egg/larvae emission (rate neglects time step)
c           * daily avg temp for devel, not point temp
c           * use named integer constants instead of integer literals 
c             in fstate/istate data manipulations
c
c     NOW: deactivate backtracking (dt<0)
c          deactivate option of initialization in user provided state 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state
      use physical_fields
      use particle_tracking
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use time_tools
      implicit none
      private                   ! default visibility



c =================================================
c ========       module data section       ========
c =================================================
      integer, parameter          :: verbose = 0    ! control output volume

      
      type(clock), pointer        :: current_time

c -------------- egg production data --------------
      integer                     :: negg_sources            ! convenience
      type(emission_box), pointer :: egg_source(:)
      integer                     :: eggID          ! counter to keep track of all (living/dead) emitted eggs 

c -------------- stochastic egg growth data --------------
c      
      real    :: hatch_begin(2)              ! days since fertilization (A,k0 > 0)
      real    :: hatch_mid(2)                ! days since fertilization (B,k1 > 0)
      real    :: hatch_end(2)                ! days since fertilization (C,k2 > 0)

c -------------- larvae production data --------------
      integer                     :: nlarv_sources            ! convenience
      type(emission_box), pointer :: larv_source(:)
      integer                     :: larvID   ! counter to keep track of all (living/dead) emitted larvae


c     the internal state descriptors of larvae emitted from
c     each larvae source
c     meaning of larvsource_state may depend on time direction
c     the meaning of larvsource_state is interpreted by larvae_production 
c     real, private, allocatable    :: larvsource_state(:,:)


c -------------- larvae growth data --------------
c    a(T)  = func(lambda(1) + lambda(2) T + lambda(3) T^2)     (T in Celcius)
c    dL/dt = a (L/L_0)^beta (1 - L/L_infty)
c
c    lambda       temperature polynomial coefficients  unit == mm/day/K^n 
c    func         temperature function                 0: func=polynomial  1: func=exp(polynomial)
c    L_0          hatch length                         unit == mm (should be fixed)
c    L_infty      Bertallanfy asymptote                unit == mm (fixed for now)   
c    beta         small larvae growth exponent         unit == 1 
c
      real            :: lambda(3)
      integer         :: tempfunc
      real            :: L_0   
      real            :: beta                           
      real, parameter :: L_infty = 218.0  


c -------------- juvenile transition data --------------
c
c    L_metamorp0 <  L_metamorp1 :  metamorphosis length interval        unit == mm  
c
c    larvae will settle in first encountered habitat when  L_metamorp0 < L < L_metamorp1
c    when L > L_metamorp1 larvae dump, irrespective of its position
c    
c    habitat_map(ix,iy) : habitat number association of this grid cell 
c                         habitat number > 0 means this grid cell is associated with a habitat
c                         habitat_map(ix,iy) == 0 means grid cell (ix,iy) is not a habitat 
c    
c    habitat_map numbes are read from file provided in habitat_map containing lines {ix, iy, habitat_number}
c    habitat_map is only loaded, if  L_metamorp0 <  L_metamorp1 
c
c    check_habitat is a logical switch for (L_metamorp0 <  L_metamorp1)
c
      real                 :: L_metamorp0,  L_metamorp1              
      integer, allocatable :: habitat_map(:,:)             ! habitat_map(nx,ny)
      logical              :: check_habitat
      integer              :: juvID  ! counter to keep track of all (living/dead) emitted eggs 
      

cccccccccccccccccccccccccc  internal state variables cccccccccccccccccccccccccc
c
c     ---------------- common istates -------------------------
c     istate(1)  == individual type  1: egg
c                                    2: larv
c                                    3: juvenile
c                                    4: settled juvenile
c
c                                  < 0: dead (value = negative of last state alive)
c
c     istate(2)  == state of origin (first istate(1) value)
c     istate(3)  == source box of origin (in any state) (0, if not originate from a box)
c     istate(4)  == Julian day of origin (in any state)
c     ----------------- eggs --------------------------
c     istate(5)  == eggID  (unique ID)
c     istate(6)  == Julian day of transition to egg state    
c     ---------------- larvae --------------------------
c     istate(7)  == larvID (unique ID)
c     istate(8)  == Julian day of transition to larval state    
c     --------------- juveniles -------------------------
c     istate(9)  == juvID  (unique ID)
c     istate(10) == Julian day of transition to juvenile state
c     istate(11) == Julian day of settlement
c     istate(12) == settlement habitat number (from habitat map)
c     --------------- death -------------------------
c     istate(13) == Julian day of death
c     ---------------- common fstates -------------------------
c     fstate(1) egg_state/larval_length    0   <= egg_state < 1     1 => hatch
c                                          L_0 <= larval_length  
c     
c     fstate(2) growth_destiny      individual relative development (growth speed seed for stochastic growth)
c     --------------------------------------------------
c     Notes:
c        * at dying, istate(1) -> -istate(1), so state prior to death is stored
c        * fields that do not apply are set to pad value 0 
c     -------------------------------------------------- 
      integer, parameter :: n_int_states  = 13
      integer, parameter :: n_real_states = 2


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc  Define public interface  cccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      public :: init_particle_state
      public :: update_tracer_state
      public :: close_particle_state
      public :: write_particle_state

c =================================================
c ========  module data manipulation part  ========
c =================================================
      contains


      subroutine init_particle_state() 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialize all data concerning higher trophic levels
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer            :: nhab, ihab, ix, iy, idum, nmaxpar
      real               :: rvec2(2)
      character*999      :: fname
c     ---------------------------------------------------------------
      write(*,*) "init_particle_state():"
      
      nmaxpar = 0 ! max number of particles emitted
      eggID   = 1  
      larvID  = 1
      juvID   = 1
      
      call create_emission_boxes("egg_source_box",    egg_source, idum)
      nmaxpar = nmaxpar + idum
      negg_sources  = size(egg_source)
      write(*,*) "found ", negg_sources, "egg_source_boxes"

      call create_emission_boxes("larvae_source_box", larv_source, idum)
      nmaxpar = nmaxpar + idum
      nlarv_sources = size(larv_source)
      write(*,*) "found ", nlarv_sources, "larvae_source_boxes"

c
c     Assume this module is the only injecting particles to 
c     particle stack. Otherwise particle stack initialization must be 
c     coordinated. 
c
      call init_particle_stack(nmaxpar, n_real_states, n_int_states) 


c --- Read egg growth parameters
      call read_control_data(ctrlfile,"egg_hatch_begin", hatch_begin)
      call read_control_data(ctrlfile,"egg_hatch_mid",   hatch_mid)
      call read_control_data(ctrlfile,"egg_hatch_end",   hatch_end)
      write(*,*) "  egg_hatch_begin =", hatch_begin
      write(*,*) "  egg_hatch_mid   =", hatch_mid
      write(*,*) "  egg_hatch_end   =", hatch_end

c --- Read larvae growth parameters
      call read_control_data(ctrlfile,"larvae_temp_coeff",  lambda)
      call read_control_data(ctrlfile,"larvae_temp_func",   tempfunc)
      call read_control_data(ctrlfile,"larvae_hatch_len",   L_0)
      call read_control_data(ctrlfile,"larvae_length_expo", beta)
      write(*,*) "  larvae_temp_coeff    =", lambda
      if     (tempfunc==0) then
          write(*,*) "  larvae temp function is polynomial"
      elseif (tempfunc==1) then
          write(*,*) "  larvae temp function is exp(polynomial)"
      else
          stop "unknown larvae temp function switch"
      endif
      write(*,*) "  larvae_hatch_len     =", L_0
      write(*,*) "  larvae_length_expo   =", beta

c --- Read juvenile transition parameters ---
      call read_control_data(ctrlfile, "larvae_metamorph_len", rvec2)
      L_metamorp0 = rvec2(1)
      L_metamorp1 = rvec2(2)
      write(*,*) "  larvae_metamorph_len =", L_metamorp0, L_metamorp1

      if (L_metamorp0 < L_metamorp1) then
         check_habitat = .true.
      else  
         check_habitat = .false.   ! also covers situation L_metamorp0 == L_metamorp1
      endif

c.....load habitat map, if L_metamorp0 < L_metamorp1 (else ignore it)
      if (check_habitat) then
          write(*,*) "  flexible habitat choice using habitat map"
          allocate( habitat_map(nx,ny) )   ! follows horizontal grid layout
          habitat_map = 0                  ! default: cell is not a habitat
          call read_control_data(ctrlfile,"habitat_map_file", fname)
          open(87, file=fname)
          nhab = 0
          do
             read(87,*,end=455) ix, iy, ihab
             nhab = nhab + 1
             if (ihab < 1) then
                write(*,337) nhab,ihab
                stop
             else
                habitat_map(ix,iy) = ihab
             endif
          enddo
 337      format("habitat",1x,i5,1x,"has invalid number",1x,i5)
 455      close(87)
          write(*,*) "  read", nhab, "habitats from file", trim(fname)
          if (nhab == 0) then
             check_habitat = .false.   
             write(*,*) "  switch off check_habitat "//
     +                  "(using L_metamorp0 only)"
          endif
      else
         write(*,*) "  no habitat check at juvenile transition point"
      endif

      end subroutine init_particle_state



      subroutine close_particle_state()     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "close_particle_state()"   
 
      nullify(current_time)  ! this module does not own this object
      if ( allocated(habitat_map)  ) deallocate( habitat_map )
      if ( associated(egg_source)  ) deallocate( egg_source  )
      if ( associated(larv_source) ) deallocate( larv_source )

      end subroutine close_particle_state



      subroutine update_tracer_state(ctime, dt)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type(clock),target :: ctime      ! must have target attribute
      real               :: dt         ! time step for update
c     ---------------------------------------------------------------
      current_time => ctime     ! store reference

      if (dt<0) stop "negative dt values currently blocked"

      call egg_devel(dt)
      call larvae_devel(dt)
      call egg_production(dt)
      call larvae_production(dt)
      call tracer_metamorphosis(dt)

      end subroutine update_tracer_state


      subroutine write_particle_state()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Dump two files:
c       1) particle state to file provided as "tracer_state_dump_file"
c       2) box emission summary to file provided as "emission_state_file"
c
c     todo: stepped dumps, using internal dump counter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*999      :: fname
      integer            :: iunit, ibox, nemit, nemax, ityp
      integer            :: negg_tot, nlarv_tot
c     ---------------------------------------------------------------
c      -------- 1) dump final simulation output --------
      call read_control_data(ctrlfile,"tracer_state_dump_file", fname)
      write(*,*) "write_particle_state: tracer dump to file ", 
     +           trim(adjustl(fname))
      call find_free_IO_unit(iunit)
      open(unit=iunit, file=fname)
      call write_tracer_data(iunit)
      close(iunit)
c      -------- 2) dump tracer emission data --------  
      negg_tot  = 0
      nlarv_tot = 0
      call read_control_data(ctrlfile,"emission_state_file", fname)
      write(*,*) "write_particle_state: tracer emission dump to file ", 
     +           trim(adjustl(fname))
      call find_free_IO_unit(iunit)
      open(unit=iunit, file=fname)
c.....egg_source loop possibly void
      ityp = 1     ! eggs
      do ibox = 1, size(egg_source)
         call get_current_emissions(egg_source(ibox), nemit, nemax)
         write(iunit,250) ibox, nemit, ityp, nemax
         negg_tot = negg_tot + nemit
      enddo
      write(*,251) negg_tot, "eggs", size(egg_source)
c.....larv_source loop possibly void  
      ityp = 2     ! larvae
      do ibox = 1, size(larv_source)
      call get_current_emissions(larv_source(ibox), nemit, nemax)
         write(iunit,250) ibox, nemit, ityp, nemax
         nlarv_tot = nlarv_tot + nemit
      enddo
      write(*,251) nlarv_tot, "larvae", size(larv_source)
c
      close(iunit)
 250  format(4(i6,1x))
 251  format("write_particle_state: ",i8,1x,a," from ",i8, " sources")     
      end subroutine write_particle_state



      subroutine egg_devel(dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Egg development and production
c
c     dt   : length of time step corresponding to this call (seconds)
c            negative for backtracking (not supported yet)
c
c     Assume module time reference (current_time) has been syncronized
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in) :: dt
      integer          :: last, ipar, ix, iy, ibot
      real             :: tm,a,b,c,g,tau
c     --------------------------------------------------------------
      write(*,*) "egg_devel:"


c------------------------------------------------------------------------------
c    development in dependence of bottom temperature
c    currently, just use cell average temperature
c   
c    tau [unit==days] is the adiabatic egg development time at 
c    this temp for this egg growth seed (fstate(2))
c  
c    Stochastic egg model based on AC/CS/HJ.. (to be published) and Buckley 1984
c    TODO: add plugs for   a + b - 2.*c == 0
c                          a-g          == 0
c------------------------------------------------------------------------------ 
      call get_last_tracer(last)

      do ipar = 1,last
         if (istate(1,ipar)==1) then         ! eggs: istate(1) == 1
            ix     = nint(positions(1,ipar))
            iy     = nint(positions(2,ipar))
            ibot   = bottom_layer(ix,iy)
            tm     = temp(ix, iy, ibot)
            a      = hatch_begin(1)*exp(-hatch_begin(2)*tm)
            b      = hatch_end(1)  *exp(-hatch_end(2)*tm)
            c      = hatch_mid(1)  *exp(-hatch_mid(2)*tm)
            g      = (a*b - c**2)/(a + b - 2.*c)  
            tau    = ((b-g)/(a-g))**fstate(2,ipar)
            tau    = g + (a-g)*tau
            fstate(1,ipar) = fstate(1,ipar) + dt/86400./tau
         endif 
      enddo
      
      end subroutine egg_devel
  



      subroutine larvae_devel(dt) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Temperature dependent larval/juvenile growth
c
c     dt   : length of time step corresponding to this call (seconds)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in) :: dt
      integer          :: last, ipar, ix, iy, iz
      real             :: tm, a, ldot, l

      write(*,*) "larvae_devel()"

c-----------------------------------------------------------------------
c     Deterministic larvae growth model based on HJ thesis and AC/CS/HJ
c     dt may be positive/negative depending on the time direction of the
c     simulation
c     currently, just use cell average temperature
c-----------------------------------------------------------------------      
      call get_last_tracer(last)

      do ipar = 1,last
         if (istate(1,ipar)>1) then         ! larv/juveniles
            ix     = nint(positions(1,ipar))
            iy     = nint(positions(2,ipar))
            iz     = nint(positions(3,ipar))
            tm     = temp(ix, iy, iz)
            a    = lambda(1) + lambda(2)*tm + lambda(3)*tm**2
            if      (tempfunc==0) then
               a    = max(a, 0.0) ! must be positive
            elseif  (tempfunc==1) then
               a    = exp(a)      ! always positive
            else
               stop "unknown tempfunc"
            endif
            l    = fstate(1,ipar)
            ldot = a * (l/L_0)**beta * (1. - l/L_infty)     ! unit == mm/day
            fstate(1,ipar) = fstate(1,ipar) + ldot*dt/86400.  ! Euler increment
         endif
      enddo
 
      end subroutine larvae_devel
      



      subroutine egg_production(dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Inject new eggs into the simulation according to input 
c
c     dt      : length of time step corresponding to this call (seconds)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in) :: dt
      integer          :: ibox, firstpar, lastpar, ipar
      integer          :: julday
c     ------------------------------------------------------------
      write(*,*) "egg_production()"

      call get_julian_day(current_time, julday)

c.....loop over egg sources (if any)

      do ibox = 1, size(egg_source)
         call activate_emission_box(current_time, egg_source(ibox), 
     +                              dt, firstpar, lastpar) 
c.....initialize new eggs
         do ipar = firstpar, lastpar
            call set_state_egg(ipar, ibox, 0.0, .true.)
         enddo   ! ipar    
         if (verbose>0) write(*,228) ibox, 1+lastpar-firstpar

      enddo      ! ibox

 228  format("new eggs in eggsource ",i4," = ",i6)

      end subroutine egg_production



      subroutine larvae_production(dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Inject new larvae into the simulation according to input 
c
c     dt      : length of time step corresponding to this call (seconds)
c               (negative for backtracking)
c           
c     backtracking note: 
c           larvae sources rates are not (currently) considered signed,
c           but implicitly positive. Currently, the code does not feature
c           larvae deletion, but only creation. The sign of dt is neglected
c           in this part.             
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in) :: dt
      integer          :: ibox, firstpar, lastpar, ipar
      integer          :: julday
c     --------------------------------------------------------------------
      write(*,*) "larvae_production()"
     
      call get_julian_day(current_time, julday)

c.....loop over larvae sources (if any)

      do ibox = 1, size(larv_source)
         call activate_emission_box(current_time, larv_source(ibox), 
     +                              dt, firstpar, lastpar) 
      
c.....initialize new larvae
         do ipar = firstpar, lastpar
            call set_state_larv(ipar, ibox, L_0, .true.)
         enddo   ! ipar

         if (verbose>0) write(*,229) ibox, 1+lastpar-firstpar
      enddo      ! ibox

 229  format("new larvae in larvsource ",i4," = ",i6)

      end subroutine larvae_production




      subroutine tracer_metamorphosis(dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Promote individuals to new types, if development 
c     thresholds has been fulfilled
c     
c     Check all particles consequetively
c     egg      -> larvae
c     larvae   -> juvenile
c     juvenile -> settled juvenile
c     juvenile -> death
c
c     Generalized scheme of metamorphosis:
c      if (L_metamorp0 < L < L_metamorp1): settle if habitat_map > 0 at current position
c      if (L > L_metamorp1):               stop, irrespective value of habitat_map at 
c                                          current position
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in) :: dt
      integer          :: ipar, julday, ix, iy, last
      integer          :: newlarv, newjuv, newsett, newdied
      logical          :: settle, die, consider_settling, must_settle 
c     --------------------------------------------------------------------
      write(*,*) "tracer_metamorphosis()"
    
c  --- initialize counters. Notice that newjuv may also appear as settlers/diers

      newlarv = 0   ! number of egg      -> larvae
      newjuv  = 0   ! number of larvae   -> juvenile
      newsett = 0   ! number of juvenile -> settled juvenile
      newdied = 0   ! number of juvenile -> death

      call get_julian_day(current_time, julday)
      call get_last_tracer(last)

      do ipar = 1,last

         if (istate(1,ipar)==4) cycle ! next particle, do nothing for settled

c        -----------------------------
c        egg -> larvae transition
c        -----------------------------
         if ((istate(1,ipar)==1).and.(fstate(1,ipar)>1.0d0)) then         
            call set_state_larv(ipar, 0, L_0, .false.)
            newlarv = newlarv + 1
            cycle               ! next particle, do not check other conditions
         endif

c        -----------------------------
c        larvae -> juvenile transition 
c        -----------------------------
c        do not cycle, but consider settling/dying immeadiately
         if ((istate(1,ipar)==2).and.(fstate(1,ipar)>L_metamorp0)) then         
            call set_state_juvenile(ipar, 0, fstate(1,ipar), .false.)
            newjuv = newjuv + 1
         endif

c        -----------------------------
c        settling/dying transition 
c        -----------------------------
         settle = .false.
         die    = .false.

         if (istate(1,ipar)==3) then 
           consider_settling = (fstate(1,ipar) > L_metamorp0)
           must_settle       = (fstate(1,ipar) > L_metamorp1)   
           if (consider_settling) then
              if (check_habitat) then
                 ix     = nint(positions(1,ipar))
                 iy     = nint(positions(2,ipar))
                 settle = (habitat_map(ix,iy) > 0) 
              else
                 settle = .true. ! we have no knowledge of habitats, mark as settle
              endif
              if (must_settle) die = .not.settle
           endif          
         endif
c        ---- service action flags ----
c             
         if (die) then
            call set_state_die(ipar)
            newdied = newdied + 1
            cycle
         endif

         if (settle) then
           call set_state_settle(ipar, .false.)
           newsett = newsett + 1
           cycle
         endif

      enddo   ! ipar
            
c --- report transformations
      if (newlarv > 0) write(*,422) "egg -> larvae", newlarv
      if (newjuv  > 0) write(*,422) "larvae   -> juvenile", newjuv
      if (newsett > 0) write(*,422) "juvenile -> settled", newsett
      if (newdied > 0) write(*,422) "juvenile -> death", newdied

 422  format("tracer_metamorphosis:",1x,a,":",1x,i5)

c...........backtracking: larvea -> eggs
c                         currently, place the egg at the bottum
cBACKTRACK            if ((dt<0).and.(ttl(zal1) .lt. L_0)) then

      end subroutine tracer_metamorphosis

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccc    specific particle transformations   cccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  set_state_X(ipar, ...) will set istate/fstate for particle 
c  ipar corresponding to a transformation/initialization into
c  particle state  X
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine set_state_egg(ipar,ibox,state,newseed)
c     -----------------------------------------------
      integer, intent(in) :: ipar,ibox
      real, intent(in)    :: state
      logical, intent(in) :: newseed
      integer             :: julday
c     ----------------------------------------------- 
      call get_julian_day(current_time, julday)
      mobility(:,ipar) = 0         ! immobilize
      istate(1,ipar)   = 1         ! current type 
      if (ibox>0) then             ! this signal a new particle
         istate(2,ipar)   = 1      ! type of origin 
         istate(3,ipar)   = ibox   ! box of origin
         istate(4,ipar)   = julday ! day of origin
      endif
      istate(5,ipar)   = eggID     ! unique ID number
      eggID = eggID + 1
      istate(6,ipar)   = julday    ! transition to egg state
      fstate(1,ipar)   = state     ! egg maturization state
      if (newseed) call random_number(fstate(2,ipar)) ! growth speed seed for stochastic growth
      end subroutine set_state_egg


      subroutine set_state_larv(ipar,ibox,state,newseed)
c     -----------------------------------------------
      integer, intent(in) :: ipar,ibox
      real, intent(in)    :: state
      logical, intent(in) :: newseed
      integer             :: julday
c     ----------------------------------------------- 
      call get_julian_day(current_time, julday)
      mobility(:,ipar) = 1         ! mobilize
      istate(1,ipar)   = 2         ! current type 
      if (ibox>0) then             ! this signal a new particle
         istate(2,ipar)   = 2      ! type of origin 
         istate(3,ipar)   = ibox   ! box of origin
         istate(4,ipar)   = julday ! day of origin
      endif
      istate(7,ipar)   = larvID     ! unique ID number
      larvID = larvID + 1
      istate(8,ipar)   = julday    ! transition to larv state
      fstate(1,ipar)   = state     ! larv maturization state
      if (newseed) call random_number(fstate(2,ipar)) ! growth speed seed for stochastic growth      
      end subroutine set_state_larv


      subroutine set_state_juvenile(ipar,ibox,state,newseed)
c     -----------------------------------------------
      integer, intent(in) :: ipar,ibox
      real, intent(in)    :: state
      logical, intent(in) :: newseed
      integer             :: julday
c     ----------------------------------------------- 
      call get_julian_day(current_time, julday)
      mobility(:,ipar) = 1         ! mobilize
      istate(1,ipar)   = 3         ! current type 
      if (ibox>0) then             ! this signal a new particle
         istate(2,ipar)   = 3      ! type of origin 
         istate(3,ipar)   = ibox   ! box of origin
         istate(4,ipar)   = julday ! day of origin
      endif
      istate(9,ipar)   = juvID     ! unique ID number
      juvID = juvID + 1
      istate(10,ipar)   = julday   ! transition to juv state
      fstate(1,ipar)    = state    ! juv maturization state
      if (newseed) call random_number(fstate(2,ipar)) ! growth speed seed for stochastic growth      
      end subroutine set_state_juvenile

      subroutine set_state_settle(ipar,newseed)
c     ------------------------------------------------------
c     Move larvae to bottom at settlement and immobilize it
c     ------------------------------------------------------
      integer, intent(in) :: ipar
      logical, intent(in) :: newseed
      integer             :: julday, ix,iy,ibot
c     ----------------------------------------------- 
      call get_julian_day(current_time, julday)
      mobility(:,ipar)  = 0        ! immobilize
      istate(1,ipar)    = 4        ! current type 
      istate(11,ipar)   = julday   ! day of settlement    
      ix     = nint(positions(1,ipar))
      iy     = nint(positions(2,ipar))
      ibot   = bottom_layer(ix,iy)
      positions(3,ipar) = ibot + 0.5  ! move tracer to sea bed
      if (check_habitat) istate(12,ipar) = habitat_map(ix,iy)
      if (newseed) call random_number(fstate(2,ipar)) ! growth speed seed for stochastic growth       
      end subroutine set_state_settle


      subroutine set_state_die(ipar)
c     -----------------------------------------------
      integer, intent(in) :: ipar
      integer             :: julday
c     ----------------------------------------------- 
      call get_julian_day(current_time, julday)
      mobility(:,ipar)  = 0                ! immobilize
      istate(1,ipar)    = -istate(1,ipar)  ! current type = negative of pre-death type
      istate(13,ipar)   = julday           ! day of death
      end subroutine set_state_die

      end module
