ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Preliminary port of sandeel biology module
c     to Phoenix version of IBMlib
c
c     
c     NOW: deactivate backtracking (dt<0)
c          removed habitat part, because it was grid specific to cmod 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser

      implicit none
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
      real        :: size                 ! egg: maturization level  - larv/juv: length in mm
      character*3 :: type                 ! "egg"/"lar"/"juv"/"die"/"non"
      real        :: growth_seed          ! represents individual variation in growth

      integer     :: spawn_day            ! Julian day
      integer     :: hatch_day            ! Julian day
      integer     :: day_of_metamorphosis ! Julian day 

      integer     :: orig_boxID ! spatial release box
      integer     :: eggID      ! negative is unset
      integer     :: larvID     ! negative is unset
      integer     :: juvID      ! negative is unset
      
      end type
      public :: state_attributes


      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes



c     =================================================
c     ========       module data section       ========
c     =================================================
      integer, parameter    :: verbose = 0    ! control output volume

c
c     -------------- stochastic egg growth model --------------
c      
      integer :: eggID    
      real    :: hatch_begin(2)              ! days since fertilization (A,k0 > 0)
      real    :: hatch_mid(2)                ! days since fertilization (B,k1 > 0)
      real    :: hatch_end(2)                ! days since fertilization (C,k2 > 0)

c     -------------- larval growth model --------------
      integer :: larvID   ! counter to keep track of all (living/dead) emitted larvae
c
c     -------------- larvae growth data --------------
c     a(T)  = func(lambda(1) + lambda(2) T + lambda(3) T^2)     (T in Celcius)
c     dL/dt = a (L/L_0)^beta (1 - L/L_infty)
c
c     lambda       temperature polynomial coefficients  unit == mm/day/K^n 
c     func         temperature function                 0: func=polynomial  1: func=exp(polynomial)
c     L_0          hatch length                         unit == mm (should be fixed)
c     L_infty      Bertallanfy asymptote                unit == mm (fixed for now)   
c     beta         small larvae growth exponent         unit == 1 
c 
      real            :: lambda(3)
      integer         :: tempfunc
      real            :: L_0   
      real            :: beta                           
      real, parameter :: L_infty = 218.0  


c     -------------- juvenile transition data --------------
c    
c     when L > L_metamorp1 larvae dump, irrespective of its position
c    
c
      real                 :: L_metamorp1             
      integer              :: juvID  ! counter to keep track of all (living/dead) emitt


c     ===============================================================
                                  contains
c     ===============================================================
     
         
 
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      write(*,*) "init_particle_state():"
      eggID   = 1    ! module stamp counters
      larvID  = 1    ! module stamp counters
      juvID   = 1    ! module stamp counters
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
      call read_control_data(ctrlfile, "larvae_metamorph_len",
     +                       L_metamorp1)
      write(*,*) "  larvae_metamorph_len =", L_metamorp1

      end subroutine 



      subroutine close_particle_state() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data
c     ---------------------------------------------------------------
      end subroutine 
      


      subroutine init_state_attributes(state_stack, space_stack,
     +                                 time_dir,             
     +                                 nstart, npar,
     +                                 initdata,emitboxID)
c     ---------------------------------------------------- 
      type(state_attributes),intent(out)     :: state_stack(:) 
      type(spatial_attributes),intent(inout) :: space_stack(:)  
      real,intent(in)                     :: time_dir            
      integer,intent(in)                  :: nstart 
      integer,intent(in)                  :: npar
      character*(*),intent(in)            :: initdata
      integer,intent(in)                  :: emitboxID
c 
      integer   :: ip
      character :: particle_type
      real      :: initsize
      integer   :: start(256), nwords ! assume less than 256 words
c     ----------------------------------------------------
c     parse initdata string
      call tokenize(initdata, start, nwords)
      if (nwords<1) stop "init_state_attributes: initdata insufficient"
      read(initdata(start(1):),*) particle_type 
      if (nwords>1) read(initdata(start(2):),*) initsize

      do ip = nstart, nstart+npar-1
         
         call set_state_new(state_stack(ip), emitboxID)

         if (particle_type=="e") then
            call set_state_egg(state_stack(ip), space_stack(ip), .true.) 
     
         elseif (particle_type=="l") then
            if (nwords<2) stop "init_state_attributes: initsize missing"
            call set_state_larv(state_stack(ip), space_stack(ip), 
     +                          initsize, .true.)

         elseif (particle_type=="j") then
            if (nwords<2) stop "init_state_attributes: initsize missing"
            call set_state_juvenile(state_stack(ip), space_stack(ip), 
     +                              initsize, .true.)
         else
            write(*,*) "init_state_attributes:unknown particle type:"//
     +                  particle_type
            stop
         endif
         
      enddo
c     ----------------------------------------------------
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine 


      subroutine update_particle_state(state, space, dt)
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"

      if (state%type == "egg")   call egg_devel(state, space, dt)

      if ((state%type == "lar").or.
     +    (state%type == "juv")) call larvae_devel(state, space, dt)

      call tracer_metamorphosis(state, space, dt) ! all types

c     --------------------------------------------------------------
      end subroutine 


      subroutine egg_devel(state, space, dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Egg development and production
c
c     dt   : length of time step corresponding to this call (seconds)
c            negative for backtracking (not supported yet)
c
c     development in dependence of bottom temperature
c     currently, just use cell average temperature
c    
c     tau [unit==days] is the adiabatic egg development time at 
c     this temp for this egg growth seed (fstate(2))
c  
c     Assume module time reference (current_time) has been syncronized
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real, intent(in)                      :: dt
      integer                               :: status
      real                                  :: tm,a,b,c,g,tau,xyz(3)
c     --------------------------------------------------------------
      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz,tm,status)
      
      a      = hatch_begin(1)*exp(-hatch_begin(2)*tm)
      b      = hatch_end(1)  *exp(-hatch_end(2)*tm)
      c      = hatch_mid(1)  *exp(-hatch_mid(2)*tm)
      g      = (a*b - c**2)/(a + b - 2.*c)  
      tau    = ((b-g)/(a-g))**state%growth_seed
      tau    = g + (a-g)*tau
      state%size = state%size + dt/86400./tau
      
      end subroutine egg_devel
  

      subroutine larvae_devel(state, space, dt) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Temperature dependent larval/juvenile growth
c
c     dt   : length of time step corresponding to this call (seconds)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real, intent(in)                      :: dt
      integer                               :: status
      real                                  :: tm, a, ldot, l,xyz(3)

c     --------------------------------------------------------------
      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz,tm,status)
     
      a    = lambda(1) + lambda(2)*tm + lambda(3)*tm**2
      if      (tempfunc==0) then
         a    = max(a, 0.0)     ! must be positive
      elseif  (tempfunc==1) then
         a    = exp(a)          ! always positive
      else
         stop "unknown tempfunc"
      endif

      l    = state%size
      ldot = a * (l/L_0)**beta * (1. - l/L_infty) ! unit == mm/day
      state%size = state%size + ldot*dt/86400. ! Euler increment
      
      end subroutine larvae_devel




      subroutine tracer_metamorphosis(state, space, dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Promote individuals to new types, if development 
c     thresholds has been fulfilled
c     
c     Check all particles consequetively
c     egg      -> larvae
c     larvae   -> juvenile
c     (juvenile -> settled juvenile -> death) currently disabled because habitat quality
c      map is not available for this module context    
c 
c     if (L > L_metamorp1):  stop, irrespective value of habitat at 
c                            current position
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real, intent(in)                        :: dt
      integer          :: ipar, julday, ix, iy, last
      integer          :: newlarv, newjuv, newsett, newdied
      logical          :: settle, die, consider_settling, must_settle 
c     --------------------------------------------------------------------
      
      if ((state%type == "die").or.(state%type == "non")) return   ! do nothing
      
      newlarv = 0
      newjuv  = 0
c     -----------------------------
c     egg -> larvae transition
c     -----------------------------
      
      if ((state%type == "egg").and.(state%size > 1.0d0)) then
         call set_state_larv(state, space, L_0, .false.)
         newlarv = newlarv + 1
         return ! no further processing for this particle
      endif


c     -----------------------------
c     larvae -> juvenile transition 
c     -----------------------------
      if ((state%type == "lar").and.(state%size > L_metamorp1)) then                
         call set_state_juvenile(state, space, state%size, .false.) ! keep size
         newjuv = newjuv + 1
      endif

c     -----------------------------
c     settling/dying transition 
c     -----------------------------
c     ... disabled       

 
c --- report transformations
      if (newlarv > 0) write(*,422) "egg    -> larvae", newlarv
      if (newjuv  > 0) write(*,422) "larvae -> juvenile", newjuv
      

 422  format("tracer_metamorphosis:",1x,a,":",1x,i5)


      end subroutine tracer_metamorphosis



      subroutine set_state_new(state,origbox)
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      integer, intent(in)                    :: origbox
c     -----------------------------------------------
      state%size = 0.0          ! 
      state%type = "non"        ! no declared
      state%growth_seed = 0.0   ! represents individual variation in growth

      state%spawn_day            = -1  ! Julian day (negative == unset)
      state%hatch_day            = -1  ! Julian day (negative == unset)
      state%day_of_metamorphosis = -1  ! Julian day (negative == unset) 

      state%orig_boxID = origbox  ! spatial release box
      state%eggID      = -1     ! negative is unset
      state%larvID     = -1     ! negative is unset
      state%juvID      = -1     ! negative is unset

      end subroutine set_state_new



      subroutine set_state_egg(state,space,newseed)
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      logical, intent(in)                    :: newseed
      integer                                :: julday
      type(clock), pointer                   :: current_time
c     ----------------------------------------------- 
      current_time => get_master_clock()
      call get_julian_day(current_time, julday)
      call set_tracer_mobility_stop(space)
      state%type      = "egg"
      state%spawn_day = julday
      state%size      = 0.0       
      state%eggID     = eggID
      eggID           = eggID + 1
      if (newseed) call random_number(state%growth_seed) ! growth speed seed for stochastic growth      
      end subroutine set_state_egg



      subroutine set_state_larv(state,space,size,newseed)
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      real, intent(in)                       :: size
      logical, intent(in)                    :: newseed
      integer                                :: julday
      type(clock), pointer                   :: current_time
c     ----------------------------------------------- 
      current_time => get_master_clock()
      call get_julian_day(current_time, julday)
      call set_tracer_mobility_free(space)
      state%type      = "lar"
      state%hatch_day = julday
      state%size      = size    
      state%larvID    = larvID
      larvID          = larvID + 1
      if (newseed) call random_number(state%growth_seed) ! growth speed seed for stochastic growth      
           
      end subroutine set_state_larv



      subroutine set_state_juvenile(state,space,size,newseed)
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      real, intent(in)                       :: size
      logical, intent(in)                    :: newseed
      integer                                :: julday  
      type(clock), pointer                   :: current_time
c     ----------------------------------------------- 
      current_time => get_master_clock()
      call get_julian_day(current_time, julday)
      call set_tracer_mobility_stop(space)
      state%type                 = "juv"
      state%day_of_metamorphosis = julday
      state%size                 = size
      state%juvID                = juvID
      juvID = juvID + 1
      if (newseed) call random_number(state%growth_seed) ! growth speed seed for stochastic growth         
      end subroutine set_state_juvenile
      


      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      write(*,*) "type ", state%type, " size = ", state%size, " mm"
      end subroutine 


      end module
