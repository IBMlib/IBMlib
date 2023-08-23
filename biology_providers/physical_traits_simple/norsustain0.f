ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Model 0 for NORSUSTAIN
c     derived from traits.f

c     * simple multi stage template
c     * allow settlement in geographical polygons defined in input
c     * backtracking enabled
c     -------------------------------------------------------------------------    
c   
c     $Rev:  $
c     $LastChangedDate:  $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state ! NB, different module name

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output
      use polygons
      
      implicit none
      private    
      
c     ---------------------------------------------
c     This type contains only specific biological attributes
c     ---------------------------------------------
      type state_attributes
c      private
         real    :: age         ! total since instantiation [days]
         integer :: source      ! emit box number
         integer :: stage       ! see below
         real    :: eggdevel    ! 0 < egg devel stage < 1 (if egg)
         real    :: length      ! current larval length [mm] (if larvae)
         real    :: age_juv_trans ! age at juvenile transition [days] where search for habitat begins
         integer :: settlement_habitat ! which habitat > 0 larvae settled in
      end type
      public :: state_attributes

c     ---- export generic particle_state interface ----
      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: delete_state_attributes
      public :: get_active_velocity
      public :: write_state_attributes
      public :: update_particle_state ! basic argument list

      public :: get_particle_version 
      public :: get_property
      public :: get_metadata_state
      interface get_property
        module procedure get_prop_state
      end interface

c     =================================================
c     ========       module data section       ========
c     =================================================
c
c     --------- egg stage  ---------
c      
      integer    :: egg_buoyancy         ! egg buoyancy model, default = 0
      real       :: sink_speed_egg       ! [m/s],    applies for buoyancy=0, positive down
      real       :: rho_egg     ! [kg/m**3] fixed particle density, applies for egg buoyancy=1
      real       :: equivalent_salinity  ! [PSU]  - see Sundby 2015
      real       :: egg_radius           ! [m]       for Stokes law, applies for egg buoyancy=1
      real       :: egg_degdays          ! Degree-days to complete egg development [degC*days]
      real       :: egg_degdays_minT ! lower temperature limit for degree day model [degC]
      
c
c     --------- epipelagic larval stage  ---------
c
      real, parameter  :: larvlength_0   = 12.0   ! first-feeding length [mm] ref Stenberg 2016
      real, parameter  :: larvlength_juv = 65.0   ! settle length   [mm] ref Stenberg 2016
      real             :: LenGrate                ! length growth rate [mm/day]      
      real             :: dLenGrate_dT            ! length growth rate T-slope [mm/day/degC]
      real             :: depthmax                ! maximum depth in epipelagic phase [m]
c
c     --------- juvenile stage  ---------
c     
      real             :: settle_period_length  ! max time since reaching settlement threshold [days] - otherwise dead
c
c     --------- stage handles  ---------
c
c     dead stages  < 0
c     alive stages >= 0
      integer, parameter :: is_unhatched       = -2 ! for backtracking
      integer, parameter :: is_dead            = -1 ! undisclosed cause, includes unborn
      
      integer, parameter :: is_egg             =  0
      integer, parameter :: is_pelagic_larvae  =  1
      integer, parameter :: is_juvenile_larvae =  2 ! searching for settlement
      integer, parameter :: is_settled_larvae  =  3
      integer, parameter :: is_unresolved      =  100 ! just track age
c
c     --------- settlement areas  ---------
c
      type(lonlat_polygon), pointer :: settlement_areas(:)   ! nullified, if no settlement areas
      
  
c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
c     ------------------------------------------------------------------------
c     ------------------------------------------------------------------------
      integer           :: ise, nse, ihit, start(5000), nwords ! buffer start must be sufficient 
      real, allocatable :: rbuf(:), nodes(:,:)
      character*5000    :: strbuf   ! should contain fairly large polygons
      character*256     :: name
c
c     --------- resolve egg parameters ---------
c
      write(*,*) "--- egg model ---"

      if (count_tags(ctrlfile, "sink_speed_egg")>0) then ! precedence, if provided

         call read_control_data(ctrlfile, "sink_speed_egg",
     +                          sink_speed_egg) 
         egg_buoyancy = 0
         write(*,390)  sink_speed_egg

      elseif (count_tags(ctrlfile, "rho_egg")>0) then 

         call read_control_data(ctrlfile,"rho_egg",rho_egg) 
         call read_control_data(ctrlfile,"egg_radius", egg_radius)  ! must be present, if rho_particle specified, no default
         egg_buoyancy = 1
         write(*,391) rho_egg, egg_radius

      elseif (count_tags(ctrlfile, "equivalent_salinity")>0) then

         call read_control_data(ctrlfile,"equivalent_salinity",
     +                          equivalent_salinity) 
         call read_control_data(ctrlfile,"egg_radius", egg_radius)  ! must be present, if rho_particle specified, no default
         egg_buoyancy = 2
         write(*,392) equivalent_salinity, egg_radius
         
      else ! default passive particle

         sink_speed_egg = 0.0
         egg_buoyancy   = 0
         write(*,399) 

      endif
c     --- egg development model
    
      call read_control_data(ctrlfile,"egg_degdays", egg_degdays)
      write(*,400)  egg_degdays
      call read_control_data(ctrlfile,"egg_degdays_minT",
     +     egg_degdays_minT)
      write(*,401)  egg_degdays_minT
      
 388  format(a20, " = ", f8.3, 1x, a)
 390  format("eggs have constant additive sinking speed = ",
     +       f8.3, " m/s")
 391  format("eggs have constant density = ", 
     +     f8.3, " kg/m**3 and radius = ", e12.3, " m")
 392  format("eggs have equiv. salinity = ", 
     +     f8.3, " kg/m**3 and radius = ", e12.3, " m")     
 399  format("eggs are neutrally buoyant")
      
 400  format("Degree-days to complete egg development = ",f8.3,
     +     " degC*days")
 401  format("Lower temperature limit for degree day model =",
     +      f8.3," degC")      
c
c     --------- resolve larval parameters ---------
c
      write(*,*) "--- larval model ---"
      call read_control_data(ctrlfile,"LenGrate", LenGrate)
      write(*,410)  LenGrate
      call read_control_data(ctrlfile,"dLenGrate_dT", dLenGrate_dT)
      write(*,411)  dLenGrate_dT
      call read_control_data(ctrlfile,"depthmax", depthmax)
      write(*,412)  depthmax
     
 410  format("Larval length growth rate = ",f8.3, " mm/day")
 411  format("Larval length growth rate temp derivative = ", f8.3,
     + " mm/day/degC")   
 412  format("Larval epipelagic depth maximum = ",f8.3, " m") 
c
c     --------- resolve juvenile parameters ---------
c     
      call read_control_data(ctrlfile,"settle_period_length",
     +     settle_period_length)
      write(*,415)  settle_period_length
 415  format("Larval max search time since juv maturation = ",
     +    f8.3, " days") 
c
c     --------- resolve settlement areas ---------
c
      nse = count_tags(ctrlfile, "settlement_area")
      if (nse == 0) then
         write(*,*) "No settlement areas defined in input file"
         nullify(settlement_areas) ! signal no areas
      else ! scan set of settlement areas
         allocate( settlement_areas(nse) )
         ihit = 1
         do ise = 1, nse
            call read_control_data(ctrlfile,"settlement_area",
     +           strbuf,ihit)
            call tokenize(strbuf, start, nwords)
            if (mod(nwords,2)==1) then
               write(*,*) "settlement area ", ise,
     +        "contain an odd number of coordinates"
               stop 492
            else
               allocate( rbuf(nwords)      )
               allocate( nodes(2,nwords/2) )
               read(strbuf,*) rbuf
               nodes = reshape(rbuf, (/2,nwords/2/))
               write(name, 371) 'settlement_area_', ise
               call init_lonlat_polygon(settlement_areas(ise),
     +              nodes, name)
               deallocate( rbuf )
               deallocate( nodes )
               write(*, 372) ise, nwords/2
            endif
            ihit = ihit + 1 ! start reading from next line
         enddo                  ! do ise =
         write(*,*) "read",nse,"settlement areas"
      endif    ! if (nse == 0)
         
 371  format(a,i6.6)
 372  format('read settlement area ',i,' with ',i,' nodes')
      
      end subroutine init_particle_state

      

      subroutine close_particle_state()  ! module operator
      end subroutine close_particle_state



      subroutine init_state_attributes(state, space, time_dir,             
     +     initdata, emitboxID)
c     ---------------------------------------------------------
c     time_dir > 0: (forward) default: newly spawned egg, 
c     otherwise provide initdata:
c                e eggdevel
c                l length
c                a age
c     time_dir < 0: (backtrack) no default - initdata must be supplied
c     ---------------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
      character         :: tag
      real              :: value
c     ---------------------------------------------------------
      state%source             = emitboxID
      state%age_juv_trans      = -1.0
      state%settlement_habitat = -1
      state%age                = 0.0
      state%stage              = is_egg
      state%eggdevel           = 0.0
      state%length             = -1.0
c
c     same parsing of initdata for forwad/backward, i.e. no split by time_dir
c
      if (trim(adjustl(initdata)) /= "") then ! overwrite default above
         
         read(initdata,*) tag,value
       
         if (tag=='e') then
            state%eggdevel  = value
         elseif (tag=='l') then
            state%stage  = is_pelagic_larvae
            state%length = value
         elseif (tag=='a') then
            state%stage = is_unresolved
            state%age   = value
         else
            write(*,*) "init_state_attributes: invalid initdata = ",
     +                 initdata
            stop 744
         endif
         
      elseif (time_dir < 0) then
         
         write(*,*) "init_state_attributes: initdata "//
     +        "must be supplied for backtracking"
         stop 293
         
      endif                     ! initdata != ""
!     
      call set_shore_BC(space, BC_reflect) 

      end subroutine init_state_attributes


      
      subroutine delete_state_attributes(state)
c     ----------------------------------------------------
      type(state_attributes), intent(inout) :: state 
      end subroutine 

      

      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Implement simple vertical velocities
c     Currently same for forward/backward dynamics
c     egg_buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity 
c     egg_buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density
c     egg_buoyancy = 2: additive terminal Stokes velocity, corresponding to fixed equivalent salinity      
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)

      real, parameter                      :: patm       = 1.01325    ! [bar] applied constant atmospheric pressure
      real, parameter                      :: rhow0      = 1027.0     ! [kg/m**3], applied constant water density for compressibility effects
      real, parameter                      :: g_x_Pa2bar = 9.81/1.e5 ! g * conversion from Pa to bar
      real, parameter                      :: g          = 9.81       ! [m/s2] gravitational constant
      real, parameter                      :: mydyn      = 1.002e-3   ! [kg/m/s] dynamic viscosity of water (at 20 oC)
      real                                 :: pressure, rho, salt, temp
      real                                 :: xyz(3), drho, rho2
      integer                              :: istat
      real, external                       :: rho_UNESCO ! defiend in ooceanography_providers/generic_elements/water_density.f
c     --------------------------------------------------------------
      v_active = 0             ! default passive particle

      if (state%stage == is_egg) then
         
         if      (egg_buoyancy == 0) then ! fixed sinking speed,
            
            v_active(3) = sink_speed_egg  ! positive down

         elseif  (egg_buoyancy == 1) then ! fixed particle density

            call get_tracer_position(space, xyz)
            call interpolate_temp(xyz,   temp,  istat)
            call interpolate_salty(xyz,  salt,  istat)
            pressure    = patm + g_x_Pa2bar*xyz(3)*rhow0 ! zeroth order approx., unit == bar
            rho         = rho_UNESCO(salt,temp,pressure) ! unit == kg/m**3
            drho        = rho_egg-rho
            v_active(3) = 2*drho*g*egg_radius**2/9/mydyn ! [m/s], Stokes terminal velocity law  (positive down)

         elseif  (egg_buoyancy == 2) then ! behave as water parcel with fixed equivalent salinity

            call get_tracer_position(space, xyz)
            call interpolate_temp(xyz,   temp,  istat)
            call interpolate_salty(xyz,  salt,  istat)
            pressure    = patm + g_x_Pa2bar*xyz(3)*rhow0 ! zeroth order approx., unit == bar
            rho         = rho_UNESCO(salt,temp,pressure) ! ambient water density {kg/m**3]
            rho2        = rho_UNESCO(equivalent_salinity,temp,pressure) !
            drho        = rho2-rho
            v_active(3) = 2*drho*g*egg_radius**2/9/mydyn ! [m/s], Stokes terminal velocity law  (positive down)
            
         else

            write(*,*) "get_active_velocity: unsupported buoyancy=",
     +              egg_buoyancy
            stop 653
         
         endif
      
      endif ! currrently no horizontal/vertical swimming for larval stages

      end subroutine get_active_velocity



      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
c     --------------------------------------------------------------
      write(*,*) "age since init = ", state%age  
c
      end subroutine write_state_attributes


      
      subroutine update_particle_state(state, space, dt)
c     --------------------------------------------------------------
c     forward/backward fanout
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds      
      if (dt>0) then
         call update_particle_state_fwd(state, space, dt) 
      else
         call update_particle_state_back(state, space, dt)
      endif
      end subroutine update_particle_state


      
      
      subroutine update_particle_state_fwd(state, space, dt)
c     --------------------------------------------------------------
c     This subroutine should integrate state FORWARD for time period dt 
c     Settlement conditions: 
c        1) settle_period_start < age < settle_period_end => set_tracer_mobility_stop
c        2) Settlement is a Poisson process with sink_rate
c     At settlement, place tracer at seabed
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real                                    :: xyz(3)
      integer                                 :: status, iset
      real                                    :: rnum, temp, cwd, t_x_t
      real                                    :: temp_clipped, dLen
      real                                    :: tsearch
      logical                                 :: inhab
c     --------------------------------------------------------------
      if (state%stage < 0) return ! do not update age, so age correspond to age of dying
      
      state%age      = state%age + dt/86400.
      call get_tracer_position(space,xyz) ! needed in all cases below
c
      select case (state%stage)
         
         case (is_egg)             ! update/promote egg
         
            call interpolate_temp(xyz,temp,status)
            temp_clipped   = max(temp,egg_degdays_minT)
            t_x_t            = (dt/86400)*temp_clipped
            state%eggdevel = state%eggdevel + t_x_t/egg_degdays
            if (state%eggdevel > 1) then ! hatch, keep eggdevel
               state%stage  = is_pelagic_larvae
               state%length = larvlength_0
            endif
            
         case (is_pelagic_larvae) ! update/promote pelagic larvae
            
            call interpolate_temp(xyz,temp,status)
            dLen = (LenGrate + dLenGrate_dT*temp)*dt/86400 ! may become negative
            dLen = max(0.0, dLen) ! never shrink
            state%length = state%length + dLen
            if (state%length > larvlength_juv) then ! hatch, keep length
               state%stage         = is_juvenile_larvae
               state%age_juv_trans = state%age ! record transition age
            endif
c           ----  check whether depthmax has been exceeded ----         
            if (xyz(3)>depthmax) then
               call interpolate_wdepth(xyz,cwd,status) ! avoid putting larvae below seabed
               xyz(3) = min(cwd, depthmax)
               call set_tracer_position(space,xyz)     ! move back to epipelagic layer
            endif
            
         case (is_juvenile_larvae) ! update/promote pelagic larvae
            
            tsearch = state%age - state%age_juv_trans
            if (tsearch > settle_period_length) then
               call set_state_dead(state, space) 
            else
               call check_if_settlement_habitat(xyz, inhab, iset)
               if (inhab) then
                  state%stage              = is_settled_larvae
                  state%settlement_habitat = iset
                  call interpolate_wdepth(xyz,cwd,status)
                  xyz(3) = cwd ! palce at bottum
                  call set_tracer_position(space,xyz)
                  call set_tracer_mobility_stop(space)
               endif
            endif
            
         case default           ! dead/unresolved
     
      end select    
c     ----------------------------------------
      end subroutine update_particle_state_fwd



      
      subroutine update_particle_state_back(state, space, dt)
c     --------------------------------------------------------------
c     This subroutine should integrate state BACKWARD for time period dt 
c     NB: dt < 0
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real                                    :: xyz(3)
      integer                                 :: status, iset
      real                                    :: rnum, temp, cwd, t_x_t
      real                                    :: temp_clipped, dLen
      real                                    :: tsearch
      logical                                 :: inhab
c     --------------------------------------------------------------
      if (state%stage < 0) return ! do not update age, so age correspond to age of dying
      
      state%age      = state%age + dt/86400. ! dt<0
      call get_tracer_position(space,xyz) ! needed in all cases below
c
      select case (state%stage)
         
         case (is_egg)             ! update/promote egg
         
            call interpolate_temp(xyz,temp,status)
            temp_clipped   = max(temp,egg_degdays_minT)
            t_x_t            = (dt/86400)*temp_clipped
            state%eggdevel = state%eggdevel + t_x_t/egg_degdays
            if (state%eggdevel < 1) then ! hatch point, freeze
               state%stage  = is_unhatched
               call set_tracer_mobility_stop(space)
            endif
            
         case (is_pelagic_larvae) ! update/promote pelagic larvae
            
            call interpolate_temp(xyz,temp,status)
            dLen = (LenGrate + dLenGrate_dT*temp)*dt/86400 ! may become positive
            dLen = min(0.0, dLen) ! never grow (= shrink backward)
            state%length = state%length + dLen
            if (state%length < larvlength_0) then ! go back to egg
               state%stage     = is_egg
               state%eggdevel  = 1.0
            endif
c           ----  check whether depthmax has been exceeded ----         
            if (xyz(3)>depthmax) then
               call interpolate_wdepth(xyz,cwd,status) ! avoid putting larvae below seabed
               xyz(3) = min(cwd, depthmax)
               call set_tracer_position(space,xyz)     ! move back to epipelagic layer
            endif
            
         case default           ! if dead Elvis already left the building
     
      end select    
c     ----------------------------------------
      end subroutine update_particle_state_back


      
      character*100 function get_particle_version()  
      end function


      
      subroutine check_if_settlement_habitat(geo, inside_any, iset)
c------------------------------------------------------------
c     Scan whether geo is inside settlement habitats
c     If settlement habitat set is empty, return inside_any = .false.
c     If inside any settlement habitat, return inside_any = .true.
c     and the first set number, where it is inside, otherwise
c     return inside_any = .false. and iset undefined
c------------------------------------------------------------
      real, intent(in)     :: geo(:)
      logical, intent(out) :: inside_any
      integer, intent(out) :: iset
c------------------------------------------------------------
      inside_any = .false.
      if (.not.associated(settlement_areas) ) then
         iset = -1
         return
      endif
c     --- if membership not found in loop then inside_any remains .false.  ---
      do iset = 1, size(settlement_areas)
         if (is_inside_polygon(settlement_areas(iset), geo)) then
            inside_any = .true.
            return
         endif   
      enddo
      
      end subroutine check_if_settlement_habitat


      
      subroutine set_state_dead(state, space)
c------------------------------------------------------------
c     Undisclosed death      
c------------------------------------------------------------  
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      state%stage              = is_dead
      call set_tracer_mobility_stop(space)
      end subroutine set_state_dead




      
      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c------------------------------------------------------------  
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in)   :: var_name
      type(variable),intent(out) :: var
      integer, intent(out)       :: status
c------------------------------------------------------------  
      end subroutine get_metadata_state


      
      
      end module
