ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Generic fish larvae biology provider 
c     (Pisces silicus)
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Implements a generic fish-larvae biology provider that is oriented towards
c     configuration at the runtime side ie via the configruation files
c
c     Parameters directly configurable from control file:
c         particle_release_scheme      ! Initial distribution in the vertical
c         growth_mdl                   ! Form of the size-based growth model
c         vertical_behaviour_model     ! Active migration in the vertical direction
c         stop_conditions              ! Criteria (age, size) to stop a particle (if any)

c     Additional parameters defining these models are also specified in the control file
c
c     Parameters configurable via the emission box initialisation argument (in order)
c         initial particle size
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output
      use random_numbers

      implicit none
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
        integer :: tracerID
        integer :: sourceBox
        real    :: size             !Can be interpreted as length or weight
        real    :: tempInt          !Temperature integral, in degree-days
        real    :: lightInt         !Integrated hours of daylight,  in days
      end type
      public :: state_attributes


      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes

      interface get_property
        module procedure get_prop_state
      end interface
      public :: get_property
      public :: get_metadata_state
      public :: get_particle_version

c     =================================================
c     ========       module parameters         ========
c     =================================================
      type(clock), pointer        :: current_time
      real                        :: time_dir       !Time direction: -1 backtrack, 1 forward
      integer                     :: tracerID
      !Vertical release scheme
      logical                     :: do_release
      real                        :: release_params(2)
      !Growth model 
      logical                     :: do_growth
      real                        :: growth_params(2) 
      !Vertical behaviour model
      logical                     :: do_vertical
      integer                     :: vert_mdl       !Vertical model number
      real                        :: vert_params(2) !Parameters for vertical behaviour models
      integer                     :: swim_mdl       !Larval swimming model
      real                        :: swim_params(1)

c     --------
      contains
c     --------
     
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
      integer :: start(256), nwords
      character*256  :: strbuf, release_scheme, growth_mdl
      character*256  :: vert_mdl_name, swim_mdl_name
      real           :: time_step
c     ---------------------------------------------------------------
c     Display version
      write (*,*) trim(get_particle_version())
      tracerID = 1  ! set counter for next particle instance
    
c --- Setup particle release scheme
      if(count_tags(ctrlfile,"release_scheme")==0) then
         do_release = .FALSE. 
         write(*,*) "release_scheme = NONE"
      else
         !Read in parameters
         do_release = .TRUE.
         call read_control_data(ctrlfile,"release_scheme",strbuf)
         strbuf=adjustl(strbuf)
         call toupper(strbuf)
         write(*,*) "release_scheme = '", trim(strbuf),"'"
         !Split variables as required
         call tokenize(strbuf, start, nwords)
         read(strbuf(start(1):),*) release_scheme
         select case(release_scheme)
         case ("DEFAULT","NONE") 
           !No release - use the default relative coordinates instead
           do_release =.FALSE. 
         case ("FIXED") 
           !Particles are released at fixed depth - corresponds to a special
           !case of the RANGE scheme
           call check_nparams("release_scheme","FIXED",1,nwords-1)
           read(strbuf(start(2):),*) release_params(1)
           release_params(2) = release_params(1)
           write(*,*) "Fixed release depth = ",release_params(1)
         case ("RANGE") 
           !Particles are released uniformly between two depths
           call check_nparams("release_scheme","RANGE",2,nwords-1)
           read(strbuf(start(2):),*) release_params(1)
           read(strbuf(start(3):),*) release_params(2)
           write(*,*) "Release range = ",release_params(1:2)
         case default
           call abort_run("init_particle_start","Release scheme '"
     +        //trim(release_scheme) //"' is unknown.")
           stop 
         end select
      endif                

c --- Setup growth model 
      if(count_tags(ctrlfile,"growth_mdl")==0) then
         do_growth = .FALSE. 
         write(*,*) "growth_mdl = NONE"
      else
         !Read in parameters
         do_growth = .TRUE.
         call read_control_data(ctrlfile,"growth_mdl",strbuf)
         strbuf=adjustl(strbuf)
         call toupper(strbuf)
         write(*,*) "growth_mdl = '", trim(strbuf),"'"
         !Split variables as required
         call tokenize(strbuf, start, nwords)
         read(strbuf(start(1):),*) growth_mdl
         select case(growth_mdl)
         case ("NONE") 
           !No growth model calculated
           do_growth =.FALSE. 
         case ("TEMP") 
           !Linear growth as a function of temperature -
           call check_nparams("growth_mdl","TEMP",2,nwords-1)
           read(strbuf(start(2):),*) growth_params(1)
           read(strbuf(start(3):),*) growth_params(2)
           write(*,*) "Temp growth coefficients",growth_params(1:2)
         case default
           call abort_run("init_particle_start","Growth model '"
     +        //trim(growth_mdl) //"' is unknown.")
           stop 
         end select
      endif                

    
c --- Read vertical behaviour scheme 
      if(count_tags(ctrlfile,"vert_behaviour")==0) then
         do_vertical = .FALSE. 
         write(*,*) "vert_behaviour= NONE"
      else
         !Read in parameters
         do_vertical = .TRUE.
         call read_control_data(ctrlfile,"vert_behaviour",strbuf)
         strbuf=adjustl(strbuf)
         call toupper(strbuf)
         write(*,*) "vert_behaviour = '", trim(strbuf),"'"
         !Split variables as required
         call tokenize(strbuf, start, nwords)
         read(strbuf(start(1):),*) vert_mdl_name
         select case(vert_mdl_name)
         case ("NONE")
           !Do nothing
           do_vertical = .FALSE.
         case ("ABOVE")
           !Active velocity is zero unless depth is deeper than
           !some arbitrary depth. Swimming velocity is determined
           !by the swimming model
           vert_mdl =1
           call check_nparams("vert_behaviour","ABOVE",1,nwords-1)
           read(strbuf(start(2):),*) vert_params(1)
           write(*,*) "Larvae minimum depth = ",vert_params(1)
         case ("DVM")
           !Diel vertical migration between night and day target
           !values. Active velocity is determined by the swimming model
           vert_mdl =2
           call check_nparams("vert_behaviour","DVM",2,nwords-1)
           read(strbuf(start(2):),*) vert_params(1)
           read(strbuf(start(3):),*) vert_params(2)
           write(*,*) "DVM day/night target depths = ",vert_params(1:2)
         case ("FIXED")
           !Particle is fixed at some arbitrary input depth
           vert_mdl =3
           call check_nparams("vert_behaviour","FIXED",1,nwords-1)
           read(strbuf(start(2):),*) vert_params(1)
           write(*,*) "Larvae fixed at depth = ",vert_params(1)
           do_vertical = .FALSE.
         case ("TARGET")
           !Particle swims towards some arbitary depth. The easy realisation
           !of this is just a DVM where the target depths are the same
           vert_mdl =2
           call check_nparams("vert_behaviour","TARGET",1,nwords-1)
           read(strbuf(start(2):),*) vert_params(1)
           write(*,*) "Larvae target depth = ",vert_params(1)
           vert_params(2)=vert_params(1)
         case default
           call abort_run("init_particle_start","Vertical behaviour "
     +        " model '" //trim(vert_mdl_name) //"' is unknown.")
           stop 
         end select
      endif

c --- Read swim model - but only if there is some vertical behaviour involved
      if(do_vertical) then  !need a swim model
         if(count_tags(ctrlfile,"swim_mdl")==0)then
            call abort_run("init_particle_state","No swim_mdl is " //
     +         "specified. This should be specified to use active " //
     +         "vertical behaviour")
         else
            !Read in parameters
            call read_control_data(ctrlfile,"swim_mdl",strbuf)
            strbuf=adjustl(strbuf)
            call toupper(strbuf)
            write(*,*) "swim_mdl = '", trim(strbuf),"'"
            !Split variables as required
            call tokenize(strbuf, start, nwords)
            read(strbuf(start(1):),*) swim_mdl_name
            select case(swim_mdl_name)
            case ("FIXED")
              !Active velocity is fixed at some constant value
              swim_mdl =1
              call check_nparams("swim_mdl","FIXED",1,nwords-1)
              read(strbuf(start(2):),*) swim_params(1)
              write(*,*) "Swiming speed = ",swim_params(1)
            case default
              call abort_run("init_particle_start","Swimming model "
     +           " '" //trim(swim_mdl_name) //"' is unknown.")
              stop 
            end select
         endif
      endif


c --- Re-read time_step from configurtation file. This is not the most elegant approach
c     to getting this information but it does work :-)
      call read_control_data(ctrlfile, "time_step", time_step)
      if(time_step < 0) then
          time_dir = -1
      else
          time_dir = 1
      endif 
      write(*,*) "Time direction (-1 backwards, 1 forwards) = ",time_dir
      end subroutine init_particle_state


      character*100 function get_particle_version()  
      get_particle_version = 
     +    "Generic larvae biology provider : $Rev$"
      end function


      subroutine close_particle_state() ! module operator
      end subroutine 
      


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ---------locals ----------
      real :: geo(3), wdepth,upper_bnd,lower_bnd,rnd
      real :: init_size
      integer :: status, start(256), nwords
c     ----------------------------------------------------
      !Make particle mobile, unless we are using the fixed vertical position scheme
      if(vert_mdl==3) then !No vertical movement
         call set_tracer_mobility(space,1,1,0)  
      else  !Otherwise free
         call set_tracer_mobility_free(space)  
      endif
      !Set states
      state%sourceBox = emitboxID       ! box of origin 
      state%tracerID  = tracerID        ! unique ID number
      state%tempInt = 0.0               ! temperature integral for calculating average temp
      state%lightInt = 0.0              ! integrated time in the light

      !Take in size as an input parameter, but only if we are dealing with growth
      !otherwise set to -1
      if(do_growth) then
         call tokenize(initdata, start, nwords)
         if (nwords<1) then 
            write(*,*) "ERROR: init_state_attributes: Need to specify "
     +      // "initialisation size for sourcebox ", emitboxID,"."
            stop
         endif
         read(initdata(start(1):),*) init_size
         state%size  = init_size
      else
        state%size  = -1 
      endif

c     Increment tracerID as last    
      tracerID = tracerID + 1    !Keeps track of the number of particles released

      !Adjust vertical position according to release_scheme user-definition
      !In situations where we want to release particles between
      !e.g. 100 and 200m, but the water depth is 60m, particles
      !will be released at the bottom
      if(do_release) then
          call get_tracer_position(space, geo)
          call interpolate_wdepth(geo,wdepth,status)
          upper_bnd = min(wdepth,release_params(1))
          lower_bnd = min(wdepth,release_params(2))
          call random_number(rnd)
          geo(3) = (lower_bnd-upper_bnd)*rnd+upper_bnd
          call set_tracer_position(space,geo)
      endif 
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     ----------------------------------            
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
c     -------locals ----------      
      real :: pos(3), swim_speed,v_vert,d_target
      type(clock),pointer :: current_time
      logical :: light
c     ----------------------------------            
      if(.NOT.do_vertical) then
        v_active=0.0
        return
      endif 

      !Get larvae swimming speed
      !Currently only one swimming speed model is considered (a fixed swimming speed)
      !However, this is the point where we would add functionality to include 
      !swimming speeds given as a function of length, for example
      swim_speed = swim_params(1) 
 
      !Choose type of active behaviour
      select case(vert_mdl)
      case (1) !----- "ABOVE" scheme ------
        !Particles that fall below given depth swim upwards 
        call get_tracer_position(space, pos)
        if(pos(3)>vert_params(1)) then
            v_vert=-swim_speed         ! m/s
        else 
            v_vert =0.0
        endif
      case (2) !----- "DVM" schemes ------
        !A diel vertical migration (DVM) scheme between the day depth
        !and night depth. 
        !Get current time and position
        current_time => get_master_clock()
        call get_tracer_position(space,pos)
        !Target depth depends on time of day
        if(is_light(current_time,pos)) then 
             d_target=vert_params(1)
        else
             d_target=vert_params(2)
        endif
        !Sign of v_vert depends on where the larvae is in relation
        !to target depth
        v_vert = sign(swim_speed,d_target-pos(3))    
      case default
        !----- Unknown schemes ------
        !Unknown scheme - throw an error
        write(*,*) "ERROR: Vertical behaviour model, '",vert_mdl,
     +     "', is unknown."
        stop
      end select

      !Direction of swimming velocity vector is context specific. According to Uffe
      !we should run the vertical component forward in time, even if we are backtracking
      !The sign of the velocity therefore needs to be adjusted accordingly to ensure that
      !this is the case. We achieve this by multipling by the time_dir variable
      v_active   = (/0.0, 0.0, v_vert*time_dir/)
      end subroutine 



      subroutine update_particle_state(state, space, time_step)
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real,intent(in)                       :: time_step
      real  :: pos(3),temp,ds
      integer :: status
      type(clock),pointer :: current_time
c     -----------------
      !Get the local temperature of the particle
      call get_tracer_position(space, pos)
      call interpolate_temp(pos,temp,status)

      !Increment temperatre integral 
      state%tempInt = state%tempInt + temp*abs(time_step)/86400. ! Temperature integral in deg.days
      
      !Integrate light integral ie feeding potential
      current_time => get_master_clock()
      if(is_light(current_time,pos)) then 
        state%lightInt = state%lightInt + abs(time_step)/86400.
      endif

      !Larval growth model
      if(do_growth) then
        ds = growth_params(1)*temp + growth_params(2)  
        state%size = state%size + time_step*ds/86400.0
      endif
      
      end subroutine 


      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      end subroutine 


      subroutine get_prop_state(state,var,bucket,status)
c     ------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in) :: var
      type(polytype), intent(out) :: bucket
      integer, intent(out) :: status
c     ------------------------------------------------------------  
      status=0  !Variable is present
      select case (get_name(var))
      case ("tracerID")
        call construct(bucket,"tracerID",state%tracerID)
      case ("sourceBox")
        call construct(bucket,"sourceBox",state%sourceBox)
      case ("size")
        call construct(bucket,"size",state%size)
      case ("tempInt")
        call construct(bucket,"tempInt",state%tempInt)
      case ("lightInt")
        call construct(bucket,"lightInt",state%lightInt)
      case default
        status=1   !Cannont find variable name
      end select
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c     ------------------------------------------------------------  
      character(*), intent(in) :: var_name
      type(variable),intent(out) :: var
      integer, intent(out) :: status
c     ------------------------------------------------------------  
      status=0 !Defaults to variable found
      select case (var_name)
       case ("tracerID")
         call construct(var,"tracerID","tracer ID number",
     +     units="-",fmt="(i6)",type="int")
       case ("sourceBox")
         call construct(var,"sourceBox","source box ID",
     +     units="-",fmt="(i6)",type="int")
       case ("size")
         call construct(var,"size","particle size (length or weight)",
     +     units="user defined",fmt="(f6.2)",type="real")
       case ("tempInt")
         call construct(var,"tempInt","temperature integral",
     +     units="deg C.day",fmt="(f6.2)",type="real")
       case ("lightInt")
         call construct(var,"lightInt","daylight hours integral",
     +     units="day",fmt="(f6.2)",type="real")
      case default
        status=1  !Cannot find variable
      end select
      end subroutine get_metadata_state


      subroutine check_nparams(tag,model,requires,supplied)
c     ------------------------------------------------------------  
      !Checks that the number of parameters available is sufficent
      !to meet the requirements of the model
      character(*), intent(in) :: tag,model
      integer, intent(in)      :: requires, supplied
      character*999 :: errstr
c     ------------------------------------------------------------  
      if(requires>supplied) then
        write(errstr,123) "'" //model //"' " //
     +     tag // " requires " , requires , " parameter(s) but only ",
     +     supplied, " was/were supplied."
123     format(a, i2, a, i2, a)
        call abort_run("init_particle_state",errstr)
      endif
      end subroutine check_nparams

      end module
