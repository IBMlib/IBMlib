ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Herring biology provider 
c     (Clupea harengus silicus)
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Implements a herring-like biology provider with user configurable growth model,
c     and vertical behavioural model
c
c     Parameters directly configurable from control file:
c         particle_release_scheme      ! Initial distribution in the vertical
c         growth_mdl                   ! Form of the length-based growth model
c         vertical_behaviour_model     ! Active migration in the vertical direction
c
c     Parameters configurable via the emission box initialisation argument (in order)
c         initial particle length      [mm]
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
        integer :: tracerID
        integer :: sourceBox
        real    :: length           !Particle length in mm
        real    :: tempInt          !Temperature integral, in degree-days
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
c     ========       module parameters         ========
c     =================================================
      type(clock), pointer        :: current_time
      integer                     :: tracerID
      character*256               :: release_scheme !Larval release scheme i.e. vertical starting point distribution
      character*256               :: growth_mdl     !Larval growth model
      character*256               :: vert_mdl       !Larval vertical behaviour scheme
      real                        :: vert_params(2) !Parameters for vertical behaviour schemes

c     --------
      contains
c     --------
     
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
      integer :: start(256), nwords
c     ---------------------------------------------------------------
c     Display version
      write (*,*) "Herring biology provider : $Rev$"
      tracerID = 1  ! set counter for next particle instance
    
c --- Read growth parameters 
      call read_control_data(ctrlfile,"growth_mdl",growth_mdl)
      growth_mdl = adjustl(growth_mdl)
      call toupper(growth_mdl)
      write(*,*) "growth_mdl = '", trim(growth_mdl),"'"

c --- Read particle release scheme 
      call read_control_data(ctrlfile,
     +    "particle_release_scheme",release_scheme)
      release_scheme=adjustl(release_scheme)
      call toupper(release_scheme)
      write(*,*) "particle_release_scheme = '", trim(release_scheme),"'"
    
c --- Read vertical behaviour scheme 
      call read_control_data(ctrlfile,"vert_behaviour",vert_mdl)
      vert_mdl=adjustl(vert_mdl)
      call toupper(vert_mdl)
      write(*,*) "vertical_behaviour_model = '", trim(vert_mdl),"'"
      if (vert_mdl(1:3)== "DVM") then
            call tokenize(vert_mdl, start, nwords)
            if (nwords/=3) then 
              write(*,*) "ERROR: init_particle_state: DVM vertical "
     +          // "behaviour needs two parameters (day/night target "
     +          // "depths) e.g. DVM 10 40"
              stop
            endif
            read(vert_mdl(start(2):),*) vert_params(1)
            read(vert_mdl(start(3):),*) vert_params(2)
            vert_mdl= "DVM"
            write(*,*) "DVM day/night target depths = ",vert_params
      endif
      end subroutine init_particle_state


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
      real :: xyz(3), wdepth,upper_bnd,lower_bnd,rnd
      real :: init_len
      integer :: status, start(256), nwords
c     ----------------------------------------------------
      !Make particle mobile
      call set_tracer_mobility_free(space)  
      !Set states
      state%sourceBox = emitboxID       ! box of origin 
      state%tracerID  = tracerID        ! unique ID number
      state%tempInt = 0.0               ! temperature integral for calculating average temp

      !Length is an input parameters
      call tokenize(initdata, start, nwords)
      if (nwords/=1) then 
        write(*,*) "ERROR: init_state_attributes: Need to specify "
     +   // "initialisation length. ", nwords," initialisation variable"
     +   // "(s) were supplied."
       stop
      endif
      read(initdata(start(1):),*) init_len
      state%length  = init_len               

c     Increment tracerID as last    
      tracerID = tracerID + 1    !Keeps track of the number of particles released

      !Adjust vertical position according to release scheme
      select case (release_scheme)
      case ("UNIFORM")
          !----- "uniform" scheme ------
          !Uniform in space within relative bounds specified in input box
          !No changes required
       case ("MIK")
          !----- "MIK" scheme ------
          !Mimics the way that MIK sampling is performed - down to
          !within 5m of the bottom, or 100m, whichever is shallower
          !We assume the distribution of larvae to be uniform within
          !this region
          call get_tracer_position(space, xyz)
          call interpolate_wdepth(xyz,wdepth,status)
          upper_bnd = 0
          lower_bnd = min(wdepth-5.0,100.0)
          call random_number(rnd)
          xyz(3) = (lower_bnd-upper_bnd)*rnd+upper_bnd
          call set_tracer_position(space,xyz)
       case ("IHLS")
          !----- "IHLS" scheme ------
          !Mimics the way that IHLS sampling is performed - as defined
          !by the NauSEA group, we release uniformly between 0 and 20m, or 
          !the bottom, whichever is shallower.
           call get_tracer_position(space, xyz)
          call interpolate_wdepth(xyz,wdepth,status)
          upper_bnd = 0
          lower_bnd = min(wdepth,20.0)
          call random_number(rnd)
          xyz(3) = (lower_bnd-upper_bnd)*rnd+upper_bnd
          call set_tracer_position(space,xyz)
       case default 
          write(*,*) "ERROR: Particle release scheme, '",
     +     trim(release_scheme), "', is unknown. "
          stop
       end select
    
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     ----------------------------------            
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
c     -------locals ----------      
      real :: xyz(3),pos(3), v_vert
      type(clock),pointer :: current_time
      logical :: light
c     ----------------------------------            
      select case(vert_mdl)
      case ("PASSIVE") 
        !----- "passive" scheme ------
        !Particles have no active behaviour
      case ("ABOVE_60M")
        !----- "Above_60m" scheme ------
        !Particles that fall below 60m swim upwards at a swimming
        !speed of 5 mm/s. 
        call get_tracer_position(space, xyz)
        if(xyz(3)>60) then
            v_vert=-5.0e-3         ! m/s
        endif
      case ("DVM")
        !----- "DVM" schemes ------
        !A diel vertical migration (DVM) scheme between the day depth
        !and night depth. Particles swim between target depths at a
        !constant 5 mm/s. 
        !Get current time and position
        current_time => get_master_clock()
        call get_tracer_position(space,pos)
        if(is_light(current_time,pos)) then
             v_vert = (vert_params(1)-pos(3))/20*5e-3    
        else
             v_vert = (vert_params(2)-pos(3))/20*5e-3    
        endif
      case default
        !----- Unknown schemes ------
        !Unknown scheme - throw an error
        write(*,*) "ERROR: Vertical behaviour model, '",trim(vert_mdl),
     +     "', is unknown."
        stop
      end select
      v_active   = (/0.0, 0.0, v_vert/)

      end subroutine 



      subroutine update_particle_state(state, space, time_step)
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real,intent(in)                       :: time_step
      real  :: xyz(3),temp,dg
      integer :: status
c     -----------------
      !Get the local temperature of the particle
      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz,temp,status)

      !Increment temperatre integral 
      state%tempInt = state%tempInt + temp*abs(time_step)/86400. ! Temperature integral in deg.days

      !Larval growth model
      !T [=] deg C, dg[=]mm/day
      select case (growth_mdl)
      case ("NONE")
          ! Do nothing
      case ("OEBERST")
          ! Growth model from Oberst et al 1999
          dg = 0.011 + 0.037 * temp
      case ("FIKSEN")
          !Growth model from Fiksen and Folkvord 1999, Eq 13
          !The original model is given in units of g/g/day - we then
          !convert it to mm/day using the length-weight relationship
          !given, SL=(ln W - 2.9)/0.21. Negative growth rates are
          !ignored and set to zero
          dg = max(0.0,(-0.024+0.012*temp)/0.21)
      case ("HEATH")
          !Growth model from Heath et al 1997, Eq 13
          !The original model is given in units of relative weight
          !gain per day - we then convert it to mm/day using the same 
          !length-weight relationship from Fiksen above. Again,
          !negative growth rates are ignored and set to zero.
          dg = max(0.0,(-0.030+0.008*temp)/0.21)
      case ("FASSLER") 
          !Growth model from Fassler 2004, equation 3 
          !The original model is given in units of weight
          !gain per day - we then convert it to mm/day using the length-weight 
          !relationship given in the same paper, W=6.1807 exp(0.263 L)
          !Negative growth rates are ignored and set to zero.
          dg = max(0.0,(-0.0191+0.0089*temp)/0.263)
      case default
            write(*,*) "ERROR: Growth model '",trim(growth_mdl),
     +       "' is unknown."
            stop 
      end select   
      !Increment length accordingly
      if(growth_mdl/="NONE") then
        state%length = state%length + time_step*dg/86400.0
      endif
      
      end subroutine 


      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      end subroutine 

      end module
