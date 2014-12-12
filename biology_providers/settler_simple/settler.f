ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Biology provider template for a simple settler
c     with embedded biology model, which is dressed with settling dynamics
c     This approach assumes settling can be separated into an 
c     independent product of a spatial map and an ontogenetic window
c     This module only advances particle state up to the point'of 
c     settlement, so particle state at the end of the simulation refers to 
c     end of the simulation/time of settlement, whatever is first
c     (embedded biology state is not aware of settlement status in this protocol)
c
c     Apply name mangling at use association for biology model
c     so that biology model just have a standard particle_state interface
c
c     Simple test of setup (using provider embedded_biology_example.f): 11 Nov 2014
c     -------------------------------------------------------------------------
c
c     The module will load an arbitrary habitat set
c     If particle is inside a habitat in the settlement window
c     particle will settle at that habitat
c   
c     $Rev: 181 $
c     $LastChangedDate: 2011-01-12 00:16:47 +0100 (Wed, 12 Jan 2011) $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use polygons             ! settlement_habitats + polygon methods
      use particle_tracking    ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output  ! access polytype for get_prop_state/get_metadata_state

c     ------- particle_state interface to be dressed -------

      use particle_state_bio, only: 
     +       state_attributes_bio       => state_attributes,  
     +       init_particle_state_bio    => init_particle_state, 
     +       close_particle_state_bio   => close_particle_state, 
     +       init_state_attributes_bio  => init_state_attributes, 
     +       get_active_velocity_bio    => get_active_velocity, 
     +       write_state_attributes_bio => write_state_attributes,
     +       update_particle_state_bio  => update_particle_state_   ! note trailing underscore 

      implicit none
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private

c     --- attributes relating to settlement dynamics
      type(clock) :: release_time         
      type(clock) :: settle_time      ! only defined if type = settled
      type(state_attributes_bio) :: bio ! embedded particle_state
      real        :: survival         ! chance of survival, 0<survival<1, only updated to settlement. survival==0 means dead 
      integer     :: type       ! 0 = free, 1 = settled, 2 = dead (missed settlement)
      integer     :: sourceBox  ! spatial release box
      integer     :: settleBox  ! spatial settlement box (settlement set possibly different from release set); only defined if type = settled
      integer     :: tracerID   ! negative is unset

      end type
      public :: state_attributes


      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes

c     --- I/O + version hand shake ---
      interface get_property
        module procedure get_prop_state
      end interface
      public :: get_particle_version
      public :: get_property
      public :: get_metadata_state


c     --- extra public particle_state methods for connectivity handshake ---
      public :: inquire_settling_success
      public :: get_settlement_habitats

      
      

c     =================================================
c     ========       module data section       ========
c     =================================================
      integer, parameter    :: verbose = 0    ! control output volume
      integer               :: tracerID       ! module stamp counters

c     --- offset must be 1, so generic offset appplies for type_meaning ---
      integer, parameter    :: particle_free    = 1 
      integer, parameter    :: particle_settled = 2
      integer, parameter    :: particle_dead    = 3
      
      character(len=7), parameter :: type_meaning(3) = (/"free   ",
     +                                                   "settled",
     +                                                   "dead   "/)


      type(lonlat_polygon),allocatable,target :: settlement_habitats(:)

c     ===============================================================
                                  contains
c     ===============================================================
     
 
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      write (*,*) trim(get_particle_version())
      tracerID   = 1    ! module stamp counters
      call load_habitat_set()
      call init_particle_state_bio()
      end subroutine init_particle_state


      subroutine close_particle_state() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data
c     ---------------------------------------------------------------
      if (allocated(settlement_habitats)) 
     +    deallocate(settlement_habitats)
      call close_particle_state_bio()
      end subroutine close_particle_state
      


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     --------- local variables ---------
c     ----------------------------------------------------
c     parse initdata string (currently no options)
      call set_state_released(state, space, emitboxID)
      call init_state_attributes_bio(state%bio, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      end subroutine init_state_attributes



      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Delegate request to embedded biology
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
c     --------------------------------------------------------------
      call get_active_velocity_bio(state%bio, space, v_active)
c     --------------------------------------------------------------
      end subroutine get_active_velocity



      subroutine update_particle_state(state, space, dt)
c     --------------------------------------------------------------
c     Only invoke update_particle_state_bio up to settlement
c
c     particle state at the end of the simulation refers to 
c     end of the simulation/time of settlement/death, whatever comes first
c     (effectively freeze particles at the point of settlement/death)
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt

      real                                    :: z
      logical                                 :: die, next  
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"
      if (state%type /= particle_free) return

      call update_particle_state_bio(state%bio,space,dt,z,die,next)
      
      state%survival = state%survival * exp(-z*dt)
      if (die) then
         call set_state_dead(state, space)  
         return
      elseif (next) then
         call try_settling(state, space)
      endif
c     --------------------------------------------------------------
      end subroutine update_particle_state




      subroutine try_settling(state, space)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test whether particle is at a suitable (static) habitat. 
c
c     If a matching habitat is found, change state to particle_settled
c     Expect particle to be in free state. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real                                    :: xyz(3)
      integer                                 :: ihab
c     --------------------------------------------------------------------
      if (state%type /= particle_free) then
         write(*,*) "try_settling: type /= particle_free"
         stop 433
      endif
c     --- scan over all habitats to see --- 
      call get_tracer_position(space, xyz)
      if (is_inside_a_habitat(xyz, ihab)) then      ! picks first matching and settle there
         call set_state_settled(state, space, ihab) ! ihab only set, if is_inside_a_habitat is .true.
         return
      endif
c     --------------------------------------------------------------------
      end subroutine try_settling
  


      subroutine set_state_released(state, space, origbox)
c     -----------------------------------------------
c     When signal start the pelagic phase
c     Do not set attributes settle_time and settleBox
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      integer, intent(in)                    :: origbox
      type(clock), pointer                   :: current_time
c     -----------------------------------------------
      state%survival  = 1.0   
      call set_tracer_mobility_free(space)   ! allow advection
      current_time    => get_master_clock()   
      call set_clock(state%release_time, current_time)
      state%type      = particle_free        !
      state%sourceBox = origbox
      state%tracerID  = tracerID
      tracerID        = tracerID + 1         ! update module counter, when checking out an ID
      end subroutine set_state_released



      subroutine set_state_settled(state, space, setBox)
c     -----------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      integer, intent(in)                    :: setBox 

      type(clock), pointer                   :: current_time
c     ----------------------------------------------- 
      current_time => get_master_clock()
      call set_clock(state%settle_time, current_time)
      call set_tracer_mobility_stop(space) 
      state%type      = particle_settled        
      state%settleBox = setBox       
      end subroutine set_state_settled
      


      subroutine set_state_dead(state, space)
c     -----------------------------------------------
c     Freeze particle at point of death
c     -----------------------------------------------
      type(state_attributes),intent(inout)    :: state
      type(spatial_attributes), intent(inout) :: space
c     -----------------------------------------------
      call set_tracer_mobility_stop(space) 
      state%type      = particle_dead   
      state%survival  = 0.0 ! death signal
      end subroutine set_state_dead




      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
c     -----------------------------------------------------  
      type(state_attributes), intent(in) :: state 
c     -----------------------------------------------------   
      write(*,*)
      write(*,*) "state_attributes instance:"
      write(*,*) "particle type = ", type_meaning(state%type)
      write(*,*) "survival      = ", state%survival
      write(*,*) "release_time  = "
      call write_clock(state%release_time)
      write(*,*) "tracerID      = ", state%tracerID
      write(*,*) "sourceBox     = ", state%sourceBox
c     --- only set, if particle_settled ---
      if (state%type == particle_settled) then
         write(*,*) "settleBox     = ", state%settleBox
         write(*,*) "settle_time   = "
         call write_clock(state%settle_time)
      endif
c     --- pass request on to embedded particle_state ---
      write(*,*) "embedded particle_state = "
      call write_state_attributes_bio(state%bio)
c
      end subroutine 

c
c ================   connectivity handshake methods   ================
c

      subroutine inquire_settling_success(state, is_settled, 
     +                                    frombox, tobox, survival)
c     ----------------------------------------------------- 
c     State polling method to inquire source/setting boxes
c     of the particle
c     Descriptor tobox is only set, if is_settled is .true.
c     frombox is always set
c     survival == 0 means particle is dead
c     ----------------------------------------------------- 
      type(state_attributes), intent(in) :: state 
      logical, intent(out)               :: is_settled
      integer, intent(out)               :: frombox    
      integer, intent(out)               :: tobox   
      real,    intent(out)               :: survival   
c     ----------------------------------------------------- 
      survival = state%survival   ! always pass on
      frombox  = state%sourceBox  ! always pass on
      if (state%type == particle_settled) then
         is_settled = .true. 
         tobox      = state%settleBox
      else
         is_settled = .false.
      endif

      end subroutine       


      function get_settlement_habitats()
c     -------------------------------------------------------------------
c     Generate a reference to settlement_habitats
c
c     If memory is allocated at for this, it should be tracked so
c     it can deallocated when invoking close_particle_state() 
c     -------------------------------------------------------------------
      type(lonlat_polygon), pointer  :: get_settlement_habitats(:)  
      if (allocated(settlement_habitats)) then
         get_settlement_habitats => settlement_habitats
      else  ! no habitats loaded, return null pointer
         nullify(get_settlement_habitats)
      endif
      end function get_settlement_habitats


c
c ================    habitat interaction methods ================
c

      subroutine load_habitat_set()
c------------------------------------------------------------ 
c     Currently only scan for square settlement habitats; later allow 
c     general polygons as well (supported by polygon type)
c     scan for:
c        square_settlement_habitat = SWlon SWlat NElon NElat
c     in control file
c------------------------------------------------------------ 
      integer  :: nhabs, ihab, inext
      real     :: corners(4), nodes(2,2) 
      character(len=32) :: habname
c------------------------------------------------------------ 
      nhabs = count_tags(ctrlfile, "square_settlement_habitat")    ! first point repeated at last
      allocate ( settlement_habitats(nhabs) )
      inext = 1
      do ihab = 1, nhabs      
         call read_control_data(ctrlfile, "square_settlement_habitat", 
     +        corners, inext) ! stepping through with inext 
         inext = inext + 1
         nodes = reshape(corners, (/2,2/))
         write(habname, 650) ihab
         call init_alligned_rectangle(settlement_habitats(ihab), 
     +                                nodes, habname) 
         if (verbose>0) write(*,*) ihab, nodes(:,1), nodes(:,2) 
      enddo     
      write(*, 639) nhabs

 638  format("load_habitat_set: read settlement habitat",1x,i8,1x,
     +       ": SW = ", 2f14.6, " NE = ", 2f14.6)
 639  format("load_habitat_set: read",1x,i8,1x,"settlement habitats")   
 650  format("settlement habitat", 1x, i8)   
      end subroutine load_habitat_set



      logical function is_inside_a_habitat(xyz, ihab)
c------------------------------------------------------------ 
c     Test whether position xyz is in loaded habitat set settlement_habitats
c     and return habitat number ihab, if inside a habitat
c     Pick first matching habitat, if there is overlap (not tested)
c
c     Return .true. if inside a habitat in settlement_habitats
c     otherwise return .false. (and ihab = number_of_habitats + 1, at loop exit)
c     As a special case return .false. (and ihab=0), if settlement_habitats is empty
c------------------------------------------------------------ 
      real, intent(in)     :: xyz(3)
      integer, intent(out) :: ihab
c------------------------------------------------------------ 
c     --- handle special case
      if (.not.allocated(settlement_habitats)) then  ! special case
         ihab            = 0
         is_inside_a_habitat = .false.  
         return
      endif
c     --- scan over habitats to see if xyź is inside one, 
c         use output ihab as loop variable    
      do ihab = 1, size(settlement_habitats)
         if (is_inside_polygon(settlement_habitats(ihab), xyz)) then
            is_inside_a_habitat = .true.  
            return
         endif
      enddo
c     --- no bounding habitat was in loop, now ihab = size(settlement_habitats) + 1
      is_inside_a_habitat = .false.   
      
      end function is_inside_a_habitat




      character*100 function get_particle_version() 
c------------------------------------------------------------ 
      get_particle_version = "simple settler" //
     +     " provider : $Rev: 181 $"
      end function




c     ----- polytype interaction/query methods ---



      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c------------------------------------------------------------  
      status=0  !Variable is present
      select case (get_name(var))
      case ("type")
         call construct(bucket,"type",     state%type)
      case ("sourceBox")
         call construct(bucket,"sourceBox",state%sourceBox)             
      case ("settleBox")
         call construct(bucket,"settleBox",state%settleBox)      
      case ("tracerID")
         call construct(bucket,"tracerID", state%tracerID)     
      case default
        status=1   !Cannont find variable name
      end select
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in) :: var_name
      type(variable),intent(out) :: var
      integer, intent(out) :: status
c------------------------------------------------------------  
      status=0 !Defaults to variable found
      select case (var_name)
      case ("type")
        call construct(var,"type", 
     +                "0 = free, 1 = settled, 2 = dead",
     +                 units="-",fmt="(i1)",type="int")  
      case ("sourceBox")
        call construct(var,"sourceBox", 
     +                "sourceBox number 1... ",
     +                 units="-",fmt="(i6)",type="int")  
      case ("settleBox")
        call construct(var,"settleBox", 
     +                "settleBox number 1... ",
     +                 units="-",fmt="(i6)",type="int")    
      case ("tracerID")
        call construct(var,"tracerID", 
     +                "ID number of tracer",
     +                 units="-",fmt="(i6)",type="int")      
      case default
        status=1  !Cannot find variable
      end select
      end subroutine get_metadata_state



      end module
