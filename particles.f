ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Manage an emsemble of particles 
c     
c     This module combines spatial and biological
c     sub classes, and redelegates overall calls to particles 
c     to spatial and biological sub classes. 
c      
c     This module exports two derived types:
c     
c       particle (combining spatial and biological sub classes)
c
c       particle_ensemble (efficient management of an array of particles)
c
c     as well as methods for initialization/update/manipulation/destruction
c     of instances of these derived types
c
c     Things we want to develop over the next years (and have data structures prepared for):
c     
c     1) particles are not a singleton stack anymore - now multiple particle stacks
c        are feasible (multi species, several ensembles)
c
c     2) parallelization over particles 
c
c     3) possible a dynamic particle stack (able to grow/shrink along the simulation)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module particles
      use particle_tracking
      use particle_state
c      use time_tools

      implicit none
      private          ! default visibility

c --------------------------------------------------------
c ---------------- module data section -------------------
c --------------------------------------------------------


c     -----------------------------------------------    
c     particle components are declared as pointer to
c     Allow spectator construction in particle_ensemble
c     -----------------------------------------------
      type particle
      private
        type(spatial_attributes), pointer :: space 
        type(state_attributes), pointer   :: state 
      end type
      public :: particle

c     -----------------------------------------------    
c     The particle_ensemble type internally has 
c     flat stacks for space/state parts of particles
c     The particles in allpart are spectator, 
c     i.e. link to space_stack/state_stack
c     The spectator construction is rather economic memory wise
c     -----------------------------------------------
      type particle_ensemble
      private
        type(particle), pointer          :: allpart(:) ! allocation may exceed active particles
        integer                          :: last       ! last active, last==0 for no active
        type(spatial_attributes),pointer :: space_stack(:) ! space data core
        type(state_attributes),pointer   :: state_stack(:) ! state data core
      end type
      public :: particle_ensemble


c     ---------------- define module public interface ------------------- 
     
      

c     .............................................................     
c     ................. module  methods .................
c     .............................................................

      public :: init_particles  ! module start up
      public :: close_particles ! module close down

c     .............................................................     
c     ................. particle  methods .................
c     .............................................................

      public :: write_particle

c     .............................................................     
c     ................. particle_ensemble methods .................
c     .............................................................

      
      public :: setup_ensemble
      public :: setup_ensemble_from_file
      public :: delete_ensemble

      interface generate_particles
         module procedure generate_particles_single
         module procedure generate_particles_vector      
      end interface
      public :: generate_particles  ! using an emission_box (overloaded subroutine)
        
      public :: get_last_particle_number
      public :: get_particle        ! Checkout a particle
      public :: get_particle_position ! delegate to particle_tracking

      public :: insert_particle     ! no copy
      public :: delete_particle
      public :: update_particles ! spatial/biological, no emit - incl respond to BC


      public :: write_ensemble 

c     .................................................   
c     .................   reexports   .................
c     .................................................

      public :: emission_box       ! reexport import from particle_tracking.f
      public :: write_emission_box ! reexport import from particle_tracking.f

c     .............................................................     
c     ................. module data section .................
c     .............................................................

      integer, parameter :: verbose = 0

c...............................................................
      contains
c...............................................................
    
      subroutine init_particles()
c     ---------------------------------------------------
c     Module start up
c     ---------------------------------------------------
      call init_particle_tracking()
      call init_particle_state()

      end subroutine init_particles


      subroutine close_particles()
c     ---------------------------------------------------
c     Module close down
c     ---------------------------------------------------
      call close_particle_state()
      call close_particle_tracking()
    
      end subroutine close_particles
     
c     .............................................................     
c     ................. particle  methods .................
c     .............................................................


      subroutine write_particle(par)
c     ---------------------------------------------------
      type(particle),intent(in) :: par
      call write_spatial_attributes(par%space)
      call write_state_attributes(par%state)
      end subroutine write_particle



c     .............................................................     
c     ................. particle_ensemble methods .................
c     .............................................................


      subroutine setup_ensemble(par_ens, npar) 
c     ---------------------------------------------------
c     Generic initialization of par_ens
c     Note that this subroutine may leak memory, if
c     par_ens is already allocated (there is no portable
c     query facilities in F90 to avoid this)
c     ---------------------------------------------------
      type(particle_ensemble), intent(out):: par_ens
      integer, intent(in)                 :: npar
      integer                             :: i
c     ---------------------------------------------------
      if (verbose>0) then
         write(*,*) "setup_particle_stack: allocating particle "
     +           // "ensemble with size", npar
      endif
      allocate( par_ens%space_stack(npar) ) 
      allocate( par_ens%state_stack(npar) ) 
      allocate( par_ens%allpart(npar)     )
      par_ens%last = 0 ! no active, initially
c
c.....setup allpart as spectator to space_stack/state_stack 
c   
      do i=1,npar
         par_ens%allpart(i)%space => par_ens%space_stack(i)
         par_ens%allpart(i)%state => par_ens%state_stack(i)
      enddo

      end subroutine setup_ensemble
      

      
      subroutine setup_ensemble_from_file(par_ens, tag, boxrefs) 
c     ---------------------------------------------------
c     Make it easy to setup an ensemble par_ens, if it is
c     associated with a single set of emission boxes
c     (identified by tag in input file). tag may 
c     appear multiple times; one emission box is created
c     for each occurence by create_emission_boxes in 
c     emission box array boxrefs (allocated in create_emission_boxes).
c     create_emission_boxes returns a nullified pointer 
c     boxref, if emission boxes are created
c     ---------------------------------------------------
      type(particle_ensemble), intent(out):: par_ens
      character*(*),intent(in)            :: tag
      type(emission_box), pointer         :: boxrefs(:)
c
      integer  :: nmaxpar
c     ---------------------------------------------------      
      call create_emission_boxes(tag, boxrefs, nmaxpar)
      call setup_ensemble(par_ens, nmaxpar)

      end subroutine setup_ensemble_from_file


      subroutine delete_ensemble(par_ens)  
c     ---------------------------------------------------
c     Assume integrity of par_ens is OK and
c     that par_ens is initialized
c     Delete reversely (from end of list) to avoid 
c     copying overhead
c     ---------------------------------------------------
      type(particle_ensemble) :: par_ens
c     ---------------------------------------------------
      do while (par_ens%last > 0) 
         call delete_particle(par_ens, par_ens%last) 
      end do

      end subroutine delete_ensemble
      

    

      subroutine generate_particles_single(par_ens, emit_box, time_dir)  
c     ---------------------------------------------------
c     Generate particles to particle ensemble par_ens from emission box emit_box
c     
c     Start by generating spatial attributes from activate_emission_box
c     directly to space_stack, then delegate initialization of 
c     corresponding state attributes of particles to init_state_attributes
c     ---------------------------------------------------
      type(particle_ensemble),intent(inout) :: par_ens
      type(emission_box),intent(inout)      :: emit_box
      real,intent(in)                       :: time_dir ! >0: forward, <0 backward
c     ...................................................
      integer                               :: npar,ip,boxID
      integer                               :: nstart
      character*999                         :: strdat
c     ---------------------------------------------------
c
c     activate_emission_box chacks that stack capacity is
c     not exceeded by new emissions
c
      nstart = par_ens%last+1 ! position where to start ghenerating new particles
      call activate_emission_box(emit_box, 
     +                           time_dir, par_ens%space_stack,
     +                           par_ens%last+1, npar)

      if (npar == 0) then    ! zero particles generated
         if (verbose>0) then
            write(*,*) "generate_particles: no new particles"
         endif
         return   ! need no further processing, last is up to date   
      endif
      
c
c     --- init bio states and assemble particle from space + state
c         as a blocked call (handle all npar particles in one short
c         to minimize parsing overhead of initialization information
c         emit_box%initialization_data)
c
      call get_initialization_data(emit_box, strdat)
      call get_emission_boxID(emit_box, boxID)  ! enable state to trace its origin
      call init_state_attributes(par_ens%state_stack, 
     +                           par_ens%space_stack,
     +                           time_dir,             
     +                           par_ens%last+1, npar, strdat, boxID) 

      par_ens%last =  par_ens%last + npar

      if (verbose>0) then
        write(*,*) "generate_particles: ", npar, "new particles"
      endif

      end subroutine generate_particles_single




      subroutine generate_particles_vector(par_ens,emit_boxes,time_dir)  
c     ---------------------------------------------------
c     Vector wrapper to generate_particles_single
c     generate_particles_vector accepts and array of emission boxes
c     whereas generate_particles_single accepts a single emission_box
c     ---------------------------------------------------
      type(particle_ensemble),intent(inout) :: par_ens
      type(emission_box),intent(inout)      :: emit_boxes(:)
      real,intent(in)                       :: time_dir ! >0: forward, <0 backward
c     ...................................................
      integer                               :: ibox
c     ---------------------------------------------------
      do ibox = 1, size(emit_boxes)
         call generate_particles_single(par_ens, emit_boxes(ibox), 
     +                                  time_dir)  
      enddo

      end subroutine generate_particles_vector





      subroutine set_last_particle(par_ens, last)
c     ---------------------------------------------------   
c     This subroutine should be private, due to
c     risk of leaking memory (particles after last
c     are not deleted by this subroutine by decreasing last)
c     Assume integrity of par_ens is OK
c     ---------------------------------------------------
      type(particle_ensemble), intent(inout) :: par_ens
      integer, intent(in)                    :: last
      if ((last>0).and.(last<=size(par_ens%allpart))) then        
         par_ens%last = last
      else
         write(*,*) "set_last_particle: illegal value for last = ",last
         stop
      endif
      end subroutine set_last_particle

      
      subroutine get_last_particle_number(par_ens, last)   
c     ---------------------------------------------------
      type(particle_ensemble), intent(in) :: par_ens
      integer, intent(out)                :: last
      last = par_ens%last
      end subroutine get_last_particle_number


      function get_particle(par_ens, ipar)
c     --------------------------------------------------- 
c     Checkout a particle from the stack
c     Most efficient usage is:
c 
c       type(particle), pointer :: a_particle
c       a_particle => get_particle(ensemb, integer)
c
c     Assume integrity of par_ens is OK.
c     No not warn, if (ipar > last)
c     Currently no check for multiple access to particles
c     ---------------------------------------------------
      type(particle), pointer  :: get_particle  
      type(particle_ensemble)  :: par_ens
      integer                  :: ipar
      if ((ipar<1).or.
     +    (ipar>size(par_ens%allpart))) then        
         write(*,*) "get_particle: illegal value ipar = ",ipar
         stop
      endif
      
      if ((verbose>0).and.(ipar>par_ens%last)) 
     +     write(*,*) "get_particle: warning, particle > last"

      get_particle => par_ens%allpart(ipar)

      end function

      
      subroutine get_particle_position(part, xyz)
c     ---------------------------------------------------
      type(particle), intent(in) :: part
      real, intent(out)          :: xyz(:)
      call get_tracer_position(part%space, xyz)
      end subroutine get_particle_position



      subroutine insert_particle(par_ens, part, pnum) 
c     ---------------------------------------------------
c     Insert particle part to ensemble par_ens at position pnum.
c     (not deleting currnt particle at pnum, but pushing
c     subsequent particle one position to the right)
c
c     If pnum is absent, add to the end of par_ens
c     Assume integrity of par_ens is OK and
c     that par_ens is initialized
c     ---------------------------------------------------
      type(particle_ensemble)        :: par_ens
      type(particle), target         :: part
      integer,intent(in),optional    :: pnum
c
      integer  :: ipos,last,ipar
c     ---------------------------------------------------
      last = par_ens%last ! short hand
      if (size(par_ens%allpart) <= last) then
         write(*,*) "add_particle: stack full", pnum
         stop
      endif
c
c     Resolve/check position ipos where to put particle part
c
      if (present(pnum)) then
         if ((pnum>0).and.(pnum < last)) then
            ipos=pnum
         else
            write(*,*) "add_particle: illegal pnum", pnum
            stop
         endif
      else  ! pnum absent, apply default (insert at end)
         ipos = last+1
      endif      
c
c     Insert part at position ipos; if this is an interior position
c     move subsequent particles one position to the right)
c
      do ipar = last, ipos, -1
         par_ens%allpart(ipar+1) = par_ens%allpart(ipar) ! copy data 
      enddo
      
      par_ens%allpart(ipar) = part    ! copy data 
      par_ens%last = par_ens%last + 1 ! range OK
      
      end subroutine insert_particle

       


      subroutine delete_particle(par_ens, pnum)  
c     ---------------------------------------------------
c     Assume integrity of par_ens is OK and
c     that par_ens is initialized
c     Invoke delete subroutine for 
c     spatial_attributes/state_attributes components of 
c     particle pnum, and remove reference to it in
c     particle_ensemble stack
c     ---------------------------------------------------
      type(particle_ensemble) :: par_ens
      integer,intent(in)      :: pnum
c     ---------------------------------------------------
      integer  :: last, ipar
c     ---------------------------------------------------
      last = par_ens%last ! shorthand
      if ((pnum<1).or.(pnum>last)) then 
         write(*,*) "delete_particle: illegal pnum", pnum
         stop
      endif
c
c     Invoke delete method for spatial/state attributes
c     to allow feeing allocated potential pointer components in attributes
c 
      call delete_spatial_attributes(par_ens%space_stack(pnum))  
      call delete_state_attributes(par_ens%state_stack(pnum))   
c
c     copy subsequent particles one position to the left
c      
      do ipar = pnum, last-1
         par_ens%allpart(ipar) = par_ens%allpart(ipar+1)
      enddo
      
      par_ens%last = par_ens%last - 1           ! range OK

      end subroutine delete_particle




      subroutine update_particles(par_ens, time_step)
c     ---------------------------------------------------
c     OK for empty particle stack
c
c     Current time retrieved from physical_fields to 
c     ensure syncronicity of time and data
c     ---------------------------------------------------
      type(particle_ensemble),intent(inout) :: par_ens
      real,intent(in)                       :: time_step
c
      type(particle), pointer   :: this
      integer                   :: ip,last
      real                      :: v_active(3)
c     ---------------------------------------------------
      do ip = 1, par_ens%last
         this => par_ens%allpart(ip)
         call get_active_velocity(this%state, this%space, v_active)
         call update_tracer_position(this%space, time_step, v_active) 
         call update_particle_state(this%state, this%space, time_step)
      enddo   
      end subroutine update_particles


      
      subroutine write_ensemble(par_ens) 
c     ---------------------------------------------------
c     ---------------------------------------------------
      type(particle_ensemble),intent(in)   :: par_ens
      type(particle), pointer              :: this
      integer :: i
      do i=1,par_ens%last
         this => get_particle(par_ens, i)
         write(*,'("------------------------------------------")')
         write(*,'("particle ",i6," spatial_attributes = ")') i
         call write_spatial_attributes(this%space)
         write(*,'("particle ",i6," state_attributes = ")') i
         call write_state_attributes(this%state) 
      enddo
      end subroutine write_ensemble 


      end module
