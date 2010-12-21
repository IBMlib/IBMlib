ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Particles tracking
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     TODO: update doc
c     old vars clean up
c     spatial_attributes access/setting
c     
c     Provide basic functionality for the spatial and temporal
c     aspect of tracer motion and emission. 
c
c     Deal with the spatial aspects of a particle which we generally refer to as 
c     a tracer (a particle generally has internal states beyond
c     space/time properties).
c     Rewrite of the tracer tracking module. New features:
c
c        1) A tracer class (spatial_attributes) representing motion state 
c           of a particle has been introduced.
c           It may (should) be used as component in derived (i.e. composite) classes
c
c        2) No tracer stack in module. All operations in the module are on single 
c           tracers now, because this module do not know about derived classes
c           that will be used in practice
c
c     ------------------------
c     Jan 27 2009: Backtracking based on local time symmetry
c                  of advective/diffusive processes added (MPA/ASC)
c     Jul 11 2009: behavior of emission_box type is being changed, 
c                  so that burst releases are supported as well.
c                  Tracer number releases now from a release schedule
c                  instead of rate*time_step (ASC)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module particle_tracking  
      use physical_fields
      use run_context, only: simulation_file
      use input_parser
      use time_tools
      use constants
      use geometry
      use random_numbers
      use output    !Provides data handling classes
      
      implicit none
      private          ! default visibility

c --------------------------------------------------------
c ---------------- module data section -------------------
c --------------------------------------------------------

      

c.....Implemented options for advec_intg_method.................
c   
c     input tag/advec_intg_method      description  
c
c      euler                    Euler forward integration      
c      rk2                      Runge-Kutta 2nd order (midpoint)  
c      rk4                      Runge-Kutta 4th order 
c...............................................................
      integer            :: advec_intg_method    
      integer, parameter :: euler = 10
      integer, parameter :: rk2   = 20
      integer, parameter :: rk4   = 30


      integer, parameter :: verbose = 0    ! control output volume of certain subroutines


      integer, parameter, public :: BC_reflect = 0 ! handle to invoke reflecting boundary conditions
      integer, parameter, public :: BC_sticky  = 1 ! handle to invoke sticky boundary conditions
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....  Neutral tracer attribute components .....................
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Can be used isolated or as component in a derived type
c     ---------------------------------------
c     position(3) is current position as (lon(deg E),lat(deg N),depth(m))
c
c     mobility(3) describes along (lon,lat,depth) (0=immobile, 1=mobile)
c                 use numeric type to allow mask operation
c
c     motion state: ashore, outofdomain, atbottom, atsurface
c                   relates to current value of position
c                   set by step checking algortihms in this module
c
c     shoreBC, domainBC, bottomBC, surfaceBC:
c     Boundary condition (BC) violation control tags:
c       BC_reflect: elastic reflection
c       BC_sticky:  sticky
c     BC handler tags are used in add_constrained_step when checking a motion step
c
c     ---------------------------------------
      type spatial_attributes
      private     ! hide internal implementation
        real              :: position(3)   ! current position as (lon,lat,depth)
        integer           :: mobility(3)   ! along lon,lat,depth
        logical           :: ashore, outofdomain, atbottom, atsurface ! motion state
        integer           :: shoreBC, domainBC, bottomBC, surfaceBC   ! BC handler tags
      end type
      public ::  spatial_attributes


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....  Tracer emitter class  ..................................
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      type emission_box 
      private                                    ! hide internal implementation
        integer           :: emission_boxID      ! unique ID stamp issued by create_emission_boxes
        type(time_period) :: emission_interval
        real              :: SW(3), NE(3)        ! grid coordinates of box corners
        integer           :: max_emissions       ! cap on number of tracers for this box
        integer           :: current_emissions   ! number of tracers emitted until now
        character*999     :: initialization_data ! unparsed data, possibly void
      end type
      public ::  emission_box 

c     Counter for issuing stamps by create_emission_boxes to emission_box instances  
c     emission_box_counter applies to next box created (not last)
c     so that total number created currently is emission_box_counter - 1

      integer :: emission_box_counter   
      
c.... Define public operator set  ..............................     

      
      public :: init_particle_tracking
      public :: close_particle_tracking

      public :: update_tracer_position
      public :: add_advection_step
      public :: add_diffusion_step
      public :: add_constrained_step

      public :: write_tracer_position !TODO: Move all write functionality to output system
      public :: write_spatial_attributes   !TODO: Move all write functionality to output system
      interface get_property
        module procedure get_prop_space
      end interface
      public :: get_property
      public :: get_metadata_space

      public :: delete_spatial_attributes
      public :: set_tracer_mobility
      public :: set_tracer_mobility_free
      public :: set_tracer_mobility_stop
      public :: get_tracer_position
      public :: set_tracer_position
c     .... currently domainBC can only be sticky
      public :: set_shore_BC  
      public :: set_bottom_BC 
      public :: set_surface_BC
      
      public :: create_emission_boxes
      public :: activate_emission_box
      public :: write_emission_box
      public :: get_current_emissions
      public :: get_emission_boxID
      public :: get_initialization_data
      public :: get_xyboxes_emissions

c...............................................................

c     ============================================================
                          contains
c     ============================================================
    
      subroutine init_particle_tracking()
c     ---------------------------------------------------
c     ---------------------------------------------------
      character(len=32)  :: intg_tag
      integer            :: idum3(3), ntt
c     ---------------------------------------------------
      call read_control_data(simulation_file, 
     +                       "advec_intg_method", intg_tag)
      intg_tag = adjustl(intg_tag)
      if     (intg_tag(1:5) == "euler") then
          advec_intg_method = euler
      elseif (intg_tag(1:3) == "rk2")   then
          advec_intg_method = rk2
      elseif (intg_tag(1:3) == "rk4")   then
          advec_intg_method = rk4
      else
         write(*,*) "init_particle_tracking: advec_intg_method " // 
     +        trim(adjustl(intg_tag)) // "not implemented"
         stop
      endif
      
      write(*,*) "init_particle_tracking: advec_intg_method ", 
     +            trim(intg_tag)

      emission_box_counter = 1   ! applies to next box created

      end subroutine init_particle_tracking


      
      subroutine close_particle_tracking()
c------------------------------------------------------------
c------------------------------------------------------------
      end subroutine close_particle_tracking


      subroutine checkstat(istat, who)
c------------------------------------------------------------
c     Local subroutine to homogenize error condition 
c     handling (module private method)
c     who is the name of the calling subroutine/function
c------------------------------------------------------------
      integer, intent(in) :: istat
      character*(*)       :: who
c------------------------------------------------------------
      if (istat /= 0) then
         write(*,*) who, "error=", istat
         stop
      endif
      end subroutine checkstat



      subroutine update_tracer_position(tracattr, h, v_active)
c------------------------------------------------------------
c     Wrapper for invoking all components of spatial update
c     for tracer described by tracer attributes tracattr
c
c     h: forward time step in seconds. h<0 for backtracking
c     v_active is the active motion velocity vector oriented (assumed constant during h)
c     in the current tangent space in meters/second   
c     
c     Return without action in case of the tracer being ashore/outofdomain/atbottom
c     Allow motion for state atsurface
c------------------------------------------------------------
      type(spatial_attributes), intent(inout) :: tracattr
      real,intent(in)                         :: h    ! time step in seconds. h<0 for backtracking
      real,intent(in)                         :: v_active(3) ! active motion velocity meters/second
      real                                    :: dR(3) ! Cartesian step      
      logical                                 :: surfenc, botenc
      logical                                 :: shoreenc, domainenc
c------------------------------------------------------------
     
c.....Allow motion for "atsurface", but ignore motion for other special states 
c     and immobile tracers
      if (tracattr%ashore .or. 
     +    tracattr%outofdomain .or.
     +    tracattr%atbottom .or.
     +    all(tracattr%mobility == 0)) return

c.....Capture sea surface fluctuation artifacts near sea bed          
      call renormalize_vertical_position(tracattr)

c.....evaluate position increments in tangent space additively into dR
      dR = h * v_active * tracattr%mobility
c      write(*,*) tracattr%position
      call add_advection_step(tracattr, h, dR) ! Cartesian step  
      call add_diffusion_step(tracattr, h, dR) ! Cartesian step   
c
c.....apply boundary consitions to position increments dR and
c     add it to tracattr
c     
      call add_constrained_step(tracattr, dR, 
     +            surfenc, botenc, shoreenc, domainenc)
      
      end subroutine update_tracer_position


      subroutine renormalize_vertical_position(tracattr)
c-------------------------------------------------------------
c     For hydrodynamic data with a free surface,
c     particles very close to the sea bed may 
c     artificially get displaced below the sea bed by 
c     sea surface fluctuations
c     Lift these cases to just above (== wlift) the sea bed
c     Do nothing, if ashore/outofdomain/atbottom or vertically frozen
c-------------------------------------------------------------
      type(spatial_attributes), intent(inout) :: tracattr
      real    :: wd
      integer :: istat
      real, parameter :: wlift = 0.01 ! [m] lift affected particles wlift above sea bed
c-------------------------------------------------------------
      if (tracattr%ashore .or. 
     +    tracattr%outofdomain .or.
     +    tracattr%atbottom .or.
     +    (tracattr%mobility(3) == 0)) return
      call interpolate_wdepth(tracattr%position,wd,istat)  
      call checkstat(istat, "renormalize_vertical_position:"//
     +               "interpolate_wdepth ")
      if (tracattr%position(3)>wd) tracattr%position(3) = wd-wlift
      end subroutine renormalize_vertical_position



      subroutine add_advection_step(tracattr, h, dR)
c-------------------------------------------------------------------------------
c     Add advection step to dR (as Cartasian step in meters**3) 
c     using Euler/RK2/RK4 of the tracer described by tracattr
c     without checking for BC violations
c
c     Predictor steps for current interpolation is based tangent space arithmics.
c     Current position of tracer is assumed valid (within bounds, but not checked)
c     h is the time step. h>0 for forward time simulation, h<0 for backtracking
c 
c     Low level method: assert tracattr is in a consistent state at a valid position
c     (not checked)
c
c     Nov 18, 2010 (asc): added test to make sure that predictor step(s) are
c                         valid wet points, else rk2/rk4 defaults to an Euler step
c-------------------------------------------------------------------------------
      type(spatial_attributes), intent(inout) :: tracattr
      real,intent(in)                         :: h     ! time step in seconds. h<0 for backtracking
      real,intent(out)                        :: dR(3) ! Cartasian step in meters**3
c-------------------------------------------------------------------------------      
      real            :: dpos(3)     ! automatic arrays (OK for last == 0)
      real            :: k1(3),k2(3),k3(3),k4(3) ! automatic array for rk1-4
      real            :: pos(3), dxyz(3)
      integer         :: istat
      logical         :: valid_predictor
      real            :: wd
c--------------------------------------------------------

c
c     resolve Cartesian position change dxyz(3) for either 
c     integration scheme Euler/RK2/RK4 
c     Apply tangent space arithmics for predictor steps
c     Assert tracattr%position is a valid wet point

      call interpolate_currents(tracattr%position,k1,istat) ! real(3) in physical units
      call checkstat(istat, "add_advection_step:interpolate_currents ")
      
c     -----------------------------------------
c.....1) Euler forward integration
c        k1  = dt*v(x)
c        x  += k1
c     -----------------------------------------
      if (advec_intg_method == euler) then

         dxyz = k1*h 


c     -----------------------------------------
c.....2) Runge-kutta 2nd order forward integration (may be optimized)
c        k1  = dt*v(x)
c        k2  = dt*v(x + 0.5*k1)
c        x  += k2
c     -----------------------------------------
      elseif (advec_intg_method == rk2) then

         pos  = tracattr%position

         dpos = 0.5*k1*h
         call dcart2dxy(pos, dpos) ! tangent space arithmics
         pos = pos + dpos          ! predictor position
         call interpolate_currents(pos, k2, istat)
         if (istat == 0) then
            dxyz  = k2*h   ! predictor OK
         else
            dxyz  = k1*h   ! predictor step not OK, default to Euler step
         endif
            
c     -----------------------------------------         
c.....3) Runge-kutta 4th order forward integration (may be optimized)
c       k1  = dt*v(x)
c       k2  = dt*v(x + 0.5*k1)
c       k3  = dt*v(x + 0.5*k2)
c       k4  = dt*v(x + k3)
c       x  += k1/6. + k2/3. + k3/3. + k4/6.
c     -----------------------------------------
      elseif (advec_intg_method == rk4) then
         
         pos  = tracattr%position
         valid_predictor = .true.

         dpos = 0.5*k1*h
         call dcart2dxy(pos, dpos) ! tangent space arithmics
         pos = pos + dpos
         call interpolate_currents(pos, k2, istat) 
         if (istat /= 0) valid_predictor = .false.

         if (valid_predictor) then
            dpos = 0.5*k2*h
            pos  = tracattr%position
            call dcart2dxy(pos, dpos) ! tangent space arithmics
            pos = pos + dpos
            call interpolate_currents(pos, k3, istat) 
            if (istat /= 0) valid_predictor = .false.
         endif

         if (valid_predictor) then
            dpos = k3*h
            pos  = tracattr%position
            call dcart2dxy(pos, dpos) ! tangent space arithmics
            pos = pos + dpos
            call interpolate_currents(pos, k4, istat) 
            if (istat /= 0) valid_predictor = .false. 
         endif
         
         if (valid_predictor) then
            dxyz = (k1/6.0d0 + k2/3.0d0 + k3/3.0d0 + k4/6.0d0)*h ! predictor OK
         else
            dxyz  = k1*h   ! predictor step not OK, default to Euler step
         endif

         
      endif
      
      dR = dR + dxyz*tracattr%mobility  ! Cartasian step in meters**3
         
      end subroutine add_advection_step

    

      subroutine add_diffusion_step(tracattr, h, dR)
c-------------------------------------------------------------------------------
c     Add a random walk step (corresponding to time step h) to dR 
c     (as Cartasian step in meters**3) corresponding
c     to physical diffusion in the vertical/horizontal diffusivity field
c     provided by module physical_fields 
c     (see subroutines interpolate_turbulence / interpolate_turbulence_derivative)
c     as described by Visser MEPS 158 pp275-281 (1997)
c     Predictor steps for turbulance interpolation is based 
c     tangent space arithmics.
c     Current position of tracer is assumed valid (within bounds, but not checked)
c     
c     TODO: 
c        1) generate R from a normal distribution rather than uniform
c        2) use a well-defined generator instead of vendor call random_number()
c        3) perform range check on "pos + uphill" before it is used+-
c
c     BACKTRACK: h>0 for forward time simulation, h<0 for backtracking 
c                for the diffusive step, the absolute value of the time step is 
c                used since a stationary random walk is time symmetric 
c
c     Low level method: assumes tracattr is in a consistent state at a valid position
c     (not checked)
c
c     Nov 18, 2010 (asc): fix over-correction bug in turbulence gradient correction to RW
c     Nov 18, 2010 (asc): added test to make sure that uphill step is
c                         valid wet points, else defaults to naive RW step
c-------------------------------------------------------------------------------
      type(spatial_attributes), intent(inout) :: tracattr
      real,intent(in)                         :: h     ! time step in seconds. h<0 for backtracking
      real,intent(out)                        :: dR(3) ! Cartasian step in meters**3
c-------------------------------------------------------------------------------         
      real            :: pos(3), uphill(3), rn3(3)
      real            :: k(3), dk(3), abs_hsca(3), abs_h
      real            :: depth, vstep, dxy(2), dxyz(3), x,y,z,znew
      integer         :: itr, istat
c--------------------------------------------------------
      if (all(tracattr%mobility == 0)) return 
      abs_h  = abs(h)

c     Assert tracattr%position is a valid wet point
      call interpolate_turbulence_deriv(tracattr%position, dk, istat) 
      call checkstat(istat, 
     +               "add_diffusion_step:interpolate_turbulence_deriv") 
c         
c.....Use tangent space approx to estimate uphill point where to 
c     evaluate diffusivity (gives higher order deviation)
c
      uphill = 0.5d0*dk*abs_h   ! step where to evaluate k
      pos    = tracattr%position 
      call dcart2dxy(pos, uphill) ! transform uphill to (lon,lat,z)
      pos    = pos + uphill     
      call get_random_number(rn3) ! avg(R)=0, Var(R)=1
c
c.....Attempt gradient corrected RW step (MEPS 158 pp275-281 (1997), Eq 6)
c     If uphill position is not valid, default to 
c     naive RW step (MEPS 158 pp275-281 (1997), Eq 2)
c
      call interpolate_turbulence(pos, k, istat) 
      if (istat == 0) then ! uphill step OK
         dxyz = dk*abs_h + rn3*sqrt(2.0d0*k*abs_h) ! factor 2 corresponding to Var(R)=1
      else                 ! uphill step not OK
         call interpolate_turbulence(tracattr%position, k, istat)
         call checkstat(istat, 
     +               "add_diffusion_step:interpolate_turbulence") 
         dxyz = rn3*sqrt(2.0d0*k*abs_h) ! factor 2 corresponding to Var(R)=1
      endif 

      dR   = dR + dxyz*tracattr%mobility      ! Cartasian step in meters**3

      end subroutine add_diffusion_step
      


      subroutine add_constrained_step(tracattr, dR, 
     +           surfenc, botenc, shoreenc, domainenc)
c--------------------------------------------------------------------------
c     Check whether proposed position increment dR (provided in 
c     Cartesian coordinates in meters**3) and add it, if it satisfies 
c     boundary conditions; otherwise manipulate dR before adding it.
c   
c     tracattr(in)  : tracer being moved
c     dR(in)        : proposed Cartesian position increment
c     surfenc(out)  : result of application of surface boundary condition on dR
c     botenc(out)   : result of application of bottom boundary condition on dR
c     shoreenc(out) : result of application of shore line boundary condition on dR
c     domainenc(out): result of application of domain violation boundary condition on dR
c     
c     There are three elements in handling boundary conditions:
c       A) Boundary condition control tags (tracattr%*BC): "what to do, if ..."
c       B) Tracer compliance with boundary condition (tracattr%(ashore|outofdomain|atbottom|atsurface)
c       C) Return flags telling whether a boundary violation has
c          occured in current step: surfenc, botenc, shoreenc, domainenc
c
c     If dR satisfies boundary conditions add it to tracattr%position, update
c     tracattr state flags and set BC flags (surfenc, botenc, shoreenc, domainenc) to false; 
c     otherwise truncate step according to applicable boundary conditions and 
c     add truncated step to tracattr, update tracattr state flags and set 
c     BC flags (surfenc, botenc, shoreenc, domainenc) accordingly.
c     Boundary condition control tags (tracattr%*BC) have the
c     following options and meaning
c
c       0: elastic reflection of step => BC state = .FALSE.
c       1: sticky tracer              => BC state = .TRUE.
c     
c     Current position (tracattr%position) is assumed valid before considering dR.
c     No tracers violating BC are removed by this subroutine
c     The algorithm below is based on separability of horizontal/vertical BC analysis
c
c     Boundary condition checks are performed in this order:
c
c     1) Horizontal domain boundary crossings: 
c
c        Response controller :
c          tracattr%domainBC = BC_reflect|BC_sticky (currently only BC_sticky is accepted)
c        Return state flags:
c          domainenc = .FALSE.: no domain boundary crossings, 
c          domainenc = .TRUE. : at least one domain boundary crossed
c        Tracer state is updated after application of BC to step:
c          tracattr%outofdomain = .TRUE.  => tracattr%mobility = 0
c          tracattr%outofdomain = .FALSE. => tracattr%mobility unchanged
c       
c
c     2) Tracers advected ashore:
c
c        Response controller:
c          tracattr%shoreBC  = BC_reflect|BC_sticky 
c        Return state flags:
c          shoreenc = .FALSE.: no coast lines crossed by dR
c          shoreenc = .TRUE. : at least one coast lines crossed
c        Tracer state is updated after application of BC to step:
c          tracattr%ashore = .TRUE.  => tracattr%mobility = 0
c          tracattr%ashore = .FALSE. => tracattr%mobility unchanged
c         
c     3) Vertical (z) boundary crossings: 
c
c        Response controller:
c          tracattr%surfaceBC = BC_reflect|BC_sticky
c        Return state flags:
c          surfenc = .FALSE.: water surface not crossed not by dR
c          surfenc = .TRUE. : water surface crossed at least once by dR
c        Tracer state is updated after application of BC to step:
c          tracattr%atsurface = .TRUE.  => tracattr%mobility = (/1,1,0/)
c                   (this allows particle to flow with surface currents)
c          tracattr%atsurface = .FALSE. => tracattr%mobility unchanged
c
c        Response controller:
c          tracattr%surfaceBC = BC_reflect|BC_sticky
c        Return state flags:
c          botenc = .FALSE.:  sea bed not crossed not by dR
c          botenc = .TRUE. :  sea bed crossed at least once by dR
c        Tracer state is updated after application of BC to step:
c          tracattr%atbottom = .TRUE.  => tracattr%mobility = 0
c          tracattr%atbottom = .FALSE. => tracattr%mobility unchanged
c          
c     Nov 19, 2010: added multiple reflection analysis
c     Nov 19, 2010: added (verbose>0) debug hooks           
c   
c     TODO: validate that particle with tracattr%atsurface = .TRUE. at entry
c           is handled consistent
c--------------------------------------------------- 
      type(spatial_attributes), intent(inout) :: tracattr
      real, intent(in)                        :: dR(3)      
      logical                                 :: surfenc, botenc
      logical                                 :: shoreenc, domainenc

 
      real    :: pos(3), virpos(3), curpos(3), xyzref(3), xyzhit(3)
      integer :: ipos, istat
      integer :: ixc,iyc,izc, ixv,iyv,izv, dimOK(3), ibot,i
      real    :: rnum(3),r1,zv,znew, depth,dxy(3)
      logical :: isOK, anycross, pingpong
c---------------------------------------------------

c     
c     Resolve the two end points of the attempted move:
c
c     curpos = current position at entry of add_constrained_step
c
c     virpos = current position with dR added. virpos is manipulated
c              iteratively below so that it finally satisfies all BC
c              At exit, virpos is assigned to tracattr%position. 
c
      if (verbose>0) then
         
         write(*,387) "add_constrained_step: begin analysis",
     +           "proposed Cartesian step dR = ", dR
 387     format(5("-"),1x,a,1x,5("-"),/,a, 3f7.2,1x,5("-"))
         write(*,*) "spatial_attributes at entry ="
         call write_spatial_attributes(tracattr)
      endif
      

      curpos = tracattr%position ! current position - assumed valid
      virpos = curpos
      dxy    = dR
      call dcart2dxy(virpos,dxy)
      virpos = virpos + dxy

      if (verbose>0) write(*,*) "initial virpos = ", virpos

c.....reset motion status flags 
      tracattr%ashore      = .FALSE. 
      tracattr%outofdomain = .FALSE. 
      tracattr%atbottom    = .FALSE. 
      tracattr%atsurface   = .FALSE.   ! ... should this be reset ?

c.....Give all step status flags a default value (no BC violations)
      surfenc   = .FALSE.
      botenc    = .FALSE.
      shoreenc  = .FALSE.
      domainenc = .FALSE.
  
      
c     -----------------------------------------------------   
c     1) Check boundary crossings 
c        If a domain violation has occured, assign 
c        virpos as current position and raise 
c        tracattr%outofdomain. Do not check other
c        boundary conditions in case a domain violation 
c        has occured. Do not try to determine an exact
c        location where the domain violation has occured,
c        but just assign first illegal position (we could also
c        apply last legal and reject virpos ...)
c        We currently do not accept domainBC == BC_reflect
c     -----------------------------------------------------
   
      isOK = horizontal_range_check(virpos)
      
      if (verbose>0) write(*,*) "horizontal_range_check =", isOK

      if (.not.isOK) then
         domainenc = .TRUE.
         tracattr%outofdomain = .TRUE.
         if (tracattr%domainBC == BC_sticky) then
            tracattr%mobility = 0        ! freeze tracer
         else
            write(*,*) "add_constrained_step:invalid tracattr%domainBC",
     +                  tracattr%domainBC
            stop
         endif
c
c........Elwis has left the building, return to sender 
c
         tracattr%position = virpos ! set illegal position

         if (verbose>0) then
            write(*,*) "exit due to horizontal range violation"
            write(*,*) "spatial_attributes at exit(1) ="
            call write_spatial_attributes(tracattr)
         endif
         return                     ! no further BC processing when outofdomain

      endif


c     -----------------------------------------------------   
c     2) Detect tracers advected ashore (horizontal component):
c
c        Response controller:
c          tracattr%shoreBC = BC_reflect|BC_sticky
c        Return state flags:
c          shoreenc = .FALSE. : no coast lines crossed by dR
c          shoreenc = .TRUE.  : at least one coast lines crossed
c        Tracer state is updated after application of BC to step:
c          tracattr%ashore = .TRUE.  => tracattr%mobility = 0
c          tracattr%ashore = .FALSE. => tracattr%mobility unchanged 
c
c        Handle the situation that both vertical and horizontal 
c        BC are violated       
c     -----------------------------------------------------   

      if (verbose>0) write(*,*) "Horizontal analysis begin"

      call multiple_reflection_path(curpos, virpos, 
     +            anycross, xyzref, xyzhit) 
      

      if (anycross) then

         if (verbose>0) write(*,*) "coast line crossing detected" 

         shoreenc = .TRUE.
         if      (tracattr%shoreBC == BC_sticky) then
            virpos    = xyzhit
            virpos(3) = 0. ! put ashore tracer sat surface
            tracattr%position = xyzhit
            tracattr%mobility = 0
            tracattr%ashore   = .TRUE.
c...........BC exit point
            
            if (verbose>0) then
               write(*,*) "particle washed ashore (shoreBC==sticky)"
               write(*,*) "spatial_attributes at exit(2) ="
               call write_spatial_attributes(tracattr)
            endif
            return         ! do not consider vertical BC, when ashore

         elseif  (tracattr%shoreBC == BC_reflect) then

            virpos = xyzref ! update virpos - proceed to vertical BC

            if (verbose>0) then
               write(*,*) "particle reflected on coast line virpos=",
     +                    virpos
            endif   

         else

            write(*,*) "add_constrained_step: invalid tracattr%shoreBC",
     +                  tracattr%shoreBC
            stop 
  
         endif
      else  ! if anycross
         if (verbose>0) write(*,*) "No coast line crossings detected" 
      endif

      if (verbose>0) write(*,*) "After horizontal analysis: virpos =", 
     +                           virpos 

c     -----------------------------------------------------   
c     3) Enforce vertical boundary conditions, if tracer is not 
c        advected ashore:
c
c        0(surface) < z <  wdepth(xyz)
c
c        At this point, the horizontal component of virpos corresponds to a wet point
c        Assume sea bed is so flat that we need not compute exactly
c        where bottom is hit, but may use the depth at virpos
c        to enforce vertical boundary conditions
c     ----------------------------------------------------- 

      call interpolate_wdepth(virpos, depth, istat)
      call checkstat(istat,"add_constrained_step:interpolate_wdepth")

      if (verbose>0) then
         write(*,*) "Vertical analysis begin: wdepth at virpos =", depth                   
      endif   

      zv    = virpos(3)
      pingpong = .FALSE.
c
c.....check for and handle surface crossing
c
      if (zv<0) then            ! a surface crossing has occured

         if (verbose>0) write(*,*) "surface crossing detected"

         surfenc              = .TRUE.
         if     (tracattr%surfaceBC == BC_sticky) then
            virpos(3)  = 0.0    ! place tracer at surface
            tracattr%mobility = (/1,1,0/) ! freeze tracer to water surface
            tracattr%atsurface = .TRUE.
         elseif (tracattr%surfaceBC == BC_reflect) then
            zv        = -zv      ! surface reflection
            virpos(3) =  zv      ! handle ping pong situation at bottom, if relevant  
         else
            write(*,*) "add_constrained_step: invalid "//
     +                 "tracattr%surfaceBC", tracattr%surfaceBC
            stop   
         endif
         if (verbose>0) write(*,*) "updated virpos = ", virpos
      endif
c
c.....check for and handle bottom crossing
c       
      if (zv>depth) then        ! a bottom crossing has occured

         if (verbose>0) write(*,*) "bottom crossing detected"

         botenc  = .TRUE.
         if     (tracattr%bottomBC == BC_sticky) then
            virpos(3)         = depth    ! place tracer at bottom
            tracattr%mobility = 0 ! freeze tracer to bottom
            tracattr%atbottom = .TRUE.
         elseif (tracattr%bottomBC == BC_reflect) then
            zv = 2.0*depth - zv ! bottom reflection
            if (surfenc.or.(zv<0)) then
               pingpong = .TRUE.
            else
               virpos(3) = zv
            endif
         else
            write(*,*) "add_constrained_step: invalid "//
     +                 "tracattr%bottomBC", tracattr%bottomBC
            stop   
         endif  

         if (verbose>0) write(*,*) "updated virpos = ", virpos

      endif
c
c.....handle a flagged ping pong situation
c         
      if (pingpong) then

         if (verbose>0) write(*,*) "pingpong condition detected"

         surfenc   = .TRUE.
         botenc    = .TRUE.
         call random_number(znew) ! compiler provided uniform RNG (0<znew<1)
         znew = znew*depth
         virpos(3) = znew  

         if (verbose>0) write(*,*) "updated virpos = ", virpos
      endif
c
c.....now virpos has been massaged so that all BC are satisfied
c     assign it here, if no previous exit point has been applied
c     
      if (verbose>0) write(*,*) "final truncated virpos = ", virpos

      tracattr%position = virpos

      if (verbose>0) then
         write(*,*) "spatial_attributes at exit(3) ="
         call write_spatial_attributes(tracattr)
         write(*,*) 
      endif


      end subroutine add_constrained_step



      subroutine multiple_reflection_path(s0, s1, anycross, sref, shit1)  ! update to particle_tracking.f 
c-------------------------------------------------------------------------
c     Compute key points of multiple horizontal coastal reflection path
c     by iterative application of coast_line_intersection primitive
c     when trying to move from (valid, wet) s0 to new position s1
c     If a coast line is crossed retrun anycross == .true.
c     Then sref is the (multiple) horizontal coastal reflected end point
c     which is valid and wet. shit1 is the point of (first) coastal intersection
c     If anycross == .false. (sref, shit1) is undefined
c-------------------------------------------------------------------------
      real,intent(in)      :: s0(3), s1(3)
      logical,intent(out)  :: anycross
      real,intent(out)     :: sref(3), shit1(3)

      real                 :: start(3), sfin(3), shit(3)
      integer              :: iter
      integer, parameter   :: max_reflections = 5 ! otherwise time step too long ...
c---------------------------------------------------------------
      if (verbose>0) write(*,*) " multiple reflection analysis begin"

      iter  = 0
      call coast_line_intersection(s0, s1, anycross, sref, shit1)

      if (verbose>0) then
         if (anycross) then 
            write(*,422) iter, s0(1:2), s1(1:2), sref(1:2), shit1(1:2)  
            write(*,*) "is_land(reflected point) = ", is_land(sref)
            write(*,*) "continue multiple reflection analysis"
         else
            write(*,423) iter, s0(1:2), s1(1:2)
            write(*,*)"wet reflection: multiple reflection analysis end"
         endif
      endif ! verbose

      if (.not.anycross) return
      
c
c     --- we hit the coast line, start multiple reflection analysis
c         (this means (sref, shit1) are defined and anycross == .true.)
c
      iter  = 1
      shit  = shit1  ! save first coastal hit
c     note: is_land = .false. at horizontal range violation
c     The loop below computes the final reflection, sref
c     In most cases we will not enter the while loop, because sref is a wet point
c
      do while (anycross .and. 
     +          is_land(sref).and.
     +          (iter <= max_reflections))
         start = shit
         sfin  = sref
         call coast_line_intersection(start, sfin, anycross, sref, shit)
         if (.not.anycross) sref=sfin ! roll back, in case sref is overwritten
         if (verbose>0) then
            if (anycross) then 
               write(*,422) iter, start(1:2), sfin(1:2), 
     +                      sref(1:2), shit(1:2)
               write(*,*) "is_land(reflected point) = ", is_land(sref)
            else
               write(*,423) iter, start(1:2), sfin(1:2)
            endif
         endif ! verbose

         iter  = iter + 1
      enddo
      if (iter > max_reflections) then
         write(*,*) "multiple_reflection_path: max_reflections exceeded"
         stop
      endif
 
      if (verbose>0) write(*,*) " multiple reflection analysis end"

 422  format("step ", i3, " :", 2f7.2, " ->", 2f7.2, " : ref=",
     +        2f7.2, " hit=", 2f7.2)    
 423  format("step ", i3, " :", 2f7.2, " ->", 2f7.2, " : no cross")
c-----------------------------------------------------------
      end subroutine
  



      subroutine write_tracer_position(tracattr)
c------------------------------------------------------------
c     write tracer position to logical unit iunit
c------------------------------------------------------------     
      type(spatial_attributes),intent(in) :: tracattr
c------------------------------------------------------------         
      write(*,829) tracattr%position
 829  format(3(f12.7))
      end subroutine write_tracer_position
      


      subroutine write_spatial_attributes(tracattr)
c------------------------------------------------------------
c     write all tracer data to logical unit iunit 
c------------------------------------------------------------  
      type(spatial_attributes),intent(in) :: tracattr
c------------------------------------------------------------  
      write(*,831) tracattr%position,
     +                  tracattr%mobility,
     +                  tracattr%ashore, 
     +                  tracattr%outofdomain, 
     +                  tracattr%atbottom, 
     +                  tracattr%atsurface,
     +                  tracattr%shoreBC, 
     +                  tracattr%domainBC, 
     +                  tracattr%bottomBC, 
     +                  tracattr%surfaceBC
 831  format(3(f12.7),2x,3(i2),2x,4(l2),2x,4(i2))

      end subroutine write_spatial_attributes


      subroutine delete_spatial_attributes(tracattr)
c------------------------------------------------------------
c     Protocol-required subroutine, currently void
c------------------------------------------------------------  
      type(spatial_attributes),intent(inout) :: tracattr
c------------------------------------------------------------  
      end subroutine delete_spatial_attributes




      subroutine set_tracer_mobility(tracattr,mobx,moby,mobz)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(inout) :: tracattr
      integer,intent(in) :: mobx,moby,mobz
      tracattr%mobility(1) = mobx
      tracattr%mobility(2) = moby
      tracattr%mobility(3) = mobz
      end subroutine set_tracer_mobility

      subroutine set_tracer_mobility_free(tracattr)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(inout) :: tracattr      
      tracattr%mobility = 1
      end subroutine set_tracer_mobility_free
  
      subroutine set_tracer_mobility_stop(tracattr)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(inout) :: tracattr
      tracattr%mobility = 0      
      end subroutine set_tracer_mobility_stop



      subroutine get_tracer_position(tracattr,xyz)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(in) :: tracattr
      real, intent(out)                   :: xyz(:)
      xyz(1:3) = tracattr%position    
      end subroutine get_tracer_position


      subroutine set_tracer_position(tracattr,xyz)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(inout) :: tracattr
      real, intent(in)                   :: xyz(:)
      tracattr%position(1:3) = xyz(1:3)
      end subroutine set_tracer_position
      

      subroutine set_shore_BC(tracattr, BChandler)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(out) :: tracattr
      integer,intent(in)                   :: BChandler
c------------------------------------------------------------ 
      if ((BChandler == BC_reflect).or.
     +    (BChandler == BC_sticky)) then
         tracattr%shoreBC = BChandler
      else
         write(*,*) "set_shore_BC: illegal BC =", BChandler
         stop
      endif
      end subroutine set_shore_BC

    

      subroutine set_bottom_BC(tracattr, BChandler)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(out) :: tracattr
      integer                              :: BChandler
c------------------------------------------------------------
      if ((BChandler == BC_reflect).or.
     +    (BChandler == BC_sticky)) then
         tracattr%bottomBC = BChandler
      else
         write(*,*) "set_bottom_BC: illegal BC =", BChandler
         stop
      endif
      end subroutine set_bottom_BC



      subroutine set_surface_BC (tracattr, BChandler)
c------------------------------------------------------------ 
      type(spatial_attributes),intent(out) :: tracattr
      integer                              :: BChandler
c------------------------------------------------------------ 
       if ((BChandler == BC_reflect).or.
     +    (BChandler == BC_sticky)) then
         tracattr%surfaceBC = BChandler
      else
         write(*,*) "set_surface_BC: illegal BC =", BChandler
         stop
      endif   
      end subroutine set_surface_BC     


c============================================================
c     Tracer emission operators
c============================================================

      subroutine calc_number_of_releases(boxref, now, timedir, newpt) 
c------------------------------------------------------------------
c     Compute how many new tracers (newpt) should be ideally 
c     released from this emission_box (boxref) corresponding 
c     to this time (now), corresponding to the settings of boxref
c     and the current number of tracers released from this box
c     Currently we support a linear release schedule within
c     a specified time interval - later we may add complementary 
c     release schedules, like normal distributions
c------------------------------------------------------------------
      type(emission_box),intent(in) :: boxref
      type(clock),intent(in)        :: now
      integer,intent(in)            :: timedir ! 1: forward, -1: backward
      integer,intent(out)           :: newpt
      integer                       :: ntarget
      real                          :: rp
c------------------------------------------------------      
      call get_relative_position(boxref%emission_interval, now, rp)
      if (boxref%current_emissions > boxref%max_emissions) then
         write(*,277) boxref%current_emissions, boxref%max_emissions
 277     format("calc_number_of_releases: "//
     +          "current_emissions=",i8,">max_emissions=",i8)
         stop
      endif
c     ..........................................
      if     (timedir ==  1) then  !  1: forward time
         rp      = max(0.0,min(rp,1.0)) 
      elseif (timedir == -1) then  ! -1: backward time
         rp      = max(0.0,min(1.0-rp,1.0)) 
      else
         write(*,*) "illegal timedir setting =", timedir
         stop
      endif
c
c     --- now 0<rp<1 ---
c     
      ntarget = nint(boxref%max_emissions*rp)
      ntarget = min(ntarget, boxref%max_emissions) ! newer exceed max_emissions
      newpt   = ntarget - boxref%current_emissions ! ideally released this number
      end subroutine calc_number_of_releases




      subroutine create_emission_boxes(identifier, boxref, nmaxpar) 
c---------------------------------------------------
c     Create emission box list corresponding to identifier
c     in the simulation file
c     format:
c
c       identifier = r(15) [bio_init_par]
c
c     where the vector r(15) means:
c
c       r(1:4)   time of emission start:             year month day sec_of_day
c       r(5:8)   time of emission end:               year month day sec_of_day
c       r(9:11)  south west corner of emission box:  lon_min lat_min z_min
c       r(12:14) north east corner of emission box:  lon_max lat_max z_max
c       r(15)    max_number_of_tracers (emitted by this box)
c     
c     z_min, z_max controls release depth; if 0 < z_min, z_max < 1 
c     this is interpreted as a relative emission depth (relative to
c     current depth at this longitude and latitude).
c     In this case tracers will be uniformly released in the
c     interal [(z_min+epsil)*current_depth, (z_max-epsil)*current_depth]
c     where epsil is a very small number (epsil=1.e-6) so that
c     z == 0 is just below the surface and
c     z == 1 is just above the bottom
c
c     if 0 > z_min, z_max this is interpreted as an absolute release depth
c     At initialization time it can not be checked that release depth
c     is above the bottom at all points in the release box - this must be done at
c     run time.
c
c     NEW: bio_init_par is the (optional) remaining part of the line, following
c     the spatio-temporal part r(15). This is just cached in the emission_box
c     and may later be used in the biological state initialization of the tracer
c
c     This initialization will abort if the corners of a
c     release box is dry (until a performance safe scheme is implemented 
c     for handling a mixed wet/dry area).
c     backward simulation: should provide start/end in normal
c     order, i.e. time(start) first, time(start)<time(end)
c
c     Return list of emission boxes in pointer boxref 
c     (data space allocated here to boxref)
c     Each box is issued a runtime unique stamp (integer) emission_boxID
c     Return nullified pointer boxref, if no sources are created
c     The maximum number of tracers that may be emitted by
c     boxes in boxref are summed in nmaxpar 
c
c     Nov 18, 2010 (asc): allow partially/fully dry corners in emission box, but issue warning
c---------------------------------------------------  
      character*(*),intent(in)    :: identifier
      type(emission_box), pointer :: boxref(:)
      integer,intent(out)         :: nmaxpar

      character*999               :: strbuf, the_rest
      integer                     :: nboxes,ihit,ibox
      integer                     :: i,start(500),nwords ! 500 arbitrary
      real                        :: r15(15), xyz(3)
      logical                     :: absspec, anydry, anyillegal
      real, parameter             :: epsil    = 1.0e-6
      real, parameter             :: oneplus  = 1 + epsil
      real, parameter             :: oneminus = 1 - epsil
c--------------------------------------------------- 
      nmaxpar = 0
      nboxes  = count_tags(simulation_file, identifier) 

      write(*,547) nboxes, trim(adjustl(identifier))
 547  format("create_emission_boxes:", i4, " boxes for identifier: ",a)

      if (nboxes > 0) then
         allocate(boxref(nboxes)) ! no check
      else
         nullify(boxref)          ! signal no boxes defined
         return 
      endif
      ihit = 1
c     ---------------------------------------------------------------
c     main loop: read and parse all boxes corresponding to identifier
c     ---------------------------------------------------------------
      do ibox = 1, nboxes
c:old         call read_control_data(simulation_file, identifier, r15, ihit)
         call read_control_data(simulation_file,identifier,strbuf,ihit)
         call tokenize(strbuf, start, nwords)
         if (nwords<15) then 
            write(*,*) "create_emission_boxes: bad value format for", 
     +                  identifier 
            write(*,*) "value=", trim(adjustl(strbuf))
            stop
         endif
         read(strbuf,*) r15   ! the remaining part ignored, if any
c
c        .... cache the remaining part of the line in the emission_box
c        .... by chopping off the 15 first words in the line value
c
         if (nwords == 15) then                    ! i.e. no initialization_data
            boxref(ibox)%initialization_data = ""  ! set an empty string
         else 
            ! chop off the 15 first words of strbuf
            ! I'm nervous about literals in the code
            boxref(ibox)%initialization_data = strbuf(start(15+1):) 
         endif
         
c        --------- parse spatio-temporal part ---------
c........time
c        r15(1:4) start ==   year month day sec_of_day
c        r15(5:8) end   ==   year month day sec_of_day
         call set_period(boxref(ibox)%emission_interval, int(r15(1:8)))

c........spatial coordinates
c        r15(9:11)      ==   lon_min lat_min z_min
c        r15(12:14)     ==   lon_max lat_max z_max
c        
         boxref(ibox)%SW = r15(9:11)
         boxref(ibox)%NE = r15(12:14)
         
c
c........basic check of box consistency
c
c        1) flag whether this is an absolute or relative depth specification
c           (and put depth specifications in standard order)
c           (do not accept a mixed absolute/relative specification)
c
         if ((boxref(ibox)%SW(3)>=0).and.(boxref(ibox)%NE(3)>=0)) then

            absspec = .FALSE.

c...........handle surface specification (z = 0 -> 0+)
            if (boxref(ibox)%SW(3) < epsil) boxref(ibox)%SW(3) = epsil
            if (boxref(ibox)%NE(3) < epsil) boxref(ibox)%NE(3) = epsil

c...........handle bottom specification (z = 1 -> 1-)
            if (abs(boxref(ibox)%SW(3)-1.0) < epsil) then
               boxref(ibox)%SW(3) = oneminus
            endif
            if (abs(boxref(ibox)%NE(3)-1.0) < epsil) then
               boxref(ibox)%NE(3) = oneminus
            endif

         elseif ((boxref(ibox)%SW(3)<0).and.(boxref(ibox)%NE(3)<0)) then
            absspec = .TRUE.
         else
             write(*,*) "create_emission_boxes: invalid specification"
             write(*,*) "z_min, z_max has unequal sign for box", ibox
             call write_emission_box(boxref(ibox))
             stop
         endif   
c
c        2) flag out-of-bound relative depth specification
c
         if ((boxref(ibox)%SW(3) > 1.0).or.    
     +       (boxref(ibox)%NE(3) > 1.0)) then
            write(*,*) "create_emission_boxes: invalid specification"
            write(*,*) "z_min or z_max > 1 for box", ibox
            call write_emission_box(boxref(ibox))
            stop
         endif 
c
c        3) flag out of domain corners 
c
         anyillegal = .FALSE.
         xyz(3) = 0.            ! dummy not used
         xyz(1:2) = boxref(ibox)%SW(1:2)
         if (.not.horizontal_range_check(xyz)) anyillegal = .TRUE.
         xyz(2)   = boxref(ibox)%NE(2)
         if (.not.horizontal_range_check(xyz)) anyillegal = .TRUE.
         xyz(1:2) = boxref(ibox)%NE(1:2)
         if (.not.horizontal_range_check(xyz)) anyillegal = .TRUE.
         xyz(2)   = boxref(ibox)%SW(2)
         if (.not.horizontal_range_check(xyz)) anyillegal = .TRUE.
         
         if (anyillegal) then
            write(*,*) "create_emission_boxes: out-of-domain"
            write(*,*) "box", ibox, "fully/partially or out-of-domain"
            call write_emission_box(boxref(ibox))
            stop
         endif          
c
c        4) Test for dry corners 
c
c           Partially/fully dry corners in emission box are accepted, but 
c           a warning is issued. Boxes with partially/fully dry corners may
c           have wet zones - this is tested at release time, and execution
c           is halted, if no wet zones are found at release time
c
         anydry = .FALSE.
         if (absspec) then
            xyz(3) = max(abs(boxref(ibox)%SW(3)),
     +                   abs(boxref(ibox)%NE(3)))
            xyz(1:2) = boxref(ibox)%SW(1:2)
            if (.not.is_wet(xyz)) anydry = .TRUE.
            xyz(2)   = boxref(ibox)%NE(2)
            if (.not.is_wet(xyz)) anydry = .TRUE.
            xyz(1)   = boxref(ibox)%NE(1)
            if (.not.is_wet(xyz)) anydry = .TRUE.
            xyz(2)   = boxref(ibox)%SW(2)
            if (.not.is_wet(xyz)) anydry = .TRUE.
         else
            xyz(3) = 0. ! dummy not used 
            xyz(1:2) = boxref(ibox)%SW(1:2)
            if (is_land(xyz)) anydry = .TRUE.
            xyz(2)   = boxref(ibox)%NE(2)
            if (is_land(xyz)) anydry = .TRUE.
            xyz(1)   = boxref(ibox)%NE(1)
            if (is_land(xyz)) anydry = .TRUE.
            xyz(2)   = boxref(ibox)%SW(2)
            if (is_land(xyz)) anydry = .TRUE.
         endif
         
         if (anydry) then
            write(*,*) "create_emission_boxes: warning: box", ibox, 
     +                 "is partially/fully dry"
         endif
c
c........initialize internal release counter and emission cap
c        r15(15)        ==   max_number_of_tracers
c
         boxref(ibox)%max_emissions     = nint(r15(15))     
         boxref(ibox)%current_emissions = 0 
         nmaxpar = nmaxpar + boxref(ibox)%max_emissions

         boxref(ibox)%emission_boxID = emission_box_counter
         emission_box_counter = emission_box_counter + 1 

         write(*,*) 
         write(*,*) "created emission_box ", ibox, ":"
         call write_emission_box(boxref(ibox))
c
         ihit = ihit + 1   ! take next box in control file
      enddo 

      end subroutine create_emission_boxes





      subroutine activate_emission_box(emit_box, 
     +                                 time_dir, parbuf, nstart, npar) 
c-----------------------------------------------------------------
c     Emit from box referenced by boxrefs, if current time (master clock in physical fields)
c     is within emission period
c     time_dir > 0 is forward time simulation, time_dir < 0 is 
c     accepted corresponding to a backtracking emission box (not absorbtion box).
c     Tracers are generated to tracer buffer parbuf (allocated externally)
c     starting at offset nstart. npar is the number of tracers 
c     generated in total by this call (and 0 if no tracers are generated)
c     So the output from this call is stored in parbuf(nstart:nstart+npar)
c 
c     activate_emission_box now does not considered whether box is active or
c     not any more, since in future release schemes (e.g. normal distributions)
c     this may not have any meaning
c
c     Jul 11 2009: time_step argument replaced with time_dir, where
c                  only the sign is important (keeps as real, for 
c                  compartability, so time_step can be used); release 
c                  is now determined from current_time in relation 
c                  to emission_interval
c
c----------------------------------------------------------------- 
      type(emission_box),intent(inout)     :: emit_box
      real,intent(in)                      :: time_dir ! >0: forward, <0
      type(spatial_attributes),intent(out) :: parbuf(:)
      integer,intent(in)                   :: nstart ! offset in parbuf
      integer,intent(out)                  :: npar   ! generated tracers
c     ..............................................
      integer, parameter     :: pos_attemps = 10**4  ! max attemps to find a wet point in box
      logical                :: pos_OK, absspec
      integer                :: iatt, nwish, maxp_mem, maxp_box,ip
      integer                :: i, forw_back, istat
      real                   :: dr(3), u3(3), abs_time_step, newpos(3)
      real                   :: curdepth
      type(clock),pointer    :: current_time
c---------------------------------------------------   
      if (time_dir>0) then
         forw_back = 1        ! forward tracking
      else
         forw_back = -1       ! backtrack
      endif

c.....compute the number of tracers to emit (nemit), not exceeding
c     total tracer buffer size and box cap
     
      maxp_mem = size(parbuf) - nstart + 1
      maxp_box = emit_box%max_emissions - emit_box%current_emissions 
      current_time => get_master_clock() 
      call calc_number_of_releases(emit_box, current_time,
     +                             forw_back, nwish)

      npar    = min(nwish, maxp_box, maxp_mem)
      
c
c.....Determine whether vertical release is specified as absolute or relative
c
      if ((emit_box%SW(3)>0).and.(emit_box%NE(3)>0)) then
            absspec = .FALSE.            
      else                ! consistency validated in create_emission_boxes
            absspec = .TRUE.
      endif 

c
c.....Initialize tracers, if any new are created (npar>0
c
      dr   = emit_box%NE - emit_box%SW    ! emission box diagonal vector

check that lon lat change is OK in emitter
check that sign change on release delimiter is OK everywhere

      do ip = nstart, nstart+npar-1         ! OK for npar==0
c
c........First horizontal: make sure tracer is initialized in a wet position
c        (eventhough basic checks are performed at setup stage of 
c         release box, box may contain dry areas at release time)
c
         iatt = 1
         do 
            call random_number(u3) ! uniform derivate, length 3
            newpos = emit_box%SW + dr*u3
c...........test if there is wet points at this horizontal position
            if (.not.is_land(newpos)) exit 
c...........terminate, if upper bound were hit (iatt==pos_attemps+1)
            if (iatt > pos_attemps) then
               write(*,*) "too many dry attempts in release box"
               call write_emission_box(emit_box)
               stop
            endif
            iatt = iatt + 1
         enddo


c........Now we know horizontal position of is newpos wet   
            
         call interpolate_wdepth(newpos, curdepth, istat) 
         call checkstat(istat, 
     +        "activate_emission_box:interpolate_wdepth")

         if (absspec) then
            newpos(3) = -newpos(3)                
         else
            newpos(3) = newpos(3)*curdepth ! scale to absolute
         endif

c........Catch other misreleases
         if (newpos(3)>curdepth) then
            write(*,*) "tracer released below bottom, pos = ", 
     +           newpos
            call write_emission_box(emit_box)
            stop
         endif
            
c........assign acceptable position newpos           
         parbuf(ip)%position     = newpos ! assign validated position
            
c........initialize other spatial attributes  
         parbuf(ip)%mobility     = 0       ! vector assignment
         parbuf(ip)%ashore       = .FALSE. ! check this ...
         parbuf(ip)%outofdomain  = .FALSE. ! check this ...
         parbuf(ip)%atbottom     = .FALSE. ! check this ...
         parbuf(ip)%atsurface    = .FALSE. ! check this ...
         parbuf(ip)%shoreBC      = BC_reflect ! customize ...
         parbuf(ip)%domainBC     = BC_sticky  ! customize ...
         parbuf(ip)%bottomBC     = BC_reflect ! customize ...
         parbuf(ip)%surfaceBC    = BC_reflect ! customize ...

      enddo

      emit_box%current_emissions = emit_box%current_emissions+npar

      if (verbose>0) then
         write(*,*) "activate_emission_box:", npar, 
     +           "new tracers emitted"
      endif

      end subroutine activate_emission_box


      subroutine write_emission_box(emit_box)
c---------------------------------------------------  
      type(emission_box), intent(in) :: emit_box
c---------------------------------------------------       
      call write_time_period(emit_box%emission_interval)
      write(*,*) "emission_boxID      =",emit_box%emission_boxID 
      write(*,*) "SW corner (x,y,z)   =",emit_box%SW
      write(*,*) "NE corner (x,y,z)   =",emit_box%NE
      write(*,*) "max_emissions       =",emit_box%max_emissions 
      write(*,*) "current_emissions   =",emit_box%current_emissions
      write(*,*) "initialization_data =", 
     +            trim(adjustl(emit_box%initialization_data))
      end subroutine write_emission_box


      subroutine get_current_emissions(emit_box, nemit, nemax)
c------------------------------------------------------------- 
      type(emission_box), intent(in) :: emit_box
      integer, intent(out)           :: nemit, nemax
c------------------------------------------------------------- 
      nemit = emit_box%current_emissions
      nemax = emit_box%max_emissions
      end subroutine get_current_emissions


      subroutine get_emission_boxID(emit_box, boxID)
c------------------------------------------------------------- 
      type(emission_box), intent(in) :: emit_box
      integer, intent(out)           :: boxID
c-------------------------------------------------------------
      boxID = emit_box%emission_boxID
      end subroutine get_emission_boxID


      subroutine get_initialization_data(emit_box, strdata)
c------------------------------------------------------------- 
      type(emission_box), intent(in) :: emit_box
      character*(*), intent(out)     :: strdata
c------------------------------------------------------------- 
      strdata = emit_box%initialization_data
      end subroutine get_initialization_data

     

      Subroutine get_xyboxes_emissions(emit_box, xy_sw, xy_ne)
c------------------------------------------------------------- 
      type(emission_box), intent(in)   :: emit_box
      real, dimension(3), intent(out)  :: xy_sw, xy_ne
c------------------------------------------------------------- 
      xy_sw = emit_box%SW
      xy_ne = emit_box%NE
      END Subroutine get_xyboxes_emissions
      
      
      subroutine get_prop_space(space,var,bucket,status)
c------------------------------------------------------------  
      type(spatial_attributes),intent(in) :: space
      type(variable), intent(inout) :: var
      type(polytype), intent(out) :: bucket
      integer, intent(out) :: status
c------------------------------------------------------------  
      status=0  !Variable is present
      select case (get_name(var))
      case ("lon")
        call construct(bucket,"lon",space%position(1))
      case ("lat")
        call construct(bucket,"lat",space%position(2))
      case("depth")
        call construct(bucket,"depth",space%position(3))
      case ("mobx")
        call construct(bucket,"mobx",space%mobility(1))
      case("moby")
        call construct(bucket,"moby",space%mobility(2))
      case("mobz")
        call construct(bucket,"mobz",space%mobility(3))
      case("ashore")
        call construct(bucket,"ashore",space%ashore)
      case("outofdomain")
        call construct(bucket,"outofdomain",space%outofdomain)
      case("atbottom")
        call construct(bucket,"atbottom",space%atbottom)
      case("atsurface")
        call construct(bucket,"atsurface",space%atsurface)
      case("shoreBC")
        call construct(bucket,"shoreBC",space%shoreBC)
      case("domainBC")
        call construct(bucket,"domainBC",space%domainBC)
      case("bottomBC")
        call construct(bucket,"bottomBC",space%bottomBC)
      case("surfaceBC")
        call construct(bucket,"surfaceBC",space%surfaceBC)
      case default
        status=1   !Cannont find variable name
      end select
      end subroutine

      subroutine get_metadata_space(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in) :: var_name
      type(variable), intent(out) :: var
      integer, intent(out) :: status
c------------------------------------------------------------  
      status=0 !Defaults to variable found
      select case (var_name)
      case ("lat")
        call construct(var,"lat","Latitude","deg N","(f6.2)","real")
      case ("lon")
        call construct(var,"lon","Longitude","deg E","(f6.2)","real")
      case ("depth")
        call construct(var,"depth","Depth positive downwards",
     +                  units="m",fmt="(f6.2)",type="real")
      case ("mobx")
        call construct(var,"mobx","Meridonal mobility",
     +                  units="T/F",fmt="(l1)",type="log")
      case ("moby")
        call construct(var,"moby","Zonal mobility",
     +                  units="T/F",fmt="(l1)",type="log")
      case ("mobz")
        call construct(var,"mobz","Vertical mobility",
     +                  units="T/F",fmt="(l1)",type="log")
      case ("ashore")
        call construct(var,"ashore","Is particle ashore?",
     +                  units="T/F",fmt="(l1)",type="log")
      case("outofdomain")
        call construct(var,"outofdomain","Is particle out of "
     +       // "the domain?", units="T/F",fmt="(l1)",type="log")
      case ("atbottom")
        call construct(var,"atbottom","Is the particle at bottom?",
     +                  units="T/F",fmt="(l1)",type="log")
      case ("atsurface")
        call construct(var,"atsurface","Is particle at the "
     +       //"surface?",units="T/F",fmt="(l1)",type="log")
      case ("shoreBC")
        call construct(var,"shoreBC","shoreBC","-","(i2)","int")
      case ("domainBC")
        call construct(var,"domainBC","domainBC","-","(i2)","int")
      case ("bottomBC")
        call construct(var,"bottomBC","bottomBC","-","(i2)","int")
      case ("surfaceBC")
        call construct(var,"surfaceBC","surfaceBC","-","(i2)","int")
      case default
        status=1  !Cannot find variable
      end select
      end subroutine
      
      
 
      end module 
      
