      program tracker
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Forward/backward particle ensemble tracking simulation
c
c     Metrics theory : for a 1D RW ensemble subject to homogeneous 
c     spatial diffusivity and no advection, the statistical properties 
c     are as follows:
c
c       n = ensemble size
c       D = diffusivity
c       T = dispersal period
c       x = real(n) = position of ensemble members
c
c     Estimators on ensemble center (c) and variance (v9:
c       ec = sum(x)/n
c       ev = sum((x-ec)**2)/n
c
c     Statistical properties of estimators (noise from
c     finite ensemble size):
c       exp(ec) = 0
c       exp(ev) = 2*D*T
c       RMS(ec) = sqrt(2.0*D*T/n)
c       RMS(ev) = sqrt(8.0/n)*D*T
c
c     The performance of the ensemble simulation, P, is tested
c     by comparing to the expected value of the estimator on
c     the scale of the estimator RMS:
c       P(c) = (ec-exp(ec))/RMS(ec)
c       P(v) = (ev-exp(ev))/RMS(ev)
c     
c     Log: Startup + R script analysis + plotting: MPA
c          metrics + bugfixes ASC: Jan 25, 2011
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use particle_tracking
      use geometry
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: release_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer      :: idum4(4), istep,last,i, ipar,iunit,dir,nsteps
      integer      :: summary,tsim,istat,tsign
      real         :: xyz(3),time_step, xyz0(3),xyz_sw(3), xyz_ne(3)
      character*999 :: output_files(2)
      character*20 :: direcname(2) 
      real         :: ec(3), exp_ec(3), diff_ec(3), RMS_ec(3), P_c(3)
      real         :: ev(3), exp_ev(3), diff_ev(3), RMS_ev(3), P_v(3)
      real         :: hdiff,vdiff,lincur(3),dxyz(3)

c     ------------   show time starts  ------------
      call init_run_context()
      
c     ------------ initialise system --------------
      call random_seed()
      call init_physical_fields()
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes

c     ------------ read parameters from file --------------
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
      call read_control_data(simulation_file,"backward_filename",
     +    output_files(1))
      call read_control_data(simulation_file,"forward_filename",
     +    output_files(2))
      call read_control_data(simulation_file, "emitbox", idum4)
      call set_clock(release_time,idum4(1),idum4(2),idum4(3),idum4(4))
      call read_control_data(simulation_file,"nsteps",nsteps)
      
      call find_free_io_unit(summary)
      open(unit=summary,file="test_summary")
      tsim = nint(abs(time_step*nsteps)) ! seconds

      call read_control_data(simulation_file,"hdiffus",hdiff)
      call read_control_data(simulation_file,"vdiffus",vdiff)
 
      direcname(1) = "backward"
      direcname(2) = "forward"

c     ------------ loop first backwards, then forwards --------------
      do dir=1,2 
         tsign = (-1)**dir
c        Setup ensemble      
         call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
         call get_xyboxes_emissions(emitboxes(1), xyz_sw, xyz_ne) ! assume only one
         xyz0 = 0.5*(xyz_sw + xyz_ne) ! center position

c        Setup output file
         call find_free_io_unit(iunit)
         open(unit=iunit,file=trim(output_files(dir)))
      
c        Set end time and run direction
         time_step = abs(time_step)*tsign
         call set_master_clock(release_time)
         call set_clock(end_time,release_time) ! assign end_time = release_time
         call add_seconds_to_clock(end_time,nint(time_step*nsteps)) ! incl sign
         current_time => get_master_clock()
      
c        =====================  main time loop =====================
         istep   = 0
         do while (compare_clocks(current_time,end_time)*tsign <= 0)
            write(*,372) dir, istep
            write(*,*) get_datetime(current_time)
            call update_physical_fields()
   
c           -------- release new tracers if needed  --------
            call generate_particles(par_ens, emitboxes, time_step)
            call get_last_particle_number(par_ens, last) 

c           -------- dump particle positions --------
            do i=1,last
               call get_particle_position(get_particle(par_ens,i),xyz)
               write(iunit,*) istep*time_step, xyz
            enddo  

c           -------- propagate tracers  --------
            call update_particles(par_ens, time_step)

c           -------- loop control       --------
            call add_seconds_to_clock(current_time, nint(time_step)) ! incl sign
            istep = istep + 1
         enddo ! while compare_clocks
         close(iunit)
c
c        ========== analyze final horizontal distribution ==========
c
c        --- evaluate ensemble properties ---
c
         ec = 0.  ! accumulate patch centroid
         ev = 0.  ! accumulate patch variance
         do i=1,last
            call get_particle_position(get_particle(par_ens,i),xyz)
            ec = ec + xyz
            ev = ev + xyz**2
         enddo  
         ec          = ec/last
         ev          = ev/last - ec**2
c
c        --- evaluate expected ensemble properties ---
c         
         exp_ec      = xyz0 
         call interpolate_currents(xyz0,lincur,istat)    ! lookup linear current vector as Cartesian
         call add_finite_step(exp_ec, lincur*tsim*tsign) ! inplace addition       
c
         exp_ev(1:2) = 2*hdiff*tsim    ! in m2
         exp_ev(3)   = 0               ! elaborate later
c
c        --- evaluate differences in Cartesian ---
c 
         diff_ec     = ec-exp_ec       ! in (lon, lat, depth)
         call dxy2dcart(ec, diff_ec)   ! now in (m,m,m)
         diff_ev     = ev              ! in (lon, lat, depth)**2 
         call dxy2dcart(ec, diff_ev) 
         call dxy2dcart(ec, diff_ev)   ! now in (m2,m2,m2)
         diff_ev     = diff_ev-exp_ev  ! in (m2,m2,m2)
c
c        --- evaluate estimator variances ---
c   
         RMS_ec(1:2) = sqrt(2.0*hdiff*tsim/last)   ! in m 
         RMS_ec(3)   = 1   ! elaborate later       ! in m 
         RMS_ev(1:2) = sqrt(8.0/last)*hdiff*tsim   ! in m2 
         RMS_ev(3)   = 1   ! elaborate later       ! in m2 
c
c        --- evaluate performance indices ---
c   
         P_c(1:2)    = diff_ec(1:2)/RMS_ec(1:2) ! elaborate later for vertical
         P_v(1:2)    = diff_ev(1:2)/RMS_ev(1:2) ! elaborate later for vertical 
         P_c(3)      = 0 ! elaborate later for vertical
         P_v(3)      = 0 ! elaborate later for vertical
      
         write(summary,811) adjustl(direcname(dir)), "centroid",  P_c
         write(summary,811) adjustl(direcname(dir)), "variance",  P_v                 
      enddo ! dir loop 
      
 811  format("direction: ",a8," patch ",a," prediction =",
     +        3f12.4, " : test OK")
 372  format(25("*"), "loop ",i2, ", step ",i5.5, " ", 25("*"))
 
c     ----------------- close down ---------------------------
     
      call close_physical_fields()
      call close_particles()
      
      write(*,*) "normal end of simulation"

      end program
