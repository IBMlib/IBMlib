      program round_trip
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     1) Test whether vertical dispersion is correct after equilibration 
c        (i.e. no vertical aggregation)
c     2) Test whether horizontal trajectory integration is correct in a gyre 
c     
c     
c     If the water depth is denoted by w, the number
c     of particles between 0 < z < d < w is given by
c     the binomial distribution B(n,p), with p = d/w and 
c        exp(m) = n*p
c        RMS(m) = sqrt( n*p*(1-p) )
c     if the vertical distribution of particles is random for 0 < z < w 
c     The performance test is then
c        P(m) = (m-exp(m))/RMS(m)
c     where m is the actual count
c  
c     Premises for vertical equilibration is that the 
c     simulation time T satisfies
c        w << sqrt(2*D*T)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use constants
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
      integer      :: summary,istat,nsurf,nbott
      real         :: geo_sw(3), geo_ne(3),avgpos(3),tsim
      real         :: vortexcen(2),dist,cspeed,difflen,wd
      real         :: equilib,time_step,errdist,traj_err
      real         :: hdiff,vdiff,lincur(3),dxyz(3),p,rms
      real         :: nexpct,noagg_surf,noagg_bott,xy(3)
      real,allocatable,target :: finpos(:,:)
      real,pointer            :: z(:)
      logical      :: equilib_ok
      character*10 :: traj_err_res, noagg_surf_res, noagg_bott_res

c     --- set criteria for accepting performance tests ---
      real, parameter :: equilib_limit  = 0.1  ! equilib << 1
      real, parameter :: traj_err_limit = 0.01 ! traj_err_limit << 1
      real, parameter :: noagg_limit    = 3    ! expect noagg ~ 1

c     ------------   show time starts  ------------
      call init_run_context()
      
c     ------------ initialise system --------------
      call init_physical_fields()
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes

c     ------------ read parameters from file --------------
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
  
      call read_control_data(simulation_file, "emitbox", idum4)
      call set_clock(release_time,idum4(1),idum4(2),idum4(3),idum4(4))
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      call set_master_clock(release_time)     
      current_time => get_master_clock()
c      
c     Compute simulation time tsim so that it corresponds to a 
c     round trip in the gyre
c
      call read_control_data(simulation_file,"lonvc",vortexcen(1)) 
      call read_control_data(simulation_file,"latvc",vortexcen(2))
      call get_xyboxes_emissions(emitboxes(1), geo_sw, geo_ne) ! assume only one
      call get_horizontal_distance(vortexcen,geo_sw,dist)
      call read_control_data(simulation_file,"cspeed",cspeed)
      write(*,*) "input time step         = ", time_step,  "sec"
      tsim      = 2*pi*dist/cspeed
      nsteps    = nint(tsim/time_step)
      time_step = tsim/nsteps ! adjust time step to make exact round trip
      write(*,*) "adjusted time step      = ", time_step,  "sec"
      write(*,*) "simulation period       = ", tsim/3600., "hours"
      write(*,*) "distance to gyre center = ", dist/1000., "km"
c
c     =====================  main time loop =====================
      do istep = 1,nsteps
         write(*,372) istep
         call update_physical_fields() ! static field
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)
      enddo
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))
                
c
c     ========== analyze final result ==========
c
      write(*,*) " ********** begin analysis ********** "

c     
c     1) analyze final horizontal distribution
c        particles should make exactly one round trip in the gyre
c        Error measure traj_err should satisfy |traj_err| << 1
c
      call get_last_particle_number(par_ens, last) 
      if (last <= 0) stop "round_trip: no particles - bailling out"
      allocate( finpos(3,last) )
      do i=1,last
         call get_particle_position(get_particle(par_ens,i),finpos(:,i))   
      enddo
      z      => finpos(3,:)
      avgpos = sum(finpos,dim=2)/last 
      call get_horizontal_distance(geo_sw, avgpos, errdist)
      traj_err = errdist/dist
     
      write(*,*)
      write(*,*) "--- horizontal closure check: ---"
      write(*,*) "initial   pos (deg E, deg N) =", geo_sw(1:2)
      write(*,*) "avg final pos (deg E, deg N) =", avgpos(1:2)
      write(*,*) "trajectory closure error (deg) =", traj_err,"<< 1 ?"
c      
      if (traj_err < traj_err_limit) then
          traj_err_res = "[OK]"
      else
          traj_err_res = "[FAILED]"
      endif
c
c     2) check vertical non aggregation
c        first check assertion wd << sqrt(2*vdiff*tsim) 
c
c        (equilibration progress) is fullfilled
c     
      call interpolate_wdepth(avgpos, wd, istat)
      call read_control_data(simulation_file,"vdiff",vdiff)
      difflen = sqrt(2*vdiff*tsim)
      equilib = wd/difflen

      write(*,*)
      write(*,*) "--- vertical non aggregation check: ---"
      write(*,*) "water depth            =", wd
      write(*,*) "diffusion length       =", difflen
      write(*,*) "equilibration progress =", equilib, "<< 1 ?"
      equilib_ok = (equilib < equilib_limit)   
c
c     the width of the surf/bott tesing region is controleld by setting the fraction p
c
      p     = 0.1
      nsurf = count(z <     p*wd) ! number of particles in slice in surf region
      nbott = count(z > (1-p)*wd) ! number of particles in slice in bott region      
      rms        = sqrt( last*p*(1-p) ) ! expected RMS on nsurf and nbott
      nexpct     = last*p               ! expected number of particles in region
c
      noagg_surf = (nsurf - nexpct)/rms ! |aggregation index| ~ 1 => no aggregation
      if (noagg_surf < noagg_limit) then
         noagg_surf_res = "[OK]"
      else
         noagg_surf_res = "[FAILED]"
      endif
c
      noagg_bott = (nbott - nexpct)/rms ! |aggregation index| ~ 1 => no aggregation
      if (noagg_bott < noagg_limit) then
         noagg_bott_res = "[OK]"
      else
         noagg_bott_res  = "[FAILED]"
      endif       
c
c     3) Finally report result of test to file "test_summary"
c
      call find_free_io_unit(summary)
      open(unit=summary, file="test_summary")
      write(summary,*) "result vertical_homogenity:"
      write(summary,568) "horizontal trajectory closure error =", 
     +                   traj_err, traj_err_res
      if (equilib_ok) then
         write(summary,568) "vertical aggregation at surface   =", 
     +                   noagg_surf, noagg_surf_res
         write(summary,568) "vertical aggregation at bottum    =", 
     +                   noagg_bott, noagg_bott_res
      else
         write(summary,*) "   vertical distribution not at equilibrium"
      endif
      close(summary)
 568  format(3x,a,1x,f9.4,1x,a)
      
c     ----------------- close down ---------------------------
     
      call close_physical_fields()
      call close_particles()
      
      write(*,*)
      write(*,*) "normal end of simulation"


      end program
