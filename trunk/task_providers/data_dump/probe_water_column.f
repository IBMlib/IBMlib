ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Template for averaging key properties of the 
c     water column at a specific position: 
c
c       depth, horizontal current shear, vertical diffusivity
c
c     ignore fluctuating sealevel at averaging
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     make ibmrun
c     ibmrun task_providers/data_dump/simpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
     
      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      real    :: time_step
      integer :: i, idum4(4), istep, start(10), nwords, nzsamp
      integer :: nstations, inext, ist, iz, istat
      character(len=256)              :: cbuf, filename
      character(len=128), allocatable :: tags(:)
      real, allocatable  :: v(:,:,:), u(:,:,:), xypos(:,:)
      real, allocatable  :: kz(:,:,:), depth(:,:), zgrid(:)
      real               :: wd, uvw(3), xyz(3), k3(3)

c     ------------   show time starts  ------------
      call init_run_context()
  
c     ------------   set clocks  ------------   
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "time_step", time_step)

c     ------------   init system  ------------   
      call init_physical_fields(start_time)
      current_time => get_master_clock()
      call update_physical_fields()

c     ------------  read zgrid   ------------  
      call read_control_data(simulation_file,"zgrid",cbuf)
      call tokenize(cbuf, start, nzsamp)
      if (nzsamp<1) then
         write(*,*) "illegal zgrid:", cbuf
         stop 731
      endif
      allocate( zgrid(nzsamp) )
      read(cbuf,*) zgrid
      write(*,*) "zgrid = ", zgrid

c     ------------  read stations  ------------  
      nstations = count_tags(simulation_file, "xystation") 
      if (nstations<1) then
         write(*,*) "number of xystation < 1:", nstations
         stop 754
      endif

      allocate( u(2,nzsamp,nstations)  )
      allocate( v(2,nzsamp,nstations)  )
      allocate( kz(2,nzsamp,nstations) )
      allocate( depth(2,nstations)     )
      allocate( xypos(2,nstations)     )
      allocate( tags(nstations)        )
      
      inext = 1
      do ist = 1, nstations
         call read_control_data(simulation_file,"xystation",cbuf,inext)
         call tokenize(cbuf, start, nwords)
         read(cbuf(start(1):start(2)-1), *) tags(ist)
         read(cbuf(start(2):), *)           xypos(:,ist)
         write(*,39) ist, nstations, trim(tags(ist)), xypos(:,ist)
         inext = inext + 1
      enddo
 39   format("read station", i4, "/", i4, ": ", a, "  xy=", 2f12.7)

      u     = 0.0
      v     = 0.0
      kz    = 0.0
      depth = 0.0
      
c     =====================  main time loop =====================
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
         call update_physical_fields()
         do ist = 1, nstations
            xyz(1:2) = xypos(:,ist)
            call interpolate_wdepth(xyz,wd,istat)
            if (istat /= 0) then
               write(*,*) "interpolate_wdepth: istat = ", istat
               stop 634
            endif
            depth(1,ist) = depth(1,ist) + wd  
            depth(2,ist) = depth(2,ist) + wd**2  
            do iz = 1, nzsamp
               xyz(3) = wd*zgrid(iz)
               call interpolate_currents(xyz,uvw,istat)
               call interpolate_turbulence(xyz,k3,istat)
               u(1,iz,ist)  = u(1,iz,ist)  + uvw(1)
               v(1,iz,ist)  = v(1,iz,ist)  + uvw(2)
               kz(1,iz,ist) = kz(1,iz,ist) + k3(3)
               u(2,iz,ist)  = u(2,iz,ist)  + uvw(1)**2  
               v(2,iz,ist)  = v(2,iz,ist)  + uvw(2)**2  
               kz(2,iz,ist) = kz(2,iz,ist) + k3(3)**2  
            enddo
         enddo
c        --- prepare for next step ---
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1
      enddo
      
      if (istep == 0) then
         write(*,*) "no steps -> no results"
         stop 461
      else
         u          = u/istep        ! average
         v          = v/istep        ! average
         kz         = kz/istep       ! average
         depth      = depth/istep    ! average
         u(2,:,:)   = sqrt(u(2,:,:)   - u(1,:,:)**2)   ! rms
         v(2,:,:)   = sqrt(v(2,:,:)   - v(1,:,:)**2)   ! rms
         kz(2,:,:)  = sqrt(kz(2,:,:)  - kz(1,:,:)**2)  ! rms
         depth(2,:) = sqrt(depth(2,:) - depth(1,:)**2) ! rms
      endif

 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))


c     -------- write results --------
      do ist = 1, nstations
         write(filename, 66) trim(tags(ist))
         open(47, file=filename) 
         write(47, 67) xypos(:,ist)
         write(47, 68) depth(:,ist)
         write(47, 69) 
         do iz = 1, nzsamp
            write(47, 70) zgrid(iz), 
     +                    u(:,iz,ist), v(:,iz,ist), kz(:,iz,ist) 
         enddo
         close(47)
      enddo
 66   format(a,".dat")
 67   format("# station lon,lat = ", 2f12.7)
 68   format("# depth   avg,rms = ", 2f12.7)
 69   format("# z  u_avg u_rms  v_avg v_rms kz_avg kz_rms")
 70   format(f8.4, 3(2x,2f12.7))
c
      write(*,*) "normal end of simulation"

      end program
