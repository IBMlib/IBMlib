ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Template for generating hydrographic indices based on some 
c     spatio-temporal averaging
c     
c     This example loads at set of sampling points and generates
c     a seasonal average of water temperature (averaged first vertically
c     and then seasonally+spatially over loaded grid) as well as 
c     vertical zooplankton maximum, averaged seasonally+spatially over loaded grid.
c
c     Notes: 
c        This example assumes points which are wet at start time remains wet
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     make ibmrun
c     ibmrun task_providers/data_dump/simpar_seasonal_indices
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
     
      implicit none
c     ------------ declarations ------------      
      type(clock),pointer  :: current_time
      type(clock)          :: start_time
      integer :: year_start, year_end, this_year
      integer :: day_start, month_start, ndays
      integer :: nspts, iday, ispt, nvpts, ivpt
      integer :: istat, nall, ndry, iunit, i
      real    :: vres, xy(2), wd, zvec(5), t
      real    :: zmax_here, temp_here
      real    :: tavg, trms, zavg, zrms         
      real, allocatable :: xyz(:,:), zmax(:), temp(:)
      character*999     :: fname
c     ------------   show time starts  ------------
      call init_run_context()
      
      call read_control_data(simulation_file, "year_start", year_start)
      call read_control_data(simulation_file, "year_end",   year_end)
      call read_control_data(simulation_file, "day_start",  day_start)
      call read_control_data(simulation_file, "month_start",month_start)
      call read_control_data(simulation_file, "days_in_period", ndays)
      call read_control_data(simulation_file,"vertical_resolution",vres)
c
c     ------------ read spatial sampling set ------------
c         set clock and update_physical_fields to ensure potential dynamic
c         topography is set for dry point detection
c
      call set_clock(start_time,year_start,month_start,day_start,0)
      call init_physical_fields(start_time)
      current_time => get_master_clock()
      call update_physical_fields()
c
      call read_control_data(simulation_file,"sampling_points", fname)
      call find_free_IO_unit(iunit) 
      open(iunit, file=fname, status='old')
c     ------- first count number of points in sampling file -------
      nspts = 0
      nall  = 0
      do 
         read(iunit,*,end=732) xy  
         nall = nall + 1 ! wet and dry
         if (is_land(xy)) ndry = ndry + 1
      enddo
 732  write(*,*) "read sampling points from file: ",trim(adjustl(fname))
      nspts = nall - ndry
      write(*,*) "total number of points in file: ",nall
      write(*,*) "number of points on land:       ",ndry
      write(*,*) "number of wet points:           ",nspts
c     ------- read only dry points, ignore rest -------
      allocate( xyz(3,nspts) )
      rewind(iunit)
      i = 1
      do 
         read(iunit,*,end=737) xy  
         if (is_land(xy)) cycle ! skip this point
         xyz(1:2,i) = xy        ! do not set element 3 (verical coordinate)
         i = i + 1
      enddo
 737  close(iunit)
c
      call read_control_data(simulation_file,"outputfile", fname)
      open(iunit, file=fname) ! unit still free
c
c     ------------ year loop ------------
c      
      allocate( zmax(nspts*ndays) )
      allocate( temp(nspts*ndays) )
      do this_year =  year_start, year_end
         call set_clock(current_time,this_year,month_start,day_start,0) ! midnight - generalize later
         zmax = 0.0
         temp = 0.0
         i    = 1 ! array counter for point in zmax/temp
c        ------ day loop ------
         do iday = 1, ndays ! loop over days in this season
            call update_physical_fields()
c           ------ horizontal loop ------
            do ispt = 1, nspts ! loop over space points in sampling set
               call interpolate_wdepth(xyz(:,ispt),wd,istat)
               nvpts = 1 + int(wd/vres) ! nvpts >= 1
               zmax_here =  0.0 ! zmax >= 0
               temp_here =  0.0
c              ------ vertical loop ------
               do ivpt = 1, nvpts
                  xyz(3,ispt) = (ivpt-1)*vres
                  call interpolate_temp(xyz(:,ispt),t,istat)
                  call interpolate_zooplankton(xyz(:,ispt),zvec,istat)
                  zmax_here = max(zmax_here, zvec(1)) 
                  temp_here = temp_here + t ! accumulate average
               enddo ! ivpt (vertical loop)
               zmax(i) = zmax_here 
               temp(i) = temp_here/nvpts  ! nvpts >= 1
               i       = i+1 ! array counter for zmax/temp
            enddo  !  ispt (horizontal loop)
            call add_seconds_to_clock(current_time, 86400)
         enddo ! iday (day loop)
c        do statistics on sampled data (zmax,temp)
         zavg = sum(zmax)/nspts/ndays
         zrms = sqrt(sum((zmax-zavg)**2)/nspts/ndays)
         tavg = sum(temp)/nspts/ndays
         trms = sqrt(sum((temp-tavg)**2)/nspts/ndays)
         write(iunit, 275) this_year, tavg, trms, zavg, zrms         
      enddo ! this_year (year loop)
 275  format(i5, 3x, 2f12.7, 3x,2e15.7)
      deallocate(zmax)
      deallocate(temp)
      deallocate(xyz)
      close(iunit)
      write(*,*) "normal end of simulation"

      end program
