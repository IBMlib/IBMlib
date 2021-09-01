ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Test template for physical fields interface
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     make ibmrun
c     ibmrun task_providers/test/testpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use constants
      use write_netcdf_lonlat_data 
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4),istat
      type(clock),target  :: start_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday,last,n,nxs,nys
      integer :: i,ix,iy,iz,dum(79),ilast,iunit
      integer :: idum(4), timeout, istep, ntracers
      real    :: time_step, now,geo(3), geo0(3)
      real    :: amp,xcen,ycen,vdiff,hdiff,res,wd
      real    :: x,y,z,dxyz_dr(3,3), xyz0(3), xyz(3)
      real    :: dxyz_dr_num(3,3), carstep,r3(3)
      real    :: x1,y1,dx,dy
      real, allocatable :: f(:,:,:)
      character*256 :: filename
      real, parameter :: dd = 0.00001
      type(netcdf_grid_2Ddata_handler) :: nch

c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
c
      call update_physical_fields()
      geo(1:2) = (/5.41,56.08/)
      call interpolate_wdepth(geo,wd,istat)
      write(*,*) wd,istat
c
c      n = 1000
c      do i=1,n
c         geo(3) = wd*i/n
c         call interpolate_turbulence(geo, r3, istat)
c         write(18,*) geo(3),r3
c     enddo

      nxs = 100
      nys = 100
      x1 = -5.
      y1 = 50.
      dx = 0.15
      dy = 0.1
      allocate( f(3,nxs,nys) )
      call init_netcdf_dataset(nch, "jj.nc",  
     +     nxs,nys,x1,y1,dx,dy, 'xyt', "")
      geo(3) = 0.5
      do ix=1,nxs
         geo(1) = x1 + (ix-1)*dx
         do iy=1,nys
            geo(2) = y1 + (iy-1)*dy
            call interpolate_turbulence(geo, f(:,ix,iy), istat)
         enddo
      enddo
      call write_netcdf_data(nch,f(1,:,:),timeval=0)
      call close_netcdf_dataset(nch)

      call close_physical_fields()
      write(*,*) "normal end of simulation"

      end program
