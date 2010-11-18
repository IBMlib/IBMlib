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
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4),istat
      type(clock),target  :: start_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      integer :: year, month, day, julday,last
      integer :: i,ix,iy,iz,dum(79),ilast,iunit
      integer :: idum(4), timeout, istep, ntracers
      real    :: time_step, now,xyz(3)
      real    :: amp,xcen,ycen,vdiff,hdiff,res
      character*256 :: filename
      
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
c
      call update_physical_fields()
      xyz(1:2) = (/6.0,56.0/)
      call interpolate_wdepth(xyz,res,istat)
      write(*,*) res,istat
c
      call close_physical_fields()
      
      write(*,*) "normal end of simulation"

      end program
