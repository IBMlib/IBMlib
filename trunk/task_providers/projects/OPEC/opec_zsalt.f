ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------------------
c     Generate a time series of depths corresponding to a specific salinity
c     at a specific position
c     -----------------------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Extraction for Mikael van Deurs
c
c     make ibmrun
c     ibmrun task_providers/projects/OPEC/simpar_zsalt
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use constants

      implicit none
c     ------------ declarations ------------      
      type(clock),pointer  :: current_time
      integer          :: start_year
      integer          :: end_year   ! incl
      integer          :: time_step  ! seconds
      character(len=999) :: ofname      
      real    :: xy(2), zpos, salttarget, t
      integer :: i, iunit, istat, year, steps_in_year
      integer :: nfaults
      real, parameter    :: zacc = 1e-3 ! meters
      integer, parameter :: seconds_in_a_year = 365*86400
c     ------------   show time starts  ------------
c
      call init_run_context()
c      
c     ------------ set clocks  ------------
c
      current_time => get_master_clock()  ! ... before init_physical_fields, hmm
      call read_control_data(simulation_file, "start_year", start_year)
      write(*,*) "start year = ", start_year
      call read_control_data(simulation_file, "end_year", end_year)
      write(*,*) "end year   = ", end_year
      call read_control_data(simulation_file, "time_step", time_step)
      write(*,*) "sampling time step = ", time_step, "secs"
      call read_control_data(simulation_file, "position", xy)
      write(*,*) "sampling position  = ", xy, "degE,degN"
      call read_control_data(simulation_file, "salinity_target", 
     +                       salttarget)
      write(*,*) "salinity_target    = ", salttarget
      call read_control_data(simulation_file, "output_file", ofname)
      call find_free_IO_unit(iunit) 
      open(iunit, file=ofname)
      write(*,*) "writing output to new file:",trim(adjustl(ofname))     
c
c      
      call set_clock(current_time, start_year, 1, 1, 0)
      call init_physical_fields()   
c
c     ------------ time loop ------------
c   
      nfaults = 0
      steps_in_year = seconds_in_a_year/time_step ! integer division
      do year = start_year, end_year
         call set_clock(current_time, year, 1, 1, 0)
         do i = 1, steps_in_year
            t = year + 1.0*(i-1)*time_step/seconds_in_a_year 
            call update_physical_fields()   
            call locate_salinity_depth(xy, salttarget, zacc, 
     +                                 zpos, istat)
            if (istat == 0) then
               write(iunit, 441) t, zpos
            else
               write(*,442) t, istat
               nfaults = nfaults + 1
            endif
            call add_seconds_to_clock(current_time, time_step)
         enddo   ! i = 
      enddo      ! year = 
 441  format(f12.4, 1x, f12.7)
 442  format("fault at t = ", f12.4, 1x, " istat = ", i2)

      close(iunit)
      call close_run_context()
      call close_physical_fields()
      write(*,*) "normal end data extraction with", nfaults, "faults"

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine locate_salinity_depth(xy, salttarget, zacc, 
     +                                 zpos, status)
c     -------------------------------------------------------------------
c     Locate at position xy = (lon,lat) the water depth zpos where
c     the salinity is salttarget, to within accuracy zacc
c     It must be satisfied that 
c        salinity(surface) < salttarget < salinity(bottom) 
c     otherwise no determination of zpos is attempted. I rare situations
c     this may miss a real solution, since water density also depends on temperature
c     Do not assess stability of the water column.
c     Currently just apply bisection because this is not a bottleneck
c     
c  
c     Return depth zpos and a status flag evaluating special conditions:
c         status =  0: successful determination of zpos
c         status =  1: salttarget < salinity(surface), return zpos = -1.0   
c         status =  2: salttarget > salinity(bottom) , return zpos = -1.0   
c         status =  3: xy is a dry position,           return zpos = -1.0   
c     -------------------------------------------------------------------
      real, intent(in)     :: xy(:)      ! (lon,lat)
      real, intent(in)     :: salttarget ! solution target where salinity(zpos) == salttarget
      real, intent(in)     :: zacc       ! absolute accuracy on solution zpos (if it exist)
      real, intent(out)    :: zpos       ! if status == 0, 0 <= zpos <= CWD, and salinity(zpos) == salttarget
      integer, intent(out) :: status     ! see above
c
      real       :: cwd, salts, saltb, saltm, geo(3)
      real       :: z0, z1, salt0, salt1
      integer    :: istat, nsteps, ist
c     -------------------------------------------------------------------
      geo(1:2) = xy(1:2) 
      zpos     = -1.0    ! set exception value

      if (is_land(xy)) then
         status =  3     ! signal dry position
         return
      endif   
      call interpolate_wdepth(geo, cwd, istat) 
      geo(3) = 0
      call interpolate_salty(geo, salts, istat) ! surface salinity
      geo(3) = cwd
      call interpolate_salty(geo, saltb, istat) ! bottum salinity

c     --- check assertion salinity(surface) < salttarget < salinity(bottom) 
      if (salttarget < salts) then
         status = 1
c     --- check for sinking to bottom (but stable water column)
      elseif (salttarget > saltb) then
         status = 2
      endif    
c
c     Now it is established that: salts < salttarget < saltb and xy is a wet point
c     Solve salinity(zpos) == salttarget by bisection
c     
      nsteps = 1 + int(log(cwd/zacc)/log(2.0))
      z0 =   0.0
      z1   = cwd
      salt0 = salts
      salt1 = saltb
      do ist = 1, nsteps
         geo(3) = 0.5*(z0+z1)
         call interpolate_salty(geo, saltm, istat) 
         if (saltm < salttarget) then ! advance left bracket
            z0   = geo(3)
            salt0 = saltm
         else                         ! decrease right bracket
            z1   = geo(3)
            salt1 = saltm
         endif
      enddo                     ! ist
c
      status = 0                ! signal normal completion
      zpos   = geo(3)
        
      end subroutine locate_salinity_depth
c
      end program
