ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Minimal example: oxygen extraction at specific times+positions
c     as provided in a file with a hard coded name
c
c 1mL dissolved O2 /L is about 1.33 mg dissolved O2 /L.
c 1 mmol O2 = 32 mg O2
c 1mL/L = 1.33 mg/L = 1.33/32 mmol/L = 1330/32 mmol/m3 = 41.5625 mmol/m3     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use geometry
      use physical_fields
      
     
      implicit none
c     ------------ declarations ------------ 
      type(clock),pointer  :: current_time
      integer       :: idum4(4),istat
      integer       :: statnum
      real          :: geo(3),wdepth, obsoxy, obott
      integer       :: iland
     
c     ------------   show time starts  ------------
      call init_run_context()
      call init_physical_fields()
      current_time => get_master_clock()
      open(34, file="stations_2005.dat")
      statnum = 0
      iland   = 0
      do
         read(34,*,end=555) geo(1:2), idum4, obsoxy
         statnum = statnum + 1  ! process station
         call set_clock(current_time,
     +        idum4(1), idum4(2), idum4(3), idum4(4))
         call update_physical_fields()
         if (.not.is_land(geo)) then
            call interpolate_wdepth(geo,wdepth,istat)
            geo(3) = max(1e-3, wdepth - 0.1)
            call interpolate_oxygen(geo,obott,istat)
            write(66,*) obsoxy, obott/41.5625 ! convert to mL/L
         else ! ignore dry points
            iland = iland + 1
         endif
      enddo                     ! station loop
c      
 555  close(34)
      close(66)
      call close_physical_fields()
      write(*,*) "processed ", statnum, " stations "
      write(*,*) iland, " dry stations were recast"
      write(*,*) "normal end of run"
    
 
      contains


      end program

      
