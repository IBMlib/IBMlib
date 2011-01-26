      program ibmrun
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test template for physical fields interface
c
c     make ibmrun
c     ibmrun task_providers/test/testpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use time_tools
      use run_context
      use physical_fields
     
      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4),istat
      type(clock),target  :: start_time
      real    :: xyz1(3),xyz2(3),geo1(3),geo2(3) 
      real    :: xyz(3),uvw(3),lld(3),dryxyz(3),drylld(3) 
      real    :: temp,salt,depth,turb(3),turb2(3),dturb(3)
      real    :: geohit(3),georef(3),xyzhit(3),xyzref(3)
      logical :: anycross
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
c
      call update_physical_fields()

      write (*,*) "------------ Coast line reflection ------------"
      xyz1=    (/90.99876,97.00860,1.0/) 
      xyz2=    (/91.00571,96.97816,1.0/) 
      call xyz2geo(xyz1,geo1) 
      call xyz2geo(xyz2,geo2) 
      geo1(3) =1
      geo2(3) =2
      write(*,*) "xyz1      : ", xyz1," wet? ",is_wet(geo1)
      write(*,*) "xyz2      : ", xyz2," wet? ",is_wet(geo2)
      call coast_line_intersection(geo1,geo2,anycross,georef,geohit)
      write(*,*) "anycross  : ", anycross
      write(*,*) "georef    : ", georef
      write(*,*) "geohit    : ", geohit
      call geo2xyz(georef,xyzref)
      call geo2xyz(geohit,xyzhit) 
      write(*,*) "xyzref    : ", xyzref
      write(*,*) "xyzhit    : ", xyzhit

      write (*,*) "----------Interpolation ----------"
      xyz=    (/33.75,104.25,10.0/)  !Coastal cell 
      xyz=    (/33.75,104.75,10.0/)  !Central cell (all neighbours wet)
      dryxyz =(/34.75,104.75,10.0/)  !On land (east of xyz)
      dryxyz= (/33.75,103.75,10.0/)  !On land (south of xyz) 
      call xyz2geo(xyz,lld) 
      call xyz2geo(dryxyz,drylld) 
      write (*,*) "xyz   : ",xyz
      write (*,*) "lld   : ",lld,"wet:",is_wet(lld)
      write (*,*) "dryxyz: ",dryxyz
      write (*,*) "drylld: ",drylld,"wet:",is_wet(drylld)
      call interpolate_wdepth(lld,depth,istat)
      write (*,*) "depth : ",istat, depth
      call interpolate_currents(lld,uvw,istat)
      write (*,*) "uvw   : ",istat,uvw
      call interpolate_temp(lld,temp,istat)
      write (*,*) "temp  : ",istat, temp
      call interpolate_salty(lld,salt,istat)
      write (*,*) "salt  : ",istat, salt

      write (*,*) "------------ Turbulence  ------------"
      call interpolate_turbulence(lld,turb,istat)
      write (*,*) "turbulence      : ", istat, turb
      lld(3) = lld(3) + 0.1
      call interpolate_turbulence(lld,turb2,istat)
      write (*,*) "turbulence + dz : ", istat, turb2
      write (*,*) "delta turb      : ", 0,turb2-turb 
      call interpolate_turbulence_deriv(lld,dturb,istat)
      write (*,*) "dturb           : ", istat, dturb


      call close_physical_fields()
      write(*,*) "normal end of simulation"
      end program



