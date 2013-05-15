ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------------
c     Core generic physics/biogeochemistry linear data extraction 
c     utility that can be used for external scripting
c     Averages phisical/biogeochemical property over a spatio-temporal range
c     Includes only wet points in statistics
c     The preprocessor flag WITH_BGC (test: #ifdef WITH_BGC) controls whether 
c     the biogeochemical interface is activated in addition to the physics-only interface
c     Currently WITH_BGC must be set locally
c
c     Report:
c        average(property),  RMS(property), minval(property), maxval(property)
c
c     over spatio-temporal range to file name given to tag output_file (or append_file). 
c     See comments in input file simpar_linear_extraction
c     for documentation on supported options + units
c     -----------------------------------------------------------------
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     make ibmrun
c     ibmrun task_providers/data_dump/simpar_seasonal_indices
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#define WITH_BGC
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
     
      implicit none
c     ------------ declarations ------------      
      type(clock),pointer  :: current_time
      type(clock)          :: start_time
      type(clock)          :: end_time
      integer              :: time_step
      real    :: lonmin, lonmax, lon_step
      real    :: latmin, latmax, lat_step
      real    :: depthmin, depthmax, ver_step
      logical :: at_bottom, at_surface, single_stratum
      character(len=99)  :: layer_tag, property, cbuf
      character(len=999) :: ofname
      character(len=4)   :: icesrec
      real    :: xyz(3), cwd, stop_depth, r3(3)
      real    :: p, pavg, pmin, pmax, pavg2, prms
      integer :: idum4(4), ic, iunit, istat
      integer :: nstatpt, nwetpt, ntimpt
c      
c     property/key word map: should be moved to template form and preprocessed
c
      integer            :: property_selector = -1
      integer, parameter :: select_temp        = 1
      integer, parameter :: select_salinity    = 2
      integer, parameter :: select_current     = 3
      integer, parameter :: select_windstress  = 4
#ifdef WITH_BGC
      integer, parameter :: select_zooplankton      = 10
      integer, parameter :: select_oxygen           = 11
      integer, parameter :: select_nh4              = 12
      integer, parameter :: select_no3              = 13
      integer, parameter :: select_po4              = 14
      integer, parameter :: select_diatoms          = 15           
      integer, parameter :: select_flagellates      = 16             
      integer, parameter :: select_cyanobacteria    = 17          
      integer, parameter :: select_organic_detritus = 18                 
      integer, parameter :: select_part_org_matter  = 19
      integer, parameter :: select_DIC              = 20
      integer, parameter :: select_alkalinity       = 21
      integer, parameter :: select_DIN              = 22
      integer, parameter :: select_chlorophyl       = 23
#endif
c
      character*(*), parameter :: temp_kword             = "temperature"
      character*(*), parameter :: salinity_kword         = "salinity"
      character*(*), parameter :: current_kword          = "currents"
      character*(*), parameter :: windstress_kword       = "windstress"
#ifdef WITH_BGC
      character*(*), parameter :: zooplankton_kword      = "zooplankton"
      character*(*), parameter :: oxygen_kword           = "oxy"
      character*(*), parameter :: nh4_kword              = "nh4"
      character*(*), parameter :: no3_kword              = "no3"
      character*(*), parameter :: po4_kword              = "po4"
      character*(*), parameter :: diatoms_kword          = "dia"         
      character*(*), parameter :: flagellates_kword      = "fla"            
      character*(*), parameter :: cyanobacteria_kword    = "cya"        
      character*(*), parameter :: organic_detritus_kword = "odt"                
      character*(*), parameter :: part_org_matter_kword  = "pom"
      character*(*), parameter :: DIC_kword              = "dic"
      character*(*), parameter :: alkalinity_kword       = "alk"
      character*(*), parameter :: DIN_kword              = "din"
      character*(*), parameter :: chlorophyl_kword       = "chl"
#endif

c
c     ------------   show time starts  ------------
c
      call init_run_context()
c      
c     ------------ set clocks  ------------
c
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      write(*,311) "start time = ", idum4
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      write(*,311) "end time = ", idum4
      call read_control_data(simulation_file, "time_step", time_step)
      write(*,312) time_step

 311  format(a12,2x,"year",1x,i4,1x,"month",1x,i2,1x,"day",1x,i2,1x,
     +       "second in day",1x,i5,1x)
 312  format("sampling time step = ", 1x,i6,1x,"seconds")
      
      call init_physical_fields(start_time)
      current_time => get_master_clock()
c
c     ------------ set horizontal ranges ------------
c
c     Either a meta area (ICES rectangle) or lon/lat range  
c     ICES_rectangle precedes longitude_range+latitude_range
c     ICES_rectangle must be accompanied by lon/lat steps
c      
      ic = count_tags(simulation_file, "ICES_rectangle") 
      if (ic > 0) then
         call read_control_data(simulation_file, "ICES_rectangle", 
     +                     cbuf)
         cbuf = adjustl(cbuf)
         icesrec = cbuf(1:4)
         write(*,322) icesrec
         call set_horizontal_range(icesrec)
         call read_control_data(simulation_file, "longitude_step",
     +                          lon_step)
         call read_control_data(simulation_file, "latitude_step",
     +                          lat_step)  
 
      else ! read an arbitrary horizontal range specification

         call read_control_data(simulation_file, "longitude_range", r3)
         lonmin   = r3(1)
         lonmax   = r3(2)
         lon_step = r3(3)
         call read_control_data(simulation_file, "latitude_range", r3)
         latmin   = r3(1)
         latmax   = r3(2)
         lat_step = r3(3) 
      endif

      write(*,321) "longitude", lonmin, lonmax, lon_step
      write(*,321) "latitude",  latmin, latmax, lat_step

 321  format("Statistics for",a10,"range = ",
     +       f9.3,1x,f9.3,1x,"with spacing",1x,f9.3)
 322  format("Statistics for ICES_rectangle", 1x,a4)
c  
c     ------------ set vertical ranges ------------   
c 
c     Either a meta layer (bottom/surface) or depth range  
c     vertical_stratum precedes vertical_range
c 
      ic = count_tags(simulation_file, "vertical_stratum") 
      at_bottom      = .false. 
      at_surface     = .false.
      single_stratum = .false.
      if (ic > 0) then
         call read_control_data(simulation_file, "vertical_stratum", 
     +                     layer_tag)
         single_stratum = .true.
         if     (trim(adjustl(layer_tag)) == "bottom") then
            at_bottom = .true.
            write(*,331) "bottom"
         elseif (trim(adjustl(layer_tag)) == "surface") then
            at_surface = .true.
            write(*,331) "surface"
         else
            write(*,*) "unknown vertical_stratum = ", 
     +                  trim(adjustl(layer_tag))
            stop
         endif

      else   ! read arbitrary vertical range specification

         call read_control_data(simulation_file, "vertical_range", r3)
         depthmin = r3(1)
         depthmax = r3(2)
         ver_step = r3(3)
         write(*,332) depthmin, depthmax, ver_step
      endif

 331  format("averaing specific vertical stratum = ",a12)
 332  format("averaing vertical range = ",f9.3,1x,f9.3,1x,
     +       "with spacing",1x,f9.3)

c  
c     ------------ set extraction property  ------------   
c       
      call read_control_data(simulation_file,"extract_data", property)
      write(*,*) "extract_data = ", trim(adjustl(property))
      call set_property_selector(property)  ! set variable property_selector
c  
c     ------------ set dump file ------------   
c      
      ic = count_tags(simulation_file, "output_file") 
      if (ic > 0) then
         call read_control_data(simulation_file,"output_file", ofname)
         call find_free_IO_unit(iunit) 
         open(iunit, file=ofname)
         write(*,*) "writing output to new file:",trim(adjustl(ofname))
      else
         ic = count_tags(simulation_file, "append_file") 
         if (ic > 0) then
            call read_control_data(simulation_file,"append_file",ofname)
            call find_free_IO_unit(iunit) 
            open(iunit, file=ofname, access='append')
            write(*,*) "appending output to file:",trim(adjustl(ofname))
         else
            write(*,*) "input error: specify output_file/append_file"
            stop
         endif
      endif
c     
c     ============== main space+time loop ==============   
c
      pavg  = 0.
      pavg2 = 0.
      pmin  =  1.e20    ! outside sensible range
      pmax  = -1.e20    ! outside sensible range
      nstatpt = 0       ! total number of sampling points
      nwetpt  = 0       ! total number of horizontal wet sampling points
      ntimpt  = 0       ! total number of time sampling points

      do while (compare_clocks(current_time, end_time) <= 0)
         call update_physical_fields()
         ntimpt = ntimpt + 1
         xyz(1) = lonmin
         do while (xyz(1) <= lonmax)
            xyz(2) = latmin
            do while (xyz(2) <= latmax)
               if (.not.is_land(xyz)) then 
                  nwetpt = nwetpt + 1
                  if (single_stratum) then
                     if (at_bottom) then
                        call interpolate_wdepth(xyz, cwd, istat) 
                        xyz(3) = cwd
                     elseif (at_surface) then
                        xyz(3) = 0.0
                     endif
                     p       = extract_data(xyz) 
                     pavg    = pavg  + p
                     pavg2   = pavg2 + p*p
                     pmin    = min(p,pmin)
                     pmax    = max(p,pmax)
                     nstatpt = nstatpt + 1
                  else          ! do vertical loop  
                     call interpolate_wdepth(xyz, cwd, istat)   
                     stop_depth = min(depthmax, cwd)
                     xyz(3) = depthmin
                     do while (xyz(3) <= stop_depth)
                        p       = extract_data(xyz) 
                        pavg    = pavg  + p
                        pavg2   = pavg2 + p*p
                        pmin    = min(p,pmin)
                        pmax    = max(p,pmax)
                        nstatpt = nstatpt + 1
                        xyz(3)  = xyz(3) + ver_step
                     enddo      ! vertical loop  
                  endif         ! if single_stratum
               endif            ! if (.not.is_land(xyz))
               xyz(2) = xyz(2) + lat_step
            enddo               ! lat loop 
            xyz(1) = xyz(1) + lon_step
         enddo                  ! lon loop             
         call add_seconds_to_clock(current_time, time_step)
      enddo                     ! time loop
c
c     ------- report basic statistics to output file -------
c      
      if (nstatpt > 0) then  
         pavg  = pavg / nstatpt
         pavg2 = pavg2 / nstatpt
         prms  = sqrt(pavg2 - pavg**2)
c        --- output file ---
         write(iunit,*) pavg,  prms, pmin, pmax
c        --- stdout ---
         write(*,610)
         write(*,620) trim(adjustl(property))
         write(*,623) pavg,  prms, pmin, pmax
         write(*,625) "total number of sampling points", nstatpt
         write(*,626) "avarage number of horizontal wet points",
     +                1.0*nwetpt/ntimpt + 1e-7
         write(*,626) "average number of vertical sampling points", 
     +                1.0*nstatpt/nwetpt + 1e-7
         write(*,610)
      else  ! undefined situation, write nothing to output
         write(*,*) "no wet sampling points in range specifications"
      endif

      close(iunit)
      call close_run_context()
      call close_physical_fields()
      write(*,*) "normal end data extraction"
 610  format(90("-"))
 620  format("spatio-temporal average of:",1x,a)     
 623  format("average =", 1x, f12.7, 3x, "RMS =",     1x, f12.7, 3x, 
     +       "min val =", 1x, f12.7, 3x, "max val =", 1x, f12.7)
 625  format(a45,1x,i12)
 626  format(a45,1x,f12.2)     
    

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine get_sw_corner(icesrec, lon_sw, lat_sw)
c     ------------------------------------------------------------------
c     Resolve the SW corner of an ICES statistical rectangle LLPP
c
c     LL = 09 -> lat =  40       latitude  step / rectangle = 0.5 degrees
c     PP = B0 -> lon = -40       longitude step / rectangle = 1.0 degrees
c     
c
c     Notes: * PP = A? is iregular so it is treated as special case   
c            * Letter in PP must be upper case
c            * This implementation makes analytical continuation onto land
c              of rectangle designation, so it will not complain about 
c              fiktive rectangles on land. All wet rectangles are correct.
c     ------------------------------------------------------------------
      character(len=4), intent(in) :: icesrec
      real, intent(out)            :: lon_sw, lat_sw
      character                    :: lonhead 
      integer                      :: latix, lonminor
c     ----------------------------------------------------------------
      read(icesrec(1:2),*) latix 
      lonhead = icesrec(3:3)
      read(icesrec(4:4),*) lonminor
      lat_sw = 40.0 + (latix-9)*0.5   ! LL = 09 -> lat(sw) = 40 ; lat step = 0.5  
      if (lonhead == "A") then        ! A? is iregular - only covers A0,A1,A2,A3
         lon_sw = -36 + 1.0*lonminor  ! lon step = 1.0
      else
         lon_sw = -40 + 10*(ichar(lonhead)-ichar("B")) + 1.0*lonminor ! lon step = 1.0
      endif
      end subroutine get_sw_corner



      subroutine set_horizontal_range(icestag)
c     ----------------------------------------------
c     Directly set lonmin, lonmax, latmin, latmax
c     in main scope
c     latitude  step / ICES rectangle = 0.5 degrees
c     longitude step / ICES rectangle = 1.0 degrees
c     ----------------------------------------------  
      character(len=4),intent(in) :: icestag
      real                        :: lon_sw, lat_sw
c     ----------------------------------------------
      call get_sw_corner(icestag, lon_sw, lat_sw)
      lonmin = lon_sw         ! in main scope
      lonmax = lon_sw + 1.0   ! in main scope
      latmin = lat_sw         ! in main scope
      latmax = lat_sw + 0.5   ! in main scope
      end subroutine set_horizontal_range



      subroutine set_property_selector(property)
c     ----------------------------------------------
c     Must be in sync with fan-out function extract_data below
c     move this code to template form and preprocess
c     ----------------------------------------------
      character*(*),intent(in) :: property
c     ----------------------------------------------
      if     (trim(adjustl(property)) == temp_kword ) then
         property_selector = select_temp
      elseif (trim(adjustl(property)) == salinity_kword) then
         property_selector = select_salinity
      elseif (trim(adjustl(property)) == current_kword) then
         property_selector = select_current
      elseif (trim(adjustl(property)) == windstress_kword ) then  
         property_selector = select_windstress
#ifdef WITH_BGC
      elseif (trim(adjustl(property)) == zooplankton_kword) then 
         property_selector = select_zooplankton
      elseif (trim(adjustl(property)) == oxygen_kword) then 
         property_selector = select_oxygen
      elseif (trim(adjustl(property)) == nh4_kword) then
         property_selector = select_nh4   
      elseif (trim(adjustl(property)) == no3_kword) then
         property_selector = select_no3
      elseif (trim(adjustl(property)) == po4_kword) then   
         property_selector = select_po4
      elseif (trim(adjustl(property)) == diatoms_kword) then   
         property_selector = select_diatoms
      elseif (trim(adjustl(property)) == flagellates_kword) then   
         property_selector = select_flagellates
      elseif (trim(adjustl(property)) == cyanobacteria_kword) then   
         property_selector = select_cyanobacteria
      elseif (trim(adjustl(property)) == organic_detritus_kword) then   
         property_selector = select_organic_detritus  
      elseif (trim(adjustl(property)) == part_org_matter_kword) then   
         property_selector = select_part_org_matter  
      elseif (trim(adjustl(property)) == DIC_kword) then   
         property_selector = select_DIC  
      elseif (trim(adjustl(property)) == alkalinity_kword) then   
         property_selector = select_alkalinity 
      elseif (trim(adjustl(property)) == DIN_kword) then   
         property_selector = select_DIN  
      elseif (trim(adjustl(property)) == chlorophyl_kword) then   
         property_selector = select_chlorophyl 
#endif 
      else
         write(*,*) "set_property_selector: unknown property",
     +               trim(adjustl(property)) 
         stop
      endif
      end subroutine set_property_selector


      real function extract_data(geo) 
c     ----------------------------------------------
c     Anonymize interpolate_X using property_selector 
c     (main scope) set at initialization
c     scalarize vector interpolations
c     waste interpolation status for now
c     move this code to template form and preprocess
c     ----------------------------------------------
      real, intent(in)    :: geo(:)
      integer, parameter  :: maxdat = 20 ! currently max 3 is used
      real                :: vecbuf(maxdat) 
      integer             :: istat 
c     ----------------------------------------------
      select case (property_selector) ! main scope
      case (select_temp) 
         call interpolate_temp(geo, extract_data, istat)
      case (select_salinity) 
         call interpolate_salty(geo, extract_data, istat)
      case (select_current) 
         call interpolate_currents(geo, vecbuf, istat)
         extract_data = sqrt(sum(vecbuf(1:3)**2))
      case (select_windstress) 
         call interpolate_wind_stress(geo, vecbuf, istat)
         extract_data = sqrt(sum(vecbuf(1:2)**2))
#ifdef WITH_BGC
      case (select_zooplankton) 
         call interpolate_zooplankton(geo, vecbuf, istat)  
         extract_data = vecbuf(1)
      case (select_oxygen) 
         call interpolate_oxygen(geo, extract_data, istat)
      case (select_nh4) 
         call interpolate_nh4(geo, extract_data, istat)  
      case (select_no3) 
         call interpolate_no3(geo, extract_data, istat)  
      case (select_po4) 
         call interpolate_po4(geo, extract_data, istat)  
      case (select_diatoms)   
         call interpolate_diatoms(geo, extract_data, istat)   
      case (select_flagellates)                
         call interpolate_flagellates(geo, extract_data, istat)       
      case (select_cyanobacteria)       
         call interpolate_cyanobacteria(geo, extract_data, istat)  
      case (select_organic_detritus)         
         call interpolate_organic_detritus(geo, extract_data, istat)  
      case (select_part_org_matter)               
         call interpolate_part_org_matter(geo, extract_data, istat)   
      case (select_DIC) 
         call interpolate_DIC(geo, extract_data, istat)  
      case (select_alkalinity)   
         call interpolate_alkalinity(geo, extract_data, istat)   
      case (select_DIN) 
         call interpolate_DIN(geo, extract_data, istat)     
      case (select_chlorophyl)  
         call interpolate_chlorophyl(geo, extract_data, istat)  
#endif  
      case default
         write(*,*) "invalid selector", property_selector
         stop 
      end select

      end function extract_data


      end program
