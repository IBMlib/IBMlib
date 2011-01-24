ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Pbi_Interface test 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Checks that a pbi provides the minimum set of functions with the appropriate
c     arguments required by the interface definition. The purpose
c     of the test is not to perform calculations - if the compiler
c     builds it, then the test is passed
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use input_parser
      use time_tools
      use run_context
      use physical_fields

      implicit none
c     ------------ declarations ------------      
      type(clock),target  :: clck
      integer ok,id
      real geo(3), geo2(3),geo3(3),geo4(3), r3(3), r2(2),r,ll(2),dt
      logical logic
      character(len=123) :: str

c     ------------ check physical fields interface ------------      
      call init_physical_fields(clck)
      call close_physical_fields()
      clck =  get_master_clock()
      call set_master_clock(clck)
      call update_physical_fields()
      call interpolate_turbulence(geo,r3,ok)
      call interpolate_turbulence_deriv(geo,r3,ok)
      call interpolate_currents(geo,r3,ok)
      call interpolate_temp(geo,r,ok) 
      call interpolate_salty(geo,r,ok)
      call interpolate_wind(geo,r2,ok)
      call interpolate_wdepth(geo,r,ok)
      logic = is_wet(ll) 
      logic = is_land(ll)
      call coast_line_intersection(geo,geo2,logic,geo3,geo4)
      logic= horizontal_range_check(ll)
      write(*,*) get_pbi_version()

      end program
