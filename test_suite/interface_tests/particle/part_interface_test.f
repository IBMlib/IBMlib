ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Interface test 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Checks that a combination of pbi and particle providers
c     provide the minimum set of functions with the appropriate
c     arguments required by the interface definition. The purpose
c     of the test is not to perform calculations - if the compiler
c     builds it, then the test is passed
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use input_parser
      use time_tools
      use run_context
c      use particle_tracking
      use particle_state
      use output

      implicit none
c     ------------ declarations ------------      
      type(state_attributes) :: state
      type(spatial_attributes) :: space
      type(variable) :: var
      type(polytype) :: bucket
      integer ok,id
      real geo(3), geo2(3),geo3(3),geo4(3), r3(3), r2(2),r,ll(2),dt
      logical logic
      character(len=123) :: str

c     ------------ check particle interface ------------      
      call init_particle_state()
      write(*,*) get_particle_version()
      call close_particle_state()
      call init_state_attributes(state,space,dt,str,id)
      call get_active_velocity(state,space,r3)
      call update_particle_state(state,space,dt)
      call delete_state_attributes(state) 
      call write_state_attributes(state)
      call get_property(state,var,bucket,ok)
      call get_metadata_state(str,var,ok)
      end program
