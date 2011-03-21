      program polcoms_proxy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    
c     compile test: 
c        ifort  -o coupled_run task_providers/online/polcoms_proxy.f libibmlib.a   -I/home/asbjorn/IBMlib/wrk9/trunk/oceanography_providers/POLCOMS+ERSEM/online -I/home/asbjorn/IBMlib/wrk9/trunk/biology_providers/passive   -I/home/asbjorn/IBMlib/wrk9/trunk/task_providers/online
c  
c     test run:   
c        coupled_run task_providers/online/simfile_POLCOMS_online
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use ibmlib
      implicit none
    
      real     :: dt=500.0 ! seconds
      integer  :: i,iz,i1,i2
      integer,parameter  :: zdim=9,londim=14,latdim=17 ! grid dimensions
      type(physical_data_bucket) :: pol_data   ! type imported from IBMlib


c     --- allocate or point to data
c         fields in pol_data are not modified (only read)  
c    
      allocate( pol_data%ucur(zdim,londim,latdim)   ) ! lon current, pos east
      allocate( pol_data%vcur(zdim,londim,latdim)   ) ! lat current, pos north
      allocate( pol_data%wcur(zdim,londim,latdim)   ) ! vertical current, pos down
      allocate( pol_data%ETW(zdim,londim,latdim)    ) ! temperature (C)
      allocate( pol_data%nuv(zdim,londim,latdim)    ) ! vertical CC diffusivity
      allocate( pol_data%depth(zdim,londim,latdim)  ) ! CC depths
      allocate( pol_data%zbnd(2,zdim,londim,latdim) ) ! face depths; IBMlib just use zbnd(2,nz,:,:)


      call init_ibmlib()

c     --- main loop ---
      do i=1,100

         write(*,281) i
 281     format(20('='), 1x," polcoms_proxy step ", i5, 1x, 20('='))
         pol_data%ucur=1.  ! emulate POLCOMS calculation
         pol_data%vcur=2.  ! emulate POLCOMS calculation
         pol_data%wcur=3.  ! emulate POLCOMS calculation
         pol_data%ETW=4.   ! emulate POLCOMS calculation
         pol_data%nuv=5.   ! emulate POLCOMS calculation
         do iz=1,zdim      ! apply static surface in this example
            pol_data%depth(iz,:,:) = -2*(iz-0.5) ! cell width = 2 meters, below water negative
         enddo
         pol_data%zbnd(2,zdim,:,:) = -2*zdim     ! cell width = 2 meters, below water negative

         call transfer_data_to_ibmlib(pol_data)
         call ibm_step(dt)
      enddo

      call close_ibmlib()

      end 
