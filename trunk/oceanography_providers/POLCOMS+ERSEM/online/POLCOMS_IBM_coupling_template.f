ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Template for coupling IBMlib with POLCOMS+ERSEM online
c
c     Makefile: produce libibm.a
c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      module pbi
      implicit none
      private
      real(kind=8), pointer :: pt2u(:,:,:) 
      integer               :: londim,latdim,zdim

      type polcoms_data_bucket
        real(kind=8),pointer :: pol_u(:,:,:)
        real(kind=8),pointer :: pol_t(:,:,:)
c       ....
      end type

      public :: polcoms_data_bucket
      public :: init_pbi, transfer_data, close_pbi

      contains 

      subroutine init_pbi(ilondim,ilatdim,izdim)
      integer :: ilondim,ilatdim,izdim,istat
      londim = ilondim
      latdim = ilatdim
      zdim   = izdim
      allocate(pt2u(londim,latdim,zdim), STAT=istat)
      write(*,*) "init_pbi: allocate istat =", istat
      pt2u = 0.
      end subroutine init_pbi

      subroutine transfer_data(pol_data)
c     POLCOMS def:  u(zdim,londim,latdim) ) 
      type(polcoms_data_bucket) :: pol_data
      pt2u = reshape(pol_data%pol_u, (/londim,latdim,zdim/),
     +               order=(/2,3,1/))
c     ...      
      call adapt_data()
      end subroutine transfer_data

      subroutine adapt_data()
      end subroutine adapt_data

      subroutine close_pbi()
      deallocate(pt2u)
      end subroutine close_pbi

      end module

      


      module IBMlib
      use pbi
      implicit none
      private
      public :: transfer_data, polcoms_data_bucket
      public :: init_IBMlib, IBM_step, close_IBMlib

      contains

      subroutine init_IBMlib(londim,latdim,zdim)
      integer :: londim,latdim,zdim
      call init_pbi(londim,latdim,zdim)
      end subroutine init_IBMlib
      
      subroutine IBM_step(dt)
      real(kind=8)             :: dt
      end subroutine IBM_step

      subroutine close_IBMlib()
      call close_pbi()
      end subroutine  close_IBMlib

      end module



      
      program polcoms
      use IBMlib
      implicit none
    
      real(kind=8)              :: dt=0.56
      integer                   :: i
      integer                   :: zdim=6,londim=4,latdim=9
      type(polcoms_data_bucket) :: pol_data   ! type from IBMlib
      allocate( pol_data%pol_u(zdim,londim,latdim) ) 

      call init_IBMlib(londim,latdim,zdim)

c     --- main loop ---
      do i=1,100
         pol_data%pol_u=37*i
         call transfer_data(pol_data)
         call IBM_step(dt)
      enddo

      call close_IBMlib()
      end 
