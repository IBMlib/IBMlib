  

      module grid_tools

      use physical_fields
      use constants

      implicit none

      contains
      
      subroutine get_cell_area(ix,iy,cellarea)
c     ---------------------------------------------------
c     compute the horizontal area of the grid cell (ix,iy) in meters**2
c     currently use tangent space approximation for simplicity
c     Then area = (jacobian(:,1) X jacobian(:,2)).(0,0,1)
c     ---------------------------------------------------
      integer, intent(in) :: ix, iy
      real, intent(out)   :: cellarea
      real                :: pos(3),jac(3,3), area 
c     ---------------------------------------------------
      pos(1) = 1.0*ix
      pos(2) = 1.0*iy
      pos(3) = 1.0   ! any vertical position
      call get_jacobian(pos, jac)
c     area is (jacobian(:,1) X jacobian(:,2)).(0,0,1)
      cellarea = -jac(1,2)*jac(2,1) + jac(1,1)*jac(2,2)
      end subroutine get_cell_area



      subroutine get_cell_volume(ix,iy,iz,cellvol)
c     ---------------------------------------------------
c     compute the volume of the grid cell (ix,iy,iz) in meters**3
c     currently use tangent space approximation for simplicity
c     ---------------------------------------------------
      integer, intent(in) :: ix, iy, iz
      real, intent(out)   :: cellvol
      real                :: pos(3),jac(3,3)
c     ---------------------------------------------------
      pos(1) = 1.0*ix
      pos(2) = 1.0*iy
      pos(3) = 1.0*iz
      call get_jacobian(pos, jac)
      cellvol = -jac(1,3)*jac(2,2)*jac(3,1) + jac(1,2)*jac(2,3)*jac(3,1) 
     +          +jac(1,3)*jac(2,1)*jac(3,2) - jac(1,1)*jac(2,3)*jac(3,2) 
     +          -jac(1,2)*jac(2,1)*jac(3,3) + jac(1,1)*jac(2,2)*jac(3,3)

      end subroutine get_cell_volume


      subroutine get_stat_moments(positions, moments1, moments2)
c     -------------------------------------------------------------------
c     Statistical centered moments of a data series positions of first 
c     order (and second, if moments2 argument is present, including covariance terms)
c     in each coordinate direction. Automatically match moment
c     calculation to the provided dimensionality of positions
c     Currently restrict calculation to 2/3 dimensions, due to metric transformation
c     Units: 
c         moments1:  grid coordinates (center of gravity)
c         moments2:  meters**2
c
c     Can easily be extended to higher order moments 
c     Do no check that moments* arrays are sufficient
c     -------------------------------------------------------------------
      real, intent(in)            :: positions(:,:)   ! (xy.., 1...maxpos)  
      real, intent(out)           :: moments1(:)      ! (x y z ...       )  
      real, intent(out), optional :: moments2(:,:)    ! (x y z ..., x y z ...)
      real                        :: metric(3)        ! restrict calculation to max 3 dimensions
      integer                     :: npos, ndim, i, j
c     -------------------------------------------------------------------
      ndim     = size(positions,1)   ! dimensionality of positions
      npos     = size(positions,2) 
     
      moments1 = sum(positions,2)/npos

c.....convert moments to local Cartesian system at center of gravity
     
      metric = 1.0d0                  ! get unit local conversion
      if (ndim==2) then
         call d_xy2d_cart(moments1, metric) ! set metric(1:2)
      elseif (ndim==3) then
         call d_xyz2d_cart(moments1, metric) ! set metric(1:3)
      else
         write(*,277) ndim
         stop
      endif


 277  format("get_stat_moments: metric not available for ndim=",1x,i3)

      if (present(moments2)) then
         do i = 1,ndim
            do j = i,ndim
               moments2(i,j) = sum((positions(i,:)-moments1(i))*
     +                             (positions(j,:)-moments1(j)))/npos
               moments2(i,j) = moments2(i,j)*metric(i)*metric(j)
               moments2(j,i) = moments2(i,j)
            enddo
         enddo
      endif
     

      end subroutine get_stat_moments



      subroutine get_rms_ellipsoid(positions, xc, yc, a, b, theta)
c     -------------------------------------------------------------------------------------
c     Determine the approximating horizontal uniform ellipsoid distribution 
c     for (grid) coordinate set positions (all data in provided array is used)
c
c     a, b:   major horizontal axes of ellipsoid (a>b) [in meters]
c             Notice that a, b is TWICE the RMS along principal ellipsoid axes
c             a, b are shape parameters of the ellipsoid distribution 
c     xc, yc: center of gravity (x,y) for positions [in grid coordinates]
c     theta:  rotation angle of larger axis a [in radians] with respect to WE tangent
c     ------------------------------------------------------------------------
      real, intent(in)  :: positions(:,:)   ! (xyz, 1...maxpos)  z not used/tested for - need not be present
      real, intent(out) :: xc, yc, a, b, theta
      real              :: xybar(2), k20, k02, k11, kappa  ! xybar/metric dimension == (2) to use d_xy2d_cart
      real              :: c2, s2, tmp, metric(2)
      real              :: xy(2,size(positions,2))  ! automatic array
      integer           :: npos
      real              :: mom1(2), mom2(2,2)    ! (
c     -----------------------------------------------------------------
      call get_stat_moments(positions(1:2,:), mom1, mom2) ! only perform xy statistics
      
      xc  = mom1(1)
      yc  = mom1(2)
      k20 = mom2(1,1)
      k02 = mom2(2,2)
      k11 = mom2(1,2)

c     -------------------------------------------------------------------------------------
c     First compute fitting ellipsoid (a,b,theta) where theta are angle between axes x and a
c     ifort atan branch def returns:  -0.25 pi < theta < 0.25 pi  (incl factor 0.5 below)
c     a-axis is defined as the principal axis, which lies in this range
c     -------------------------------------------------------------------------------------
      kappa = (2.*k11)/(k20-k02)
      theta = 0.5*atan(kappa)  ! just pick atan implementation branch
      c2    = cos(theta)**2
      s2    = sin(theta)**2
      a     = sqrt( 4.0d0*( c2*k20 - s2*k02) / (c2-s2) )
      b     = sqrt( 4.0d0*(-s2*k20 + c2*k02) / (c2-s2) )


c     -------------------------------------------------------------------------------------
c     Flip axes, if nescessary, so that a>b, and correspondingly theta
c     The solution follows branch shift symmetries:
c       I)   (a, b, theta) -> (b, a, theta+pi/2)
c       II)  (a, b, theta) -> (a, b, theta+pi)              (implied by symmetry I)
c     Apply symmetries to pick solution safisfying 
c       1)   a > b
c       2)  -0.5 pi < theta < 0.5 pi
c     -------------------------------------------------------------------------------------
      if (a<b) then                                ! apply symmetry I
         tmp   = b
         b     = a
         a     = tmp
         theta = theta + pi/2.0d0
      endif
      theta = theta - pi*floor(theta/pi + 0.5d0)   ! apply symmetry II
     
c.....output result      
c      write(*,333)  a, b, theta
c 333  format(" a=",f12.7," b=",f12.7," theta=",f12.7)

      end subroutine get_rms_ellipsoid



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Evaluate concentration in boxes in units of tracers/volume
c     
c     Only provide computation on physical grid (as exported by physical_fields)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_xyz_concentration(pos, conc)
c     ---------------------------------------------------
c     Core updater of 3D concentration array
c     Discard positions outside valid grid range
c     Return conc = 0, if all pos violates valid grid range
c     ---------------------------------------------------
      real, intent(in)  :: pos(:,:)     ! (xyz, 1...maxpos) 
      real, intent(out) :: conc(:,:,:)  ! (nx,ny,nz)
c     ---------------------------------------------------
      integer :: ix, iy, iz, ipos, n_ok
      real    :: cellvol
c     ---------------------------------------------------
      conc = 0.
      n_ok = 0.
      do ipos = 1, size(pos,2)
         ix   = nint(pos(1,ipos))
         if ( (ix>nx) .or. (ix<1) ) cycle
         iy   = nint(pos(2,ipos))
         if ( (iy>ny) .or. (iy<1) ) cycle
         iz   = nint(pos(3,ipos))
         if ( (iz>nz) .or. (iz<1) ) cycle
c........this position is within valid range, accept it
         n_ok = n_ok + 1
         conc(ix,iy,iz) = conc(ix,iy,iz) + 1.0d0
      enddo
c.....divide by box volumes  
      do ix=1,nx
         do iy=1,ny
            do iz=1,nz
               call get_cell_volume(ix,iy,iz,cellvol)
               conc(ix,iy,iz) = conc(ix,iy,iz)/cellvol
            enddo
         enddo
      enddo

      end subroutine get_xyz_concentration


      subroutine get_xy_concentration(pos, conc)
c     ---------------------------------------------------
      real, intent(in)  :: pos(:,:)          ! (xyz, 1...maxpos) 
      real, intent(out) :: conc(:,:)         ! (nx,ny)
      real              :: xyzconc(nx,ny,nz) ! automatic aux array
c     ---------------------------------------------------
      call get_xyz_concentration(pos, xyzconc)
      conc = sum(xyzconc, 3)         ! sum over nz
      end subroutine get_xy_concentration

      subroutine get_z_concentration(pos, conc)
c     ---------------------------------------------------
      real, intent(in)  :: pos(:,:)          ! (xyz, 1...maxpos) 
      real, intent(out) :: conc(:)           ! (nz)
      real              :: xyzconc(nx,ny,nz) ! automatic aux array
c     ---------------------------------------------------
      call get_xyz_concentration(pos, xyzconc)
      conc = sum(sum(xyzconc,1),1)   ! sum over nx,ny
      end subroutine get_z_concentration

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Concentration writers   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
  
      subroutine write_xyz_concentration(filename, pos)
c     ---------------------------------------------------
      character*(*), intent(in) :: filename
      real, intent(in)          :: pos(:,:)  ! (xyz, 1...maxpos)
      integer                   :: iunit,ix,iy,iz
      real                      :: xyzconc(nx,ny,nz) ! automatic aux array
c     ---------------------------------------------------
      call find_free_IO_unit(iunit)
      open(iunit, file=filename)
      call get_xyz_concentration(pos, xyzconc)
      do ix=1,nx
         do iy=1,ny
            do iz=1,nz
               write(iunit,203) ix, iy, iz, xyzconc(ix,iy,iz)
            enddo
         enddo
      enddo     
      close(iunit)
 203  format(3i4, 1x, e12.6)     
      end subroutine write_xyz_concentration


      subroutine write_xy_concentration(filename, pos)
c     ---------------------------------------------------
      character*(*), intent(in) :: filename
      real, intent(in)          :: pos(:,:)  ! (xyz, 1...maxpos)
      integer                   :: iunit,ix,iy
      real                      :: xyconc(nx,ny) ! automatic aux array
c     ---------------------------------------------------
      call find_free_IO_unit(iunit)
      open(iunit, file=filename)
      call get_xy_concentration(pos, xyconc)
      do ix=1,nx
         do iy=1,ny
            write(iunit,202) ix, iy, xyconc(ix,iy)
         enddo  
      enddo     
      close(iunit)
 202  format(2i4, 1x, e12.6)
      end subroutine write_xy_concentration


      subroutine write_z_concentration(filename, pos)
c     ---------------------------------------------------
      character*(*), intent(in) :: filename
      real, intent(in)          :: pos(:,:)  ! (xyz, 1...maxpos)
      integer                   :: iunit,iz
      real                      :: zconc(nz) ! automatic aux array
c     ---------------------------------------------------
      call find_free_IO_unit(iunit)
      open(iunit, file=filename)
      call get_z_concentration(pos, zconc)
      do iz=1,nz
         write(iunit,201) iz, zconc(iz)      
      enddo     
      close(iunit)
 201  format(1i4, 1x, e12.6)
      end subroutine write_z_concentration



      end module

      
c$$$      program jj
c$$$c.....ifort -e90 grid_tools.f  input_parser.o
c$$$      use grid_tools
c$$$      use physical_fields
c$$$      real    :: pos(3,2000)
c$$$c     ---------------------------------------
c$$$      call random_number(pos)
c$$$      pos(1,:) = pos(1,:)   + 2.5
c$$$      pos(2,:) = 2*pos(2,:) + 3.5
c$$$      pos(3,:) = 3*pos(3,:) + 4.5
c$$$      nx = 8
c$$$      ny = 9
c$$$      nz = 10
c$$$      call write_xyz_concentration("jjxyz", pos)
c$$$      call write_xy_concentration("jjxy", pos)
c$$$      call write_z_concentration("jjz", pos)
c$$$      end program
