ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     Interpolators
c     -----------------------------------------------------------
c     $Rev: 231 $
c     $LastChangedDate: 2011-01-26 11:07:27 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: mpay $ 
c
c     This component of the NORWECOM oceanography provider deals with
c     accessor functions, speicifically with returning a given physical
c     parameter at a given point in space. It relies heavilty on the
c     functions contained in the NORWECOM_grid.f file. This code is
c     intended to be incorporated into the NORWECOM_main.f file using
c     an include statement.
c
c     Dealing with data points that are close to land is a bit tricky
c     as it essentially involves an extrapolation in the region between
c     the centre of the last cell and the coast. Here we solve the
c     problem in a pragmatic manner by switching from linear
c     interpolation to neearest-neighbour in cases where the
c     particle is in this region (i.e. within 0.5 of a cell of the
c     coast). Alternative schemes can be imagined, but this is the
c     simplest to implement in the meantime
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine interpolator_component_version()
c     Writes the SVN revision number out for logging purposes
      write(*,*) "Interpolators : $Rev: 231 $"
      end subroutine interpolator_component_version

c########################################################################
c########################################################################
c#                    Public Accessor Methods                           # 
c########################################################################

      subroutine interpolate_currents(geo, vels, status)
c     ----------------------------------------------------------------
c     Perform linear two-face interpolation in current fields 
c     position, geo, is suppliered as (lon, lat, depth)   shape == real(3+)
c     vels are interpolated currents in m/s, shape == real(3+), parallel
c     to the lon, lat and depth axis (increasing towards the east, north
c     and downwards)
c     Interpolation is limited to the data bounds, and will return an
c     error if we try and go outside those bounds.
c     All velocities on land are assumed to be padded with an
c     appropriate value - in this case zero. The effect of land on the
c     interpolation is therefore quite different to the other properties
c     (e.g. temperature, salinity) and can essentially be ignored.
c     ------------- conventions / grid layout ------------------------
c     (x,y) positive along grid coordinates, z positive downwards
c     cell (ix,iy,iz) centered at xyz = (ix-0.5,iy-0.5,iz-0.5)
c     u,v positive along lat,lon coordinate system in units m/s
c     w   positive upward                          in units m/s
c     uvw is a vector holding these velocities
c    
c     u(ix,iy,iz) at grid position (ix-1  , iy    , iz)     (i.e. western  cell face)
c     v(ix,iy,iz) in grid position (ix    , iy-1.0, iz)     (i.e. southern cell face)
c     w(ix,iy,iz) in grid position (ix    , iy    , iz-1.0) (i.e. upper    cell face)
c     The boundary condition:
c           w( ix, iy, nz ) = 0
c     is implicit 
c
c     Implied interpolation ranges:
c           0.0 < x < nx-1
c           0.0 < y < ny-1
c           0.0 < z < nz (due to bottom BC)
c     ----------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: vels(:)
      integer, intent(out) :: status
      real              :: uvw(3) 
      integer           :: ix,iy,iz,nexteast,nextnorth,nextlow,ibot
      integer           :: nextrp, nintp, idum(3),iunit
      real              :: sx,sy,sz,xyz(3)
c     ----------------------------------------------------------------
c     Check bounds - if out of bounds return to calling function
      if (horizontal_range_check(geo) ) then
         status=0 !In range 
      else
         status=1  !Trying to interpolate out of bounds
         vels= NaN
         return    !so just return nothing here
      endif
c     No further bounds checks are performed beyond this point
c     Now we now that we are within the range of values
c     and can use direct lookups

c.....determine integer cell associations of point, constrained to 1 <= i <= n
c     i.e. the indices of the cell in which position is to be found 
      call geo2xyz(geo,xyz)
      call cell_index(xyz,ix,iy,iz)
      ibot = bottom_layer(ix,iy)       ! For NORWECOM, this always equals nz, as it is terrain following

c.....determine relative intracell coordinate s, constrained to 0 < s < 1 
c     we have already checked whether the variable is in or out of range
c     so safety mechanism for when point exceed grid range is removed
      sx = xyz(1)-real(ix)+1.0d0
      sy = xyz(2)-real(iy)+1.0d0           
      sz = xyz(3)-real(iz)+1.0d0

c.....face-to-face interpolation in physical units
c     cell indices satisfy 1 <= (ix,iy,iz) <= (nx-1,ny-1,ibot)
      uvw(1) = (1.-sx)*u(ix, iy, iz)    + sx*u(ix+1, iy, iz)
      uvw(2) = (1.-sy)*v(ix, iy,iz)     + sy*v(ix, iy+1, iz)
      if (iz==ibot) then        !If in the bottom cell
         uvw(3) = (1.-sz)*w(ix, iy, iz) ! w=0 at sea bed, iz corresponds to upper cell face
      else
         uvw(3) = (1.-sz)*w(ix, iy, iz) + sz*w(ix, iy, iz+1) 
      endif
        
c.....resolve vector components from being  oriented 
c     along the NORWECOM grid to oriented along lon, lat
      call resolve_vector(xyz,uvw,vels)

      end subroutine interpolate_currents 


      subroutine interpolate_turbulence(geo, k, status)
c     ----------------------------------------------------------------
c     Interpolate module turbulence data 
c
c                 hdiffus(nx,ny,nz) -> k(1:2)
c                 vdiffus(nx,ny,nz+1) -> k(3)
c 
c     on position (in grid coordinates)
c     hdiffus/vdiffus are assumed up to date (no update invoked)
c     The output buffer k(:) is assumed sufficiently large (not checked)  
c     All hard work is done by the turbulence interpolator
c     ----------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: k(:)
      integer, intent(out) :: status
c     ----------------------------------------------------------------------------
      call turbulence_interpolator(geo, k,.FALSE., status)
      end subroutine interpolate_turbulence


      subroutine interpolate_turbulence_deriv(geo, dk, status)
c     ----------------------------------------------------------------
c     Interpolation derivatives of module turbulence data with respect to 
c     scaled coordinates (NOT Cartesian coordinates)
c
c              (d/dx, d/dy) hdiffus(nx,ny,nz)   -> dk(1:2)
c              (d/dz)       vdiffus(nx,ny,nz+1) -> dk(3)
c 
c     at position  (in grid coordinates)
c     hdiffus/vdiffus are assumed up to date (no update invoked)
c     All hard work is done by the turbulence interpolator
c     ----------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: dk(:)
      integer, intent(out) :: status
c     ----------------------------------------------------------------------------
      call turbulence_interpolator(geo,dk,.TRUE., status)
      end subroutine interpolate_turbulence_deriv 


      subroutine interpolate_temp (geo, r, status) 
c     ----------------------------------------------------------------------------
c     Temperature is measured at the centre of the grid cell i.e.
c     the temperature of the ith, jth cell is at (i-0.5,j-0.5)
c     Function is wrapper to the cell-centered interpolation routine
c     ----------------------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      real, pointer        :: dat(:,:,:)
c     ----------------------------------------------------------------------------
      dat => temp
      call interpolate_CC_data(geo, dat,r, status)
      end subroutine 


      subroutine interpolate_salty(geo, r, status)
c     ----------------------------------------------------------------------------
c     Salinity is measured on the centre of the grid cell i.e.
c     the salinity of the ith, jth cell is at (i-0.5,j-0.5)
c     Function is wrapper to the cell-centered interpolation routine
c     ----------------------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      real, pointer        :: dat(:,:,:)
c     ----------------------------------------------------------------------------
      dat => salinity
      call interpolate_CC_data(geo,dat, r,status)
      end subroutine 


      subroutine interpolate_wind(ll, r2, status)
c     ----------------------------------------------------------------------------
c     NORWECOM does not provide wind. Therefore throw error if called
c     ----------------------------------------------------------------------------
      real, intent(in)     :: ll(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
c     ----------------------------------------------------------------------------
      status=1
      r2=NaN
      call abort("interpolate_wind","No wind provided in NORWECOM")
      end subroutine 


      subroutine interpolate_wdepth(ll, r, status) 
c     ----------------------------------------------------------------------------
c     Interpolation of total water depth r at geo position (lon,lat) corresponding
c     to z coordinates. Replaces direct array access, hides coordinate transform
c     return result in meters
c     NORWECOM bottom is treated as being jagged, therefore no interpolation
c     ----------------------------------------------------------------------------
      real, intent(in)     :: ll(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      real                 :: xy(2)
      integer              :: ix,iy
c     --------------------------------------------------      
c     Check data range - return NaN if out of range
      call lonlat2xy(ll,xy)
      if(data_range_check(xy(1:2))) then
         status=0
      else
         status=1   !Out of range
         r=NaN
         return
      endif
c     Check wet/dry - if cell is dry, retur NaN
c     Cell is assumed to be wet from hereon 
      call cell_index(xy,ix,iy)
      if (land(ix,iy)) then
         status=2  !In range but Dry point
         r=NaN
         return    !so just return NaN here
      endif
      r  = bath(ix,iy)  
      end subroutine 


      subroutine interpolate_CC_data(geo,dat,res,status)
c     ----------------------------------------------------------------
c     Interpolate a cell centered 3D array dat at position,geo 
c     (lon,lat, depth), using appropriate interpolator. This interpolation 
c     applies to vertically flat cell boundaries - NORWECOM grid space is
c     orthogonal and regular (in grid space, not in real space) and
c     therefore we can get away with such tricks.
c     Assume dimensions dat(1:nx, 1:ny, 1:nz) (not checked)
c     Surface/bottom breaking and nearby land cells is handled by the
c     data extraction subroutine, get_data_octet
c     ----------------------------------------------------------------
      real, intent(in)  :: geo(3) ! assumed shape real(3+)
      real, intent(in)  :: dat(:,:,:) 
      real, intent(out) :: res
      integer, intent(out) :: status 
      integer           :: ixyz(3),cnr000(3), cnr111(3),i,j
      real              :: vc(8),xyz(3),sxyz(3)
c     ----------------------------------------------------------------
c     Check bounds - if out of bounds return to calling function
c     Not further bounds checks are performed beyond this point
      if (horizontal_range_check(geo)) then
         status=0 !In range 
      else
         status=1  !Trying to interpolate out of bounds
         res=NaN
         return    !so just return NaN here
      endif
c     Now we now that we are within the range of values
c     and can use direct lookups
c     Check wet/dry - if cell is dry, return
      call geo2xyz(geo,xyz)
      call cell_index(xyz,ixyz(1),ixyz(2),ixyz(3))
      if (land(ixyz(1),ixyz(2))) then
         status=2  !In range but Dry point
         res=NaN
         return    !so just return NaN here
      endif
c     Cell is assumed to be wet from hereon 
c     retrieve corner positions of cube enclosing xyz 
      call get_CC_cube(xyz,cnr000,cnr111,sxyz)
c     retrieve data points at those corners      
      call get_data_octet(dat,ixyz,cnr000, cnr111,vc)
c     Farm out the actual interpolation to other functions
      call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,0,res) 
      end subroutine interpolate_CC_data


      subroutine turbulence_interpolator(geo, r, deriv, status)
c     ----------------------------------------------------------------
c     Interpolate module turbulence data 
c
c                 hdiffus(nx,ny,nz) 
c                 vdiffus(nx,ny,nz+1) 
c 
c     based on position (in geo coordinates)
c     Because of similarities in the way the data is setup between
c     calculating the value and the derivatives, we have combined
c     the two potential functions into one here
c     hdiffus/vdiffus are assumed up to date (no update invoked)
c     The output buffer r(:) is assumed sufficiently large (not checked)  
c     
c     The interpolation scheme follows that employed in
c     interpolate_CC_data, especially with regard to surface/ottom crossing
c     and interaction with land
c     ----------------------------------------------------------------
      real, intent(in)     :: geo(:) ! assumed shape real(3+)
      real, intent(out)    :: r(:)
      integer, intent(out) :: status
      logical, intent(in)  :: deriv
c     ------locals------
      integer           :: cnr000(3), cnr111(3),ixyz(3)
      real              :: sxyz(3),xyz(3),k,vc(8)
      real              :: invjacob(2,2),trans_invjacob(2,2),dDdgeo(2)
      real              :: invdistmat(2,2)
      real, pointer     :: dat(:,:,:)
c     ----------------------------------------------------------------
c     Check bounds - if out of bounds return to calling function
c     No further bounds checks are performed beyond this point
      if (horizontal_range_check(geo)) then
         status=0 !In range 
      else
         status=1  !Trying to interpolate out of bounds
         r =NaN    
         return    !so just return NaN here
      endif
c     Now we now that we are within the range of values
c     and can use direct lookups
c     Check wet/dry - if cell is dry, return
      call geo2xyz(geo,xyz)
      call cell_index(xyz,ixyz(1),ixyz(2),ixyz(3))
      if (land(ixyz(1),ixyz(2))) then
         status=2  !In range but Dry point
         r =NaN    
         return    !so just return NaN here
      endif
c     Cell is assumed to be wet from hereon 
c     
c.....Determine horizontal diffusivity interpolation cube
c     hdiffus is cell centered
      call get_CC_cube(xyz,cnr000,cnr111,sxyz)
c     Extract data from hdiffusivity     
      dat => hdiffus
      call get_data_octet(dat,ixyz,cnr000,cnr111,vc)
c     Farm out the interpolatioins
      if(deriv) then
        call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,1,r(1)) 
        call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,2,r(2)) 
      else 
        call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,0,k)
        r(1:2) = k
      endif
c
c.....Determine vertical diffusivity interpolation cube
c     vdiffus is face centered at upper face, so z index + offset 
c     are different from horizontal diffusivity array, but x,y same
      cnr000(3) = ceiling(xyz(3)) ! shallower z face
      cnr111(3) = cnr000(3) +1     ! deeper face 
      sxyz(3)  = xyz(3)-real(cnr000(3))+1.0 
c     Extract data from vdiffusivity     
      dat => vdiffus
      call get_data_octet(dat,ixyz,cnr000,cnr111,vc)
c     Farm out the interpolatioins
      if(deriv) then
        call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,3,r(3)) 
      else 
        call interp_3Dbox_data(sxyz(1),sxyz(2),sxyz(3),vc,0,r(3))
      endif

c.....Transform and rotate grids
c     No transformation is required for the normal turbulences. However,
c     the derivatives are given as m^2/s / grid unit oriented along the
c     grid axes, xyz, whereas we need to return them in units of m/s
c     oriented along the spherical geographic coordinate system.
c     The transformation is performed as described in the NORWECOM pbi 
c     documentation using matrix multiplication
c           dDdst = invdistmat %*% jacob^-T %*% dDdxy
c     to give the derivates oriented along the geographical coordinate
c     system in the desired units m^2/s / m
      if(deriv) then
         call get_inv_jacobian(xyz,invjacob)
         trans_invjacob=transpose(invjacob)
         dDdgeo = matmul(trans_invjacob,r(1:2)) !Intermediate variable
         invdistmat(1,1)=1/(earth_radius*cos(geo(2)*deg2rad))
         invdistmat(1,2)=0 
         invdistmat(2,1)=geo(1)*deg2rad*tan(geo(2)*deg2rad)/earth_radius
         invdistmat(2,2)=1/earth_radius
         r(1:2) = matmul(invdistmat,dDdgeo)
         r(3) = r(3)/layer_width(ixyz(1),ixyz(2),ixyz(3)) 
      endif
      end subroutine turbulence_interpolator


      subroutine get_data_octet(dat,ixyz,i0,i1,octet)
c     ----------------------------------------------------------------
c     Extracts the data points at the corners of the cube defined by i0
c     and i1 and places them into a format suitable for use by 
c     an interpolattion function
c     As this function looks up the data directly in the array, surface/bottom
c     crossing and land issues must be handled here. 
c      * Surface/bottom crossing is dealt with by flattening the 3D cube
c        into a flat 2D slab
c      * If one or more of the corner cells is land, we can strike all
c        sorts of problems with interpolation. In this case, we take the
c        easiest approach by enforcing a nearest neighbour interpolation 
c        in the horizontal (but not the vertical)
c     ----------------------------------------------------------------
      integer, intent(inout) :: ixyz(3),i0(3), i1(3)
      real, intent(in)    :: dat(:,:,:)
      real, intent(out)   :: octet(8)
      logical cnr_on_land(0:1,0:1)
c     ----------------------------------------------------------------
c     Surface/bottom crossing is relatively straightforward to deal with
c     If cube has broken the surface or bottom, collapse the 
c     two horizontal layers into one
      i0(3) = max(1,min(i0(3) ,nz))
      i1(3) = max(1,min(i1(3) ,nz))
c     Now check whether we are in the near-coastal region ie
c     one of the cube corners is on land
      cnr_on_land(0,0) = land(i0(1),i0(2)) 
      cnr_on_land(1,0) = land(i1(1),i0(2))
      cnr_on_land(0,1) = land(i0(1),i1(2)) 
      cnr_on_land(1,1) = land(i1(1),i1(2)) 
      if(any(cnr_on_land)) then !use nearest neighbour in the horizontal
        i0(1:2) = ixyz(1:2)
        i1(1:2) = ixyz(1:2)
      endif
c     retrieve data points at the corners      
      octet(1) = dat(i0(1), i0(2), i0(3))
      octet(2) = dat(i0(1), i0(2), i1(3))
      octet(3) = dat(i0(1), i1(2), i0(3))
      octet(4) = dat(i0(1), i1(2), i1(3))
      octet(5) = dat(i1(1), i0(2), i0(3))
      octet(6) = dat(i1(1), i0(2), i1(3))
      octet(7) = dat(i1(1), i1(2), i0(3))
      octet(8) = dat(i1(1), i1(2), i1(3))
      end subroutine 

