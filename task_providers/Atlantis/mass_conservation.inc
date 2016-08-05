cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Include functionality for module atlantis_grid
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_mass_conservation(hflux, vflux)
c     ----------------------------------------------------------------------------------
c     Diagnostic subroutine to assess relative mass conservation violation rate (1/V)dV/dt
c     for horizontal fluxes hflux and vertical fluxes vflux
c     hflux positive from right to left box of faces
c     Positive direction for vflux across faces is upward;
c     vflux(i,j) is the flux from layer i+1 to layer i (offset 1) in box j (offset 0)
c
c     Write result to stdout - report with IBMlib layer indexing
c     ----------------------------------------------------------------------------------
      real, intent(in)            :: hflux(:,0:) ! [m3/s] shape = (1:nlayers, 0:nfaces-1)
      real, intent(in)            :: vflux(:,0:) ! [m3/s] shape = (1:nlayers, 0:nboxes-1)
      integer                     :: ibox, ilay, ipt, i
      integer                     :: nboxes, nlayers, nfaces
      integer                     :: istat, ifc, nscan
      integer                     :: lbox, rbox
      type(AtlantisFace), pointer :: face
      type(AtlantisBox), pointer  :: left, right, box
      real,allocatable            :: fluxsum(:,:), absflxs(:,:)
      real                        :: area, h, uvw(3), geo(3)
c
      integer, parameter          :: probefaces(3) = (/88,89,94/)
c     ----------------------------------------------------------------------------------
      nboxes  = size(Boxes)
      nlayers = size(zgrid_CC) 
      nfaces  = size(Faces)
      allocate( fluxsum(1:nlayers, 0:nboxes-1) )
      allocate( absflxs(1:nlayers, 0:nboxes-1) )
      fluxsum  = 0.     ! set pad value for dry layers
      absflxs  = 1.e-6  ! avoid zero division later 

c
c     ---------- First scan horizontal fluxes ----------
c
      do ifc = 0, nfaces-1
         face  => Faces(ifc)
         lbox  = face%left_right(1)
         rbox  = face%left_right(2)
         left  => Boxes(lbox)
         right => Boxes(rbox)  
         nscan  = min(left%nwet, right%nwet)
c        ...... fluxsum(ilay,ibox) is flux into (ilay,ibox) ......
         do ilay = 1, nscan  ! top-down
            fluxsum(ilay,lbox) = fluxsum(ilay,lbox) + hflux(ilay,ifc)  ! flux positive from right to left
            fluxsum(ilay,rbox) = fluxsum(ilay,rbox) - hflux(ilay,ifc)  ! flux positive from right to left
            absflxs(ilay,lbox) = absflxs(ilay,lbox) + 
     +                           abs(hflux(ilay,ifc))
            absflxs(ilay,rbox) = absflxs(ilay,rbox) + 
     +                           abs(hflux(ilay,ifc))
         enddo ! ilay
      enddo    ! ifc
c
c     ---------- Then scan vertical fluxes ----------
c    
      do ibox = 0, nboxes-1
         box => Boxes(ibox)
         do ilay = 1, box%nwet  ! top-down, vflux(i,j) is the flux from layer i+1 to layer i
             fluxsum(ilay,ibox) = fluxsum(ilay,ibox) + vflux(ilay,ibox)
             absflxs(ilay,ibox) = absflxs(ilay,ibox) + 
     +                            abs(vflux(ilay,ibox))
             if (ilay < box%nwet) then
                 fluxsum(ilay+1,ibox) = fluxsum(ilay+1,ibox) 
     +                                  - vflux(ilay,ibox)
                 absflxs(ilay+1,ibox) = absflxs(ilay+1,ibox) 
     +                                  + abs(vflux(ilay,ibox))    
             endif
         enddo ! ilay
      enddo    ! ibox
c
c     ---------- scale fluxsum by cell volumes ----------
c    
      do ibox = 0, nboxes-1
         area = sum(Boxes(ibox)%da)
         do ilay = 1, nlayers
            h = zgrid_bounds(ilay+1) - zgrid_bounds(ilay) 
            absflxs(ilay,ibox) = abs(fluxsum(ilay,ibox))  ! evaluate accuracy before rescaling
     +                           / absflxs(ilay,ibox)
            fluxsum(ilay,ibox) = fluxsum(ilay,ibox)/area/h
         enddo 
      enddo
c
c     ---------- write fluxsum (also for dry layers, in case) ----------
c    
      write(*,871) "ibox","ilay","dV/dt/V", "errorlev"
      do ibox = 0, nboxes-1
         do ilay = 1, Boxes(ibox)%nwet
             write(*,872) ibox, ilay, fluxsum(ilay,ibox), 
     +                   absflxs(ilay,ibox)
         enddo 
      enddo
 871  format(a6,1x,a6,1x,a12,  1x,a12)
 872  format(i6,1x,i6,1x,f12.7,1x,f12.7)
c
c     --------- clean up ----------
c
      deallocate( fluxsum )
      deallocate( absflxs )
c
      end subroutine check_mass_conservation






      subroutine apply_top_down_mass_cons(hflux, vflux)
c     ----------------------------------------------------------------------------------
c     Enforce mass conservation by defining vertical fluxes from sums of 
c     horizontal fluxes into each (box,layer), starting from the surface and working
c     down toward the bottom. The vertical fluxes through the bottom of each column
c     are NOT redefined (so there may be an residual in the bottom box, indicating
c     sea level fluctuations). Also this preserves fluxes between open bottoms and box0
c     Fluxes as well as integrated fluxes can be supplied
c
c     hflux is positive from right to left box of faces
c     Positive direction for vflux across faces is upward;
c     vflux(i,j) is the flux from layer i+1 to layer i (offset 1) in box j (offset 0)
c     ----------------------------------------------------------------------------------
      real, intent(in)            :: hflux(:,0:) ! [m3/s] shape = (1:nlayers, 0:nfaces-1)
      real, intent(inout)         :: vflux(:,0:) ! [m3/s] shape = (1:nlayers, 0:nboxes-1)

      integer                     :: ibox, ilay, ipt, i
      integer                     :: nboxes, nlayers, nfaces
      integer                     :: istat, ifc, nscan
      integer                     :: lbox, rbox
      type(AtlantisFace), pointer :: face
      type(AtlantisBox), pointer  :: left, right, box
      real,allocatable            :: fluxsum(:,:)
c     ----------------------------------------------------------------------------------
      nboxes  = size(Boxes)
      nlayers = size(zgrid_CC) 
      nfaces  = size(Faces)
      allocate( fluxsum(1:nlayers, 0:nboxes-1) )
      fluxsum  = 0.     ! set pad value for dry layers
c
c     ---------- First scan horizontal fluxes ----------
c                fluxsum(ilay,ibox) is mass transport into (ilay,ibox)
c
      do ifc = 0, nfaces-1
         face  => Faces(ifc)
         lbox  = face%left_right(1)
         rbox  = face%left_right(2)
         left  => Boxes(lbox)
         right => Boxes(rbox)  
         nscan  = min(left%nwet, right%nwet)
c        ...... fluxsum(ilay,ibox) is flux into (ilay,ibox) ......
         do ilay = 1, nscan  ! top-down
            fluxsum(ilay,lbox) = fluxsum(ilay,lbox) + hflux(ilay,ifc)  ! flux positive from right to left
            fluxsum(ilay,rbox) = fluxsum(ilay,rbox) - hflux(ilay,ifc)  ! flux positive from right to left
         enddo ! ilay
      enddo    ! ifc
c
c     ---------- Define vertical fluxes top-down in each box column ----------
c                Do not redefied flux through floor of bottom layer
c
      do ibox = 0, nboxes-1
         box => Boxes(ibox)
         do ilay = 1, box%nwet-1  ! top-down, vflux(i,j) is the flux from layer i+1 to layer i, skip floor of bottom layer
            vflux(ilay,ibox)     = -fluxsum(ilay,ibox) 
            fluxsum(ilay+1,ibox) = fluxsum(ilay+1,ibox)-vflux(ilay,ibox)
         enddo ! ilay
      enddo    ! ibox
c
c     --------- clean up ----------
c
      deallocate( fluxsum )
c
      end subroutine apply_top_down_mass_cons
