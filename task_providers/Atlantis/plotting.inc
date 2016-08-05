cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Include functionality for module atlantis_grid
c
c     xmgrace ../task_providers/Atlantis/tmp/box?.xy ../task_providers/Atlantis/tmp/box??.xy -settype  xyvmap  hflux_plot.dat
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine plot_horizontal_fluxes(hflux, layer, iunit, maxlen)
c     ----------------------------------------------------------------------------------
c     Plot horizontal flux vectors of layer as (lon, lat, flux_u, flux_v) where (lon, lat) 
c     from flux map hflux. Fluxes are palced on the mid point of the face, and flux is 
c     multipled onto the face normal. Data is written to logical unit iunit, and fluxes
c     are scaled by some scale factor so that the largest flux has length maxlen (as a lon-lat step)
c     for plotting convenience. iunit is not closed at the end
c     nor bound to a file in this subroutine - this allows overlays, if desired.
c     ----------------------------------------------------------------------------------
      real, intent(in)            :: hflux(:,0:) ! [m3/s] shape = (1:nlayers, 0:nfaces-1)
      integer, intent(in)         :: layer, iunit
      real, intent(in)            :: maxlen
c
      integer                     :: nfaces, ifc
      type(AtlantisFace), pointer :: face
      real                        :: xy(2), uv(2), scale, fmax

c     ----------------------------------------------------------------------------------
      nfaces  = size(Faces)
c
c     ---------- scan horizontal fluxes ----------
c
      fmax = maxval(abs(hflux(layer,:)))
      if (fmax<1e-12) then
         write(*,*) "plot_horizontal_fluxes: fluxes are zero"
         stop 755
      else
         scale = maxlen/fmax
      endif
c
      do ifc = 0, nfaces-1
         face  => Faces(ifc)
         xy = 0.5*(face%nodes(:,1) + face%nodes(:,2))   ! face mid point
         uv = hflux(layer, ifc)*face%normal(1:2)        ! flux along face normal
         write(iunit,528) xy, uv*scale
      enddo    ! ifc

 528  format(2f12.7, 2x, 2f12.7)
 
      end subroutine plot_horizontal_fluxes