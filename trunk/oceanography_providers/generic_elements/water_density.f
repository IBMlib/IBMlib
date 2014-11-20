      module water_density
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Plugin to evaluate water density from depth, salinity and temperature,
c     based on UNESCO relation. Relies on directly on mesh_grid arrays 
c     Can be used directly around physical_fields. Based on pre-IBMlib version
c     NSParticleTracking/water_pressure.f
c     
c     Optimize automatic update overhead by monitoring data revision tag provided
c     by mesh_grid
c
c     Values of pressure and rhow in dry points are not guarantied
c
c     12Nov2014: rho_UNESCO moved to separate file, to allow usage independent link without need for build this module build
c                currently placed in physical_oceanography/rho_UNESCO.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh_grid
      implicit none
      private

      real,allocatable  :: rhow(:,:,:)      ! nx,ny,nz [kg/m**3] density at center of cell 
      real,allocatable  :: pressure(:,:,:)  ! nx,ny,nz [bar]     pressure at center of cell 

      integer(kind=selected_int_kind(16)) :: last_update = -1      ! track data revision tag in mesh_grid; last_update<0 means unset
      logical                             :: auto_update = .false. ! enable automatic update before interpolation

      public init_water_density
      public close_water_density
      public interpolate_rhow               ! with automitic update, if enabled
      public update_water_density           ! force update
      

      contains


      subroutine init_water_density(autoupd)
c     ---------------------------------------------------
c     nx,ny,nz from mesh_grid
c     ---------------------------------------------------
      logical, intent(in), optional :: autoupd
c     ---------------------------------------------------
      write(*,*) "init_water_pressure: allocating pressure + rhow"
      if (present(autoupd)) then
         auto_update = autoupd
      else
         auto_update = .false. ! default
      endif
      
      if (auto_update) then
         write(*,*) "init_water_pressure: automatic update before "//
     +              "interpolation"
      else
         write(*,*) "init_water_pressure: manual update mode"
      endif
      if (.not.allocated(pressure)) allocate( pressure(nx,ny,nz) )
      if (.not.allocated(rhow))     allocate( rhow(nx,ny,nz)     )
      last_update = -1 
      end subroutine init_water_density



      subroutine close_water_density()
c     ---------------------------------------------------
c     ---------------------------------------------------
      write(*,*) "close_turbulenc(): deallocating module data"
      if (allocated(pressure)) deallocate(pressure)
      if (allocated(rhow))     deallocate(rhow)
      last_update = -1 
      end subroutine close_water_density



      subroutine interpolate_rhow  (geo, r, status) 
c     ------------------------------------------ 
c     Do not accept data revision tag < 0 from mesh_grid
c     in auto_update mode - most likely this means that
c     update_data_rev_tag_mesh_grid has not been invoked 
c     by mesh_grid client
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)    ! (lon,lat,depth)
      real, intent(out)    :: r
      integer, intent(out) :: status

      integer(kind=selected_int_kind(16)) :: data_rev 
c     ------------------------------------------ 
      if (auto_update) then
         call get_data_rev_tag_mesh_grid(data_rev)
         if (data_rev < 0) then
            write(*,*) "data revision tag from mesh_grid < 0"
            stop
         endif
         if (data_rev /= last_update) call update_water_density() ! syncronizes last_update to data_rev
      endif
c     finally set (r,status) using interpolate_cc_3Dgrid_data(geo,array,deriv,padval,result,status)
c
      call interpolate_cc_3Dgrid_data(geo, rhow, 0, 0.0, r, status) 
c     ------------------------------------------ 
      end subroutine interpolate_rhow  




      subroutine update_water_density()
c-----------------------------------------------------------------------
c
c     Updates module arrays (rhow, pressure) iteratively to selfconsistency  
c     Water density and pressure are determined by two coupled, non-linear
c     relations: 
c          1) dp   = g rhow(p) dz
c          2) rhow = rho_UNESCO(s,T,p)
c
c     (rhow, p=pressure) are cell-centered and considered cell-wise constants and therefore
c     eqs. 1,2 reads (iz=2,...)
c
c          3)  p(iz)    = p(iz-1) + (g/2) * (h(iz-1)*rhow(iz-1) + h(iz)*rhow(iz) )
c          4)  rhow(iz) = rho_UNESCO(s,T,p(iz))
c
c     3,4 are iterated to selfconsistency for each layer, after which it is solved 
c     for next layer below, iz+1. The iterative loop usually converges to numerical 
c     precision within 3-4 iteratative cycles of eqs 3,4
c 
c     When arrays (rhow, pressure) are updated, set the data revision tag last_update correspondingly,
c     but do not compare prior to update (i.e. forceful update)
c
c     Notes:
c       1 bar = 100 000 pascals (Pa) 
c       Assume atmospheric pressure is constant patm = 101.325 kPa = 1.01325 bar
c     
c-----------------------------------------------------------------------
      real            :: h(nz)        ! h(iz) = layer width of layer iz in meters, incl. DSLM  (automatic array)
      real, parameter :: g     =  9.81d00    ! [meter/sec**2] added, asbjorn
      real, parameter :: gfac  =  0.5*g/1.e5 ! g/2 * conversion from Pa to bar
      real, parameter :: patm  =  1.01325d0  ! [bar]          constant atmospheric pressure
      real, parameter :: rhow0 =  1.027d3    ! [kg/m**3]
      real, parameter :: dplim    =  1.e-12  ! [bar] break condition for selfconsistency loop
      real, parameter :: dpaccept =  1.e-5   ! [bar] tolerance limit: due to finite float rep, SC loop may very rarely get trapped  
      integer, parameter :: nscmax = 10      ! max number of iterations per layer for eqs 3+4 above
      logical, parameter :: debug  = .false. ! write state variables to stdout 
      integer         :: ix,iy,iz,im
      real            :: dp, plast, s, t
c-----------------------------------------------------------------------    
      pressure = patm  ! initial condition on all cells     
      rhow     = rhow0 ! initial condition on all cells
     
      if (debug) write(*,458)
      do iy = 1, ny
         do ix = 1, nx
            if (bottom_layer(ix,iy) == 0) cycle ! this point is dry, go to next
            h = acc_width(ix,iy,2:nz+1) - acc_width(ix,iy,1:nz) ! h(iz) = layer width of layer iz
            iz = 1
            dp = 2*dplim         ! do at least one iteration
            s  = salinity(ix,iy,iz) ! bottom_layer(ix,iy) >= 1
            t  = temp(ix,iy,iz)  ! bottom_layer(ix,iy) >= 1

            do im=1,nscmax
               plast             = pressure(ix,iy,iz)
               pressure(ix,iy,iz)= patm+gfac*h(iz)*rhow(ix,iy,iz) ! gfac = g/2
               rhow(ix,iy,iz)    = rho_UNESCO(s,t,pressure(ix,iy,iz))
               dp                = abs(pressure(ix,iy,iz)-plast)
               if (debug) write(*,459) ix,iy,iz,im,
     +                      rhow(ix,iy,iz),pressure(ix,iy,iz)  
               if (dp<dplim) exit
            enddo
            if (dp>dpaccept) then
               write(*,*) "ix,iy,iz = ",ix,iy,iz
               write(*,*) "dp = ", dp, " < dpaccept = ",dpaccept
               stop "SC rho/p iterations not converged"
            endif

            do iz = 2, bottom_layer(ix,iy)  ! downward scan; loop void if bottom_layer(ix,iy) <= 1
               dp = 2*dplim         ! do at least one iteration
               s  = salinity(ix,iy,iz) ! 
               t  = temp(ix,iy,iz)  ! 
 
               do im=1,nscmax
                  plast             = pressure(ix,iy,iz)
                  pressure(ix,iy,iz)= pressure(ix,iy,iz-1) +
     +                 gfac*(h(iz-1)*rhow(ix,iy,iz-1) +    ! gfac = g/2
     +                       h(iz  )*rhow(ix,iy,iz))
                  rhow(ix,iy,iz)    = rho_UNESCO(s,t,pressure(ix,iy,iz))
                  dp                = abs(pressure(ix,iy,iz)-plast)
                  if (debug) write(*,459) ix,iy,iz,im,
     +                         rhow(ix,iy,iz),pressure(ix,iy,iz) 
                  if (dp<dplim) exit
               enddo
               if (dp>dpaccept) then
                  write(*,*) "ix,iy,iz = ",ix,iy,iz
                  write(*,*) "dp = ", dp, " < dpaccept = ",dpaccept
                  stop "SC rho/p iterations not converged"
               endif

            enddo               ! iz = 2, ...
         enddo                  ! ix = 1, ...
      enddo                     ! iy = 1, ...

      call get_data_rev_tag_mesh_grid(last_update)  ! syncronize data revision tag
      
c      ix=60
c      iy=80
c      do iz=1,bottom_layer(ix,iy)
c      write(81,*)iz,salinity(ix,iy,iz),temp(ix,iy,iz),pressure(ix,iy,iz),
c     +           rhow(ix,iy,iz)
c      enddo
c      STOP 777


 458  format(2x,"ix",2x,"iy",2x,"iz",2x,"im",3x,
     +       "rhow[kg/m**3]",1x,"pressure[bar]")
 459  format(4i4, 2x, 2f12.5)
      end subroutine update_water_density

      end module 


c$$$      program test_module
c$$$c     -----------------------------------------------
c$$$      use water_pressure
c$$$      implicit none
c$$$      integer, parameter :: nx=2, ny=2, nz=20
c$$$      integer   :: botlayer(nx,ny)
c$$$      real      :: salty(nx,ny,nz), temp(nx,ny,nz), h(nx,ny,nz)
c$$$c     -----------------------------------------------
c$$$      call init_water_pressure(nx,ny,nz)
c$$$      botlayer(1,1) = 0
c$$$      botlayer(1,2) = 1
c$$$      botlayer(2,1) = 10
c$$$      botlayer(2,2) = nz
c$$$      call random_number(salty)
c$$$      salty =  30.d0 + 10.d0 * salty
c$$$      call random_number(temp)
c$$$      temp  = 15.d0 * temp
c$$$      call random_number(h(:,:,1))
c$$$      h(:,:,1)     = 4.d0 + 4.d0 * h(:,:,1)
c$$$      h(:,:,2:4)   = 8.d0 
c$$$      h(:,:,5:10)  = 1.6d1
c$$$      h(:,:,11:15) = 5.0d1 
c$$$      h(:,:,16:)   = 1.0d2 
c$$$      call update_rhow(salty, temp, botlayer, h)
c$$$      call close_water_pressure()
c$$$
c$$$      end program
