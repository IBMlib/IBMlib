ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Updates the arrays
c
c        vertical diffusivity    vdiffus(:,:,:)   ! ix,iy,iz   [m**2/s]
c        horizontal diffusivity  hdiffus(:,:,:)   ! ix,iy,iz   [m**2/s]
c     
c     exported to physical_fields. Standalone verion that can be applied with regular lon-lat data grids
c     for posterior evaluation of vertical/horizontal diffusivities (undeveloped module from NSParticletracking)
c
c     Do not dublicate array dimensions (nx,ny,nz) as explicit module data, but
c     perform query at the entry of relevant methods
c     The module current also holds and updates ocean surface stress taus(nx,ny) 
c     is private attribute (may be moved to another module later)
c
c     Provides:
c       init_particle_state:                    parameters from input file
c       update_turbulence(...):                 arguments specific to each provider
c       interpolate_turbulence(pos, k)          pos = scaled coor    k( 3) == turbulence at pos
c       interpolate_turbulence_deriv(pos, dk)   pos = scaled coor    dk(3) == turbulence deriv at pos
c       close_turbulence:                       civilized module close-down 
c
c     TODO:
c        alternative posteori turbulence algorithms from GOTM  (www.gotm.net)
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module turbulence
      use input_parser
      use run_context
      use constants
      implicit none
      private                   ! default visibility

c --------------------------------------------------------
c ---------------- module data section -------------------
c --------------------------------------------------------
c
c vdiffus(:,:,2) is situated between layers 1,2 etc, shape(vdiffus) = (nx,ny,nz+1)
c vdiffus(:,:,:) is face centered on upper cell face
c i.e. vdiffus(ix,iy,iz) corresponds to (x,y,z) = (ix,iy,iz-0.5)
c hdiffus is cell centered, i.e. hdiffus(ix,iy,iz) corresponds to (x,y,z) = (ix,iy,iz)
c interior data ranges:
c   vdiffus(2:nx,   1:ny-1, 1:botlay+1) 
c   hdiffus(2:nx-1, 2:ny-1, 1:botlay) 

      real, allocatable, public    :: vdiffus(:,:,:) ! nx,ny,nz+1  face centered vertical   diffusivity [m**2/s]                               
      real, allocatable, public    :: hdiffus(:,:,:) ! nx,ny,nz    cell centered horizontal diffusivity [m**2/s]
      real, allocatable            :: taus(:,:)      ! nx,ny       cell centered surface stress         [N/m**2]

      character*4                  :: vertical_diffusivity_update     ! "none" / "cons" / "ppmx" / "pohl" / "bshm"
      character*4                  :: horizontal_diffusivity_update   ! "none" / "cons" / "smag"
      real                         :: vertical_diffusivity_constant   ! (only defined for vertical_diffusivity_update   == "cons")
      real                         :: horizontal_diffusivity_constant ! (only defined for horizontal_diffusivity_update == "cons")
      character*4                  :: vertical_diffusivity_BC         ! "diri" (vdiffus=0) / "nman" (dvdiffus/dz=0) 
c     
c     set compile time defaults
c
      real                         :: smagorinsky_constant     =  0.31    ! 0.31 (effective small viscosity DMI setting)
      real                         :: smagorinsky_schmidtnum   =  1.0     ! 1.0 
      real                         :: smagorinsky_eddyvisc_max =  0.1     ! unit == m2/sec (suggested NS value 0.1 m2/sec)
c
c     notation as J. Phys. Ocean. 11(1981) p1450, but in SI units
c
      real                         :: ppmix_vb = 1.0e-4   ! m2/sec        (background eddy viscousity) 
      real                         :: ppmix_kb = 1.0e-5   ! m2/sec        (background eddy diffusivity)
      real                         :: ppmix_n  = 2.0      ! dimensionless (denominator power)
      real                         :: ppmix_a  = 5.0      ! dimensionless (Richardson number prefactor)
      real                         :: ppmix_v0 = 50.0e-4  ! m2/sec        (eddy viscousity numerator)
c.....set public calling interface

      public init_turbulence
      public close_turbulence
      public update_turbulence
      public interpolate_turbulence
      public interpolate_turbulence_deriv

c.....ONLY temp public      
      public taus
      public update_taus


      contains
 

      subroutine init_turbulence(nx,ny,nz)
c     ---------------------------------------------------
c     nx,ny,nz from physical_fields are not accessible via 
c     use association - do not dublicate to module data
c     
c     look for input tage: 
c         vertical_diffusivity_scheme
c         horizontal_diffusivity_scheme
c         (vertical_diffusivity_constant)
c         (horizontal_diffusivity_constant)
c     
c     v/h_diffusivity_scheme == none:     no update of v/h diffusivity
c     v/h_diffusivity_scheme == constant: set v/h diffusivity to v/h_diffusivity_constant
c
c     v/h_diffusivity_scheme is not present: use intrinsic update (eddyvis_v/h)
c     ---------------------------------------------------   
      integer, intent(in) :: nx,ny,nz
      character*256       :: chr_option
      real                :: diffus_val
      integer             :: ntags
c     ---------------------------------------------------  
      write(*,*) "init_turbulence: allocating"   
      allocate( vdiffus(nx,ny,nz+1) ) 
      allocate( hdiffus(nx,ny,nz) ) 
      allocate( taus(nx,ny)       ) 
c
c     ----------------------------------------
c     ------  resolve vertical update   ------
c     ----------------------------------------
c
c     set default, if vertical_diffusivity_scheme is not provided
c
      ntags = count_tags(simulation_file, "vertical_diffusivity_scheme")
      if (ntags == 1) then     
         call read_control_data(simulation_file, 
     +                       "vertical_diffusivity_scheme", chr_option)
         chr_option = adjustl(chr_option)
         write(*,229) "init_turbulence: vertical_diffusivity_scheme = ",
     +                trim(chr_option)
         if     (chr_option(1:4) == "none") then
            vertical_diffusivity_update = "none"

         elseif (chr_option(1:8) == "constant") then
            vertical_diffusivity_update = "cons"
            call read_control_data(simulation_file, 
     +           "vertical_diffusivity_constant", diffus_val)
            vertical_diffusivity_constant = diffus_val
            write(*,231) "init_turbulence: vdiffus =", diffus_val

         elseif (chr_option(1:8) == "pohlmann") then
            vertical_diffusivity_update = "pohl" 
            
         elseif (chr_option(1:6) == "bshmod") then
            vertical_diffusivity_update = "bshm" 
         
         elseif (chr_option(1:5) == "ppmix") then
            vertical_diffusivity_update = "ppmx" 
           
         else
            write(*,229) "init_turbulence: unknown " //
     +                   "vertical_diffusivity_scheme = ",
     +                   trim(chr_option)
            stop ! fatal
         endif

      else   ! vertical_diffusivity_scheme not provided, set default
             ! Pacanowski/Philander scheme
         vertical_diffusivity_update = "ppmx"
         write(*,229) "init_turbulence: using "//
     +                "intrinsic vertical diffusivity scheme", 
     +                 vertical_diffusivity_update
      endif 
c
c     ----------------------------------------
c     ------ resolve horizontal update  ------
c     ----------------------------------------
c
c     set default, if horizontal_diffusivity_scheme is not provided
c
      ntags=count_tags(simulation_file,"horizontal_diffusivity_scheme")
      if (ntags == 1) then     
         call read_control_data(simulation_file, 
     +                     "horizontal_diffusivity_scheme", chr_option)
         chr_option = adjustl(chr_option)
         write(*,229) "init_turbulence: "//
     +                "horizontal_diffusivity_scheme = ",
     +                trim(chr_option)

         if     (chr_option(1:4) == "none") then
            horizontal_diffusivity_update = "none"

         elseif (chr_option(1:8) == "constant") then
            horizontal_diffusivity_update = "cons"
            call read_control_data(simulation_file, 
     +           "horizontal_diffusivity_constant", diffus_val)
            horizontal_diffusivity_constant = diffus_val
            write(*,231) "init_turbulence: hdiffus =", diffus_val

         elseif (chr_option(1:12) == "smagorinsky") then
            horizontal_diffusivity_update = "smag"

            call read_control_data(simulation_file, 
     +         "smagorinsky_constant", smagorinsky_constant)
            write(*,231) "smag constant", smagorinsky_constant

            call read_control_data(simulation_file, 
     +         "smagorinsky_schmidtnum", smagorinsky_schmidtnum)
            write(*,231) "smag schmidtnum", smagorinsky_schmidtnum

            call read_control_data(simulation_file, 
     +         "smagorinsky_eddyvisc_max", smagorinsky_eddyvisc_max)
            write(*,231) "smag eddyvisc_max", smagorinsky_eddyvisc_max

         else
            write(*,229) "init_turbulence: unknown " //
     +                   "horizontal_diffusivity_scheme = ",
     +                   trim(chr_option)
            stop   ! fatal
         endif

      else  ! horizontal_diffusivity_scheme not provided, set default
         horizontal_diffusivity_update = "smag"
         write(*,230) "init_turbulence: using "//
     +                "intrinsic horizontal diffusivity scheme=",
     +                 horizontal_diffusivity_update 
c        --- set compile time smagorinsky defaults (ignore potential apperances in input)
         write(*,231) "smag constant     =", smagorinsky_constant
         write(*,231) "smag schmidtnum   =", smagorinsky_schmidtnum
         write(*,231) "smag eddyvisc_max =", smagorinsky_eddyvisc_max

      endif 
c
c     ----------------------------------------
c     ---- resolve turbulence BC options  ----
c     ---------------------------------------- 
      
      ntags=count_tags(simulation_file,"vertical_diffusivity_BC")
      if (ntags == 1) then  
         call read_control_data(simulation_file, 
     +                     "vertical_diffusivity_BC", chr_option)
         chr_option = adjustl(chr_option)
         write(*,229) "init_turbulence: "//
     +                "vertical_diffusivity_BC = ",
     +                trim(chr_option)

         if     (chr_option(1:9) == "dirichlet") then
            vertical_diffusivity_BC  = "diri"

         elseif (chr_option(1:7) == "neumann") then
            vertical_diffusivity_BC  = "nman"
         else
            write(*,229) "init_turbulence: unknown " //
     +                   "vertical_diffusivity_BC = ",
     +                   trim(chr_option)
            stop   ! fatal
         endif

      else  ! vertical_diffusivity_BC not provided, set default
         vertical_diffusivity_BC = "diri"
         write(*,229) "init_turbulence: using "//
     +                "default vertical diffusivity BC = ",
     +                 vertical_diffusivity_BC   
      endif 
      
 229  format(a,1x,a)
 230  format(a)
 231  format(a,1x,f12.7)

      end subroutine init_turbulence



      subroutine close_turbulence()
c     ---------------------------------------------------
      write(*,*) "close_turbulence: deallocating"
      if (allocated(vdiffus)) deallocate(vdiffus)
      if (allocated(hdiffus)) deallocate(hdiffus)
      if (allocated(taus))    deallocate(taus)
      end subroutine close_turbulence


      subroutine update_turbulence(u, v, uswind, vswind, rho, h, depth,
     +                             botlayer, dth, xfacelen, yfacelen)   
c     -----------------------------------------------------------------------------
c     Update horizontal/vertical turbulence fields from present physical fields
c     
c     Note:
c         In principle, depth(ix,iy) = sum( h(ix,iy,iz), iz=1, botlayer(ix,iy) )
c         This redundancy in arguments has been maintained for convenience
c     -----------------------------------------------------------------------------
      real, intent(in)    :: u(:,:,:)      ! W:E water current (assumed shape array)
      real, intent(in)    :: v(:,:,:)      ! S:N water current (assumed shape array)
      real, intent(in)    :: uswind(:,:)   ! W:E surface wind  (assumed shape array)
      real, intent(in)    :: vswind(:,:)   ! S:N surface wind  (assumed shape array) 
      real, intent(in)    :: rho(:,:,:)    ! water density     (assumed shape array)
      real, intent(in)    :: h(:,:,:)      ! layer thickness   (assumed shape array)
      real, intent(in)    :: depth(:,:)    ! local water depth (assumed shape array)
      integer, intent(in) :: botlayer(:,:) ! bottom layer      (assumed shape array)
      real, intent(in)    :: dth           ! hydrodynamic time step
      real, intent(in)    :: xfacelen(:)   ! WE cell face length (assumed shape array)
      real, intent(in)    :: yfacelen      ! SN cell face length
c     ---------------------------------------------------
      write(*,*) "update_turbulence: vertica/horizontal turbulence"
      call update_taus(uswind, vswind)

c
c     ------ select vertical update ------
c
      if     (vertical_diffusivity_update == "none") then
         write(*,*) "update_turbulence: no vdiffus update"
      elseif (vertical_diffusivity_update == "cons") then
         write(*,*) "update_turbulence: set constant vdiffus"
         vdiffus = vertical_diffusivity_constant
      elseif (vertical_diffusivity_update == "ppmx") then
         write(*,*) "update_turbulence: call ppmix_vertical_viscosity"
         call ppmix_vertical_viscosity(u, v, rho, h, botlayer)
      elseif (vertical_diffusivity_update == "pohl") then
         write(*,*) "update_turbulence: call hamsom_vertical_viscosity"
         call hamsom_vertical_viscosity(u, v, rho, h, botlayer)
      elseif (vertical_diffusivity_update == "bshm") then
         write(*,*) "update_turbulence: call eddyvis_vertical"
         call eddyvis_vertical(u, v, rho, h, depth, botlayer, dth)         
      else
         write(*,*) "update_turbulence: " //
     +              "unhandled vertical_diffusivity_update = ",
     +              vertical_diffusivity_update
         stop
      endif

      call apply_vdiffus_BC(botlayer)
c
c     ------ select horizontal update ------
c
      if     (horizontal_diffusivity_update == "none") then
         write(*,*) "update_turbulence: no hdiffus update"
      elseif (horizontal_diffusivity_update == "cons") then
         write(*,*) "update_turbulence: set constant hdiffus"
         hdiffus = horizontal_diffusivity_constant
      elseif (horizontal_diffusivity_update == "smag") then
         write(*,*) "update_turbulence: call eddyvis_horizontal"
         call eddyvis_horizontal(u,v, botlayer, xfacelen,yfacelen)
      else
         write(*,*) "update_turbulence: " //
     +              "unhandled horizontal_diffusivity_update",
     +              horizontal_diffusivity_update
         stop
      endif
      
      call apply_hdiffus_BC()
c
c
      end subroutine update_turbulence

      
   
      subroutine update_taus(u_wind, v_wind)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update surface wind stress taus(nx,ny) (currently module data)
c     from current surface wind speed
c     taus unit = N/m**2
c     Update based on heuristic expression from DMI (Dec. 13, 2006)
c     Assign sea/air wind stress irrespective of whether point is wet or dry (i.e. no wet check)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in)  :: u_wind(:,:)          ! ix,iy,  surface W:E wind speed  [meter/second]
      real, intent(in)  :: v_wind(:,:)          ! ix,iy,  surface S:N wind speed  [meter/second]
c     ------------------------------------------------------------------
      real,allocatable  :: ws2(:,:) ! ix,iy,  wind speed, automatic array
      real, parameter   :: rhoair   =  1.2d0    ! [kg/m3]
c     ------------------------------------------------------------------
      allocate( ws2(size(u_wind,1),size(u_wind,2)) ) ! same shape as u_wind
      ws2  = u_wind**2 + v_wind**2                           
      taus = 1.d-3 * rhoair * (0.5 + 0.07*sqrt(ws2)) * ws2 ! vectorized expression
      deallocate(ws2)
      end subroutine update_taus

      
      subroutine apply_vdiffus_BC(botlayer)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Apply boundary conditions on vertical diffusivity vdiffus
c
c     Currently only BC in vertical direction, because horizontal edges are open
c     Options:
c           "diri" Dirichlet BC    (vdiffus=0)
c           "nman" Neumann   BC    (d(vdiffus)/dz=0) 
c
c     Assume interior part of vdiffus has been updated
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer, intent(in) :: botlayer(:,:) ! bottom layer      (assumed shape array)
      integer             :: ix,iy,iz,kb,nx,ny
c     --------------------------------------------------------------------------------
      nx      = size(botlayer,1)     
      ny      = size(botlayer,2)     

      if (vertical_diffusivity_BC == "diri") then
         do ix = 1, nx
            do iy = 1, ny  
               kb = botlayer(ix,iy)
               if (kb>0) then   ! i.e. wet point
                  vdiffus(ix,iy,1)    = 0.0d0 ! OK for single layer (kb=1)
                  vdiffus(ix,iy,kb+1) = 0.0d0 ! OK for single layer (kb=1)
               endif
            enddo
         enddo

      elseif (vertical_diffusivity_BC == "nman") then
         do ix = 1, nx
            do iy = 1, ny  
               kb = botlayer(ix,iy)
               if (kb>0) then   ! i.e. wet point             
                  vdiffus(ix,iy,1)    = vdiffus(ix,iy,2)  ! OK for single layer (kb=1)
                  vdiffus(ix,iy,kb+1) = vdiffus(ix,iy,kb) ! OK for single layer (kb=1)               
               endif
            enddo
         enddo

      else
         write(*,197) vertical_diffusivity_BC
         stop
      endif

 197  format("unknown vertical_diffusivity_BC option: ", a)

      end subroutine apply_vdiffus_BC



      subroutine apply_hdiffus_BC()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Apply boundary conditions on horizontal diffusivity hdiffus
c
c     Assume interior part of hdiffus has been updated
c     currently void, but keep for future elaborations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end subroutine apply_hdiffus_BC




      subroutine ppmix_vertical_viscosity(u, v, rho, h, botlayer)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Provide simplified vertical eddy viscosity/diffusivity scheme of 
c     Pacanowski and Philander (J. Phys. Ocean. 11(1981) p1450)
c
c     Output: set module data vdiffus(ix,iy,iz) == vertical diffusivity on tracers [meter**2/second] 
c
c     vdiffus(:,:,2) is situated between layers 1,2 etc, 
c     i.e. vdiffus(:,:,iz) corresponds to z = iz-0.5
c     vdiffus(:,:,1) and vdiffus(bottom_layer+1) are set by boundary conditions
c
c     For single-layer cells (botlayer=1), vdiffus(:,:,1:2) = vdiffus_shallow
c     where vdiffus_shallow is a constant (for now)    
c         
c
c     Units: SI     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

      real, intent(in)    :: u(:,:,:)      ! nx,ny,nz  W:E current   [meter/second]
      real, intent(in)    :: v(:,:,:)      ! nx,ny,nz  S:N current   [meter/second]
      real, intent(in)    :: rho(:,:,:)    ! nx,ny,nz  water density [kg/m3]
      real, intent(in)    :: h(:,:,:)      ! nx,ny,nz  actual layer thickness (including surface elevation) [meter]
      integer, intent(in) :: botlayer(:,:) ! nx,ny     bottom layer thickness in this water column
      
      real, parameter   :: vdiffus_shallow   = 1.45d0   ! diffusion length per day == 500 m
 
      integer       :: ix, iy, iz, nx, ny, nz, kb, izu, izl,last
      real          :: du_dz2, dv_dz2, drho_dz, hcen, arg 
      real          :: h_umx, h_lmx, BVfreq2, ricut, ppdenom, rhomid
      real          :: rich(size(u,3))       ! Richadrson numbers for water column
      real          :: eddyvis(size(u,3))    ! eddy viscosity [meter**2/second                    
c     ------------------------------
      vdiffus = 0.            ! set pad value (for undefined cells)
      nx      = size(u,1)     ! assume u,v matches
      ny      = size(u,2)     ! assume u,v matches
      nz      = size(u,3)     ! assume u,v matches

      do ix = 2, nx
         do iy = 1, ny-1            
            rich    = 0.0d0
c           ---------------------------------------------------
c           Assign vdiffus(ix,iy,:) for column (ix,iy)
c           --------------------------------------------------- 
            kb = botlayer(ix,iy)      

            if     (kb == 0) then                     ! do nothing for kb==0 (this is a dry point, vdiffus=0)
               cycle
            elseif (kb == 1) then                     ! shallow area limit
               vdiffus(ix,iy,1:2) = vdiffus_shallow   !
            elseif (kb  > 1) then  
c              ------------------------------------------------------------------          
c              loop over interfaces in water column with at least two layers 
c              ------------------------------------------------------------------   
               eddyvis = 0.0
               do iz = 2, kb                          ! 
                  hcen    = 0.5*(h(ix,iy,iz) + h(ix,iy,iz-1)) ! distance between layer centers
                  drho_dz = (rho(ix,iy,iz) - rho(ix,iy,iz-1)) / hcen    ! positive direction down  
                  rhomid  = 0.5*(rho(ix,iy,iz) + rho(ix,iy,iz-1))
                  du_dz2  = ((u(ix,iy,iz)  + u(ix-1,iy,iz) 
     +                    -  u(ix,iy,iz-1) - u(ix-1,iy,iz-1))
     +                       /(2.d0*hcen))**2
                  dv_dz2  = ((v(ix,iy,iz)  + v(ix,iy+1,iz) 
     +                    -  v(ix,iy,iz-1) - v(ix,iy+1,iz-1))
     +                       /(2.d0*hcen))**2
                  BVfreq2  =  gdrho*drho_dz   ! frequency positive for a stable water column
                  rich(iz) = BVfreq2 / (du_dz2 + dv_dz2 + 1.e-15) 
                  ricut    = max(-0.5/ppmix_a, rich(iz)) ! avoid singular denominator if water column unstable (rich<0)                
                  ppdenom  = (1.d0 + ppmix_a*ricut)        ! guarantied > 0
                  eddyvis(iz) = ppmix_v0/(ppdenom**ppmix_n) + ppmix_vb  ! Pacanowski and Philander Eq 1
                  vdiffus(ix,iy,iz) = eddyvis(iz)/ppdenom + ppmix_kb    ! Pacanowski and Philander Eq 2    
c
c                  write(23,63)ix,iy,iz,rich(iz),rhomid,vdiffus(ix,iy,iz) 
c 63              format(3i4,2f12.4,e12.5)
c                  
               enddo
               
c              ------------ print test block: begin------------
c               write(23,*) ix,iy,eddyvis(2)*1.0d4  ! cm**2/s
c               write(23,*) ix,iy,rich(2)
c               write(23,*) ix,iy,vforce(2)
c               write(81,*) ix,iy, sum(eddyvis(2:kb))/(kb-1)*1.0d4  ! cm**2/s
c               write(*,*) ix,iy, eddyvis(2:kb)
c               if ((rich(2) < 0.3).and.(kb>2)) write(81,*) ix,iy 
c               write(81,*) ix,iy, eddyvis(3)*1.0d4
c               if ((ix==40).and.(iy==40)) then
c                  do iz=1,kb
c               write(81,*)iz,eddyvis(iz)*1.0d4,rich(iz)  !,vdiffus(ix,iy,iz)
c                      write(82,*)iz,rich(iz)
c                     write(81,*)iz,rho(ix,iy,iz)
c                  enddo
c               endif
c              ------------ print test block: end------------

            endif

         enddo   ! do iy = 1, ny-1 
      enddo      ! do ix = 2, nx

c      do iz = 2, kb      
c      write(*,*)ix,iy,iz,rich(iz),eddyvis(iz),schmidt(iz),rho(ix,iy,iz)
c      enddo

      end subroutine ppmix_vertical_viscosity







      subroutine hamsom_vertical_viscosity(u, v, rho, h, botlayer)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Provide simplified vertical k-epsilon eddy viscosity/diffusivity
c     as implemented in HAMSOM as reported in 
c
c     Thomas Pohlmann, Cont. Shelf Res. 16(2), pp. 131-146, 1996
c     Thomas Pohlmann, Cont. Shelf Res. 16(2), pp. 147-161, 1996 
c     Corinna Schrum,  Cont. Shelf Res. 17(6), pp. 689-716, 1997
c
c     Output: set module data vdiffus(ix,iy,iz) == vertical diffusivity on tracers [meter**2/second] 
c
c     vdiffus(:,:,2) is situated between layers 1,2 etc, 
c     i.e. vdiffus(:,:,iz) corresponds to z = iz-0.5
c     vdiffus(:,:,1) and vdiffus(bottom_layer+1) are set by boundary conditions
c
c     For single-layer cells (botlayer=1), vdiffus(:,:,1:2) = vdiffus_shallow
c     where vdiffus_shallow is a constant (for now)    
c     
c     The layer structure of the water column is assumed to be 
c
c           [upper mixed layer] | [middle laminar layer] | [lower mixed layer]
c
c     where any of these may be absent. If middle laminar layer, upper/lower mixed
c     are merged into one upper mixed layer. The layer thickness is only resolved as
c     integer layers (i.e. the intralayer position of interface is not assessed by
c     e.g. interpolation (as Pohlmann) in this implementation). The middle laminar layer
c     is assigned molecular viscosity (1.34d-6 m**2/s)
c     Diffusivity of large molecules in quiescent fluids < 10^-10 m2/sec
c
c     Units: SI     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

      real, intent(in)    :: u(:,:,:)      ! nx,ny,nz  W:E current   [meter/second]
      real, intent(in)    :: v(:,:,:)      ! nx,ny,nz  S:N current   [meter/second]
      real, intent(in)    :: rho(:,:,:)    ! nx,ny,nz  water density [kg/m3]
      real, intent(in)    :: h(:,:,:)      ! nx,ny,nz  actual layer thickness (including surface elevation) [meter]
      integer, intent(in) :: botlayer(:,:) ! nx,ny     bottom layer thickness in this water column
      
      real, parameter   :: vdiffus_shallow   = 1.45d0   ! diffusion length per day == 500 m
      real, parameter   :: eddyvis_molecular = 1.34d-6  ! laminar water viscosity [meter**2/second]
      real, parameter   :: eddyvis_max       = 1000.e-4 ! heuristic upper limit from HAMSOM (1000 cm**2/s)
      real, parameter   :: rich_crit         = 0.22d0   ! Richardson number for laminar stability  (ri > ri_crit -> laminar) 
      real, parameter   :: c_ml              = 0.05d0  ! scale factor. Schrum calibration: 0.025  Pohlmann: 0.05
      real, parameter   :: ricut_min         = 3.e-3    ! heuristic min limit on Richardson number
      real, parameter   :: ricut_max         = 1.e6

      integer       :: ix, iy, iz, nx, ny, nz, kb, izu, izl,last
      real          :: du_dz2, dv_dz2, drho_dz, hcen, arg 
      real          :: h_umx, h_lmx, BVfreq2, ricut
      real          :: rich(size(u,3)), schmidt(size(u,3)) ! Richadrson and Schmidt numbers for water column
      real          :: eddyvis(size(u,3))                  ! eddy viscosity [meter**2/second]
      real          :: vforce(size(u,3))                           
c     ------------------------------
      vdiffus = 0.            ! set pad value (for undefined cells)
      nx      = size(u,1)     ! assume u,v matches
      ny      = size(u,2)     ! assume u,v matches
      nz      = size(u,3)     ! assume u,v matches

      do ix = 2, nx
         do iy = 1, ny-1            
            rich    = 0.0d0
            schmidt = 1.0d0
c           ---------------------------------------------------
c           Assign vdiffus(ix,iy,:) for column (ix,iy)
c           --------------------------------------------------- 
            kb = botlayer(ix,iy)      

            if     (kb == 0) then                     ! do nothing for kb==0 (this is a dry point, vdiffus=0)
               cycle
            elseif (kb == 1) then                     ! shallow area limit
               vdiffus(ix,iy,1:2) = vdiffus_shallow   !
            elseif (kb  > 1) then           
c              ------------------------------------------------------------------          
c              loop over interfaces in water column with at least two layers 
c              ------------------------------------------------------------------      
               do iz = 2, kb                          ! 
                  hcen    = 0.5*(h(ix,iy,iz) + h(ix,iy,iz-1)) ! distance between layer centers
                  drho_dz = (rho(ix,iy,iz) - rho(ix,iy,iz-1)) / hcen    ! positive direction down             
                  du_dz2  = ((u(ix,iy,iz)  + u(ix-1,iy,iz) 
     +                    -  u(ix,iy,iz-1) - u(ix-1,iy,iz-1))
     +                       /(2.d0*hcen))**2
                  dv_dz2  = ((v(ix,iy,iz)  + v(ix,iy+1,iz) 
     +                    -  v(ix,iy,iz-1) - v(ix,iy+1,iz-1))
     +                       /(2.d0*hcen))**2
                  BVfreq2  =  gdrho*drho_dz   ! frequency positive for a stable water column
                  rich(iz) = BVfreq2 / (du_dz2 + dv_dz2 + 1.e-15)  
                  ricut    = max(ricut_max, min(rich(iz), ricut_min)) ! limit before evaluating schmidt,vforce
                  arg      = ricut**2 - 0.316*ricut + 0.0346  ! arg >= 0.009636 
                  schmidt(iz) = ricut/(ricut + 0.186 - sqrt(arg))/0.725 

c                  write(23,*) BVfreq2
c 
c                 negative BVfreq2 (unstable stratification) increases viscosity
c
                  arg = max(1.e-15, du_dz2+dv_dz2 - BVfreq2/schmidt(iz))
                  vforce(iz) = sqrt(arg) !  vforce >= 0

               enddo
               
c              ------------------------------------------------------------------
c              Now the water column is characterized by hydrodynamic 
c              numbers (rich,schmidt). Now deduce the layer structure by the
c              rich(iz) >< rich_crit criterium
c              ------------------------------------------------------------------
               eddyvis = eddyvis_molecular       ! laminar middle (if any)
c
c              -- successively increase upper mixed layer downward, until rich(izu) > rich_crit
c                 at loop exit, izu is the first laminar interface
c              
               h_umx = 0.
               do izu = 2,kb
                  if (rich(izu) >  rich_crit) exit
                  h_umx = h_umx + h(ix,iy,izu-1) ! add a mixed layer
               enddo
               if (izu > 2) then                 ! i.e. non zero upper mixed layer
                  eddyvis(2:izu-1) = (c_ml*h_umx)**2*vforce(2:izu-1)
               endif
c
c              -- successively increase lower mixed layer upward, until rich(izl) > rich_crit
c                 at loop exit, izu is the first laminar interface
c                             
               h_lmx = 0.
               do izl = kb, izu, -1                ! izu >= 2
                  if (rich(izl) >  rich_crit) exit
                  h_lmx = h_lmx + h(ix,iy,izl) ! add a mixed layer
               enddo
               if (izl < kb) then              ! i.e. non zero lower mixed layer
                  eddyvis(izl+1:kb)=(c_ml*h_lmx)**2*vforce(izl+1:kb)
               endif
c
c              -- Apply heuristic upper/lower bounds on eddy viscosity
c
               where (eddyvis < eddyvis_molecular) 
                  eddyvis = eddyvis_molecular
               end where 
               where (eddyvis > eddyvis_max) 
                  eddyvis = eddyvis_max
               end where 

c
c              -- Now eddyvis are defined for the water column (ix,iy,layers >= 2) 
c                 finally generate eddy diffusivity by division with schmidt number
c
               vdiffus(ix,iy, 2:kb) = eddyvis(2:kb) / schmidt(2:kb)

c              ------------ print test block: begin------------
c               write(23,*) ix,iy,eddyvis(2)*1.0d4  ! cm**2/s
c               write(23,*) ix,iy,rich(2)
c               write(23,*) ix,iy,vforce(2)
c               write(81,*) ix,iy, sum(eddyvis(2:kb))/(kb-1)*1.0d4  ! cm**2/s
c               write(*,*) ix,iy, eddyvis(2:kb)
c               if ((rich(2) < 0.3).and.(kb>2)) write(81,*) ix,iy 
c               write(81,*) ix,iy, eddyvis(3)*1.0d4
c               if ((ix==40).and.(iy==40)) then
c                  do iz=1,kb
c               write(81,*)iz,eddyvis(iz)*1.0d4,rich(iz)  !,vdiffus(ix,iy,iz)
c                      write(82,*)iz,rich(iz)
c                     write(81,*)iz,rho(ix,iy,iz)
c                  enddo
c               endif
c              ------------ print test block: end------------


            endif

         enddo   ! do iy = 1, ny-1 
      enddo      ! do ix = 2, nx

c      do iz = 2, kb      
c      write(*,*)ix,iy,iz,rich(iz),eddyvis(iz),schmidt(iz),rho(ix,iy,iz)
c      enddo

      end subroutine hamsom_vertical_viscosity


      subroutine eddyvis_vertical(u, v, rho, h, depth, botlayer, dt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Port of eddyvis.f90, provided by Per Berg from DMI (per@dmi.dk) December, 2006
c     
c     Output: set module data vdiffus(ix,iy,iz) == vertical diffusivity on tracers
c         
c     vdiffus(:,:,2) is situated between layers 1,2 etc, 
c     i.e. vdiffus(:,:,iz) corresponds to z = iz-0.5
c     vdiffus(:,:,1) and vdiffus(bottom_layer+1) are set by boundary conditions
c
c     For single-layer cells (botlayer=1), vdiffus(:,:,1:2) = vdiffus_shallow
c     where vdiffus_shallow is a constant (for now)    
c     
c     Units: SI 
c     Changes (Asbjorn, Jan 2007): 
c         * modified DOS EOL->linux EOL
c         * shift from fortran free form to fixed form source
c         * for offline mode, only use coldstart update (corresponding to lll=1) in original code, 
c           i.e. no avv mixing in time
c         * removed ice coverage functionality, so setup is implicitly ices free (casus==.true. )
c           (deleted casus, cdrag, cdragrhoe, eisco, ueis, veis, rhoe, rhoi)
c         * inserted parameters:
c              real(8), parameter :: alpha0 = 10.0d00    (alpha0 removed from argument list)
c              real(8), parameter :: r      = 0.0021d00  (rdt removed from argument list)
c              real(8), parameter :: g      =    9.81d00
c              real(8), parameter :: rhow   = 1027.00d00
c              real(8), parameter :: vkarm  =    0.40d00
c         * deleted use associations (to modules dmi_omp, constants)
c         * removed "!FPsave" lines to enhance clarity
c         * externalize calculation of taus = sqrt(tauxs**2 + tauys**2), give taus(:,:) as 
c           subroutine argument. Deleted coefficient 0.5e0 in parameter hdrhovk2 (which belongs
c           to the cell averaging of taus) and renamed parameter odrhovk2 (One Divided by Rho and Von Karman^2)
c           (correspondance with per@dmi.dk Thu 1/18/2007) 
c         * removed unused variables
c         * renamed output variable avt -> vdiffus
c         * changed real(8) -> real to conform rest of code
c         * replaced tiefe(mi1) -> depth(ix,iy,1)
c           (  tiefe(mi1) er vanddybden i pkt (i,j), mens h(mi1) er lagtykkelsen af øverste lag.
c              tiefe(:) er et wetpoint-compressed 2D array
c              h(:) er et wetpoint-compressed 3D array  )
c         * restrict horizontal loop to
c              do iy = 1, ny-1
c              do ix = 2, nx
c           to avoid illegal indexation
c
c     Time step ambiguity (dt):  er vort time step, som for update af eddy visc. er et eller andet multiplum
c           af  time step for momentum equations. Jeg kan ikke svare på, hvad disse har
c           været i alle versioner gennem hele modellens operationelle levetid, men dt er
c           sikkert af størrelsen 3 min for update af vertikal eddy.
c
c     Implicitly assume w=0 ?
c     -------------------------------------------------------------------
c     Notes on DMI variables:
c        h == layer thickness  h(1) == first, h(kb) == last wet      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c - modules -------------------------------------------------------------------

c - directives ----------------------------------------------------------------
      implicit none

c - arguments -----------------------------------------------------------------
   
      real, intent(in)    :: u(:,:,:)      ! nx,ny,nz  W:E current   [meter/second]
      real, intent(in)    :: v(:,:,:) ! nx,ny,nz  S:N current   [meter/second]
      real, intent(in)    :: rho(:,:,:)    ! nx,ny,nz  water density [kg/m3]
      real, intent(in)    :: h(:,:,:)      ! nx,ny,nz  actual layer thickness (including surface elevation) [meter]
      real, intent(in)    :: depth(:,:)    ! nx,ny     actual depth of this water column (incl surface elevation) [meter]
      integer, intent(in) :: botlayer(:,:) ! nx,ny     bottom layer thickness in this water column
      real, intent(in)    :: dt            ! hydrodynamic time step [seconds]

c .... important references to module data:
c     
c        vdiffus(ix,iy,iz)  == vertical diffusivity on tracers  [meter**2/second] 
c        taus(ix,iy)        == cell centered surface stress     [N/meter**2] 
c        g (gravity constant) moved to constants.f (original value = 9.81d0)

c - local constants + vars ---------------------------------------------------------
      real, parameter         :: alpha0 = 10.0d00 ! locked as parameter, value from per@dmi.dk
      real, parameter         :: r      = 0.0021d00     ! added, asbjorn
      real, parameter         :: vkarm  = 0.40d00    ! added, asbjorn
      real, parameter         :: rhom   = rhow0
      real, parameter         :: odrhovk2 = 1.0e0/(rhom*vkarm*vkarm)
      real, parameter         :: dseven = 1.e0/7.e0

c.....added constants, Asbjorn
      real, parameter         :: vdiffus_shallow = 1.45d0 ! diffusion length per day == 500 m

      integer       :: k, kb, ki, ix, iy, nx, ny, nz
      real          :: wu, eu, sv, nv, ell0_reci, qdt
      real          :: ell0_k, fac0, fac, vau_k, ell_k, rrdvk2, avv_k
      real          :: hho, hhu, hh, rhoo, rhou, drhodz, Rg, alri
      real          :: ratio, turblm, turbgm, rdt, rr, fac12
      real          :: denom, richf_k
      real, dimension(1:size(u,3))           :: alpha, richg, rsmpr
      real, dimension(1:size(u,3))           :: redf, turblmo, turbgmo
      real, dimension(1:max(size(u,3)+1, 3)) :: duvdzq, duvdz, qmrf
      real, dimension(1:max(size(u,3)+1, 3)) :: turblmu, turbgmu
      
c -----------------------------------------------------------------------------
c  Parameters for calculation of reciproc Schmidt-Prandtl number Pr and flux 
c  Richardson number Rf from gradient Richardson number Rg using the 
c  Mellor & Yamada scheme:
c 
c    1/Pr = D/(Rg + B + sqrt(Rg*(Rg - A) + E))
c    Rf   = Rg/Pr
c 
      real, parameter :: A = 0.316e0
      real, parameter :: B = 0.186e0
      real, parameter :: C = 2.0e0*B + A
      real, parameter :: D = 0.725e0*C     ! = 0.4988
      real, parameter :: E = 0.0346e0
      real, parameter :: F = 0.5e0*D       ! = Rf,max = 0.2494
      real, parameter :: Rgcut = 100.e0*A  ! for expanding sqrt
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     - some initialisations:

      vdiffus = 0.            ! set pad value (for undefined cells)
      nx      = size(u,1)
      ny      = size(u,2)
      nz      = size(u,3)
      rdt     = r*dt                  ! rdt is bed friction parameter
      qdt     = 0.25e0*dt
      rrdvk2  = rdt/(dt*vkarm*vkarm)

c     - do the computations: (horizontal loop)  
      
      do iy = 1, ny-1
         do ix = 2, nx
            kb = botlayer(ix,iy)     ! do nothing for kb==0
            if     (kb == 1) then                     ! <-- NEW 
               vdiffus(ix,iy,1:2) = vdiffus_shallow   ! <-- NEW
            elseif (kb  > 1) then

c        -----------------------------------------------------------------
c        - Obtain (du/dz)**2, gradient Richardson number and alpha
c
c        asbjorn 19-01-2007: change in data layout from DMI setup, generally    
c                 mi1 -> ix,iy,1
c                 mi  -> ix,iy,k
c                 mib -> ix,iy,kb
               hhu      = h(ix,iy,1) ! top layer thickness
               rhou     = rho(ix,iy,1)   
               alpha(1) = 0.e0  
               do k=2,kb
                  hho    = hhu
                  hhu    = h(ix,iy,k)
                  hh     = 0.5e0*(hho+hhu)
                  rhoo   = rhou
                  rhou   = rho(ix,iy,k)
                  drhodz = (rhoo-rhou)/hh             ! positive direction up
               

c        asbjorn 19-01-2007: change in data layout from DMI setup  
c               wu = (u(m(i,j-1,k-1))-u(m(i,j-1,k)))**2
c               eu = (u(m(i,j  ,k-1))-u(mi        ))**2
c               sv = (v(m(i  ,j,k-1))-v(mi        ))**2
c               nv = (v(m(i-1,j,k-1))-v(m(i-1,j,k)))**2
               
                  wu = ( u(ix-1,iy,   k) - u(ix-1,iy,   k-1) )**2
                  eu = ( u(ix,  iy,   k) - u(ix,  iy,   k-1) )**2
                  sv = ( v(ix,  iy,   k) - v(ix,  iy,   k-1) )**2
                  nv = ( v(ix,  iy+1, k) - v(ix,  iy+1, k-1) )**2

                  duvdzq(k) = 0.5e0*((wu+eu) + (sv+nv))/(hh*hh) + 1.e-8
                  
                  richg(k) = -gdrho*drhodz/duvdzq(k)  ! positive direction up
                  alpha(k) = alpha0
               enddo

c           - account for wind stress (this part of NS ice free, removed ice parts, asbjorn):

               duvdzq(1) = odrhovk2*taus(ix,iy)/h(ix,iy,1)**2 + 1.e-8

c           - account for bedfriction:  mib = index of bottom layer
c            wu = u(m(i,j-1,kb)) 
c            eu = u(mib)          
c            sv = v(mib)            
c            nv = v(m(i-1,j,kb))    

               rr = rrdvk2/h(ix,iy,kb)**2
               wu = u(ix-1,iy,   kb)
               eu = u(ix,  iy,   kb) 
               sv = v(ix,  iy,   kb) 
               nv = v(ix,  iy+1, kb)
               
               duvdzq(kb+1) = 0.25e0*rr*((wu+eu)**2+(sv+nv)**2) + 1.e-8

c           -----------------------------------------------------------------
c           obtain |du/dz| from (du/dz)**2
               do k=1,kb+1
                  duvdz(k) = sqrt(duvdzq(k))
               enddo

c           -----------------------------------------------------------------
c           - Obtain reciproc Schmidt-Prandl number and flux Richardson No
c             using assymptotic expansion to maintain sensible FP precision.
c             We store 0.25-Rf in qmrf for use in later loop calculating L.
c             NB: remove the comments called !FPsave calculate the 
c                 the expressions safely in FP precision.
               qmrf(1)    = 0.25e0
               qmrf(kb+1) = 0.25e0
               rsmpr(1)   = 0.e0
               do k=2,kb
                  Rg = richg(k)
c                 - normal range:
                  denom    = Rg + B + sqrt(Rg*(Rg-A) + E)
                  rsmpr(k) = D/denom ! 0.002494 < 1/Pr < 1.45
                  richf_k  = Rg*rsmpr(k) ! -145 < Rf < 0.2494
                  qmrf(k)  = 0.25e0 - richf_k ! 0.0006 < qmrf < 145.25
               enddo


c           -----------------------------------------------------------------
c           - Obtain length scale:
c               (with the above expansion, 0.25-Rf >= 0.25-0.2494 = 0.0006
c                which should be sensibly handled by FP operations)
c
               fac0 = 5.e0/depth(ix,iy)     ! replaced tiefe(mi1) -> depth(ix,iy)
c
               do  k=2,kb
                  fac = h(ix,iy,k-1) + h(ix,iy,k)
                  
                  ell0_reci = abs( duvdz(k-1)*qmrf(k-1)             
     +                            -duvdz(k+1)*qmrf(k+1) )              
     +                           / ( duvdz(k)*qmrf(k)*vkarm*fac )
                  ell0_k    = 1.e0/(ell0_reci+fac0) + 1.e-2*fac
                  ratio     = 0.5e0*fac*vkarm/ell0_k
c     redf(k)   = exp(-7.e0*ratio/(7.e0+ratio))
                  redf(k)   = exp(-ratio/(1.e0+dseven*ratio))
               enddo

               turblmo(1) = 0.e0
               turbgmo(1) = 0.e0
               do k=2,kb
                  fac = vkarm*h(ix,iy,k-1)
                  turblmo(k) = redf(k)*(turblmo(k-1)    +2.e0*fac)
                  turbgmo(k) = redf(k)*(turbgmo(k-1)+duvdz(k)*fac)
               enddo
               turblmu(kb+1)=0.e0
               turbgmu(kb+1)=0.e0
               do ki=2,kb
                  k   = kb-ki+2
                  fac = vkarm*h(ix,iy,k)
                  turblmu(k) = redf(k)*(turblmu(k+1)    +2.e0*fac)
                  turbgmu(k) = redf(k)*(turbgmu(k+1)+duvdz(k)*fac)
               enddo

c           -----------------------------------------------------------------
c           - Obtain vertical eddy viscosity and diffusivity:
               do k=2,kb
                  turblm = min(turblmo(k),turblmu(k))
                  turbgm = max(turbgmo(k),turbgmu(k))
                  ell_k  = turblm + 1.e-2*( h(ix,iy,k-1) + h(ix,iy,k) )
                  vau_k  = turbgm
c..............use cold start expression (avv=0 && rlx=0)
c              original expression: 
c                   fac    = qdt*duvdz(k)
c                   avv_k  = (avv(mi) + fac*ell_k*vau_k)/(rlx+fac)
                  avv_k  = ell_k*vau_k
                  alri = alpha(k)*richg(k)
                  if (alri > 0.e0) then
                     avv_k = avv_k/sqrt(1.e0+alri)
                  endif

                  vdiffus(ix,iy,k) = rsmpr(k)*avv_k + 1.e-7   ! previous name == avt(mi)
c                 avv(mi) =          avv_k + 1.e-5

c                  if (k==2) write(23,*) ix,iy,1.0/rsmpr(k)  !avv_k*1.0d4
c                  if (k==2) write(23,*) ix,iy,avv_k*1.0d4

               enddo  ! k=2,kb
               
            endif     ! if (kb > 1) then                                      
         enddo        ! do iy=1,ny
      enddo           ! do ix=1,nx

      end subroutine eddyvis_vertical



      subroutine eddyvis_horizontal(u,v,botlayer,xfacelen,yfacelen)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Update horizontal diffusivity  hdiffus(2:nx-1, 2:ny-1, wetpts)   ! ix,iy,iz   [m**2/s]
c  
c  hdiffus is cell centered, i.e. hdiffus(ix,iy,iz) corresponds to (x,y,z) = (ix,iy,iz)
c  ---------------------------------------------------------------------------
c  Horizontalt benytter vi Smagorinsky parametriseringen i momentligningerne;
c  første side i vedhængte EddyViscosity.doc angiver det basale i denne (glem
c  side 2 i denne sammenhæng). For tracers, benytter vi så en effektiv
c  horizontal diffusivitet D_eff beregnet udfra den horizontale eddy viskositet
c  Eh som
c  
c  D_eff = 0.25*Eh/(0.1 + Eh*dt/L**2)
c  
c  hvor L er længdeskalaen i Smagorinsky modellen
c  
c  L**2 = Csmag**2 * dx * dy,    Csmag = 0.2
c  
c  Damping envelope:
c   L**2/dt ~ 2000^2/90. ~ 44444.
c   f(x) = 0.25*x/(0.1 + x/44444.)
c   Plot[{x,0.25*x/(0.1 + x/44444.)}, {x,0,500}]
c
c  Shear opeartors ranges (dudx, dudy, dvdx, dvdy) have a joint definition 
c  range of [2:nx-1, 2:ny-1, :] 
c  Entries of hdiffus outside valid range (and dry points) are assigned zero
c
c  Note: DMI damping has been dropped, use constant limit enstead, independent 
c  on tracer time step
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real, intent(in)    :: u(:,:,:)       ! ix,iy,iz  W:E current, on cell east face   [meter/second]
      real, intent(in)    :: v(:,:,:)       ! ix,iy,iz  S:N current, on cell south face  [meter/second]
      integer, intent(in) :: botlayer(:,:)  ! last wet layer (0 for dry points) ix,iy
      real, intent(in)    :: xfacelen(:)    ! xfacelen(j) is W:E length of cell j south face [meters]
                                            ! array size = xfacelen(1:ny+1) 
      real, intent(in)    :: yfacelen       ! S:N cell length [meters]
c      real, intent(in)    :: dt             ! hydrodynamic time step
c --- locals ---------
      real, parameter     :: Csmag2 = 0.2**2    ! Smagorinsky prefactor squared
      real                :: L2                 ! Smagorinsky squared length scale
      real                :: dudx, dudy, dvdx, dvdy, deform  ! v shear components, cell centered
      real                :: xf0,xf1,xfm, Eh
      integer             :: nx,ny,ix,iy,iz
c     ---------------------------------------------------------------------------
      hdiffus = 0.         ! set pad value (for undefined cells)
      nx      = size(u,1)  ! retrieve grid dimension
      ny      = size(u,2)  ! retrieve grid dimension
c   
c     --- update diffus(2:nx-1, 2:ny-1, 1:nz) 
c
      do iy=2,ny-1
         xf0 = xfacelen(iy)
         xf1 = xfacelen(iy+1)  ! ny+1 slot is defined
         xfm = 0.5*(xf0+xf1) 
         L2  = (smagorinsky_constant**2)*yfacelen*xfm

         do ix=2,nx-1
            do iz=1,botlayer(ix,iy)

               dudx = (u(ix,iy,iz) - u(ix-1,iy,iz))/xfm
               dudy = 0.25d0*(u(ix,  iy+1,iz)-u(ix,  iy-1,iz) 
     +              +         u(ix-1,iy+1,iz)-u(ix-1,iy-1,iz))/yfacelen
               dvdx = 0.25d0*((v(ix+1,iy+1,iz)-v(ix-1,iy+1,iz))/xf1
     +              +         (v(ix+1,iy  ,iz)-v(ix-1,iy  ,iz))/xf0)
               dvdy = (v(ix,iy+1,iz) - v(ix,iy,iz))/yfacelen

               deform = sqrt((dudx-dvdy)**2 + (dudy+dvdx)**2)

c..............DMI procedure: (hdiffus <- damped Eh) 
c               Eh   = L2*deform               
c               hdiffus(ix,iy,iz) = 0.25*Eh/(0.1 + Eh*dt/L2)
c
               hdiffus(ix,iy,iz) = L2*deform 
            enddo
         enddo
      enddo

c.....Impose upper limit on viscosity
      where(hdiffus > smagorinsky_eddyvisc_max) 
         hdiffus = smagorinsky_eddyvisc_max
      end where

c.....Divide by Schmidt number to obtain the diffusivity
      hdiffus = hdiffus/smagorinsky_schmidtnum


      end subroutine eddyvis_horizontal
      

      logical function range_check(x,y,z,nx,ny,nz)
c     ---------------------------------------------
c     Assess 0.5 <= x,y,z <= nx+0.5,ny+0.5,nz+0.5
c     which is the relevant interpolation ranges 
c     for cell centered data (it is not checked that
c     position is in a wet point)
c     ---------------------------------------------
      real       ::  x, y, z
      integer    :: nx,ny,nz
c     ---------------------------------------------
      range_check = (((0.5d0 < x).and.(x < (nx+0.5d0))).and.
     +               ((0.5d0 < y).and.(y < (ny+0.5d0))).and.
     +               ((0.5d0 < z).and.(z < (nz+0.5d0))))
      end function range_check



      subroutine interpolate_turbulence(position, k) 
c     -----------------------------------------------------------------------------
c     Interpolate module turbulence data 
c
c                 hdiffus(nx,ny,nz) -> k(1:2)
c                 vdiffus(nx,ny,nz+1) -> k(3)
c 
c     on position (in grid coordinates)
c     hdiffus/vdiffus are assumed up to date (no update invoked)
c     The output buffer k(:) is assumed sufficiently large (not checked)  
c     
c     The interpolation does not check for bottom crossings, but assumes diffusivities 
c     are set with appropriate padding values
c     -----------------------------------------------------------------------------
      real, intent(in)     :: position(:) ! assumed shape real(3+)
      real, intent(out)    :: k(:)        ! assumed shape real(3+)
      
      integer           :: ix0,iy0,iz0,ix1,iy1,iz1
      integer           :: nx,ny,nz,nzp1
      integer           :: nextrp, nintp
      real              :: x,y,z, sx,sy,sz
      real              :: k000,k100,k010,k110, k001,k101,k011,k111
c     -----------------------------------------------------------------------------
      nextrp = 0    ! extrapolation counter
      nintp  = 0    ! interpolation counter
      nx     = size(vdiffus, 1) ! shape(vdiffus) = (nx,ny,nz+1)
      ny     = size(vdiffus, 2)
      nzp1   = size(vdiffus, 3)
      nz     = nzp1 - 1

      x = position(1)
      y = position(2)
      z = position(3)
      
      if ( range_check(x,y,z,nx,ny,nz) ) then
         nintp  = nintp + 1
      else
         nextrp = nextrp + 1    ! do not distinguish array bound violations
      endif
c      
c.....Perform trilinear interpolation of hdiffus/vdiffus at actual position (x,y,z)
c     At grid edges, interpolation cube deflates to a square/line/point
c     Interpolation is robust against this defaltion
c 
c.....Determine horizontal diffusivity interpolation cube
c     This definition extrapolation is performed with zero slopes and continuous at edges
c     hdiffus is cell centered
c
      ix0 = max(1, min(int(x),   nx)) ! lower x face
      iy0 = max(1, min(int(y),   ny)) ! lower y face
      iz0 = max(1, min(int(z),   nz)) ! lower z face
      ix1 = max(1, min(int(x+1), nx)) ! upper x face
      iy1 = max(1, min(int(y+1), ny)) ! upper y face
      iz1 = max(1, min(int(z+1), nz)) ! upper z face
      sx = x - ix0
      sy = y - iy0
      sz = z - iz0 

c.....cube degenerate at extraplation      
      k000 = hdiffus(ix0,iy0,iz0)
      k100 = hdiffus(ix1,iy0,iz0)
      k010 = hdiffus(ix0,iy1,iz0)
      k110 = hdiffus(ix1,iy1,iz0)
      k001 = hdiffus(ix0,iy0,iz1)
      k101 = hdiffus(ix1,iy0,iz1)
      k011 = hdiffus(ix0,iy1,iz1)
      k111 = hdiffus(ix1,iy1,iz1)
      

      k(1:2) = (1-sz)*((1-sy)*(k000*(1-sx) + sx*k100)  + 
     &                  sy *(k010*(1-sx) + sx*k110)) + 
     &                  sz *((1-sy)*(k001*(1-sx) + sx*k101)  + 
     &                       sy *(k011*(1-sx) + sx*k111))
         
c
c.....Determine vertical diffusivity interpolation cube
c     This definition extrapolation is performed with zero slopes and continuous at edges
c     vdiffus is face centered at upper face, so z index + offset 
c     are different from horizontal diffusivity array, but x,y same
c
      iz0 = max(1, min(int(z+0.5d0), nzp1)) ! lower z face
      iz1 = max(1, min(int(z+1.5d0), nzp1)) ! upper z face
      sz = z - iz0 + 0.5d0
 
c.....cube degenerate at extraplation     
      k000 = vdiffus(ix0,iy0,iz0)
      k100 = vdiffus(ix1,iy0,iz0)
      k010 = vdiffus(ix0,iy1,iz0)
      k110 = vdiffus(ix1,iy1,iz0)
      k001 = vdiffus(ix0,iy0,iz1)
      k101 = vdiffus(ix1,iy0,iz1)
      k011 = vdiffus(ix0,iy1,iz1)
      k111 = vdiffus(ix1,iy1,iz1)
         
      k(3) = (1-sz)*((1-sy)*(k000*(1-sx) + sx*k100)  + 
     &                          sy *(k010*(1-sx) + sx*k110)) + 
     &                  sz *((1-sy)*(k001*(1-sx) + sx*k101)  + 
     &                          sy *(k011*(1-sx) + sx*k111))

      
c      write(*,*) "interpolate_turbulence : ", nextrp, "extrapolations"
c      write(*,*) "interpolate_turbulence : ", nintp,  "interpolations"

      end subroutine interpolate_turbulence



      subroutine interpolate_turbulence_deriv(position, dk) 
c     -----------------------------------------------------------------------------
c     Interpolation derivatives of module turbulence data with respect to 
c     scaled coordinates (NOT Cartesian coordinates)
c
c              (d/dx, d/dy) hdiffus(nx,ny,nz)   -> dk(1:2)
c              (d/dz)       vdiffus(nx,ny,nz+1) -> dk(3)
c 
c     at position  (in grid coordinates)
c     hdiffus/vdiffus are assumed up to date (no update invoked)
c     Expressions are generated by Mathematica as :
c     k = (1-sz) ((1-sy) (k000 (1-sx) + sx k100) + sy (k010 (1-sx) + sx k110)) + 
c            sz  ((1-sy) (k001 (1-sx) + sx k101) + sy (k011 (1-sx) + sx k111))
c     D[k, sx] -> dk(1)  where (k000, .. k111) comes from hdiffus
c     D[k, sy] -> dk(2)  where (k000, .. k111) comes from hdiffus
c     D[k, sz] -> dk(3)  where (k000, .. k111) comes from vdiffus
c
c     The output buffer k(:) is assumed sufficiently large (not checked) 
c     The derivative scheme is consistent with scheme used in interpolate_turbulence
c     (implementation validated by numerical differentation)
c
c     The interpolation does not check for bottom crossings, but assumes 
c     diffusivities are set with appropriate padding values
c     
c     -----------------------------------------------------------------------------
      real, intent(in)     :: position(:)  ! assumed shape real(3+)
      real, intent(out)    :: dk(:)        ! assumed shape real(3+)
      
      integer           :: ix0,iy0,iz0,ix1,iy1,iz1
      integer           :: nx,ny,nz,nzp1
      integer           :: nextrp, nintp
      real              :: x,y,z, sx,sy,sz
      real              :: k000,k100,k010,k110, k001,k101,k011,k111
c     -----------------------------------------------------------------------------     
      nextrp = 0    ! extrapolation counter
      nintp  = 0    ! interpolation counter
      nx     = size(vdiffus, 1)
      ny     = size(vdiffus, 2)
      nzp1   = size(vdiffus, 3)
      nz     = nzp1 - 1

      x = position(1)
      y = position(2)
      z = position(3)

      if ( range_check(x,y,z,nx,ny,nz) ) then
         nintp  = nintp + 1
      else
         nextrp = nextrp + 1    ! do not distinguish array bound violations
      endif
c     
c.....Perform trilinear interpolation of hdiffus/vdiffus at actual position (x,y,z)
c     At grid edges, interpolation cube deflates to a square/line/point
c     Interpolation is robust against this defaltion
c 
c.....Determine horizontal diffusivity interpolation cube
c     This definition extrapolation is performed with zero slopes and continuous at edges
c     hdiffus is cell centered
c
      ix0 = max(1, min(int(x),   nx)) ! lower x face
      iy0 = max(1, min(int(y),   ny)) ! lower y face
      iz0 = max(1, min(int(z),   nz)) ! lower z face
      ix1 = max(1, min(int(x+1), nx)) ! upper x face
      iy1 = max(1, min(int(y+1), ny)) ! upper y face
      iz1 = max(1, min(int(z+1), nz)) ! upper z face
      sx = x - ix0
      sy = y - iy0
      sz = z - iz0 
 
c.....cube degenerate at extraplation     
      k000 = hdiffus(ix0,iy0,iz0)
      k100 = hdiffus(ix1,iy0,iz0)
      k010 = hdiffus(ix0,iy1,iz0)
      k110 = hdiffus(ix1,iy1,iz0)
      k001 = hdiffus(ix0,iy0,iz1)
      k101 = hdiffus(ix1,iy0,iz1)
      k011 = hdiffus(ix0,iy1,iz1)
      k111 = hdiffus(ix1,iy1,iz1)

         
      dk(1) = ((-k000 + k100)*(1 - sy) +
     &                 (-k010 + k110)*sy)*(1 - sz) +
     &                ((-k001 + k101)*(1 - sy)     + 
     &                 (-k011 + k111)*sy)*sz

      dk(2) = (-(k000*(1 - sx)) + k010*(1 - sx) - 
     &                   k100*sx + k110*sx)*(1 - sz) +
     &                (-(k001*(1 - sx)) + k011*(1 - sx) - 
     &                   k101*sx + k111*sx)*sz

c
c.....Determine vertical diffusivity interpolation cube
c     This definition extrapolation is performed with zero slopes and continuous at edges
c     vdiffus is face centered at upper face, so z index + offset 
c     are different from horizontal diffusivity array, but x,y same
c
      iz0 = max(1, min(int(z+0.5d0), nzp1)) ! lower z face
      iz1 = max(1, min(int(z+1.5d0), nzp1)) ! upper z face
      sz = z - iz0 + 0.5d0
 
c.....cube degenerate at extraplation      
      k000 = vdiffus(ix0,iy0,iz0)
      k100 = vdiffus(ix1,iy0,iz0)
      k010 = vdiffus(ix0,iy1,iz0)
      k110 = vdiffus(ix1,iy1,iz0)
      k001 = vdiffus(ix0,iy0,iz1)
      k101 = vdiffus(ix1,iy0,iz1)
      k011 = vdiffus(ix0,iy1,iz1)
      k111 = vdiffus(ix1,iy1,iz1)
   
      
      dk(3) = -((k000*(1 - sx) + k100*sx)*(1 - sy)) +
     &                  (k001*(1 - sx) + k101*sx)*(1 - sy) - 
     &                  (k010*(1 - sx) + k110*sx)*sy +
     &                  (k011*(1 - sx) + k111*sx)*sy


c      write(*,*) "interpolate_turbulence_deriv : ",
c     +            nextrp, "extrapolations"
c      write(*,*) "interpolate_turbulence_deriv : ", 
c     +            nintp,  "interpolations"

      end subroutine interpolate_turbulence_deriv
    

      end module 
      
     
    
c$$$      program test_module
c$$$c     -----------------------------------------------
c$$$c     ifort -e90 turbulence.f
c$$$c     tests:
c$$$c       1) taus magnitude
c$$$c       2) rotate frame 90 degrees
c$$$c     -----------------------------------------------
c$$$      use turbulence
c$$$      implicit none
c$$$      integer, parameter :: nx=10, ny=10, nz=10, npos=2
c$$$      real, parameter    :: mudyn20 = 1.e-3 ! [Pa*sec] dynamic viscosity of water at 20 oC
c$$$      real, parameter    :: dd = 1.e-3
c$$$      real    :: u(nx,ny,nz)      ! W:E water current 
c$$$      real    :: v(nx,ny,nz)      ! S:N water current 
c$$$      real    :: uswind(nx,ny)    ! W:E surface wind  
c$$$      real    :: vswind(nx,ny)    ! S:N surface wind  
c$$$      real    :: rho(nx,ny,nz)    ! water density     
c$$$      real    :: h(nx,ny,nz)      ! layer thickness   
c$$$      real    :: lcd(nx,ny,nz)    ! layer center depth
c$$$      real    :: depth(nx,ny)     ! current depth at this point
c$$$      integer :: botlayer(nx,ny)  ! bottom layer
c$$$      real    :: dth              ! hydrodynamic time step
c$$$      real    :: xfacelen(ny+1)   ! WE cell face length
c$$$      real    :: yfacelen         ! SN cell face length
c$$$      real    :: pos(3,npos),k(2,npos),dk(3,npos)
c$$$      real    :: kup,kdw
c$$$c
c$$$      real    :: ddep, swfac, u0, v0
c$$$      integer :: ix,iy,iz,ib,ik,ipos
c$$$c     -----------------------------------------------
c$$$      call init_turbulence(nx,ny,nz)
c$$$    
c$$$      u        = 0.
c$$$      v        = 0.
c$$$      uswind   = 9.
c$$$      vswind   = 0.
c$$$      
c$$$      ib       = nz
c$$$      rho      = 1027.00d0
c$$$      botlayer = ib
c$$$      dth      = 90.d0
c$$$      yfacelen = 1.0d4
c$$$      xfacelen = 1.0d4
c$$$
c$$$c ------------ set h ------------
c$$$c
c$$$c      call random_number(h(1,1,:))
c$$$c      do iz = 1,ib
c$$$c         h(:,:,iz) = 3 + 5*h(1,1,iz)
c$$$c      enddo
c$$$      h = 40./nz
c$$$c      
c$$$c ------------ set depth + lcd ------------
c$$$      do ix = 1, nx
c$$$         do iy = 1, ny
c$$$c           .......................................
c$$$            depth(ix,iy) = 0. 
c$$$            do iz = 1, botlayer(ix,iy) 
c$$$               depth(ix,iy) = depth(ix,iy) + h(ix,iy,iz) 
c$$$            enddo
c$$$c           .......................................
c$$$            lcd(ix,iy,1) = 0.5* h(ix,iy,1)     
c$$$            do iz = 2, botlayer(ix,iy) 
c$$$               ddep = 0.5*( h(ix,iy,iz-1) + h(ix,iy,iz) )
c$$$               lcd(ix,iy,iz) = lcd(ix,iy,iz-1) + ddep
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$c ------------ setup rho
c$$$      do ix = 1, nx
c$$$         do iy = 1, ny            
c$$$            do iz = 1, botlayer(ix,iy) 
c$$$               ddep = lcd(ix,iy,iz)
c$$$               swfac = 1./(1+exp(-(ddep-15)/5.))
c$$$               rho(ix,iy,iz) =  1027.00d0 + swfac*15.
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$c
c$$$c ------------ setup a shear flow corresponding 
c$$$c              to the surface stress taus
c$$$c
c$$$
c$$$c      --- z shear of u ---
c$$$      u0 = 0.5
c$$$      do iz = 1,ib
c$$$         u(:,:,iz) = u0*(1. - float(iz)/(ib+1))
c$$$      enddo
c$$$      
c$$$c$$$c      --- y shear of u ---
c$$$c$$$      u0 = 0.5
c$$$c$$$      do iy = 1,ny
c$$$c$$$         u(:,iy,:) = u0*float(iy)/ny
c$$$c$$$      enddo
c$$$
c$$$c$$$c      --- x shear of v ---
c$$$c$$$      v0 = 0.5
c$$$c$$$      do ix = 1,nx
c$$$c$$$         v(ix,:,:) = v0*float(ix)/nx
c$$$c$$$      enddo
c$$$     
c$$$
c$$$      call update_turbulence(u, v, uswind, vswind, rho, h, depth,
c$$$     +                       botlayer, dth, xfacelen, yfacelen) 
c$$$       
c$$$
c$$$      do iz=1,nz
c$$$c         write(33,*) iz,vdiffus(2,2,iz)
c$$$c         write(33,*) lcd(2,2,iz), vdiffus(2,2,iz)
c$$$c         write(*,*) lcd(2,2,iz), hdiffus(4,4,iz), hdiffus(4,6,iz) 
c$$$c         write(33,*)  lcd(1,1,iz),rho(1,1,iz) 
c$$$      enddo
c$$$
c$$$c$$$c.....test derivative interpolation consistency.......................
c$$$c$$$
c$$$c$$$      call random_number(vdiffus)
c$$$c$$$      call random_number(hdiffus)
c$$$c$$$      call random_number(pos)
c$$$c$$$      pos=pos+2.0
c$$$c$$$ 
c$$$c$$$cccc.....line scan
c$$$c$$$ccc      do ix=1,1000
c$$$c$$$ccc         pos(1,1)=pos(1,1)+0.005
c$$$c$$$ccc         call interpolate_turbulence(pos, k)
c$$$c$$$ccc         write(*,*) pos(1,1), k(:,1), "xx"
c$$$c$$$ccc      enddo
c$$$c$$$ccc      stop
c$$$c$$$
c$$$c$$$      call interpolate_turbulence_deriv(pos, dk)
c$$$c$$$      do ipos=1, npos
c$$$c$$$         write(*,*) "pos=",pos(:,ipos)
c$$$c$$$c
c$$$c$$$         pos(1,ipos) = pos(1,ipos)+dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kup=k(1,ipos)
c$$$c$$$         pos(1,ipos) = pos(1,ipos)-2*dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kdw=k(1,ipos)
c$$$c$$$         pos(1,ipos) = pos(1,ipos)+dd
c$$$c$$$         write(*,*) "dk1(ana)=", dk(1,ipos),"dk1(num)=",(kup-kdw)/2/dd
c$$$c$$$c
c$$$c$$$         pos(2,ipos) = pos(2,ipos)+dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kup=k(1,ipos)
c$$$c$$$         pos(2,ipos) = pos(2,ipos)-2*dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kdw=k(1,ipos)
c$$$c$$$         pos(2,ipos) = pos(2,ipos)+dd
c$$$c$$$         write(*,*) "dk2(ana)=", dk(2,ipos),"dk2(num)=",(kup-kdw)/2/dd
c$$$c$$$c
c$$$c$$$         pos(3,ipos) = pos(3,ipos)+dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kup=k(2,ipos)
c$$$c$$$         pos(3,ipos) = pos(3,ipos)-2*dd
c$$$c$$$         call interpolate_turbulence(pos, k)
c$$$c$$$         kdw=k(2,ipos)
c$$$c$$$         pos(3,ipos) = pos(3,ipos)+dd
c$$$c$$$         write(*,*) "dk3(ana)=", dk(3,ipos),"dk3(num)=",(kup-kdw)/2/dd
c$$$c$$$c
c$$$c$$$      enddo
c$$$      
c$$$      call close_turbulence()
c$$$
c$$$      end program

c$$$      program test_module
c$$$c     -----------------------------------------------
c$$$c     ifort -e90 turbulence.f input_parser.o run_context.o constants.o
c$$$c     -----------------------------------------------
c$$$      use turbulence
c$$$      implicit none
c$$$      integer, parameter :: nx=10, ny=10, nz=10
c$$$      real, parameter    :: dd = 1.e-3
c$$$      real    :: pos(3),k(3),dk(3),z
c$$$      integer :: ix,iy,iz
c$$$c     -----------------------------------------------
c$$$   
c$$$      allocate( vdiffus(nx,ny,nz+1) )
c$$$      allocate( hdiffus(nx,ny,nz)   )
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            do iz = 1, nz+1
c$$$               vdiffus(ix,iy,iz) = 7.0d0*(iz-0.5)         ! 7z 
c$$$               if (iz<=nz) hdiffus(ix,iy,iz) = 5.0d0*iz + 3*ix + 4*iy ! 5z + 3x + 4y
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$
c$$$      pos(1:2) = 2.4
c$$$      do z = 0.5, nz+0.5, 0.01
c$$$         pos(3)=z
c$$$         call interpolate_turbulence(pos, k)
c$$$         call interpolate_turbulence_deriv(pos, dk)
c$$$         write(*,*) z,k(3),dk(3)
c$$$      enddo
c$$$ 
c$$$
c$$$      end program
