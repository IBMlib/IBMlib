      module eulerian_production_rate 
c     ------------------------------------------------------------------
c     Implement logistic growth rate for interface eulerian_production_rate
c     along with a horizontally heterogeneous grazing rate
c
c     Currently, only one field type is represented, but in 
c     the future, multiple (possibly overlapping) field may be
c     represented
c     ------------------------------------------------------------------
      use physical_fields 
      use mesh_grid, only: nz, ccdepth
      use horizontal_grid_transformations, only:
     +           nx,ny,get_horiz_grid_coordinates

      use input_parser           ! parameter I/O
      use run_context            ! parameter I/O
      
      implicit none
      private                    ! default scope
      
      public :: init_eulerian_production_rate   ! mandatory
      public :: close_eulerian_production_rate  ! mandatory
      public :: update_eulerian_production_rate ! mandatory
      public :: point_production_rate           ! mandatory
        

c     ----------- module data (currently only one type field) -----------

      integer              :: verbose = 1        ! module output switch

      real                 :: rate0              ! Eulerian production rate   as X/m3/sec     (static input)
      real                 :: ccap               ! Eulerian equilibrium concentartion as X/m3 (dynamic)
      real                 :: my0                ! background grazing rate as 1/sec           (dynamic)
      real,parameter       :: ccap_min = 1.0e-8  ! imposed lower limit on ccap X/m3           (static)
      real,parameter       :: my0_min  = 1.0e-10 ! imposed lower limit on my0 1/sec           (static)

      real                 :: mixfac            ! mixing factor for correction vector of (ccap and my0)
      integer              :: iter              ! running average counter
      integer              :: mix_start         ! start correction (ccap and my0) after this iter
      real                 :: intg_time_step    ! in sec - local copy, for mixing algortihms
      real                 :: mass_pa_avg       ! running average mass per area as X/m2           (in grazing area)  (static control var)
      real                 :: prod_pa_avg       ! running average production per area as X/m2/sec (in grazing area)  (static control var) 

      real                 :: consump_pa_target ! prescribed consumption as X/m2/sec         (in grazing area)  (static control var)
      real                 :: mass_pa_target    ! prescribed mass per area as X/m2           (in grazing area)  (static control var)
      real                 :: prod_pa_target    ! prescribed production per area as X/m2/sec (in grazing area)  (static control var) 

      real, allocatable    :: grazing_mask(:,:) ! (nx,ny) : vertically dynamic
      integer, allocatable :: habitat_occ(:,:)  ! (nx,ny) : vertically constant (not screened for dry overlap etc.)
      real, allocatable    :: bioactive(:)      ! (nz)    : 0 < bioactive < 1: weight of vertical biological activity


      contains

      
      subroutine init_eulerian_production_rate()
c     ------------------------------------------------
c     read production parameters 
c     ------------------------------------------------
      integer :: minmax(2)
c     ------------------------------------------------      
      call read_control_data(simulation_file, 
     +                       "eulerian_production_rate", rate0)
      call read_control_data(simulation_file, 
     +                       "eulerian_production_ccap", ccap)   
      call read_control_data(simulation_file,
     +         "eulerian_production_background_grazing",  my0)
      call read_control_data(simulation_file,
     +         "eulerian_production_mixfac",  mixfac)
      call read_control_data(simulation_file,
     +         "eulerian_production_mix_start",  mix_start)
      call read_control_data(simulation_file,
     +   "eulerian_production_consumption_target",consump_pa_target)
      call read_control_data(simulation_file,
     +   "eulerian_production_mass_target",       mass_pa_target)
      call read_control_data(simulation_file,
     +   "eulerian_production_prod_target",       prod_pa_target)
      call read_control_data(simulation_file,
     +   "eulerian_production_bioactive_minmax",  minmax)
c     ---- (re)read integration time step (needed for mixing algorithms)
      call read_control_data(simulation_file,"time_step",intg_time_step)
c
      allocate( bioactive(nz) )
      bioactive = 0.0 ! default: inactive
      bioactive( max(1,minmax(1)): min(nz,minmax(2)) ) = 1.0 ! activate slice
      iter        = 0 ! averaging counter
      mass_pa_avg = 0.0      
      prod_pa_avg = 0.0 
c
c     --- report parsed input ---
c
      write(*,181) "eulerian_production_rate",                  rate0
      write(*,181) "eulerian_production_ccap",                   ccap
      write(*,181) "eulerian_production_background_grazing",      my0
      write(*,184) "eulerian_production_mixfac",               mixfac
      write(*,183) "eulerian_production_mix_start",         mix_start
      write(*,181) "eulerian_production_consumption_target", 
     +                                              consump_pa_target
      write(*,181) "eulerian_production_mass_target",  mass_pa_target
      write(*,181) "eulerian_production_prod_target",  prod_pa_target
      write(*,182) "eulerian_production_bioactive",            minmax

 181  format("init_eulerian_production_rate",a," = ",999f10.5)
 182  format("init_eulerian_production_rate",a," = ",999f7.4)
 183  format("init_eulerian_production_rate",a," = ",999i8)
 184  format("init_eulerian_production_rate",a," = ",999e10.4)     
c
      call setup_grazing_representation()
c
      end subroutine init_eulerian_production_rate
      

      subroutine close_eulerian_production_rate()
c     ---------------------------------------------------------- 
      if ( allocated(bioactive)  ) deallocate( bioactive  )
      call close_grazing_representation()
      end subroutine close_eulerian_production_rate
      


      subroutine setup_grazing_representation()
c     ---------------------------------------------------------- 
c     Setup vertically constant grazing mask from file
c     giving in input file file_name under tag
c
c           eulerian_grazing_distribution = file_name
c
c     The file file_name contains lines with point observations
c
c          lon[deg E]  lat[deg N]  grazing[1/s]
c
c     Each point observation will represent the grid box
c     containing the point. If several points are provided for a 
c     grid box, the last point will apply (no averaging attempted)
c     Grid boxes are interpreted cell-centered, so that cell centers
c     are at integer grid coordinates and cell faces at half integers
c     Flag range violation if grid coordinates exceeds
c     (0.5,0.5) <= (x,y) <= (nx+0.5,ny+0.5)
c
c     One easy way to set grazing = 0 is to let file_name be empty
c     ---------------------------------------------------------- 
      character(len=999) :: fname
      integer            :: iunit, nread, ix,iy
      real               :: geo(2), grazing, x,y
c     ---------------------------------------------------------- 
      allocate( grazing_mask(nx,ny) ) ! vertically dynamic
      allocate( habitat_occ(nx,ny)  ) ! vertically constant
      grazing_mask = 0.0              ! initialize default
      habitat_occ  = 0                ! set default

      call read_control_data(simulation_file, 
     +                       "eulerian_grazing_distribution", fname)  
      call find_free_IO_unit(iunit) ! runtime_tools.o
      open(iunit, file=fname, status='old', action='read')
c     ---- read grazing distribution ----
      nread = 0
      do while (.true.)   ! let end clause in read() terminate loop
         read(iunit,*,end=344) geo, grazing 
         nread = nread+1 ! increment when read
         call get_horiz_grid_coordinates(geo, x, y)
         ix = int(x+0.5)  ! cell cast
         iy = int(y+0.5)  ! cell cast
         if ((ix<1).or.(ix>nx)) then
            write(*,423) nread, geo, nx,ny
            stop 
         elseif ((iy<1).or.(iy>ny)) then
            write(*,423) nread, geo, nx,ny
            stop
         else
            habitat_occ(ix,iy) = 1          ! point represents this cell
         endif
      enddo
 344  close(iunit)
      write(*,455) nread,trim(adjustl(fname))

 423  format("setup_grazing_representation: range violation (point",
     +        i6, "=", 2f12.7,") nx,ny=",2i5)
 455  format("setup_grazing_representation: read ",i6," points from ",a)
      

      end subroutine setup_grazing_representation


      subroutine close_grazing_representation()
c     ----------------------------------------------------------   
c     ----------------------------------------------------------   
      if ( allocated(grazing_mask) ) deallocate( grazing_mask )
      if ( allocated(habitat_occ)  ) deallocate( habitat_occ  )
      end subroutine close_grazing_representation





      subroutine update_eulerian_production_rate() 
c     -----------------------------------------------------------------
c     Update auxillary data needed for point by point
c     calls to point_production_rate, to avoid
c     repetitive operations. 
c     (currently void)
c     -----------------------------------------------------------------  
      end subroutine update_eulerian_production_rate       



      subroutine point_production_rate(dc_dt, c, incl,
     +                                 ixoff, iyoff, nxs, nys, nzs,
     +                                 surf_area, cell_volume)
c     -----------------------------------------------------------
c     Evaluate local Eulerian net production dc_dt in a sub box 
c     of size (nxs, nys, nzs) starting at (ixoff, iyoff) on parent grid,
c     including both growth and decay processes.
c     dc_dt, c, incl is only in interior part of the sub grid, so that index range
c     is (1:nxs, 1:nys, 1:nzs)
c     Time is assumed syncroneous with physical_fields
c     -----------------------------------------------------------
      real, intent(out)   :: dc_dt(:,:,:)  ! Eulerian production    as X/m3/sec (sub grid indexing)
      real, intent(in)    :: c    (:,:,:)  ! Eulerian concentartion as X/m3     (sub grid indexing)  
      logical, intent(in) :: incl (:,:,:)  ! exclude dry cells etc.             (sub grid indexing)  
      integer, intent(in) :: ixoff, iyoff     ! sub box offset
      integer, intent(in) :: nxs, nys, nzs    ! sub box size
      real, intent(in)    :: surf_area(:,:)     ! parent grid indexing
      real, intent(in)    :: cell_volume(:,:,:) ! parent grid indexing
c
      integer             :: i,j,k 
      real                :: my, mass, prod, sarea, cv, hvol
      real                :: mass_pa, prod_pa, vcomass, ccd
      real                :: mulfac, mass_dev, prod_dev, dccap, dmy0
      real                :: grazing_pv, prod_pv, deff, totmass
      real                :: inhab_sgz_frac, net_sgz_frac
      real                :: jac(2,2), invjac(2,2)
      real                :: det, my_star, ccap_star
      real,parameter      :: maxfac = 1.1
c     ------------------------------------------------
c
c     calculate mass (mass) and production (prod) on habitats 
c     in interior part of sub grid
c      
      mass    = 0.0    ! total mass       of food in habitat vol
      prod    = 0.0    ! total production of food in habitat vol
      sarea   = 0.0    ! total surface area of habitat area
      hvol    = 0.0    ! total volume of habitat (debug variable)
      vcomass = 0.0    ! vertical center of mass in habitat
      totmass = 0.0    ! total mass on sub grid
c
      do k=1,nzs
         do j=1,nys
            do i=1,nxs
c              ------------  habitat clause  ------------
               if ( (habitat_occ(i+ixoff, j+iyoff) > 0 ) .and.
     +               incl(i, j, k) ) then                  
                  cv      = cell_volume(i+ixoff,j+iyoff,k) ! look up on parent grid - unit = m3
                  hvol    = hvol + cv          ! unit = m3
                  mass    = mass + c(i,j,k)*cv ! unit = X
                  ccd     = ccdepth(i+ixoff,j+iyoff,k) ! unit = m
                  vcomass = vcomass + c(i,j,k)*cv*ccd  ! unit = X*m
                  prod_pv = rate0 * bioactive(k) * c(i,j,k) * ! logistic 
     +                       (1.0 - c(i,j,k)/ccap)
                  prod    = prod  + prod_pv*cv ! unit = X/sec 
                  if (k==1) sarea = sarea + surf_area(i+ixoff,j+iyoff) ! unit = m2  
               endif
c              ------------  sub grid clause  ------------
               if ( incl(i, j, k) ) then               
                  cv      = cell_volume(i+ixoff,j+iyoff,k) ! look up on parent grid - unit = m3
                  totmass = totmass + c(i,j,k)*cv ! unit = X                
               endif
c              
            enddo
         enddo
      enddo


c      
      vcomass      = vcomass/mass
      mass_pa      = mass/sarea                   ! normalize mass to X/m2
      prod_pa      = prod/sarea                   ! normalize production to X/m2/sec 
      my           = consump_pa_target*sarea/mass ! compute my>0 that implies consump_pa_target
      grazing_mask = my0 + habitat_occ*my         ! vector operation
c
c     ---- update running averages ----   
c
      iter        = iter + 1   ! initial value = 0
      mass_pa_avg = mass_pa_avg + (mass_pa - mass_pa_avg)/iter
      prod_pa_avg = prod_pa_avg + (prod_pa - prod_pa_avg)/iter
c      
c     --------- diagnostic dump section ---------
c
c      write(*,*) "deff =", hvol/sarea
c      sample point in habitat: (i,j,k) = (14, 19, 1:11)
c           
c      write(44,*) c(14, 19, 1)
c      write(44,*) mass_pa, prod_pa*86400, my*86400
c      write(44,*) mass_pa_avg, prod_pa_avg
c
c
c     --------- dynamic adjustments of dynamic parameters --------- 
c    
c     implied logistic parameters (assuming quasi equilibrium)   
c
      my_star   = prod_pa_avg/mass_pa_avg
      ccap_star = mass_pa_avg/(1 - prod_pa_avg/rate0/mass_pa_avg)
      deff      = vcomass*2  ! length scale for conversion from per-area to per-volume
c
c     compute correction vector (dccap,dmy0)
c     jac(:,:) is the response matrix to change in (ccap, my0)
c
      jac(1,1)  = 1 - my_star/rate0              ! dc_dccap
      jac(1,2)  = - ccap_star/rate0              ! dc_dmy
      jac(2,1)  = my_star*jac(1,1)               ! dp_dccap
      jac(2,2)  = mass_pa_avg + my_star*jac(1,2) ! dp_dmy
      det         = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
      invjac(1,1) =  jac(2,2)/det
      invjac(1,2) = -jac(1,2)/det
      invjac(2,1) = -jac(2,1)/det
      invjac(2,2) =  jac(1,1)/det
      mass_dev  = mass_pa_target - mass_pa_avg
      prod_dev  = prod_pa_target - prod_pa_avg
      dccap = invjac(1,1)*mass_dev + invjac(1,2)*prod_dev
      dmy0  = invjac(2,1)*mass_dev + invjac(2,2)*prod_dev
      dccap = dccap/deff         ! X/m2 -> X/m3
c
c     mix in a fraction of the correction vector (to avoid instability)
c
      if (iter > mix_start) then
c         write(*,*) mixfac*dccap/ccap,  mixfac*dmy0/my0
         ccap = ccap + mixfac*dccap            
         my0  = my0  + mixfac*dmy0
         ccap = max(ccap, ccap_min) ! impose lower limit to avoid over compensation
         my0  = max(my0,  my0_min)  ! impose lower limit to avoid over compensation
      endif

c     
      if (verbose > 0) then
c          write(44,*) iter, ccap, my0, 
c     +                abs(mass_dev/mass_pa_target), 
c     +                abs(prod_dev/prod_pa_target)
c          write(44,*) iter*intg_time_step/86400, 
c     +                abs(mass_dev/mass_pa_target), 
c     +                abs(prod_dev/prod_pa_target)
c          write(45,*) iter*intg_time_step/86400, ccap, my0
c          write(*,*) "PARADJST", ccap, my0
c          write(*,*) "PARADJST", mass/hvol, mass_pa_rdev,  prod_pa
c          write(45,*) iter*intg_time_step/86400, deff
c
         inhab_sgz_frac = my / (my0+my)
         net_sgz_frac   = my*mass / (my0*totmass + my*mass)
         write(48,*) iter*intg_time_step/86400, 
     +               inhab_sgz_frac, net_sgz_frac
c
      endif
c      
c     --- assign interior elements (1:nxs, 1:nys, 1:nzs) of sub box ---
c 
      do k=1,nzs
         do j=1,nys
            do i=1,nxs
               prod_pv      = rate0 * bioactive(k) * c(i,j,k) *  ! logistic 
     +                       (1.0 - c(i,j,k)/ccap)
               grazing_pv   = c(i,j,k)*grazing_mask(i+ixoff, j+iyoff)   
               dc_dt(i,j,k) = prod_pv - grazing_pv   
            enddo
         enddo
      enddo

      end subroutine point_production_rate

      end module 
