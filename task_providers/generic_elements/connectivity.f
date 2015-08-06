c=======================================================================
c     Generate a connectivity_matrix from a vector of 
c     state_attributes accepting a method
c        inquire_settling_success(state, is_settled, frombox, tobox, success)
c     This tool is not concerned with the layout of release/settlement 
c     habitats, but assumes each form a ordered set of habitats and that 
c     (frombox, tobox, success) refers to this ordering.
c     Offer diagnostic + file dump methods
c=======================================================================
      module connectivity

      use particle_state
      use netcdf

      implicit none
      private

c ----------------------------------------------------------------------
c ---------------------- module data section ---------------------------
c ----------------------------------------------------------------------

      public :: get_connectivity_matrix  ! generate connectivity matrix from sampling set
      public :: tmat_basic_diagnostics   ! basic diagnostics to stdout
      public :: tmatnc_start             ! initialize netCDF write of connectivity matrix
      public :: tmatnc_add_tmat          ! add a connectivity matrix to file
      public :: tmatnc_close             ! close connectivity matrix write

c     ----------------------------------------------
c     class tmatnc_file_handler is a local type
c     to ease I/O with a connectivity netCDF set
c     Cache various netCDF set handlers
c     ----------------------------------------------
      type tmatnc_file_handler
        integer :: ncid        ! netCDF set ID
        integer :: dimid_unlim ! netCDF unlimited dim ID
        integer :: yearid      ! netCDF var ID for time vector
        integer :: tmatid      ! netCDF var ID for conectivity matrices
      end type
      public :: tmatnc_file_handler

c ----------------------------------------------------------------------
c ---------------------- module methods section ------------------------
c ----------------------------------------------------------------------
      contains
c ----------------------------------------------------------------------


      subroutine get_connectivity_matrix(state_stack, nemit,
     +                                   tmat_prior, tmat)
c ----------------------------------------------------------------------
c     Generate connectivity matrix tmat(nsrc, ndest) from the sampling set 
c     represented by state_stack; destination habitats needs not cover all
c     possible destination, so sums over destinations may be less than one
c     even without mortality
c
c     nemit(1:nsrc): particles emitted by each source (for normalization)
c     tmat_prior(1:nsrc) is a background value applied to each tmat(isrc, :)
c     as finite sampling smoothing. tmat must have been allocaetd prior to
c     calling this subroutine
c ----------------------------------------------------------------------
      type(state_attributes), intent(in) :: state_stack(:) ! all active
      integer, intent(in)                :: nemit(:)
      real, intent(in)                   :: tmat_prior(:)
      real, intent(out)                  :: tmat(:,:)
      logical :: is_settled
      integer :: ist, frombox, tobox, isrc, nsrc, ndest
      real    :: survival, nofac
      integer :: nactive, nfree, nsett, ndead
c     ------------------------------------------------------------------
      nsrc  = size(tmat, 1)
      ndest = size(tmat, 2)
      nactive = 0
      nfree   = 0
      nsett   = 0
      ndead   = 0
c
c     ------ apply prior set to tmat ------
c
      do isrc = 1, nsrc    
         tmat(isrc,:) = tmat_prior(isrc)
      enddo
c
c     ------ add sampling set to tmat ------
c
      do ist = 1, size(state_stack)
         call inquire_settling_success(state_stack(ist), is_settled, 
     +                                 frombox, tobox, survival)
         if (is_settled) then
               tmat(frombox, tobox) = tmat(frombox, tobox) + survival
         endif
c        --- basic particle statistics ---
         nactive = nactive + 1
         if (is_settled) then
            nsett = nsett + 1
         elseif (survival == 0) then ! set to 0 for dead 
            ndead = ndead + 1
         else
            nfree = nfree + 1
         endif
      enddo   
c     
c     ------ normalize transport matrix by emissions and priors ------
c     
      do isrc = 1, nsrc    
         nofac = nemit(isrc) + ndest*tmat_prior(isrc)
         tmat(isrc,:) = tmat(isrc,:)/nofac
      enddo
c     
c     ------ basic report ------
c
      write(*,239) 
      write(*,240) nactive
      write(*,241) "free=",    nfree, 100.*nfree/nactive
      write(*,241) "settled=", nsett, 100.*nsett/nactive
      write(*,241) "dead=",    ndead, 100.*ndead/nactive
 239  format("get_connectivity_matrix: particle count summary")
 240  format("get_connectivity_matrix: active=", 1x, i12)
 241  format("get_connectivity_matrix:", a8, 1x, i12, "(", f12.7, " %)")
c
      end subroutine get_connectivity_matrix



      subroutine tmat_basic_diagnostics(tmat)
c ----------------------------------------------------------------------
c     Write basic diagnostics on connectivity matrix to stdout
c     retention is only defined if source set = destination set
c     in which case retention is tmat(i,i)
c ----------------------------------------------------------------------
      real, intent(in)       :: tmat(:,:)
     
      write(*,344) minval(sum(tmat,dim=2)), maxval( sum(tmat,dim=2))
      write(*,345) minval(sum(tmat,dim=1)), maxval( sum(tmat,dim=1))
 344  format(f12.7, " < surv < ", f12.7)
 345  format(f12.7, " < sink < ", f12.7)      

      end subroutine tmat_basic_diagnostics



cccccccccccccccc NetCDF related (prefix tmatnc_ for public methods) cccccccccccccccc



      subroutine nccheck(status)
c     -------------------------------------------------------------
c     Check status of netCDF operation and stop at error conditions
c     Apply as:
c         call nccheck( < netcdf_lib_call > )
c     Keep private
c     -------------------------------------------------------------
      integer, intent(in) :: status
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
      end subroutine nccheck  




      subroutine tmatnc_start(FILE_NAME, nfrom, nto, ncfh)
c     -------------------------------------------------------------
c     Prepare netcdf for storing year vector and a vector of tmats
c
c     The dimids array is used to pass the IDs of the dimensions of
c     the variables. Note that in fortran arrays are stored in
c     column-major format. Time dimension is assigned to the 
c     unlimited axis
c     -------------------------------------------------------------
      character*(*), intent(in)               :: FILE_NAME
      integer, intent(in)                     :: nfrom, nto
      type(tmatnc_file_handler), intent(out)  :: ncfh
      integer                                 :: dimids1(1), dimids3(3)
      integer                                 :: dimid_to, dimid_from ! keep local, only used in def mode
c     -------------------------------------------------------------      
      call nccheck( nf90_create(FILE_NAME, NF90_CLOBBER, ncfh%ncid)) 
 
c     --- create dims
      if ((nto < 1).or.(nfrom < 1)) then
         write(*,*) "tmatnc_start: dimensions < 1 for connect matrix"
         write(*,*) "ndest   = ", nto
         write(*,*) "nsource = ", nfrom
         stop 382
      endif
      call nccheck( nf90_def_dim(ncfh%ncid, "ndest",   nto,   
     +                           dimid_to) )
      call nccheck( nf90_def_dim(ncfh%ncid, "nsource", nfrom, 
     +                           dimid_from) )
      call nccheck( nf90_def_dim(ncfh%ncid, "time",    nf90_unlimited, 
     +                           ncfh%dimid_unlim) )

c     --- create vars
      dimids1 =  (/ ncfh%dimid_unlim /) 
      call nccheck(nf90_def_var(ncfh%ncid,"time", NF90_FLOAT, dimids1, 
     +                          ncfh%yearid))
      dimids3 =  (/ dimid_from, dimid_to, ncfh%dimid_unlim /) ! same order as declared in fortran
      call nccheck(nf90_def_var(ncfh%ncid,"tmat", NF90_FLOAT, dimids3, 
     +                          ncfh%tmatid))

c     --- end def mode
      call nccheck( nf90_enddef(ncfh%ncid) )

      end subroutine tmatnc_start




      subroutine tmatnc_add_tmat(ncfh, tmat, time)
c     -------------------------------------------------------------
c     Store tmat corresponding to time to open netCDF set 
c     with handler ncfh as last element along unlimited dimension
c     -------------------------------------------------------------
      type(tmatnc_file_handler), intent(in) :: ncfh
      real,                 intent(in) :: tmat(:,:)
      real,                 intent(in) :: time
      integer                          :: ioff1(1), ioff3(3), last
      character(len=NF90_MAX_NAME)     :: dimname
c     -------------------------------------------------------------
c     First probe current length of the unlimited dimension
      call nccheck( nf90_inquire_dimension(ncfh%ncid, ncfh%dimid_unlim, 
     +              dimname, last))
      ioff1(1) = last + 1
      call nccheck( nf90_put_var(ncfh%ncid, ncfh%yearid, time, ioff1) )
c
      ioff3    = 1
      ioff3(3) = last + 1
      call nccheck( nf90_put_var(ncfh%ncid, ncfh%tmatid, tmat, ioff3) )
c
      end subroutine tmatnc_add_tmat



      subroutine tmatnc_close(ncfh)
c     -------------------------------------------------------------
c     Close connectivity matrix write
c     -------------------------------------------------------------
      type(tmatnc_file_handler), intent(in) :: ncfh
      call nccheck( nf90_close(ncfh%ncid) )
      end subroutine tmatnc_close


      end module

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      Module test section
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Uncomment here below to test module connectivity
c     
c     module particle_state is a proxy test module, illustrating 
c     which parts of particle_state are used by connectivity
c     Build demo:
c        1) copy particle_state below to particle_state.f
c        2) ifort -I/usr/local/include -L/usr/local/include  -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz  -o testprogram  particle_state.f connectivity.f 
c     Run demo: 
c           testprogram; ncdump test.nc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      module particle_state
c$$$c     ---------------------------------------------------------------------------------------
c$$$c     proxy test module
c$$$c     public scope for simplicity, even though canonical interface has provate scope
c$$$c     ---------------------------------------------------------------------------------------
c$$$      implicit none
c$$$      type state_attributes
c$$$         integer :: frombox    
c$$$         integer :: tobox  
c$$$         real    :: success
c$$$         logical :: settled 
c$$$      end type
c$$$      contains
c$$$      subroutine inquire_settling_success(state, settled, 
c$$$     +                                 frombox, tobox, success)
c$$$      type(state_attributes), intent(in) :: state 
c$$$      logical, intent(out)               :: settled
c$$$      integer, intent(out)               :: frombox    
c$$$      integer, intent(out)               :: tobox   
c$$$      real,    intent(out)               :: success  
c$$$      frombox = state%frombox
c$$$      tobox   = state%tobox
c$$$      success = state%success
c$$$      settled = state%settled
c$$$      end subroutine inquire_settling_success
c$$$      end module

c$$$      program test_module
c$$$c     -------------------------------------------------------------
c$$$c     module demo
c$$$c     -------------------------------------------------------------
c$$$      use particle_state
c$$$      use connectivity
c$$$      implicit none
c$$$c     -------------------------------------------------------------
c$$$      integer, parameter :: nsrc  = 2
c$$$      integer, parameter :: ndest = 3
c$$$      integer, parameter :: npar  = 38
c$$$      integer, parameter :: npar_per_src = 14
c$$$      integer                :: i, isrc, idest
c$$$      real                   :: tmat(nsrc, ndest), tmat_prior(nsrc)
c$$$      integer                :: nemit(nsrc)
c$$$      type(state_attributes) :: pstates(npar) 
c$$$      type(tmatnc_file_handler) :: ncfh
c$$$c     -------------------------------------------------------------
c$$$      nemit = 0
c$$$      do i = 1, npar
c$$$         isrc  = min(1 + int(i / npar_per_src), nsrc)
c$$$         idest = 1 + mod(i, ndest)
c$$$c         idest = 2
c$$$         nemit(isrc) = nemit(isrc) + 1
c$$$         pstates(i)%frombox = isrc  
c$$$         pstates(i)%tobox   = idest
c$$$         pstates(i)%success = 1.0
c$$$         pstates(i)%settled = .true. 
c$$$      enddo
c$$$      do isrc = 1, nsrc
c$$$         tmat_prior(isrc) = 1.0/ndest/max(1, nemit(isrc))
c$$$      enddo
c$$$c     --- generate and write connectivity_matrix
c$$$      call get_connectivity_matrix(pstates, nemit, tmat_prior, tmat)
c$$$      call tmatnc_start("test.nc", nsrc, ndest, ncfh)   ! initialize netCDF write of connectivity matrix
c$$$      call tmatnc_add_tmat(ncfh, tmat, 2012.0)
c$$$c     --- generate and add new connectivity_matrix
c$$$      tmat_prior = 0.
c$$$      call get_connectivity_matrix(pstates, nemit, tmat_prior, tmat)
c$$$      call tmatnc_add_tmat(ncfh, tmat, 2013.0)
c$$$c     ---------------
c$$$      call tmatnc_close(ncfh)
c$$$      do isrc = 1, nsrc
c$$$         write(*,*) tmat(isrc,:) 
c$$$      enddo
c$$$     
c$$$      end program
