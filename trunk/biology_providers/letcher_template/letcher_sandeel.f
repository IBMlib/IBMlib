ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Preliminary skeleton for Letcher type larval bioenergetics
c     with daily growth increments
c
c     Phoenix version of IBMlib
c
c     asc:July2010 - made rough skeleton for Letcher model
c     zegy:Aug2010 - writes parameters and equations of growth and survival of larvae
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods 
      use physical_fields
      use run_context          ! only: ctrlfile => simulation_file
      use input_parser
      use output  ! THIS

      implicit none            
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
      real        :: length               ! length in mm
      real        :: weight               
      real        :: survival             ! 0 < survival < 1 of this larvae
      real        :: max_weight           ! previous maximum weight of this larvae
      real        :: accum_time           ! update counter of bioenergetics at a daily step

      logical     :: alive
      integer     :: hatch_day            ! Julian day
      integer     :: death_day            ! Julian day

      integer     :: orig_boxID           ! spatial release box
      integer     :: larvID               ! ensemble ID - negative is unset
      
      end type
      public :: state_attributes          ! make type state_attributes visible outside


c.....All other subroutines than these below are private to the module
      public :: init_particle_state       ! module operator
      public :: close_particle_state      ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes

      interface get_property               !THIS:1
      module procedure get_prop_state
      end interface
      public :: get_particle_version
      public :: get_property
      public :: get_metadata_state          !THIS:2
      

c     ==============================================================
c     ========             module data section              ========
c     ==============================================================
      integer, parameter    :: verbose = 0   ! control output volume

      integer               :: larv_counter  ! for generating IDs for larvae

      
c     here one should define data and parameters relating to the
c     species, not the individual, e.g. hatch length

c     0.8*ve=0.8*predsize*predvss=0.8*50*3
c     ve: swimming speed of predator in mm/second
c     predvss=3: predator's swimming speed multiplier

      real, parameter :: assimmax = 0.8      ! maximum assimilation efficiency
      real, parameter :: assimsh  = 0.002    ! assimilation efficiency shape parameter
      real, parameter :: metabnum = 4500     ! routine metabolism numerator 
      real, parameter :: metabden = 45000    ! routine metabolism denominator
      real, parameter :: metabexp = 0.97     ! routine metabolism exponent
      real, parameter :: actmetab = 2.5      ! activity metabolism multiplier
      real, parameter :: light    = 0.5417   ! proportion of day available for feeding
      real, parameter :: cmaxint  = 2.8275   ! cmax function intercept
      real, parameter :: cmaxexp  = 0.8496   ! cmax function exponent
      real, parameter :: winf     = 31940000 ! adult weight limit in microgram from Macer
      real, parameter :: linf     = 218      ! adult length limit in mm from Macer
      real, parameter :: ge       = 3.068    ! weight-length exponent from Macer
      real, parameter :: prop     = 0.5      ! proportion of reactive area used for feeding
      real, parameter :: pi       = 3.14159  ! pi value
      real, parameter :: htint    = 0.264    ! handling time intercept
      real, parameter :: htsl     = 7.0151   ! handling time slope
      real, parameter :: thresh   = 0.75     ! starvation threshold
      real, parameter :: m50sl    = 0.801    ! slope of the 50 % mortality function
      real, parameter :: ssint    = 0.776    ! swimming speed intercept
      real, parameter :: ssexp    = 1.07     ! swimming speed exponent
      real, parameter :: predden  = 1e-18    ! predator density in number/mm**3
      real, parameter :: predsize = 50       ! predator length in mm
      real, parameter :: predrdm  = 120      ! predator's reactive distance multiplier      
      real, parameter :: captexp  = 2.28     ! predator's capture success exponent
      real, parameter :: captnum  = 3.37     ! predator's capture success numerator
      real, parameter :: captden  = 4.476    ! predator's capture success denominator
      real, parameter :: propatak = 0.5      ! proportion of encounters predator attacks
      real, parameter :: ve       = 150      ! swimming speed of predator in mm/second
      integer, parameter :: nprey = 4        ! number of prey classes
      real, parameter :: prey_len(nprey) = (/0.16, 0.31, 0.71, 1.34/)     ! average prey length in mm
      real, parameter :: prey_mass(nprey) = (/0.25, 1.44, 15.8, 57.5/)    ! average prey biomass in microgram
      real, parameter :: csnum(nprey) = (/0.95, 0.9, 0.7, 0.9/)           ! capture success numerator
      real, parameter :: csden(nprey) = (/10, 750, 50000000, 500000000/)  ! capture success denominator
      real, parameter :: hsizel   = 6.3      ! average hatch length of larvae in mm
      real, parameter :: hsizew   = 605.8    ! average hatch weight in microgram from Macer model and = 195.4 from Letcher model
c     real, parameter :: predens(nprey) = 
c    +                  (/0.2265, 0.0906, 0.00604, 0.000604/)             ! prey densities from nominal IBM in number/ml

c     ================================================================
                                  contains
c     ================================================================
     
 
      subroutine init_particle_state()  ! module operator
c     ----------------------------------------------------------------
c     This subroutine should initialize this module 
c     ----------------------------------------------------------------
      write(*,*) "init_particle_state():"
      larv_counter = 1 ! ID for next larvae

c     you may read parameters like this:
c
c     call read_control_data(ctrlfile,"my_parameter_name", value)
    

      end subroutine init_particle_state



      subroutine close_particle_state() ! module operator
c     ----------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ----------------------------------------------------------------
      write(*,*) "close_particle_state():"
      larv_counter               = 0
      end subroutine close_particle_state


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------------------
c     This subroutine should initialize state
c
c     Their space part has already been initialized so that e.g. 
c     local physical conditions can be assessed using their position.
c    
c     Currently initdata is unused (in this subroutine)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ----------------------------------------------------------------
      if (time_dir < 0) then
         write(*,*) "init_state_attributes: time_dir<0 not implemented" 
         stop
      endif
        
c     parse initdata string, if needed

c     ---- initialization loop over individual particles -------------

      state%orig_boxID = emitboxID        ! store releasing box
      state%larvID     = larv_counter    
      state%accum_time = 0
      larv_counter     = larv_counter + 1 ! module data

      call set_state_hatched(state, space)


c     ----------------------------------------------------------------
      end subroutine init_state_attributes


      subroutine get_active_velocity(state, space, v_active)
c     ----------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c     ----------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     ----------------------------------------------------------------
      end subroutine get_active_velocity


   
      subroutine update_particle_state(state, space, dt)
c     ----------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt (need not be small)
c     Currently only dt == 1 day is accepted
c     ----------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space  
      real, intent(in)                        :: dt       ! in seconds
      real                                    :: xyz(3), zbiomass(1),   
     +                                           predens(nprey), 
     +                                           ingestion
      integer                                 :: status
c     ----------------------------------------------------------------
      
c     ----------------------------------------------------------------
      if (dt < 0) stop "negative dt values currently blocked"
c      if (abs(dt-86400.) > 1.e-1) then
c         write(*,*) "update_particle_state: currently only daily steps" 
c         stop 
c      endif
      write(*,*) "update_particle_state:begin"
      if (state%accum_time < 86400) then
         state%accum_time = state%accum_time + dt
         write(*,*) "update_particle_state:cached time"
         return
      else
         state%accum_time = 0
         write(*,*) "update_particle_state:1"
         call get_tracer_position(space, xyz)
         write(*,*) "update_particle_state:2"
         call interpolate_zooplankton(xyz, zbiomass, status) 
         write(*,*) "update_particle_state:3"
c        assess food availability (zbiomass -> predens)    
         call calculate_predens(zbiomass, predens)
         write(*,*) "update_particle_state:4"
         call calculate_ingestion(state, predens, ingestion)

         write(95,*) xyz, zbiomass, predens, ingestion
         write(*,*) "update_particle_state:5"
         call grow_larvae(ingestion, state)

         write(96,*) xyz, state%length, state%weight, state%max_weight 
         write(*,*) "update_particle_state:6"
         call update_survival_chance(state, space)

         write(97,*) xyz, state%survival, state%alive
         write(98,*) xyz, state%hatch_day, state%death_day

      endif
      write(*,*) "update_particle_state:end"
c     ----------------------------------------------------------------
      end subroutine update_particle_state



      subroutine delete_state_attributes(state)
c     ---------------------------------------------------------------- 
      type(state_attributes), intent(inout) :: state 
c     ----------------------------------------------------------------
      stop "delete_state_attributes not implemented"
      end subroutine delete_state_attributes



      subroutine write_state_attributes(state)
c     ----------------------------------------------------------------
c     log output of state
c     ----------------------------------------------------------------
      type(state_attributes), intent(in) :: state 
c     ----------------------------------------------------------------
      end subroutine write_state_attributes


c     ================================================================
c     ========      internal (non public) subroutine        ==========
c     ================================================================


      subroutine length_to_weight(length, weight)
c     ---------------------------------------------------
c     Convert larval length to larval weight
c     ---------------------------------------------------
c     Letcher model uses weight = lwint*(length**lwexp)
c     where lwint = 0.1674 and
c     lwexp = 3.837
c     ---------------------------------------------------
      real, intent(in)  :: length
      real, intent(out) :: weight
c     ---------------------------------------------------
      weight = winf*((length/linf)**ge)     ! Macer model
c     ---------------------------------------------------      
      end subroutine length_to_weight



      subroutine weight_to_length(weight, length_p)
c     ---------------------------------------------------
c     Convert weight to length 
c     ---------------------------------------------------
c     Letcher model uses length = (weight/lwint)**(1/lwexp)
c     where lwint = 0.1674 and
c     lwexp = 3.837
c     ---------------------------------------------------
      real, intent(in)  :: weight
      real, intent(out) :: length_p
c     ---------------------------------------------------
      length_p = linf*((weight/winf)**(1/ge)) ! Macer model
c     ---------------------------------------------------            
      end subroutine weight_to_length



      subroutine set_state_hatched(state, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      type(clock)                             :: current_time
c     ---------------------------------------------------      
      call set_tracer_mobility_free(space) 
      state%length     = hsizel          
      state%weight     = hsizew                      
      state%survival   = 1.0                  ! fresh start         
      state%max_weight = state%weight
      state%alive      = .true.
      call get_julian_day(current_time, state%hatch_day)
c     ---------------------------------------------------
      end subroutine set_state_hatched



      subroutine set_state_dead(state, space)
c     ---------------------------------------------------
c     Set larval state to dead
c     leave most variable at value just before larvae died
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      type(clock)                             :: current_time 
c     ---------------------------------------------------     
      call set_tracer_mobility_stop(space)           
      state%survival   = 0.0    
      state%alive      = .false.
      call get_julian_day(current_time, state%death_day)
c     ---------------------------------------------------
      end subroutine set_state_dead



      subroutine calculate_predens(zbiomass, predens)
c     ---------------------------------------------------
c     Calculate abundances in each prey class in number/ml
c     bulk zooplankton must be estimated down to zooplankton in sandeel diet
c     conversion from mmol-N/m**3 to total individuals is applied
c     0.19875 => microgram-dw/ml
c     0.26575 => av. weight for ind. in microgram-dw
c     ... x0.19875  / 0.26575 = ... x0.75 = ... number/ml in total
c     proportion of stages to total zooplankton in number is as in Letcher
c     0.713 for prey 1; 0.267 for prey 2; 0.0178 for prey 3;
c     0.00178 for prey 4
c     final multiplication factor to give number/ml for each prey type
c     0.53475 for prey 1; 0.20025 for prey 2; 0.01335 for prey 3;
c     0.0.001335 for prey 4
c     ---------------------------------------------------
      real, intent(in)  :: zbiomass(1)
      real, intent(out) :: predens(nprey)
      predens = 
     + (/zbiomass*0.53, zbiomass*0.2, zbiomass*0.013, zbiomass*0.0013/)
c     ---------------------------------------------------
      end subroutine calculate_predens



      subroutine calculate_ingestion(state, predens, ingestion)
c     ---------------------------------------------------
c     Letcher Eqs 3-13
c     day of first feeding ff=4.09-ffslope*length ffslope=0.237 first feeding function slope
c     compute encounter rates for each prey class
c     compute capture success for each prey class
c     compute handling time for each prey class
c     rank prey -> diet composition -> ingestion
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state     
      real, intent(in)                :: predens(nprey) 
      real, intent(out)               :: ingestion
      real                            :: alpha
      real                            :: length
      real                            :: weight
      real                            :: swim_speed 
      real                            :: react_dist
      real                            :: react_area
      real                            :: sear_vol
      real                            :: sear_vol_day
      real                            :: encount_rate(nprey)
      real                            :: encount_rate_day
      real                            :: encount_rate_day_poiss(nprey)
      real                            :: capt_success(nprey)
      real                            :: handl_time(nprey)
      real                            :: prey_rank(nprey) 
      real                            :: profit(nprey)
      real                            :: profit_sum(nprey)
      real                            :: profit_den(nprey)
      real                            :: profit_den_sum(nprey)
      real                            :: profitability(nprey)
      real                            :: encount_rate_day_poiss_f(nprey)
      real                            :: prey_bin
      real                            :: prey_bin_mass(nprey)
      real                            :: profitability_diff(nprey-1)
      real                            :: prey_bin_tot_mass
      real                            :: cons_max
      integer                         :: ip
      integer                         :: rr
      integer                         :: rank(nprey) 
      integer                         :: KFLAG = -1
      integer                         :: IER
      integer, external               :: ignbin
      integer, external               :: ignpoi
c     ---------------------------------------------------
      length = state%length
      weight = state%weight
      alpha = 0.000291*exp(9.14-(2.4*log(length))+
     +       (0.229*(log(length))**2))                                   ! in arc minutes converted to radians
      write(*,*) "update_particle_state:7" 
      swim_speed = ssint*(length**ssexp)                                 ! average swimming speed of larvae in mm/second
      write(*,*) "update_particle_state:8" 
      do ip = 1, nprey
         react_dist = prey_len(ip)/(2*(tan(alpha/2)))                    ! reactive distance in mm
      write(*,*) "update_particle_state:9" 
         react_area = (react_dist**2)*pi*prop                            ! reactive area in mm**2
      write(*,*) "update_particle_state:11" 
         sear_vol = swim_speed*react_area*0.001                          ! search volume in ml/second
      write(*,*) "update_particle_state:12" 
         sear_vol_day = sear_vol*60*60*24*light                          ! daily search volume in ml/day
      write(*,*) "update_particle_state:13"  
         encount_rate(ip) = sear_vol*predens(ip)                         ! number of prey items/second
      write(*,*) "update_particle_state:14"  
         encount_rate_day = sear_vol_day*predens(ip)                     ! number of prey items/day
      write(*,*) "update_particle_state:15"
      write(*,*) "XTRR", encount_rate_day
         encount_rate_day_poiss(ip) = ignpoi(encount_rate_day)               ! number of prey items/day from poisson distribution sampling
      write(*,*) "update_particle_state:16" 
         capt_success(ip) = (csnum(ip)*(weight**2))/(csden(ip)+
     +                      (weight**2))                                 ! capture success
      write(*,*) "update_particle_state:17" 
         handl_time(ip) = exp(htint*(10**(htsl*(prey_len(ip)/length))))  ! handling time in seconds
      write(*,*) "update_particle_state:18"
         prey_rank(ip) = (prey_mass(ip)*capt_success(ip))/handl_time(ip) ! rank of prey types 
      write(*,*) "update_particle_state:19"
      enddo
      write(*,*) encount_rate(1:nprey)
      write(*,*) encount_rate_day_poiss(1:nprey)
c     ---------------------------------------------------
c     optimal foraging
c     --------------------------------------------------- 
      call DPSORT(prey_rank, nprey, rank, KFLAG, IER) 
      do ip = 1, nprey
         profit(ip) = (prey_mass(rank(ip))*encount_rate(rank(ip)))*
     +             capt_success(rank(ip))                                ! benefit-cost ratio
      enddo   
      call cumsum(profit, profit_sum)
      do ip = 1, nprey 
         profit_den(ip) = encount_rate(rank(ip))*handl_time(rank(ip))
      enddo    
      call cumsum(profit_den, profit_den_sum)
      do ip = 1, nprey
         profitability(ip) = profit_sum(rank(ip))/
     +                      (1+profit_den_sum(rank(ip)))
         encount_rate_day_poiss_f = aint(encount_rate_day_poiss)
         prey_bin = ignbin(encount_rate_day_poiss_f(rank(ip)),
     +                     capt_success(rank(ip)))                       ! realized number of prey/day
         prey_bin_mass(ip) = prey_bin*prey_mass(rank(ip))                ! prey mass in microgram/day
      enddo
      profitability_diff = profitability(2:nprey)-
     +                     profitability(1:(nprey-1))                    ! decreasing profitability        
      call imaxloc(profitability_diff, rr)
      prey_bin_tot_mass = sum(prey_bin_mass(1:nprey), dim = 1)           ! total prey mass in microgram/day
      if ( rr /= 0 ) then
         prey_bin_tot_mass = sum(prey_bin_mass(1:rr), dim = 1)
      endif
      cons_max = cmaxint*(weight**cmaxexp)                               ! maximum consumption in microgram
      ingestion = min(prey_bin_tot_mass, cons_max)                       ! ingestion in microgram      
c     ---------------------------------------------------     
      end subroutine calculate_ingestion



      SUBROUTINE DPSORT(DX, N, IPERM, KFLAG, IER)
C***BEGIN PROLOGUE  DPSORT
C***PURPOSE  Return the permutation vector generated by sorting a given
C            array and, optionally, rearrange the elements of the array.
C            The array may be sorted in increasing or decreasing order.
C            A slightly modified quicksort algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A1B, N6A2B
C***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
C***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
C***AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DPSORT returns the permutation vector IPERM generated by sorting
C   the array DX and, optionally, rearranges the values in DX.  DX may
C   be sorted in increasing or decreasing order.  A slightly modified
C   quicksort algorithm is used.
C
C   IPERM is such that DX(IPERM(I)) is the Ith value in the
C   rearrangement of DX.  IPERM may be applied to another array by
C   calling IPPERM, SPPERM, DPPERM or HPPERM.
C
C   The main difference between DPSORT and its active sorting equivalent
C   DSORT is that the data are referenced indirectly rather than
C   directly.  Therefore, DPSORT should require approximately twice as
C   long to execute as DSORT.  However, DPSORT is more general.
C
C   Description of Parameters
C      DX - input/output -- double precision array of values to be
C           sorted.  If ABS(KFLAG) = 2, then the values in DX will be
C           rearranged on output; otherwise, they are unchanged.
C      N  - input -- number of values in array DX to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the value in the original order of the
C              DX array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting DX in increasing order and sort DX also.
C            =  1  means return the permutation vector resulting from
C                  sorting DX in increasing order and do not sort DX.
C            = -1  means return the permutation vector resulting from
C                  sorting DX in decreasing order and do not sort DX.
C            = -2  means return the permutation vector resulting from
C                  sorting DX in decreasing order and sort DX also.
C      IER - output -- error indicator:
C          =  0  if no error,
C          =  1  if N is zero or negative,
C          =  2  if KFLAG is not 2, 1, -1, or -2.
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified by John A. Wisniewski to use the Singleton
C           quicksort algorithm.
C   870423  Modified by Gregory S. Rhoads for passive sorting with the
C           option for the rearrangement of the original data.
C   890619  Double precision version of SPSORT created by D. W. Lozier.
C   890620  Algorithm for rearranging the data vector corrected by R.
C           Boisvert.
C   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
C   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
C   920507  Modified by M. McClain to revise prologue text.
C   920818  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
C***END PROLOGUE  DPSORT
C     .. Scalar Arguments ..
      INTEGER IER, KFLAG, N
C     .. Array Arguments ..
      real DX(*)
      INTEGER IPERM(*)
C     .. Local Scalars ..
      real R, TEMP
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  DPSORT
      write(*,*) DX(1:N)
      write(*,*) N
      write(*,*) KFLAG
c      stop 66
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         IER = 1
         stop "point 1"
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         IER = 2
         stop "point 1"
         RETURN
      ENDIF
C
C     Initialize permutation vector
C
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE
C
C     Return if only one value is to be sorted
C
      IF (NN .EQ. 1) RETURN
C
C     Alter array DX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 20 I=1,NN
            DX(I) = -DX(I)
   20    CONTINUE
      ENDIF
C
C     Sort DX only
C
      M = 1
      I = 1
      J = NN
      R = .375D0
C
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
   40 K = I
C
C     Select a central element of the array and save it in location L
C
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
C
C     If first element of array is greater than LM, interchange with LM
C
      IF (DX(IPERM(I)) .GT. DX(LM)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than LM, interchange with LM
C
      IF (DX(IPERM(J)) .LT. DX(LM)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
C
C        If first element of array is greater than LM, interchange
C        with LM
C
         IF (DX(IPERM(I)) .GT. DX(LM)) THEN
            IPERM(IJ) = IPERM(I)
            IPERM(I) = LM
            LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
C
C     Find an element in the second half of the array which is smaller
C     than LM
C
   60 L = L-1
      IF (DX(IPERM(L)) .GT. DX(LM)) GO TO 60
C
C     Find an element in the first half of the array which is greater
C     than LM
C
   70 K = K+1
      IF (DX(IPERM(K)) .LT. DX(LM)) GO TO 70
C
C     Interchange these elements
C
      IF (K .LE. L) GO TO 50
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 90
C
C     Begin again on another portion of the unsorted array
C
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
C
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
C
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (DX(IPERM(I)) .LE. DX(LM)) GO TO 100
      K = I
C
  110 IPERM(K+1) = IPERM(K)
      K = K-1
      IF (DX(LM) .LT. DX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
C
C     Clean up
C
  120 IF (KFLAG .LE. -1) THEN
         DO 130 I=1,NN
            DX(I) = -DX(I)
  130    CONTINUE
      ENDIF
C
C     Rearrange the values of DX if desired
C
      IF (KK .EQ. 2) THEN
C
C        Use the IPERM vector as a flag.
C        If IPERM(I) < 0, then the I-th value is in correct location
C
         DO 150 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0) THEN
               INDX = ISTRT
               INDX0 = INDX
               TEMP = DX(ISTRT)
  140          IF (IPERM(INDX) .GT. 0) THEN
                  DX(INDX) = DX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 140
               ENDIF
               DX(INDX0) = TEMP
            ENDIF
  150    CONTINUE
C
C        Revert the signs of the IPERM values
C
         DO 160 I=1,NN
            IPERM(I) = -IPERM(I)
  160    CONTINUE
C
      ENDIF
C
      RETURN
      END SUBROUTINE DPSORT


      
      subroutine cumsum(arr, cum_sum)
c     --------------------------------------------------- 
c     Cumulative sum(n) of array(n)
c     ---------------------------------------------------
      real, dimension(:), intent(in) :: arr
      real, dimension(size(arr)), intent(out) :: cum_sum
      integer :: n, j
c     --------------------------------------------------- 
      n = size(arr)
      cum_sum(1) = arr(1)
      do j = 2, n
         cum_sum(j) = cum_sum(j-1)+arr(j)
      enddo
c     --------------------------------------------------- 
      end subroutine cumsum



      subroutine imaxloc(arr, imax_loc)
c     --------------------------------------------------- 
c     Location of array maximum as an integer
c     ---------------------------------------------------
      real,  intent(in)     :: arr(:)
      integer, intent(out)  :: imax_loc
      integer, dimension(1) :: imax
c     --------------------------------------------------- 
      imax = maxloc(arr(:))
      imax_loc = imax(1)
c     --------------------------------------------------- 
      end subroutine imaxloc



      subroutine grow_larvae(ingestion, state)
c     -----------------------------------------------------------
c     Letcher Eqs 14-17
c     compute assimilation efficiency
c     growth = ingestion*assimilation efficiency - total cost
c     update state%weight 
c     update state%length
c     update state%max_weight
c     -----------------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      real, intent(in)                       :: ingestion
      real                                   :: tot_cost
      real                                   :: growth
      real                                   :: assim_eff
      real                                   :: length
      real                                   :: length_p
      real                                   :: weight     
c     -----------------------------------------------------------
      length = state%length
      weight = state%weight      
      assim_eff = assimmax*(1-(0.25*exp(-1*assimsh*(weight-10)))) ! assimilation efficiency
      call calculate_tot_cost(state, ingestion, tot_cost)         
      growth = (ingestion*assim_eff)-tot_cost                     ! growth in microgram
      state%weight = weight+growth                                ! update weight
      call weight_to_length(weight, length_p)
      state%length = max(state%length, length_p)                  ! update length
      state%max_weight = max(state%max_weight, state%weight)      ! update max_weight
c     -----------------------------------------------------------     
      end subroutine grow_larvae



      subroutine calculate_tot_cost(state, ingestion, tot_cost)
c     -----------------------------------------------------------
c     Letcher Eqs 16-17
c     -----------------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      real, intent(in)                       :: ingestion
      real                                   :: length
      real                                   :: weight
      real                                   :: rout_metab
      real                                   :: assim_eff
      real                                   :: maint
      real                                   :: beta
      real                                   :: coef
      real                                   :: sub_maint_metab
      real                                   :: tot_cost 
c     -----------------------------------------------------------     
      length = state%length
      weight = state%weight    
      rout_metab = (metabnum/metabden)*(weight**metabexp)         ! routine metabolism in microgram
c     -----------------------------------------------------------
c     Starvation metabolism check
c     -----------------------------------------------------------
      assim_eff = assimmax*(1-(0.25*exp(-1*assimsh*(weight-10)))) 
      maint = (rout_metab*(1+(actmetab*light)))/(assim_eff-0.1)
      beta = thresh**(1/(m50sl*length))
      coef = log(maint/((1-beta)*weight))
      if ( ingestion < maint ) then
         sub_maint_metab = ((1-beta)*weight)+
     +                     (maint*(1-(exp(-1*coef*ingestion/maint))))
         tot_cost = sub_maint_metab
      elseif ( ingestion >= maint ) then
         tot_cost = rout_metab+actmetab*light*rout_metab+ingestion*0.3
      endif
c     -----------------------------------------------------------      
      end subroutine calculate_tot_cost



      subroutine update_survival_chance(state, space)
c     ---------------------------------------------------
c     Update the survival counter for this larvae 
c     corresponding to life this day 
c     Letcher Eqs 18-22
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      type(spatial_attributes), intent(inout)  :: space
      real                                   :: starv_surv
      real                                   :: pred_surv
      real                                   :: tot_surv 
c     ---------------------------------------------------
      if ( state%weight < state%max_weight*0.75 ) then  ! death from starvation          
         call set_state_dead(state, space) 
      elseif ( state%weight >= state%max_weight*0.75 ) then
         starv_surv = state%survival      
         call calculate_predation_survival(state, space, pred_surv)  
         tot_surv = starv_surv*pred_surv                ! for this day
         state%survival = state%survival*tot_surv       ! for life time
         if ( state%survival == 0 ) then
            call set_state_dead(state, space)
         endif
      endif
c     ---------------------------------------------------
      end subroutine update_survival_chance
      


      subroutine calculate_predation_survival(state, space, pred_surv)
c     --------------------------------------------------- 
c     Letcher Eqs 19-22
c     predation survival probability
c     At some point we may also overlay predator maps (space available)
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      type(spatial_attributes), intent(inout)  :: space
      real, intent(out)                      :: pred_surv
      real                                   :: length
      real                                   :: r_l
      real                                   :: predrd
      real                                   :: encount_rate_w_pred
      real                                   :: swim_speed
      real                                   :: pred_prey_capt_success
      real                                   :: lambda
      real                                   :: ze
      real                                   :: vuln(5)
      real                                   :: vuln_s
      integer                                :: ip
c     ---------------------------------------------------     
      length = state%length
      swim_speed = ssint*(length**ssexp)
      r_l = (2*length)/(pi**2)
      predrd = predrdm+r_l                                               ! encounter radius of predator in  mm
      encount_rate_w_pred = pi*(predrd**2)*(((swim_speed**2)*3*(ve**2))/
     +                      3*ve)*predden                                ! encounter rate with predator in number/second
      pred_prey_capt_success = 1-((((predsize/length)+captnum)/
     +                         captden)**(-1*captexp))                   ! summary predator-prey capture success function
      lambda = encount_rate_w_pred*86400*light                           ! daily mean number of encounters
      ze = propatak*pred_prey_capt_success                               ! probability of larva being eaten once encountered
      do ip = 1, 5
         vuln = ((((lambda)**ip)*exp(-1*lambda))/fact(ip))*
     +          ((1-ze)**(ip-1))*ze
      enddo
      vuln_s = sum(vuln, dim = 1)                                        ! vulnerability of larva
      pred_surv = 1-vuln_s
c     ---------------------------------------------------      
      end subroutine calculate_predation_survival      



      function fact(n)
c     --------------------------------------------------- 
c     Factorial calculation
c     ---------------------------------------------------
      integer :: fact, i, n
c     --------------------------------------------------- 
      fact = 1
      do i = 2, n
         fact = i*fact
      enddo
c     --------------------------------------------------- 
      end function fact     


      character*100 function get_particle_version() ! FROM HERE
c------------------------------------------------------------ 
      get_particle_version = "Zeren model for sandeel larval growth" //
     +     " provider : $Rev: 181 $"
      end function



      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in) :: var
      type(polytype), intent(out) :: bucket
      integer, intent(out) :: status
c------------------------------------------------------------  
      stop "not implemented yet"
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in) :: var_name
      type(variable),intent(out) :: var
      integer, intent(out) :: status
c------------------------------------------------------------ 
      stop "not implemented yet"
      end subroutine get_metadata_state    ! TOHERE



      end module
