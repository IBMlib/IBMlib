      module larval_properties
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module implements the specific parameterization 
c     of the bioenergetics of a species
c
c     References for cod parameters:  
c         Daewel_2011:    Daewel et al, CJFAS 68, 2011, pp 426-443
c         Lough_2005:     Lough  et al, Fish. Oceanogr 14(4), 2005, pp 241-262
c         Ione_1996:      Ione et al, J . Morphology 227, 1996, pp 1-35,37-50
c         MacKenzie_1995: MacKenzie and Kiorboe, Limnol. Oceanogr 40(7), 1995, pp 1278-1289
c         Peck_2006:      Peck et al, Environ. Biol. Fishes 75(4), 2006, pp 419-429
c         Peck_2007:      Peck and Daewel,MEPS 347, 2007, pp 171-183       
c
c     Minimal testing guide:
c       test value/derivative consistency of CS, HT, SV by numerical differencing
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods 
      use physical_fields
      use particle_state_base

      implicit none
      private

c     --------------------------------------------------------------------
c     Define class physiology which are incorporated 
c     into state_attributes
c
c     NB: public scope. state_attributes must be able to access length+weight
c     --------------------------------------------------------------------

      type larval_physiology
         real    :: length             ! mandatory entry: standard length in mm 
         real    :: weight             ! mandatory entry: dry weight in myg <= nominel weight at length
         integer :: hatch_day          ! mandatory entry: Julian day
         real    :: stomach_content    ! myg

      end type

c     =================== biological parameters ===================

      real, parameter :: length0  = 5.5         ! [mm] first exogeneous feeding, ref Ione_1996
      real, parameter :: w2lxpo   = 0.247       ! ref Lough_2005
      real, parameter :: l2wxpo   = 1.0/w2lxpo  ! implied
      real, parameter :: weight0  = 68.670     ! [myg] (length0/1.935)**l2wxpo - ref Daewel_2011)
     
      real, parameter :: length_meta  = 30.0   ! [mm] end of larval phase, ref Ione_1996
  
      real, parameter :: min_condition_number = 0.4 ! ...

c
c     =================== define public scope ===================
c
      public :: init_larval_properties   ! module initialization
      public :: close_larval_properties  ! module close down

      public :: larval_physiology        ! component for state_attributes
      public :: init_larval_physiology   ! class constructor
      public :: close_larval_physiology  ! class destructor
      public :: set_larvae_hatched       ! alternative class constructor

      public :: length_to_nominal_weight ! length-weight key
      public :: weight_to_nominal_length ! corresponding weight-length key
             
      public :: capture_sucess           ! evaluate capture success of encounter 
      public :: handling_time            ! evaluate handling time 
      public :: search_volume            ! evaluate search volume            
      public :: grow_larvae
      public :: inquire_stage_change
      public :: evaluate_mortality_rate
      
c     ===============================================================
                             contains
c     ===============================================================

      subroutine init_larval_properties()
      end subroutine init_larval_properties
      
      subroutine close_larval_properties()
      end subroutine close_larval_properties


      subroutine init_larval_physiology(self,len,wgt)
c     ------------------------------------------
      type(larval_physiology), intent(out) :: self
      real, intent(in)                     :: len,wgt
c     ------------------------------------------
      self%length       = -1   ! signal unset
      self%weight       = -1   ! signal unset
      self%hatch_day    = -1   ! signal unset 
      self%stomach_content = 0.0
      self%length          = len          
      self%weight          = wgt 
      end subroutine init_larval_physiology


      subroutine close_larval_physiology(self)
      type(larval_physiology), intent(inout) :: self
      end subroutine close_larval_physiology


      subroutine set_larvae_hatched(self, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(larval_physiology), intent(inout) :: self
      type(spatial_attributes),intent(inout) :: space
      type(clock), pointer                   :: current_time
c     ---------------------------------------------------
      call init_larval_physiology(self,length0,weight0)                      
      current_time    => get_master_clock()
      call get_julian_day(current_time, self%hatch_day)

      end subroutine set_larvae_hatched



      subroutine length_to_nominal_weight(length, weight)
c     ---------------------------------------------------
c     Convert larval length to nominel larval weight
c     ---------------------------------------------------
      real, intent(in)  :: length
      real, intent(out) :: weight
      weight = weight0*(length/length0)**l2wxpo
      end subroutine length_to_nominal_weight



      subroutine weight_to_nominal_length(weight, length)
c     ---------------------------------------------------
c     Convert larval nominel weight to length 
c     ---------------------------------------------------
      real, intent(in)  :: weight
      real, intent(out) :: length
      length = length0*(weight/weight0)**w2lxpo
      end subroutine weight_to_nominal_length
      


      subroutine capture_sucess(lp, llarv, local_env, csuc, dcsuc_dlp)
c     -------------------------------------------------------------------- 
c     prey capture success and optionally its derivative wrt. lp: dcsuc_dlp
c     
c     ref: Daewel_2011
c     TODO: parameterize turbulence impact on capture_sucess
c           along with search_volume (see MacKenzie 1994, Evidence ...)  
c     -------------------------------------------------------------------- 
      REAL,intent(in)                    :: lp    ! prey length [mm]
      REAL,intent(in)                    :: llarv ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: csuc      ! capture_sucess
      REAL,intent(out),optional          :: dcsuc_dlp ! (d/dlp) capture_sucess 

      REAL                               :: lpmax ! max size an llarv can eat

      real,parameter :: CS_0         =  1.1
      real,parameter :: lpmax_inf    =  1.2936 ! [mm] 
      real,parameter :: llarv_scale  =  1.2758 ! [mm]
      real,parameter :: llarv_trans  =  2.8579 ! [mm]    
c     -------------------------------------------------------------------- 
      lpmax = lpmax_inf/(1.0 + exp(-(llarv-llarv_trans)/llarv_scale))

      if (lp>lpmax) then
         csuc = 0.0 
         if ( present(dcsuc_dlp) ) dcsuc_dlp = 0.0
         return
      endif
      
      csuc = CS_0*(1.0 - lp/lpmax)

      if ( present(dcsuc_dlp) ) dcsuc_dlp = -CS_0/lpmax

      end subroutine capture_sucess




      subroutine handling_time(lp, llarv, local_env, ht, dht_dlp)
c     -------------------------------------------------------------------- 
c     prey handling time [seconds] and optionally its derivative wrt. lp: dht_dlp
c
c     ref: Daewel_2011
c     numerical overflow protection:
c         have (d/dlp) handling_time > 0 to avoid artifacts
c     -------------------------------------------------------------------- 
      REAL,intent(in)                    :: lp    ! prey length [mm]
      REAL,intent(in)                    :: llarv ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions 
      REAL,intent(out)                   :: ht    ! handling time [seconds]
      REAL,intent(out),optional          :: dht_dlp ! (d/dlp) handling time [seconds]
c    
      REAL               :: w,dw_dlp 
      REAL,parameter     :: log10 = log(10.0)
      REAL,parameter     :: wmax = 5.0 ! numerical overflow limit

      real,parameter :: HT_fac1 = 0.264
      real,parameter :: HT_fac2 = 7.0151
c     ------------------------------------------------------------------------
      w      = HT_fac2*log10*(lp/llarv)
      dw_dlp = HT_fac2*log10/llarv
      if (w<wmax) then 
         ht = exp(HT_fac1 * exp(w))
         if ( present(dht_dlp) ) then
             dht_dlp = ht * HT_fac1 * exp(w) * dw_dlp
         endif
      else   ! avoid numerical overflow; have (d/dlp) handling_time > 0
         ht = exp(HT_fac1 * exp(wmax)) + (w-wmax)
         if ( present(dht_dlp) ) then
             dht_dlp = dw_dlp
         endif 
      endif

      end subroutine handling_time




      subroutine search_volume(lp, llarv, local_env, svol, dsvol_dlp)
c     -------------------------------------------------------------------- 
c     Search volume (unit = mm3/sec) of larvae with length llarv
c     exhibiting a pause-travel search pattern,
c     wrt. prey of length lp and optionally its derivative wrt. lp: dsvol_dlp)
c     
c     ref: Daewel_2011 and MacKenzie_1995
c
c     deriv = false: return search volume (unit = mm3/sec)   
c     deriv = true : return (d/dlp) search volume (currently not used) 
c
c     SSlarv = larval swim speed (assuming active)              [mm/sec]
c     alpha  = angle of visual acuity                           [radians]
c     SSprey = prey swim speed (before persuit?)                [mm/sec]
c     Vpause  = velocity component perpendicular to RD in pauses mode [mm/sec]
c     Vcruise = velocity component perpendicular to RD in cruise mode [mm/sec]
c     pause_freq = pause_frequency [1/sec] 
c     pause_dur  = pause_duratio    [sec] 
c
c     Note that ref Daewel_2011 contains several misprints on search_volume
c
c     TODO: * parameterize wturb
c           * FO.17:5.p333.2008 and JFB.66.p882.2005 uses search_volume_relSSprey = 3
c             this is probably based on escape velocity, which should not enter 
c             the encounter part (pre-pursuit)
c           * check prefactor 2 in V = sqrt(... + 2*wturb**2)
c             (depends of frame of reference for wturb)
c     -------------------------------------------------------------------- 
      REAL,intent(in)                    :: lp    ! prey length [mm]
      REAL,intent(in)                    :: llarv ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: svol  ! search volume [mm3/sec]
      REAL,intent(out),optional          :: dsvol_dlp ! (d/dlp) search volume 

      REAL,parameter     :: pi      = 4.0*atan(1.0)
      REAL,parameter     :: deg2rad = pi/180.0
      REAL,parameter     :: wturb = 0.0 ! current setting, later move to local_env
      REAL               :: alpha, RD, SSlarv, SSprey, Vpause
      REAL               :: dVpause_dlp
      REAL               :: dRD_dlp, dSSprey_dlp

c     encounter process parameters
      real, parameter :: SV_relSSprey = 3.0     ! prey swim speed in body lengths per sec, ref Daewel_2011
      real, parameter :: pause_freq   = 32./60. ! [1/sec] pause_frequency, ref Daewel_2011
      real, parameter :: pause_dur    = 1.4     ! [sec] pause_duration, ref Daewel_2011
c     ------------------------------------------------------------------------
c
c     --- determine reactive distance RD       
c
      alpha = 0.0167*exp(9.14 - 2.4*log(llarv) + 0.229*log(llarv)**2 ) ! degrees
      alpha = alpha*deg2rad
      RD    = lp/2.0/tan(alpha/2.0) !
c 
c     --- determine velocity component V perpendicular to RD 
c
      SSlarv = 0.261*llarv**(1.552*llarv**0.92) - 5.289/llarv ! ref Peck_2006/Daewel_2011
      SSprey = lp*SV_relSSprey         
c      Vcruise = sqrt(SSlarv**2 + SSprey**2 + 2*wturb**2)
      Vpause = sqrt(SSprey**2 + 2*wturb**2)
c
c     evaluate value/derivative of search_volume (MacKenzie_1995, Eq 1b)
c
      svol =  2.0*pi*RD**3*pause_freq/3. +
     +            pi*RD**2*Vpause*pause_freq*pause_dur

c
      if ( present(dsvol_dlp) ) then
         dRD_dlp     = 1.0/2.0/tan(alpha/2.0) ! independent of lp
         dSSprey_dlp = SV_relSSprey ! independent of lp
         dVpause_dlp = dSSprey_dlp * SSprey / Vpause
         dsvol_dlp   =  2.0*pi*RD**2*pause_freq*dRD_dlp +
     +            2.0*pi*RD*dRD_dlp *Vpause *    pause_freq*pause_dur + 
     +                pi*RD**2      *dVpause_dlp*pause_freq*pause_dur
      endif
   
      end subroutine search_volume



      subroutine grow_larvae(self, local_env, irate, dt)
c     ---------------------------------------------------
c     Update larval mass and length  corresponding to interval dt
c     subject to (maximal) ingestion rate.
c
c     ingestion/stomach evacuation dynamics is based on exact
c     ODE solution for time interval dt:
c
c      d(stomach) 
c      -----------  = [irate - evac_rate*gut_capacity] * H(0 < stomach < gut_capacity)
c          dt
c
c     where H(true/false) = 1/0. Ingestion rate is capped to evac_rate*gut_capacity, 
c     if  stomach >= gut_capacity
c     ---------------------------------------------------
      type(larval_physiology),intent(inout) :: self 
      type(local_environment),intent(in)    :: local_env ! ambient conditions
      real, intent(in)                      :: irate     ! (maximal) ingestion rate [myg/sec] 
      real, intent(in)                      :: dt        ! seconds 
c      
      real            :: growth, max_growth
      real            :: ingestion, AE
      real            :: RS, k, dl, RMcosts,w
      real            :: gut_capacity,tempfac,gut_evac_rate
      real            :: dt1,empty_space,net_rate
      real            :: b,c

      real, parameter :: hour12  = 12*3600.0
      real, parameter :: onehour = 3600.0 ! in seconds

      real, parameter :: metab_conv    = 0.00463*227.0 ! muLO2 -> myg DW conversion
      real, parameter :: act_metab_fac = 2.5  ! activity multiplier
      real, parameter :: AEinfty = 0.7   ! asymptotic assimilation efficiency
      real, parameter :: AEreduc = 0.4   ! reduction factor in AE for small larvae
      real, parameter :: AEdecay = 0.003 ! myg^-1

      real, parameter :: SDA = 0.35      ! [1]     ref: Daewel_2011 
      real, parameter :: Q10_GER = 3.0   ! [1/sec] ref: Daewel_2011 
c
c     ---------------------------------------------------      
      w = self%weight 
      
c     Compute metabolic costs RMcosts corresponding to time interval dt 
c     at ambient conditions local_env
c     FO.17:5.p333.2008, table 1, eqs 4+5
c     Q10 def: R2 = R1*Q10**((T2-T1)/10) 
c
c     --- RS = standard metabolic rate RS as myg/hour at night
 

      b = 1.029  - 0.00774*log(w)
      c = 0.1072 - 0.0032*log(w)
      RS = metab_conv * 0.00114 * w**b * exp(c*local_env%temp) ! myg/hour (ref Lough_2005)
c
c     activity (==light) correction
c      
      if (local_env%light) then
         if (self%length > 5.5) then
            k = act_metab_fac ! large larv limit
         else
            k = 1.4           ! small larv limit
         endif
      else
         k=1.0  ! night
      endif

      RMcosts = k*RS*(dt/onehour)  ! [myg] for time interval dt 

c      
c     Compute actual ingestion <= irate*dt and update stomach_content 
c     The analysis splits on whether stomach_content exceeds gut_capacity
c     during time interval dt. gut_evac_rate is linear evacuation as fraction of gut_capacity/sec
c
      gut_capacity  = 3.24 + 0.064*w   ! [myg] ref Lough_2005
      tempfac       = Q10_GER**((local_env%temp-12.0)/10.0) ! 
      gut_evac_rate = 1.792 * (self%length)**(-0.828) * tempfac  ! [1/sec] ref: Peck_2007 

      net_rate    = irate - gut_evac_rate*gut_capacity
      empty_space = gut_capacity - self%stomach_content
      if (net_rate*dt > empty_space) then ! we'll hit the roof within dt
         dt1 = empty_space/net_rate ! time where we hit the roof 
         dt1 = max(0.0, min(dt1,dt)) ! 0/0 protextion - ensure 0 < dt1 < dt
         ingestion    = irate*dt1 + gut_evac_rate*gut_capacity*(dt-dt1) ! < irate*dt 
         self%stomach_content = gut_capacity ! end with full stomach
      else
         ingestion    = irate*dt ! there is room in stomach, eat at max rate
         self%stomach_content = self%stomach_content + net_rate*dt
         self%stomach_content = max(0.0, self%stomach_content) ! never go negative
      endif

c     ---- assess assimilation efficiency AE
      AE  = AEinfty*(1.0 - AEreduc*exp(-AEdecay*(w-weight0)))

c     ---- assess standard dynamic action SDA (no variability with size/temp ??)
      
c     ---- combine pieces into a potential growth ----
      growth = ingestion*AE*(1.0 - SDA) - RMcosts   ! [myg] for time interval dt 

c     potentially impose a max growth here as well

c     ---- update larval weight for time interval dt      
      self%weight = self%weight + growth   
      call update_larval_length(self) ! sync length/weight

      end subroutine grow_larvae



      subroutine update_larval_length(self)    
c     ---------------------------------------------------
c     Syncronize length to current weight
c     Enforce that larval length can not decrease
c     This is a private method
c     ---------------------------------------------------
      type(larval_physiology),intent(inout) :: self  
      real :: lnom
      call weight_to_nominal_length(self%weight, lnom)
      if (lnom > self%length) self%length = lnom
      end subroutine update_larval_length



      subroutine inquire_stage_change(self, next) 
c     -------------------------------------------------------------------------------------------
c     Test whether to advance (or regress) ontogenetic stage 
c     Return next = .true. to advance (or regress) to ontogenetic stage (depending on sign of dt)
c     Do not consider sign of simulation time arrow
c     assumes weight/length are appropriately updated
c     -------------------------------------------------------------------------------------------
      type(larval_physiology),intent(in) :: self 
      logical,intent(out)                :: next
      real,parameter                     :: epsil = 1.e-4 ! avoid flagging if at stage boundary
c     ---------------------------------------------------   
      if ((self%length < (length0-epsil)).or.
     +    (self%length > (length_meta+epsil))) then
         next = .true.
      else
         next = .false.
      endif
      end subroutine inquire_stage_change



      subroutine evaluate_mortality_rate(self, local_env, 
     +                                   mortality_rate, die)
c     ------------------------------------------------------------------
c     Based on the current state and in the current environment
c     evaluate the current mortality_rate and wheter the larvae 
c     should die absolutely (avoid setting mortality_rate = infinity)
c     ------------------------------------------------------------------ 
      type(larval_physiology),intent(in) :: self 
      type(local_environment),intent(in) :: local_env
      real,intent(out)                   :: mortality_rate ! [1/sec]
      logical,intent(out)                :: die
      real                               :: curr_cond, nomweight 
c     ------------------------------------------------------------------
      mortality_rate = 0.0 ! elaborate later

      call length_to_nominal_weight(self%length, nomweight)
      curr_cond = self%weight / nomweight
      if (curr_cond  < min_condition_number) then
         die = .true.
      else
         die = .false.
      endif

      end subroutine evaluate_mortality_rate



      end ! module
