      module species_properties
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This module implements the specific parameterization 
c     of the bioenergetics of a species
c
c     References for sprat parameters:  
c         FO.17:5.p333.2008 = U Daewel etal, Fish Ocean 17:5.p333 (2008)
c         JFB.66.p882.2005  = MA Peck  etal, J.Fish Biol. 66,p882 (2005)
c   
c     Biology synopsis: http://www.ices.dk/marineworld/fishmap/ices/pdf/sprat.pdf 
c 
c     Russell, F.S. 1976. The eggs and planktonic stages of British marine fishes. 
c     Academic Press : metamorphosis, at 32-41mm
c     Onset of schooling: l=15 mm
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

         real    :: stomach_content      ! myg, currently just an ingestion cache
         real    :: RMcosts              ! myg
      end type

c     =================== biological parameters ===================

      real, parameter :: length0  = 6.2    ! [mm] first exogeneous feeding, FO.17:5.p333.2008
      real, parameter :: l2wxpo   = 5.022  ! ref FO.17:5.p333.2008  eq. 6
      real, parameter :: weight0  = 9.240225 ! [myg] (length0/3.982)**l2wxpo)
      real, parameter :: w2lxpo   = 1.0/l2wxpo
      
      real, parameter :: length_meta  = 32.0 ! [mm] lower size at metamorphosis 

      real,parameter :: CS_0        =  1.1
      real,parameter :: lpmax_inf   =  1.5249 ! [mm] (suspect misprint in FO.17:5.p333.2008)
      real,parameter :: lpmax_scale = 14.0    ! [mm]
      real,parameter :: lpmax_expo  = -1.63   ! should be negative

      real,parameter :: HT_fac1 = 0.264
      real,parameter :: HT_fac2 = 7.0151
               
c     avoid evaluating functions at lp=0      
      real, parameter :: lpmin = 1.e-6       ! [mm] left bracket for rank search   

c     .... metabolic parameters from  FO.17:5.p333.2008    
      real, parameter :: std_metab_fac     = 0.00272*0.00463*227.0 
      real, parameter :: std_metab_massexp = 0.80
      real, parameter :: std_metab_q10     = 2.57  ! deg celcius
      real, parameter :: std_metab_Tref    = 8.0   ! deg celcius

c     FO.17:5.p333.2008: k = 2.5 @ l = 15 mm
c                        k = 1.9 @ l = length0
      real, parameter :: act_metab_fac     = 2.5
      real, parameter :: act_metab_len     = 15.0 ! mm
      real, parameter :: act_metab_dfac    = 0.6
     
      real, parameter :: AEinfty = 0.7   ! asymptotic assimilation efficiency
      real, parameter :: AEreduc = 0.3   ! reduction factor in AE for small larvae
      real, parameter :: AEdecay = 0.003 ! myg^-1

c     FO.17:5.p333.2008: SDA(l = length0)             = 0.11
c                        SDA(l = length_meta = 32 mm) = 0.14 
      real, parameter :: SDAmin  = 0.11
      real, parameter :: SDAmeta = 0.14  
      real, parameter :: SDAslope=(SDAmeta-SDAmin)/(length_meta-length0)

c     cmax parameterized for myg/12 h
      real, parameter :: cmax_fac     = 1.315
      real, parameter :: cmax_massexp = 0.83
      real, parameter :: cmax_q10     = 2.872  ! deg celcius
      real, parameter :: cmax_Tref    = 15.0   ! deg celcius

c     encounter process parameters
      real, parameter :: SV_relSSprey = 3.0       ! prey swim speed in body lengths per sec
      real, parameter :: SV_visual_fraction = 0.5 ! fraction of encounter disc where prey are detected
      real, parameter :: min_condition_number = 0.4 ! FO.17:5.p333.2008
c
c     =================== define public scope ===================
c
      public :: larval_physiology        ! component for state_attributes
      public :: init_larval_physiology   ! constructor
      public :: close_larval_physiology  ! destructor
      public :: length_to_nominal_weight ! length-weight key
      public :: weight_to_nominal_length ! corresponding weight-length key
      public :: set_larvae_hatched       
      public :: CS                          ! evaluate capture success of encounter 
      public :: HT                          ! evaluate handling time 
      public :: SV                          ! evaluate search volume
      public :: add_routine_metabolic_costs 
      public :: add_ingestion               
      public :: grow_larvae
      public :: min_condition_number        ! hard limit for survival


c     ===============================================================
                             contains
c     ===============================================================

      subroutine init_larval_physiology(self)
      type(larval_physiology), intent(out) :: self
      self%length       = -1   ! signal unset
      self%weight       = -1   ! signal unset
      self%hatch_day    = -1   ! signal unset 
      self%stomach_content = 0.0
      self%RMcosts         = 0.0
      end subroutine init_larval_physiology


      subroutine close_larval_physiology(self)
      type(larval_physiology), intent(inout) :: self
      end subroutine close_larval_physiology



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



      subroutine set_larvae_hatched(self, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(larval_physiology), intent(inout) :: self
      type(spatial_attributes),intent(inout) :: space
      type(clock), pointer                   :: current_time
c     ---------------------------------------------------
      self%length     = length0          
      self%weight     = weight0                          
      current_time    => get_master_clock()
      call get_julian_day(current_time, self%hatch_day)

      end subroutine set_larvae_hatched
      


      real function CS(lp,llarv,deriv)
c     ---------------------------------- 
c     prey capture success  
c
c     TODO: parameterize turbulence impact on CS
c           along with SV (see MacKenzie 1994, Evidence ...)
c
c     deriv = false: return capture success 
c     deriv = true : return (d/dlp) capture success 
c     ---------------------------------- 
      REAL,intent(in)    :: lp    ! prey length [mm]
      REAL,intent(in)    :: llarv ! larval length [mm]
      logical,intent(in) :: deriv
      REAL               :: lpmax ! max size an llarv can eat
c     ---------------------------------- 
      lpmax = lpmax_inf/(1.0 + (llarv/lpmax_scale)**lpmax_expo)
      if (lp>lpmax) then
         CS = 0.0 ! both for deriv == value+slope 
      else
         if (deriv) then
            CS = -CS_0/lpmax
         else
            CS =  CS_0*(1.0 - lp/lpmax)
         endif   
      endif  ! lp>lpmax
c     write(*,*) "cs",lp,llarv,deriv,CS
      end function CS



      real function HT(lp,llarv,deriv)
c     ---------------------------------- 
c     prey handling time [seconds]
c
c     deriv = false: return handling time   
c     deriv = true : return (d/dlp) handling time   
c
c     numerical overflow protection:
c      
c     ---------------------------------- 
      REAL,intent(in)    :: lp  ! prey length [mm]
      REAL,intent(in)    :: llarv ! larval length [mm]
      logical,intent(in) :: deriv
      REAL               :: w,arg
      REAL,parameter     :: log10 = log(10.0)
      REAL,parameter     :: wmax = 5.0 ! numerical overflow limit
c     --------------------------------------
      w   = HT_fac2*log10*(lp/llarv)
      if (w<wmax) then  
         HT = exp(HT_fac1 * exp(w))
         if (deriv) then
            HT = HT * HT_fac1* exp(w) * w /lp
         endif   
      else   ! avoid numerical overflow 
         if (deriv) then
            HT = 0.0 + HT_fac2*log10/llarv
         else
            HT = exp(HT_fac1 * exp(wmax)) + (w-wmax)
         endif
      endif
c      write(*,*) "ht",lp,llarv,deriv,HT
      end function HT




      real function SV(lp,llarv,local_env,deriv)
c     ---------------------------------- 
c     Search volume (unit = mm3/sec) of larvae with length llarv
c     wrt. prey of length lp 
c     
c     deriv = false: return search volume (unit = mm3/sec)   
c     deriv = true : return (d/dlp) search volume (currently not used) 
c
c     SSlarv = larval swim speed (assuming active) [mm/sec]
c     alpha  = angle of visual acuity              [radians]
c     SSprey = prey swim speed (before persuit?)   [mm/sec]
c     V      = velocity component perpendicular to RD   [mm/sec]
c     
c     TODO: * parameterize wturb
c           * FO.17:5.p333.2008 and JFB.66.p882.2005 uses SV_relSSprey = 3
c             this is probably based on escape velocity, which should not enter 
c             the encounter part (pre-pursuit)
c           * check prefactor 2 in V = sqrt(... + 2*wturb**2)
c             (depends of frame of reference for wturb)
c           * implement deriv== .true.
c     ---------------------------------- 
      REAL,intent(in)    :: lp    ! prey length [mm]
      REAL,intent(in)    :: llarv ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      logical,intent(in) :: deriv
              
      REAL,parameter     :: pi      = 4.0*atan(1.0)
      REAL,parameter     :: deg2rad = pi/180.0
      REAL,parameter     :: wturb = 0.0 ! current setting, later move to local_env
      REAL               :: alpha, RD, SSlarv, SSprey, V
      REAL               :: dRD,dSSprey,dV
c     --------------------------------------
c
c     --- determine reactive distance RD       
c
      alpha = 0.0167*exp(9.14 - 2.4*log(llarv) + 0.229*log(llarv)**2 ) ! degrees
      alpha = alpha*deg2rad
      RD    = lp/2.0/tan(alpha/2.0) !
c 
c     --- determine velocity component V perpendicular to RD 
c
      SSlarv = 181.15/(1.0+exp(-(llarv-29.52)/5.57)) ! Munk (1992)
      SSprey = lp*SV_relSSprey
      V      = sqrt(SSlarv**2 + SSprey**2 + 2*wturb**2)
c
c     evaluate value/derivative of SV
c
      if (deriv) then

         dRD     = 1.0/2.0/tan(alpha/2.0)
         dSSprey = SV_relSSprey
         dV      = dSSprey*SSprey/V

         SV      = SV_visual_fraction * pi * (2*RD*dRD*V + RD**2*dV)

      else ! deriv == value
  
         SV      = SV_visual_fraction * pi * RD**2 * V

      endif
      end function SV



      subroutine add_routine_metabolic_costs(self,local_env,dt)
c     ---------------------------------------------------
c     Add standard metabolic costs to metabolic buffer self%RMcosts
c     corresponding to time interval dt at ambient conditions local_env
c
c     FO.17:5.p333.2008, table 1, eqs 4+5
c
c     Q10 def: R2 = R1*Q10**((T2-T1)/10) 
c     ---------------------------------------------------
      type(larval_physiology),intent(inout)   :: self
      type(local_environment),intent(in)      :: local_env ! ambient conditions
      real,intent(in)                         :: dt        ! in seconds
c     -------------------------------
      real,parameter                          :: onehour = 3600.0 ! in seconds
      real                                    :: RS, k, dl
c     ---------------------------------------------------
c
c     compute standard metabolic rate RS as myg/hour at night
c 
      RS = std_metab_fac * (self%weight)**std_metab_massexp * 
     +     std_metab_q10**((local_env%temp - std_metab_Tref)/10.0)
c
c     activity (==light) correction
c      
      if (local_env%light) then
         k = act_metab_fac
         if (self%length < act_metab_len) then
            dl = act_metab_len - self%length  ! dl > 0
            k = k - act_metab_dfac*dl/(act_metab_len-length0)
         endif
      else
         k=1.0  ! night
      endif

      self%RMcosts = self%RMcosts + k*RS*(dt/onehour)

      end subroutine add_routine_metabolic_costs



      subroutine add_ingestion(self, ingestion, dt)
c     ---------------------------------------------------
c     Add ingested mass[myg] to stomach_content corresponding
c     to time interval dt
c     Possibly apply an ingestion rate limit
c     ---------------------------------------------------
      type(larval_physiology),intent(inout)  :: self
      real,intent(in)                        :: ingestion  ! unit = myg  
      real,intent(in)                        :: dt         ! unit = sec 
c     ---------------------------------------------------
      self%stomach_content = self%stomach_content + ingestion

      end subroutine add_ingestion



      subroutine grow_larvae(self, tempavg, dt)
c     ---------------------------------------------------
c     Update larval mass (but not length) corresponding to ingestion
c     at average temperature tempavg and time interval dt
c     ---------------------------------------------------
      type(larval_physiology),intent(inout) :: self 
      real, intent(in)                      :: tempavg   ! deg celcius
      real, intent(in)                      :: dt        ! seconds
      real                                  :: growth, max_growth
      real                                  :: ingestion, AE, SDA, cmax
      real, parameter                       :: hour12 = 12*3600.0
c     ---------------------------------------------------      
      AE  = AEinfty*(1.0 - AEreduc*exp(-AEdecay*(self%weight-weight0)))
      SDA = min(1.0, SDAmin + SDAslope*(self%length - length0))
    
c     ---- impose a max consumption limit (for now using avg temp) ----
      cmax = cmax_fac * (self%weight)**cmax_massexp * 
     +                   cmax_q10**((tempavg - cmax_Tref)/10.0) ! myg/12 h
      
      ingestion = min(cmax*dt/hour12, self%stomach_content)

c     ---- combine pieces into a potential growth ----
      growth = ingestion*AE*(1.0 - SDA) - self%RMcosts

c     potentially impose a max growth here as well

c     ---- update larval weight      
      self%weight = self%weight + growth

c     ---- clear metabolic buffers
      self%stomach_content = 0.0
      self%RMcosts         = 0.0

      end subroutine grow_larvae


      end ! module
