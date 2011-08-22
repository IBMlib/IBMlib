cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Provide the prey size spectrum from environmental 
c     descriptors and prey properties
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module prey_community

      use particle_state_base
      use constants
      
      implicit none            
      private     

c     -------- Prey DW m-l key: --------
c 
c      Letcher model 4 class fit: 
c
c               m_prey[myg] = 7.4427 * l_prey[mm]^2.3348 
c 
c      real, parameter :: prey_len0     = 1.0    ! [mm]
c      real, parameter :: prey_expo     = 2.3348 
c      real, parameter :: prey_mass0    = 7.4427 ! [myg] at prey_len0 
c    
c      Daewel 2008, JPR 30(1),pp1-5, eq. 2
c
c               m_prey[myg] = 21.6197 * l_prey[mm]^2.772
c
c      where 21.6197 = 1000^2.772 * 10^-7.476 / 0.32
c
      real, parameter :: prey_len0     = 1.0    ! [mm]
      real, parameter :: prey_expo     = 2.772 
      real, parameter :: prey_mass0    = 21.6197 ! [myg] at prey_len0     
c
c     Food distribution window
c
      real, parameter :: lp_SSmin = 0.02 ! [mm] lower Microplankton size limit
      real, parameter :: lp_SSmax = 2.0  ! [mm] upper Mesoplankton size limit
c      
c     --- parameters from Hufnagl & Peck (2011):
c         
      real, parameter  :: SS_slope_avg = -1.85 
      real, parameter  :: SS_slope_amp =  0.7
      real, parameter  :: SS_dayoffset = -175+365/2. ! cosine offset
      real, parameter  :: day2rad      =  2*pi/365.


c     =================== define public scope ===================

      public :: prey_mass                ! length-weight key for prey
      public :: prey_spectrum_density    ! evaluate size spectrum 


c     ===============================================================
                                  contains
c     ===============================================================
      

      subroutine prey_mass(lp, pm, dpm_dlp)
c     -------------------------------------------------------------------- 
c     prey length to prey mass function [myg]   
c     -------------------------------------------------------------------- 
      REAL,intent(in)           :: lp      ! prey length [mm]
      REAL,intent(out)          :: pm      ! prey mass   [myg] 
      REAL,intent(out),optional :: dpm_dlp ! (d/dlp) prey mass [myg/mm] 
c     --------------------------------------------------------------------
      pm = prey_mass0*(lp/prey_len0)**prey_expo
      if ( present(dpm_dlp) ) dpm_dlp = pm * prey_expo / lp
      
      end subroutine prey_mass




      subroutine prey_spectrum_density(lp, local_env, psd, dpsd_dlp)
c     ---------------------------------------------------- 
c     Lookup function for prey density psd [unit == individuals/mm3/mm] 
c     per length interval, i.e. prey_density*dL is the
c     concentration prey in a length interval [lp - 0.5*dL; lp + 0.5*dL]
c     in unit individuals/mm3
c
c     Currently the biomass is distributed on the size interval 
c     for lp : [0.02; 2.0] mm which covers
c     the commonly accepted size associations:
c
c        Microplankton: (0.02 - 0.2 mm) 
c        Mesoplankton:	(0.2  - 2.0 mm)
c
c     using a power law ansatz with a seasonal slope:
c
c         prey_density [lp < lp_SSmin]             = 0
c         prey_density [lp_SSmin < lp < lp_SSmax]  A * (lp/prey_len0)**slope(time)
c         prey_density [lp > lp_SSmax]             = 0
c
c     One principal issue with this parameterization is that
c     either lp_SSmin or lp_SSmax are very influential on
c     normalization factor A. Normally slope < 0 
c     ----------------------------------------------------   
      REAL,intent(in)                    :: lp        ! prey length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: psd       ! prey_spectrum_density
      REAL,intent(out),optional          :: dpsd_dlp  ! (d/dlp) prey_spectrum_density
c
      REAL                               :: slope,A,phase,xpo,s1,s0
      REAL                               :: bmass
c     ---------------------------------------------------- 
      if ((lp<lp_SSmin).or.(lp>lp_SSmax)) then
         psd = 0.0 ! both deriv == value/slope
         if ( present(dpsd_dlp) ) dpsd_dlp = 0.0
      else 
c
c     --- now we know lp_SSmin < lp < lp_SSmax ---
c         Determine the prey_density normalization factor A
c         local_env%zbiomass is in unit kg/m3
c         conversion factor to local unit (myg/mm3) is 1
c
         phase = day2rad*(local_env%julday - SS_dayoffset)
         slope = SS_slope_avg + SS_slope_amp*cos(phase) ! slope < 0
         xpo   = 1.0 + slope + prey_expo

         s1  = lp_SSmax/prey_len0 
         s0  = lp_SSmin/prey_len0 
         A   = local_env%zbiomass  * xpo 
     +          / prey_mass0 / prey_len0 / (s1**xpo - s0**xpo)
c  
         psd = A * (lp/prey_len0)**slope  
         if ( present(dpsd_dlp) ) dpsd_dlp = psd * slope / lp
         
      endif
      end subroutine prey_spectrum_density


      
      end  ! module
c     
