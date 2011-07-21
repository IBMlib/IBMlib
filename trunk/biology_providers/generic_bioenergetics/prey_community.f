cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Provide the prey size spectrum from environmental 
c     descriptors and prey properties
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module prey_community

      use particle_state_base

      implicit none            
      private     

c     -------- Prey DW m-l key: --------
c 
c      Letcher model 4 class fit: 
c
c               m_prey[myg] = 7.4427 * l_prey[mm]^2.3348 
c 
c      Daewel 2008, JPR 30(1),pp1-5, eq. 2
c
c               m_prey[myg] = 21.6197 * l_prey[mm]^2.772
c
c      where 21.6197 = 1000^2.772 * 10^-7.476 / 0.32
c
      real, parameter :: prey_len0     = 1.0    ! [mm]
      real, parameter :: prey_expo     = 2.3348 
      real, parameter :: prey_mass0    = 7.4427 ! [myg] at prey_len0     
c
c     Food distribution window
c
      real, parameter :: lp_SSmin = 0.02 ! mm, lower Microplankton size limit
      real, parameter :: lp_SSmax = 2.0  ! mm, upper Mesoplankton size limit
c      
c     --- parameters from Hufnagl & Peck (2011):
c         slope = -1.85 + 0.39*SIN(2*pi*(day+94.75)/365.)
c
      real, parameter  :: SS_slope_avg = -1.85 
      real, parameter  :: SS_slope_amp =  0.39
      real, parameter  :: SS_dayoffset = 94.75

c     =================== define public scope ===================

      public :: prey_mass                 ! length-weight key for prey
      public :: spectrum_prey_density     ! evaluate size spectrum 


c     ===============================================================
                                  contains
c     ===============================================================
      

      real function prey_mass(lp,deriv)
c     ---------------------------------- 
c     prey length to prey mass function [myg] 
c
c     deriv = false: return prey_mass 
c     deriv = true : return (d/dlp) prey_mass 
c     ---------------------------------- 
      REAL,intent(in)    :: lp ! prey length [mm]
      logical,intent(in) :: deriv
c     --------------------------------------
      if (deriv) then
         prey_mass = (prey_mass0*prey_expo/prey_len0)*
     +               (lp/prey_len0)**(prey_expo-1.0) 
      else
         prey_mass = prey_mass0*(lp/prey_len0)**prey_expo
      endif
      end function prey_mass




      real function spectrum_prey_density(lp,local_env,deriv)
c     ---------------------------------------------------- 
c     Lookup function for prey density [unit == individuals/mm3/mm] 
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
c         prey_density [lp_SSmin < lp < lp_SSmax]  A * lp**slope(time)
c         prey_density [lp > lp_SSmax]             = 0
c
c     One principal issue with this parameterization is that
c     either lp_SSmin or lp_SSmax are very influential on
c     normalization factor A. Normally slope < 0
c
c     deriv = false: return prey density 
c     deriv = true : return (d/dlp) prey density  
c
c     TODO: Currently beta is constant
c     ----------------------------------------------------   
      REAL,intent(in)                    :: lp        ! prey length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      logical,intent(in)                 :: deriv
      REAL              :: slope,A,phase,xpo
      REAL,parameter    :: pi = 4.0*atan(1.0)
c     ----------------------------------------------------  
      if ((lp<lp_SSmin).or.(lp>lp_SSmax)) then
         spectrum_prey_density = 0.0 ! both deriv == value/slope
         return
      endif
c
c     --- now we know lp_SSmin < lp < lp_SSmax ---
c         Determine the prey_density normalization factor A
      phase = 2*pi*(local_env%julday + SS_dayoffset)/365.
      slope = SS_slope_avg + SS_slope_amp*SIN(phase) ! slope < 0
      xpo   = 1.0 + slope + prey_expo

      A = local_env%zbiomass * xpo * (prey_len0**prey_expo) 
     +    / prey_mass0 / (lp_SSmax**xpo - lp_SSmin**xpo)
  
      if (deriv) then
         spectrum_prey_density = A * slope * lp**(slope-1.0)
      else
         spectrum_prey_density = A * lp**slope
      endif   

      end function spectrum_prey_density

      

      
      end  ! module
c     
