ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Minimalistic Plastic module addressing 3 major traits :
c           1) Plastic fragmentation
c           2) Bouyancy
c           3) Biofouling
c
c     Dynamics:   
c             horizontal passive particle with fixed sink_speed (buoyancy = 0) OR 
c             fixed object density (buoyancy = 1) in the latter case, the Stokes velocity is applied 
c             completely passive particle obtained by sink_speed = 0 and buoyancy = 0
c
c     Settlement: 
c             1) settle_period_start < age < settle_period_end
c             2) set sticky BC along shore, and poll if ashore
c     -------------------------------------------------------------------------    
c   
c     $Rev:  $
c     $LastChangedDate:  $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state_bio ! NB, different module name

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
     
      implicit none
c      private
      
c     ---------------------------------------------
c     This type contains only specific biological attributes
c     Generic aspects, like accumulated survival, origin etc
c     are maintained by embedding type
c     ---------------------------------------------
      type state_attributes
c      private
          real    :: age,distance            ! since instantiation [days]
          real    :: nb_algae       ! number of algae on the surface

          end type
      public :: state_attributes

c     ---- export generic particle_state interface ----
      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: write_state_attributes
      public :: update_particle_state_ ! extended argument list
      public :: write_state_attributes_tofile


c     =================================================
c     ========       module data section       ========
c     =================================================
      integer, parameter    :: verbose = 0    ! control output volume      

      real       :: settle_period_start  ! since instantiation [days]
      real       :: settle_period_end    ! since instantiation [days]



c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity (default = 0.0)
c     buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density (no default) 

      integer    :: buoyancy             ! buoyancy model, default = 0
      real       :: sink_speed           ! [m/s],    applies for buoyancy=0, positive down
      real       :: rho_particle         ! [kg/m**3] fixed particle density, applies for buoyancy=1
      real       :: radius_particle      ! [m]       consider the radius of the DR for Dietrich Law
      real       :: length_particle      ! [m]       for Dietrich law
      real       :: PIV                  ! particle initial volume
      real       :: PIM                  ! Particle initial mass
      real       :: P_surf               ! Particle surface
      integer    :: file =79
      real       :: set_law                 !law used for calculating active vertical velocity. 0 = stokes, 1 = Dietrich

c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
c     ------------------------------------------------------------------------
      call read_control_data(ctrlfile,"settle_period_start",
     +                       settle_period_start)
      call read_control_data(ctrlfile,"settle_period_end",
     +                       settle_period_end)
    
      write(*,388) "settle_period_start",settle_period_start, "days"
      write(*,388) "settle_period_end"  ,settle_period_end,   "days"
c
c     Resolve buoyancy mode
c
c     fixed sinking speed takes precedence, if present (ignore rho_particle/radius_particle if present)
c
      if (count_tags(ctrlfile, "sink_speed")>0) then ! precedence, if provided

         call read_control_data(ctrlfile,"sink_speed",sink_speed) 
         buoyancy = 0
         write(*,390)  sink_speed

      elseif (count_tags(ctrlfile, "rho_particle")>0) then 

         call read_control_data(ctrlfile,"rho_particle",rho_particle) 
         call read_control_data(ctrlfile,"radius_particle",
     +                                   radius_particle)  ! must be present, if rho_particle specified, no default
         call read_control_data(ctrlfile,"length_particle",
     +                                   length_particle)
        call read_control_data(ctrlfile,"settling_law",
     +                                   set_law)
         buoyancy = 1
         PIV      = 3.14 * length_particle* radius_particle**2   !Particle initial volume
         PIM      =  PIV * rho_particle  !Particle initial mass
         P_surf        = 2*3.14*radius_particle*length_particle
         write(*,391) rho_particle, radius_particle,
     +                            length_particle, set_law

      else ! default passive particle

         sink_speed = 0.0
         buoyancy   = 0 
         write(*,392) "particles are neutrally buoyant"

      endif

 388  format(a20, " = ", f8.3, 1x, a)
 390  format("particles have constant additive sinking speed = ",
     +       f8.3, " m/s")
 391  format("particles has constant density = ", 
     +       f8.3, " kg/m**3 and radius = ", e12.3, " m")
 392  format("particles are neutrally buoyant")

      end subroutine init_particle_state


      subroutine close_particle_state()  ! module operator
      end subroutine close_particle_state



      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ---------------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ---------------------------------------------------------
      state%age           = 0.0
      state%nb_algae      = 0.0
      state%distance     =0.0
c      state%ESD     =(3 * length_particle* radius_particle**2
c     +              /4)**(0.3333)! Equivalent spherical diameter
      !state%ESD     =0! Equivalent spherical diameter
c      state%P_mass =PIM


      call set_shore_BC(space, BC_reflect)
      call set_bottom_BC(space, BC_sticky)
      end subroutine init_state_attributes




      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Implement simple vertical velocities
c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity 
c     buoyancy = 1: additive terminal Dietrich velocity, corresponding to constant object density
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)

      real, parameter                      :: patm       = 1.01325    ! [bar] applied constant atmospheric pressure
      real, parameter                      :: rhow0      = 1027.0     ! [kg/m**3], applied constant water density for compressibility effects
      real, parameter                      :: g_x_Pa2bar = 9.81/1.e5 ! g * conversion from Pa to bar
      real, parameter                      :: g          = 9.81       ! [m/s2] gravitational constant
      real, parameter                      :: mydyn      = 1.002e-3   ! [kg/m/s] dynamic viscosity of water (at 20 oC)
      real                                 :: pressure, rho, salt, temp
      real                                 :: xyz(3),wind(2),hor(2)
      integer                              :: istat, status
      real, external                       :: rho_UNESCO ! defiend in ooceanography_providers/generic_elements/water_density.f
      real                                 :: V_kin,  V_set_c,V_set_d                ! kinematic viscosity of seawater
      real                                 :: W_star =0, correction               ! dimentionmess settling velocity
      real                                 :: D_star  =0              ! dimentionmess particle diameter
      real                                 :: V_set, V_set_b,V_stokes                ! All variables to get the kinematik viscosity
      real                                 :: pi=3.1415
      real                                 :: biofouling_m, biofouling_v
      real                                 :: P_volume, P_mass, rho_tot
      real                                 :: Alg_v=2E-16  !!! [m3]
      real                                 :: rho_b=1388 !kg/m3
      real                                 :: ESD,r2(2), g_prim
      real                                 :: depth, V_ekman(2)
      real                                 :: sign_of(2),ws_khat


c     --------------------------------------------------------------
      v_active(1:2) = 0             ! no horizontal component
      if      (buoyancy == 0) then  ! fixed sinking speed,

         v_active(3) = sink_speed   ! positive down

      elseif  (buoyancy == 1) then  ! fixed particle density
c ---------- find density parapeters ----
         biofouling_v   =state%nb_algae*Alg_v*P_surf
         biofouling_m   =biofouling_v*rho_b
         P_volume       =PIV+biofouling_v
         P_mass         =PIM+biofouling_m
         rho_tot        =P_mass/P_volume
         ESD            =(6*P_volume/pi)**0.333
c         write(*,*)  rho_tot
         call get_tracer_position(space, xyz)
         call interpolate_temp(xyz,   temp,  istat)
         call interpolate_salty(xyz,  salt,  istat)
         pressure    = patm + g_x_Pa2bar*xyz(3)*rhow0   ! zeroth order approx., unit == bar
         rho         = rho_UNESCO(salt,temp,pressure)   ! unit == kg/m**3
         V_ekman=Ekman_transport(xyz,rho)
         v_active(1)=V_ekman(1)
         v_active(2)=V_ekman(2)
c         write(*,*) V_ekman
        !  write(*,*) "settling velocity", set_law
            if (set_law == 0) then
          !   write(*,*) "settling velocitystokes"
            v_active(3)=stokes(ESD,rho,rho_tot,mydyn)   !! compile stokes velocity using: ESD,rho_tot,mydyn
            elseif (set_law == 1) then
         !    write(*,*) "settling velocity Dietrich"
            v_active(3)=Dietrich(ESD,rho, rho_tot,mydyn)!Dietrich(ESD,rho,rho_tot,mydyn)

            elseif (set_law == 2) then
                V_kin       =  mydyn/rho
c            correction=2.*radius_particle*1000.*length_particle*1000./
c     +                 (52.749*length_particle*1000.+12.691)
c                v_active_mm=pi*g*(rho_tot-rho)*correction/
c     +                      (1000.*rho*2.*V_kin)
             g_prim=g*(rho_tot-rho)/(rho*1000.0);
             correction=2*radius_particle*1000.0*length_particle*1000.0/
     +       ((52.749*length_particle*1000.0)+11.442)
             ws_khat=pi*g_prim*correction/(2*V_kin)


c                g_prim=g*(rho_tot-rho)/(rho*1000.0)
c                correction=2.0*3.e-4*1000.0*0.1*1000.0/
c     +              ((52.749*0.1*1000.0)+11.442)
c                ws_khat=pi*g_prim*correction/(2.0*V_kin)
             v_active(3)=ws_khat/1000
c            write(*,*) ws_khat
            else
            write(*,*) "Settling law, not entered in Simpar"
            endif

      else
         write(*,*) "get_active_velocity: unsupported buoyancy=",
     +              buoyancy
         stop 653
      endif
      !   write(*,*) v_active(3)
c      call interpolate_wind_stress(geo, r2, status)
      end subroutine get_active_velocity

      function stokes(ESD,rho,rho_tot,mydyn)
c     --------------------------------------------------------------
c     Implement simple stokes
c     --------------------------------------------------------------
       real:: stokes
       real:: ESD
       real:: rho_tot
       real:: rho
       real:: g=9.81
       real:: mydyn
      stokes = 2.*(rho_tot-rho)*g*
     +                 (ESD/2.)**2./9./mydyn  ! [m/s], Stokes terminal velocity law  (positive down) //// we dont want stokes
        return
      end function stokes


      function Dietrich(ESD,rho,rho_tot,mydyn)
c     --------------------------------------------------------------
c     Implement simple Dietrich
c     --------------------------------------------------------------

      real:: Dietrich
      real:: ESD
      real:: rho_tot
      real:: rho
      real:: g=9.81
      real:: mydyn, diff
      real:: V_kin, W_star,  D_star, V_set
      diff=rho_tot-rho
      V_kin       =  mydyn/rho
      D_star      = (rho_tot-rho)*g*
     +                (ESD**3.0)/(rho*(V_kin**2.0))
       if (D_star  <  0.05) then ! use Stokes for Da <= 0.05, from Dietrich, 1982 paper
        W_star=D_star**2.0 /5832.
c      V_set=((rho_tot-rho)*9.81*W_star*V_kin/rho)**0.3333
        V_set=sign(1.0,diff)*(abs(rho_tot-rho)*9.81*W_star*
     +        V_kin/rho)**0.3333

c      V_set=1.74E-4**0.3333 *ESD**2.0 *(rho_tot-rho)
c     +                                      *g/(rho*V_kin)  ! V_set works but the other give a value of 1 (strangely...)
c              write(*,*) "small d star", D_star, state%ESD
      elseif (D_star<= 5.0E9) then
        W_star=exp(-3.76715 + 1.92944 * log(D_star)
     +                -0.09815*log(D_star)**2.0 - 0.00575*log(D_star)**3.0+
     +                 0.00056*log(D_star)**4.0)

c       V_set=((rho_tot-rho)*9.81*W_star*V_kin/rho)**0.3333
        V_set=sign(1.0,diff)*(abs(rho_tot-rho)*9.81*W_star*
     +        V_kin/rho)**0.3333

c             write(*,*) "execute dstar big too",W_star, D_star
      else
        V_set=0.0
      end if
      Dietrich=V_set
      return
      end function Dietrich

      function Ekman_transport(xyz,rho)
c     --------------------------------------------------------------
c     Compile Ekmann transport from  wind stress
c     --------------------------------------------------------------
        real      :: xyz(3), hor(2), V_ekman(2),rho
        real      ::  wind(2), depth,Ekman_transport(2)
        real      :: Eddy_v  ![m2/day]
        real      :: f_cor ![m2/day]
        real      :: coef !decrease with depth
        real      :: pi=3.14

        integer   :: status
        hor(1)=xyz(1)
        hor(2)=xyz(2)
        f_cor=2.* 7.292e-5
        Eddy_v=20./(24.*60.*60.)
        coef=sqrt(f_cor/(2.*Eddy_v))
        call interpolate_wind_stress(hor, wind, status)
        depth = xyz(3)
        V_ekman(1)=sign(1.,wind(1))*sqrt(abs(wind(1))/rho)
     +              *exp(-coef*depth)*cos(pi/4+coef*depth)
        V_ekman(2)=sign(1.,wind(2))*sqrt(abs(wind(2))/rho)
     +              *exp(-coef*depth)*sin((pi/4)+(coef*depth))
        Ekman_transport=V_ekman
c        write(*,*) wind, coef, Eddy_v, f_cor, Ekman_transport
c          wind_10(1)=sign(1.,wind(1))*sqrt(abs(wind(1))/0.0158)
c          wind_10(2)=sign(1.,wind(2))*sqrt(abs(wind(2))/0.0158)
      return
      end function Ekman_transport

      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
c     --------------------------------------------------------------
      write(78,*)  state%age,  state%nb_algae,
     +  state%distance
      end subroutine write_state_attributes

      subroutine write_state_attributes_tofile(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
c      integer                              :: file
c     --------------------------------------------------------------

      write(79,*)  state%age,  state%nb_algae,
     +   state%distance
      end subroutine write_state_attributes_tofile


      subroutine update_particle_state_(state, space, dt, 
     +                       mortality_rate, die, next)
c     --------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt (extended argument list)
c     Also set auxillary parameters:
c        mortality_rate: current net mortality rate [1/sec]
c        die           : if .true. kill this entity (avoid having mortality_rate=infinity)
c        next          : if .true. advance (or regress) to ontogenetic stage (depending on sign of dt)
c
c     Settlement conditions: 
c          1) settle_period_start < age < settle_period_end
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real                                    :: V_settling
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next 
      real                                    :: xyz(3)
      integer                                 :: stat,istat2
      integer                                 :: statud
      real                                    :: Mp   ! Particule mass
      real                                    :: Vp !Particule volume
      real                                    :: rho_b=1388. !kg/m3
      real                                    :: rnd(5)
      integer, parameter                      :: nmax = 4
      real                                    :: size(2**nmax)
      integer                                 :: icut, n, i
      real                                    :: l1,l2
      real                                    :: Alg_v=2.E-16  !!! [m3]
      real                                    :: Alg_m  !!! [kg]
      real                                    :: alg_death=4.5139E-6   ! death rate of algaes [1/s]
      real                                    :: phyto ! [Algae nb/m3] phytoplankton concentration that will be replaced by the actual concentration from copernicus
      real, parameter                         :: mydyn      = 1.002e-3   ! [kg/m/s] dynamic viscosity of water (at 20 oC)
      real                                    :: temp !, pressure, rho, salt
      real                                    :: K_boltz =1.38064E-23   ! Boltzmann constant [m2 kg/s2 K]
      real                                    :: diff_P, diff_A ! diffusivity of plastic particle and algae [m2/s]
      real                                    :: shear = 2 ! [1/s]
      real                                    :: encounter
      real                                    :: growth
      real                                    :: grazing, P_volume
      real                                    :: respiration
      real                                    :: r_a! phytoplancton radius  [m]
      real                                    :: B_A! encouter kernel rate = collision frequency [m3/s]
      real                                    :: r_tot !!! radius of particle [m]
      integer                                 :: istat, julday
      real                                    :: K_cap !! carrying capacity on the plastic
      real                                    :: Chl_a! =5E-6 ! [mg/m3] phytoplankton concentration that will be replaced by the actual concentration from
      real                                    :: v_active(3)
      real                                    :: f_CA=2726E-9! [kg carbon/ cell] carbon/cell ration
      real                                    :: pi=3.1415
      real                                    :: Z,t,p ! depth
      real                                    :: k=0.2 ! exctinction coeff of water
      real                                    :: teta, I_o , I_z ! depth
      real                                    :: alpha=0.12 ! intitial slope
      real                                    :: mu_max =1.85 ! maximal growth rate of algae [per day]
      real                                    :: Mu_tot,Mu_opt, day ! growth rate (total and optimal) (depend on the light and temperature [per day]
      real                                    :: depth,dist,I_max
      real                                    :: B_Aa, B_Ab, resp
      real                                    :: carbon_tot,k1(3)
      real                                    :: f_Chla_C ! =0.03!ratio of the total organic carbon during exponential growth
      real                                    :: geo(3),xy2(2),xy1(2)
c      real                                    :: chla_m= 893.51e-6![kg_chla/mmol]
      real                                    :: phyto_C,v(3)
      real                                    :: I_opt=1.8e13
      type(clock),pointer                     :: current_time
      integer                                  :: isec

c ------------    get traveled distance --------------------------

c       xy2(1)=xyz(1)
c       xy2(2)=xyz(2)
      call get_tracer_position(space, xyz)
c       xy1(1)=xyz(1)
c       xy1(2)=xyz(2)
c      call get_horizontal_distance(xy1, xy2, dist)
c      state%distance=state%distance+dist
      call interpolate_chlorophyl(xyz, Chl_a, stat)
      call interpolate_temp(xyz,   temp,  istat)
      call get_active_velocity(state, space, v_active)
      depth            = xyz(3)
      state%age      = state%age + dt/86400.
c   ----- light ---
      current_time => get_master_clock()
      call get_julian_day(current_time, julday)
      day            = 2.*pi*julday/365.-pi/2.
      I_max          = 4.6*60.*60.*24.*((sin(day)+1.)*375.+ 100.)
      p=(sin(day)+1.)*3./24.+ 9/24.
      call get_second_in_day(current_time,isec)
      t              = isec/(24.*60.*60.)
      if (t<=p) then
       I_o=I_max*sin(pi*t/p)
      else
       I_o=0
      endif
      I_z            = I_o*exp(k*depth)
      mortality_rate = 0.0! for the particle. not the biofouling organisms
c     --------------- algae concentration parameters ---------------
      r_a            = (Alg_v*3/(4*3.14))**0.333
      Alg_m          = Alg_v*rho_b
c      phyto_C        = Chl_a*1.e-6*117.*12./(2.*14.)! use of chl because chla is declared as (diatoms+flagelates+cyanobacteria )*2
c      carbon_tot     = Chl_a*chla_m/f_Chla_C! total organic carbon (carbon during exponential growth)
      f_Chla_C       =0.003+1.0154*exp(0.050*temp)*
     +                exp(-0.059*I_z*1.0e-6)
c      write(*,*)  I_z,temp,f_Chla_C
      carbon_tot     = Chl_a/f_Chla_C! total organic carbon (kg/m^3)
      phyto          = carbon_tot/f_CA ![cell/m3]
c       phyto2          =phyto_C/f_CA
c        write(*,*)  I_o
c     --------------- if settlement of particles     ---------------
      if (state%age < settle_period_start) then
      die            = .false.
      next           = .false.
      elseif ((state%age > settle_period_start).and.
     +        (state%age < settle_period_end)) then
      die            = .false.
      next           = .true.   ! try settling
      else ! missed settlement
      die            = .true.
      next           = .false.
      endif
c     ---------------- calculate the biofilm growth ------------

       K_cap         = 15./(pi*r_a**2.) ! maximum nb of algae per m^2
       V_settling   = v_active(3)
c ---------------- encounter -------------------------------
       P_volume=state%nb_algae*Alg_v*P_surf+PIV
      r_tot=(3.*P_volume/(4.*pi))**0.333
      B_Aa  = 3.14 * r_tot**2. *abs(V_settling) /2.!settling collision frequency
      B_Ab  =1.3 * shear * (r_tot + r_a)**3.
      B_A   =B_Aa+B_Ab !kernel encounter rate= brownian + settling + shear
c      write(*,*) "diff_p, B_A", diff_p, B_A, mydyn, r_tot
      encounter  = B_A * phyto
c       write(*,*) encounter
c        write(*,*)   Chl_a,carbon_tot, phyto_C, phyto1, phyto2,
c     +encounter
c ----------------- Growth ---------------------------------
c   ----- temp ---
      teta       = 1.5**((temp-20.)/10.)
      Mu_opt     = mu_max*I_z* alpha/(alpha*I_z+mu_max*
     +              ((I_z/I_opt)-1.)**2.)
      Mu_tot     = Mu_opt*teta/(24.*60.*60.)
      growth     = state%nb_algae * Mu_tot * (1-state%nb_algae/K_cap)!* G_max! growth on the plastic is proportional to the ambiant primary production
      grazing    = alg_death * state%nb_algae
      resp       =(0.1*state%nb_algae* 2.**((temp-20.)/10.))
     +          /(24.*60.*60.)
      state%nb_algae = state%nb_algae + (encounter/P_surf + growth
     +                  - grazing-resp)*dt!
c       write(*,*)  encounter, growth, grazing, resp,
c     +          state%nb_algae
c -------------get distance traveled -----------------------
c      call interpolate_currents(xyz,k1,istat2) ! real(3) in physical units
c      v=k1+v_active
c      Dist_h=sqrt(v(1)**2.0+v(2)**2.0)*dt
c      Dist_v=v(3)*dt
c      state%distance=state%distance+Dist_h
c      write(*,*) v_active,k1, v,dt, Dist_h, state%distance
c     ----------------------------------------
      end subroutine update_particle_state_

      end module


