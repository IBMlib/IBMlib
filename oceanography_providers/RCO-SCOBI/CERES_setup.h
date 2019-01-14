c     ------------------------------------------------------------
c     Data descriptors in relation to RCO-SCOBI setup
c     provided by SMHI in relation to the CERES project
c     Data entry mapping, see /home/data/RCO-Scobi_test/modvars.h
c     ------------------------------------------------------------
      real, parameter    :: layer_width = 3.0 ! meters, same for all layers
      integer, parameter :: nz_expect   = 83  !
c
      integer, parameter :: levnum_dslm = 1
      integer, parameter :: levnum_t    = 4
      integer, parameter :: levnum_s    = 5
      integer, parameter :: levnum_u    = 8
      integer, parameter :: levnum_v    = 9
c       
      integer, parameter :: levnum_zoo  = 10   
      integer, parameter :: levnum_diat = 11         
      integer, parameter :: levnum_flag = 12 
      integer, parameter :: levnum_cyan = 13
      integer, parameter :: levnum_orgd = 14
      integer, parameter :: levnum_ammo = 15
      integer, parameter :: levnum_nit  = 16
      integer, parameter :: levnum_phos = 17
      integer, parameter :: levnum_oxy  = 18
c
