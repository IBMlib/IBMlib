      real function rho_UNESCO(ss,tt,pp)
c-----------------------------------------------------------------------
c
c     Function that calculates the density of sea water at a given point
c     using the UNESCO equation of state (from Per Berg, per@dmi.dk)
c
c     Output:
c       local water density in [kg/m**3]
c     
c     Input:
c       ss = salinity    [PSU]
c       tt = temperature [Celcius]
c       pp = pressure    [bar] 1 bar = 100 000 pascals (Pa) 
c
c     Comments
c       Re-written from the original F77 function in cmodlib
c       to F90 such that it is ready for in-lining if we choose to do so
c       converted from free format to fixed form
c       replaced real(8) -> real
c
c-----------------------------------------------------------------------
      
      real, intent(in) :: ss
      real, intent(in) :: tt
      real             :: pp
      real :: s1, s2, s3, s15, t1, t2, t3, t4, t5, rhonul, sbm
c**********      coefficients   **************************************
      real, parameter :: a0 = 999.842594d00                            
      real, parameter :: a1 =   6.793952d-2                            
      real, parameter :: a2 =  -9.095290d-3                            
      real, parameter :: a3 =   1.001685d-4                            
      real, parameter :: a4 =  -1.120083d-6                            
      real, parameter :: a5 =   6.536332d-9                            
      real, parameter :: b0 =  +8.24493d-1                             
      real, parameter :: b1 =  -4.0899d-3                              
      real, parameter :: b2 =  +7.6438d-5                              
      real, parameter :: b3 =  -8.2467d-7                              
      real, parameter :: b4 =  +5.3875d-9                              
      real, parameter :: c0 =  -5.72466d-3                             
      real, parameter :: c1 =  +1.0227d-4                              
      real, parameter :: c2 =  -1.6546d-6                              
      real, parameter :: d0 =  +4.8314d-4                              
      real, parameter :: e0 = 19652.21d00                              
      real, parameter :: e1 =   148.4206d00                            
      real, parameter :: e2 =    -2.327105d00                          
      real, parameter :: e3 =     1.360477d-2                          
      real, parameter :: e4 =    -5.155288d-5

      real, parameter :: f0 =    54.6746d00                            
      real, parameter :: f1 =    -0.603459d00                          
      real, parameter :: f2 =     1.09987d-2                           
      real, parameter :: f3 =    -6.1670d-5                            
      real, parameter :: g0 =     7.944d-2                             
      real, parameter :: g1 =     1.6483d-2                            
      real, parameter :: g2 =    -5.3009d-4                            
      real, parameter :: h0 =  3.239908d00                             
      real, parameter :: h1 =  1.43713d-3                              
      real, parameter :: h2 =  1.16092d-4                              
      real, parameter :: h3 = -5.77905d-7                              
      real, parameter :: i0 =  2.2838d-3                               
      real, parameter :: i1 = -1.0981d-5                               
      real, parameter :: i2 = -1.6078d-6                               
      real, parameter :: j0 =  1.91075d-4                              
      real, parameter :: k0 =  8.50935d-5                              
      real, parameter :: k1 = -6.12293d-6                              
      real, parameter :: k2 =  5.2787d-8                               
      real, parameter :: l0 = -9.9348d-7                               
      real, parameter :: l1 =  2.0816d-8                               
      real, parameter :: l2 =  9.1697d-10

c*********************************************************************
      s1 = ss
      s2 = ss*ss
      s3 = s2*ss
      s15= sqrt(s3)
      t1 = tt
      t2 = tt*tt
      t3 = t2*tt
      t4 = t2*t2
      t5 = t4*tt
 
      rhonul = a0+a1*t1+a2*t2+a3*t3+a4*t4+a5*t5                           
     +     + (b0+b1*t1+b2*t2+b3*t3+b4*t4)*s1                             
     +     +(c0+c1*t1+c2*t2)*s15                                        
     +     + d0*s2

      sbm    =   e0+e1*t1+e2*t2+e3*t3+e4*t4                                 
     +     +(f0+f1*t1+f2*t2+f3*t3)*s1                                   
     +     +(g0+g1*t1+g2*t2)*s15                                        
     +     +(  h0+h1*t1+h2*t2+h3*t3                                     
     +     +(i0+i1*t1+i2*t2)*s1                                       
     +     + j0*s15                  ) * pp                          
     +     +(  k0+k1*t1+k2*t2                                           
     +     +(l0+l1*t1+l2*t2)*s1      ) * pp * pp

      rho_UNESCO = rhonul/(1.e0-pp/sbm)

      end  function rho_UNESCO
