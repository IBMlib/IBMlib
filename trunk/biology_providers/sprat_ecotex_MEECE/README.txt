===========================================================================
Parameterization of copper impact on early life stages (egg+yolksac) of sprat, 
based on experiments of Tanja St. John within EU FP7 MEECE

Parameterization covers:
   1) The decrease in egg survival (hatch fraction)
   2) The decrease in larval hatch length
including temperature dependence. Experiments only cover egg+yolksac stages
and only addresses egg survival and hatch length (no other impact mechanisms)
===========================================================================

1) The decrease in egg survival (hatch fraction) in the meta parameterization is represented as 

          S(c,T) =  S0(T) exp(- c / K2(T))

          S0  = A0*(1 - exp(-A1*(T-A2)))
          K2  = 50 yg/l 
	  A0  = 0.231794
          A1  = 0.286788
          A2  = 5.68092

with T as degC and copper concentration c as yg/l.
I have skipped copper-temperature interaction, because variability of K2(T)
from experiments seem spurious, and I fixed K2  = 50 yg/l

This corresponds to a mortality formulation as 
           S(c,T) =  exp( -Z(c,T))
---------------------------------------------------------------------------
2) The decrease in larval hatch length in the meta parameterization is represented as 

 Lh(c,T) = L0 exp(- c / K1(T)) + Linf

            L0   = 7 mm
            Linf = 5 mm
            K1   = 0.00049 * T^5.95 yg/l

with T as degC and copper concentration c as yg/l

As you note in your draft, toxiticity increases with lower temperature
You also may want to note in your paper that if you extrapolate 
temperature-copper interaction below experimental temperatures,
current Baltic average copper concentration (c = 0.1-0.6 yg/l)
becomes toxic below Tx = (0.6/0.00049)**(1.0/5.95) = 3.3 deg Celcius.
This extrapolation is of course sensitive to the power law assumption
when fitting K1(T), but it may nevertheless be interesting in the discussion to
play with the existence of such a toxicity temperature Tx at current copper levels.

