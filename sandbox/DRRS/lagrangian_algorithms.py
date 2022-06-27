#!/usr/bin/env python
#####################################################################################
#  Assess Lagrangian algorithms for the 1D pollution plume problem with constant
#  diffusivity D, advection velocity v and sink rate lam:
#  Compare DRRS (dynamical renormalization resampling scheme) to conventional Lagrangian
#
#   dP/dt = -dJ/dx - lam*P
#       J = -D*dP/dx + v*P
#      J0 = s              (flux boundary condition at x=0)
#
#  Implicit interior domain 0 < x < 1
#
#  Staionary solution (down wind, k>0):
#    P = s*exp(-k*x) / (v+D*k)
#    k = -v/(2D) + sqrt[(v/(2D))^2 + lam/D]
#
#  If boundary condition at x=1 are enforced, the upwind branch (k<0) gets excited as well
#  ----------------------------------------------------------------------------------
#  Points to stress:
#    1) tau does not affect the dynamics of DRRS; it is a posterior plot-time parameter
#    
#####################################################################################

from numpy    import *   
from numpy.linalg import *
from numpy.random import *
from scipy.special import erfc


def err_metric(y,yref):
     return sqrt(mean((y-yref)**2))

 
def equilibrium_solution(x, s, lam, v, dh):
    # -----------------------------------------------------------
    # Evaluate analytical equilibrium solution on array x
    # -----------------------------------------------------------
    if dh==0:
        if v==0:
            raise ValueError("singular limit dh==0 and v==0 not defined")
        else: # dh==0 and v>0
            k = lam/v
    else: # dh>0
        z = v/2/dh
        k = -z + sqrt(z**2 + lam/dh)
    #
    A = s/(v+dh*k)
    return A*exp(-k*x)




# -------- integral of unnormalized probability density erfc((mu-x)/sigma/sqrt(2))  --------

px    = lambda z,m,s: erfc((m-z)/s/sqrt(2))/2
intpx = lambda z,m,s: ((-exp(-(1 + m)**2/(2.*s**2)) + exp(-(m - z)**2/(2.*s**2)))*sqrt(2/pi)*s + (1 + m)*erfc((1 + m)/(sqrt(2)*s)) + (-m + z)*erfc((m - z)/(sqrt(2)*s)))/2.


# -----------------------------------------------------------------------------------
# Avoid slow scipy.stats.rv_continuous methods, using vectorized analytic CDF and PDF
# and Newton-Raphson for solving r==CDF(x)
# Distribution has support on [-1,0]
# -----------------------------------------------------------------------------------
class backtail_fast:
    _NRsingular = 1e-12 # Stabilize Newton-Raphson cycles for flat tails in CDF
    _zmin       = -1
    _zmax       =  0
    def __init__(self, mu, sigma):
        self.mu    = mu
        self.sigma = sigma
        self.nofac = intpx(0, mu, sigma)
    def _pdf(self, x):        
        return px(x, self.mu, self.sigma)/self.nofac
    def _cdf(self, x):
        return intpx(x, self.mu, self.sigma)/self.nofac
    def draw(self, n, tol=1e-12, itmax=30):
        # ----------------------------------------------------------
        # Solve random(n) == CDF(x) by Newton-Raphson cycles
        # Stabilize Newton-Raphson cycles for flat tails in CDF
        # ----------------------------------------------------------
        r     = random(n)
        z     = 0.5*(self._zmin+self._zmax)*ones(n, float) # mid point in support interval
        it    = 0
        dzmax = 1e20
        while dzmax > tol:
            dz    = (r-self._cdf(z))/(self._pdf(z)+self._NRsingular) 
            dzmax = amax(abs(dz))
            z    += dz
            z     = where(z>self._zmin, z, self._zmin) # enforce support interval
            z     = where(z<self._zmax, z, self._zmax) # enforce support interval
            it   += 1
            #print("iter %d : dzmax = %f" % (it, dzmax))
            if it>itmax:
                raise RuntimeError("backtail_fast:itmax exceeded")
        return z
    

    
class standard_IBM:
    # ------------------------------------------------------------------------------
    # Standard IBM without resampling of 1D down stream pollution plume
    #
    # Currently no provision for right-end influx (from background pollution),
    # which may become important in the diffusive limit v->0.
    # For simplicity do not explicitly remove particles exiting the interior domain 0 < x < 1
    # 
    # np: total number of particles on stack
    # Pint0: initial sum of plastic items dispersed on interior interval(rest withheld for source emission)
    # ppp: fixed conversion ratio: plastic items per particle
    # ------------------------------------------------------------------------------
    def __init__(self, s, lam, v, dh, np, Pint0, ppp):
        self.t   = 0.
        self.s   = s    # source strength (plastic items emitted per unit time)
        self.lam = lam  # removal rate
        self.v   = v    # advective velocity > 0
        self.dh  = dh   # horizontal diffusivity
        self.ppp = ppp  # fixed conversion ratio: plastic items per particle
        n0 = round(Pint0/ppp)  # corresponding initial number of particles
        assert np >= n0 > 0
        self.x   = ma.array(zeros(np, float), mask=np*[True]) 
        self.x.mask[:n0] = False           # activate initial distribution
        self.x.data[:n0] = uniform(0,1,n0) # scatter n0 over inner domain
        self.next        = n0              # index of next tracer to be activated

    def spatial_step(self, dt):
        # -------------------------------------------------------
        # spatial propagation of tracers incl timer increment    
        # -------------------------------------------------------
        np = len(self.x) # all active       
        self.x += self.v*dt                                       # advection, only unmasked elements modified        
        self.x += normal(0, 1, np)*sqrt(2*self.dh*dt)             # diffusion, only unmasked elements modified       
        #self.x.data[:] = where(self.x.data > 0, self.x.data, 0.0) # left end clip BC (pile-up)
        self.x.data[:] = abs(self.x.data)                         # reflective left end clip BC
        self.t += dt
        
    def forward_step(self, dt):
        self.spatial_step(dt)
        # ------ release new particles from source ------
        nemit = round(self.s*dt/self.ppp)
        assert self.next+nemit <= len(self.x) # check stack is not exhausted
        self.x.mask[self.next:self.next+nemit] = False
        self.next += nemit
        # ------ deactivate lost particles according to rate lam ------
        iactive = flatnonzero(self.x.mask[:self.next]==False)  # indices of particles than can be removed
        remfrac = 1-exp(-self.lam*dt)                          # finite step expression
        unlucky = flatnonzero(random(len(iactive)) < remfrac)  # indices of particles in iactive list to remove
        remove  = take(iactive, unlucky)                       # indices in full list of particles to remove
        put(self.x.mask, remove, True)                         # mask these elements
        
        
    def eulerian_cast(self, edges, *args):
        # compute concentration of plastic items on provided grid edges
        p = histogram(self.x.data[self.x.mask==False], edges)[0]
        dx = edges[1:]-edges[:-1]
        return self.ppp * p/dx



    
class DRRS_IBM(standard_IBM):
    # -----------------------------------------------------------------------------------------------------
    # DRRS (dynamical renormalization resampling scheme) resampling IBM for 1D down stream pollution plume
    #
    # Currently no provision for right-end influx (from background pollution),
    # which may become important in the diffusive limit v->0.
    # 
    # np: total number of particles on stack
    # tau: smoothing time scale for computing plastic per particle ratio ppp
    # -----------------------------------------------------------------------------------------------------
    def __init__(self, s, lam, v, dh, np):
        self.t   = 0.
        self.s   = s    # source strength (plastic items emitted per unit time)
        self.lam = lam  # removal rate
        self.v   = v    # advective velocity > 0
        self.dh  = dh   # horizontal diffusivity
        self.x   = ma.array(uniform(0,1,np), mask=False) # scatter np over inner domain
        self.resample = {"time": [1.0*self.t], "nemit":[]} # just add starting time
        
    def forward_step(self, dt):
        self.spatial_step(dt)
        # ---- resample particles removed or leaving interior
        interior = flatnonzero(self.x.data < 1)                  # indices of interior particles than can be removed
        remfrac  = 1-exp(-self.lam*dt)                           # finite step expression
        unlucky  = flatnonzero(random(len(interior)) < remfrac)  # indices of particles in interior list to remove
        remove   = take(interior, unlucky)                       # indices in full list of particles to remove
        exterior = flatnonzero(self.x.data >= 1)                 # indices of particles that left interior, all resampled
        nemit    = len(unlucky) + len(exterior)
        self.resample["time"].append(1.0*self.t)
        self.resample["nemit"].append(nemit)
        put(self.x.data, remove, 0.0)                            # emit resampled at source point
        put(self.x.data, exterior, 0.0)                          # emit resampled at source point
        #
        
    def eulerian_cast(self, edges, tau):
        # compute concentration of plastic items on provided grid edges
        # set ppp and invoke parent cast
        # tau: smoothing time scale for computing plastic per particle ratio ppp
        if len(self.resample)==0:
            if hasattr(self, "ppp"):
                print("eulerian_cast: using preset ppp = %f" % self.ppp)
            else:
                raise RuntimeError("resample statistics void, but no plastic-per-particle set ")
        else: # compute conversion ratio plastic-per-particle 
            t     = array(self.resample["time"])
            nemit = array(self.resample["nemit"])
            dt    = t[1:]-t[:-1] # elapsed time since last emission event, same length as nemit
            w     = exp(t[1:]/tau) # exclude starting time, which does not have emission
            w    /= sum(w)
            avg_par_emit_rate = sum(w*nemit/dt)
            self.ppp = self.s/avg_par_emit_rate
        return standard_IBM.eulerian_cast(self, edges)


class DRRS_IBM_backdiffcorr(DRRS_IBM):
    # -----------------------------------------------------------------------------------------------------
    # DRRS (dynamical renormalization resampling scheme) resampling IBM for 1D down stream pollution plume,
    # with back-diffusion correction at right boundary x=1
    # time step dt must bew fixed, to avoid recompiling back-diffusion generator
    # -----------------------------------------------------------------------------------------------------
    def __init__(self, s, lam, v, dh, np, dt):
        self.dt       = dt
        self.backdiff = backtail_fast(v*dt, sqrt(2*dh*dt)) # compile once
        DRRS_IBM.__init__(self, s, lam, v, dh, np)
        
    def forward_step(self, dt):
        assert dt == self.dt # allow ignored argument for algorithmic comparison
        self.spatial_step(dt)
        # ---- resample particles removed or leaving interior
        
        interior = flatnonzero(self.x.data < 1)                  # indices of interior particles than can be removed
        remfrac  = 1-exp(-self.lam*dt)                           # finite step expression
        unlucky  = flatnonzero(random(len(interior)) < remfrac)  # indices of particles in interior list to remove
        remove   = take(interior, unlucky)                       # indices in full list of particles to remove
        exterior = flatnonzero(self.x.data >= 1)                 # indices of particles that left interior, all resampled
        leaving  = concatenate((remove, exterior))               # particles leaving
        n0       = len(exterior)/(self.v*self.dt)                # estimator on particle concentration at right boundary
        nback    = round(float(0.5*n0*self.backdiff.nofac))      # cast to float (since round(numpy.float64) gives numpy.float64)
        #
        assert nback <= len(leaving)
        nemit    = len(leaving) - nback                          # at left source
        self.resample["time"].append(1.0*self.t)
        self.resample["nemit"].append(nemit)
        #
        xback = 1 + self.backdiff.draw(nback)                    # 1: add offset, since backdiff referred to left bound
        xback = where(xback>0, xback, 0.0)                       # enforce left boundary
        put(self.x.data, leaving[:nback], xback)                 # back diffuse nback from right side
        put(self.x.data, leaving[nback:], 0.0)                   # emit rest at source point
        #
    
#########################################################################
if __name__ == "__main__":
    s   = 200000 # plastic items per time unit
    lam = 0.5
    v   = 0.2
    dh  = 0.01
    ppp = 200
    bins = linspace(0,1,100)
    cpts = 0.5*(bins[1:]+bins[:-1])
    np   = 10**6
    nt   = 1000
    tau  = 0.2
    dt   = 0.1
    #ens  = standard_IBM(s, lam, v, dh, np, s/lam, ppp)
    #ens  = DRRS_IBM(s, lam, v, dh, np)
    ens  = DRRS_IBM_backdiffcorr(s, lam, v, dh, np, dt)
    for i in range(nt):
        ens.forward_step(dt)
        #f = open("jjj%d" % i, "w")
        #for xy in zip(cpts, ens.eulerian_cast(bins, tau)):
        #    f.write("%f %f\n" % xy)
        print(i,sum(ens.eulerian_cast(bins, tau)))
        #f.close()
    f = open("jjjlast", "w")
    conc = ens.eulerian_cast(bins, tau)
    for (x,y) in zip(cpts, conc):     
        ya = equilibrium_solution(x, s, lam, v, dh)
        f.write("%f %f %f\n" % (x,y,ya))
    f.close()
    
