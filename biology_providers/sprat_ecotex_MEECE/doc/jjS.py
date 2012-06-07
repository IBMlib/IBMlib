#!/usr/bin/env python
from numpy import *

A0  = 0.231794
A1  = 0.286788
A2  = 5.68092
K2  = 50
T   = 12
L0 = 7
Li = 5

for c in arange(0,20,0.1):
    S =  A0*(1 - exp(-A1*(T-A2)))*exp(- c / K2)
    K1 = 0.00049 * T**5.95 
    Lh = L0 * exp(- c / K1) + Li
    #print c,S
    print c,Lh
