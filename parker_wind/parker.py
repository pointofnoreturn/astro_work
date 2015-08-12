#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as mpl


# little Newton-Raphson solver to solve
# the transcendental equation for v/c_s
def parker_f(vbar,rbar,C):
    return vbar**2 - np.log(vbar**2) - 4.0*np.log(rbar) \
        - 4.0/rbar - C

def parker_dfdvbar(vbar):
    return 2.0*vbar - 2.0 / vbar

def parker_err(vbar,rbar,C):
    return np.abs( parker_f(vbar,rbar,C) / (4.0*np.log(rbar) \
        + 4.0/rbar + C))

def parker(rbar,C,vbar_guess):
    tol = 1.0e-10
    # handle bifurcation point at r = r_c
    if rbar > 1.0:
        vbar = np.fmax(vbar_guess,1.0+tol)
    else:
        # also can't have vbar = 0 (log (0) = bad)
        vbar = np.fmax(vbar_guess,tol)
    print vbar
    it = 0
    while parker_err(vbar,rbar,C) > tol:

        dvbar = - parker_f(vbar,rbar,C) / \
               parker_dfdvbar(vbar)

        # Limit changes of vbar to be no larger
        # than 20% per iteration step.
        # This turns out to be neccessary to
        # keep the solution from wandering off into
        # bad territory (like vbar < 0).
        fac1 = np.fmin(0.2/np.abs(dvbar/vbar),1.0)
    
        vbar = vbar + fac1 * dvbar

# debug output
#        print it, fac1, parker_f(vbar,rbar,C), parker_err(vbar,rbar,C), vbar
        it = it+1

    return vbar


# integration constant that gives us the Parker V solution
C=-3

# make a grid of r/r_c of 1000 points
N = 1000
rad = np.linspace(0.01,5,N)
vel = np.zeros(N)

# fill in the points with solutions for
# v/c_s as a function of r/r_c
vguess = 1.0e-10
for i in range(N):
    vel[i] = parker(rad[i],C,vguess)
    vguess = vel[i]

mpl.plot(rad,vel)
mpl.show()


