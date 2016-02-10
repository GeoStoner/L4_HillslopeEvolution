# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:28:01 2016

@author: ryanstoner

for modeling class lab 3
"""

"""
Initialize
"""
import numpy as np # importing math functions
from numpy import *

import matplotlib.pyplot as plt # importing plotting function
plt.close()

P0 = 0.0001 # soil prod. rate @ H=0, m/yr
kappa = 0.05 # diffusivity, m**2/yr
hstar = 0.5 # length scale of soil production, m
denrck = 2700. # rock density, kg/m**3
denreg = 2000. # regolith density, kg/m**3
elow = 0.00000002 # m, yr lowering rate of streams @ boundary condition
slip = 0.0003 # m/yr, slip rate of fault

xmax = 50 # boundaries
 # m
dx = 1 # m
N = int(((xmax*2)/dx)+1) # number of intervals

x = linspace(-xmax,xmax,N)

tmax = 1*10**6. # yr
dt = 10**2 # yr

t = arange(0,tmax,dt) # setting up time array
Nplots = 100 # number of plots
tplot = tmax/Nplots # interval between plots

z = np.zeros(len(x))

H = 5*np.ones(len(z)) # regolith thickness, m
b = z - H # basement, m

dzdx = diff(z)/dx # finding slope 

#print dzdx
q = np.zeros(len(dzdx)) # 
dqdx = np.zeros(len(dzdx-1)) # initial flux is 0


"""
Run
"""
# w, production rate of soil, m/y
# b, change in basement
# z, total height of slope
fig = plt.figure()

line1, = plt.plot(x,z,label='initial state')

for i in range(len(t)):
    plt.ion()
    
    if (t[i+1] % tplot)==0: # clear plot, t[i+1] so it won't clear last plot
        line2.remove()
    
    w = P0*np.exp(-H/hstar)
    b -= w*dt
    b[0:len(z)/2] -= slip*dt # moving left part of fault downwards
    b[len(z)/2:len(z)-1] += slip*dt # moving right part upwards
    dzdx = diff(z)/dx
    q = - kappa * dzdx
    dqdx = diff(q)/dx
    
    H[1:N-1] += (denrck/denreg)*w[1:N-1]*dt - (1/denreg)*dqdx*dt
    z = b + H 
   
    z[0] -= elow*dt*i 
    
    z[len(z)-1] += slip*dt*i-elow*dt*i
    
    if (t[i] % tplot) == 0:
        line2, = plt.plot(x,z,'r-',label='final state') # plotting
        plt.xlabel('distance (m)')
        plt.ylabel('height (m)')
        plt.xlim((-xmax, xmax))
        plt.ylim((-300, 300))
        plt.title('Evolution of hillslope over time')
        
        
        plt.legend(handles=[line1, line2])
    
        plt.pause(0.001)
        plt.grid(True)
        

"""
Close
"""

# Add the analytical solution to compare

#z2 = (elow/2*kappa)*(xmax**2-x**2)-elow*dt*i
#
#plt.plot(x,z2,'g--')

    
    
    