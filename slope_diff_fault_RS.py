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

P0 = 0.0001 # soil prod. rate @ H=0, m/yr
kappa = 0.1 # diffusivity, m**2/yr
hstar = 0.5 # length scale of soil production, m
denrck = 2700. # rock density, kg/m**3
denreg = 2000. # regolith density, kg/m**3
elow = 0.0001 # m,yr lowering rate of streams @ boundary condition
slip = 0.03 # m/yr, slip rate of fault

xmax = 40 # boundaries
 # m
dx = 1 # m
N = int(((xmax*2)/dx)+1) # number of intervals

x = linspace(-xmax,xmax,N)

tmax = 10**6. # yr
dt = 10**3 # yr

t = arange(0,tmax,dt) # setting up time array

z = np.copy(x)
z[0:len(z)/2] *= 0 # setting up initial slope
z[len(z)/2:len(z)] *= 0
H = 10*np.ones(len(z)) # regolith thickness, m
b = z - H # basement, m

dzdx = diff(z)/dx # finding slope 

#print dzdx
q = np.zeros(len(dzdx)) # 
dqdx = np.zeros(len(dzdx-1)) # initial flux is 0

plt.close()

"""
Run
"""
# w, production rate of soil, m/y
# b, change in basement
# z, total height of slope
fig = plt.figure()

line1, = plt.plot(x,z,label='initial state')

for index, i in enumerate(t):
    
    w = P0*np.exp(-H/hstar)
    b -= w*dt
    
    
    dzdx = diff(z)/dx
    q = - kappa * dzdx
    dqdx = diff(q)/dx
    
    H[1:N-1] += (denrck/denreg)*w[1:N-1]*dt - (1/denreg)*dqdx*dt
    z = b + H 
    z[0] -= elow*dt 
    z[len(z)-1] -= elow*dt
    z[0:len(z)/2] -= slip*dt # moving fault upwards
    
line2, = plt.plot(x,z,label='final state') # plotting
plt.xlabel('distance (m)')
plt.ylabel('height (m)')
plt.title('Evolution of hillslope over time')

plt.legend(handles=[line1, line2])

plt.show()
    
    
    
    