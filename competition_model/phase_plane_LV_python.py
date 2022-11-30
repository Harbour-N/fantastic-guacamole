# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 17:48:04 2022

@author: nicho
"""

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

# all the parateres taken from Wikipidia page on LV equations
alpha = 1.1  # prey growth rate 
beta = 0.4    # prey death rate 
gamma = 0.4 # predator death rate
delta = 0.1 # predator growth rate


# define system in terms of a Numpy array
def Sys(X, t=0):
    # here X[0] = x and x[1] = y    
    return np.array([ alpha*X[0] - beta*X[0]*X[1],
                     delta*X[0]*X[1] - gamma*X[1] ])


# generate 1000 linearly spaced numbers for x-axes
t = np.linspace(0, 100, num = 10000)
# initial values: x0 = 10, y0 = 10
Sys0 = np.array([10, 10])

# 1st figure plot evolution in time and in phase space side by side
det_solution = scipy.integrate.odeint(Sys, Sys0, t)
x = det_solution[:,0]
y = det_solution[:,1]

fig = plt.figure(figsize=(15,5))
fig.subplots_adjust(wspace = 0.5, hspace = 0.3)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.plot(t,x, '-', label='prey')
ax1.plot(t,y, '--',label='predator')
ax1.set_title("Dynamics in time")
ax1.set_xlabel("time")
ax1.grid()
ax1.legend(loc='best')

ax2.plot(x, y, color="blue")
ax2.set_xlabel("x")
ax2.set_ylabel("y")  
ax2.set_title("Phase space")
ax2.grid()
plt.show()

# 2nd figure plot phase space as quiver plot
x = np.linspace(0, 30, 20)
y = np.linspace(0, 10, 20)

X1 , Y1  = np.meshgrid(x, y)                    # create a grid
DX1, DY1 = Sys([X1, Y1])                        # compute growth rate on the grid
M = (np.hypot(DX1, DY1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M

plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid')

plt.legend()
plt.grid()
plt.show()

# 3rd figure plot phase space as stream plot showing trajectories more clearly
plt.streamplot(X1, Y1, DX1, DY1)
plt.show()





