#!/usr/bin/env python

#Author: Craig Lage, NYU; 
#Date: 11/17/09

'''
This program manages the running of the euler.c 2D fluid flow program, now with added gravitational forces.

'''

from pylab import *
from subprocess import *
import time
import matplotlib.pyplot as plt   # For plotting graphs.
matplotlib.use('Agg')
import subprocess                 # For issuing commands to the OS.
import os
import sys                        # For determining the Python version.


def cploty(file,plotcount):

    infile = open(file)

    x=zeros([NBx,NBy])
    y=zeros([NBx,NBy])
    rho=zeros([NBx,NBy])    
    vx=zeros([NBx,NBy])    
    vy=zeros([NBx,NBy])    
    P=zeros([NBx,NBy])
    phi=zeros([NBx,NBy])

    for line in infile:

        strval = line.strip().split()
        xtemp = float(strval[0])
        ytemp = float(strval[1])
        i = round((xtemp-Bxmin)/dx)
        j = round((ytemp-Bymin)/dx)
        x[i,j] = max(min(xtemp,Bxmax-dx/2),Bxmin+dx/2)
        y[i,j] = max(min(ytemp,Bymax-dx/2),Bymin+dx/2)
        rho[i,j]= float(strval[2])
        vx[i,j] = float(strval[3])
        vy[i,j] = float(strval[4])
        P[i,j] = float(strval[5])
        phi[i,j] = float(strval[6])

    plt.contourf(x,y,rho,levels=arange(0,20,1))
    filename = str('%03d' % plotcount) + '.png'
    plt.savefig(filename, dpi=100)
    print 'Wrote file', filename
    plt.clf()
    plotcount=plotcount+1

    return plotcount

#****************MAIN PROGRAM*****************
k=1
g=0.0
NBx=100 # Array size
NBy=100 # Array size
NPx=2 # Number of processors in X-direction
NPy=2 # Number of processors in Y-direction
Nsteps=5000 # of time steps
Outsteps=100 # Will save output when time step counter is a multiple of this
dx=.02
dy=.02
Bxmin=-1.0
Bymin=-1.0
gamma=1.4
theta=1.5
CFL=0.5
alpha=1.8
G=-10.0

BC=0 # 0 - Periodic, 1 - No backflow, 2 - Reflective

Bxmax = Bxmin+NBx*dx
Bymax = Bymin+NBx*dx

plotcount=1

plotcount=cploty('input.dat',plotcount)

for timestep in range(100, 30000, 100):

    plotcount=cploty('output'+str(timestep)+'.dat',plotcount)


command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'rho.avi')

#os.spawnvp(os.P_WAIT, 'mencoder', command)

print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)

print "\n\n The movie was written to 'output.avi'"

print "\n\n You may want to delete *.png now.\n\n"
print 'Initializing data set...'   # Let the user know what's happening.


#************END MAIN PROGRAM*************************

