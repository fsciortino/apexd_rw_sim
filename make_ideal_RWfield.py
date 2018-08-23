from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from bbr_particle import * 
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
from matplotlib.pyplot import cm
import subprocess
import scipy
from shutil import copyfile
import os
import math

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 18

# Copy settings from a SIMION file
version=10
path_0='./simion_Efields/ver.1'
path = './simion_Efields/ver.' + str(version)

if not os.path.exists(path):
    os.makedirs(path)
if not os.path.exists(path+'/e3d_hd.txt'):
    copyfile(path_0+'/e3d_hd.txt',path+'/e3d_hd.txt')

# Now read
file=open(path+'/e3d_hd.txt','rb')

xmin,xmax=file.readline().split()
ymin,ymax=file.readline().split()
zmin,zmax=file.readline().split()
nx,ny,nz=file.readline().split()
dx,dy,dz=file.readline().split()

xmin=float(xmin.replace('D','E'))
xmax=float(xmax.replace('D','E'))
ymin=float(ymin.replace('D','E'))
ymax=float(ymax.replace('D','E'))
zmin=float(zmin.replace('D','E'))
zmax=float(zmax.replace('D','E'))

nx=int(nx.replace('D','E'))
ny=int(ny.replace('D','E'))
nz=int(nz.replace('D','E'))
dx=float(dx.replace('D','E'))
dy=float(dy.replace('D','E'))
dz=float(dz.replace('D','E'))

tot_vals=nx*ny
xx = np.linspace(xmin,xmax,nx)
yy = np.linspace(ymin,ymax,ny)
X,Y=np.meshgrid(xx,yy)

# Create 2D toroidal E-field
amp=1.0
f_RW = 1.0e5
Phrw=1.0
V=1.0
num_electrodes=4
Vel=np.zeros(num_electrodes)
mm=1

# split one cycle into 4 times
time = np.concatenate(([0.0,],np.linspace(1.0/(f_RW*4),(1.0/f_RW) ,4)))[:-1]
Ex=np.zeros((nx,ny,len(time)))
Ey=np.zeros((nx,ny,len(time)))

for ti in range(len(time)):
    t=time[ti]

    for iele in range(num_electrodes):
        # note phasing at the end of the following expression:
        Vel[iele] = amp*np.sin(- 2.0*np.pi*f_RW*t + mm*Phrw*iele*np.pi*(2.0/num_electrodes))
        
    for i in range(nx):        
        for j in range(ny): #int(ny/2),ny):
            
            rr=np.sqrt(xx[i]**2+yy[j]**2)+1.0e-10
            #theta = np.arccos(xx[i]/rr)
            theta = math.atan2(yy[j],xx[i])
            
            for iele in range(num_electrodes):

                #Ex[i,j,ti]= - (V)*np.sin(theta)
                #Ey[i,j,ti]= + (V)*np.cos(theta)
            
                Ex[i,j,ti]+= - (mm*Vel[iele]/rr)*np.sin(theta)#* np.sin( - 2*np.pi*f_RW*t)
                Ey[i,j,ti]+= + (mm*Vel[iele]/rr)*np.cos(theta)#* np.cos(theta - 2*np.pi*f_RW*t)

# complicated trick to get vortex vector map...
#Ex[:,:int(ny/2),:] = - Ex[:,int(ny/2):,:][:,::-1,:]
#Ey[:,:int(ny/2),:] =  Ey[:,int(ny/2):,:][:,::-1,:]

# sub-sample
ss_factor = 2
skip = 1 # to avoid crossing arrows at edges
Exx=Ex[skip:-skip:ss_factor,skip:-skip:ss_factor,:]
Eyy=Ey[skip:-skip:ss_factor,skip:-skip:ss_factor,:]
X=X[skip:-skip:ss_factor,skip:-skip:ss_factor]
Y=Y[skip:-skip:ss_factor,skip:-skip:ss_factor]

# Plot            
fig2=plt.figure(figsize=(10,12))

ax11 = plt.subplot(221)
ax11.quiver(Y,X,Exx[:,:,0],Eyy[:,:,0], color='r')
ax11.set_ylabel('y [m]')
ax11.set_title(r'$t=0 \mu s$')
ax11.axis('equal')
ax11.set_xticks(np.arange(0.0,X.max(), 0.1))

ax22 = plt.subplot(222,sharex=ax11)
ax22.quiver(Y,X,Exx[:,:,1],Eyy[:,:,1], color='r')
#ax22.set_xlabel('x [m]')
#ax22.set_ylabel('y [m]')
ax22.set_title(r'$t=50 \mu s$')
ax22.axis('equal')
ax22.set_xticks(np.arange(0.0,X.max(), 0.1))

ax33 = plt.subplot(223,sharex=ax11)
ax33.quiver(Y,X,Exx[:,:,2],Eyy[:,:,2], color='r')
ax33.set_xlabel('x [m]')
ax33.set_ylabel('y [m]')
ax33.set_title(r'$t=100 \mu s$')
ax33.axis('equal')
ax33.set_xticks(np.arange(0.0,X.max(), 0.1))

ax44 = plt.subplot(224,sharex=ax11)
ax44.quiver(Y,X,Exx[:,:,3],Eyy[:,:,3], color='r')
ax44.set_xlabel('x [m]')
#ax44.set_ylabel('y [m]')
ax44.set_title(r'$t=150 \mu s$')
ax44.axis('equal')
ax44.set_xticks(np.arange(0.0,X.max(), 0.1))

#plt.savefig('/home/sciortino/APEXD/obt_1/report_figs/RWfield_v%i.png'%version)
