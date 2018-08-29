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

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 18

version=1
path='/home/sciortino/APEXD/obt_1/simion_Efields/ver.%i/'%version
file=open(path+'/e3d_hd.txt','rb')
if version>4:
    num_electrodes=8
else:
    num_electrodes=4
    
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

tot_vals=nx*ny*nz

Ex=np.zeros((nx+1,ny+1,nz+1,num_electrodes))
Ey=np.zeros((nx+1,ny+1,nz+1,num_electrodes))
Ez=np.zeros((nx+1,ny+1,nz+1,num_electrodes))

for electrode in range(1,num_electrodes+1):
    with open(path+'/e3d%i.txt'%electrode,'rb') as f:
        content=f.readlines()

    for line in content:
        linevars = line.split() 
        idx=int(linevars[3])
        idy=int(linevars[4])
        idz=int(linevars[5])
        Ex[idx,idy,idz,electrode-1]+=float(linevars[0].replace('D','E'))
        Ey[idx,idy,idz,electrode-1]+=float(linevars[1].replace('D','E'))
        Ez[idx,idy,idz,electrode-1]+=float(linevars[2].replace('D','E'))

xx=[dx*i for i in range(nx+1)]
yy=[dy*i for i in range(ny+1)]
zz=[dz*i for i in range(nz+1)]
X,Y=np.meshgrid(xx,yy)


#### plot electric field from each electrode individually ####
fig1=plt.figure(1,figsize=(10,8))
if num_electrodes==4: sbp=221
elif num_electrodes==8: sbp=421

ax1 = plt.subplot(sbp)
ax1.quiver(X,Y,Ex[:,:,25,0],Ey[:,:,25,0], color='b')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.axis('equal')

for ii in range(num_electrodes-1):
    ax2 = plt.subplot(sbp+ii+1,sharex=ax1)
    ax2.quiver(X,Y,Ex[:,:,25,ii+1],Ey[:,:,25,ii+1], color='b')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.axis('equal')


#### Next, obtain combination field from multiple electrodes ####
amp=1
Fel=1.0e5
Phrw=1.0

time=[0.0,2.5e-6,5.0e-6,7.5e-6]

Vel=np.zeros(num_electrodes)
Exx=np.zeros((nx+1,ny+1,len(time)))
Eyy=np.zeros((nx+1,ny+1,len(time)))

for t in range(len(time)):
    tp=time[t]
    for iele in range(num_electrodes):
        Vel[iele] = amp*np.sin( 2.0*np.pi*Fel*tp + Phrw*iele*np.pi/(num_electrodes/2.0))
        
    for nnx in range(nx+1):
        for nny in range(ny+1):
            for iele in range(num_electrodes):
                # field must be cylindrically symmetric and equal along z
                Exx[nnx,nny,t] += Vel[iele]*Ex[nnx,nny,25,iele]
                Eyy[nnx,nny,t] += Vel[iele]*Ey[nnx,nny,25,iele]

# sub-sample
ss_factor = 3
skip = 3 # to avoid crossing arrows at edges
Exx=Exx[skip:-skip:ss_factor,skip:-skip:ss_factor,:]
Eyy=Eyy[skip:-skip:ss_factor,skip:-skip:ss_factor,:]
X=X[skip:-skip:ss_factor,skip:-skip:ss_factor]
Y=Y[skip:-skip:ss_factor,skip:-skip:ss_factor]
                           
# plot rotating fields
fig2=plt.figure(figsize=(10,12))


ax11 = plt.subplot(221)
ax11.quiver(X,Y,Exx[:,:,0],Eyy[:,:,0], color='r')
#ax11.set_xlabel('x [m]')
ax11.set_ylabel('y [m]')
ax11.set_title(r'$t=0 \mu s$')
ax11.axis('equal')
ax11.set_xticks(np.arange(0.0,X.max(), 0.1))

ax22 = plt.subplot(222,sharex=ax11)
ax22.quiver(X,Y,Exx[:,:,1],Eyy[:,:,1], color='r')
#ax22.set_xlabel('x [m]')
#ax22.set_ylabel('y [m]')
ax22.set_title(r'$t=50 \mu s$')
ax22.axis('equal')
ax22.set_xticks(np.arange(0.0,X.max(), 0.1))

ax33 = plt.subplot(223,sharex=ax11)
ax33.quiver(X,Y,Exx[:,:,2],Eyy[:,:,2], color='r')
ax33.set_xlabel('x [m]')
ax33.set_ylabel('y [m]')
ax33.set_title(r'$t=100 \mu s$')
ax33.axis('equal')
ax33.set_xticks(np.arange(0.0,X.max(), 0.1))

ax44 = plt.subplot(224,sharex=ax11)
ax44.quiver(X,Y,Exx[:,:,3],Eyy[:,:,3], color='r')
ax44.set_xlabel('x [m]')
#ax44.set_ylabel('y [m]')
ax44.set_title(r'$t=150 \mu s$')
ax44.axis('equal')
ax44.set_xticks(np.arange(0.0,X.max(), 0.1))

plt.savefig('/home/sciortino/APEXD/obt_1/report_figs/RWfield_v%i.png'%version)
