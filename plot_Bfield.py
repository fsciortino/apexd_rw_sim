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

import warnings

with warnings.catch_warnings():
    import scipy.interpolate as interp

with open('./b.txt','rb') as f:
    content=f.readlines()
        
xi=[]; zi=[]; Psi=[]; Bs=[]; Bx=[]; By=[]; Bz=[];
for line in content:
    es = line.split()
    elems=[e.replace('D','E') for e in es]
            
    try:
        Xxi,Xzi,XPsi,XBs,XBx,XBy,XBz = elems
        # modify Fortran format
        xxi=float(Xxi.replace('D','E'))
        xzi=float(Xzi.replace('D','E'))
        xPsi=float(XPsi.replace('D','E'))
        xBs=float(XBs.replace('D','E'))
        xBx=float(XBx.replace('D','E'))
        xBy=float(XBy.replace('D','E'))
        xBz=float(XBz.replace('D','E'))
            
        # save
        xi.append(xxi); zi.append(xzi)
        Psi.append(xPsi); Bs.append(xBs)
        Bx.append(xBx); By.append(xBy); Bz.append(xBz)
                
    except:
        pass

xi=np.asarray(xi)
zi=np.asarray(zi)
Psi=np.asarray(Psi)
Bs=np.asarray(Bs)
Bx=np.asarray(Bx)
By=np.asarray(By)
Bz=np.asarray(Bz)

xi_u = np.unique(xi)
zi_u = np.unique(zi)

Psi_m=np.zeros((len(xi_u),len(zi_u)))
Bs_m=np.zeros((len(xi_u),len(zi_u)))
Bx_m=np.zeros((len(xi_u),len(zi_u)))
By_m=np.zeros((len(xi_u),len(zi_u)))
Bz_m=np.zeros((len(xi_u),len(zi_u)))

for i in range(len(xi)):
    for ii in range(len(xi_u)):
        for jj in range(len(zi_u)):
            if xi[i]==xi_u[ii] and zi[i]==zi_u[jj]:
                Psi_m[ii,jj]=Psi[i]
                Bs_m[ii,jj]=Bs[i]
                Bx_m[ii,jj]=Bx[i]
                By_m[ii,jj]=By[i]
                Bz_m[ii,jj]=Bz[i]


#fig1=plt.figure(1,figsize=(10,8))
#ax1 = plt.subplot(221)

#xmin=ax3.get_xlim()[0]; xmax=ax3.get_xlim()[1]
#ymin=ax3.get_ylim()[0]; ymax=ax3.get_ylim()[1]
#zmin=ax3.get_zlim()[0]; zmax=ax3.get_zlim()[1]

#mask_1=np.ones((len(xi_u),len(zi_u)))

#for i in range(len(xi_u)):
#    for j in range(len(zi_u)):
#        if xi_u[i]>0 and xi_u[i]<xmax and zi_u[j]>zmin and zi_u[j]<zmax:
#            mask_1[i,j]=0

#Bs_masked=np.ma.masked_array(Bs_m,mask=mask_1)
#X_masked=np.ma.masked_array(X,mask=mask_1)
#Z_masked=np.ma.masked_array(Z,mask=mask_1)

#ax3.contour(X_masked,Z_masked,Bs_masked,zdir='y', offset=ax3.get_ylim()[1])

# Create dense grid:
xmd=np.linspace(min(xi_u),max(xi_u),900)
zmd=np.linspace(min(zi_u),max(zi_u),1000)
Xd,Zd = np.meshgrid(xmd,zmd)

# Interpolate fields on dense grid
Bxmd = interp.interp2d(xi_u,zi_u,Bx_m,kind='cubic')(xmd,zmd)
Bymd = interp.interp2d(xi_u,zi_u,By_m,kind='cubic')(xmd,zmd)
Bzmd = interp.interp2d(xi_u,zi_u,Bz_m,kind='cubic')(xmd,zmd)
Bsd = interp.interp2d(xi_u,zi_u,Bs,kind='cubic')(xmd,zmd)

plt.figure()
#plt.contour(Z_masked,X_masked, Bs_masked)
plt.contourf(Xd,Zd, Bsd)
plt.xlabel('x [m]'); plt.ylabel('z [m]')
plt.colorbar()

# gyrofrequency
eom= 1.7588200236e11
w_c = eom*Bsd
plt.figure()
plt.contourf(Xd,Zd,np.log10(w_c))
cbar=plt.colorbar(label='log-gyrofrequency')
plt.xlabel('x [m]'); plt.ylabel('z [m]')
