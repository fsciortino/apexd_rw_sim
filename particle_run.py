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
from helper_functions import *

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 14

import os
import multiprocessing
import pdb
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from scipy import optimize
    from scipy.stats import maxwell  
    from scipy.signal import savgol_filter
    from scipy.interpolate import interp1d
    import scipy.fftpack
    from scipy.signal import argrelextrema

# Run in virtualenv 'aenv' to use emcee
from emcee.interruptible_pool import InterruptiblePool
me=9.10938356e-31
qe=1.6e-19
np.random.seed(433) # fix for reproducibility

run_f90=False
plot_distr=False
parallel_run=True
if run_f90:
    command='./bbr.exe'
    subprocess.call(command)

# =============================================
# Initial parameters randomization
num_particles=1
colors=cm.rainbow(np.linspace(0,1,num_particles))

# distribution in (x,z), with y=0
# coil of radius 0.1. Cannot run particles on the current loop!
rad=1e-10 #0.001
norm_mu=np.asarray([0.2,0.0])
cov=np.asarray([[rad,0.0],[0.0,rad]])
ttt=np.random.multivariate_normal(norm_mu,cov,num_particles)
pmx=ttt[:,0]; pmz=ttt[:,1]
pmy=[0.0,]*len(pmx)

# Gaussian distribution for parallel energy
#Tfl=5.0; Tpar=1.0e-5; Tperp=1.0e-5;
Tfl=5.0; Tpar=1.0; Tperp=1.0;
pmk1 = np.random.normal(Tfl,Tpar,num_particles)
#pmk2 = np.random.exponential(Tperp,num_particles)
pmk2 = maxwell.rvs(loc=0.0,scale=Tperp,size=num_particles)
pma = 2*np.pi*np.random.random(num_particles)

# temporary
#pmk1=[5.0,]
#pmk2=[4.0e-5,]
#pmx=[0.2,]
#pmz=[0.0,]

# set iend/imbk < 10,000
dt = 1.0e-12
# Very short, very detailed path:
#iend = 1.0e6
#imbk = 1.0e2
#imbk2 = 1.0e2
# Short, detailed path:
#iend = 5.0e6 # 1.0e7 gives 10us
#imbk = 1.0e3
#imbk2 = 1.0e3
# Medium-long, not-so-detailed path:
iend = 1.0e9 #1.0e7#1.0e9   #15000000
imbk = 5e4#1.0e3 #5 #1.0e7
imbk2 = 5e4#1.0e3 #5 #1.0e7

print "Evolving until dt*iend = ", round(dt*iend,8), "s"
time_unit= 1.0e6 # to get us
time_unit2= 1.0e6
time_label='time $[\mu s]$'
# Turn off rotating wall E-field at tE1:
tE1=4.0 #40.0e-6 #0.0e6# 5000e-6 #4.0e-6

# ================================================

if plot_distr:
    # Gaussian
    distr=np.random.multivariate_normal(norm_mu,cov,1000)
    
    plt.figure()
    plt.scatter(distr[:,0],distr[:,1],c='k',label='Gaussian')
    plt.scatter(ttt[:,0],ttt[:,1],c='r',label='sample')
    plt.xlabel('x [m]',fontsize=14)
    plt.ylabel('z [m]',fontsize=14)
    plt.axis('equal')

    # Parallel energy
    plt.figure()
    plt.hist(np.random.normal(Tfl,Tpar,size=10000),bins=1000)
    axes=plt.gca()
    for en in range(len(pmk1)):
        par_en=pmk1[en]
        plt.plot([par_en,par_en],axes.get_ylim(),c=colors[en],linestyle='--')
    plt.title('Parallel energy distribution and sampling')

    # Perpendicular energy
    plt.figure()
    plt.hist(maxwell.rvs(loc=0.0,scale=Tperp,size=10000),bins=1000)
    axes=plt.gca()
    for en in range(len(pmk2)):
        perp_en=pmk2[en]
        plt.plot([perp_en,perp_en],axes.get_ylim(),c=colors[en],linestyle='--')
    plt.plot([perp_en,perp_en],axes.get_ylim(),'r--')
    plt.title('Perpendicular energy distribution and sampling')
    
# ==========================
If=5.0e2 #3
Rc=0.1
pmrwa=[0.0,]*len(pmx)
pmrwf=[1.0e5,]*len(pmx) #1.0e5
pmrwfg=[5.0e10,]*len(pmx) #5.0e9
pmrwd=[1.0,]*len(pmx)
fpm1=[0.0,]*len(pmx)
fpm2=[0.0,]*len(pmx)
fpm3=[0.0,]*len(pmx)
fpm4=[0.0,]*len(pmx)
fpm5=[0.0,]*len(pmx)
fpm6=[0.0,]*len(pmx)

# Choose which electric field file to load:
Efield_file=1

# =================================

if parallel_run:
    # Run in parallel
    def f_par(idx):
        print "Running particle # ", idx+1
        pid =  os.getpid() #multiprocessing.current_process().pid
        
        res = bbr_particle(pmx[idx],pmy[idx],pmz[idx],pmk1[idx],pmk2[idx],pma[idx],pmrwa[idx],pmrwf[idx],pmrwfg[idx],pmrwd[idx],fpm1[idx],fpm2[idx],fpm3[idx],fpm4[idx],fpm5[idx],fpm6[idx],dt,iend,imbk,imbk2,tE1,If,Rc,pid, Efield_file)
        
        return (pid,res)

    n_cpu = multiprocessing.cpu_count()
    pool = InterruptiblePool(processes=n_cpu if num_particles>n_cpu else num_particles)
    out = pool.map(f_par,range(num_particles))

    pids = [out[i][0] for i in range(num_particles)]
    ress = [out[i][1] for i in range(num_particles)]
    
else:
    # Run serially:
    res = bbr_particle(pmx,pmy,pmz,pmk1,pmk2,pma,pmrwa,pmrwf,pmrwfg,pmrwd,fpm1,fpm2,fpm3,fpm4,fpm5,fpm6,dt,iend,imbk,imbk2,tE1,If,Rc,os.getpid(),Efield_file)
    
    stt,stx,sty,stz,strr,stgx,stgy,stgz,stgr,stK,stKpr,stKpp,stps1,stps2,stmu,strt = res

    pids = np.asarray([os.getpid()]*num_particles)
    
########## First adiabatic moment -- mu ###########

mu_vars_all=OrderedDict()

for nop in range(1,num_particles+1):
    mu_vars=OrderedDict()
    mu_vars['time']=[]; mu_vars['mu']=[]; mu_vars['psi']=[]
    
    with open('./orb/mu{:04}_{:05}bbr.txt'.format(nop if parallel_run==False else 1,pids[nop-1]),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(mu_vars)): 
            mu_vars[mu_vars.keys()[i]].append(float(elems[i]))
    
    mu_vars['time']=np.asarray(mu_vars['time'])*time_unit2
    mu_vars_all['%s'%(str(nop))] = mu_vars
    
# plot gyro-averaged mu
plot_gyro_averaged_mu(mu_vars_all)

############ Second adiabatic moment -- J ###########

J_vars_all=OrderedDict()

for nop in range(1,num_particles+1):

    J_vars=OrderedDict()
    J_vars['time']=[]; J_vars['J']=[]; J_vars['psi']=[]
    
    with open('./orb/jp{:04}_{:05}bbr.txt'.format(nop if parallel_run==False else 1, pids[nop-1]),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(J_vars)): 
            J_vars[J_vars.keys()[i]].append(float(elems[i]))
    
    J_vars['time']=np.asarray(J_vars['time'])*time_unit2
    J_vars_all['%s'%(str(nop))] = J_vars
    
# plot J and psi
plot_bounce_averaged_J(J_vars_all)


############# orbits ####################

psi_gc_list = []
time_list = []

# Collect info for all particles
vars_all=OrderedDict()

# loop over particles
for nop in range(1,num_particles+1):
    vars=OrderedDict()
    vars['x']=[];vars['y']=[];vars['z']=[]
    vars['r']=[];vars['vx']=[];vars['vy']=[]
    vars['vz']=[];vars['time']=[];vars['Kttl']=[]
    vars['Kpara']=[];vars['Kperp']=[];vars['Bx']=[]
    vars['By']=[];vars['Bz']=[];vars['B']=[]
    vars['Ex']=[];vars['Ey']=[];vars['Ez']=[]
    vars['Psi']=[];vars['xgc']=[];vars['ygc']=[]
    vars['zgc']=[];vars['rgc']=[];vars['Psi_g']=[]
    vars['mu']=[];vars['nrot']=[]

    with open('./orb/orb{:04}_{:05}bbr.txt'.format(nop if parallel_run==False else 1, pids[nop-1]),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(vars)): 
            vars[vars.keys()[i]].append(float(elems[i]))
        
    vars['time']=np.asarray(vars['time'])*time_unit 

    # change into an ndarray
    for i in range(len(vars)): 
        vars[str(vars.keys()[i])]=np.asarray(vars[str(vars.keys()[i])])

    psi_gc_list.append(vars['Psi_g'])
    time_list.append(vars['time'])

    # save particle orbit info
    vars_all['%s'%(str(nop))] = vars



############
t_transport=[]; x_transport=[]; vx_transport=[]
 
mod_y=vars['y']; mod_y[0]=mod_y[1] # prevent change of sign detection

# find sign of y
asign=np.sign(mod_y)

# detect whether it changed sign wrt its neighboring points
signchange=((np.roll(asign,1)-asign) !=0).astype(bool)

# choose quadrant where x>0
pos_x=vars['x']>0

# save values
t_transport.append(vars['time'][signchange*pos_x])
x_transport.append(vars['x'][signchange*pos_x])
vx_transport.append(vars['x'][signchange*pos_x])

t_transport=np.asarray(t_transport)
x_transport=np.asarray(x_transport)
vx_transport=np.asarray(vx_transport)

plt.figure()
plt.scatter(x_transport,vx_transport)
#plt.scatter(vars['x'],vars['vx'])
plt.xlabel('x [m]')
plt.ylabel('v_x [m]')
plt.title('Poincare map')



stoppp












    
# plot 3D orbits:
plot_3D_orbits(vars_all)

# Plot radius of gyrocenter vs. time:
plot_rgc(vars_all)

    


#####  Summary plots  #####
summary_plots=True
if summary_plots:
    plot_summary_1(vars_all)
    plot_summary_2(vars_all, pmrwa, pmrwf, pmrwfg)

stop

#### Test gyro-center approximation ####

gyrocenter_approx_metrics=False
if gyrocenter_approx_metrics:
    # Check gyrocenter approximation metrics, evaluated at the particle's position:
    get_gyrocenter_metrics(vars)

    # Check gyrocenter approximation metrics, averaged over gyro-motion:
    get_gyrocenter_averaged_metrics(vars)

    
####### Estimate radial transport  #####
hist_transport=False
if hist_transport:
    get_hist_transport(vars, t_transport,x_transport)







# =======================================
# check if all particles have equal time vectors:
max_len=max([len(time_list[i]) for i in range(num_particles)])
max_len_idx=np.argmax([len(time_list[i]) for i in range(num_particles)])
psi_gc_all = np.ma.empty((num_particles,max_len))
psi_gc_all.mask=True
for pp in range(num_particles):
    psi_gc_all[pp,:len(psi_gc_list[pp])] = psi_gc_list[pp]

psi_gc_avg = np.median(psi_gc_all.data,axis=0)

time_list=np.asarray(time_list)

#plt.figure(15)
#plt.plot(time_list[max_len_idx],psi_gc_avg*1e4,c='k',linewidth=0.1,marker='*',label='Particle median')

# Subtract mean
#if np.mod(len(psi_gc_avg)/100,2)!=0:
#    if len(psi_gc_avg)/100>1:
#        window_length=len(psi_gc_avg)/100
#    else:
#        window_length=3
#else:
#    window_length=len(psi_gc_avg)/100+3
#psi_smoothed= savgol_filter(np.ndarray.tolist(psi_gc_avg),window_length,1)
#plt.plot(time_list[max_len_idx],psi_smoothed*1e4,'b-')
#psi_oscillation=psi_gc_avg-psi_smoothed

psi_oscillation=scipy.signal.detrend(psi_gc_avg, type='constant')

#plt.figure()
#plt.plot(time_list[max_len_idx],psi_oscillation,'g-')

yf=scipy.fftpack.fft(psi_oscillation)
N=len(psi_oscillation)
DeltaT=np.diff(time_list[max_len_idx]/time_unit)[0]
xf=np.linspace(0.0,1.0/(2.0*DeltaT),N/2.0)

fig,ax=plt.subplots()
Sf=2.0/N * np.abs(yf[:N//2])
ax.plot(xf,Sf)
ax.set_xlabel('f [Hz]',fontsize=14)
ax.set_ylabel('S(f)',fontsize=14)

#Y=np.fft.fft(psi_oscillation)
#freq=np.fft.fftfreq(len(psi_oscillation),DeltaT)
#fig,axxxx=plt.subplots()
#axxxx.plot(freq,np.abs(Y))

freq_osc=xf[np.argmax(Sf)]
print "Frequency of radial oscillaton: ", freq_osc/1.0e3, " kHz"

### Find peak frequency over different section of the time trace
num_segments=10

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

psi_chunks=np.asarray(list(split(psi_oscillation,num_segments)))
psigc_chunks=np.asarray(list(split(psi_gc_avg,num_segments)))
time_chunks=np.asarray(list(split(time_list[max_len_idx],num_segments)))

fig,(ax1,ax2,ax3)=plt.subplots(nrows=1,ncols=3,figsize=(18,8))
#ax2.plot(time_list[max_len_idx],psi_gc_avg*1e4,c='k',linewidth=0.1,marker='*',label='Particle median')

cols3=cm.rainbow(np.linspace(0,1,num_segments))
for ii in range(num_segments):
    psii=psi_chunks[ii]; timei=time_chunks[ii]
    psigci=psigc_chunks[ii]
    maxt_idx=len(timei)

    yf=scipy.fftpack.fft(psii)
    N=len(psii)
    DeltaT=np.diff(timei/time_unit)[0]
    xf=np.linspace(0.0,1.0/(2.0*DeltaT),N/2.0)

    Sf=2.0/N * np.abs(yf[:N//2])
    ax1.plot(xf,Sf,c=cols3[ii])
    ax2.plot(timei,psigci,c=cols3[ii]) #[ax2.get_ylim()[0],ax2.get_ylim()[1]],'k-')

ax1.set_xlabel('f [Hz]',fontsize=14)
ax1.set_ylabel('S(f)',fontsize=14)
ax1.set_xlim([0,5e5])
ax2.set_xlabel(r'%s'%time_label,fontsize=14)
ax2.set_ylabel(r'$\Psi_{gc}\times 10^4 [Tm^2]$',fontsize=14)

f_RW=pmrwf[0]+pmrwfg[0]*time_list[max_len_idx]*1.0e-6
ax3.plot(time_list[max_len_idx], f_RW/1.0e3, 'k-',linewidth=1.0)
ax3.set_xlabel(r'%s'%time_label,fontsize=14)
ax3.set_ylabel(r'$f_{RW}$ [kHz]',fontsize=14)

#######################################################
#
# ======================================================
#
#########################################################

# plot to compare to other RW conditions
fig7=plt.figure(7)
ax7=plt.subplot(111)
ax7.plot(time_list[max_len_idx],psi_gc_avg*1e4,linewidth=6,marker='*',
         label='Particle median,A={} V,f={} kHz'.format(pmrwa[0],pmrwf[0]/1e3))
ax7.set_xlabel(r'%s'%time_label,fontsize=14)
ax7.set_ylabel(r'$\Psi_{gc} \times 10^4 [Tm^{-2}]$',fontsize=14)
leg7=plt.legend()
leg7.draggable()

plt.figure(10)
ax10=plt.subplot(111)
survived_particles=np.array([np.sum(~psi_gc_all[:,i].mask) for i in range(psi_gc_all.shape[1])])
ax10.plot(time_list[max_len_idx], survived_particles,'b*')
ax10.set_xlabel(r'%s'%time_label)
ax10.set_ylabel('Number of surviving particles')




############## short orbit section ################
#ffig3 = plt.figure(figsize=(12,6))
#aax3 = plt.subplot(111, projection='3d')

start_idx= int( len(vars['x'])*4/5.0) #4500
end_idx= len(vars['x'])-1

plot_3D_orbits(vars, start_idx, end_idx)
# plot starting and ending points
#aax3.scatter(vars['x'][start_idx], vars['y'][start_idx], vars['z'][start_idx],c=colors[nop-1],marker='o',s=13)
#aax3.scatter(vars['x'][end_idx], vars['y'][end_idx], vars['z'][end_idx],c=colors[nop-1],marker='s',s=13)
    
# plot entire orbits
#aax3.plot(vars['x'][start_idx:end_idx], vars['y'][start_idx:end_idx], vars['z'][start_idx:end_idx],c=colors[nop-1],linestyle='-')
#aax3.plot(vars['xgc'][start_idx:end_idx], vars['ygc'][start_idx:end_idx], vars['zgc'][start_idx:end_idx],c=colors[nop-1],linestyle='--')
#aax3.set_xlabel('x [m]')
#aax3.set_ylabel('y [m]')
#aax3.set_zlabel('z [m]')
        



