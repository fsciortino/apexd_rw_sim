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

import os
import multiprocessing
import pdb

import matplotlib as mpl
mpl.rcParams['xtick.labelsize']=14
mpl.rcParams['ytick.labelsize']=14

# Run in virtualenv 'aenv' to use emcee
from emcee.interruptible_pool import InterruptiblePool

np.random.seed(13) # fix for reproducibility

run_f90=False
plot_distr=False
parallel_run=True
if run_f90:
    command='./bbr.exe'
    subprocess.call(command)

# =============================================
# Initial parameters randomization
num_particles=8

# distribution in (x,z), with y=0
# coil of radius 0.1. Cannot run particles on the current loop!
rad=0.001
mu=np.asarray([0.2,0.0])
cov=np.asarray([[rad,0.0],[0.0,rad]])
ttt=np.random.multivariate_normal(mu,cov,num_particles)
pmx=ttt[:,0]; pmz=ttt[:,1]
pmy=[0.0,]*len(pmx)

# Gaussian distribution for parallel energy
Tfl=5.0; Tpar=1.0; Tperp=1.0;
pmk1 = np.random.normal(Tfl,Tpar,num_particles)
pmk2 = np.random.exponential(Tperp,num_particles)
pma = 2*np.pi*np.random.random(num_particles)

# set iend/imbk < 10,000
dt = 1.0e-12
iend = 15000000 
imbk = 20000
imbk2 = 20000
print "Evolving until dt*iend = ", dt*iend
time_unit= 1.0e6 # to get us
time_unit2= 1.0e6

# ================================================

if plot_distr:
    # Gaussian
    plt.figure(100)
    plt.scatter(ttt[:,0],ttt[:,1],label='Gaussian')
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')
    plt.axis('equal')
    
# ==========================

scan_amp=[4.0,]*6 #[3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0]
scan_freq=np.logspace(np.log10(9e4),np.log10(2.0e5),6)
# Turn off rotating wall E-field at tE1:
tE1=6.0e-6
    
fig7=plt.figure(7)
ax7=plt.subplot(111)
fig10=plt.figure(10)
ax10=plt.subplot(111)
    
for i in range(len(scan_amp)):
    pmrwa=[scan_amp[i],]*len(pmx)
    pmrwf=[scan_freq[i],]*len(pmx)
    pmrwd=[1.0,]*len(pmx)
    fpm1=[0.0,]*len(pmx)
    fpm2=[0.0,]*len(pmx)
    fpm3=[0.0,]*len(pmx)
    fpm4=[0.0,]*len(pmx)
    fpm5=[0.0,]*len(pmx)
    fpm6=[0.0,]*len(pmx)

    # =================================

    if parallel_run:
        # Run in parallel
        def f_par(idx):
            print "Running particle # ", idx+1
            pid =  os.getpid() #multiprocessing.current_process().pid
        
            res = bbr_particle(pmx[idx],pmy[idx],pmz[idx],pmk1[idx],pmk2[idx],pma[idx],pmrwa[idx],pmrwf[idx],pmrwd[idx],fpm1[idx],fpm2[idx],fpm3[idx],fpm4[idx],fpm5[idx],fpm6[idx],dt,iend,imbk,imbk2,tE1,pid)
        
            return (pid,res)

        n_cpu = multiprocessing.cpu_count()
        pool = InterruptiblePool(processes=n_cpu if num_particles>n_cpu else num_particles)
        out = pool.map(f_par,range(num_particles))

        pids = [out[i][0] for i in range(num_particles)]
        ress = [out[i][1] for i in range(num_particles)]
    
    else:
        # Run serially:
        res = bbr_particle(pmx,pmy,pmz,pmk1,pmk2,pma,pmrwa,pmrwf,pmrwd,fpm1,fpm2,fpm3,fpm4,fpm5,fpm6,dt,iend,imbk,imbk2,tE1,os.getpid())
    
        stt,stx,sty,stz,strr,stgx,stgy,stgz,stgr,stK,stKpr,stKpp,stps1,stps2,stmu,strt = res

        pids = np.asarray([os.getpid()]*num_particles)
    
    #######################
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    psi_gc_list = []
    time_list = []
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
        
        vars['time']=[float(tt)*time_unit for tt in vars['time']]  #imbk gives sampling frequency  

        # change into an ndarray
        for i in range(len(vars)): 
            vars[str(vars.keys()[i])]=np.asarray(vars[str(vars.keys()[i])])
            
        psi_gc_list.append(vars['Psi_g'])
        time_list.append(vars['time'])
 

    # particle statistics
    max_len=max([len(time_list[i]) for i in range(num_particles)])
    max_len_idx=np.argmax([len(time_list[i]) for i in range(num_particles)])
    psi_gc_all = np.ma.empty((num_particles,max_len))
    psi_gc_all.mask=True
    for pp in range(num_particles):
        psi_gc_all[pp,:len(psi_gc_list[pp])] = psi_gc_list[pp]

    psi_gc_avg = np.median(psi_gc_all,axis=0)

    # plot to compare to other RW conditions
    fig7=plt.figure(7)
    ax7.plot(time_list[max_len_idx],psi_gc_avg*1e4,linewidth=6,marker='*',label='Particle median,A={:.1f} V,f={:.1f} kHz'.format(pmrwa[0],pmrwf[0]/1e3))
    ax7.set_xlabel(r'time $[\mu s]$',fontsize=14)
    ax7.set_ylabel(r'$\Psi_{gc} \times 10^4 [Tm^{-2}]$',fontsize=14)
    ax7.plot([tE1, tE1],[ax7.get_ylim()[0],ax7.get_ylim()[1]],'k-')
    leg7=plt.legend(loc='best',fontsize=12)
    leg7.draggable()

    fig10=plt.figure(10)
    survived_particles=np.array([np.sum(~psi_gc_all[:,i].mask) for i in range(psi_gc_all.shape[1])])
    ax10.plot(time_list[max_len_idx], survived_particles,'*',label='Particle median,A={:.1f} V,f={:.1f} kHz'.format(pmrwa[0],pmrwf[0]/1e3))
    ax10.set_xlabel(r'time $[\mu s]$',fontsize=14)
    ax10.set_ylabel('Number of surviving particles',fontsize=14)


