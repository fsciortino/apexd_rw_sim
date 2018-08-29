from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
from matplotlib.pyplot import cm
import subprocess
import scipy
from helper_functions import *
import math
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
    from scipy.stats import maxwell  
    from scipy.signal import savgol_filter
    import scipy.fftpack

# emcee package is used to create interruptible pool for parallelization
from emcee.interruptible_pool import InterruptiblePool

# APEX-D Particle Dynamics Simulation library
import apds

me=9.10938356e-31
qe=1.6e-19
np.random.seed(433) # fix for reproducibility

plot_distr=False

# ==================== Simulation settings ========================= #
num_particles=8

# Spatial distribution in (x,y=0,z). Coil normally of radius 0.1.
# NB: Cannot run particles on the current loop itself! 
rad= 0.001
norm_mu=np.asarray([0.2,0.0])
cov=np.asarray([[rad,0.0],[0.0,rad]])

# Parameters for Gaussian distribution for parallel energy:
Kpar_mean=5.0;
#Tpar=1.0; Tperp=1.0;    # NEPOMUC Open port beam
Tpar  = 0.015; Tperp = 0.020;    # Buffer-gas trap expected output

# Time step:
dt = 1.0e-12

# Choose type of run: 0: 1us, 1: 20us, 2: 100us, 3: 1ms, 4: 10ms
run_type = 2

If=1.0e3 #50.0e3      # Current \propto magnetic field intensity
Rc=0.1         # radius of coil
rwa= 4.5      # RW field amplitude
rwf=0.75e5 #1.0e4 #1.0e5    # RW starting frequency
rwfg=10.0e9   # RW frequency gradient
rwd=1.0      # RW direction: +1 -> p+ drift direction, -1 -> e- drift direction

# RW E-field option: (0) Use E=0 always, everywhere;
#                    (1) Use ideal azimuthal E-field;
#                    (2) Read E-field file from  "./simion_Efields/ver.#". 
#                    (3) Locally-rotating E-field vortex (just a test)
iexe = 2

Efile=1 # Choose which electric field file to load (only used if iexe!=2)
tE1=1#8.0e-6 # Turn off rotating wall E-field at tE1 [s]
mm=1 # mode number for azimuthal field case

# B-field option: (0): Biot-Savart calculation (standard)
#                 (1): set Bz=0.1 everywhere
#                 (2): use external B field file (not well tested)
iexb = 0

# ============= Commands interpretation =============== #
colors=cm.rainbow(np.linspace(0,1,num_particles))

ttt=np.random.multivariate_normal(norm_mu,cov,num_particles)
pmx=ttt[:,0]; pmz=ttt[:,1]
pmy=[0.0,]*len(pmx)

Kpar = np.random.normal(Kpar_mean,Tpar,num_particles)
Kperp = np.random.exponential(scale=Tperp,size=num_particles)
#maxwell.rvs(loc=0.0,scale=Tperp,size=num_particles)
pma = 2*np.pi*np.random.random(num_particles)

# Set Kpar, Kperp and pma arbitrarily (for tests)
#Kpar = np.linspace(0.0,5, num_particles) #[2.9,]*len(pmx)
#Kperp = [1.7,]*len(pmx)
#pma=[0.5,]*num_particles #np.linspace(0.0,2*np.pi,num_particles)
              
# Get pitch angle for each particle:
theta0 = []
for comb in zip(Kperp,Kpar):
    vperp,vpar=comb
    theta0.append(math.atan2(vperp,vpar))

print "Kpar, Kperp, pma and theta_0 for each particle:"
for i in range(num_particles):
    print " ---> ", Kpar[i],Kperp[i],pma[i],theta0[i]
    
# Depending on "run_type", set different lengths and detail of simulation
if run_type==0:
    iend=iend = 1.0e6
    downsampling = 1.0e1 
elif run_type==1:
    iend = 2.0e7
    downsampling = 1.0e3
elif run_type==2:
    iend = 1.0e8
    downsampling = 1.0e4
elif run_type==3:
    iend = 1.0e9
    downsampling = 1.0e5
elif run_type==4:
    iend = 1.0e10
    downsampling = 1.0e6

print " Evolving until dt*iend = ", round(dt*iend,8), "s"
time_unit= 1.0e6 # to get us
time_label='time $[\mu s]$'

# ====================== Plot space and energy distributions ====================== #

if plot_distr:
    # Gaussian in space
    distr=np.random.multivariate_normal(norm_mu,cov,1000)
    
    plt.figure()
    plt.scatter(distr[:,0],distr[:,1],c='k',label='Gaussian')
    plt.scatter(ttt[:,0],ttt[:,1],c='r',label='sample')
    plt.xlabel('x [m]',fontsize=14)
    plt.ylabel('z [m]',fontsize=14)
    plt.axis('equal')

    # Parallel energy
    plt.figure()
    plt.hist(np.random.normal(Kpar_mean,Tpar,size=10000),bins=1000)
    axes=plt.gca()
    for en in range(len(Kpar)):
        par_en=Kpar[en]
        plt.plot([par_en,par_en],axes.get_ylim(),c=colors[en],linestyle='--')
    plt.xlabel(r'v$_{||} [eV]$')
    plt.ylabel('Counts [a.u.]')

    # Perpendicular energy
    plt.figure()
    plt.hist(maxwell.rvs(loc=0.0,scale=Tperp,size=10000),bins=1000)
    axes=plt.gca()
    for en in range(len(Kperp)):
        perp_en=Kperp[en]
        plt.plot([perp_en,perp_en],axes.get_ylim(),c=colors[en],linestyle='--')
    plt.plot([perp_en,perp_en],axes.get_ylim(),'r--')
    plt.xlabel(r'v$_{\perp} [eV]$')
    plt.ylabel('Counts [a.u.]')


# =================================
pid =  os.getpid() #multiprocessing.current_process().pid

# Run in parallel
def f_par(idx):
    #print "Running particle # ", idx+1
        
    res = apds.bbr_single_sim(pmx[idx],pmy[idx],pmz[idx],Kpar[idx],Kperp[idx],pma[idx],
             rwa,rwf,rwfg,rwd,dt,iend,downsampling,tE1,If,Rc,pid+idx,iexe,mm,Efile,iexb)
        
    return (pid,res)

n_cpu = multiprocessing.cpu_count()
pool = InterruptiblePool(processes=n_cpu if num_particles>n_cpu else num_particles)
out = pool.map(f_par,range(num_particles))

# Collect outputs in 'ress'. In the current version, there is none, but future updates may include some. 
pids = [out[i][0] for i in range(num_particles)]
ress = [out[i][1] for i in range(num_particles)]
    
    
########## First adiabatic moment -- mu ###########

mu_vars_all=OrderedDict()

for nop in range(1,num_particles+1):
    mu_vars=OrderedDict()
    mu_vars['time']=[]; mu_vars['mu']=[]; mu_vars['psi']=[]
    
    with open('./orb/mu_{:05}bbr.txt'.format(pids[nop-1]+nop-1),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(mu_vars)): 
            mu_vars[mu_vars.keys()[i]].append(float(elems[i]))

    for i in range(len(mu_vars)):
        mu_vars[mu_vars.keys()[i]] = np.asarray(mu_vars[mu_vars.keys()[i]])
        
    mu_vars['time']=np.asarray(mu_vars['time'])*time_unit
    mu_vars_all['%s'%(str(nop))] = mu_vars
    
# plot gyro-averaged mu
plot_gyro_averaged_mu(mu_vars_all)

############ Second adiabatic moment -- J ###########

J_vars_all=OrderedDict()

for nop in range(1,num_particles+1):

    J_vars=OrderedDict()
    J_vars['time']=[]; J_vars['J']=[]; J_vars['psi']=[]
    
    with open('./orb/jp_{:05}bbr.txt'.format(pids[nop-1]+nop-1),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(J_vars)): 
            J_vars[J_vars.keys()[i]].append(float(elems[i]))

    for i in range(len(J_vars)):         
        J_vars[J_vars.keys()[i]] = np.asarray(J_vars[J_vars.keys()[i]])
            
    J_vars['time']=np.asarray(J_vars['time'])*time_unit
    J_vars_all['%s'%(str(nop))] = J_vars
    
# plot J and psi
plot_bounce_averaged_J(J_vars_all)

# Plot Psi_gc, mu and J next to each other:
try:
    plot_adiabatic_conservation(mu_vars_all, J_vars_all)
except:
    pass

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

    with open('./orb/orb_{:05}bbr.txt'.format(pids[nop-1]+nop-1),'rb') as f:
        content=f.readlines()

    for line in content:
        es = line.split()
        elems=[e.replace('D','E') for e in es]

        for i in range(len(vars)): 
            vars[vars.keys()[i]].append(float(elems[i]))
            #print float(elems[i])
        
    vars['time']=np.asarray(vars['time'])*time_unit 

    # change into an ndarray
    for i in range(len(vars)): 
        vars[str(vars.keys()[i])]=np.asarray(vars[str(vars.keys()[i])])

    psi_gc_list.append(vars['Psi_g'])
    time_list.append(vars['time'])

    # save particle orbit info
    vars_all['%s'%(str(nop))] = vars



############
# Get Poincare' plots for a specific particle:
plot_poincare_tor(vars_all['1'])

# plot 3D orbits:
plot_3D_orbits(vars_all)

# ... and 2D orbit projection on x-y plane:
plot_2D_orbits(vars_all)
               
# Plot radius of gyrocenter vs. time:
plot_rgc(vars_all)

# Check phases of E-field wrt \Psi_gc    
plot_Efield_path(vars) # latest particle

# Compare \Psi_gc with drift-rotation frequency
plot_Psigc_rot(vars_all,rwa,rwf,rwfg)
                
#####  Summary plots  #####
plot_summary_1(vars_all,dt,downsampling,time_unit)
plot_summary_2(vars_all, rwa, rwf, rwfg)

stop

#### Test gyro-center approximation ####
# Check gyrocenter approximation metrics, evaluated at the particle's position:
get_gyrocenter_metrics(vars)

# Check gyrocenter approximation metrics, averaged over gyro-motion:
get_gyrocenter_averaged_metrics(vars)

    
####### Estimate radial transport  #####
hist_transport=False
if hist_transport:
    get_hist_transport(vars, t_transport,x_transport)



#################### tests

Emag=np.sqrt(vars['Ex']**2+vars['Ey']**2)
rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
dt_s=np.diff(vars['time'])[0]/time_unit
rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

# Compute E_{theta} and E_r
E_r = vars['Ex']*np.cos(rot_angle) + vars['Ey']*np.sin(rot_angle)
E_theta = - vars['Ex']*np.sin(rot_angle) + vars['Ey']*np.cos(rot_angle)
vExB_r=E_theta*vars['Bz']
vExB_tor= - E_r * vars['Bz']
    
# Get average rotation rate 
hist,bin_edges = np.histogram(vars['nrot'],int(max(vars['nrot']))+1)

# digitize
inds = np.digitize(vars['nrot'], bin_edges)

Etheta_means = [E_theta[inds == i].mean() for i in range(min(inds), max(inds))]
vExB_r_means = [vExB_r[inds == i].mean() for i in range(min(inds), max(inds))]

# find index of when the last complete turn begun
endturn_idx = np.argmin(np.abs(vars['nrot'] - bin_edges[-1]))
nrot = vars['nrot'][:endturn_idx]
rot_time = vars['time'][:endturn_idx]

# exclude last (incomplete) rotation
hist_s=hist[1:-1]; bin_edges_s=bin_edges[1:-1]

time_per_rotation = hist_s* dt * downsampling

# average "nrot" per rotation should constan
rate = 1.0/time_per_rotation

# time at the start of each rotation:
cum_rot_time = np.cumsum(time_per_rotation)*time_unit

cum_rot_time0 = np.concatenate(([0,],cum_rot_time))

plt.figure()
plt.plot(cum_rot_time, rate/1.0e3,'*-')
plt.xlabel(r'$time [\mu s]$')
plt.ylabel(r'$\langle f_{tor} \rangle$ [kHz]')

plt.figure();
plt.plot(cum_rot_time0, Etheta_means[:-1], 's-')
plt.xlabel(r'$time [\mu s]$')
plt.ylabel(r'$\langle E_{\theta} \rangle$ [V/m]')

plt.figure();
plt.plot(cum_rot_time0, vExB_r_means[:-1], 's-')
plt.xlabel(r'$time [\mu s]$')
plt.ylabel(r'$\langle v_{E\times B, r} \rangle$ [m/s]')


# Second method: low-pass filter of ftor
rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
dt_s=np.diff(vars['time'])[0]/time_unit
rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

window_length = 101 #must be odd
f_ave = savgol_filter(rot_freq,window_length,1)

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
start_idx= int( len(vars['x'])*4/5.0) #4500
end_idx= len(vars['x'])-1
plot_3D_orbits(vars, start_idx, end_idx)

        



