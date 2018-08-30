'''
Helper functions for particle pushing simulations. Most of these functions are related to plotting of quantities that are obtained from simulations via bbr_particle.f90 and read via particle_run.py. 

F. Sciortino, August 2018
'''

import numpy as np
import pdb
import warnings
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['axes.labelsize'] = 18

m_p = 0.910938356e-30 # e/p mass
q_p=1.6e-19

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from scipy.signal import argrelextrema
    from scipy.stats import binned_statistic
    from scipy.signal import savgol_filter


##########################
def get_pol_drift(vars, time_unit=1.0e6):
    ''' Plot the polarization drift, using Eq. (92) of Fitzpatrick's notes:
    http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node17.html
    General formulation.
    '''
    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))
        
    # Polarization drift
    dt_s=np.diff(vars['time'])[0]/time_unit
    omega_c = vars['B']*qe/me

    # ExB drift cartesian components
    vExB_x = (vars['Ey']*vars['Bz'] - vars['Ez']*vars['By'])/(vars['B']**2)
    vExB_y = (-vars['Ex']*vars['Bz'] + vars['Ez']*vars['Bx'])/(vars['B']**2)
    vExB_z = (vars['Ex']*vars['By'] - vars['Ey']*vars['Bx'])/(vars['B']**2)

    # time derivatives of ExB drifts along cartesian directions
    dvExB_x_dt = np.gradient(vExB_x,dt_s)
    dvExB_y_dt = np.gradient(vExB_y,dt_s)
    dvExB_z_dt = np.gradient(vExB_z,dt_s)

    # polarization drift
    vp_x = (1.0/(omega_c * vars['B'])) * (vars['By']*dvExB_z_dt - vars['Bz']*dvExB_y_dt)
    vp_y = (1.0/(omega_c * vars['B'])) * (-vars['Bx']*dvExB_z_dt + vars['Bz']*dvExB_x_dt)
    vp_z = (1.0/(omega_c * vars['B'])) * (vars['Bx']*dvExB_y_dt - vars['By']*dvExB_x_dt)

    plt.figure()
    plt.plot(vars['time'],vp_x/1.0e3,c='r',label=r'$v_{p,x}$')
    plt.plot(vars['time'],vp_y/1.0e3,c='b',label=r'$v_{p,y}$')
    plt.plot(vars['time'],vp_z/1.0e3,c='g',label=r'$v_{p,z}$')
    plt.xlabel(r'%s'%time_label)
    plt.ylabel(r'$v_p$ [km/s]')
    plt.legend(fontsize=20, loc='best').draggable()

    
###########################
def plot_gyro_averaged_mu(mu_vars_all, time_unit=1.0e6):

    num_particles=len(mu_vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))
        
    fig0, ax0 = plt.subplots()
    
    for nop in range(1,num_particles+1):
        mu_vars = mu_vars_all['%s'%(str(nop))]
        
        ax0.plot(mu_vars['time'], mu_vars['mu'],c=colors[nop-1],marker='o')
        ax0.set_xlabel(r'%s'%time_label)
        ax0.set_ylabel(r'$\mu [m^2/s^2 T]$')

        
#########################
def plot_bounce_averaged_J(J_vars_all,time_unit=1.0e6):
    ''' Plot bounce-averaged 2nd adiabatic invariant.
    Every value of J corresponds to an integration between time when 
    a particle crosses the z=0 axis. 
    '''

    num_particles=len(J_vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    fig1,host = plt.subplots(1,figsize=(10,7))
    ax2=host.twinx()

    for nop in range(1,num_particles+1):
      
        J_vars = J_vars_all['%s'%(str(nop))]

        p1, =host.plot(J_vars['time'],J_vars['J']/1.0e4,c=colors[nop-1],marker='o')
        p2, =ax2.plot(J_vars['time'],J_vars['psi'],'k--', linewidth=2.0) #color=colors[nop-1],linestyle='--') 

    host.set_xlabel(r'%s'%time_label)
    host.set_ylabel(r'$J [10^4 m^2/s]$')
    ax2.set_ylabel(r'$\Psi$')

    host.yaxis.label.set_color(p1.get_color())
    ax2.yaxis.label.set_color(p2.get_color())
    tkw=dict(size=4,width=1.5)
    host.tick_params(axis='y',colors=p1.get_color(),**tkw)
    ax2.tick_params(axis='y',colors=p2.get_color(),**tkw)
    host.tick_params(axis='x',**tkw)


########################
def plot_adiabatic_conservation(mu_vars_all, J_vars_all, time_unit=1.0e6):
    ''' Plot mu, J and \psi_gc over time of simulation to assess conservation
    of adiabatic invariants. 
    '''
    num_particles=len(mu_vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))
        
    figg=plt.figure(figsize=(20,7))
    
    for nop in range(1,num_particles+1):
        mu_vars = mu_vars_all['%s'%(str(nop))]
        J_vars = J_vars_all['%s'%(str(nop))]

        mu_min = np.min(mu_vars['mu'])
        mu_mean = np.mean(mu_vars['mu'])
        mu_max = np.max(mu_vars['mu'])
        
        ax0=plt.subplot(131)
        ax0.plot(J_vars['time'], J_vars['psi']/1.0e4,c=colors[nop-1])
        ax0.set_xlabel(r'%s'%time_label)
        ax0.set_ylabel(r'$\Psi_{gc} [10^4 Tm^2]$')
        ax0.grid()
        
        ax1=plt.subplot(132)
        ax1.plot(mu_vars['time'], mu_vars['mu'],c=colors[nop-1])
        ax1.set_xlabel(r'%s'%time_label)
        ax1.set_ylabel(r'$\mu [m^2/s^2 T]$')
        ax1.grid()
        ax1.set_ylim([mu_mean-3*(mu_mean-mu_min), mu_mean+3*(mu_max-mu_mean)])
        
        ax2=plt.subplot(133)
        ax2.plot(J_vars['time'], J_vars['J']/1.0e4,c=colors[nop-1])
        ax2.set_xlabel(r'%s'%time_label)
        ax2.set_ylabel(r'$J [10^4 Tm^2]$')
        ax2.grid()
        

    
########################

def plot_rgc(vars_all,time_unit=1.0e6):
    ''' Plot radial position of gyrocenter vs. time. 
    Two figures are produced, with radial position expressed as distance from 
    the center of the simulation domain and in terms of the flux label. 

    '''
    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))
    
    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    fig4, ax4 = plt.subplots()
    for nop in range(1,num_particles+1):
        vars=vars_all['%s'%(str(nop))]
        
        # Analyze guiding centers:
        ax4.plot(vars['time'],vars['rgc'],c=colors[nop-1],marker='o')
        ax4.set_xlabel(r'%s'%time_label)
        ax4.set_ylabel(r'$r_{gc}$ [m]')

    fig5, ax5 = plt.subplots()
    for nop in range(1,num_particles+1):
        vars=vars_all['%s'%(str(nop))]
        
        # plot flux surface of guiding centers
        ax5.plot(vars['time'],vars['Psi_g']*1e4,c=colors[nop-1],marker='o')
        ax5.set_xlabel(r'%s'%time_label)
        ax5.set_ylabel(r'$\Psi_{gc} \times 10^4 [Tm^{-2}]$')
        
##################

def plot_Psigc_rot(vars_all,rwa,rwf,rwfg,time_unit=1.0e6):
    ''' Plot radial position of gyrocenter (in terms of the flux function of the gyrocenter)
     and the toroidal rotation frequency vs. time. 

    '''
    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))
    
    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    figg=plt.figure(figsize=(16,8))
    
    for nop in range(1,num_particles+1):
        vars = vars_all['%s'%(str(nop))]
        
        Emag=np.sqrt(vars['Ex']**2+vars['Ey']**2)
        rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
        dt_s=np.diff(vars['time'])[0]/time_unit
        rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

        # Compute E_{theta} and E_r
        E_r = vars['Ex']*np.cos(rot_angle) + vars['Ey']*np.sin(rot_angle)
        E_theta = - vars['Ex']*np.sin(rot_angle) + vars['Ey']*np.cos(rot_angle)
        vExB_r=E_theta*vars['Bz']
        vExB_tor= - E_r * vars['Bz']

        #vpar=2.0*vars['Kpara']*qe/me
        #vperp=vpar=2.0*vars['Kperp']*qe/me

        # Summary plot (style #2):
        ax001=plt.subplot(121)
        ax001.plot(vars['time'],vars['Psi_g']*1.0e4,c=colors[nop-1])
        ax001.set_xlabel(r'%s'%time_label)
        ax001.set_ylabel(r'$\Psi_{gc}\times 10^4 [Tm^2]$')
        ax001.grid()

        f_RW=rwf+rwfg*vars['time']/time_unit
        ax002 = plt.subplot(122, sharex=ax001)
        ax002.set_xlabel(r'%s'%time_label)
        ax002.set_ylabel(r'$f_{tor} [kHz]$')
        ax002.plot(vars['time'],rot_freq/1.0e3,c=colors[nop-1])
        if rwa!=0:
            ax002.plot(vars['time'],f_RW/1.0e3,c='k',linewidth=4.0,label='$f_{RW}$')
        if nop==1:
            leg=plt.legend(loc='best', fontsize=16).draggable()
        ax002.grid()


####################
def plot_2D_orbits(vars_all, start_idx=0, end_idx=0):
    ''' Plot 2D particle orbits in the x-y plane, given data collected 
    from one or multiple simulations. 
    '''
    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    fig4 = plt.figure(figsize=(8,6))
    ax4 = plt.subplot(111)
    
    for nop in range(1,num_particles+1):
        vars = vars_all['%s'%(str(nop))]

        # if requested end is beyond length of simulation, plot until the end
        if end_idx==0 or end_idx>len(vars['x'])-1:
            eid = len(vars['x'])-1

        if start_idx>len(vars['x'])-1:
            raise ValueError('Cannot plot 3D trajectory! Requested start index is longer than particle history.')
        
        # if an end-index is not given, plot the entire trajectory
        if end_idx==0: end_idx=len(vars['x'])-1
    
        # plot starting and ending points
        ax4.scatter(vars['x'][start_idx],vars['y'][start_idx],c=colors[nop-1],marker='o',s=13)
        ax4.scatter(vars['x'][eid],vars['y'][eid],c=colors[nop-1],marker='s',s=13)
    
        # plot requested section of the orbits
        ax4.plot(vars['x'][start_idx:eid],vars['y'][start_idx:eid],c=colors[nop-1],linestyle='-')
        #ax4.plot(vars['xgc'][start_idx:eid],vars['ygc'][start_idx:eid],c=colors[nop-1],linestyle='--')
        
        ax4.set_xlabel('x [m]')
        ax4.set_ylabel('y [m]')

        ax4.axis('equal')





###################
def plot_3D_orbits(vars_all, start_idx=0, end_idx=0):
    ''' Plot 3D particle orbits, given data collected from one or multiple 
    simulations. 
    '''
    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    fig3 = plt.figure(figsize=(12,6))
    ax3 = plt.subplot(111, projection='3d')
    
    for nop in range(1,num_particles+1):
        vars = vars_all['%s'%(str(nop))]

        # if requested end is beyond length of simulation, plot until the end
        if end_idx==0 or end_idx>len(vars['x'])-1:
            eid = len(vars['x'])-1
        else:
            eid = end_idx
            
        if start_idx>len(vars['x'])-1:
            raise ValueError('Cannot plot 3D trajectory! Requested start index is longer than particle history.')
    
        # plot starting and ending points
        ax3.scatter(vars['x'][start_idx],vars['y'][start_idx],vars['z'][start_idx],
                    c=colors[nop-1],marker='o',s=13)
        ax3.scatter(vars['x'][eid],vars['y'][eid],vars['z'][eid],
                    c=colors[nop-1],marker='s',s=13)
    
        # plot requested section of the orbits
        ax3.plot(vars['x'][start_idx:eid],vars['y'][start_idx:eid],vars['z'][start_idx:eid],
                 c=colors[nop-1],linestyle='-')
        ax3.plot(vars['xgc'][start_idx:eid],vars['ygc'][start_idx:eid],vars['zgc'][start_idx:eid],
                 c=colors[nop-1],linestyle='--')
        
        ax3.set_xlabel('x [m]')
        ax3.set_ylabel('y [m]')
        ax3.set_zlabel('z [m]')


        
#####################
def plot_poincare_tor(vars):
    '''
    Obtain Poincare' maps for the particle motion. Here, we select positions at which 
    the particle crosses the z=0 and the y=0 axis (giving 2 periodic motions).
    '''
    ###### y=0 crossing:
    t=[]; x=[]; vx=[]
 
    mod_y=vars['y']; mod_y[0]=mod_y[1] # prevent change of sign detection at start from 0

    # find sign of y
    asign=np.sign(mod_y)

    # detect whether it changed sign wrt its neighboring points
    signchange=((np.roll(asign,1)-asign) !=0).astype(bool)

    # choose quadrant where x>0
    pos_x=vars['x']>0

    # save values
    t.append(vars['time'][signchange*pos_x])
    x.append(vars['x'][signchange*pos_x])
    vx.append(vars['vx'][signchange*pos_x])

    t=np.asarray(t)
    x=np.asarray(x)
    vx=np.asarray(vx)

    # Poincare' plot:
    plt.figure()
    plt.scatter(x,vx/1.0e3)
    plt.xlabel('x [m]')
    plt.ylabel(r'$v_x$ [km/s]')

    ###### z=0 crossing:
    t=[]; x=[]; vx=[]
    mod_z=vars['z']; mod_z[0]=mod_z[1] # prevent change of sign detection

    # find sign of z
    asign=np.sign(mod_z)

    # detect whether it changed sign wrt its neighboring points
    signchange=((np.roll(asign,1)-asign) !=0).astype(bool)

    # choose quadrant where z>0
    pos_z=vars['z']>0

    # save values
    t.append(vars['time'][signchange*pos_z])
    x.append(vars['x'][signchange*pos_z])
    vx.append(vars['vx'][signchange*pos_z])

    t=np.asarray(t); x=np.asarray(x); vx=np.asarray(vx)

    # Poincare' plot:
    plt.figure()
    plt.scatter(x,vx/1.0e3)
    plt.xlabel('x [m]')
    plt.ylabel(r'$v_x$ [km/s]')
  


#######################################    
def plot_summary_1(vars_all, dt, downsampling, time_unit=1.0e6):
    ''' Particle simulation overview plots'''

    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles))

    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    figg=plt.figure(figsize=(18,10))
    
    for nop in range(1,num_particles+1):
        vars = vars_all['%s'%(str(nop))]

        # eliminate last 2 points to avoid visualization issues for 'dead' particles
        for i in range(len(vars)): 
            vars[vars.keys()[i]]=vars[vars.keys()[i]][:-2]
        
        Emag=np.sqrt(vars['Ex']**2+vars['Ey']**2)
        rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
        dt_s=np.diff(vars['time'])[0]/time_unit
        rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

        # Compute E_{theta} and E_r
        E_r = vars['Ex']*np.cos(rot_angle) + vars['Ey']*np.sin(rot_angle)
        E_theta = - vars['Ex']*np.sin(rot_angle) + vars['Ey']*np.cos(rot_angle)
        vExB_r=E_theta*vars['Bz']
        vExB_tor= - E_r * vars['Bz']

        # Obtain frequency of complete azimuthal turns from nrot
        hist,bin_edges = np.histogram(vars['nrot'],int(max(vars['nrot']))+1)

        # find index of when the last complete turn begun
        endturn_idx = np.argmin(np.abs(vars['nrot'] - bin_edges[-1]))
        nrot = vars['nrot'][:endturn_idx]
        rot_time = vars['time'][:endturn_idx]

        # exclude last (incomplete) rotation
        hist_s=hist[1:-1]; bin_edges_s=bin_edges[1:-1]
        
        time_per_rotation = hist_s* dt * downsampling
        f_tor_rot = 1.0/time_per_rotation

        # time at the start of each rotation:
        cum_rot_time = np.cumsum(time_per_rotation)*time_unit

        # Summary plots
        ax001=plt.subplot(331)
        ax001.plot(vars['time'],vars['Kpara'],c=colors[nop-1])
        ax001.set_ylabel(r'$K_{||} [eV]$')
        ax001.grid()

        ax002=plt.subplot(332,sharex=ax001)
        ax002.plot(vars['time'],vars['Kperp'],c=colors[nop-1])
        ax002.set_ylabel(r'$K_{\perp} [eV]$')
        ax002.grid()

        ax003=plt.subplot(333,sharex=ax001)
        ax003.plot(vars['time'],vars['Psi_g']*1.0e4,c=colors[nop-1])
        ax003.set_ylabel(r'$\Psi_{gc}\times 10^4 [Tm^2]$')
        ax003.grid()

        ax004=plt.subplot(334,sharex=ax001)
        ax004.plot(vars['time'],Emag,c=colors[nop-1])
        ax004.set_ylabel(r'|E| [Vm]')
        ax004.grid()

        ax005=plt.subplot(335,sharex=ax001)
        ax005.plot(vars['time'],vars['B']*1.0e3,c=colors[nop-1])
        ax005.set_ylabel(r'$B_{tot} [mT]$')
        ax005.grid()

        ax006=plt.subplot(336,sharex=ax001)
        ax006.plot(vars['time'],vars['rgc'],c=colors[nop-1])
        ax006.set_ylabel(r'$r_{gc}$ [m]')
        ax006.grid()

        ax007=plt.subplot(337,sharex=ax001)
        ax007.plot(vars['time'],vars['mu'],c=colors[nop-1])
        ax007.set_xlabel(r'%s'%time_label)
        ax007.set_ylabel(r'$\mu [Am^2]$')
        ax007.grid()

        ax008=plt.subplot(338,sharex=ax001)
        ax008.plot(cum_rot_time,f_tor_rot/1.0e3,c=colors[nop-1],marker='*')
        ax008.set_xlabel(r'%s'%time_label)
        ax008.set_ylabel(r'$\langle f_{tor}\rangle [kHz]$')
        ax008.grid()

        ax009 = plt.subplot(339,sharex=ax001)
        ax009.set_xlabel(r'%s'%time_label)
        ax009.set_ylabel(r'$v_{E\times B,tor} [m/s^{-1}]$')
        ax009.plot(vars['time'],vExB_tor,c=colors[nop-1])
        ax009.grid()

    
########################
def plot_summary_2(vars_all,rwa,rwf,rwfg, time_unit=1.0e6):
    ''' Particle simulation overview plots. 
    Inputs: vars (dictionary containing many quantities at the particle position,
            rwa: RW field amplitude
            rwf: RW starting frequency
            rwfg: RW frequency gradient
    '''

    num_particles=len(vars_all)
    colors=cm.rainbow(np.linspace(0,1,num_particles)) 

    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    figg=plt.figure(figsize=(18,10))

    for nop in range(1,num_particles+1):
        vars = vars_all['%s'%(str(nop))]

        # eliminate last 2 points to avoid visualization issues for 'dead' particles
        for i in range(len(vars)): 
            vars[vars.keys()[i]]=vars[vars.keys()[i]][:-2]
            
        Emag=np.sqrt(vars['Ex']**2+vars['Ey']**2)
        rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
        dt_s=np.diff(vars['time'])[0]/time_unit
        rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

        # low-pass filter of rotation frequency to average over bouncing & gyration
        window_length =int(len(vars['time'])/32.0)

        if np.mod(window_length,2)==0:
            window_length+=1
            
        #window_length = 101 #must be odd
        
        f_ave = savgol_filter(rot_freq,window_length,1)

        # Compute E_{theta} and E_r
        E_r = vars['Ex']*np.cos(rot_angle) + vars['Ey']*np.sin(rot_angle)
        E_theta = - vars['Ex']*np.sin(rot_angle) + vars['Ey']*np.cos(rot_angle)
        vExB_r=E_theta*vars['Bz']
        vExB_tor= - E_r * vars['Bz']

        #vpar=2.0*vars['Kpara']*qe/me
        #vperp=vpar=2.0*vars['Kperp']*qe/me

        # Summary plot (style #2):
        ax001=plt.subplot(421)
        ax001.plot(vars['time'],vars['Psi_g']*1.0e4,c=colors[nop-1])
        ax001.set_ylabel(r'$\Psi_{gc}\times 10^4 [Tm^2]$')
        ax001.grid()

        ax002=plt.subplot(422, sharex=ax001)
        ax002.plot(vars['time'],vars['mu'],c=colors[nop-1])
        ax002.set_ylabel(r'$\mu [Am^2]$')
        ax002.grid()
        
        ax003=plt.subplot(423, sharex=ax001)
        ax003.plot(vars['time'],vars['Kpara'], c=colors[nop-1])
        ax003.set_ylabel(r'$K_{||} [eV]$')
        ax003.grid()

        ax004=plt.subplot(424, sharex=ax001)
        ax004.plot(vars['time'],vars['Kperp'], c=colors[nop-1]) 
        ax004.set_ylabel(r'$K_{\perp} [eV]$')
        ax004.grid()

        ax005=plt.subplot(425, sharex=ax001)
        ax005.plot(vars['time'],vars['rgc'],c=colors[nop-1])
        ax005.set_ylabel(r'$r_{gc} [m]$')
        ax005.grid()

        ax006=plt.subplot(426, sharex=ax001)
        ax006.plot(vars['time'],vars['B']*1.0e3,c=colors[nop-1])
        ax006.set_ylabel(r'$B_{tot} [mT]$')
        ax006.grid()

        if num_particles==1:
            ccc='r'
        else:
            ccc=colors[nop-1]
        
        f_RW=rwf+rwfg*vars['time']/time_unit
        ax007 = plt.subplot(427, sharex=ax001)
        ax007.set_xlabel(r'%s'%time_label)
        ax007.set_ylabel(r'$f_{tor} [kHz]$')
        ax007.plot(vars['time'],rot_freq/1.0e3,c=colors[nop-1])
        ax007.plot(vars['time'],f_ave/1.0e3,c=ccc,linewidth=3.0)
        if rwa!=0:
            ax007.plot(vars['time'],f_RW/1.0e3,c='k',linewidth=4.0,label='$f_{RW}$')
        if nop==1 and rwa!=0:
            plt.legend(loc='best').draggable()
        ax007.grid()

        ax008 = plt.subplot(428, sharex=ax001)
        ax008.set_xlabel(r'%s'%time_label)
        ax008.set_ylabel(r'$v_{E\times B,r} [m/s^{-1}]$')
        ax008.plot(vars['time'],vExB_r,c=colors[nop-1])
        ax008.grid()



    
def plot_Efield_path(vars, time_unit=1.0e6):
    ''' Plot electric field components and total magnitude 
    at the particle position during a simulation 
    '''
    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))
    
    Emag=np.sqrt(vars['Ex']**2+vars['Ey']**2)
    rot_angle=np.arccos(vars['xgc']/vars['rgc']) # in radians
    dt_s=np.diff(vars['time'])[0]/time_unit
    rot_freq=np.abs(np.gradient(rot_angle,dt_s ))/(2*np.pi) # in Hz

    # Compute E_{theta} and E_r
    E_r = vars['Ex']*np.cos(rot_angle) + vars['Ey']*np.sin(rot_angle)
    E_theta = - vars['Ex']*np.sin(rot_angle) + vars['Ey']*np.cos(rot_angle)
    minE = np.min([np.min(E_r),np.min(E_theta)])
    maxE = np.max([np.max(E_r),np.max(E_theta)])
    
    vExB_r=E_theta*vars['Bz']
    vExB_tor= - E_r * vars['Bz']
    minvE = np.min([np.min(vExB_r),np.min(vExB_tor)])
    maxvE = np.max([np.max(vExB_r),np.max(vExB_tor)])

    #vpar=2.0*vars['Kpara']*qe/me
    #vperp=vpar=2.0*vars['Kperp']*qe/me

    colors=cm.rainbow(np.linspace(0,1,5))

    ##### Check role of E-field phases
    figgg=plt.figure(figsize=(10,10))
    ax0001=plt.subplot(511)
    ax0001.plot(vars['time'],vars['Psi_g']*1.0e4,c=colors[0])
    ax0001.set_ylabel(r'$\Psi_{gc}\times 10^4$ [T$m^2$]')
    ax0001.grid()
    
    ax0006=plt.subplot(512,sharex=ax0001)
    ax0006.plot(vars['time'],E_r,c=colors[1])
    ax0006.set_ylabel(r'$E_r$ [Vm]')
    ax0006.set_ylim([minE,maxE])
    ax0006.grid()
    
    ax0007=plt.subplot(513,sharex=ax0001)
    ax0007.plot(vars['time'],E_theta,c=colors[2])
    ax0007.set_ylabel(r'$E_{\theta}$ [Vm]')
    ax0007.set_ylim([minE,maxE])
    ax0007.grid()
    
    ax0008=plt.subplot(514,sharex=ax0001)
    ax0008.plot(vars['time'],vExB_r,c=colors[3])
    ax0008.set_ylabel(r'$v_{E\times B, r}$ [m/s]')
    ax0008.set_ylim([minvE,maxvE])
    ax0008.grid()
    
    ax0009=plt.subplot(515,sharex=ax0001)
    ax0009.plot(vars['time'],vExB_tor,c=colors[4])
    ax0009.set_xlabel(r'%s'%time_label)
    ax0009.set_ylabel(r'$v_{E\times B,\theta}$ [m/s]')
    ax0009.set_ylim([minvE,maxvE])
    ax0009.grid()

def plot_Bfield_path(vars,time_unit=1.0e6):
    ''' Plot magnetic field components and total magnitude 
    at the particle position during a simulation 
    '''
    if time_unit==1.0e6:
        time_label='time $[\mu s]$'
    else:
        time_label='time $[10^{%f} s]$'%(np.log10(1.0/time_unit))

    colors=cm.rainbow(np.linspace(0,1,4))
    
    figgg=plt.figure(figsize=(18,10))
    ax0006=plt.subplot(411)
    ax0006.plot(vars['time'],vars['Bx'],c=colors[0])
    ax0006.set_ylabel(r'$B_x$ [T]')
    ax0006.grid()
    
    ax0007=plt.subplot(412,sharex=ax0006)
    ax0007.plot(vars['time'],vars['By'],c=colors[1])
    ax0007.set_ylabel(r'$B_y$ [T]')
    ax0007.grid()
    
    ax0008=plt.subplot(413,sharex=ax0006)
    ax0008.plot(vars['time'],vars['Bz'],c=colors[2])
    ax0008.set_ylabel(r'$B_z$ [T]')
    ax0008.grid()
    
    ax0009=plt.subplot(414,sharex=ax0006)
    ax0009.plot(vars['time'],vars['B'],c=colors[3])
    ax0009.set_xlabel(r'%s'%time_label)
    ax0009.set_ylabel(r'$B$ [T]')
    ax0009.grid()



    
######################
def get_hist_transport(vars):
    ''' Produce histograms to show radial transport   '''

    t_transport=[]; x_transport=[]
    
    ####### detect passing through y=0  #####
    mod_y=vars['y']; mod_y[0]=mod_y[1] # prevent change of sign detection
    asign=np.sign(mod_y)
    signchange=((np.roll(asign,1)-asign) !=0).astype(bool)
    pos_x=vars['x']>0
    t_transport.append(vars['time'][signchange*pos_x])
    x_transport.append(vars['x'][signchange*pos_x])

    #####
    x_vals=np.asarray([]); t_vals=np.asarray([])
    time_max=10.0
    for particle in range(len(x_transport)):

        sel=tuple([t_transport[particle]<time_max])
        x_vals=np.concatenate((x_vals,x_transport[particle][sel]))
        t_vals=np.concatenate((t_vals,t_transport[particle][sel]))

    # get bins
    histt,bin_edges=np.histogram(t_vals,bins=6)

    # split x values in bins
    xx=[]
    for tbin_idx in range(len(bin_edges)-1):
        xx.append([x_vals[i] for i in range(x_vals.shape[0])
                   if t_vals[i]>bin_edges[tbin_idx] and t_vals[i]<bin_edges[tbin_idx+1]])

    plt.figure()
    plt.hist(xx[0],bins=50)
    for bi in range(len(bin_edges)-1):
        plt.hist(xx[bi],bins=10,label=r'Time bin: [{:.2f},{:.2f}] $\mu s$'.format(bin_edges[bi],bin_edges[bi+1]))
    plt.legend().draggable()

    plt.figure()
    bi=0
    plt.hist(xx[bi],bins=10,label=r'Time bin: [{:.2f},{:.2f}] $\mu s$'.format(bin_edges[bi],bin_edges[bi+1]))
    plt.legend().draggable()
    
    #plt.figure()
    #plt.plot(t_vals,x_vals,'bo')
    #plt.plot(vars['time'],vars['x'],'r')
        


#######################                       
def get_gyrocenter_metrics(vars):
    '''Test gyro-center approximation by plotting
    \rho/L_B = \rho \nabla B / B
    and 
    1/(\omega_c \tau_b) = (B/\omega_c) [\partial B/\partial t]^{-1}

    Here, all quantities are evaluated for each time step of the simulation.
    Gradients are evaluated between two adjacent positions of the particle.
    '''

    # gyrofrequency:
    w_gyro = q_p * vars['B'] / m_p
    
    plt.figure()
    plt.plot(vars['time'],w_gyro/1.0e9, 'b-')
    plt.xlabel(r'time $[\mu s]$')
    plt.ylabel(r'$\omega_c [GHz]$')
    
    # compute \rho/L_B and 1/(\omega_c \tau_B) for all gyro-orbits
    B_mag = np.zeros(len(vars['time']))
    rho_over_L_B = np.zeros(len(vars['time']))
    one_over_omegac_tau_B = np.zeros(len(vars['time']))

    for pp in range(1,len(vars['time'])):
        
        # get v-parallel (convert parallel energy in Joules first) 
        vpar = np.sqrt(2*vars['Kpara'][pp]*1.6e-19/m_p)
        
        rho = abs(vars['rgc'][pp]-vars['r'][pp])

        # Use magnitude of B
        #dB_dx = (vars['B'][pp] - vars['B'][pp-1])/(vars['x'][pp] - vars['x'][pp-1])
        #dB_dy = (vars['B'][pp] - vars['B'][pp-1])/(vars['y'][pp] - vars['y'][pp-1])
        #dB_dz = (vars['B'][pp] - vars['B'][pp-1])/(vars['z'][pp] - vars['z'][pp-1])

        # or better to use dB_dr?
        dB_dr = (vars['B'][pp] - vars['B'][pp-1])/(vars['r'][pp] - vars['r'][pp-1])
        
        #grad_B_max = max([dB_dx,dB_dy,dB_dz])
        grad_B_max = abs(dB_dr)
        
        # L_B:
        L_B = vars['B'][pp] / grad_B_max
        
        # tau_B:
        dB_dt = grad_B_max * vpar
        tau_B = (dB_dt/vars['B'][pp])**(-1)

        # metric
        rho_over_L_B[pp] = rho /L_B 
        one_over_omegac_tau_B[pp] =  1.0/(w_gyro[pp] * tau_B)


    # Plot gyro-center approximation metrics

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,7))#, sharey=True)
    ax1.plot(vars['time'][2:], rho_over_L_B[2:], 'ro-')
    ax1.set_xlabel(r'time $[\mu s]$')
    ax1.set_ylabel(r'$\rho/L_B$')
    ax2.plot(vars['time'][2:], one_over_omegac_tau_B[2:], 'bo-')
    ax2.set_xlabel(r'time $[\mu s]$')
    ax2.set_ylabel(r'$1/(\omega_c \tau_B)$')
    ax1.set_title('No gyro-averaging performed')

    # Repeat plot, but bin values and take median within those bins
    vals1, bin_edges1, bin_num1 = binned_statistic(vars['time'][2:], rho_over_L_B[2:], 'median', bins=300)
    bin_centers1 = bin_edges1[:-1] + np.diff(bin_edges1)[0]
    vals2, bin_edges2, bin_num2 = binned_statistic(vars['time'][2:], one_over_omegac_tau_B[2:], 'median', bins=300)
    bin_centers2 = bin_edges2[:-1] + np.diff(bin_edges2)[0]
    
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,7))#, sharey=True)
    ax1.plot(bin_centers1, vals1, 'ro-', label='Median values, no gyro-averaging')
    ax1.set_xlabel(r'time $[\mu s]$')
    ax1.set_ylabel(r'$\rho/L_B$')
    ax2.plot(bin_centers2, vals2, 'bo-', label='Median values, no gyro-averaging')
    ax2.set_xlabel(r'time $[\mu s]$')
    ax2.set_ylabel(r'$1/(\omega_c \tau_B)$')
    ax2.legend().draggable()
    #ax1.set_title('No gyro-averaging performed')
    #ax2.set_title('Median values within bins')
    
################################################
def get_gyrocenter_averaged_metrics(vars):
    '''Test gyro-center approximation by plotting
    \rho/L_B = \rho \nabla B / B
    and 
    1/(\omega_c \tau_b) = (B/\omega_c) [\partial B/\partial t]^{-1}

    Both quantities are evaluated as an average over a gyro-orbit.

    NOT COMPLETED. REQUIRES MORE ATTENTION. 
    '''
    
    # gyrofrequency:
    w_gyro = q_p * vars['B'] / m_p

    plt.figure()
    plt.plot(vars['time'],w_gyro/1.0e9, 'b-')
    plt.xlabel(r'time $[\mu s]$')
    plt.ylabel(r'$\omega_c [GHz]$')

    # times over which to average:
    tau_gyro=1.0/w_gyro
    
    # compute \rho/L_B and 1/(\omega_c \tau_B) for all gyro-orbits
    gyro_time_steps=[0,]
    rho_mean = []
    B_gyro_mean = []
    grad_B_gyro_max = []
    L_B = []
    rho_over_L_B = []
    tau_B = []
    dB_dt = []
    one_over_omegac_tau_B = []
    
    for pp in range(1,len(tau_gyro)):
        gyro_time_steps.append(gyro_time_steps[-1]+tau_gyro[pp])
        idxs = (vars['time']<gyro_time_steps[-1])*(vars['time']>=gyro_time_steps[-2])

        Bx_mean = np.mean(vars['Bx'][idxs])
        By_mean = np.mean(vars['By'][idxs])
        Bz_mean = np.mean(vars['Bz'][idxs])

        Bmag_min = np.min(vars['B'][idxs])
        Bmag_max = np.max(vars['B'][idxs])
        
        # get v-parallel (convert parallel energy in Joules first)
        Kpara_Joules = np.mean(vars['Kpara'][idxs]) * 1.6e-19
        vpar_mean = np.sqrt(2*Kpara_Joules/m_p)
        
        rho_mean.append(np.mean(abs(vars['rgc'][idxs]-vars['r'][idxs])))
    
        B_gyro_mean.append(np.sqrt(Bx_mean**2+By_mean**2+Bz_mean**2))
        grad_B_gyro_max.append(np.max(abs(Bmag_max-Bmag_min))/rho_mean[-1])
            
        #pdb.set_trace()
        L_B.append(B_gyro_mean[-1] / grad_B_gyro_max[-1])
        
        # Now tau_B:
        dB_dt.append(grad_B_gyro_max[-1] * vpar_mean)
        tau_B.append((dB_dt[-1]/B_gyro_mean[-1])**(-1))
        
        # gyro-center theory averaged metrics
        rho_over_L_B.append(rho_mean[-1]/L_B[-1])
        one_over_omegac_tau_B.append( 1.0/(w_gyro[pp] * tau_B[-1]))


    # transform into numpy arrays for convenience
    gyro_time_steps=np.asarray(gyro_time_steps)
    rho_mean=np.asarray(rho_mean)
    B_gyro_mean = np.asarray(B_gyro_mean)
    grad_B_gyro_max = np.asarray(grad_B_gyro_max)
    L_B = np.asarray(L_B)
    rho_over_L_B = np.asarray(rho_over_L_B)
    one_over_omegac_tau_B = np.asarray(one_over_omegac_tau_B)

    # Plot gyro-center approximation metrics

    #pdb.set_trace()
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,7))#, sharey=True)
    ax1.plot(gyro_time_steps[1:], rho_over_L_B, 'ro-', label='gyro-averaged values')
    ax1.set_xlabel(r'time $[\mu s]$')
    ax1.set_ylabel(r'$\rho/L_B$')
    ax2.plot(gyro_time_steps[1:], one_over_omegac_tau_B, 'bo-', label='gyro-averaged values')
    ax2.set_xlabel(r'time $[\mu s]$')
    ax2.set_ylabel(r'$1/(\omega_c \tau_B)$')
    ax1.legend().draggable()
