import numpy as np
#We model a polystyrene bead of diameter 1 micrometer
# Base units: pg, μm, ms
# Energy unit: zJ (zeptojoules)
# Force unit: fN (femtonewtons)
zeta=6*np.pi*890*0.5 #drag coefficient in pg/ms, calculated using Stokes' law 
kB=0.0138 #zJ/K
T=295 #K
kBT=kB*T 
def F_drag(v):
    return -zeta*v
def F_harm(x):
    k=100 #spring constant of harmonic potential in zJ/(micrometer)^2
    return - k*x
def F_r(dt): #dt = timestep of simulation
    z = np.random.standard_normal()
    return z*(2*kBT*zeta/dt)**0.5

dt=0.000000001  #dt = timestep of simulation in ms
v_0=10 #initial velocity, in micrometer/ms= mm/s
xlist=[0] #list of positions in micrometers with initial position x_0=0
vlist=[v_0] # list of velocities 
tlist=[0] #list of time steps
m=0.55 #mass of particle in pg
steps=1000000

F_r1=F_r(dt) #initial random kick
for i in range(0,steps):
    F_r2=F_r(dt)
    x_new= xlist[-1]+vlist[-1]*dt + dt**2*(F_drag(vlist[-1])+F_harm(xlist[-1])+F_r1)/(2*m)
    v_new=(vlist[-1]+dt*(F_r1+F_r2 +F_drag(vlist[-1])+F_harm(x_new)+F_harm(xlist[-1]))/(2*m))/(1+dt*zeta/(2*m))
    xlist.append(x_new)
    vlist.append(v_new)
    tlist.append(tlist[-1]+dt)
    F_r1=F_r2

import matplotlib.pyplot as plt
fig, axs = plt.subplots(2, 1, figsize=(12,10))
axs[0].plot(tlist, xlist)
axs[0].set_title("Position wrt Time")
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('Position (μm)')
#axs[0].axvline(x=zeta/F_harm(-1), color='magenta', linestyle='--', linewidth=2, label=r'$\tau_x$ (Position Relaxation Time)')
axs[1].plot(tlist, vlist)
axs[1].set_title("Velocity wrt Time")
axs[1].set_xlabel('Time (ms)')
axs[1].set_ylabel('Velocity (mm/s)')
axs[1].axvline(x=m/zeta, color='red', linestyle='--', linewidth=2, label=r'$\tau_v$ (Velocity Relaxation Time)')
plt.legend()
plt.tight_layout()
plt.show()


                                                  
                                                  
                                                  

