import numpy as np
#We model a polystyrene bead of diameter 1 micrometer
# Base units: pg, μm, ms
# Energy unit: zJ (zeptojoules)
# Force unit: fN (femtonewtons)
eta=890 #viscosity in micro-Pa*s, 890 for water.
zeta=6*np.pi*eta*0.5 #drag coefficient in pg/ms, calculated using Stokes' law 
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

dt=0.00000005  #dt = timestep of simulation in ms 
v_0=10 #initial velocity, in micrometer/ms= mm/s
xlist=[50000] #list of positions in micrometers with initial position x_0=0
vlist=[v_0] # list of velocities 
tlist=[0] #list of time steps
m=0.55 #mass of particle in pg
steps=10000
numberofparticles=2000
# Initialize 1000 independent particles at equlibrium using equipartition theorem.
v_0list = np.random.normal(loc=0.0, scale=np.sqrt(kBT/m), size=numberofparticles)
x_0list = np.zeros(numberofparticles)
vv_0=[]
for k in range(0,steps):
    tlist.append(tlist[-1]+dt)
for j in range(0,len(x_0list)):
    xlist=[x_0list[j]] #list of positions 
    vlist=[v_0list[j]] # list of velocities 

    F_r1=F_r(dt) #initial random kick
    for i in range(0,steps):
        F_r2=F_r(dt)
        x_new= xlist[-1]+vlist[-1]*dt + dt**2*(F_drag(vlist[-1])+F_harm(xlist[-1])+F_r1)/(2*m)
        v_new=(vlist[-1]+dt*(F_r1+F_r2 +F_drag(vlist[-1])+F_harm(x_new)+F_harm(xlist[-1]))/(2*m))/(1+dt*zeta/(2*m))
        xlist.append(x_new)
        vlist.append(v_new)
        F_r1=F_r2
    vv_0.append(np.array(vlist)*v_0list[j])

VAC = np.sum(vv_0, axis=0) / numberofparticles

import matplotlib.pyplot as plt
plt.plot(tlist, VAC, label='1000 particle ensemble')
tau_v=m/zeta
plt.axvline(x=tau_v, color='red', linestyle='--', linewidth=2, label=r'$\tau_v$ (Velocity Relaxation Time)')
plt.suptitle("Velocity Autocorrelation for a Free Particle")
plt.xlabel('Time (ms)')
plt.ylabel(r'$\langle v(t)v(0)\rangle$ ($\mu\text{m}^2/\text{ms}^2$)')
plt.plot(tlist,(kBT/m)*np.exp(-np.array(tlist)/tau_v),label='Theoretical result:'+ r'$\frac{K_BT}{m}\exp\left(-\frac{|t|}{\tau_v}\right)$')
plt.legend()
plt.show()
                                                  
                                 
