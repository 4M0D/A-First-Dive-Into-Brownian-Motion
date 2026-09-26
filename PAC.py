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

dt=0.001  #dt = timestep of simulation in ms
v_0=10 #initial velocity, in micrometer/ms= mm/s
xlist=[0] #list of positions in micrometers with initial position x_0=0
vlist=[v_0] # list of velocities 
tlist=[0] #list of time steps
m=0.55 #mass of particle in pg
numberofparticles = 1000
steps = 600000
# Initialize 1000 independent particles at equlibrium using Boltzmann distribution.
v_0list = np.random.normal(loc=0.0, scale=np.sqrt(kBT/m), size=numberofparticles)
x_0list = np.random.normal(loc=0.0, scale=np.sqrt(kBT/F_harm(-1)), size=numberofparticles)
tlist = np.arange(steps + 1) * dt
xx_0sum=np.zeros(len(tlist))
for j in range(0,len(x_0list)):#iterating over each particle
    x=x_0list[j] #latest position and velocity of particle
    v=v_0list[j] 
    xx_0sum[0]+=x*x_0list[j] #adding intial contribution to correlation from the particle

    F_r1=F_r(dt) #initial random kick
    for i in range(0,steps):
        F_r2=F_r(dt)
        x_new= x+v*dt + dt**2*(F_drag(v)+F_harm(x)+F_r1)/(2*m)
        v_new=(v+dt*(F_r1+F_r2 +F_drag(v)+F_harm(x_new)+F_harm(x))/(2*m))/(1+dt*zeta/(2*m))
        x=x_new
        v=v_new
        F_r1=F_r2
        xx_0sum[i+1]+=x*x_0list[j] #adding contribution at (i+1)th time iteration

PAC = xx_0sum / numberofparticles

import matplotlib.pyplot as plt
plt.plot(tlist, PAC, label='1000 particle ensemble')
tau=zeta/F_harm(-1)
plt.axvline(x=tau, color='red', linestyle='--', linewidth=2, label=r'$\tau=\zeta/k$ (Trap Relaxation Time)')
plt.suptitle("Position Autocorrelation in an Harmonic Potential")
plt.xlabel('Time (ms)')
plt.ylabel(r'$\langle x(t)x(0)\rangle$ ($\mu\text{m}^2$)')
plt.plot(tlist,(kBT/F_harm(-1))*np.exp(-np.array(tlist)/tau),label='Theoretical result:'+ r'$\frac{K_BT}{k}\exp\left(-\frac{|t|}{\tau}\right)$')
plt.legend()
plt.show()
                                                  
                                 
