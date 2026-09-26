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

m=0.55 #mass of particle in pg
steps=100000
numberofparticles=1000
dt=0.000001 #dt = timestep of simulation in ms
# Initialize 1000 independent particles at equlibrium using Boltzmann distribution.
v_0list = np.random.normal(loc=0.0, scale=np.sqrt(kBT/m), size=numberofparticles)
x_0list = np.random.normal(loc=0.0, scale=np.sqrt(kBT/F_harm(-1)), size=numberofparticles)
tlist = np.arange(steps + 1) * dt


displacementlist=[]
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
    displacementlist.append(np.array(xlist)-x_0list[j])

MSD = np.sum(np.array(displacementlist)**2, axis=0) / numberofparticles

import matplotlib.pyplot as plt

plt.plot(tlist, MSD, label='1000 particle ensemble')
plt.plot(tlist,2*kBT*np.array(tlist)/zeta,label='2D|t|')
plt.suptitle("Verifying Einstein's Relation")
plt.xlabel('Time (ms)')
plt.ylabel(r'$\langle (x(t)-x(0))^2\rangle$ ($\mu\text{m}^2$)')
plt.legend()
plt.show()
                                                  
                                 
