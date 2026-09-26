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
k=F_harm(-1)
dt=0.002 
m=0.55 #mass of particle in pg
numberofparticles = 10000
steps =80000 
velocities=[]
positions=[]
for j in range(1,numberofparticles+1):#iterating over each particle
    x=0#initializing each particle at rest at x=0
    v= np.random.normal(loc=0.0, scale=np.sqrt(kBT/m))
    F_r1=F_r(dt) #initial random kick
    for i in range(0,steps):
        F_r2=F_r(dt)
        x_new= x+v*dt + dt**2*(F_drag(v)+F_harm(x)+F_r1)/(2*m)
        v_new=(v+dt*(F_r1+F_r2 +F_drag(v)+F_harm(x_new)+F_harm(x))/(2*m))/(1+dt*zeta/(2*m))
        x=x_new
        v=v_new
        F_r1=F_r2
        if (i+1)==steps: 
            positions.append(x)
            velocities.append(v)
            
import matplotlib.pyplot as plt
def Pos_Boltzmann(x):
    return np.sqrt(k/(2*np.pi*kBT))*np.exp(-k*x**2/(2*kBT))
def Vel_Boltzmann(v):
    return np.sqrt(m/(2*np.pi*kBT))*np.exp(-m*v**2/(2*kBT))
tau_v=m/zeta
#sigma_sq=kBT*(2*tau_v/dt - (2*tau_v**2)*(1-np.exp(-dt/tau_v))/dt**2)/m
sigma_sq=(kBT/m)/(1+(dt/(2*tau_v)))
def V_avg_Boltz(v):
    return np.sqrt(1/(2*np.pi*sigma_sq))*np.exp(-v**2/(2*sigma_sq))
vlist = np.linspace(min(velocities), max(velocities), 300)
xlist=np.linspace(min(positions), max(positions), 300)
plt.figure()
plt.hist(velocities, bins=30, density=True, alpha=0.6, color='mediumpurple', edgecolor='black', label='Equilibrium Velocity Distribution')
plt.plot(vlist, Vel_Boltzmann(vlist), 'r-', linewidth=2, label=r'Theory: $P(v) \propto \exp\left(-\frac{m v^2}{2 k_B T}\right)$')      
plt.plot(vlist, V_avg_Boltz(vlist), 'b-', linewidth=2, label=r'$P_{v_\text{sim}}(v)$')
plt.xlabel('Velocity (mm/s)')  
plt.ylabel('Probability Density')
plt.legend()
plt.figure()
plt.hist(positions, bins=30, density=True, alpha=0.6, color='lightseagreen', edgecolor='black', label='Equilibrium Position Distribution')
plt.plot(xlist, Pos_Boltzmann(xlist), 'r-', linewidth=2, label=r'Theory: $P(x) \propto \exp\left(-\frac{k x^2}{2 k_B T}\right)$')
plt.xlabel(r'Position ($\mu\text{m}$)')  
plt.ylabel('Probability Density')
plt.legend()
plt.show()
                                                  
                                 
