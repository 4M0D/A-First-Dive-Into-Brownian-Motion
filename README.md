## Preface

Cells, the units of life, are sophisticated microscopic machines made of organelles and macromolecules. These components of the cell undergo collisions with molecules of the cytoplasmic fluid (cytosol), leading to erratic motion. Motivated by their significance to us as living beings ourselves, this text is a stepping stone into such dynamics as an exploration of Brownian motion in one dimension. We start our discussion with how position and velocity correlations evolve with time for a Brownian particle and how the ‘Langevin equation’ helps us model its dynamics. We then briefly review some concepts from stochastic calculus and use them to derive the Fokker–Planck equation, a second-order linear PDE describing the evolution of the probability distribution of the position and velocity of such a particle. Finally, we test some results derived in the text via numerical simulation.


---

## Python scripts for numerical simulation discussed:

The Python scripts implement the numerical Langevin integration schemes detailed in **Section 4: Simulating Langevin Dynamics**:

* **`basic_brownian.py`** *(Sec 4.1)* — Implements the basic underdamped second-order Langevin integrator to output single-particle position $x(t)$ and velocity $v(t)$ trajectories.
* **`Einstein_Relation.py`** *(Sec 4.2)* — Verifies the Einstein diffusion relation ($S(t) \propto t$) at short timescales ($t \ll \tau$) inside a harmonic well.
* **`DAC.py`** *(Sec 4.2)* — Plots Displacement Autocorrelation (DAC) relaxation to demonstrate exponential saturation due to confinement inside an optical trap.
* **`PAC.py`** *(Sec 4.2)* — Plots Position Autocorrelation (PAC) decay validating the overdamped approximation for $t \gg \tau_v$.
* **`VAC.py`** *(Sec 4.2)* — Resolves Velocity Autocorrelation (VAC) exponential decay over inertial timescales ($\tau_v = m/\zeta$).
* **`Boltzmanns.py`** *(Sec 4.3)* — Demonstrates spatial and velocity equilibration to stationary Maxwell-Boltzmann distributions.
