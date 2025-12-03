import numpy as np
import matplotlib.pyplot as plt

R = 1e-9  # Radius of the Brownian particle [m].
eta = 1e-3  # Viscosity of the medium.
gamma = 6 * np.pi * R * eta  # Drag coefficient of the medium. 
rho = 2e3  # Density of the particle [kg/m^3]
m = 4 * np.pi / 3 * rho * R ** 3  # Mass of the particle [kg].

tau = m / gamma  # Momentum relaxation time.

dt = 1e-12  # Time step [s].
duration = 10e-6  # Total time [s].



x0 = 0  # Initial position [m].
v0 = 0  # Initial velocity [m/s].

def spring_force(x):
    k = 0.0002
    y = x-8e-9
    return -y*k

def binding_force(x):
    kBT = 4.11e-21
    sigma = 1e-9/2
    xpos = 16e-9
    x = x- xpos
    U = 16*kBT
    return -U*x/(sigma**2)*np.exp(-x**2/ (2*sigma**2))

#x_grid = np.linspace(-5*1e-9/2, 50*1e-9/2, 500)

#plt.plot(x_grid, binding_force(x_grid)+spring_force(x_grid))
#plt.show()

def evolution_viscous(x0, gamma, dt, duration):
    """
    Function to generate the solution for the Langevin equation with 
    inertia.
    
    Parameters
    ==========
    x0 : Initial position of the oscillator [m].
    gamma : Friction coefficient [N*s/m].    
    dt : Time step for the numerical solution [s].
    duration : Total time for which the solution is computed [s].
    """
    
    kBT = 4.11e-21  # kB*T at room temperature [J].
    
    D = kBT / gamma  # Diffusion constant [m^2 / s].
    
    # Coefficients for the finite difference solution.
    c_noise = np.sqrt(2 * D * dt)

    N = int(np.ceil(duration / dt))  # Number of time steps.

    x = np.zeros(N)
    rn = np.random.normal(0, 1, N - 1)
    
    x[0] = x0
    x_eq = 16e-9
    stable_time = 5e-9
    eps = 3e-10

    stable_steps_required = int(stable_time / dt)
    stable_counter = 0
    
    for i in range(N - 1):
        f = spring_force(x[i]) + binding_force(x[i])
        x[i + 1] = x[i] + c_noise * rn[i] + f*dt/gamma
        if i % N*0.01 == 0:
            print ('max time elapsed')
            print (i/N)
        if abs(x[i+1] - x_eq) < eps:
            stable_counter += 1
        else:
            stable_counter = 0
        
        # Stop when condition met
        if stable_counter >= stable_steps_required:
            print(f"Stopped early at t = {i*dt:.3e} s.")
            return x[:i+2], i
    return x, i
