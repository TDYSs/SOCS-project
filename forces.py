import numpy as np
import matplotlib.pyplot as plt

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

x_grid = np.linspace(-5*1e-9/2, 50*1e-9/2, 500)

plt.plot(x_grid, binding_force(x_grid)+spring_force(x_grid))
plt.show()