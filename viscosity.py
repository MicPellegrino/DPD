import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def read_xvg(file_name):
    x = []
    y = []
    with open(file_name, 'r') as input_xvg:
        for line in input_xvg:
            cols = line.split()
            if not( cols[0][0] == '#' or cols[0][0] == '@' ):
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)

# Read input parameters for forcing amplitude and box size
A = 0.0
k = 0.0
Lz = 0.0
m = 0.0
with open('parameters.txt', 'r') as input_par:
    for line in input_par:
        cols = line.split()
        if cols[0] == 'amplitude_x':
            A = float(cols[2])
        if cols[0] == 'wave_number':
            k = float(cols[2])
        if cols[0] == 'Lz':
            Lz = float(cols[2])
        if cols[0] == 'mass':
            m = float(cols[2])

# Read binned velocity profile output
vx, z = read_xvg('velocity.xvg')
plt.plot(z, vx, 'rx', label='MD')

# Perform the best cosine fit
# fit_cosine = lambda t, a, v0 : a*np.cos(k*t) + v0
fit_cosine = lambda t, a : a*np.cos(k*t)
popt, pcov = opt.curve_fit(fit_cosine, z, vx)
zfit = np.linspace(0.0, Lz, 100)
vfit = fit_cosine(zfit, *popt)
plt.plot(zfit, vfit, 'k-', label='fit')

# Estimate viscosity
# ...

# Plot
plt.xlabel('z')
plt.ylabel('vx')
plt.legend()
plt.show()
