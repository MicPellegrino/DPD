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
A = 0.2
Lx = 7.5
Ly = 7.5
Lz = 7.5
m = 1.0
N = 300
with open('parameters.txt', 'r') as input_par:
    for line in input_par:
        cols = line.split()
        if cols[0] == 'amplitude_x':
            A = float(cols[2])
        if cols[0] == 'Lx':
            Lx = float(cols[2])
        if cols[0] == 'Ly':
            Ly = float(cols[2])
        if cols[0] == 'Lz':
            Lz = float(cols[2])
        if cols[0] == 'mass':
            m = float(cols[2])
        if cols[0] == 'n_part':
            N = int(cols[2])
k = 2.0*np.pi/Lz

# Read binned velocity profile output
vx, z = read_xvg('velocity.xvg')
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
ax1.plot(z, vx, 'rx', markersize=12.5, markeredgewidth=2.5, label='MD')

# Perform the best cosine fit
# fit_cosine = lambda t, a, v0 : a*np.cos(k*t) + v0
fit_cosine = lambda t, V : V*np.cos(k*t+0.5*np.pi)
popt, pcov = opt.curve_fit(fit_cosine, z, vx)
zfit = np.linspace(0.0, Lz, 100)
vfit = fit_cosine(zfit, *popt)
ax1.plot(zfit, vfit, 'k:', linewidth=3.0, label='LS fit')

# Estimate viscosity
V = popt[0]
eta = (A/V) * (N*m*Lz) / (4.0*Lx*Ly*np.pi**2)
eta_err_p = (A/V-A/(V+pcov[0])) * (N*m*Lz) / (4.0*Lx*Ly*np.pi**2)
eta_err_m = (A/V-A/(V-pcov[0])) * (N*m*Lz) / (4.0*Lx*Ly*np.pi**2)
print("Viscosity: eta = "+str(eta))
print("err + = "+str(eta_err_p[0]))
print("err - = "+str(eta_err_m[0]))

# Plot
ax1.set_xlabel(r'$z$', fontsize=25.0)
ax1.set_ylabel(r'$v_x$', fontsize=25.0)
ax1.legend(fontsize=20.0)
ax1.xaxis.set_tick_params(labelsize=17.5)
ax1.yaxis.set_tick_params(labelsize=17.5)
ax1.set_xlim([0.0, Lz])

# Full picture

visc = dict()
std_p = dict()
std_m = dict()

a = [0.5, 0.4, 0.3, 0.2, 0.1]

visc[1.0] = [0.17209059465994817, 0.2066474834433838, 0.20072765410549262, 0.19275525943669286, 0.15165505986023037]
std_p[1.0] = [0.008482977921073217, 0.012825457117872937, 0.004018844250183898, 0.0021207423680370315, 0.004852264855495643]
std_m[1.0] = [0.009410759799668779, 0.014643086601155457, 0.004186482549941916, 0.002168458225392576, 0.005183992735668083]

visc[0.8] = [0.17365916497164896, 0.15142649680353293, 0.16843736403247203, 0.14602431546171366, 0.24742711527004302]
std_p[0.8] = [0.0008043300781139267, 0.0013972038078165243, 0.005682120558331204, 0.0006496325891672714, 0.006749785970966174]
std_m[0.8] = [0.0008118505079464486, 0.0014234724082411486, 0.006093222005097702, 0.0006554646479585324, 0.007139304964802186]

visc[0.6] = [0.1381726119903816, 0.14206248406194294, 0.16070579519272907, 0.12936618911267728, 0.16326194132733654]
std_p[0.6] = [0.007018040831694271, 0.0006508368267659457, 0.0024421222308707215, 0.002062028752136271, 0.0011281357713311846]
std_m[0.6] = [0.007811568683416718, 0.0006568553851111004, 0.0025186708335837088, 0.002129928627681886, 0.0011439450322132116]

visc[0.4] = [0.1187629232670612, 0.12365187108592596, 0.12697874355575314, 0.11923469144582394, 0.08179111216125565]
std_p[0.4] = [0.003013190852452409, 0.0017062401677317878, 0.0003182601762714305, 0.002584821645821012, 0.0016310914824502626]
std_m[0.4] = [0.003174262302534058, 0.0017546644901075274, 0.0003198635917762567, 0.002701970641921434, 0.001698848917865413]

visc[0.2] = [0.10569829450475558, 0.09373019607315802, 0.12699278848363185, 0.1139155373959572, 0.08648786763427933, ]
std_p[0.2] = [0.0012045213627230734, 0.0013401399626724745, 0.0009316307897298681, 0.0008221882533405261, 0.002143338998831068]
std_m[0.2] = [0.0012326147352330478, 0.0013795903036620538, 0.0009455033895845582, 0.0008342304112478213, 0.002255111154057006]

ax2.errorbar(a, visc[1.0], yerr=[std_p[1.0],std_m[1.0]], fmt='mD', elinewidth=3.0, capsize=7.0, \
        capthick=3.0, ms=12.5, label=r'$f=1.0$')
ax2.errorbar(a, visc[0.8], yerr=[std_p[0.8],std_m[0.8]], fmt='ro', elinewidth=3.0, capsize=7.0, \
        capthick=3.0, ms=12.5, label=r'$f=0.8$')
ax2.errorbar(a, visc[0.6], yerr=[std_p[0.6],std_m[0.6]], fmt='bx', elinewidth=3.0, capsize=7.0, \
        capthick=3.0, ms=12.5, mew=3.0, label=r'$f=0.6$')
ax2.errorbar(a, visc[0.4], yerr=[std_p[0.4],std_m[0.4]], fmt='gs', elinewidth=3.0, capsize=7.0, \
        capthick=3.0, ms=12.5, label=r'$f=0.4$')
ax2.errorbar(a, visc[0.2], yerr=[std_p[0.2],std_m[0.2]], fmt='kH', elinewidth=3.0, capsize=7.0, \
        capthick=3.0, ms=12.5, label=r'$f=0.2$')
ax2.legend(fontsize=17.5)
ax2.set_xlabel(r'$\mathcal{A}$', fontsize=25.0)
ax2.set_ylabel(r'$\eta$', fontsize=25.0)
ax2.xaxis.set_tick_params(labelsize=17.5)
ax2.yaxis.set_tick_params(labelsize=17.5)
plt.show()
