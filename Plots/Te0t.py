import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fn = 'Simulation.nc'
ds = nc.Dataset(fn)
x  = ds['x']
t  = ds['t']
n  = ds['T0_0']
l  = ds['LT_0']
w  = ds['WT_0']

xx = np.asarray(x)
tt = np.asarray(t)
nn = np.asarray(n)
ll = np.asarray(l)
ww = np.asarray(w)

fig = plt.figure (figsize=(12.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (1, 2, 1)

plt.xlim (0, tt[-1])
plt.ylim (1.05*nn[-1], 0)

plt.plot (tt, nn, color = 'red', linewidth = 1.5, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize="20")
plt.ylabel(r'$\hat{T}_e$', fontsize="20")

plt.subplot (1, 2, 2)

plt.xlim (0, tt[-1])
plt.ylim (0, 1.05*ww[-1])

plt.plot (tt, ll, color = 'red', linewidth = 1.5, linestyle = 'solid')
plt.plot (tt, ww, color = 'blue', linewidth = 1.5, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize="20")
plt.ylabel(r'$\hat{L}$', fontsize="20")

plt.tight_layout(pad=0.5)

plt.show ()
