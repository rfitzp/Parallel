import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fn    = 'Simulation.nc'
ds    = nc.Dataset(fn)
x     = ds['x']
t     = ds['t']
n     = ds['Qe_0']

xx = np.asarray(x)
tt = np.asarray(t)
nn = np.asarray(n)

fig, ax = plt.subplots()
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.xlim (xx[0], xx[-1]);
plt.ylim (1.05*np.min(nn), 1.05*np.max(nn))
plt.axhline (0., color='black', linewidth=2., linestyle='dotted')
plt.axvline (0., color='black', linewidth=2., linestyle='dotted')

plt.xlabel(r'$\hat{x}$', fontsize="20")
plt.ylabel(r'$\hat{Q}_e$', fontsize="20")

line, = ax.plot([], [], linewidth='2.0', color = 'black')

def init ():
    return line,

def update(i):
    line.set_data(xx, nn[:,i])
    return line,

ani = animation.FuncAnimation(fig, update, frames = range(tt.size), init_func = init, interval=100, blit=True)

plt.show ()
