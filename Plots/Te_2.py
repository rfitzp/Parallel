#matplotlib widget
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fn = 'Simulation.nc'
ds = nc.Dataset(fn)
x  = ds['x']
t  = ds['t']
n  = ds['Te_2']

xx = np.asarray(x)
tt = np.asarray(t)
nn1 = np.asarray(n)
b = (1-0.5)/2
c = 1/np.sqrt(2)*100
nn2 = tt
nn2 = nn2.reshape(1,501)

fig, ax = plt.subplots()
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.xlim (-0.01, 0.01);
plt.ylim (1.05*np.min(nn1), 1.05*np.max(nn1))
plt.axhline (0., color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0., color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel (r'$\hat{x}$',         fontsize = 17)
plt.ylabel (r'$\delta\hat{n}_e$', fontsize = 17)

line1, = ax.plot ([], [], linewidth = 1.0, color = 'red', label='Numerical')
line2, = ax.plot([], [], linewidth=1.0, color='blue', label='Theoretical')

ax.legend()

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2

def update(i):
    print(f"Frame {i}: Updating line with data")
    line1.set_data(xx, nn1[:, i])  
    line2.set_data(xx, nn2[:, i])  
    return line1, line2
    
ani = animation.FuncAnimation (fig, update, frames = range(tt.size), init_func = init, interval = 100, blit = True, repeat = False)

plt.show ()
