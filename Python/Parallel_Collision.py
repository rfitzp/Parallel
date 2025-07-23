#!/usr/bin/env python
# coding: utf-8

# In[1]:Define the function


import matplotlib.pyplot as plt
import numpy as np
from plasmapy.dispersion import plasma_dispersion_func, plasma_dispersion_func_deriv


def G_xi(xi): # This is G_xi.
    Z_0 = -xi * plasma_dispersion_func(xi)
    numerator = (xi ** 2 - 1) - (xi ** 2 - 3 / 2) * Z_0
    denominator = 2 * xi ** 2 - (2 * xi ** 2 - 1) * Z_0

    return numerator / denominator

def tan_alpha(w, L):
    k = np.pi / L
    xi = 1j * (1 - 1j * w) / k
    Gi = G_xi(xi).imag
    Gr = G_xi(xi).real
    numerator = w / 2 - w * Gr + Gi
    denominator = - Gr - w * Gi
    return numerator / denominator

def Wc_Wc0(w, L):
    k = np.pi / L
    tan_a = tan_alpha(w, L)
    cos_a = np.cos(np.arctan(tan_a))
    numerator = 3 / 4 * k ** 2 * cos_a
    xi = 1j * (1 - 1j * w) / k
    Gi = G_xi(xi).imag
    Gr = G_xi(xi).real
    denominator = - Gr - w * Gi
    return (numerator / denominator) ** (1/4)


# In[2]:Plots for Testing


L = np.linspace(0.001, 10, 500)
w = 0.1

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
ax1.plot(L, tan_alpha(w, L), label=fr'$\tan\alpha$ - $L$ at $\omega = {w}$')
ax1.set_xlabel(r'$L$ (wavelength)')
ax1.set_ylabel(r'$\tan \alpha$')
ax1.set_title(r'$\tan\alpha$ as a Function of $L$ for Fixed $\omega$ (General)')
ax1.grid(True)
ax1.legend()

ax2.plot(L, Wc_Wc0(w,L), label=fr'$W_c/W_{{c0}}$ - $L$ at $\omega = {w}$')
ax2.set_xlabel(r'$L$ (wavelength)')
ax2.set_ylabel(r'$W_c/W_{c0}$')
ax2.set_title(r'$W_c/W_{c0}$ as a Function of $L$ for Fixed $\omega$ (General)')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig("Fix_w_1.png", dpi=300)
plt.show()







# In[3]:Plots of Short-Free-Path Limit


plt.figure(figsize=(16, 10))

plt.subplot(2,2,1)
w = np.linspace(10, 50, 500)
L = 1
plt.plot(w, tan_alpha(w, L), label=fr'$\tan\alpha$ - $\omega$ at $L = {L}$')
plt.plot(w, 2/3 * w * (1+w**2) / (np.pi)**2 * L **2, label = r'$ \frac{2}{3}\frac{\omega(1+\omega^2)}{k_0^2} $')
plt.xlabel(r'$\omega$', fontsize=12)
plt.ylabel(r'$\tan \alpha$', fontsize=12)
plt.title(r'$\tan\alpha$ as a Function of $\omega$ for Fixed $L$', fontsize=14)
plt.grid(True)
plt.legend()

plt.subplot(2,2,2)
w = 10
L = np.linspace(1,10,500)
plt.plot(L, tan_alpha(w, L), label=fr'$\tan\alpha$ - $L$ at $\omega = {w}$')
plt.plot(L, 2/3 * w * (1+w**2) / (np.pi)**2 * L **2, label = r'$ \frac{2}{3}\frac{\omega(1+\omega^2)}{k_0^2} $')
plt.xlabel(r'$L$', fontsize=12)
plt.ylabel(r'$\tan \alpha$', fontsize=12)
plt.title(r'$\tan\alpha$ as a Function of $L$ for Fixed $\omega$', fontsize=14)
plt.grid(True)
plt.legend()

plt.subplot(2,2,3)
w = np.linspace(0.1, 50, 500)
L = 1
plt.plot(w, Wc_Wc0(w, L), label=fr'$W_c/W_{{c0}}$ - $\omega$ at $L = {L}$')
plt.plot(w, (3/2*(np.pi/L)**2/w)**0.25, label = r'$ (\frac{3k_0^2}{2\omega})^{1/4} $')
plt.xlabel(r'$\omega$', fontsize=12)
plt.ylabel(r'$W_c/W_{{c0}}$', fontsize=12)
plt.title(r'$W_c/W_{{c0}}$ as a Function of $\omega$ for Fixed $L$', fontsize=14)
plt.grid(True)
plt.legend()

plt.subplot(2,2,4)
w = 10
L = np.linspace(1,20,500)
plt.plot(L, Wc_Wc0(w, L), label=fr'$W_c/W_{{c0}}$ - $L$ at $\omega = {w}$')
plt.plot(L,(3/2*(np.pi/L)**2/w)**0.25, label = r'$ (\frac{3k_0^2}{2\omega})^{1/4} $')
plt.xlabel(r'$L$', fontsize=12)
plt.ylabel(r'$W_c/W_{{c0}}$', fontsize=12)
plt.title(r'$W_c/W_{{c0}}$ as a Function of $L$ for Fixed $\omega$', fontsize=14)
plt.grid(True)
plt.legend()
plt.suptitle('Short Mean-Free-Path Limit',fontsize=20)
plt.tight_layout()
plt.savefig("Short_mfp1.png", dpi=300)
plt.show()


# In[4]:Plots of Long Mean-Free-Path Limit


plt.figure(figsize=(16, 10))

plt.subplot(2,2,1)
L = 0.1
w = np.linspace(0.001, 0.1, 500)
k0 = np.pi / L
plt.plot(w, tan_alpha(w, L), label=fr'$\tan\alpha$ - $\omega$ at $L = {L}$')
plt.plot(w, np.sqrt(np.pi) * w / k0 *(4/np.pi - 1), label = r'$ \frac{\pi^{1/2}\omega}{k_0}(4/\pi-1)$')
plt.xlabel(r'$\omega$', fontsize=12)
plt.ylabel(r'$\tan \alpha$', fontsize=12)
plt.title(r'$\tan\alpha$ as a Function of $\omega$ for Fixed $L$', fontsize=14)
plt.grid(True)
plt.legend()

plt.subplot(2,2,2)
L = np.linspace(0.001, 0.1, 500)
w = 0.1
k0 = np.pi / L
plt.plot(L, tan_alpha(w, L), label=fr'$\tan\alpha$ - $L$ at $\omega = {w}$')
plt.plot(L, np.sqrt(np.pi) * w / k0 *(4/np.pi - 1), label = r'$ \frac{\pi^{1/2}\omega}{k_0}(4/\pi-1)$')
plt.xlabel(r'$L$', fontsize=12)
plt.ylabel(r'$\tan \alpha$', fontsize=12)
plt.title(r'$\tan\alpha$ as a Function of $L$ for Fixed $\omega$ ', fontsize=14)
plt.grid(True)
plt.legend()

plt.subplot(2,2,3)
L = 0.1
w = np.linspace(0.1, 1, 500)
k0 = np.pi / L
plt.plot(w, Wc_Wc0(w, L), label=fr'$W_c/W_{{c0}}$ - $\omega$ at $L = {L}$')
plt.plot(w, (3/4*((np.pi)/L)**2 / (np.sqrt(np.pi)/L) *w/w)**0.25, label = r'$(\frac{3\pi^{1/2}k_0}{4})^{1/4}$')
plt.xlabel(r'$\omega$', fontsize=12)
plt.ylabel(r'$W_c/W_{{c0}}$', fontsize=12)
plt.ylim(1,3)
plt.title(r'$W_c/W_{{c0}}$ as a Function of $\omega$ for Fixed $L$', fontsize=14)
plt.grid(True)
plt.legend()


plt.subplot(2,2,4)
L = np.linspace(0.1, 1, 500)
w = 0.1
k0 = np.pi / L
plt.plot(L, Wc_Wc0(w, L), label=fr'$W_c/W_{{c0}}$ - $L$ at $\omega = {w}$')
plt.plot(L, (3/4*((np.pi)/L)**2 / (np.sqrt(np.pi)/L) *w/w)**0.25, label = r'$(\frac{3\pi^{1/2}k_0}{4})^{1/4}$')
plt.xlabel(r'$L$', fontsize=12)
plt.ylabel(r'$W_c/W_{{c0}}$', fontsize=12)
plt.title(r'$W_c/W_{{c0}}$ as a Function of $L$ for Fixed $\omega$', fontsize=14)
plt.grid(True)
plt.legend()

plt.suptitle('Long Mean-Free-Path Limit',fontsize=20)
plt.tight_layout()
plt.savefig("Long_mfp1.png", dpi=300)
plt.show()


# In[5]:Contour of tan_alpha


# Grid
num = 2
x = np.linspace(0, num, 400)
y = np.linspace(0, num, 400)
w, L = np.meshgrid(x, y)
Z = tan_alpha(w, L)

# Contour
plt.figure(figsize=(8, 6))
contour = plt.contourf(w, L, Z, 360, cmap='bwr')
cbar = plt.colorbar(contour, ticks=np.arange(0, 1.4, 0.2))  # 从 1.6 到 2.0，步长 0.2
cbar.set_label(r'$\tan\alpha$')
plt.contour (w, L, Z, levels = 20, colors = 'black', linestyles = 'solid', linewidths = 0.5)
plt.title(r'Contour plot of $\tan\alpha$')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$L$')
plt.grid(True)
plt.savefig("tan_Contour.png", dpi=300)
plt.show()


# In[6]:Contour of W_c/W_{c0}


# Grid
start1 = 0
end1 = 5
start2 = 0.1
end2 = 5
x = np.linspace(start1, end1, 400)
y = np.linspace(start2, end2, 400)
w, L = np.meshgrid(x, y)
Z = Wc_Wc0(w, L)

# Contour
plt.figure(figsize=(8, 6))

contour = plt.contourf(w, L, Z, 360, cmap='bwr')
cbar = plt.colorbar(contour, ticks=np.arange(0.6, 2.6, 0.2))  # 从 1.6 到 2.0，步长 0.2
cbar.set_label(r'$W_c/W_{c0}$')
plt.contour (w, L, Z, levels = 20, colors = 'black', linestyles = 'solid', linewidths = 0.5)
plt.title(r'Contour plot of $W_c/W_{c0}$')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$L$')
plt.grid(True)
plt.savefig("W0_Contour.png", dpi=300)
plt.show()


# In[7]:Zero-Frequency's Plot of W_c/W_{c0} 


L = np.linspace(0.1, 20, 500)
w = 0
k0 = np.pi / L
plt.plot(L, Wc_Wc0(w, L), label=fr'$W_c/W_{{c0}}$ - $L$ at $\omega = {w}$')
plt.plot(L, (3/4*((np.pi)/L)**2 / (np.sqrt(np.pi)/L) *(w+1)/(w+1))**0.25, label = r'$(\frac{3\pi^{1/2}k_0}{4})^{1/4}$')
plt.xlabel(r'$L$', fontsize=12)
plt.ylabel(r'$W_c/W_{{c0}}$', fontsize=12)
plt.title(r'$W_c/W_{{c0}}$ as a Function of $L$ for $\omega = 0$', fontsize=14)
plt.grid(True)
plt.legend()
plt.savefig("W0_wis0.png", dpi=300)
plt.show()




