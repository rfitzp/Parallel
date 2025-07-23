# Island.py

# ######################################
# Script to investigate island transport
# ######################################

import math
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

sqpi = math.sqrt(math.pi)

# ##############
# Function Z_bar
# ##############
def Z_bar (xi):

    res = 1j * sqpi * sp.wofz(xi)

    if (xi.imag < 0.):
        res -= 2. * 1j * sqpi * math.exp(-xi*xi)

    return res

# ############
# Function Z_0
# ############
def Z_0 (xi):

    return - xi * Z_bar (xi)

# ##########
# Function G
# ##########
def G (xi):

    z0  = Z_0 (xi)
    xi2 = xi*xi

    return ((xi2 - 1.) - (xi2 - 1.5) * z0) /(2.*xi2 - (2.*xi2 - 1.) * z0)

# #############################
# Function to return tan(alpha)
# #############################
def tan_alpha (L, w):

    k0 = math.pi /L
    xi = w/k0 + 1j/k0

    GG = G(xi)
    Gr = GG.real
    Gi = GG.imag

    return (w/2. - w*Gr + Gi) /(-Gr - w* Gi)

# ###################################################################
# Function to return short mean-free-path approximation to tan(alpha)
# ###################################################################
def tan_alpha_smfp (L, w):

    k0 = math.pi /L

    return (2./3.) * w * (1. + w*w) /k0/k0

# ###################################################################
# Function to return short long-free-path approximation to tan(alpha)
# ###################################################################
def tan_alpha_lmfp (L, w):

    k0 = math.pi /L

    return sqpi * (4./math.pi - 1.) * w /k0

# ###########################
# Function to return W_c/W_c0
# ###########################
def W_c (L, w):

    k0 = math.pi /L
    xi = w/k0 + 1j/k0

    talpha = tan_alpha (L, w)
    calpha = 1. /math.sqrt (1. + talpha*talpha)

    GG = G(xi)
    Gr = GG.real
    Gi = GG.imag

    return (0.75 * k0*k0 * calpha /(- Gr - w*Gi))**0.25

# #################################################################
# Function to return short mean-free-path approximation to W_c/W_c0
# #################################################################
def W_c_smfp (L, w):

    k0 = math.pi /L

    return (1.5 * k0*k0)**0.25 /(w*w + 2.25 * k0**4 /(1. + w*w)**2)**0.125

# ################################################################
# Function to return long mean-free-path approximation to W_c/W_c0
# ################################################################
def W_c_lmfp (L, w):

    k0 = math.pi /L

    return (0.75 * k0*k0)**0.25 /((4./math.pi - 1.)**2 * w*w + (k0/sqpi + 4./math.pi - 0.75)**2)**0.125

# ###############################################################################
# Function to return long mean-free-path zero-frequency approximation to W_c/W_c0
# ###############################################################################
def W_c_lmfp_z (L, w):

    k0 = math.pi /L

    #return (0.75*k0*k0 /(k0/sqpi + 4./math.pi - 0.75))**0.25
    return (0.75*k0*k0 /(k0/sqpi))**0.25

# #########################
# Plot zero frequency limit
# #########################

L_min = 0.01
L_max = 100.
L_num = 10000

LL = np.linspace (L_min, L_max, L_num)

xa = []
xb = []
for L in LL:
    xa.append(W_c       (L, 0.))
    xb.append(W_c_lmfp_z (L, 0.))

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font)

plt.xlim (L_min, L_max)

plt.xscale('log') 
plt.plot    (LL, xa, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$W_c/W_{c\,0}$')
plt.plot    (LL, xb, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$W_{c\,3}/W_{c\,0}$')
plt.axhline (1.,     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{L}$',fontsize = font)
plt.legend (fontsize = font)

#plt.show()
plt.savefig ("zero.pdf")

# ###############################
# Plot short-mean-free path limit
# ###############################

w_min = 10.
w_max = 50.
w_num = 1000

ww = np.linspace (w_min, w_max, w_num)

y1 = []
y2 = []
x1 = []
x2 = []
for w in ww:
    y1.append(tan_alpha      (1., w))
    y2.append(tan_alpha_smfp (1., w))
    x1.append(W_c            (1., w))
    x2.append(W_c_smfp       (1., w))

L_min = 1.
L_max = 20.
L_num = 1000

LL = np.linspace (L_min, L_max, L_num)

y3 = []
y4 = []
x3 = []
x4 = []
for L in LL:
    y3.append(tan_alpha      (L, 1.))
    y4.append(tan_alpha_smfp (L, 1.))
    x3.append(W_c            (L, 1.))
    x4.append(W_c_smfp       (L, 1.))
  
font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (2, 2, 1)

plt.xlim (w_min, w_max)
 
plt.plot    (ww, y1, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$\tan(\alpha)$')
plt.plot    (ww, y2, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$\tan(\alpha_{smfp})$')
plt.axhline (0.,     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{\omega}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 2)

plt.xlim (L_min, L_max)
 
plt.plot    (LL, y3, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$\tan(\alpha)$')
plt.plot    (LL, y4, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$\tan(\alpha_{smfp})$')
plt.axhline (0.,     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{L}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 3)

plt.xlim (w_min, w_max)
 
plt.plot    (ww, x1, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$W_c/W_{c\,0}$')
plt.plot    (ww, x2, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$W_{c\,smfp}/W_{c\,0}$')
plt.axhline (0.,     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{\omega}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 4)

plt.xlim (L_min, L_max)
 
plt.plot    (LL, x3, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$W_c/W_{c\,0}$')
plt.plot    (LL, x4, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$W_{c\,smfp}/W_{c\,0}$')
plt.axhline (0.,     color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{L}$',fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()
plt.savefig ("smfp.pdf")

# #########################
# Long-mean-free path limit
# #########################

w_min = 0.
w_max = 0.1
w_num = 1000

ww = np.linspace (w_min, w_max, w_num)

yy1 = []
yy2 = []
xx1 = []
xx2 = []
for w in ww:
    yy1.append(tan_alpha      (0.1, w))
    yy2.append(tan_alpha_lmfp (0.1, w))
    xx1.append(W_c            (0.1, w))
    xx2.append(W_c_lmfp       (0.1, w))

L_min = 1.e-4
L_max = 0.1
L_num = 1000

LL = np.linspace (L_min, L_max, L_num)

yy3 = []
yy4 = []
xx3 = []
xx4 = []
for L in LL:
    yy3.append(tan_alpha      (L, 0.1))
    yy4.append(tan_alpha_lmfp (L, 0.1))
    xx3.append(W_c            (L, 0.1))
    xx4.append(W_c_lmfp       (L, 0.1))

font = 20
fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = font) 
plt.rc ('ytick', labelsize = font) 

plt.subplot (2, 2, 1)

plt.xlim (w_min, w_max)
 
plt.plot    (ww, yy1, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$\tan(\alpha)$')
plt.plot    (ww, yy2, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$\tan(\alpha_{lmfp})$')
plt.axhline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{\omega}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 2)

plt.xlim (L_min, L_max)
 
plt.plot    (LL, yy3, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$\tan(\alpha)$')
plt.plot    (LL, yy4, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$\tan(\alpha_{smfp})$')
plt.axhline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{L}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 3)

plt.xlim (w_min, w_max)
 
plt.plot    (ww, xx1, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$W_c/W_{c\,0}$')
plt.plot    (ww, xx2, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$W_{c\,lmfp}/W_{c\,0}$')
plt.axhline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{\omega}$',fontsize = font)
plt.legend (fontsize = font)

plt.subplot (2, 2, 4)

plt.xlim (L_min, L_max)
 
plt.plot    (LL, xx3, color = 'blue',  linewidth = 2,  linestyle = 'solid',  label = r'$W_c/W_{c\,0}$')
plt.plot    (LL, xx4, color = 'red',   linewidth = 2,  linestyle = 'dashed', label = r'$W_{c\,lmfp}/W_{c\,0}$')
plt.axhline (0.,      color = 'black', linewidth = 1., linestyle = 'dotted')

plt.xlabel (r'$\hat{L}$',fontsize = font)
plt.legend (fontsize = font)

plt.tight_layout ()
plt.savefig ("lmfp.pdf")
