import numpy as np
#from plasmapy.dispersion import plasma_dispersion_func, plasma_dispersion_func_deriv
import scipy.special as sp
import scipy.constants as sc
import math

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
def G_xi (xi):

    z0  = Z_0 (xi)
    xi2 = xi*xi

    return ((xi2 - 1.) - (xi2 - 1.5) * z0) /(2.*xi2 - (2.*xi2 - 1.) * z0)

"""
def G_xi(xi): # This is G_xi - g_1 / 2.
    Z_0 = -xi * plasma_dispersion_func(xi)
    numerator = (xi ** 2 - 1) - (xi ** 2 - 3 / 2) * Z_0
    denominator = 2 * xi ** 2 - (2 * xi ** 2 - 1) * Z_0

    return numerator / denominator
"""

# #########################
# Function to solve problem
# #########################

def solve_x_and_kpar(
        *,
        n        = 1,     # toroidal mode number
        s_s      = 0.7,   # magnetic shear
        R_0      = 1.67,  # major radius (m)
        r_s      = 0.3,   # minor radius (m)
        n_e      = 5e19,  # electron number density (m^-3)
        T_e      = 2.5,   # electron temperarature (keV)
        chi_perp = 1.0,   # perpendicular diffusivity (m^2/s)
        G_xi     = G_xi,
        x0       = 2.0,
        max_iter = 50,
        tol      = 1e-8,
        verbose  = False,
):

    # This function returns 3 different results:
    # 1. x, which is the ratio of W_c and r_s;
    # 2. chi_parallel;
    # 3. the ratio of chi_parallel and chi_perpendicular as well as their string type.
    
    #  G_xi：The complicated function we got
    
    if G_xi is None:
        G_xi = lambda z: 1.0

    # Derived quantities
    Te_raw     = T_e * sc.e * 1.e3 # Convert electron temperature to energy units
    vt_e       = (2.*Te_raw /sc.m_e)**0.5
    tau_e      = 6.*2.**0.5 * math.pi**1.5 * sc.epsilon_0**2 * sc.m_e**0.5 * Te_raw**1.5 /15. /sc.e**4 /n_e
    l_e        = vt_e * tau_e
    kappa_perp = n_e * chi_perp
    epsilon_s  = r_s / R_0
    alpha      = n * s_s * l_e / (2.0 * np.pi * R_0)
    beta       = (2.0 * np.pi * alpha) / (n * epsilon_s * s_s) * np.sqrt(kappa_perp / (n_e * vt_e * l_e))

    # Iteration
    x = float(x0)
    history = [x]
    for it in range(int(max_iter)):

        # rhs = beta * |G(i/(alpha*x))|^{-1/2}

        val_G = abs(G_xi(1j / (alpha * x)))
        rhs   = beta * (val_G ** (-0.5))

        if verbose:
            print(f"iter {it:2d}: x = {rhs:.6e}   |G| = {val_G:.6e}")

        # Converge or not
        if abs(rhs - x) < tol * max(1.0, abs(x)):
            x = rhs
            history.append(x)
            iters = it + 1
            break

        x = rhs
        history.append(x)
    else:
        iters = max_iter  
   
    kappa_parallel = ((beta ** 2) / (alpha ** 2) / (x ** 4)) * n_e * vt_e * l_e
    chi_parallel   = ((beta ** 2) / (alpha ** 2) / (x ** 4)) * vt_e * l_e

    result = dict(
        x              = x,
        W_c            = x * r_s,
        kappa_parallel = kappa_parallel,
        alpha          = alpha,
        beta           = beta,
        eps            = epsilon_s,
        iters          = iters,
        history        = history,
    )
    return (x, chi_parallel, chi_parallel/chi_perp), \
           (f"{x:.3e}", f"{chi_parallel:.3e}", f"{chi_parallel/chi_perp:.3e}")

print(solve_x_and_kpar(
    n        = 1,
    s_s      = 1.0,
    R_0      = 6.2,
    r_s      = 1.2,
    n_e      = 1e20,
    T_e      = 7.0,
    chi_perp = 1.0,
    G_xi     = G_xi,
    x0       = 2.0,
    max_iter = 50,
    tol      = 1e-8,
    verbose  = False,))
