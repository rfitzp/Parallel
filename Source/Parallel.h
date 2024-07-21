// Parallel.h

// #################################################################
// Class to solve time-dependent parallel electron transport problem
// in plasma of arbitrary collisionality
// #################################################################

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex>

#include <blitz/array.h>
#include <gsl/gsl_spline.h>
#include <netcdf.h>

#include "Faddeeva.h"

using namespace std;
using namespace blitz;
using namespace Faddeeva;

// Namelist reading function
extern "C" void NameListRead (double* mue, double* lambdaD, double* sigma,
			      double* xmax, int* Nx,
			      double* gmax, int* Ng,
			      double* tmax, int* Nt,
			      double* gamma, double* ksmax);

// ############
// Class header
// ############
class Parallel
{
 private:

  // ..................
  // Physics parameters
  // ..................
  double mue;     // Electron-ion collision frequency divided by total electron collision frequency (read from namelist)
  double lambdaD; // Electron Debye length divided by electron mean-free-path between collisions (read from namelist)
  double sigma;   // Width of particle and energy source deposition profiles (in units of the electron mean-free-path) (read from namelist)

  // ......................
  // Calculation parameters
  // ......................
  double xmax;    // Simulation extends from x = 0 to x = xmax (read from namelist)
  int    Nx;      // Number of equally spaced points on x and k grids (read from namelist)
  double gmax;    // Bromwich contour runs from g_i = 0 to g_i = g_max (read from namelist)
  int    Ng;      // Number of unequally spaced grid points on Bromwich contour (read from namelist)
  double tmax;    // Simulation extends from t = 0 to t = tmax (read from namelist)
  int    Nt;      // Number of equally spaced times between 0 and tmax (read from namelist)
  double gamma;   // Bromwich contour evaluated at g_r = gamma (read from namelist)
  double ksmax;   // Maximum value of k * sigma (read from namelist)

  // ................
  // Calculation data
  // ................
  Array<double,1>          xx;         // x grid points
  Array<double,1>          kk;         // k grid points
  Array<double,1>          gg;         // g grid points
  Array<double,1>          tt;         // t grid points

  Array<complex<double>,2> Fn0;        // Fn0 (x, g)
  Array<complex<double>,2> FT0;        // FT0 (x, g)
  Array<complex<double>,2> Fn2;        // Fn2 (x, g)
  Array<complex<double>,2> FT2;        // FT2 (x, g)

  Array<double,2>          ne0;        // delta n_e (x, t) for particle source
  Array<double,2>          Te0;        // delta T_e (x, t) for particle source
  Array<double,2>          ne2;        // delta n_e (x, t) for energy source
  Array<double,2>          Te2;        // delta T_e (x, t) for energy source

  Array<double,2>          ne0t;       // delta d/dt n_e (x, t) for particle source
  Array<double,2>          Te0t;       // delta d/dt T_e (x, t) for particle source
  Array<double,2>          ne2t;       // delta d/dt n_e (x, t) for energy source
  Array<double,2>          Te2t;       // delta d/dt T_e (x, t) for energy source
  
  Array<double,2>          Ve0;        // V_e (x, t) for particle source
  Array<double,2>          Qe0;        // Q_e (x, t) for particle source
  Array<double,2>          qe0;        // q_e (x, t) for particle source
  Array<double,2>          Ve2;        // V_e (x, t) for energy source
  Array<double,2>          Qe2;        // Q_e (x, t) for particle source
  Array<double,2>          qe2;        // q_e (x, t) for energy source

  Array<double,1>          Lne0;       // Spatial half-width of n_e (t, x) for particle source
  Array<double,1>          LTe0;       // Spatial half-width of T_e (t, x) for particle source
  Array<double,1>          Lne2;       // Spatial half-width of n_e (t, x) for energy source
  Array<double,1>          LTe2;       // Spatial half-width of T_e (t, x) for energy source

  Array<double,1>          Wne0;       // Spatial 90%-width of n_e (t, x) for particle source
  Array<double,1>          WTe0;       // Spatial 90%-width of T_e (t, x) for particle source
  Array<double,1>          Wne2;       // Spatial 90%-width of n_e (t, x) for energy source
  Array<double,1>          WTe2;       // Spatial 90%-width of T_e (t, x) for energy source
  
  Array<double,1>          ne00;       // n_e (t, 0) for particle source
  Array<double,1>          Te00;       // T_e (t, 0) for particle source
  Array<double,1>          ne20;       // n_e (t, 0) for energy source
  Array<double,1>          Te20;       // T_e (t, 0) for energy source
  
  gsl_interp_accel*        acc_n0r;    // Accelerator for Fn0_r
  gsl_interp_accel*        acc_n0i;    // Accelerator for Fn0_i
  gsl_interp_accel*        acc_T0r;    // Accelerator for Fn0_r
  gsl_interp_accel*        acc_T0i;    // Accelerator for Fn0_i
  gsl_interp_accel*        acc_n2r;    // Accelerator for Fn0_r
  gsl_interp_accel*        acc_n2i;    // Accelerator for Fn0_i
  gsl_interp_accel*        acc_T2r;    // Accelerator for Fn0_r
  gsl_interp_accel*        acc_T2i;    // Accelerator for Fn0_i

  gsl_spline*              spline_n0r; // Spline interpolator for Fn2_r
  gsl_spline*              spline_n0i; // Spline interpolator for Fn2_i
  gsl_spline*              spline_T0r; // Spline interpolator for FT2_r
  gsl_spline*              spline_T0i; // Spline interpolator for FT2_i
  gsl_spline*              spline_n2r; // Spline interpolator for Fn2_r
  gsl_spline*              spline_n2i; // Spline interpolator for Fn2_i
  gsl_spline*              spline_T2i; // Spline interpolator for FT2_i
  gsl_spline*              spline_T2r; // Spline interpolator for FT2_r
  
  // -------------------------------
  // Adaptive integration parameters
  // -------------------------------
  double acc;       // Integration accuracy
  double h0;        // Initial step-length
  double hmin;      // Minimum step-length
  double hmax;      // Maximum step-length
  int    maxrept;   // Maximum number of step recalculations
  int    flag;      // Integration error calcualation flag
  
  // ------------------
  // RK4/RK5 parameters
  // ------------------
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ----
  // Misc
  // ----
  int    count;
  double T;

public:

  // Constructor
  Parallel ();
  // Destructor
  ~Parallel ();

  // Solve problem
  void Solve ();

private:

  // Calculate inverse Laplace transform target functions
  void CalcF ();

  // Evaluate Z0
  complex<double> GetZ0 (complex<double> xi);
  // Evaluate G
  complex<double> GetG (complex<double> xi);
  // Evaluate LambdaD
  double GetLambdaD (double k);
  // Function to evaluate Kn0
  complex<double> GetKn0 (complex<double> g, double k);
  // Function to evaluate KT0
  complex<double> GetKT0 (complex<double> g, double k);
  // Function to evaluate Kn2
  complex<double> GetKn2 (complex<double> g, double k);
  // Function to evaluate KT2
  complex<double> GetKT2 (complex<double> g, double k);
   
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, Array<double,1> y, Array<double,1> dydx);
  // Adaptive step length RK4/RK5 integration routine
  void RK4RK5Adaptive (double& x, Array<double,1> y, double& h, 
		       double& t_err, double acc, double S, double T, int& rept,
		       int maxrept, double h_min, double h_max, int flag, 
		       int diag, FILE* file);
  // Fixed step length RK4/RK5 integration routine
  void RK4RK5Fixed (double& x, Array<double,1> y, Array<double,1> err, double h);
  
  // Open new file for writing
  FILE* OpenFilew (char* filename);
};
