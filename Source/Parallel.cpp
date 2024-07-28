// Parallel.cpp

#include "Parallel.h"

// ###########
// Constructor
// ###########
Parallel::Parallel ()
{
  // .............
  // Read namelist
  // .............
  NameListRead (&mue, &lambdaD, &sigma, &xmax, &Nx, &gmax, &Ng, &tmax, &Nt, &gamma, &ksmax, &flg);

  printf ("\nmue = %10.3e  lambdaD = %10.3e  sigma = %10.3e\n", mue, lambdaD, sigma);
  printf ("\nxmax = %10.3e  Nx = %4d\n", xmax, Nx);
  printf ("\ngmax = %10.3e  Ng = %4d\n", gmax, Ng);
  printf ("\ntmax = %10.3e  Nt = %4d\n", tmax, Nt);
  printf ("\ngamma = %10.3e  ksmax = %10.3e\n", gamma, ksmax);
  printf ("\nflg = %1d\n", flg);
 
  // -----------------------------------
  // Set adaptive integration parameters
  // -----------------------------------
  acc     = 1.e-10;
  h0      = 1.e-2;
  hmin    = 1.e-10;
  hmax    = 1.e-1;
  maxrept = 50;
  flag    = 2;

  // ----------------------
  // Set RK4/RK5 parameters
  // ----------------------
  aa1 = 0.;
  aa2 = 1./5.;
  aa3 = 3./10.;
  aa4 = 3./5.;
  aa5 = 1.;
  aa6 = 7./8.;

  cc1 =  37./378.;
  cc3 = 250./621.;
  cc4 = 125./594.;
  cc6 = 512./1771.;

  ca1 = cc1 -  2825./27648.;
  ca3 = cc3 - 18575./48384.;
  ca4 = cc4 - 13525./55296.;
  ca5 =     -   277./14336.;
  ca6 = cc6 -     1./4.;

  bb21 = 1./5.;

  bb31 = 3./40.;
  bb32 = 9./40.;

  bb41 =   3./10.;
  bb42 = - 9./10.;
  bb43 =   6./5.;

  bb51 = - 11./54.;
  bb52 =    5./2.;
  bb53 = - 70./27.;
  bb54 =   35./27.;

  bb61 =  1631./55296.;
  bb62 =   175./512.;
  bb63 =   575./13824.;
  bb64 = 44275./110592.;
  bb65 =   253./4096.;
}

// ##########
// Destructor
// ##########
Parallel::~Parallel ()
{
}

// #########################
// Function to solve problem
// #########################
void Parallel::Solve ()
{
  // ...............
  // Allocate memory
  // ...............
  xx.resize (Nx+1);
  kk.resize (Nx+1);
  gg.resize (Ng+1);
  tt.resize (Nt+1);

  KKn0.resize (Nx+1, Ng+1);
  KKT0.resize (Nx+1, Ng+1);
  KKn2.resize (Nx+1, Ng+1);
  KKT2.resize (Nx+1, Ng+1);
  
  Fn0.resize (Nx+1, Ng+1);
  FT0.resize (Nx+1, Ng+1);
  Fn2.resize (Nx+1, Ng+1);
  FT2.resize (Nx+1, Ng+1);

  ne0.resize (Nx+1, Nt+1);
  Te0.resize (Nx+1, Nt+1);
  ne2.resize (Nx+1, Nt+1);
  Te2.resize (Nx+1, Nt+1);

  ne0t.resize (Nx+1, Nt+1);
  Te0t.resize (Nx+1, Nt+1);
  ne2t.resize (Nx+1, Nt+1);
  Te2t.resize (Nx+1, Nt+1);

  Ve0.resize (Nx+1, Nt+1);
  Qe0.resize (Nx+1, Nt+1);
  qe0.resize (Nx+1, Nt+1);
  Ve2.resize (Nx+1, Nt+1);
  Qe2.resize (Nx+1, Nt+1);
  qe2.resize (Nx+1, Nt+1);

  Lne0.resize (Nt+1);
  LTe0.resize (Nt+1);
  Lne2.resize (Nt+1);
  LTe2.resize (Nt+1);

  Wne0.resize (Nt+1);
  WTe0.resize (Nt+1);
  Wne2.resize (Nt+1);
  WTe2.resize (Nt+1);

  ne00.resize (Nt+1);
  Te00.resize (Nt+1);
  ne20.resize (Nt+1);
  Te20.resize (Nt+1);

  acc_n0r = gsl_interp_accel_alloc ();
  acc_n0i = gsl_interp_accel_alloc ();
  acc_T0r = gsl_interp_accel_alloc ();
  acc_T0i = gsl_interp_accel_alloc ();
  acc_n2r = gsl_interp_accel_alloc ();
  acc_n2i = gsl_interp_accel_alloc ();
  acc_T2r = gsl_interp_accel_alloc ();
  acc_T2i = gsl_interp_accel_alloc ();

  spline_n0r = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_n0i = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_T0r = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_T0i = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_n2r = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_n2i = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_T2r = gsl_spline_alloc (gsl_interp_cspline, Ng+1);
  spline_T2i = gsl_spline_alloc (gsl_interp_cspline, Ng+1);

  // ...........................
  // Set up x, k, g, and t grids
  // ...........................
  for (int j = 0; j <= Nx; j++)
    {
      double fx = double (j) /double (Nx);
      double x  = xmax * fx;

      xx(j) = x;
    }
  for (int j = 0; j <= Nx; j++)
    {
      double fk   = double (j) /double (Nx);
      double kmax = ksmax /sigma;
      double k    = kmax * fk;

      kk(j) = k;
    }
  for (int j = 0; j <= Ng; j++)
    {
      double fg = double (j) /double (Ng);
      double g  = gmax * fg * fg * fg;

      gg(j) = g;
    }
  for (int j = 0; j <= Nt; j++)
    {
      double ft = double (j) /double (Nt);
      double t  = tmax * ft;

      tt(j) = t;
    }

  // ................................
  // Calculate Kn0, KT0, Kn2, and KT2
  // ................................
  for (int i = 0; i <= Nx; i++)
    for (int j = 0; j <= Ng; j++)
      {
	if (i == 0)
	  {
	    KKn0 (i, j) = GetKn0 (complex<double> (gamma, gg(j)), kk(1));
	    KKT0 (i, j) = GetKT0 (complex<double> (gamma, gg(j)), kk(1));
	    KKn2 (i, j) = GetKn2 (complex<double> (gamma, gg(j)), kk(1));
	    KKT2 (i, j) = GetKT2 (complex<double> (gamma, gg(j)), kk(1));
	  }
	else
	  {
	    KKn0 (i, j) = GetKn0 (complex<double> (gamma, gg(j)), kk(i));
	    KKT0 (i, j) = GetKT0 (complex<double> (gamma, gg(j)), kk(i));
	    KKn2 (i, j) = GetKn2 (complex<double> (gamma, gg(j)), kk(i));
	    KKT2 (i, j) = GetKT2 (complex<double> (gamma, gg(j)), kk(i));
	  }
      }

  // ....................................................
  // Calculate inverse Laplace transform target functions
  // ....................................................
  CalcF ();

  // ..................................
  // Perform inverse Laplace transforms
  // ..................................
  printf ("Calculating inverse Laplace transforms:\n");
  for (int i = 0; i <= Nx; i++)
    {
      double x = xx(i);
      
      // Initialize interpolation
      double* gx = new double[Ng+1];
      double* Fr = new double[Ng+1];
      double* Fi = new double[Ng+1];
      
      for (int l = 0; l <= Ng; l++)
	{
	  gx[l] = gg(l);
	  if (flg)
	    {
	      Fr[l] = real (KKn0(i, l));
	      Fi[l] = imag (KKn0(i, l));
	    }
	  else
	    {
	      Fr[l] = real (Fn0(i, l));
	      Fi[l] = imag (Fn0(i, l));
	    }
	}
      gsl_spline_init (spline_n0r, gx, Fr, Ng+1);
      gsl_spline_init (spline_n0i, gx, Fi, Ng+1);

      for (int l = 0; l <= Ng; l++)
	{
	  if (flg)
	    {
	      Fr[l] = real (KKT0(i, l));
	      Fi[l] = imag (KKT0(i, l));
	    }
	  else
	    {
	      Fr[l] = real (FT0(i, l));
	      Fi[l] = imag (FT0(i, l));
	    }
	}
      gsl_spline_init (spline_T0r, gx, Fr, Ng+1);
      gsl_spline_init (spline_T0i, gx, Fi, Ng+1);

      for (int l = 0; l <= Ng; l++)
	{
	  if (flg)
	    {
	      Fr[l] = real (KKn2(i, l));
	      Fi[l] = imag (KKn2(i, l));
	    }
	  else
	    {
	      Fr[l] = real (Fn2(i, l));
	      Fi[l] = imag (Fn2(i, l));
	    }
	}
      gsl_spline_init (spline_n2r, gx, Fr, Ng+1);
      gsl_spline_init (spline_n2i, gx, Fi, Ng+1);

      for (int l = 0; l <= Ng; l++)
	{
	  if (flg)
	    {
	      Fr[l] = real (KKT2(i, l));
	      Fi[l] = imag (KKT2(i, l));
	    }
	  else
	    {
	      Fr[l] = real (FT2(i, l));
	      Fi[l] = imag (FT2(i, l));
	    }
	}
      gsl_spline_init (spline_T2r, gx, Fr, Ng+1);
      gsl_spline_init (spline_T2i, gx, Fi, Ng+1);

      delete[] gx; delete[] Fr; delete[] Fi;
      
      for (int j = 0; j <= Nt; j++)
	{
	  if (i%10 == 0 && j%10 == 0)
	    printf (".");
	  fflush (stdout);
    
	  T = tt(j);
	  
	  double          g, h, t_err;
	  int             rept; 
	  Array<double,1> y  (4);
	  Array<double,1> err(4);
	  
	  g     = 0.;
	  h     = h0;
	  count = 0;
	  y(0)  = 0.;
	  y(1)  = 0.;
	  y(2)  = 0.;
	  y(3)  = 0.;
	  	  
	  do
	    {
	      RK4RK5Adaptive (g, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag, 0, NULL);
	    }
	  while (g + h < gmax);
	  RK4RK5Fixed (g, y, err, gmax - g);

	  ne0(i, j) = y(0);
	  Te0(i, j) = y(1);
	  ne2(i, j) = y(2);
	  Te2(i, j) = y(3);
	}
      if (i%10 == 0)
	printf ("%4d\n", i);
    }
 
  // ..........................
  // Calculate time derivatives
  // ..........................
  double h = tt(1) - tt(0);
  for (int i = 0; i <= Nx; i++)
    {
      ne0t(i, 0) = (ne0(i, 1) - ne0(i, 0)) /h;
      Te0t(i, 0) = (Te0(i, 1) - Te0(i, 0)) /h;
      ne2t(i, 0) = (ne2(i, 1) - ne2(i, 0)) /h;
      Te2t(i, 0) = (Te2(i, 1) - Te2(i, 0)) /h;
      
      for (int j = 1; j < Nt; j++)
	{
	  ne0t(i, j) = (ne0(i, j+1) - ne0(j, j-1)) /2./h;
	  Te0t(i, j) = (Te0(i, j+1) - Te0(j, j-1)) /2./h;
	  ne2t(i, j) = (ne2(i, j+1) - ne2(j, j-1)) /2./h;
	  Te2t(i, j) = (Te2(i, j+1) - Te2(j, j-1)) /2./h;
	}

      ne0t(i, Nt) = (ne0(i, Nt) - ne0(i, Nt-1)) /h;
      Te0t(i, Nt) = (Te0(i, Nt) - Te0(i, Nt-1)) /h;
      ne2t(i, Nt) = (ne2(i, Nt) - ne2(i, Nt-1)) /h;
      Te2t(i, Nt) = (Te2(i, Nt) - Te2(i, Nt-1)) /h;
    }

  // ................
  // Calculate fluxes
  // ................
  for (int j = 0; j <= Nt; j++)
    {
      Ve0(0, j) = 0.;
      Qe0(0, j) = 0.;
      qe0(0, j) = 0.;
      Ve2(0, j) = 0.;
      Qe2(0, j) = 0.;
      qe2(0, j) = 0.;

      h = xx(1) - xx(0);
      for (int i = 1; i <= Nx; i++)
	{
	  double SC0 = 0.5  * (  exp (- xx(i-1)*xx(i-1) /2./sigma/sigma) /sqrt(2.*M_PI) /sigma - ne0t(i-1, j)
			       + exp (- xx(i  )*xx(i  ) /2./sigma/sigma) /sqrt(2.*M_PI) /sigma - ne0t(i,   j));
	  double SC2 = 0.5  * ( - ne2t(i-1, j) - ne2t(i, j));
	  double SE0 = 0.25 * ( - Te0t(i-1, j) - Te0t(i, j));
	  double SE2 = 0.25 * (  exp (- xx(i-1)*xx(i-1) /2./sigma/sigma) /sqrt(2.*M_PI) /sigma - Te2t(i-1, j)
			       + exp (- xx(i  )*xx(i  ) /2./sigma/sigma) /sqrt(2.*M_PI) /sigma - Te2t(i,   j));

	  Ve0(i, j) = Ve0(i-1, j) + h * SC0;
	  Qe0(i, j) = Qe0(i-1, j) + h * SE0;
	  qe0(i, j) = Qe0(i,   j) - Ve0(i, j);
	  Ve2(i, j) = Ve2(i-1, j) + h * SC2;
	  Qe2(i, j) = Qe2(i-1, j) + h * SE2;
	  qe2(i, j) = Qe2(i,   j) - Ve2(i, j);
	}
    }

  // ........................
  // Calculate central values
  // ........................
  for (int j = 0; j <= Nt; j++)
    {
      ne00(j) = ne0(0, j);
      Te00(j) = Te0(0, j);
      ne20(j) = ne2(0, j);
      Te20(j) = Te2(0, j);
    }
  
  // .............................
  // Calculate spatial half-widths
  // .............................
  Lne0(0) = 0.;
  LTe0(0) = 0.;
  Lne2(0) = 0.;
  LTe2(0) = 0.;
  for (int j = 1; j <= Nt; j++)
    {
      int I = 0; double xm1, x0;
      for (int i = 1; i <= Nx; i++)
	if (ne0(i, j) /ne0(0, j) < 0.5 && I == 0)
	  I = i;

      xm1 = 2. * ne0(I-1, j) /ne0(0, j) - 1.;
      x0  = 2. * ne0(I,   j) /ne0(0, j) - 1.;
      Lne0(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);

      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (Te0(i, j) /Te0(0, j) < 0.5 && I == 0)
	  I = i;

      xm1 = 2. * Te0(I-1, j) /Te0(0, j) - 1.;
      x0  = 2. * Te0(I,   j) /Te0(0, j) - 1.;
      LTe0(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
   
      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (ne2(i, j) /ne2(0, j) < 0.5 && I == 0)
	  I = i;

      xm1 = 2. * ne2(I-1, j) /ne2(0, j) - 1.;
      x0  = 2. * ne2(I,   j) /ne2(0, j) - 1.;
      Lne2(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
        
      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (Te2(i, j) /Te2(0, j) < 0.5 && I == 0)
	  I = i;

      xm1 = 2. * Te2(I-1, j) /Te2(0, j) - 1.;
      x0  = 2. * Te2(I,   j) /Te2(0, j) - 1.;
      LTe2(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
    }

  // ............................
  // Calculate spatial 90%-widths
  // ............................
  Wne0(0) = 0.;
  WTe0(0) = 0.;
  Wne2(0) = 0.;
  WTe2(0) = 0.;
  for (int j = 1; j <= Nt; j++)
    {
      int I = 0; double xm1, x0;
      for (int i = 1; i <= Nx; i++)
	if (ne0(i, j) /ne0(0, j) < 0.1 && I == 0)
	  I = i;

      xm1 = 10. * ne0(I-1, j) /ne0(0, j) - 1.;
      x0  = 10. * ne0(I,   j) /ne0(0, j) - 1.;
      Wne0(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);

      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (Te0(i, j) /Te0(0, j) < 0.1 && I == 0)
	  I = i;

      xm1 = 10. * Te0(I-1, j) /Te0(0, j) - 1.;
      x0  = 10. * Te0(I,   j) /Te0(0, j) - 1.;
      WTe0(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
   
      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (ne2(i, j) /ne2(0, j) < 0.1 && I == 0)
	  I = i;

      xm1 = 10. * ne2(I-1, j) /ne2(0, j) - 1.;
      x0  = 10. * ne2(I,   j) /ne2(0, j) - 1.;
      Wne2(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
        
      I = 0;
      for (int i = 1; i <= Nx; i++)
	if (Te2(i, j) /Te2(0, j) < 0.1 && I == 0)
	  I = i;

      xm1 = 10. * Te2(I-1, j) /Te2(0, j) - 1.;
      x0  = 10. * Te2(I,   j) /Te2(0, j) - 1.;
      WTe2(j) = (xx(I-1) * x0 - xx(I) * xm1) / (x0 - xm1);
    }
  
  // ..................
  // Output netcdf file
  // ..................

  int err = 0, dataFile;
  err = nc_create ("Plots/Simulation.nc", NC_CLOBBER, &dataFile);
  
  if (err != 0)
    {
      printf ("Error opening Plots/Simulation.nc\n");
      exit (1);
    }

  int x_d, t_d;
  err += nc_def_dim (dataFile, "Nx", 2*Nx+1, &x_d);
  err += nc_def_dim (dataFile, "Nt", Nt+1, &t_d);

  double* xxx    = new double[2*Nx+1];
  double* ttt    = new double[Nt+1];
  double* data1  = new double[(2*Nx+1)*(Nt+1)];
  double* data2  = new double[(2*Nx+1)*(Nt+1)];
  double* data3  = new double[(2*Nx+1)*(Nt+1)];
  double* data4  = new double[(2*Nx+1)*(Nt+1)];
  double* data5  = new double[(2*Nx+1)*(Nt+1)];
  double* data6  = new double[(2*Nx+1)*(Nt+1)];
  double* data7  = new double[(2*Nx+1)*(Nt+1)];
  double* data8  = new double[(2*Nx+1)*(Nt+1)];
  double* data9  = new double[(2*Nx+1)*(Nt+1)];
  double* data10 = new double[(2*Nx+1)*(Nt+1)];
  double* data11 = new double[Nt+1];
  double* data12 = new double[Nt+1];
  double* data13 = new double[Nt+1];
  double* data14 = new double[Nt+1];
  double* data15 = new double[Nt+1];
  double* data16 = new double[Nt+1];
  double* data17 = new double[Nt+1];
  double* data18 = new double[Nt+1];
  double* data19 = new double[Nt+1];
  double* data20 = new double[Nt+1];
  double* data21 = new double[Nt+1];
  double* data22 = new double[Nt+1];

  for (int i = 0; i <= Nt; i++)
    ttt[i] = tt(i);

  for (int i = 0; i < Nx; i++)
    xxx[i] = - xx(Nx-i);
  xxx[Nx] = xx(0);
  for (int i = 1; i <= Nx; i++)
    xxx[Nx+i] = xx(i);

  for (int i = 0; i < Nx; i++)
    for (int j = 0; j <= Nt; j++)
      {
	data1 [j + i*(Nt+1)] =   ne0(Nx-i, j);
	data2 [j + i*(Nt+1)] =   Te0(Nx-i, j);
	data3 [j + i*(Nt+1)] =   ne2(Nx-i, j);
	data4 [j + i*(Nt+1)] =   Te2(Nx-i, j);
	data5 [j + i*(Nt+1)] = - Ve0(Nx-i, j);
	data6 [j + i*(Nt+1)] = - Qe0(Nx-i, j);
	data7 [j + i*(Nt+1)] = - qe0(Nx-i, j);
	data8 [j + i*(Nt+1)] = - Ve2(Nx-i, j);
	data9 [j + i*(Nt+1)] = - Qe2(Nx-i, j);
	data10[j + i*(Nt+1)] = - qe2(Nx-i, j);
      }
  for (int j = 0; j <= Nt; j++)
    {
      data1 [j + Nx*(Nt+1)] = ne0(0, j);
      data2 [j + Nx*(Nt+1)] = Te0(0, j);
      data3 [j + Nx*(Nt+1)] = ne2(0, j);
      data4 [j + Nx*(Nt+1)] = Te2(0, j);
      data5 [j + Nx*(Nt+1)] = Ve0(0, j);
      data6 [j + Nx*(Nt+1)] = Qe0(0, j);
      data7 [j + Nx*(Nt+1)] = qe0(0, j);
      data8 [j + Nx*(Nt+1)] = Ve2(0, j);
      data9 [j + Nx*(Nt+1)] = Qe2(0, j);
      data10[j + Nx*(Nt+1)] = qe2(0, j);
    }
  for (int i = 1; i <= Nx; i++)
    for (int j = 0; j <= Nt; j++)
      {
	data1 [j + (Nx+i)*(Nt+1)] = ne0(i, j);
	data2 [j + (Nx+i)*(Nt+1)] = Te0(i, j);
	data3 [j + (Nx+i)*(Nt+1)] = ne2(i, j);
	data4 [j + (Nx+i)*(Nt+1)] = Te2(i, j);
	data5 [j + (Nx+i)*(Nt+1)] = Ve0(i, j);
	data6 [j + (Nx+i)*(Nt+1)] = Qe0(i, j);
	data7 [j + (Nx+i)*(Nt+1)] = qe0(i, j);
	data8 [j + (Nx+i)*(Nt+1)] = Ve2(i, j);
	data9 [j + (Nx+i)*(Nt+1)] = Qe2(i, j);
	data10[j + (Nx+i)*(Nt+1)] = qe2(i, j);
      }
  for (int i = 0; i <= Nt; i++)
    {
      data11[i] = Lne0(i);
      data12[i] = LTe0(i);
      data13[i] = Lne2(i);
      data14[i] = LTe2(i);
      data15[i] = ne00(i);
      data16[i] = Te00(i);
      data17[i] = ne20(i);
      data18[i] = Te20(i);
      data19[i] = Wne0(i);
      data20[i] = WTe0(i);
      data21[i] = Wne2(i);
      data22[i] = WTe2(i);
    }

  int xx_d, tt_d;
  err += nc_def_var (dataFile, "x", NC_DOUBLE, 1, &x_d, &xx_d);
  err += nc_def_var (dataFile, "t", NC_DOUBLE, 1, &t_d, &tt_d);

  int data_d[2];
  data_d[0] = x_d;
  data_d[1] = t_d;

  int ne0_d, Te0_d, ne2_d, Te2_d, Ve0_d, Qe0_d, qe0_d, Ve2_d, Qe2_d, qe2_d;
  err += nc_def_var (dataFile, "ne_0", NC_DOUBLE, 2, data_d, &ne0_d);
  err += nc_def_var (dataFile, "Te_0", NC_DOUBLE, 2, data_d, &Te0_d);
  err += nc_def_var (dataFile, "ne_2", NC_DOUBLE, 2, data_d, &ne2_d);
  err += nc_def_var (dataFile, "Te_2", NC_DOUBLE, 2, data_d, &Te2_d);
  err += nc_def_var (dataFile, "Ve_0", NC_DOUBLE, 2, data_d, &Ve0_d);
  err += nc_def_var (dataFile, "Qe_0", NC_DOUBLE, 2, data_d, &Qe0_d);
  err += nc_def_var (dataFile, "qe_0", NC_DOUBLE, 2, data_d, &qe0_d);
  err += nc_def_var (dataFile, "Ve_2", NC_DOUBLE, 2, data_d, &Ve2_d);
  err += nc_def_var (dataFile, "Qe_2", NC_DOUBLE, 2, data_d, &Qe2_d);
  err += nc_def_var (dataFile, "qe_2", NC_DOUBLE, 2, data_d, &qe2_d);

  int Ln0_d, LT0_d, Ln2_d, LT2_d;
  err += nc_def_var (dataFile, "Ln_0", NC_DOUBLE, 1, &t_d, &Ln0_d);
  err += nc_def_var (dataFile, "LT_0", NC_DOUBLE, 1, &t_d, &LT0_d);
  err += nc_def_var (dataFile, "Ln_2", NC_DOUBLE, 1, &t_d, &Ln2_d);
  err += nc_def_var (dataFile, "LT_2", NC_DOUBLE, 1, &t_d, &LT2_d);

  int n00_d, T00_d, n20_d, T20_d;
  err += nc_def_var (dataFile, "n0_0", NC_DOUBLE, 1, &t_d, &n00_d);
  err += nc_def_var (dataFile, "T0_0", NC_DOUBLE, 1, &t_d, &T00_d);
  err += nc_def_var (dataFile, "n0_2", NC_DOUBLE, 1, &t_d, &n20_d);
  err += nc_def_var (dataFile, "T0_2", NC_DOUBLE, 1, &t_d, &T20_d);

  int Wn0_d, WT0_d, Wn2_d, WT2_d;
  err += nc_def_var (dataFile, "Wn_0", NC_DOUBLE, 1, &t_d, &Wn0_d);
  err += nc_def_var (dataFile, "WT_0", NC_DOUBLE, 1, &t_d, &WT0_d);
  err += nc_def_var (dataFile, "Wn_2", NC_DOUBLE, 1, &t_d, &Wn2_d);
  err += nc_def_var (dataFile, "WT_2", NC_DOUBLE, 1, &t_d, &WT2_d);

  err += nc_enddef (dataFile);

  if (err != 0)
    {
      printf ("Error defining variables in Plots/Simulation.nc\n");
      exit (1);
    }

  err += nc_put_var_double (dataFile, xx_d,  xxx);
  err += nc_put_var_double (dataFile, tt_d,  ttt);
  err += nc_put_var_double (dataFile, ne0_d, data1);
  err += nc_put_var_double (dataFile, Te0_d, data2);
  err += nc_put_var_double (dataFile, ne2_d, data3);
  err += nc_put_var_double (dataFile, Te2_d, data4);
  err += nc_put_var_double (dataFile, Ve0_d, data5);
  err += nc_put_var_double (dataFile, Qe0_d, data6);
  err += nc_put_var_double (dataFile, qe0_d, data7);
  err += nc_put_var_double (dataFile, Ve2_d, data8);
  err += nc_put_var_double (dataFile, Qe2_d, data9);
  err += nc_put_var_double (dataFile, qe2_d, data10);
  err += nc_put_var_double (dataFile, Ln0_d, data11);
  err += nc_put_var_double (dataFile, LT0_d, data12);
  err += nc_put_var_double (dataFile, Ln2_d, data13);
  err += nc_put_var_double (dataFile, LT2_d, data14);
  err += nc_put_var_double (dataFile, n00_d, data15);
  err += nc_put_var_double (dataFile, T00_d, data16);
  err += nc_put_var_double (dataFile, n20_d, data17);
  err += nc_put_var_double (dataFile, T20_d, data18);
  err += nc_put_var_double (dataFile, Wn0_d, data19);
  err += nc_put_var_double (dataFile, WT0_d, data20);
  err += nc_put_var_double (dataFile, Wn2_d, data21);
  err += nc_put_var_double (dataFile, WT2_d, data22);

  if (err != 0)
    {
      printf ("Error writing Plots/Simulation.nc\n");
      exit (1);
    }
  
   err += nc_close (dataFile);

  if (err != 0)
    {
      printf ("Error closing Plots/Simulation.nc\n");
      exit (1);
    }
  
  // ........
  // Clean up
  // ........
  delete[] xxx;    delete[] ttt;
  delete[] data1;  delete[] data2;  delete[] data3;  delete[] data4;
  delete[] data5;  delete[] data6;  delete[] data7;  delete[] data8;
  delete[] data9;  delete[] data10; delete[] data11; delete[] data12;
  delete[] data13; delete[] data14; delete[] data15; delete[] data16;
  delete[] data17; delete[] data18; delete[] data19; delete[] data20;
  delete[] data21; delete[] data22;
  
  gsl_spline_free (spline_n0r); gsl_spline_free (spline_n0i);
  gsl_spline_free (spline_T0r); gsl_spline_free (spline_T0i);
  gsl_spline_free (spline_n2r); gsl_spline_free (spline_n2i);
  gsl_spline_free (spline_T2r); gsl_spline_free (spline_T2i);
  
  gsl_interp_accel_free (acc_n0r); gsl_interp_accel_free (acc_n0i);
  gsl_interp_accel_free (acc_T0r); gsl_interp_accel_free (acc_T0i);
  gsl_interp_accel_free (acc_n2r); gsl_interp_accel_free (acc_n2i);
  gsl_interp_accel_free (acc_T2r); gsl_interp_accel_free (acc_T2i);
}

// ###############################################################
// Function to evaluate inverse Laplace transform target functions
// ###############################################################
void Parallel::CalcF ()
{
  // .................
  // Output integrands
  // .................
  FILE* file1 = OpenFilew ("Plots/Kn0.out");
  FILE* file2 = OpenFilew ("Plots/KT0.out");
  FILE* file3 = OpenFilew ("Plots/Kn2.out");
  FILE* file4 = OpenFilew ("Plots/KT2.out");

  for (int i = 1; i <= Nx; i++)
    {
      double k     = kk(i);
      double ks    = k * sigma;
      double gauss = exp (- ks*ks /2.) /M_PI;

      fprintf (file1, "%11.4e ", k);
      fprintf (file2, "%11.4e ", k);
      fprintf (file3, "%11.4e ", k);
      fprintf (file4, "%11.4e ", k);

      for (int j = 0; j <= Ng; j = j + Ng/10)
	{
	  complex<double> g  = complex<double> (gamma, gg(j));

	  complex<double> kn0 = GetKn0 (g, k) * gauss;
	  complex<double> kT0 = GetKT0 (g, k) * gauss;
	  complex<double> kn2 = GetKn2 (g, k) * gauss;
	  complex<double> kT2 = GetKT2 (g, k) * gauss;

	  fprintf (file1, "%11.4e %11.4e ", real (kn0), imag (kn0));
	  fprintf (file2, "%11.4e %11.4e ", real (kT0), imag (kT0));
	  fprintf (file3, "%11.4e %11.4e ", real (kn2), imag (kn2));
	  fprintf (file4, "%11.4e %11.4e ", real (kT2), imag (kT2));
	}
      fprintf (file1, "\n"); fprintf (file2, "\n"); fprintf (file3, "\n"); fprintf (file4, "\n"); 
    }

  fclose (file1); fclose (file2); fclose (file3); fclose (file4);

  // ...................
  // Set Simpson weights
  // ...................
  double* Weight = new double[Nx+1];
  double  h      = ksmax /sigma /double (Nx);
  for (int j = 0; j <= Nx; j++)
    {
      Weight[j] = 0.;
    }
  for (int j = 0; j < Nx-1; j += 2)
    {
      Weight[j]   +=      h /3.;
      Weight[j+1] += 4. * h /3.;
      Weight[j+2] +=      h /3.;
    }

  // ..............................................
  // Calculate Laplace transformed target functions
  // ..............................................
  printf ("\nCalculating inverse Laplace transformed target functions:\n");
  for (int i = 0; i <= Nx; i++)
    {
      if (i%10 == 0)
	printf (".");
      fflush (stdout);
      
      double x = xx(i);

      for (int j = 0; j <= Ng; j++)
	{
	  complex<double> g  = complex<double> (gamma, gg(j));

	  complex<double> fn0 = complex<double> (0., 0.);
	  complex<double> fT0 = complex<double> (0., 0.);
	  complex<double> fn2 = complex<double> (0., 0.);
	  complex<double> fT2 = complex<double> (0., 0.);

	  for (int l = 1; l <= Nx; l++)
	    {
	      double k     = kk (l);
	      double ks    = k * sigma;
	      double gauss = exp (- ks*ks /2.) /M_PI;
	      double ckx   = gauss * cos (k * x);

	      if (l == 1)
		{
		  fn0 += (Weight[0] /M_PI + Weight[1] * ckx) * GetKn0 (g, k);
		  fT0 += (Weight[0] /M_PI + Weight[1] * ckx) * GetKT0 (g, k);
		  fn2 += (Weight[0] /M_PI + Weight[1] * ckx) * GetKn2 (g, k);
		  fT2 += (Weight[0] /M_PI + Weight[1] * ckx) * GetKT2 (g, k);
		}
	      else
		{
		  fn0 += Weight[l] * GetKn0 (g, k) * ckx;
		  fT0 += Weight[l] * GetKT0 (g, k) * ckx;
		  fn2 += Weight[l] * GetKn2 (g, k) * ckx;
		  fT2 += Weight[l] * GetKT2 (g, k) * ckx;
		}
	    }

	  Fn0 (i, j) = fn0;
	  FT0 (i, j) = fT0;
	  Fn2 (i, j) = fn2;
	  FT2 (i, j) = fT2;
	}
    }
  printf ("\n");
  delete[] Weight;

  // ...........................................
  // Output Laplace transformed target functions
  // ...........................................
  file1 = OpenFilew ("Plots/Fn0.out");
  file2 = OpenFilew ("Plots/FT0.out");
  file3 = OpenFilew ("Plots/Fn2.out");
  file4 = OpenFilew ("Plots/FT2.out");

  for (int j = 0; j <= Ng; j++)
    {
      double g_i = gg(j);

      fprintf (file1, "%11.4e ", g_i);
      fprintf (file2, "%11.4e ", g_i);
      fprintf (file3, "%11.4e ", g_i);
      fprintf (file4, "%11.4e ", g_i);

      for (int i = 0; i <= Nx; i = i + Nx/10)
	{
	  fprintf (file1, "%11.4e %11.4e ", real (Fn0 (i, j)), imag (Fn0 (i, j)));
	  fprintf (file2, "%11.4e %11.4e ", real (FT0 (i, j)), imag (FT0 (i, j)));
	  fprintf (file3, "%11.4e %11.4e ", real (Fn2 (i, j)), imag (Fn2 (i, j)));
	  fprintf (file4, "%11.4e %11.4e ", real (FT2 (i, j)), imag (FT2 (i, j)));
	}
      fprintf (file1, "\n"); fprintf (file2, "\n"); fprintf (file3, "\n"); fprintf (file4, "\n"); 
    }

  fclose (file1); fclose (file2); fclose (file3); fclose (file4);
}

// #######################
// Function to evaluate Z0
// #######################
complex<double> Parallel::GetZ0 (complex<double> xi)
{
  complex<double> Z0;
  complex<double> fad = w (xi);
  complex<double> Z   = complex<double> (0., 1.) * sqrt (M_PI) * fad;

  if (imag (xi) < 0.)
    Z -= 2. * complex<double> (0., 1.) * sqrt (M_PI) * exp (-xi*xi);
  
  Z0 = - xi * Z;
 
  if (isnan (real (Z0)) || isnan (imag (Z0)) || isinf (real (Z0)) || isinf (imag (Z0)))
    {
      printf ("GetZ0: Error: xi = (%11.4e, %11.4e)  w = (%11.4e, %11.4e)  Z0 = (%11.4e, %11.4e)\n",
	      real(xi), imag(xi), real(fad), imag(fad), real(Z0), imag(Z0));
      exit (1); 
    }
  else
    return Z0;
}

// ######################
// Function to evaluate G
// ######################
complex<double> Parallel::GetG (complex<double> xi)
{
  complex<double> G;

  if (abs (xi) > 100.)
    {
      G = 3./4./xi/xi + 3./2./xi/xi/xi/xi;

      if (isnan (real (G)) || isnan (imag (G)) || isinf (real (G)) || isinf (imag (G)))
	{
	  printf ("GetG: Error: xi = (%11.4e, %11.4e) G = (%11.4e, %11.4e)\n",
		  real(xi), imag(xi), real(G), imag(G));
	  exit (1);
	}
      else
	return G;
    }
  else
    {
      complex<double> Z0  = GetZ0 (xi);
      complex<double> top = (xi*xi - 1.) - (xi*xi - 1.5) * Z0;
      complex<double> bot = 2. * xi*xi - (2.*xi*xi - 1.) * Z0;
      
      G = top /bot;
      
      if (isnan (real (G)) || isnan (imag (G)) || isinf (real (G)) || isinf (imag (G)))
	{
	  printf ("GetG: Error: xi = (%11.4e, %11.4e) Z0 = (%11.4e, %11.4e) top = (%11.4e, %11.4e) bot = (%11.4e, %11.4e) G = (%11.4e, %11.4e)\n",
		  real(xi), imag(xi), real(Z0), imag(Z0), real(top), imag(top), real(bot), imag(bot), real(G), imag(G));
	  exit (1); 
	}
      else
	return G;
    }
}

// ############################
// Function to evaluate LambdaD
// ############################
double Parallel::GetLambdaD (double k)
{
  double kl  = k * lambdaD;
  double kl2 = kl*kl;

  return (1. + kl2) /kl2;
}

// ########################
// Function to evaluate Kn0
// ########################
complex<double> Parallel::GetKn0 (complex<double> g, double k)
{
  complex<double> Kn0;
  complex<double> xi      = complex<double> (0., 1.) * (1. + g) /k;
  complex<double> g1      = g /(1. + g);
  complex<double> g2      = 2. * (1. - mue + g) * xi*xi /(1. + g);
  complex<double> G       = GetG (xi);
  double          LambdaD = GetLambdaD (k);
  complex<double> top     = - 1. - g2 * (G - g1/2.);
  complex<double> bot     = g * (1. + g) * ((G - g1/2.) * (LambdaD - g1 * g2) - g1);

  Kn0 =  top /bot;

  if (isnan (real (Kn0)) || isnan (imag (Kn0)))
    {
      printf ("GetKn0: Error: xi = (%11.4e, %11.4e) G = (%11.4e, %11.4e) top = (%11.4e, %11.4e) bot = (%11.4e, %11.4e) Kn0 = (%11.4e, %11.4e)\n",
	      real(xi), imag(xi), real(G), imag(G), real(top), imag(top), real(bot), imag(bot), real(Kn0), imag(Kn0));
      exit (1); 
    }
  else
    return Kn0;
}

// ########################
// Function to evaluate KT0
// ########################
complex<double> Parallel::GetKT0 (complex<double> g, double k)
{
  complex<double> KT0;
  complex<double> xi      = complex<double> (0., 1.) * (1. + g) /k;
  complex<double> g1      = g /(1. + g);
  complex<double> g2      = 2. * (1. - mue + g) * xi*xi /(1. + g);
  complex<double> G       = GetG (xi);
  double          LambdaD = GetLambdaD (k);
  complex<double> top     = LambdaD;
  complex<double> bot     = g * (1. + g) * ((G - g1/2.) * (LambdaD - g1 * g2) - g1);

  KT0 =  top /bot;

  if (isnan (real (KT0)) || isnan (imag (KT0)))
    {
      printf ("GetKT0: Error: xi = (%11.4e, %11.4e)  KT0 = (%11.4e, %11.4e)\n",
	      real(xi), imag(xi), real(KT0), imag(KT0));
      exit (1); 
    }
  else
    return KT0;
}

// ########################
// Function to evaluate Kn2
// ########################
complex<double> Parallel::GetKn2 (complex<double> g, double k)
{
  complex<double> Kn2;
  complex<double> xi      = complex<double> (0., 1.) * (1. + g) /k;
  complex<double> g1      = g /(1. + g);
  complex<double> g2      = 2. * (1. - mue + g) * xi*xi /(1. + g);
  complex<double> G       = GetG (xi);
  double          LambdaD = GetLambdaD (k);
  complex<double> top     = 1.;
  complex<double> bot     = 2. * g * (1. + g) * ((G - g1/2.) * (LambdaD - g1 * g2) - g1);

  Kn2 = top /bot;

  if (isnan (real (Kn2)) || isnan (imag (Kn2)))
    {
      printf ("GetKn2: Error: xi = (%11.4e, %11.4e)  Kn2 = (%11.4e, %11.4e)\n",
	      real(xi), imag(xi), real(Kn2), imag(Kn2));
      exit (1); 
    }
  else
    return Kn2;
}

// ########################
// Function to evaluate KT2
// ########################
complex<double> Parallel::GetKT2 (complex<double> g, double k)
{
  complex<double> KT2;
  complex<double> xi      = complex<double> (0., 1.) * (1. + g) /k;
  complex<double> g1      = g /(1. + g);
  complex<double> g2      = 2. * (1. - mue + g) * xi*xi /(1. + g);
  complex<double> G       = GetG (xi);
  double          LambdaD = GetLambdaD (k);
  complex<double> top     = - LambdaD + g1 * g2;
  complex<double> bot     = 2. * g * (1. + g) * ((G - g1/2.) * (LambdaD - g1 * g2) - g1);

  KT2 = top /bot;

  if (isnan (real (KT2)) || isnan (imag (KT2)))
    {
      printf ("GetKn2: Error: xi = (%11.4e, %11.4e)  KT2 = (%11.4e, %11.4e)\n",
	      real(xi), imag(xi), real(KT2), imag(KT2));
      exit (1); 
    }
  else
    return KT2;
}

// ###############################################################
// Function to evaluate right-hand sides of differential equations
// ###############################################################
void Parallel::Rhs (double x, Array<double,1> y, Array<double,1> dydx)
{
  double egt = exp (gamma * T);

  dydx(0) = egt * (gsl_spline_eval (spline_n0r, x, acc_n0r) * cos (x * T) - gsl_spline_eval (spline_n0i, x, acc_n0i) * sin (x * T)) /M_PI;
  dydx(1) = egt * (gsl_spline_eval (spline_T0r, x, acc_T0r) * cos (x * T) - gsl_spline_eval (spline_T0i, x, acc_T0i) * sin (x * T)) /M_PI;
  dydx(2) = egt * (gsl_spline_eval (spline_n2r, x, acc_n2r) * cos (x * T) - gsl_spline_eval (spline_n2i, x, acc_n2i) * sin (x * T)) /M_PI;
  dydx(3) = egt * (gsl_spline_eval (spline_T2r, x, acc_T2r) * cos (x * T) - gsl_spline_eval (spline_T2i, x, acc_T2i) * sin (x * T)) /M_PI;
}

// #######################################################################
//  Function to advance set of coupled first-order o.d.e.s by single step
//  using adaptive step length fourth-order/fifth-order Runge-Kutta scheme
//
//     x       ... independent variable
//     y       ... array of dependent variables
//     h       ... step length
//     t_err   ... actual truncation error per step 
//     acc     ... desired truncation error per step
//     S       ... safety factor
//     T       ... step length cannot change by more than this factor from step to step
//     rept    ... number of step recalculations		  
//     maxrept ... maximum allowable number of step recalculations		  
//     h_min   ... minimum allowable step length
//     h_max   ... maximum allowable step length
//     flag    ... controls manner in which truncation error is calculated	
//
//  Function advances equations by single step while attempting to maintain 
//  constant truncation error per step of acc:
//
//    flag = 0 ... error is absolute
//    flag = 1 ... error is relative
//    flag = 2 ... error is mixed
//
// #######################################################################
void Parallel::RK4RK5Adaptive (double& x, Array<double,1> y, double& h, 
			       double& t_err, double acc, double S, double T, int& rept,
			       int maxrept, double h_min, double h_max, int flag, 
			       int diag, FILE* file)
{
  int neqns = y.extent(0);

  Array<double,1> y0  (neqns);
  Array<double,1> Err (neqns);
  double          hin = h;

  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0(i) = y(i);

  // Take RK4/RK5 step 
  RK4RK5Fixed (x, y, Err, h);

  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
        {
          err   = fabs (Err(i) /y(i));
          t_err = (err > t_err) ? err : t_err;
        }
    }
  else 
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
        {
          err1  = fabs (Err(i) /y(i));
	  err2  = fabs (Err(i));
          err   = (err1 < err2) ? err1 : err2;
          t_err = (err > t_err) ? err  : t_err;
        }
    }

  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;

  // Calculate new step length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow (fabs (acc /t_err), 0.20);
  else
    h_est = S * h * pow (fabs (acc /t_err), 0.25);

  // Prevent step length from changing by more than factor T
  if (h_est /h > T)
    h *= T;
  else if (h_est /h < 1./T)
    h /= T;
  else
    h = h_est;

  // Prevent step length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h /fabs(h) : h;

  // Prevent step length from falling below h_min
  if (fabs(h) < h_min)
    { 
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }

  // Diagnose step
  if (diag) 
    fprintf (file, "x = %11.4e hin = %11.4e err = %11.4e acc = %11.4e hout = %11.4e count = %3d\n", 
	     x, hin, t_err, acc, h, count);

  // Check if truncation error acceptable
  if ((t_err <= acc) || (count >= maxrept))
    {
      // If truncation error acceptable take step 
      rept  = count;
      count = 0;
    }
  else 
    {
      // If truncation error unacceptable repeat step 
      count++;
      x = x0;
      for (int i = 0; i < neqns; i++)
	y(i) = y0(i);
      RK4RK5Adaptive (x, y, h, t_err, acc, S, T, rept, 
		      maxrept, h_min, h_max, flag, diag, file);
    }
}

// #####################################################################
// Function to advance set of coupled first-order o.d.e.s by single step
// using fixed step length fourth-order/fifth-order Runge-Kutta scheme
//
//     x       ... independent variable
//     y       ... array of dependent variables 
//     err     ... array of errors
//     h       ... step length
//     
// #####################################################################
void Parallel::RK4RK5Fixed (double& x, Array<double,1> y, Array<double,1> err, double h)
{
  int neqns = y.extent(0);
  
  Array<double,1> dydx (neqns);
  Array<double,1> k1   (neqns);
  Array<double,1> k2   (neqns);
  Array<double,1> k3   (neqns);
  Array<double,1> k4   (neqns);
  Array<double,1> k5   (neqns);
  Array<double,1> k6   (neqns);
  Array<double,1> f    (neqns);

  // First stage
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1(i) = h * dydx(i);
      f (i) = y(i) + bb21 * k1(i);
    }

  // Second stage
  Rhs (x + aa2 * h, f, dydx);  
  for (int i = 0; i < neqns; i++)
    {
      k2(i) = h * dydx(i);
      f (i) = y(i) + bb31 * k1(i) + bb32 * k2(i);
    }

  // Third stage
  Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3(i) = h * dydx(i);
      f (i) = y(i) + bb41 * k1(i) + bb42 * k2(i) + bb43 * k3(i);
    }

  // Fourth stage
  Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4(i) = h * dydx(i);
      f (i) = y(i) + bb51 * k1(i) + bb52 * k2(i) + bb53 * k3(i) + bb54 * k4(i);
    }

  // Fifth stage
  Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5(i) = h * dydx(i);
      f (i) = y(i) + bb61 * k1(i) + bb62 * k2(i) + bb63 * k3(i) + bb64 * k4(i) + bb65 * k5(i);
    }

  // Sixth stage
  Rhs (x + aa6 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k6(i) = h * dydx(i);
    }

  // Actual step 
  for (int i = 0; i < neqns; i++)
    {
      y  (i) = y(i) + cc1 * k1(i) + cc3 * k3(i) + cc4 * k4(i)               + cc6 * k6(i);
      err(i) =        ca1 * k1(i) + ca3 * k3(i) + ca4 * k4(i) + ca5 * k5(i) + ca6 * k6(i);
    }
  x += h;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Parallel::OpenFilew (char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL) 
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}

