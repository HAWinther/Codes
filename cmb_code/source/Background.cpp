
//=================================================================
//	General methods for computing the background evolution
//=================================================================

#include "Background.h"

// Calculate the conformal time
void Background::calc_eta(){

#ifdef VERBOSE
  cout << endl;
  cout <<  "**********************************************" << endl;
  cout <<	 "*** Calculate background expansion history ***" << endl;
  cout <<	 "**********************************************" << endl;
#endif

  // Integration variable
  VecDoub y(1);

  // Initial value. Approx value assuming LCDM in the past with Omega_relativistic = 8.3e-5
  eta[0] = 109.76*exp(x_main[0]);
  y[0]   = eta[0];

  double atol = 1.0e-10;
  double rtol = 1.0e-10;
  double h1 = 1.0e-5;
  double hmin = 0.0;

  // Set up ODE
  Output out;
  derivs_eta deta(this);
  Odeint<StepperBS<derivs_eta> > ode_eta(y,x_main[0],x_main[1],atol,rtol,h1,hmin,out,deta);

#ifdef VERBOSE
  cout << endl << "---> Calculate conformal time" << endl;
#endif

  // Integrate
  for (int i = 1; i<nmain; i++) {
    ode_eta.x1 = x_main[i-1];
    ode_eta.x2 = x_main[i];
    ode_eta.integrate();
    eta[i] = y[0];
  }

  // Spline eta
  eta_spline = Spline(x_main, eta, 1.0e30, 1.0e30);
  //eta_spline.set_linear_x();
  eta_spline.set_three_regions_x(np1, np2);

#ifdef OUTPUT_BACKGROUND

#ifdef VERBOSE
  cout << endl;
  cout << "---> Output Omega_i(x)" << endl;
  cout <<	"     Filename: '" << cosm_simu->outprefix << "background_omega.dat'" << endl;
#endif

  ofstream omega_out;
  string omega_filename = cosm_simu->outprefix + "background_omega.dat";
  omega_out.open(omega_filename.c_str());
  double x, H2;
  for (int i = 0; i<nmain; i++) {
    x = x_main[i];
    H2 = H(x)*H(x);
    omega_out << x << " " << cosm_simu->Omega_m*exp(-3.0*x)/H2 << " " << cosm_simu->Omega_b*exp(-3.0*x)/H2 << " " << cosm_simu->Omega_r*exp(-4.0*x)/H2 << " " << cosm_simu->Omega_nu*exp(-4.0*x)/H2 << " " << cosm_simu->Omega_lambda/H2 << " " << get_eta(x) << endl;
  }
#endif

}

double Background::get_eta(double x){
  return eta_spline.get_spline(x);
}

void Background::get_all_H_c(const double x, double &Hc, double &dHHc, double &ddHHc){
  double Hx   = H(x);
  double dHx  = dH(x);
  double ddHx = ddH(x);

  Hc    = Hx*exp(x);
  dHHc  = 1.0 + dHx/Hx;
  ddHHc = (1.0 + 2.0*dHx/Hx + ddHx/Hx);
}
