
//==================================================
// General methods to compute the recombination
// history of the universe.
//==================================================

#include "Recombination.h"

// Calculate the full recombiantion history
void Recombination::calc_rec_history(){	

#ifdef VERBOSE
  cout  << endl;
  cout  << "**********************************************" << endl;
  cout	<< "***    Calculate recombination history     ***" << endl;
  cout  << "**********************************************" << endl;
#endif

  // Start time
  cosm_simu->time_now = clock();

  // Calculate electron number density n_e(x)
#ifdef VERBOSE
  cout  << endl	<< "---> Calculate electron density" << endl;
  calc_n_e();

  cout  << endl	<< "---> Calculate optical debth" << endl;
  calc_tau();

  cout  << endl	<< "---> Calculate visibility function" << endl;
  calc_g();
#else
  calc_n_e();
  calc_tau();
  calc_g();
#endif

#ifdef OUTPUT_RECOMBINATION

#ifdef VERBOSE
  cout << endl;
  cout << "---> Output tau(x), tau'(x) and tau''(x)" << endl; 
  cout <<	"     Filename: '" << cosm_simu->outprefix << "recombination_tau.dat'" << endl;
#endif

  ofstream out_tau;
  string out_tau_filename = cosm_simu->outprefix + "recombination_tau.dat";
  out_tau.open(out_tau_filename.c_str());
  double xx;
  int nn = 10000;
  for (int i = 0; i<=nn; i++) {
    xx = x_start_rec*double(i)/double(nn);
    out_tau << xx << " " << get_tau(xx) << " " << get_dtau(xx) << " " << get_ddtau(xx) << endl;
  }
  out_tau.close();

#ifdef VERBOSE
  cout << endl;
  cout << "---> Output g(x), g'(x) and g''(x)" << endl;
  cout <<	"     Filename: '" << cosm_simu->outprefix << "recombination_g.dat'" << endl;
#endif

  ofstream out_g;
  string out_g_filename = cosm_simu->outprefix + "recombination_g.dat";
  out_g.open(out_g_filename.c_str());
  for (int i = 0; i<=nn; i++) {
    xx = x_start_rec*double(i)/double(nn);
    out_g << xx << " " << get_g(xx) << " " << get_dg(xx)/10.0 << " " << get_ddg(xx)/100.0 << endl;
  }
  out_g.close();

#ifdef VERBOSE
  cout << endl;
  cout << "---> Output X_e" << endl;
  cout <<	"     Filename: '" << cosm_simu->outprefix << "recombination_xe.dat'" << endl;
#endif

  ofstream out_xe;
  string out_xe_filename = cosm_simu->outprefix + "recombination_xe.dat";
  out_xe.open(out_xe_filename.c_str());
  for (int i = 0; i<=nn; i++) {
    xx = x_start_rec*double(i)/double(nn);
    out_xe << xx << " " << get_x_e(xx) << endl;
  }
  out_xe.close();
#endif

  // End time
  cosm_simu->time_now = clock() - cosm_simu->time_now;

#ifdef VERBOSE
  cout << endl;
  cout << "---> Calculated values" << endl;
  cout <<	"     Electron density today     : " << get_n_e(0.0) << endl;
  cout <<	"     Optical debth today        : " << get_tau(0.0) << endl;
  if (cosm_simu->z_reion>0.0)
    cout <<	"     Optical debth at reion     : " << get_tau(-log(1.0+cosm_simu->z_reion + 5.0*cosm_simu->dz_reion)) << endl;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now/1000000.0 << endl;
#endif
}

// Compute the visibility function g(x)
void Recombination::calc_g(){
  double x;

  // Initialize
  g  = VecDoub(nmain);
  g2 = VecDoub(nmain);

  for (int i = 0; i<nmain; i++) {
    x = BGSUB->x_main[i];
    g[i] = -get_dtau(x)*exp(-get_tau(x));
  }

  // Spline g
  g_spline = Spline(BGSUB->x_main, g, 1.0e30, 1.0e30);
  //g_spline.set_linear_x();
  g_spline.set_three_regions_x(np1, np2);

  // Spline g''
  g_spline.extract_y2(g2);
  g2_spline = Spline(BGSUB->x_main, g2, 1.0e30, 1.0e30);
  //g2_spline.set_linear_x();
  g2_spline.set_three_regions_x(np1, np2);
}

// Calculate the optical debth tau(x)
void Recombination::calc_tau(){

  // Integration variables
  const double rtol = 1.0e-10;
  const double atol = 1.0e-10;
  const double h1   = 1.0e-5;

  // Intialize
  tau  = VecDoub(nmain);
  tau2 = VecDoub(nmain);

  // Set up the tau equation
  derivs_tau dtau(this);
  Output out;
  VecDoub y(1);
  y[0] = 0.0;

  // Integrate tau equation
  Odeint<StepperBS<derivs_tau> > ode_tau(y,BGSUB->x_main[nmain-1],BGSUB->x_main[nmain-2],atol,rtol,h1,0.0,out,dtau);
  for (int i = nmain-1; i>=1; i--) {
    ode_tau.x1 = BGSUB->x_main[i];
    ode_tau.x2 = BGSUB->x_main[i-1];
    ode_tau.integrate();
    tau[i-1] = y[0];
  }
  tau[nmain-1] = 0.9*tau[nmain-2]; // Change to allow for spline of log(tau). Do not affect physics

  // Change tau to log(tau)
  for (int i = 0; i<nmain; i++) {
    tau[i] = log(abs(tau[i]));
  }

  // Spline tau
  tau_spline = Spline(BGSUB->x_main, tau, 1.0e30, 1.0e30);
  //tau_spline.set_linear_x();
  tau_spline.set_three_regions_x(np1, np2);

  // Spline tau''
  tau_spline.extract_y2(tau2);
  tau2_spline = Spline(BGSUB->x_main, tau2, 1.0e30, 1.0e30);
  //tau2_spline.set_linear_x();
  tau2_spline.set_three_regions_x(np1, np2);
}

// The electron density including the effect of reionization
double Recombination::get_n_e(double x){
  double reion_fact = 0.0;
  double n_b		  = 0.0;

  if (cosm_simu->z_reion > 0.0) {
    reion_fact  = 0.5 * (1.0 + tanh((pow(1.0+cosm_simu->z_reion,1.5)-exp(-1.5*x))/1.5/sqrt(1.0+cosm_simu->z_reion)/cosm_simu->dz_reion));
    n_b			= (1.0-cosm_simu->Y_p)*(cosm_simu->Omega_b*cosm_simu->rhoc)*exp(-3.0*x)/m_H;
  }
  return exp(n_e_spline.get_spline(x))*(1.0-reion_fact) + reion_fact*n_b;
}

// The electron fraction including the effect of reionization
double Recombination::get_x_e(double x){
  double reion_fact = 0.0;
  double n_b	      = (1.0-cosm_simu->Y_p)*(cosm_simu->Omega_b*cosm_simu->rhoc)*exp(-3.0*x)/m_H;

  if (cosm_simu->z_reion > 0.0) {
    reion_fact  = 0.5 * (1.0 + tanh((pow(1.0+cosm_simu->z_reion,1.5)-exp(-1.5*x))/1.5/sqrt(1.0+cosm_simu->z_reion)/cosm_simu->dz_reion));
    n_b			= (1.0-cosm_simu->Y_p)*(cosm_simu->Omega_b*cosm_simu->rhoc)*exp(-3.0*x)/m_H;
  }
  return exp(n_e_spline.get_spline(x))*(1.0-reion_fact)/n_b + reion_fact;
}

// The optical debth tau(x)
double Recombination::get_tau(double x){
  return exp( tau_spline.get_spline(x) );
}

// The derivative of the optical debth dtau(x)
double Recombination::get_dtau(double x){
  return -c*sigma_T*get_n_e(x)/(BGSUB->H(x)*cosm_simu->H_0);
}

// The second derivative of the optical debth dtau(x)
double Recombination::get_ddtau(double x){
  double dTau = get_dtau(x);
  double Tau  = get_tau(x);
  return Tau < 1.0e-10 ? 0.0 : Tau*(tau2_spline.get_spline(x)+dTau*dTau/(Tau*Tau));
}

// Get both tau, tau', tau''
void Recombination::get_all_tau(double x, double &Tau, double &dTau, double &ddTau){
  Tau   = get_tau(x);
  dTau  = get_dtau(x),
        ddTau = Tau*(tau2_spline.get_spline(x)+dTau*dTau/(Tau*Tau));
}

// The visibility function g(x)
double Recombination::get_g(double x){
  return -get_dtau(x)*exp(-get_tau(x));
}

// The derivative of the visibility function g(x)
double Recombination::get_dg(double x){
  double dTau = get_dtau(x);
  return get_g(x)*(get_ddtau(x)/dTau - dTau);
}

// The second derivative of the visibility function g(x)
double Recombination::get_ddg(double x){
  return g2_spline.get_spline(x);
}

// Get both g, dg, ddg
void Recombination::get_all_g(double x, double &G, double &dG, double &ddG){
  double Tau, dTau, ddTau;
  get_all_tau(x, Tau, dTau, ddTau);
  G     = -dTau*exp(-Tau);
  dG    = G*(ddTau/dTau - dTau);
  ddG   = g2_spline.get_spline(x);
}

// Get both tau, dtau, ddtau, g, dg, ddg
void Recombination::get_all_tau_g(double x, double &Tau, double &dTau, double &ddTau, double &G, double &dG, double &ddG){
  Tau   = get_tau(x);
  dTau  = get_dtau(x),
        ddTau = Tau*(tau2_spline.get_spline(x)+dTau*dTau/(Tau*Tau));
  G     = -dTau*exp(-Tau);
  dG    = G*(ddTau/dTau - dTau);
  ddG   = g2_spline.get_spline(x);
}
