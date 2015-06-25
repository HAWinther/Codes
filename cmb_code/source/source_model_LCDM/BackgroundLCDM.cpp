
//==================================================================
// Implementation of the standard LCDM background
//==================================================================

#include "BackgroundLCDM.h"

// Constructor of the background object
BackgroundLCDM::BackgroundLCDM(ParamSimu *param){
  cosm_simu = param;

  eta    = VecDoub(nmain);
  eta2   = VecDoub(nmain);
  x_main = VecDoub(nmain);

  // Initialize eta arrays
  for (int i = 0; i<nmain; i++) {
    eta[i]  = 0.0;
    eta2[i] = 0.0;
  }

  // Set up x-array
  if (cosm_simu->z_reion > 0.0){
    // Add extra points close to reionization period
    for (int i = 0; i<nmain; i++) {
      if (i < np1) {
        x_main[i] = x_start_main + (cosm_simu->x_start_reion - x_start_main)*i/double(np1);
      } else if (i >= np1 && i <= np1 + np2) {
        x_main[i] = cosm_simu->x_start_reion + (cosm_simu->x_end_reion - cosm_simu->x_start_reion)*(i-np1)/double(np2);
      } else {
        x_main[i] = cosm_simu->x_end_reion + (0.0 - cosm_simu->x_end_reion)*(i-np1-np2)/double(nmain-1-np1-np2);
      }
    }
  } else {
    for (int i = 0; i<nmain; i++) {
      x_main[i] = x_start_main + (0.0 - x_start_main)*double(i)/double(nmain-1);
    }
  }
}

// Compute and spline quantities related to background evolution
void BackgroundLCDM::calc_expansion_history(){
  cosm_simu->time_now = clock();

  // Calculate conformal time eta(x)
  calc_eta();

#ifdef VERBOSE
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout << endl << "---> Calculated values" << endl;
  cout <<			"     Conformal time today       : " << get_eta(0.0) << endl;
  cout <<			"     Time used (sec)            : " << cosm_simu->time_now/1000000.0 << endl;
#endif
}

//==================================================================
// Implementation of the virtual background functions in the
// background object
//==================================================================

double BackgroundLCDM::H(double x){
  return sqrt((cosm_simu->Omega_m+cosm_simu->Omega_b)*exp(-3.0*x) + (cosm_simu->Omega_r+cosm_simu->Omega_nu)*exp(-4.0*x) + cosm_simu->Omega_lambda);
}

double BackgroundLCDM::dH(double x){
  return -(3.0*(cosm_simu->Omega_m+cosm_simu->Omega_b)*exp(-3.0*x) + 4.0*(cosm_simu->Omega_r+cosm_simu->Omega_nu)*exp(-4.0*x))/(2.0*H(x));
}

double BackgroundLCDM::ddH(double x){
  double Hx  = H(x);
  double dHx = dH(x);
  return (9.0*(cosm_simu->Omega_m+cosm_simu->Omega_b)*exp(-3.0*x) + 16.0*(cosm_simu->Omega_r+cosm_simu->Omega_nu)*exp(-4.0*x))/(2.0*Hx) - dHx*dHx/Hx;
}

double BackgroundLCDM::H_c(double x){
  return exp(x)*H(x);
}

double BackgroundLCDM::dH_c(double x){
  return exp(x)*(H(x) + dH(x));
}

double BackgroundLCDM::ddH_c(double x){
  return exp(x)*(H(x) + 2.0*dH(x) + ddH(x));
}

