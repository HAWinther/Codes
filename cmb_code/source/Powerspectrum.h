
//==================================================
// General class for computing power-sepectra
//==================================================

#pragma once

#include "Parameters.h"
#include "Background.h"
#include "Recombination.h"
#include "Perturbations.h"
#include "Matrix.h"    
#include "Spline.h"
#include "Odeint.h"

class Powerspectrum {
  public:
    Background *BG;
    Recombination *REC;
    Perturbations *P;

    // Parameters
    ParamSimu *cosm_simu;

    // Bessel function arrays and splines
    MatDoub j_l, j_l2;
    VecDoub z_bessel;
    Spline bessel_spline;

    // Theta_l arrays and splines
    MatDoub Theta_l_T, Theta_l_E;
    Spline thetal_T_spline, thetal_E_spline;

    // Temperature source function arrays and splines
    MatDoub S_high_T_k;
    Spline low_res_spline_T;
    Spline high_res_spline_T;

    // E-mode source function and splines
    MatDoub S_high_E_k;
    Spline low_res_spline_E;
    Spline high_res_spline_E;

    // TT powerspectum and spline
    VecDoub cls_TT, cls_EE, cls_TE;
    Spline cl_TT_spline, cl_EE_spline, cl_TE_spline;

    // High resolution x and k arrays
    VecDoub x_high, k_high;

    // l arrays
    VecDoub lsd;

    // Matter powerspectrum array and splines
    VecDoub matterpow;
    Spline matterpow_spline;

    // Global bookkeeping variable
    int k_current;

    Powerspectrum(Background *BG, Recombination *REC, Perturbations *P, ParamSimu *param);
    ~Powerspectrum(){
      BG  = NULL;
      REC = NULL;
      P   = NULL;
      cosm_simu = NULL;
    };

    void calc_high_res_source_E();
    void calc_high_res_source_T();
    void calc_scalar_powerspectrums();
    void calc_bessel();
    void calc_matter_powerspectrum();
    void calc_sigma8();	
    void integrate_thetal_T();
    void integrate_thetal_E();
    void integrate_cl_TT();	
    void integrate_cl_EE();
    void integrate_cl_TE();
    void integrate_cl_all();
    void find_cl_max(int &lmax, double &clmax);

    double get_bessel(double x, int l);
    double get_high_source_T(double x, int k);
    double get_high_source_E(double x, int k);
    double get_thetal_T(double x, int l);
    double get_thetal_E(double x, int l);
    double get_cl_TT(double l);
    double get_cl_EE(double l);
    double get_cl_TE(double l);
    double get_primordial_power_spectrum(double k);
    double get_matter_power_spectrum(double k);

};

// Theta_l_T ODE Structure
struct derivs_thetal_T {
  Powerspectrum *P;
  int lnumber;

  derivs_thetal_T(Powerspectrum *PP, int lin) : P(PP) {
    lnumber = lin;
  }
  ~derivs_thetal_T(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double source_T, arg;

    source_T = P->get_high_source_T(x, P->k_current);
    arg		 = P->k_high[P->k_current]*(P->BG->get_eta(0.0) - P->BG->get_eta(x));

    for (int l = 0; l<lnumber; l++) {
      dydx[l] = source_T * P->get_bessel(arg, l);
    }
  }
};

// Theta_l_E ODE Structure
struct derivs_thetal_E {
  Powerspectrum *P;
  VecDoub lfact;

  derivs_thetal_E(Powerspectrum *PP) : P(PP) {
    lfact = VecDoub(l_num);
    for(int l = 0;l<l_num;l++)
      lfact(l) = sqrt( (ls[l]+2.0)*(ls[l]+1.0)*ls[l]*(ls[l]-1.0) );
  }
  ~derivs_thetal_E(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double source_E, arg;

    source_E = P->get_high_source_E(x, P->k_current);
    arg		 = P->k_high[P->k_current]*(P->BG->get_eta(0.0) - P->BG->get_eta(x));

    for (int l = 0; l<l_num; l++) {
      dydx[l] = source_E * P->get_bessel(arg, l) * lfact(l);
    }
  }
};

// Combined Theta_l_T and Theta_l_E ODE Structure
struct derivs_thetal_TE {
  Powerspectrum *P;

  derivs_thetal_TE(Powerspectrum *PP) : P(PP) {}
  ~derivs_thetal_TE(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double bessel, source_E, source_T, arg, fact;
    int know = P->k_current;

    source_E = P->get_high_source_E(x, know);
    source_T = P->get_high_source_T(x, know);
    arg		 = P->k_high[know]*(P->BG->get_eta(0.0) - P->BG->get_eta(x));

    for (int l = 0; l<l_num; l++) {
      bessel			= P->get_bessel(arg, l);
      dydx[l]			= source_T * bessel;
      dydx[l_num+l]	= source_E * bessel * sqrt( (ls[l]+2.0)*(ls[l]+1.0)*ls[l]*(ls[l]-1.0) );
    }
  }
};

// C_l_TT ODE Structure
struct derivs_cl_TT {
  Powerspectrum *P;
  double units;

  derivs_cl_TT(Powerspectrum *PP) : P(PP) {
    units  = 4.0 * pi * pow(1.0e6*P->cosm_simu->T_0,2);
  }
  ~derivs_cl_TT(){ P = NULL; } 

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double tl, power;

    power = P->get_primordial_power_spectrum(x)/x * units;

    for (int l = 0; l<l_num; l++) {
      tl =  P->get_thetal_T(x,l);
      dydx[l] =  tl*tl*power;
    }
  }
};

// C_l_EE ODE Structure
struct derivs_cl_EE {
  Powerspectrum *P;
  double units;

  derivs_cl_EE(Powerspectrum *PP) : P(PP) {
    units  = 4.0 * pi * pow(1.0e6*P->cosm_simu->T_0,2);
  }
  ~derivs_cl_EE(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double el, power;

    // Last factor is to get the result in (muK)^2
    power = P->get_primordial_power_spectrum(x)/x * units;

    for (int l = 0; l<l_num; l++) {
      el =  P->get_thetal_E(x,l);
      dydx[l] =  el*el*power;
    }
  }
};

// C_l_TE ODE Structure
struct derivs_cl_TE {
  Powerspectrum *P;
  double units;

  derivs_cl_TE(Powerspectrum *PP) : P(PP) { 
    units  = 4.0 * pi * pow(1.0e6*P->cosm_simu->T_0,2);
  }
  ~derivs_cl_TE(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double tl, el;

    // Last factor is to get the result in (muK)^2
    double power = P->get_primordial_power_spectrum(x)/x * units;

    for (int l = 0; l<l_num; l++) {
      tl =  P->get_thetal_E(x,l);
      el =  P->get_thetal_T(x,l);
      dydx[l] =  el*tl*power;
    }
  }
};

// Combined C_l TT, EE and TE ODE Structure
struct derivs_cl_all {
  Powerspectrum *P;
  double units, factor;

  derivs_cl_all(Powerspectrum *PP) : P(PP) { 
    units  = 4.0 * pi * pow(1.0e6*P->cosm_simu->T_0,2);
  }
  ~derivs_cl_all(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double tl, el, power;

    // Last factor is to get the result in (muK)^2
    power = P->get_primordial_power_spectrum(x)/x * units;

    for (int l = 0; l<l_num; l++) {
      tl =  P->get_thetal_T(x,l);
      el =  P->get_thetal_E(x,l);
      dydx[l]			= tl*tl*power;
      dydx[l+l_num]	= tl*el*power;
      dydx[l+2*l_num]	= el*el*power;
    }
  }
};

// sigma8 ODE Structure
struct derivs_sigma8 {
  Powerspectrum *P;

  derivs_sigma8(Powerspectrum *PP) : P(PP) {}
  ~derivs_sigma8(){ P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    double window, r8, power, kR8;

    power = P->get_matter_power_spectrum(x) / x * pow(x / 2998.0,3)/ (2.0*pi*pi);
    kR8 = 8.0/2998.0 * x;

    if (kR8<0.01){
      window = 1.0 - 0.1*kR8*kR8;
    }else{
      window = 3.0/(kR8*kR8*kR8) * (sin(kR8) - kR8*cos(kR8));
    }

    dydx[0] = power * window * window;
  }
};

