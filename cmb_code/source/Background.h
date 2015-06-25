
//========================================================
// Base class which defines functions for computing 
// the background evolution for a general background
//========================================================

#pragma once

#include "Parameters.h"
#include "Matrix.h"    
#include "Spline.h"
#include "Odeint.h"

class Background {
  public:
    ParamSimu *cosm_simu;

    // Global variables common to all backgrounds
    VecDoub eta, eta2, x_main;
    Spline eta_spline;

    Background(){};
    virtual ~Background(){};

    // Global functions common to all backgrounds
    void calc_eta();
    void get_all_H_c(const double x, double &Hc, double &dHHc, double &ddHc);
    double get_eta(double x);

    // Virtual functions which needs a background dependent implementation
    virtual double H(const double x) = 0;
    virtual double dH(const double x) = 0;
    virtual double ddH(const double x) = 0;
    virtual double H_c(const double x) = 0;
    virtual double dH_c(const double x) = 0;
    virtual double ddH_c(const double x) = 0;
    virtual void calc_expansion_history() = 0;

    void output_omega();

};

// General ODE structure for the conformal time eta(x)
struct derivs_eta {
  Background *B;
  derivs_eta(Background *A){
    B = A;
  }
  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    dydx[0] = 1.0/(B->H_c(x));
  }
};

