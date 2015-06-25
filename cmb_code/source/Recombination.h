
//================================================
// Abstract class containing general methods to 
// compute the recombination history in any
// given background. There is also the possibilty
// to change the physics during recombination.
//================================================

#pragma once

#include "Parameters.h"
#include "Background.h"
#include "Matrix.h"    
#include "Spline.h"
#include "Odeint.h"

class Recombination {
  public:

    Background *BGSUB;

    ParamSimu *cosm_simu;

    VecDoub X_e, n_e;
    VecDoub tau, tau2;
    VecDoub g, g2;

    Spline n_e_spline, tau_spline, tau2_spline, g_spline, g2_spline;

    Recombination() {}
    virtual ~Recombination() { BGSUB = NULL; cosm_simu = NULL; }

    void calc_rec_history();
    void calc_tau();
    void calc_g();
    double get_tau(double x);
    double get_dtau(double x);
    double get_ddtau(double x);
    double get_g(double x);
    double get_dg(double x);
    double get_ddg(double x);
    double get_n_e(double x);
    double get_x_e(double x);
    void get_all_tau(double x, double &Tau, double &dTau, double &ddTau);
    void get_all_g(double x, double &G, double &dG, double &ddG);
    void get_all_tau_g(double x, double &Tau, double &dTau, double &ddTau, double &G, double &dG, double &ddG);

    virtual void calc_n_e() = 0;
    virtual double get_derivs_peebles(double x, VecDoub_I &y) = 0;
};

// Peebles ODE Structure
struct derivs_peebles {
  Recombination *R;

  derivs_peebles(Recombination *RR) : R(RR) {}
  ~derivs_peebles(){ R = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub &dydx){
    dydx[0] = R->get_derivs_peebles(x,y);
  }
};

// Tau ODE Structure
struct derivs_tau {
  Recombination *R;

  derivs_tau(Recombination *RR) : R(RR) {}
  ~derivs_tau(){ R = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx){
    dydx[0] = R->get_dtau(x);
  }
};
