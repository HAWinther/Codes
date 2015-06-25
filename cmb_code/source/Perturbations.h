
//===========================================================
// Abstract class which contains all the methods to
// integrate perturbations in a general model
//===========================================================

#pragma once

#include "Parameters.h"
#include "Background.h"
#include "Recombination.h"
#include "Matrix.h"    
#include "Spline.h"
#include "Odeint.h"

class Perturbations {
  public:
    Background *BGSUB;
    Recombination *RECSUB;

    ParamSimu *cosm_simu;

    MatDoub Phi, deltaM;

    VecDoub x_rec;
    MatDoub S_low, SE_low;
    VecDoub ks;

    // Global book-keeping variables
    double k_current, x_current;

    Perturbations() {}
    virtual ~Perturbations() { BGSUB = NULL; RECSUB = NULL; cosm_simu = NULL; }

    void integrate_perturbations();
    double get_tight_coupling_time(double k);

    virtual void set_initial_conditions_tight_coupling(double x, VecDoub_IO &y) = 0;
    virtual void set_initial_conditions_full(double x, VecDoub_I &ytc, VecDoub_IO &y) = 0;

    virtual void get_derivs_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dydx) = 0;
    virtual void get_derivs_full_equation(const double x, VecDoub_I &y, VecDoub_O &dydx) = 0;

    virtual void get_jacobian_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) = 0;
    virtual void get_jacobian_full_equation(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) = 0;

    void init_store_values();
    void store_values(double x, double k, VecDoub *y, bool tight_coupling); 
    void analyze_store_values();
};

// ODE structure for tight coupling equations
struct derivs_tight_coupling {
  Perturbations *P;

  derivs_tight_coupling(Perturbations *PP) : P(PP) {}
  ~derivs_tight_coupling() { P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx){
    P->get_derivs_tight_coupling(x, y, dydx);		
  }

  void jacobian(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy){
    P->get_jacobian_tight_coupling(x, y, dfdx, dfdy);
  }
};

// ODE structure for full perturbations equations
struct derivs_full {
  Perturbations *P;

  derivs_full(Perturbations *PP) : P(PP) {}
  ~derivs_full() { P = NULL; }

  void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx){
    P->get_derivs_full_equation(x, y, dydx);
  }

  void jacobian(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy){
    P->get_jacobian_full_equation(x, y, dfdx, dfdy);
  }
};

