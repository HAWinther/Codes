#pragma once

#include "../Parameters.h"
#include "BackgroundLCDM.h"
#include "../source_recmodel_STANDARD/RecombinationSTANDARD.h"
#include "../Perturbations.h"

class PerturbationsLCDM : public Perturbations {
  public:
    Background *BG;
    Recombination *REC;

    PerturbationsLCDM() {}
    PerturbationsLCDM(Background *BB, Recombination *RR, ParamSimu *param){
      BG = BB; REC = RR; BGSUB = BG; RECSUB = REC; cosm_simu = param;
    }
    ~PerturbationsLCDM() { BG = NULL; REC = NULL; cosm_simu = NULL; }

    void set_initial_conditions_tight_coupling(double x, VecDoub_IO &y);
    void set_initial_conditions_full(double x, VecDoub_I &ytc, VecDoub_IO &y);

    void get_derivs_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dydx);
    void get_derivs_full_equation(const double x, VecDoub_I &y, VecDoub_O &dydx);

    void get_jacobian_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy);
    void get_jacobian_full_equation(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy);

};
