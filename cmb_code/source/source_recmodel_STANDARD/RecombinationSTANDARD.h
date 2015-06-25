#pragma once

#include "../Background.h"
#include "../Recombination.h"

class RecombinationSTANDARD : public Recombination {
  public:
    Background *BG;

    RecombinationSTANDARD(){}
    RecombinationSTANDARD(Background *BB, ParamSimu *param) { BG = BB; BGSUB = BG; cosm_simu = param; }
    ~RecombinationSTANDARD() { BG = NULL; cosm_simu = NULL; }

    void calc_n_e();
    double get_derivs_peebles(double x, VecDoub_I &y);

};
