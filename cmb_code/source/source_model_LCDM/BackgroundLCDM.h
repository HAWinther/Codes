#pragma once

#include "../Parameters.h"
#include "../Background.h"

class BackgroundLCDM : public Background {
  public:

    BackgroundLCDM(ParamSimu *param);
    ~BackgroundLCDM(){ cosm_simu = NULL; };

    // Implementation of virtual Background methods
    double H(const double x);
    double dH(const double x);
    double ddH(const double x);
    double H_c(const double x);
    double dH_c(const double x);
    double ddH_c(const double x);
    void calc_expansion_history();
};

