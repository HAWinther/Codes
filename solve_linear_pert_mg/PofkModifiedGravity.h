#ifndef _LINEARPERTSOLVER_
#define _LINEARPERTSOLVER_
#include "OdeSolver.h"
#include "Spline.h"

//========================================
// Global book-keeping variable
//========================================

struct CosmologyParameters {
  // Cosmology
  double Omegal, Omegam;

  // Modified gravity f(R) params
  double fofr0, n_fofr;

  // Current value of k/H0
  double know;
} params;

class PofkModifiedGravity{

  private:

    //==================================================
    //
    // Compute P_MG / P_LCDM(k,a) where P_MG (P_LCDM) is 
    // the linear power-spectrum in modified gravity (LCDM)
    //
    // k is in units of 1/H0
    // Time-variable is x = log(a)
    //
    // Hans A. Winther (2015) (hans.a.winther@gmail.com)
    //
    //==================================================

    double kmin;                            // kmin we integrate pert.eq. for
    double kmax;                            // kmax we integrate pert.eq. for
    int nk;                                 // Number of k-points we solve the pert.eq. for
    double aini;                            // Integrate from a = aini
    double aend;                            // Integrate till a = aend
    int nint;                               // Number of a-points we store pofk ratio in
    std::vector<double> kk;                 // The k-array
    bool solved;                            // Have we solved already?
    DSpline pofkratio_spline;               // The spline containing the pofk ratio
    DSpline delta_lcdm_x_spline;            // Spline of delta_m(a) for LCDM (this is k-independent)
    std::vector<DSpline> delta_mg_x_spline; // Spline of delta_m(a,k) for modified gravity

    // The right hand side of the perturbation equation for LCDM and MG
    void (*derivs_lcdm)(const double &, const std::vector<double>&, std::vector<double>&);
    void (*derivs_mg)(const double &, const std::vector<double>&, std::vector<double>&);

  public:

    PofkModifiedGravity() :
                        kmin(0.0),
                        kmax(0.0),
                        nk(0),
                        aini(0.0),
                        aend(0.0),
                        nint(0),
                        kk(0),
                        solved(false),
                        derivs_mg(NULL),
                        derivs_lcdm(NULL) {}

    PofkModifiedGravity(double _kmin, double _kmax, int _nk, 
                        double _aini, double _aend, int _nint, 
                        void (*_derivs_lcdm)(const double &, const std::vector<double>&, std::vector<double>&),  
                        void (*_derivs_mg)(const double &, const std::vector<double>&, std::vector<double>&)) : 
                        kmin(_kmin),
                        kmax(_kmax),
                        nk(_nk),
                        aini(_aini),
                        aend(_aend),
                        nint(_nint),
                        kk(std::vector<double>(nk, 0.0)),
                        solved(false),
                        delta_mg_x_spline(std::vector<DSpline>(nk)),
                        derivs_mg(_derivs_mg),
                        derivs_lcdm(_derivs_lcdm) {}

    //========================================
    // Solve the LCDM and MG ODE's
    //========================================
    void solve(bool verbose);
    
    //========================================
    // P_MG/P_LCDM(a,k), k in units of h/Mpc
    // If a is not given then a=aend is used
    //========================================
    double get_pofk_ratio(double k, double a);
    double get_pofk_ratio(double k);
};

void PofkModifiedGravity::solve(bool verbose){
  double delta_lcdm_0;
  double xmin, xmax;
  std::vector<double> x, y, ic;
  std::vector<double> pofk(nk, 0.0);

  // Set the initial conditions
  // At early times delta = a = exp(x) ic[0] = ic[1] = 1.0
  xmin = log(aini);
  xmax = log(aend);
  ic   = std::vector<double>(2, 1.0);

  // Set up solver and solve the ODEs
  OdeSolver lcdm(nint, 2, derivs_lcdm);
  lcdm.set_ic(xmin, xmax, ic);
  lcdm.solve(false);

  // The values of delta today
  y = lcdm.get_y(0);
  x = lcdm.get_x();
  delta_lcdm_0 = y[nint-1];

  // DSpline delta_lcdm(x)
  delta_lcdm_x_spline = DSpline(x, y, 1e30, 1e30, "delta_LCDM");

  // Loop over k'values and caluclate delta_m(a,k)
  for(int i = 0; i < nk; i++){
    // Set params k-values
    params.know = exp(log(kmin) + log(kmax/kmin)*i/double(nk-1));

    // Set up solver and solve ODE
    OdeSolver mg(nint, 2, derivs_mg);
    mg.set_ic(xmin, xmax, ic);
    mg.solve(false);

    // The values of delta(x)
    x = lcdm.get_x();
    y = mg.get_y(0);
    
    // Spline it up
    delta_mg_x_spline[i] = DSpline(x, y, 1e30, 1e30, "delta_mg");

    if(verbose && i % 20 == 0) {
      // The value delta_m(k,a=1) / delta_m_LCDM(k)
      std::cout << "k: " << params.know << " Delta/Delta_LCDM: " << y[nint-1]/delta_lcdm_0 << std::endl;
    }

    // Store values
    kk[i] = params.know;
    pofk[i] = pow(y[nint-1]/delta_lcdm_0, 2);
  }

  // Spline up the ratio
  pofkratio_spline = DSpline(kk, pofk, 1e30, 1e30, "P(k) / P(k)_LCDM ratio");
  solved = true;
}
    
double PofkModifiedGravity::get_pofk_ratio(double k, double a){
  double x = log(a), delta_lcdm, ratio, delta_mg_low, delta_mg_high;
  int indlow, indhigh;

  // Get LCDM value (note: same value for all k)
  delta_lcdm = delta_lcdm_x_spline.f(x);

  // Get index in k-array
  indlow = int(log(k/kmin)/log(kmax/kmin)*(nk-1) + 0.5);

  // Bounds check
  if(indlow >= nk || indlow < 0){
    std::cout << "Error: out of bounds in get_pofk_ratio_of_a_k" << std::endl;
    std::cout << "kmin < k < kmax : " << kmin << " " << k << " " << kmax << std::endl;
    exit(1);
  }
  indhigh = indlow + 1;

  // Spline lookup + linear interpolation
  if(indhigh == nk){
    delta_mg_low  = delta_mg_x_spline[indlow].f(x);
    ratio = pow(delta_mg_low/delta_lcdm, 2);	
  } else {
    delta_mg_low  = delta_mg_x_spline[indlow].f(x);
    delta_mg_high = delta_mg_x_spline[indhigh].f(x);
    ratio  = pow(delta_mg_low/delta_lcdm, 2);
    ratio += (k - kk[indlow])/(kk[indhigh] - kk[indlow]) * (pow(delta_mg_high/delta_lcdm,2) - pow(delta_mg_low/delta_lcdm,2)); 
  }

  return ratio;
}

double PofkModifiedGravity::get_pofk_ratio(double k){
  if(!solved){
    std::cout << "Error in PofkModifiedGravity. Cannot get pofk ratio before solving. Returning 1.0" << std::endl;
    return 1.0;
  }
  return pofkratio_spline.f(k);
}

#endif
