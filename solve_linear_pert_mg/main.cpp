#include <iostream>
#include <fstream>
#include "linear_pert_solver.h"

/////////////////////////////////////////////////////
// Specify Geff(a,k) below
//
// Code units:
//  k is in units of h/Mpc
//  x = log(a)
//
// Outputs the ratio of the linear P(k)
// in modified gravity to that of LCDM
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
/////////////////////////////////////////////////////

//////////////////////////////////////////
// G_effective/G for modified gravity
// k is k_physical/H0
// Here for Hu-Sawicky f(R)
//////////////////////////////////////////

realT Geff(realT a, realT k){
  realT OlOm =  global.Omegal/global.Omegam, fac;

  // This is [m^2/H0^2] * a^2 * [H0/k]^2 and for f(R) we have m^2 = 1/[3f_RR]
  fac = pow(a/k/2998.0,2)/((global.n_fofr+1.0)*global.fofr0) * global.Omegam * (1. + 4.*OlOm)*pow( (1.0/(a*a*a) + 4.*OlOm) / (1.0 + 4.*OlOm), global.n_fofr + 2.0);

  return 1.0 + 0.333333/(1.0+fac);
}

///////////////////////////////////////////////////////
// LCDM equation for delta(a)
// x = log(a)
// The equation and variables are:
// 0: ddelta[x]/dx = q/(Ha^2)
// 1: dq/dx = 1.5 Om/(aH) delta
///////////////////////////////////////////////////////

void dydx_lcdm_eq(realT x, realT *y, realT *dydx){
  realT H = sqrt(global.Omegam*exp(-3.0*x) + global.Omegal);
  dydx[0] = y[1]/(H * exp(2.0*x));
  dydx[1] = 1.5*global.Omegam/(H * exp(x))*y[0];
}

///////////////////////////////////////////////////////
// MG equation for delta(k,a)
// x = log(a)
// k in units of h/Mpc is set globally in know
// The equation and variables are:
// 0: ddelta[x]/dx = q/(Ha^2)
// 1: dq/dx = 1.5 Om/(aH) delta Geff(a,k)
///////////////////////////////////////////////////////

void dydx_mg_eq(realT x, realT *y, realT *dydx){
  realT H = sqrt(global.Omegam*exp(-3.0*x) + global.Omegal);
  dydx[0] = y[1]/(H * exp(2.0*x));
  dydx[1] = 1.5*global.Omegam/(H * exp(x))*y[0]*Geff(exp(x),global.know);
}

//////////////////////////////////////////
// Solve for the ratio of the linear 
// power-spectra in MG and LCDM:
// P(k)_MG / P(k)_LCDM
//////////////////////////////////////////

void solve_linear_pert(){
  realT aini, aend, kmin, kmax, k;
  int n, nk;
  bool verbose = true;

  // Make solver
  PofkModifiedGravity *psolver;

  // Number of time and k-points
  n = 1000;
  nk = 100;

  // Set cosmological and modified gravity parameters parameters
  global.Omegam = 0.3;
  global.Omegal = 0.7;
  global.fofr0  = 1e-6;
  global.n_fofr = 1.0;

  // Initial and final scale factor (the one we output P/P_LCDM at)
  aini = 0.002;
  aend = 1.0;

  // Min max k-values in h/Mpc
  kmin = 0.001;
  kmax = 100.0;

  // Solve for the ratio P_mg(k,aend) / P_LCDM(k,aend)
  psolver = new PofkModifiedGravity(kmin, kmax, nk, aini, aend, n, dydx_lcdm_eq,  dydx_mg_eq);
  psolver->solve(verbose);

  // Output to file
  std::ofstream fp("ratio_pofk.txt");
  fp << "#k (h/Mpc)    vs  P(k)_linear_f(R) / P(k)_linear_LCDM   at   a = " << aend << std::endl;
  for(int i=0;i<nk;i++){
    k = exp(log(kmin) + log(kmax/kmin)*i/realT(nk-1));
    fp << k << " " << psolver->pofk_ratio(k) << std::endl;
  }
  fp.close();

  delete psolver;
}

int main(int argv, char **argc){
  solve_linear_pert();
}
