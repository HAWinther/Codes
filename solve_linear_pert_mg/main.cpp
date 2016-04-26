#include <iostream>
#include <fstream>
#include "PofkModifiedGravity.h"

//==================================================
//
// Code to compute the ratio P_MG(k,a)/P_LCDM(k,a)
// in linear perturbation theory
//
// The modified gravity model is specified by
// defining Geff(a,k) and H(a) below
//
// Code units:
//  k is in units of h/Mpc
//  x = log(a)
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//
//==================================================

//===========================================
// G_effective/G for Hu-Sawicky f(R) gravity
// Below fac = m(a)^2 / [a^2k^2]
// For f(R) we have m(a)^2 = 1/[3f_RR(a)]
//===========================================

double Geff(double a, double kH0){
  double OlOm =  params.Omegal/params.Omegam, fac;
  fac  = pow( a / (kH0 * 2998.0), 2) / ( (params.n_fofr + 1.0) * params.fofr0 ) * params.Omegam;
  fac *= (1.0 + 4.0*OlOm)*pow( (1.0/(a*a*a) + 4.0*OlOm) / (1.0 + 4.0*OlOm), params.n_fofr + 2.0);
  return 1.0 + (1.0/3.0)/(1.0 + fac);
}

//===========================================
// Hubble function H/H0(x) in LCDM and MG
//===========================================

double hubble_lcdm(double x){
  return sqrt(params.Omegam * exp(-3.0*x) + params.Omegal);
}
double hubble_mg(double x){
  return sqrt(params.Omegam * exp(-3.0*x) + params.Omegal);
}

//===========================================
// LCDM equation for delta(a)
// x = log(a)
// The equation and variables are:
// 0: ddelta[x]/dx = q/(Ha^2)
// 1: dq/dx = 1.5 Om/(aH) delta
//===========================================

void dydx_lcdm_eq(const double& x, const std::vector<double>& y, std::vector<double>& dydx){
  double H = hubble_lcdm(x);
  dydx[0] = y[1]/(H * exp(2.0*x));
  dydx[1] = 1.5*params.Omegam/(H * exp(x))*y[0];
}

//===========================================
// MG equation for delta(k,a)
// x = log(a)
// k in units of h/Mpc is set globally in know
// The equation and variables are:
// 0: ddelta[x]/dx = q/(Ha^2)
// 1: dq/dx = 1.5 Om/(aH) delta Geff(a,k)
//===========================================

void dydx_mg_eq(const double& x, const std::vector<double>& y, std::vector<double>& dydx){
  double H = hubble_mg(x);
  dydx[0] = y[1]/(H * exp(2.0*x));
  dydx[1] = 1.5*params.Omegam/(H * exp(x))*y[0]*Geff(exp(x),params.know);
}

//===========================================
// Solve for the ratio of the linear 
// power-spectra in MG and LCDM:
// P(k)_MG / P(k)_LCDM
//===========================================

void solve_linear_pert(){
  double aini, aend, kmin, kmax, k;
  int n, nk;
  bool verbose = true;
  PofkModifiedGravity solver;

  // Number of time and k-points
  n = 1000;
  nk = 100;

  // Set cosmological and modified gravity parameters parameters
  params.Omegam = 0.3;
  params.Omegal = 0.7;
  params.fofr0  = 1e-6;
  params.n_fofr = 1.0;

  // Initial and final scale factor (the one we output P/P_LCDM at)
  aini = 0.002;
  aend = 1.0;

  // Min max k-values in h/Mpc
  kmin = 0.001;
  kmax = 100.0;

  // Solve for the ratio P_mg(k,aend) / P_LCDM(k,aend)
  solver = PofkModifiedGravity(kmin, kmax, nk, aini, aend, n, dydx_lcdm_eq,  dydx_mg_eq);
  solver.solve(verbose);

  // Output to file
  std::ofstream fp("ratio_pofk.txt");
  fp << "#k (h/Mpc)    vs  P(k)_linear_f(R) / P(k)_linear_LCDM   at   a = " << aend << std::endl;
  for(int i = 0; i < nk; i++){
    k = exp(log(kmin) + log(kmax/kmin)*i/double(nk-1));
    fp << k << " " << solver.get_pofk_ratio(k) << std::endl;
  }
  fp.close();

}

int main(int argv, char **argc){
  solve_linear_pert();
}
