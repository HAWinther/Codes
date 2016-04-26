#include <iomanip>
#include <vector>
#include "OdeSolver.h"

//=======================================================
// 
// Demonstration of simple class to solve first order 
// coupled ODEs.
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//
//=======================================================

//=======================================================
// The single ODE {dy/dx = - 2xy} which has the 
// analytical solution:
// 
//     y(x) = y(0) * Exp(-x^2)
//=======================================================
void ode_single_rhs(const double &x, const std::vector<double> &y, std::vector<double> &dydx){
  dydx[0] = -2.0 * x * y[0];
}
double ode_single_analytical(double x, std::vector<double> &ic){
  return ic[0] * exp(-x * x);
}

//=======================================================
// The coupled ODE {dy_1/dx = y_2, dy_2/dx = 1} which 
// has the analytical solution:
//
//     y_1(x) = y_1(0) + y_2(0)x + x^2/2
//     y_2(x) = y_2(0) + x
//=======================================================
void ode_coupled_rhs(const double &x, const std::vector<double> &y, std::vector<double> &dydx){
  dydx[0] = y[1];
  dydx[1] = 1.0;
}
double ode_coupled_analytical(double x, std::vector<double> &ic, int i){
  if(i==0) return ic[0] + ic[1] * x + x * x / 2.0;
  return ic[1] + x;
}

//=======================================================
// Solve a single ODE
//=======================================================

void solve_single_ode(){
  std::vector<double> x, y, ic;
  double xmin = 0.0, xmax = 1.0;
  int n, neq;
  bool verbose = false;

  std::cout << "=========================================" << std::endl;
  std::cout << "Solve single ODE:     dy/dx = -2xy       " << std::endl;
  std::cout << "=========================================" << std::endl;
  
  // Number of points between xmin and xmax to store the solution in
  n = 20;

  // Number of equations
  neq = 1;

  // Initial conditions for ODE: y(0) = 1
  ic = std::vector<double>(neq, 1.0);

  OdeSolver myode(n, neq, ode_single_rhs);
  myode.set_ic(xmin, xmax, ic);
  myode.solve(verbose);

  // Extract solution
  x = myode.get_x();
  y = myode.get_y();

  // Print data
  std::cout << std::endl << "Solution:" << std::endl;
  for(int i = 0; i < n; i++){
    std::cout << std::setw(2) << i << " / " << n << "  x: " << std::setw(12) << x[i];
    std::cout << "  y: " << std::setw(12) << y[i] << " Error: " <<  std::setw(12) << y[i] - ode_single_analytical(x[i], ic) << std::endl; 
  }
  std::cout << std::endl;
}

//=======================================================
// Solve a coupled ODE system
//=======================================================

void solve_coupled_ode(){
  std::vector<double> x, y1, y2, ic;
  double xmin = 0.0, xmax = 1.0;
  int n, neq;
  bool verbose = true;

  std::cout << "=========================================" << std::endl;
  std::cout << "Solve coupled ODE:   dy1/dx=y2, dy2/dx=1 " << std::endl;
  std::cout << "=========================================" << std::endl;
  
  // Number of points between xmin and xmax to store the solution in
  n = 10;

  // Number of equations
  neq = 2;

  // Initial conditions for ODEs: y1(0) = y2(0) = 1
  ic = std::vector<double>(neq, 1.0);

  OdeSolver myode(n, neq, ode_coupled_rhs);
  myode.set_ic(xmin, xmax, ic);
  
  // Set precision goal and performance params: [epsilon, h_start, hmin]
  myode.set_precision(1e-20, 1e-12, 0.0);   
  
  myode.solve(verbose);                     

  // Extract solution
  x  = myode.get_x();
  y1 = myode.get_y(0);
  y2 = myode.get_y(1);

  // Print data
  std::cout << std::endl << "Solution:" << std::endl;
  for(int i = 0; i < n; i++){
    std::cout <<  std::setw(2) << i+1 << " / " << n << "  x: " <<  std::setw(12) << x[i];
    std::cout << "  y[1]: " <<  std::setw(12) << y1[i] << " Error[1]: " <<  std::setw(12) << y1[i] - ode_coupled_analytical(x[i], ic, 0);
    std::cout << "  y[2]: " <<  std::setw(12) << y2[i] << " Error[2]: " <<  std::setw(12) << y2[i] - ode_coupled_analytical(x[i], ic, 1) << std::endl; 
  }
  std::cout << std::endl;
}

int main(int argv, char **argc){
  solve_single_ode();
  solve_coupled_ode();
}

