#include<iomanip>
#include "Spline.h"

void output_result(double xnow, double ynow, double dynow);
void test_spline();

//=========================================================================
// Code to make simple qubic splines
// Use demonstrated for the function y = x^2
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//=========================================================================

void test_spline(){
  double xmin, xmax, xnow, ynow, dynow, dydx1, dydxn;
  Spline<double> func, *gfunc, hfunc;
  std::vector<double> x, y;
  int n, splinetype;

  // Calculate y = x2 for x in [-10,10]
  n    = 15;
  xmin = -10.0; xmax =  10.0;
  for(int i=0;i<n;i++){
    x.push_back(xmin + (xmax-xmin)*i/double(n-1));
    y.push_back(x[i]*x[i]);
  }

  //================================================================
  // Example 1: 
  //================================================================

  // If the x-array spacing is not regular (or known)
  splinetype = 0;

  // Natural spline boundary conditions
  dydx1 = dydxn = 1e30;

  // Contruct spline
  func = Spline<double>(x,y,n,dydx1,dydxn,splinetype,"Test spline 1");

  std::cout << "Test 1 of spline for y = x2" << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/double(n-1);
    ynow  = func(xnow); 
    dynow = func.dfdx(xnow); 
    output_result(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  //================================================================
  // Example 2: Using new + direct lookup
  //================================================================

  // The x-array is lineary spaced so use direct lookup to speed it up
  splinetype = 1;

  // Natural spline boundary conditions
  dydx1 = dydxn = 1e30;

  // Contruct spline
  gfunc = new Spline<double>(x,y,n,dydx1,dydxn,splinetype,"Test spline 2");

  std::cout << "Test 2 of spline for y = x2 | delta_y = y-x2" << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/double(n-1);
    ynow  = gfunc[0](xnow);
    dynow = gfunc->dfdx(xnow);
    output_result(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  //================================================================
  // Example 3: Use boundary conditions for dy/dx
  //================================================================

  // The x-array is lineary spaced so use direct lookup to speed it up
  splinetype = 1;

  // Analytical boundary conditions
  dydx1 = 2*xmin;
  dydxn = 2*xmax;

  // Contruct spline
  hfunc = Spline<double>(x,y,n,dydx1,dydxn,splinetype,"Test spline 3");

  std::cout << "Test 3 of spline for y = x2 with analytical BC (y'=2x) | delta_y = y-x2 " << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/double(n-1)+1.0;
    ynow  = hfunc.f(xnow);
    dynow = hfunc.dfdx(xnow);
    output_result(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  // Clean up memory
  delete gfunc;
}

void output_result(double xnow, double ynow, double dynow){
  std::cout << "  x = " << std::setw(12) << xnow;
  std::cout << "  y = " << std::setw(12) << ynow;
  std::cout << "  |y - x2|    = " << std::setw(12) << fabs(ynow - xnow*xnow);
  std::cout << "  |dydx - 2x| = " << std::setw(12) << fabs(dynow-2.0*xnow);
  std::cout << std::endl; 
}

int main(int argc, char **argv){
  test_spline();
}
