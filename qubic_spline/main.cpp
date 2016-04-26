#include <iomanip>
#include <vector>
#include "Spline.h"

//==================================================================
// Example use of the QubicSpline class
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//==================================================================

void output(double xx, double yy, double dy);
double f(double x){ return x * x; }
double dfdx(double x){ return 2.0 * x; }

int main(int argc, char** argv){
  double xmin, xmax, xx, yy, dy;
  float *xp, *yp;
  int n, m;
  std::vector<double> x, y;
  DSpline f_spline, g_spline;

  //================================================================
  // Set (x,y) to satisfy y = x^2 for x in [1,10]
  //================================================================
  n     = 10,  m = 2*n;
  xmin  = 1.0, xmax  = 10.0;
  x = y = std::vector<double>(n, 0.0);
  for(int i = 0; i < n; i++){
    x[i] = xmin + (xmax - xmin) * i/double(n-1);
    y[i] = f(x[i]);
  }
  
  //================================================================
  // Example 1: Natural spline
  //================================================================

  f_spline = DSpline(x, y, "Test spline 1");

  std::cout << "[Test 1 of spline for f(x) = x^2]" << std::endl;
  for(int i = 0; i < m; i++){
    xx = xmin + (xmax - xmin) * i/double(m-1);
    yy = f_spline(xx); 
    dy = f_spline.dfdx(xx); 
    output(xx, yy, dy);
  }
  std::cout << std::endl;

  //================================================================
  // Example 2: We supply boundary conditions
  //================================================================

  g_spline = DSpline(x, y, dfdx(xmin), dfdx(xmax), "Test spline 2");

  std::cout << "[Test 2 of spline for f(x) = x^2 with analytical BC (dfdx = 2x)]" << std::endl;
  for(int i = 0; i < m; i++){
    xx = xmin + (xmax - xmin) * i/double(m-1);
    yy = g_spline.f(xx);
    dy = g_spline.dfdx(xx);
    output(xx, yy, dy);
  }
  std::cout << std::endl;

  //================================================================
  // Example 3: Init from *, using float and show some features
  //================================================================

  xp = new float[n];
  yp = new float[n];
  for(int i = 0; i < n; i++){
    xp[i] = x[i];
    yp[i] = y[i];
  }
  FSpline h_spline = FSpline(xp, yp, n, "Test spline 3");

  // Set to show error message if x is out of bounds
  h_spline.show_warning(true);

  std::cout << "[Test 3 with x out of bounds and info about spline]" << std::endl;
  h_spline(xmin - 0.5);
  h_spline(xmax + 0.5);

  // Output info about the spline
  h_spline.info();

  delete[] xp;
  delete[] yp;
}

//================================================================
// Dump to screen
//================================================================
void output(double xx, double yy, double dy){
  std::cout.precision(3);
  std::cout << " x = " << std::setw(8) << xx;
  std::cout << " f(x) = " << std::setw(8) << yy;
  std::cout << " Error f = " << std::setw(8) << fabs( yy / f(xx) - 1.0 );
  std::cout << " Error df/dx = " << std::setw(8) << fabs( dy / dfdx(xx) - 1.0 );
  std::cout << std::endl; 
  std::cout.precision(0);
}

