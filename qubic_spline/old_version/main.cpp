#include<iomanip>
#include "spline.h"

//=========================================================================
// Code to make simple qubic splines
//
// As an example we show how to make and use a spline for the function y = x^2
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//=========================================================================

void print(realT xnow, realT ynow, realT dynow){
  std::cout << "   x = " << std::setw(12) << xnow << "   y = " << std::setw(12) << ynow;
  std::cout << "   |y - x^2| = " << std::setw(12) << fabs(ynow - xnow*xnow);
  std::cout << "   |dy/dx - 2x| = " << std::setw(12) << fabs(dynow-2.0*xnow) << std::endl; 
}

void make_spline(){
  Spline func, *gfunc, hfunc;
  realT *x, *y;
  realT xmin, xmax, xnow, ynow, dynow;
  int n;

  // Number of points in array
  n = 25;

  // Initialize arrays 
  x = new realT[n];
  y = new realT[n];

  // Make y[i] = x[i]^2 for x[i] in [-10,10]
  xmin = -10.0;
  xmax =  10.0;
  for(int i=0;i<n;i++){
    x[i] = xmin + (xmax-xmin)*i/realT(n-1);
    y[i] = x[i]*x[i];
  }

  //================================================================
  // Example 1:
  //================================================================

  func.create_spline(x,y,n,1e30,1e30,1,"Test spline 1");

  std::cout << "Test of spline object for y = x^2" << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/realT(n-1);
    ynow  = func(xnow); 
    dynow = func.dfdx(xnow); 
    print(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  //================================================================
  // Example 2: As pointer
  //================================================================

  gfunc = new Spline(x,y,n,1e30,1e30,1,"Test spline 2");

  std::cout << "Test of spline pointer for y = x^2 | delta_y = y-x^2" << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/realT(n-1);
    ynow  = gfunc[0](xnow);
    dynow = gfunc->dfdx(xnow);
    print(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  //================================================================
  // Example 3: If boundary condition for dy/dx is known we can 
  // make a more accurate spline
  //================================================================

  hfunc.create_spline(x,y,n,2*(-10.0),2*(10.0),1,"Test spline 1");

  std::cout << "Test of spline object for y = x^2 with the correct BC (y'=2x) | delta_y = y-x^2 " << std::endl;
  for(int i=0;i<n-1;i++){
    xnow  = xmin + (xmax-xmin)*(i+0.5)/realT(n-1)+1.0;
    ynow  = hfunc(xnow);
    dynow = hfunc.dfdx(xnow);
    print(xnow,ynow,dynow);
  }
  std::cout << std::endl;

  // Clean up
  delete gfunc;
  delete[] x,y;
}

int main(int argc, char **argv){
  make_spline();
}
