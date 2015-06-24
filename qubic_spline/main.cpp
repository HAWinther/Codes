#include<iomanip>
#include "spline.h"

///////////////////////////////////
// Make and test the spline for
// the function y = x^2
///////////////////////////////////

void make_spline(){
  Spline func;
  realT *x, *y;
  realT xmin, xmax, xnow;
  int n;

  std::cout << "==================" << std::endl;
  std::cout << "   Make Spline   " << std::endl;
  std::cout << "==================" << std::endl;

  // Number of points in array
  n = 25;

  // Initialize arrays 
  x = new realT[n];
  y = new realT[n];

  xmin = -10.0;
  xmax =  10.0;

  // Make y[i] = x[i]^2 for x[i] in [-10,10]
  for(int i=0;i<n;i++){
    x[i] = xmin + (xmax-xmin)*i/realT(n-1);
    y[i] = x[i]*x[i];
  }

  // Create the spline
  func.make_spline(x,y,n,1e30,1e30,1,"Test spline");

  // Test to see how accurate the spline is
  for(int i=0;i<n-1;i++){
    xnow = xmin + (xmax-xmin)*(i+0.5)/realT(n-1);
    std::cout << " x: " << std::setw(12) << xnow << " y: " << std::setw(12) << func.get_spline(xnow);
    std::cout << " delta_y: " << std::setw(12) << func.get_spline(xnow) - xnow*xnow << std::endl;
  }
}

int main(int argc, char **argv){
  make_spline();
}
