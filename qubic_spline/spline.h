#ifndef SPLINEHEADER_INC
#define SPLINEHEADER_INC
#include <iostream> 
#include <stdlib.h> 
#include <vector>   
#include <math.h>

//====================================================
//
// Class to construct a Qubic Spline of a function
//
// n is number of datapoints
//
// Assumes monotone x-array: 
//    x_start = x[0] < x[1] < ... < x[n-1] = x_end
//
// If x-array is regulary spaced we use direct lookup
// to save time:
//   type = 0 : arbritrary (binary search used)
//   type = 1 : linear spacing
//   type = 2 : logaritmic spacing
//
// Boundary conditions:
//   dydx1, dydxn are dy/dx at the two boundary points
//   Use >0.99e30 to get the so-called natural spline
//
//=====================================================

template <class T> 
class Spline {
  private:
    std::string name;
    std::vector<T> y, x, y2;
    T x_start, x_end;
    T dydx1, dydxn;
    int n, type;   

  public:

    Spline(): n(0), type(0), x_start(0.0), x_end(0.0), name("") {}
    ~Spline(){}

    // Initialized x,y from T-array
    Spline(T *x, T *y, int n, T dydx1, T dydxn, int type, std::string name) : 
      n(n), type(type), name(name), x_start(x[0]), x_end(x[n-1]), dydx1(dydx1), dydxn(dydxn),
      x(std::vector<T>(x,x+n)), y(std::vector<T>(y,y+n)), y2(std::vector<T>(n,0.0)){
        create_spline();
    }

    // Initialize x,y from T-vector
    Spline(std::vector<T> &x, std::vector<T> &y, int n, T dydx1, T dydxn, int type, std::string name) : 
      n(n), type(type), name(name), x_start(x[0]), x_end(x[n-1]), dydx1(dydx1), dydxn(dydxn),
      x(x), y(y), y2(std::vector<T>(n,0.0)){
        create_spline();
    }

    // Assignment operator to allow for 'myspline(x)' useage
    T operator()(const T& x){
      return f(x);
    }

    //=====================================================
    // Calculate index klo such that x[klo] <= x0 < x[klo+1]
    //=====================================================

    inline int lower_index(double x0){
      int klo, khi, k;
      if (type == 1) {
        klo = std::min(int((x0 - x_start)/(x_end-x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else if (type == 2) {
        klo = std::min(int((log(x0/x_start))/log(x_end/x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else {
        klo = 0;
        khi = n-1;
        while (khi-klo>1) {
          k = (khi+klo) >> 1;
          if (x[k]>x0) {
            khi = k;
          } else {
            klo = k;
          }
        }
      }
      return klo;
    }

    //=====================================================
    // Calculate the spline
    //=====================================================

    void create_spline(){
      T sig, p, un;
      std::vector<T> u(n);
      int i;

      // Boundary conditions for the spline at left end
      if (dydx1 > 0.99e30){
        y2[0] = u[0]  = 0.0;
      } else {
        y2[0] = -0.5;
        u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-dydx1);
      }

      // Create spline by solving recurence relation
      for (i=1;i<n-1;i++){
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p = sig*y2[i-1]+2.0;
        y2[i] = (sig-1.0)/p;
        u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
      }

      // Boundary condition for the spline at right end
      if (dydxn > 0.99e30){
        y2[n-1] = 0.0;
      } else {
        un = (3.0/(x[n-1]-x[n-2]))*(dydxn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
        y2[n-1] = (un-0.5*u[n-2])/(0.5*y2[n-2]+1.0);
      }

      // Calculate y''
      for (i = n-2; i>=0; i--) y2[i] = y2[i]*y2[i+1]+u[i];
    }

    //=====================================================
    // Extract function value from the spline.
    // If x0 is outside range return the closest value.
    //=====================================================

    T f(T x0){
      int klo, khi;
      T h, b, a, result;

      // Calculate x[klo] < x0 < x[khi]
      klo = lower_index(x0);
      khi = std::min(klo+1,n-1);
      h = x[khi]-x[klo];

      // Check for error
      if (h == 0.0){
        std::cout << "Error in Spline<" << name << "> f(x); h = 0; x-values must be distict!" << std::endl;
        exit(1);
      }

      // Calculate interpolation value
      a = (x[khi]-x0)/h;
      b = (x0-x[klo])/h;
      result = (a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
      return result;
    }
    
    //=====================================================
    // Extract the derivative of the splined function.
    // If x0 is outside range return the closest value.
    //=====================================================

    T dfdx(T x0){
      int klo, khi;
      T h, b, a, result;

      // Calculate x[klo] < x0 < x[khi]
      klo = lower_index(x0);
      khi = std::min(klo+1,n-1);
      h = x[khi]-x[klo];

      // Check for error
      if (h == 0.0){
        std::cout << "Error in Spline<" << name << "> dfdx(x); h = 0; x-values must be distict!" << std::endl;
        exit(1);
      }

      // Calculate interpolation value
      a = (x[khi]-x0)/h;
      b = (x0-x[klo])/h;
      result = (y[khi]-y[klo])/h + h/6.0*(-(3*a*a-1)*y2[klo] + (3*b*b-1)*y2[khi]);
      return result;
    }

    //=====================================================
    // Clean up memory if needed
    //=====================================================

    void clean(){
      n = type = 0; name = "";
      x.clear(); y.clear(); y2.clear();
    }
};

#endif
