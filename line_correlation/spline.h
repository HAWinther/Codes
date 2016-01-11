#ifndef  SPLINEHEADER_INC
#define SPLINEHEADER_INC
#include <iostream>
#include <stdlib.h>
#include <math.h>
#define MIN(x,y) ((x)>(y) ? (y) : (x))

//=======================================================
// Spline class
//=======================================================

class Spline {
  private:
		std::string name;
    double *y, *x, *y2;
    double x_start, x_end;
    int n;
    int type;  

  public:
    
    //===============================
    // Type of x-array (for faster
    // look-up)
    // type = 0 : arbritrary x array
    // type = 1 : x is linear
    // type = 2 : x is logaritmic
    //===============================

		Spline() : n(0) {}

		Spline(double *xin, double *yin, int nin, double yp1, double ypn, int typein, std::string namein){
			make_spline(xin,yin,nin,yp1,ypn,typein,namein);
		}

		~Spline(){free_spline();}

		// Operator overloading to allow spline(x) to get the value of the spline
		double operator()(double x0){
			return get_y(x0);
		}

		// Make a spline
    void make_spline(double *x, double *y, int n, double yp1, double ypn, int type, std::string name){
      double sig, p, *u, un;
      int i;

      // Set constants
      this->name = name;
      this->type = type;
      this->n    = n;
      x_start    = x[0];
      x_end      = x[n-1];

      // Allocate memory
      u        = new double[n];
      y2       = new double[n];
      this->y  = new double[n];
      this->x  = new double[n];
      for (i = 0; i<n; i++) {
        this->y[i]  = y[i];
        this->x[i]  = x[i];
        y2[i] = 0.0;
      }

      // Boundary conditions for the spline at left end
      if (yp1 > 0.99e30){
        y2[0] = u[0]  = 0.0;
      } else {
        y2[0] = -0.5;
        u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
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
      if (ypn > 0.99e30){
        y2[n-1] = 0.0;
      } else {
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
        y2[n-1] = (un-0.5*u[n-2])/(0.5*y2[n-2]+1.0);
      }

      // Calculate y''
      for (i = n-2; i>=0; i--) y2[i] = y2[i]*y2[i+1]+u[i];

      delete[] u;
    }
  
    // Get the value from a spline
    double get_y(double x0){
      int klo, khi, k;
      double h, b, a, result;

			if(x0>x_end)
				std::cout << "Spline error: x0 = " << x0 << ">" << x_end << std::endl;
			if(x0<x_start)
				std::cout << "Spline error: x0 = " << x0 << "<" << x_start << std::endl;

      // Find closest x-value to x0
      if (type == 1) {
        klo = MIN(int((x0 - x_start)/(x_end-x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else if (type == 2) {
        klo = MIN(int((log(x0/x_start))/log(x_end/x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else {
        klo = 0;
        khi = n-1;
        while (khi-klo>1) {
          k = (khi+klo) >> 1;
          if (x[k]>x0) {
            khi = k;
          }else {
            klo = k;
          }
        }
      }
      h = x[khi]-x[klo];
      if (h == 0.0){
				std::cout << "Spline error: h=0 in " << name << ". The x-values must be distict. Exiting." << std::endl;
        exit(1);
      }
      a = (x[khi]-x0)/h;
      b = (x0-x[klo])/h;
      result = (a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
      return result;
    }

    // Get the derivative of the quantity splines
		double get_dy(double x0){
      int klo, khi, k;
      double h, b, a, result;

			if(x0>x_end)
				std::cout << "Spline error: x0 = " << x0 << ">" << x_end << std::endl;
			if(x0<x_start)
				std::cout << "Spline error: x0 = " << x0 << "<" << x_start << std::endl;

      // Find closest x-value to x0
      if (type == 1) {
        klo = MIN(int((x0 - x_start)/(x_end-x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else if (type == 2) {
        klo = MIN(int((log(x0/x_start))/log(x_end/x_start)*(n-1)),n-2);
        khi = klo + 1;
      } else {
      	// Binary search
        klo = 0;
        khi = n-1;
        while (khi-klo>1) {
          k = (khi+klo) >> 1;
          if (x[k]>x0) {
            khi = k;
          }else {
            klo = k;
          }
        }
      }
      h = x[khi]-x[klo];
      if (h == 0.0){
				std::cout << "Spline error: h=0 in " << name << ". The x-values must be distict. Exiting." << std::endl;
        exit(1);
      }
      a = (x[khi]-x0)/h;
      b = (x0-x[klo])/h;
      result = (y[khi]-y[klo])/h + h/6.0*(-(3*a*a-1)*y2[klo] + (3*b*b-1)*y2[khi]);
      return result;
    }

    // Delete all arrays in the spline
    void free_spline(){
    	if(x!=NULL){
      	delete[] x;
      	delete[] y;
      	delete[] y2;
			}
    }
};

#endif
