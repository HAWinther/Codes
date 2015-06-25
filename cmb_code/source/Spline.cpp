
//=======================================================
// Qubic spline methods
//=======================================================

#include "Spline.h"

Spline::Spline() : n(0), linear_x(false), quadratic_x(false), three_regions_x(false), log_x(false), four_regions_x(false) {}

// Constructor for 1D spline from array
Spline::Spline(VecDoub &x1, VecDoub &y1, double yp1_in, double ypn_in) : n(x1.n_elem), 
  yp1(yp1_in), ypn(ypn_in), linear_x(false), quadratic_x(false), three_regions_x(false), log_x(false), four_regions_x(false) {

    x  = x1;
    y  = y1;
    y2 = VecDoub(n);

    make_spline();
  }

// Constructor for spline from matrix
Spline::Spline(VecDoub &x1, MatDoub &yy1, double yp1_in, double ypn_in, bool first1) : 
  first(first1), yp1(yp1_in), ypn(ypn_in), linear_x(false), quadratic_x(false), three_regions_x(false), log_x(false), four_regions_x(false) {
    nx = yy1.nx;
    ny = yy1.ny;
    x  = x1;
    yy = yy1;
    yy2  = MatDoub(nx,ny);

    // Make spline of each row (collumn) in the matrix
    if (first) {
      n = nx;

      for (int i = 0; i<ny ; i++) {
        make_spline_matrix(i);
      }

    } else {
      n = ny;

      for (int i = 0; i<nx ; i++) {
        make_spline_matrix(i);
      }
    }
  }

// Let spline object know that the x-array is lineary spaced
void Spline::set_linear_x(){
  if (n>0) {
    linear_x = true;
    x_start = x[0];
    x_end = x[n-1];
  }
}

// Let spline object know that the x-array is quadratic spaced
void Spline::set_quadratic_x(){
  if (n>0) {
    quadratic_x = true;
    x_start = x[0];
    x_end = x[n-1];
  }
}

// Let spline object know that the x-array is linear in three regions (with different points in each)
void Spline::set_three_regions_x(int np_1, int np_2){
  if (n>0) {
    three_regions_x = true;
    np1		  = np_1;           // Points in region 1
    np2		  = np_2;           // Points in region 2
    np3		  = n - np1 - np2;  // Points in region 3
    x_start   = x[0];         // Start value of region 1
    x_start_reg2 = x[np1];    // Start value of region 2
    x_start_reg3 = x[np1+np2];// Start value of region 3
    x_end		 = x[n-1];        // End value
  }
}

// Let spline object know that the x-array is linear in three regions (with different points in each)
void Spline::set_four_regions_x(int np_1, int np_2, int np_3){
  if (n>0) {
    four_regions_x = true;
    np1		  = np_1;					        // Points in region 1
    np2		  = np_2;					        // Points in region 2
    np3	    = np_3;					        // Points in region 3
    np4		  = n - np1 - np2 - np3;  // Points in region 4
    x_start = x[0];                 // Start value of region 1
    x_start_reg2 = x[np1];          // Start value of region 2
    x_start_reg3 = x[np1+np2];      // Start value of region 3
    x_start_reg4 = x[np1+np2+np3];  // Start value of region 4
    x_end		     = x[n-1];          // End value
  }
}

// Log spaced
void Spline::set_log_x(){
  if (n>0) {
    log_x   = true;
    x_start = x[0];
    x_end   = x[n-1];
  }
}

// Make the spline
void Spline::make_spline(){
  VecDoub u(n);
  double sig, p;

  // Natural spline with y'' = 0 at endpoint 1 if yp1>=1.0e30
  if (yp1 > 0.99e30){
    y2[0] = 0.0;
    u[0]  = 0.0;
  } else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  // Tridiagonal matrix decomposition loop
  for (int i=1;i<n-1;i++){
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  // Natural spline with y'' = 0 at endpoint n if ypn>=1.0e30
  if (ypn > 0.99e30){
    y2[n-1] = 0.0;
  } else {
    double un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    y2[n-1] = (un-0.5*u[n-2])/(0.5*y2[n-2]+1.0);
  }

  // Backsubstitution loop of tridiagonal algorithm
  for (int i = n-2; i>=0; i--) {
    y2[i] = y2[i]*y2[i+1]+u[i];
  }
};

void Spline::make_spline_matrix(int j){
  double sig, p;

  if (first) {

    // Spline over first index : y[0...nx,j]

    VecDoub u(n);

    // Natural spline with y'' = 0 at endpoint 1 if yp1>=1.0e30
    if (yp1 > 0.99e30){
      yy2(0,j) = 0.0;
      u[0]  = 0.0;
    } else {
      yy2(0,j) = -0.5;
      u[0]  = (3.0/(x[1]-x[0]))*((yy(1,j)-yy(0,j))/(x[1]-x[0])-yp1);
    }

    // Tridiagonal matrix decomposition loop
    for (int i=1;i<n-1;i++){
      sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p = sig*yy2(i-1,j)+2.0;
      yy2(i,j) = (sig-1.0)/p;
      u[i] = (yy(i+1,j)-yy(i,j))/(x[i+1]-x[i]) - (yy(i,j)-yy(i-1,j))/(x[i]-x[i-1]);
      u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    // Natural spline with y'' = 0 at endpoint n if ypn>=1.0e30
    if (ypn > 0.99e30){
      yy2(n-1,j) = 0.0;
    } else {
      double un = (3.0/(x[n-1]-x[n-2]))*(ypn-(yy(n-1,j)-yy(n-2,j))/(x[n-1]-x[n-2]));
      yy2(n-1,j) = (un-0.5*u[n-2])/(0.5*yy2(n-2,j)+1.0);
    }

    // Backsubstitution loop of tridiagonal algorithm
    for (int i = n-2; i>=0; i--) {
      yy2(i,j) = yy2(i,j)*yy2(i+1,j)+u[i];
    }

  } else {

    // Spline over second index : y[j, 0...nx]

    VecDoub u(n);

    // Natural spline with y'' = 0 at endpoint 1 if yp1>=1.0e30
    if (yp1 > 0.99e30){
      yy2(0,j) = 0.0;
      u[0]  = 0.0;
    } else {
      yy2(0,j) = -0.5;
      u[0]  = (3.0/(x[1]-x[0]))*((yy(j,1)-yy(j,0))/(x[1]-x[0])-yp1);
    }

    // Tridiagonal matrix decomposition loop
    for (int i=1;i<n-1;i++){
      sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p = sig*yy2(j,i-1)+2.0;
      yy2(j,i) = (sig-1.0)/p;
      u[i] = (yy(j,i+1)-yy(j,i))/(x[i+1]-x[i]) - (yy(j,i)-yy(j,i-1))/(x[i]-x[i-1]);
      u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }

    // Natural spline with y'' = 0 at endpoint n if ypn>=1.0e30
    if (ypn > 0.99e30){
      yy2(j,n-1) = 0.0;
    } else {
      double un = (3.0/(x[n-1]-x[n-2]))*(ypn-(yy(j,n-1)-yy(j,n-2))/(x[n-1]-x[n-2]));
      yy2(j,n-1) = (un-0.5*u[n-2])/(0.5*yy2(j,n-2)+1.0);
    }

    // Backsubstitution loop of tridiagonal algorithm
    for (int i = n-2; i>=0; i--) {
      yy2(j,i) = yy2(j,i)*yy2(j,i+1)+u[i];
    }
  }
};


double Spline::get_spline(double x0){
  int klo, khi, k;
  double h,b,a;

  if (linear_x) {
    klo = min(int((x0 - x_start)/(x_end-x_start)*(n-1)),n-2);
  } else if (three_regions_x) {

    if (x0 < x_start_reg2) {
      klo = min(int((x0 - x_start)/(x_start_reg2-x_start)*(np1-1)),np1-2);
    } else if (x0 < x_start_reg3) {
      klo = np1 + min(int((x0 - x_start_reg2)/(x_start_reg3-x_start_reg2)*(np2-1)),np2-2);
    } else {
      klo = np1 + np2 + min(int((x0 - x_start_reg3)/(x_end-x_start_reg3)*(np3-1)),np3-2);
    }
  } else if (four_regions_x) {

    if (x0 < x_start_reg2) {
      klo = min(int((x0 - x_start)/(x_start_reg2-x_start)*(np1-1)),np1-2);
    } else if (x0 < x_start_reg3) {
      klo = np1 + min(int((x0 - x_start_reg2)/(x_start_reg3-x_start_reg2)*(np2-1)),np2-2);
    } else if (x0 < x_start_reg4) {
      klo = np1 + np2 + min(int((x0 - x_start_reg3)/(x_start_reg4-x_start_reg3)*(np3-1)),np3-2);
    } else {
      klo = np1 + np2 + np3 + min(int((x0 - x_start_reg4)/(x_end-x_start_reg4)*(np4-1)),np4-2);
    }

  } else if (log_x) {
    klo = min(int(log(x0/x_start)/log(x_end/x_start) * (n-1)),n-2);
  } else if (quadratic_x) {
    klo = min(int(sqrt((x0 - x_start)/(x_end-x_start))*(n-1)),n-2);
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
  khi = klo + 1;

  h = x[khi]-x[klo];
  if (h == 0.0){
    cout << "The x-values must be distict. Exiting." << endl;
    exit(1);
  }
  a = (x[khi]-x0)/h;
  b = (x0-x[klo])/h;

  return a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;
};

double Spline::get_spline_matrix(double x0, int j){
  int klo, khi, k;
  double h,b,a;

  if (linear_x) {
    klo = min(int((x0 - x_start)/(x_end-x_start) * (n-1)),n-2);
  } else if (log_x) {
    klo = min(int(log(x0/x_start)/log(x_end/x_start) * (n-1)),n-2);
  } else if (three_regions_x) {
    if (x0 < x_start_reg2) {
      klo = min(int((x0 - x_start)/(x_start_reg2-x_start)*(np1-1)),np1-2);
    } else if (x0 < x_start_reg3) {
      klo = np1 + min(int((x0 - x_start_reg2)/(x_start_reg3-x_start_reg2)*(np2-1)),np2-2);
    } else {
      klo = np1 + np2 + min(int((x0 - x_start_reg3)/(x_end-x_start_reg3)*(np3-1)),np3-2);
    }
  } else if (four_regions_x) {
    if (x0 < x_start_reg2) {
      klo = min(int((x0 - x_start)/(x_start_reg2-x_start)*(np1-1)),np1-2);
    } else if (x0 < x_start_reg3) {
      klo = np1 + min(int((x0 - x_start_reg2)/(x_start_reg3-x_start_reg2)*(np2-1)),np2-2);
    } else if (x0 < x_start_reg4) {
      klo = np1 + np2 + min(int((x0 - x_start_reg3)/(x_start_reg4-x_start_reg3)*(np3-1)),np3-2);
    } else {
      klo = np1 + np2 + np3 + min(int((x0 - x_start_reg4)/(x_end-x_start_reg4)*(np4-1)),np4-2);
    }
  } else if (quadratic_x) {
    klo = min(int(sqrt((x0 - x_start)/(x_end-x_start)) * (n-1)),n-2);
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
  khi = klo + 1;

  if (first) {
    h = x[khi]-x[klo];
    if (h == 0.0){
      cout << "The x-values must be distict. Exiting." << endl;
      exit(1);
    }
    a = (x[khi]-x0)/h;
    b = (x0-x[klo])/h;

    return a*yy(klo,j)+b*yy(khi,j)+((a*a*a-a)*yy2(klo,j)+(b*b*b-b)*yy2(khi,j))*(h*h)/6.0;
  } else {		
    h = x[khi]-x[klo];
    if (h == 0.0){
      cout << "The x-values must be distict. Exiting." << endl;
      exit(1);
    }
    a = (x[khi]-x0)/h;
    b = (x0-x[klo])/h;

    return a*yy(j,klo)+b*yy(j,khi)+((a*a*a-a)*yy2(j,klo)+(b*b*b-b)*yy2(j,khi))*(h*h)/6.0;
  }

};

void Spline::extract_y2(VecDoub_O &y1){
  if (n != y1.n_elem)
    cout << "Error in Spline::extract_y2 : Size of arrays do not match." << endl;
  y1 = y2;
}

void Spline::extract_yy2(MatDoub_O &yy1){
  if (nx != yy1.nx && ny != yy1.ny)
    cout << "Error in Spline::extract_yy2 : Size of arrays do not match." << endl;
  yy1 = yy2;
}

void Spline::set_external_spline(VecDoub_I &x_in, MatDoub_I &yy_in, MatDoub_I &yy2_in, bool first1){
  first = first1;
  yy2 = yy2_in;
  yy  = yy_in;
  x   = x_in;
  nx  = yy_in.nx;
  ny  = yy_in.ny;
  n   = first ? nx : ny;
}
