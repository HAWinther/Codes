#ifndef ODESOLVER_INC
#define ODESOLVER_INC

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

//=======================================================
// Ode Solver Class
//=======================================================

class OdeSolver{
  private:

    // Constants
    static const int MAXSTP = 1000000;

    static const double SAFETY = 0.9;
    static const double PGROW  = -0.2;
    static const double PSHRNK = -0.25;
    static const double ERRCON = 1.89e-4;

    static const double a2  = (0.2);
    static const double a3  = (0.3);
    static const double a4  = (0.6);
    static const double a5  = (1.0);
    static const double a6  = (0.8750);
    static const double c1  = (37.0/378.0);
    static const double c3  = (250.0/621.0);
    static const double c4  = (125.0/594.0);
    static const double c6  = (512.0/1771.0);
    static const double dc1 = ((37.0/378.0)  -  2825.0/27648.0);
    static const double dc3 = ((250.0/621.0) -  18575.0/48384.0);
    static const double dc4 = ((125.0/594.0) -  13525.0/55296.0);
    static const double dc5 = (-277.0/14336.0);
    static const double dc6 = ((512.0/1771.0) -  0.250);
    static const double b21 = (0.2);
    static const double b31 = (3.0/40.0);
    static const double b32 = (9.0/40.0);
    static const double b41 = (0.3);
    static const double b42 = (-0.9); 
    static const double b43 = (1.2);
    static const double b51 = (-11.0/54.0); 
    static const double b52 = (2.5);
    static const double b53 = (-70.0/27.0); 
    static const double b54 = (35.0/27.0);
    static const double b61 = (1631.0/55296.0);
    static const double b62 = (175.0/512.0);
    static const double b63 = (575.0/13824.0); 
    static const double b64 = (44275.0/110592.0); 
    static const double b65 = (253.0/4096.0);

    // Storeage vectors
    std::vector<double> y, xstore;
    std::vector< std::vector<double> > ystore;
 
    int nok, nbad, n, inow, neq;
    double h1, hmin, eps, xmin, xmax, xold;
    bool is_ic_set;

    void (* dydx_func)(double, std::vector<double> &, std::vector<double> &);
  public:

    OdeSolver(int n, int neq, void (*derivs)(double, std::vector<double> &, std::vector<double> &)){

      // Number of steps
      this->n = n;

      // Number of equations
      this->neq = neq;

      // Derivative function
      this->dydx_func = derivs;

      // Initialize
      y = std::vector<double> (n,0.0);
      ystore = std::vector< std::vector<double> >(neq);
      for(int i=0;i<neq;i++)
        ystore[i] = std::vector<double>(n,0.0);
      xstore = std::vector<double>(n,0.0);

      // Standard accuracy parameters
      nok  = nbad = 0;
      hmin = 0.0;
      eps  = 1e-15;
      h1   = 1e-5;

      // IC are not set yet
      is_ic_set = false;
    }

    ~OdeSolver(){}

    // Return the x-values of the solution
    std::vector<double> x_array(){
      return xstore;
    }

    // Return the solution
    std::vector<double> y_array(int i){
      return ystore[i];
    }

    // Set initial conditions from an array
    void set_initial_conditions(double xmin, double xmax, std::vector<double> &ic){
      this->xmin = xmin;
      this->xmax = xmax;
      xold = xstore[0] = xmin;
      for(int i=0;i<neq;i++)
        y[i] = ystore[i][0] = ic[i];
      inow = 1;
      is_ic_set = true;
    }

    // Solve the whole system by taking steps until we are done
    void solve(bool verbose){
      while(take_step(verbose));
    }

    // Take one single RK step
    bool take_step(bool verbose){
      double x;

      // Check if IC are set of we are done
      if(is_ic_set == false){
        std::cout << "Error in ode-solver. IC not set!" << std::endl;
        exit(1);
      }	
      if(inow > n-1){
        std::cout << "Warning in ode-solver - i > n" << std::endl;
        return false;
      }

      // Current x-value
      x = xmin + (xmax-xmin)*inow/double(n-1);

      // Integrate
      rk_odeint(y, neq, xold, x, eps, h1, hmin, &nok, &nbad, dydx_func);

      // Verbose
      if(verbose){
        std::cout << "i = " << inow << " / " << n << " x: " << x << " y: ";
        for(int i=0;i<neq;i++) std::cout << y[i] << " ";
        std::cout << std::endl;
      }

      // Store all values
      xstore[inow] = x;
      for(int i=0;i<neq;i++)
        ystore[i][inow] = y[i];

      // Update old values
      xold = x;
      inow++;

      return (inow <= n-1);
    }

    // List of solver methods
    void rkck(std::vector<double> &y, std::vector<double> &dydx, int n, double x, double h, std::vector<double> &yout, std::vector<double> &yerr, void (*derivs)(double, std::vector<double> &, std::vector<double> &));

    void rkqs(std::vector<double> &y, std::vector<double> &dydx, int n, double *x, double htry, double eps, std::vector<double> &yscal, double *hdid, double *hnext, void (*derivs)(double, std::vector<double> &, std::vector<double> &));
    
    void rk_odeint(std::vector<double> &ystart, int nvar, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, void (*derivs)(double, std::vector<double> &, std::vector<double> &));

};

//=======================================================
//  Runge Kutta Methods
//=======================================================

void OdeSolver::rkqs(std::vector<double> &y, std::vector<double> &dydx, int n, double *x, double htry, double eps, std::vector<double> &yscal, double *hdid, double *hnext, void (*derivs)(double, std::vector<double> &, std::vector<double> &)){
  double errmax, h, htemp, xnew;
  std::vector<double> yerr(n,0.0);
  std::vector<double> ytemp(n,0.0);
  int i;

  // Check for NaN
  if(y[1] != y[1]){
    std::cout << "NAN found in rkqs n = " << n << std::endl;
    exit(1);
  }

  // Main loop
  h=htry;
  for (;;) {

    rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);

    errmax = 0.0;

    for (i=0;i<n;i++){
      if (errmax < fabs(yerr[i]/yscal[i]))
        errmax = fabs(yerr[i]/yscal[i]);
    }
    errmax /= eps;

    // Check for error
    if (errmax <= 1.0) break;

    // New stepsize
    htemp = SAFETY*h*pow(errmax,PSHRNK);
    if(h>=0.0){
      if(htemp>0.1*h) h = htemp;
      else h = 0.1*h;
    } else {
      if(htemp<0.1*h) h = htemp;
      else h = 0.1*h;
    }

    // Check for stepsize underflow
    xnew=(*x)+h;
    if (xnew == *x){
      std::cout << "Stepsize underflow in rkqs" << std::endl;
      exit(1);
    }
  }

  // Next stepsize
  if(errmax > ERRCON)
    *hnext=SAFETY*h*pow(errmax,PGROW);
  else
    *hnext=5.0*h;
  *x += (*hdid=h);

  // Copy over solution
  for(i=0;i<n;i++)
    y[i] = ytemp[i];
}

//=======================================================
//  Runge Kutta methods
//=======================================================

void OdeSolver::rkck(std::vector<double> &y, std::vector<double> &dydx, int n, double x, double h, std::vector<double> &yout, std::vector<double> &yerr, void (*derivs)(double, std::vector<double> &, std::vector<double> &)){
  std::vector<double> ak2, ak3, ak4, ak5, ak6, ytemp;
  int i;

  ak2   = std::vector<double>(n,0.0);
  ak3   = std::vector<double>(n,0.0);
  ak4   = std::vector<double>(n,0.0);
  ak5   = std::vector<double>(n,0.0);
  ak6   = std::vector<double>(n,0.0);
  ytemp = std::vector<double>(n,0.0);

  // Runge-Kutta steps
  for (i=0;i<n;i++)
    ytemp[i] = y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=0;i<n;i++)
    yout[i] = y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<n;i++)
    yerr[i] = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

//=======================================================
//  Driver for solving the system
//=======================================================

void OdeSolver::rk_odeint(std::vector<double> &ystart, int nvar, double x1, double x2, double eps, double h1, 
    double hmin, int *nok, int *nbad, void (*derivs)(double, std::vector<double> &, std::vector<double> &)){
  double xsav,x,hnext,hdid,h,dxsav;
  std::vector<double> yscal, y, dydx;
  int nstp,i; 

  dydx  = std::vector<double>(nvar,0.0);
  yscal = std::vector<double>(nvar,0.0);
  y     = ystart;

  // Fixed constants...
  const double TINY = 1.0e-30;

  // Initial x and stepsize
  x = x1;
  h = (h1*(x2-x1) > 0.0) ? h1 : -h1;
  *nok = (*nbad) = 0;

  // Take steps
  for(nstp=1;nstp<=MAXSTP;nstp++) {

    // Get derivatives
    (*derivs)(x,y,dydx); 

    // Take a steps
    for(i=0;i<nvar;i++) yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;

    if((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;

    // Runge-Kutta steps
    rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);

    // Check if step was successful
    if(hdid == h) 
      ++(*nok);
    else 
      ++(*nbad);

    // Check if we are done
    if ((x-x2)*(x2-x1) >= 0.0) {
      // Copy over solution
      for (i=0;i<nvar;i++)
        ystart[i] = y[i];
      return;
    }

    // Check if we have reached precision limit in stepsize
    if (fabs(hnext) <= hmin) {
      std::cout << "Step size too small in odeint" << std::endl;
      exit(1);
    }
    h = hnext;
  } 

  // If we are here then steps > maxsteps
  std::cout << "Too many steps in routine odeint" << std::endl;
  exit(1);
}

#endif
