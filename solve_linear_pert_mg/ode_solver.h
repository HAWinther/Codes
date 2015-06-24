#ifndef ODESOLVER_INC
#define ODESOLVER_INC

#include <stdlib.h>
#include <math.h>
#include <iostream>

// Type to be used
#if defined(DOUBLE)
typedef double realT;
#elif defined(LONGDOUBLE)
typedef long double realT;
#elif defined(FLOAT)
typedef float realT;
#else
typedef double realT;
#endif

// Some constants needed for the solver
#define a2      ((realT) 0.2e0)
#define a3      ((realT) 0.3e0)
#define a4      ((realT) 0.6e0)
#define a5      ((realT) 1.e0)
#define a6      ((realT) 0.875e0)
#define c1      ((realT) 37.e0/378.e0)
#define c3      ((realT) 250.e0/621.e0)
#define c4      ((realT) 125.e0/594.e0)
#define c6      ((realT) 512.e0/1771.e0)
#define dc1     (c1 -  (realT) 2825.e0/27648.e0)
#define dc3     (c3 -  (realT) 18575.e0/48384.e0)
#define dc4     (c4 -  (realT) 13525.e0/55296.e0)
#define dc5     ((realT) -277.e0/14336.e0)
#define dc6     (c6 -  (realT) 0.25e0)
#define MAXSTP 1000000
#define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkck(realT *y, realT *dydx, int n, realT x, realT h, realT *yout, realT *yerr, void (*derivs)(realT, realT *, realT *));

///////////////////////////////////////////////////
//  Runge Kutta Method
///////////////////////////////////////////////////

void rkqs(realT *y, realT *dydx, int n, realT *x, realT htry, realT eps, realT *yscal, realT *hdid, realT *hnext, void (*derivs)(realT, realT *, realT *)){
  int i;
  realT errmax,h,htemp,xnew,*yerr,*ytemp;

  // Check for NaN
  if(y[1] != y[1]){
    std::cout << "NAN found in rkqs n = " << n << std::endl;
    exit(1);
  }

  // Allocate
  yerr  = new realT[n];
  ytemp = new realT[n];

  // Initialize
  for(i=0;i<n;i++){
    yerr[i] = 0.0;
    ytemp[i] = 0.0;
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

  // Clean up
  delete[] ytemp;
  delete[] yerr;
}

///////////////////////////////////////////////////
//  Runge Kutta method
///////////////////////////////////////////////////

void rkck(realT *y, realT *dydx, int n, realT x, realT h, realT *yout, realT *yerr, void (*derivs)(realT, realT *, realT *)){
  realT *ak2, *ak3, *ak4, *ak5, *ak6, *ytemp;
  realT b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9, b43=1.2, 
        b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0, b54=35.0/27.0, b61=1631.0/55296.0, 
        b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0;
  int i;

  // Allocate memory
  ak2 = new realT[n];
  ak3 = new realT[n];
  ak4 = new realT[n];
  ak5 = new realT[n];
  ak6 = new realT[n];
  ytemp = new realT[n];

  // Runge-Kutta steps
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
  for (i=0;i<n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

  // Clean up
  delete[] ytemp;
  delete[] ak6;
  delete[] ak5;
  delete[] ak4;
  delete[] ak3;
  delete[] ak2;
}


///////////////////////////////////////////////////
//  Integrate coupled ODE using RK4
///////////////////////////////////////////////////

void rk_odeint(realT *ystart, int nvar, realT x1, realT x2, realT eps, realT h1, 
    realT hmin, int *nok, int *nbad, void (*derivs)(realT, realT *, realT *)){
  int nstp,i; 
  realT xsav,x,hnext,hdid,h; 
  realT *yscal, *y, *dydx;
  realT dxsav;

  // Allocate memory
  yscal = new realT[nvar];
  y     = new realT[nvar];
  dydx  = new realT[nvar];

  // Initial x
  x = x1;

  // Stepsize
  if(h1*(x2-x1) > 0.0)
    h = h1;
  else
    h = -h1;
  *nok = (*nbad) = 0;

  // Make copy of y
  for(i=0;i<nvar;i++)
    y[i] = ystart[i];

  // Take steps
  for(nstp=1;nstp<=MAXSTP;nstp++) {

    // Get derivatives
    (*derivs)(x,y,dydx); 

    // Take a steps
    for(i=0;i<nvar;i++)
      yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;

    if((x+h-x2)*(x+h-x1) > 0.0)
      h = x2-x;

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
      // Clean up
      delete[] dydx;
      delete[] y;
      delete[] yscal;
      return;
    }

    // Check if we have reached precision limit in stepsize
    if (fabs(hnext) <= hmin) {
      std::cout << "Step size too small in odeint" << std::endl;
      exit(1);
    }
    h = hnext;
  } 

  // Clean up
  delete[] dydx;
  delete[] y;
  delete[] yscal;

  // If we are here then steps > maxsteps
  std::cout << "Too many steps in routine odeint" << std::endl;
  exit(1);
}

///////////////////////////////////////
// Ode Solver Class
///////////////////////////////////////

class OdeSolver{
  private:
    realT *y, **ystore, *xstore;
    int nok, nbad, n, inow, neq;
    realT h1, hmin, eps, xmin, xmax, xold;
    void (* dydx_func)(realT, realT *, realT *);
    bool ic_set;
  public:

    OdeSolver(int n, int neq, void (*derivs)(realT, realT *, realT *)){

      // Number of steps
      this->n = n;

      // Number of equations
      this->neq = neq;

      // Derivative function
      this->dydx_func = derivs;

      // Allocate
      y = new realT[neq];
      ystore = new realT*[neq];
      for(int i=0;i<neq;i++)
        ystore[i] = new realT[n];
      xstore = new realT[n];

      // Standard parameters
      nok  = nbad = 0;
      hmin = 0.0;
      eps  = 1e-15;
      h1   = 1e-5;

      // IC are not set yet
      ic_set = false;
    }

    ~OdeSolver(){
      delete[] y;
      for(int i=0;i<neq;i++)
        delete[] ystore[i];
      delete[] ystore;
      delete[] xstore;
      y = NULL;
      ystore = NULL;
      xstore = NULL;
    }

    realT *x_array(){return xstore;}
    realT *y_array(int i){return ystore[i];}

    void set_initial_conditions(realT xmin, realT xmax, realT *ic){
      this->xmin = xmin;
      this->xmax = xmax;
      xold = xstore[0] = xmin;
      for(int i=0;i<neq;i++)
        y[i] = ystore[i][0] = ic[i];
      inow = 1;
      ic_set = true;
    }

    void solve(bool verbose){
      while(take_step(verbose));
    }

    bool take_step(bool verbose){
      realT x;

      // Check if IC are set of we are done
      if(ic_set == false){
        std::cout << "Error in ode-solver. IC not set!" << std::endl;
        exit(1);
      }	
      if(inow > n-1){
        std::cout << "Warning in ode-solver - i > n" << std::endl;
        return false;
      }

      // Current x-value
      x = xmin + (xmax-xmin)*inow/realT(n-1);

      // Integrate
      rk_odeint(y, neq, xold, x, eps, h1, hmin, &nok, &nbad, dydx_func);

      // Verbose
      if(verbose){
        std::cout << "i = " << inow << " / " << n << " x: " << x << " y: ";
        for(int i=0;i<neq;i++) std::cout << y[i] << " ";
        std::cout << std::endl;
      }

      // Store values
      xstore[inow] = x;
      for(int i=0;i<neq;i++)
        ystore[i][inow] = y[i];

      // Update old values
      xold = x;
      inow++;

      return (inow <= n-1);
    }
};
#endif
