#ifndef ODESOLVER_INC
#define ODESOLVER_INC
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>


class OdeSolver{

  //=======================================================
  // 
  // Code to solve coupled systems of ODEs
  // Methods described in Numerical Recipies
  //
  // We solve the ODE y_vec' = f_vec(x,y_vec) where y is
  // a [neq] long vector of variables to solve for and [f_vec] 
  // is the RHS of the system provided by the user in a 
  // function.
  //
  // The ODE is solved from [xmin] to [xmax] and stored at 
  // [n] lineary spaced points
  //
  // When done the solution (x_k, y_j(x_k)) is given as
  //       y_j(x_k) => y_array(j)[k]
  //       x_k = x_array[k]
  //
  //       where x_k = xmin + k * (xmax-xmin)/(n-1)
  //
  //=======================================================

  private:

    // Hardcoded constants
    static constexpr int    MAXSTP = 1000000;
    static constexpr double TINY   = (1.0e-30);
    static constexpr double SAFETY = (0.9);
    static constexpr double PGROW  = (-0.2);
    static constexpr double PSHRNK = (-0.25);
    static constexpr double ERRCON = (1.89e-4);
    static constexpr double a2     = (0.2);
    static constexpr double a3     = (0.3);
    static constexpr double a4     = (0.6);
    static constexpr double a5     = (1.0);
    static constexpr double a6     = (0.8750);
    static constexpr double c1     = (37.0/378.0);
    static constexpr double c3     = (250.0/621.0);
    static constexpr double c4     = (125.0/594.0);
    static constexpr double c6     = (512.0/1771.0);
    static constexpr double dc1    = (37.0/378.0  - 2825.0/27648.0);
    static constexpr double dc3    = (250.0/621.0 - 18575.0/48384.0);
    static constexpr double dc4    = (125.0/594.0 - 13525.0/55296.0);
    static constexpr double dc5    = (-277.0/14336.0);
    static constexpr double dc6    = ((512.0/1771.0) - 0.250);
    static constexpr double b21    = (0.2);
    static constexpr double b31    = (3.0/40.0);
    static constexpr double b32    = (9.0/40.0);
    static constexpr double b41    = (0.3);
    static constexpr double b42    = (-0.9); 
    static constexpr double b43    = (1.2);
    static constexpr double b51    = (-11.0/54.0); 
    static constexpr double b52    = (2.5);
    static constexpr double b53    = (-70.0/27.0); 
    static constexpr double b54    = (35.0/27.0);
    static constexpr double b61    = (1631.0/55296.0);
    static constexpr double b62    = (175.0/512.0);
    static constexpr double b63    = (575.0/13824.0); 
    static constexpr double b64    = (44275.0/110592.0); 
    static constexpr double b65    = (253.0/4096.0);

    int n;    // Number of x-points     
    int neq;  // Number of equations

    std::vector<double> y;
    std::vector<double> dydx;
    std::vector<double> xstore;
    std::vector<std::vector<double> > ystore;

    int nok;        // Number of accepted steps
    int nbad;       // Number of rejected steps

    double hmin;    // Minimum step-size
    double eps;     // Accuracy control
    double h1;      // First step-size to try
    double hdid;    // Step-size of old step

    bool is_ic_set; // Is IC set or not?

    // The right hand side of the ODE system: y_vec' = f_vec(x,y_vec)
    void (* dydx_func)(const double &, const std::vector<double>&, std::vector<double>&);
    
    double xmin;    // Startvalue of x
    double xmax;    // Endvalue of x
    double xold;    // Current startvalue of x
    int inow;       // Current index in x-array we are at
    
  public:

    OdeSolver() : 
          n(0), 
          neq(0),
          y(0),
          dydx(0),
          xstore(0),
          ystore(0),
          nok(0),
          nbad(0),
          hmin(0.0),
          eps(1e-8),
          h1(1e-4),
          hdid(0.0),
          is_ic_set(false) {}
    
    OdeSolver(int _n, int _neq, void (*rhs_func)(const double &, const std::vector<double> &, std::vector<double> &)) : 
          n(_n), 
          neq(_neq), 
          y(std::vector<double> (n, 0.0)), 
          dydx(std::vector<double> (n, 0.0)), 
          xstore(std::vector<double> (n, 0.0)), 
          ystore(std::vector< std::vector<double> >(neq, std::vector<double>(n, 0.0))), 
          nok(0), 
          nbad(0), 
          hmin(0.0), 
          eps(1e-8), 
          h1(1e-4), 
          hdid(0.0),
          is_ic_set(false),
          dydx_func(rhs_func) {}

    ~OdeSolver(){}

    //=======================================================
    // Change the precision parameters
    //=======================================================
    void set_precision(double _eps, double _h1, double _hmin){
      eps  = _eps;
      hmin = _hmin;
      h1   = _h1;
    }

    //=======================================================
    // Extract data
    //=======================================================
    std::vector<double> get_x(){ return xstore; }
    std::vector<double> get_y(int i = 0){ return ystore[i]; }

    //=======================================================
    // Set initial conditions from an array
    //=======================================================
    void set_ic(double _xmin, double _xmax, std::vector<double>& ic){
      xmin = _xmin;
      xmax = _xmax;
      xold = xstore[0] = _xmin;
      for(int i = 0; i < neq; i++)
        y[i] = ystore[i][0] = ic[i];
      inow = 1;
      is_ic_set = true;
    }

    //=======================================================
    // Solve the whole system by taking steps until we are done
    //=======================================================
    void solve(bool verbose){
      while(solve_one(verbose));
    }

    //=======================================================
    // Propagate the ODE from current xvalue to the next
    //=======================================================
    bool solve_one(bool verbose){
      double x;

      // Check if IC are set of we are done
      if(is_ic_set == false)
        report_error("Error in ode-solver. IC not set!",1);
      if(inow > n-1){
        report_error("Warning in ode-solver - i > n",0);
        return false;
      }

      // Current x-value
      x = xmin + (xmax-xmin)*inow/double(n-1);

      // Integrate the ODE from xold to x
      rk_driver(xold, x);

      // Verbose
      if(verbose){
        std::cout << "Tried a total of " << nok+nbad << " steps. Accepted steps: " << nok << " Rejected steps: " << nbad << std::endl;
        std::cout << "i = " << std::setw(8) << inow << " / " << n << " x: " << std::setw(12) << x << " y: ";
        for(int i = 0; i < neq; i++) std::cout << std::setw(12) << y[i] << " ";
        std::cout << std::endl;
      }

      // Store solution
      xstore[inow] = x;
      for(int i = 0; i < neq; i++)
        ystore[i][inow] = y[i];

      // Update old values
      xold = x;
      return (++inow < n);
    }

    //=======================================================
    // Handling of error messages
    //=======================================================
    void report_error(std::string errormessage, int errorcode){
      std::cout << errormessage << std::endl;
      if(errorcode != 0) exit(0);
    }

    //=======================================================
    // List of solver methods
    //=======================================================
    void rk_step(std::vector<double> &y, double x, double h, std::vector<double> &yout, std::vector<double> &yerr);
    double take_and_adjust_step(std::vector<double> &y, double *x, double htry);
    void rk_driver(double x1, double x2);

};

//=======================================================
// Adjust stepsize until we have an acceptable step
//=======================================================

double OdeSolver::take_and_adjust_step(std::vector<double> &y, double *x, double htry){
  double errmax, h, htemp, xnew;
  std::vector<double> yerr (neq,0.0);
  std::vector<double> ytemp(neq,0.0);
  std::vector<double> yscal(neq,0.0);

  // Sum of the absolute value of the two sides of the ODE. Used as normalization when estimating the error
  for(int i = 0; i < neq; i++) yscal[i] = fabs(y[i]) + fabs(dydx[i]*htry) + TINY;

  // Check for NaN in solution
  if(y[1] != y[1]) report_error("NAN found in take_and_adjust_step",1);

  // Initial stepsize
  h = htry;

  // Find acceptable stepsize
  while(1) {

    // Try taking a step with stepsize h
    rk_step(y,*x,h,ytemp,yerr);

    // Estimate error
    errmax = 0.0;
    for (int i = 0; i < neq; i++){
      if (errmax < fabs(yerr[i]/yscal[i]))
        errmax = fabs(yerr[i]/yscal[i]);
    }
    errmax /= eps;

    // If error is small enough accept step
    if (errmax <= 1.0) break;

    // ...otherwise try a smaller stepsize
    htemp = SAFETY*h*pow(errmax,PSHRNK);
    if(h>=0.0){
      h = (htemp>0.1*h) ? htemp : 0.1*h;
    } else {
      h = (htemp<0.1*h) ? htemp : 0.1*h;
    }

    // Check for stepsize underflow
    xnew = (*x) + h;
    if (xnew == *x) report_error("Stepsize underflow in take_and_adjust_step",1);
  }

  // Increment x
  *x += (hdid = h);

  // Copy over solution
  y = ytemp;

  // Return next stepsize
  return (errmax > ERRCON) ? SAFETY*h*pow(errmax,PGROW) : 5.0*h;
}

//=======================================================
// Propagate y forward one step using a RK5 method
//=======================================================

void OdeSolver::rk_step(std::vector<double> &y, double x, double h, std::vector<double> &yout, std::vector<double> &yerr){
  std::vector<double> ak2, ak3, ak4, ak5, ak6, ytemp;
  ak2   = std::vector<double>(neq,0.0);
  ak3   = std::vector<double>(neq,0.0);
  ak4   = std::vector<double>(neq,0.0);
  ak5   = std::vector<double>(neq,0.0);
  ak6   = std::vector<double>(neq,0.0);
  ytemp = std::vector<double>(neq,0.0);

  // Runge-Kutta step 1
  for (int i = 0; i < neq; i++)
    ytemp[i] = y[i]+b21*h*dydx[i];

  // Runge-Kutta step 2
  dydx_func(x + a2*h, ytemp, ak2);
  for (int i = 0; i < neq; i++)
    ytemp[i] = y[i]+h*(b31*dydx[i]+b32*ak2[i]);

  // Runge-Kutta step 3
  dydx_func(x + a3*h, ytemp, ak3);
  for (int i = 0; i < neq; i++)
    ytemp[i] = y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);

  // Runge-Kutta step 4
  dydx_func(x + a4*h, ytemp, ak4);
  for (int i = 0; i < neq; i++)
    ytemp[i] = y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

  // Runge-Kutta step 5
  dydx_func(x + a5*h, ytemp, ak5);
  for (int i = 0; i < neq; i++)
    ytemp[i] = y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);

  // Runge-Kutta step 6
  dydx_func(x + a6*h, ytemp, ak6);
  for (int i = 0; i < neq; i++)
    yout[i] = y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);

  // Estimate error
  for (int i = 0; i < neq; i++)
    yerr[i] = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

//=======================================================
//  Driver for solving the system y'=f(x,y) on [x1,x2]
//=======================================================

void OdeSolver::rk_driver(double x1, double x2){
  double x, hnext, h;

  // Initialize x, stepsize and num. acc./rej. steps
  x   = x1;
  h   = (h1*(x2-x1) > 0.0) ? h1 : -h1;
  nok = nbad = 0;

  // Take steps
  for(int nstp = 0; nstp < MAXSTP; nstp++) {
    hdid = 0.0;

    // Calculate the right hand side dydx
    dydx_func(x,y,dydx); 

    // Initial guess for the stepsize
    if((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;

    // Make a step and get new stepsize
    hnext = take_and_adjust_step(y,&x,h);

    // Check if step was successful
    (hdid == h) ? ++nok : ++nbad;

    // Check if we are done
    if ((x-x2)*(x2-x1) >= 0.0) return;

    // Check if we have reached precision limit in stepsize
    if (fabs(hnext) <= hmin)
      report_error("Step size too small in odeint",1);
    h = hnext;
  } 

  // If we are here then steps > maxsteps
  report_error("Too many steps in routine odeint",1);
}

#endif
