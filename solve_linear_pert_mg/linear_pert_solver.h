#ifndef LINEARPERTSOLVER_HEADER
#define LINEARPERTSOLVER_HEADER
#include "ode_solver.h"
#include "spline.h"

// Spline up P(k,a)/P(k,a)_LCDM
#define FULLEVOSPLINES

//////////////////////////////////////////
// Global book-keeping variable
//////////////////////////////////////////

struct params{
  // Cosmology
  double Omegal, Omegam;

  // Modified gravity
  double fofr0, n_fofr;

  // k/H0
  double know;
} global;

//////////////////////////////////////////
// P(k) Solver Class
//////////////////////////////////////////

class PofkModifiedGravity{
  private:
    double kmin, kmax;
    int nk;

    double aini, aend;
    int nint;

    double *kk, *pofk;

    void (*derivs_lcdm)(realT, realT *, realT *);
    void (*derivs_mg)(realT, realT *, realT *);

    Spline *pofkratio_spline;

#ifdef FULLEVOSPLINES
    // Splines for the full P(k,a)/P_LCDM(k,a)
    Spline *delta_lcdm_x_spline, **delta_mg_x_spline;
#endif

    bool solved;
  public:

    PofkModifiedGravity(){
      solved = false;
    }

    PofkModifiedGravity(double kmin, double kmax, int nk, double aini, double aend, int nint, 
        void (*derivs_lcdm)(realT, realT *, realT *),  void (*derivs_mg)(realT, realT *, realT *)){
      this->kmin = kmin;
      this->kmax = kmax;
      this->nk = nk;

      this->aini = aini;
      this->aend = aend;
      this->nint = nint;

      this->derivs_mg = derivs_mg;
      this->derivs_lcdm = derivs_lcdm;

      kk = new double[nk];
      pofk = new double[nk];
      solved = false;
    }

    ~PofkModifiedGravity(){

      if(kk !=NULL)
        delete[] kk;
      if(pofk !=NULL)
        delete[] pofk;
      kk = pofk = NULL;
      if(pofkratio_spline != NULL)
        delete pofkratio_spline;
      pofkratio_spline = NULL;
      solved = false;
#ifdef FULLEVOSPLINES
      if(delta_lcdm_x_spline != NULL)
        delete delta_lcdm_x_spline;
      if(delta_mg_x_spline != NULL){
        for(int i=0;i<nk;i++)
          delete delta_mg_x_spline[i];
        delete[] delta_mg_x_spline;
      }
#endif
    }
#ifdef FULLEVOSPLINES
    // k in units of h/Mpc
    double pofk_ratio_of_a_k(double a, double k){
      double x = log(a);
      int indlow, indhigh;

      // Get LCDM value
      double delta_lcdm = delta_lcdm_x_spline->get_spline(x);

      // Get index in k-array
      indlow = int(log(k/kmin)/log(kmax/kmin)*(nk-1) + 0.5);

      // Bounds check
      if(indlow > nk-1 || indlow < 0){
        std::cout << "Error: out of bounds in pofk_ratio_of_a_k" << std::endl;
        std::cout << "kmin < k < kmax : " << kmin << " " << k << " " << kmax << std::endl;
        exit(1);
      }
      indhigh = indlow+1;

      double ratio;
      if(indhigh > nk-1){
        double delta_mg_low  = delta_mg_x_spline[indlow]->get_spline(x);
        ratio = pow(delta_mg_low/delta_lcdm,2);	
      } else {
        double delta_mg_low  = delta_mg_x_spline[indlow]->get_spline(x);
        double delta_mg_high = delta_mg_x_spline[indhigh]->get_spline(x);

        // Linear interpolation
        ratio = pow(delta_mg_low/delta_lcdm,2) + (k - kk[indlow])/(kk[indhigh]-kk[indlow]) * (pow(delta_mg_high/delta_lcdm,2) - pow(delta_mg_low/delta_lcdm,2)); 
      }

      return ratio;
    }
#endif

    // k in units of h/Mpc
    double pofk_ratio(double k){
      if(!solved){
        std::cout << "Error in PofkModifiedGravity. Cannot get ratio before solving" << std::endl;
        exit(1);
      }
      return pofkratio_spline->get_spline(k);
    }

    void solve(bool verbose){
      double *ic, *y, *x, delta_lcdm_0;
      double xmin, xmax;

      // Set the initial conditions
      xmin = log(aini);
      xmax = log(aend);

      // At early times delta = a = exp(x) ic[0] = ic[1]
      ic = new double[2];
      ic[0] = 1.0;
      ic[1] = 1.0;

      // Set up solver
      OdeSolver lcdm(nint, 2, derivs_lcdm);

      // Set initial conditions
      lcdm.set_initial_conditions(xmin, xmax, ic);

      // Solve
      lcdm.solve(false);

      // The values of delta today
      y = lcdm.y_array(0);
      x = lcdm.x_array();
      delta_lcdm_0 = y[nint-1];

#ifdef FULLEVOSPLINES
      // Spline delta_lcdm(x)
      delta_lcdm_x_spline = new Spline();
      delta_lcdm_x_spline->make_spline(x,y,nint,1e30,1e30,0,"delta_LCDM");

      // Allocate memory for spline
      delta_mg_x_spline = new Spline*[nk];
#endif

      // Loop over k'values and caluclate delta_m(a,k)
      for(int i=0;i<nk;i++){
        // Set global k-values
        global.know = exp(log(kmin) + log(kmax/kmin)*i/double(nk-1));

        // Set up solver
        OdeSolver mg(nint, 2, derivs_mg);

        // Set initial conditions
        mg.set_initial_conditions(xmin, xmax, ic);

        // Solve
        mg.solve(false);

        // The values of delta(x)
        y = mg.y_array(0);
        x = lcdm.x_array();

#ifdef FULLEVOSPLINES
        // Spline it up
        delta_mg_x_spline[i] = new Spline();
        delta_mg_x_spline[i]->make_spline(x,y,nint,1e30,1e30,0,"delta_mg");
#endif

        // The value delta_m(k,a=1) / delta_m_LCDM(k)
        if(verbose && i % 20 == 0) std::cout << "k: " << global.know << " Delta/Delta_LCDM: " << y[nint-1]/delta_lcdm_0 << std::endl;

        // Store values
        kk[i] = global.know;
        pofk[i] = pow(y[nint-1]/delta_lcdm_0,2);
      }

      // Spline
      pofkratio_spline = new Spline();
      pofkratio_spline->make_spline(kk,pofk,nk,1e30,1e30,2,"P(k) / P(k)_LCDM ratio");
      delete[] ic;
      solved = true;
    }
};

#endif
