#ifndef SOLVEPOISSON_INC
#define SOLVEPOISSON_INC
#include <stdlib.h>
#include <iostream>
#include <fftw3.h>
#include <math.h>

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

//===========================================================
// Class for solving Poisson-like equations 
// on the form D^2 Phi - m^2 Phi = S'
//
// The mass is assumes to be 0.0 be default
// Define this in set_massterm(realT)
//===========================================================

class PoissonSolver{
  private:
    fftw_complex *in, *out;
    fftw_plan p_forward, p_backward;
    bool weallocatein, weallocateout, madeplans;

    // 1D grid-size
    int n;

    // For testing
    realT nvec[3], norm;

    // Mass term in code units
    realT mass;
  public:

    PoissonSolver(){
      in=out=NULL;
      weallocatein=weallocateout=madeplans=false;
      mass=0.0;
      n=0;
    }

    //===========================================================
    // Only gridsize given so allocate in and out ourself
    //===========================================================
    PoissonSolver(int n){
      this->n = n;

      in  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n*n*n);
      out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n*n*n);

      p_forward  = fftw_plan_dft_3d(n, n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      p_backward = fftw_plan_dft_3d(n, n, n, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);

      weallocatein=weallocateout=true;
      madeplans=true;
      mass = 0.0;
    }

    //===========================================================
    // In and out arrays given so don't allocate anything here
    // After in will have the solution and out will have the
    // FFTed source multiplied by the laplacian term 1/k^2+m^2
    //===========================================================
    PoissonSolver(fftw_complex *in, fftw_complex *out, int n){
      this->n = n;

      this->in  = in;
      this->out = out;

      p_forward  = fftw_plan_dft_3d(n, n, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
      p_backward = fftw_plan_dft_3d(n, n, n, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);

      weallocatein=weallocateout=false;
      madeplans=true;
      mass = 0.0;
    }

    //===========================================================
    // In array given so just use this. In will be overwritten with the solution
    //===========================================================
    PoissonSolver(fftw_complex *in, int n){
      this->n = n;

      this->in  = in;
      this->out = in;

      p_forward  = fftw_plan_dft_3d(n, n, n, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
      p_backward = fftw_plan_dft_3d(n, n, n, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);

      weallocatein=weallocateout=false;
      madeplans=true;
      mass = 0.0;
    }

    //===========================================================
    // Real in array is given. Allocate internal in array and copy over data
    // In the end the solution will be in the internal in[] only
    //===========================================================
    PoissonSolver(realT *in, int n){
      this->n = n;

      this->in  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n*n*n);
      this->out = this->in;

      for(int i=0;i<n*n*n;i++){
        this->in[i][0] = in[i];
        this->in[i][1] = in[0];
      }

      p_forward  = fftw_plan_dft_3d(n, n, n, this->in, this->in, FFTW_FORWARD, FFTW_ESTIMATE);
      p_backward = fftw_plan_dft_3d(n, n, n, this->in, this->in, FFTW_BACKWARD, FFTW_ESTIMATE);

      weallocatein=true;
      weallocateout=false;
      madeplans=true;
      mass = 0.0;
    }

    ~PoissonSolver(){
      clearall();
    }

    //===========================================================
    // This is the term 'm' in 'D^2 Phi - m^2 Phi = S'
    // We have m_code = (m_real * Boxsize)/(2pi) in units of c=hbar=1
    //===========================================================
    void set_massterm(realT mass){
      this->mass = mass;
    }

    //===========================================================
    // Does everything: FFT, *= -1/(k^2+m^2) and the FFT back with
    // correct normalization
    //===========================================================
    void execute(){
      fftw_execute(p_forward);
      divide_by_laplacian();
      fftw_execute(p_backward);
    }

    //===========================================================
    // Delete everything
    //===========================================================
    void clearall(){
      if(weallocatein)
        fftw_free(in);
      if(weallocateout)
        fftw_free(out);
      if(madeplans){
        fftw_destroy_plan(p_forward);
        fftw_destroy_plan(p_backward);
      }
      n=0;
      mass=0.0;
      weallocatein=false;
      weallocateout=false;
      madeplans=false;
      in=out=NULL;
    }

    //===========================================================
    // This divides the FFTed source by the term 
    // k^2 + m^2 and takes care of units
    //===========================================================
    void divide_by_laplacian(){
      realT norm_poisson, norm_fftw, fac;
      int n2 = n*n, nover2 = n/2;

      norm_poisson  = -1.0/(4.0*acos(-1)*acos(-1));
      norm_fftw     = 1.0/realT(n * n * n);
      fac           = norm_poisson * norm_fftw;

      for(int i=0;i<n;i++){
        int ii = (i < nover2 ? i: i-n);
        for(int j=0;j<n;j++){
          int jj = (j < nover2 ? j: j-n);
          for(int k=0;k<n;k++){
            int kk = (k < nover2 ? k: k-n);
            int ind = i + j*n + k*n2;
            if(ind > 0){
              realT facnow = fac/(ii*ii+jj*jj+kk*kk + mass*mass);
              out[ind][0] *=  facnow;
              out[ind][1] *=  facnow;
            }
          }
        }
      }
      out[0][0] = 0.0;
      out[0][1] = 0.0;
    }

    //===========================================================
    // The function source Source = A sin(2pi nx x) sin(2pi ny y) sin(2pi nz z)   
    // has the solution Phi = B Source for a constant B. 
    // We use this as a test-case where nx,ny,nz are randomly generate
    //===========================================================
    realT test_func(int i, bool before){
      static realT twopi = 2.0*acos(-1);
      realT x,y,z, func;

      x = ((i % n) + 0.5)/realT(n);
      y = ((i/n % n) + 0.5)/realT(n);
      z = ((i/(n*n) % n) + 0.5)/realT(n);

      func = sin(twopi * nvec[0] * x) * sin(twopi * nvec[1] * y) * sin(twopi * nvec[2] * z);

      // Make sure we give the correct norm factor (before is 'A/B', after is 1)
      if(before)
        return func*norm;
      else
        return func;
    }

    //===========================================================
    // Run tests with the solver
    //===========================================================
    void run_test(){
      realT twopi = 2.0*acos(-1);

      std::cout << "=============================" << std::endl;
      std::cout << "    Running standard test      " << std::endl;
      std::cout << "     for Poisson solver        " << std::endl;
      std::cout << "=============================" << std::endl;

      // Generate random integer vector
      nvec[0] = 1.*(rand() % 10 - 5);
      nvec[1] = 1.*(rand() % 10 - 5);
      nvec[2] = 1.*(rand() % 10 - 5);

      std::cout << "Generate random integer vector (n1,n2,n3) = ( ";
      std::cout << nvec[0] << " , " << nvec[1] << " , " << nvec[2] << " ) " << std::endl;
      std::cout << "Set source S corresponding to the solution:" << std::endl;
      std::cout << "Phi(x,y,z) = sin(twopi * n1 * x) * sin(twopi * n2 * y) * sin(twopi * n3 * z)" << std::endl;

      // Norm for test_func
      norm = -(twopi*twopi)*(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2] + mass*mass);

      init_to_test_function();
      print_test(20,true);
      execute();
      print_test(20,false);
    }

    //===========================================================
    // Initialize source to test_func(x)
    //===========================================================
    void init_to_test_function(){
      for(int i=0;i<n*n*n;i++){
        in[i][0] = test_func(i,true);
        in[i][1] = 0.0;
      }
    }

    //===========================================================
    // Outputs results from test-runs
    //===========================================================
    void print_test(int r, bool before){
      std::cout << std::endl;
      if(before){
        std::cout << "Printing the source for every " << r << " elements of in: " << std::endl;
      }else{
        std::cout << "Printing the solution for every " << r << " elements of in: " << std::endl;
      }
      for(int i=0;i<n*n*n;i++){
        if(i % ((n*n*n)/r +1) == 0)
          if(before)
            printf("S = % 10.5e  S_imag = % 10.5e\n",in[i][0],in[i][1]);
          else
            printf("y = % 10.5e  (ytrue = % 10.5e)  yimag = % 10.5e  (yimagtrue = % 10.5e)\n",out[i][0],test_func(i,before),out[i][1],0.0);
      }
      std::cout << std::endl;
    }

    //===========================================================
    // Copy the results from the internal arrays
    //===========================================================
    void copy_over_array(realT *res, bool copy_in_array, bool realpart){
      int ii = 0;
      if(!realpart) ii = 1;
      if(copy_in_array){
        for(int i=0;i<n*n*n;i++)
          res[i] = in[i][ii];
      } else {
        for(int i=0;i<n*n*n;i++)
          res[i] = out[i][ii];
      }
    }
};
#endif
