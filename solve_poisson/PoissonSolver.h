#ifndef _POISSONSOLVERHEADER_
#define _POISSONSOLVERHEADER_
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>

int n_to_ndim(int n, int ndim);

template<int NDIM>
class PoissonSolver{

  //===========================================================
  // 
  // Class for solving 1D,2D,3D Poisson equations on the form 
  //               D^2 Phi - m^2 Phi = S
  //
  // The mass is assumes to be 0.0 be default
  // Define this in set_mass(mass)
  //
  // If the box is not-square then after rescaling coordinates the
  // most general equation is
  //
  //    [A d/dx^2 + B d/dy^2 + C d/dz^2 - m^2] Phi = S
  // 
  // The factors A,B,C is 1 by standard but can be changed by
  // calling set_scalefactor(A,B,C)
  //
  //===========================================================

  private:

    unsigned int n;          // 1D gridsize
    unsigned int ntot;       // Total number of gridnodes

    fftw_complex *in;        // The fftw arrays
    fftw_complex *out;       

    bool allocate_in;        // Do we allocate the in array or just use pointer?
    bool allocate_out;       // Do we allocate the out array or just use pointer?

    double mass;             // Mass term in code units

    double scalefactor[NDIM]; // For unequal axes. Standard = {1,1,1,...}
    
    double nvec[NDIM], norm; // For testing
  public:

    PoissonSolver() : 
          n(0),
          ntot(0),
          in(NULL), 
          out(NULL),
          allocate_in(false), 
          allocate_out(false), 
          mass(0.0) {
            set_scalefactor();
          }

    //===========================================================
    // Only gridsize given so allocate in and out ourself
    //===========================================================
    PoissonSolver(unsigned int _n) : 
          n(_n), 
          ntot(n_to_ndim(n,NDIM)),
          in((fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ntot)),
          out((fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ntot)),
          allocate_in(true), 
          allocate_out(true), 
          mass(0.0) {
            set_scalefactor();
          }

    //===========================================================
    // In and out arrays given so don't allocate anything here
    // After in will have the solution and out will have the
    // FFTed source multiplied by the laplacian term 1/k^2+m^2
    //===========================================================
    PoissonSolver(fftw_complex *_in, fftw_complex *_out, unsigned int _n) : 
          n(_n), 
          ntot(n_to_ndim(n,NDIM)),
          in(_in),
          out(_out),
          allocate_in(false),
          allocate_out(false),
          mass(0.0) {
            set_scalefactor();
          }

    //===========================================================
    // In array given so just use this. In will be overwritten with the solution
    //===========================================================
    PoissonSolver(fftw_complex *_in, unsigned int _n) : 
          n(_n), 
          ntot(n_to_ndim(n,NDIM)),
          in(_in), 
          out(_in),
          allocate_in(false),
          allocate_out(false),
          mass(0.0) {
            set_scalefactor();
          }

    //===========================================================
    // Real in array is given. Allocate internal in array and copy over data
    // In the end the solution will be in the internal in[] only
    //===========================================================
    PoissonSolver(double *_in, unsigned int _n) : 
          n(_n),
          ntot(n_to_ndim(n,NDIM)),
          in((fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ntot)),
          out(in),
          allocate_in(true),
          allocate_out(false),
          mass(0.0) {
            for(int i = 0; i < ntot; i++){
              in[i][0] = _in[i];
              in[i][1] = 0.0;
            }
            set_scalefactor();
          }
    
    //===========================================================
    // Copy constructor
    //===========================================================
    PoissonSolver(PoissonSolver& p){
      n            = p.n;
      ntot         = p.ntot;
      if(p.allocate_in){
        in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * p.ntot);
        std::memcpy(in, p.in, p.ntot * sizeof(fftw_complex));
      } else {
        in = p.in;
      }
      if(p.allocate_out){
        out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * p.ntot);
        std::memcpy(out, p.out, p.ntot * sizeof(fftw_complex));
      } else {
        out = p.out;
      }
      allocate_in  = p.allocate_in;
      allocate_out = p.allocate_out;
      mass         = p.mass;
      std::memcpy(scalefactor, p.scalefactor, NDIM * sizeof(double));
    }
    
    //===========================================================
    // Assignment operator
    //===========================================================
    PoissonSolver operator=(PoissonSolver& rhs){
      PoissonSolver tmp( rhs );
      std::swap(n, tmp.n);
      std::swap(ntot, tmp.ntot);
      std::swap(in, tmp.in);
      std::swap(out, tmp.out);
      std::swap(allocate_in, tmp.allocate_in);
      std::swap(allocate_out, tmp.allocate_out);
      std::swap(mass, tmp.mass);
      std::swap(scalefactor, tmp.scalefactor);
      return *this;
    }

    ~PoissonSolver() { clear(); }
    
    //===========================================================
    // Extract data
    //===========================================================
    fftw_complex* get_in(){ return in; }
    fftw_complex* get_out(){ return out; }
    double get_mass(){ return mass; }
    int get_n(){ return n; }
    int get_ndim(){ return NDIM; }

    //===========================================================
    // This is the term 'm' in 'D^2 Phi - m^2 Phi = S'
    // We have m_code = (m_real * Boxsize)/(2pi) in units of c=hbar=1
    //===========================================================
    void set_mass(double _mass) { mass = _mass; }
  
    //===========================================================
    // In the case of a non-square box the Laplace operator can
    // be written [A_x d/dx^2 + A_y d/dy^2 + A_z d/dz^2 - m^2]Phi = S
    // Set the factors A_x, A_y, A_z here
    //===========================================================
    double set_scalefactor(double Ax = 1.0, double Ay = 1.0, double Az = 1.0){
      scalefactor[0] = Ax;
      if(NDIM >= 2)
        scalefactor[1] = Ay;
      if(NDIM == 3)
        scalefactor[2] = Az;
    }

    //===========================================================
    // Does everything: FFT, *= -1/(k^2+m^2) and the FFT back with
    // correct normalization
    //===========================================================
    void execute(){
      fftw_plan p_forward, p_backward;

      // Make plans
      if(NDIM == 3){
        p_forward  = fftw_plan_dft_3d(n, n, n, in,  out, FFTW_FORWARD,  FFTW_ESTIMATE);
        p_backward = fftw_plan_dft_3d(n, n, n, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
      } else if (NDIM == 2){
        p_forward  = fftw_plan_dft_2d(n, n,    in,  out, FFTW_FORWARD , FFTW_ESTIMATE);
        p_backward = fftw_plan_dft_2d(n, n,    out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
      } else if (NDIM == 1){
        p_forward  = fftw_plan_dft_1d(n,       in,  out, FFTW_FORWARD,  FFTW_ESTIMATE);
        p_backward = fftw_plan_dft_1d(n,       out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
      } else {
        std::cout << "Error [SolvePoisson]->[execute] ndim is not valid" << std::endl;
        exit(1);
      }

      // Transform
      fftw_execute(p_forward);
      divide_by_laplacian();
      fftw_execute(p_backward);

      // Clean up
      fftw_destroy_plan(p_forward);
      fftw_destroy_plan(p_backward);
    }

    //===========================================================
    // Delete everything
    //===========================================================
    void clear(){
      if(allocate_in){
        fftw_free(in);
        allocate_in = false;
      }
      if(allocate_out){
        fftw_free(out);
        allocate_out = false;
      }
      n = ntot = 0;
    }

    //===========================================================
    // This divides the FFTed source by the term 
    // k^2 + m^2 and takes care of units
    //===========================================================
    void divide_by_laplacian(){
      double norm_poisson, norm_fftw, fac;
      int n2 = n*n, nover2 = n/2, ix, ll;

      // Normalization factors
      norm_poisson  = -1.0/(4.0*M_PI*M_PI);
      norm_fftw     = 1.0/double(ntot);
      fac           = norm_poisson * norm_fftw;

      // ind = ix + n*iy + n^2*iz and facnow = m^2 + kx^2 + ky^2 + kz^2
      for(int i = 0; i < ntot; i++){
        int ind = 0;
        double facnow = mass * mass;
        for(int j = 0, nn = 1; j < NDIM; j++, nn *= n){
          ix = (i/nn % n);
          ll = ix < nover2 ? ix : ix - n;
          ind += ix * nn;
          facnow += ll * ll * scalefactor[j];
        }
        facnow = fac / facnow;
        if(ind > 0){
          out[ind][0] *= facnow;
          out[ind][1] *= facnow;
        }
      }

      // Set zero-mode to zero
      out[0][0] = 0.0;
      out[0][1] = 0.0;
    }

    //===========================================================
    // The function source Source = A Prod_i sin(2pi n_i x_i)
    // has the solution Phi = B Source for a constant B. 
    // We use this as a test-case where n_i is randomly generated
    //===========================================================
    double test_func(int i, bool before){
      double x[NDIM], func = 1.0;

      for(int j = 0, nn = 1; j < NDIM; j++, nn *= n)
        x[j] = ((i/nn % n) + 0.5)/double(n);
     
      for(int j = 0; j < NDIM; j++)
        func *= sin(2.0 * M_PI * nvec[j] * x[j]);

      // Make sure we give the correct norm-factor (before is 'A/B', after is 1)
      return (before ? func * norm : func);
    }

    //===========================================================
    // Run tests with the solver
    //===========================================================
    void run_test(){
      std::cout << "=======================================================" << std::endl;
      std::cout << "  Running standard test for Poisson solver NDIM = " << NDIM << std::endl;
      std::cout << "=======================================================" << std::endl;

      // Generate random integer vector
      for(int i = 0; i < NDIM; i++)
        nvec[i] = 1.0*(rand() % 10 + 1);

      std::cout << "Generate random integer vector n = ( ";
      std::cout << nvec[0];
      for(int i = 1; i < NDIM; i++)
        std::cout << " , " << nvec[i];
      std::cout << " )" << std::endl;
      std::cout << "Set source S corresponding to the solution:" << std::endl;
      std::cout << "Phi = Prod_i sin(2pi * ni * x_i)" << std::endl;

      // Norm for test_func
      norm = mass * mass;
      for(int i = 0; i < NDIM; i++)
        norm += nvec[i] * nvec[i] * scalefactor[i];
      norm *= -4.0 * M_PI * M_PI;

      // Run the test
      init_source_to_testfunction();
      print_test(20, true);
      execute();
      print_test(20, false);
    }

    //===========================================================
    // Initialize source to test_func(x)
    //===========================================================
    void init_source_to_testfunction(){
      for(int i = 0; i < ntot; i++){
        in[i][0] = test_func(i, true);
        in[i][1] = 0.0;
      }
    }

    //===========================================================
    // Outputs results from test-runs
    //===========================================================
    void print_test(int r, bool before){
      double avg_err = 0.0, max_err = 0.0, err;

      std::cout << std::endl;
      if(before){
        std::cout << "Printing the source for every " << r << " elements of in: " << std::endl;
      }else{
        std::cout << "Printing the solution for every " << r << " elements of in: " << std::endl;
      }
      for(int i = 0; i < ntot; i++){
        if(!before){
          err = test_func(i, before) - out[i][0];
          avg_err += fabs(err);
          if(err > max_err) max_err = err;
        }
        if(i % (ntot/r + 1) == 0)
          if(before)
            std::cout << "Re[S] = " << std::setw(12) << in[i][0] << "  Im[S] = " << std::setw(12) << in[i][1] << std::endl;
          else {
            std::cout << "Re[y] = " << std::setw(12) << out[i][0] << " (Re[y - y_exact] = " << std::setw(12) << test_func(i, before) - out[i][0];
            std::cout << " Im[y] = " << std::setw(12) << out[i][1] << std::endl;
          }
      }
      if(!before){
        avg_err /= double(ntot);
        std::cout << std::endl;
        std::cout << "Avg error: " << std::setw(12) << avg_err << "    Max error: " << std::setw(12) << max_err << std::endl;
        std::cout << std::endl;
      }
    }

    //===========================================================
    // Copy the results from the internal arrays
    //===========================================================
    void copy_over_array(double *res, bool copy_in_array, bool realpart){
      int ii = 0;
      if(!realpart) ii = 1;
      if(copy_in_array){
        for(int i = 0; i < ntot; i++)
          res[i] = in[i][ii];
      } else {
        for(int i = 0; i < ntot; i++)
          res[i] = out[i][ii];
      }
    }
};

// Compute n^NDIM
int n_to_ndim(int n, int ndim){
  if(ndim <= 0 || ndim > 3){
    std::cout << "Error [PoissonSolver]: not valid ndim = " << ndim << std::endl;
    exit(1);
  }
  int ntondim = 1;
  for(int i = 0; i < ndim; i++)
    ntondim *= n;
  return ntondim;
}

#endif
