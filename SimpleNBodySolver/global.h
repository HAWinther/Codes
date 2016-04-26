#ifndef _GLOBALINC
#define _GLOBALINC
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <assert.h>
#ifdef OMP
#include <omp.h>
#endif
#include "Spline.h"
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

using namespace std;

class MultiGrid {

  private:

  //=========================================================
  // A container for all the grids needed for the multigrid
  // solver
  //=========================================================

  public:

    int N, Ntot;
    int Nlevel;
    int Nmin;
    bool verbose;
    double rms_res, rms_res_i, rms_res_old;
    double **f, **res, **rho, **source;

    MultiGrid(double *ff, double *rres, double *rrho, double *ssource, int ngrid, bool vverbose){
      N      = ngrid;
      Ntot   = pow3(ngrid);
      Nmin   = 3;
      Nlevel = int(log2(N)-Nmin+1);
      verbose = vverbose;
      assert(Nlevel >= 2);
      assert(ff != NULL);
      assert(rrho != NULL);
      assert(ssource != NULL);
      assert(rres != NULL);

      // Allocate all levels
      f      = new double*[Nlevel];
      res    = new double*[Nlevel];
      rho    = new double*[Nlevel];
      source = new double*[Nlevel];

      // Pointers to the domain-grid arrays
      f[0]      = ff;
      res[0]    = rres;
      rho[0]    = rrho;
      source[0] = ssource;

      // Allocate level for level
      for (int i = 1, Nnow = N/2; i < Nlevel; i++, Nnow /= 2) {
        int NnowTot = pow3(Nnow);
        f[i]      = new double[NnowTot];
        rho[i]    = new double[NnowTot];
        res[i]    = new double[NnowTot];
        source[i] = new double[NnowTot];
      }
    }

    ~MultiGrid(){
      for(int i = 1;i<Nlevel;i++){
        delete[] f[i];
        delete[] rho[i];
        delete[] res[i];
        delete[] source[i];
      }
      f[0] = rho[0] = res[0] = source[0] = NULL;
      delete[] f;
      delete[] rho;
      delete[] res;
      delete[] source;
      f = rho = res = source = NULL;
    }
};

//=========================================================
// A particle container
//=========================================================

class Particle{
  public:
    double _x[3], _v[3];

    Particle(){}
    Particle(double *x, double *v){
      for(int i = 0; i < 3; i++){
        _x[i] = x[i];
        _v[i] = v[i];
      }
    }
    ~Particle(){}
};

//=========================================================
// Global simulations struct containing all the arrays
// and variables for the simulation
//=========================================================

struct Simulation{
  // Grid related variables
  int  _ngrid, _ngrid_tot;

  // Cosmology
  double _aexp_ini, _aexp;
  double _BoxSize, _omega_m, _omega_lambda;

  // Time-step control
  double _courant_fac;

  // Newtonian Gravity variables
  double *_Rho;
  fftw_complex *_Phi;       
  fftw_complex *_Force[3]; 

  // FFT plans
  fftw_plan _plan_phi_forward, _plan_phi_backward;
  fftw_plan _plan_force_backward[3];

  // Modified gravity variables
  bool _modified_gravity;
  double *_ScalarField;
  double *_FifthForce[3];
  double *_MGSource, *_MGResidual;
  MultiGrid *_mg;

  // f(R) gravity
  double _nfofr, _fofr0;

  // Symmetron
  double _beta_symm, _lambda_symm, _assb_symm;

  // Convergence criterion and number of sweeps over the grid per step in V-cycle
  double _EPS_CONVERGE, _EPS_RELDIR;
  int _NGS_FINE, _NGS_COARSE;
  bool _verboseMG;

  // Spline for tt(a) and a(tt)
  DSpline _ttofa_spline;
  DSpline _aoftt_spline;

  // Particle data
  int  _npart;
  string _particlefile;
  Particle *_Particle;

  // Other stuff
  bool _static;  
  bool _analytic_matter_density;

} sim;

// Multigrid methods
void make_initial_guess(double *Phi, double aexp);
double l_operator(int level, int i, double kinetic, bool makesource);
double dl_operator(int level, int i, double dkinetic);
void solve_phi_with_multigrid();
double calc_gs_change(int level, int i, double kinetic, double dkinetic);
double residual(int level, int i, double kinetic);
double calculate_residual(int level);
bool check_convergence(int istep);
void prolonge_up(int to_level, double* Old, double* New);
void Gauss_Seidel(int level, int redblack);
void solve_current_level(int level);
void correct_solution(int level);
void recursive_go_up(int to_level);
void make_prolongation_array(int level);
void make_new_source(int level);
void recursive_go_down(int from_level);
void restrict_down_density_all();
void restrict_down_single(int from_level, double *Top, double *Bottom);
void calc_min_max_scalarfield();
void calculate_fifth_force();
void calculate_gradient(double *phi, double *gradient, int dim, int ngrid);

// Methods in main.cpp
double Hubble(double aexp);
double calc_vrms(double aexp);
void init(int argc, char **argv);
void allocate_memory();
void deallocate_memory();
void read_particles_from_file(string particlefile);
void main_timestep_loop();
void assign_particles_to_grid(double *rho, int ngrid);
void leapfrog_drift_particles(double dt);
void leapfrog_kick_particles(double dt);
void enforce_periodic_bc_particles();
void compute_gravitational_acceleration(double aexp);
void output_particles(double aexp);
void compute_time_scalefactor_spline();
void output_grid(fftw_complex *grid, int ngrid, string filename);
void output_grid(double *grid, int ngrid, string filename);
void set_analytical_density_field(double *rho, int ngrid);
void static_test(double aexp);
pair<double,double> calc_energy();
double calc_energy_conservation_constant();
double calc_new_timestep();
inline void calc_forces_grid(Particle *p, double *f);
inline double cic_kernel(double *xg, double *xp);
inline double tsc_kernel(double *xg, double *xp);

#endif
