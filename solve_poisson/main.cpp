#include<limits>
#include<iomanip>
#include "PoissonSolver.h"
#define NDIM 3

//======================================================
//
// Code to solve 1D, 2D and 3D Poisson equations
//
//              D^2 Phi = m^2 Phi + S
//
// where D^2 is the NDIM-dimensional Laplacian in a box 
// with periodic boundary conditions.
// 
// The mass 'm' is zero by default. Change it by calling
// set_mass(m).
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//
//======================================================

void run_standard_test(int ngrid){

  //======================================================
  // Run built-in test of the solver to check that it works
  // Set source so that the solution will
  // be Phi(x,y,z) = Prod_i Sin(2pi*n1*x_i)
  //======================================================

  // Run standard test for ndim = 1
  PoissonSolver<1> psolver_1(ngrid);
  psolver_1.run_test();

  // Run standard test for ndim = 2
  PoissonSolver<2> psolver_2(ngrid);
  psolver_2.run_test();
  
  // Run standard test for ndim = 3
  PoissonSolver<3> psolver_3(ngrid);
  psolver_3.run_test();
}

void solve_poisson(int ngrid){
  int ntot = n_to_ndim(ngrid, NDIM);
  fftw_complex *source, *solution;
  double avg = 0.0, x[NDIM];

  //======================================================
  // Make a random source and solve the equation
  // For this case we keep the source and solution in 
  // different arrays. Other alternatives availiable.
  //======================================================

  std::cout << "=======================================================" << std::endl;
  std::cout << "       Solve Poisson equation with random source       " << std::endl;
  std::cout << "=======================================================" << std::endl;

  // Allocate memory
  source   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ntot);
  solution = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ntot);

  // Fill source with random data
  for(int i = 0; i < ntot; i++){
    source[i][0] = (rand() % std::numeric_limits<int>::max()) / double(std::numeric_limits<int>::max()) - 0.5;
    source[i][1] = 0.0;
    solution[i][0] = solution[i][1] = 0.0;
    avg += source[i][0];
  }

  // For this test to be consistent enforce <S> = 0 as is required by the equation
  avg /= double(ntot);
  for(int i = 0; i < ntot; i++)
    source[i][0] -= avg;

  // Set up solver
  PoissonSolver<NDIM> mysolver(source, solution, ngrid);

  // Solve the equation
  mysolver.execute();

  // Print the solution
  for(int i = 0; i < ntot; i++){

    // Print 20 of the elements
    if(i % (ntot/20 + 1) == 0){

      // How x,y,z in the grid are related to i
      for(int j = 0, nn = 1; j < NDIM; j++, nn *= ngrid)
        x[j] = (i/nn % ngrid) / double(ngrid);

      std::cout << " x = " << std::setw(12) << x[0] << "  Re[Phi] = " << std::setw(12) << solution[i][0];
      std::cout << " Im[Phi] = " << std::setw(12) << solution[i][1] << std::endl;
    }
  }
}

int main(int argc, char **argv){
  srand(time(NULL));
  run_standard_test(32);
  solve_poisson(32);
}

