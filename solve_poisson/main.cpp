#include<limits>
#include<iomanip>
#include "solve_poisson.h"

///////////////////////////////////////
// Run built-in test of the
// solver to check that it works
//
// Set source so that the solution will
// be Phi(x,y,z) = Sin(2pi*n1*x) * 
//                 Sin(2pi*n1*x) * 
//                 Sin(2pi*n1*z)
///////////////////////////////////////

void run_standard_test(){
  int ngrid = 128;

  // Set up solver
  PoissonSolver mysolver(ngrid);

  // Run standard test
  mysolver.run_test();
}

///////////////////////////////////////
// Make a random source and solve
// the equation.
//
// For this case we keep the source 
// and solution in different arrays
// Other alternatives availiable:
// see solve_poisson.h
//////////////////////////////////////

void solve_poisson(){
  int ngrid = 32;
  fftw_complex *source, *solution;
  realT avg = 0.0;

  std::cout << "=============================" << std::endl;
  std::cout << "   Solve Poisson equation    " << std::endl;
  std::cout << "     with random source      " << std::endl;
  std::cout << "=============================" << std::endl;

  // Allocate memory
  source   = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ngrid*ngrid*ngrid);
  solution = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ngrid*ngrid*ngrid);

  // Fill source with random data
  for(int i = 0;i<ngrid*ngrid*ngrid;i++){
    source[i][0] = (rand() % std::numeric_limits<int>::max()) / realT(std::numeric_limits<int>::max()) - 0.5;
    source[i][1] = 0.0;
    solution[i][0] = solution[i][1] = 0.0;
    avg += source[i][0];
  }

  // For this test to be consistent enforce <S> = 0 as is required by the equation
  avg /= realT(ngrid*ngrid*ngrid);
  for(int i=0;i<ngrid*ngrid*ngrid;i++)
    source[i][0] -= avg;

  // Set up solver
  PoissonSolver mysolver(source, solution, ngrid);

  // Solve the equation
  mysolver.execute();

  // Print the solution
  for(int i=0;i<ngrid*ngrid*ngrid;i++){

    // Only print every 1000 of the elements
    if(rand() % 1000 <= 1){
      // This illustrates how x,y,z in the grid are related to i
      realT x = (i % ngrid) / realT(ngrid);
      realT y = i/ngrid % ngrid / realT(ngrid);
      realT z = i/(ngrid*ngrid) / realT(ngrid);

      std::cout << " x = " << std::setw(12) << x << "  Phi = " << std::setw(12) << solution[i][0];
      std::cout << " Phi_imag = " << std::setw(12) << solution[i][1] << std::endl;
    }
  }
}

int main(int argc, char **argv){
  run_standard_test();
  solve_poisson();
}

