#include <iomanip>
#include "ode_solver.h"

using namespace std;

///////////////////////////////////////
// The ODE dy/dx = -2xy
// ===> y[x] = y[0] * Exp[-x^2]
///////////////////////////////////////

void ode1(realT x, realT *y, realT *dydx){
	dydx[0] = -2.0*x*y[0];
}

///////////////////////////////////////
// The ODE {dy_1/dx = y_2, dy_2/dx = 1}
// ==> y_2[x] = y_2[0] + x
//     y_1[x] = y_1[0] + y_2[0]x + x^2/2
///////////////////////////////////////

void ode2(realT x, realT *y, realT *dydx){
	dydx[0] = y[1];
	dydx[1] = 1.0;
}

///////////////////////////////////
// Solve ODE system 1
///////////////////////////////////

void solve_ode1(){
	realT *x, *y, *ic;
	realT xmin, xmax;
	realT yini;
  int n, neq;

  std::cout << "==================" << std::endl;
  std::cout << "   Solve ODE 1    " << std::endl;
  std::cout << "==================" << std::endl;

	// Number of equations
	neq = 1;

	// Initial conditions for ODE1
	ic = new realT[neq];
	xmin = 0.0, xmax = 1.0;
	ic[0] = yini = 1.0;

  // Number of points between xmin and xmax
  // to store the solution in
  n = 20;
	
	// Set up solver for ODE1
	OdeSolver myode(n, neq, ode1);

	// Set initial conditions
	myode.set_initial_conditions(xmin, xmax, ic);

	// Solve
	myode.solve(false);

	// Get pointers to solution
	x = myode.x_array();
	y = myode.y_array(0);

	// Print data
	for(int i=0;i<n;i++){
		std::cout << setw(2) << i << " / " << n << "  x: " << setw(12) << x[i];
		std::cout << "  y: " <<  setw(12) << y[i] << " delta_y: " <<  setw(12) << y[i] - exp(-x[i]*x[i]) << std::endl; 
  }

	delete[] ic;
}

///////////////////////////////////
// Solve ODE system 2
///////////////////////////////////

void solve_ode2(){
	realT *x, *y1, *y2, *ic;
	realT xmin, xmax;
	realT y1ini, y2ini;
  int n, neq;

  std::cout << "==================" << std::endl;
  std::cout << "   Solve ODE 2    " << std::endl;
  std::cout << "==================" << std::endl;

	// Number of equations
	neq = 2;

	// Initial conditions for ODE1
	ic = new realT[neq];
	xmin = 0.0, xmax = 1.0;
  ic[0] = y1ini = 1.0;
	ic[1] = y2ini = 1.0;

  // Number of points between xmin and xmax
  // to store the solution in
  n = 20;
	
	// Set up solver for ODE2
	OdeSolver myode(n, neq, ode2);

	// Set initial conditions
	myode.set_initial_conditions(xmin, xmax, ic);

	// Solve
	myode.solve(false);

	// Get pointers to solution
	x = myode.x_array();
	y1 = myode.y_array(0);
	y2 = myode.y_array(1);

	// Print data
	for(int i=0;i<n;i++){
		std::cout <<  setw(2) << i+1 << " / " << n << "  x: " <<  setw(12) <<x[i];
	  std::cout << "  y1: " <<  setw(12) << y1[i] << " delta_y1: " <<  setw(12) << y1[i] - 1.0 - x[i] - x[i]*x[i]/2.0;
    std::cout << "  y2: " <<  setw(12) << y2[i] << " delta_y2: " <<  setw(12) << y2[i] - 1.0 - x[i] << std::endl; 
  }

	delete[] ic;
}

int main(){
  solve_ode1();
  solve_ode2();
}

