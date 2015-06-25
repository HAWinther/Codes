
//=========================================================
// Methods to compute standard recombination history
// Standard means no change to the physics of rec. so
// same equations as in LCDM (but H(a) etc. are of course
// that of the implemented model)
//=========================================================

#include "RecombinationSTANDARD.h"

// Calculate electron number density by integrating the Peebles ODE
void RecombinationSTANDARD::calc_n_e(){
  double a, T_b, n_b, fac, z_now, reion_fact, factor, x1, x2, xH, f_e, f_e_old, fac_1, fac_2, fac_H, Y_p, T_0;
  bool saha = true;

  Y_p = cosm_simu->Y_p;
  T_0 = cosm_simu->T_0;

  // Integration variables
  const double rtol = 1.0e-10;
  const double atol = 1.0e-10;
  const double h1   = 1.0e-5;

  // Initialize arrays
  X_e  = VecDoub(nmain);
  n_e  = VecDoub(nmain);

  // Initialize Peebles ODE
  derivs_peebles dpeebles(this);
  VecDoub y(1);
  Output out;

  // Calculate recombination history
  for (int i = 0; i<nmain; i++) {

    // Use Saha equation
    if (saha) {
      a = exp(BG->x_main[i]);
      T_b = k_b*T_0/a;
      n_b = (cosm_simu->Omega_b*cosm_simu->rhoc)/(m_H*a*a*a);
      factor = pow(m_e*T_b/(hbar*hbar*2.0*pi),1.5)/n_b;

      // Start by assuming f_e = 1.0
      f_e     = 1.0;
      f_e_old = 0.0;

      // Iterate until f_e is found to a good accuracy
      while (abs(f_e-f_e_old) > 1.0e-10) {

        // R.h.s of the Saha equations for He++, He+ and H+ respectivily
        fac_1 = 2.0*factor*exp(-xhi0/T_b)/f_e;
        fac_2 = 4.0*factor*exp(-xhi1/T_b)/f_e;
        fac_H = factor*exp(-epsilon_0/T_b)/f_e;

        // Abundances of He++, He+ and H+ respectivily
        x1 = fac_1/(1.0 + fac_1 + fac_1*fac_2);
        x2 = fac_2 * x1;
        xH = fac_H/(1.0 + fac_H);

        // Calculate new f_e
        f_e_old = f_e;
        f_e = (2.0*x2+x1)*Y_p/4.0 + (1.0-Y_p)*xH;
      }

      // Electron fraction and number density
      X_e[i] = f_e/(1.0-Y_p);
      n_e[i] = f_e*n_b;

      // Test if Saha approximation is still valid
      if (X_e[i] < saha_limit){
        saha = false;
        y[0] = X_e[i];
      }
    }else{
      // Use Peebles equation
      a = exp(BG->x_main[i]);
      n_b = (cosm_simu->Omega_b*cosm_simu->rhoc)/(m_H*a*a*a);

      // Integrate Peebles equation
      Odeint<StepperBS<derivs_peebles> > ode_peebles(y, BG->x_main[i-1], BG->x_main[i], atol, rtol, h1, 0.0,out, dpeebles);
      ode_peebles.integrate();

      // Electron fraction and number density
      X_e[i] = y[0];
      n_e[i] = X_e[i]*(1.0-Y_p)*n_b;
    }
  }

  // Change n_e to log(n_e)
  for (int i = 0; i<nmain; i++) {
    if (n_e[i] > 0.0) {
      n_e[i] = log(n_e[i]);
    }
  }

  // Spline n_e
  n_e_spline = Spline(BG->x_main, n_e, 1.0e30, 1.0e30);
  //n_e_spline.set_linear_x();
  n_e_spline.set_three_regions_x(np1, np2);

}

// Right hand side of the Peebles ODE dX/dx = RHS
double RecombinationSTANDARD::get_derivs_peebles(double x, VecDoub_I &y){
  double T_b, n_b, phi_2, alpha_2, beta, beta_2, n1s, lambda_a, lambda_2s1s, c_r, delta, H1, Y_p;

  Y_p = cosm_simu->Y_p;
  T_b = k_b*cosm_simu->T_0*exp(-x);
  n_b = (cosm_simu->Omega_b*cosm_simu->rhoc)*exp(-3.0*x)/m_H;
  H1  = cosm_simu->H_0*BG->H(x);

  phi_2 = 0.448*log(epsilon_0/T_b);
  alpha_2 = phi_2*(8.0/sqrt(3.0*pi))*sigma_T*sqrt(epsilon_0/T_b);				// Dim metre^2
  beta = alpha_2*c*pow(m_e*T_b/(2.0*pi*hbar*hbar),1.5)*exp(-epsilon_0/T_b);	// Dim 1/sec

  delta = ((1.0-y[0])*beta/(c*n_b*alpha_2*y[0]*y[0]));

  if (delta > 1.0e-20) {
    beta_2 = beta*exp(0.75*epsilon_0/T_b);									// Dim 1/sec
    n1s = (1.0-y[0])*(1.0-Y_p)*n_b;											// Dim 1/metre^3
    lambda_a = H1*pow(3.0*epsilon_0/(hbar*c),3)*pow(8*pi,-2)/n1s;			// Dim 1/sec
    lambda_2s1s = 8.227;													// Dim 1/sec
    c_r = (lambda_2s1s+lambda_a)/(lambda_2s1s+lambda_a+beta_2);			
    return c_r/H1*(beta*(1.0-y[0])-c*n_b*(1.0-Y_p)*alpha_2*y[0]*y[0]);
  }else {
    return -1.0/H1*c*n_b*(1.0-Y_p)*alpha_2*y[0]*y[0];
  }

}
