
//==========================================================
// Implementation of perturbations equations for LCDM
//==========================================================

#include "PerturbationsLCDM.h"

// Set initial conditions for the full system
void PerturbationsLCDM::set_initial_conditions_full(double x, VecDoub_I &ytc, VecDoub_IO &y){
  double ck_H = k_current/(BG->H_c(x));
  double ck_Htau = ck_H/REC->get_dtau(x);

  // delta, delta_b, v, v_b, Phi, Theta_0, Theta_1
  for (int i = 0; i<7; i++) {
    y[i] = ytc[i];
  }

  // Theta_2
  y[7] = -8.0/15.0*ck_Htau*ytc[6];

  // Theta_l
  for (int l = 3; l<=lmax_int; l++) {
    y[5+l] = -double(l)/double(2*l+1)*ck_Htau*y[4+l];
  }

  // Theta_p_l
  y[6+lmax_int] = 1.25*y[7];
  y[7+lmax_int] = -0.25*ck_Htau*y[7];
  y[8+lmax_int] = 0.25*y[7];
  for (int l = 3; l<=lmax_int; l++) {
    y[6 + lmax_int + l] = -double(l)/double(2*l+1)*ck_Htau*y[5 + lmax_int + l];
  }

  // Nu_l
  for (int l = 0; l<=lmax_nu; l++) {
    y[7+2*lmax_int+l] = ytc[7+l];
  }

}

// Set initial conditions for the tight coupling regime
void PerturbationsLCDM::set_initial_conditions_tight_coupling(double x, VecDoub_IO &y){
  double ck_H = k_current/(BG->H_c(x));
  double lambda = (1.0+2.0/5.0*cosm_simu->f_nu);

  y[0] = 1.5/lambda;			// delta
  y[1] = y[0];				// delta_b
  y[2] = ck_H/2.0/lambda;		// v
  y[3] = y[2];				// v_b
  y[4] = 1.0;					// Phi
  y[5] = 0.5/lambda;			// Theta_0
  y[6] = -ck_H/6.0/lambda;	// Theta_1
  y[7] = 0.5/lambda;			// Nu_0
  y[8] = -ck_H/6.0/lambda;	// Nu_1
  y[9] = -k_current*k_current*exp(2.0*x)/12.0/cosm_simu->Omega_nu*(1.0-1.0/lambda); // Nu_2

  for (int l =3; l<=lmax_nu; l++) {
    y[7+l] = ck_H/double(2*l+1)*y[6+l]; // Nu_l
  }
}

// The RHS of the ODE system dx_i/deta = RHS_i  in the tight coupling regime
void PerturbationsLCDM::get_derivs_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dydx){
  double Hx, dHH, ddHH, tau, dtau, ddtau, ck_H, ck_Htau, Theta_2, R, q, Psii;

  BG->get_all_H_c(x, Hx, dHH, ddHH);
  REC->get_all_tau(x, tau, dtau, ddtau);

  ck_H = k_current/Hx;
  ck_Htau = ck_H/dtau;
  Theta_2 = -8.0/15.0*ck_Htau*y[6];
  Psii = -y[4]-12.0*exp(-2.0*x)/(k_current*k_current)*(cosm_simu->Omega_r*Theta_2+cosm_simu->Omega_nu*y[9]);

  // The variables are...
  // [0 = delta], [1 = delta_b], [2 = v], [3 = v_b], [4 = Phi]

  // Phi
  dydx[4] = Psii - ck_H*ck_H/3.0*y[4] + 0.5/(Hx*Hx)*(cosm_simu->Omega_m*y[0]*exp(-x)+cosm_simu->Omega_b*y[1]*exp(-x)+4.0*exp(-2.0*x)*(cosm_simu->Omega_r*y[5]+cosm_simu->Omega_nu*y[7]));

  R = cosm_simu->R0*exp(-x);
  q = (-((1.0-R)*dtau+(1.0+R)*ddtau)*(3.0*y[6]+y[3])-ck_H*Psii+(1.0-dHH)*ck_H*(-y[5]+2.0*Theta_2)-ck_H*(-ck_H*y[6]-dydx[4]))/((1.0+R)*dtau+dHH-1.0);

  // delta
  dydx[0] = ck_H*y[2]-3.0*dydx[4];

  // delta_b
  dydx[1] = ck_H*y[3]-3.0*dydx[4];

  // v
  dydx[2] = -y[2]-ck_H*Psii;

  // v_b
  dydx[3] = ( -y[3] - ck_H*Psii + R*(q+ck_H*(-y[5]+2.0*Theta_2) - ck_H*Psii) )/(1.0+R);

  // Theta_0
  dydx[5] = -ck_H*y[6] - dydx[4];

  // Theta 1
  dydx[6] = (q-dydx[3])/3.0;

  // Nu 0
  dydx[7] = -ck_H*y[8] - dydx[4];

  // Nu 1
  dydx[8] = ck_H/3.0*(y[7]-2.0*y[9]+Psii);

  // Nu l
  for (int l = 2; l < lmax_nu; l++) {
    dydx[7+l] = ck_H/double(2*l+1)*(double(l)*y[6+l]-double(l+1)*y[8+l]);
  }

  // Nu lmax_nu
  dydx[7+lmax_nu] = ck_H*y[6+lmax_nu] - double(lmax_nu+1)*y[7+lmax_nu]/(BG->get_eta(x)*Hx);
}

// The Jacobian matrix for the tight coupling ODE
void PerturbationsLCDM::get_jacobian_tight_coupling(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy){
  double expx, expx2, Hx, tau, dtau, ddtau, ck_H, ck_Htau, dHH, R, Theta_2, t1, ddHH;
  double dpsi_4, dpsi_6, dpsi_9;
  double dq_3, dq_4, dq_5, dq_6, dq_9;

  BG->get_all_H_c(x, Hx, dHH, ddHH);
  REC->get_all_tau(x, tau, dtau, ddtau);

  expx     = exp(-x);
  expx2    = expx*expx;
  ck_H     = k_current/Hx;
  ck_Htau  = ck_H/dtau;
  R        = cosm_simu->R0*expx;
  Theta_2  = -8.0/15.0*ck_Htau;

  // Initialize matrix
  for (int i = 0; i<ntight; i++) {
    for (int j = 0; j<ntight; j++) {
      dfdy(i,j) = 0.0;
    }
  }

  // dPsi/dy(j)
  t1 = -12.0*expx2/(k_current*k_current);
  dpsi_4 = -1.0;
  dpsi_6 = t1*cosm_simu->Omega_r*Theta_2;
  dpsi_9 = t1*cosm_simu->Omega_nu;

  // dq / dy(j) without the dPhi/dx terms
  t1 = -1.0/((1.0+R)*dtau+dHH-1.0);
  dq_3 = t1*(1.0+R)*ddtau;
  dq_4 = t1*dpsi_4;
  dq_5 = -t1*(1.0-dHH)*ck_H;
  dq_6 = 3.0*dq_3 + t1*(1.0-dHH)*ck_H*2.0*Theta_2 + t1*ck_H*ck_H + t1*dpsi_6;
  dq_9 = t1*dpsi_9;

  // The (4,j) components
  t1 = 0.5/(Hx*Hx); 
  dfdy(4,0) = t1*cosm_simu->Omega_m*expx;
  dfdy(4,1) = t1*cosm_simu->Omega_b*expx;
  dfdy(4,4) = -ck_H*ck_H/3.0 + dpsi_4;
  dfdy(4,5) = t1*expx2*cosm_simu->Omega_r;
  dfdy(4,6) = dpsi_6;
  dfdy(4,7) = t1*expx2*cosm_simu->Omega_nu;
  dfdy(4,9) = dpsi_9;

  // The (0 - 1 - 5,j) components
  dfdy(0,2) = ck_H;
  dfdy(1,3) = ck_H;
  dfdy(5,6) = -ck_H;
  for (int j = 0; j<ntight; j++) {
    dfdy(0,j) += -3.0*dfdy(4,j);
    dfdy(1,j) += -3.0*dfdy(4,j);
    dfdy(5,j) += -dfdy(4,j);
  }

  // The (2,j) components
  dfdy(2,2) = -1.0;
  dfdy(2,4) = -ck_H*dpsi_4;
  dfdy(2,6) = -ck_H*dpsi_6;
  dfdy(2,9) = -ck_H*dpsi_9;

  // The (3,j) components
  t1 = -1.0/((1.0+R)*dtau+dHH-1.0);
  for (int j = 0; j<ntight; j++) {
    dfdy(3,j) += -R*t1*ck_H*dfdy(4,j)/(1.0+R);
  }
  dfdy(3,3) += -1.0/(1.0+R) + R*dq_3/(1.0+R);
  dfdy(3,4) += -ck_H*dpsi_4/(1.0+R) + R*dq_4/(1.0+R); 
  dfdy(3,5) += -ck_H*R/(1.0+R) + R*dq_5/(1.0+R);
  dfdy(3,6) += 2.0*R/(1.0+R)*Theta_2-ck_H*dpsi_6/(1.0+R) + R*dq_6/(1.0+R);
  dfdy(3,9) += -ck_H*dpsi_9/(1.0+R) + R*dq_9/(1.0+R); 

  // The (6,j) components
  for (int j = 0; j<ntight; j++) {
    dfdy(6,j) += -t1*ck_H*dfdy(4,j)/3.0 - dfdy(3,j)/3.0;
  }
  dfdy(6,3) += dq_3/3.0;
  dfdy(6,4) += dq_4/3.0;
  dfdy(6,5) += dq_5/3.0;
  dfdy(6,6) += dq_6/3.0;
  dfdy(6,9) += dq_9/3.0;

  // The (7,j) components
  for (int j = 0; j<ntight; j++) {
    dfdy(7,j) = -dfdy(4,j);
  }
  dfdy(7,8) += -ck_H;

  // The (8,j) components		
  dfdy(8,4) = ck_H/3.0*dpsi_4;
  dfdy(8,6) = ck_H/3.0*dpsi_6;
  dfdy(8,7) = ck_H/3.0;
  dfdy(8,9) = -2.0*ck_H/3.0+ck_H/3.0*dpsi_9;

  // The (9+l,j) components
  for (int l = 2; l<lmax_nu; l++) {
    dfdy(7+l,6+l) = ck_H/double(2*l+1)*double(l);
    dfdy(7+l,8+l) = -ck_H/double(2*l+1)*double(l+1);
  }

  // The (npar-1, j) components
  dfdy(ntight-1,ntight-2) = ck_H;
  dfdy(ntight-1,ntight-1) = -double(lmax_nu+1)/(BG->get_eta(x)*Hx);
}

// The RHS of the ODE system dx_i/deta = RHS_i 
void PerturbationsLCDM::get_derivs_full_equation(const double x, VecDoub_I &y, VecDoub_O &dydx){
  double expx, Hx, dHH, ddHH, tau, dtau, ddtau, ck_H, ck_Htau, LargePi, Psii, R, etaH;

  BG->get_all_H_c(x, Hx, dHH, ddHH);
  REC->get_all_tau(x, tau, dtau, ddtau);

  expx    = exp(-x);
  ck_H    = k_current/Hx;
  ck_Htau = ck_H/dtau;
  LargePi = y[7] + y[6+lmax_int] + y[8+lmax_int];
  Psii    = -y[4]-12.0*expx*expx/(k_current*k_current)*(cosm_simu->Omega_r*y[7]+cosm_simu->Omega_nu*y[9+2*lmax_int]);
  R       = cosm_simu->R0*expx;
  etaH    = BG->get_eta(x)*Hx;

  // Phi
  dydx[4] = Psii - ck_H*ck_H/3.0*y[4] + 0.5/(Hx*Hx)*((cosm_simu->Omega_m*y[0]+cosm_simu->Omega_b*y[1])*expx+4.0*expx*expx*(cosm_simu->Omega_r*y[5]+cosm_simu->Omega_nu*y[7+2*lmax_int]));

  // delta
  dydx[0] = ck_H*y[2] - 3.0*dydx[4];

  // delta_b
  dydx[1] = ck_H*y[3] - 3.0*dydx[4];

  // v
  dydx[2] = -y[2] - ck_H*Psii;

  // v_b
  dydx[3] = -y[3] - ck_H*Psii + dtau*R*(3.0*y[6]+y[3]);

  // Theta 0
  dydx[5] = -ck_H*y[6] - dydx[4];

  // Theta 1
  dydx[6] = ck_H/3.0*(y[5] - 2.0*y[7] + Psii) + dtau*(y[6]+y[3]/3.0);

  // Theta 2
  dydx[7] = ck_H/5.0*(2.0*y[6]-3.0*y[8]) + dtau*(y[7]-0.1*LargePi);

  // Theta l
  for (int l = 3; l<lmax_int; l++) {
    dydx[5 + l] = ck_H/double(2*l+1) * (double(l)*y[4+l] - double(l+1)*y[6+l]) + dtau*y[5+l];
  }

  // Theta lmax_int
  dydx[5 + lmax_int] = ck_H*y[4+lmax_int] - double(lmax_int+1)*y[5+lmax_int]/etaH + dtau*y[5+lmax_int]; 

  // Theta_p 0
  dydx[6 + lmax_int] = -ck_H*y[7+lmax_int] + dtau*(y[6+lmax_int]-0.5*LargePi);

  // Theta p 1
  dydx[7 + lmax_int] = ck_H/3.0*(y[6+lmax_int]-2.0*y[8+lmax_int]) + dtau*y[7+lmax_int];

  // Theta p 2
  dydx[8 + lmax_int] = ck_H/5.0*(2.0*y[7+lmax_int]-3.0*y[9+lmax_int]) + dtau*(y[8+lmax_int]-0.1*LargePi);

  // Theta p l
  for (int l = 3; l<lmax_int; l++) {
    dydx[6 + lmax_int + l] = ck_H/double(2*l+1)*(double(l)*y[5+lmax_int+l]-double(l+1)*y[7+lmax_int+l])+dtau*y[6+lmax_int+l];
  }

  // Theta p lmax_int
  dydx[6 + 2*lmax_int] = ck_H*y[5+2*lmax_int] -double(lmax_int+1)*y[6+2*lmax_int]/etaH;

  // nu 0
  dydx[7 + 2*lmax_int] = -ck_H*y[8+2*lmax_int]-dydx[4];

  // nu 1
  dydx[8 + 2*lmax_int] = ck_H/3.0*(y[7+2*lmax_int]-2.0*y[9+2*lmax_int]+Psii);

  // nu l
  for (int l = 2; l<lmax_nu; l++) {
    dydx[7 + 2*lmax_int + l] = ck_H/double(2*l+1)*(double(l)*y[6+2*lmax_int+l]-double(l+1)*y[8+2*lmax_int+l]);
  }

  // nu lmax_nu
  dydx[7+2*lmax_int+lmax_nu] = ck_H*y[6+2*lmax_int+lmax_nu]-double(lmax_nu+1)*y[7+2*lmax_int+lmax_nu]/etaH;		
}

// Jacobian matrix for the full perturbation ODE system
void PerturbationsLCDM::get_jacobian_full_equation(const double x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy){

  double expx, Hx, ck_H, dtau, R, dpsi_4, dpsi_7, dpsi_nu2, t1, t2, etaH;

  expx  = exp(-x);
  Hx    = BG->H_c(x);
  ck_H  = k_current/Hx;
  dtau  = REC->get_dtau(x);
  etaH  = BG->get_eta(x)*Hx;

  // Initialize array
  for (int i = 0; i<npar; i++) {
    for (int j = 0; j<npar; j++) {
      dfdy(i,j) = 0.0;
    }
  }

  // The d(Phi) / dyj components
  t1 = -12.0*expx*expx/(k_current*k_current);
  dpsi_4   = -1.0;
  dpsi_7   = t1*cosm_simu->Omega_r;
  dpsi_nu2 = t1*cosm_simu->Omega_nu;

  t2 = 0.5/(Hx*Hx)*expx;
  dfdy(4,0)            = t2*cosm_simu->Omega_m;	// delta
  dfdy(4,1)            = t2*cosm_simu->Omega_b;	// delta_b
  dfdy(4,4)            = dpsi_4 - ck_H*ck_H/3.0;		// Phi
  dfdy(4,7)            = dpsi_7;						// Theta_2
  dfdy(4,9+2*lmax_int) = dpsi_nu2;					// Nu_2
  dfdy(4,5)            = 4.0*expx*t2*cosm_simu->Omega_r;			// Theta_0
  dfdy(4,7+2*lmax_int) = 4.0*expx*t2*cosm_simu->Omega_nu;		// Nu_0

  // The d(delta) / dyj components
  for (int j = 0; j<10+2*lmax_int; j++) { // Loop only as long as we need
    dfdy(0,j) = -3.0*dfdy(4,j);						// dPhi
  }
  dfdy(0,2) += ck_H;									// v

  // The d(deltab) / dyj  components
  for (int j = 0; j<10+2*lmax_int; j++) { // Loop only as long as we need
    dfdy(1,j) = -3.0*dfdy(4,j);						// dPhi
  }
  dfdy(1,3) += ck_H;									// v_b

  // The d(v) / dyj components
  dfdy(2,2)            = -1.0 + 3.0*R*dtau;			// v
  dfdy(2,4)            = -ck_H*dpsi_4;				// Phi      (from Psi)
  dfdy(2,7)            = -ck_H*dpsi_7;				// Theta_2  (from Psi)
  dfdy(2,9+2*lmax_int) = -ck_H*dpsi_nu2;				// Nu_2		(from Psi)

  // The d(v_b) / dyj  components
  dfdy(3,3)            = -1.0;						// v_b
  dfdy(3,4)            = -ck_H*dpsi_4;				// Phi      (from Psi)
  dfdy(3,6)            = 3.0*R*dtau;					// Theta_1
  dfdy(3,7)            = -ck_H*dpsi_7;				// Theta_2  (from Psi)
  dfdy(3,9+2*lmax_int) = -ck_H*dpsi_nu2;				// Nu_2		(from Psi)

  // The d(Theta_0) / dyj components
  for (int j = 0; j<=9+2*lmax_int; j++) {
    dfdy(5,j) = -dfdy(4,j);							// dPhi
  }
  dfdy(5,6) += -ck_H;									// Theta_1

  // The d(Theta_1) / dyj components
  dfdy(6,3)            = dtau/3.0;						// v_b
  dfdy(6,4)            = ck_H*dpsi_4/3.0;					// Phi        (from Psi)
  dfdy(6,5)            = ck_H/3.0;						// Theta_0
  dfdy(6,6)            = dtau;							// Theta_1
  dfdy(6,7)            = ck_H*dpsi_7/3.0-2.0*ck_H/3.0;	// Theta_2  (+ from Psi)
  dfdy(6,9+2*lmax_int) = ck_H*dpsi_nu2/3.0;				// Nu_2       (from Psi)

  // The d(Theta_2) / dyj components
  dfdy(7,6)            = 2.0*ck_H/5.0;	// Theta_1
  dfdy(7,7)            = dtau - 0.1;		// Theta_2	  (+ from LargePi)
  dfdy(7,8)            = -3.0*ck_H/5.0;	// Theta_2
  dfdy(7,6+lmax_int)   = - 0.1;			// Theta_p_0  (from LargePi)
  dfdy(7,8+lmax_int)   = - 0.1;			// Theta_p_2  (from LargePi)

  // The d(Theta_l) / dyj components
  for (int l = 3; l<lmax_int; l++) {
    dfdy(5 + l,4+l) = ck_H/double(2*l+1)*double(l);		// Theta_(l-1)
    dfdy(5 + l,5+l) = dtau;								// Theta_l
    dfdy(5 + l,6+l) = -ck_H/double(2*l+1)*double(l+1);	// Theta_(l+1)
  }

  // The  d(Theta_lmax) / dyj components
  dfdy(5+lmax_int,4+lmax_int) = ck_H;								// Theta_(lmax-1)
  dfdy(5+lmax_int,5+lmax_int) = - double(lmax_int+1)/etaH + dtau;	// Theta_(lmax)


  // The  d(Theta_p_0) / dyj components
  dfdy(6+lmax_int,7)          = -0.5*dtau;
  dfdy(6+lmax_int,6+lmax_int) = 0.5*dtau; 
  dfdy(6+lmax_int,7+lmax_int) = -ck_H; 
  dfdy(6+lmax_int,8+lmax_int) = -0.5*dtau;

  // The  d(Theta_p_1) / dyj components
  dfdy(7+lmax_int,6+lmax_int) = ck_H/3.0;;
  dfdy(7+lmax_int,7+lmax_int) = dtau; 
  dfdy(7+lmax_int,8+lmax_int) = -2.0*ck_H/3.0;

  // The  d(Theta_p_2) / dyj components
  dfdy(8+lmax_int,7)          = -0.1*dtau;
  dfdy(8+lmax_int,6+lmax_int) = -0.1*dtau;  
  dfdy(8+lmax_int,7+lmax_int) = 2.0*ck_H/5.0;
  dfdy(8+lmax_int,8+lmax_int) = 0.9*dtau;
  dfdy(8+lmax_int,9+lmax_int) = -3.0*ck_H/5.0;

  // The  d(Theta_p_l) / dyj components
  for (int l = 3; l<lmax_int; l++) {
    dfdy(6+lmax_int + l,5+lmax_int + l) = ck_H/double(2*l+1)*double(l);
    dfdy(6+lmax_int + l,6+lmax_int + l) = dtau;
    dfdy(6+lmax_int + l,7+lmax_int + l) = -ck_H/double(2*l+1)*double(l+1);
  }

  // The  d(Theta_p_lmax) / dyj components
  dfdy(6+2*lmax_int,5+2*lmax_int) = ck_H;
  dfdy(6+2*lmax_int,6+2*lmax_int) = -double(lmax_int+1)/etaH;

  // The  d(Nu_0) / dyj components
  for (int j = 0; j<= 9+2*lmax_int; j++) {
    dfdy(7+2*lmax_int,j) = -dfdy(4,j); 
  }
  dfdy(7+2*lmax_int,8+2*lmax_int) += -ck_H;

  // The (8+2*lmax_int,j) components
  dfdy(8+2*lmax_int,4)            = ck_H/3.0*dpsi_4;
  dfdy(8+2*lmax_int,7)            = ck_H/3.0*dpsi_7;
  dfdy(8+2*lmax_int,9+2*lmax_int) = ck_H/3.0*(dpsi_nu2 - 2.0);
  dfdy(8+2*lmax_int,7+2*lmax_int) = ck_H/3.0;

  // The (8+2*lmax_int+l,j) components
  for (int l = 2; l<lmax_nu; l++) {
    dfdy(7 + 2*lmax_int + l,6+2*lmax_int+l) =  ck_H/double(2*l+1)*double(l);
    dfdy(7 + 2*lmax_int + l,8+2*lmax_int+l) = -ck_H/double(2*l+1)*double(l+1);
  }

  // The (8+2*lmax_int+lmax_nu,j) components
  dfdy(7+2*lmax_int+lmax_nu,6+2*lmax_int+lmax_nu) = ck_H;
  dfdy(7+2*lmax_int+lmax_nu,7+2*lmax_int+lmax_nu) = -double(lmax_nu+1)/etaH;

}
