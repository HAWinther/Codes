
//=================================================================
// Methods to compute power-spectra P(k), C_l, etc.
//=================================================================

#include "Powerspectrum.h"
#include "Sphbessel.h"

// Initialize the power-spectrum object
Powerspectrum::Powerspectrum(Background *BG, Recombination *REC, Perturbations *P, ParamSimu *param){

#ifdef VERBOSE
  cout << endl;
  cout << "**********************************************" << endl;
  cout <<	"*** Initialize Power-spectrum calculation  ***" << endl;
  cout <<	"**********************************************" << endl;
#endif

  // Initialize pointers to background, recombination etc.
  this->P   = P;
  this->BG  = BG;
  this->REC = REC;
  cosm_simu = param;

  // Initialize high resolution k grid
  k_high = VecDoub(n_k_high);

  for (int k = 0; k<n_k_high; k++) {
    k_high[k] = k_min*pow(k_max/k_min, k/double(n_k_high-1));
  }

  // Initialize high resolution x grid
  x_high = VecDoub(ntot_high);
  if (cosm_simu->z_reion > 0.0) {
    // Add extra points during reionization period if z_reion > 0
    for (int i = 0; i<ntot_high; i++) {
      if (i < n1_high) {
        x_high[i] = x_start_rec + (x_end_rec-x_start_rec) * double(i)/double(n1_high);
      }else if(i < n1_high + n2_high) {
        x_high[i] = x_end_rec + (cosm_simu->x_start_reion - x_end_rec) * double(i-n1_high)/double(n2_high);
      }else if(i < n1_high + n2_high + n3_high) {
        x_high[i] = cosm_simu->x_start_reion + (cosm_simu->x_end_reion - cosm_simu->x_start_reion) * double(i-n1_high-n2_high)/double(n3_high);
      }else {
        x_high[i] = cosm_simu->x_end_reion + (0.0 - cosm_simu->x_end_reion) * double(i-n1_high-n2_high-n3_high)/double(n4_high-1);
      }
    }
  } else {
    for (int i = 0; i<ntot_high; i++) {
      if (i < n1_high + n2_high) {
        x_high[i] = x_start_rec + (x_end_rec-x_start_rec) * double(i)/double(n1_high + n2_high);
      }else{
        x_high[i] = x_end_rec + (0.0 - x_end_rec) * double(i-n1_high-n2_high)/double(ntot_high - 1 - n1_high - n2_high);
      }
    }
  }

  // Calculate Bessel functions
  calc_bessel();	
}

// Calculate the matter power-spectrum
void Powerspectrum::calc_matter_powerspectrum(){
  double transfer, know, max, sigma8;
  matterpow = VecDoub(n_k);

  // Late time growth-factor! xxx NB: Only valid for LCDM xxx
  // growth = (P->deltaM(n_k-1,ntot-1) - P->deltaM(n_k-1,ntot-2))/(P->x_rec[ntot-1]-P->x_rec[ntot-2]) / P->deltaM(n_k-1,ntot-1);

  // Compute the matter powerspectrum
  for (int k = 0; k<n_k; k++) {
    know = P->ks[k]; 

    // P(k) = <O_m*delta_m + O_b*delta_b / O_m+ O_b>^2 Primordial(k) in units (h/Mpc)^3
    matterpow[k] = 8.0 * pi * pi / 9.0 * know * pow(P->Phi(k,ntot-1),2)/pow(cosm_simu->Omega_m + cosm_simu->Omega_b,2) * get_primordial_power_spectrum(know) * pow(2998.0,3);
  }

  // Spline matter power-spectrum
  matterpow_spline = Spline(P->ks, matterpow, 1.0e30, 1.0e30);
  matterpow_spline.set_log_x();

  // Sigma8 normalization
  if (cosm_simu->norm_cl == 1) {

    sigma8 = cosm_simu->sigma8;
    calc_sigma8(); 
    cosm_simu->norm = pow(sigma8 / cosm_simu->sigma8,2);
    cosm_simu->sigma8 = sigma8;

#ifdef VERBOSE
    cout << endl << "     Applying sigma8 normalization" << endl;
    cout <<			"     sigma8                     : " << cosm_simu->sigma8 << endl;
    cout <<			"     A_s Norm                   : " << cosm_simu->norm << endl;
#endif
  }
}

double Powerspectrum::get_matter_power_spectrum(double k){
  return (cosm_simu->norm*matterpow_spline.get_spline(k));
}

void Powerspectrum::integrate_thetal_T(){
  double a_int, r_int, h1_int;
  double arg_low, arg_mid, eta0 = BG->get_eta(0.0);
  int llow = 17, lmid = 45; // Three regimes
  Output out;

  // Integration variablesx
  a_int  = 1.0e-6; // Increase to 8-10 to get better accuracy for l = 1000-2000
  r_int  = 1.0e-6; 
  h1_int = 1.0e-3;

  derivs_thetal_T dthetal_T_low(this,llow);
  derivs_thetal_T dthetal_T_mid(this,lmid);
  derivs_thetal_T dthetal_T_all(this,l_num);

  VecDoub ylow = VecDoub(llow);
  VecDoub ymid = VecDoub(lmid);
  VecDoub yall = VecDoub(l_num);

  Theta_l_T = MatDoub(l_num,n_k_high);
  for (int l = 0; l<l_num; l++)
    for (int k = 0; k<n_k_high; k++)
      Theta_l_T(l,k) = 0.0;

  for (int l = 0; l<l_num; l++) yall(l) = 0.0;
  Odeint<StepperBS<derivs_thetal_T> > ode_thetal_T_all(yall,x_start_rec,0.0,a_int,r_int,h1_int,0.0,out,dthetal_T_all);

  for (int k = 0; k<n_k_high; k++) {
    k_current = k;

    arg_low = k_high[k]*eta0 / ls[llow];
    arg_mid = k_high[k]*eta0 / ls[lmid];

    if (arg_low < 1.0) {
      // Only integrate for l<llow. Theta_l = 0 for l > llow
      for (int l = 0; l<llow; l++)  ylow(l) = 0.0;
      Odeint<StepperBS<derivs_thetal_T> > ode_thetal_T_low(ylow,x_start_rec,0.0,a_int,r_int,h1_int,0.0,out,dthetal_T_low);
      ode_thetal_T_low.integrate();
      for (int l = 0; l<llow; l++)  Theta_l_T(l,k) = ylow(l);
    } else if (arg_mid < 1.0) {
      // Only integrate for l<lmid. Theta_l = 0 for l > lmid
      for (int l = 0; l<lmid; l++)  ymid(l) = 0.0;
      Odeint<StepperBS<derivs_thetal_T> > ode_thetal_T_mid(ymid,x_start_rec,0.0,a_int,r_int,h1_int,0.0,out,dthetal_T_mid);
      ode_thetal_T_mid.integrate();
      for (int l = 0; l<lmid; l++)  Theta_l_T(l,k) = ymid(l);
    } else {
      // Integrate all Theta_l
      for (int l = 0; l<l_num; l++) yall(l) = 0.0;			
      Odeint<StepperBS<derivs_thetal_T> > ode_thetal_T_all(yall,x_start_rec,0.0,a_int,r_int,h1_int,0.0,out,dthetal_T_all);
      ode_thetal_T_all.integrate();			
      for (int l = 0; l<l_num; l++) Theta_l_T(l,k) = yall(l);
    }
  }

  // Spline Theta_l_T
  thetal_T_spline = Spline(k_high, Theta_l_T, 1.0e30, 1.0e30, false);
  thetal_T_spline.set_log_x();
}

void Powerspectrum::integrate_thetal_E(){
  double a_int, r_int, h1_int;
  Output out;

  // Integration variables
  a_int  = 1.0e-6;
  r_int  = 1.0e-6;
  h1_int = 1.0e-3;

  derivs_thetal_E dthetal_E(this);

  // Make integration and Theta_l_E array
  VecDoub yint = VecDoub(l_num);
  Theta_l_E = MatDoub(l_num,n_k_high);

  // Initialize Theta_l_E
  for (int l = 0; l<l_num; l++) {
    yint[l] = 0.0;
  }

  for (int k = 0; k<n_k_high; k++) {
    k_current = k;

    // Integrate dTheta_l_E/dx
    Odeint<StepperBS<derivs_thetal_E> > ode_thetal_E(yint,x_start_rec,0.0,a_int,r_int,h1_int,0.0,out,dthetal_E);
    ode_thetal_E.integrate();

    // Save Theta_l_E
    for (int l = 0; l<l_num; l++) {
      Theta_l_E(l,k) = yint[l];
      yint[l] = 0.0;
    }
  }

  // Spline Theta_l_E
  thetal_E_spline	= Spline(k_high, Theta_l_E, 1.0e30, 1.0e30, false);
  thetal_E_spline.set_log_x();
}

void Powerspectrum::integrate_cl_TT(){
  double a_int, r_int, h1_int;
  Output out;

  // Integration variables
  a_int  = 1.0e-10;
  r_int  = 1.0e-10;
  h1_int = 1.0e-8;

  derivs_cl_TT dcl_TT(this);

  // Initialize Cl
  VecDoub yint = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    yint[l] = 0.0;
  }

  // Integrate dCl/dk
  Odeint<StepperBS<derivs_cl_TT> > ode_cl(yint,k_min,k_max,a_int,r_int,h1_int,0.0,out,dcl_TT);
  ode_cl.integrate();	

  // Save Cl
  cls_TT = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    cls_TT[l] = ls[l] * (1.0+ls[l]) / (2.0*pi) * yint[l];
  }

  // Spline C_l
  cl_TT_spline = Spline(lsd, cls_TT, 1.0e30, 1.0e30);
}

void Powerspectrum::integrate_cl_EE(){
  double a_int, r_int, h1_int;
  Output out;

  // Integration variables
  a_int  = 1.0e-10;
  r_int  = 1.0e-10;
  h1_int = 1.0e-8;

  derivs_cl_EE dcl_EE(this);

  // Initialize Cl
  VecDoub yint = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    yint[l] = 0.0;
  }

  // Integrate dCl/dk
  Odeint<StepperBS<derivs_cl_EE> > ode_cl(yint,k_min,k_max,a_int,r_int,h1_int,0.0,out,dcl_EE);
  ode_cl.integrate();	

  // Save Cl
  cls_EE = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    cls_EE[l] = ls[l] * (1.0+ls[l]) / (2.0*pi) * yint[l];
  }

  // Spline C_l
  cl_EE_spline = Spline(lsd, cls_EE, 1.0e30, 1.0e30);
}

void Powerspectrum::integrate_cl_TE(){
  double a_int, r_int, h1_int;
  Output out;

  // Integration variables
  a_int  = 1.0e-10;
  r_int  = 1.0e-10;
  h1_int = 1.0e-8;

  derivs_cl_TE dcl_TE(this);

  // Initialize Cl
  VecDoub yint = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    yint[l] = 0.0;
  }

  // Integrate dCl/dk
  Odeint<StepperBS<derivs_cl_TE> > ode_cl(yint,k_min,k_max,a_int,r_int,h1_int,0.0,out,dcl_TE);
  ode_cl.integrate();	

  // Save Cl
  cls_TE = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    cls_TE[l] = ls[l] * (1.0+ls[l]) / (2.0*pi) * yint[l];
  }

  // Spline C_l
  cl_TE_spline = Spline(lsd, cls_TE, 1.0e30, 1.0e30);

}

void Powerspectrum::integrate_cl_all(){
  double a_int, r_int, h1_int, lfact;
  Output out;

  // Integration variables
  a_int  = 1.0e-10;
  r_int  = 1.0e-10;
  h1_int = 1.0e-8;

  derivs_cl_all dcl_all(this);

  // Initialize Cl
  VecDoub yint = VecDoub(3*l_num);
  for (int l = 0; l<3*l_num; l++) {
    yint[l] = 0.0;
  }

  // Integrate dCl/dk
  Odeint<StepperBS<derivs_cl_all> > ode_cl(yint,k_min,k_max,a_int,r_int,h1_int,0.0,out,dcl_all);
  ode_cl.integrate();	

  // Save Cl
  cls_TT = VecDoub(l_num);
  cls_TE = VecDoub(l_num);
  cls_EE = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    lfact = ls[l] * (1.0+ls[l]) / (2.0*pi);
    cls_TT[l] = lfact * yint[l];
    cls_TE[l] = lfact * yint[l+l_num];
    cls_EE[l] = lfact * yint[l+2*l_num];
  }

  // Spline C_l
  cl_TT_spline = Spline(lsd, cls_TT, 1.0e30, 1.0e30);
  cl_TE_spline = Spline(lsd, cls_TE, 1.0e30, 1.0e30);
  cl_EE_spline = Spline(lsd, cls_EE, 1.0e30, 1.0e30);

}

// Calculate all scalar power-spectrums
void Powerspectrum::calc_scalar_powerspectrums(){
  int lmax, kkk; 
  double clmax;

#ifdef VERBOSE
  cout << endl; 
  cout << "**********************************************" << endl;
  cout <<	"***   Calculate Scalar Power-spectrums     ***" << endl;
  cout <<	"**********************************************" << endl;
#endif

  // Initialize arrays of l's and Cl's for Cl evaluation
  lsd = VecDoub(l_num);
  for (int l = 0; l<l_num; l++) {
    lsd[l] = double(ls[l]);
  }

#ifdef VERBOSE
  cosm_simu->time_now = clock();
  cout << endl << "---> Calculate matter power-spectrum" << endl;
  calc_matter_powerspectrum();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<			"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#else
  calc_matter_powerspectrum();
#endif


#ifdef VERBOSE
  cosm_simu->time_now = clock();
  cout << endl;
  cout << "---> Calculate High resolution T-source" << endl;
  calc_high_res_source_T();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;

  cosm_simu->time_now = clock();
  cout << endl;
  cout << "---> Calculate High resolution E-source" << endl;
  calc_high_res_source_E();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#else
  calc_high_res_source_T();
  calc_high_res_source_E();
#endif

  // Source at a given k-value (here for the middle point)
  kkk = n_k_high/2;
  ofstream outsource;
  string outsource_filename = cosm_simu->outprefix + "powerspectrum_ThetaTTintegrand.dat";
  outsource.open(outsource_filename.c_str());
  for (int i = 0; i<ntot_high; i++) {
    outsource << x_high[i] << " " << get_high_source_T(x_high[i], kkk)*get_bessel(k_high[kkk]*(BG->get_eta(0.0)-BG->get_eta(x_high[i])), 17) << endl;
  }
  outsource.close();

#ifdef VERBOSE
  cosm_simu->time_now = clock();
  cout << endl;
  cout << "---> Calculate Theta_l_T" << endl;
  integrate_thetal_T();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;

  cosm_simu->time_now = clock();
  cout << endl;
  cout << "---> Calculate Theta_l_E" << endl;
  integrate_thetal_E();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#else
  integrate_thetal_T();
  integrate_thetal_E();
#endif

  ofstream outtheta;
  string outtheta_filename = cosm_simu->outprefix + "powerspectrum_ThetaTT.dat";
  outtheta.open(outtheta_filename.c_str());
  for (int l = 0; l<44; l++) {
    for (int k = 0; k<n_k_high; k++) {
      if (k_high[k] < 1000.0) {
        outtheta << k_high[k] << " " <<  get_thetal_T(k_high[k], l) << endl;
      }
    }
    outtheta << endl;
  }
  outtheta.close();

#ifdef VERBOSE
  cosm_simu->time_now = clock();
  cout << endl;
  cout << "---> Calculate Cl_TT, Cl_EE and CL_TE" << endl;
  integrate_cl_all();
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#else
  integrate_cl_TT();
  integrate_cl_EE();
  integrate_cl_TE();
#endif

  // First-peak normalization
  lmax = 200;
  find_cl_max(lmax, clmax);	
  if (cosm_simu->norm_cl == 0) {
    cosm_simu->norm = cosm_simu->first_peak / clmax;
  } else {
    cosm_simu->first_peak = clmax;
  }

#ifdef VERBOSE
  cout << endl << "---> After normalization we have" << endl;
  cout <<			"     A_s                        : " << cosm_simu->A_s*cosm_simu->norm << endl;
  calc_sigma8();
  cout <<			"     sigma8                     : " << cosm_simu->sigma8 << endl;
  cout <<			"     C_l(220)                   : " << get_cl_TT(220.0) << endl;
#endif

#ifdef VERBOSE
  cout << endl << "---> Output C_l" << endl;
  cout <<			"     Filename: '" << cosm_simu->outprefix << "powerspectrum_cl.dat'" << endl;
#endif

#ifndef NO_OUTPUT_CL
  // Output normalized C_l togeather with WMAP and error
  ofstream out_cl;
  string out_cl_filename = cosm_simu->outprefix + "powerspectrum_cl.dat";
  out_cl.open(out_cl_filename.c_str());
  for (int l = 2; l<ls[l_num-1]; l++) {
    out_cl << l << " " << get_cl_TT(l) << " " << get_cl_EE(l) << " " << get_cl_TE(l) << endl;
  }
  out_cl.close();
#endif

#ifndef NO_OUTPUT_MATTERPOW
#ifdef VERBOSE
  cout << endl;
  cout << "---> Output Matter Power Spectrum" << endl;
  cout <<	"     Filename: '" << cosm_simu->outprefix << "powerspectrum_pofk.dat'" << endl;
#endif

  // Output normalized matter powerspectrum
  double know, normf;
  ofstream out_mpow;
  string out_mpow_filename = cosm_simu->outprefix + "powerspectrum_pofk.dat";
  out_mpow.open(out_mpow_filename.c_str());
  for (int k = 0; k<n_k; k++) {
    know = P->ks[k]/2998.0;
    normf = get_matter_power_spectrum(P->ks[n_k-1])/(pow(P->deltaM(n_k-1,ntot-1),2)/pow(P->ks[n_k-1]/2998.0,3)*get_primordial_power_spectrum(P->ks[n_k-1]));
    out_mpow << know << " " << get_matter_power_spectrum(P->ks[k]) << " " << normf*pow(P->deltaM(k,ntot-1),2)/pow(know,3)*get_primordial_power_spectrum(P->ks[k]) << endl;
  }
  out_mpow.close();
#endif

}

double Powerspectrum::get_high_source_T(double x, int k){
  return high_res_spline_T.get_spline_matrix(x, k);
} 

double Powerspectrum::get_high_source_E(double x, int k){
  return high_res_spline_E.get_spline_matrix(x, k);
}

double Powerspectrum::get_thetal_T(double x, int l){
  return thetal_T_spline.get_spline_matrix(x, l);
}

double Powerspectrum::get_thetal_E(double x, int l){
  return thetal_E_spline.get_spline_matrix(x, l);
}

double Powerspectrum::get_cl_TT(double l){
  return (cosm_simu->norm * cl_TT_spline.get_spline(l));
}

double Powerspectrum::get_cl_EE(double l){
  return (cosm_simu->norm * cl_EE_spline.get_spline(l));
}

double Powerspectrum::get_cl_TE(double l){
  return (cosm_simu->norm * cl_TE_spline.get_spline(l));
}

// Find maximum of C_l (for first peak normalization)
void Powerspectrum::find_cl_max(int &lmax, double &clmax){
  bool found = false;
  double clnow;

  clmax = get_cl_TT(double(lmax));

  while (!found) {
    lmax++;
    clnow = get_cl_TT(1.0*lmax);
    if (clnow > clmax){
      clmax = clnow;
    } else {
      found = true;
    }
    if (lmax > 300) {
      lmax = 220;
      cout << "     Error in find_cl_max: C_l has maximum outside 200 < l < 300" << endl;
      cout << "     Proceding by putting maximum to l = 220" << endl;
      found = true;
    }
  }
}

// Using the definition P(k) = k^3/2pi^2 P'(k) = A_s (k/k_pivot)^(n_s - 1)
double Powerspectrum::get_primordial_power_spectrum(double k){
  //return (cosm_simu->A_s * pow(k * cosm_simu->hubb / k_pivot , cosm_simu->n_s-1.0));
  return (cosm_simu->A_s * pow(k, cosm_simu->n_s-1.0));
}

// Calculate the high res source function(s)
void Powerspectrum::calc_high_res_source_T(){

  // Initialize high resolution source function
  S_high_T_k   = MatDoub(n_k_high,ntot);

  // Spline across k
  low_res_spline_T = Spline(P->ks, P->S_low, 1.0e30, 1.0e30, true);
  low_res_spline_T.set_log_x();

  // Calculate high resolution grid using low resolution spline
  for (int k = 0; k<n_k_high; k++) {
    for (int i = 0; i<ntot; i++) {
      S_high_T_k(k,i) = low_res_spline_T.get_spline_matrix(k_high[k], i);
    }
  }

  // Spline across x
  high_res_spline_T = Spline(P->x_rec, S_high_T_k, 1.0e30, 1.0e30, false);
  if (cosm_simu->z_reion>0.0) {
    high_res_spline_T.set_four_regions_x(n1, n2, n3);
  }
}

void Powerspectrum::calc_high_res_source_E(){

  // Initialize E-mode high resolution source function
  S_high_E_k = MatDoub(n_k_high,ntot);

  // Spline across k
  low_res_spline_E = Spline(P->ks, P->SE_low, 1.0e30, 1.0e30, true);
  low_res_spline_E.set_log_x();

  // Calculate high resolution grid using low resolution spline
  for (int k = 0; k<n_k_high; k++) {
    for (int i = 0; i<ntot; i++) {
      S_high_E_k(k,i) = low_res_spline_E.get_spline_matrix(k_high[k], i);
    }
  }

  // Spline across x
  high_res_spline_E = Spline(P->x_rec, S_high_E_k, 1.0e30, 1.0e30, false);
  if (cosm_simu->z_reion>0.0) {
    high_res_spline_E.set_four_regions_x(n1, n2, n3);
  }
}

// Calculate spherical bessel-functions
void Powerspectrum::calc_bessel(){

#ifdef VERBOSE
  cout << endl << "---> Calculate Bessel-functions" << endl;	
  cosm_simu->time_now = clock();
#endif

  int lold, nold, ll, ii;
  double jlx;

  // Initialize j_l
  j_l  = MatDoub(l_num, n_bessel);

  // Initialize grid for Bessel function evaluation
  z_bessel = VecDoub(n_bessel);
  double besselmax  = 4.0*k_max; // NB: We need besselmax >= kmax*eta0. Check for this!
  for (int i = 0; i<n_bessel; i++) {
    z_bessel[i] = besselmax*double(i)/double(n_bessel-1);
  }

  // Open previous made file containing Bessel functions
  ifstream bessel_file;
  string bessel_filename = cosm_simu->outprefix + "precomputed_bessel.dat";
  bessel_file.open(bessel_filename.c_str());

  if (!bessel_file) {

#ifdef VERBOSE
    cout << "     Besselfunctions not found on disc" << endl;
    cout << "     Compute Besselfunctions spline" << endl;
#endif

    // Spherical besselfunction object
    Bessel Bess;

    // Evaluate Besselfunctions on a grid
    for (int l = 0; l<l_num; l++) {
      j_l(l,0) = 0.0;
      for (int i = 1; i<n_bessel; i++) {
        jlx = Bess.sphbesj(ls[l],z_bessel[i]);
        j_l(l,i) = abs(jlx) > 1.0e-10 ? jlx : 0.0; 
      }
    }

    // Spline Bessel functions over z
    bessel_spline = Spline(z_bessel, j_l, 1.0e30, 1.0e30, false);
    bessel_spline.set_linear_x();

    // Write to file
    ofstream save_bessel;
    save_bessel.open(bessel_filename.c_str());
    save_bessel << l_num << " " << n_bessel << " " << k_max << endl;
    for (int l = 0; l<l_num; l++) {
      for (int i = 0; i<n_bessel; i++) {
        save_bessel << l << "  " << i << " " << j_l(l,i) << endl;
      }
    }
    save_bessel.close();

#ifdef VERBOSE
    cout << endl;
    cout << "---> Output bessel functions" << endl;
    cout <<	"     Filename: '" << bessel_filename << "'" << endl;
#endif

  } else {

#ifdef VERBOSE
    cout << "     Besselfunction file found so will read from this" << endl;
    cout << "     NB: Assuming that l_num, ls[l], n_bessel and besselmax is the same as in this run" << endl;
    cout <<	"     Filename: '" << cosm_simu->outprefix << "precomputed_bessel.dat'" << endl;
#endif

    int l_num_file, n_bessel_file, k_max_file;
    bessel_file >> l_num_file;
    bessel_file >> n_bessel_file;
    bessel_file >> k_max_file;
    if(!(l_num_file == l_num && n_bessel_file == n_bessel && k_max_file == k_max)){
      cout << "Error: l_num, n_bessel and/or k_max in the precomputed besselfile do not agree with the ones we use" << endl;
      cout << "l_num    = " << l_num << " found " << l_num_file << endl;
      cout << "n_bessel = " << n_bessel << " found " << n_bessel_file << endl;
      cout << "k_max    = " << k_max << " found " << k_max_file << endl;
      exit(1);
    }

    while (!bessel_file.eof()){
      bessel_file >> ll; 
      bessel_file >> ii; 
      bessel_file >> j_l(ll,ii);
    }
    bessel_file.close();
    bessel_spline = Spline(z_bessel, j_l, 1.0e30, 1.0e30, false);
    bessel_spline.set_linear_x();
  }

#ifdef VERBOSE
  cosm_simu->time_now = clock() - cosm_simu->time_now;
  cout <<			"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#endif

}

double Powerspectrum::get_bessel(double x, int l){

  // The value of l
  int lvalue = ls[l];

  // Analytical values for l=2
  if (lvalue==2) {
    if (x > 0.1) {
      return sin(x)*(3.0/(x*x*x)-1.0/x)-3.0*cos(x)/(x*x);
    } else {
      return x*x/15.0;
    }
  }

  // If x smaller than sensitivity put j_l = 0.0 
  if (lvalue >= 20) {
    if (x < 0.2*lvalue) {
      return 0.0;
    }
  } else if (lvalue >= 100) {
    if (x < 0.7*lvalue) {
      return 0.0;
    }
  }
  return bessel_spline.get_spline_matrix(x, l);
}

// Calculate sigma8 = sigma(R=8Mpc/h)
void Powerspectrum::calc_sigma8(){
  double atol, rtol, h1;
  VecDoub yint = VecDoub(1);
  Output out;

  // Integration variables
  atol = 1.0e-10;
  rtol = 1.0e-10;
  h1   = 1.0e-8;
  yint[0] = 0.0;

  derivs_sigma8 dsigma8(this);

  Odeint<StepperBS<derivs_sigma8> > ode_sigma8(yint, k_min, k_max, atol, rtol, h1, 0.0, out, dsigma8);
  ode_sigma8.integrate();

  cosm_simu->sigma8 = sqrt(yint[0]);
}

