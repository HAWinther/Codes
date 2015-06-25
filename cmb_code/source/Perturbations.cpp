
//===========================================================
// General methods to integrate perturbations equations
//===========================================================

#include "Perturbations.h"

// Integrate the perturbations
void Perturbations::integrate_perturbations(){	
  double g, dg, ddg, tau, dtau, ddtau, H, dHH, ddHH, LPi, dLPi, ddLPi; 
  double Psi, dPmP, x_old, x_tight_coupling, ck_H, etafactor;

#ifdef VERBOSE
  cout << endl;
  cout << "**********************************************" << endl;
  cout <<	"***      Integrate the perturbations       ***" << endl;
  cout <<	"**********************************************" << endl;
  cosm_simu->time_now = clock();
#endif

  // Initialize k-array
  ks = VecDoub(n_k);
  for (int k = 0; k<n_k; k++) {
    ks[k] = k_min*pow(k_max/k_min, k/double(n_k-1));
  }

  // Initialize x-arrays
  x_rec = VecDoub(ntot);

  if (cosm_simu->z_reion > 0.0) {
    // Add extra points during reionization period if z_reion > 0
    for (int i = 0; i<ntot; i++) {
      if (i < n1) {
        x_rec[i] = x_start_rec + (x_end_rec-x_start_rec) * double(i)/double(n1);
      }else if(i < n1+ n2) {
        x_rec[i] = x_end_rec + (cosm_simu->x_start_reion - x_end_rec) * double(i-n1)/double(n2);
      }else if(i < n1 + n2 + n3) {
        x_rec[i] = cosm_simu->x_start_reion + (cosm_simu->x_end_reion - cosm_simu->x_start_reion) * double(i-n1-n2)/double(n3);
      }else {
        x_rec[i] = cosm_simu->x_end_reion + (0.0 - cosm_simu->x_end_reion) * double(i-n1-n2-n3)/double(n4-1);
      }
    }
  } else {
    for (int i = 0; i<ntot; i++) {
      if (i < n1+n2) {
        x_rec[i] = x_start_rec + (x_end_rec-x_start_rec) * double(i)/double(n1+n2);
      } else {
        x_rec[i] = x_end_rec + (0.0 - x_end_rec) * double(i-n1-n2)/double(ntot - 1 - n1 - n2);
      }
    }
  }

  // Initialize temperature source function
  S_low = MatDoub(n_k, ntot);

  // Initialize E-field source function
  SE_low = MatDoub(n_k, ntot);

  // Initialize gravitational potential
  Phi = MatDoub(n_k, ntot);

  // Initialize matter perturbations
  deltaM = MatDoub(n_k,ntot);

  // Initialize arrays that hold the integration variables
  VecDoub y_tight_coupling(ntight);
  VecDoub y(npar);
  VecDoub dydx(npar);

  // The ODE structures
  derivs_tight_coupling dtc(this);
  derivs_full dfull(this), get_derivs_full(this);

  // Integration variable
  Output out;

#ifdef OUTPUT_PERTURBATIONS

#ifdef VERBOSE
  cout << endl;
  cout << "---> Output perturbation quantites" << endl;
  cout <<	"     Filename: '" << cosm_simu->outprefix << "pert_QUANTITY.dat'" << endl;
#endif

  // Output
  ofstream out_phi, out_delta, out_v, out_theta, out_nu, out_thetap, out_source;
  string out_phi_file    = cosm_simu->outprefix + "perturbations_phi.dat";
  string out_delta_file  = cosm_simu->outprefix + "perturbations_delta.dat";
  string out_v_file      = cosm_simu->outprefix + "perturbations_v.dat";
  string out_theta_file  = cosm_simu->outprefix + "perturbations_theta.dat";
  string out_nu_file     = cosm_simu->outprefix + "perturbations_nu.dat";
  string out_thetap_file = cosm_simu->outprefix + "perturbations_thetap.dat";

  out_phi.open   (out_phi_file.c_str());
  out_delta.open (out_delta_file.c_str());
  out_v.open     (out_v_file.c_str());
  out_theta.open (out_theta_file.c_str());
  out_nu.open    (out_nu_file.c_str());
  out_thetap.open(out_thetap_file.c_str());

#endif

  //=====================================================================
  // Method in case we want to do store any do something with quantity(k,x)
  // We init here, store values while integrating and analyzing afterwards
  // Implement: init_store_values, store_values and analyze_store_values
  //=====================================================================
  init_store_values();

  // For debugging: Read low resoulution source from disc if needed
  if (cosm_simu->use_precomputed_pert){

    // Check that files exist
    ifstream infile(cosm_simu->precomputed_pert_sourcefile.c_str());
    if(!infile.good()) cosm_simu->use_precomputed_pert = false;;
    infile.close();
    infile.open(cosm_simu->precomputed_pert_sourcefile.c_str());
    if(!infile.good()) cosm_simu->use_precomputed_pert = false;
    infile.close();

    cout << endl;
    cout << "---> Skip integrating and read perturbations from file" << endl;

    if(cosm_simu->use_precomputed_pert){
      goto end_perturbations;
    } else {
      cout << "     Warning: use_precomputed_pert = true but no files found! Let us integrate then..."  << endl;
    }
  }

  // OpenMP paralellization: Make private y, y_tight_coupling, x_tight_coupling, k_current (put in threadprivate variable)
  // x_current, dydx, ck_H, LPi, dLPi, ddLPi, Psi, dPmP, etafactor

#ifdef OPENMP
#pragma omp parallel for private(y,y_tight_coupling,x_tight_coupling,dydx,ck_H,LPi,dLPi,ddLPi,Psi,dPmP,etafactor)
#endif
  // Loop over wavenumbers 'k' and calculate the perturbations
  for (int k = 0; k<n_k; k++) {

    // Set current wavenumber
    k_current = ks[k];

    // Set initial conditions
    set_initial_conditions_tight_coupling(x_init, y_tight_coupling);

    // Get time where tight coupling ends
    x_tight_coupling = get_tight_coupling_time(k_current);

    // Integrate the system until tight coupling is no longer valid
    if (k<ktight) {
      //Odeint<StepperBS<derivs_tight_coupling> >  ode_tc (y_tight_coupling,x_init,x_tight_coupling,a_bs, r_bs, h1_bs, 0.0,out,dtc);
      //ode_tc.integrate();

      int ttt = 200;
      for (int s = 0; s<ttt; s++) {
        double ddx = (x_tight_coupling - x_init)/double(ttt);
        double xx = x_init + (x_tight_coupling - x_init)*s/double(ttt);
        Odeint<StepperBS<derivs_tight_coupling> >  ode_tc (y_tight_coupling,xx,xx+ddx,a_bs, r_bs, h1_bs, 0.0,out,dtc);
        ode_tc.integrate();

        ////////////////////////////////////////////////////////////////////
        // Here we can store quantity(k,t) at times t < t_tight_coupling
        ////////////////////////////////////////////////////////////////////
        store_values(xx+ddx, k_current, &y_tight_coupling, true);

#ifdef OUTPUT_PERTURBATIONS
        // Output to file for 5 values of k
        if (int(100*k/double(n_k)) % 20 == 1) {
          x_current = xx+ddx;
          out_delta  << x_current << " " << y_tight_coupling[0]	<< " " << y_tight_coupling[1]	<< endl;
          out_v      << x_current << " " << y_tight_coupling[2]	<< " " << y_tight_coupling[3]	<< endl;
          out_phi    << x_current << " " << y_tight_coupling[4]	<< " " << 0.0	<< endl;
          out_theta  << x_current << " " << y_tight_coupling[5] << " " << y_tight_coupling[6]	<< " " << y_tight_coupling[7]	<< endl;
          out_thetap << x_current << " " << 0.0		<< " " << 0.0		<< " " << 0.0		<< endl;
          out_nu     << x_current << " " << y_tight_coupling[7]	<< " " << y_tight_coupling[8]	<< " " << 0.0	<< endl;
        }
#endif
      }
    } else {
      //Odeint<StepperSie<derivs_tight_coupling> > ode_tc (y_tight_coupling,x_init,x_tight_coupling,a_sie,r_sie,h1_sie,0.0,out,dtc);
      //ode_tc.integrate();

      int ttt = 200;
      for (int s = 0; s<ttt; s++) {
        double ddx = (x_tight_coupling - x_init)/double(ttt);
        double xx = x_init + (x_tight_coupling - x_init)*s/double(ttt);
        Odeint<StepperBS<derivs_tight_coupling> >  ode_tc (y_tight_coupling,xx,xx+ddx,a_bs, r_bs, h1_bs, 0.0,out,dtc);
        ode_tc.integrate();

        ////////////////////////////////////////////////////////////////////
        // Here we can store quantity(k,t) at times t < t_tight_coupling
        ////////////////////////////////////////////////////////////////////
        store_values(xx+ddx, k_current, &y_tight_coupling, true);

#ifdef OUTPUT_PERTURBATIONS
        // Output to file for 5 values of k
        if (int(100*k/double(n_k)) % 20 == 1) {
          x_current = xx+ddx;
          out_delta  << x_current << " " << y_tight_coupling[0]				<< " " << y_tight_coupling[1]	<< endl;
          out_v      << x_current << " " << y_tight_coupling[2]				<< " " << y_tight_coupling[3]	<< endl;
          out_phi    << x_current << " " << y_tight_coupling[4]				<< " " << 0.0	<< endl;
          out_theta  << x_current << " " << y_tight_coupling[5]				<< " " << y_tight_coupling[6]				<< " " << y_tight_coupling[7]				<< endl;
          out_thetap << x_current << " " << 0.0		<< " " << 0.0		<< " " << 0.0		<< endl;
          out_nu     << x_current << " " << y_tight_coupling[7]	<< " " << y_tight_coupling[8]	<< " " << 0.0	<< endl;
        }
#endif
      }
    }

    // Set initial conditions for the full system
    x_current = x_tight_coupling;
    set_initial_conditions_full(x_current, y_tight_coupling, y);

    // Integrate until recombination
    if (x_current<= x_start_rec) {
      Odeint<StepperSie<derivs_full> > ode_full_br(y,x_current,x_start_rec,a_sie,r_sie,h1_sie,0.0,out,dfull);
      ode_full_br.integrate();
    }

    // Initalize ODE objects
    Odeint<StepperBS<derivs_full > > ode_full (y, x_rec[0], x_rec[1], a_bs,  r_bs,  h1_bs,  0.0, out, dfull);
    Odeint<StepperSie<derivs_full> > ode1_full(y, x_rec[0], x_rec[1], a_sie, r_sie, h1_sie, 0.0, out, dfull);

    // Make arrays
    VecDoub yn(npar);
    VecDoub dfdxn(npar);
    MatDoub dfdyn(npar,npar);

    // Integrate from recombination and until today
    for (int i = 1; i<ntot; i++) {

      // Set current x
      x_current = x_rec[i];

      // Integrate system one step
      if (k<kfull) {
        ode_full.x1 = x_rec[i-1];
        ode_full.x2 = x_rec[i];
        ode_full.integrate();
      } else {
        ode1_full.x1 = x_rec[i-1];
        ode1_full.x2 = x_rec[i];
        ode1_full.integrate();
      }

      // Get derivatives
      get_derivs_full(x_current, y, dydx);

      // Calculate the source function (H = H(x)*exp(x))
      RECSUB->get_all_tau_g(x_current, tau, dtau, ddtau, g, dg, ddg);
      BGSUB->get_all_H_c(x_current, H, dHH, ddHH);

      // Modify: make calc_source function and implement this model by model!
      // set_source(int i, int k){

      ck_H  = k_current/H;
      LPi   =    y[7] +    y[6+lmax_int] +    y[8+lmax_int];
      dLPi  = dydx[7] + dydx[6+lmax_int] + dydx[8+lmax_int];
      ddLPi = 0.1*ck_H*(  4.0*(-dHH*y[6] + dydx[6]) - 6.0*(-dHH*(y[8]+y[7+lmax_int]+y[9+lmax_int]) + dydx[8]+dydx[7+lmax_int]+dydx[9+lmax_int])) + 0.3*(ddtau*LPi+dtau*dLPi);
      Psi   = -y[4] - 12.0*exp(-2.0*x_current)/(k_current*k_current)*(cosm_simu->Omega_r*y[7]+cosm_simu->Omega_nu*y[9+2*lmax_int]);
      dPmP  = -2.0*dydx[4]-12.0/(k_current*k_current)*exp(-2.0*x_current)*(cosm_simu->Omega_r*(-2.0*y[7] + dydx[7])+cosm_simu->Omega_nu*(-2.0*y[9+2*lmax_int] + dydx[9+2*lmax_int]));

      // Store matter perturbations (to test P(k) = delta_m^2)
      deltaM(k,i) = (y[0]*cosm_simu->Omega_m+y[1]*cosm_simu->Omega_b)/(cosm_simu->Omega_m+cosm_simu->Omega_b);

      // Store the T-source function
      S_low(k,i)  = g*(y[5] + Psi + 0.25*LPi) + exp(-tau)*dPmP - (1.0/ck_H)*((dHH*g+dg)*y[3] + g*dydx[3]) + (3.0/(4.0*ck_H*ck_H))*((dHH*dHH+ddHH)*g*LPi + 3.0*dHH*(dg*LPi+g*dLPi) + ddg*LPi + 2.0*dg*dLPi + g*ddLPi);

      // Store the E-source function
      etafactor = k_current*(BGSUB->get_eta(0.0)-BGSUB->get_eta(x_current));
      etafactor *= etafactor;
      SE_low(k,i) = (etafactor > 1.0e-20) ? g*0.75*LPi/etafactor : 0.0;

      // Store the gravitational potential
      Phi(k,i) = y[4];

      ////////////////////////////////////////////////////////////////////
      // Here we can store quantity(k,t) at times t_today > t > t_tight_coupling
      ////////////////////////////////////////////////////////////////////
      store_values(x_current, k_current, &y, false);

#ifdef OUTPUT_PERTURBATIONS
      // Output to file for 5 values of k
      if (int(100*k/double(n_k)) % 20 == 1) {
        out_delta  << x_current << " " << y[0]				<< " " << y[1]	<< endl;
        out_v      << x_current << " " << y[2]				<< " " << y[3]	<< endl;
        out_phi    << x_current << " " << y[4]				<< " " << Psi	<< endl;
        out_theta  << x_current << " " << y[5]				<< " " << y[6]				<< " " << y[7]				<< endl;
        out_thetap << x_current << " " << y[6+lmax_int]		<< " " << y[7+lmax_int]		<< " " << y[8+lmax_int]		<< endl;
        out_nu     << x_current << " " << y[7+2*lmax_int]	<< " " << y[8+2*lmax_int]	<< " " << y[8+2*lmax_int]	<< endl;
      }
#endif
    }

#ifdef VERBOSE
    if (k == 0) {
      fprintf(stdout,"\n---> Progress meter:\n     ");
    }

    if( (((100*(k+1))/n_k % 10 == 0 ) && (100*k)/n_k % 10 != 0)) 
      fprintf(stdout,"%i%% ",(100*(k+1))/n_k);

    if (k ==  n_k - 1) {
      fprintf(stdout,"\n");
    }

    fflush(stdout);
#endif

    }

    // General method to do something with data stored using store_perturbations_values()
    analyze_store_values();

    // Fill in the last values by hand for smooth spline at endpoints
    for (int k = 0; k<n_k; k++) {
      S_low(k,0)	= S_low(k,1);
      SE_low(k,0) = SE_low(k,1);
      Phi(k,0)	= Phi(k,1);
      deltaM(k,0) = deltaM(k,1);
    }

    // Here we arrive if use_precomputed_pert = true
end_perturbations:

    if (!cosm_simu->use_precomputed_pert) {

#ifdef VERBOSE
      cout << endl;
      cout << "---> Output S_low_TT(k,x) and S_low_E(k,x)" << endl;
      cout <<	"     Filename: '" << cosm_simu->outprefix << "precomputed_source.dat'"<< endl;
#endif	

      ofstream outsource;
      string outsource_filename = cosm_simu->outprefix + "precomputed_source.dat";
      outsource.open(outsource_filename.c_str());
      outsource << ntot << " " << n_k << endl;
      for (int i = 0; i<ntot; i++) {
        for (int k = 0; k<n_k; k++) {
          outsource << i << " " << k << " " << S_low(k,i) << " " << SE_low(k,i) << endl;
        }
      }
      outsource.close();

#ifdef VERBOSE
      cout << endl << "---> Output Phi(k,x)"<< endl;
      cout <<	"     Filename: '" << cosm_simu->outprefix << "precomputed_phi.dat'"<< endl;
#endif

      ofstream outphi;
      string outphi_filename = cosm_simu->outprefix + "precomputed_phi.dat";
      outphi.open(outphi_filename.c_str());
      outphi << ntot << " " << n_k << endl;
      for (int i = 0; i<ntot; i++) {
        for (int k = 0; k<n_k; k++) {
          outphi << i << " " << k << " " << Phi(k,i) << endl;
        }
      }
      outphi.close();
    } else {

#ifdef VERBOSE
      cout << endl;
      cout << "---> Reading precomputed S_low_TT(k,x) and S_low_E(k,x)"<< endl;
      cout <<	"     Filename: '" << cosm_simu->precomputed_pert_sourcefile << "'" << endl;
#endif

      int iin, kin;
      double sin;

      ifstream readsource;
      readsource.open(cosm_simu->precomputed_pert_sourcefile.c_str());
      readsource >> iin;
      readsource >> kin;
      if(!(iin == ntot && kin == n_k)){
        cout << "Error: ntot and n_k in sourcefile does not match the ones we use now" << endl;
        cout << "ntot = " << ntot << " found " << iin << endl;
        cout << "n_k  = " << n_k  << " found " << kin << endl;
        exit(1);
      }
      while (readsource) {
        readsource >> iin;
        readsource >> kin;
        readsource >> sin;
        S_low(kin,iin) = sin;
        readsource >> sin;
        SE_low(kin,iin) = sin;
      }
      readsource.close();

#ifdef VERBOSE
      cout << endl;
      cout << "---> Reading precomputed Phi(k,x)"<< endl;
      cout <<	"     Filename: '" << cosm_simu->precomputed_pert_phifile << "'" << endl;
#endif
      ifstream readphi;
      readphi.open(cosm_simu->precomputed_pert_phifile.c_str());
      readphi >> iin;
      readphi >> kin;
      if(!(iin == ntot && kin == n_k)){
        cout << "Error: ntot and n_k in phifile does not match the ones we use now" << endl;
        cout << "ntot = " << ntot << " found " << iin << endl;
        cout << "n_k  = " << n_k  << " found " << kin << endl;
        exit(1);
      }
      while (readphi) {
        readphi >> iin;
        readphi >> kin;
        readphi >> Phi(kin,iin);
      }
      readphi.close();
    }

#ifdef OUTPUT_PERTURBATIONS
    out_delta.close();
    out_v.close();
    out_phi.close();
    out_thetap.close();
    out_theta.close();
    out_nu.close();
#endif

#ifdef VERBOSE
    cosm_simu->time_now = clock() - cosm_simu->time_now;
    cout << endl;
    cout << "---> Time consumption" << endl;
    cout <<	"     Time used (sec)            : " << cosm_simu->time_now*1.0e-6 << endl;
#endif
  }

  // Calculate the time when tight-coupling approximation ends
  double Perturbations::get_tight_coupling_time(double k){
    int    n = 500;
    double x_start_search = -10.0;
    double x = x_start_search;
    double dtau;

    // Tight coupling is valid as long as dtau < 10 and ck/(H*dtau) > 0.1 and x < x_start_rec
    while (true) {
      x += (x_start_rec-x_start_search)/double(n);
      dtau = abs(RECSUB->get_dtau(x));
      if (x >= x_start_rec) {
        return x_start_rec;
      }else if(dtau < 10.0 || 10.0*k/(BGSUB->H_c(x)*dtau) > 1.0){
        return x;
      }
    }
  }

  //===========================================================
  // General methods to store and analyze perturbation 
  // quantities Q(x,t). Init. Store while integrating
  // and analyze afterwards
  //===========================================================

  void Perturbations::init_store_values(){


  }

  void Perturbations::store_values(double x, double k, VecDoub *y, bool tight_coupling){
    // Store quantity(k,x)...

    if(tight_coupling){

    } else {

    }

  }

  void Perturbations::analyze_store_values(){

    // Do something with the data stored...

    // Deallocate array...

  }


