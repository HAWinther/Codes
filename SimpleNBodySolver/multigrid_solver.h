#ifndef MGSOLVER_INC
#define MGSOLVER_INC
#include <math.h>

//================================================
//
// General multigrid solver
// 
// All model-specific implementation inside
// -DMODEL define tags
//
//================================================

//================================================
// The cosmological value of the field we are
// solving for
//================================================
double phi_backgroud(double aexp){
  double phi;

#if defined(FOFR)
  // This if phi in |f_R| = a^{-2} e^{phi}
  phi = (sim._nfofr+1.0) * log( (1.0 + 4.0*sim._omega_lambda/sim._omega_m) * pow3(sim._aexp) * 
        pow(sim._fofr0*sim._aexp,1.0/(sim._nfofr+1.0)) / (1.0 + 4.0*pow3(sim._aexp)*sim._omega_lambda/sim._omega_m) );
#elif defined(SYMMETRON)
  phi = aexp > sim._assb_symm ? sqrt(1.0 - pow3(sim._assb_symm / aexp)) : 0.001;
#elif defined(GR)
  phi = 0.0;
#else
  .................
  ...Not defined...
  .................
#endif

  return phi;
}

//================================================
// The initial guess for the solver
//================================================
void make_initial_guess(double *phi, double aexp){
  int Ntot = sim._mg->Ntot;
  double phi_guess = phi_backgroud(aexp);

  std::cout << "Taking phi = " << phi_guess << " as starting guess for solver" << std::endl;

#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0;i < Ntot; i++){
    phi[i] = phi_guess;
  }
}

//================================================
// The L operator (EOM is written as L(phi) = 0 )
//================================================
double l_operator(int level, int i, double kinetic, bool makesource){
  double rho, h, source, phi;
  double fac1, fac2;

  // Gridspacing
  h = 1.0/double(sim._mg->N) * (1 << level);

  // Density and scalar field in node
  rho = sim._mg->rho[level][i];
  phi = sim._mg->f[level][i];

#if defined(FOFR)
  // Equation of motion for f(R) gravity L = kinetic - source = D[ e^phi Dphi ] - source(phi,rho,a)
  fac1   = 1.0 + 4.0*pow3(sim._aexp)*sim._omega_lambda/sim._omega_m;
  fac2   = (1.0 + 4.0*sim._omega_lambda/sim._omega_m)*pow3(sim._aexp)*pow(sim._fofr0*sim._aexp,1.0/(sim._nfofr+1.0));
  source = sim._omega_m * sim._aexp * pow2(sim._BoxSize/2998.0) * (rho - 1.0 + fac1 - fac2 * exp(-phi/(sim._nfofr+1.0)));
#elif defined(SYMMETRON)
  source = 0.5*pow2(sim._aexp * sim._BoxSize / sim._lambda_symm) * ( pow3(sim._assb_symm/sim._aexp) * rho * phi  + pow3(phi)- phi );
#elif defined(GR)
  // The Poisson equation for the gravitational potential
  source = 1.5 * sim._omega_m * sim._aexp * (rho - 1.0);
#else
  .................
  ...Not defined...
  .................
#endif

  // Add the source term arising from restricting the equation down to the lower level
  if( (level > 0) && !makesource ){
    source += sim._mg->source[level][i];
  }

  return (kinetic/(h*h) - source);
}


//================================================
// The derivative of the L-operator
//================================================
double dl_operator(int level, int i, double dkinetic){
  double dl, h, dsource, phi, rho, fac2;

  // Gridspacing
  h = 1.0/double(sim._mg->N) * (1 << level);

  // Density and scalar field in node
  rho = sim._mg->rho[level][i];
  phi = sim._mg->f[level][i];

#if defined(FOFR)
  // dL/dphi for Hu-Sawicky f(R) gravity
  fac2    = (1.0 + 4.0*sim._omega_lambda/sim._omega_m)*pow3(sim._aexp)*pow(sim._fofr0*sim._aexp,1.0/(sim._nfofr+1.0));
  dsource = sim._omega_m * sim._aexp * pow2(sim._BoxSize/2998.0) * ( fac2 * exp(-phi/(sim._nfofr+1.0)))/(1.0+sim._nfofr);
#elif defined(SYMMETRON)
  dsource = 0.5*pow2(sim._aexp * sim._BoxSize / sim._lambda_symm) * ( pow3(sim._assb_symm /sim._aexp) * rho + 3.0*phi*phi - 1.0);
#elif defined(GR)
  // Poisson equation
  dsource = 0.0;
#else
  .................
  ...Not defined...
  .................
#endif

  dl = dkinetic/(h*h) - dsource;
  if(fabs(dl) < 1e-10){
    std::cout << "Error: dl close to 0" << std::endl;
    exit(1);
  }

  return dl;
}

//================================================
// Calculate fifth-force F_fifth which is the term
// in the geodesic equation:
// x'' = -DPhi - F_fifth
//================================================

void calculate_fifth_force(){
  double *fx = sim._FifthForce[0];
  double *fy = sim._FifthForce[1];
  double *fz = sim._FifthForce[2];

  // First calculate the gradient of phi
  for(int k = 0; k < 3; k++)
    calculate_gradient(sim._ScalarField, sim._FifthForce[k], k, sim._ngrid);

#ifdef OMP
#pragma omp parallel for
#endif
  // Translate to force in geodesic equation
  for(int i = 0; i < sim._ngrid_tot; i++){

#if defined(FOFR)
    double norm = 0.5/pow2(sim._BoxSize/2998.0) * exp(sim._ScalarField[i]);
#elif defined(SYMMETRON)
    double norm = 6.0*sim._omega_m*pow2(sim._beta_symm * sim._lambda_symm * sim._aexp)/pow3(sim._assb_symm) * sim._ScalarField[i];
#elif defined(GR)
    double norm = 1.0;
#else
  .................
  ...Not defined...
  .................
#endif
    fx[i] *= norm;
    fy[i] *= norm;
    fz[i] *= norm;
  }
}

//================================================
// Solve phi(x,y,z) using multigrid
//================================================
void solve_phi_with_multigrid(){
  int maxsteps = 1000;
  int vcycles;

  std::cout << "======> Solve for scalar field" << std::endl;
  
  // Initialize arrays
  for(int i = 0; i < sim._ngrid_tot; i++){
    sim._MGSource[i] = sim._MGResidual[i] = 0.0;
  }

#if defined(SYMMETRON)
  // Do not solve if a < assb as the field is zero here anyway
  if(sim._aexp < sim._assb_symm){
    std::cout << "       aexp < assb so no need to solve for symmetron" << std::endl;
    return;
  }
#endif

  // Make initial guess (only needed first step as we can use old solution as the initial guess after that)
  // make_initial_guess(sim._ScalarField, sim._aexp);

  // Set density on all levels
  restrict_down_density_all();

  // Pre-solve on domaingrid
  solve_current_level(0);

  // Define the initial residual before the multigrid
  sim._mg->rms_res_i = sim._mg->rms_res;
  sim._mg->rms_res_old = 0.0;

  // Check if already have convergence
  if(check_convergence(0)) return;

  // The V-cycle
  std::cout << "======> Initial residual = " << sim._mg->rms_res_i << std::endl;
  vcycles = 1;
  while (1) {

    // Go down to the bottom
    recursive_go_down(0);

    // Go up to the top
    recursive_go_up(sim._mg->Nlevel-2);

    // Check for convergence
    if(check_convergence(vcycles)) break;

    std::cout << "======> V-cycle i = " << vcycles << " res = " << sim._mg->rms_res << std::endl;
    calc_min_max_scalarfield();

    // Break after NSTEP to avoid infinite loop...
    if(vcycles >= maxsteps){
      std::cout << "WARNING: Multigrid failed to converge..." << std::endl;
      std::cout << " res = " << sim._mg->rms_res << " res_old = " << sim._mg->rms_res_old;
      std::cout << " res_i = " << sim._mg->rms_res_i << std::endl;
      break;
    }
    vcycles++;
  }
}

//================================================
// Calculate max/min/mean values [ we compute f_R = e^{phi}/a^2 ]
//================================================
void calc_min_max_scalarfield(){
  double phi_avg, phi_max = -1e30, phi_min = 1e30, phi_back, phinow;

  // Background field-value
#if defined(FOFR)
  phi_back = exp(phi_backgroud(sim._aexp))/pow2(sim._aexp);
#elif defined(SYMMETRON)
  phi_back = phi_backgroud(sim._aexp);
#elif defined(GR)
  phi_back = 1.0;
#else
  .................
  ...Not defined...
  .................
#endif

  phi_avg = 0.0;
  for(int i = 0;i < sim._ngrid_tot; i++){
#if defined(FOFR)
    // Translate from phi_code to |f_R|
    phinow = exp(sim._ScalarField[i])/pow2(sim._aexp);
#elif defined(SYMMETRON) || defined(GR)
    phinow = sim._ScalarField[i];
#else
  .................
  ...Not defined...
  .................
#endif
    if(phinow > phi_max) phi_max = phinow;
    if(phinow < phi_min) phi_min = phinow;
    phi_avg += phinow;
  }
  phi_avg /= double(sim._ngrid_tot);
  cout << "        Scalar field Min/Avg/Max : " << phi_min << " " << phi_avg << " " << phi_max << "( phi_b = " << phi_back << " )" << endl;
}

//================================================
// GS change for node
//================================================
double calc_gs_change(int level, int i, double kinetic, double dkinetic){
  double l, dl;
  l  = l_operator(level, i, kinetic, false);
  dl = dl_operator(level, i, dkinetic);
  return (l/dl);
}

//================================================
// Calculates the residual of the operator equation for one node
// and stores minus the residual
//================================================
double residual(int level, int i, double kinetic){
  return (-l_operator(level, i, kinetic, false));
}

//================================================
// Calculates the residual on a grid
//================================================
double calculate_residual(int level){
  int Ntot,N,N2,ixP,ixM,iyP,iyM,izP,izM,ix,iy,iz;
  double kinetic, totres, b0, bp, bm;
  double *res, *ff;

  N    = sim._mg->N / (1 << level);
  N2   = N*N;
  Ntot = N*N*N;

  // Pointers to the grid variables
  ff  = sim._mg->f[level];
  res = sim._mg->res[level];

  // Calculate residual in each node
#ifdef OMP
#pragma omp parallel for private(ix,iy,iz,kinetic,ixP,ixM,iyP,iyM,izP,izM,b0,bp,bm)
#endif
  for (int i = 0; i<Ntot; i++) {
    ix = i % N;
    iy = (i/N) % N;
    iz = i/N2;

    ixP = (ix == N-1)  ? 1-N : +1;
    ixM = (ix == 0)    ? N-1 : -1;
    iyP = ((iy == N-1) ? 1-N : +1)*N;
    iyM = ((iy == 0)   ? N-1 : -1)*N;
    izP = ((iz == N-1) ? 1-N : +1)*N2;
    izM = ((iz == 0)   ? N-1 : -1)*N2;

#if defined(FOFR)
    // Non-canonical kinetic term D[ b Dphi ] with b = e^phi
    b0 = exp(ff[i]);
    bp = 0.5 * (exp(ff[i + ixP]) + b0);
    bm = 0.5 * (exp(ff[i + ixM]) + b0);
    kinetic  = bp *(ff[i + ixP]-ff[i]) + bm*(ff[i + ixM]-ff[i]);
    bp = 0.5 * (exp(ff[i + iyP]) + b0);
    bm = 0.5 * (exp(ff[i + iyM]) + b0);
    kinetic  += bp*(ff[i + iyP]-ff[i]) + bm*(ff[i + iyM]-ff[i]);
    bp = 0.5 * (exp(ff[i + izP]) + b0);
    bm = 0.5 * (exp(ff[i + izM]) + b0);
    kinetic  += bp*(ff[i + izP]-ff[i]) + bm*(ff[i + izM]-ff[i]);
#elif defined(SYMMETRON) || defined(GR)
    // Canonical kinetic term D^2 phi
    kinetic = 0.0;
    kinetic += ff[i + ixP] + ff[i + ixM] - 2.0*ff[i];
    kinetic += ff[i + iyP] + ff[i + iyM] - 2.0*ff[i];
    kinetic += ff[i + izP] + ff[i + izM] - 2.0*ff[i];
#else
    .................
    ...Not defined...
    .................
#endif

    res[i] = residual(level,i,kinetic);
  }

  // Calculate RMS residual
  totres = 0.0;
#ifdef OMP
#pragma omp parallel for reduction(+:totres)
#endif
  for (int i = 0; i < Ntot; i++)
    totres += res[i]*res[i];
  totres = sqrt(totres/double(Ntot));
  return totres;
}

//================================================
// Criterion for convergence here
// Standard ways are: based on residual or
// the ration of the residual to the initial
// residual (err)
//================================================
bool check_convergence(int istep){
  double err, rel_diff, res, res_i, res_o;
  res   = sim._mg->rms_res;
  res_i = sim._mg->rms_res_i;
  res_o = sim._mg->rms_res_old;
  err   = res_i != 0.0 ? res/res_i : 1.0;
  sim._mg->rms_res_old = sim._mg->rms_res;

  // Print out the residual and error information
  if(sim._mg->verbose){
    std::cout << "Checking for convergence... step = " << istep << std::endl;
    std::cout << "Residual = " << res << "  Residual_old = " <<  res_o << "  Residual_i = " << res_i << std::endl;
    std::cout << "Error = " <<  err << std::endl;
  }

  // Another option is [err < ...]
  if(res < sim._EPS_CONVERGE){
    std::cout << "======> The solution has converged (res=" << res << ") (err=" << err << ") istep = " << istep << std::endl;
    return true;
  }

  return false;
}

//================================================
// Prolonge up solution phi
//================================================
void prolonge_up(int to_level, double* Old, double* New){
  int ix, iy, iz, N, N2, Ntot, Nold, Nold2, ix_old, iy_old, iz_old;
  int dx_P, dy_P, dz_P, i_old;
  double jfac, kfac, ifac, tfac;

  if(sim._mg->verbose)
    std::cout << "Prolonge dphi up to level = " << to_level << std::endl;

  // Number of grid nodes on top and bottom
  N     = sim._mg->N / (1 << to_level);
  N2    = N*N;
  Ntot  = N*N*N;
  Nold  = N/2;
  Nold2 = Nold*Nold;

  // Trilinear prolongation
#ifdef OMP
#pragma omp parallel for private(ix,iy,iz,ix_old,iy_old,iz_old,dx_P,dy_P,dz_P,ifac,jfac,kfac,tfac,i_old)
#endif
  for (int i = 0; i<Ntot; i++) {
    // Gather indices in new and old array and distances
    ix     = i % N;
    iy     = i/N % N;
    iz     = i/N2;

    ix_old = ix/2;
    iy_old = iy/2;
    iz_old = iz/2;

    dx_P   = ((ix_old == Nold-1) ? 1-Nold :  1);
    dy_P   = ((iy_old == Nold-1) ? 1-Nold :  1)*Nold;
    dz_P   = ((iz_old == Nold-1) ? 1-Nold :  1)*Nold2;

    ifac   = ix%2==0 ? 0.0 : 1.0;
    jfac   = iy%2==0 ? 0.0 : 1.0;
    kfac   = iz%2==0 ? 0.0 : 1.0;

    tfac = 1.0/double((ifac+1)*(jfac+1)*(kfac+1));
    i_old = ix_old + Nold*iy_old + Nold2*iz_old;

    // Prolongate
    New[i] = tfac* ( Old[i_old] +
        ifac *  Old[i_old + dx_P] +
        jfac *  Old[i_old + dy_P] +
        kfac *  Old[i_old + dz_P] +
        ifac * jfac * Old[i_old + dx_P + dy_P] +
        jfac * kfac * Old[i_old + dy_P + dz_P] +
        kfac * ifac * Old[i_old + dz_P + dx_P] +
        ifac * jfac * kfac * Old[i_old + dx_P + dy_P + dz_P]);
  }
}

//================================================
// The Gauss-Seidel Sweeps
//================================================
void Gauss_Seidel(int level, int redblack){
  double *ff, kinetic, dkinetic, bp, bm, b0;
  int ixP,ixM,iyP,iyM,izP,izM,ix,iy,iz, Ntot,N,N2;

  ff   = sim._mg->f[level];
  N    = sim._mg->N / (1 << level);
  N2   = N*N;
  Ntot = N*N*N;

  // Loop over all cells
#ifdef OMP
#pragma omp parallel for private(ix,iy,iz,ixP,ixM,iyP,iyM,izP,izM,kinetic,dkinetic,b0,bp,bm)
#endif
  for (int i=0; i<Ntot; i++) {
    ix = i % N;
    iy = i/N % N;
    iz = i/N2;

    // The color function here [redblack is chessboard ordering]
    if( ((ix+iy+iz) % 2) == redblack){
      ixP = ((ix == N-1) ? 1-N : +1);
      ixM = ((ix == 0)   ? N-1 : -1);
      iyP = ((iy == N-1) ? 1-N : +1)*N;
      iyM = ((iy == 0)   ? N-1 : -1)*N;
      izP = ((iz == N-1) ? 1-N : +1)*N2;
      izM = ((iz == 0)   ? N-1 : -1)*N2;

#if defined(FOFR)
      // Calculate non-canonical kinetic term  D[ b Dphi ] with b = e^phi
      b0 = exp(ff[i]);
      bp = 0.5 * (exp(ff[i + ixP]) + b0);
      bm = 0.5 * (exp(ff[i + ixM]) + b0);
      kinetic   = bp* (ff[i + ixP]-ff[i]) + bm*(ff[i + ixM]-ff[i]);
      dkinetic  = -6.0*ff[i]*b0;
      dkinetic += -(bp+bm) + (ff[i + ixP] + ff[i + ixM])*b0;
      bp = 0.5 * (exp(ff[i + iyP]) + b0);
      bm = 0.5 * (exp(ff[i + iyM]) + b0);
      kinetic  += bp*(ff[i + iyP]-ff[i]) + bm*(ff[i + iyM]-ff[i]);
      dkinetic += -(bp+bm) + (ff[i + iyP] + ff[i + iyM])*b0;
      bp = 0.5 * (exp(ff[i + izP]) + b0);
      bm = 0.5 * (exp(ff[i + izM]) + b0);
      kinetic  += bp*(ff[i + izP]-ff[i]) + bm*(ff[i + izM]-ff[i]);
      dkinetic += -(bp+bm) + (ff[i + izP] + ff[i + izM])*b0;
#elif defined(SYMMETRON) || defined(GR)
      // Canonical kinetic term D^2 phi and d/dphi [D^2 phi]
      kinetic  = 0.0;
      kinetic += ff[i + ixP] + ff[i + ixM] - 2.0*ff[i];
      kinetic += ff[i + iyP] + ff[i + iyM] - 2.0*ff[i];
      kinetic += ff[i + izP] + ff[i + izM] - 2.0*ff[i];
      dkinetic = -6.0; 
#else
      .................
      ...Not defined...
      .................
#endif
      // Update the solution phi = phi - L/dL
      ff[i] -= calc_gs_change(level, i, kinetic, dkinetic);
    }
  }
}

//================================================
// Solve on the current level
//================================================
void solve_current_level(int level){
  int ngs_sweeps;

  // Number of sweeps we do
  if(level == 0)
    ngs_sweeps = sim._NGS_FINE;
  else
    ngs_sweeps = sim._NGS_COARSE;

  // Do N Gauss-Seidel Sweeps
  for (int i=0; i<ngs_sweeps; i++) {
   
    Gauss_Seidel(level,0); // Red sweep
    Gauss_Seidel(level,1); // Black sweep

    // For debug. Gives residual every 3 sweeps for all levels
    if(i % 3 == 0 && sim._mg->verbose){
      cout << "        level = " << setw(5) << level << " NGS Sweep = " << setw(5) << i;
      cout << " Residual = " << setw(10) << calculate_residual(level) << endl;
    }
  }

  // Store domaingrid residual
  if (level == 0){
    sim._mg->rms_res = calculate_residual(0);
  }
}

//================================================
// Correct phi on level = level
//================================================
void correct_solution(int level){
  int Ntot;
  double *f, *df, max_corr = -1e30;
  f    = sim._mg->f[level];
  df   = sim._mg->res[level];
  Ntot = pow3(sim._mg->N / (1 << level));
#ifdef OMP
#pragma omp parallel for
#endif
  for (int i=0; i<Ntot; i++){
    f[i] += df[i];
    if(fabs(df[i]) > max_corr) max_corr = fabs(df[i]);
  }
  //cout << "        Max corr: " << max_corr << endl;
}

//================================================
// V-cycle go all the way up
//================================================
void recursive_go_up(int to_level){
  int from_level = to_level + 1;

  // Restrict down R[f] and store in res
  restrict_down_single(to_level, sim._mg->f[to_level], sim._mg->res[from_level]);

  // Make prolongation array ready at from-level
  make_prolongation_array(from_level);

  // Prolonge up solution from-level to to-level (stored in res)
  prolonge_up(to_level, sim._mg->res[to_level+1], sim._mg->res[to_level]);

  // Correct solution at to-level
  correct_solution(to_level);

  // Calculate new residual
  calculate_residual(to_level);

  // Solve on current level
  solve_current_level(to_level);

  // Go up again
  if(to_level > 0)
    recursive_go_up(to_level-1);
  else {
    return;
  }
}

//================================================
// Make the array we are going to prolonge up
//================================================
void make_prolongation_array(int level){
  int Ntot;
  Ntot = pow3(sim._mg->N / (1 << level));

#ifdef OMP
#pragma omp parallel for
#endif
  for (int i = 0; i < Ntot; i++)
    sim._mg->res[level][i] = (sim._mg->f[level][i] - sim._mg->res[level][i]);
}

//================================================
// Make new source and store in res
//================================================
void make_new_source(int level){
  int ix, iy, iz, N, N2, Ntot, ixP, ixM, iyP, iyM, izP, izM;
  double kinetic, res, *ff, b0,bp,bm;

  ff   = sim._mg->f[level];
  N    = sim._mg->N / (1 << level);
  N2   = N*N;
  Ntot = N*N*N;

  // Calculate the new source
#ifdef OMP
#pragma omp parallel for private(ix,iy,iz,ixP,ixM,iyP,iyM,izP,izM,kinetic,res,b0,bp,bm)
#endif
  for(int i = 0; i < Ntot; i++){
    ix = i % N;
    iy = i/N % N;
    iz = i/N2;

    ixP = ((ix == N-1) ? 1-N : +1);
    ixM = ((ix == 0)   ? N-1 : -1);
    iyP = ((iy == N-1) ? 1-N : +1)*N;
    iyM = ((iy == 0)   ? N-1 : -1)*N;
    izP = ((iz == N-1) ? 1-N : +1)*N2;
    izM = ((iz == 0)   ? N-1 : -1)*N2;

#if defined(FOFR)
    // Calculate non-canonical kinetic term: D[ b D phi ]
    b0 = exp(ff[i]);
    bp = 0.5 * (exp(ff[i + ixP]) + b0);
    bm = 0.5 * (exp(ff[i + ixM]) + b0);
    kinetic  = bp *(ff[i + ixP]-ff[i]) + bm*(ff[i + ixM]-ff[i]);
    bp = 0.5 * (exp(ff[i + iyP]) + b0);
    bm = 0.5 * (exp(ff[i + iyM]) + b0);
    kinetic += bp *(ff[i + iyP]-ff[i]) + bm*(ff[i + iyM]-ff[i]);
    bp = 0.5 * (exp(ff[i + izP]) + b0);
    bm = 0.5 * (exp(ff[i + izM]) + b0);
    kinetic += bp *(ff[i + izP]-ff[i]) + bm*(ff[i + izM]-ff[i]);
#elif defined(SYMMETRON) || defined(GR)
    // The canonical kinetic term D^2 phi
    kinetic = 0.0;
    kinetic += ff[i + ixP] + ff[i + ixM] - 2.0*ff[i];
    kinetic += ff[i + iyP] + ff[i + iyM] - 2.0*ff[i];
    kinetic += ff[i + izP] + ff[i + izM] - 2.0*ff[i];
#else
    .................
    ...Not defined...
    .................
#endif

    res = l_operator(level, i, kinetic, true);

    // Calculate the new source
    sim._mg->source[level][i] = sim._mg->res[level][i] + res;
  }
}

//================================================
// V-cycle go all the way down
//================================================
void recursive_go_down(int from_level){
  int to_level = from_level+1;

  // Check if we are at the bottom
  if(from_level+1 >= sim._mg->Nlevel) return;

  // Restrict residual
  restrict_down_single(from_level, sim._mg->res[from_level], sim._mg->res[to_level]);

  // Restrict solution
  restrict_down_single(from_level, sim._mg->f[from_level], sim._mg->f[to_level]);

  // Make new source and store in res
  make_new_source(to_level);

  // Solve on current level
  solve_current_level(from_level+1);

  // Recursive call
  recursive_go_down(from_level+1);
}

//================================================
// Restrict down grid Top and store in Bottom
//================================================

void restrict_down_single(int from_level, double *Top, double *Bottom){
  int N, N2, Ntot, N_bottom, N_bottom2, N_bottom_tot;
  int ix, iy, iz, i_bottom;
  double fac;

  // Nodes on top and bottom level
  N         = sim._mg->N / (1 << from_level);
  N2        = N*N;
  Ntot      = N*N*N;
  N_bottom  = N/2;
  N_bottom2 = N_bottom*N_bottom;
  N_bottom_tot = pow3(N_bottom);

  if(sim._mg->verbose) 
    cout << "Restrict down single: N = " << N << " Nb = " << N_bottom << " from_level = " << from_level << endl; 

  // Clear bottom array
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < N_bottom_tot; i++)
    Bottom[i] = 0.0;

  // Loop over top grid
  fac = 1.0/double(1 << 3);
  for (int i = 0; i < Ntot; i++) {
    ix = i % N;
    iy = i/N % N;
    iz = i/N2;

    // Bottom array index
    i_bottom = (ix/2) + (iy/2)*N_bottom + (iz/2)*N_bottom2;

    // Restrict down field (Note: OMP CRITICAL slowed it down)
    Bottom[i_bottom] += Top[i]*fac;
  }
}

//================================================
// Restricts down all grids (used for density only)
//================================================
void restrict_down_density_all(){
  for(int i = 1; i < sim._mg->Nlevel; i++){

    // Assign particles directly to grid / or use analytical formula instead of restriction
    if(sim._analytic_matter_density)
      set_analytical_density_field(sim._mg->rho[i],sim._ngrid / (1 << i));
    else
      assign_particles_to_grid(sim._mg->rho[i],sim._ngrid / (1 << i));
    
    // The alternative: to restrict it down can lead to solver braking down
    // restrict_down_single(i, sim._mg->rho[i-1], sim._mg->rho[i]);
  }
}

//================================================
// Calculate gradient of field
//================================================

void calculate_gradient(double *phi, double *gradient, int dim, int ngrid){
  int Ntot, N, N2, ix,iy,iz;
  int dx_P,dx_PP,dy_P,dy_PP,dz_P,dz_PP,dx_M,dx_MM,dy_M,dy_MM,dz_M,dz_MM;
  double w[4], dx;

  N      = ngrid;
  N2     = N*N;
  Ntot   = N*N*N;

  // Weights for 5-point stencil of the derivative
  dx   = 1.0/double(N);
  w[0] = 8.0/12.0/dx;
  w[1] = -1.0/12.0/dx;
  w[2] = -8.0/12.0/dx;
  w[3] = 1.0/12.0/dx;

  // Loop over all cells and calculate the derivative
#ifdef OMP
#pragma omp parallel for private(ix,iy,iz,dx_P,dx_PP,dy_P,dy_PP,dz_P,dz_PP,dx_M,dx_MM,dy_M,dy_MM,dz_M,dz_MM)
#endif
  for (int i = 0; i < Ntot; i++) {

    // Gradient in x direction
    if(dim == 0){
      ix = i % N;
      dx_P = 1; dx_PP = 2; dx_MM = -2; dx_M = -1;
      if (ix < 2 || ix > N-3) {
        if(ix==0)	  {dx_MM = N-2; dx_M  = N-1;}
        if(ix==1)	  {dx_MM = N-2;}
        if(ix==N-2)	{dx_PP = 2-N;}
        if(ix==N-1)	{dx_P  = 1-N; dx_PP = 2-N;}
      }
      gradient[i] = w[0]*phi[i + dx_P] + w[1]*phi[i + dx_PP] + w[2]*phi[i + dx_M] + w[3]*phi[i + dx_MM];

    // Gradient in y direction
    } else if(dim ==1){
      iy = i/N % N;
      dy_P = 1; dy_PP = 2; dy_MM = -2; dy_M = -1;
      if (iy < 2 || iy > N-3) {
        if(iy==0)	  {dy_MM = N-2; dy_M  = N-1;}
        if(iy==1)	  {dy_MM = N-2;}
        if(iy==N-2)	{dy_PP = 2-N;}
        if(iy==N-1)	{dy_P  = 1-N; dy_PP = 2-N;}
      }
      dy_M *= N; dy_P *= N; dy_MM *= N; dy_PP *= N;
      gradient[i] = w[0]*phi[i + dy_P] + w[1]*phi[i + dy_PP] + w[2]*phi[i + dy_M] + w[3]*phi[i + dy_MM];

    // Gradient in z direction
    } else {
      iz = i/N2;
      dz_P = 1; dz_PP = 2; dz_MM = -2; dz_M = -1;
      if (iz < 2 || iz > N-3) {
        if(iz==0)	  {dz_MM = N-2; dz_M  = N-1;}
        if(iz==1)	  {dz_MM = N-2;}
        if(iz==N-2)	{dz_PP = 2-N;}
        if(iz==N-1)	{dz_P  = 1-N; dz_PP = 2-N;}
      }
      dz_M *= N2; dz_P *= N2; dz_MM *= N2; dz_PP *= N2;
      gradient[i] = w[0]*phi[i + dz_P] + w[1]*phi[i + dz_PP] + w[2]*phi[i + dz_M] + w[3]*phi[i + dz_MM];
    }
  }
}
#endif
