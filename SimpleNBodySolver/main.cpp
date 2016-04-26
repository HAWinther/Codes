#include "global.h"
#include "multigrid_solver.h"
#include "readlua.h"
#include <random>
#include <algorithm>

//=========================================================
//
// Very simple N-body code for Newtonian and 
// Modified gravity simulations
// Hans Winther [hans.a.winther@gmail.com]
//
//=========================================================

//=========================================================
// Generate IC from a power-spectrum
//=========================================================

void create_ic(Particle *p, int npart){
  int npoints = 1000;
  double *kk, *pofk;
  
  // Read in and store P(k,z=0)
  kk   = new double[npoints];
  pofk = new double[npoints];
  ifstream fp("pofk.txt");
  for(int i = 0;i < npoints; i++){
    fp >> kk[i];
    fp >> pofk[i];
  }
  DSpline pofk_spline = DSpline(kk, pofk, npoints, "P(k) spline");

  // Integrate growth-factor
  // ... not implemented

  // Generate random numbers
  random_device rd;
  int seed = rd();
  mt19937 rg(seed);
  uniform_real_distribution<double> unidist(0, 1);
  int nrandom = npart * 3 * 2;
  double *randoms = new double[nrandom];
  for(int i = 0; i < nrandom; i++){
    randoms[i] = unidist(rg);
  }

  // Make realisation in Fourier space
  // ... not implemented

  // FFT back to get displacement field
  // ... not implemented

  // Assign particles and move them
  // ... not implemented

  delete[] randoms;
}

//=========================================================
// Hubble factor E[a] = H(a) / H_0
//=========================================================

double Hubble(double aexp){
  return sqrt( sim._omega_m / pow3(aexp) + sim._omega_lambda);
}

//=========================================================
// Compute the mapping from scale factor a to the super
// co-moving time variable
// d[tt] = d[t] H_0/ a^2 => tt(a) = Int_ai^a[ dx/[ E[x]x^3 ] ]
//=========================================================

void compute_time_scalefactor_spline(){
  double aexp_ini, aexp_end, da, integrand;
  double *tt, *aa;
  int nsteps;

  // Initialize
  nsteps   = 5000;
  aexp_ini = sim._aexp_ini/2.0;
  aexp_end = 2.0;
  da       = (aexp_end - aexp_ini)/double(nsteps-1);
  tt       = new double[nsteps];
  aa       = new double[nsteps];
  for(int i = 0; i < nsteps; i++) {
    tt[i] = 0.0;
    aa[i] = aexp_ini + da*i;
  }

  // Integrate up
  for(int i = 1; i < nsteps; i++){
    integrand = da/( Hubble(aa[i]) * pow3(aa[i]) );
    tt[i] = tt[i-1] + integrand;
  }

  // Shift to tt=0 for a = 0 (just for consistency)
  int aunityind = int( (1.0 - aexp_ini) / da + 0.5 );
  double t0 = tt[aunityind];
  for(int i = 0; i < nsteps; i++)
    tt[i] -= t0;

  // Make a spline
  sim._ttofa_spline = DSpline(aa, tt, nsteps, "The tt(a) function");
  sim._aoftt_spline = DSpline(tt, aa, nsteps, "The a(tt) function");

  delete[] aa;
  delete[] tt;
}

//=========================================================
// Output grid
//=========================================================

void output_grid(double *grid, int ngrid, string filename){
  int ngrid2    = pow2(ngrid);
  int ngrid_tot = pow3(ngrid);

  ofstream fp(filename.c_str());
  for(int i = 0;i < ngrid_tot; i++){
    double x, y, z;

    // Position of node
    x = (i % ngrid + 0.5)/double(ngrid);
    y = (i/ngrid % ngrid + 0.5)/double(ngrid);
    z = (i/ngrid2 % ngrid + 0.5)/double(ngrid);

    // Write grid-data to file
    fp << x << "  " << y << "  " << "  " << z << "  " << grid[i] << endl;
  }
}

void output_grid(fftw_complex *grid, int ngrid, string filename){
  int ngrid2    = pow2(ngrid);
  int ngrid_tot = pow3(ngrid);

  ofstream fp(filename.c_str());
  for(int i = 0;i < ngrid_tot; i++){
    double x, y, z;
    
    // Position of node
    x = (i % ngrid + 0.5)/double(ngrid);
    y = (i/ngrid % ngrid + 0.5)/double(ngrid);
    z = (i/ngrid2 % ngrid + 0.5)/double(ngrid);

    // Write grid-data to file (real-part of fft-grid)
    fp << x << "  " << y << "  " << "  " << z << "  " << grid[i][0] << endl;
  }
}

//=========================================================
// Output particles to ascii. Outputname is [output_a0.XX]
// where 0.XX = int[aexp * 100]/100
//=========================================================

void output_particles(double aexp){
  Particle *p = sim._Particle;
  double xfac = sim._BoxSize;
  double vfac = 100*sim._BoxSize/aexp;
  ostringstream strs;
  string filename;
  ofstream fp;

  // Make filename
  strs << double(int(aexp*100)/100.0);
  filename = "output_a" + strs.str();
  cout << "====> Dumping to file [" << filename << "]" << endl;
  
  // Write to file
  fp.open(filename.c_str());
  fp << sim._npart << endl;
  for(int i = 0;i < sim._npart; i++){
    fp << p[i]._x[0]*xfac << "  " << p[i]._x[1]*xfac << "  " << p[i]._x[2]*xfac << " ";
    fp << p[i]._v[0]*vfac << "  " << p[i]._v[1]*vfac << "  " << p[i]._v[2]*vfac << endl;
  }
}

//=========================================================
// Compute the gravitational force DPhi
// We first compute Phi(k) by FFT of Poissons equation
// D^2 Phi = 1.5 Omegam a (rho-1)
// Then we compute DPhi by inverse FFT
//
// Phi(k)   = -F[ 1.5 Omegam a (rho-1) ] / k^2
// (DPhi)_x = F^{-1}[ ik_x Phi(k) ]
//=========================================================

void compute_gravitational_acceleration(double aexp){
  int ngrid = sim._ngrid, ngrid_tot = sim._ngrid_tot, nover2 = ngrid/2;
  double phi_fac, laplacian_fac;

  cout << "====> Computing gravitational acceleration" << endl;

  // Pointers to force, phi and rho arrays
  fftw_complex *fx = sim._Force[0];
  fftw_complex *fy = sim._Force[1];
  fftw_complex *fz = sim._Force[2];
  fftw_complex *phi = sim._Phi;
  double *rho = sim._Rho;

  // Make Phi source D^2 Phi = 1.5 Om a (rho - 1)
  phi_fac = 1.5 * sim._omega_m * aexp; 
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i = 0; i < ngrid_tot; i++){
    phi[i][0] = phi_fac * (rho[i]-1.0);
    phi[i][1] = 0.0;
  }

  // Fourier transform the density-field
  fftw_execute(sim._plan_phi_forward);

  // Set zero mode to 0
  phi[0][0] = phi[0][1] = 0.0;
  fx[0][0]  = fx[0][1]  = 0.0;
  fy[0][0]  = fy[0][1]  = 0.0;
  fz[0][0]  = fz[0][1]  = 0.0;

  // Divide by the Laplacian
  laplacian_fac = -1.0/pow2(2.0*M_PI)/double(ngrid_tot);
#ifdef OMP
#pragma omp parallel for
#endif
  for(int i=0;i<ngrid;i++){
    int ii = (i < nover2 ? i: i-ngrid);
    for(int j=0;j<ngrid;j++){
      int jj = (j < nover2 ? j: j-ngrid);
      for(int k=0;k<ngrid;k++){
        int kk = (k < nover2 ? k: k-ngrid);
        int ind = i + ngrid*(j + k*ngrid);

        if(ind > 0){
          double norm = laplacian_fac / double(ii*ii + jj*jj + kk*kk);
          phi[ind][0] *= norm; 
          phi[ind][1] *= norm; 

          // Compute the fourier force: F(k) = ik Phi(k)
          norm = 2.0 * M_PI;
          fx[ind][0] = -norm * ii * phi[ind][1]; 
          fy[ind][0] = -norm * jj * phi[ind][1]; 
          fz[ind][0] = -norm * kk * phi[ind][1]; 

          fx[ind][1] = norm * ii * phi[ind][0]; 
          fy[ind][1] = norm * jj * phi[ind][0]; 
          fz[ind][1] = norm * kk * phi[ind][0]; 
        } 
      }
    }
  }

  // Fourier transform to get the force in real-space
  fftw_execute(sim._plan_force_backward[0]);
  fftw_execute(sim._plan_force_backward[1]);
  fftw_execute(sim._plan_force_backward[2]);
  
  // Calculate Phi (only needed for energy conservation check)
  fftw_execute(sim._plan_phi_backward);
}

//=========================================================
// Set analytical density field
//=========================================================

void set_analytical_density_field(double *rho, int ngrid){
  int ngrid2    = pow2(ngrid);
  int ngrid_tot = pow3(ngrid);
  double rho_avg, sig2, rho0;
  
  // Size and amplitude of profile
  sig2 = 0.01; rho0 = 10.0;

  rho_avg = 0.0;
  for(int i = 0; i < ngrid_tot; i++){
    double x, y, z, r2;

    // Position of node
    x  = (i % ngrid + 0.5)/double(ngrid);
    y  = (i/ngrid % ngrid + 0.5)/double(ngrid);
    z  = (i/ngrid2 % ngrid + 0.5)/double(ngrid);
    r2 = pow2(x-0.5) +  pow2(y-0.5) +  pow2(z-0.5);

    // Mexican top-hat at center of the box with size 'sig' and rho(0) = rho0
    rho[i] = 1.0 + rho0*(1.0 - r2/(3.0*sig2))*exp(-r2/(2.0*sig2));
    rho_avg += rho[i];
  }
  rho_avg /= double(ngrid_tot);

  cout << endl;
  cout << "Analytical density field. Average density in the box: " << rho_avg << endl;
}

//=========================================================
// Static test of the solvers
//=========================================================

void static_test(double aexp){

  // Assign density to grid [or use analytical formula]
  if(sim._analytic_matter_density){
    set_analytical_density_field(sim._Rho, sim._ngrid);
  } else {
    assign_particles_to_grid(sim._Rho, sim._ngrid);
  }

  // Solve for gravitational force
  compute_gravitational_acceleration(aexp);

  if(sim._modified_gravity){
    // Initial guess for modified gravity solver
    make_initial_guess(sim._ScalarField, aexp);

    // Solve for scalar field  
    solve_phi_with_multigrid();

    // Compute fifth-force
    calculate_fifth_force();
  }

  // Dump grids
  if(sim._modified_gravity)
    output_grid(sim._ScalarField, sim._ngrid, "scalarfield_grid.txt");
  output_grid(sim._Rho, sim._ngrid, "rho_grid.txt");
  output_grid(sim._Phi, sim._ngrid, "rho_grid.txt");
}

//=========================================================
// The main time-stepping loop. Everything is run from here
//=========================================================

void main_timestep_loop(){
  double aexp, aold, amid, t, told, tmid, dt, econs = 0.0;
  int nstep = 0;

  // Initialize time variables
  aexp = aold = sim._aexp = sim._aexp_ini;
  t    = sim._ttofa_spline.f(aexp);

  // Verbose
  cout << endl;
  cout << "=============================================================================" << endl;
  cout << "Starting time-integration aexp: " << setw(10) << aexp << endl;
  cout << "=============================================================================" << endl;

  // Set initial guess for modified gravity solver
  if(sim._modified_gravity)
    make_initial_guess(sim._ScalarField, aexp);

  // Main time-step loop
  while(aexp < 1.0){

    // Update time and translate between scale-factor and code time
    dt   = calc_new_timestep();
    if(t+dt > 0.0) dt = -t + 0.001;
    told = t;
    aold = sim._aoftt_spline.f(told);
    tmid = t + dt/2.0;
    amid = sim._aoftt_spline.f(tmid);
    t    = t + dt;
    aexp = sim._aoftt_spline.f(t);

    // Verbose
    cout << endl;
    cout << "=============================================================================" << endl;
    cout << "Timestep     nstep : " << setw(10) << nstep << "  aexp : " << setw(10) << aexp << endl;
    cout << "                 t : " << setw(10) << t     << "    dt : " << setw(10) << dt   << endl; 
    cout << "              vrms : " << setw(10) << calc_vrms(aexp);
    cout << " econs : " << setw(10) << econs << endl;
    cout << "=============================================================================" << endl;
    cout << endl;

    // Update positions of particles - 1/2 time-step
    leapfrog_drift_particles(dt/2.0); sim._aexp = amid;

    // Assign density to grid
    assign_particles_to_grid(sim._Rho, sim._ngrid);

    // Solve for gravitational force
    compute_gravitational_acceleration(amid);
    
    // Solve for the scalar field and compute fifth-force
    if(sim._modified_gravity){
      solve_phi_with_multigrid();
      calculate_fifth_force();
    }

    // Update velocity of particles - 1 time-step
    leapfrog_kick_particles(dt);

    // Update positions of particles - 1/2 time-step
    leapfrog_drift_particles(dt/2.0); sim._aexp = aexp;

    // Check energy conservation
    econs = calc_energy_conservation_constant();

    // Output when a = 0.1, 0.2, ...., 1.0
    if( (int(aexp*100) % 10 == 0) && ( int(aold*100) % 10 != 0) )
      output_particles(aexp);

    ++nstep;
  }
}

//=========================================================
// Make sure all particles are in [0,1] after we have 
// moved them
//=========================================================

void enforce_periodic_bc_particles(){
  Particle *curpart;
  int npart = sim._npart;
#ifdef OMP
#pragma omp parallel for private(curpart)
#endif
  for(int i = 0; i < npart; i++){
    curpart = &sim._Particle[i];
    for(int k = 0; k < 3; k++){
      if(curpart->_x[k] < 0.0) curpart->_x[k] += 1.0;
      if(curpart->_x[k] > 1.0) curpart->_x[k] -= 1.0;
    }
  }
}

//=========================================================
// Calculate the rms velocity Sqrt[ <v^2> ]
//=========================================================

double calc_vrms(double aexp){
  double vrms = 0.0;
#ifdef OMP
#pragma omp parallel for reduction(+:vrms)
#endif
  for(int i = 0; i < sim._npart; i++){
    for(int k = 0; k < 3; k++)
      vrms += pow2(sim._Particle[i]._v[k]);
  }
  return sqrt(vrms/double(sim._npart)) * (100.0 * sim._BoxSize) / aexp;
}

//=========================================================
// Calculate kinetic and potential (for GR only) energy
//=========================================================

pair<double,double> calc_energy(){
  int ngrid = sim._ngrid;
  double kinetic_energy = 0.0, potential_energy = 0.0;
#ifdef OMP
#pragma omp parallel for reduction(+:kinetic_energy) reduction(+:potential_energy)
#endif
  for(int i = 0; i < sim._npart; i++){
    int ix[3], ind;
    for(int k = 0; k < 3; k++){
      ix[k] = int(sim._Particle[i]._x[k] * ngrid) % ngrid;
      kinetic_energy += pow2(sim._Particle[i]._v[k]);
    }
    ind = ix[0] + ngrid*(ix[1]+ngrid*ix[2]);
  }
  kinetic_energy /= double(sim._npart);
 
  // Compute potential energy in box
  for(int i = 0; i < sim._ngrid_tot; i++)
    potential_energy += sim._Phi[i][0] * (sim._Rho[i] - 1.0);
  potential_energy /= double(sim._ngrid_tot);

  // Equivalent way (by integration by parts) to calculate potential energy
  //for(int i = 0; i < sim._ngrid_tot; i++)
  //  potential_energy += pow2(sim._Force[0][i][0]) + pow2(sim._Force[1][i][0]) + pow2(sim._Force[2][i][0]);
  //potential_energy /= - double(sim._ngrid_tot) * 1.5 * sim._omega_m * sim._aexp;

  return make_pair(kinetic_energy, potential_energy);
}

//=========================================================
// Check Layzer-Irvine energy conservation (for GR only)
// d/dt[T + U] - U/a = 0
//=========================================================

double calc_energy_conservation_constant(){
  double kinetic_energy, potential_energy, econst;
  pair<double,double> energy;

  static double epot_int     = 0.0;
  static double constant     = 0.0;
  static double epot_old     = 0.0;
  static double aexpold      = 0.0;
 
  // Calculate energy
  energy = calc_energy();
  kinetic_energy   = energy.first;
  potential_energy = energy.second;

  if(constant == 0.0){
    // Set initial energy
    constant  = kinetic_energy + potential_energy;
    aexpold   = sim._aexp;
    epot_old  = potential_energy;
    econst    = 0.0;
  } else {
    // Integrate up U/a
    epot_int += (potential_energy + epot_old) * 0.5 * log(sim._aexp/aexpold);
    epot_old  = potential_energy;
    aexpold   = sim._aexp;
    econst    = ( kinetic_energy + potential_energy - epot_int - constant) /
      (kinetic_energy +  fabs(potential_energy) + fabs(epot_int)); 
  }

  // Output energy file
  ofstream fp("energy.txt", std::ofstream::app);
  fp << sim._aexp << " " << econst << endl;
  fp.close();

  return econst;
}

//=========================================================
// Update positions of particles
// x' = v => x_{n+1} = x_n + v_{n+1/2} dt_code
//=========================================================

void leapfrog_drift_particles(double dt){
  Particle *curpart;
  double dx, dx_max = -1e30;

#ifdef OMP
#pragma omp parallel for private(curpart,dx) reduction(max : dx_max)
#endif
  for(int i = 0; i < sim._npart; i++){
    curpart = &sim._Particle[i];
    for(int k = 0; k < 3; k++){
      dx = curpart->_v[k] * dt;
      if(fabs(dx) > dx_max) dx_max = fabs(dx);
   
      // Move particle
      curpart->_x[k] += dx;
    }
  }
  cout << "====> Drift Particles | Maximum dx / dx_cell = " << dx_max*sim._ngrid << endl;
  
  // Particles that have moves outside the box must be moved
  enforce_periodic_bc_particles();

  // If we move across the whole box then something is wrong so just exit
  assert(dx_max < 1.0);
}

//=========================================================
// Move velocity of particles
// v' = F => v_{n+1/2} = v_{n-1/2} + F_{n} dt_code
//=========================================================

void leapfrog_kick_particles(double dt){
  Particle *curpart;
  double force[3];
  double dv, dv_max = -1e30;

  // Loop over all particles
#ifdef OMP
#pragma omp parallel for private(curpart,force,dv) reduction(max:dv_max)
#endif
  for(int i = 0; i < sim._npart; i++){
    curpart = &sim._Particle[i];
 
    // Calculates the force ( + fifth-force if modified gravity)
    // at the position of the particle
    calc_forces_grid(curpart, force);
    for(int k = 0; k < 3; k++){
      dv = -force[k] * dt;
      if(fabs(dv) > dv_max) dv_max = fabs(dv);
      
      // Kick velocity
      curpart->_v[k] += dv;
    }
  }
  cout << "====> Kick Particles  | Maximum dv (km/s) = " << dv_max * (100.0 * sim._BoxSize) / sim._aexp << endl;
}

//=========================================================
// Calculate new time-step based on a few conditions
//=========================================================

double calc_new_timestep(){
  double dt1, dt2, dt3, dt, cfac = 0.1;
  double rhomax = 0.0, vmax = 0.0;

  // Maximim density
#ifdef OMP
#pragma omp parallel for reduction(max:rhomax)
#endif
  for(int i = 0;i < sim._ngrid_tot; i++)
    if(sim._Rho[i] > rhomax) 
      rhomax = sim._Rho[i];
  
  // Maximum velocity
#ifdef OMP
#pragma omp parallel for reduction(max:vmax)
#endif
  for(int i = 0;i < sim._npart; i++)
    for(int k = 0; k < 3; k++)
      if(fabs(sim._Particle[i]._v[k]) > vmax) 
        vmax = fabs(sim._Particle[i]._v[k]);

  // dt based gravity free fall time
  dt1 = cfac * sqrt(pow2(M_PI)/(4.0 * sim._omega_m * sim._aexp * rhomax));

  // dt based on background expansion
  dt2 = 0.1/(sim._aexp*Hubble(sim._aexp));

  // dt based on courant condition for max velocity (do not allow travel more than 1/2 gridcell)
  dt3 = cfac/(double(sim._ngrid) * vmax);

  // Take the minimum time-step
  dt = min(dt1,dt2);
  dt = min(dt,dt3);

  return dt;
}

//=========================================================
// Extract newtonian force at the position of particle
// using NGP, CIC or TSC
//=========================================================

inline void calc_forces_grid(Particle *p, double *f){
  int ix[3], ind, ngrid = sim._ngrid;

#if defined(CIC)
  int sign[3], ixnbor[3]; 
  double xnbor[3], xp[3], weight;
  f[0] = f[1] = f[2] = 0.0;

  // Position of particle / dx and index of nearest grid-point
  for(int k=0;k<3;k++){
    xp[k] = p->_x[k] * ngrid;
    ix[k] = int(xp[k]);
  }

  // Sign tells us if to go left or right in every direction
  for(int k = 0; k < 3; k++)
    sign[k] = xp[k] - (ix[k]+0.5) > 0.0 ? +1 : -1;

  // Loop over 8 nodes containing particle cloud
  for(int ii=0;ii<8;ii++){
    for(int k=0,n=1;k<3;k++,n*=2){
      ixnbor[k] = ix[k] + sign[k] * (ii/n % 2);
      if(ixnbor[k] < 0)      ixnbor[k] += ngrid;
      if(ixnbor[k] >= ngrid) ixnbor[k] -= ngrid;
      xnbor[k] = ixnbor[k] + 0.5;
    }
    ind = ixnbor[0] + ngrid*(ixnbor[1] + ngrid*ixnbor[2]);

    // Compute CIC weight
    weight = cic_kernel(xnbor, xp);

    // Add up force
    for(int k = 0; k < 3; k++)
      f[k] += sim._Force[k][ind][0] * weight;

    // Add fifth-force
    if(sim._modified_gravity){
      for(int k = 0; k < 3; k++)
        f[k] += sim._FifthForce[k][ind] * weight;
    }
  }
#elif defined(TSC)
  int ixnbor[3]; 
  double xnbor[3], xp[3], weight;
  f[0] = f[1] = f[2] = 0.0;

  // Position of particle / dx and index of nearest grid-point
  for(int k=0;k<3;k++){
    xp[k] = p->_x[k] * ngrid;
    ix[k] = int(xp[k]);
  }

  // Loop over all 27 neighbor grid-points
  for(int ii=0;ii<27;ii++){
    for(int k=0,n=1;k<3;k++,n*=3){
      ixnbor[k] = ix[k] - 1 + (ii/n % 3);
      if(ixnbor[k] < 0)      ixnbor[k] += ngrid;
      if(ixnbor[k] >= ngrid) ixnbor[k] -= ngrid;
      xnbor[k] = ixnbor[k] + 0.5;
    }
    ind = ixnbor[0] + ngrid*(ixnbor[1] + ngrid*ixnbor[2]);

    // Compute TSC weight
    weight = tsc_kernel(xnbor, xp);

    // Add up force
    for(int k = 0; k < 3; k++)
      f[k] += sim._Force[k][ind][0] * weight;

    // Add fifth-force
    if(sim._modified_gravity){
      for(int k = 0; k < 3; k++)
        f[k] += sim._FifthForce[k][ind] * weight;
    }
  }
#else
  // Calculate index of nearest grid-point
  for(int k = 0; k < 3; k++)
    ix[k] = int(p->_x[k] * sim._ngrid) % sim._ngrid;
  ind = ix[0] + sim._ngrid*(ix[1] + sim._ngrid * ix[2]);

  // Add up force
  for(int k = 0; k < 3; k++)
    f[k] = sim._Force[k][ind][0];

  // Add fifth-force
  if(sim._modified_gravity){
    for(int k = 0; k < 3; k++)
      f[k] += sim._FifthForce[k][ind];
  }
#endif
}

//=========================================================
// CIC and TSC kernels
// yg position of node
// xp positions of the interpolation point
//=========================================================

inline double cic_kernel(double *xg, double *xp){
  double w = 1.0;
  for(int k=0;k<3;k++){
    double dx = fabs(xp[k] - xg[k]);
    w *= (dx < 1.0 ? 1.0 - dx: 0.0);
  }
  return w;
}

inline double tsc_kernel(double *xg, double *xp){
  double w = 1.0;
  for(int k=0;k<3;k++){
    double dx = fabs(xp[k] - xg[k]);
    w  *= (dx < 0.5) ? 0.75 - pow2(dx) : ( dx < 1.5 ? 0.5 * pow2(1.5-dx) : 0.0 );
  }
  return w;
}

//=========================================================
// Simple nearest grid-point (NGP) or cloud-in-cell (CIC)
// or triangular-shaped-cloud (TSC) density assignment. 
// Bins particle positions to get Rho(x,y,z)
// Density is normalized to the background density
//
// Compile with -DCIC to use CIC or -DTSC to use TSC
//=========================================================

void assign_particles_to_grid(double *rho, int ngrid){
  int ngrid_tot = pow3(ngrid), npart = sim._npart, ix[3], ind;
  double rho_norm_fac, rho_avg;

  // Reset density array
#ifdef OMP
#pragma omp parallel for 
#endif
  for(int i = 0; i < ngrid_tot; i++) 
    rho[i] = 0.0;

  // Loop over all particles and assign to grid
  for(int i = 0; i < npart; i++){
#if defined(CIC)
    int ixnbor[3], sign[3];
    double xnbor[3], xp[3];

    // Position of particle / dx and index of nearest grid-point
    for(int k=0;k<3;k++){
      xp[k] = sim._Particle[i]._x[k] * ngrid;
      ix[k] = int(xp[k]);
    }

    // Sign tells us if cloud is located left or right of node
    for(int k = 0; k < 3; k++)
      sign[k] = xp[k] - (ix[k]+0.5) > 0.0 ? +1 : -1;

    // Loop over all 8 neighbor nodes containing particle cloud
    for(int ii=0;ii<8;ii++){
      for(int k=0,n=1;k<3;k++,n*=2){
        ixnbor[k] = ix[k] + sign[k] * (ii/n % 2);
        if(ixnbor[k] < 0)      ixnbor[k] += ngrid;
        if(ixnbor[k] >= ngrid) ixnbor[k] -= ngrid;
        xnbor[k] = ixnbor[k] + 0.5;
      }
      ind = ixnbor[0] + ngrid*(ixnbor[1] + ngrid*ixnbor[2]);
      rho[ind] += cic_kernel(xnbor, xp);
    }
#elif defined(TSC)
    int ixnbor[3]; 
    double xnbor[3], xp[3];

    // Position of particle / dx and index of nearest grid-point
    for(int k=0;k<3;k++){
      xp[k] = sim._Particle[i]._x[k] * ngrid;
      ix[k] = int(xp[k]);
    }

    // Loop over all 27 neighbor grid-points
    for(int ii=0;ii<27;ii++){
      for(int k=0,n=1;k<3;k++,n*=3){
        ixnbor[k] = ix[k] - 1 + (ii/n % 3);
        if(ixnbor[k] < 0)      ixnbor[k] += ngrid;
        if(ixnbor[k] >= ngrid) ixnbor[k] -= ngrid;
        xnbor[k] = ixnbor[k] + 0.5;
      }
      ind = ixnbor[0] + ngrid*(ixnbor[1] + ngrid*ixnbor[2]);
      rho[ind] += tsc_kernel(xnbor, xp);
    }
#else
    // Compute index of nearest grid-point
    for(int k = 0; k < 3; k++){
      ix[k] = int(sim._Particle[i]._x[k] * ngrid) % ngrid;
    }
    ind = ix[0] + ngrid*(ix[1] + ngrid*ix[2]);
    rho[ind] += 1.0;
#endif
  }

  // Factor to translate from number of particles in cell to rho/rho_background
  rho_norm_fac  = 1.0/double(npart) * pow3(ngrid);

  // Normalize to get rho/rho_background
  rho_avg = 0.0;
#ifdef OMP
#pragma omp parallel for reduction(+:rho_avg)
#endif
  for(int i = 0; i < ngrid_tot; i++){
    rho[i] *= rho_norm_fac;
    rho_avg += rho[i];
  }
  rho_avg /= double(ngrid_tot);

  cout << "====> Assigning particles to grid. Average density-to-mean-density in the box is: " << rho_avg << endl;
}

//=========================================================
// Read a parameterfile
//=========================================================

void ReadFileWithLua::read(){
  bool status;

  // Extract data from LUA parameterfile
  sim._aexp = sim._aexp_ini    = read_double("aexp_ini",         &status, true);
  sim._BoxSize                 = read_double("BoxSize",          &status, true);
  sim._courant_fac             = read_double("courant_fac",      &status, true);
  sim._omega_m                 = read_double("omega_m",          &status, true);
  sim._ngrid                   = read_int   ("ngrid",            &status, true);
  sim._NGS_FINE                = read_int   ("NGS_FINE",         &status, false);
  if(!status) sim._NGS_FINE    = 5;
  sim._NGS_COARSE              = read_int   ("NGS_COARSE",       &status, false);
  if(!status) sim._NGS_COARSE  = 5;
  sim._EPS_CONVERGE            = read_int   ("EPS_CONVERGE",     &status, false);
  if(!status) sim._EPS_CONVERGE= 1e-5;
  sim._static                  = read_bool  ("static",           &status, true);
  if(!status) sim._static      = false;
  sim._verboseMG               = read_bool  ("verboseMG",        &status, true);
  sim._analytic_matter_density = read_bool  ("verboseMG",        &status, true);
  sim._modified_gravity        = read_bool  ("modified_gravity", &status, false);
  if(!status) sim._modified_gravity = false;
#if defined(FOFR)
  sim._nfofr                   = read_double("nfofr",            &status, true);
  sim._fofr0                   = read_double("fofr0",            &status, true);
#elif defined(SYMMETRON)
  sim._assb_symm               = read_double("assb_symm",        &status, true);
  sim._beta_symm               = read_double("beta_symm",        &status, true);
  sim._lambda_symm             = read_double("lambda_symm",      &status, true);
#endif
  sim._particlefile            = read_string("particlefile",     &status, false);

  // If no particlefile assume static-mode
  if(!status) sim._static      = sim._analytic_matter_density = true;

  // Derived parameters and checks
  sim._omega_lambda = 1.0 - sim._omega_m;
  sim._ngrid_tot = pow3( sim._ngrid );
  assert( sim._ngrid == (1 << int(log2(sim._ngrid))) );
  if(sim._analytic_matter_density) sim._npart = 0;
}

//=========================================================
// Initialize the simulation
// For a restart just change particle-file to output_a0.XX
// file and change aexp_init = 0.XX
//=========================================================

void init(int argc, char **argv){

  // Read parameterfile or use standard parameters
  if(argc == 1){
    cout << endl;
    cout << "=============================================================================" << endl;
    cout << "Initializing sim using standard values" << endl; 
    cout << "=============================================================================" << endl;

    sim._aexp_ini = sim._aexp = 0.1;        // Starting scalefactor
    sim._courant_fac = 0.5;                 // Time-step control (lower-values = smaller timestep)
    sim._BoxSize = 200.0;                   // Boxsize in Mpc/h
    sim._omega_m = 0.3;                     // Cosmological parameters
    sim._omega_lambda = 1.0 - sim._omega_m;
    sim._particlefile = "ic/part.ascii";    // Particle file
    sim._ngrid     = 64;                    // Grid-size. Must be a power of 2
    sim._ngrid_tot = pow3( sim._ngrid );
    assert( sim._ngrid == (1 << int(log2(sim._ngrid))) );
    sim._modified_gravity = false;          // Modified gravity sim?
    sim._static = false;                    // Solve at a single time
    sim._analytic_matter_density = false;   // Set density analytically. 
                                            // Used for static_test only (do not require particles)
#if defined(FOFR)
    sim._nfofr = 1.0;                       // f(R) parameters
    sim._fofr0 = 1.0e-6;                    
#elif defined(SYMMETRON)
    sim._assb_symm   = 0.5;                 // Symmetron parameters
    sim._lambda_symm = 5.0;                 // Range in Mpc/h 
    sim._beta_symm   = 1.0;                 // Coupling strength
#endif

                                            // Multigrid solver parameters
    sim._EPS_CONVERGE = 1.0e-7;             // Convergence criterion
    sim._NGS_FINE     = 20;                 // Number of sweeps per V-cycle step (finest level)
    sim._NGS_COARSE   = 10;                 // (coarser levels)
    sim._verboseMG    = false;              // Add (alot) of verbose for the mg solver
  } else {
    cout << endl;
    cout << "=============================================================================" << endl;
    cout << "Initializing sim from parameterfile" << endl; 
    cout << "=============================================================================" << endl;
    cout << "Reading parameters from: " << argv[1] << endl;  

    ReadFileWithLua paramfile;
    paramfile.open_read_close(argv[1]);
  } 
  if(!sim._static) sim._analytic_matter_density = false;
 
  // Read top of particle-file to get npart
  if(!sim._analytic_matter_density){
    ifstream fp(sim._particlefile.c_str());
    fp >> sim._npart;
    fp.close();
  } else {
    sim._npart = 0;
  }

  // Allocate all memory
  allocate_memory();

  // Read particles from file
  if(sim._npart > 0)
    read_particles_from_file(sim._particlefile);
  
  // Create spline(s) to translate between code time and scale-factor
  compute_time_scalefactor_spline();

  // Verbose
  cout << endl;
  cout << "Running with the following parameters: " << endl;
  cout << "Npart     = " << setw(10) << sim._npart << endl;
  cout << "Ngrid     = " << setw(10) << sim._ngrid << endl;
  cout << "Omega_m   = " << setw(10) << sim._omega_m << endl;
  cout << "Omega_l   = " << setw(10) << sim._omega_lambda << endl;
  cout << "aexp_ini  = " << setw(10) << sim._aexp_ini << endl;
  cout << "BoxSize   = " << setw(10) << sim._BoxSize << " Mpc/h" << endl;
  cout << "Courantfac= " << setw(10) << sim._courant_fac << endl;
  if(sim._static)
    cout << "Will run in static mode" << endl;
  if(sim._modified_gravity){
    cout << endl;
#if defined(FOFR)
    cout << "Modified gravity simulation - Hu-Sawkicy f(R) gravity:" << endl;
    cout << "fofr0     = " << sim._fofr0 << endl;
    cout << "n         = " << sim._nfofr << endl;
#elif defined(SYMMETRON)
    cout << "Modified gravity simulation - The symmetron model:" << endl;
    cout << "assb      = " << sim._assb_symm   << endl;
    cout << "beta      = " << sim._beta_symm   << endl;
    cout << "lambda    = " << sim._lambda_symm << endl;
#elif defined(GR)
    cout << "Modified gravity simulation - GR as MG!" << endl;
#else
    .................
    ...Not defined...
    .................
#endif
    cout << endl;
    cout << "Multigrid parameters:" << endl;
    cout << "epsilon   = " << sim._EPS_CONVERGE << endl;
    cout << "ngs_fine  = " << sim._NGS_FINE << endl;
    cout << "ngs_corse = " << sim._NGS_COARSE << endl;
  }

}

//=========================================================
// Allocate all memory
//=========================================================

void allocate_memory(){
  double mem_newtonian, mem_scalar, mem_particles, mem_multigrid, mem_tot;

  // Calculate the amount of memory needed
  mem_newtonian = sim._ngrid_tot * sizeof(fftw_complex) * 5 / pow2(1024.);
  mem_particles = sim._npart * sizeof(double) * 6 / pow2(1024.);
  if(sim._modified_gravity){
    mem_scalar    = sim._ngrid_tot * sizeof(double) * 6 / pow2(1024.);
    mem_multigrid = sim._ngrid_tot * sizeof(double) * 4 / pow2(1024.);
  } else {
    mem_scalar = mem_multigrid = 0.0;
  }
  mem_tot = mem_newtonian + mem_scalar + mem_particles + mem_multigrid;

  // Verbose
  cout << endl;
  cout << "=============================================================================" << endl;
  cout << "Memory needed for simulation:" << endl;
  cout << "=============================================================================" << endl;
  cout << "Newtonian        : " << setw(10) << mem_newtonian << " Mb" << endl;
  cout << "Modified Gravity : " << setw(10) << mem_scalar    << " Mb" << endl;
  cout << "MultiGrid Arrays : " << setw(10) << mem_multigrid << " Mb" << endl;
  cout << "Total            : " << setw(10) << mem_tot       << " Mb" << endl;

#ifdef OMP_FFT
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  // FFT related variables needed for solving for Newtonian force
  sim._Rho      = new double[sim._ngrid_tot];
  sim._Phi      = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * sim._ngrid_tot );
  sim._Force[0] = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * sim._ngrid_tot );
  sim._Force[1] = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * sim._ngrid_tot );
  sim._Force[2] = (fftw_complex *) fftw_malloc( sizeof(fftw_complex) * sim._ngrid_tot );

  // Make FFT plans
  sim._plan_phi_forward  = fftw_plan_dft_3d(sim._ngrid, sim._ngrid, sim._ngrid, sim._Phi, sim._Phi, FFTW_FORWARD , FFTW_MEASURE); 
  sim._plan_phi_backward = fftw_plan_dft_3d(sim._ngrid, sim._ngrid, sim._ngrid, sim._Phi, sim._Phi, FFTW_BACKWARD, FFTW_MEASURE); 
  for(int k = 0; k < 3; k++)
    sim._plan_force_backward[k] = fftw_plan_dft_3d(sim._ngrid, sim._ngrid, sim._ngrid, sim._Force[k], sim._Force[k], FFTW_BACKWARD, FFTW_MEASURE); 

  // Allocate for modified gravity
  if(sim._modified_gravity){
    sim._ScalarField   = new double[sim._ngrid_tot];
    sim._MGSource      = new double[sim._ngrid_tot];
    sim._MGResidual    = new double[sim._ngrid_tot];
    sim._FifthForce[0] = new double[sim._ngrid_tot];
    sim._FifthForce[1] = new double[sim._ngrid_tot];
    sim._FifthForce[2] = new double[sim._ngrid_tot];

    // Initialize grid
    for(int i = 0; i < sim._ngrid_tot; i++)
      sim._ScalarField[i] = sim._MGSource[i] = sim._MGResidual[i] = 0.0;

    // Init MultiGrid solver
    sim._mg = new MultiGrid( sim._ScalarField, sim._MGResidual, sim._Rho, sim._MGSource, sim._ngrid , sim._verboseMG );
  }

  // Allocate memory for particles
  if(sim._npart > 0)
    sim._Particle = new Particle[sim._npart];
}

//=========================================================
// Free up all the memory
//=========================================================

void deallocate_memory()
{
  // Deallocate arrays
  if(sim._modified_gravity){
    delete[] sim._ScalarField;
    delete[] sim._MGSource;
    delete[] sim._MGResidual;
    delete[] sim._FifthForce[0];
    delete[] sim._FifthForce[1];
    delete[] sim._FifthForce[2];
  }
  delete[] sim._Rho;
  delete[] sim._Particle;
  delete   sim._mg;
  
  // Deallocate FFT arrays
  fftw_free(sim._Phi);
  fftw_free(sim._Force[0]);
  fftw_free(sim._Force[1]);
  fftw_free(sim._Force[2]);

  // Destroy FFT plans
  fftw_destroy_plan(sim._plan_phi_forward);
  fftw_destroy_plan(sim._plan_phi_backward);
  for(int k=0;k<3;k++)
    fftw_destroy_plan(sim._plan_force_backward[k]);
}

//=========================================================
// Read particles from file. Assuming here a simple 
// ascii-file with format:
//
// npart
// x  y  z  vx  vy  vz
// x  y  z  vx  vy  vz
// ...
//
// [x] in units of [Mpc/h] 
// [v] in units of [km/s]
//=========================================================

void read_particles_from_file(string particlefile){
  ifstream fp(particlefile.c_str());
  double max_pos, min_pos, max_vel, rms_vel, vfac, xfac, vel2;
  double x[3], v[3];

  // Abort if file do not exist
  assert(fp != NULL);

  // Factor to transform v from km/s to a * v / (H0B0)
  vfac = sim._aexp_ini / (100.0 * sim._BoxSize);

  // Factor to transform x from Mpc/h to [x/BoxSize]
  xfac = 1.0/sim._BoxSize;

  // Initialize max/min
  max_vel = -1e30;
  max_pos = -1e30;
  min_pos = 1e30;
  rms_vel = 0.0;

  // Verbose
  cout << endl;
  cout << "=============================================================================" << endl;
  cout << "Reading particles from file: " << particlefile << endl;
  cout << "=============================================================================" << endl;

  fp >> sim._npart;
  cout << "We have npart = " << sim._npart << endl;

  // Read all particles
  for(int i = 0; i < sim._npart; i++){
    for(int k = 0; k < 3 ; k++) fp >> x[k];
    for(int k = 0; k < 3 ; k++) fp >> v[k];

    // Verbose
    if(i < 3 || i > sim._npart - 4){
      cout << "Particle: " << setw(10) << i << " Pos: " << setw(10) << x[0] << " ";
      cout << setw(10) << x[1] << " " << setw(10) << x[2] << endl;
      if(i==2) cout << "..." << endl;
    }

    // Transform to code-units
    for(int k = 0; k < 3 ; k++){
      x[k] *= xfac;
      v[k] *= vfac;

      // Enforce periodic boundary conditions
      if(x[k] >= 1.0) x[k] -= 1.0;
    }

    // Compute some quantities
    vel2 = 0.0;
    for(int k = 0; k < 3 ; k++){
      if(x[k] > max_pos) max_pos = x[k];
      if(x[k] < min_pos) min_pos = x[k];
      vel2 += pow2(v[k]);
    }
    rms_vel += vel2;
    if(vel2 > max_vel) max_vel = vel2;

    // Make new particle
    sim._Particle[i] = Particle(x, v);
  }
  rms_vel = sqrt(rms_vel / double(sim._npart)) / vfac;
  max_vel = sqrt(max_vel) / vfac;

  // Output some info about the particles read
  cout << endl;
  cout << "Done reading " << sim._npart << " particles" << endl;
  cout << "Min x: " << setw(10) << min_pos << " Max x : " << setw(10) << max_pos <<  endl;
  cout << "Max v: " << setw(10) << max_vel << " RMS(v): " << setw(10) << rms_vel << " km/s" << endl;

  // Check for errors
  assert(min_pos >= 0.0);
  assert(max_pos < 1.0);
}

//=========================================================
// Main is main...
//=========================================================

int main(int argc, char **argv){

  // Initialize the simulation
  init(argc, argv);

  // Static test
  if(sim._static){
    static_test(sim._aexp);
    deallocate_memory();
    exit(1);
  }
   
  // Perform a full simulation
  main_timestep_loop();
  deallocate_memory();
}

