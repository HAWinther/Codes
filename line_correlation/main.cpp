#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <complex>
#include <limits.h>
#include <unistd.h>
#include <algorithm>
#include <vector>
#include <iomanip> 
#include <ctime>
#include "spline.h"
#ifdef CXXRANDOM
  #include <random>
#endif
#ifdef WITHMPI
  #include <mpi.h>
#endif

using namespace std;

#ifdef CXXRANDOM
mt19937 rg;
random_device rd;
#endif

//============================================================
// Code description and units
//============================================================
/*
 Code to compute line correction functions from N-body data.
 Can also compute P(k) and zeta(r)
 
 The gridsize used for the density field is hardcoded below.
 Grid-assignment used is NGP.

 For P(k) we deconvolve the window function and subtract shot-noise

 The internal code units are:
   Positions in units where boxsize = 1 (r_code = r/Boxsize)
   k is usually integer wavenumber k_phys = 2pi k / Boxsize
   In some places k = 2pi * k_int i.e. in units of boxsize = 1
   In outputs and boxsize etc. we always use Mpc/h!

 Memory needed per CPU: ~ [2(NGRID/512)^3 + 3(NPART/512)^3 * O(1)/NFILES ] GB
 where NFILES is the number of RAMSES input files. The RAMSES files is read
 on the fly and particles added to the grid (to save memory).
 After particles are read then this memory is freed.

 Line correlation function:
   For large r the double sum in l(r) is calculated exactly.
   This is a O(kmax^6) = O(1/r^6) algorithm so very slow for small 'r'.

   For small r we use Monte-Carlo integration. The stopping
   criterion is set in EPSILON_MC_INTEGRATE. The maximum number 
   of samples (per CPU) is set in NSAMPLEMAX_HARDCODED

 To save computation time we precompute j0(2pi Sqrt(x)) on x=[0,5]
 with NPOINTS_J0_PRECOMP points

 Use proper random number generator (Mersenne twister) with CXXRANDOM 

 Define WITHMPI to run in parallel. 
 NB: All CPUs have their own copy of the grid.

 Compile as:
 mpicc main.cpp -I/../include -L/../local/lib -lfftw3 -DWITHMPI -O3 -o lcorr_mpi
   icc main.cpp -I/../include -L/../local/lib -lfftw3 -O3 -o lcorr_nompi

 Run as:
 ./lcorr_nompi /path/to/folder/containing/output_xxxxx/ xxxxx numfiles npart boxsize nbins outputprefix
 mpirun -np N ./lcorr_mpi /path/to/folder/containing/output_xxxxx/ xxxxx numfiles npart boxsize nbins outputprefix

 Time per run:
  For NGRID=512 it usually takes ~30min on 10 CPUs to compute the line-correlation function [plus P(k) / zeta(r)]
  with nbins=20 and accuracy EPSILON = 0.5%. Most of the time is used on the low-r values (r < Boxsize * 10/NGRID).

*/

//============================================================
// Hardcoded constants
//============================================================

// Grid size to use in the analysis
#define NGRID_HARDCODED      512

// Max number of samples in MC integration
#define NSAMPLEMAX_HARDCODED 1000000000000

// Check for convergence every N samples in MC integration
#define NSAMPLECONVCHECK     100000000

// Convergence criterion for MC integration "Est. error / mean value < epsilon"
#define EPSILON_MC_INTEGRATE 0.005

// Use precomputed lookup-function for j0. Much faster
#define USEBESSELLOOKUP
#define BESSELXMAX 5.0

// Sinc is precomputed to avoid computing billions of sin(x)
#define NPOINTS_J0_PRECOMP   int(1e6)

// Density assignment (for deconvolving the window function) 
// NGS=1 CIC=2 TSC=3. Only NGP is implemented.
#define DENSITYASSIGNMENT 1

#define TWOPI                6.283185307
#define pow2(x)              ((x)*(x))
#define pow3(x)              ((x)*(x)*(x))

//============================================================
// All methods in this code...
//============================================================

void   init(int argv, char **argc);
void   myexit(int i, string err);

//-- Read methods
void   read_and_bin_particles(int nparttot, string filepath, int filenum, int nfiles);
int    read_int       (FILE* fp);
void   read_int_vec   (FILE* fp, int *buffer,    int n);
void   read_double_vec(FILE* fp, double *buffer, int n);
inline string ramses_num_as_string(int i);

//-- FFT methods
void   compute_phase_and_power( bool returnreal = false);
inline void add_to_pofk(int index, fftw_complex* grid);

//-- Line correlation methods
double calc_line_corr_single(double r);
void   calc_line_corr_all();

//-- Monte-Carlo integration of l(r) double integral
double monte_carlo_integrate_line_corr(double r);
double estimate_progress_monte_carlo(double x, double y, bool reset = false);

//-- Power-spectrum methods
void   compute_correlation_function();

//-- Spherical besselfunction j0
void   precompute_j0(int npoints, double xmax);
inline double j0_func(double x);

//-- Output to file
void output(bool output_l, bool output_pofk, bool output_zeta);

//============================================================
// Read methods
//============================================================

// Reads one single integer
int read_int(FILE* fp){
  int tmp, skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(&tmp, sizeof(int), 1, fp);
  fread(&skip, sizeof(int), 1, fp);
  return tmp;
}

// Reads an integer vector
void read_int_vec(FILE* fp, int *buffer, int n){
  int skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(buffer, sizeof(int), n, fp);
  fread(&skip, sizeof(int), 1, fp);
}

// Reads a double vector
void read_double_vec(FILE* fp, double *buffer, int n){
  int skip;
  fread(&skip, sizeof(int), 1, fp);
  fread(buffer, sizeof(double), n, fp);
  fread(&skip, sizeof(int), 1, fp);
}

// Converts a number X to a 5 digit string '0000X'
inline string ramses_num_as_string(int i){
  char *filenameout = new char[200];
  if(i<10)
    sprintf(filenameout,"0000%i",i);
  else if(i<100)
    sprintf(filenameout,"000%i",i);
  else if(i<1000)
    sprintf(filenameout,"00%i",i);
  else if(i<10000)
    sprintf(filenameout,"0%i",i);
  else{
    cout << "ramses_num_as_string not implemented for 10000 <= " << i << endl;
    myexit(0,"ramsesnum");
  }
  return string(filenameout);
}

// Exit from MPI
void myexit(int i, string err){
#ifdef WITHMPI
  printf("%s\n",err.c_str());
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD,i);
#endif
  exit(i);
}

//============================================================
// All global information. Handled all allocation and dealloc
// of memory apart from temp read buffers
//============================================================

class Sim{
  public:
  // Grid related data
  fftw_complex *grid;
  int ngrid, ngridtot;
  bool alloc_grid;

  // Power-spectrum and 2-point correlation function
  int npofk;
  double *pofk, *pofk_count, *pofk_k;
  bool alloc_pofk;
  int nzeta;
  double *zeta, *rzeta;

  // For RAMSES input files
  string filepath;
  int filenum;
  int nfiles;

  // Output filename
  string outputprefix;

  // Boxsize
  double boxsize;

  // Binning parameters
  double *r, *l;
  double rmin, rmax;
  int nbins;
  bool alloc_bins;

  // MPI
  int myid, ncpu;

  // Precomputed j0
  double *precomp_j0;
  double  precomp_j0_fac;
  int     precomp_j0_npoints;
  bool    alloc_precomp;

  Sim(){ alloc_grid = alloc_precomp = alloc_bins = alloc_pofk = false; }
  ~Sim(){
    deallocate_grid();
    deallocate_precomp();
    deallocate_bins();
    deallocate_pofk();
  }

  // Allocate/deallocate memory for binning
  void allocate_bins(int n, double r1, double r2){
    // Sanity check
    if(n > 1e3 || n <= 0){
      cout << "Do we really want " << n << " bins?" << endl;
      myexit(1,"allocbins");
    }
    nbins = n;
    rmin = r1;
    rmax = r2;
    l = new double[n];
    r = new double[n];
    for(int i=0;i<n;i++) l[i] = r[i] = 0.0;
    alloc_bins = true;
  }
  void deallocate_bins(){
    if(alloc_bins) delete[] r,l;
    alloc_bins = false;
  }

  // Allocate/deallocate memory for fftw grid
  void allocate_grid(int n){
    // Sanity check
    if(n > 1024 || n < 32){
      cout << "Do we really want " << n << "^3 grid points?" << endl;
      myexit(2,"allocgrid");
    }
    ngrid = n;
    ngridtot = pow3(n);
    grid = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * ngridtot);
    alloc_grid= true;
    for(int i=0;i<n*n*n;i++) grid[i][0] = grid[i][1] = 0.0;
  }
  void deallocate_grid(){
    if(alloc_grid) fftw_free(grid);
    alloc_grid = false;
  }

  // Allocate/deallocate memory for fftw grid
  void allocate_precomp(int n, double fac){
    // Sanity check
    if(n > 1e9 && n < 100){
      cout << "Do we really want " << n << " precomp points?" << endl;
      myexit(3,"allocprecomp");
    }
    precomp_j0_npoints = n;
    precomp_j0_fac = fac;
    precomp_j0 = new double[n];
    alloc_precomp = true;
    for(int i=0;i<n;i++) precomp_j0[i] = 0.0;
  }
  void deallocate_precomp(){
    if(alloc_precomp) delete[] precomp_j0; 
    alloc_precomp = false;
  }

  // Allocate bins for pofk
  void allocate_pofk(int npower, int ntwoptcorr){
    // Sanity check
    if(npower>1e9 || npower <= 0 || ntwoptcorr>1e9 || ntwoptcorr<= 0){
      cout << "Do we really want " << npower << " / " << ntwoptcorr << " pofk/zeta points?" << endl;
      myexit(4,"allocpofk");
    }
    pofk       = new double[npower];
    pofk_k     = new double[npower];
    pofk_count = new double[npower];
    zeta       = new double[ntwoptcorr];
    rzeta      = new double[ntwoptcorr];
    for(int i=0;i<npower;i++){
      pofk[i]   = pofk_count[i] = 0.0;
      pofk_k[i] = TWOPI * i;
    }
    for(int i=0;i<ntwoptcorr;i++)
      zeta[i] = rzeta[i] = 0.0;
    npofk = npower;
    nzeta = ntwoptcorr;
    alloc_pofk = true;
  }
  void deallocate_pofk(){
    if(alloc_pofk) delete[] pofk, pofk_count, pofk_k, zeta, rzeta;
    alloc_pofk = false;
  }

} global;

//============================================================
// Initialize
//
// INPUT ARGUMENTS:
// argc  1   filepath
// argc  2   filenum
// argc  3   nfiles
// argc  4   nparticles 
// argc  5   boxsize Mpc/h
// argc  6   nbins
// argc  7   outputfilename
//
//============================================================

void init(int argv, char **argc){
  char cwd[1024];
  FILE *fp;
  long seed, *allseeds;

  // Initialize MPI
#ifdef WITHMPI
  MPI_Init(&argv,&argc);
  MPI_Comm_rank(MPI_COMM_WORLD, &global.myid);
  MPI_Comm_size(MPI_COMM_WORLD, &global.ncpu);
#else
  global.myid = 0;
  global.ncpu = 1;
#endif

  // Check we have enough input
  if(argv < 7){
    cout << "Error in init. Not enough data: " << endl;;
    for(int i=0;i<argv;i++) cout << i << " : " << argc[i] << endl;
    myexit(5,"init");
  }

  // Initialize random number generator. Set random seed using urandom or by time + id if not availiable
#ifdef CXXRANDOM
  seed = rd();
  rg.seed(seed);

  // Warm-up
  uniform_real_distribution<double> uni_dist(0, 1);
  for(int i=0;i<10000;i++)
    uni_dist(rg);
#else
  if((fp = fopen("/dev/urandom","r")) != NULL){
    fread(&seed, sizeof(int), 1, fp);
    fclose(fp);
  } else {
    seed = time(NULL) + 10 * (global.myid + 1);
  }
  srand(seed);
#endif

  // Fetch all random seeds to allow for reproducibillity of the results
  if(global.myid==0) allseeds = new long[global.ncpu];
#ifdef WITHMPI
  MPI_Gather(&seed,1,MPI_LONG,allseeds,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif
  if(global.myid==0){
    cout << " -- Random seeds used by different CPUs: " << endl;
    for(int i=0;i<global.ncpu;i++) cout << " -- id = " << i << " seed = " << allseeds[i] << endl;
    cout << endl;
    delete[] allseeds;
  }

  // Read filepath and ramses output nr
  global.filepath = argc[1];
  global.filenum  = atoi(argc[2]);
  global.nfiles   = atoi(argc[3]);
  if(global.filepath == ".") global.filepath = getcwd(cwd, sizeof(cwd));
 
  // Allocate and initialize grid for density field (and later epsilon = F{delta}/|F{delta}|)
  global.allocate_grid(NGRID_HARDCODED);

  // Boxsize and rmax
  global.boxsize = atoi(argc[5]);

  // Make binning for l(r): npoints from rmin to rmax
  global.allocate_bins(atoi(argc[6]), 1.0/double(global.ngrid/2), 0.5);

  // Output filename
  global.outputprefix = argc[7];

  // Verbose
  if(global.myid==0){
    cout << "Initializing" << endl;
    cout << " -- Running with    " << global.ncpu                << " CPUs" << endl;
    cout << " -- RAMSES files    " << global.filepath            << " " << global.filenum << " " << global.nfiles << endl;
    cout << " -- Gridsize        " << global.ngrid               << endl;
    cout << " -- Nparticles      " << atoi(argc[4])              << endl;
    cout << " -- Rmin (Mpc/h)    " << global.boxsize*global.rmin << endl;
    cout << " -- Rmax (Mpc/h)    " << global.boxsize*global.rmax << endl;
    cout << " -- Box  (Mpc/h)    " << global.boxsize             << endl;
    cout << " -- Nbins           " << global.nbins               << endl;
    cout << " -- Outputprefix    " << global.outputprefix         << endl;
  }

  // Read and bin particles
  read_and_bin_particles(atoi(argc[4]), global.filepath, global.filenum, global.nfiles);

  // Precompute j0(2pi Sqrt(x)). If not called we compute it directly
#ifdef USEBESSELLOOKUP
  precompute_j0(NPOINTS_J0_PRECOMP, BESSELXMAX);
#endif

  // Allocate memory for pofk and zeta(r) binning. 
  // If not called it will not be computed!
  int kmax = int(NGRID_HARDCODED/2 * sqrt(3.0)) + 1;
  int ncorr_ptr = 50;
  global.allocate_pofk(kmax, ncorr_ptr);
}

//============================================================
// Read halo catalog (here for AHF)
//============================================================

void read_and_bin_halos(string filename){
  int nhalos, nparttot, ncolhalofile = 43;
  int ngrid, ngridtot;
  double tmp[ncolhalofile], oneoverboxsizekpc, norm;
  double avg, avg2;
  clock_t begin = clock();
  string line;

  // Initialize
  ngrid = global.ngrid;
  ngridtot = global.ngridtot;
  oneoverboxsizekpc = global.boxsize*1000.0;

  // Open file
  ifstream fp(filename.c_str());
  
  // Count number of lines in file
  nhalos = std::count(std::istreambuf_iterator<char>(fp),
                      std::istreambuf_iterator<char>(), '\n') - 1;
  
  // Verbose
  if(global.myid==0){
    cout << "Read and bin halos from " << filename << endl;
    cout << " -- Found " << nhalos << " halos" << endl;
  }

  // Rewind and add
  fp.clear();
  fp.seekg(0);
  getline(fp,line);

  // Read all halos 
  nparttot = 0;
  for(int i=0;i<nhalos;i++){
    double x[3], rvir, mass;
    int npart, ix[3], index;

    // Read line
    for(int j=0;j<ncolhalofile;j++) fp >> tmp[j];    

    // Extract values 
    mass  = tmp[3];
    npart = int(tmp[4]);
    x[0]  = tmp[5]  * oneoverboxsizekpc;
    x[1]  = tmp[6]  * oneoverboxsizekpc;
    x[2]  = tmp[7]  * oneoverboxsizekpc;
    rvir  = tmp[11] * oneoverboxsizekpc;

    // Add to grid (Add CIC smoothing?) 
    for(int j=0;j<3;j++){ 
      ix[j] = int(x[j]*ngrid);
      if(ix[j] >= ngrid){ 
        if(ix[j] == ngrid){
          ix[j] = 0;
        } else {
          // Error particle position not in [0,1] ( x >= 1 + 1/ngrid )
          cout << "Error in read particles. Index too large: ix = " << ix[j] << endl;
          myexit(6,"readpart");
        } 
      }
    }
    index = ix[0] + ngrid*(ix[1] + ngrid*ix[2]);
    if(index < 0 || index >= ngridtot) {
      cout << "Too low / high index: ind = " << index << endl;
      myexit(7,"readpartindex");
    }
    global.grid[index][0] += npart;
    nparttot += npart;
  }

  // Normalize and compute useful info
  norm = 1.0/double(nparttot) * ngridtot;
  avg = avg2 = 0.0;
  for(int i=0;i<ngridtot;i++){
    global.grid[i][0] = global.grid[i][0]*norm-1.0;
    global.grid[i][1] = 0.0;
    avg  += global.grid[i][0];
    avg2 += pow2(global.grid[i][0]);
  }
  avg /= double(ngridtot);
  avg2 = sqrt(avg2/double(ngridtot));

  // Verbose
  clock_t end = clock();
  if(global.myid==0){
    cout << " -- <Density> = " << avg << " <Density^2> = " << avg2 << endl;
    cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;
  }
}

//============================================================
// Read N-body data: RAMSES data
//============================================================

void read_and_bin_particles(int nparttot, string filepath, int filenum, int nfiles){
  double **buffer, norm, avg, avg2;
  int ngrid, ngridtot, nread, nbuffer;
  clock_t begin = clock();

  // Initialize
  ngrid = global.ngrid;
  ngridtot = global.ngridtot;
  nread = 0;

  // Verbose
  if(global.myid==0) cout << " -- Read RAMSES particle files" << endl;

  // Sanity check
  if(nparttot > pow3(1024) || nparttot <= pow3(32)){
    cout << "Do we really want to run with " << nparttot << " particles?" << endl;
    myexit(8,"nparttotoff");
  }

  // Read header of RAMSES files to get max particle number we need to allocate
  int npart_max = 0;
  for(int i=1;i<=nfiles;i++){
    string filename = filepath + "output_" + ramses_num_as_string(filenum) + "/part_" + ramses_num_as_string(filenum) + ".out" + ramses_num_as_string(i);
    FILE* fp;
    if((fp = fopen(filename.c_str(),"r")) == NULL){
      cout << "Error cannot read file..." << endl;
      myexit(9,"cannotreadfile");
    }
    int ncpu_loc, ndim_loc, npart_loc;
    ncpu_loc      = read_int(fp);
    ndim_loc      = read_int(fp);
    npart_loc     = read_int(fp);
    if(npart_loc > npart_max) npart_max = npart_loc;
    fclose(fp);
  }
  
  // Allocate read buffer 
  nbuffer = int( npart_max );
  buffer = new double*[3];
  for(int i=0;i<3;i++)
    buffer[i] = new double[nbuffer];

  // Loop over input-files and read one by one. nread keeps track of how many particles we have read so far.
  for(int i=1;i<=nfiles;i++){
    string filename = filepath + "output_" + ramses_num_as_string(filenum) + "/part_" + ramses_num_as_string(filenum) + ".out" + ramses_num_as_string(i);
    FILE* fp;

    // Check for error
    if((fp = fopen(filename.c_str(),"r")) == NULL){
      cout << "Error cannot read file..." << endl;
      myexit(10,"cannotreadfile");
    }

    // Verbose
    if(global.myid==0) cout << " ---- Reading from... " << filename << " ( " << nread << " ) " << endl;

    // Read RAMSES header
    int ncpu_loc, ndim_loc, npart_loc, localseed_loc[4], nstar_tot_loc, mstar_tot_loc[2], mstar_lost_loc[2], nsink_loc;
    ncpu_loc      = read_int(fp);
    ndim_loc      = read_int(fp);
    npart_loc     = read_int(fp);
    read_int_vec(fp, localseed_loc, 4);
    nstar_tot_loc = read_int(fp);
    read_int_vec(fp, mstar_tot_loc, 2);
    read_int_vec(fp, mstar_lost_loc, 2);
    nsink_loc     = read_int(fp);

    // Check if buffer is large enough!
    if(npart_loc > nbuffer){
      cout << "Error read buffer not large enough: " << nbuffer << " vs " << npart_loc << " Resizing it..." << endl;
      myexit(11,"buffersizeread");
    }

    // Read particle positions
    for(int j=0;j<ndim_loc;j++)
      read_double_vec(fp, buffer[j], npart_loc);
    fclose(fp);

    // Assign particle data to grid
    for(int p=0;p<npart_loc;p++){
      int ix[3], index;
      for(int j=0;j<3;j++){ 
        ix[j] = int(buffer[j][p]*ngrid);
        if(ix[j] >= ngrid){ 
          if(ix[j] == ngrid){
            ix[j] = 0;
          } else {
            // Error particle position not in [0,1] ( x >= 1 + 1/ngrid )
            cout << "Error in read particles. Index too large: ix = " << ix[j] << endl;
            myexit(12,"readpart");
          } 
        }
      }
      index = ix[0] + ngrid*(ix[1] + ngrid*ix[2]);
      if(index < 0 || index >= ngridtot) {
        cout << "Too low / high index: ind = " << index << " i = " << p << endl;
        myexit(13,"indexoffread");
      }
      global.grid[index][0] += 1.0;
    }

    // Update nread. NB: must be done after reading particle positions!
    nread += npart_loc;
  }

  // Check we got all particles
  if(nread != nparttot){
    cout << "Error did not get all particles " << nread << " != " << nparttot << endl; 
    myexit(14,"notallpart");
  }

  // Verbose
  if(global.myid==0) cout << " ---- Done reading " << nread << " particles" << endl;

  // Normalize grid to get delta_m and compute some useful info
  norm = 1.0/double(nparttot) * ngridtot;
  avg = avg2 = 0.0;
  for(int i=0;i<ngridtot;i++){
    global.grid[i][0] = global.grid[i][0]*norm-1.0;
    global.grid[i][1] = 0.0;
    avg  += global.grid[i][0];
    avg2 += pow2( global.grid[i][0]);
  }
  avg /= double(ngridtot);
  avg2 = sqrt(avg2/double(ngridtot));

  // Verbose
  clock_t end = clock();
  if(global.myid==0) {
    cout << " ---- <Density> = " << avg << " <Density^2> = " << avg2 << endl;
    cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;
  }

  // Cleanup 
  for(int i=0;i<3;i++)
    delete[] buffer[i];
  delete[] buffer;
}

// Calculate power and mode-norm from a given index in a fft grid
inline void add_to_pofk(int index, fftw_complex* grid){
    double power, window, arg, know, weight;
    int k[3], kind, ngrid = global.ngrid;

    // Calculate k-mode
    for(int j=0,n=1;j<3;j++,n *= ngrid){
      k[j] = index/n % ngrid;
      if(k[j] >= ngrid/2) k[j] -= ngrid;
    }
    know = sqrt(pow2(k[0])+pow2(k[1])+pow2(k[2]));
    kind = floor(know);

    // Calculate power (and deconvolve the window function DENSITYASSIGNMENT is 1 (NGP), 2 (CIC), 3 (TSC))
    window = 1.0;
    for(int j=0;j<3;j++){
      arg = pow2(k[j]/(2.0*ngrid));
      window *= j0_func(arg);
    }
    window = pow(window,DENSITYASSIGNMENT);
    power = ( pow2(grid[index][0]) + pow2(grid[index][1]) ) / pow2(window) / pow2(double(global.ngridtot));
    
    // Add to bin. Use linear interpolation between neighboring bins for k>3
    if( kind <= 3){  
      global.pofk[kind] += power;
      global.pofk_count[kind] += 1.0;
    } else if(kind < global.npofk-1) {
      weight = 1.0 + double(kind) - know;
      global.pofk[kind] += power * weight;
      global.pofk_count[kind] += weight;
      global.pofk[kind+1] += power * (1.0-weight);
      global.pofk_count[kind+1] += 1.0-weight;
    }
} 

//============================================================
// Compute the phase of the density field. At the same time
// compute the power-spectrum P(k) 
//
// returnreal = true => global.grid = F^{-1}[d_k/|d_k|] = eps(r)
// is the real-space whitened field
//============================================================

void compute_phase_and_power(bool returnreal){
  double norm,avg,avg2;
  int ngrid, ngridtot;
  fftw_complex *grid;
  fftw_plan p,q;
  clock_t begin = clock();

  // Initialize
  grid = global.grid;
  ngrid = global.ngrid;
  ngridtot = global.ngridtot; 

  // Verbose
  if(global.myid==0){
    cout << "Perform whitening of the density field n = " << ngrid << " " << ngridtot << endl;
    if(global.alloc_pofk) cout << " -- Will also compute P(k)" << endl;
  }

  // Calculate some useful quantities
  avg = avg2 = 0.0;
  for(int i=0;i<ngridtot;i++){
    avg += grid[i][0];
    avg2 += pow2(grid[i][0]);
  }
  avg /= double(ngridtot);
  avg2 = sqrt(avg2/double(ngridtot));

  // Verbose
  if(global.myid==0) cout << " -- Before whitening <field> = " << avg << " <field^2> = " << avg2 << endl;

  // Transform to fourier space...
  p = fftw_plan_dft_3d(ngrid, ngrid, ngrid, grid, grid, FFTW_FORWARD,  FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  for(int i=1;i<ngridtot;i++){
    // Bin up power-spectrum (because why not)
    if(global.alloc_pofk) add_to_pofk(i, grid);

    // Whitening...
    norm = sqrt(pow2(grid[i][0]) + pow2(grid[i][1]));
    grid[i][0] /= norm;
    grid[i][1] /= norm;
  }
  grid[0][0] = grid[0][1] = 0.0;

  // Normalize pofk and subtract shotnoise
  if(global.alloc_pofk){
    for(int i=0;i<global.npofk;i++)
      if(global.pofk_count[i] > 0.0) global.pofk[i] = global.pofk[i]/global.pofk_count[i] - 1.0/double(ngridtot);
  }

  // Option to transform back to real space or leave it as it is...
  if(returnreal){

    // Transform back to real space
    q = fftw_plan_dft_3d(ngrid, ngrid, ngrid, grid, grid, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(q);
    fftw_destroy_plan(q);

    // Calculate some useful quantities
    avg = avg2 = 0.0;
    norm = double(ngridtot);
    for(int i=0;i<ngridtot;i++){
      grid[i][0] = grid[i][0] * norm;
      avg += grid[i][0];
      avg2 += pow2(grid[i][0]);
    }
    avg /= double(ngridtot);
    avg2 = sqrt(avg2/double(ngridtot));

    // Verbose
    if(global.myid==0) cout << " -- After whitening <field> = " << avg << " <field^2> = " << avg2 << endl;
  }

  // Verbose
  clock_t end = clock();
  if(global.myid==0) cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;
}

//============================================================
// Spherical bessel-function [ NB: Sinc[ 2*pi * Sqrt[x] ] ]
// Precomputing or direct calculation. The first option is much faster
//============================================================

// Precompute Sinc[2pi Sqrt[x]] and store in global
void precompute_j0(int npoints, double xmax){

  // Verbose
  if(global.myid==0) cout << "\nPrecomputing Sinc(2pi Sqrt[x])" << endl;

  // Allocate memory and set global params
  global.allocate_precomp(npoints, double(npoints) / xmax);

  // Calculate
  for(int i=0;i<npoints;i++){
    double x = (i+0.5)/global.precomp_j0_fac;
    if(x>1e-3){
      global.precomp_j0[i] = sin(TWOPI * sqrt(x)) / (TWOPI * sqrt(x));
    } else {
      global.precomp_j0[i] = 1.0 - pow2(TWOPI)/6.0*x;
    }

    // Verbose
    if(global.myid==0 && (i % (npoints/10) == 0)) cout << " -- x =" << setw(10) << x << "  j0[2 * pi * Sqrt(x)] =  " << setw(12) << global.precomp_j0[i] <<  " Test lookup (0.0 = OK): " << j0_func(x) - global.precomp_j0[i] << endl;
  } 
}

#ifdef USEBESSELLOOKUPTBLE

// Lookup of precomputed j0
inline double j0_func(double x){
  int i = int(x * global.precomp_j0_fac);
  // Bounds check
  if( i >= global.precomp_j0_npoints || i < 0){ cout << "Error in j0_func (lookup): x = " << x << " i = " << i << endl; myexit(15,"j0lookup");} 
  return global.precomp_j0[i];
}

#else

// Expensive, but guaranteed accurate, j0-calculation
inline double j0_func(double x){
  double y = sqrt(x)*TWOPI;
  if(y<1e-3 && y>-1e-3)
    return 1.0 - y*y/6.0;
  return sin(y)/y;
}
#endif

//============================================================
// Compute the real space correlation function from a P(k) 
// obtained from a FT grid with gridsize ngrid (which
// determines rmin, rmax)
//============================================================

void compute_correlation_function(){
  double logkmin, dlogk, kmin, kmax, know, integral;
  double rmin, rmax, rnow, dlogr, logrmin, oldval;
  int nint, ngrid = global.ngrid; 
  clock_t begin = clock();

  // If P(k) is not computed then no point of doing this...
  if(!global.alloc_pofk) return;

  // Verbose
  if(global.myid==0) cout << "Calculate 2-point correlation function from P(k)" << endl;

  // Fill in any zero-values in P(k) array. Just for safety!
  oldval = global.pofk[1];
  for(int i=0;i<global.npofk;i++){
    if(global.pofk[i] == 0.0) global.pofk[i] = oldval;
    else oldval = global.pofk[i];
  }

  // Spline [k, P(k)] with k in units of 1/B0 and P(k) in units of B0^3
  Spline pofk_spline(global.pofk_k,global.pofk,global.npofk,1e30,1e30,1,"pofk_spline");

  // Integral limits [kmin/kmax is set according to what modes we have availiable]
  nint = 10000;
  kmin = TWOPI;
  kmax = TWOPI * double(ngrid/2);
  logkmin = log(kmin);
  dlogk = log(kmax/kmin)/double(nint);

  // Correlation function limits [NB: zeta(r) will not be accurate for all these r's]
  rmin  = 1.0/double(ngrid/2.0);
  rmax  = 0.5;
  logrmin = log(rmin);
  dlogr = log(rmax/rmin)/double(global.nzeta);

  // Loop over all r-values
  for(int k=0;k<global.nzeta;k++){
    rnow  = exp(logrmin + dlogr*k);

    // Integrate zeta(r) = 1/2pi^2 * Int_{kmin}^{kmax} dk/k [k^3 P(k)] sin(kr)/(kr)
    integral = 0.0;
    for(int i=0;i<nint;i++){
      know = exp(logkmin + dlogk*i);
      integral += pow3(know) * pofk_spline(know) * sin(know*rnow)/(know*rnow);
    }
    integral *= 2.0/pow2(TWOPI) * dlogk;

    // Store values
    global.rzeta[k] = rnow;    
    global.zeta[k]  = integral;
 
    // Verbose 
    if(global.myid==0)
      cout << " -- r = " << global.rzeta[k]*global.boxsize << " zeta = " << global.zeta[k] << endl; 
  }

  clock_t end = clock();
  if(global.myid==0) cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;
}

// Modulo function for translating integer wavenumber to grid-index
inline int mmod(int i){
  // Bounds check
  if(i>=NGRID_HARDCODED) {cout << "Error in mmod: " << i << endl; myexit(16,"mmod");}
  return i>=0 ? i : i + NGRID_HARDCODED;
}

//============================================================
// Calculate the line correlation function for a single 'r'
// r is in units of the boxsize
// kmax is the maximum integer wavenumber we probe
//============================================================

double calc_line_corr_single(double r){
  double integrand,r2,c1,c2,c3,s1,s2,s3,window,integrand_tot;
  int kmax, kmax2, kx_low_local, kx_high_local;
  int nchunk,extra,myid,index;
  fftw_complex *grid = global.grid;
  clock_t begin = clock();

  // Max integer wavenumber [kphys = kint * 2pi/B0 <= 2.0*pi/rphys]
  kmax      = int(1.0/r);
  r2        = r*r;
  integrand = 0.0;

  // Should not happen if global.rmin is set appropriately, but check anyway!
  if(kmax >= NGRID_HARDCODED/2){
    cout << " -- Warning changing kmax = " << kmax << " -> " << NGRID_HARDCODED/2 << endl; 
    kmax = NGRID_HARDCODED/2;
  }
  kmax2 = kmax*kmax;

  // Verbose
  if(global.myid==0){
    cout << " -- Calculate line correation function for r = " << r*global.boxsize << " Mpc/h" << endl;
    cout << " ---- kmax = " << kmax << " ngrid/2 = " << NGRID_HARDCODED/2 << endl;
  }

  // If kmax is too large use Monte-Carlo integration instead
  // since otherwise it will take forever to compute the sums
  if(kmax > 25) return monte_carlo_integrate_line_corr(r);

#ifdef WITHMPI

  // Calculate MPI workload
  nchunk  = (2*kmax+1)/global.ncpu;
  extra   = (2*kmax+1)-nchunk*global.ncpu;
  myid    = global.myid;

  // Divide up MPI tasks
  if(2*kmax+1 <= global.ncpu){
    kx_low_local  = -kmax + myid;
    kx_high_local = kx_low_local + 1;

    // If more cpu's than tasks let some cpu's rest...
    if(myid>=2*kmax){
      // For debugging:
      // cout << "---- Cpu: " << global.myid << " laying low" << endl; 
      goto endloop; 
    }
  } else {
    // The first few CPUs take some extra work
    if(myid<extra){
      kx_low_local  = -kmax + (nchunk+1)*myid;
      kx_high_local = -kmax + (nchunk+1)*(myid+1);
    } else {
      kx_low_local  = -kmax + extra + (nchunk)*myid;
      kx_high_local = -kmax + extra + (nchunk)*(myid+1);
    }
  }
  if(kx_high_local > kmax) kx_high_local = kmax;

  // For debugging:
  // cout << "---- Cpu: " << global.myid << " Doing " << kx_low_local << " " << kx_high_local-1 << endl;
  
#else
  kx_low_local = -kmax;
  kx_high_local = kmax;
#endif

  // Loop over all k-vectors and integrate up...
  for(int kx=kx_low_local;kx<kx_high_local;kx++){
    if(global.myid==0) cout << " ---- " << global.myid << " : " <<  kx-kx_low_local+1 <<  " / " << kx_high_local-kx_low_local << endl;
    int kylim = int(sqrt(kmax*kmax-kx*kx)+0.5); if(kylim > kmax) kylim = kmax;
    for(int ky=-kylim;ky<kylim;ky++){
      int kzlim = int(sqrt(kmax*kmax-kx*kx-ky*ky)+0.5); if(kzlim > kmax) kzlim = kmax;
      for(int kz=-kzlim;kz<kzlim;kz++){

        // Check if we are within the sphere |k| <= kmax
        if(pow2(kx) + pow2(ky) + pow2(kz) <= kmax2){
          
          // Compute c1, s1
          index = mmod(kx)  + NGRID_HARDCODED*(mmod(ky) + NGRID_HARDCODED*mmod(kz));
          c1 = grid[index][0];
          s1 = grid[index][1];

          // Loop over q-vectors
          for(int qx=-kmax;qx<kmax;qx++){
            int qylim = int(sqrt(kmax*kmax-qx*qx)+0.5); if(qylim > kmax) qylim = kmax;
            for(int qy=-qylim;qy<qylim;qy++){
              int qzlim = int(sqrt(kmax*kmax-qx*qx-qy*qy)+0.5); if(qzlim > kmax) qzlim = kmax;
              for(int qz=-qzlim;qz<qzlim;qz++){

                // Check if we are within the sphere |q| <= kmax and |k+q| <= kmax
                if(pow2(qx) + pow2(qy) + pow2(qz) <= kmax2 && pow2(kx+qx) + pow2(ky+qy) + pow2(kz+qz) <= kmax2){

                  // Compute c2, s2
                  index = mmod(qx)  + NGRID_HARDCODED*(mmod(qy) + NGRID_HARDCODED*mmod(qz));
                  c2 = grid[index][0];
                  s2 = grid[index][1];
               
                  // Compute c2, s3 
                  index = mmod(-qx-kx)  + NGRID_HARDCODED*(mmod(-qy-ky) + NGRID_HARDCODED*mmod(-qz-kz));
                  c3 = grid[index][0];
                  s3 = grid[index][1];
  
                  // Compute window-function
                  window = j0_func( r2 * (pow2(kx-qx) + pow2(ky-qy) + pow2(kz-qz)) );

                  // Sum up the integrand
                  integrand += ( c1*(c2*c3 - s2*s3) - s1*(s2*c3 + s3*c2) ) * window;
                }
              }
            }
          }
        }
      }
    }
  }

#ifdef WITHMPI
  endloop:

  // For debugging:
  // cout << " ---- " << global.myid << " is done  r = " << r*global.boxsize  << " (Mpc/h)  l_local = " << integrand << endl; 

  MPI_Allreduce(&integrand, &integrand_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  integrand = integrand_tot;
#endif

  // Normalize
  integrand *= pow(r,4.5);

  clock_t end = clock();
  if(global.myid==0) cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;

  return integrand;
}

//============================================================
// Calculate l(r) for a range of 'r' using Monte-Carlo
// integration with [n x ncpu] samples.
//
// [delta l(r) / l(r)] < epsilon is the convergence criterion.
// Does a maximum of [n x ncpu] samples
//============================================================

double monte_carlo_integrate_line_corr(double r){
  double c1,c2,c3,s1,s2,s3,window,avg,var,integrand,integrand2, curvalue, sigma, epsilon;
  double r2, volratio, norm;
  int kmax,twokmax,kmax2,kx,ky,kz,qx,qy,qz,index;
  long long nmaxsample, n_between_checks;
  fftw_complex *grid;
  clock_t begin = clock();

  // Verbose
  if(global.myid==0) cout << " ---- Using Monte-Carlo integration to calculate l(r). Maxsamples = " << double(NSAMPLEMAX_HARDCODED) << " Check convergence every " << double(NSAMPLECONVCHECK) << " samples" << endl;

  // Initialize
  grid       = global.grid;
  nmaxsample = NSAMPLEMAX_HARDCODED;
  epsilon    = EPSILON_MC_INTEGRATE;
  n_between_checks = NSAMPLECONVCHECK/global.ncpu;
  kmax       = int(1.0/r);
  kmax2      = pow2(kmax);
  twokmax    = 2*kmax;
  r2         = r*r;

  // Random number generator
#ifdef CXXRANDOM
  uniform_int_distribution<int> uni_dist(-kmax, kmax-1);
#endif

  // The volume inside the region |q|,|k|,|k-q|<= 1 relative to the volume of the cube |k_i|,|q_i| <= 1  times  (kmax)^6
  volratio = 5.0 * pow2(TWOPI) / 3.0 / pow(2.0,3);
  volratio *= pow(kmax,6);

  // Reset timing
  estimate_progress_monte_carlo(0.0, 0.0, true);

  // Do brute-force Monte-Carlo integration
  integrand  = integrand2 = 0.0;
  for(long long i=0;i<nmaxsample;i++){

    // Generate random k-vector
    while(1){
#ifdef CXXRANDOM
      kx = uni_dist(rg);
      ky = uni_dist(rg);
      kz = uni_dist(rg);
#else
      kx = (rand() % twokmax)-kmax;
      ky = (rand() % twokmax)-kmax;
      kz = (rand() % twokmax)-kmax;
#endif
      if(kx*kx + ky*ky + kz*kz < kmax2) break;
    }

    // Generate random q-vector
    while(1){
#ifdef CXXRANDOM
      qx = uni_dist(rg);
      qy = uni_dist(rg);
      qz = uni_dist(rg);
#else
      qx = (rand() % twokmax)-kmax;
      qy = (rand() % twokmax)-kmax;
      qz = (rand() % twokmax)-kmax;
#endif
      if(qx*qx + qy*qy + qz*qz < kmax2 && pow2(qx+kx) + pow2(qy+ky) + pow2(qz+kz) < kmax2 ) break;
    }
      
    // The e(k) term
    index = mmod(kx) + NGRID_HARDCODED*(mmod(ky) + NGRID_HARDCODED*mmod(kz));
    c1 = grid[index][0];
    s1 = grid[index][1];

    // The e(q) term
    index = mmod(qx) + NGRID_HARDCODED*(mmod(qy) + NGRID_HARDCODED*mmod(qz));
    c2 = grid[index][0];
    s2 = grid[index][1];

    // The e(-k-q) term
    index = mmod(-qx-kx) + NGRID_HARDCODED*(mmod(-qy-ky) + NGRID_HARDCODED*mmod(-qz-kz));
    c3 = grid[index][0];
    s3 = grid[index][1];

    // Window function
    window = j0_func( r2 * (pow2(kx-qx) + pow2(ky-qy) + pow2(kz-qz)) );

    // Integrate up [e(k)e(q)e(-k-q) * j0] and the L2 norm of the same integrand
    curvalue  = ( c1*(c2*c3 - s2*s3) - s1*(s2*c3 + s3*c2) ) * window;
    integrand  += curvalue;
    integrand2 += pow2(curvalue);

    // Calculate average and variance over all CPU's
    if(i % n_between_checks == 0 && i>0){
      // Make comm arrays
      double avg_local[global.ncpu], avg_all[global.ncpu], var_local[global.ncpu], var_all[global.ncpu];

      // Total number of samples done so far
      double nsamples = double(i+1) * global.ncpu;

      // Initialize. Store local value in [myid]'th position
      for(int j=0;j<global.ncpu;j++) 
        avg_local[j] = var_local[j] = 0.0;
      avg_local[global.myid] = integrand;
      var_local[global.myid] = integrand2;

#ifdef WITHMPI
      // Collect samples from all the CPU's and copy back to _local arrays
      MPI_Allreduce(avg_local, avg_all, global.ncpu, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(var_local, var_all, global.ncpu, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for(int j=0;j<global.ncpu;j++){ 
        avg_local[j] = avg_all[j];
        var_local[j] = var_all[j];
      }
#endif

      // Calculate average over all samples
      avg = 0.0;
      for(int j=0;j<global.ncpu;j++) 
        avg += avg_local[j];
      avg /= nsamples;

      // Calculate variance
      var = 0.0;
      for(int j=0;j<global.ncpu;j++)
        var += var_local[j];
      var /= nsamples;  
      var  = fabs(var - pow2(avg));

      // Estimate for error
      sigma = sqrt(var) / sqrt(nsamples);

      // Set the correct normalization
      norm = pow(r,4.5) * volratio;
      avg   *= norm;
      sigma *= norm;
      
      // Verbose
      if(global.myid == 0){ 
        cout << " ---- r = " << r*global.boxsize << " <I> = " << avg;
        cout << " Estimated Progress: " << estimate_progress_monte_carlo(log(nsamples), log(fabs(sigma/avg))) << "%" << endl;
        cout << "    Done with nsamples = 10^{" << log10(nsamples) << "} <I> = " << avg << " Sigma(I) / <I> = " << fabs(sigma/avg) << endl;
      }

      // Check for convergence
      if(fabs(sigma/avg) < epsilon){
        if(global.myid==0)  cout << " ---- Converged after 10^{" << log10(nsamples) << "} samples" << endl;
        integrand = avg;
        break;
      }
    }

    // If we did not get convergence... Should increase maxsteps!
    if(i+1==nmaxsample){
      if(global.myid==0)  cout << " ---- WARNING: We have done a total of NSAMPLEMAX x NCPU without convergence so defining convergence" << endl;
#ifdef WITHMPI
      MPI_Allreduce(&integrand, &integrand2, global.ncpu, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      integrand = integrand2/double(global.ncpu);
#endif
      norm = 1.0/double(i+1) * pow(r,4.5) * volratio;
      integrand *= norm;
    }
  }

  // Verbose 
  clock_t end = clock();
  if(global.myid==0) cout << " ---- Elapsed time (sec): " << double(end - begin) / CLOCKS_PER_SEC << endl;

  return integrand;
}

//============================================================
// Calculate progress in MC integration...
// We expect log(sigma) = a  - 1/2 log(nsamples) so perform running 
// least square fit to log(sigma) = a + b log(nsamples) and output
// estimated percentage of work done.
//
// Variables: [y_i = log(eps_i)  x_i = log(nsamples_i)]
//============================================================

double estimate_progress_monte_carlo(double x, double y, bool reset){
  static double sumx  = 0.0;
  static double sumy  = 0.0;
  static double sumxx = 0.0;
  static double sumyy = 0.0;
  static double sumxy = 0.0;
  static double n     = 0.0;
  double dom, a, b;

  // Reset counters...
  if(reset) {
    sumx = sumy = sumxx = sumyy = sumxy = n = 0.0;
    return 0.0;
  }

  // Update counters
  sumx  += x;
  sumy  += y;
  sumxx += x*x;
  sumyy += y*y;
  sumxy += x*y;
  n     += 1.0;

  if(n==1.0) return 0.0;

  // Least square fit
  dom = n * sumxx - sumx * sumx;
  a = (sumy * sumxx - sumx * sumxy) / dom;
  b = (n * sumxy - sumx * sumy) / dom;

  // Calculate percentage n / n_expected * 100 
  return exp(x - (log(EPSILON_MC_INTEGRATE) - a)/b) * 100.0;
}

//============================================================
// Calculate l(r) for a range of 'r' values and output to file
//============================================================

void calc_line_corr_all(){
  double rmin, rmax;
  int nbins;
  string filename;
  fftw_complex *grid;

  // Initialize
  grid = global.grid;
  nbins = global.nbins;
  rmin = global.rmin;
  rmax = global.rmax;

  // Verbose
  if(global.myid==0) cout << "Calculate line-correlation function from " << rmin*global.boxsize << " ->  " << rmax*global.boxsize << " (Mpc/h) with " << nbins << " points." << endl;   

  // Calculate it up starting with large 'r'
  for(int i=nbins-1;i>=0;i--){
    global.r[i] = exp(log(rmin) + log(rmax/rmin)*(i+0.5)/double(nbins-1));
    global.l[i] = calc_line_corr_single(global.r[i]);

    // Verbose
    if(global.myid==0) cout << " -- r =  " << global.r[i]*global.boxsize << "   l(r) = " << global.l[i] << "\n" << endl;
  }
}

//============================================================
// Output to file
//============================================================

void output(bool output_l, bool output_pofk, bool output_zeta){
  // Output l(r) to file
  if(global.myid==0 && output_l){
    string filename = global.outputprefix + "_line.txt";
    cout << " -- Output l(r) to " << filename << endl; 
    ofstream fp(filename.c_str());
    for(int i=0;i<global.nbins;i++)
      fp << global.r[i]*global.boxsize << "  " << global.l[i] << endl;
    fp.close();
  }

  // Output P(k) to file
  if(global.myid==0 && output_pofk && global.alloc_pofk){
    string filename = global.outputprefix + "_pofk.txt";
    cout << " -- Output P(k) to " << filename << endl; 
    ofstream fp(filename.c_str());
    for(int i=0;i<global.npofk;i++)
      // Only output non-empty bins and only up to the Nyquist frequency
      if(global.pofk[i] != 0.0 && global.pofk_k[i]/TWOPI <= global.ngrid/2) 
        fp << global.pofk_k[i]/global.boxsize << "  " << global.pofk[i] * pow3(global.boxsize) << endl;
    fp.close();
  }

  // Output zeta(r) to file
  if(global.myid==0 && output_zeta && global.alloc_pofk){
    string filename = global.outputprefix + "_zeta.txt";
    cout << " -- Output zeta(r) to " << filename << endl; 
    ofstream fp(filename.c_str());
    for(int i=0;i<global.nzeta;i++)
      fp << global.rzeta[i]*global.boxsize << "  " << global.zeta[i] << endl;
    fp.close();
  }
}

//============================================================
// Main is main...
//============================================================

int main(int argv, char **argc){
  clock_t begin = clock();

  // Initialize
  init(argv,argc);
  if(global.myid==0) cout << endl;

  // Perform whitening and calculate P(k). After this the grid contains the phase the density field
  compute_phase_and_power(false);
  output(false,true,false);
  if(global.myid==0) cout << endl;

  // Compute the real space correlation function by integrating up P(k) 
  compute_correlation_function();
  output(false,false,true);
  if(global.myid==0) cout << endl;

  // Calculate the line correlation function
  calc_line_corr_all();
  output(true,false,false);
  if(global.myid==0) cout << endl;

#ifdef WITHMPI
  MPI_Finalize();
#endif

  // Verbose
  clock_t end = clock();
  if(global.myid==0) cout << "Finished... Whole run took " << double(end - begin) / CLOCKS_PER_SEC << " sec" << endl;
}

