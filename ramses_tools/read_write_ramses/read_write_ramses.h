#ifndef READWRITERAMSES_HEADER
#define READWRITERAMSES_HEADER
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <iomanip>

#define NDIM 3      // Number of dimensions in simulation
#define VELOCITY    // If not defined we skip reading v-data from files

class Particle{
  public:
    double x[NDIM];
#ifdef VELOCITY
    double v[NDIM];
#endif
    Particle(){}
    Particle(double *_x) {
      std::memcpy(x, _x, NDIM * sizeof(double));
    }
#ifdef VELOCITY
    Particle(double *_x, double *_v) {
      std::memcpy(x, _x, NDIM * sizeof(double));
      std::memcpy(v, _v, NDIM * sizeof(double));
    }
#endif
    ~Particle(){}
};

class ReadRamsesSnapshot{

  //===========================================
  // 
  // Read (and write) files related to RAMSES 
  // snapshots
  //
  // So far implemented read of:
  //  info-file
  //  particle-file -> reads and stores particle
  //                   data to Particle *p
  //
  // AMR and HYDRO files not implemented yet
  //
  // Implemented write of:
  //  ic_deltab IC file
  //  RAMSES data as ASCII
  //
  //===========================================

  private:
    
    static const int MAXLEVEL = 50;
 
    // File description
    std::string filepath;
    int outputnr;

    // Data obtained from info-file
    int ncpu;
    int npart;
    int levelmin;
    int levelmax;
    int ngridmax;
    int nstep_coarse;
    double aexp;
    double time;
    double boxlen;
    double h0;
    double omega_m;
    double omega_l;
    double omega_k;
    double omega_b;
    double unit_l;
    double unit_d;
    double unit_t;
    double boxlen_ini;
    std::vector<int> npart_file;

    // Data obtained from particle-files
    Particle *p;

    // Book-keeping variables
    bool infofileread;
    bool partfilesread;
    int npart_read;
    int nsink_read;

  public:

    ReadRamsesSnapshot() : 
      filepath(""),
      outputnr(0),
      p(NULL),
      infofileread(false),
      partfilesread(false),
      npart_read(0),
      nsink_read(0) {}

    ReadRamsesSnapshot(std::string _filepath, int _outputnr) : 
      filepath(_filepath), 
      outputnr(_outputnr), 
      p(NULL),
      infofileread(false),
      partfilesread(false),
      npart_read(0),
      nsink_read(0) {}

    ReadRamsesSnapshot(ReadRamsesSnapshot& rhs){
      std::cout << "Object not meant to be copied" << std::endl;
      exit(0);
    }
    
    ReadRamsesSnapshot operator=(const ReadRamsesSnapshot& rhs){
      std::cout << "Object not meant to be assigned" << std::endl;
      exit(0);
    }

    ~ReadRamsesSnapshot(){
      if(p != NULL)
        delete[] p;
    }

    //====================================================
    // Read all the particles and store pos/vel if wanted
    //====================================================

    void read_particles(){
      if(!infofileread) read_info();

      std::cout << std::endl; 
      std::cout << "==================================" << std::endl;
      std::cout << "Read particle files:              " << std::endl;
      std::cout << "==================================" << std::endl;

      std::cout << "Allocate memory for particles... ";  
      p = new Particle[npart];
      std::cout << "Done..." << std::endl;

      // Loop and read all particle files
      npart_file = std::vector<int>(ncpu, 0);
      for(int i = 0; i < ncpu; i++)
        read_particle_file(i);

      std::cout << "Done reading n = " << npart_read << " particles" << std::endl;
    }

    //====================================================
    // Integer to ramses string
    //====================================================

    std::string int_to_ramses_string(int i){
      std::stringstream rnum;
      rnum << std::setfill('0') << std::setw(5) << i;
      return rnum.str();
    }

    //====================================================
    // Read binary methods. The skips's are due to
    // files to be read are written by fortran code
    //====================================================

    int read_int(FILE* fp){
      int tmp, skip;
      fread(&skip, sizeof(int), 1, fp);
      fread(&tmp,  sizeof(int), 1, fp);
      fread(&skip, sizeof(int), 1, fp);
      return tmp;
    }

    void read_int_vec(FILE* fp, int *buffer, int n){
      int skip;
      fread(&skip,  sizeof(int), 1, fp);
      fread(buffer, sizeof(int), n, fp);
      fread(&skip,  sizeof(int), 1, fp);
    }

    void read_double_vec(FILE* fp, double *buffer, int n){
      int skip;
      fread(&skip,  sizeof(int),    1, fp);
      fread(buffer, sizeof(double), n, fp);
      fread(&skip,  sizeof(int),    1, fp);
    }

    //====================================================
    // Read a ramses info file
    //====================================================

    void read_info(){
      int ndim_loc;
      std::string numbers = int_to_ramses_string(outputnr);
      std::string infofile = filepath + "/output_" + numbers + "/info_" + numbers + ".txt";
      FILE *fp;

      if( (fp = fopen(infofile.c_str(), "r")) == NULL){
        std::cout << "Error opening info file " << infofile << std::endl;
        exit(0);
      }

      // Read the info-file
      fscanf(fp, "ncpu        =  %d\n", &ncpu);
      fscanf(fp, "ndim        =  %d\n", &ndim_loc);
      fscanf(fp, "levelmin    =  %d\n", &levelmin);
      fscanf(fp, "levelmax    =  %d\n", &levelmax);
      fscanf(fp, "ngridmax    =  %d\n", &ngridmax);
      fscanf(fp, "nstep_coarse=  %d\n", &nstep_coarse);
      fscanf(fp, "\n");
      fscanf(fp, "boxlen      =  %lf\n",&boxlen);
      fscanf(fp, "time        =  %lf\n",&time);
      fscanf(fp, "aexp        =  %lf\n",&aexp);
      fscanf(fp, "H0          =  %lf\n",&h0);
      fscanf(fp, "omega_m     =  %lf\n",&omega_m);
      fscanf(fp, "omega_l     =  %lf\n",&omega_l);
      fscanf(fp, "omega_k     =  %lf\n",&omega_k);
      fscanf(fp, "omega_b     =  %lf\n",&omega_b);
      fscanf(fp, "unit_l      =  %lf\n",&unit_l);
      fscanf(fp, "unit_d      =  %lf\n",&unit_d);
      fscanf(fp, "unit_t      =  %lf\n",&unit_t);
      fclose(fp);

      // Sanity check
      if(ndim_loc != NDIM){
        std::cout << "Error: Incompatible ndim = " << ndim_loc << " != " << NDIM << ". Change NDIM, recompile and run again!" << std::endl;
        exit(0);
      }

      // Calculate boxsize
      boxlen_ini = unit_l*h0/100.0/aexp/3.08e24;

      // Calculate number of particles [n = (2^levelmin) ^ ndim]
      npart = 1;
      for(int j = 0; j < NDIM; j++)
        npart = npart << levelmin;

      // Verbose
      std::cout << std::endl; 
      std::cout << "==================================" << std::endl;
      std::cout << "Read info-file:                   " << std::endl;
      std::cout << "==================================" << std::endl;
      std::cout << "Infofile     = "      << infofile   << std::endl;
      std::cout << "Box (Mpc/h)  = "      << boxlen_ini << std::endl;
      std::cout << "ncpu         = "      << ncpu       << std::endl;
      std::cout << "npart        = "      << npart      << std::endl;
      std::cout << "ndim         = "      << NDIM       << std::endl;
      std::cout << "aexp         = "      << aexp       << std::endl;
      std::cout << "H0           = "      << h0         << std::endl;
      std::cout << "omega_m      = "      << omega_m    << std::endl;

      infofileread = true;
    }

    //====================================================
    // Store the positions 
    //====================================================

    void store_positions(double* pos, const int dim, const int npart){
      Particle *pp = &p[npart_read];
      for(int i = 0; i < npart; i++)
        pp[i].x[dim] = pos[i];
    }

#ifdef VELOCITY
    //====================================================
    // Store the velocities
    //====================================================
    void store_velocity(double *vel, const int dim, const int npart){
      Particle *pp = &p[npart_read];
      double velfac = 100.0 * boxlen_ini / aexp;
      for(int i = 0; i < npart; i++)
        pp[i].v[dim] = vel[i] * velfac;
    }
#endif

    //====================================================
    // Read a single particle file
    //====================================================

    void read_particle_file(const int i){
      std::string numberfolder = int_to_ramses_string(outputnr);
      std::string numberfile = int_to_ramses_string(i+1);
      std::string partfile = filepath + "/output_" + numberfolder + "/part_" + numberfolder + ".out" + numberfile;
      double *buffer;
      FILE *fp;

      // Local variables used to read into
      int ncpu_loc, ndim_loc, npart_loc, localseed_loc[4], nstar_tot_loc, mstar_tot_loc[2], mstar_lost_loc[2], nsink_loc;

      std::cout << "Reading " << partfile << std::endl;

      // Open file
      if( (fp = fopen(partfile.c_str(), "r")) == NULL){
        std::cout << "Error opening particle file " << partfile << std::endl;
        exit(0);
      }

      // Read header
      ncpu_loc      = read_int(fp);
      ndim_loc      = read_int(fp);
      npart_loc     = read_int(fp);
      read_int_vec(fp, localseed_loc, 4);
      nstar_tot_loc = read_int(fp);
      read_int_vec(fp, mstar_tot_loc, 2);
      read_int_vec(fp, mstar_lost_loc, 2);
      nsink_loc     = read_int(fp);

      // Sanity check
      if(ndim_loc != NDIM){
        std::cout << "Error: Incompatible ndim = " << ndim_loc << " != " << NDIM << ". Change NDIM, recompile and run again!" << std::endl;
        exit(0);
      }

      // Store npcu globally
      if(!infofileread)
        ncpu = ncpu_loc;

      npart_file[i] = npart_loc;

      // Allocate memory for buffer
      buffer = new double[npart_loc];

      // Verbose
      std::cout << "npart = " << npart_loc << std::endl;

      // Read positions:
      for(int j = 0; j < ndim_loc; j++){
        read_double_vec(fp, buffer, npart_loc);
        store_positions(buffer, j, npart_loc);
      }

#ifdef VELOCITY
      // Read velocities:
      for(int j=0;j<ndim_loc;j++){
        read_double_vec(fp, buffer, npart_loc);
        store_velocity(buffer, j, npart_loc);
      }
#endif

      // Update global variables
      nsink_read += nsink_loc;
      npart_read += npart_loc;

      fclose(fp);
      delete[] buffer;
    } 

    //====================================================
    // Write a ramses ascii file. Assumes all the particle
    // pos and vel are already defined. Assumes the
    // mass is the same for all particles.
    //
    // Format:
    //
    // npart
    // [x1 y1 z1] [vx1 vy1 vz1] [mass]
    // ...
    // [xn yn zn] [vxn vyn vzn] [mass]
    //
    // If !writevel or VELOCITY is not defined then v-block
    // is ommited
    // If !writemass then mass-block is ommited
    //====================================================

    void write_ramses_ascii(const std::string filename, bool writemass = true, bool writevel = true){  
      float mass = 1.0/float(npart);
      std::ofstream fp;
      fp.open(filename.c_str());

      if(!fp){
        std::cout << "Error opening " << filename << std::endl;
        exit(0);
      }

      fp << npart << std::endl;
      for(int i = 0; i < npart; i++){
        for(int j = 0; j < NDIM; j++)
          fp << p[i].x[j]-0.5 << "  ";
#ifdef VELOCITY
        if(writevel)
          for(int j = 0; j < NDIM; j++)
            fp << p[i].v[j]-0.5 << "  ";
#endif
        if(writemass)
          fp << mass;
        fp << std::endl;
      }
    }

    //====================================================
    // Write a ramses ic_deltab file for a given cosmology
    // and levelmin. The variables used here are completely
    // independent from the onee in the header of this
    // class apart from the global define NDIM which is used
    //====================================================

    void write_icdeltab(const std::string filename, const float _astart, const float _omega_m, const float _omega_l, 
                             const float _boxlen_ini, const float _h0, const int _levelmin){
      FILE *fp;
      int tmp, n1, n2, n3;
      float xoff1, xoff2, xoff3, dx;

      // Assume zero offset
      xoff1 = 0.0;
      xoff2 = 0.0;
      xoff3 = 0.0;

      // Assumes same in all directions
      n1 = 1 << _levelmin;
      n2 = 1 << _levelmin;
      n3 = 1 << _levelmin;

      // dx in Mpc (CHECK THIS)
      dx  = _boxlen_ini / float(n1) * 100.0/_h0;    

      // Number of floats and ints we write
      tmp = 8*sizeof(float) + 3*sizeof(int);

      // Open file
      if( (fp = fopen(filename.c_str(), "w")) == NULL){
        std::cout << "Error opening " << filename << std::endl;
        exit(0);
      }

      // Write file
      fwrite(&tmp,      sizeof(int),   1, fp);
      fwrite(&n1,       sizeof(int),   1, fp);
      if(NDIM>1)
        fwrite(&n2,     sizeof(int),   1, fp);
      if(NDIM>2)
        fwrite(&n3,     sizeof(int),   1, fp);
      fwrite(&dx,       sizeof(float), 1, fp);
      fwrite(&xoff1,    sizeof(float), 1, fp);
      if(NDIM>1)
        fwrite(&xoff2,  sizeof(float), 1, fp);
      if(NDIM>2)
        fwrite(&xoff3,  sizeof(float), 1, fp);
      fwrite(&_astart,  sizeof(float), 1, fp);
      fwrite(&_omega_m, sizeof(float), 1, fp);
      fwrite(&_omega_l, sizeof(float), 1, fp);
      fwrite(&_h0,      sizeof(float), 1, fp);
      fwrite(&tmp,      sizeof(int),   1, fp);
      fclose(fp);
    }
};

#endif

