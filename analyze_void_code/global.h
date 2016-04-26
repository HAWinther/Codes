#ifndef GLOBAL_H_INC
#define GLOBAL_H_INC
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct Simulation {

  //====================================================
  //
  // Contains general informaton about the simulation
  //
  //====================================================

  int nthreads;		          // For OpenMP parallelisation

  double boxsize;           // Boxsize in Mpc/h
  int np;			              // Number of particles in simulation
  int nv;			              // Number of voids

  // Input files
  char partfile[256];       // Particle file  [/path/to/gadget.]
  char volfile[256];        // Volume file    [Vol.dat]
  char voidfile[256];       // ZobovVoid list [Voidlist.dat]
  char partinzonefile[256]; // ZobovP file    [PartInZone.dat]
  char zoneinvoidfile[256]; // ZobovZ file    [ZoneInVoid.dat]
  char adjfile[256];        // Adj file       [Adj.dat]
  char voidtreefile[256];   // Voidtree       [voidtree.dat]
                            // Voidtree file does not need to exist
} global;
typedef struct Simulation SIM;

struct ZobovVoid {

  //====================================================
  //
  // Contains all information about a ZobovVoid/Zone
  // Positions are in units of the boxsize
  //
  //====================================================

  int used;       		      // Book-keeping variable to determine if void has been used
  int isinvoidstack;  		  // Book-keeping variable to determine if the void is part of the stack 

  int n;				            // Index of the void in [allvoid] 
  int	CoreParticle;		      // Particle number of coreparticle in [allpart] 
  double CoreDensity;		    // Density of coreparticle in units of mean 
  int nCentralZone;		      // Number of particles in central zone of the void 
  double VolumeCentalZone;  // Volume of central zone of void in units of mean density particle 
  double VoidDensContrast;  // Density contrast of void 
  double VoidProb;		      // Probabillity of void due to Poisson noise
  double VoidVolume;		    // Volume of the void, in units of the volume occupied by a mean-density particle 
  double reff;			        // Effective radius of the void 
  double x,y,z;			        // Position of the void 

  // Zones 
  int nZones;			          // Number of zones in the void 
  int *zonelist;			      // Inception list: List of voids within the void 
  // Note that the void itself is included in zonelist[0] 

  // Particles 
  int np;				            // Total number of particles in void as said by voidfile (includes subvoids) 
  int npl;			            // Number of particles in this void (not including subvoids) as given by particle list 
  int *particlelist;		    // List of particles in the void 

  int children;			        // Number of subvoids 
  int ischild;			        // Book-keeping variable to determine if the void is a subvoid itself

  // Void tree
  struct ZobovVoid **childlist; // Pointers to subvoids in the tree 
  struct ZobovVoid *Down;   // Pointers to voids in the void tree
  struct ZobovVoid *Up;     // NULL means there is no void below/up/etc.  
  struct ZobovVoid *Next;
  struct ZobovVoid *Prev;
  
  // Other
  double InertiaTensor[6];  // Moment of inertia tensor for all particles in void

} *allvoid, rootvoid;
typedef struct ZobovVoid VOID;

struct SingleParticle {

  //====================================================
  //
  // Contains all information about a single particles
  // x,y,z are in units of the boxsize, i.e. in [0,1]
  // Volumes are in units of the boxsize^3
  // Velocities are in units of km/s
  //
  //====================================================

  int n;                    // Index of particle in array 
  double x, y, z;	          // Position of particle	
  double volume;	          // Volume of particle	

#ifdef ADJS
  int nadj;	                // Number of adjacent particles
  int *adj;	                // List of the index (in allpart) to adjecent particles 
#endif

#ifdef VELOCITY
  double vx, vy, vz;        // Velocity of particle 
#endif
} *allpart;
typedef struct SingleParticle PARTICLE;

typedef int (*comp)(const void*,const void*);

struct treecalc {

  //====================================================
  // Thread-private book-keeping struct for computing void
  // quantities like the barycenter etc. which is 
  // done via recursion over sub-sub-sub-...-voids
  //====================================================

  // Barycenter
  double x, y, z;

  // Pos of central part for void in question
  double xcp, ycp, zcp;

  // Volume and particles
  double volume;
  int npart;

  // Book-keeping array for all voids [which void is used]
  // Must be allocated/deallocated in OMP region if needed
  int *used;

  // Moment of inertia
  double Ixx, Iyy, Izz;
  double Ixy, Iyz, Izx;
} q;
typedef struct treecalc BOOKKEEP;

#ifdef USEOPENMP
#pragma omp threadprivate(q)
#endif

#endif
