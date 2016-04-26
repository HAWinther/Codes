#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "init.h"
#include "analysis.h"

//=================================================================================
// 
// Code to read in particle data and output from Zobov and organize everything
// in simple datastructures.
//
// SIM:      Contains general information about the simulations like boxsize,
//           number of particles and voids.
//
// PARTICLE: Position [x,y,z], index in allpart array and tesselation volume
//           Positions are in units of the boxsize, i.e. in [0,1]
// 
// VOID:     A Zobov Void/Zone. Contains alot of info about the void including a
//           list of subvoids and particles in the void (NB: to get all particles
//           one must loop over subvoids)
//
// Defines:
//           VERBOSE : Output some info as we go along
//           SUPERVERBOSE : For debugging. Outputs ALOT
//           VELOCITY : Read velocity and store it with particles (in units of km/s)
//           ADJS : Read adj file and store tesselation particle adjecencies
//
// Structs are defined in global.h, the read-file routines are in init.h and the
// routines for computing stuff are in analysis.h
//
//=================================================================================

int main(int argc, char **argv) {

  //============================================
  // 
  // Read data from file
  //
  //============================================

  printf("\n*** Read inputfile           ***\n");
  init(argc,argv);

  printf("\n*** Read particle positions  ***\n");
  init_particle_position();

  printf("\n*** Read particle volumes    ***\n");
  init_particle_volume();

  printf("\n*** Read particle adjecents  ***\n");
  init_particle_adj();

  printf("\n*** Read voids               ***\n");
  init_void();

  printf("\n*** Read zones in voids      ***\n");
  init_zones_in_voids();

  printf("\n*** Read particles in zones  ***\n");
  init_particles_in_zones();

  //============================================
  // 
  // Analysis
  //
  //============================================

  printf("\n*** Compute void radius      ***\n");
  calc_effective_void_radius();
 
  printf("\n*** Compute void tree        ***\n");
  calc_void_tree();

  printf("\n*** Compute barycenter all   ***\n");
  calc_barycenter_all();

  printf("\n*** Walk tree and find voids ***\n");
  find_usable_voids_in_tree(); 

  printf("\n*** Free up memory and exit  ***\n");
  free_all_memory();
}


