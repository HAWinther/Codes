#include <sstream>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#ifdef OPENMP
#include <omp.h>
#endif
#include "global.h"
#include "io.h"
#include "analysis.h"

//=======================================================
//
// Code to read halos and profiles from AHF and bin them
//
// All positions in the code is in units of the boxsize
// i.e. in [0,1]
//
// All masses are in units of Msun/h
//
// All io related stuff is in io.h
// Output names are defined as prefix_name_suffix
// where prefix/suffix is defined in init() in io.h
//
// All of the analysis is done in analysis.h
//
// This file and global.h contain general methods
// and classes.
//
// Halo-file is assumed to be in the AHF format.
// Particles are assumed to be in ascii format.
// These things are easily changed in io.h
//
// The working precision can be set by changing realT
//
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//
//=======================================================

using namespace std;

std::pair<realT,int> calc_nhalos_in_massrange(realT mmin, realT mmax);

//=======================================================
// Calculate nhalos and <rvir> over all halos in the
// mass-range mmin < M < mmax
//=======================================================

pair<realT,int> calc_nhalos_in_massrange(realT mmin, realT mmax){
  int admissablehalos, nhalos = global.nhalos;
  Halo *halos = global.halos;
  realT rviravg;

  rviravg = 0.0;
  admissablehalos = 0;
  for(int i=0;i<nhalos;i++){
    if(halos[i].mass > mmin && halos[i].mass < mmax){
      rviravg += halos[i].rvir;
      admissablehalos++;
    }
  }
  rviravg /= realT(admissablehalos);

  return make_pair(rviravg,admissablehalos);
}

//=======================================================
// Calculate max-value of rvir in a mass range
//=======================================================

realT max_rvir_in_massrange(realT mmin, realT mmax){
  int nhalos = global.nhalos;
  Halo *halos = global.halos;
  realT rmax = 0.0;

  for(int i=0;i<nhalos;i++)
    if(halos[i].mass > mmin && halos[i].mass < mmax){
      if(halos[i].rvir > rmax) rmax = halos[i].rvir;
    }

  return rmax;
}

//=======================================================
// Bin the AHF profile files over a mass-range
//=======================================================

void bin_ahf_profile(realT mmin, realT mmax, int nbins, string filename){
  int admissablehalos, ncolprofilefile = NCOLPROFILEFILE, nprofile = ncolprofilefile-1;
  Halo *halos = global.halos;
  realT rviravg, rmin, rmax;
  ofstream fp;

  // Count number of halos and get rviravg
  pair<realT,int> rvir_and_nhalo = calc_nhalos_in_massrange(mmin, mmax);
  rviravg = rvir_and_nhalo.first;
  admissablehalos = rvir_and_nhalo.second;

  // XXX ADD BIN IN PHYSICAL UNITS XXX

  cout << "========================================" << endl;
  cout << "Average profile of AHF halos... M = " << mmin << " - " << mmax << endl;
  cout << "========================================" << endl;

  // Return if no halos are found
  if(admissablehalos == 0){
    cout << "No halos found..." << endl;
    return;
  } else {
    cout << "Taking average over nhalos = " << admissablehalos << endl;
  }

  // Define the bins (in terms of rvir)
  rmin  = RVIRMINFAC;
  rmax  = 1.0;

  // Allocate
  realT **y = new realT*[nprofile];
  realT **sigma = new realT*[nprofile];
  realT *r  = new realT[nbins];
  for(int k=0;k<nprofile;k++){
    y[k] = new realT[nbins];
    sigma[k] = new realT[nbins];
  }
  for(int j=0;j<nbins;j++){
    r[j] = exp(log(rmin) + log(rmax/rmin)*(j+0.5)/realT(nbins));
    for(int k=0;k<nprofile;k++){
      y[k][j] = sigma[k][j] = 0.0;
    }
  }

  // Loop over halos and bin
  for(int i=0;i<global.nhalos;i++){
    if(halos[i].mass > mmin && halos[i].mass < mmax){
      for(int k=0;k<nprofile;k++){
        for(int j=0;j<nbins;j++){
          y[k][j] += halos[i].get_profile(k, halos[i].rvir * r[j] * global.boxsizeinkpc);
        }
      } 
    }
  }

  // Make average
  for(int k=0;k<nprofile;k++){
    for(int j=0;j<nbins;j++){
      y[k][j] /= realT(admissablehalos);
    }
  }

  // Calculate standard deviation
  for(int i=0;i<global.nhalos;i++){
    if(halos[i].mass > mmin && halos[i].mass < mmax){
      for(int k=0;k<nprofile;k++){
        for(int j=0;j<nbins;j++){
          sigma[k][j] += pow2(y[k][j] - halos[i].get_profile(k, halos[i].rvir * r[j] * global.boxsizeinkpc));
        }
      } 
    }
  }

  // Make average
  for(int k=0;k<nprofile;k++){
    for(int j=0;j<nbins;j++){
      sigma[k][j] = sqrt(sigma[k][j] / realT(admissablehalos)); 

      // Standard error
      sigma[k][j] = sigma[k][j]  / sqrt(realT(admissablehalos)); 
    }
  }

  // Output
  fp.open(filename.c_str());
  fp << "#Stacked ahf profile with nhalos = " << admissablehalos << " of mass " << mmin << " to " << mmax << endl; 
  for(int j=0;j<nbins;j++){
    realT rfac = rviravg * global.boxsizeinkpc;
    rfac = 1.0;
    fp << r[j] * rfac << " ";
    for(int k=0;k<nprofile;k++){
      fp << y[k][j] << " " << sigma[k][j] << " ";
    }
    fp << endl;
  }
  fp.close();

  // Clean up
  for(int k=0;k<nprofile;k++)
    delete[] y[k];
  delete[] r;
  delete[] y;
}

//=======================================================
// Main is main! Allocate memory and do some work!
//=======================================================

int main(int argv, char **argc){

  //==============================================
  // Read parameters from file / set them here 
  //==============================================
  init(argv, argc);

  //==============================================
  // Analysis
  //==============================================
  analysis();

  return 0;
}

