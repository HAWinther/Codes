#ifndef GLOBALHEADER_INC
#define GLOBALHEADER_INC
#include <string.h>

// Read input file with lua
#include "readlua.h"

// Add spline class
#include "spline.h"

// Add halo class
#include "halo.h"

// Useful method
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

// Working precision
#if defined(DOUBLE)
typedef double realT;
#elif defined(LONGDOUBLE)
typedef long double realT;
#elif defined(FLOAT)
typedef float realT;
#else
typedef double realT;
#endif

//=======================================================
// Global variables / pointers
//=======================================================

struct GlobalInfo{
  Halo *halos;
  int nhalos;
  realT boxsizeinkpc;
  realT aexp, omegam;
  int nbins;

  std::string halofile;
  std::string haloprofilefile;

  std::string outfileprefix;
  std::string outfilesuffix;

  int nprofiles;
  double Mmin[100], Mmax[100];
} global;

//=======================================================
// Methods needed globally
//=======================================================

void bin_ahf_profile(realT mmin, realT mmax, int nbins, std::string filename);
#endif
