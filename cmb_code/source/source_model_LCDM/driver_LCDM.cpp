
//==================================================
// Main method where everything runs from
//
// Here we go step by step
// though and calculate the CMB spectra
// by first doing background the recombiantion
// then perturbations and finally line of
// sight integration to get power-spectra
//=================================================

#include "../Parameters.h"
#include "BackgroundLCDM.h"
#include "../source_recmodel_STANDARD/RecombinationSTANDARD.h"
#include "PerturbationsLCDM.h"
#include "../Powerspectrum.h"

int main(int argv, char **argc){
  string filename;

  // Read paramfile from input
  if(argv <= 1){
    cout << "Run as ./run paramfile" << endl;
    exit(1);
  } else {
    filename = argc[1];
  }

#ifdef VERBOSE
  cout <<	endl;
  cout <<  "**********************************************" << endl;
  cout <<  "***         STARTING CM-FUCKING-B          ***" << endl;
  cout <<	 "**********************************************" << endl;
#endif

  // Pointers to (general) objects
  ParamSimu *cosm_simu;
  Background *B;
  RecombinationSTANDARD *R;
  Perturbations *P;
  Powerspectrum *PS;

  // Initialize parameters
  cosm_simu = new ParamSimu(filename);

  // Calculate bakcground expansion for LCDM
  B = new BackgroundLCDM(cosm_simu);
  B->calc_expansion_history();

  // Calculate recombination history 
  R = new RecombinationSTANDARD(B, cosm_simu);
  R->calc_rec_history();

  // Integrate perturbations for LCDM
  P = new PerturbationsLCDM(B, R, cosm_simu);
  P->integrate_perturbations();

  // Calculate power-spectrum(s)
  PS = new Powerspectrum(B, R, P, cosm_simu);
  PS->calc_scalar_powerspectrums();

#ifdef VERBOSE
  cosm_simu->time_start = clock() - cosm_simu->time_start;
  cout << endl << "---> Exiting CMFB" << endl;
  cout <<			"     Total Time used (sec)      : " << cosm_simu->time_start*1.0e-6 << endl;
  cout <<			"**********************************************" << endl;
#endif
};
