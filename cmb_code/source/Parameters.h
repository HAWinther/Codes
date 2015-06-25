
//========================================================
// Here all the parameters used by the code is set
//========================================================

#pragma once

#include <cmath>
#include "Matrix.h"

// Units and constants
static const double Mpc     = 3.085680e22;			    // Megaparsec in m 
static const double eV      = 1.60217e-19;			    // Electronvolt in J
static const double pi      = 3.141592653;          // Just Pi
static const double k_b     = 1.38065e-23;			    // Bolzmanns constant in J/K
static const double m_e     = 9.10938e-31;			    // Mass of electron in kg
static const double m_H     = 1.67352e-27;				  // Mass of hydrogen atom in kg
static const double c       = 2.9979245e8;				  // Speed of light m/s
static const double G       = 6.67258e-11;				  // Gravitational constant in N*m^2/kg^2
static const double hbar    = 1.05457e-34;		      // Reduced Plancks constant in Js
static const double sigma_T = 6.65246e-29;			    // Thomas scattering cross-section in m^2

// Derived constants
static const double epsilon_0 = 13.6056917*eV;	    // Ionization energy for the ground state of hydrogen in J
static const double xhi0      =	24.5874*eV;			    // Ionization energy for neutral Helium in J
static const double xhi1      = 4.0*epsilon_0;		  // Ionization energy for singly ionized Helium in J
static const double Mpl       = 1.0/sqrt(8.0*pi*G); // Planck mass in units of hbar=c=1

// Parameters
static const double k_pivot		   = 0.002;		        // Pivot scale for the power-spectrum at which A_s is defined (in units of H_0/c) Not used in this version!

// Constants we should not have to change
static const double x_init       = -15.0;					  // Initial x for integration of the perturbations	
static const double x_start_rec  = -log(1631.0);		// x at start of recombination
static const double x_end_rec    = -log(615.0); 		// x at end of recombination
static const double saha_limit   = 0.99;		   	    // X_e for when Saha equation needs to be replaced by the Peebles equation
static const double x_start_main = log(1.0e-10);	  // Initial x for integrations of eta / tau / X_e / g

// Specification for x_main array 
static const int nmain     = 5000;                  // Number of points for integration (and spline) of eta / tau / X_e / g		
static const int np1       = int(4*nmain/10.0);     // Before reionization
static const int np2       = int(2*nmain/10.0);     // During reionization
static const int np3       = nmain - np1 - np2;     // After  reionization

// Max and min wavenumbers and number of points in k-array
static const int n_k       = 500;					          // Number of points in k-space for integration of perturbations
static const int n_k_high  = 8000;					        // Points in high resolution k grid
static const double k_min  = 0.0001 * 2998.0;     	  // Minimum wavenumber in k integration (final value in units of H_0/c)
static const double k_max  = 1.0  * 2998.0;			  // Maximum wavenumber in k integration

// Number of ODE's we have to solve. Cutoff of l-hiearachy ofr temp. and polarization
static const int lmax_int  = 8;					            // Highest l-moment in the determination of Theta and Theta_p
static const int lmax_nu   = 10;					          // Highest l-moment in the determination of Nu
static const int npar      = 8+2*lmax_int+lmax_nu;	// Total number of coupled ODE's for the full perturbations equations
static const int ntight    = 8+lmax_nu;				      // Total number of coupled ODE's during tight coupling

static const int n_bessel  = 20000;					        // Number of points for Bessel function evaluation

// Number of points in time
static const int n1        = 300;							      // Points during recombination
static const int n2        = 150;							      // Points between recombination and reionization
static const int n3        = 100;							      // Points during reionization
static const int n4        = 200;							      // Points from end of reionization and until today
static const int ntot      = n1+n2+n3+n4;				    // Total points between rec and today
static const int ntot_high = 5000;					        // Total points in high resolution x grid
static const int n1_high   = (n1*ntot_high)/ntot;   // High resolution equivalents...
static const int n2_high   = (n2*ntot_high)/ntot;
static const int n3_high   = (n3*ntot_high)/ntot;
static const int n4_high   = ntot_high-n1_high-n2_high-n3_high;

// k value when to switch from BS to SIE method in integration of perturbations
static const int ktight    = int(0.6*double(n_k)); // Should be derived from k_max
static const int kfull     = int(0.6*double(n_k));

// Integrations paraemeters for perturbations ODE (StepperBS)
static const double r_bs   = 1.0e-7;
static const double a_bs   = 1.0e-7;
static const double h1_bs  = 1.0e-3;

// Integrations paraemeters for perturbations ODE (StepperSie) (8,8,6 high precision, 4,4,3 fast)
static const double r_sie  = 1.0e-5;
static const double a_sie  = 1.0e-5;
static const double h1_sie = 1.0e-3;

// Array of l's for C_l evaluation
static const int ls[]		= 
{	2,		3,		4,		6,		8,		10,		12,		15,		20,		30,		40,
  50,		60,		70,		80,		90,		100,	120,	140,	160,	180,	200,
  225,	250,	275,	300,	350,	400,	450,	500,	550,	600,	650,
  700,	750,	800,	850,	900,	950,	1000,	1050,	1100,	1150,	1200,
  1250,	1300,	1350,	1400,	1450,	1500,	1550,	1600,	1650,	1700,	1750, 
  1800,	1850,	1900,	1950,	2000};
static const int l_num = sizeof (ls) / sizeof (ls[0]);

// Parameter container
struct ParamSimu {
  ifstream infile;

  ParamSimu(){
    init("input_param.txt");
  }

  ParamSimu(string filename){
    init(filename);
  }

  void init(string filename){

    // Start timer
    time_start = clock();

    // Read parameters from file
    infile.open(filename.c_str());

    // Check if file exist and read parameters
    if (infile.is_open()) {
      init_param();
    } else {
      int ii;
      cout << endl;
      cout << "======================================" << endl;
      cout << "Error:   Parameterfile do not exit    " << endl;
      cout << "======================================" << endl;
      cout << "Enter 1 to proceed with standard param" << endl;
      cin  >> ii;
      cout << "======================================" << endl;
      if(ii==1)
        init_param_standard();
      else
        exit(1);
    }
  }

  // Filenames
  string outprefix;

  // Use precomputed data if possible?
  bool use_precomputed_pert;
  bool use_precomputed_bessel;
  string precomputed_pert_phifile;
  string precomputed_pert_sourcefile;
  string precomputed_besselfile;

  // Time control
  unsigned time_start, time_now;

  /*********************************************************************************************/
  // Cosmological parameters
  /*********************************************************************************************/
  double Omega_m;			  // Dark matter density parameter today
  double Omega_b;			  // Baryonic matter density parameter today
  double Omega_r;			  // Radiation density parameter today
  double Omega_nu;		  // Neutrino density parameter today
  double Omega_lambda;	// Dark energy density parameter today
  double hubb;			    // Hubble factor hubb = H0 / (100km/sMpc)
  double n_nu_eff;		  // Effective number of relativistic neutrion species 				
  double z_reion;			  // Reionization redshift
  double dz_reion;		  // Reionization redshift width
  double n_s;				    // Scalar spectral index

  double A_s;				    // Initial amplitude
  double sigma8;			  // RMS fluctuations on scale 8 Mpc/h
  double first_peak;		// Amplitude of first peak in muK^2
  /*********************************************************************************************/

  int norm_cl;			    // Describes type of normalization - 0: first-peak, 1: sigma8, 2; A_s
  double norm;			    // Factor used to normalize power-spectrums after final normalization

  double T_0;
  double Y_p;

  // Derived contants
  double H_0 ;			// Hubble constant today in 1/s
  double R0;				// 1.33 * Omega_r/Omega_b
  double f_nu;			// Effective number of relativistic neutrinos
  double rhoc;			// Critical density today
  double x_start_reion;	// Start of reionization array
  double x_end_reion;		// End of reionization array

  // Initialize parameters from file
  void init_param(){
    infile >> norm_cl;
    infile >> Omega_m;
    infile >> Omega_b;
    infile >> Omega_lambda;
    infile >> hubb;
    infile >> n_nu_eff;
    infile >> z_reion;
    infile >> dz_reion;
    infile >> n_s;
    infile >> A_s;
    infile >> sigma8;
    infile >> first_peak;

    infile >> T_0;
    infile >> Y_p;

    infile >> outprefix;

    // Use precomputed data?
    int tmp;
    infile >> tmp;
    use_precomputed_pert = (tmp==1);
    infile >> precomputed_pert_phifile;
    infile >> precomputed_pert_sourcefile;
    infile >> tmp;
    use_precomputed_bessel = (tmp==1);
    infile >> precomputed_besselfile;

    Omega_r  = 2.469e-5/(hubb*hubb) * pow(T_0/2.725,4.0);
    Omega_nu = 0.2271 * n_nu_eff * Omega_r;

    H_0		= hubb * 100.0*1.0e3/Mpc;

    R0		= 4.0*Omega_r/(3.0*Omega_b);
    f_nu	= 1.0/(8.0/(7.0*n_nu_eff)*pow(11.0/4.0,4.0/3.0)+1.0);
    rhoc	= 3.0*H_0*H_0*Mpl*Mpl;

    x_start_reion	= -log(1.0+z_reion + 6.0*dz_reion);
    x_end_reion		= -log(1.0+z_reion - 6.0*dz_reion);

#ifdef VERBOSE
    cout <<	endl;
    cout << "**********************************************"	<< endl;
    cout <<	"***       Reading parameters from file     ***"	<< endl;
    cout <<	"**********************************************"	<< endl;
#endif

    print_ic();
  }

  void print_ic(){
#ifdef VERBOSE
    cout <<     "Outputprefix                    : " << outprefix << endl;
    cout <<			"Omega_m                         : " << Omega_m		<< endl;
    cout <<			"Omega_b                         : " << Omega_b		<< endl;
    cout <<			"Omega_r                         : " << Omega_r		<< endl;
    cout <<			"Omega_nu                        : " << Omega_nu	<< endl;
    cout <<			"Omega_L                         : " << Omega_lambda << endl;
    cout <<			"h                               : " << hubb		<< endl;
    cout <<			"n_s                             : " << n_s			<< endl;

    if (z_reion>0.0) {
      cout <<		"z_reion                         : " << z_reion		<< endl;
      cout <<		"dz_reion                        : " << dz_reion	<< endl;
    } else {
      cout <<	endl << "NOT INCLUDING REIONIZATION EFFECTS" << endl;
    }

    cout <<			"Y_p                             : " << Y_p			<< endl;
    cout <<			"T_0                             : " << T_0			<< endl;
    cout << endl;

    switch (norm_cl) {
      case 0:
        cout << "Using maximum normalization" << endl;
        cout << "Max Amplitude (muK^2)           : " << first_peak << endl;
        sigma8		= 0.0;
        first_peak	= 5775.0;
        break;
      case 1:
        cout << "Using sigma8 normalization" << endl;
        cout << "sigma8                          : " << sigma8 << endl;
        first_peak = 5775.0;
        break;
      case 2:
        cout << "Using initial P(k) amplitude" << endl;
        cout << "A_s                             : " << A_s << endl;
        sigma8		= 0.0;
        first_peak	= 5775.0;
        break;
      default:
        cout << "Normalization method does not exist" << endl;
        exit(1);
        break;
    }
    norm = 1.0;

    cout << endl;
    if (use_precomputed_pert) {
      cout << "Use precomputed perturbation data" << endl;
      cout << "--> " << precomputed_pert_phifile << endl;
      cout << "--> " << precomputed_pert_sourcefile << endl;
    }
    cout << endl;
    if (use_precomputed_bessel) {
      cout << "Use precomputed bessel functions" << endl;
      cout << "--> " << precomputed_besselfile <<  endl;
    }

    cout <<			"**********************************************" << endl;
#endif
  }

  // Initialize to standard values
  void init_param_standard(){
    norm_cl		 = 0;
    hubb		   = 0.72;
    n_nu_eff	 = 3.04;
    z_reion		 = 10.5;
    dz_reion	 = 0.5;
    n_s			   = 0.96;
    A_s			   = 3.14e-9;
    sigma8		 = 0.8;
    first_peak = 5775.12;
    T_0			   = 2.725;
    Y_p			   = 0.24;

    outprefix  = "output/cmfb_";

    use_precomputed_pert = false;
    precomputed_pert_phifile = "";
    precomputed_pert_sourcefile = "";
    use_precomputed_bessel = false;
    precomputed_besselfile = "";

    Omega_m		 = 0.224;
    Omega_b		 = 0.046;
    Omega_r		 = 2.469e-5/(hubb*hubb) * pow(T_0/2.725,4.0);
    Omega_nu	 = 0.2271 * n_nu_eff * Omega_r;
    Omega_lambda = 1.0 - Omega_m - Omega_b - Omega_r - Omega_nu;

    // Derived constants
    H_0		= hubb * 100.0*1.0e3/Mpc;
    R0		= 4.0*Omega_r/(3.0*Omega_b);
    f_nu	= 1.0/(8.0/(7.0*n_nu_eff)*pow(11.0/4.0,4.0/3.0)+1.0);
    rhoc	= 3.0*H_0*H_0*Mpl*Mpl;

    x_start_reion	= -log(1.0+z_reion + 6.0*dz_reion);
    x_end_reion		= -log(1.0+z_reion - 6.0*dz_reion);

#ifdef VERBOSE
    cout <<	endl;
    cout <<     "**********************************************" << endl;
    cout <<			"***  Using WMAP LCDM best fit parameters   ***" << endl;
    cout <<			"**********************************************" << endl;
#endif

    print_ic();
  }
};

