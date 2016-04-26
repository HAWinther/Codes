#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<iostream>

//===========================================
// Write to file a RAMSES ic_deltab file
// Input: levelmin (gridsize = [2^levelmin])
//        astart    (start aexp)
//        boxlen_ini (boxsize in Mpc/h)
//        outdir (string with)
// Output: ic_deltab
//===========================================

void write_ramses_header(float astart, int levelmin, float omega_m, float omega_l, float h0, float boxlen_ini, std::string outfile){
	FILE *fp;
	int tmp, n1, n2, n3;
	float xoff1, xoff2, xoff3, dx;

	// Assume zero offset
	xoff1 = xoff2 = xoff3 = 0.0;

	// Assumes the same in all directions n = 2^levelmin
	n1 = n2 = n3 = 1 << levelmin;

	// dx in Mpc
	dx  = boxlen_ini / float(n1) * 100.0/h0;    

	// Number of floats and ints we write
	tmp = 8*sizeof(float) + 3*sizeof(int);

	// Open file
	if( (fp = fopen(outfile.c_str(), "w")) == NULL){
		std::cout << "Error opening file: " << outfile << std::endl;
		exit(0);
	}

	// Write the file
	fwrite(&tmp,     sizeof(int),   1, fp);
	fwrite(&n1,      sizeof(int),   1, fp);
	fwrite(&n2,      sizeof(int),   1, fp);
	fwrite(&n3,      sizeof(int),   1, fp);
	fwrite(&dx,      sizeof(float), 1, fp);
	fwrite(&xoff1,   sizeof(float), 1, fp);
	fwrite(&xoff2,   sizeof(float), 1, fp);
	fwrite(&xoff3,   sizeof(float), 1, fp);
	fwrite(&astart,  sizeof(float), 1, fp);
	fwrite(&omega_m, sizeof(float), 1, fp);
	fwrite(&omega_l, sizeof(float), 1, fp);
	fwrite(&h0,      sizeof(float), 1, fp);
	fwrite(&tmp,     sizeof(int),   1, fp);
	fclose(fp);
}

//===========================================
// Read parameters from file and write a
// ic_deltab file
//===========================================

void generate_ic_deltab(int argv, char ** argc){
	int levelmin;
	float astart, omegam, omegal, h0, boxsize;
  std::string outfile;

	// Check if we have enough input
	if(argv < 7){
		std::cout << "Run as ./genicdeltab astart levelmin omegam omegal h0 boxsize outputdir" << std::endl;
		std::cout << "E.g ./genicdeltab 0.5 7 0.3 0.7 70.0 200.0 /usr/hansw/ic_files" << std::endl;
		exit(1);
	} else {
		// Get input from standard input
		astart   = atof(argc[1]);
		// levelmin is in 3D related to particle number as npart_tot = 2^(3*levelmin)
		levelmin = atoi(argc[2]); 
		omegam   = atof(argc[3]);
		omegal   = atof(argc[4]);
		h0       = atof(argc[5]);
		boxsize  = atof(argc[6]);
		if(argv==7)
		  outfile = "ic_deltab";
    else 
      outfile = std::string(argc[7]) + "/ic_deltab";
	}

	// Verbose
	std::cout << "=================================" << std::endl;
  std::cout << "Generate ic_deltab with param:   " << std::endl;
  std::cout << "=================================" << std::endl;
	std::cout << "astart          = "    << astart   << std::endl;
	std::cout << "levelmin        = "    << levelmin << std::endl;
	std::cout << "Omega_m         = "    << omegam   << std::endl;
	std::cout << "Omega_l         = "    << omegal   << std::endl;
	std::cout << "h               = "    << h0       << std::endl;
	std::cout << "Boxsize (Mpc/h) = "    << boxsize  << std::endl;
	std::cout << "Output file     = "    << outfile   << std::endl;

	// Write the file
	write_ramses_header(astart, levelmin, omegam, omegal, h0, boxsize, outfile);
}

int main(int argv, char ** argc){
  generate_ic_deltab(argv, argc);
}
