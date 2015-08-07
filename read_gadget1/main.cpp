#include "read_gadget.h"
#include <string>
#include <sstream>

//////////////////////////////////////////////////////
// Code to read GADGET1 files
// Hans A. Winther (2015) (hans.a.winther@gmail.com)
//////////////////////////////////////////////////////

//////////////////////////////////////
// A copy of to_string since its not in std
//////////////////////////////////////

template <typename T> std::string to_string(T value){
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

//////////////////////////////////////
// Store the positions of the particles
//////////////////////////////////////

void store_position(float *f, float boxsize, int npart_now, int npart_read){
  double x,y,z;

  // Verbose
  std::cout << "Read positions (First particle x/B = ";
  std::cout << f[0]/boxsize << " y/B = " << f[1]/boxsize << " z/B = " << f[2]/boxsize <<  std::endl;

  for(int i=0;i<npart_now;i++){
    // The current particle is number 'npart_read + i'

    // Position (x,y,z)
    x = f[3*i];
    y = f[3*i+1];
    z = f[3*i+2];

    // Store the data
    // ...
  }
}

//////////////////////////////////////
// Store the velocities of the particles
//////////////////////////////////////

void store_velocity(float *f, int npart_now, int npart_read){
  double vx,vy,vz;

  // Verbose
  std::cout << "Read velocities (First particle vx [km/s] = ";
  std::cout << f[0] << " vy [km/s] = " << f[1] << " vz [km/s] = " << f[2] <<  std::endl;

  for(int i=0;i<npart_now;i++){
    // The current particle is number 'npart_read + i'

    // Velocity (vx,vy,vz)
    vx = f[3*i];
    vy = f[3*i+1];
    vz = f[3*i+2];

    // Store the data
    // ...
  }
}

//////////////////////////////////////
// Read contents of a single GADGET1 file
//////////////////////////////////////

void read_single_gadget_file(std::string fileprefix, int filenum, int *numfiles, int *npart_read){
  std::string filename = fileprefix + to_string(filenum);
  FILE *fp;
  char *buffer;
  int npart_now;
  float *f;
  bool verbose = false;

  // Open file
  fp = fopen(filename.c_str(),"r");

  // Read header and print it for the first file only
  if(filenum==0)
    read_gadget_header(fp, true);
  else
    read_gadget_header(fp, false);

  // Get number of files
  if(filenum==0)
    *numfiles = header.num_files;

  // Allocate memory for read buffer
  npart_now = header.npart[1];;
  buffer = new char[sizeof(float) * npart_now * 3];
  f = (float *) buffer;

  // Verbose
  std::cout << std::endl;
  std::cout << "Reading file " << filenum+1 << " / " << *numfiles << std::endl;
  std::cout << "Npart_file = " << npart_now << " Numfiles = " << *numfiles << " Npart_read = " << *npart_read << std::endl;

  // Read positions
  read_gadget_float_vector(fp, buffer, 3, npart_now, verbose);
  store_position(f, header.BoxSize, npart_now, *npart_read);

  // Read velocity
  read_gadget_float_vector(fp, buffer, 3, npart_now, verbose);
  store_velocity(f, npart_now, *npart_read);

  // Update how many particles in total we have read
  *npart_read += npart_now;

  fclose(fp);
  delete[] buffer;
}

//////////////////////////////////////
// Read all GADGET1 files
//////////////////////////////////////

void read_all_gadget_files(std::string fileprefix){
  int nfiles, npart_read = 0;

  // Read first file to get the number of output files
  read_single_gadget_file(fileprefix, 0, &nfiles, &npart_read);

  // Loop over all the other files
  for(int i=1;i<nfiles;i++)
    read_single_gadget_file(fileprefix, i, &nfiles, &npart_read);

  // Verbose
  std::cout << std::endl;
  std::cout << "Done reading. We have read a total of " << npart_read << " particles" << std::endl;
}

int main(int argc, char **argv){
  read_all_gadget_files("/Users/winther/mg_comp/lcdm_z0.0/gadget.");
}
