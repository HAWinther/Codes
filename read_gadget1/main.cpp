#include "read_gadget.h"
#include <string>
#include <sstream>

//////////////////////////////////////
// A copy of to_string since its not in std
//////////////////////////////////////

template <typename T> std::string to_string(T value){
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

//////////////////////////////////////
// Read contents of a GADGET1 file
//////////////////////////////////////

void read_single_gadget_file(int filenum){
  std::string filename = "/Users/winther/ahf_test_data/gadget." + to_string(filenum);
  FILE *fp;
  char *buffer;
  int npart_now, NDIM = 3;
  float *f;

  // Open file
  fp = fopen(filename.c_str(),"r");

  // Read header
  read_gadget_header(fp, true);

  // Allocate memory
  npart_now = header.npart[1];;
  buffer = new char[sizeof(float) * npart_now * NDIM];
  f = (float *) buffer;

  std::cout << "Npart = " << npart_now << std::endl;

  // Read positions
  std::cout << "Read positions: " << std::endl;
  read_gadget_float_vector(fp, buffer, NDIM, npart_now, true);
  for(int i=0;i<npart_now;i++)
    for(int k=0;k<NDIM;k++){
      // Do something with the data
      // Position x[k] of particle 'i' is f[3*i+k]
    }

  // Read velocity
  std::cout << "Read velocity: " << std::endl;
  read_gadget_float_vector(fp, buffer, NDIM, npart_now, true);
  for(int i=0;i<npart_now;i++)
    for(int k=0;k<NDIM;k++){
      // Do something with the data
      // Velocity v[k] of particle 'i' is f[3*i+k]
    }

  fclose(fp);
  delete[] buffer;
}

int main(int argc, char **argv){
  read_single_gadget_file(0);
}
