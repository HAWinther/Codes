#ifndef IO_GADGET
#define IO_GADGET
#include <stdio.h>
#include <stdlib.h>

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void read_gadget_float_vector(FILE *fd, void *buffer, int dim, int npart);
void read_gadget_int_vector(FILE *fd, void *buffer, int dim, int npart);
void read_gadget_header(FILE *fd, int verbose);

//=======================================================
// Read binary method
//=======================================================

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream){
  size_t nread;
  if((nread = fread(ptr, size, nmemb, stream)) != nmemb){
    printf("I/O error (fread)");
    exit(1);
  }
  return nread;
}

//=======================================================
// Gadget header
//=======================================================

struct io_header{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char fill[60];
} header;

//=======================================================
// Print contents of gadget header
//=======================================================

void print_header(){
  printf("\n======================================\n");
  printf("Gadget Header:\n");
  printf("======================================\n");
  printf("N = %i %i %i %i %i %i\n",header.npart[0],header.npart[1],header.npart[2],header.npart[3],header.npart[4],header.npart[5]);
  printf("M = %e %e %e %e %e %e\n",header.mass[0],header.mass[1],header.mass[2],header.mass[3],header.mass[4],header.mass[5]);
  printf("a               = %f\n",header.time);
  printf("z               = %f\n",header.redshift);
  printf("BoxSize (Mpc/h) = %f\n",header.BoxSize/1000.0);
  printf("OmegaM          = %f\n",header.Omega0);
  printf("OmegaL          = %f\n",header.OmegaLambda);
  printf("HubbleParam     = %f\n",header.HubbleParam);
  printf("======================================\n\n");
}

//=======================================================
// Read gadget header
//=======================================================

void read_gadget_header(FILE *fd, int verbose){
  int tmp;
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(&header, sizeof(header), 1, fd);
  my_fread(&tmp,sizeof(int),1,fd);
  if(verbose) print_header();
}

//=======================================================
// Read gadget vector of floats of size npart * dim
//=======================================================

void read_gadget_float_vector(FILE *fd, void *buffer, int dim, int npart){
  int tmp;
  float *val = (float *)buffer;

  // Read
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(float), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);
}

//=======================================================
// Read gadget vector of ints of size npart * dim
//=======================================================

void read_gadget_int_vector(FILE *fd, void *buffer, int dim, int npart){
  int tmp;
  int *val = (int *)buffer;

  // Read
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(int), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);
}

/*
void read_and_bin_particles_gadget(char* fileprefix, int filenum, int *npart_tot, double *readbuffer, int *nbuffer){
  char filename[512];
  sprintf(filename,"%s%i",filenum);
  FILE *fp;
  int npart_now;
  char *buffer = (char *) readbuffer;

  // Open file
  fp = fopen(filename,"r");

  // Read header and print it for the first file only
  if(filenum==0)
    read_gadget_header(fp, 1);
  else
    read_gadget_header(fp, 0);

  npart_now = header.npart[1];;

  // Check that buffer is large enough
  if(3 * npart_now > *nbuffer){
    exit(1); 
  }
   
  // Verbose

  // Read positions
  read_gadget_float_vector(fp, readbuffer, 3, npart_now);

  process_particle_data((double *) readbuffer, npart_now, double(header.BoxSize));

  // Update how many particles in total we have read
  *npart_tot += npart_now;

  fclose(fp);
}
*/

#endif
