#ifndef READGADGET_INC
#define READGADGET_INC
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void read_gadget_float_vector(FILE *fd, void *buffer, int dim, int npart, bool print);
void read_gadget_int_vector(FILE *fd, void *buffer, int dim, int npart, bool print);
void read_gadget_header(FILE *fd, bool verbose);

//=======================================================
// Read binary method
//=======================================================

size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream){
  size_t nread;
  if((nread = fread(ptr, size, nmemb, stream)) != nmemb){
    std::cout << "I/O error (fread) ";
    std::cout << "nread=" << nread << " size=" << size << " nmemb=" << nmemb << std::endl;
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
  printf("\n");
  printf("======================================\n");
  printf("Gadget Header:\n");
  printf("======================================\n");
  printf("N=%i %i %i %i %i %i\n",header.npart[0],header.npart[1],header.npart[2],header.npart[3],header.npart[4],header.npart[5]);
  printf("M=%e %e %e %e %e %e\n",header.mass[0],header.mass[1],header.mass[2],header.mass[3],header.mass[4],header.mass[5]);
  printf("a=%e\n",header.time);
  printf("z=%e\n",header.redshift);
  printf("BoxSize=%e\n",header.BoxSize);
  printf("OmegaM=%e\n",header.Omega0);
  printf("OmegaL=%e\n",header.OmegaLambda);
  printf("HubbleParam=%e\n",header.HubbleParam);
  printf("======================================\n");
  printf("\n");
}

//=======================================================
// Read gadget header
//=======================================================

void read_gadget_header(FILE *fd, bool verbose = true){
  int tmp;
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(&header, sizeof(header), 1, fd);
  my_fread(&tmp,sizeof(int),1,fd);
  if(verbose) print_header();
}

//=======================================================
// Read gadget vector of floats of size npart * dim
//=======================================================

void read_gadget_float_vector(FILE *fd, void *buffer, int dim, int npart, bool verbose = true){
  int tmp;
  float *val = (float *)buffer;

  // Read
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(float), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);

  // Verbose
  if(verbose){
    std::cout << "***********************" << std::endl;
    for(int i=0;i<npart;i++){
      if(i>=npart-10 || i<=10){
        std::cout << "i = " << i << " ";
        for(int k=0;k<dim;k++){
          std::cout << val[dim*i+k] << " ";
        }
        std::cout << std::endl;
      }
    }
    std::cout << "***********************" << std::endl;
    std::cout << std::endl;
  }
}

//=======================================================
// Read gadget vector of ints of size npart * dim
//=======================================================

void read_gadget_int_vector(FILE *fd, void *buffer, int dim, int npart, bool verbose = true){
  int tmp;
  int *val = (int *)buffer;

  // Read
  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(buffer, dim*sizeof(int), npart, fd);
  my_fread(&tmp,sizeof(int),1,fd);

  // Verbose
  if(verbose){
    std::cout << "***********************" << std::endl;
    for(int i=0;i<npart;i++){
      if(i>=npart-3 || i<=3){
        for(int k=0;k<dim;k++){
          std::cout << val[dim*i+k] << " ";
        }
        std::cout << std::endl;
      }
    }
    std::cout << "***********************" << std::endl;
    std::cout << std::endl;
  }
}

//=======================================================
// Read a gadget2 label[4]
//=======================================================

int read_gadget2_label(FILE *fd, char *labelin){
  int nextblock, tmp;
  char label[4];

  my_fread(&tmp,sizeof(int),1,fd);
  my_fread(&label, sizeof(char), 4, fd);
  my_fread(&nextblock, sizeof(int), 1, fd);
  my_fread(&tmp,sizeof(int),1,fd);

  my_fread(&tmp,sizeof(int),1,fd);

  printf("Reading header => '%c%c%c%c' (%d byte)\n", label[0], label[1], label[2], label[3], nextblock);

  for(int j=0;j<4;j++) labelin[j] = label[j];
  return nextblock;
}

#endif
