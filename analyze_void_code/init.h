#ifndef IO_H_INC
#define IO_H_INC
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "io_gadget.h"

void init(int argc, char** argv);
void init_particle_position();
void init_particle_volume();
void init_void();
void init_particles_in_zones();
void init_zones_in_voids();
void init_particle_adj();
int  free_all_memory();

void init(int argc, char** argv){
  FILE *fp;

  //====================================================
  //
  // Read the inputfile and initializes some variables
  //
  // Inputfile format (ascii):
  // 
  // PARTFILE		  (particle pos and vel) 
  // VOLFILE		  (volume of particles, ZOBOV)
  // VOIDFILE		  (list of voids, ZOBOV)
  // PFILE		    (particles in voids, ZOBOV)
  // ZFILE		    (subvoids in voids, ZOBOV)
  // ADJFILE		  (particle adjs, ZOBOV)
  // VOIDTREE     (filename to store voidtree as
  //               if this file exists it will be read
  //               if one tries to compute the void tree)
  //
  //====================================================

  // Initialize OpenMP
#ifdef USEOPENMP
#pragma omp parallel
  {
    if(omp_get_thread_num() == 0){
      global.nthreads = omp_get_num_threads();
      printf("OpenMP in use. We have %i threads availiable\n",global.nthreads);
    }
  }
#else
  global.nthreads = 1;
#endif

  // Validate input
  if(argc == 1){
    printf("Run as ./analyze input.txt\n");
    exit(1);
  }

  // Read parameters from inputfile
  if((fp = fopen(argv[1],"r")) != NULL){
    printf("Read parameters from file %s\n", argv[1]);
    fscanf(fp,"%s\n",&global.partfile);
    fscanf(fp,"%s\n",&global.volfile);
    fscanf(fp,"%s\n",&global.voidfile);
    fscanf(fp,"%s\n",&global.partinzonefile);
    fscanf(fp,"%s\n",&global.zoneinvoidfile);
    fscanf(fp,"%s\n",&global.adjfile);
    fscanf(fp,"%s\n",&global.voidtreefile);
    fclose(fp);
  } else {
    printf("Error opening parameterfile %s\n", argv[1]);
    exit(1);
  }

#ifdef VERBOSE
  printf("Partfile = %s\n",global.partfile);
  printf("Volfile  = %s\n",global.volfile);
  printf("Voidfile = %s\n",global.voidfile);
  printf("Pfile    = %s\n",global.partinzonefile);
  printf("Zfile    = %s\n",global.zoneinvoidfile);
  printf("Adjfile  = %s\n",global.adjfile);
  printf("Treefile = %s\n",global.voidtreefile);
#endif
}

void init_particle_position(){
  int nfiles, np, i, j, nbuffer;
  float *buffer;

  //========================================================
  //
  // Reads GADGET1 format file (no labels).
  // Assumes boxsize is in [kpc/h] when storing boxsize
  // Allocates [allpart] array
  // Stores positions in units of the boxsize, i.e. in [0,1]
  // If VELOCITY is defined it reads and stored velocity
  // in units of km/s
  //
  //========================================================

#ifdef VERBOSE
  printf("Read gadgetfiles from %s\n", global.partfile);
#endif

  // Sum up number of particles and find max particles in file
  np = nbuffer = 0;
  for(i=0;;i++){
    FILE *fp;
    char filename[512];

    sprintf(filename,"%s%i", global.partfile, i);
    fp = fopen(filename,"r");
    read_gadget_header(fp, 0);
    fclose(fp);

    np += header.npart[1];
    if(header.npart[1] > nbuffer) nbuffer = header.npart[1];
    if(i == 0) {
      nfiles = header.num_files;
      global.boxsize = header.BoxSize/1000.; // Assumes boxsize is given in [kpc/h] 
    }
    if(i == nfiles-1) break;
  }

  // Allocate particles and buffer
  buffer = (float *)malloc(3*nbuffer*sizeof(float));
  allpart = (PARTICLE *)malloc(np*sizeof(PARTICLE));

#ifdef VERBOSE
  printf("We found %i files with np = %i (nbuffer = %i). Boxsize: %lf Mpc/h\n", nfiles, np, nbuffer, global.boxsize);
#ifdef VELOCITY
  printf("Memory needed to store particles (x,y,z,vol,index) not counting adjs: %lf GB\n", (7*8+1*4)*np/1024./1024./1024.);
#else
  printf("Memory needed to store particles (x,y,z,vol,index) not counting adjs: %lf GB\n", (4*8+1*4)*np/1024./1024./1024.);
#endif
#endif

  // Loop over files and read positions
  np = 0;
  for(i=0;i<nfiles;i++){
    FILE *fp;
    int curnp;
    double boxfac;
    char filename[512];

    // Make filename
    sprintf(filename,"%s%i", global.partfile, i);
#ifdef SUPERVERBOSE
    printf("Opening file %s\n",filename);
#endif
    fp = fopen(filename,"r");

    // Read header and get number of particles
    read_gadget_header(fp, 0);
    curnp = header.npart[1];
    boxfac = 1.0/header.BoxSize;

    // Read positions
#ifdef SUPERVERBOSE
    printf("[%i] Reading np = %i\n", i, curnp);
#endif
    read_gadget_float_vector(fp, buffer, 3, curnp);

    // Process positions
    for(j=0;j<curnp;j++){
      allpart[np+j].x = (double) buffer[3*j+0]*boxfac;
      allpart[np+j].y = (double) buffer[3*j+1]*boxfac;
      allpart[np+j].z = (double) buffer[3*j+2]*boxfac;

      // Periodic boundary conditions
      if(allpart[np+j].x >= 1.0) allpart[np+j].x -= 1.0;
      if(allpart[np+j].y >= 1.0) allpart[np+j].y -= 1.0;
      if(allpart[np+j].z >= 1.0) allpart[np+j].z -= 1.0;
    }

    // Read velocities
#ifdef VELOCITY
    read_gadget_float_vector(fp, buffer, 3, curnp);
    for(j=0;j<curnp;j++){
      allpart[np+j].vx = (double) buffer[3*j+0];
      allpart[np+j].vy = (double) buffer[3*j+1];
      allpart[np+j].vz = (double) buffer[3*j+2];
    }
#endif

    // Add up number of particles
    np += curnp;
    fclose(fp);
  }

  // Initialize the rest which we don't know yet
  global.np = np;
  for (i = 0; i<np; i++){
    allpart[i].n      = i;
    allpart[i].volume = 0.0;
#ifdef ADJS
    allpart[i].nadj   = 0;
    allpart[i].adj    = NULL;
#endif
  }

#ifdef VERBOSE
  printf("Done reading %i particles\n", global.np);
  printf("x0   x1 ...  xn: %lf %lf ... %lf\n", allpart[0].x, allpart[1].x, allpart[np-1].x);
  printf("y0   y1 ...  yn: %lf %lf ... %lf\n", allpart[0].y, allpart[1].y, allpart[np-1].y);
  printf("z0   z1 ...  zn: %lf %lf ... %lf\n", allpart[0].z, allpart[1].z, allpart[np-1].z);
#ifdef VELOCITY
  printf("vx0 vx1 ... vxn: %lf %lf ... %lf\n", allpart[0].vx, allpart[1].vx, allpart[np-1].vx);
  printf("vy0 vy1 ... vyn: %lf %lf ... %lf\n", allpart[0].vy, allpart[1].vy, allpart[np-1].vy);
  printf("vz0 vz1 ... vzn: %lf %lf ... %lf\n", allpart[0].vz, allpart[1].vz, allpart[np-1].vz);
#endif
#endif

  // Free temp memory
  free(buffer);
}

void init_particle_volume(){
  FILE *fp;
  double *buffer;
  double minv = 1e30, maxv, avgv;
  int np, i;
  maxv = avgv = 0.0;

  //========================================================
  //
  // Reads a Zobov volume file. The volumes are in units of
  // boxsize^3. The volume of all particles should sum to 1.
  //
  //========================================================

  // Open volume file
  if ((fp = fopen(global.volfile,"r")) == NULL) {
    printf("Error in init_particle_volume cannot open %s\n",global.volfile);
    exit(0);
  }

  // Read in number of particles
  fread(&np, 1, sizeof(int),fp);

  // Allocate buffer memory
  buffer = (double *) malloc(np*sizeof(double));

  // Check if consistent
  if (np != global.np){
    printf("\n\n*** Particle number in volumefile (%i) does not match (%i) ***\n\n",np, global.np);
    exit(0);
  }

  // Read in volumes
  fread(buffer,np,sizeof(double),fp);
  for (i = 0; i<np; i++){
    allpart[i].volume = buffer[i];
    if (buffer[i] < minv) minv = buffer[i];
    if (buffer[i] > maxv) maxv = buffer[i];
    avgv += buffer[i];
  }
  avgv /= (double)np;

#ifdef VERBOSE
  printf("Min Avg Max Volume: %lf  %lf  %lf\n",minv,avgv,maxv);
  printf("Vol1, Vol2, ... Voln:   %lf  %lf ... %lf\n",buffer[0],buffer[1],buffer[np-1]);
#endif

  // Close file and free temp memory
  fclose(fp);
  free(buffer);
}

void init_void(){	
  FILE *fp;
  int np, nv, i, j;
  VOID *curvoid;

  //========================================================
  // 
  // Reads a Zobov void-list
  // Allocates the [allvoid] array
  //
  //========================================================

  // The quanitites we input from voidfile
  int itemp, voidn, corep, zonep, nzone, vpart; 
  double cored, zonev, vvol, vdc, vprob;

  // Make sure second line is not too long, if so increase number
  char ctemp[256];

  // Open void file (ascii)
  if ((fp = fopen(global.voidfile,"r")) == NULL){
    printf("Error in init_void cannot open %s\n", global.voidfile);
    exit(0);
  }
  fscanf(fp,"%i %i\n",&np,&nv);

#ifdef VERBOSE
  printf("Reading voidfile from %s. Number of voids: %i\n", global.voidfile, nv);
#endif

  // Read second line of file which is just text
  fgets(ctemp,sizeof(ctemp),fp);

  /* Check if consistent */
  if (np != global.np){
    printf("\n\n*** Particle number in voidfile (%i) does not match (%i) ***\n\n",np, global.np);
    exit(0);
  }

  // Allocate memory 
  allvoid = (VOID *)malloc(nv*sizeof(VOID));

  // Initialize important simulation details
  global.nv = nv;

  // Read void by void
  for (i = 0; i<nv; i++) {
    fscanf(fp,"%i %i %i %lf %lf %i %i %lf %i %lf %lf\n",&itemp, &voidn, &corep, &cored, &zonev, 
        &zonep, &nzone, &vvol, &vpart, &vdc, &vprob);

    // Assign values to the current void
    curvoid                   = &allvoid[voidn];
    curvoid->n                = voidn;
    curvoid->CoreParticle     = corep;
    curvoid->CoreDensity      = cored;
    curvoid->VolumeCentalZone = zonev;
    curvoid->nCentralZone     = zonep;
    curvoid->nZones           = nzone;
    curvoid->VoidVolume       = vvol;
    curvoid->np               = vpart;
    curvoid->VoidDensContrast = vdc;
    curvoid->VoidProb         = vprob;
    curvoid->used             = 0;
    curvoid->childlist		     = (VOID **)malloc(nzone*sizeof(VOID *));
    for(j=0;j<nzone;j++) curvoid->childlist[j] = NULL;
    curvoid->children         = 0;
    curvoid->ischild          = 0;
    curvoid->isinvoidstack    = 0;
    curvoid->Up = curvoid->Down = curvoid->Next = curvoid->Prev = NULL;
    
#ifdef SUPERVERBOSE
    if(i % 1000 == 0)
      printf("Vold  %i | %i %i %lf %lf %i %i %lf %i %lf %lf\n",i, voidn, corep, cored, zonev, zonep, nzone, vvol, vpart, vdc, vprob);
#endif
    // Allocate memory for zonelist
    curvoid->zonelist = (int *)malloc(curvoid->nZones*sizeof(int));

    // Init some values to be computed later
    curvoid->reff = curvoid->x = curvoid->y = curvoid->z = 0.0;
    for(j=0;j<6;j++) curvoid->InertiaTensor[j] = 0.0;
  }

  // Init the rootvoid containing all subvoids
  rootvoid.n = -1;
  rootvoid.nZones = global.nv;
  rootvoid.childlist = (VOID **)malloc(global.nv*sizeof(VOID *));
  for(j=0;j<global.nv;j++)
    rootvoid.childlist[j] = NULL;
  rootvoid.children = rootvoid.ischild = rootvoid.isinvoidstack = 0;
  rootvoid.Up = rootvoid.Down = rootvoid.Next = rootvoid.Prev = NULL;

  fclose(fp);
}

void init_particles_in_zones(){
  FILE *fp;
  int np, nv, nc, i, j, total, totalv;
  int *buffer;
  VOID *curvoid;

  //========================================================
  // 
  // Reads a Zobov PartInZones file
  // Allocates [particlelist] in Void struct
  // Must be run after init_void has been run!
  //
  //========================================================

  // Open particle-in-zones file
  if ((fp = fopen(global.partinzonefile,"r")) == NULL){
    printf("Error in init_particles_in_zones cannot open %s\n",global.partinzonefile);
    exit(0);
  }

#ifdef VERBOSE
  printf("Reading particles in zones from %s\n", global.partinzonefile);
  printf("Memory needed to store particlelist: %lf GB\n", 4*np/1024./1024./1024.);
#endif

  // Read in total particle number and total void number
  fread(&np,1,sizeof(int),fp);
  fread(&nv,1,sizeof(int),fp);

  // Check if consistent
  if (np != global.np || nv != global.nv){
    printf("\n\n*** Particle number in voidfile (%i) does not match (%i) ***\n\n",np, global.np);
    printf("\n\n*** Void number in voidfile (%i) does not match (%i) ***\n\n",nv, global.nv);
    exit(0);
  }

  // Allocate buffer
  buffer = (int *)malloc(np*sizeof(int));

  // Loop over all voids
  total = totalv = 0;
  for (i=0; i<nv; i++) {
    curvoid = &allvoid[i];

    // Number of particles in void
    fread(&nc,1,sizeof(int),fp);
    curvoid->npl = nc;
    curvoid->particlelist = (int *)malloc(nc*sizeof(int));

    // Update total number of particles in all voids and in voids+subvoids
    total  += nc;
    totalv += curvoid->np;

    // Read in particle numbers
    fread(buffer,nc,sizeof(int),fp);

    // Initialize particle list
    for (j = 0; j<nc; j++) curvoid->particlelist[j] = buffer[j];

#ifdef SUPERVERBOSE
    if(i % 1000 == 0)
      printf("Void %i has %i particles\n", i, nc);
#endif
  }

#ifdef VERBOSE
  printf("Total particles in list: %i   (%i)\n",total,global.np);
  printf("Total particles in voids: %i (void-in-void leads to double counting)\n",totalv);
#endif

  // Close file and free up temp memory
  fclose(fp);
  free(buffer);
}

void init_zones_in_voids(){
  FILE *fp;
  int nv, nz, i, j, *buffer;
  VOID *curvoid;

  //========================================================
  //
  // Reads a Zobov ZoneInVoid file
  // Allocates [zonelist] in Void struct
  // Must be run after init_void has been run!
  //
  //========================================================

  // Open zones-in-void file
  if ((fp = fopen(global.zoneinvoidfile,"r")) == NULL){
    printf("Error in init_zones_in_voids cannot open %s\n",global.zoneinvoidfile);
    exit(0);
  }

#ifdef VERBOSE
  printf("Reading zones in voids from %s\n", global.zoneinvoidfile);
#endif

  // Read in total void number
  fread(&nv,1,sizeof(int),fp);

  // Check for consistency
  if (global.nv != nv) {
    printf("\n\n*** Number of voids does not match global void: %i != %i ***\n\n",i,nv,global.nv);
    exit(0);
  }

  // Allocate temp memory
  buffer = (int *)malloc(nv*sizeof(int));

  // Read void by void
  for (i=0; i<nv; i++) {
    curvoid = &allvoid[i];

    fread(&nz,1,sizeof(int),fp);
    fread(buffer,nz,sizeof(int),fp);

    // Check for consistency
    if (curvoid->nZones != nz) {
      printf("\n\n*** Number of zones in void %i does not match: %i != %i ***\n\n",i,nz,curvoid->nZones);
      exit(0);
    }

    // Allocate memory for zonelist
    curvoid->zonelist = (int *)malloc(nz*sizeof(int));

    for (j=0; j<nz; j++) curvoid->zonelist[j] = buffer[j];
  }

  // Close file and free temp memory
  fclose(fp);
  free(buffer);
}

#ifdef ADJS
void init_particle_adj(){
  FILE *fp;
  int np, i, j, nout, ind;
  int *buffer;

  //========================================================
  //
  // Reads Zobov adjs file
  // Allocates [adj] in Particle struct
  // Must be run after init_particle_position has been run!
  //
  //========================================================

  // Open zones-in-void file
  if ((fp = fopen(global.adjfile,"r")) == NULL){
    printf("Error in init_particle_adj cannot open %s\n",global.adjfile);
    exit(0);
  }

  // Read in total number of particles
  fread(&np,1,sizeof(int),fp);

  // Check for consistency
  if (global.np != np) {
    printf("\n\n*** Number of particles do not match: %i != %i ***\n\n",i,np,global.np);
    exit(0);
  }

  // Allocate temp memory
  buffer = (int *)malloc(np*sizeof(int));
  for(i=0;i<np;i++) buffer[i] = 0;

  // Read number of adj per particle
  double avg = 0.0;
  for (i=0;i<np;i++){
    fread(&allpart[i].nadj,1,sizeof(int),fp);

    // Allocate memory
    allpart[i].adj = (int *)malloc(allpart[i].nadj*sizeof(int));

    avg += allpart[i].nadj;
  }
  avg /= (double) np;

#ifdef VERBOSE
  printf("Average adjs per particle: %lf Memory needed to store it %lf GB\n",avg, (1.0+avg)*np*4/1024./1024./1024.);
#endif

  // Read adjs
  for (i=0;i<np;i++){
    if (allpart[i].nadj > 0) {

      // Read adjs and set adjs (file written to avoid double-counting)
      fread(&nout,1,sizeof(int),fp);
      for(j=0;j<nout;j++){
        fread(&ind,1,sizeof(int),fp);
        allpart[i].adj[buffer[i]++] = ind;
        allpart[ind].adj[buffer[ind]++] = i;
      }
    }
  }

  // Consistency check
  for (i=0;i<np;i++){
    if(buffer[i] != allpart[i].nadj)
      printf("Error in readadjs for particle %i: %i != %i\n", i, buffer[i], allpart[i].nadj); 
  }

#ifdef VERBOSE
  printf("Particle 0 has adjs:");
  for(i=0;i<allpart[0].nadj;i++){
    printf(" %i", allpart[0].adj[i]);
  }
  printf("\n");
  j=allpart[0].adj[0];

  printf("Particle %i has adjs:",j);
  for(i=0;i<allpart[j].nadj;i++){
    printf(" %i", allpart[j].adj[i]);
  }
  printf("\n");
#endif

  // Close file and clear up temp memory
  fclose(fp);
  free(buffer);
}
#else
  void init_particle_adj(){
    printf("Particle adjs will not be read since ADJS is not defined\n");
  }
#endif

int free_all_memory(){
  int i;

  //========================================================
  //
  // Frees up all memory allocated above
  // Assumes we have allocated everything (does not check)
  //
  //========================================================

  // Free particle data
#ifdef ADJS
  for(i=0;i<global.np;i++)
    free(allpart[i].adj);
#endif
  free(allpart);

  // Free void data
  for(i=0;i<global.nv;i++){
    free(allvoid[i].zonelist);
    free(allvoid[i].childlist);
    free(allvoid[i].particlelist);
  }
  free(allvoid);
}
#endif
