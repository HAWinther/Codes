#ifndef ANALYSIS_H_INC
#define ANALYSIS_H_INC
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "global.h"

int  sort_void(VOID *a, VOID *b);
int  sort_voidp(const void *aa, const void *bb);
int  check_if_subvoid(VOID* large, VOID* small);
int  check_if_void_is_good(VOID *curvoid);
void set_used(VOID *curvoid);
void calc_effective_void_radius();
void calc_barycenter_void(VOID *curvoid);
void calc_barycenter_loop_subvoids(VOID *v);
void calc_void_tree();
void calc_barycenter_all();
void calc_barycenter_void(VOID *curvoid);
void calc_barycenter_loop_subvoids(VOID *curvoid);
void find_usable_voids_in_tree_serial(VOID *thisvoid);
void find_usable_voids_in_tree_parallel();
void find_usable_voids_in_tree();
void calc_eigenvalues_3x3_matrix(double *M, double *Eig);
void compute_void_shapes_single(VOID *curvoid);

void calc_effective_void_radius(){
  VOID *curvoid;
  double avgpartvol, voidvol, maxrad, minrad;
  int i;
 
  //========================================================
  //
  // Computes the effective radius of the void from the
  // Zobov void-volume using VolumeVoid = 4pi/3 Reff^3
  //
  // Reff is stored in units of the boxsize
  //
  //========================================================

  // Average particle volume in units of boxsize^3
  avgpartvol = 1.0/(double)global.np;

  // Loop over all voids
  maxrad = -1e30; minrad = 1e30;
  for (i=0; i<global.nv; i++) {
    curvoid = &allvoid[i];

    // Volume of void
    voidvol = curvoid->VoidVolume * avgpartvol;

    // Reff in units of the boxsize
    curvoid->reff = pow(voidvol*3.0/(4.0*M_PI),0.3333);

    if(curvoid->reff > maxrad) maxrad = curvoid->reff;
    if(curvoid->reff < minrad) minrad = curvoid->reff;

#ifdef SUPERVERBOSE
    if(i % 1000 == 0)
      printf("Effective radius of void %i in Mpc/h:  %lf\n", i, curvoid->reff*global.boxsize);
#endif

    if(curvoid->reff > 1.0) {
      printf("Error effective radius of void %i has radius %lf larger than boxsize", i, curvoid->reff);
      exit(1);
    }
  }

#ifdef VERBOSE
  printf("The min/max void radius we have is %lf  - %lf Mpc/h\n", minrad*global.boxsize, maxrad*global.boxsize);
#endif
}

int sort_void(VOID *a, VOID *b){

  //========================================================
  //
  // Sorting function for qsort
  // Sort voids by number of subvoids
  //
  //========================================================
  
  return (a->nZones < b->nZones);
}

int sort_voidp(const void *aa, const void *bb){

  //========================================================
  //
  // Sorting function for qsort
  // Sort an array of pointers to voids by zones without
  // sorting the acctual allvoid array
  //
  //========================================================

  const VOID *a = *(VOID **)aa;	
  const VOID *b = *(VOID **)bb;	
  return (a->nZones < b->nZones);
}

int check_if_subvoid(VOID* large, VOID* small){
  int nl, ns, nlcur, nscur, i, j;	

  //========================================================
  //
  // Checks if void "small" is a subvoid of void "large"
  // by checking if all of the zones in "small" is also
  // in "large"
  //
  //========================================================

  nl = large->nZones;
  ns = small->nZones;

  if(nl == 1){
    printf("Error in check_if_subvoid: Large void has only 1 zone:  %i >= %i\n", nl, ns);
    return 0;
  }

  if(nl < ns){
    printf("Error in check_if_subvoid: Large void is smaller than small void\n");
    return 0;
  }

  // Loop over small's zonelist and check if they are all in larges zonelist 
  for(i=0;i<ns;i++){	
    nscur = small->zonelist[i];
    for(j=1;j<nl;j++){
      nlcur = large->zonelist[j];
      if(nlcur == nscur) j = 2*nl;
    }

    // The current small-zone is not in large's zonelist
    if(j != 2*nl+1) return 0;
  }

#ifdef SUPERVERBOSE
  printf("%i is a subvoid of %i\n",small->n,large->n);
#endif

  // All zones in small is in large	
  return 1;
}

void calc_void_tree(){
  int i, j, k, assigned, ncur, n, children, ischild, vid;
  VOID *curvoid, **v;
  FILE *fp, *fpout;
#ifdef USEOPENMP
  int id, nthreads = global.nthreads;
  int assignedtemp[nthreads];
#endif
  n = global.nv;
  v = (VOID **)malloc(n*sizeof(VOID *));

  //========================================================
  // 
  // This routine makes a void tree by sorting voids
  // according to their subvoids. This can take a while
  // to compute to we have the option to read it from file
  //
  //========================================================

  // Check if the void tree exists in file
  if((fp = fopen(global.voidtreefile,"r")) != NULL){
#ifdef VERBOSE
    printf("Reading void-tree from file\n");
#endif
    fscanf(fp,"%i\n", &n);
    if(n != global.nv){
      printf("Error in void-tree [%s] The number of voids do not match %i !=%i\n",global.voidtreefile, n, global.nv);
      exit(1);
    }
    for(i=0;i<global.nv;i++){
      curvoid = &allvoid[i];
      fscanf(fp,"%i %i ", &ischild, &children);
      curvoid->ischild = ischild;
      curvoid->children = children;
      if(children > 0)
        for(j=0;j<children;j++){
          fscanf(fp,"%i ",&vid);
          curvoid->childlist[j] = &allvoid[vid];
        }
      fscanf(fp,"\n");
    }
    fclose(fp);
  } else {
    assigned = 0;

    // Initialize void list 
    for(i=0;i<global.nv;i++)
      v[i] = &allvoid[i];

    // Sort all voids by number of zones (nZones)
    qsort(v,n,sizeof(v), (comp) sort_voidp);	

    // Loop from smallest to largest void (in terms of zones)
#ifdef USEOPENMP
    for(id=0;id<nthreads;id++) assignedtemp[id] = 0;
#pragma omp parallel for default(shared) private(curvoid,ncur,j,k,id)
#endif
    for(i=n-1;i>0;i--){
#ifdef USEOPENMP
      id = omp_get_thread_num();
#endif

      // Take out void
      curvoid = v[i];
      ncur = curvoid->nZones;

#ifdef SUPERVERBOSE
      printf("Doing void %i (nZones = %i)\n", i, curvoid->nZones);
      fflush(stdout);
#endif
      // Loop over all the other voids and assign its parent
      for(j=i-1;j>=0;j--){
        if(v[j]->nZones > ncur){
          if(check_if_subvoid(v[j],curvoid) == 1){
            // Parent found, set it up
            curvoid->ischild = 1;
#ifdef USEOPENMP
            assignedtemp[id]++;
#else
            assigned++;
#endif
            // Make sure 'v' isn't written to simultainiously
#pragma omp critical
            {
              v[j]->childlist[v[j]->children++] = curvoid;
            }

            // Parent found..abort search
            j = -1;
          }
        }
      }
    }

    // Output voidtree to file...
#ifdef VERBOSE
    printf("Outputing voidtree to file\n");
#endif
    fpout = fopen(global.voidtreefile,"w");
    fprintf(fpout,"%i\n",global.nv);
    for(i=0;i<global.nv;i++){
      curvoid = &allvoid[i];
      children = curvoid->children;

      fprintf(fpout,"%i %i ",curvoid->ischild,children);
      if(children > 0)
        for(j=0;j<children;j++) fprintf(fpout, "%i ", (curvoid->childlist[j])->n);
      fprintf(fpout,"\n");
    }
    fclose(fpout);

#ifdef USEOPENMP
    /* Sum up assigned voids for parallel version... */
    for(id=0;id<nthreads;id++) assigned += assignedtemp[id];
#endif

#ifdef SUPERVERBOSE
    printf("We assigned %i of total %i\n",assigned,global.nv);
    fflush(stdout);
#endif
  }

  // Loop over all the voids and add remaining voids to rootvoids childlist
  k = 0;
  for(i=0;i<global.nv;i++){
    if(allvoid[i].ischild == 0){
      rootvoid.childlist[rootvoid.children] = &allvoid[i];
      rootvoid.children++;
      allvoid[i].ischild = 1;
      k += allvoid[i].children;
    }
  }

#ifdef VERBOSE
  printf("The rootvoid has %i children\n",rootvoid.children);
  printf("The next level has %i children\n",k);
#endif

  // Consistency check: No voids can have more children than subvoids!!
  for(i=0;i<n;i++)
    if(allvoid[i].children >= allvoid[i].nZones){
      printf("Error: Void %i has %i children and %i zones\n",i, allvoid[i].children,allvoid[i].nZones);
      exit(1);
    }

  // Make the tree link pointers
  for(i=0;i<=n;i++){
    if(i<n)
      curvoid = &allvoid[i];
    else
      curvoid = &rootvoid;
    children = curvoid->children;

    if(children >= 2){
      // Set next for children
      for(j=0;j<children-1;j++)
        curvoid->childlist[j]->Next = curvoid->childlist[j+1];

      // Set prev for children
      for(j=1;j<children;j++)
        curvoid->childlist[j]->Prev = curvoid->childlist[j-1];

      // Set down
      curvoid->Down = curvoid->childlist[0];

      // Set up for children
      for(j=0;j<children;j++)
        curvoid->childlist[j]->Up = curvoid;

    } else if(children == 1){
      // Set up for child
      curvoid->childlist[0]->Up = curvoid;

      // Set down
      curvoid->Down = curvoid->childlist[0];
    }
  }

  // Clear up temp memory
  free(v);
}

void set_used(VOID *curvoid){
  int i;

  //========================================================
  //
  // Set book-keeping variable "used" in Void stuct to "1"
  // to signal that the void has been selected. Does this
  // for all subvoids also
  //
  //========================================================
  
  curvoid->used = 1;
  if(curvoid->nZones == 1) return;

  for(i=1;i<curvoid->nZones;i++)
    set_used(&allvoid[curvoid->zonelist[i]]);
}

int check_if_void_is_good(VOID *curvoid){

  //========================================================
  //
  // Check if a single void satisfy criterion and if so 
  // mark that it is good and return that we don't need to
  // continue down the tree
  //
  // Requires that we have computed the void-tree
  // 
  //========================================================

  // Condition based on rmin
  if(curvoid->reff > 10.0/global.boxsize){
    // Condition based on rmax
    if(curvoid->reff < 50.0/global.boxsize){ 
      // Condition based on probabillity that void is due to Poisson noise
      if(curvoid->VoidProb <= 0.0027){
#ifdef VERBOSE
        // Check if a void (or subvoid) is reused in the analysis
        if(curvoid->used == 1) 
          printf("WARNING: Void used again...should not happen!\n");
        printf("# New void found: nZones=%6i    reff=%12.4lf   Prob=%12.4lf  (VoidDens=%12.4lf)\n",
            curvoid->nZones,curvoid->reff*global.boxsize,curvoid->VoidProb,curvoid->VoidDensContrast);
        fflush(stdout);
#endif
        // Set used so we don't reuse the void (or it's subvoids) again
        set_used(curvoid);
        curvoid->isinvoidstack = 1;

        return 0;
      }
    }
  }
  return 1;
}

void calc_barycenter_all(){
  int i;

  //========================================================
  // 
  // Computes the barycenter for all the voids 
  //
  //========================================================

#ifdef USEOPENMP
#pragma omp parallel 
  {
    q.used = (int *) malloc(global.nv*sizeof(int)); 
#pragma omp for
    for(i=0;i<global.nv;i++)
      calc_barycenter_void(&allvoid[i]);
    free(q.used);
  }
#else
  for(i=0;i<global.nv;i++)
    calc_barycenter_void(&allvoid[i]);
#endif

#ifdef VERBOSE
  printf("Center x0 y0 z0: %lf %lf %lf\n",allvoid[0].x, allvoid[0].y, allvoid[0].z);
  printf("Center x1 y1 z1: %lf %lf %lf\n",allvoid[1].x, allvoid[1].y, allvoid[1].z);
  printf("Center xn yn zn: %lf %lf %lf\n",allvoid[global.nv-1].x, allvoid[global.nv-1].y, allvoid[global.nv-1].z);
#endif
}

void calc_barycenter_void(VOID *curvoid){
  double xc,yc,zc;
  int i;

  //========================================================
  // 
  // Calculate the barycenter of the void [curvoid]
  // by looping over itself and all subvoids and summing over
  // all particles belonging to those subvoids
  // x = Sum x_i vol_i / Sum vol_i
  //
  // Also computes the inertia-tensor
  //
  //========================================================

  // Initialize all voids
  for(i=0;i<global.nv;i++) q.used[i] = 0;

  // Initialize thread-private struct that is used
  // to store the data we compute
  q.npart = 0;
  q.x = q.y = q.z = q.volume = 0.0;
  q.xcp = allpart[curvoid->CoreParticle].x;
  q.ycp = allpart[curvoid->CoreParticle].y;
  q.zcp = allpart[curvoid->CoreParticle].z;
  
  // Moment of inertia
  q.Ixx = q.Iyy = q.Izz = q.Ixy = q.Iyz = q.Izx = 0.0;

  // Loop over all subvoids
  calc_barycenter_loop_subvoids(curvoid);

  // Calculate center
  xc = q.x/q.volume;
  yc = q.y/q.volume;
  zc = q.z/q.volume;

  // Periodic boundary conditions
  if(xc <  0.0) xc += 1.0;
  if(xc >= 1.0) xc -= 1.0;
  if(yc <  0.0) yc += 1.0;
  if(yc >= 1.0) yc -= 1.0;		
  if(zc <  0.0) zc += 1.0;
  if(zc >= 1.0) zc -= 1.0;

  // Store barycenter
  curvoid->x = xc;
  curvoid->y = yc;
  curvoid->z = zc;

  // Store inertia-tensor 
  curvoid->InertiaTensor[0] = q.Ixx / (double) q.npart;
  curvoid->InertiaTensor[1] = q.Iyy / (double) q.npart;
  curvoid->InertiaTensor[2] = q.Izz / (double) q.npart;
  curvoid->InertiaTensor[3] = q.Ixy / (double) q.npart;
  curvoid->InertiaTensor[4] = q.Iyz / (double) q.npart;
  curvoid->InertiaTensor[5] = q.Izx / (double) q.npart;

  // Compute void shape
  // compute_void_shapes_single(curvoid);

  // Check if voidvolume matches particle volume
  if(fabs(1.0 - q.volume/curvoid->VoidVolume) > 1.0e-3){
    printf("Volume does not match with particle volume: %lf != %lf\n", q.volume, curvoid->VoidVolume);
    exit(1);
  }
}

void compute_void_shapes_single(VOID *curvoid){
  double Eig[3];
  double a, b, c, eps, rell;

  //========================================================
  //
  // Compute void shapes (ellipticity)
  // Assumes the inertia-tensor has been computed, i.e.
  // that calc_barycenter_void has been run
  //
  // For formulas see "First Structure Formation: A 
  // Simulation of Small-Scale Structure at High Redshift"
  // by Hannah Jang-Condell and Lars Hernquist
  //
  //========================================================

  // Calculate eigenvalues
  calc_eigenvalues_3x3_matrix(curvoid->InertiaTensor,Eig);

  // Calculate the shapes
  a = sqrt(5.0/2.0*(Eig[1] + Eig[2] - Eig[0]));
  b = sqrt(5.0/2.0*(Eig[2] + Eig[0] - Eig[1]));
  c = sqrt(5.0/2.0*(Eig[0] + Eig[1] - Eig[2]));

  // Check for zero or complex eigenvalues
  if(a != a || b != b || c != c || a == 0.0 || b == 0.0 || c == 0.0){
    printf("Warning Comp-shape: eigenvalue is 0 or imaginary in void %i: %lf %lf %lf\n",curvoid->n, a, b, c);
    printf("%12.6e  %12.6e  %12.6e\n",curvoid->InertiaTensor[0],curvoid->InertiaTensor[3],curvoid->InertiaTensor[5]);
    printf("%12.6e  %12.6e  %12.6e\n",curvoid->InertiaTensor[3],curvoid->InertiaTensor[1],curvoid->InertiaTensor[4]);
    printf("%12.6e  %12.6e  %12.6e\n",curvoid->InertiaTensor[5],curvoid->InertiaTensor[4],curvoid->InertiaTensor[2]);
    return;
  }
  eps = 1.0 - sqrt(c/a);

  /*
  // Normalize a, b, c
  a /= (double) q.npart;
  b /= (double) q.npart;
  c /= (double) q.npart;
  */

  // Radius
  rell = pow(a*b*c,0.33333);

  //printf("Reff %12.6lf   Rell %12.6lf   Eps %12.6lf\n",curvoid->reff*128.0,rell/curvoid->reff, eps);
}

void calc_barycenter_loop_subvoids(VOID *curvoid){
  int i, j, nv;
  VOID *newvoid;
 
  //========================================================
  // 
  // Recursive nightmare of a function that loops over
  // subvoids and subvoids of subvoids etc. to compute
  // the barycenter. We can also use this routine
  // to compute other stuff (inertia tensor etc.)
  //
  //========================================================

  // Calculate quanity for current void
  if (q.used[curvoid->n] != 1) {
    PARTICLE *curpart;
    double xx, yy, zz, vol, dr2;

    // Set used
    q.used[curvoid->n] = 1;
  
    // Loop over all particles in void
    for (i = 0; i < curvoid->npl; i++) {
      curpart = &allpart[curvoid->particlelist[i]];
      vol = curpart->volume;

      // xcp,yxp,zcp is location of the core particle (need something to measure againgst)
      xx = curpart->x-q.xcp;
      yy = curpart->y-q.ycp;
      zz = curpart->z-q.zcp;

      // Periodic boundary conditions
      if(xx > 0.5) xx = xx-1.0;
      if(yy > 0.5) yy = yy-1.0;
      if(zz > 0.5) zz = zz-1.0;

      if(xx < -0.5) xx = 1.0+xx;
      if(yy < -0.5) yy = 1.0+yy;
      if(zz < -0.5) zz = 1.0+zz;

      // Update positions, volume and number of particles added
      q.x += vol*xx + vol*q.xcp;
      q.y += vol*yy + vol*q.ycp;
      q.z += vol*zz + vol*q.zcp;
  
      // For moment of inertia comp [Instead of particles a points do particles as vol-elements, i.e. multiply by rho = 1/vol ?!]
      dr2 = xx*xx+yy*yy+zz*zz;
      q.Ixx += dr2-xx*xx;
      q.Iyy += dr2-yy*yy;
      q.Izz += dr2-zz*zz;
      q.Ixy -= xx*yy;
      q.Iyz -= yy*zz;
      q.Izx -= zz*xx;

      // Update volume and number of particles added
      q.volume += vol;
      q.npart++;
    }
  }

  // Number of subvoids
  nv = curvoid->nZones;

  // Return if no subvoids
  if(nv == 1) return;

  // Loop over subvoids and recursively call this method
  for (i = 1; i<nv; i++) {
    newvoid = &allvoid[curvoid->zonelist[i]];
    calc_barycenter_loop_subvoids(newvoid);
  }
  return;
}

void find_usable_voids_in_tree_serial(VOID *thisvoid){
  int continue_down, i;
  VOID *curvoid;

  //========================================================
  //
  // Transverses the tree and searches for voids that 
  // satisfy the conditions [rmin,rmax,voidprob] etc.
  // we have defined.
  //
  //========================================================

  curvoid = thisvoid;
  while(1) {
    continue_down = check_if_void_is_good(curvoid);
    if(curvoid->Down != NULL && continue_down == 1)
      find_usable_voids_in_tree_serial(curvoid->Down);
    if(curvoid->Next == NULL)
      break;
    else
      curvoid = curvoid->Next;
  }
}

void find_usable_voids_in_tree_parallel(){
  int nchildren, i, continue_down;
  VOID *curvoid, *startvoid = &rootvoid;

  //========================================================
  //
  // Transverses the tree and searches for voids that 
  // satisfy the conditions [rmin,rmax,voidprob] etc.
  // we have defined
  //
  // This method is a parallelised version of 
  // find_usable_voids_in_tree_serial. Works by giving the 
  // children of the root node of the tree to the threads.
  // This makes sure we get a speedup and don't mess up 
  // anything! IMPORTANT: 'q' needs to be threadprivate!
  //
  //========================================================

  // First get how many children the rootvoid has
  nchildren = startvoid->children;
#ifdef VERBOSE
  printf("Rootvoid has %i children\n",nchildren);
#endif

  // Move down from rootvoid until we have more than one child to allow for paralellisation 
  // Note: We do not check if this void is in the stack (as then there would be only 1 void in the stack)
  while(nchildren == 1){
    startvoid = startvoid->Down;
    nchildren = startvoid->children;
#ifdef VERBOSE
    printf("Next layer has %i children...\n", nchildren);
#endif
  }
  printf("We start parallel version of code. Layer has %i children\n",nchildren);

#ifdef USEOPENMP
#pragma omp parallel for private(curvoid,continue_down)
#endif
  for(i=0;i<nchildren;i++){

    // Continue with the children of the startvoid. Make sure this part is sequential!
    curvoid = startvoid->childlist[i];
    continue_down = check_if_void_is_good(curvoid);

    if(continue_down == 1 && curvoid->Down != NULL)
      find_usable_voids_in_tree_serial(curvoid->Down);
  }
}

void find_usable_voids_in_tree(){
  int i;

  //========================================================
  //
  // After this is done all the viable voids have
  // isinvoidstack = 1 and are ready to be analyzed
  //
  //========================================================

  // Initialize
  for(i=0;i<global.nv;i++)
    allvoid[i].isinvoidstack = 0;

#ifdef USEOPENMP
  find_usable_voids_in_tree_parallel();
#else
  find_usable_voids_in_tree_serial(&rootvoid);
#endif

  // Count how many void we found
  int ntot = 0;
  for(i=0;i<global.nv;i++){
    ntot += allvoid[i].isinvoidstack;
  }
  printf("We found %i (non-overlapping) voids that satisfy the requirements\n", ntot);
}

void calc_eigenvalues_3x3_matrix(double *M, double *Eig){
  double A11, A12, A13;
  double A22, A23, A33;
  double e1, e2, e3;
  double r, phi, p, qq;

  //========================================================
  //
  // Computes the eigenvalues of a 3x3 matrix
  //
  //========================================================

  // Get matrix coefficients
  A11 = M[0]; A22 = M[1];
  A33 = M[2]; A12 = M[3];
  A23 = M[4]; A13 = M[5];

  // Two cases, diagonal or not...
  p = A12*A12 + A13*A13 + A23*A23;
  if(p == 0.0){
    e1 = A11;
    e2 = A22;
    e3 = A33;
  } else {
    qq = (A11+A22+A33)/3.0;
    p = (A11-qq)*(A11-qq) + (A22-qq)*(A22-qq) + (A33-qq)*(A33-qq) + 2.0*p;
    p = sqrt(p/6.0);
    A11 = (A11 - qq)/p;
    A22 = (A22 - qq)/p;
    A33 = (A33 - qq)/p;
    A12 /= p;
    A13 /= p;
    A23 /= p;

    r = (A11*(A22*A33 - A23*A23) - A12*(A12*A33 - A13*A23) + A13*(A12*A23-A22*A13))/2.0;
    if(r <= -1.0){
      phi = M_PI/3.0;
    } else if(r>=1.0){
      phi = 0.0;
    } else {
      phi = acos(r)/3.0;
    }

    // The eigenvalues satisfy e3 <= e2 <= e1
    e1 = qq + 2.0*p*cos(phi);
    e2 = qq + 2.0*p*cos(phi + 2.0*M_PI/3.0);
    e3 = 3.0*qq - e1 - e2;
  }

  // Sort from smallest to largest
  if(e2>e1){
    p = e2;
    e2 = e1;
    e1 = p;
  }
  if(e3>e1){
    p = e3;
    e3 = e1;
    e1 = p;
  }
  if(e3>e2){
    p = e3;
    e3 = e2;
    e2 = p;
  }

  // Send back the eigenvalues sorted from smallest to largest
  Eig[0] = e3;
  Eig[1] = e2;
  Eig[2] = e1;
}

#endif
