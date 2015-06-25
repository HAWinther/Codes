#ifndef HALOHEADER_INC
#define HALOHEADER_INC
#include "spline.h"

#if defined(DOUBLE)
typedef double realT;
#elif defined(LONGDOUBLE)
typedef long realT realT;
#elif defined(FLOAT)
typedef float realT;
#else
typedef double realT;
#endif

//=======================================================
// A copy of to_string since its not in std
//=======================================================

template <typename T> std::string to_string(T value){
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

//=======================================================
// Halo class
//=======================================================

class Halo{
  public:
    realT x[NDIM], v[NDIM];
    realT rvir;
    realT mass;
    int npart;
    int id;

    // For halo profiles if they are given
    realT *r, **y;
    int nprofile, nbins;
    Spline *yspline;

    Halo(){
      for(int i=0;i<NDIM;i++)
        x[i] = v[i] = 0.0;
      rvir = mass = 0.0;
      npart = id = 0;
      r=NULL; y=NULL;
      nprofile=nbins=0;
    }

    Halo(realT *x, realT rvir, realT mass, int npart, int id){
      for(int i=0;i<NDIM;i++){
        this->x[i] = x[i];
        this->v[i] = 0.0;
      }
      this->rvir = rvir;
      this->mass = mass;
      this->npart = npart;
      this->id = id;
      r=NULL; y=NULL;
      nprofile=nbins=0;
    }

    Halo(realT *x, realT *v, realT rvir, realT mass, int npart, int id){
      for(int i=0;i<NDIM;i++){
        this->x[i] = x[i];
        this->v[i] = v[i];
      }
      this->rvir = rvir;
      this->mass = mass;
      this->npart = npart;
      this->id = id;
      r=NULL; y=NULL;
      nprofile=nbins=0;
    }

    ~Halo(){
      if(nprofile>0){
        delete[] r;
        for(int i=0;i<nprofile;i++){
          delete[] y[i];
        }
        delete[] y;
      }
    }

    // Dump some basic halo information to standard output
    void print_info(){
      std::cout << "============================" << std::endl;
      std::cout << "Halo information id= " << id  << std::endl;
      std::cout << "============================" << std::endl;
      std::cout << "Pos   = "; for(int ii=0;ii<NDIM;ii++) std::cout << " " << x[ii]; std::cout << std::endl;
      std::cout << "Vel   = "; for(int ii=0;ii<NDIM;ii++) std::cout << " " << v[ii]; std::cout << std::endl;
      std::cout << "Mass  = " << mass << std::endl;
      std::cout << "Rvir  = " << rvir << std::endl;
      std::cout << "npart = " << npart << std::endl;
      std::cout << "============================" << std::endl;
      std::cout << std::endl;
    }

    // If profiles are given lets allocate some memory for them
    void allocate_memory_profile(int nbins, int nprofile){
      this->nprofile = nprofile;
      this->nbins = nbins;
      r = new realT[nbins];
      y = new realT*[nprofile];
      for(int i=0;i<nprofile;i++){
        y[i] = new realT[nbins];
      }
    }

    // Store the profiles given to us
    void store_profile_data(realT *data, int j){
      r[j] = fabs(data[0]); // We might want to rescale wrt the virial radius here!
      for(int i=0;i<nprofile;i++){
        y[i][j] = data[i+1];
      }
    }

    // Spline the profiles to make them easier to work with
    void spline_halo_profiles(){
      yspline = new Spline[nprofile];

      for(int i=0;i<nprofile;i++){
        yspline[i].make_spline(r,y[i],nbins,1e30,1e30,0,"Halospline halo#"+to_string(id)+" profile #"+to_string(i));
      }
    }

    // Get the value of profile i at r=rr
    realT get_profile(int i,realT rr){
      return yspline[i].get_spline(rr);
    }
};

#endif
