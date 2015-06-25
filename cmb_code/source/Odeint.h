
//=================================================
// Code is modified "Numerical Recipies" code
// Methods to integrate ODE's
//=================================================

#pragma once

#include "Matrix.h"

#ifdef USELAPACK
#include <f2c.h>
#include <clapack.h>
#else
#ifdef LUSAFE
#include "LUdcmp-safe.h"
#else
#include "LUdcmp.h"
#endif
#endif

struct Output {
  int kmax;
  int nvar;
  int nsave;
  bool dense;
  int count;
  double x1,x2,xout,dxout;
  VecDoub xsave;
  MatDoub ysave;
  Output() : kmax(-1),dense(false),count(0) {}
  Output(const int nsavee) : kmax(500),nsave(nsavee),count(0),xsave(kmax) {
    dense = nsave > 0 ? true : false;
  }
  void init(const int neqn, const double xlo, const double xhi) {
    nvar=neqn;
    if (kmax == -1) return;
    ysave.resize(nvar,kmax);
    if (dense) {
      x1=xlo;
      x2=xhi;
      xout=x1;
      dxout=(x2-x1)/nsave;
    }
  }
  void resize() {
    int kold=kmax;
    kmax *= 2;
    VecDoub tempvec(xsave);
    xsave.resize(kmax);
    for (int k=0; k<kold; k++)
      xsave(k)=tempvec(k);
    MatDoub tempmat(ysave);
    ysave.resize(nvar,kmax);
    for (int i=0; i<nvar; i++)
      for (int k=0; k<kold; k++)
        ysave(i,k)=tempmat(i,k);
  }
  template <class Stepper>
    void save_dense(Stepper &s, const double xout, const double h) {
      if (count == kmax) resize();
      for (int i=0;i<nvar;i++)
        ysave(i,count)=s.dense_out(i,xout,h);
      xsave(count++)=xout;
    }
  void save(const double x, VecDoub_I &y) {
    if (kmax <= 0) return;
    if (count == kmax) resize();
    for (int i=0;i<nvar;i++)
      ysave(i,count)=y(i);
    xsave(count++)=x;
  }
  template <class Stepper>
    void out(const int nstp,const double x,VecDoub_I &y,Stepper &s,const double h) {
      if (!dense)
        throw("dense output not set in Output!");
      if (nstp == -1) {
        save(x,y);
        xout += dxout;
      } else {
        while ((x-xout)*(x2-x1) > 0.0) {
          save_dense(s,xout,h);
          xout += dxout;
        }
      }
    }
};

struct StepperBase {
  double &x;
  double xold;
  VecDoub &y,&dydx;
  double atol,rtol;
  bool dense;
  double hdid;
  double hnext;
  double EPS;
  int n,neqn;
  VecDoub yout,yerr;
  StepperBase(VecDoub_IO &yy, VecDoub_IO &dydxx, double &xx, const double atoll,
      const double rtoll, bool dens) : x(xx),y(yy), dydx(dydxx), atol(atoll),
  rtol(rtoll), dense(dens), n(y.n_elem), neqn(n), yout(n), yerr(n) {}
};

template<class Stepper>
struct Odeint {
  static const int MAXSTP=50000;
  double EPS;
  int nok;
  int nbad;
  int nvar;
  double x1,x2,hmin;
  bool dense;
  VecDoub y,dydx;
  VecDoub &ystart;
  Output &out;
  typename Stepper::Dtype &derivs;
  Stepper s;
  int nstp;
  double x,h;
  Odeint(VecDoub_IO &ystartt,const double xx1,const double xx2,
      const double atol,const double rtol,const double h1,
      const double hminn,Output &outt,typename Stepper::Dtype &derivss);
  void integrate();
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const double xx1, const double xx2,
    const double atol, const double rtol, const double h1, const double hminn,
    Output &outt,typename Stepper::Dtype &derivss) : nvar(ystartt.n_elem),
  y(nvar),dydx(nvar),ystart(ystartt),x(xx1),nok(0),nbad(0),
  x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
  s(y,dydx,x,atol,rtol,dense) {
    EPS=numeric_limits<double>::epsilon();
    h=SIGN(h1,x2-x1);
    for (int i=0;i<nvar;i++) y(i)=ystart(i);
    out.init(s.neqn,x1,x2);
  }

template<class Stepper>
void Odeint<Stepper>::integrate() {
  derivs(x,y,dydx);
  if (dense)
    out.out(-1,x,y,s,h);
  else
    out.save(x,y);
  for (nstp=0;nstp<MAXSTP;nstp++) {
    if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
      h=x2-x;
    s.step(h,derivs);
    if (s.hdid == h) ++nok; else ++nbad;
    if (dense)
      out.out(nstp,x,y,s,s.hdid);
    else
      out.save(x,y);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (int i=0;i<nvar;i++) ystart(i)=y(i);
      if (out.kmax > 0 && abs(out.xsave(out.count-1)-x2) > 100.0*abs(x2)*EPS)
        out.save(x,y);
      return;
    }
    if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
    h=s.hnext;
  }
  throw("Too many steps in routine Odeint");
}

template <class D>
struct StepperSie : StepperBase {
  typedef D Dtype;
  static const int KMAXX=12,IMAXX=KMAXX+1;
  int k_targ;
  VecInt nseq;
  VecDoub cost;
  MatDoub table;
  MatDoub dfdy;
  VecDoub dfdx;
  double jac_redo;
  bool calcjac;
  double theta;
  MatDoub a;
  int kright;
  MatDoub coeff;
  MatDoub fsave;
  VecDoub dens;
  VecDoub factrl;

  // LAPACK Parameters
  VecInt pivBlas;

  StepperSie(VecDoub &yy, VecDoub &dydxx, double &xx, const double atol,
      const double rtol, bool dens);
  void step(const double htry,D &derivs);
  bool dy(VecDoub &y, const double htot, const int k, VecDoub &yend,
      int &ipt,VecDoub &scale,D &derivs);
  void polyextr(const int k, MatDoub &table, VecDoub &last);
  void prepare_dense(const double h,VecDoub &ysav,VecDoub &scale,
      const int k, double &error);
  double dense_out(const int i,const double x,const double h);
  void dense_interp(const int n, VecDoub &y, const int imit);

  void backsubs(MatDoub_I &atemp, VecDoub_IO &xtemp, VecDoub_I &btemp);
};

  template<class D>
StepperSie<D>::StepperSie(VecDoub &yy, VecDoub &dydxx, double &xx,
    const double atoll,const double rtoll, bool dens)
: StepperBase(yy,dydxx,xx,atoll,rtoll,dens),nseq(IMAXX),cost(IMAXX),
  dfdy(n,n),dfdx(n),calcjac(false),
  a(n,n),coeff(IMAXX,IMAXX),pivBlas(n),
  fsave((IMAXX-1)*(IMAXX+1)/2+2,n),dens((IMAXX+2)*n),factrl(IMAXX) {

    table = MatDoub(KMAXX,n);

    static const double costfunc=1.0,costjac=5.0,costlu=1.0,costsolve=1.0;
    EPS=numeric_limits<double>::epsilon();
    jac_redo=MIN(1.0e-4,rtol);
    theta=2.0*jac_redo;
    nseq(0)=2;
    nseq(1)=3;
    for (int i=2;i<IMAXX;i++)
      nseq(i)=2*nseq(i-2);
    cost(0)=costjac+costlu+nseq(0)*(costfunc+costsolve);
    for (int k=0;k<KMAXX;k++)
      cost(k+1)=cost(k)+(nseq(k+1)-1)*(costfunc+costsolve)+costlu;
    hnext=-1.0e99;
    double logfact=-log10(rtol+atol)*0.6+0.5;
    k_targ=MAX(1,MIN(KMAXX-1,int(logfact)));
    for (int k=0; k<IMAXX; k++) {
      for (int l=0; l<k; l++) {
        double ratio=double(nseq(k))/nseq(l);
        coeff(k,l)=1.0/(ratio-1.0);
      }
    }
    factrl(0)=1.0;
    for (int k=0; k<IMAXX-1; k++)
      factrl(k+1)=(k+1)*factrl(k);
  }
template <class D>
void StepperSie<D>::step(const double htry,D &derivs) {
  const double STEPFAC1=0.6,STEPFAC2=0.93,STEPFAC3=0.1,STEPFAC4=4.0,
        STEPFAC5=0.5,KFAC1=0.7,KFAC2=0.9;
  static bool first_step=true,last_step=false;
  static bool forward,reject=false,prev_reject=false;
  static double errold;
  int i,k;
  double fac,h,hnew,err;
  bool firstk;
  VecDoub hopt(IMAXX),work(IMAXX);
  VecDoub ysav(n),yseq(n);
  VecDoub ymid(n),scale(n);
  work(0)=1.e30;
  h=htry;
  forward = h>0 ? true : false;
  for (i=0;i<n;i++) ysav(i)=y(i);
  if (h != hnext && !first_step) {
    last_step=true;
  }
  if (reject) {
    prev_reject=true;
    last_step=false;
    theta=2.0*jac_redo;
  }
  for (i=0;i<n;i++)
    scale(i)=atol+rtol*abs(y(i));
  reject=false;
  firstk=true;
  hnew=abs(h);
compute_jac:
  if (theta > jac_redo && !calcjac) {
    derivs.jacobian(x,y,dfdx,dfdy);
    calcjac=true;
  }
  while (firstk || reject) {
    h = forward ? hnew : -hnew;
    firstk=false;
    reject=false;
    if (abs(h) <= abs(x)*EPS)
      throw("step size underflow in StepperSie");
    int ipt=-1;
    for (k=0; k<=k_targ+1;k++) {
      bool success=dy(ysav,h,k,yseq,ipt,scale,derivs);
      if (!success) {
        reject=true;
        hnew=abs(h)*STEPFAC5;
        break;
      }
      if (k == 0)
        y=yseq;
      else
        for (i=0;i<n;i++)
          table(k-1,i)=yseq(i);
      if (k != 0) {
        polyextr(k,table,y);
        err=0.0;
        for (i=0;i<n;i++) {
          scale(i)=atol+rtol*abs(ysav(i));
          err+=SQR((y(i)-table(0,i))/scale(i));
        }
        err=sqrt(err/n);
        if (err > 1.0/EPS || (k > 1 && err >= errold)) {
          reject=true;
          hnew=abs(h)*STEPFAC5;
          break;
        }
        errold=max(4.0*err,1.0);
        double expo=1.0/(k+1);
        double facmin=pow(STEPFAC3,expo);
        if (err == 0.0)
          fac=1.0/facmin;
        else {
          fac=STEPFAC2/pow(err/STEPFAC1,expo);
          fac=MAX(facmin/STEPFAC4,MIN(1.0/facmin,fac));
        }
        hopt(k)=abs(h*fac);
        work(k)=cost(k)/hopt(k);
        if ((first_step || last_step) && err <= 1.0)
          break;
        if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
          if (err <= 1.0)
            break;
          else if (err>nseq(k_targ)*nseq(k_targ+1)*4.0) {
            reject=true;
            k_targ=k;
            if (k_targ>1 && work(k-1)<KFAC1*work(k))
              k_targ--;
            hnew=hopt(k_targ);
            break;
          }
        }
        if (k == k_targ) {
          if (err <= 1.0)
            break;
          else if (err>nseq(k+1)*2.0) {
            reject=true;
            if (k_targ>1 && work(k-1)<KFAC1*work(k))
              k_targ--;
            hnew=hopt(k_targ);
            break;
          }
        }
        if (k == k_targ+1) {
          if (err > 1.0) {
            reject=true;
            if (k_targ>1 && work(k_targ-1)<KFAC1*work(k_targ))
              k_targ--;
            hnew=hopt(k_targ);
          }
          break;
        }
      }
    }
    if (reject) {
      prev_reject=true;
      if (!calcjac) {
        theta=2.0*jac_redo;
        goto compute_jac;
      }
    }
  }
  calcjac=false;
  if (dense)
    prepare_dense(h,ysav,scale,k,err);
  xold=x;
  x+=h;
  hdid=h;
  first_step=false;
  int kopt;
  if (k == 1)
    kopt=2;
  else if (k <= k_targ) {
    kopt=k;
    if (work(k-1) < KFAC1*work(k))
      kopt=k-1;
    else if (work(k) < KFAC2*work(k-1))
      kopt=MIN(k+1,KMAXX-1);
  } else {
    kopt=k-1;
    if (k > 2 && work(k-2) < KFAC1*work(k-1))
      kopt=k-2;
    if (work(k) < KFAC2*work(kopt))
      kopt=MIN(k,KMAXX-1);
  }
  if (prev_reject) {
    k_targ=MIN(kopt,k);
    hnew=MIN(abs(h),hopt(k_targ));
    prev_reject=false;
  }
  else {
    if (kopt <= k)
      hnew=hopt(kopt);
    else {
      if (k<k_targ && work(k)<KFAC2*work(k-1))
        hnew=hopt(k)*cost(kopt+1)/cost(k);
      else
        hnew=hopt(k)*cost(kopt)/cost(k);
    }
    k_targ=kopt;
  }
  if (forward)
    hnext=hnew;
  else
    hnext=-hnew;	
}

//==============================================================
// TEST METHOD FOR BACKSUBSITUTÌON AFTER FOUND LU DECOMPOSITION
//==============================================================
template <class D>
void StepperSie<D>::backsubs(MatDoub_I &atemp, VecDoub_IO &xtemp, VecDoub_I &btemp){
  // Backsubstitution algorithm (FORTRAN matrix order (j,i) instead of (i,j))

  int ntemp = xtemp.n_elem;

  VecDoub ytemp(ntemp);

  // Ly = P^-1 b
  ytemp(0) = btemp(pivBlas[0]-1);
  for (int i = 1; i<ntemp; i++) {
    ytemp(i) = btemp(pivBlas[i]-1);
    for (int j = i-1; j>=0; j--) {
      ytemp(i) -= atemp(j,i)*ytemp(j); 
    }
  }

  // Ux = y
  xtemp(n-1) = ytemp(ntemp-1)/atemp(ntemp-1,ntemp-1);
  for (int i = ntemp-2; i>=0; i--) {
    xtemp(i) = ytemp(i);
    for (int j = i+1; j<ntemp-1; j++) {
      xtemp(i) -= atemp(j,i)*xtemp(j); 
    }
    xtemp(i) /= atemp(i,i);
  }
}

template <class D>
bool StepperSie<D>::dy(VecDoub &y,const double htot,const int k,VecDoub &yend, int &ipt,VecDoub &scale,D &derivs) {

  VecDoub del(n), ytemp(n), dytemp(n);

  int nstep = nseq(k);
  double h = htot/nstep;
  char typea = 'T';
  int nrhs = 1, info;

  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++){
      a(i,j) = -dfdy(i,j);
    }
    a(i,i) += 1.0/h;
  }

#ifndef USELAPACK
  LUdcmp alu(a);
#endif

  double xnew=x+h;
  derivs(xnew,y,del);
  for (int i=0;i<n;i++)
    ytemp(i)=y(i);

#ifdef USELAPACK
  dgetrf_(&n, &n, &a.M[0], &n, &pivBlas[0], &info );
  dgetrs_(&typea, &n, &nrhs, &a.M[0], &n, &pivBlas[0], &del[0], &n, &info );
#else
  alu.solve(del,del);
#endif

  // Compute LU:
  // dgetrf_(&N, &N, &a.M[0], &N, &pivBlas[0], &info );
  // Solve system:
  // dgetrs_(&typea[0], &N, &nrhs, &a.M[0], &N, &pivBlas[0], B-array->X-array, &nrhs, &info );

  if (dense && nstep==k+1) {
    ipt++;
    for (int i=0;i<n;i++)
      fsave(ipt,i)=del(i);
  }
  for (int nn=1;nn<nstep;nn++) {
    for (int i=0;i<n;i++)
      ytemp(i) += del(i);
    xnew += h;
    derivs(xnew,ytemp,yend);
    if (nn ==1 && k<=1) {
      double del1=0.0;
      for (int i=0;i<n;i++)
        del1 += SQR(del(i)/scale(i));
      del1=sqrt(del1);
      derivs(x+h,ytemp,dytemp);
      for (int i=0;i<n;i++)
        del(i)=dytemp(i)-del(i)/h;

#ifdef USELAPACK
      dgetrs_(&typea, &n, &nrhs, &a.M[0], &n, &pivBlas[0], &del[0], &n, &info );
#else
      alu.solve(del,del);
#endif

      double del2=0.0;
      for (int i=0;i<n;i++)
        del2 += SQR(del(i)/scale(i));
      del2=sqrt(del2);
      theta=del2/MAX(1.0,del1);
      if (theta > 1.0)
        return false;
    }

#ifdef USELAPACK
    for (int i = 0; i<n; i++) {
      del[i] = yend[i];
    }
    dgetrs_(&typea, &n, &nrhs, &a.M[0], &n, &pivBlas[0], &del[0], &n, &info );
#else
    alu.solve(yend,del);
#endif

    if (dense && nn >= nstep-k-1) {
      ipt++;
      for (int i=0;i<n;i++)
        fsave(ipt,i)=del(i);
    }
  }
  for (int i=0;i<n;i++)
    yend(i)=ytemp(i)+del(i);

  return true;
}
template <class D>
void StepperSie<D>::polyextr(const int k,MatDoub &table,VecDoub &last) {
  int l=last.n_elem;
  for (int j=k-1; j>0; j--)
    for (int i=0; i<l; i++)
      table(j-1,i)=table(j,i)+coeff(k,j)*(table(j,i)-table(j-1,i));
  for (int i=0; i<l; i++)
    last(i)=table(0,i)+coeff(k,0)*(table(0,i)-last(i));
}
template <class D>
void StepperSie<D>::prepare_dense(const double h,VecDoub &ysav,VecDoub &scale, const int k,double &error) {
  kright=k;
  for (int i=0; i<n; i++) {
    dens(i)=ysav(i);
    dens(n+i)=y(i);
  }
  for (int klr=0; klr < kright; klr++) {
    if (klr >= 1) {
      for (int kk=klr; kk<=k; kk++) {
        int lbeg=((kk+3)*kk)/2;
        int lend=lbeg-kk+1;
        for (int l=lbeg; l>=lend; l--)
          for (int i=0; i<n; i++)
            fsave(l,i)=fsave(l,i)-fsave(l-1,i);
      }
    }
    for (int kk=klr; kk<=k; kk++) {
      double facnj=nseq(kk);
      facnj=pow(facnj,klr+1)/factrl(klr+1);
      int ipt=((kk+3)*kk)/2;
      int krn=(kk+2)*n;
      for (int i=0; i<n; i++) {
        dens(krn+i)=fsave(ipt,i)*facnj;
      }
    }
    for (int j=klr+1; j<=k; j++) {
      double dblenj=nseq(j);
      for (int l=j; l>=klr+1; l--) {
        double factor=dblenj/nseq(l-1)-1.0;
        for (int i=0; i<n; i++) {
          int krn=(l+2)*n+i;
          dens(krn-n)=dens(krn)+(dens(krn)-dens(krn-n))/factor;
        }
      }
    }
  }
  for (int in=0; in<n; in++) {
    for (int j=1; j<=kright+1; j++) {
      int ii=n*j+in;
      dens(ii)=dens(ii)-dens(ii-n);
    }
  }
}
template <class D>
double StepperSie<D>::dense_out(const int i,const double x,const double h) {
  double theta=(x-xold)/h;
  int k=kright;
  double yinterp=dens((k+1)*n+i);
  for (int j=1; j<=k; j++)
    yinterp=dens((k+1-j)*n+i)+yinterp*(theta-1.0);
  return dens(i)+yinterp*theta;
}

/*
 *  stepperbs.h
 */

template <class D>
struct StepperBS : StepperBase {
  typedef D Dtype;
  static const int KMAXX=8,IMAXX=KMAXX+1;
  int k_targ;
  VecInt nseq;
  VecInt cost;
  MatDoub table;
  VecDoub dydxnew;
  int mu;
  MatDoub coeff;
  VecDoub errfac;
  MatDoub ysave;
  MatDoub fsave;
  VecInt ipoint;
  VecDoub dens;
  StepperBS(VecDoub_IO &yy, VecDoub_IO &dydxx, double &xx, const double atol,
      const double rtol, bool dens);
  void step(const double htry,D &derivs);
  virtual void dy(VecDoub_I &y, const double htot, const int k, VecDoub_O &yend,
      int &ipt, D &derivs);
  void polyextr(const int k, MatDoub_IO &table, VecDoub_IO &last);
  virtual void prepare_dense(const double h,VecDoub_I &dydxnew, VecDoub_I &ysav,
      VecDoub_I &scale, const int k, double &error);
  virtual double dense_out(const int i,const double x,const double h);
  virtual void dense_interp(const int n, VecDoub_IO &y, const int imit);
};

template <class D>
StepperBS<D>::StepperBS(VecDoub_IO &yy,VecDoub_IO &dydxx,double &xx,
    const double atoll,const double rtoll, bool dens) :
  StepperBase(yy,dydxx,xx,atoll,rtoll,dens),nseq(IMAXX),cost(IMAXX),table(KMAXX,n),coeff(IMAXX,IMAXX),
  dydxnew(n),errfac(2*IMAXX+2),ysave(IMAXX,n),
  fsave(IMAXX*(2*IMAXX+1),n),ipoint(IMAXX+1),dens((2*IMAXX+5)*n) {

    EPS=numeric_limits<double>::epsilon();
    if (dense)
      for (int i=0;i<IMAXX;i++)
        nseq(i)=4*i+2;
    else
      for (int i=0;i<IMAXX;i++)
        nseq(i)=2*(i+1);
    cost(0)=nseq(0)+1;
    for (int k=0;k<KMAXX;k++) cost(k+1)=cost(k)+nseq(k+1);
    hnext=-1.0e99;
    double logfact=-log10(MAX(1.0e-12,rtol))*0.6+0.5;
    k_targ=MAX(1,MIN(KMAXX-1,int(logfact)));

    for (int k = 0; k<IMAXX; k++) {
      for (int l=0; l<k; l++) {
        double ratio=double(nseq(k))/nseq(l);
        coeff(k,l)=1.0/(ratio*ratio-1.0);
      }
    }
    for (int i=0; i<2*IMAXX+1; i++) {
      int ip5=i+5;
      errfac(i)=1.0/(ip5*ip5);
      double e = 0.5*sqrt(double(i+1)/ip5);
      for (int j=0; j<=i; j++) {
        errfac(i) *= e/(j+1);
      }
    }
    ipoint(0)=0;
    for (int i=1; i<=IMAXX; i++) {
      int njadd=4*i-2;
      if (nseq(i-1) > njadd) njadd++;
      ipoint(i)=ipoint(i-1)+njadd;
    }
  }

template <class D>
double StepperBS<D>::dense_out(const int i,const double x,const double h) {
  double theta=(x-xold)/h;
  double theta1=1.0-theta;
  double yinterp=dens(i)+theta*(dens(n+i)+theta1*(dens(2*n+i)*theta
        +dens(3*n+i)*theta1));
  if (mu<0)
    return yinterp;
  double theta05=theta-0.5;
  double t4=SQR(theta*theta1);
  double c=dens(n*(mu+4)+i);
  for (int j=mu;j>0; j--)
    c=dens(n*(j+3)+i)+c*theta05/j;
  yinterp += t4*c;
  return yinterp;
}

template <class D>
void StepperBS<D>::dense_interp(const int n, VecDoub_IO &y, const int imit) {
  double y0,y1,yp0,yp1,ydiff,aspl,bspl,ph0,ph1,ph2,ph3,fac1,fac2;
  VecDoub a(31);
  for (int i=0; i<n; i++) {
    y0=y(i);
    y1=y(2*n+i);
    yp0=y(n+i);
    yp1=y(3*n+i);
    ydiff=y1-y0;
    aspl=-yp1+ydiff;
    bspl=yp0-ydiff;
    y(n+i)=ydiff;
    y(2*n+i)=aspl;
    y(3*n+i)=bspl;
    if (imit < 0) continue;
    ph0=(y0+y1)*0.5+0.125*(aspl+bspl);
    ph1=ydiff+(aspl-bspl)*0.25;
    ph2=-(yp0-yp1);
    ph3=6.0*(bspl-aspl);
    if (imit >= 1) {
      a(1)=16.0*(y(5*n+i)-ph1);
      if (imit >= 3) {
        a(3)=16.0*(y(7*n+i)-ph3+3*a(1));
        for (int im=5; im <=imit; im+=2) {
          fac1=im*(im-1)/2.0;
          fac2=fac1*(im-2)*(im-3)*2.0;
          a(im)=16.0*(y((im+4)*n+i)+fac1*a(im-2)-fac2*a(im-4));
        }
      }
    }
    a(0)=(y(4*n+i)-ph0)*16.0;
    if (imit >= 2) {
      a(2)=(y(n*6+i)-ph2+a(0))*16.0;
      for (int im=4; im <=imit; im+=2) {
        fac1=im*(im-1)/2.0;
        fac2=im*(im-1)*(im-2)*(im-3);
        a(im)=(y(n*(im+4)+i)+a(im-2)*fac1-a(im-4)*fac2)*16.0;
      }
    }
    for (int im=0; im<=imit; im++)
      y(n*(im+4)+i)=a(im);
  }
}


template <class D>
void StepperBS<D>::dy(VecDoub_I &y,const double htot,const int k,VecDoub_O &yend,
    int &ipt,D &derivs) {
  VecDoub ym(n),yn(n);
  int nstep=nseq(k);
  double h=htot/nstep;
  for (int i=0;i<n;i++) {
    ym(i)=y(i);
    yn(i)=y(i)+h*dydx(i);
  }
  double xnew=x+h;
  derivs(xnew,yn,yend);
  double h2=2.0*h;
  for (int nn=1;nn<nstep;nn++) {
    if (dense && nn == nstep/2) {
      for (int i=0;i<n;i++)
        ysave(k,i)=yn(i);
    }
    if (dense && abs(nn-nstep/2) <= 2*k+1) {
      ipt++;
      for (int i=0;i<n;i++)
        fsave(ipt,i)=yend(i);
    }
    for (int i=0;i<n;i++) {
      double swap=ym(i)+h2*yend(i);
      ym(i)=yn(i);
      yn(i)=swap;
    }
    xnew += h;
    derivs(xnew,yn,yend);
  }
  if (dense && nstep/2 <= 2*k+1) {
    ipt++;
    for (int i=0;i<n;i++)
      fsave(ipt,i)=yend(i);
  }
  for (int i=0;i<n;i++)
    yend(i)=0.5*(ym(i)+yn(i)+h*yend(i));
}


template <class D>
void StepperBS<D>::polyextr(const int k, MatDoub_IO &table, VecDoub_IO &last) {
  int l=last.n_elem;
  for (int j=k-1; j>0; j--)
    for (int i=0; i<l; i++)
      table(j-1,i)=table(j,i)+coeff(k,j)*(table(j,i)-table(j-1,i));
  for (int i=0; i<l; i++)
    last(i)=table(0,i)+coeff(k,0)*(table(0,i)-last(i));
}


template <class D>
void StepperBS<D>::prepare_dense(const double h,VecDoub_I &dydxnew,
    VecDoub_I &ysav,VecDoub_I &scale,const int k,double &error) {
  mu=2*k-1;
  for (int i=0; i<n; i++) {
    dens(i)=ysav(i);
    dens(n+i)=h*dydx(i);
    dens(2*n+i)=y(i);
    dens(3*n+i)=dydxnew(i)*h;
  }
  for (int j=1; j<=k; j++) {
    double dblenj=nseq(j);
    for (int l=j; l>=1; l--) {
      double factor=SQR(dblenj/nseq(l-1))-1.0;
      for (int i=0; i<n; i++)
        ysave(l-1,i)=ysave(l,i)+(ysave(l,i)-ysave(l-1,i))/factor;
    }
  }
  for (int i=0; i<n; i++)
    dens(4*n+i)=ysave(0,i);
  for (int kmi=1; kmi<=mu; kmi++) {
    int kbeg=(kmi-1)/2;
    for (int kk=kbeg; kk<=k; kk++) {
      double facnj=pow(nseq(kk)/2.0,kmi-1);
      int ipt=ipoint(kk+1)-2*kk+kmi-3;
      for (int i=0; i<n; i++)
        ysave(kk,i)=fsave(ipt,i)*facnj;
    }
    for (int j=kbeg+1; j<=k; j++) {
      double dblenj=nseq(j);
      for (int l=j; l>=kbeg+1; l--) {
        double factor=SQR(dblenj/nseq(l-1))-1.0;
        for (int i=0; i<n; i++)
          ysave(l-1,i)=ysave(l,i)+
            (ysave(l,i)-ysave(l-1,i))/factor;
      }
    }
    for (int i=0; i<n; i++)
      dens((kmi+4)*n+i)=ysave(kbeg,i)*h;
    if (kmi == mu) continue;
    for (int kk=kmi/2; kk<=k; kk++) {
      int lbeg=ipoint(kk+1)-1;
      int lend=ipoint(kk)+kmi;
      if (kmi == 1) lend += 2;
      for (int l=lbeg; l>=lend; l-=2)
        for (int i=0; i<n; i++)
          fsave(l,i)=fsave(l,i)-fsave(l-2,i);
      if (kmi == 1) {
        int l=lend-2;
        for (int i=0; i<n; i++)
          fsave(l,i)=fsave(l,i)-dydx(i);
      }
    }
    for (int kk=kmi/2; kk<=k; kk++) {
      int lbeg=ipoint(kk+1)-2;
      int lend=ipoint(kk)+kmi+1;
      for (int l=lbeg; l>=lend; l-=2)
        for (int i=0; i<n; i++)
          fsave(l,i)=fsave(l,i)-fsave(l-2,i);
    }
  }
  dense_interp(n,dens,mu);
  error=0.0;
  if (mu >= 1) {
    for (int i=0; i<n; i++)
      error += SQR(dens((mu+4)*n+i)/scale(i));
    error=sqrt(error/n)*errfac(mu-1);
  }
}


template <class D>
void StepperBS<D>::step(const double htry,D &derivs) {
  const double STEPFAC1=0.65,STEPFAC2=0.94,STEPFAC3=0.02,STEPFAC4=4.0,
        KFAC1=0.8,KFAC2=0.9;
  static bool first_step=true,last_step=false;
  static bool forward,reject=false,prev_reject=false;
  int i,k;
  double fac,h,hnew,hopt_int,err;
  bool firstk;
  VecDoub hopt(IMAXX),work(IMAXX);
  VecDoub ysav(n),yseq(n);
  VecDoub ymid(n),scale(n);
  work(0)=0;
  h=htry;
  forward = h>0 ? true : false;
  for (i=0;i<n;i++) ysav(i)=y(i);
  if (h != hnext && !first_step) {
    last_step=true;
  }
  if (reject) {
    prev_reject=true;
    last_step=false;
  }
  reject=false;
  firstk=true;
  hnew=abs(h);
interp_error:
  while (firstk || reject) {
    h = forward ? hnew : -hnew;
    firstk=false;
    reject=false;
    if (abs(h) <= abs(x)*EPS)
      throw("step size underflow in StepperBS");
    int ipt=-1;
    for (k=0; k<=k_targ+1;k++) {
      dy(ysav,h,k,yseq,ipt,derivs);
      if (k == 0)
        y=yseq;
      else
        for (i=0;i<n;i++)
          table(k-1,i)=yseq(i);
      if (k != 0) {
        polyextr(k,table,y);
        err=0.0;
        for (i=0;i<n;i++) {
          scale(i)=atol+rtol*MAX(abs(ysav(i)),abs(y(i)));
          err+=SQR((y(i)-table(0,i))/scale(i));
        }
        err=sqrt(err/n);
        double expo=1.0/(2*k+1);
        double facmin=pow(STEPFAC3,expo);
        if (err == 0.0)
          fac=1.0/facmin;
        else {
          fac=STEPFAC2/pow(err/STEPFAC1,expo);
          fac=MAX(facmin/STEPFAC4,MIN(1.0/facmin,fac));
        }
        hopt(k)=abs(h*fac);
        work(k)=cost(k)/hopt(k);
        if ((first_step || last_step) && err <= 1.0)
          break;
        if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
          if (err <= 1.0)
            break;
          else if (err>SQR(nseq(k_targ)*nseq(k_targ+1)/(nseq(0)*nseq(0)))) {
            reject=true;
            k_targ=k;
            if (k_targ>1 && work(k-1)<KFAC1*work(k))
              k_targ--;
            hnew=hopt(k_targ);
            break;
          }
        }
        if (k == k_targ) {
          if (err <= 1.0)
            break;
          else if (err>SQR(nseq(k+1)/nseq(0))) {
            reject=true;
            if (k_targ>1 && work(k-1)<KFAC1*work(k))
              k_targ--;
            hnew=hopt(k_targ);
            break;
          }
        }
        if (k == k_targ+1) {
          if (err > 1.0) {
            reject=true;
            if (k_targ>1 && work(k_targ-1)<KFAC1*work(k_targ))
              k_targ--;
            hnew=hopt(k_targ);
          }
          break;
        }
      }
    }
    if (reject)
      prev_reject=true;
  }
  derivs(x+h,y,dydxnew);
  if (dense) {
    prepare_dense(h,dydxnew,ysav,scale,k,err);
    hopt_int=h/MAX(pow(err,1.0/(2*k+3)),0.01);
    if (err > 10.0) {
      hnew=abs(hopt_int);
      reject=true;
      prev_reject=true;
      goto interp_error;
    }
  }
  dydx=dydxnew;
  xold=x;
  x+=h;
  hdid=h;
  first_step=false;
  int kopt;
  if (k == 1)
    kopt=2;
  else if (k <= k_targ) {
    kopt=k;
    if (work(k-1) < KFAC1*work(k))
      kopt=k-1;
    else if (work(k) < KFAC2*work(k-1))
      kopt=MIN(k+1,KMAXX-1);
  } else {
    kopt=k-1;
    if (k > 2 && work(k-2) < KFAC1*work(k-1))
      kopt=k-2;
    if (work(k) < KFAC2*work(kopt))
      kopt=MIN(k,KMAXX-1);
  }
  if (prev_reject) {
    k_targ=MIN(kopt,k);
    hnew=MIN(abs(h),hopt(k_targ));
    prev_reject=false;
  }
  else {
    if (kopt <= k)
      hnew=hopt(kopt);
    else {
      if (k<k_targ && work(k)<KFAC2*work(k-1))
        hnew=hopt(k)*cost(kopt+1)/cost(k);
      else
        hnew=hopt(k)*cost(kopt)/cost(k);
    }
    k_targ=kopt;
  }
  if (dense)
    hnew=MIN(hnew,abs(hopt_int));
  if (forward)
    hnext=hnew;
  else
    hnext=-hnew;	
}
