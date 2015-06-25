
//===============================================
// Methods to calculate Bessel functions
// Code modified from 'Numerical Recipies'
//===============================================

#ifdef BESSELLIB
#include<complex_bessel.h>
#endif

struct Bessel {
  static const double cc1[7], cc2[8];

  // Spherical Bessel function J_(n+1/2)	
  double sphbesj(const int n, const double x){
    const double RTPIO2 = 1.253314137315500251;
    double te;

    // Assign value if n = 0
    if (x == 0.0) {
      if (n==0) {
        return 1.0;
      } else {
        return 0.0;
      }
    }

    // Check for bad arguments
    if (n<0 || x < 0.0) throw("Bad arguments in sphbesj: n<0 or x<0");

    // For small values of x use asymptotic form
    if (x < 0.1*sqrt(n)) {
      if (n < 100){
        double gamma = 1.0;
        for (int i = 0; i<=n; i++) {
          gamma *= double(2*i+1);
        }
        return pow(x,n)/gamma;
      }
      return 0.0;
    } else if (n > 200 && x < double(n)/2.0) {
      return 0.0;
    }

#ifdef BESSELLIB
    std::complex<double> w(x,0.0);
    std::complex<double> Jnu = sp_bessel::besselJ(n+0.5, w);
    te = RTPIO2/sqrt(x) * Jnu.real();
#else
    // Calculate function value
    te = RTPIO2/sqrt(x) * besseljy(n+0.5,x);
#endif

    if (te != te) {
      cout << "NaN" << endl;
      cout << "n = " << n << endl;
      cout << "x = " << x << endl;
      exit(1);
    }
    return te;
  }

  // Bessel function J_nu
  double besseljy(const double nu, const double x){
    const int    MAXIT = 1000; //10000
    const double EPS = numeric_limits<double>::epsilon();
    const double FPMIN = numeric_limits<double>::min()/EPS;
    const double XMIN = 2.0;
    const double PI = 3.141592653589793;

    // Bookeeping variables
    double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact;
    double fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2;
    double q,r,rjl,rjl1,rjmu,rjpl,rjtemp;
    double sum,sum1,temp,w,x2,xi,xi2,xmu,xx;
    int i,isign,l,nl;

    // Test if arguments are within bounds
    if (x < 0.0 || nu < 0.0) throw("Bad arguments in besseljy: n < 0 or x<0");

    // Apply (modified) numerical recipies algorithm
    nl   = (x < XMIN ? int(nu+0.5) : MAX(0,int(nu-x+1.5)));
    xmu  = nu-nl;
    xi   = 1.0/x;
    xi2  = 2.0*xi;
    w    = xi2/PI;
    isign = 1;
    h = nu*xi;
    if (h < FPMIN) h=FPMIN;
    b = xi2*nu;
    d = 0.0;
    c = h;

    for (i=0;i<MAXIT;i++) {
      b += xi2;

      d  = b-d;
      if (abs(d) < FPMIN) d=FPMIN;

      c = b-1.0/c;
      if (abs(c) < FPMIN) c=FPMIN;

      d = 1.0/d;
      del = c*d;
      h = del*h;
      if (d < 0.0) isign = -isign;
      if (abs(del-1.0) <= EPS) break;
    }

    char *buf;
    if (i >= MAXIT){
      // xxx Quick fix - does not solve problem xxx
      //buf = new char[200];
      //sprintf(buf, "x too large in besseljy: try asymptotic expansion  x=%g  nu=%g",x,nu);
      //throw(buf);
      return  1.0/x * cos(x - nu*acos(-1.0)/2.0 - acos(-1.0)/4.0);
    }
    rjl  = isign*FPMIN;
    rjpl = h*rjl;
    rjl1 = rjl;
    fact = nu*xi;

    for (l=nl-1;l>=0;l--) {
      rjtemp = fact*rjl+rjpl;
      fact  -= xi;
      rjpl   = fact*rjtemp-rjl;
      rjl    = rjtemp;
    }

    if (rjl == 0.0) rjl=EPS;
    f = rjpl/rjl;

    if (x < XMIN) {
      x2 = 0.5*x;
      pimu = PI*xmu;
      fact = (abs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
      d = -log(x2);
      e = xmu*d;
      fact2 = (abs(e) < EPS ? 1.0 : sinh(e)/e);
      xx = 8.0*SQR(xmu)-1.0;
      gam1 = chebev(cc1,7,xx);
      gam2 = chebev(cc2,8,xx);
      gampl = gam2-xmu*gam1;
      gammi = gam2+xmu*gam1;
      ff = 2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
      e = exp(e);
      p = e/(gampl*PI);
      q = 1.0/(e*PI*gammi);
      pimu2 = 0.5*pimu;
      fact3 = (abs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
      r = PI*pimu2*fact3*fact3;
      d = -x2*x2;
      c = 1.0;
      sum = ff+r*q;
      sum1 = p;
      for (i=1;i<=MAXIT;i++) {
        ff=(i*ff+p+q)/(i*i-xmu*xmu);
        c *= (d/i);
        p /= (i-xmu);
        q /= (i+xmu);
        del =c*(ff+r*q);
        sum += del;
        del1 = c*p-i*del;
        sum1 += del1;
        if (abs(del) < (1.0+abs(sum))*EPS) break;
      }

      if (i > MAXIT) throw("Error: Bessy series failed to converge");
      rjmu=w/(-xmu*xi*sum+sum1*xi2+f*sum);
    } else {
      a=0.25-xmu*xmu;
      p = -0.5*xi;
      q=1.0;
      br=2.0*x;
      bi=2.0;

      fact=a*xi/(p*p+q*q);
      cr=br+q*fact;
      ci=bi+p*fact;

      den=br*br+bi*bi;
      dr=br/den;
      di = -bi/den;
      dlr=cr*dr-ci*di;
      dli=cr*di+ci*dr;

      temp=p*dlr-q*dli;
      q=p*dli+q*dlr;
      p=temp;

      for (i=1;i<MAXIT;i++) {
        a += 2*i;
        bi += 2.0;
        dr=a*dr+br;
        di=a*di+bi;
        if (abs(dr)+abs(di) < FPMIN) dr=FPMIN;

        fact=a/(cr*cr+ci*ci);
        cr=br+cr*fact;
        ci=bi-ci*fact;
        if (abs(cr)+abs(ci) < FPMIN) cr=FPMIN;

        den=dr*dr+di*di;
        dr /= den;
        di /= -den;

        dlr=cr*dr-ci*di;
        dli=cr*di+ci*dr;

        temp=p*dlr-q*dli;
        q=p*dli+q*dlr;
        p=temp;
        if (abs(dlr-1.0)+abs(dli) <= EPS) break;
      }

      if (i >= MAXIT) throw("Error in besseljy: cf2 failed");
      rjmu=sqrt(w/((p-f)*(p-f)/q+q)); 
      rjmu=SIGN(rjmu,rjl);        
    }
    return rjl1*(rjmu/rjl);
  }

  inline double chebev(const double *c, const int m, const double x){
    double d = 0.0;
    double dd = 0.0;
    double sv;
    int j;
    for (j=m-1; j>0; j--) {
      sv = d;
      d = 2.0*x*d-dd+c[j];
      dd = sv;
    }
    return x*d-dd+0.5*c[0];
  }

};

const double Bessel::cc1[] = {-1.142022680371168,6.5165112670737e-3, 3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,3.67795e-11,-1.356e-13};
const double Bessel::cc2[] = {1.843740587300905,-7.68528408447867e-2, 1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8, 2.423096e-10,-1.702e-13,-1.49e-15};


