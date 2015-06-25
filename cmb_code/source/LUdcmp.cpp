
//=======================================================
// Code is modified "Numerical Recipies" code
// LU decomposition of a matrix (safe version)
//=======================================================

#include "LUdcmp.h"

LUdcmp::LUdcmp(MatDoub &A) : a(A), n(A.nx){
  ludcmp();
}

void LUdcmp::ludcmp(){
  int i, j, k, nm1;
  double sum, diag;

  diag = 1.0/a(0,0);
  for (i=1; i<n; i++) a(0,i) *= diag;

  nm1 = n - 1;
  for (j=1; j<nm1; j++) {

    for (i=j; i<n; i++) {
      sum = 0.0;
      for (k=0; k<j; k++) sum += a(i,k)*a(k,j);
      a(i,j) -= sum;
    }
    diag = 1.0/a(j,j);
    for (k=j+1; k<n; k++) {
      sum = 0.0;
      for (i=0; i<j; i++) sum += a(j,i)*a(i,k);
      a(j,k) = (a(j,k)-sum)*diag;
    }
  }

  sum = 0.0;
  for (k=0; k<nm1; k++) sum += a(nm1,k)*a(k,nm1);
  a(nm1,nm1) -= sum;
}

void LUdcmp::solve(VecDoub &b, VecDoub &x){
  int i,j;
  double sum;

  for (i=0; i<n; i++) {
    x[i] = b[i];
  }

  x[0] /= a(0,0);
  for (i=1; i<n; i++) {
    sum = 0.0;
    for (j=0; j<i; j++) sum += a(i,j)*x[j];
    x[i] = (x[i]-sum)/a(i,i);
  }

  for (i=n-2; i>=0; i--) {
    sum = 0.0;
    for (j=i+1; j<n; j++) sum += a(i,j) * x[j];
    x[i] -= sum;
  }
}

void LUdcmp::solve(VecDoub &x){
  int i,j;
  double sum;

  x[0] /= a(0,0);
  for (i=1; i<n; i++) {
    sum = 0.0;
    for (j=0; j<i; j++) sum += a(i,j)*x[j];
    x[i] = (x[i]-sum)/a(i,i);
  }

  for (i=n-2; i>=0; i--) {
    sum = 0.0;
    for (j=i+1; j<n; j++) sum += a(i,j) * x[j];
    x[i] -= sum;
  }
}
