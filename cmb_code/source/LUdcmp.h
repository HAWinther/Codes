#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "Matrix.h"

struct LUdcmp{
  MatDoub a;
  int n;

  LUdcmp(MatDoub &A);

  void ludcmp();
  void solve(VecDoub &b, VecDoub &x);
  void solve(VecDoub &x);

};

