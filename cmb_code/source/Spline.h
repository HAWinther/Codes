
//=========================================
// Qubic spline class
//=========================================

#pragma once

#include "Matrix.h"

class Spline {
  private:
    VecDoub x;

    // 1D spline from array
    int n;
    VecDoub y, y2;

    // 1D splines from matrix
    bool first;
    double nx, ny;
    MatDoub yy, yy2;

    // Optional parameters for faster look-up in the get_spline method
    bool linear_x, quadratic_x, three_regions_x, log_x, four_regions_x;
    double x_start, x_end;

    int np1, np2, np3, np4;
    double x_start_reg2, x_start_reg3, x_start_reg4;

    double yp1, ypn;

  public:
    // Constructors
    Spline();
    Spline(VecDoub &x1, VecDoub &y1, double yp1_in, double ypn_in);
    Spline(VecDoub &x1, MatDoub &yy1, double yp1_in, double ypn_in, bool first1);
    ~Spline(){};

    // Methods for making spline
    void make_spline();
    void make_spline_matrix(int j);

    // Method for spline interpolation
    double get_spline(double x0);
    double get_spline_matrix(double x0, int j);

    // Test function
    double get_spline_matrix_first_linear(double x0, int j);

    // Let spline object know if the x-array is lineary spaced
    void set_linear_x();

    // Let spline object know if the x-array is quadratic spaced
    void set_quadratic_x();

    // Log spaced
    void set_log_x();

    // Let spline object know if the x-array is linear in three different regions
    void set_three_regions_x(int np_1, int np_2);

    // Let spline object know if the x-array is linear in four different regions
    void set_four_regions_x(int np_1, int np_2, int np_3);

    // Extract arrays
    void extract_y2(VecDoub_O &y1);
    void extract_yy2(MatDoub_O &yy1);

    // Set the yy and yy2 arrays externally (from file)
    void set_external_spline(VecDoub_I &x, MatDoub_I &yy_in, MatDoub_I &yy2_in, bool first1);
};

