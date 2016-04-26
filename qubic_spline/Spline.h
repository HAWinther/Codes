#ifndef _SPLINEHEADER_
#define _SPLINEHEADER_
#include <iostream> 
#include <iomanip> 
#include <vector>   
#include <math.h>   

enum ArraySpacing{
  _LINEAR_,
  _LOG_,
  _ARBITRARY_
};

const std::string arrayspacing_type[] = {"linear", "log", "arbitrary"};

template<class T>
class QubicSpline {

  //====================================================
  //
  // Class to construct a Qubic Spline of a set of 
  // distinct (x,y) points.
  //
  // Intended for T = {float, double} 
  //
  // Assumes monotone x-array: 
  //    x_min = x[0] < x[1] < ... < x[n-1] = x_max
  // This is checked in check_arrays()
  //
  // If x-array is regulary spaced we can use direct lookup
  // to save time. Implemented checks for:
  //   type = 0 (_ARBITRARY_) : arbritrary (binary search used)
  //   type = 1 (_LINEAR_)    : linear spacing
  //   type = 2 (_LOG_)       : logaritmic spacing
  // Automatically checked and set in check_arrays()
  //
  // Boundary conditions:
  //   dfdx1, dfdxn are dy/dx at the two boundary points
  //   Use >= 1e30 to get the so-called natural spline
  //   Natural spline is assumed if the bc is not provided
  //
  // Typedefs:
  // QubicSpline<double> == DSpline
  // QubicSpline<float>  == FSpline
  //
  // If x is not in [xmin,xmax] then spline will return
  // the closest value (i.e. f(xmin) / f(xmax)) instead
  // To get error-messages if this happens call
  // show_warning(true)
  //
  //=====================================================

  private:
 
    int n;              // Number of points in arrays
    int type;           // x-array type (see above)

    std::string name;   // Description of spline 

    T x_min;            // x-value at left endpoint
    T x_max;            // x-value at right endpoint
    T dfdx1;            // Derivative at left endpoint
    T dfdxn;            // Derivative at right endpoint
    
    bool range_warning; // Display warning when using spline 
                        // if x is not in range [xmin,xmax]

    std::vector<T> x;   // x values
    std::vector<T> y;   // y values
    std::vector<T> y2;  // Second derivatives

  public:

    //=====================================================
    // Default constructor
    //=====================================================
    QubicSpline() : 
          n(0), 
          name(""),
          x_min(0.0), 
          x_max(0.0), 
          dfdx1(0.0),
          dfdxn(0.0),
          range_warning(false),
          x(0),
          y(0),
          y2(0) {}

    //=====================================================
    // Initialized x,y from T-array
    //=====================================================
    QubicSpline(T *_x, T *_y, int _n, T _dfdx1, T _dfdxn, std::string _name) : 
          n(_n), 
          type(_ARBITRARY_), 
          name(_name), 
          x_min(_x[0]), 
          x_max(_x[n-1]), 
          dfdx1(_dfdx1), 
          dfdxn(_dfdxn),
          range_warning(false),
          x(std::vector<T>(_x, _x + _n)), 
          y(std::vector<T>(_y, _y + _n)), 
          y2(std::vector<T>(_n, 0.0)){
            check_arrays();
            create_spline();
          }

    //=====================================================
    // Initialize x,y from T-vector
    //=====================================================
    QubicSpline(std::vector<T> &_x, std::vector<T> &_y, T _dfdx1, T _dfdxn, std::string _name) : 
          n(_x.size()), 
          type(_ARBITRARY_), 
          name(_name), 
          x_min(_x[0]), 
          x_max(_x[n-1]), 
          dfdx1(_dfdx1), 
          dfdxn(_dfdxn),
          range_warning(false),
          x(_x), 
          y(_y), 
          y2(std::vector<T>(n, 0.0)){
            check_arrays();
            create_spline();
          }
    
    //=====================================================
    // ...if no BC is specified assume natural spline
    //=====================================================
    QubicSpline(T *_x, T *_y, int _n, std::string _name) : QubicSpline(_x, _y, _n, 1e30, 1e30, _name) {} 
    QubicSpline(std::vector<T> &_x, std::vector<T> &_y, std::string _name) : QubicSpline(_x, _y, 1e30, 1e30, _name) {}
    
    //=====================================================
    // Extract data
    //=====================================================
    int get_n() { return n; }
    std::string get_name() { return name; }
    std::string get_type() { return arrayspacing_type[type]; }
    std::vector<T>& get_x()  { return x;  }
    std::vector<T>& get_y()  { return y;  }
    std::vector<T>& get_y2() { return y2; }

    //=====================================================
    // Checks that x-array is stricktly increasing and
    // that x and y has the correct size. Determines if
    // x is log or linear spaced
    //=====================================================
    void check_arrays();
   
    //=====================================================
    // Turn on/off showing errors when x is out of range
    // when using spline
    //=====================================================
    void show_warning(bool v) { range_warning = v; }

    //=====================================================
    // Assignment operator to allow for 'myspline(x)' useage
    //=====================================================
    T operator()(const T& _x){ return f(_x); }
    
    //=====================================================
    // Calculate index klo such that x[klo] <= x0 < x[klo+1]
    //=====================================================
    int lower_index(const T& _x);

    //=====================================================
    // Calculate the y2 array
    //=====================================================
    void create_spline();

    //=====================================================
    // Extract function value from the spline.
    // If _x is outside range return the closest value.
    //=====================================================
    T f(const T& _x);

    //=====================================================
    // Extract the derivative of the splined function.
    // If x0 is outside range return the closest value.
    //=====================================================
    T dfdx(const T& _x);

    //=====================================================
    // Checks if x is in [xmin,xmax]. Outputs error if not
    //=====================================================
    bool is_x_in_range(const T& _x);
    
    //=====================================================
    // Output info about the spline
    //=====================================================
    void info();

    //=====================================================
    // Delete containts of arrays and reset parameters
    //=====================================================
    void clean();
};

//=====================================================
// Some typedefs for easier use
//=====================================================
typedef class QubicSpline<double>  DSpline;
typedef class QubicSpline<float>   FSpline;

#endif

