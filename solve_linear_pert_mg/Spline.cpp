#include "Spline.h" 

//=====================================================
// Calculate index klow such that x[klow] <= _x < x[klow+1]
//=====================================================
template<class T>
int QubicSpline<T>::lower_index(const T& _x){
  int klow, khigh, k;
  if (type == _LINEAR_) {
    klow = std::min(int((_x - x_min)/(x_max-x_min)*(n-1)),n-2);
    khigh = klow + 1;
  } else if (type == _LOG_) {
    klow = std::min(int((log(_x/x_min))/log(x_max/x_min)*(n-1)),n-2);
    khigh = klow + 1;
  } else {
    klow = 0;
    khigh = n-1;
    while (khigh-klow>1) {
      k = (khigh+klow) >> 1;
      if (x[k] > _x) {
        khigh = k;
      } else {
        klow = k;
      }
    }
  }
  return klow;
}

//=====================================================
// Calculate the spline
//=====================================================
template<class T>
void QubicSpline<T>::create_spline(){
  T t, p;
  std::vector<T> u(n+1, 0.0);

  // Boundary conditions for the spline at left end
  if (dfdx1 < 1e30){
    u[0]  = (3.0/(x[1] - x[0]))*((y[1] - y[0])/(x[1] - x[0]) - dfdx1);
    y2[0] = -0.5;
  }
  
  // Create spline by solving recurence relation (see Wiki::spline_interpolation)
  for (int i = 1; i < n-1; i++){
    t    = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
    p    = t * y2[i-1] + 2.0;
    y2[i]= (t - 1.0)/p;
    u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
    u[i] = (6.0 * u[i]/(x[i+1] - x[i-1]) - t * u[i-1])/p;
  }
  
  // Boundary condition for the spline at right end
  if (dfdxn < 1e30){
    u[n]    = (3.0/(x[n-1] - x[n-2]))*(dfdxn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    y2[n-1] = (u[n] - 0.5 * u[n-2])/(0.5 * y2[n-2] + 1.0);
  }

  // Calculate y''
  for (int i = n-2; i >= 0; i--) 
    y2[i] = y2[i] * y2[i+1] + u[i];
}

//=====================================================
// Extract function value from the spline.
// If _x is outside range return the closest value.
//=====================================================
template<class T> 
T QubicSpline<T>::f(const T& _x){
  int klow, khigh;
  T h, b, a, result;

  // Check and show error if x is outside range
  if(range_warning)
    is_x_in_range(_x);
  
  // If outside of range return closest value
  if(_x <= x_min) return y[0];
  if(_x >= x_max) return y[n-1];

  // Calculate x[klow] < _x < x[khigh]
  klow = lower_index(_x);
  khigh = std::min(klow+1,n-1);
  h = x[khigh]-x[klow];

  // Calculate interpolation value
  a = (x[khigh] - _x)/h;
  b = (_x - x[klow])/h;
  result = (a * y[klow] + b * y[khigh] + ((a * a * a - a) * y2[klow] + (b * b * b - b) * y2[khigh]) * (h * h)/6.0);
  return result;
}

//=====================================================
// Extract the derivative of the splined function.
// If _x is outside range return the closest value.
//=====================================================
template<class T>
T QubicSpline<T>::dfdx(const T& _x){
  int klow, khigh;
  T h, b, a, result;

  // Check and show error if x is outside range
  if(range_warning)
    is_x_in_range(_x);

  // If outside of range return closest value
  if(_x <= x_min) return dfdx1;
  if(_x >= x_max) return dfdxn;

  // Calculate x[klow] < _x < x[khigh]
  klow = lower_index(_x);
  khigh = std::min(klow + 1, n - 1);
  h = x[khigh] - x[klow];

  // Calculate interpolation value
  a = (x[khigh] - _x)/h;
  b = (_x - x[klow])/h;
  result = (y[khigh] - y[klow])/h + h/6.0 * (-(3.0 * a * a - 1.0) * y2[klow] + (3.0 * b * b - 1.0) * y2[khigh]);
  return result;
}

//=====================================================
// Clean up everything
//=====================================================
template<class T> 
void QubicSpline<T>::clean(){
  n = 0;
  type = _ARBITRARY_;
  name = "";
  x_min = x_max = 0.0;
  dfdx1 = dfdxn = 0.0;
  x.clear();
  y.clear();
  y2.clear();
}

//=====================================================
// Checks that x-array is stricktly increasing and
// that x and y has the correct size
//=====================================================
template<class T> 
void QubicSpline<T>::check_arrays(){
  // Allowed fractional difference error in checks of being log/linear spaced
  const T eps = 1e-6;
  
  // Check size is correct
  if(x.size() != y.size()){
    std::cout << "Error [QubicSpline]->[check_arrays]: the size of x and y do not match ";
    std::cout << x.size() << " != " << y.size() << " Name: [" << name << "]" << std::endl;
  }

  // Check that x is monotone
  for(int i = 1; i < n; i++)
    if(x[i] <= x[i-1])
      std::cout << "Error [QubicSpline]->[check_arrays]: the x-array is not stricktly increasing. Name: [" << name << "]" << std::endl;
  
  // Check if x is log-spaced
  bool islog = true;
  if(x_min <= 0.0 && x_max >= 0.0) {
    islog = false;
  } else {
    T ratio = x[1]/x[0];
    for(int i = 2; i < n; i++){
      T r = fabs(x[i]/x[i-1]);
      if(fabs(r - ratio) > fabs(ratio)*eps)
        islog = false;
    }
  }

  // Check that x is lineary-spaced
  bool islinear = true;
  T diff = x[1]-x[0];
  for(int i = 2; i < n; i++){
    T d = x[i]-x[i-1];
    if(fabs(d - diff) > fabs(diff)*eps)
      islinear = false;
  }

  // Set array type
  if(islog) type = _LOG_;
  if(islinear) type = _LINEAR_;
}

//=====================================================
// Check that x is in [xmin,xmax]
//=====================================================
template<class T> 
bool QubicSpline<T>::is_x_in_range(const T& _x){
  if(_x < x_min || _x > x_max){
    std::cout.precision(3);
    std::cout << "Warning [QubicSpline]->[is_x_in_range]: x = " << std::setw(8) << _x;
    std::cout << " is not in range [" << std::setw(8) << x_min << "," << std::setw(8) << x_max << "] Name: [" << name << "] ";
    if(fabs(_x - x_min) < fabs(_x-x_max))
      std::cout << "Returning f(xmin) = " << y[0] << std::endl;
    else  
      std::cout << "Returning f(xmax) = " << y[n-1] << std::endl;
    return false;
    std::cout.precision(0);
  }
  return true;
}

//=====================================================
// Output info about the spline
//=====================================================
template<class T> 
void QubicSpline<T>::info(){
  std::string T_type = "unknown";
  if(std::is_same<T, double>::value) T_type = "double";
  if(std::is_same<T, float>::value)  T_type = "float";

  std::cout.precision(3);
  std::cout << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << "QubicSpline [" << name << "] of type [" << T_type << "] has n = " << n << " points" << std::endl;
  std::cout << "The x-array spacing is [" << get_type() << "]" << std::endl;
  std::cout << "The x-range is       xmin = " << std::setw(8) << x_min << " xmax  = " << std::setw(8) << x_max << std::endl;
  std::cout << "Boundary conditions dfdx1 = " << std::setw(8) << dfdx1 << " dfdxn = " << std::setw(8) << dfdxn << std::endl;
  if(dfdx1 >= 1e30 && dfdxn >= 1e30) 
    std::cout << "=> This correponds to the natural spline" << std::endl;
  if(range_warning)
    std::cout << "Warning is shown if x is out of bounds" << std::endl;
  else
    std::cout << "Warning is not shown if x is out of bounds" << std::endl;
  std::cout << "If x is out of bounds then closest value in [xmin,xmax] is used" << std::endl;
  std::cout << "===============================================================" << std::endl;
  std::cout << std::endl;
  std::cout.precision(0);
}

//=====================================================
// Explicit template specialisation
//=====================================================
template class QubicSpline<double>;
template class QubicSpline<float>;

