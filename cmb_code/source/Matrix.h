
//=================================
// General array and matrix class
//=================================

#pragma once

#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <String>

using namespace std;

//==============================
// The vector class
//==============================

template<class T>
class SimpleArray {
  private:
    T *x;
  public:
    // Number of elements
    int n_elem;

    // Constructors
    SimpleArray();										// Default constructor
    explicit SimpleArray(int n);						// Zero-based
    SimpleArray(const SimpleArray<T> &m);				// Copy constructor
    ~SimpleArray();										// Destructor
    SimpleArray &operator=(const SimpleArray<T> &m);	// Assignment

    // Make the value of T availible externally
    typedef T value_type;

    // Get the i'th element
    inline T &operator()(int i);
    inline T  operator()(int i) const;
    inline T &operator[](int i);
    inline T  operator[](int i) const;

    // Resize (destroys array)
    void resize(int n);
};

//==============================
// Inline func. for vector class
//==============================

// Constructors
template<class T>
SimpleArray<T>::SimpleArray() : n_elem(0), x(NULL) {}

template<class T>
SimpleArray<T>::SimpleArray(int n) : n_elem(n), x(n>0 ? new T[n] : NULL) {}

template<class T>
SimpleArray<T>::SimpleArray(const SimpleArray<T> &m) : n_elem(m.n_elem), x(n_elem>0 ? new T[n_elem] : NULL) {
  for (int i = 0; i<n_elem; i++) x[i] = m[i];
}

// Assignment
template<class T>
SimpleArray<T> & SimpleArray<T>::operator=(const SimpleArray<T> &m){

  // Check if m refers to the same object
  if (this != &m) {

    // Check if the size matches
    if (n_elem != m.n_elem) {
      if (x != NULL) delete[] (x);
      n_elem = m.n_elem;
      x = n_elem > 0 ? new T[n_elem] : NULL;
    }

    // Copy the contents of m
    for (int i = 0; i<n_elem; i++) {
      x[i] = m[i];
    }
  }
  return *this;
}

// Get i'th element by x(i) : NO BOUND CHECKING
template<class T>
inline
T & SimpleArray<T>::operator()(int i){
  return x[i];
}

template<class T>
inline
T SimpleArray<T>::operator()(int i) const{
  return x[i];
}


// Get i'th element by x[i] : NO BOUND CHECKING
template<class T>
inline
T & SimpleArray<T>::operator[](int i){
  return x[i];
}

template<class T>
inline
T SimpleArray<T>::operator[](int i) const{
  return x[i];
}

// Resize
template<class T>
void SimpleArray<T>::resize(int n){

  // Check if resize is needed
  if (n != n_elem) {
    if (x != NULL) delete[] (x);
    n_elem = n;
    x = n > 0 ? new T[n] : NULL;
  }
}

template<class T>
SimpleArray<T>::~SimpleArray(){
  if (x != NULL) delete[] (x);
}

//==============================
// The matrix class
//==============================

template<class T>
class Matrix {
  public:
    T *M;
    int nx, ny;

    // Constructors
    Matrix();
    Matrix(int n1, int n2);				// Zero based array
    Matrix(const Matrix<T> &m);			// Copy constructor
    Matrix &operator=(const Matrix &m); // Assignment
    ~Matrix();	

    // Make T externally availible
    typedef T value_type;

    // Operators
    inline T &operator()(int i, int j);
    inline T  operator()(int i, int j) const;
    inline T &operator[](int i);
    inline T  operator[](int i) const;

    Matrix &operator=(T const m){set(m);}
    Matrix &operator+(T const m){add(m);}
    Matrix &operator+=(T const m){add(m);}
    Matrix &operator-(T const m){subtract(m);}
    Matrix &operator-=(T const m){subtract(m);}
    Matrix &operator*(T const m){multiply(m);}
    Matrix &operator*(Matrix<T> const m){multiply(m);}
    Matrix &operator*=(T const m){multiply(m);}
    Matrix &operator/(T const m){divide(m);}
    Matrix &operator/=(T const m){divide(m);}

    // Methods
    void set(T const m);
    void add(T const m);
    void subtract(T const m);
    void multiply(T const m);
    void multiply(Matrix<T> const &m);
    void divide(T const m);
    void resize(int newx, int newy);
};

//==============================
// Inline func. for matrix class
//==============================

// Constructors
template<class T>
Matrix<T>::Matrix() : nx(0), ny(0), M(NULL) {}

template<class T>
Matrix<T>::Matrix(int n1, int n2) : nx(n1), ny(n2),  M(n1*n2 > 0 ? new T[n1*n2] : NULL) {}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &m) : nx(m.nx), ny(m.ny), M(m.nx*m.ny > 0 ? new T[m.nx*m.ny] : NULL) {
  for (int i = 0; i<nx*ny; i++){
    M[i] = m[i];
  }
}

// Assignment
template<class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &m){

  // Check if m refers to the same object
  if (this != &m) {

    // Check if the size matches
    if (nx != m.nx || ny != m.ny) {
      if (M != NULL) delete[] (M);
      nx = m.nx;
      ny = m.ny;
      M = nx*ny > 0 ? new T[nx*ny] : NULL;
    }

    // Copy the contents of m
    for (int i = 0; i<nx*ny; i++) {
      M[i] = m[i]; 
    }
  }
  return *this;
}

// Resize (destoys containts)
template <class T>
void Matrix<T>::resize(int newx, int newy){
  if (newx != nx || newy != ny) {
    if (M != NULL) delete[] M;
    nx = newx;
    ny = newy;
    M  = nx*ny > 0 ? new T[nx*ny]: NULL;
  }
}

// Get i'th element (no bounds check)
template<class T>
inline
T & Matrix<T>::operator()(int i, int j){
  return M[i*ny + j];
}

template<class T>
inline
T Matrix<T>::operator()(int i, int j) const{
  return M[i*ny + j];
}

template<class T>
inline
T & Matrix<T>::operator[](int i){
  return M[i];
}

template<class T>
inline
T Matrix<T>::operator[](int i) const{
  return M[i];
}

// Assign with constant value
template<class T>
inline void Matrix<T>::set(T const m){
  for (int i = 0; i<nx; i++) {
    for (int j = 0; j<ny; j++) {
      M[i*ny + j] = m;
    }
  }
}

// Add constant matrix
template<class T>
inline void Matrix<T>::add(T const m){
  for (int i = 0; i<nx*ny; i++) {
    M[i] += m;
  }
}

// Subtract constant matrix
template<class T>
inline void Matrix<T>::subtract(T const m){
  for (int i = 0; i<nx*ny; i++) {
    M[i] -= m;
  }
}

// Multiply with constant
template<class T>
inline void Matrix<T>::multiply(T const m){
  for (int i = 0; i<nx*ny; i++) {
    M[i] *= m;
  }
}

// Divide with constant
template<class T>
inline void Matrix<T>::divide(T const m){
  for (int i = 0; i<nx*ny; i++) {
    M[i] /= m;
  }
}

// Matrix multiplication
template<class T>
inline void Matrix<T>::multiply(Matrix<T> const &m){
  Matrix<T> temp(nx,ny);	
  for (int i = 0; i<nx; i++) {
    for (int j = 0; j<ny; j++) {
      temp(i,j) = 0.0;
      for (int k = 0; k<ny; k++) {
        temp(i,j) += M[i*nx+k]*m[k*nx+j];
      }
    }
  }
  for (int i = 0; i<nx; i++) {
    M[i] = temp[i];
  }
}

// Destructor
template<class T>
Matrix<T>::~Matrix(){
  if (M != NULL) delete[] (M);
}

//==============================
// Typedefs
//==============================

typedef SimpleArray<double> VecDoub, VecDoub_O, VecDoub_IO;
typedef const SimpleArray<double> VecDoub_I;
typedef SimpleArray<int> VecInt;
typedef Matrix<double> MatDoub, MatDoub_O, MatDoub_IO;
typedef const Matrix<double> MatDoub_I;

//==============================
// General inline functions
//==============================

template<class T>
inline T SQR(const T a) {return a*a;}

  template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
{return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
{return b > a ? float(b) : (a);}

  template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
{return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
{return b < a ? float(b) : (a);}

  template<class T>
inline T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;}
