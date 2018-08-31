#ifndef _FUNCTION_HH
#define _FUNCTION_HH

#include <vector>
#include <algorithm>
#include <gsl/gsl_spline.h>
#include <functional>
#include <cmath>

using namespace std;

class Spline
{
  size_t _size;
  double * _x, * _y;
  double _small_y, _large_y;
  gsl_interp_accel * _acc;
  gsl_spline * _spline;
public:
  Spline(const size_t & count, double x[], double y[], double small_y, double large_y) :
    _size(count),
    _x(new double[_size]),
    _y(new double[_size]),
    _small_y(small_y),
    _large_y(large_y)
  {
    copy(x, x + _size, _x);
    copy(y, y + _size, _y);
    _acc = gsl_interp_accel_alloc();
    _spline = gsl_spline_alloc(gsl_interp_linear, _size);
    gsl_spline_init(_spline, _x, _y, _size);
  }
  Spline(const vector<double> & x, const vector<double> & y, double small_y, double large_y) :
    _size(x.size()),
    _x(new double[_size]),
    _y(new double[_size]),
    _small_y(small_y),
    _large_y(large_y)
  {
    copy(x.begin(), x.end(), _x);
    copy(y.begin(), y.end(), _y);
    _acc = gsl_interp_accel_alloc();
    _spline = gsl_spline_alloc(gsl_interp_linear, _size);
    gsl_spline_init(_spline, _x, _y, _size);
  }
  Spline(const Spline & function) :
    _size(function._size),
    _x(new double[_size]),
    _y(new double[_size]),
    _small_y(function._small_y),
    _large_y(function._large_y)
  {
    copy(function._x, function._x + _size, _x);
    copy(function._y, function._y + _size, _y);
    _acc = gsl_interp_accel_alloc();
    _spline = gsl_spline_alloc(gsl_interp_linear, _size);
    gsl_spline_init(_spline, _x, _y, _size);
  }
  Spline & operator= (const Spline & function)
  {
    if(&function != this)
      {
	dealloc();
	_size = function._size;
	_x = new double[_size];
	_y = new double[_size];
	_small_y = function._small_y;
	_large_y = function._large_y;
	copy(function._x, function._x + _size, _x);
	copy(function._y, function._y + _size, _y);
	_acc = gsl_interp_accel_alloc();
	_spline = gsl_spline_alloc(gsl_interp_linear, _size);
	gsl_spline_init(_spline, _x, _y, _size);
      }
    return * this;
  }
  ~Spline() 
  {
    dealloc();
  }
  double operator() (const double & x) const 
  { 
    return eval(x);
  }
  double eval(const double & x) const 
  { 
    if(x < _x[0]) return _small_y;
    else if(x > _x[_size - 1]) return _large_y;
    else return gsl_spline_eval(_spline, x, _acc);
  }
  double deriv(const double & x) const
  {
    return gsl_spline_eval_deriv(_spline, x, _acc);
  }
  double deriv2(const double & x) const
  {
    return gsl_spline_eval_deriv2(_spline, x, _acc);
  }
  double integ(const double & a, const double & b) const
  {
    return gsl_spline_eval_integ(_spline, a, b, _acc);
  }
  double minX() const { return _x[0]; }
  double maxX() const { return _x[_size - 1]; }
  size_t size() const { return _size; }

//   Spline inverse() const 
//   {
//     return Spline(_size, _y, _x);    
//   }

  vector<double> getGrid() const
  {
    vector<double> grid(_size);
    copy(_x, _x + _size, grid.begin());
    return grid;
  }
private:
  void dealloc() 
  {  
    gsl_spline_free(_spline);
    gsl_interp_accel_free(_acc);
    delete [] _y;
    delete [] _x; 
  }
};

// class LinearFunction
// {
//   size_t _size;
//   double * _x, * _y;
//   gsl_interp_accel * _acc;
//   gsl_spline * _spline;
// public:
//   LinearFunction(const size_t & count, double x[], double y[]) :
//     _size(count),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     copy(x, x + _size, _x);
//     copy(y, y + _size, _y);
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_linear, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   LinearFunction(const vector<double> & x, const vector<double> & y) :
//     _size(x.size()),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     copy(x.begin(), x.end(), _x);
//     copy(y.begin(), y.end(), _y);
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_linear, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   LinearFunction(const LinearFunction & function) :
//     _size(function._size),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     copy(function._x, function._x + _size, _x);
//     copy(function._y, function._y + _size, _y);
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_linear, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   LinearFunction & operator= (const LinearFunction & function)
//   {
//     if(&function != this)
//       {
// 	dealloc();
// 	_size = function._size;
// 	_x = new double[_size];
// 	_y = new double[_size];
// 	copy(function._x, function._x + _size, _x);
// 	copy(function._y, function._y + _size, _y);
// 	_acc = gsl_interp_accel_alloc();
// 	_spline = gsl_spline_alloc(gsl_interp_linear, _size);
// 	gsl_spline_init(_spline, _x, _y, _size);
//       }
//     return * this;
//   }
//   ~LinearFunction() 
//   {
//     dealloc();
//   }
//   double operator() (const double & x) const 
//   { 
//     return gsl_spline_eval(_spline, x, _acc);
//   }
//   double eval(const double & x) const 
//   { 
//     return gsl_spline_eval(_spline, x, _acc);
//   }
//   double deriv(const double & x) const
//   {
//     return gsl_spline_eval_deriv(_spline, x, _acc);
//   }
//   double deriv2(const double & x) const
//   {
//     return gsl_spline_eval_deriv2(_spline, x, _acc);
//   } 
//   double integ(const double & a, const double & b) const
//   {
//     return gsl_spline_eval_integ(_spline, a, b, _acc);
//   }
//   double minX() const { return _x[0]; }
//   double maxX() const { return _x[_size - 1]; }
//   size_t size() const { return _size; }

//   LinearFunction inverse() const 
//   {
//     return LinearFunction(_size, _y, _x);    
//   }
// private:
//   void dealloc() 
//   {  
//     gsl_spline_free(_spline);
//     gsl_interp_accel_free(_acc);
//     delete [] _y;
//     delete [] _x; 
//   }
// };

// class PowerLaw
// {
//   size_t _size;
//   double * _x, * _y;
//   gsl_interp_accel * _acc;
//   gsl_spline * _spline;
// public:
//   PowerLaw(const size_t & count, double x[], double y[]) :
//     _size(count),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     transform(x, x + _size, _x, ptr_fun<double, double> (&log));
//     transform(y, y + _size, _y, ptr_fun<double, double> (&log));
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_cspline, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   PowerLaw(const vector<double> & x, const vector<double> & y) :
//     _size(x.size()),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     transform(x.begin(), x.end(), _x, ptr_fun<double, double> (&log));
//     transform(y.begin(), y.end(), _y, ptr_fun<double, double> (&log));
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_cspline, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   PowerLaw(const PowerLaw & function) :
//     _size(function._size),
//     _x(new double[_size]),
//     _y(new double[_size])
//   {
//     copy(function._x, function._x + _size, _x);
//     copy(function._y, function._y + _size, _y);
//     _acc = gsl_interp_accel_alloc();
//     _spline = gsl_spline_alloc(gsl_interp_cspline, _size);
//     gsl_spline_init(_spline, _x, _y, _size);
//   }
//   PowerLaw & operator= (const PowerLaw & function)
//   {
//     if(&function != this)
//       {
// 	dealloc();
// 	_size = function._size;
// 	_x = new double[_size];
// 	_y = new double[_size];
// 	copy(function._x, function._x + _size, _x);
// 	copy(function._y, function._y + _size, _y);
// 	_acc = gsl_interp_accel_alloc();
// 	_spline = gsl_spline_alloc(gsl_interp_cspline, _size);
// 	gsl_spline_init(_spline, _x, _y, _size);
//       }
//     return * this;
//   }
//   ~PowerLaw() 
//   {
//     dealloc();
//   }
//   double operator() (const double & x) const 
//   { 
//     return exp(gsl_spline_eval(_spline, log(x), _acc));
//   }
//   double eval(const double & x) const 
//   { 
//     return exp(gsl_spline_eval(_spline, log(x), _acc));
//   }
// //   double deriv(const double & x) const
// //   {
// //     return gsl_spline_eval_deriv(_spline, x, _acc);
// //   }
// //   double deriv2(const double & x) const
// //   {
// //     return gsl_spline_eval_deriv2(_spline, x, _acc);
// //   }
// //   double integ(const double & a, const double & b) const
// //   {
// //     return gsl_spline_eval_integ(_spline, a, b, _acc);
// //   }
//   double minX() const { return exp(_x[0]); }
//   double maxX() const { return exp(_x[_size - 1]); }
//   size_t size() const { return _size; }

//   PowerLaw inverse() const 
//   {
//     vector<double> temp_x(_size), temp_y(_size);

//     for(size_t i = 0; i < _size; i ++)
//       {
// 	temp_x[_size - i - 1] = exp(_x[i]);
// 	temp_y[_size - i - 1] = exp(_y[i]);
//       }
//     return PowerLaw(temp_y, temp_x);    
//   }

//   vector<double> getGrid() const
//   {
//     vector<double> grid(_size);
//     transform(_x, _x + _size, grid.begin(), ptr_fun<double, double> (&exp));
//     return grid;
//   }
// private:
//   void dealloc() 
//   {  
//     gsl_spline_free(_spline);
//     gsl_interp_accel_free(_acc);
//     delete [] _y;
//     delete [] _x; 
//   }
// };

#endif






