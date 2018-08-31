

#include "logspace.h"
#include "data.h"
#include <cmath>

vector<double> logspace(double min, double max, double steps)
{

  double d, f;
  vector<double> space;
  double step = 1.0/steps;
  
  for (d = min; d <= 1.01*max; d += step)
  {
    f = pow(10.0, d);
    space.push_back(f);
  }

  return space;
}








