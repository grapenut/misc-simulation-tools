

#include "logspace.h"
#include "data.h"
#include <cmath>

vector<double> logspace(double min, double max, double steps)
{

  double d;
  vector<double> space;
  double step = 1.0/steps;
  
  for (d = min; d <= max; d += step)
  {
    space.push_back(pow(10.0, d));
  }

  return space;
}








