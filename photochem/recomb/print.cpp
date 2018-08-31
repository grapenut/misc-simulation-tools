
#include <cmath>
#include "print.h"

using std::abs;

void print_data(FILE *fout, data &d)
{
  fprintf(fout, "%15g ", d.temp);
  fprintf(fout, "%15g ", d.metal);
  fprintf(fout, "%15g ", d.thermal_time);
  fprintf(fout, "%15g\n", d.recomb_rate);
}

void debug_data(data &d)
{
  printf("%15g ", d.temp);
  printf("%15g ", d.metal);
  printf("%15g ", d.thermal_time);
  printf("%15g\n", d.recomb_rate);
}

