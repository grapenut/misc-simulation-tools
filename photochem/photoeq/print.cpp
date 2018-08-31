
#include <cmath>
#include "print.h"

void print_data(FILE *fout, data &d)
{
  fprintf(fout, "%15g ", d.flux);
  fprintf(fout, "%15g ", d.metal);
  fprintf(fout, "%15g ", d.temp);
  fprintf(fout, "%15g ", d.eden);
  fprintf(fout, "%15g ", d.h);
  fprintf(fout, "%15g ", d.hplus);
  fprintf(fout, "%15g ", d.he);
  fprintf(fout, "%15g ", d.heplus);
  fprintf(fout, "%15g ", d.heplusplus);
  fprintf(fout, "%15g ", d.c);
  fprintf(fout, "%15g ", d.cplus);
  fprintf(fout, "%15g ", d.cplusplus);
  fprintf(fout, "%15g ", d.thermal_time);
  fprintf(fout, "%15g ", d.recomb_time);
  fprintf(fout, "%15g\n", d.xmetal);
}

void debug_data(data &d)
{
  printf("%15g ", d.flux);
  printf("%15g ", d.metal);
  printf("%15g ", d.temp);
  printf("%15g ", d.eden);
  printf("%15g ", d.h);
  printf("%15g ", d.hplus);
  printf("%15g ", d.he);
  printf("%15g ", d.heplus);
  printf("%15g ", d.heplusplus);
  printf("%15g ", d.c);
  printf("%15g ", d.cplus);
  printf("%15g ", d.cplusplus);
  printf("%15g ", d.thermal_time);
  printf("%15g ", d.recomb_time);
  printf("%15g\n", d.xmetal);
}

