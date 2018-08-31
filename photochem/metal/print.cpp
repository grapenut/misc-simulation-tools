
#include <cmath>
#include "print.h"

void print_data(FILE *fout, data &d)
{
  fprintf(fout, "%15g ", d.dens);
  fprintf(fout, "%15g ", d.metal);
  fprintf(fout, "%15g ", d.temp);
  fprintf(fout, "%15g ", d.eden);
  fprintf(fout, "%15g ", d.cooling);
  fprintf(fout, "%15g ", d.heating);
  fprintf(fout, "%15g ", d.h);
  fprintf(fout, "%15g ", d.hplus);
  fprintf(fout, "%15g ", d.he);
  fprintf(fout, "%15g ", d.heplus);
  fprintf(fout, "%15g ", d.heplusplus);
  fprintf(fout, "%15g ", d.c);
  fprintf(fout, "%15g ", d.cplus);
  fprintf(fout, "%15g\n", d.cplusplus);
}

void debug_data(data &d)
{
}

