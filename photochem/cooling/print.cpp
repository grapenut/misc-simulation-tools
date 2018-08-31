
#include <iostream>

#include "data.h"

void print_data(data& d)
{
  fprintf(fout, "%15g ", d.dens);
  fprintf(fout, "%15g ", d.temp);
  fprintf(fout, "%15g ", d.metal);
  fprintf(fout, "%15g ", d.flux);
  fprintf(fout, "%15g ", d.eden);
  fprintf(fout, "%15g ", d.cooling);
  fprintf(fout, "%15g ", d.heating);
  fprintf(fout, "%15g ", d.);
  fprintf(fout, "%15g ", d.);
  fprintf(fout, "%15g ", d.);
  fprintf(fout, "%15g ", d.);
  fprintf(fout, "%15g ", d.);
  fprintf(fout, "%15g ", d.);
  %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n", dens, metal, temp, flux, eden, cooling, heating, h, hplus, he, heplus, heplusplus);
}








