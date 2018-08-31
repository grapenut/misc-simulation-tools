
// generate data tables for the PhotoChem variant of FLASH chemistry

#include "cddefines.h"
#include "cddrive.h"

#include <iostream>
#include <vector>

//#include "cpu.h"
//#include "hmi.h"
//#include "mole.h"
#include "dense.h"

#include "data.h"
#include "logspace.h"

using std::vector;

int main(int argc, char **argv)
{
  vector<double> temp_list, metal_list, dens_list, flux_list;
  int num_temp, num_metal, num_dens, num_flux;
  
  double dens, temp, flux, metal;
  int n, t, f, z;
  
  char filename[256];
  FILE *fout;
  
  data result;

  temp_list = logspace(TEMP_MIN, TEMP_MAX, TEMP_PER_DECADE);
  metal_list = logspace(METAL_MIN, METAL_MAX, METAL_PER_DECADE);
  flux_list = logspace(FLUX_MIN, FLUX_MAX, FLUX_PER_DECADE);
  dens_list = logspace(DENS_MIN, DENS_MAX, DENS_PER_DECADE);
  
  num_temp = temp_list.size();
  num_flux = flux_list.size();
  num_metal = metal_list.size();
  num_dens = dens_list.size();
  
  //printf("n = %d \t z = %d \t t = %d \t f = %d\n", num_dens, num_metal, num_temp, num_flux);
  //printf("Done with init. num_total = %d\n", num_dens*num_metal*num_temp*num_flux);

  sprintf(filename, "output.photoeq.dat");
  if ((fout = fopen(filename, "w")) == NULL)
  {
    printf("ERROR: Could not open %s for writing.\n", filename);
    cdEXIT(EXIT_FAILURE);
  }

  result.age = 1.0e6;
  
  for (z = 0; z < num_metal; z++)
  {
    result.metal = metal_list[z];

    for (n = 0; n < num_dens; n++)
    {
      result.dens = dens_list[n];
      
      for (t = 0; t < num_temp; t++)
      {
        result.temp = temp_list[t];
        
        for (f = 0; f < num_flux; f++)
        {
          result.flux = flux_list[f];
          
          run_cloudy(&result);
          print_data(&result);
        }
        
        fclose(fout);
      }
    }
  }

  printf("Done.\n");

}

