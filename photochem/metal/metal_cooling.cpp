
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
#include "cloudy.h"
#include "print.h"

// number of data points per decade in log-space
#define DENS_PER_DECADE		3
#define TEMP_PER_DECADE		5
#define METAL_PER_DECADE	4

// log10 minima/maxima for parameter space
#define TEMP_MIN	1.0
#define TEMP_MAX	9.0

#define DENS_MIN	-4.0
#define DENS_MAX	8.0

#define METAL_MIN	-6.0
#define METAL_MAX	-0.75

using std::vector;

int main(int argc, char **argv)
{
  vector<double> metal_list, dens_list, temp_list;
  int num_metal, num_dens, num_temp;
  
  int d, t, z;
  
  const int len = 256;
  char filename[len];
  FILE *fout;
  
  data result;

  double helium_frac = 0.08;
  double carbon_frac = 2.45e-4;
  double oxygen_frac = 4.90e-4;
  double metal_scale = (1.0+4.0*helium_frac) / (14.0*carbon_frac + 16.0*oxygen_frac);
  
  // disable output buffering
  setbuf(stdout, NULL);

  metal_list = logspace(METAL_MIN, METAL_MAX, METAL_PER_DECADE);
  temp_list = logspace(TEMP_MIN, TEMP_MAX, TEMP_PER_DECADE);
  dens_list = logspace(DENS_MIN, DENS_MAX, DENS_PER_DECADE);
  
  num_dens = dens_list.size();
  num_temp = temp_list.size();
  num_metal = metal_list.size();
  
  //printf("n = %d \t z = %d \t t = %d \t f = %d\n", num_dens, num_metal, num_temp, num_flux);
  //printf("Done with init. num_total = %d\n", num_dens*num_metal*num_temp*num_flux);

  snprintf(filename, len, "output.metal_cooling.txt");
  if ((fout = fopen(filename, "w")) == NULL)
  {
    printf("ERROR: Could not open %s for writing.\n", filename);
    cdEXIT(EXIT_FAILURE);
  }

  result.age = 1.0e6;
  
  for (z = 0; z < num_metal; z++)
  {
    result.metal = metal_list[z];
    result.xmetal = result.metal / (1.0 - result.metal) * metal_scale;

    for (d = 0; d < num_dens; d++)
    {
      result.dens = dens_list[d];
      
      for (t = 0; t < num_temp; t++)
      {
        result.temp = temp_list[t];
        
        printf("PROC: %10g %10g %10g %10g\n", result.dens, result.temp, result.metal, result.xmetal);
        run_cloudy(result);

        //debug_data(result);
        print_data(fout, result);
      }
    }
  }

  fclose(fout);
  
  cdEXIT(EXIT_SUCCESS);
}

