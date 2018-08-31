
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
#define DENS_PER_DECADE		1
#define TEMP_PER_DECADE		10
#define METAL_PER_DECADE	1
#define FLUX_PER_DECADE		1

// log10 minima/maxima for parameter space
#define TEMP_MIN	0.0
#define TEMP_MAX	8.0

#define DENS_MIN	-4.0
#define DENS_MAX	8.0

#define METAL_MIN	-3.0
#define METAL_MAX	3.0

#define FLUX_MIN	-10.0
#define FLUX_MAX	20.0

using std::vector;

int main(int argc, char **argv)
{
  vector<double> metal_list, temp_list;
  int num_metal, num_temp;
  
  int t, z;
  
  char filename[256];
  FILE *fout;
  
  data result;
  
  // disable output buffering
  setbuf(stdout, NULL);

  //metal_list = logspace(METAL_MIN, METAL_MAX, METAL_PER_DECADE);
  temp_list = logspace(TEMP_MIN, TEMP_MAX, TEMP_PER_DECADE);
  
  num_temp = temp_list.size();
  //num_metal = metal_list.size();
  
  //printf("n = %d \t z = %d \t t = %d \t f = %d\n", num_dens, num_metal, num_temp, num_flux);
  //printf("Done with init. num_total = %d\n", num_dens*num_metal*num_temp*num_flux);

  if (argc > 1)
  {
    strcpy(filename, argv[1]);
  } else {
    sprintf(filename, "output.recomb.dat");
  }
  
  if ((fout = fopen(filename, "w")) == NULL)
  {
    printf("ERROR: Could not open %s for writing.\n", filename);
    cdEXIT(EXIT_FAILURE);
  }

  result.age = 1.0e6;
  result.dens = 1.0;
  result.metal = 1.0;
  
  //for (z = 0; z < num_metal; z++)
  //{
  //  result.metal = metal_list[z];

    for (t = 0; t < num_temp; t++)
    {
      result.temp = temp_list[t];
      
      printf("RUN: %10g\n", result.temp);
      run_cloudy(result);
      debug_data(result);
      print_data(fout, result);
    }
  //}
  
  fclose(fout);
  
  //cdEXIT(EXIT_SUCCESS);
  
  return 0;
}

