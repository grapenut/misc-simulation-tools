#ifndef __DATA_H
#define __DATA_H

// number of data points per decade in log-space
#define DENS_PER_DECADE		1
#define TEMP_PER_DECADE		1
#define METAL_PER_DECADE	1
#define FLUX_PER_DECADE		1

// log10 minima/maxima for parameter space
#define TEMP_MIN	1.0
#define TEMP_MAX	8.0

#define DENS_MIN	-4.0
#define DENS_MAX	8.0

#define METAL_MIN	-6.0
#define METAL_MAX	3.0

#define FLUX_MIN	-15.0
#define FLUX_MAX	10.0



struct data
{
  // phase space coordinates
  double dens, temp, metal, flux;
  
  // rates
  double heating, cooling;
  
  // ions
  double eden, h, hplus, he, heplus, heplusplus, c, cplus;
};

#endif

