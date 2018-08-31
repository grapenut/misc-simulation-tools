#ifndef __DATA_H
#define __DATA_H

struct data
{
  // phase space coordinates
  double dens, temp, metal, flux;
  
  // ions
  double eden, h, hplus, he, heplus, heplusplus, c, cplus, cplusplus;
  
  // rates
  double heating, cooling;
  
  // timescales
  double age, thermal_time, recomb_rate;
};

#endif

