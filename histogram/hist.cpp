// C++ program file parts.cpp

/*
  By Nathan C. Hearn
     December 11, 2006

  Demonstration of particle-related functions for the QuickFlash library.




  loop through the particles and create a Point object for each one
  
  put the Point objects in a vector
  
  call KDtree(vector<Points>)
  
  


*/

#include <iostream>
#include <hdf5.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <map>
#include <functional>
#include <numeric>
#include <sstream>
#include <iterator>
#include <signal.h>
#include <execinfo.h>


#include "quickflash_file_siminfo.hpp"
#include "quickflash_file_datafile.hpp"
#include "quickflash_block_blockdata.hpp"
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_file_meshinfo.hpp"
#include "quickflash_particles_data.hpp"

#include "nr3.h"
#include "pointbox.h"
#include "kdtree.h"

#define PC_TO_CM	3.08567758e18
#define CM_TO_PC	3.2407793e-19
#define MSUN		2.0e33
#define mH		1.67e-24
#define PI		3.14159265359

static double dens_min = 1.0e-3;
static double dens_max = 1.0e5;
static double log_range = 0.0;
static size_t num_bins = 17;
static size_t num_sn = 7;

size_t binIndex(double val)
{
  if (val < dens_min) 
    return 0; 
  else if (val >= dens_max) 
    return num_bins - 1; 
  else
    return size_t(std::floor(double(num_bins - 1) * std::log(val / dens_min) / log_range) + 1.0);
}

double binCenter(size_t i)
{
  if(i == 0) 
    return 0.5 * dens_min;
  else
    return dens_min * std::exp((double(i) - 0.5) * log_range / double(num_bins - 1));
}

double binLeft(size_t i)
{
  if(i == 0) 
    return 0.0;
  else
    return dens_min * std::exp((double(i) - 1.0) * log_range / double(num_bins - 1));
}

double binRight(size_t i)
{
  return dens_min * std::exp((double(i)) * log_range / double(num_bins - 1));
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

int main(int argc, char * argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;

  double SAMPLE_RADIUS_SQ;
  string inputfile;

  double part_mass[num_sn][num_bins], fluid_metal_mass[num_bins], fluid_total_mass[num_bins], part_total_mass[num_bins];
  unsigned int part_count[num_sn][num_bins], fluid_count[num_bins];
  
  vector<double> cell_x, cell_y, cell_z, cell_dens, cell_vol, cell_mass, cell_metal, cell_metal_mass, cell_parts[num_sn];
  
  // read data from argv[1]
  if (argc > 2)
  {
    double r = atof(argv[1]);
    cout << "Using radius " << r << " parsecs" << endl;
    SAMPLE_RADIUS_SQ = std::pow(r * PC_TO_CM, 2.0);
    inputfile = argv[2];
  }
  else if (argc > 1)
  {
    SAMPLE_RADIUS_SQ = 1e99;
    inputfile = argv[1];
  }
  else
  {
    cerr << "Usage: " << argv[0] << " [<distance>] <file>" << endl;
    return 0;
  }
  ifstream din(inputfile.c_str());
  
  cout << "Reading data from: " << inputfile << endl;
  
  while (din.good())
  {
    double cx, cy, cz, cdens, cvol, cmetal, cparts[num_sn], r2;
  
    din >> cx >> cy >> cz >> cdens >> cvol >> cmetal;
    
    //if (cdens > dens_max) dens_max = cdens;
    //if (cdens < dens_min) dens_min = cdens;

    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      din >> cparts[sn];
    }
    
    r2 = cx*cx + cy*cy + cz*cz;
    if (r2 > SAMPLE_RADIUS_SQ) continue;

    cell_x.push_back(cx);
    cell_y.push_back(cy);
    cell_z.push_back(cz);
    cell_dens.push_back(cdens/mH);
    cell_vol.push_back(cvol);
    cell_mass.push_back(cdens*cvol);
    cell_metal.push_back(cmetal);
    cell_metal_mass.push_back(cdens*cvol*cmetal);

    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      cell_parts[sn].push_back(cparts[sn]);
    }
    
  }
  
  //dens_min *= 0.999;
  //dens_max *= 1.001;
  
  log_range = std::log(dens_max / dens_min);
  
  unsigned int num_cells = cell_dens.size();
  
  cout << "Done reading file: " << num_cells << " cells" << endl;
  
  // bin data
  // add up particle mass for each supernova in bins and keep count of the cells that contribute
  // add up fluid mass and keep count of cells
  
  // init bins/counts
  for (unsigned int bin = 0; bin < num_bins; bin++)
  {
    fluid_metal_mass[bin] = 0.0;
    fluid_total_mass[bin] = 0.0;
    fluid_count[bin] = 0;
    
    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      part_mass[sn][bin] = 0.0;
      part_count[sn][bin] = 0;
    }
  }

  // bin data
  for (unsigned int cell = 0; cell < num_cells; cell++)
  {
    unsigned int bin = binIndex(cell_dens[cell]);
    
    fluid_metal_mass[bin] += cell_metal_mass[cell];
    fluid_total_mass[bin] += cell_mass[cell];
    fluid_count[bin]++;
  
    double total_mass = 0.0;
    
    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      total_mass += cell_parts[sn][cell];
      part_mass[sn][bin] += cell_parts[sn][cell];
      part_count[sn][bin]++;
    }
    
    part_total_mass[bin] += total_mass;
  }

  double average_fluid_metallicity[num_bins];
  double variance_fluid_metallicity[num_bins];

  double average_part_metallicity[num_sn][num_bins];
  
  for (unsigned int bin = 0; bin < num_bins; bin++)
  {
    average_fluid_metallicity[bin] = 0.0;
    variance_fluid_metallicity[bin] = 0.0;
    
    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      average_part_metallicity[sn][bin] = 0.0;
    }
  }
  
  // compute averages for fluid and particles
  for (unsigned int bin = 0; bin < num_bins; bin++)
  {
    if (fluid_total_mass[bin] > 0.0)
      average_fluid_metallicity[bin] = fluid_metal_mass[bin] / fluid_total_mass[bin];
    else
      average_fluid_metallicity[bin] = 0.0;
    
    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      if (fluid_total_mass[bin] > 0.0)
        average_part_metallicity[sn][bin] = part_mass[sn][bin] / fluid_total_mass[bin]; // * average_fluid_metallicity[bin];
      else
        average_part_metallicity[sn][bin] = 0.0;
    }
  }

  // compute variances for fluid bins
  for (unsigned int cell = 0; cell < num_cells; cell++)
  {
    unsigned int bin = binIndex(cell_dens[cell]);
    
    if (fluid_total_mass[bin] > 0.0)
      variance_fluid_metallicity[bin] += cell_mass[cell] * std::pow(cell_metal[cell] - average_fluid_metallicity[bin], 2.0) / fluid_total_mass[bin];
  }
  
  // output binned data, averages, and variances
  stringstream fname;
  fname << "hist_" << argv[1] << ends;
  
  ofstream dout(fname.str().c_str());
  dout << '#';
  dout << setw(19) << "bin";
  dout << setw(20) << "binCenter";
  dout << setw(20) << "fluid_total_mass";
  dout << setw(20) << "fluid_metal_mass";
  dout << setw(20) << "avg_fluid_metal";
  dout << setw(20) << "var_fluid_metal";

  dout << setw(20) << "part_total_mass";
  for (unsigned int sn = 0; sn < num_sn; sn++)
  {
    stringstream mass_field, avg_field;
    mass_field << "part_mass" << sn;
    avg_field << "avg_part_metal" << sn;
    
    dout << setw(20) << mass_field.str();
    dout << setw(20) << avg_field.str();
  }
  
  dout << endl;

  for (unsigned int bin = 0; bin < num_bins; bin++)
  {
    dout << setw(20) << bin;
    dout << setw(20) << binCenter(bin);
    dout << setw(20) << fluid_total_mass[bin];
    dout << setw(20) << fluid_metal_mass[bin];
    dout << setw(20) << average_fluid_metallicity[bin];
    dout << setw(20) << variance_fluid_metallicity[bin];
    dout << setw(20) << part_total_mass[bin];

    for (unsigned int sn = 0; sn < num_sn; sn++)
    {
      dout << setw(20) << part_mass[sn][bin];
      dout << setw(20) << average_part_metallicity[sn][bin];
    }
    
    dout << endl;
  }
  dout.close();

  
  return 0;
  


  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////


  
}


