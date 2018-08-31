/*

  find the mass of a gas cloud around maximum gas density

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


#define PC_TO_CM	3.08567758e18
#define CM_TO_PC	3.2407793e-19
#define MSUN		2.0e33
#define mH		1.67e-24
#define PI		3.14159265359

typedef struct
{
  double x, y, z, dens, mass, vol, dist, temp, metal;
  
} cdata;

bool operator<(cdata const &a, cdata const &b)
{
  if (a.dist == b.dist) return (a.dens > b.dens);
  return (a.dist < b.dist);
}


void traceback(int sig)
{
  void *array[10];
  size_t size;
    
  // get void*'s for all entries on the stack
  size = backtrace(array, 10);
        
  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main(int argc, char * argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;

  double SAMPLE_RADIUS = 1000.0 * PC_TO_CM;
  double SAMPLE_RADIUS_SQ = SAMPLE_RADIUS*SAMPLE_RADIUS;
  double FOURTHIRDSPI = 4.0 * PI / 3.0;

  std::time_t start, stop;
  
  signal(SIGSEGV, traceback);

  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " filename" << endl << endl;
    return 1;
  }

  int argptr = 1;

  const std::string file_name = argv[argptr++];

  // Open the file

  QuickFlash::File::DataFile dfile(file_name);

  const QuickFlash::File::SimInfo & siminfo = dfile.get_sim_info();
  
  double scale = siminfo.get_real_scalar("scalefactor");
  double scale2 = scale*scale;
  double scale3 = scale*scale*scale;
  //double redshift = (1.0/scale) - 1.0;
  
  double sim_star_x, sim_star_y, sim_star_z;
  
  if (siminfo.get_real_scalar("sim_star_created") > 0.0)
  {
    sim_star_x = siminfo.get_real_scalar("sim_star_x");
    sim_star_y = siminfo.get_real_scalar("sim_star_y");
    sim_star_z = siminfo.get_real_scalar("sim_star_z");
  } else {
    sim_star_x = 0.0;
    sim_star_y = 0.0;
    sim_star_z = 0.0;
  }
  
  cerr << "Star particle: " << sim_star_x << " " << sim_star_y << " " << sim_star_z << endl;
  
  //////////////////////////////////////////////////////////////
  
  // load the grid cells into an array
  const unsigned int buffer_size = 512;
  unsigned int cache_size = 0;

  // Get mesh info
  const QuickFlash::File::MeshInfo & meshinfo = dfile.get_mesh_info();

  const unsigned int num_blocks = meshinfo.get_num_blocks();

  //cerr << "Num blocks [ " << num_blocks << " ]" << endl;

  // Open a dataset
  const QuickFlash::File::Dataset & dset_dens
    = dfile.get_dataset("dens", buffer_size);

  dset_dens.set_cache_size(cache_size);
  dset_dens.set_report_stats();

  const QuickFlash::File::Dataset & dset_temp
    = dfile.get_dataset("temp", buffer_size);

  dset_temp.set_cache_size(cache_size);
  dset_temp.set_report_stats();

  const QuickFlash::File::Dataset & dset_z
    = dfile.get_dataset("z", buffer_size);

  dset_z.set_cache_size(cache_size);
  dset_z.set_report_stats();

  //cerr << "Dataset loaded" << endl;

  // Set up temporary storage for cell data
  QuickFlash::Block::BlockData<double> data_dens, data_temp, data_z;
  std::vector<double> cell_center;
  std::vector<cdata> cell_data;
  
  unsigned int cell_count = 0;
  
  /////////////////////////////////////////////////////////////////////
  // capture cell data within SAMPLE_RADIUS
  
  // Run through the leaf nodes
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
  {
    if(meshinfo.is_leaf(block_index))
    {
      // Get the block
      const QuickFlash::Block::BlockInfo & block_info 
	= meshinfo.get_block_info(block_index);

      dset_dens.get_block_data(block_index, data_dens);
      dset_temp.get_block_data(block_index, data_temp);
      dset_z.get_block_data(block_index, data_z);

      const unsigned int num_cells = block_info.get_num_cells();
      const double cell_vol = block_info.get_cell_volume();

      for (unsigned int index = 0; index < num_cells; index++)
      {
	block_info.get_cell_center(index, cell_center);
	
	double dx = scale * (cell_center[0] - sim_star_x);
	double dy = scale * (cell_center[1] - sim_star_y);
	double dz = scale * (cell_center[2] - sim_star_z);
	double r2 = dx*dx + dy*dy + dz*dz;
	
	if (r2 > SAMPLE_RADIUS_SQ) continue;
        
        cell_count++;
        
        cdata c;
        
        c.x = cell_center[0];
        c.y = cell_center[1];
        c.z = cell_center[2];
        
        c.vol = cell_vol * scale3;
	c.dens = data_dens[index] / scale3;
	c.mass = data_dens[index] * cell_vol;
	c.temp = data_temp[index] * scale2;
	c.metal = data_z[index];
        c.dist = sqrt(r2);

	cell_data.push_back(c);
	
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  // radially bin the gas

  // lets just sort the arrays by cell_dist and read off the items between r and r+dr
  std::sort(cell_data.begin(), cell_data.end());
  
  // bin mass radially between bin_min and bin_max
  double bin_min = 3.0e17;
  double bin_max = 3.0e21;
  
  int nbins = 100;
  
  std::vector<double> mass, metal, volume, temperature;
  mass.resize(nbins);
  metal.resize(nbins);
  volume.resize(nbins);
  temperature.resize(nbins);
  
  double log_max_over_min = std::log10(bin_max / bin_min);
  
  for (int i = 0; i < cell_count; i++)
  {
    int bin;

    if (cell_data[i].dist <= bin_min)
    {
      bin = 0;
    } else {
      bin = std::floor(double(nbins-2) * std::log10(cell_data[i].dist / bin_min) / log_max_over_min) + 1.0;
    }
    
    mass[bin] += cell_data[i].mass;
    metal[bin] += cell_data[i].metal * cell_data[i].mass;
    volume[bin] += cell_data[i].vol;
    temperature[bin] += cell_data[i].mass * cell_data[i].temp;
  }
  
  
  double mass_sum = 0.0;

  // output binned data
  for (int i = 0; i < nbins; i++)
  {
    double inner, outer, vol;

    if (i == 0)
    {
      inner = 0.0;
      outer = bin_min;
    } else {
      inner = bin_min * std::pow((bin_max/bin_min), double(i)/double(nbins-1));
      outer = bin_min * std::pow((bin_max/bin_min), double(i+1)/double(nbins-1));
    }
    
    //vol = 4.0*M_PI/3.0 * (outer*outer*outer - inner*inner*inner);
    
    if (volume[i] > 0.0)
    {
      mass_sum += mass[i];
      cout << i << " " << inner << " " << outer << " " << mass_sum << " " << mass[i] << " " << temperature[i]/mass[i] << " " << metal[i]/mass[i] << " " << volume[i] << " " << mass[i] / volume[i] << endl;
    }
  }
  
/*
  ofstream fout("bins.txt");
  fout << setprecision(12);

  for (int i = 0; i < cell_count; i++)
  {
    fout << setw(20) << cell_data[i].vol 
         << setw(20) << cell_data[i].dens 
         << setw(20) << cell_data[i].temp 
         << setw(20) << cell_data[i].ion 
         << endl;
  }
  
  fout.close();
*/
  
}
