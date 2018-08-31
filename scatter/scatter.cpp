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
#define METAL_FLOOR	1.0e-10

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
  
  int checkpoint = siminfo.get_integer_scalar("checkpointfilenumber");
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

  const QuickFlash::File::Dataset & dset_gpot
    = dfile.get_dataset("gpot", buffer_size);

  dset_gpot.set_cache_size(cache_size);
  dset_gpot.set_report_stats();

  //cerr << "Dataset loaded" << endl;

  // Set up temporary storage for cell data
  QuickFlash::Block::BlockData<double> data_dens, data_temp, data_z, data_gpot;
  std::vector<double> cell_center;
  std::vector<cdata> cell_data;
  
  unsigned int cell_count = 0;
  
  /////////////////////////////////////////////////////////////////////
  // capture cell data within SAMPLE_RADIUS
  
  double max_dens = 0.0;
  double min_dens = 1.0e99;
  double max_metal = 0.0;
  double min_metal = 1.0e99;
  double max_dist = 0.0;
  double min_dist = 1.0e99;

  double max_gpot = 0.0;
  double halo_x = 0.0;
  double halo_y = 0.0;
  double halo_z = 0.0;
  
  double min_size = 1.0e99;
  
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
      dset_gpot.get_block_data(block_index, data_gpot);

      const unsigned int num_cells = block_info.get_num_cells();
      const double cell_vol = block_info.get_cell_volume();
      
      double len = scale*cbrt(cell_vol);
      if (len < min_size) min_size = len;

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
        
        if (c.metal < METAL_FLOOR) c.metal = METAL_FLOOR;
        
        if (c.metal > max_metal) max_metal = c.metal;
        if (c.metal < min_metal) min_metal = c.metal;
        if (c.dens > max_dens) max_dens = c.dens;
        if (c.dens < min_dens) min_dens = c.dens;
        
        if (data_gpot[index] < max_gpot)
        {
          max_gpot = data_gpot[index];
          halo_x = c.x;
          halo_y = c.y;
          halo_z = c.z;
        }

	cell_data.push_back(c);
      }
    }
  }
  
  for (int i = 0; i < cell_count; i++)
  {
    double dx = scale * (cell_data[i].x - halo_x);
    double dy = scale * (cell_data[i].y - halo_y);
    double dz = scale * (cell_data[i].z - halo_z);
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (r < min_size) r = min_size;
    cell_data[i].dist = r;

    if (cell_data[i].dist > max_dist) max_dist = cell_data[i].dist;
    if (cell_data[i].dist < min_dist) min_dist = cell_data[i].dist;
  }
  
  cerr << "Metal limits: " << min_metal << " " << max_metal << endl;
  cerr << "Density limits: " << min_dens << " " << max_dens << endl;
  cerr << "Distance limits: " << min_dist << " " << max_dist << endl;
  cerr << "Highest dark matter density: " << halo_x << " " << halo_y << " " << halo_z << endl;
  
  double off_x = halo_x - sim_star_x;
  double off_y = halo_y - sim_star_y;
  double off_z = halo_z - sim_star_z;
  double offset = sqrt(off_x*off_x + off_y*off_y + off_z*off_z);
  
  cerr << "Distance: " << CM_TO_PC*offset*scale << endl;
  
  /////////////////////////////////////////////////////////////////////
  // radially bin the gas

  // lets just sort the arrays by cell_dist and read off the items between r and r+dr
  std::sort(cell_data.begin(), cell_data.end());
  
  // bin mass radially between bin_min and bin_max
  int nbins = 100;
  
  double mass_sum[nbins][nbins];
  double dist_sum[nbins][nbins];

  double log_metal_max_over_min = std::log10(max_metal / min_metal);
  double log_dens_max_over_min = std::log10(max_dens / min_dens);
  double log_dist_max_over_min = std::log10(max_dist / min_dist);
  
  for (int i = 0; i < nbins; i++) 
  {
    for (int j = 0; j < nbins; j++)
    {
      mass_sum[i][j] = 1.0e-99;
      dist_sum[i][j] = 1.0e-99;
    }
  }
  
  for (int i = 0; i < cell_count; i++)
  {
    int mbin, nbin, rbin;

    mbin = std::floor(double(nbins-1) * std::log10(cell_data[i].metal / min_metal) / log_metal_max_over_min);
    nbin = std::floor(double(nbins-1) * std::log10(cell_data[i].dens / min_dens) / log_dens_max_over_min);
    rbin = std::floor(double(nbins-1) * std::log10(cell_data[i].dist / min_dist) / log_dist_max_over_min);
    
    mass_sum[mbin][nbin] += cell_data[i].mass;
    dist_sum[rbin][nbin] += cell_data[i].mass;
  }
  

  // output binned data
  double mdist, ndist, rdist;
  
  stringstream metalfile;
  metalfile << "metal_" << setw(4) << setfill('0') << checkpoint << ".dat";

  ofstream dout(metalfile.str().c_str());
  for (int i = 0; i < nbins; i++)
  {
    mdist = min_metal * std::pow((max_metal/min_metal), (double(i)+0.5)/double(nbins-1));
    
    for (int j = 0; j < nbins; j++)
    {
      ndist = min_dens * std::pow((max_dens/min_dens), (double(j)+0.5)/double(nbins-1));
      
      if (mass_sum[i][j] > 1.0e-90)
      {
        dout << i << " " << j << " " << ndist << " " << mdist << " " << mass_sum[i][j] << endl;
      }
    }
  }
  dout.close();

  stringstream distfile;
  distfile << "dist_" << setw(4) << setfill('0') << checkpoint << ".dat";

  ofstream fout(distfile.str().c_str());
  for (int i = 0; i < nbins; i++)
  {
    rdist = min_dist * std::pow((max_dist/min_dist), (double(i)+0.5)/double(nbins-1));
    
    for (int j = 0; j < nbins; j++)
    {
      ndist = min_dens * std::pow((max_dens/min_dens), (double(j)+0.5)/double(nbins-1));
      
      if (dist_sum[i][j] > 1.0e-90)
      {
        fout << i << " " << j << " " << rdist << " " << ndist << " " << dist_sum[i][j] << endl;
      }
    }
  }
  fout.close();
  
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
