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
#include "pointbox.h"
#include "kdtree.h"

#define PC_TO_CM	3.08567758e18
#define CM_TO_PC	3.2407793e-19
#define MSUN		2.0e33
#define mH		1.67e-24
#define PI		3.14159265359
#define G		6.673e-8
#define kB		1.380658e-16


typedef struct
{
  double x, y, z, dens, temp, mass, vol, dist, vx, vy, vz, pde, pden;
  
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
  
  if (siminfo.get_integer_scalar("timestep_index") > 61)
  {
    sim_star_x = siminfo.get_real_scalar("sim_star_x");
    sim_star_y = siminfo.get_real_scalar("sim_star_y");
    sim_star_z = siminfo.get_real_scalar("sim_star_z");
  } else {
    sim_star_x = 0.0;
    sim_star_y = 0.0;
    sim_star_z = 0.0;
  }
  
  cerr << "Star particle : " << sim_star_x << " " << sim_star_y;
  cerr << " " << sim_star_z << endl;
  
  //cout << "Sim time : " << siminfo.get_sim_time() << endl;
  
  //////////////////////////////////////////////////////////////
  
  // load the grid cells into a points array
  const unsigned int buffer_size = 512;
  unsigned int cache_size = 0;

  // Get mesh info
  const QuickFlash::File::MeshInfo & meshinfo = dfile.get_mesh_info();

  const unsigned int num_blocks = meshinfo.get_num_blocks();

  cerr << "Num blocks [ " << num_blocks << " ]" << endl;

  // Open a dataset
  const QuickFlash::File::Dataset & dset_dens
    = dfile.get_dataset("dens", buffer_size);

  dset_dens.set_cache_size(cache_size);
  dset_dens.set_report_stats();

  const QuickFlash::File::Dataset & dset_pde
    = dfile.get_dataset("pde", buffer_size);

  dset_pde.set_cache_size(cache_size);
  dset_pde.set_report_stats();

  const QuickFlash::File::Dataset & dset_pden
    = dfile.get_dataset("pden", buffer_size);

  dset_pden.set_cache_size(cache_size);
  dset_pden.set_report_stats();

  const QuickFlash::File::Dataset & dset_temp
    = dfile.get_dataset("temp", buffer_size);

  dset_temp.set_cache_size(cache_size);
  dset_temp.set_report_stats();

  const QuickFlash::File::Dataset & dset_vx
    = dfile.get_dataset("velx", buffer_size);

  dset_vx.set_cache_size(cache_size);
  dset_vx.set_report_stats();

  const QuickFlash::File::Dataset & dset_vy
    = dfile.get_dataset("vely", buffer_size);

  dset_vy.set_cache_size(cache_size);
  dset_vy.set_report_stats();

  const QuickFlash::File::Dataset & dset_vz
    = dfile.get_dataset("velz", buffer_size);

  dset_vz.set_cache_size(cache_size);
  dset_vz.set_report_stats();

  cerr << "Dataset loaded" << endl;

  // Set up temporary storage for cell data
  QuickFlash::Block::BlockData<double> data_dens;
  QuickFlash::Block::BlockData<double> data_pde;
  QuickFlash::Block::BlockData<double> data_pden;
  QuickFlash::Block::BlockData<double> data_temp;
  QuickFlash::Block::BlockData<double> data_vx;
  QuickFlash::Block::BlockData<double> data_vy;
  QuickFlash::Block::BlockData<double> data_vz;
  std::vector<double> cell_center;
  std::vector<cdata> cell_data;
  
  unsigned int cell_count = 0;
  
  /////////////////////////////////////////////////////////////////////
  // Find the density maximum location
  double dens_max = 0.0;
  double pde_max = 0.0;
  double pden_max = 0.0;
  double dens_max_x, dens_max_y, dens_max_z, dens_max_dist;
  double pde_max_x, pde_max_y, pde_max_z, pde_max_dist;
  double pden_max_x, pden_max_y, pden_max_z, pden_max_dist;
  
  double dx, dy, dz, r2;

  // Run through the leaf nodes
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
  {
    if(meshinfo.is_leaf(block_index))
    {
      // Get the block
      const QuickFlash::Block::BlockInfo & block_info 
	= meshinfo.get_block_info(block_index);

      dset_dens.get_block_data(block_index, data_dens);
      dset_pde.get_block_data(block_index, data_pde);
      dset_pden.get_block_data(block_index, data_pden);
      dset_temp.get_block_data(block_index, data_temp);
      dset_vx.get_block_data(block_index, data_vx);
      dset_vy.get_block_data(block_index, data_vy);
      dset_vz.get_block_data(block_index, data_vz);

      const unsigned int num_cells = block_info.get_num_cells();
      const double cell_vol = block_info.get_cell_volume();
      
      for (unsigned int index = 0; index < num_cells; index++)
      {
	block_info.get_cell_center(index, cell_center);
	
	dx = scale * (cell_center[0] - sim_star_x);
	dy = scale * (cell_center[1] - sim_star_y);
	dz = scale * (cell_center[2] - sim_star_z);
	r2 = dx*dx + dy*dy + dz*dz;
	
	if (r2 > SAMPLE_RADIUS_SQ) continue;
        
        cell_count++;
        
        cdata c;
        
        c.x = cell_center[0];
        c.y = cell_center[1];
        c.z = cell_center[2];
        
        c.vol = cell_vol * scale3;
	c.dens = data_dens[index] / scale3;
	c.pde = data_pde[index] / scale3;
	c.pden = data_pden[index] / scale3;
	c.temp = data_temp[index] * scale2;
	c.mass = data_dens[index] * cell_vol;
	
	c.vx = data_vx[index] * scale;
	c.vy = data_vy[index] * scale;
	c.vz = data_vz[index] * scale;
	
	cell_data.push_back(c);
	
	if (c.dens > dens_max)
	{
	  dens_max = c.dens;
	  dens_max_x = cell_center[0];
	  dens_max_y = cell_center[1];
	  dens_max_z = cell_center[2];
	  dens_max_dist = sqrt(r2);
        }

	if (c.pde > pde_max)
	{
	  pde_max = c.pde;
	  pde_max_x = cell_center[0];
	  pde_max_y = cell_center[1];
	  pde_max_z = cell_center[2];
	  pde_max_dist = sqrt(r2);
        }

	if (c.pden > pden_max)
	{
	  pden_max = c.pden;
	  pden_max_x = cell_center[0];
	  pden_max_y = cell_center[1];
	  pden_max_z = cell_center[2];
	  pden_max_dist = sqrt(r2);
        }
      }
    }
  }
  
  double vx, vy, vz, mass_sum = 0.0;
  
  for (unsigned int index = 0; index < cell_count; index++)
  {
    vx += cell_data[index].mass * cell_data[index].vx;
    vy += cell_data[index].mass * cell_data[index].vy;
    vz += cell_data[index].mass * cell_data[index].vz;
    mass_sum += cell_data[index].mass;
  }
  
  vx = vx / mass_sum;
  vy = vy / mass_sum;
  vz = vz / mass_sum;
  
  dx = scale * (pde_max_x - dens_max_x);
  dy = scale * (pde_max_y - dens_max_y);
  dz = scale * (pde_max_z - dens_max_z);
  r2 = dx*dx + dy*dy + dz*dz;
  pde_max_dist = sqrt(r2);

  dx = scale * (pden_max_x - dens_max_x);
  dy = scale * (pden_max_y - dens_max_y);
  dz = scale * (pden_max_z - dens_max_z);
  r2 = dx*dx + dy*dy + dz*dz;
  pden_max_dist = sqrt(r2);

  cerr << "Dataset cell count [ " << cell_count << " ]" << endl;
  cerr << "Density max: " << dens_max << endl;
  cerr << "  peak distance from star: " << dens_max_dist*CM_TO_PC << endl;

  cerr << "PDE max: " << pde_max << endl;
  cerr << "  peak distance from dens_max: " << pde_max_dist*CM_TO_PC << endl;

  cerr << "PDEN max: " << pden_max << endl;
  cerr << "  peak distance from dens_max: " << pden_max_dist*CM_TO_PC << endl;

  cerr << "Bulk velocity: " << sqrt(vx*vx + vy*vy + vz*vz)/1.0e5 << " km/s" << endl;

  /////////////////////////////////////////////////////////////////////
  // radially bin the gas around the density maximum

  // lets just sort the arrays by cell_dist and read off the items between r and r+dr

  // first calculate distance to the gas density maximum
  for (unsigned int index = 0; index < cell_count; index++)
  {
    dx = scale * (cell_data[index].x - dens_max_x);
    dy = scale * (cell_data[index].y - dens_max_y);
    dz = scale * (cell_data[index].z - dens_max_z);
    
    r2 = dx*dx + dy*dy + dz*dz;
    
    cell_data[index].dist = sqrt(r2);
  }
  
  // now sort by distance
  std::sort(cell_data.begin(), cell_data.end());
  
  // bin mass radially between bin_min and bin_max
  double bin_min = 3.0e17;
  double bin_max = 3.0e21;
  
  int nbins = 100;
  
  std::vector<double> mass, dens, pde, pden, temp, velocity;
  mass.resize(nbins);
  dens.resize(nbins);
  pde.resize(nbins);
  pden.resize(nbins);
  temp.resize(nbins);
  velocity.resize(nbins);
  
  double log_max_over_min = std::log(bin_max / bin_min);
  
  for (unsigned int i = 0; i < cell_count; i++)
  {
    int bin;

    if (cell_data[i].dist <= bin_min)
    {
      bin = 0;
    } else {
      bin = std::floor(double(nbins) * std::log(cell_data[i].dist / bin_min) / log_max_over_min) - 1;
    }
    
    mass[bin] += cell_data[i].mass;
    dens[bin] += cell_data[i].dens;
    pde[bin] += cell_data[i].pde;
    pden[bin] += cell_data[i].pden;
    temp[bin] += cell_data[i].temp * cell_data[i].mass;
    velocity[bin] += cell_data[i].mass * (std::pow((cell_data[i].vx - vx), 2.0) + std::pow((cell_data[i].vy - vy), 2.0) + std::pow((cell_data[i].vz - vz), 2.0));
    //vely[bin] += cell_data[i].mass * std::pow((cell_data[i].vy - vy), 2.0);
    //velz[bin] += cell_data[i].mass * std::pow((cell_data[i].vz - vz), 2.0);
    
  }

  double total_mass = 0.0;
  double total_temp = 0.0;
  double total_vol = 0.0;
  double total_velocity = 0.0;
  double avg_soundspeed = 0.0;
  double alt_jeans_mass = 0.0;

  // output binned data
  for (int i = 0; i < nbins; i++)
  {
    double inner, outer, vol, rms_velocity, jeans_mass, turb_jeans_mass;

    if (i == 0)
    {
      inner = 0.0;
      outer = bin_min;
    } else {
      inner = bin_min * std::pow((bin_max/bin_min), double(i)/double(nbins-1));
      outer = bin_min * std::pow((bin_max/bin_min), double(i+1)/double(nbins-1));
    }
    
    vol = 4.0*M_PI/3.0 * (outer*outer*outer - inner*inner*inner);
    
    //if (outer <= dens_max_dist)
    //{
      total_mass += mass[i];
      total_temp += temp[i];
      total_vol += vol;
      total_velocity += velocity[i];
      rms_velocity = total_velocity / total_mass; // < v^2 >
      avg_soundspeed = (kB * total_temp/total_mass) / mH; // c_s^2
    //}
    
    jeans_mass = 500.0 * std::pow((total_temp/total_mass/200.0),1.5) * std::pow((total_mass/total_vol/mH/1.0e4),-0.5) * MSUN;
    alt_jeans_mass = std::pow((avg_soundspeed), 1.5) * std::pow(G*G*G*(total_mass/total_vol),-0.5);
    turb_jeans_mass = std::pow((avg_soundspeed + rms_velocity/3.0), 1.5) * std::pow(G*G*G*(total_mass/total_vol),-0.5);
    
    cout << i << " " << inner << " " << outer << " " << mass[i] << " " << mass[i] / vol << " " << total_mass << " " 
         << temp[i]/mass[i] << " " << total_temp/total_mass << " " << jeans_mass << " " << sqrt(rms_velocity) << " " 
         << sqrt(avg_soundspeed) << " " << turb_jeans_mass << " " << alt_jeans_mass << " " << pde[i] << " " << pden[i] << endl;
         
    //cout << i << " " << outer << " " << total_mass << " " << jeans_mass << endl;
  }
  
  cerr << "Total Mass within distance from density peak to star:" << total_mass / 2e33 << endl;
}
