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

  double SAMPLE_RADIUS = 5000.0 * 3.0e18;
  double SAMPLE_RADIUS_SQ = SAMPLE_RADIUS*SAMPLE_RADIUS;
  double FOURTHIRDSPI = 4.0 * PI / 3.0;

  std::time_t start, stop;
  
  signal(SIGSEGV, traceback);

  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " filename" << endl
	 << endl;
    return 1;
  }

  int argptr = 1;

  const std::string file_name = argv[argptr++];

  // Open the file

  QuickFlash::File::DataFile dfile(file_name);

  if (!(dfile.particles_present()))
  {
    cerr << "Error: No particle data present" << endl;
    return -2;
  }

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
    sim_star_x = 1.6116045491e+24;
    sim_star_y = 1.4655038168e+24;
    sim_star_z = 1.6378772557e+24;
  }
  
  cout << "Star particle : " << sim_star_x << " " << sim_star_y;
  cout << " " << sim_star_z << endl;
  
  cout << "Sim time : " << siminfo.get_sim_time() << endl;
  
  std::vector<double> ejecta_mass;
  ejecta_mass.push_back(2.5);
  ejecta_mass.push_back(2.5);
  ejecta_mass.push_back(2.5);
  ejecta_mass.push_back(2.5);
  ejecta_mass.push_back(2.5);
  ejecta_mass.push_back(2.35);
  ejecta_mass.push_back(1.85);
  
  //double total_ejecta_mass = 16.7;

  // Get the particle data object
  const QuickFlash::Particles::PartData & partdata = dfile.get_part_data();

  // Output some basic information

  cout << "Particles : " << partdata.get_num_particles() << endl;

  //cout << "Variables : " << partdata.get_num_variables() << endl;
  //const std::vector<std::string> & var_names = partdata.get_variable_names();
  //const unsigned int num_variables = var_names.size();
  //for (unsigned int index = 0; index < num_variables; index++)
  //  cout << "  " << index << " : " << var_names[index] << endl;

  //////////////////////////////////////////////////////////////
  
  // load the grid cells into a points array
  const unsigned int buffer_size = 512;
  unsigned int cache_size = 0;

  // Open the file
  //QuickFlash::File::DataFile dfile(file_name);

  // Get mesh info
  const QuickFlash::File::MeshInfo & meshinfo = dfile.get_mesh_info();

  const unsigned int num_blocks = meshinfo.get_num_blocks();
  //const unsigned int dims = meshinfo.get_dims();

  cout << "Num blocks [ " << num_blocks << " ]" << endl;

  // Get information
  //const std::vector<unsigned int> & block_dims = meshinfo.get_block_dims();

  //const std::vector<unsigned int> & base_block_dims 
  //  = meshinfo.get_base_block_dims();

  // Open a dataset
  const QuickFlash::File::Dataset & dset_dens
    = dfile.get_dataset("dens", buffer_size);

  const QuickFlash::File::Dataset & dset_metal
    = dfile.get_dataset("z", buffer_size);

  const QuickFlash::File::Dataset & dset_temp
    = dfile.get_dataset("temp", buffer_size);

  dset_dens.set_cache_size(cache_size);
  dset_dens.set_report_stats();

  dset_metal.set_cache_size(cache_size);
  dset_metal.set_report_stats();

  dset_temp.set_cache_size(cache_size);
  dset_temp.set_report_stats();

  cout << "Dataset loaded" << endl;

  // Set up temporary storage for cell data
  unsigned int num_sn = 7;
  std::vector< Point<3> > amrpoints;
  QuickFlash::Block::BlockData<double> data_dens, data_metal, data_temp;
  std::vector<double> cell_center, cell_x, cell_y, cell_z, cell_density, 
                      cell_mass, cell_metal, cell_temp, cell_volume, cell_parts[num_sn];
  
  unsigned int cell_count = 0;
  
  double dens_min = 1.0e99;
  double dens_max = 0.0;

  // Run through the leaf nodes
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
  {
    if(meshinfo.is_leaf(block_index))
    {
      // Get the block
      const QuickFlash::Block::BlockInfo & block_info 
	= meshinfo.get_block_info(block_index);

      dset_dens.get_block_data(block_index, data_dens);
      dset_metal.get_block_data(block_index, data_metal);
      //dset_temp.get_block_data(block_index, data_temp);

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
        
        cell_x.push_back(dx);
        cell_y.push_back(dy);
        cell_z.push_back(dz);
        
        cell_volume.push_back(cell_vol * scale3);
	cell_density.push_back(data_dens[index] / scale3);
	cell_mass.push_back(data_dens[index] * cell_vol);
	cell_metal.push_back(data_metal[index]);
	//cell_temp.push_back(data_temp[index] * scale2);
	
	if (cell_density.back() > dens_max) dens_max = cell_density.back();
	if (cell_density.back() < dens_min) dens_min = cell_density.back();
	
	for (unsigned int which = 0; which < num_sn; which++)
	{
	  cell_parts[which].push_back(0.0);
	}

        amrpoints.push_back( Point<3>(dx, dy, dz) );

      }
    }
  }

  cout << "Dataset cell count [ " << cell_count << " ]" << endl;
  
  KDtree<3> *amrtree = new KDtree<3>(amrpoints);
  
  cout << "AMR data loaded into KDtree" << endl;
  
  cout << "Density min/max: " << dens_min << " / " << dens_max << endl;


  ///////////////////////////////////////

  // load particle data into points vectors
  std::vector< Point<3> > points[num_sn], allpoints;
  unsigned int num_particles = partdata.get_num_particles();
  int count = 0;
  int dm_count = 0;
  int sn_count[num_sn];
  
  for (unsigned int index = 0; index < num_sn; index++)
  {
    sn_count[index] = 0;
  }
  
  std::vector<double> types, tags, posx, posy, posz;
  
  start = std::time(0);
  cout << "Loading data... " << ctime(&start);
  
  partdata.get_variable_data("type", types);
  partdata.get_variable_data("tag", tags);
  partdata.get_variable_data("posx", posx);
  partdata.get_variable_data("posy", posy);
  partdata.get_variable_data("posz", posz);

  stop = std::time(0);
  double duration = difftime(stop, start);
  cout << "Done loading data... " << ctime(&stop);
  cout << "  duration.......... " << duration << endl;

  // Load point vectors for each supernova
  for (unsigned int index = 0; index < num_particles; index++)
  {
    if (types[index] == 1.0)
    {
      double dx = scale * (posx[index] - sim_star_x);
      double dy = scale * (posy[index] - sim_star_y);
      double dz = scale * (posz[index] - sim_star_z);
      double r2 = dx*dx + dy*dy + dz*dz;
      unsigned int tag = (unsigned int)(tags[index]) - 1;
      
      points[tag].push_back( Point<3>(dx, dy, dz) );
      allpoints.push_back( Point<3>(dx, dy, dz) );

      if (r2 <= SAMPLE_RADIUS_SQ)
      {
        int c = amrtree->nearest( Point<3>(dx, dy, dz) );
        cell_parts[tag][c] += 1.0;

        count++;
        sn_count[tag]++;
        
        if (!(count == 0) && !(count & (count-1)))
          cout << "count: " << count << ' ' << c << endl;
      }
    } else {
      dm_count++;
    }
  }
  
  for (unsigned int index = 0; index < num_sn; index++)
  {
    ejecta_mass[index] *= MSUN / (double)(points[index].size());
    cout << "SN counts: " << sn_count[index] << " with total mass " << sn_count[index] * ejecta_mass[index] << " / " << points[index].size() * ejecta_mass[index] << endl;
  }
  
  for (unsigned int index = 0; index < amrpoints.size(); index++)
  {
    for (unsigned int which = 0; which < num_sn; which++)
    {
      cell_parts[which][index] *= ejecta_mass[which];
    }
  }
  
  ofstream dout("data_output.txt");
  for (unsigned int index = 0; index < amrpoints.size(); index++)
  {
    dout << cell_x[index] << ' ' << cell_y[index] << ' ' << cell_z[index] << ' ';
    dout << cell_density[index] << ' ' << cell_volume[index] << ' ' << cell_metal[index];;
  
    for (unsigned int which = 0; which < num_sn; which++)
    {
      dout << ' ' << cell_parts[which][index];
    }
    dout << endl;
  }
  dout.close();
  
  return 0;
  

}
