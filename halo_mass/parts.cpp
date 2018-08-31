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

  double SAMPLE_RADIUS = 6.0e21;
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
  
  double sim_star_x = siminfo.get_real_scalar("sim_star_x");
  double sim_star_y = siminfo.get_real_scalar("sim_star_y");
  double sim_star_z = siminfo.get_real_scalar("sim_star_z");
  
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
  
#if 0

  // Set up temporary storage for cell data
  std::vector< Point<3> > amrpoints;
  QuickFlash::Block::BlockData<double> data_dens, data_metal, data_temp;
  std::vector<double> cell_center, cell_x, cell_y, cell_z, cell_density, cell_mass, cell_metal, cell_temp, cell_volume;
  
  unsigned int cell_count = 0;

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
      dset_temp.get_block_data(block_index, data_temp);

      const unsigned int num_cells = block_info.get_num_cells();
      const double cell_vol = block_info.get_cell_volume();

      for (unsigned int index = 0; index < num_cells; index++)
      {
	block_info.get_cell_center(index, cell_center);
	
	double dx = scale * (cell_center[0] - sim_star_x);
	double dy = scale * (cell_center[1] - sim_star_y);
	double dz = scale * (cell_center[2] - sim_star_z);
	
	double r = std::sqrt(dx*dx + dy*dy + dz*dz);
	
	if (r > SAMPLE_RADIUS) continue;
        
        cell_count++;
        
        cell_x.push_back(dx);
        cell_y.push_back(dy);
        cell_z.push_back(dz);
        
        cell_volume.push_back(cell_vol * scale3);
	cell_density.push_back(data_dens[index] / scale3);
	cell_mass.push_back(data_dens[index] * cell_vol);
	cell_metal.push_back(data_metal[index]);
	cell_temp.push_back(data_temp[index] * scale2);

        amrpoints.push_back( Point<3>(cell_center[0], cell_center[1], cell_center[2]) );

      }
    }
  }

  cout << "Dataset cell count [ " << cell_count << " ]" << endl;
  
  KDtree<3> *amrtree = new KDtree<3>(amrpoints);
  
  cout << "AMR data loaded into KDtree" << endl;


  ///////////////////////////////////////
#endif

  // load particle data into points vectors
  unsigned int num_sn = 7;
  std::vector< Point<3> > points[num_sn], allpoints;
  unsigned int num_particles = partdata.get_num_particles();
  int count = 0;
  int dm_count = 0;
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
      count++;

      unsigned int tag = (unsigned int)(tags[index]) - 1;

      points[tag].push_back( Point<3>(posx[index], posy[index], posz[index]) );
      allpoints.push_back( Point<3>(posx[index], posy[index], posz[index]) );

    } else {
      dm_count++;
    }
  }
  
  cout << "Particle counts: DM " << dm_count << " and metal " << count << endl;
  
  // create KDtrees for each supernova
  KDtree<3> *tree[num_sn], *alltree;

  alltree = new KDtree<3>(allpoints);

  for (unsigned int index = 0; index < num_sn; index++)
  {
    cout << "SN" << index+1 << " " << points[index].size() << endl;
    
    tree[index] = new KDtree<3>(points[index]);
  }
  
  ///////////////////////////////////////
  
#if 0  
  // calculate densities at amr points
  std::vector<double> density[num_sn], alldensity;
  
  ofstream fout("output_data.txt");
  
  int nsamples = 10000;
  vector<unsigned int> random_index;
  for (int i = 0; i < nsamples; i++)
  {
    double rand = drand48() * amrtree->npts;
    random_index.push_back(rand);
  }
  
  double metal_mass_cutoff = 1.0e-6 * MSUN;
  
  vector<double> mass_cuts;
  
  mass_cuts.push_back(1.0e-6 * MSUN);
  mass_cuts.push_back(1.0e-5 * MSUN);
  mass_cuts.push_back(1.0e-4 * MSUN);
  mass_cuts.push_back(1.0e-3 * MSUN);
  mass_cuts.push_back(1.0e-2 * MSUN);
  mass_cuts.push_back(1.0e-1 * MSUN);
  
  int mass_counts[6];
  for (int i = 0; i < 6; i++)
  {
    mass_counts[i] = 0;
  }
  
  for (int i = 0; i < amrtree->npts; i++)
  {
    double metal_mass = cell_mass[i] * cell_metal[i];
    for (unsigned int j = 0; j < mass_cuts.size(); j++)
    {
      if (metal_mass > mass_cuts[j])
      {
        mass_counts[j]++;
      }
    }
  }
  
  for (int i = 0; i < 6; i++)
  {
    cout << "Mass cut " << mass_cuts[i] << " has " << mass_counts[i] << " cells." << endl;
  }
  
  for (int i = 0; i < amrtree->npts; i++)
//  for (int index = 0; index < nsamples; index++)
  {
  
//    const unsigned int i = random_index[index];
    int n = 10;//, neghmax = 1000;
    int negh[n];//, neghlist[neghmax];
    double d[n];
    double dst[num_sn];
    double max_dst = 0.0;
    double sum = 0.0;
    
    if (!(i % (amrtree->npts)))
    {
      cout << "Progress " << i << " out of " << amrtree->npts << endl;
    }
    
    double local_metal_mass = cell_mass[i] * cell_metal[i];
    
    if (local_metal_mass > metal_mass_cutoff)
    {
      fout << setw(20) << i;
      fout << setw(20) << cell_x[i];
      fout << setw(20) << cell_y[i];
      fout << setw(20) << cell_z[i];
      fout << setw(20) << cell_volume[i];
      fout << setw(20) << cell_density[i];
      fout << setw(20) << cell_mass[i];
      fout << setw(20) << cell_metal[i];
      fout << setw(20) << cell_temp[i];

      for (unsigned int sn = 0; sn < num_sn; sn++)
      {
        tree[sn]->nnearest(amrpoints[i], negh, d, n);
        
        dst[sn] = d[0];

        for (int k = 1; k < n; k++)
        {
          if (d[k] > dst[sn])
            dst[sn] = d[k];
        }
        
        if (dst[sn] > max_dst)
          max_dst = dst[sn];
        
        max_dst *= scale;
        
        double dens = (double)(n) * MSUN * ejecta_mass[sn] / (FOURTHIRDSPI * max_dst * max_dst * max_dst) / (double)(tree[sn]->npts);
        //double npart = tree[sn]->locatenear(amrtree->ptss[i], max_dst, neghlist, neghmax);

        density[sn].push_back(dens);
        sum += dens;
        
        fout << setw(20) << dens;
      }
      
      fout << setw(20) << sum << endl;
      
      
      alldensity.push_back(sum);
    }
  }
  
  fout.close();
  
  cout << "Wrote data to file..." << endl;
  
  return 0;
  
  // variance and mean
  double mean[num_sn];
  double variance[num_sn];
  double allmean = 0.0;
  double allvar = 0.0;
  double dmin = 1.0e99;
  double dmax = 0.0;
  
  for (unsigned int i = 0; i < num_sn; i++)
  {
    mean[i] = 0.0;
    variance[i] = 0.0;
  }
  
  for (int i = 0; i < amrtree->npts; i++)
  {
    for (unsigned int j = 0; j < num_sn; j++)
    {
      if (density[j][i] > dmax)
        dmax = density[j][i];
      if (density[j][i] < dmin)
        dmin = density[j][i];
      
      mean[j] += density[j][i];
      allmean += density[j][i];
    }
  }
  
  allmean /= (double)(alltree->npts);
  cout << "allmean = " << allmean << endl;
  
  for (unsigned int i = 0; i < num_sn; i++)
  {
    mean[i] /= (double)(tree[i]->npts);
    cout << "avg[" << i << "] = " << mean[i] << endl;
  }
  
  for (int i = 0; i < amrtree->npts; i++)
  {
    for (unsigned int j = 0; j < num_sn; j++)
    {
      variance[j] += (density[j][i] - mean[j]) * (density[j][i] - mean[j]);
      allvar += (density[j][i] - allmean) * (density[j][i] - allmean);
    }
  }
  
  allvar /= (double)(alltree->npts);
  cout << "allvar = " << allvar << endl;
  
  for (unsigned int i = 0; i < num_sn; i++)
  {
    variance[i] /= (double)(tree[i]->npts);
    cout << "var[" << i << "] = " << variance[i] << endl;
  }
  

  cout << "min = " << dmin << endl;
  cout << "max = " << dmax << endl;
#endif
  return 0;
}
