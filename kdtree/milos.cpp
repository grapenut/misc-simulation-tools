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

inline double sq(double x) { return x * x; }
inline double cub(double x) { return x * x * x; }

inline double Epanechnikov(double x, double kernel_radius) 
{
  if(x >= kernel_radius) return 0;
  else return 1.0 - sq(x / kernel_radius); 
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

struct chiHistogram : vector<double> 
{
  double chi_min, chi_max, log_range;

  chiHistogram(double _chi_min, double _chi_max, size_t _chi_bins) :
    vector<double> (_chi_bins, 0.0), chi_min(_chi_min), chi_max(_chi_max), log_range(std::log(chi_max / chi_min)) { }

  size_t bin(double chi) const
  {
    if (chi < chi_min) 
      return 0; 
    else if (chi >= chi_max) 
      return size() - 1; 
    else
      return size_t(std::floor(double(size() - 1) * std::log(chi / chi_min) / log_range) + 1.0);
  }

  void addMass(double chi, double mass) 
  {
    at(bin(chi)) += mass;
  }

  void normalize() 
  {
    const double total_mass = std::accumulate(begin(), end(), 0.0);

    if(total_mass > 0.0) std::transform(begin(), end(), begin(), std::bind2nd(std::divides<double> (), total_mass));
  }

  bool isEmpty() const { return find_if(begin(), end(), bind2nd(greater<double> (), 0.0)) == end(); }

  double binCenter(size_t i) const
  {
    if(i == 0) 
      return 0.5 * chi_min;
    else
      return chi_min * std::exp((double(i) - 0.5) * log_range / double(size() - 1));
  }
  double binLeft(size_t i) const
  {
    if(i == 0) 
      return 0.0;
    else
      return chi_min * std::exp((double(i) - 1.0) * log_range / double(size() - 1));
  }
  double binRight(size_t i) const
  {
    return chi_min * std::exp((double(i)) * log_range / double(size() - 1));
  }
};

ostream & operator<< (ostream & out, const chiHistogram & hist)
{
  for(size_t i = 0; i < hist.size(); i ++)
    out << hist.binCenter(i) << " " << hist.binLeft(i) << " " << hist.binRight(i) << " " << hist[i] << endl;

  return out;
}

inline bool is_power_2(size_t x)
{
  return x != 0 && !(x & (x - 1));
}

struct Neighbor
{
  size_t tag;
  double distance, mass;
  Neighbor(size_t _tag, double _distance, double _mass) : tag(_tag), distance(_distance), mass(_mass) { }
};

bool operator< (const Neighbor & a, const Neighbor & b)
{
  return a.distance < b.distance;
}

double massWithinKernel(const vector<Neighbor> & neighbors, double kernel_radius)
{
  double mass = 0.0;

  for(vector<Neighbor>::const_iterator i = neighbors.begin(); i != neighbors.end(); i ++)
    {
      if(i->distance >= kernel_radius) break;

      mass += i->mass * Epanechnikov(i->distance, kernel_radius);
    }

  return mass;
}

double kernelRadius(const vector<Neighbor> & neighbors, double target_mass)
{
  const double tolerance = 1.0e-3;

  double rmin = neighbors.front().distance;
  double rmax = neighbors.back().distance;
 
  while(rmax > (1.0 + tolerance) * rmin)
    {
      const double rmid = 0.5 * (rmin + rmax);
      
      const double mass = massWithinKernel(neighbors, rmid);

      if(mass < target_mass) rmin = rmid;
      else rmax = rmid;
    }

  return 0.5 * (rmin + rmax);
}

int main(int argc, char * argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  
  const double kpc = 3.0e21;
  const double SAMPLE_RADIUS = 5.0 * kpc;
  const double SAMPLE_RADIUS_SQ = SAMPLE_RADIUS * SAMPLE_RADIUS;
  const double FOURTHIRDSPI = 4.0 * M_PI / 3.0;

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

  //  const QuickFlash::File::Dataset & dset_temp
  //    = dfile.get_dataset("temp", buffer_size);

  if(false)
    {

      dset_dens.set_cache_size(cache_size);
      dset_dens.set_report_stats();

      dset_metal.set_cache_size(cache_size);
      dset_metal.set_report_stats();

      //  dset_temp.set_cache_size(cache_size);
      //  dset_temp.set_report_stats();

      cout << "Dataset loaded" << endl;

      // Set up temporary storage for cell data
      std::vector< Point<3> > amrpoints;
      QuickFlash::Block::BlockData<double> data_dens, data_metal; //, data_temp;
      std::vector<double> cell_center, cell_x, cell_y, cell_z, cell_density, cell_mass, cell_metal, cell_volume; // cell_temp;
  
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
	      // dset_temp.get_block_data(block_index, data_temp);

	      const unsigned int num_cells = block_info.get_num_cells();
	      const double cell_vol = block_info.get_cell_volume();

	      for (unsigned int index = 0; index < num_cells; index++)
		{
		  block_info.get_cell_center(index, cell_center);
	
		  const double dx = scale * (cell_center[0] - sim_star_x);
		  const double dy = scale * (cell_center[1] - sim_star_y);
		  const double dz = scale * (cell_center[2] - sim_star_z);
	
		  const double r_sq = dx*dx + dy*dy + dz*dz;
	
		  if (r_sq > SAMPLE_RADIUS_SQ) continue;

		  cell_count++;
        
		  cell_x.push_back(dx);
		  cell_y.push_back(dy);
		  cell_z.push_back(dz);
        
		  cell_volume.push_back(cell_vol * scale3);
		  cell_density.push_back(data_dens[index] / scale3);
		  cell_mass.push_back(data_dens[index] * cell_vol);
		  cell_metal.push_back(data_metal[index]);
		  //	cell_temp.push_back(data_temp[index] * scale2);

		  amrpoints.push_back( Point<3>(dx, dy, dz) );

		}
	    }
	}

      cout << "Dataset cell count [ " << cell_count << " ]" << endl;
  
      KDtree<3> * amrtree = new KDtree<3> (amrpoints);
  
      cout << "AMR data loaded into KDtree" << endl;

    }


  ///////////////////////////////////////

  // load particle data into points vectors 
  
  const unsigned int num_sn = 7;
  const double ejecta_mass[num_sn] = { 2.5, 2.5, 2.5, 2.5, 2.5, 2.35, 1.85 };
  
  //double total_ejecta_mass = 16.7;
 
  std::vector< Point<3> > points[num_sn], allpoints;
  std::vector<size_t> alltags;
  
  // Get the particle data object
  const QuickFlash::Particles::PartData & partdata = dfile.get_part_data();

  // Output some basic information

  cout << "Particles : " << partdata.get_num_particles() << endl;

  const unsigned int num_particles = partdata.get_num_particles();

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
  const double duration = difftime(stop, start);
  cout << "Done loading data... " << ctime(&stop);
  cout << "  duration.......... " << duration << endl;

  // Load point vectors for each supernova
  for (unsigned int index = 0; index < num_particles; index++)
    {
      if (types[index] == 1.0) // metal tracer particles
	{
	  count++;

	  const unsigned int tag = (unsigned int)(tags[index]) - 1;

	  const double dx = scale * (posx[index] - sim_star_x);
	  const double dy = scale * (posy[index] - sim_star_y);
	  const double dz = scale * (posz[index] - sim_star_z);

	  allpoints.push_back( Point<3> (dx, dy, dz) );
	  alltags.push_back(tag);
	  points[tag].push_back( allpoints.back() );

	} else {
	dm_count++;
      }
    }
  
  double particle_mass[num_sn];
  for(size_t tag = 0; tag < num_sn; tag ++)
    {
      particle_mass[tag] = ejecta_mass[tag] / double(points[tag].size());
    }

  cout << "Particle counts: DM " << dm_count << " and metal " << count << endl;

  cout << "Metal particle masses = " << endl;
  std::copy(particle_mass, particle_mass + num_sn, std::ostream_iterator<double> (cout, " "));
  cout << endl;
  
  // create KDtrees for each supernova
  KDtree<3> *tree[num_sn], *alltree;

  alltree = new KDtree<3>(allpoints);

  for (unsigned int index = 0; index < num_sn; index++)
    {
      cout << "SN" << index+1 << " " << points[index].size() << endl;
    
      tree[index] = new KDtree<3>(points[index]);
    }
  
  cout << "Constructed KD trees" << endl;

  ///////////////////////////////////////
  
  const double chi_min = 1.0e-2;
  const double chi_max = 1.0;
  const size_t chi_bins = 40;

  chiHistogram * hist[num_sn][num_sn];
  double sum_chi_mass[num_sn][num_sn];
  double sum_chi_sq_mass[num_sn][num_sn];
  double sum_mass[num_sn][num_sn];

  for(size_t i = 0; i < num_sn; i ++)
    for(size_t j = 0; j < num_sn; j ++)
      if(i != j)
	{
	  hist[i][j] = new chiHistogram(chi_min, chi_max, chi_bins);
	  sum_chi_mass[i][j] = 0.0;
	  sum_chi_sq_mass[i][j] = 0.0;
	  sum_mass[i][j] = 0.0;
	}

  const double reference_metal_mass = 4.0e-5; // solar masses

  const size_t sn_max = 2;
  const size_t sn_stride = 1;

  const double safety = 15.0;

  for(size_t sn_1 = 0; sn_1 < sn_max; sn_1 += sn_stride)
    {
      const size_t num_neigh_1 = size_t(std::ceil(safety * reference_metal_mass / particle_mass[sn_1]));

      cout << "num_neigh[" << sn_1 + 1 << "] = " << num_neigh_1 << endl;

      int * neighbors_1 = new int [num_neigh_1];
      double * distances_1 = new double [num_neigh_1];

      int count_point = 0;
     
      for(vector<Point<3> >::const_iterator p = points[sn_1].begin(); p != points[sn_1].end(); p ++)
	{
	  if(is_power_2(++ count_point)) cout << count_point << " " << std::flush;

	  tree[sn_1] -> nnearest(* p, neighbors_1, distances_1, num_neigh_1);

	  std::vector<Neighbor> all_neighbors_1;
	  all_neighbors_1.reserve(num_neigh_1);
	  
	  for(size_t i = 0; i < num_neigh_1; i ++)
	    all_neighbors_1.push_back(Neighbor(sn_1, distances_1[i], particle_mass[sn_1]));

	  std::sort(all_neighbors_1.begin(), all_neighbors_1.end());
	 
	  const double mass_within_max_kernel_1 = massWithinKernel(all_neighbors_1, all_neighbors_1.back().distance);
	  if(mass_within_max_kernel_1 < reference_metal_mass)
	    {	 
	      cout << "Insufficent neighbors, increase safety, reference_metal mass = " << reference_metal_mass << ",  mass_within_max_kernel_1 = " << mass_within_max_kernel_1 << endl;
	      std::abort();
	    }

	 
	  for(size_t sn_2 = 0; sn_2 < sn_max; sn_2 += sn_stride)
	    if(sn_1 != sn_2)
	      {
		const size_t num_neigh_2 = size_t(std::ceil(safety * reference_metal_mass / particle_mass[sn_2]));
		int * neighbors_2 = new int [num_neigh_2];
		double * distances_2 = new double [num_neigh_2];
	  
		tree[sn_2] -> nnearest(* p, neighbors_2, distances_2, num_neigh_2);

		std::vector<Neighbor> all_neighbors_2;
		all_neighbors_2.reserve(num_neigh_2);
		
	        for(size_t i = 0; i < num_neigh_2; i ++)
		  all_neighbors_2.push_back(Neighbor(sn_2, distances_2[i], particle_mass[sn_2]));

		delete [] neighbors_2;
		delete [] distances_2;

		std::sort(all_neighbors_2.begin(), all_neighbors_2.end());

		const double mass_within_max_kernel_2 = massWithinKernel(all_neighbors_2, all_neighbors_2.back().distance);
		if(mass_within_max_kernel_2 < reference_metal_mass)
		  {	 
		    cout << "Insufficent neighbors, increase safety, reference_metal mass = " << reference_metal_mass << ",  mass_within_max_kernel_2 = " << mass_within_max_kernel_2 << endl;
		    std::abort();
		  }

		vector<Neighbor> all_neighbors(all_neighbors_1);
		all_neighbors.insert(all_neighbors.end(), all_neighbors_2.begin(), all_neighbors_2.end()); 
		
		std::sort(all_neighbors.begin(), all_neighbors.end());
	      
		const double kernel_radius = kernelRadius(all_neighbors, reference_metal_mass);
	
		const double mass_1 = massWithinKernel(all_neighbors_1, kernel_radius);
		const double mass_2 = massWithinKernel(all_neighbors_2, kernel_radius);
	
		const double chi_12 = mass_1 / (mass_1 + mass_2);
		const double chi_21 = 1.0 - chi_12;

		hist[sn_1][sn_2]->addMass(chi_12, particle_mass[sn_1]);
		hist[sn_2][sn_1]->addMass(chi_21, particle_mass[sn_1]);

		sum_chi_mass[sn_1][sn_2] += particle_mass[sn_1] * chi_12;
		sum_chi_sq_mass[sn_1][sn_2] += particle_mass[sn_1] * sq(chi_12);
		sum_mass[sn_1][sn_2] += particle_mass[sn_1];	

		sum_chi_mass[sn_2][sn_1] += particle_mass[sn_1] * chi_21;
		sum_chi_sq_mass[sn_2][sn_1] += particle_mass[sn_1] * sq(chi_21);
		sum_mass[sn_2][sn_1] += particle_mass[sn_1];
	      }
	}
      cout << endl;
	  
      delete [] neighbors_1;
      delete [] distances_1;
    }


  for(size_t sn_1 = 0; sn_1 < num_sn; sn_1 ++)
    for(size_t sn_2 = 0; sn_2 < num_sn; sn_2 ++)
      if(sn_1 != sn_2) 
	if(!hist[sn_1][sn_2]->isEmpty())
	  {
	    hist[sn_1][sn_2]->normalize();

	    ostringstream file_name;
	    file_name << "hist_" << sn_1 + 1 << "_" << sn_2 + 1 << ".dat" << ends;
	    ofstream out(file_name.str().c_str());
	    out << * hist[sn_1][sn_2];

	    const double chi_average = sum_chi_mass[sn_1][sn_2] / sum_mass[sn_1][sn_2];
	    const double chi_sq_average = sum_chi_sq_mass[sn_1][sn_2] / sum_mass[sn_1][sn_2];
	    const double chi_sigma = std::sqrt(chi_sq_average - chi_average * chi_average);
	    
	    const double fully_mixed_average = ejecta_mass[sn_1] / (ejecta_mass[sn_1] + ejecta_mass[sn_2]);
	    
	    cout << "statistics: " << sn_1 << " " << sn_2 << " " << fully_mixed_average << " " << chi_average << " " << chi_sigma << endl;
	  }

  return 0;
}
