// C++ program file halo.cpp
/*
    1) determine halo center
        find the maximum DM density within 100 pc of the star particle
    2) sort the DM particles by distance to the halo center
    3) iterate over the DM particles, calculating enclosed mass density
    4) determine where the enclosed density falls below 200*rho_critical
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
#define YR_TO_SEC	3.15569e7
#define SEC_TO_YR	3.16888e-8

#define MSUN		2.0e33
#define mH		1.67e-24
#define PI		3.14159265359
//#define SN_START_TIME	6.0158037347e15
#define SN_START_TIME	6.02855e+15

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

bool compare_pairs(const std::pair< double, double > &left, const std::pair< double, double > &right)
{
  return left.first < right.first;
}

int main(int argc, char * argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;

  double SAMPLE_RADIUS = 100.0 * 3.0e18;
  double SAMPLE_RADIUS_SQ = SAMPLE_RADIUS*SAMPLE_RADIUS;
  double FOURTHIRDSPI = 4.0 * PI / 3.0;

  std::time_t start, stop;


  // Initialize traceback handler for segfaults
  signal(SIGSEGV, traceback);


  // Parse checkpoint name from command line
  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " filename" << endl
	 << endl;
    return 1;
  }

  const std::string file_name = argv[1];


  // Open the checkpoint file
  QuickFlash::File::DataFile dfile(file_name);

  if (!(dfile.particles_present()))
  {
    cerr << "Error: No particle data present" << endl;
    return -2;
  }


  // Read the parameters for this checkpoint (redshift, scale, time, etc)
  const QuickFlash::File::SimInfo & siminfo = dfile.get_sim_info();
  
  double scale = siminfo.get_real_scalar("scalefactor");
  double scale2 = scale*scale;
  double scale3 = scale*scale*scale;
  double redshift = (1.0/scale) - 1.0;
  
  SAMPLE_RADIUS_SQ /= scale2;
  
  double sim_star_x, sim_star_y, sim_star_z;
  
  int checkpoint = siminfo.get_integer_scalar("checkpointfilenumber");
  if (checkpoint > 12)
  {
    sim_star_x = siminfo.get_real_scalar("sim_star_x");
    sim_star_y = siminfo.get_real_scalar("sim_star_y");
    sim_star_z = siminfo.get_real_scalar("sim_star_z");
  } else {
    //sim_star_x = 1.6116045491e+24;
    //sim_star_y = 1.4655038168e+24;
    //sim_star_z = 1.6378772557e+24;
    sim_star_x = 1.60934e+24;
    sim_star_y = 1.46397e+24;
    sim_star_z = 1.63851e+24;
  }
  
  //cout << "Star particle : " << sim_star_x << " " << sim_star_y << " " << sim_star_z << endl;
  //cout << "Sim time : " << siminfo.get_sim_time() << endl;


  ////////////////////////////////////////////////////////
  // Find maximum DM density within 100pc of star particle

  // Read the particle data object
  const QuickFlash::Particles::PartData & partdata = dfile.get_part_data();

  // Load particle data into points vectors
  std::vector< Point<3> > points;
  std::vector< double > masses;
  unsigned int num_particles = partdata.get_num_particles();
  unsigned int dm_count = 0;  
  double total_mass = 0.0;

  std::vector<double> ptype, pmass, px, py, pz;
  partdata.get_variable_data("type", ptype);
  partdata.get_variable_data("mass", pmass);
  partdata.get_variable_data("posx", px);
  partdata.get_variable_data("posy", py);
  partdata.get_variable_data("posz", pz);

  // Create points vector
  for (unsigned int index = 0; index < num_particles; index++)
  {
    if (ptype[index] == 2.0)
    {
      dm_count++;

      points.push_back( Point<3>(px[index], py[index], pz[index]) );
      masses.push_back(pmass[index]);
    
      total_mass += pmass[index];
    }
  }
  
  // total_mass / mpc^3
  double rho_average = total_mass / (1e18*PC_TO_CM*PC_TO_CM*PC_TO_CM);
  double rho_200 = 200.0 * rho_average;
  
  //cout << "Total DM particles: " << dm_count << endl;
  //cout << "Total Mass: " << total_mass << endl;
  //cout << "Rho Critical: " << rho_average << endl;
  
  // Create KDtree from points
  KDtree<3> *tree = new KDtree<3>(points);

  double max_dens = 0.0;
  unsigned int max_index;
  int n = 10;
  int negh[n];
  double d[n];
  
  // Loop through points, estimate density at each, find the maximum
  for (unsigned int index = 0; index < dm_count; index++)
  {
    double dx = points[index].x[0] - sim_star_x;
    double dy = points[index].x[1] - sim_star_y;
    double dz = points[index].x[2] - sim_star_z;
    double r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 > SAMPLE_RADIUS_SQ) continue;
    
    tree->nnearest(points[index], negh, d, n);
    
    double max_dist = d[0];
    double sum_mass = masses[index] + masses[negh[0]];

    for (int k = 1; k < n; k++)
    {
      sum_mass += masses[negh[k]];
      if (d[k] > max_dist)
        max_dist = d[k];
    }
    
    double density = sum_mass / (max_dist*max_dist*max_dist);
    if (density < max_dens) continue;
    
    max_dens = density;
    max_index = index;
  }
  
  double max_x = points[max_index].x[0];
  double max_y = points[max_index].x[1];
  double max_z = points[max_index].x[2];
  //cout << "Maximum DM density found: " << max_dens << " at " << max_x << " " << max_y << " " << max_z << endl;

  ///////////////////////////////////////

  // Create mass/radius pairs
  std::vector< std::pair< double, double > > parts;
  
  for (unsigned int index = 0; index < dm_count; index++)
  {
    double dx = points[index].x[0] - max_x;
    double dy = points[index].x[1] - max_y;
    double dz = points[index].x[2] - max_z;
    double r2 = dx*dx + dy*dy + dz*dz;
    
    parts.push_back(make_pair(r2, masses[index]));
  }

  std::sort(parts.begin(), parts.end(), compare_pairs);
  
  //cout << "First particle: " << std::sqrt(parts.front().first)*scale*CM_TO_PC << " and " << parts.front().second/MSUN << endl;
  //cout << " Last particle: " << std::sqrt(parts.back().first)*scale*CM_TO_PC << " and " << parts.back().second/MSUN << endl;
  
  double m200, r200, rho;
  int pcount = 1;
  
  m200 = parts[0].second;
  for (unsigned int index = 1; index < dm_count; index++)
  {
    pcount++;
    m200 += parts[index].second;
    r200 = std::sqrt(parts[index].first);
    rho = m200 / (FOURTHIRDSPI*r200*r200*r200);
    
    if (rho < rho_200) break;
  }

  double back_m200, back_r200;
  int back_pcount = dm_count;
  
  back_m200 = total_mass + parts[dm_count-1].second;
  for (int index = dm_count-1; index > -1; index--)
  {
    back_pcount--;
    back_m200 -= parts[index].second;
    back_r200 = std::sqrt(parts[index].first);

    rho = back_m200 / (FOURTHIRDSPI*back_r200*back_r200*back_r200);
    
    if (rho > rho_200) break;
  }
  
  // open the output file and begin appending data

  //string outfile = "profile.dat";
  
  stringstream outfile;
  outfile << "profile_" << setw(4) << setfill('0') << checkpoint << ".dat";

  ofstream dout(outfile.str().c_str());
  
  dout << "# z=" << redshift;
  dout << ", time=" << (siminfo.get_sim_time() - SN_START_TIME) * SEC_TO_YR / 1e6;
  dout << ", r200=" << r200 * scale * CM_TO_PC << ", m200=" << m200 / MSUN << ", numparts=" << pcount;
  //dout << ", pos=(" << max_x << ", " << max_y << ", " << max_z << ")" << endl;
  printf(", pos=(%20.13e,%20.13e,%20.13e)\n", max_x, max_y, max_z);

  cout << "# z=" << redshift;
  cout << ", time=" << (siminfo.get_sim_time() - SN_START_TIME) * SEC_TO_YR / 1e6;
  cout << ", r200=" << r200 * scale * CM_TO_PC << ", m200=" << m200 / MSUN << ", numparts=" << pcount;
  cout << ", pos=(" << max_x << ", " << max_y << ", " << max_z << ")" << endl;
  
  double m, r;
  m = parts[0].second;
  for (unsigned int index = 1; index < dm_count; index++)
  {
    m += parts[index].second;
    r = std::sqrt(parts[index].first);
    
    if (r > 10*r200) break;
    
    rho = m / (FOURTHIRDSPI*r*r*r);
    
    dout << r * scale * CM_TO_PC << " " << rho / scale3  << " " << m / MSUN << endl;
  }
  dout.close();

  return 0;
}
