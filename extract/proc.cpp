

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#define PI	3.14159
#define MSUN	1.99e33
#define PC	3.0856776e18
#define H0	70.0
#define G	6.67e-8

using namespace std;

typedef struct {
  double px, py, pz, vx, vy, vz, mass, dens, len, temp;
} particle;

int main(int argc, void *argv)
{
  // location of the densest particle from gadget extraction
  // could find this automatically, but why not just write it in for now
  // 9.732840049e+23 7.518643153e+23 1.552571738e+23
  // 9.73276063088405e+23 7.51858330296794e+23 1.5527222517297e+23
  // 9.73272685053183e+23 7.51855542928802e+23 1.55278391978447e+23
  // 9.73274e+23 7.51856e+23 1.55289e+23
  
  // here i've iterated before using d = 100 pc
  // 9.73267558078785e+23 7.51856417433053e+23 1.55281636997704e+23
  double dcx = 9.73267558078785e23;
  double dcy = 7.51856417433053e23;
  double dcz = 1.55281636997704e23;
  
  double cx, cy, cz;
  
  double d = 20.0 * PC;
  int ntimes = 10;
  
  double d2 = d*d;
  double r2, oldr2;
  double dx, dy, dz, dr;
  
  double dens_max;
  double x_max;
  double y_max;
  double z_max;
  
  int i, n;
  
  // 1 Mpc box @ z=10 (physical)
  double width = 1000000.0 * PC / 11.0;
  double halfwidth = 0.5 * width;
  
  double scale = 1.0 / 11.0;
  double velocity_fix = 1.0;
  
  ifstream gasfile("baryon.txt");
  ifstream dmfile("darkmatter.txt");
  
  const char *gas_output_file = "gas_cutout.txt";
  const char *dm_output_file = "dm_cutout.txt";

  vector<particle> gasdata, dmdata;

  double px, py, pz, vx, vy, vz, mass, dens, len, temp;
  particle p;
  
  double total_gmass, total_gmvx, total_gmvy, total_gmvz;
  int total_count;

  double total_mass, total_mvx, total_mvy, total_mvz;
  int total_gcount;

  FILE *fp;
  
  
  dens_max = 0.0;
  x_max = 0.0;
  y_max = 0.0;
  z_max = 0.0;
  
  // initialize center of momentum aggregators
  total_gmass = 0.0;
  total_gmvx = 0.0;
  total_gmvy = 0.0;
  total_gmvz = 0.0;
  total_gcount = 0;

  total_mass = 0.0;
  total_mvx = 0.0;
  total_mvy = 0.0;
  total_mvz = 0.0;
  total_count = 0;
  
  // load gas particles
  while (gasfile >> px >> py >> pz >> vx >> vy >> vz >> mass >> dens >> len >> temp)
  {
    p.px = px;
    p.py = py;
    p.pz = pz;

    p.vx = vx * velocity_fix;
    p.vy = vy * velocity_fix;
    p.vz = vz * velocity_fix;
    
    p.dens = dens;
    p.len = len;
    p.temp = temp;
    
    p.mass = mass;
    
    dx = px - dcx;
    dy = py - dcy;
    dz = pz - dcz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    if (dens > dens_max) {
      dens_max = dens;
      x_max = p.px;
      y_max = p.py;
      z_max = p.pz;
    }
    
    if (r2 < d2) {
      total_gmvx += mass * p.px;
      total_gmvy += mass * p.py;
      total_gmvz += mass * p.pz;
      total_gmass += mass;
      total_gcount++;
    }
    
    gasdata.push_back(p);
  }

  ////////////////////////////////////////////////////////////////////////////
  
  // max density particle
  dx = x_max - dcx;
  dy = y_max - dcy;
  dz = z_max - dcz;
  r2 = dx*dx + dy*dy + dz*dz;
  
  cout << "Density:" << std::setprecision(15) << x_max << " " 
       << std::setprecision(15) << y_max << " " 
       << std::setprecision(15) << z_max << " " << sqrt(r2) / PC << endl;

  // center of momentum velocity vector
  total_gmvx = total_gmvx / total_gmass;
  total_gmvy = total_gmvy / total_gmass;
  total_gmvz = total_gmvz / total_gmass;
  
  dx = total_gmvx - dcx;
  dy = total_gmvy - dcy;
  dz = total_gmvz - dcz;
  r2 = dx*dx + dy*dy + dz*dz;
  
  cout << "GasCoM: " << total_gcount << " " << total_gmass/MSUN << " " 
       << std::setprecision(15) << total_gmvx << " " 
       << std::setprecision(15) << total_gmvy << " "
       << std::setprecision(15) << total_gmvz << " " << sqrt(r2) / PC << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  // load dark matter particles
  while (dmfile >> px >> py >> pz >> vx >> vy >> vz >> mass)
  {
    p.px = px;
    p.py = py;
    p.pz = pz;
    
    p.vx = vx * velocity_fix;
    p.vy = vy * velocity_fix;
    p.vz = vz * velocity_fix;
    
    p.dens = 0.0;
    p.len = 0.0;
    p.temp = 0.0;
    
    p.mass = mass;
    
    dx = px - dcx;
    dy = py - dcy;
    dz = pz - dcz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 < d2) {
      total_mvx += mass * p.px;
      total_mvy += mass * p.py;
      total_mvz += mass * p.pz;
      total_mass += mass;
      total_count++;
    }
    
    dmdata.push_back(p);
  }
  
  ////////////////////////////////////////////////////////////////////////////
  
  // center of momentum velocity vector
  total_mvx = total_mvx / total_mass;
  total_mvy = total_mvy / total_mass;
  total_mvz = total_mvz / total_mass;
  
  dx = total_mvx-dcx;
  dy = total_mvy-dcy;
  dz = total_mvz-dcz;
  r2 = dx*dx + dy*dy + dz*dz;
  
  cout << "DM CoM: " << total_count << " " << total_mass/MSUN << " "
       << std::setprecision(15) << total_mvx << " "
       << std::setprecision(15) << total_mvy << " "
       << std::setprecision(15) << total_mvz << " " << sqrt(r2) / PC << endl;
  
  // center of momentum velocity vector
  cx = (total_mvx*total_mass + total_gmvx*total_gmass) / (total_mass + total_gmass);
  cy = (total_mvy*total_mass + total_gmvy*total_gmass) / (total_mass + total_gmass);
  cz = (total_mvz*total_mass + total_gmvz*total_gmass) / (total_mass + total_gmass);
  
  dx = cx-dcx;
  dy = cy-dcy;
  dz = cz-dcz;
  r2 = dx*dx + dy*dy + dz*dz;
  
  cout << "TotCoM: " << total_count+total_gcount << " " << (total_mass+total_gmass)/MSUN << " "
       << std::setprecision(15) << cx << " "
       << std::setprecision(15) << cy << " "
       << std::setprecision(15) << cz << " " << sqrt(r2) / PC << endl;

  
  cout << "Loaded " << gasdata.size() << " gas particles." << endl;
  cout << "Loaded " << dmdata.size() << " dark matter particles." << endl;

  ////////////////////////////////////////////////////////////////////////////
  // iterate on the center of mass until it converges
  oldr2 = r2;
  for (n = 0; n < ntimes; ++n) {
    cout << "-- " << n << " -----------------------------------------------------------" << endl;
    total_gmass = 0.0;
    total_gmvx = 0.0;
    total_gmvy = 0.0;
    total_gmvz = 0.0;
    total_gcount = 0;
  
    total_mass = 0.0;
    total_mvx = 0.0;
    total_mvy = 0.0;
    total_mvz = 0.0;
    total_count = 0;
  
    for (i = 0; i < gasdata.size(); ++i) {
      p = gasdata[i];
      dx = p.px - cx;
      dy = p.py - cy;
      dz = p.pz - cz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < d2) {
        total_gmass += p.mass;
        total_gmvx += p.mass * p.px;
        total_gmvy += p.mass * p.py;
        total_gmvz += p.mass * p.pz;
        total_gcount++;
      }
    }
    
    for (i = 0; i < dmdata.size(); ++i) {
      p = dmdata[i];
      dx = p.px - cx;
      dy = p.py - cy;
      dz = p.pz - cz;
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < d2) {
        total_mass += p.mass;
        total_mvx += p.mass * p.px;
        total_mvy += p.mass * p.py;
        total_mvz += p.mass * p.pz;
        total_count++;
      }
    }

    total_gmvx = total_gmvx / total_gmass;
    total_gmvy = total_gmvy / total_gmass;
    total_gmvz = total_gmvz / total_gmass;

    dx = total_gmvx - cx;
    dy = total_gmvy - cy;
    dz = total_gmvz - cz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    cout << "GasCoM: " << total_gcount << " " << total_gmass/MSUN << " " 
         << std::setprecision(15) << total_gmvx << " " 
         << std::setprecision(15) << total_gmvy << " "
         << std::setprecision(15) << total_gmvz << " " << sqrt(r2) / PC << endl;
  
    total_mvx = total_mvx / total_mass;
    total_mvy = total_mvy / total_mass;
    total_mvz = total_mvz / total_mass;
    
    dx = total_mvx - cx;
    dy = total_mvy - cy;
    dz = total_mvz - cz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    cout << "DM CoM: " << total_count << " " << total_mass/MSUN << " "
         << std::setprecision(15) << total_mvx << " "
         << std::setprecision(15) << total_mvy << " "
         << std::setprecision(15) << total_mvz << " " << sqrt(r2) / PC << endl;

    px = (total_mvx*total_mass + total_gmvx*total_gmass) / (total_mass + total_gmass);
    py = (total_mvy*total_mass + total_gmvy*total_gmass) / (total_mass + total_gmass);
    pz = (total_mvz*total_mass + total_gmvz*total_gmass) / (total_mass + total_gmass);
    
    dx = px - cx;
    dy = py - cy;
    dz = pz - cz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    cout << "TotCoM: " << total_count+total_gcount << " " << (total_mass+total_gmass)/MSUN << " "
         << std::setprecision(15) << px << " "
         << std::setprecision(15) << py << " "
         << std::setprecision(15) << pz << " " << sqrt(r2) / PC << endl;
    
    cx = px;
    cy = py;
    cz = pz;
    
    dr = sqrt(r2) - sqrt(oldr2);
    cout << "Relative change: " << dr / PC << endl;
    
    if (dr == 0.0) {
      cout << "Stopping iteration." << endl;
      break;
    }
    
    dr = fabs(dr / sqrt(oldr2));
    if (dr < 0.01) {
      cout << "Stopping iteration." << endl;
      break;
    }

    oldr2 = r2;
  }


  ////////////////////////////////////////////////////////////////////////////
  // find bulk velocity of the center of mass sphere
  
  total_mvx = 0.0;
  total_mvy = 0.0;
  total_mvz = 0.0;
  total_mass = 0.0;
  total_count = 0;
  
  for (i = 0; i < gasdata.size(); ++i)
  {
    p = gasdata[i];
    dx = p.px - cx;
    dy = p.py - cy;
    dz = p.pz - cz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 < d2) {
      total_count++;
      total_mass += p.mass;
      total_mvx += p.mass * p.vx;
      total_mvy += p.mass * p.vy;
      total_mvz += p.mass * p.vz;
    }
  }

  for (i = 0; i < dmdata.size(); ++i)
  {
    p = dmdata[i];
    dx = p.px - cx;
    dy = p.py - cy;
    dz = p.pz - cz;
    r2 = dx*dx + dy*dy + dz*dz;
    
    if (r2 < d2) {
      total_count++;
      total_mass += p.mass;
      total_mvx += p.mass * p.vx;
      total_mvy += p.mass * p.vy;
      total_mvz += p.mass * p.vz;
    }
  }
  
  total_mvx = total_mvx / total_mass;
  total_mvy = total_mvy / total_mass;
  total_mvz = total_mvz / total_mass;
  
  r2 = total_mvx*total_mvx + total_mvy*total_mvy + total_mvz*total_mvz;
  
  cout << "Velocity shifts: " << total_mvx / 1.0e5 << " " << total_mvy / 1.0e5 << " " << total_mvz / 1.0e5 << " = " << sqrt(r2) << endl;
  
  
  ////////////////////////////////////////////////////////////////////////////

  fp = fopen(gas_output_file, "w");
  if (!fp)
  {
    printf("ERROR unable to open gas output file\n");
    exit(0);
  }
  
  total_gcount = 0;
  for (i = 0; i < gasdata.size(); ++i)
  {
    p = gasdata[i];
    
    p.px = p.px - cx + halfwidth;
    p.py = p.py - cy + halfwidth;
    p.pz = p.pz - cz + halfwidth;

    if (p.px >= width || p.px <= 0.0) continue;
    if (p.py >= width || p.py <= 0.0) continue;
    if (p.pz >= width || p.pz <= 0.0) continue;
    
    total_gcount++;
    
    p.vx -= total_mvx;
    p.vy -= total_mvy;
    p.vz -= total_mvz;

    fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", p.px, p.py, p.pz, p.vx, p.vy, p.vz, p.mass, p.dens, p.len, p.temp);
  }
  
  fclose(fp);
  
  fp = fopen(dm_output_file, "w");
  if (!fp)
  {
    printf("ERROR unable to open dark matter output file\n");
    exit(0);
  }
  
  total_count = 0;
  for (i = 0; i < dmdata.size(); ++i)
  {
    p = dmdata[i];
    
    p.px = p.px - cx + halfwidth;
    p.py = p.py - cy + halfwidth;
    p.pz = p.pz - cz + halfwidth;
  
    if (p.px >= width || p.px <= 0.0) continue;
    if (p.py >= width || p.py <= 0.0) continue;
    if (p.pz >= width || p.pz <= 0.0) continue;
    
    total_count++;
    
    p.vx -= total_mvx;
    p.vy -= total_mvy;
    p.vz -= total_mvz;

    fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", p.px, p.py, p.pz, p.vx, p.vy, p.vz, p.mass);
  }
  
  fclose(fp);
  
  cout << "Wrote " << total_gcount << " gas particles." << endl;
  cout << "Wrote " << total_count << " dark matter particles." << endl;

  return 0;
}

