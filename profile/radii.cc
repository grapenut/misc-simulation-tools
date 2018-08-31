#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <iterator>

#include "function.hh"

using namespace std;

const double M_sun = 2.0e33;
const double pc = 3.0e18;
const double Myr = 1.0e6 * 3.14e7;

inline double cub(double x) { return x * x * x; }

vector<double> enclosed_mass(const vector<double> & volume, const vector<double> & rho)
{
  const double safety = 1.0 + 1.0e-5;
  const double safety_tiny = 1.0e-100;

  vector<double> M(volume.size());

  M[0] = 0.0;

  for(size_t i = 1; i < volume.size(); i ++)
    {
      M[i] = M[i - 1] * safety + volume[i] * rho[i] + safety_tiny;
    }

  return M;
}


void make_plot_file(const char *infile, const char *outfile, size_t num_mass_points)
{
 ofstream plot("plot.gp");
  
  plot << "set yrange [0.1:1e4]" << endl
       << "set log x" << endl
       << "set log y" << endl
       << "set xlabel \"time\"" << endl
       << "set ylabel \"radius\"" << endl
       << "set term postscript landscape" << endl
       << "set out \"" << outfile << "\"" << endl
       << "set nokey" << endl
       << "plot ";

  for(size_t i = 0; i < num_mass_points; i ++)
    {
      if(i != 0) plot << ", ";

      plot << "\"" << infile << "\" using (($1-6.01568e15)/3.14e13):($" << i + 3 << "/3e18) with lines ";
    }
  
  plot << endl;
  
  plot.close();

  system("gnuplot plot.gp");


}




int main()
{


  ofstream out_t("t.out");
  ofstream out_b("b.out");
  ofstream out_m("m.out");

  ofstream star_out_t("star_t.out");
  ofstream star_out_b("star_b.out");
  ofstream star_out_m("star_m.out");
  
  ofstream part_out("part.out");
  ofstream star_part_out("star_part.out");


  out_t.setf(ios::scientific, ios::floatfield);
  out_b.setf(ios::scientific, ios::floatfield);
  out_m.setf(ios::scientific, ios::floatfield);
  out_t.precision(10);
  out_b.precision(10);
  out_m.precision(10);

  star_out_t.setf(ios::scientific, ios::floatfield);
  star_out_b.setf(ios::scientific, ios::floatfield);
  star_out_m.setf(ios::scientific, ios::floatfield);
  star_out_t.precision(10);
  star_out_b.precision(10);
  star_out_m.precision(10);


  part_out.setf(ios::scientific, ios::floatfield);
  part_out.precision(10);

  star_part_out.setf(ios::scientific, ios::floatfield);
  star_part_out.precision(10);

  ofstream metallicity_out("z.out");
  metallicity_out.setf(ios::scientific, ios::floatfield);
  metallicity_out.precision(10);

  ofstream mdot_out("mdot.out");
  mdot_out.setf(ios::scientific, ios::floatfield);
  mdot_out.precision(10);

  const double M_min = M_sun / 1024.0;
  const size_t num_mass_points = 40;

  const size_t index_begin = 60;
  const size_t index_count = 320;
  const size_t index_step = 5;

  vector<double> M(num_mass_points);
  //star_out_b << 0.0 << " " << 0.0;
  for(size_t i = 0; i < num_mass_points; i ++)
    {
      M[i] = M_min * pow(2.0, double(i));
  //    star_out_b << " " << M[i];
    }
//    star_out_b << endl;
  
  // vector<vector<double> > radii_b_all(M.size()), radii_m_all(M.size());




  const size_t coarse_number = 20;
  const double coarse_r_min = 0.1 * pc;
  const double coarse_r_max = 1000.0 * pc;

  vector<double> coarse_r(coarse_number);
  for(size_t i = 0; i < coarse_number; i ++) 
    coarse_r[i] = coarse_r_min * exp(double(i) * log(coarse_r_max / coarse_r_min) / double(coarse_number));

  for(size_t index = 0; index < index_count; index+=index_step)
    {

      ostrstream file_name;
      file_name << "profile_" << setfill('0') << setw(4) << index_begin + index << ".txt" << ends;
      ifstream in(file_name.str());

      if(!in) { cout << "file open fail " << file_name.str() << endl; abort(); }

      double t, z;

      vector<double> r(1, 0.0), volume(1, 0.0), rho_t(1, 0.0), rho_b(1, 0.0), rho_m(1, 0.0), v_b(1, 0.0), v_m(1, 0.0), star_volume(1, 0.0), star_rho_t(1, 0.0), star_rho_b(1, 0.0), star_rho_m(1, 0.0), star_v_b(1, 0.0), star_v_m(1, 0.0), part_volume(1, 0.0), part_rho(1, 0.0), part_v(1, 0.0), star_part_rho(1, 0.0), star_part_v(1, 0.0);

      double _r, _volume, _rho_t, _rho_b, _rho_m, _v_b, _v_m, _star_volume, _star_rho_t, _star_rho_b, _star_rho_m, _star_v_b, _star_v_m, _part_vol, _part_rho, _part_v, _star_part_rho, _star_part_v;

      in >> t >> z;

      while(in >> _r >> _volume >> _rho_t >> _rho_b >> _rho_m >> _v_b >> _v_m >> _star_volume >> _star_rho_t >> _star_rho_b >> _star_rho_m >> _star_v_b >> _star_v_m >> _part_vol >> _part_rho >> _part_v >> _star_part_rho >> _star_part_v)
	{
	  r.push_back(_r);
	  volume.push_back(_volume);
	  rho_t.push_back(_rho_t);
	  rho_b.push_back(_rho_b);
	  rho_m.push_back(_rho_m);
	  v_b.push_back(_v_b);
	  v_m.push_back(_v_m);
	  
	  star_volume.push_back(_star_volume);
	  star_rho_t.push_back(_star_rho_t);
	  star_rho_b.push_back(_star_rho_b);
	  star_rho_m.push_back(_star_rho_m);
	  star_v_b.push_back(_star_v_b);
	  star_v_m.push_back(_star_v_m);
	  
	  part_volume.push_back(_part_vol);

	  part_rho.push_back(_part_rho);
	  part_v.push_back(_part_v);

	  star_part_rho.push_back(_star_part_rho);
	  star_part_v.push_back(_star_part_v);
	  
	}

      const vector<double> M_t = enclosed_mass(volume, rho_t);
      const vector<double> M_b = enclosed_mass(volume, rho_b);
      const vector<double> M_m = enclosed_mass(volume, rho_m);

      const vector<double> star_M_t = enclosed_mass(star_volume, star_rho_t);
      const vector<double> star_M_b = enclosed_mass(star_volume, star_rho_b);
      const vector<double> star_M_m = enclosed_mass(star_volume, star_rho_m);

      const vector<double> part_M = enclosed_mass(part_volume, part_rho);
      const vector<double> star_part_M = enclosed_mass(part_volume, star_part_rho);

      Spline M_part_spline(r, star_part_M, 1.0 * M_sun, 1.0e9 * M_sun);
      Spline M_b_spline(r, star_M_b, 1.0 * M_sun, 1.0e9 * M_sun);
      
      mdot_out << t << " " << z << " " << M_b_spline(20.0 * pc) << " " << M_part_spline(20.0 * pc) << endl;


      vector<double> coarse_M_part(coarse_number), coarse_M_b(coarse_number);
      
      transform(coarse_r.begin(), coarse_r.end(), coarse_M_part.begin(), M_part_spline);
      transform(coarse_r.begin(), coarse_r.end(), coarse_M_b.begin(), M_b_spline);

      vector<double> r_metallicity(coarse_number - 1), z_metallicity(coarse_number - 1);

      for(size_t i = 0; i < (coarse_number - 1); i ++)
	{
	  r_metallicity[i] = sqrt(coarse_r[i] * coarse_r[i + 1]);
          
          if (index == 0) 
            metallicity_out << r_metallicity[i] << " ";

	  z_metallicity[i] = (coarse_M_part[i + 1] - coarse_M_part[i]) / (coarse_M_b[i + 1] - coarse_M_b[i]);
	}

      if (index == 0)
        metallicity_out << endl;
      
      for(size_t i = 1; i < r_metallicity.size(); i ++)
	{
	  metallicity_out << z_metallicity[i] << " ";
	}
      metallicity_out << endl;



      cout << M_m.back() / M_sun << " " << star_M_m.back() / M_sun << " " << part_M.back() / M_sun << " " << star_part_M.back() / M_sun  << endl;

      //  for(size_t i = 0; i < r.size(); i ++)
      //     {
      //       cout << r[i] / pc << " " << M_m[i] / M_sun << endl;
      //     }
      //   return 0;


      const double small_r = 1.0;
      const double large_r = 1.0e100;

      Spline r_of_M_t(M_t, r, small_r, large_r);
      Spline r_of_M_b(M_b, r, small_r, large_r);
      Spline r_of_M_m(M_m, r, small_r, large_r);

      Spline star_r_of_M_t(star_M_t, r, small_r, large_r);
      Spline star_r_of_M_b(star_M_b, r, small_r, large_r);
      Spline star_r_of_M_m(star_M_m, r, small_r, large_r);
  
      Spline part_r_of_M(part_M, r, small_r, large_r);
      Spline star_part_r_of_M(star_part_M, r, small_r, large_r);




      vector<double> radii_t(num_mass_points);
      vector<double> radii_b(num_mass_points);
      vector<double> radii_m(num_mass_points);

      vector<double> star_radii_t(num_mass_points);
      vector<double> star_radii_b(num_mass_points);
      vector<double> star_radii_m(num_mass_points);

      vector<double> part_radii(num_mass_points);
      vector<double> star_part_radii(num_mass_points);




      transform(M.begin(), M.end(), radii_t.begin(), r_of_M_t);
      transform(M.begin(), M.end(), radii_b.begin(), r_of_M_b);
      transform(M.begin(), M.end(), radii_m.begin(), r_of_M_m);

      transform(M.begin(), M.end(), star_radii_t.begin(), star_r_of_M_t);
      transform(M.begin(), M.end(), star_radii_b.begin(), star_r_of_M_b);
      transform(M.begin(), M.end(), star_radii_m.begin(), star_r_of_M_m);

      transform(M.begin(), M.end(), part_radii.begin(), part_r_of_M);
      transform(M.begin(), M.end(), star_part_radii.begin(), star_part_r_of_M);

   //    for(size_t i = 9; i < M.size(); i ++)
// 	{
// 	  radii_b_all[i].push_back(radii_b[i]);
// 	  radii_m_all[i].push_back(radii_m[i]);
// 	}
      
      //  for(size_t i = 0; i < num_mass_points; i ++)
      //     {
      //       cout << M[i] / M_sun << " " << radii_b[i] / pc << " " << radii_m[i] / pc << endl;
      //     }

      
      
      out_t << t << " " << z << " "; copy(radii_t.begin(), radii_t.end(), ostream_iterator<double> (out_t, " ")); out_t << endl;
      out_b << t << " " << z << " "; copy(radii_b.begin(), radii_b.end(), ostream_iterator<double> (out_b, " ")); out_b << endl;
      out_m << t << " " << z << " "; copy(radii_m.begin(), radii_m.end(), ostream_iterator<double> (out_m, " ")); out_m << endl;

      star_out_t << t << " " << z << " "; copy(star_radii_t.begin(), star_radii_t.end(), ostream_iterator<double> (star_out_t, " ")); star_out_t << endl;
      star_out_b << t << " " << z << " "; copy(star_radii_b.begin(), star_radii_b.end(), ostream_iterator<double> (star_out_b, " ")); star_out_b << endl;
      star_out_m << t << " " << z << " "; copy(star_radii_m.begin(), star_radii_m.end(), ostream_iterator<double> (star_out_m, " ")); star_out_m << endl;

      part_out << t << " " << z << " "; copy(part_radii.begin(), part_radii.end(), ostream_iterator<double> (part_out, " ")); part_out << endl;
      star_part_out << t << " " << z << " "; copy(star_part_radii.begin(), star_part_radii.end(), ostream_iterator<double> (star_part_out, " ")); star_part_out << endl;
      


    }



  out_t.close();
  out_b.close();
  out_m.close();

  star_out_t.close();
  star_out_b.close();
  star_out_m.close();
  
  part_out.close();
  star_part_out.close();

  metallicity_out.close();
  mdot_out.close();
  
  make_plot_file("t.out", "t.ps", num_mass_points);
  make_plot_file("b.out", "b.ps", num_mass_points);
  make_plot_file("m.out", "m.ps", num_mass_points);

  make_plot_file("star_t.out", "star_t.ps", num_mass_points);
  make_plot_file("star_b.out", "star_b.ps", num_mass_points);
  make_plot_file("star_m.out", "star_m.ps", num_mass_points);

  make_plot_file("part.out", "part.ps", num_mass_points);
  make_plot_file("star_part.out", "star_part.ps", num_mass_points);

 
  
  return 0;
}
