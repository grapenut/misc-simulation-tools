#include <vector>
#include <iostream> 
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>

#include "function.hh"

using namespace std;

const double c = 3.0e10;
const double Omega_m = 0.3;
const double Omega_b = 0.05;
const double Omega_Lambda = 0.7;
const double h = 0.7;
const double km = 1.0e5;
const double pc = 3.0e18;
const double Mpc = 1.0e6 * pc;
const double H_0 = h * 100 * km / Mpc;
const double G = 6.67e-8;
const double M_sun = 2.0e33;

const double f_df = 1.0;  // fixme
const double ln_Lambda = 2.0;

const double tolerance_exact = 1.0e-10;
const double tolerance_ode = 1.0e-6;
const double stepsize_ode = 1.0e-6;

const double delta_c_factor = 3.0 * pow(12.0 * M_PI, 2.0 / 3.0) / 20.0;
const double Delta_vir_factor = 18.0 * M_PI * M_PI;

const string sigma_file_name("sigma_dsigma.mass");

inline double sq(double x) { return x * x; }

inline double H(double z) 
{
  const double one_plus_z = 1.0 + z;

  return H_0 * sqrt(Omega_Lambda + one_plus_z * one_plus_z * one_plus_z * Omega_m);
}

double distanceIntegrand(double z, void * param)
{
  return c / H(z);
}

double angularDiameterDistance(double z)
{
  static const size_t work_size = 1000;
  static gsl_integration_workspace * w = gsl_integration_workspace_alloc(work_size);
  
  gsl_function f;
  f.function = &distanceIntegrand;
  f.params = 0;

  double integral;

  const double eps_abs = 0.0;
  const double eps_rel = tolerance_exact;
  double abserr;
  gsl_integration_qag(&f, 0, z, eps_abs, eps_rel, work_size, 6, w, &integral, &abserr);
  const double one_plus_z = 1.0 + z;

  return integral / one_plus_z;
}

double GrowthFactorIntegrand(double z, void * param)
{
  const double one_plus_z = 1.0 + z;
  
  return one_plus_z * pow(Omega_Lambda + one_plus_z * one_plus_z * one_plus_z * Omega_m, -1.5);
}

double GrowthNormalization()
{
  static const size_t work_size = 1000;
  static gsl_integration_workspace * w = gsl_integration_workspace_alloc(work_size);
  
  gsl_function f;
  f.function = &GrowthFactorIntegrand;
  f.params = 0;

  double integral;

  const double eps_abs = 0.0;
  const double eps_rel = tolerance_exact;
  double abserr;
  gsl_integration_qagiu(&f, 0.0, eps_abs, eps_rel, work_size, w, &integral, &abserr);

  return 1.0 / integral;
}

double GrowthFactor(double z)
{
  static const size_t work_size = 1000;
  static gsl_integration_workspace * w = gsl_integration_workspace_alloc(work_size);
  static const double normalization = GrowthNormalization();

  gsl_function f;
  f.function = &GrowthFactorIntegrand;
  f.params = 0;

  double integral;

  const double eps_abs = 0.0;
  const double eps_rel = tolerance_exact;
  double abserr;
  gsl_integration_qagiu(&f, z, eps_abs, eps_rel, work_size, w, &integral, &abserr);
  const double one_plus_z = 1.0 + z;

  return normalization * sqrt(Omega_Lambda + one_plus_z * one_plus_z * one_plus_z * Omega_m) * integral;
}

double GrowthFactorDeriv(double z)
{
  static const size_t work_size = 1000;
  static gsl_integration_workspace * w = gsl_integration_workspace_alloc(work_size);
  static const double normalization = GrowthNormalization();

  gsl_function f;
  f.function = &GrowthFactorIntegrand;
  f.params = 0;

  double integral;

  const double eps_abs = 0.0;
  const double eps_rel = tolerance_exact;
  double abserr;
  gsl_integration_qagiu(&f, z, eps_abs, eps_rel, work_size, w, &integral, &abserr);
  const double one_plus_z = 1.0 + z;

  return 
    normalization * (0.5 / sqrt(Omega_Lambda + one_plus_z * one_plus_z * one_plus_z * Omega_m)) * 
    (Omega_m * 3.0 * one_plus_z * one_plus_z) * integral -
    normalization * one_plus_z / (Omega_Lambda + one_plus_z * one_plus_z * one_plus_z * Omega_m);
}


double delta_c(double z) 
{
  return delta_c_factor / GrowthFactor(z);
}

double delta_c_deriv(double z)
{
  const double f = GrowthFactor(z);
  const double df = GrowthFactorDeriv(z);
  return - delta_c_factor * df / (f * f);
}


int main(int argc, char *argv[])
{
  if(argc < 4)
    {
      cout << "usage: " << argv[0] << " <Mpc> <mass> <redshift>" << endl;
      return 0;
    }

  double L = atof(argv[1]);
  double M = atof(argv[2]);
  double z = atof(argv[3]);
  double sigma_coeff;
  
  if (argc < 5)
    sigma_coeff = 1.0;
  else
    sigma_coeff = atof(argv[4]);

  cout << "size = " << L << endl;
  cout << "mass = " << M << endl;
  cout << "   z = " << z << endl;
  cout << "sigm = " << sigma_coeff << endl;

  ifstream in("sigma_m.dat");
  if(!in)
    {
      cout << "cannot open file " << endl;
      abort();
    }
  double m, s, ds;
  vector<double> _mass, _sigma;
  while(in >> m >> s)
    {
      _mass.push_back(m);
      _sigma.push_back(sigma_coeff * s);
    }

  Function sigma(_mass, _sigma);

  static const double factor = 3.0 * H_0 * H_0 / (8.0 * M_PI * G);
  
  const double one_plus_z = 1.0 + z;
  double scale = 1.0 / (1.0 + z);
        
  const double density = Omega_m * factor * one_plus_z * one_plus_z * one_plus_z;  

  vector<double> _dn;
  double val;
  
  double coeff = sqrt(2.0 / M_PI);
  
  double dc;
  
  size_t count = _mass.size();
  
  double m_max = 1.0e13;
  double m_min = 1.0e2;
  
  double volume = L*Mpc * L*Mpc * L*Mpc * scale * scale * scale;
  //volume = 4.0 * M_PI / 3.0 * (L*Mpc*scale) * (L*Mpc*scale) * (L*Mpc*scale);

  double tol = 1.0e-2;
  
  double z_start = 40.0;
  double z_end = 10.0;
  double dz = -0.5;
  
  double sig_start = 1.0;
  double sig_end = 2.0;
  double sig_delta = 0.1;

  double m_bot, m_top;


  dc = delta_c(z);
  s = sigma(M);
  
  cout << "Extremness = " << dc / s << endl;


  Function *dn;
  ofstream fout("most_massive.dat");

  
  for (z = z_start; z > z_end; z += dz)
  {
  cout.precision(10);
  cout.setf(ios::fixed,ios::scientific);
    cout << z;
    fout << z;

    dc = delta_c(z);
    
    
    for (sigma_coeff = sig_start; sigma_coeff <= sig_end + sig_delta/2.0; sigma_coeff += sig_delta)
    {
      _dn.clear();
      for (size_t i = 0; i < _mass.size(); i++)
      {
        m = _mass[i];
        
        s = sigma_coeff * sigma(m);
        ds = sigma_coeff * fabs(sigma.deriv(m));
      
        val = coeff * density / m / M_sun * dc / s / s * ds * exp(- dc * dc / 2.0 / s / s);
        
        _dn.push_back(val);
      }
      dn = new Function(_mass, _dn);
      

      m_top = m_max;
      m_bot = m_min;
      m = m_top / 2.0;
      while (1)
      {
      
        if (m == m_max)
        {
          cout << "m = m_max! " << m << endl;
          exit(0);
        }
      
        val = volume * dn->integ(m, m_max);
        
        if ((fabs(1.0 - val) < tol))
        {
cout.width(10);
  cout.precision(5);
  cout.setf(ios::fixed,ios::scientific);
          cout << " " << m;
          fout << " " << m;
          break;
        }
        
        if (val > 1.0)
        {
          m_bot = m;
          m = 0.5*(m + m_top);
        } else {
          m_top = m;
          m = 0.5*(m + m_bot);
        }
        
        //cout << m << " " << val << endl;
      }
    }
    
    cout << endl;
    fout << endl;
  }
  
  //double normal = 1.0 / val;

  //ofstream out("dn_dm.dat");
  //for (size_t i = 0; i < count; i++)
  //{
  //  m = m_min * exp(double(i) * log(m_max / m_min) / double(count - 1));
  // 
  //  out << m << " " << sigma(m) << " " << dn(m) << endl;
  //}
  
  //cout << "Total mass is " << volume * density << " = " << volume * density / M_sun << endl;
  
  //cout << "The value is " << volume * val << endl;

  return 0;
}
