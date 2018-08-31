/*
 * flash2gadget
 *   written by Jeremy Ritter (2011)
 *
 * Usage: flash2gadget inputfile outputfile
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <math.h>
#include <string.h>

#define STRING_LEN 80
#define SMALL_STRING_LEN 24

#define CM_TO_MPC 3.24077885e-25
#define CM_TO_KM 1.0e-5
#define SEC_TO_GYR 3.16887646e-17
#define H0_TO_h 3.08568025e17
#define GRAM_TO_OMEGA (1.0 / 1.989e43)

typedef struct
{
  double scalefactor;
  double boxsize;
  int npart;
  double Omega0;
  double OmegaBaryon;
  double OmegaLambda;
  double redshift;
  double time;
  double hubbleparam;
  int map_x;
  int map_y;
  int map_z;
  int map_mass;
  int map_vx;
  int map_vy;
  int map_vz;
} hdata;

typedef struct
{
  char name[STRING_LEN];
  double value;
} real_list_t;
    
void read_data(const char *, double ***, hdata *);
double read_scalar(hid_t, const char *, const char *);
int read_particle_map(hid_t, const char *);

void write_data(const char *, double **, hdata *);

void write_attr_array_d(hid_t, const char *, double *, hsize_t);
void write_attr_array_i(hid_t, const char *, int *, hsize_t);
void write_attr_array_ui(hid_t, const char *, unsigned int *, hsize_t);

void write_attr_i(hid_t, const char *, int);
void write_attr_d(hid_t, const char *, double);
void write_attr_f(hid_t, const char *, float);
void write_attr_s(hid_t, const char *, const char *, hsize_t);
void write_header(hid_t, hdata *);
void write_constants(hid_t);
void write_units(hid_t);
void write_parameters(hid_t);

int main(int argc, char **argv)
{

  double **flash_data;
  int npart;
  char infile_name[200], outfile_name[200];
  int i, j;
  hdata header;
  
  if (argc != 3)
  {
    printf("Usage: %s inputfile outputfile\n", argv[0]);
    exit(1);
  }
  
  sprintf(infile_name, argv[1]);
  sprintf(outfile_name, argv[2]);

  read_data(infile_name, &flash_data, &header);
  write_data(outfile_name, flash_data, &header);

  free(flash_data[0]);
  free(flash_data);

  printf("Done.\n");
  return 0;
}

///////////////////////////////////////////////////////////////////////////

void read_data(const char *filename, double ***indata, hdata *header)
{
  hid_t file, space, dset;
  herr_t status;
  hsize_t dims[2];
  int ndims, sizex, sizey, i, j;
  double **data;
  
  printf("Opening input file: %s\n", filename);
  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
  {
    printf("Could not open input file!\n");
    exit(1);
  }
  
  printf("  Reading particle dataset.\n");
  dset = H5Dopen(file, "tracer particles");
  if (dset < 0)
  {
    printf("Could not open input dataset!\n");
    exit(1);
  }
  
  space = H5Dget_space(dset);
  ndims = H5Sget_simple_extent_dims(space, dims, NULL);
  sizex = dims[0];
  sizey = dims[1];

  data = (double **) malloc(sizex * sizeof(double *));
  data[0] = (double *) malloc(sizex * sizey * sizeof(double));
  for (i=1; i<sizex; i++)
    data[i] = data[0] + i * sizey;

  status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data[0]);
  status = H5Sclose(space);
  status = H5Dclose(dset);

  printf("  Reading header data.\n");
  header->npart = sizex;
  header->scalefactor = read_scalar(file, "real scalars", "scalefactor");
  header->redshift = (1.0 / header->scalefactor) - 1.0;
  header->time = read_scalar(file, "real scalars", "time") * SEC_TO_GYR;

  header->Omega0 = read_scalar(file, "real runtime parameters", "omegamatter");
  header->OmegaBaryon = read_scalar(file, "real runtime parameters", "omegabaryon");
  header->OmegaLambda = read_scalar(file, "real runtime parameters", "cosmologicalconstant");
  header->boxsize = read_scalar(file, "real runtime parameters", "xmax") * CM_TO_MPC;
  header->hubbleparam = read_scalar(file, "real runtime parameters", "hubbleconstant") * H0_TO_h;

  printf("  Reading particle map.\n");
  header->map_x = read_particle_map(file, "posx");
  header->map_y = read_particle_map(file, "posy");
  header->map_z = read_particle_map(file, "posz");
  header->map_mass = read_particle_map(file, "mass");
  header->map_vx = read_particle_map(file, "velx");
  header->map_vy = read_particle_map(file, "vely");
  header->map_vz = read_particle_map(file, "velz");

  status = H5Fclose(file);

  *indata = data;

  return;
}

////////////////////////////////////////////////////////////////////////////////////

double read_scalar(hid_t file, const char *dataset, const char *which)
{
  hid_t dset, space, memspace;
  hid_t string_type, real_list_type;
  hsize_t dims;
  herr_t status;
  int i;
  int found = 0;
  double retval = 0.0;
  real_list_t *real_list;
  
  char nameval[STRING_LEN];
  
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, STRING_LEN);
  
  dset = H5Dopen(file, dataset);
  if (dset < 0)
  {
    printf("Could not open real parameter list: %s\n", dataset);
    exit(1);
  }
  
  space = H5Dget_space(dset);
  H5Sget_simple_extent_dims(space, &dims, NULL);

  real_list = (real_list_t *) malloc(dims * sizeof(real_list_t));

  memspace = H5Screate_simple(1, &dims, NULL);

  real_list_type = H5Tcreate(H5T_COMPOUND, sizeof(real_list_t));

  H5Tinsert(real_list_type, 
          "name", 
          HOFFSET(real_list_t, name),
          string_type);

  H5Tinsert(real_list_type, 
          "value", 
          HOFFSET(real_list_t, value),
          H5T_IEEE_F64LE);
  
  status = H5Dread(dset, real_list_type, memspace, space, H5P_DEFAULT, real_list);
  if (status < 0) {
    printf("Could not read real parameters list: %s\n", dataset);
    exit(1);
  }

  for (i = 0; i < dims; i++) {
    strncpy(nameval, real_list[i].name, STRING_LEN);
    if (!strncasecmp(nameval, which, strlen(which)))
    {
      found = 1;
      retval = real_list[i].value;
    }
  }
  
  if (found == 0)
  {
    printf("Could not find real parameter: %s\n", which);
    exit(1);
  }

  free(real_list);
  H5Tclose(string_type);
  H5Tclose(real_list_type);
  H5Sclose(memspace);
  H5Sclose(space);
  H5Dclose(dset);

  return retval;
  
}

////////////////////////////////////////////////////////////////////////////////////

int read_particle_map(hid_t file, const char *which)
{
  hid_t dset, space, memspace;
  hid_t string_type;
  hsize_t dims[2];
  herr_t status;
  int i;
  int found = 0;
  int retval = 0.0;
  char **string_list;
  
  char nameval[SMALL_STRING_LEN];
  
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, SMALL_STRING_LEN);
  
  dset = H5Dopen(file, "particle names");
  if (dset < 0)
  {
    printf("Could not open particle names!\n");
    exit(1);
  }
  
  space = H5Dget_space(dset);
  H5Sget_simple_extent_dims(space, dims, NULL);

  string_list = (char **) malloc(dims[0] * sizeof(char *));
  string_list[0] = (char *) malloc(dims[0] * SMALL_STRING_LEN); 
  for(i=1;i<dims[0];i++)
    string_list[i] = string_list[0] + i * SMALL_STRING_LEN;

  memspace = H5Screate_simple(2, dims, NULL);

  status = H5Dread(dset, string_type, memspace, space, H5P_DEFAULT, string_list[0]);
  if (status < 0) {
    printf("Could not read particle names!\n");
    exit(1);
  }

  for (i = 0; i < dims[0]; i++) {
    strncpy(nameval, string_list[i], SMALL_STRING_LEN);
    if (!strncasecmp(nameval, which, strlen(which)))
    {
      printf("    Found mapping variable(%d): %s - %s\n", i, which, nameval);
      found = 1;
      retval = i;
    }
  }
  
  if (found == 0)
  {
    printf("Could not find particle name: %s\n", which);
    exit(1);
  }

  free(string_list[0]);
  free(string_list);
  H5Tclose(string_type);
  H5Sclose(memspace);
  H5Sclose(space);
  H5Dclose(dset);

  return retval;
  
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

void write_data(const char *filename, double **data, hdata *header)
{
  hid_t file, space, group, dset, aspace, attr;
  hsize_t dims[2];
  herr_t status;
  double attrdata;
  int i, j, npart;
  float **coords, *mass_data, **vel_data;
  unsigned int *ids;
  double Omega0, OmegaB, del_omega, h, a, sqrt_a;

  npart = header->npart;
  Omega0 = header->Omega0;
  OmegaB = header->OmegaBaryon;
  del_omega = Omega0 / (Omega0 - OmegaB);
  h = header->hubbleparam;
  a = header->scalefactor;
  sqrt_a = sqrt(a);
  
  // Build arrays to hold gadget data
  coords = (float **) malloc(npart * sizeof(float *));
  coords[0] = (float *) malloc(npart * 3 * sizeof(float));
  for (i=0; i<npart; i++)
    coords[i] = coords[0] + i * 3;
  
  mass_data = (float *) malloc(npart * sizeof(float));
  ids = (unsigned int *) malloc(npart * sizeof(unsigned int));

  vel_data = (float **) malloc(npart * sizeof(float *));
  vel_data[0] = (float *) malloc(npart * 3 * sizeof(float));
  for (i=0; i<npart; i++)
    vel_data[i] = vel_data[0] + i * 3;
    
  for (i=0; i < npart; i++)
  {
    coords[i][0] = (float)( data[i][header->map_x] * CM_TO_MPC * h );
    coords[i][1] = (float)( data[i][header->map_y] * CM_TO_MPC * h );
    coords[i][2] = (float)( data[i][header->map_z] * CM_TO_MPC * h );

    vel_data[i][0] = (float)( data[i][header->map_vx] * CM_TO_KM * sqrt_a );
    vel_data[i][1] = (float)( data[i][header->map_vy] * CM_TO_KM * sqrt_a );
    vel_data[i][2] = (float)( data[i][header->map_vz] * CM_TO_KM * sqrt_a );

    mass_data[i] = (float)( data[i][header->map_mass] * GRAM_TO_OMEGA * del_omega * h );
     
    ids[i] = i+1;
  }
  
  printf("Opening output file: %s\n", filename);
  file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0)
  {
    printf("Could not open output file!\n");
    exit(1);
  }


  // Write dark matter particles
  printf("  Creating output group PartType1.\n");
  group = H5Gcreate(file, "PartType1", 0);
  if (group < 0)
  {
    printf("Could not open PartType1 group!\n");
    exit(1);
  }

  // Write out Coordinates dataset
  printf("  Writing Coordinates dataset.\n");
  dims[0] = npart;
  dims[1] = 3;
  space = H5Screate_simple(2, dims, NULL);
  dset = H5Dcreate(group, "Coordinates", H5T_IEEE_F32LE, space, H5P_DEFAULT);
  status = H5Dwrite(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *coords);
  status = H5Sclose(space);
  write_attr_d(dset, "CGSConversionFactor", 1.0 / CM_TO_MPC);
  write_attr_f(dset, "aexp-scale-exponent", 1.0f);
  write_attr_f(dset, "h-scale-exponent", -1.0f);
  status = H5Dclose(dset);

  // Write out Velocity dataset
  printf("  Writing Velocity dataset.\n");
  dims[0] = npart;
  dims[1] = 3;
  space = H5Screate_simple(2, dims, NULL);
  dset = H5Dcreate(group, "Velocity", H5T_IEEE_F32LE, space, H5P_DEFAULT);
  status = H5Dwrite(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *vel_data);
  status = H5Sclose(space);
  write_attr_d(dset, "CGSConversionFactor", 1.0 / CM_TO_KM);
  write_attr_f(dset, "aexp-scale-exponent", 0.5f);
  write_attr_f(dset, "h-scale-exponent", 0.0f);
  status = H5Dclose(dset);

  // Write out Mass dataset
  printf("  Writing Mass dataset.\n");
  dims[0] = npart;
  dims[1] = 1;
  space = H5Screate_simple(1, dims, NULL);
  dset = H5Dcreate(group, "Mass", H5T_IEEE_F32LE, space, H5P_DEFAULT);
  status = H5Dwrite(dset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass_data);
  status = H5Sclose(space);
  write_attr_d(dset, "CGSConversionFactor", 1.0 / GRAM_TO_OMEGA);
  write_attr_f(dset, "aexp-scale-exponent", 0.0f);
  write_attr_f(dset, "h-scale-exponent", -1.0f);
  status = H5Dclose(dset);

  // Write out IDs dataset
  printf("  Writing ParticleIDs dataset.\n");
  dims[0] = npart;
  dims[1] = 1;
  space = H5Screate_simple(1, dims, NULL);
  dset = H5Dcreate(group, "ParticleIDs", H5T_STD_U32LE, space, H5P_DEFAULT);
  status = H5Dwrite(dset, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
  status = H5Sclose(space);
  write_attr_d(dset, "CGSConversionFactor", 1.0);
  write_attr_f(dset, "aexp-scale-exponent", 0.0f);
  write_attr_f(dset, "h-scale-exponent", 0.0f);
  status = H5Dclose(dset);

  status = H5Gclose(group);

  // Write out the rest of the metadata
  printf("  Writing header data.\n");
  write_header(file, header);
  write_constants(file);
  write_units(file);
  write_parameters(file);

  status = H5Fclose(file);
  
  free(coords[0]);
  free(coords);
  free(mass_data);
  free(ids);
  free(vel_data[0]);
  free(vel_data);
  
  return;
}

//////////////////////////////////////////////////////////////////////////

void write_header(hid_t file, hdata *header)
{
  hid_t group;
  herr_t status;
  
  int numparts[6] = {0, 0, 0, 0, 0, 0};
  unsigned int unumparts[6] = {0, 0, 0, 0, 0, 0};
  unsigned int zeros[6] = {0, 0, 0, 0, 0, 0};
  double dzeros[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  numparts[1] = header->npart;
  unumparts[1] = header->npart;
  
  group = H5Gcreate(file, "Header", 0);
  if (group < 0)
  {
    printf("Could not open header group!\n");
    exit(1);
  }
  
  write_attr_d(group, "BoxSize", header->boxsize * header->hubbleparam);
  write_attr_d(group, "ExpansionFactor", header->scalefactor);
  write_attr_d(group, "HubbleParam", header->hubbleparam);
  write_attr_d(group, "Omega0", header->Omega0);
  write_attr_d(group, "OmegaBaryon", header->OmegaBaryon);
  write_attr_d(group, "OmegaLambda", header->OmegaLambda);
  write_attr_d(group, "Redshift", header->redshift);
  write_attr_d(group, "Time_GYR", header->time);

  write_attr_i(group, "Flag_Cooling", 0);
  write_attr_i(group, "Flag_Feedback", 0);
  write_attr_i(group, "Flag_Metals", 0);
  write_attr_i(group, "Flag_Sfr", 0);
  write_attr_i(group, "Flag_StellarAge", 0);
  write_attr_i(group, "NumFilesPerSnapshot", 1);
  
  write_attr_array_i(group, "NumPart_ThisFile", numparts, 6);
  write_attr_array_ui(group, "NumPart_Total", unumparts, 6);
  write_attr_array_ui(group, "NumPart_Total_HighWord", zeros, 6);
  write_attr_array_d(group, "MassTable", dzeros, 6);

  status = H5Gclose(group);
  
  return;
}

//////////////////////////////////////////////////////////////////////////

void write_constants(hid_t file)
{
  hid_t group;
  herr_t status;
  
  group = H5Gcreate(file, "Constants", 0);
  if (group < 0)
  {
    printf("Could not open header group!\n");
    exit(1);
  }
  
  write_attr_d(group, "AVOGADRO", 6.0222E23);
  write_attr_d(group, "BOLTZMANN", 1.38066E-16);
  write_attr_d(group, "C", 2.9979E10);
  write_attr_d(group, "CM_PER_MPC", 3.085678E24);
  write_attr_d(group, "ELECTRONCHARGE", 4.8032E-10);
  write_attr_d(group, "ELECTRONMASS", 9.10953E-28);
  write_attr_d(group, "EV_TO_ERG", 1.60217646E-12);
  write_attr_d(group, "GAMMA", 1.6666666666666667);
  write_attr_d(group, "GAS_CONST", 8.31425E7);
  write_attr_d(group, "GRAVITY", 6.672E-8);
  write_attr_d(group, "HUBBLE", 3.2407789E-18);
  write_attr_d(group, "PI", 3.141592653589793);
  write_attr_d(group, "PLANCK", 6.6262E-27);
  write_attr_d(group, "PROTONMASS", 1.6726E-24);
  write_attr_d(group, "RAD_CONST", 7.565E-15);
  write_attr_d(group, "SEC_PER_MEGAYEAR", 3.155E13);
  write_attr_d(group, "SEC_PER_YEAR", 3.155E7);
  write_attr_d(group, "SOLAR_LUM", 3.826E33);
  write_attr_d(group, "SOLAR_MASS", 1.989E33);
  write_attr_d(group, "STEFAN", 7.5657E-15);
  write_attr_d(group, "THOMPSON", 6.6524587E-25);
  write_attr_d(group, "T_CMB0", 2.728);
  write_attr_d(group, "Z_Solar", 0.012663729);

  status = H5Gclose(group);
  
  return;
}

//////////////////////////////////////////////////////////////////////////

void write_units(hid_t file)
{
  hid_t group;
  herr_t status;
  
  group = H5Gcreate(file, "Units", 0);
  if (group < 0)
  {
    printf("Could not open header group!\n");
    exit(1);
  }
  
  write_attr_d(group, "UnitDensity_in_cgs", 6.769911178294543E-31);
  write_attr_d(group, "UnitEnergy_in_cgs", 1.989E53);
  write_attr_d(group, "UnitLength_in_cm", 3.085678E24);
  write_attr_d(group, "UnitMass_in_g", 1.989E43);
  write_attr_d(group, "UnitPressure_in_cgs", 6.769911178294542E-21);
  write_attr_d(group, "UnitTime_in_s", 3.085678E19);
  write_attr_d(group, "UnitVelocity_in_cm_per_s", 100000.0);

  status = H5Gclose(group);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////

void write_parameters(hid_t file)
{
  hid_t group;
  herr_t status;
  
  group = H5Gcreate(file, "Parameters", 0);
  if (group < 0)
  {
    printf("Could not open header group!\n");
    exit(1);
  }
  
  status = H5Gclose(group);
  
  return;
}

//////////////////////////////////////////////////////////////////////////////

void write_attr_array_d(hid_t dset, const char *name, double *data, hsize_t size)
{
  hid_t attr, space;
  hsize_t dims;
  herr_t status;
  
  space = H5Screate_simple(1, &size, NULL);
  attr = H5Acreate(dset, name, H5T_IEEE_F64LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_IEEE_F64LE, data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_array_i(hid_t dset, const char *name, int *data, hsize_t size)
{
  hid_t attr, space;
  hsize_t dims;
  herr_t status;
  
  space = H5Screate_simple(1, &size, NULL);
  attr = H5Acreate(dset, name, H5T_STD_I32LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_STD_I32LE, data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_array_ui(hid_t dset, const char *name, unsigned int *data, hsize_t size)
{
  hid_t attr, space;
  hsize_t dims;
  herr_t status;
  
  space = H5Screate_simple(1, &size, NULL);
  attr = H5Acreate(dset, name, H5T_STD_U32LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_STD_U32LE, data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_i(hid_t dset, const char *name, int data)
{
  hid_t attr, space;
  herr_t status;

  space = H5Screate(H5S_SCALAR);
  attr = H5Acreate(dset, name, H5T_STD_I32LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_STD_I32LE, &data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_d(hid_t dset, const char *name, double data)
{
  hid_t attr, space;
  herr_t status;

  space = H5Screate(H5S_SCALAR);
  attr = H5Acreate(dset, name, H5T_IEEE_F64LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_IEEE_F64LE, &data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_f(hid_t dset, const char *name, float data)
{
  hid_t attr, space;
  herr_t status;

  space = H5Screate(H5S_SCALAR);
  attr = H5Acreate(dset, name, H5T_IEEE_F32LE, space, H5P_DEFAULT);
  status = H5Awrite(attr, H5T_IEEE_F32LE, &data);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}

void write_attr_s(hid_t dset, const char *name, const char *data, hsize_t size)
{
  hid_t attr, space, string_type;
  herr_t status;

  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, size);
  space = H5Screate(H5S_SCALAR);
  attr = H5Acreate(dset, name, string_type, space, H5P_DEFAULT);
  status = H5Awrite(attr, string_type, &data);
  status = H5Tclose(string_type);
  status = H5Sclose(space);
  status = H5Aclose(attr);
}



