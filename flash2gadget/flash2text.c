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
  double redshift;
  double time;
  double x, y, z, r;
  int map_x;
  int map_y;
  int map_z;
  int map_mass;
  int map_vx;
  int map_vy;
  int map_vz;
  int map_dens;
  int map_tag;
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

int main(int argc, char **argv)
{

  double **flash_data;
  int npart;
  char infile_name[200], outfile_name[200];
  int i, j;
  hdata header;
  
  if (argc != 7)
  {
    printf("Usage: %s inputfile outputfile x y z r\n", argv[0]);
    exit(1);
  }
  
  sprintf(infile_name, argv[1]);
  sprintf(outfile_name, argv[2]);
  
  header.x = atof(argv[3]);
  header.y = atof(argv[4]);
  header.z = atof(argv[5]);
  header.r = atof(argv[6]);

  printf("%g %g %g %g\n", header.x, header.y, header.z, header.r);

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
  header->time = read_scalar(file, "real scalars", "time");
  header->boxsize = read_scalar(file, "real runtime parameters", "xmax");

  printf("  Reading particle map.\n");
  header->map_x = read_particle_map(file, "posx");
  header->map_y = read_particle_map(file, "posy");
  header->map_z = read_particle_map(file, "posz");
  header->map_mass = read_particle_map(file, "mass");
  header->map_vx = read_particle_map(file, "velx");
  header->map_vy = read_particle_map(file, "vely");
  header->map_vz = read_particle_map(file, "velz");
  header->map_dens = read_particle_map(file, "grid_dens");
  header->map_tag = read_particle_map(file, "tag");

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
  int i, npart;
  
  double x, y, z, dx, dy, dz, vx, vy, vz, mass, dens, tag, dist;
  
  FILE *fp;

  npart = header->npart;

  if((fp=fopen(filename, "w")) == NULL) {
    printf("Cannot open output file.\n");
    exit(1);
  }

  for (i=0; i < npart; i++)
  {
    x = (float)( data[i][header->map_x]);
    y = (float)( data[i][header->map_y]);
    z = (float)( data[i][header->map_z]);

    vx = (float)( data[i][header->map_vx]);
    vy = (float)( data[i][header->map_vy]);
    vz = (float)( data[i][header->map_vz]);

    mass = (float)( data[i][header->map_mass]);
    
    dens = (float)( data[i][header->map_dens]);
    
    tag = (float)( data[i][header->map_tag]);
     
    dx = x - header->x;
    dy = y - header->y;
    dz = z - header->z;
    
    dist = sqrt(dx*dx + dy*dy + dz*dz);
     
    if (dist <= header->r)
    {
      fprintf(fp, "%17.8g %17.8g %17.8g %17.8g\n", tag, x, y, z);
    }
  }
  
  fclose(fp);

  return;
}

//////////////////////////////////////////////////////////////////////////

