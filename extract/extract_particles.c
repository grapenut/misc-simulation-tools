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

#define SAMPLE_FACTOR 0.02

#define STRING_LEN 80
#define SMALL_STRING_LEN 24

#define CM_TO_MPC 3.24077885e-25
#define MPC_TO_CM 3.085678e24
#define CM_TO_KM 1.0e-5
#define SEC_TO_GYR 3.16887646e-17
#define H0_TO_h 3.08568025e17
#define GRAM_TO_OMEGA (1.0 / 1.989e43)

typedef struct
{
  double scale;
  double hubble;
  double boxsize;
  
  char snapshot[200];
} hdata;

typedef struct
{
  int num_dm;
  int num_star;
  int num_baryon;
  
  double **dm_pos;
  double **dm_vel;
  double *dm_mass;
  
  double **star_pos;
  double **star_vel;
  double *star_mass;
  double *star_dens;
  double *star_len;
  
  double **bar_pos;
  double **bar_vel;
  double *bar_mass;
  double *bar_dens;
  double *bar_len;
  double *bar_temp;
  
} pdata;


///////////////////////////////////////////////////////////////////


double *make_array(int);
double **make_2darray(int, int);

hid_t open_file(const char *);
void close_file(hid_t);

hid_t open_group(hid_t, const char *);
void close_group(hid_t);

void read_header(hid_t, hdata *);
void read_snapshot(hdata *, pdata *);

void read_dataset(hid_t, const char *, hdata *, double **, int *);
void read_2d_dataset(hid_t, const char *, hdata *, double ***, int *, int *);
double read_attr(hid_t, const char *);

void write_particles(hdata *, pdata *);


///////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{

  pdata particles;
  hdata header;
  
  sprintf(header.snapshot, argv[1]);
  printf("Reading file %s\n", header.snapshot);
  
  read_snapshot(&header, &particles);

  write_particles(&header, &particles);

  free(particles.bar_pos[0]);
  free(particles.bar_pos);
  free(particles.bar_vel[0]);
  free(particles.bar_vel);
  free(particles.bar_mass);
  free(particles.bar_dens);
  free(particles.bar_temp);
  free(particles.bar_len);

  free(particles.dm_pos[0]);
  free(particles.dm_pos);
  free(particles.dm_vel[0]);
  free(particles.dm_vel);
  free(particles.dm_mass);

  printf("Done.\n");
  return 0;
}

///////////////////////////////////////////////////////////////////////////

void read_snapshot(hdata *header, pdata *particles)
{
  hid_t file, space, grp;
  int cols, rows;

  // Populate particle data structure
  file = open_file(header->snapshot);
  
  read_header(file, header);

  grp = open_group(file, "PartType0");
  read_2d_dataset(grp, "Coordinates", header, &(particles->bar_pos), &rows, &cols);
  read_2d_dataset(grp, "Velocity", header, &(particles->bar_vel), &rows, &cols);
  read_dataset(grp, "Mass", header, &(particles->bar_mass), &rows);
  read_dataset(grp, "Density", header, &(particles->bar_dens), &rows);
  read_dataset(grp, "Temperature", header, &(particles->bar_temp), &rows);
  read_dataset(grp, "SmoothingLength", header, &(particles->bar_len), &rows);
  particles->num_baryon = rows;
  close_group(grp);
  
  grp = open_group(file, "PartType1");
  read_2d_dataset(grp, "Coordinates", header, &(particles->dm_pos), &rows, &cols);
  read_2d_dataset(grp, "Velocity", header, &(particles->dm_vel), &rows, &cols);
  read_dataset(grp, "Mass", header, &(particles->dm_mass), &rows);
  //read_dataset(grp, "SmoothingLength", header, &(particles->dm_len), &rows);
  particles->num_dm = rows;
  close_group(grp);
  
  printf("DEBUG: READING STARS\n");
  
  grp = open_group(file, "PartType4");
  read_2d_dataset(grp, "Coordinates", header, &(particles->star_pos), &rows, &cols);
  read_2d_dataset(grp, "Velocity", header, &(particles->star_vel), &rows, &cols);
  read_dataset(grp, "Mass", header, &(particles->star_mass), &rows);
  read_dataset(grp, "Density", header, &(particles->star_dens), &rows);
  read_dataset(grp, "SmoothingLength", header, &(particles->star_len), &rows);
  particles->num_star = rows;
  close_group(grp);
  
  close_file(file);

}
    
///////////////////////////////////////////////////////////////////////////

void read_header(hid_t file, hdata *header)
{
  hid_t grp;
  
  grp = open_group(file, "Header");
  
  header->scale = read_attr(grp, "ExpansionFactor");
  header->boxsize = read_attr(grp, "BoxSize");
  header->hubble = read_attr(grp, "HubbleParam");
  
  close_group(grp);
  
  printf("Header summary: a=%g, L=%g, H=%g\n", header->scale, header->boxsize, header->hubble);
  
}

///////////////////////////////////////////////////////////////////////////

void read_dataset(hid_t file, const char *name, hdata *header, double **outdata, int *dsize)
{
  hid_t dset, space, attr;
  herr_t status;
  hsize_t dims[2];
  double cgs_factor;
  double scale_exp, scale_factor;
  double h_exp, h_factor;
  double *data;
  int sizex;
  int i;
  
  if (!dsize)
  {
    printf("No return variable for dataset array size!\n");
    exit(1);
  }
  
  dset = H5Dopen(file, name);
  if (dset < 0)
  {
    printf("Could not open input dataset!\n");
    exit(1);
  }
  
  cgs_factor = read_attr(dset, "CGSConversionFactor");
  scale_exp = read_attr(dset, "aexp-scale-exponent");
  h_exp = read_attr(dset, "h-scale-exponent");
  
  scale_factor = pow(header->scale, scale_exp);
  h_factor = pow(header->hubble, h_exp);

  printf("DEBUG: cgs=%g, scale_exp=%g, h_exp=%g, a=%g, h=%g\n", cgs_factor, scale_exp, h_exp, scale_factor, h_factor);
  
  space = H5Dget_space(dset);
               
  status = H5Sget_simple_extent_dims(space, dims, NULL);
  sizex = dims[0];
  *dsize = sizex;
  
  printf("Reading 1d dataset %s, rows=%d\n", name, *dsize);

  data = make_array(sizex);

  status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Sclose(space);
  status = H5Dclose(dset);
  
  for (i = 0; i < sizex; i++)
  {
    data[i] = data[i] * cgs_factor * scale_factor * h_factor;
  }

  *outdata = data;
}

///////////////////////////////////////////////////////////////////////////

void read_2d_dataset(hid_t file, const char *name, hdata *header, double ***outdata, int *rows, int *cols)
{
  hid_t dset, space, attr;
  herr_t status;
  hsize_t dims[2];
  double cgs_factor;
  double scale_exp, scale_factor;
  double h_exp, h_factor;
  double dmin, dmax;
  double **data;
  int sizex, sizey;
  int i, j;
  
  if (!cols || !rows)
  {
    printf("No return variable for dataset array size!\n");
    exit(1);
  }
  
  dset = H5Dopen(file, name);
  if (dset < 0)
  {
    printf("Could not open input dataset!\n");
    exit(1);
  }
  
  cgs_factor = read_attr(dset, "CGSConversionFactor");
  scale_exp = read_attr(dset, "aexp-scale-exponent");
  h_exp = read_attr(dset, "h-scale-exponent");
  
  scale_factor = pow(header->scale, scale_exp);
  h_factor = pow(header->hubble, h_exp);
  
  printf("DEBUG: cgs=%g, scale_exp=%g, h_exp=%g, a=%g, h=%g\n", cgs_factor, scale_exp, h_exp, scale_factor, h_factor);
  
  space = H5Dget_space(dset);
               
  status = H5Sget_simple_extent_dims(space, dims, NULL);
  sizex = dims[0];
  sizey = dims[1];
  *rows = sizex;
  *cols = sizey;
  
  printf("Reading 2d dataset %s, rows=%d, cols=%d\n", name, *rows, *cols);
  
  data = make_2darray(sizex, sizey);

  status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data[0]);
  status = H5Sclose(space);
  status = H5Dclose(dset);
  
  dmin = 1.0e100;
  dmax = -1.0e100;
  
  for (i = 0; i < sizex; i++)
  {
    for (j = 0; j < sizey; j++)
    {
      data[i][j] = data[i][j] * cgs_factor * scale_factor * h_factor;
      
      if (data[i][j] > dmax) dmax = data[i][j];
      
      if (data[i][j] < dmin) dmin = data[i][j];

      //if (i == 0 && j == 0) printf("DEBUG: data[0][0] = %g\n" , data[i][j]);
    }
  }
  
  printf("DEBUG: dmin = %g, dmax = %g\n", dmin, dmax);

  *outdata = data;
}

///////////////////////////////////////////////////////////////////////////

hid_t open_file(const char *name)
{
  hid_t file;
  
  file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
  {
    printf("Could not open input file!\n");
    exit(1);
  }

  return file;
}

void close_file(hid_t file)
{
  herr_t status;
  
  status = H5Fclose(file);
}

///////////////////////////////////////////////////////////////////////////

hid_t open_group(hid_t file, const char *name)
{
  hid_t grp;
  
  grp = H5Gopen(file, name);
  if (grp < 0)
  {
    printf("Could not open input file!\n");
    exit(1);
  }

  return grp;
}

void close_group(hid_t grp)
{
  herr_t status;
  
  status = H5Gclose(grp);
}

///////////////////////////////////////////////////////////////////////////

double read_attr(hid_t file, const char *name)
{
  hid_t attr;
  herr_t status;
  double ret;
  
  attr = H5Aopen_name(file, name);
  if (attr < 0)
  {
    printf("Could not open attribute!\n");
    return 1.0;
  }

  status = H5Aread(attr, H5T_IEEE_F64LE, &ret);
  status = H5Aclose(attr);
  
  return ret;
}

///////////////////////////////////////////////////////////////////////////

void write_particles(hdata *header, pdata *particles)
{
  FILE *fp;
  int i;
  
  //printf("HALO: %g %g %g (%g %g %g) %g %g %g\n", halo->x, halo->y, halo->z, halo->cm_x, halo->cm_y, halo->cm_z, halo->vx, halo->vy, halo->vz);
  
  fp = fopen("baryon.txt", "w");
  if (!fp)
  {
    printf("Unable to open baryon output file!\n");
    exit(1);
  }

  for (i = 0; i < particles->num_baryon; i++)
  {
    fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", particles->bar_pos[i][0], particles->bar_pos[i][1], particles->bar_pos[i][2],
        particles->bar_vel[i][0], particles->bar_vel[i][1], particles->bar_vel[i][2], particles->bar_mass[i],
        particles->bar_dens[i], particles->bar_len[i], particles->bar_temp[i]);
  }

  for (i = 0; i < particles->num_star; i++)
  {
    //fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", particles->star_pos[i][0], particles->star_pos[i][1], particles->star_pos[i][2],
    //    particles->star_vel[i][0], particles->star_vel[i][1], particles->star_vel[i][2], particles->star_mass[i]);
    fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", particles->star_pos[i][0], particles->star_pos[i][1], particles->star_pos[i][2],
        particles->star_vel[i][0], particles->star_vel[i][1], particles->star_vel[i][2], particles->star_mass[i],
        particles->star_dens[i], particles->star_len[i], 1.0e4);
  }

  fclose(fp);


  fp = fopen("darkmatter.txt", "w");
  if (!fp)
  {
    printf("Unable to open output file!\n");
    exit(1);
  }

  for (i = 0; i < particles->num_dm; i++)
  {
    fprintf(fp, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", particles->dm_pos[i][0], particles->dm_pos[i][1], particles->dm_pos[i][2],
        particles->dm_vel[i][0], particles->dm_vel[i][1], particles->dm_vel[i][2], particles->dm_mass[i]);
  }
  
  fclose(fp);
}

///////////////////////////////////////////////////////////////////////////

double *make_array(int rows)
{
  double *data;
  
  data = (double *) malloc(rows * sizeof(double));
  
  return data;
}

///////////////////////////////////////////////////////////////////////////

double **make_2darray(int rows, int cols)
{
  int i;
  double **data;
  
  data = (double **) malloc(rows * sizeof(double *));
  data[0] = (double *) malloc(rows * cols * sizeof(double));
  for (i = 1; i < rows; i++)
    data[i] = data[0] + i*cols;
  
  return data;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

