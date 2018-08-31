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
  
  char subfind_base[200];
  char snapshot_base[200];
  char output_base[200];
  int num_files;
} hdata;

typedef struct
{
  int num_dm;
  int num_baryon;
  
  double **dm_pos;
  double **dm_vel;
  double *dm_mass;
  
  double **bar_pos;
  double **bar_vel;
  double *bar_mass;
  double *bar_dens;
  double *bar_len;
  double *bar_temp;
  
} pdata;

typedef struct
{
  int file_num;
  int halo_num;
  double mass;
  double x, y, z;
  double cm_x, cm_y, cm_z;
  double vx, vy, vz;
} sdata;


///////////////////////////////////////////////////////////////////


double *make_array(int);
double **make_2darray(int, int);

hid_t open_file(const char *);
void close_file(hid_t);

hid_t open_group(hid_t, const char *);
void close_group(hid_t);

void read_header(hid_t, hdata *);
void read_subfind(hdata *, sdata *);
void read_snapshot(hdata *, sdata *, pdata *);

void read_dataset(hid_t, const char *, hdata *, double **, int *);
void read_2d_dataset(hid_t, const char *, hdata *, double ***, int *, int *);
double read_attr(hid_t, const char *);

void write_particles(hdata *, sdata *, pdata *);


///////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{

  pdata particles;
  sdata halo;
  hdata header;
  
  if (argc != 5)
  {
    printf("Usage: %s subfind_base snapshot_base num_files output_file\n", argv[0]);
    exit(1);
  }
  
  sprintf(header.subfind_base, argv[1]);
  sprintf(header.snapshot_base, argv[2]);
  
  header.num_files = atoi(argv[3]);
  
  sprintf(header.output_base, argv[4]);

  read_subfind(&header, &halo);
  read_snapshot(&header, &halo, &particles);

  write_particles(&header, &halo, &particles);

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

void read_subfind(hdata *header, sdata *halo)
{
  hid_t file, grp;
  int sizex;
  double *data;
  char filename[200];
  
  int n, n_max, i, i_max;
  double m_max;
    
  m_max = 0.0;

  for (n = 0; n < header->num_files; n++)
  {
  
    sprintf(filename, "%s.%d.hdf5", header->subfind_base, n);
    
    file = open_file(filename);
    grp = open_group(file, "SUBFIND");
    
    read_header(file, header);

    read_dataset(grp, "Mass", header, &data, &sizex);
    
    for (i=0; i<sizex; i++)
    {
      if (data[i] > m_max)
      {
        m_max = data[i];
        n_max = n;
        i_max = i;
      }
    }
    
    free(data);
    
    close_group(grp);
    close_file(file);
  }
    
  halo->mass = m_max;
  halo->file_num = n_max;
  halo->halo_num = i_max;

  sprintf(filename, "%s.%d.hdf5", header->subfind_base, n_max);
  
  file = open_file(filename);
  grp = open_group(file, "SUBFIND");
  

  read_dataset(grp, "Position", header, &data, &sizex);
  
  halo->x = data[3 * i_max];
  halo->y = data[3 * i_max + 1];
  halo->z = data[3 * i_max + 2];

  free(data);


  read_dataset(grp, "CenterOfMass", header, &data, &sizex);
  
  halo->cm_x = data[3 * i_max];
  halo->cm_y = data[3 * i_max + 1];
  halo->cm_z = data[3 * i_max + 2];

  free(data);


  read_dataset(grp, "CenterOfMassVelocity", header, &data, &sizex);
  
  halo->vx = data[3 * i_max];
  halo->vy = data[3 * i_max + 1];
  halo->vz = data[3 * i_max + 2];
  
  free(data);
  

  close_group(grp);
  close_file(file);
    

}

///////////////////////////////////////////////////////////////////////////

void read_snapshot(hdata *header, sdata *halo, pdata *particles)
{
  hid_t file, space, grp;
  int cols, rows, i, j;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double **pos, **vel;
  double *mass, *dens, *len, *temp;
  char filename[200];
  int n;
  int num_dm, num_dm_total, num_baryon, num_baryon_total;
  int cur_dm, cur_bar;
  double factor;
  double xoffset, yoffset, zoffset;
  
  xoffset = 1.0e100;
  yoffset = 1.0e100;
  zoffset = 1.0e100;

  num_dm = 0;
  num_dm_total = 0;
  
  num_baryon = 0;
  num_baryon_total = 0;
  
  cur_bar = 0;
  cur_dm = 0;

  factor = SAMPLE_FACTOR * 0.5 * MPC_TO_CM * header->scale;

  xmin = halo->x - factor;
  xmax = halo->x + factor;
  ymin = halo->y - factor;
  ymax = halo->y + factor;
  zmin = halo->z - factor;
  zmax = halo->z + factor;

  for (n = 0; n < header->num_files; n++)
  {
  
    sprintf(filename, "%s.%d.hdf5", header->snapshot_base, n);
    file = open_file(filename);

    grp = open_group(file, "PartType0");
    read_2d_dataset(grp, "Coordinates", header, &pos, &cols, &rows);
    close_group(grp);
    
    for (i=0; i<cols; i++)
    {
      num_baryon_total++;
      
      if (pos[i][0] > xmin && pos[i][0] < xmax &&
          pos[i][1] > ymin && pos[i][1] < ymax &&
          pos[i][2] > zmin && pos[i][2] < zmax)
      {
        if (pos[i][0] < xoffset) xoffset = pos[i][0];
        if (pos[i][1] < yoffset) yoffset = pos[i][1];
        if (pos[i][2] < zoffset) zoffset = pos[i][2];

        num_baryon++;
      }
    }
    
    free(pos[0]);
    free(pos);
    

    grp = open_group(file, "PartType1");
    read_2d_dataset(grp, "Coordinates", header, &pos, &cols, &rows);
    close_group(grp);
    
    for (i=0; i<cols; i++)
    {
      num_dm_total++;
      
      if (pos[i][0] > xmin && pos[i][0] < xmax &&
          pos[i][1] > ymin && pos[i][1] < ymax &&
          pos[i][2] > zmin && pos[i][2] < zmax)
      {
        if (pos[i][0] < xoffset) xoffset = pos[i][0];
        if (pos[i][1] < yoffset) yoffset = pos[i][1];
        if (pos[i][2] < zoffset) zoffset = pos[i][2];

        num_dm++;
      }
    }
    
    free(pos[0]);
    free(pos);
    
    close_file(file);
  }
  
  halo->x -= xoffset;
  halo->cm_x -= xoffset;
  halo->y -= yoffset;
  halo->cm_y -= yoffset;
  halo->z -= zoffset;
  halo->cm_z -= zoffset;
  
  printf("Found %d out of %d baryons.\n", num_baryon, num_baryon_total);
  printf("Found %d out of %d dark matter.\n", num_dm, num_dm_total);
  
  
  // Initialize particle data structure
  particles->num_baryon = num_baryon;
  particles->num_dm = num_dm;
  
  particles->dm_pos = make_2darray(num_dm, 3);
  particles->dm_vel = make_2darray(num_dm, 3);
  particles->dm_mass = make_array(num_dm);

  particles->bar_pos = make_2darray(num_baryon, 3);
  particles->bar_vel = make_2darray(num_baryon, 3);
  particles->bar_mass = make_array(num_baryon);
  particles->bar_dens = make_array(num_baryon);
  particles->bar_len = make_array(num_baryon);
  particles->bar_temp = make_array(num_baryon);

  // Populate particle data structure
  for (n = 0; n < header->num_files; n++)
  {
  
    sprintf(filename, "%s.%d.hdf5", header->snapshot_base, n);
    file = open_file(filename);

    grp = open_group(file, "PartType0");
    read_2d_dataset(grp, "Coordinates", header, &pos, &cols, &rows);
    read_2d_dataset(grp, "Velocity", header, &vel, &cols, &rows);
    read_dataset(grp, "Mass", header, &mass, &cols);
    read_dataset(grp, "Density", header, &dens, &cols);
    read_dataset(grp, "Temperature", header, &temp, &cols);
    read_dataset(grp, "SmoothingLength", header, &len, &cols);
    close_group(grp);
    
    for (i=0; i<cols; i++)
    {
      if (pos[i][0] > xmin && pos[i][0] < xmax &&
          pos[i][1] > ymin && pos[i][1] < ymax &&
          pos[i][2] > zmin && pos[i][2] < zmax)
      {
        particles->bar_pos[cur_bar][0] = pos[i][0] - xoffset;
        particles->bar_pos[cur_bar][1] = pos[i][1] - yoffset;
        particles->bar_pos[cur_bar][2] = pos[i][2] - zoffset;

        particles->bar_vel[cur_bar][0] = vel[i][0] - halo->vx;
        particles->bar_vel[cur_bar][1] = vel[i][1] - halo->vy;
        particles->bar_vel[cur_bar][2] = vel[i][2] - halo->vz;
        
        particles->bar_mass[cur_bar] = mass[i];
        particles->bar_dens[cur_bar] = dens[i];
        particles->bar_temp[cur_bar] = temp[i];
        particles->bar_len[cur_bar] = len[i];
      
        cur_bar++;
      }
    }
    
    free(pos[0]);
    free(pos);
    free(vel[0]);
    free(vel);
    free(mass);
    free(dens);
    free(temp);
    free(len);
    

    grp = open_group(file, "PartType1");
    read_2d_dataset(grp, "Coordinates", header, &pos, &cols, &rows);
    read_2d_dataset(grp, "Velocity", header, &vel, &cols, &rows);
    read_dataset(grp, "Mass", header, &mass, &cols);
    close_group(grp);
    
    for (i=0; i<cols; i++)
    {
      if (pos[i][0] > xmin && pos[i][0] < xmax &&
          pos[i][1] > ymin && pos[i][1] < ymax &&
          pos[i][2] > zmin && pos[i][2] < zmax)
      {
        particles->dm_pos[cur_dm][0] = pos[i][0] - xoffset;
        particles->dm_pos[cur_dm][1] = pos[i][1] - yoffset;
        particles->dm_pos[cur_dm][2] = pos[i][2] - zoffset;

        particles->dm_vel[cur_dm][0] = vel[i][0] - halo->vx;
        particles->dm_vel[cur_dm][1] = vel[i][1] - halo->vy;
        particles->dm_vel[cur_dm][2] = vel[i][2] - halo->vz;
        
        particles->dm_mass[cur_dm] = mass[i];

        cur_dm++;
      }
    }
    
    free(pos[0]);
    free(pos);
    free(vel[0]);
    free(vel);
    free(mass);
    
    close_file(file);
  }
  
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
  
  space = H5Dget_space(dset);
               
  status = H5Sget_simple_extent_dims(space, dims, NULL);
  sizex = dims[0];
  *dsize = sizex;
  
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

void read_2d_dataset(hid_t file, const char *name, hdata *header, double ***outdata, int *cols, int *rows)
{
  hid_t dset, space, attr;
  herr_t status;
  hsize_t dims[2];
  double cgs_factor;
  double scale_exp, scale_factor;
  double h_exp, h_factor;
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
  
  space = H5Dget_space(dset);
               
  status = H5Sget_simple_extent_dims(space, dims, NULL);
  sizex = dims[0];
  sizey = dims[1];
  *cols = sizex;
  *rows = sizey;
  
  data = make_2darray(sizex, sizey);

  status = H5Dread(dset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data[0]);
  status = H5Sclose(space);
  status = H5Dclose(dset);
  
  for (i = 0; i < sizex; i++)
  {
    for (j = 0; j < sizey; j++)
    {
      data[i][j] = data[i][j] * cgs_factor * scale_factor * h_factor;
    }
  }

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

void write_particles(hdata *header, sdata *halo, pdata *particles)
{
  char filename[200];
  FILE *fp;
  int i;
  
  printf("HALO: %g %g %g (%g %g %g) %g %g %g\n", halo->x, halo->y, halo->z, halo->cm_x, halo->cm_y, halo->cm_z, halo->vx, halo->vy, halo->vz);
  
  sprintf(filename, "%s_baryon.txt", header->output_base);
  fp = fopen(filename, "w");
  if (!fp)
  {
    printf("Unable to open output file!\n");
    exit(1);
  }

  for (i = 0; i < particles->num_baryon; i++)
  {
    fprintf(fp, "%g %g %g %g %g %g %g %g %g %g\n", particles->bar_pos[i][0], particles->bar_pos[i][1], particles->bar_pos[i][2],
        particles->bar_vel[i][0], particles->bar_vel[i][1], particles->bar_vel[i][2], particles->bar_mass[i],
        particles->bar_dens[i], particles->bar_len[i], particles->bar_temp[i]);
  }

  fclose(fp);


  sprintf(filename, "%s_dm.txt", header->output_base);
  fp = fopen(filename, "w");
  if (!fp)
  {
    printf("Unable to open output file!\n");
    exit(1);
  }

  for (i = 0; i < particles->num_dm; i++)
  {
    fprintf(fp, "%g %g %g %g %g %g %g\n", particles->dm_pos[i][0], particles->dm_pos[i][1], particles->dm_pos[i][2],
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

