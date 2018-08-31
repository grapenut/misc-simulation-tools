
// generate data tables for the PhotoChem variant of FLASH chemistry

#include "cddefines.h"
#include "cddrive.h"

#include <iostream>

//#include "cpu.h"
//#include "hmi.h"
//#include "mole.h"
#include "dense.h"

// number of data points per decade in log-space
#define NH_PER_DECADE		1
#define TEMP_PER_DECADE		1
#define METAL_PER_DECADE	1
#define FLUX_PER_DECADE		1

// log10 minima/maxima for parameter space
#define TEMP_MIN	1.0
#define TEMP_MAX	8.0

#define NH_MIN		-4.0
#define NH_MAX		8.0

#define METAL_MIN	-6.0
#define METAL_MAX	3.0

#define FLUX_MIN	-15.0
#define FLUX_MAX	10.0

using std::vector;

int main(int argc, char **argv)
{
  vector<double> temp_list, metal_list, nh_list, flux_list;
  
  double step, temp, flux, metal, heating, cooling, nh, eden, h, hplus, he, heplus, heplusplus;
  
  int num_temp, num_metal, num_nh, num_flux;
  
  int t, z, n, f;
  
  char filename[256];

  int ierr;

  FILE *fout;
  char buff[100];

  if (argc < 4) {
    printf("ERROR: need to supply n, z, and t indices.\n");
  }
  
  n = atoi(argv[1]);
  z = atoi(argv[2]);
  t = atoi(argv[3]);
  
  step = 1.0/TEMP_PER_DECADE;
  for (temp = TEMP_MIN; temp <= TEMP_MAX; temp += step)
  {
    temp_list.push_back(pow(10.0, temp));
    //printf("TEMP %f\n", temp_list.back());
  }
  num_temp = temp_list.size();

  step = 1.0/FLUX_PER_DECADE;
  for (flux = FLUX_MIN; flux <= FLUX_MAX; flux += step)
  {
    flux_list.push_back(pow(10.0, flux));
    //printf("FLUX %f\n", flux_list.back());
  }
  num_flux = flux_list.size();

  step = 1.0/METAL_PER_DECADE;
  for (metal = METAL_MIN; metal <= METAL_MAX; metal += step)
  {
    metal_list.push_back(pow(10.0, metal));
    //printf("METAL %f\n", metal_list.back());
  }
  num_metal = metal_list.size();

  step = 1.0/NH_PER_DECADE;
  for (nh = NH_MIN; nh <= NH_MAX; nh += step)
  {
    nh_list.push_back(pow(10.0, nh));
    //printf("nH %f\n", nh_list.back());
  }
  num_nh = nh_list.size();
  
  printf("n = %d \t z = %d \t t = %d \t f = %d\n", num_nh, num_metal, num_temp, num_flux);

  //printf("Done with init. num_total = %d\n", num_nh*num_metal*num_temp*num_flux);

  // ready to run cloudy
  cdOutput("summary.out" );
  temp = temp_list[t];
  metal = metal_list[z];
  nh = nh_list[n];

  /* file to output short form of calculation's results */
  sprintf(filename, "output/photochem_%d_%d_%d.txt", n, z, t);
  printf("Outputting file %s.\n", filename);

  if ((fout = fopen(filename, "w")) == NULL)
  {
    printf("ERROR: Could not open %s for writing.\n", filename);
    cdEXIT(EXIT_FAILURE);
  }

  for (f = 0; f < num_flux; f++)
  {
    flux = flux_list[f];

    cdInit();

    cdTalk(true);
    
    cdRead("element helium abundance 0.08 linear");
    cdRead("element carbon scale -3.0 log");
    cdRead("element oxygen scale -3.0 log");
    
    cdRead("element nitrogen off");
    cdRead("element neon off");
    cdRead("element sodium off");
    cdRead("element magnesium off");
    cdRead("element aluminum off");
    cdRead("element silicon off");
    cdRead("element sulphur off");
    cdRead("element argon off");
    cdRead("element calcium off");
    cdRead("element iron off");
    cdRead("element nickel off");
    cdRead("element Lithium off");
    cdRead("element Beryllium off");
    cdRead("element Boron off");
    cdRead("element Fluorine off");
    cdRead("element Phosphor off");
    cdRead("element Chlorine off");
    cdRead("element Potassium off");
    cdRead("element Scandium off");
    cdRead("element Titanium off");
    cdRead("element Vanadium off");
    cdRead("element Chromium off");
    cdRead("element Manganese off");
    cdRead("element Cobalt off");
    cdRead("element Copper off");
    cdRead("element Zinc off");

    sprintf(buff, "metals %f log", log10(metal));
    cdRead(buff);

    cdRead("cmb 10 "  );
    cdRead("no molecules " );
    cdRead("cosmic ray background -10 log" );
    //cdRead("no photoionization " );
    //cdRead("no compton effect " );
    cdRead("age 6 log years " );

    cdRead("stop zone 1 "  );
    cdRead("set dr 0 "  );
    
    sprintf(buff, "hden %f ", log10(nh));
    cdRead(buff);

    sprintf(buff, "constant temper %f log", log10(temp));
    cdRead(buff);
    
    //sprintf(buff, "ionization parameter %f", log10(flux));
    //cdRead(buff);

    //cdRead("punch heating \"summary.het\" last no hash no clobber "  );
    //cdRead("punch cooling \"summary.col\" last no hash no clobber "  );
    //cdRead("punch overview \"summary.ovr\" last no hash no clobber "  );
    
    ierr = cdDrive();
    //printf("===== %d,%d,%d,%d ===== problems!!\n", n, z, t, f);
    
    eden = cdEDEN_last();

    cooling = cdCooling_last() / nh / nh;
    
    heating = cdHeating_last() / nh / nh;
    
    /* molecular fraction of hydrogen */
    //h2frac = MAX2(SMALLFLOAT,2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN]);

    /* carbon monoxide molecular fraction of CO */
    //cofrac = SDIV(dense.gas_phase[ipCARBON]);
    //cofrac = findspecies("CO")->hevmol/cofrac;
    //cofrac = MAX2(SMALLFLOAT, cofrac );
    
    h = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN]);
    hplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]);

    he = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][0]/dense.gas_phase[ipHELIUM]);
    heplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][1]/dense.gas_phase[ipHELIUM]);
    heplusplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][2]/dense.gas_phase[ipHELIUM]);

    //printf("%-10g %-10g %-10g %-10g = ne %-10g cl %-10g ht %-10g h %-10g / %-10g he %-10g / %-10g / %-10g\n", nh, metal, temp, flux, eden, cooling, heating, h, hplus, he, heplus, heplusplus);
    
    fprintf(fout, "%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n", nh, metal, temp, flux, eden, cooling, heating, h, hplus, he, heplus, heplusplus);
  }

  fclose(fout);

  printf("Done.\n");

}

