
//#define USE_BLACKBODY

#include "cddefines.h"
#include "cddrive.h"
#include "timesc.h"
#include "dense.h"

#include "cloudy.h"

#define LIGHTSPEED	2.99792458e10

void run_cloudy(data &result)
{
  int ierr;
  char buff[256];

  cdOutput("summary.out");

  cdInit();

  cdTalk(true);

  ///////////////////////////////////////////////////////////////////
  // chemical abundances, primordial + carbon/oxygen
  cdRead("element helium abundance 0.08 linear");
  cdRead("element carbon abundance 2.45e-4 linear");
  cdRead("element oxygen abundance 4.90e-4 linear");

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
  cdRead("no molecules");

  // metal scaling must come after chemical abundances have been set
  sprintf(buff, "metals %f log", log10(result.xmetal));
  cdRead(buff);

  ///////////////////////////////////////////////////////////////////
  // cloud thermodynamic/radiative state
  sprintf(buff, "hden %f log", log10(result.dens));
  cdRead(buff);

  //sprintf(buff, "ionization parameter %f log", log10(result.flux/LIGHTSPEED));
  //cdRead(buff);

  // set constant temperature for reading off cooling rates
  sprintf(buff, "constant temper %f log", log10(result.temp));
  cdRead(buff);

  // cloud age (only for metal cooling rates?)
  sprintf(buff, "age %f log years", log10(result.age));
  cdRead(buff);

  ///////////////////////////////////////////////////////////////////
  // source spectrum
#ifdef USE_BLACKBODY
  // single temperature blackbody (from 2012, 2014, 2016 papers)
  //cdRead("blackbody 4.62325 log");
#else
  // starburst99 generic SED for M=10^5, z=0.001
  //sprintf(buff, "table stars \"z10burst.mod\" age %f log", log10(result.age));
  //cdRead(buff);
#endif

  // might be needed for molecular chemistry not to crash
  // some basic minimal radiation field settings
  cdRead("cosmic ray background -10 log");
  cdRead("cmb 10");
  cdRead("no photoionization");
  cdRead("no compton effect");

  // geometry
  // we are just reading off 1-zone values so we can take them from the first zone
  //cdRead("iterate to convergence");
  cdRead("stop zone 1");
  cdRead("set dr 0");

  // run cloudy with the above parameters
  ierr = cdDrive();
  
  //result.temp = cdTemp_last();;

  result.eden = cdEDEN_last();
  
  result.cooling = cdCooling_last() / result.dens / result.dens;
  result.heating = cdHeating_last() / result.dens / result.dens;
  
  result.thermal_time = timesc.time_therm_long;
  result.recomb_time = timesc.time_Hrecom_long;
  
  result.h = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN]);
  result.hplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]);

  result.he = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][0]/dense.gas_phase[ipHELIUM]);
  result.heplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][1]/dense.gas_phase[ipHELIUM]);
  result.heplusplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][2]/dense.gas_phase[ipHELIUM]);

  result.c = MAX2(SMALLFLOAT,dense.xIonDense[ipCARBON][0]/dense.gas_phase[ipCARBON]);
  result.cplus = MAX2(SMALLFLOAT,dense.xIonDense[ipCARBON][1]/dense.gas_phase[ipCARBON]);
  result.cplusplus = MAX2(SMALLFLOAT,dense.xIonDense[ipCARBON][2]/dense.gas_phase[ipCARBON]);

}


