

#include "cddefines.h"
#include "cddrive.h"

#include "dense.h"

#include "cloudy.h"

void run_cloudy(data& result)
{

  int ierr;
  char buff[256];

  cdOutput("summary.out");

  cdInit();

  cdTalk(true);

  // chemical abundances, primordial + carbon/oxygen
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
  cdRead("no molecules");


  // input parameters
  // metal scaling must come after chemical abundances have been set
  sprintf(buff, "metals %f log", log10(metal));
  cdRead(buff);

  sprintf(buff, "hden %f log", log10(dens));
  cdRead(buff);

  sprintf(buff, "constant temper %f log", log10(temp));
  cdRead(buff);

  sprintf(buff, "ionization parameter %f log", log10(flux));
  cdRead(buff);

  sprintf(buff, "age %f log years", log10(age));
  cdRead(buff);

  // source spectrum
  
  // single temperature blackbody (from 2012, 2014, 2016 papers)
  //cdRead("blackbody 4.970 log");

  // might be needed for molecular chemistry not to crash
  //cdRead("cosmic ray background -10 log");

  // starburst99 generic SED for M=10^5, z=0.001
  sprintf(buff, "table stars \"z10galaxy.mod\" age %f log", log10(age));
  cdRead(buff);

  // geometry
  // we are just reading off 1-zone values so we can take them from the first zone
  cdRead("stop zone 1");
  cdRead("set dr 0");

  // run cloudy with the above parameters
  ierr = cdDrive();

  result.eden = cdEDEN_last();

  result.cooling = cdCooling_last() / dens / dens;

  result.heating = cdHeating_last() / dens / dens;

  result.h = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN]);
  result.hplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]);

  result.he = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][0]/dense.gas_phase[ipHELIUM]);
  result.heplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][1]/dense.gas_phase[ipHELIUM]);
  result.heplusplus = MAX2(SMALLFLOAT,dense.xIonDense[ipHELIUM][2]/dense.gas_phase[ipHELIUM]);

  result.c = MAX2(SMALLFLOAT,dense.xIonDense[ipCARBON][0]/dense.gas_phase[ipCARBON]);
  result.cplus = MAX2(SMALLFLOAT,dense.xIonDense[ipCARBON][1]/dense.gas_phase[ipCARBON]);



}


