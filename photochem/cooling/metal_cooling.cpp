/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* runs pure collisional models at range of temperatures, prints cooling */
#include "cddefines.h"
#include "cddrive.h"

#include "cpu.h"
#include "hmi.h"
#include "mole.h"
#include "dense.h"

// production values
// dens 3 per decade +1
// temp 5 per decade +1
// metal 3 per decade +1

#define DENS_PER_DECADE		3
#define TEMP_PER_DECADE		5
#define METAL_PER_DECADE	3
#define FLUX_PER_DECADE		5

#define TEMP_MIN	1.0
#define TEMP_MAX	9.0

#define HDEN_MIN	-4.0
#define HDEN_MAX	8.0

#define METAL_MIN	-6.0
#define METAL_MAX	-0.66666667

#define FLUX_MIN	0.0
#define FLUX_MAX	6.0

/*int main( int argc, char *argv[] )*/
int main( void )
{

  double temp_list[NTEMP];
  double metal_list[NMETAL];
  double hden_list[NHDEN];
  
  double cool_list[NTEMP][NMETAL][NHDEN];
  double heat_list[NTEMP][NMETAL][NHDEN];
  double eden_list[NTEMP][NMETAL][NHDEN];
  double h_list[NTEMP][NMETAL][NHDEN];
  double hplus_list[NTEMP][NMETAL][NHDEN];
  double he_list[NTEMP][NMETAL][NHDEN];
  double heplus_list[NTEMP][NMETAL][NHDEN];
  double heplusplus_list[NTEMP][NMETAL][NHDEN];
  
  double metal, telog, heating , cooling, hden , eden, h, hplus, he, heplus, heplusplus;
  
  double helium_frac = 0.08;
  double carbon_frac = 2.45e-4;
  double oxygen_frac = 4.90e-4;
  double metal_scale = (1.0+4.0*helium_frac) / (14.0*carbon_frac + 16.0*oxygen_frac);
  double xmetal;
  
  int i;
  int j;
  int k;
  
  char filename[256];
  

	int exit_status = EXIT_FAILURE;

	DEBUG_ENTRY( "main()" );

	try {
		int lgBAD ;

		FILE *ioRES;
		char chLine[100];

		/* this will be limit to the number of command chLines we can still put in */
		long int nleft;


		/* we do not want to generate any output */
		cdOutput( "metal_cooling.out" );

                for (i=0;i<NTEMP;i++)
                {
                  temp_list[i] = TEMP_MIN + ((i-0.0)/(NTEMP-1.0))*(TEMP_MAX-TEMP_MIN);
                  printf("TEMP[%d] == %f\n", i, pow(10.0,temp_list[i]));
                }

                for (j=0;j<NMETAL;j++)
                {
                  metal_list[j] = METAL_MIN + ((j-0.0)/(NMETAL-1.0))*(METAL_MAX-METAL_MIN);
                  printf("METAL[%d] == %f\n", j, pow(10.0,metal_list[j]));
                }

                for (k=0;k<NHDEN;k++)
                {
                  hden_list[k] = HDEN_MIN + ((k-0.0)/(NHDEN-1.0))*(HDEN_MAX-HDEN_MIN);
                  printf("HDEN[%d] == %f\n", k, pow(10.0,hden_list[k]));
                }
                
                for (i=0;i<NTEMP;i++)
                {
                  telog = temp_list[i];
                  
                  //printf("TEMP: %g | ", telog);

                  for (j=0;j<NMETAL;j++)
                  {
                    metal = pow(10.0, metal_list[j]);
                    xmetal = metal / (1.0 - metal) * metal_scale;
                    
                    //printf("METAL: %g | ", metal);

                    for (k=0;k<NHDEN;k++)
                    {
                        hden = hden_list[k];
                      
			/* initialize the code for this run */
			//cdInit();

			//cdTalk( true );

			/* if this is uncommented the calculation will not be done,
			 * but all parameters will be generated, as a quick way to see
			 * that grid is set up properly */
			/*cdNoExec( );*/

			/* cosmic background, microwave and hard parts */
			//nleft = cdRead( "background 0 "  );

			/* cosmic ray background, this is included because drives
			 * the chemistry in molecular gas */
			//nleft = cdRead( "cosmic ray background "  );

			/* this is a pure collisional model to turn off photoionization 
			   nleft = cdRead( "no photoionization "  );*/

			/* set a continuum shape, even though not used 
			   nleft = cdRead( "table agn "  );
			   nleft = cdRead( "phi(h) -4 "  );*/

			/* do only one zone */
			//nleft = cdRead( "stop zone 1 "  );

			/* set the hydrogen density */
			//sprintf(chLine,"hden %f ",hden);
			//nleft = cdRead( chLine  );

			/* sets the gas kinetic temperature */
			//sprintf(chLine,"constant temper %f ",telog);
			//nleft = cdRead( chLine  );

			/* identify sources of heating and cooling */
			//nleft = cdRead( "punch heating \"hazy_coolingcurve.het\" last no hash no clobber "  );
			//nleft = cdRead( "punch cooling \"hazy_coolingcurve.col\" last no hash no clobber "  );

			/* actually call the code */
			//lgBAD = cdDrive();

			/* get cooling for last zone */
			//cooling = cdCooling_last();

			/* want to print cooling over density squared */
			//cooling = cooling / pow(10.,hden*hden);
                  
                  ////////////////////////////////////////////////////////////
			/* initialize the code for this run */
			cdInit();

			cdTalk( true );
			
			nleft = cdRead("element helium abundance 0.08 linear");
			nleft = cdRead("element carbon abundance 2.45e-4 linear");
			nleft = cdRead("element oxygen abundance 4.90e-4 linear");
			
			nleft = cdRead("element nitrogen off");
			nleft = cdRead("element neon off");
			nleft = cdRead("element sodium off");
			nleft = cdRead("element magnesium off");
			nleft = cdRead("element aluminum off");
			nleft = cdRead("element silicon off");
			nleft = cdRead("element sulphur off");
			nleft = cdRead("element argon off");
			nleft = cdRead("element calcium off");
			nleft = cdRead("element iron off");
			nleft = cdRead("element nickel off");
			nleft = cdRead("element Lithium off");
			nleft = cdRead("element Beryllium off");
			nleft = cdRead("element Boron  off");
			nleft = cdRead("element Fluorine  off");
			nleft = cdRead("element Phosphor off");
			nleft = cdRead("element Chlorine off");
			nleft = cdRead("element Potassium off");
			nleft = cdRead("element Scandium  off");
			nleft = cdRead("element Titanium off");
			nleft = cdRead("element Vanadium off");
			nleft = cdRead("element Chromium off");
			nleft = cdRead("element Manganese off");
			nleft = cdRead("element Cobalt off");
			nleft = cdRead("element Copper off");
			nleft = cdRead("element Zinc  off");

			sprintf(chLine, "metals %f log", log10(xmetal));
			nleft = cdRead(chLine);

			nleft = cdRead( "cmb 18 "  );
			nleft = cdRead( "no molecules " );
			nleft = cdRead( "cosmic ray background -10 log" );
                        nleft = cdRead( "no photoionization " );
                        nleft = cdRead( "no compton effect " );
                        nleft = cdRead( "age 6 log years " );

			nleft = cdRead( "stop zone 1 "  );
			nleft = cdRead( "set dr 0 "  );
			
			sprintf(chLine,"hden %f ",hden);
			nleft = cdRead( chLine  );

			sprintf(chLine,"constant temper %f log",telog);
			nleft = cdRead( chLine  );

			nleft = cdRead( "punch heating \"metal_cooling.het\" last no hash no clobber "  );
			nleft = cdRead( "punch cooling \"metal_cooling.col\" last no hash no clobber "  );
			nleft = cdRead( "punch overview \"metal_cooling.ovr\" last no hash no clobber "  );
			
			lgBAD = cdDrive();
			if( lgBAD )
			{
				fprintf(stderr,"===== %d,%d,%d ===== problems!!\n", i, j, k);
			}
  			
  			
			eden = cdEDEN_last();

			cooling = cdCooling_last() / pow(10.0,hden+hden);
			
			heating = cdHeating_last() / pow(10.0,hden+hden);
                        
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

                        printf("%-10g | %-10g | %-10g = %-10g :Cool %-10g :Heat %-10g :H %-10g / %-10g :He %-10g / %-10g / %-10g\n", hden, metal, telog, eden, cooling, heating, h, hplus, he, heplus, heplusplus);
			
			cool_list[i][j][k] = cooling;
			heat_list[i][j][k] = heating;
			eden_list[i][j][k] = eden;

			h_list[i][j][k] = h;
			hplus_list[i][j][k] = hplus;
			he_list[i][j][k] = he;
			heplus_list[i][j][k] = heplus;
			heplusplus_list[i][j][k] = heplusplus;

                    }
		  }
		}

                for (k=0;k<NHDEN;k++)
                {
                  for (j=0;j<NMETAL;j++)
                  {
                    /* file to output short form of calculation's results */
                    sprintf(filename, "metal_cooling_%d_%d.txt", k+1, j+1);
                    printf("Outputting file %s.\n", filename);
                    
                    if( (ioRES = fopen(filename,"w")) == NULL )
                    {
                           printf(" could not open %s for writing.\n", filename);
                           cdEXIT(EXIT_FAILURE);
                    }
                    
                    for (i=0;i<NTEMP;i++)
                    {
                      fprintf(ioRES, "%10g %10g %10g %10g %10g %10g %10g %10g %10g %10g %10g\n", pow(10.0,hden_list[k]), pow(10.0,metal_list[j]), pow(10.0,temp_list[i]), eden_list[i][j][k], cool_list[i][j][k], heat_list[i][j][k], h_list[i][j][k], hplus_list[i][j][k], he_list[i][j][k], heplus_list[i][j][k], heplusplus_list[i][j][k]);
                    }

                    fclose(ioRES);
                  }
                }

		cdEXIT(lgBAD);
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Bailing out...\n" );
		cdPrepareExit();
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		cdPrepareExit();
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line(), e.comment() );
		cdPrepareExit();
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT )
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
			else if( e.sig() == SIGTERM )
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
			else if( e.sig() == SIGILL )
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
			else if( e.sig() == SIGFPE )
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
			else if( e.sig() == SIGSEGV )
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
#			endif
			else
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );

		}
		cdPrepareExit();
	}
#endif
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		cdPrepareExit();
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		cdPrepareExit();
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		cdPrepareExit();
	}



	return exit_status;
}

