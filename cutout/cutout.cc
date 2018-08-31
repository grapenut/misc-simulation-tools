#include <string>
#include <cstdlib>
#include <iostream>
#include <strstream>
#include <iomanip>


using namespace std;

int main()
{
  const size_t varn = 5;

  const string prefix = "";
  const string suffix = "_unformatted.dat";
  const string variables[varn] = { "dens", "temp", "h", "hplu", "htwo" };

  const size_t index_start = 0;
  const size_t index_end = 0;

  const string vdf_file = "cutout.vdf";

  const size_t grid_size = 512;

  const size_t level = 4;

  ostrstream vdfcreate;

  vdfcreate << "vdfcreate" 
	    << " -dimension " << grid_size << "x" << grid_size << "x" << grid_size
	    << " -numts " << index_end - index_start + 1
	    << " -level " << level 
	    << " -vars3d ";

  for(size_t v = 0; v < varn; v ++)
    {
      if(v != 0) vdfcreate << ":";
      vdfcreate << variables[v];
    }
  
  vdfcreate << " " << vdf_file << ends;

  cout << "executing: " << vdfcreate.str() << endl;

  system(vdfcreate.str());

  const bool checkpoints = false;

  const size_t checkpoint_digits = 2;

  for(size_t i = index_start; i <= index_end; i ++)
    {
      for(size_t v = 0; v < varn; v ++)
	{
	  ostrstream file_name;
	  file_name << prefix 
		    << variables[v];

	  if(checkpoints) 
	    {
	      file_name << "_" 
			<< setw(checkpoint_digits) << setfill('0') << i;
	    }

	  file_name << suffix
		    << ends;
	  
	  ostrstream raw2vdf;

	  raw2vdf << "raw2vdf"
		  << " -dbl"
		  << " -ts " << i - index_start
		  << " -varname " << variables[v]
		  << " " << vdf_file
		  << " " << file_name.str()
		  << ends;

	  cout << "executing: " << raw2vdf.str() << endl;

	  system(raw2vdf.str());
	}
    }

  return 0;
}
