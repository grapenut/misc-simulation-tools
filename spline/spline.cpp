#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "function.hh"

using namespace std;

int main(int argc, char **argv)
{
	int i, n;
	double X, Y;

        double x_max;
        double x_min;
	double y_cutoff = 1.0e-10;

        int num_col;
        int row_count = 100;

	// post-spline output data
        vector<double> X_array(row_count);
        map<int, vector<double> > Y_grid;

	// input data
        vector<double> tmp_X_array;
        map< int, vector<double> > tmp_Y_grid;

	if (argc < 4) {
		cerr << "usage: spline min max cols" << endl;
		return 0;
	}

	x_min = atof(argv[1]);
	x_max = atof(argv[2]);
	num_col = atoi(argv[3]);

        // load non-uniform data in X and Y
	// read at the end of the loop to set EOF on blank ending line
	// pre-read it before the loop because of this
	cin >> X;
        while( !cin.eof() )
        {       
                tmp_X_array.push_back(X);

                for(i = 0; i < num_col; i++)
                {
                        cin >> Y;
                        tmp_Y_grid[i].push_back(Y);
                }
		
		// read the next X value, or set EOF
                cin >> X;
        }
        
        // build the uniform X array
        for(i = 0; i < row_count; i++)
        {
                X_array[i] = x_min * exp(double(i) * log(x_max / x_min) / double(row_count - 1));
        }

	// build the uniform Y grid
        for (i = 0; i < num_col; i++)
        {
                const Spline splined_Y_value(tmp_X_array, tmp_Y_grid[i]);
       
                for(n = 0; n < row_count; n++)
                {
			if (splined_Y_value(X_array[n]) > y_cutoff)
			{
	                        Y_grid[i].push_back(splined_Y_value(X_array[n]));
			} else {
	                        Y_grid[i].push_back(0.0);
			}
                }
        }
       
        // output the uniform table
        for(i = 0; i < row_count; i++)
        {
		cout << setprecision(5) << scientific << X_array[i];  // output X index

                for(n = 0; n < num_col; n++)
                {
                        cout << " " << Y_grid[n][i];
                }
                
                cout << endl;
        } 
}
