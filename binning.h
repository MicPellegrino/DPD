#ifndef BINNING_H
#define BINNING_H

#include "ensemble.h"

#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cassert>

class Binning
{

private:
	// Convention: binning along z
	int Nz;
	double Lz;
	double hz;
	int frames;
	int idx;
	std::vector<double> density;
	std::vector<double> velocity;

public:
	Binning(int n, double l): 
		Nz(n), Lz(l), hz(l/(double)n), frames(0),
       		density(n, 0.0), velocity(n, 0.0) { }
	
	void add_frame(Ensemble& ens)
	{
		int np = ens.n_particles();
		for (int i = 0; i < np; ++i)
		{
			idx = (int)(ens.pz[i]/hz);
			assert( (idx<Nz) && "Atom position out of binning grid!" );
			density[idx] += 1.0;
			velocity[idx] += ens.vx[i];
		}
		frames++;
	}

	void average_over_frames(void)
	{
		for (int j = 0; j < Nz; j++)
		{
			density[j] /= frames;
			velocity[j] /= frames;
		}
	}

	void output_density(std::string file_name)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# density (1) position\n";
			output_file << "@    title \"Density profile\"\n";
			output_file << "@    xaxis  label \"density [1]\"\n";
			output_file << "@    yaxis  label \"position [1]\"\n";
			output_file << "@TYPE xy\n";		
			for (int i = 0; i<Nz; i++)
			{
				sprintf(c, "%.3f %.3f\n", density[i], hz*((double)i+0.5) );
				output_file << c;
			}		
			output_file.close();
  		}
	}

	void output_velocity(std::string file_name)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# velocity (1) position\n";
			output_file << "@    title \"Velocity profile\"\n";
			output_file << "@    xaxis  label \"velocity [1]\"\n";
			output_file << "@    yaxis  label \"position [1]\"\n";
			output_file << "@TYPE xy\n";		
			for (int i = 0; i<Nz; i++)
			{
				sprintf(c, "%.3f %.3f\n", velocity[i], hz*((double)i+0.5) );
				output_file << c;
			}		
			output_file.close();
  		}
	}

};

#endif /* BINNING_H */
