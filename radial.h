#ifndef RADIAL_H
#define RADIAL_H

#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include <stdio.h>

class RadialDistribution
{

private:
	int n_bins;
	double delta_r;
	std::vector<int> bins;
	std::vector<double> d;

public:

	RadialDistribution(int nb, double r_max):
		n_bins(nb), delta_r(r_max/nb), bins(nb, 0), d(nb, 0.0) 
		{
			std::cout << "r.d.f. bin size = " << delta_r << std::endl;
		}

	void add_to_bin(double r)
	{
		int idx = (int)(r/delta_r);
		assert( idx<n_bins && "Binning error: outside maximal minimal-image distance" );
		bins[idx] += 1;
	}
	
	void reset(void)
	{
		for (int i=0; i<n_bins; ++i)
			bins[i] = 0.0;
	}

	void avg_frame(int nf)
	{
		for (int i=0; i<n_bins; ++i)
			d[i] += bins[i]/nf;
		reset();
	}

	void avg_distribution(int nf)
	{
		for (int i=0; i<n_bins; ++i)
			d[i] /= nf;
	}

	void output_distribution(std::string file_name)
	{
		double r, s=0.0;
		for (int i=0; i<n_bins; ++i)
		{
			r = delta_r*(i+0.5);
			d[i] /= (r*r);
			s+=d[i];
		}
		s *= delta_r;
		for (int i=0; i<n_bins; ++i)
			d[i] /= s;
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# r (1) dist\n";
			output_file << "@    title \"Radial distribution function\"\n";
			output_file << "@    xaxis  label \"r [1]\"\n";
			output_file << "@    yaxis  label \"d [1]\"\n";
			output_file << "@TYPE xy\n";		
			for (int i = 0; i<n_bins; i++)
			{
				sprintf(c, "%.5f %.5f\n",  delta_r*(i+0.5), d[i]);
				output_file << c;
			}		
			output_file.close();
  		}
	}

};

#endif /* RADIAL_H */
