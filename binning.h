#ifndef BINNING_H
#define BINNING_H

#include "ensemble.h"

#include <vector>

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
			density[idx] += 1.0;
			velocity[idx] += ens.vz[i];
		}
	}

}

#endif /* BINNING_H */
