#ifndef NONEQUILIBRIUM_H
#define NONEQUILIBRIUM_H

#include "ensemble.h"

#include <cmath>

class Forcing
{

public:
	virtual void apply(Ensemble& ens) = 0;

};

class ConstantForcing : public Forcing
{

private:
	double ampx;
	double ampy;
	double ampz;

public:
	ConstantForcing(double ax, double ay, double az): 
		ampx(ax), ampy(ay), ampz(az) { }
	void apply(Ensemble& ens) override
	{
		int np = ens.n_particles();
		for (int i = 0; i < np; ++i)
		{
			ens.fx[i] += ampx;
			ens.fy[i] += ampy;
			ens.fz[i] += ampz;
		}
	}

};

class CosineForcing : public Forcing
{

private:
	// Convention: acceleration along X as funcion of Z
	double amp;
	double wn;

public:
	CosineForcing(double a, double w): amp(a), wn(w) { }
	void apply(Ensemble& ens) override
	{
		int np = ens.n_particles();
		for (int i = 0; i < np; ++i)
			ens.fx[i] += amp*std::cos(wn*ens.pz[i]);
	}

};

#endif /* NONEQUILIBRIUM_H */
