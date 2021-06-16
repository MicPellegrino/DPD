#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "rng.h"
#include "ensemble.h"

#include <vector>
#include <cmath>

class DPD
{
private:
	double f0;
	double rc;
	int n_coll;	// Number of collisions
	int n_part;
	double m;
	double dt;
	double T0;

	double r=0;
	double f=0;
	double g=0;
	double v_relx=0;
	double v_rely=0;
	double v_relz=0;
	double dvx_temp = 0;
	double dvy_temp = 0;
	double dvz_temp = 0;
	std::vector<double> dvx;
	std::vector<double> dvy;
	std::vector<double> dvz;

public:
	DPD(double friction, double cutoff,  int nc, int np, double mass, double time_step, double temp):
		f0(friction), rc(cutoff), n_coll(nc), n_part(np), m(mass), dt(time_step), T0(temp), 
		dvx(np, 0.0), dvy(np, 0.0), dvz(np, 0.0) { }

	void dpd_step(Ensemble& ens, EngineWrapper& rng, double& e)
	{
		for (int i = 0; i<n_part; i++)
		{
			ens.vx[i] += dt*ens.fx[i]/m;
			ens.vy[i] += dt*ens.fy[i]/m;
			ens.vz[i] += dt*ens.fz[i]/m;
		}
		for (int i = 0; i<n_coll; i++)
		{
			int j, k;
			j = rng.random_index();
			do { k = rng.random_index(); } while(k==j);
			v_relx = ens.vx[j] - ens.vx[k];
			v_rely = ens.vy[j] - ens.vy[k];
			v_relz = ens.vz[j] - ens.vz[k];
			r = ens.pbc_dist(j, k);
			f = r<rc ? f0*(1.0-r/rc) : 0.0;
			g = sqrt(0.5*f*(2.0-f)*T0/m);
			dvx_temp = -f*v_relx + g*rng.gaussian(1.0);
			dvy_temp = -f*v_rely + g*rng.gaussian(1.0);
			dvz_temp = -f*v_relz + g*rng.gaussian(1.0);
			dvx[j] += dvx_temp;
			dvy[j] += dvy_temp;
			dvz[j] += dvz_temp;
			dvx[k] -= dvx_temp;
			dvy[k] -= dvy_temp;
			dvz[k] -= dvz_temp;
		}
		for (int i = 0; i<n_part; i++)
		{
			ens.px[i] += dt*(ens.vx[i]+0.5*dvx[i]);
			ens.py[i] += dt*(ens.vy[i]+0.5*dvy[i]);
			ens.pz[i] += dt*(ens.vz[i]+0.5*dvz[i]);
			ens.vx[i] += dvx[i];
			ens.vy[i] += dvy[i];
			ens.vz[i] += dvz[i];
			dvx[i] = 0.0;
			dvy[i] = 0.0;
			dvz[i] = 0.0;
			e += 0.5*m*( ens.vx[i]*ens.vx[i]+ens.vy[i]*ens.vy[i]+ens.vz[i]*ens.vz[i] );
		}
	}

};

#endif /* THERMOSTAT_H */