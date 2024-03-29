#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "rng.h"
#include "ensemble.h"

#include <vector>
#include <cmath>

class Andersen
{

private:
	double T0;
	double m;
	int n_coll;
	int n_part;
	double sdev;

public:
	Andersen(double temp, double mass, int nc, int np):
		T0(temp), m(mass), n_coll(nc), n_part(np)
		{
			sdev = sqrt(T0/m);
		}

	void step(Ensemble& ens, EngineWrapper& rng, double& dek, double& dep)
	{
		for (int i = 0; i < std::ceil(0.5*n_coll*n_part); ++i)
		{
			int j = rng.random_index();
			ens.vx[j] = rng.gaussian(sdev);
			ens.vy[j] = rng.gaussian(sdev);
			ens.vz[j] = rng.gaussian(sdev);
		}
	}

};

class DPD
{

private:
	double f0;
	double rc;
	int n_coll;	// Number of collisions (now PER-PARTICLE)
	int n_part;
	double m;
	double dt;
	double T0;

	double r=0.0;
	double ex=0.0;
	double ey=0.0;
	double ez=0.0;
	double f=0.0;
	double g=0.0;
	double xi=0.0;
	double v_relx=0.0;
	double v_rely=0.0;
	double v_relz=0.0;
	double v_par=0.0;
	double dvx_temp = 0.0;
	double dvy_temp = 0.0;
	double dvz_temp = 0.0;

public:
	DPD(double friction, double cutoff,  int nc, int np, double mass, double time_step, double temp):
		f0(friction), rc(cutoff), n_coll(nc), n_part(np), m(mass), dt(time_step), T0(temp) { }

	void step(Ensemble& ens, EngineWrapper& rng, double& dek, double& dep)
	{
		for (int i = 0; i<std::ceil(0.5*n_coll*n_part); i++)
		{
			int j, k;
			j = rng.random_index();
			do { k = rng.random_index(); } while(k==j);
			v_relx = ens.vx[j] - ens.vx[k];
			v_rely = ens.vy[j] - ens.vy[k];
			v_relz = ens.vz[j] - ens.vz[k];
			ens.pbc_dist(j, k, ex, ey, ez, r);
			ex /= r;
			ey /= r;
			ez /= r;
			v_par = ( v_relx*ex + v_rely*ey + v_relz*ez );
			/* Let's try some different formula */
			f = r < rc ? f0*(1.0-r/rc) : 0.0;
			// f = f0 / ( 1.0 + ( r / (0.1*rc) ) );
			// f = r < rc ? f0 : 0.0;
			// f = f0;
			g = sqrt(2.0*f*(2.0-f)*T0/m);
			xi = -f*v_par + g*rng.gaussian(1.0);
			dvx_temp = xi*ex;
			dvy_temp = xi*ey;
			dvz_temp = xi*ez;
			// ??? Conserved kinetic energy ???
			dek += 0.25*( dvx_temp*dvx_temp + dvy_temp*dvy_temp + dvz_temp*dvz_temp ) +
				0.5*( dvx_temp*v_relx + dvy_temp*v_rely + dvz_temp*v_relz );
			// ??? Conserved potential energy ???
			dep += 0.25*dt*( dvx_temp*(ens.fx[j]-ens.fx[k]) 
					+ dvy_temp*(ens.fy[j]-ens.fy[k]) 
					+ dvz_temp*(ens.fz[j]-ens.fz[k]) );
			ens.vx[j] += 0.5*dvx_temp;
			ens.vy[j] += 0.5*dvy_temp;
			ens.vz[j] += 0.5*dvz_temp;
			ens.vx[k] -= 0.5*dvx_temp;
			ens.vy[k] -= 0.5*dvy_temp;
			ens.vz[k] -= 0.5*dvz_temp;
			ens.px[j] += 0.25*dt*dvx_temp;
			ens.py[j] += 0.25*dt*dvy_temp;
			ens.pz[j] += 0.25*dt*dvz_temp;
			ens.px[k] -= 0.25*dt*dvx_temp;
			ens.py[k] -= 0.25*dt*dvy_temp;
			ens.pz[k] -= 0.25*dt*dvz_temp;
		}
	}

	void half_step(Ensemble& ens, EngineWrapper& rng, double& dek)
	{
		for (int i = 0; i<std::ceil(0.25*n_coll*n_part); i++)
		{
			int j, k;
			j = rng.random_index();
			do { k = rng.random_index(); } while(k==j);
			v_relx = ens.vx[j] - ens.vx[k];
			v_rely = ens.vy[j] - ens.vy[k];
			v_relz = ens.vz[j] - ens.vz[k];
			ens.pbc_dist(j, k, ex, ey, ez, r);
			ex /= r;
			ey /= r;
			ez /= r;
			v_par = ( v_relx*ex + v_rely*ey + v_relz*ez );
			/* Let's try some different formula */
			f = r < rc ? f0*(1.0-r/rc) : 0.0;
			// f = f0 / ( 1.0 + ( r / (0.1*rc) ) );
			// f = r < rc ? f0 : 0.0;
			// f = f0;
			g = sqrt(2.0*f*(2.0-f)*T0/m);
			xi = -f*v_par + g*rng.gaussian(1.0);
			dvx_temp = xi*ex;
			dvy_temp = xi*ey;
			dvz_temp = xi*ez;
			// ??? Conserved kinetic energy ???
			dek += 0.25*( dvx_temp*dvx_temp + dvy_temp*dvy_temp + dvz_temp*dvz_temp ) +
				0.5*( dvx_temp*v_relx + dvy_temp*v_rely + dvz_temp*v_relz );
			ens.vx[j] += 0.5*dvx_temp;
			ens.vy[j] += 0.5*dvy_temp;
			ens.vz[j] += 0.5*dvz_temp;
			ens.vx[k] -= 0.5*dvx_temp;
			ens.vy[k] -= 0.5*dvy_temp;
			ens.vz[k] -= 0.5*dvz_temp;
		}
	}

};

#endif /* THERMOSTAT_H */
