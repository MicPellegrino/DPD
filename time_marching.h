#ifndef TIME_MARCHING_H
#define TIME_MARCHING_H

#include "rng.h"

class LeapFrog
{

private:
	double m;
	double dt;

public:
	LeapFrog(double mass, double time_step): m(mass), dt(time_step) { }

	void half_step_back(double& force, double& v, double& e)
	{
		v -= 0.5*dt*force/m;
		e += 0.5*m*v*v;
	}

	void advect(double& force, double& x, double& v, double& e)
	{
		v += 0.5*dt*force/m;
		e += 0.5*m*v*v;
		v += 0.5*dt*force/m;
		x += dt*v;
	}

};

class VelocityVerlet
{

private:
	double m;
	double dt;

public:
	VelocityVerlet(double mass, double time_step): m(mass), dt(time_step) { }

	void half_step_velocity(double& force, double& v)
	{
		v += 0.5*dt*force/m;
	}

	void full_step_position(double& x, double& v, double& e)
	{
		x += dt*v;
		e += 0.5*m*v*v;
	}

};

#endif /* TIME_MARCHING_H */
