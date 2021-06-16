#include "rng.h"
#include "forces.h"
#include "ensemble.h"
#include "time_marching.h"
#include "input_parameters.h"
#include "energy.h"
#include "thermostat.h"

#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>

// Just for dome static preliminary tests

#ifndef N_TEST_RNG
#define N_TEST_RNG 100
#endif

#ifndef SEED
#define SEED 1234
#endif

#ifndef N_PART
#define N_PART 20
#endif

void init_ensemble
(Ensemble&, EngineWrapper&, 
 double, double, double, double, double);

void update_forces
(LennardJones&, Ensemble&, double&);

/* ------------------------------------------------------------------- */
/* ----- MAIN -------------------------------------------------------- */
/* ------------------------------------------------------------------- */

int main ()
{

// Initialization

InputParameters ip("parameters.txt");
ip.print_info();

double Lx = ip.get_Lx();
double Ly = ip.get_Ly();
double Lz = ip.get_Lz();
double T0 = ip.get_T0();
int np = ip.get_n_part();
EngineWrapper rng(ip.get_seed(), np);
LennardJones forces(ip.get_epsilon(), ip.get_sigma());
Ensemble atoms(np, Lx, Ly, Lz);
init_ensemble(atoms, rng, T0, Lx, Ly, Lz, 1.122462048309373*ip.get_sigma());

// Testing output
std::string file_name("conf/conf00000.gro");
atoms.dump(file_name, Lx, Ly, Lz);

// Testing time marching
int N = ip.get_n_steps();
double dt = ip.get_dt();
double m = ip.get_mass();
LeapFrog integrator(m, dt);
DPD thermostat(ip.get_friction(), ip.get_cutoff(), ip.get_coll_num(), np, m, dt, T0);
Energy ener(dt);
int n_zeros = 5;
std::string label;
double Epot=0;
double Ekin=0;
double Etot=0;
int n_steps_collisions = ip.get_coll_num();
int n_dump_energy = ip.get_n_dump_energy();
int n_dump_trajec = ip.get_n_dump_trajec();

system("rm traj.gro");

update_forces(forces, atoms, Epot);
for (int j = 0; j<np; j++)
{
	integrator.half_step_back(atoms.fx[j], atoms.vx[j]);
	integrator.half_step_back(atoms.fy[j], atoms.vy[j]);
	integrator.half_step_back(atoms.fz[j], atoms.vz[j]);
}
atoms.apply_pbc();

for (int i = 0; i<N; i++)
{
	if (i%n_steps_collisions==0)
	{
		thermostat.dpd_step(atoms, rng, Ekin);
	}
	else
	{
		for (int j = 0; j<np; j++)
		{
			integrator.advect(atoms.fx[j], atoms.px[j], atoms.vx[j], Ekin);
			integrator.advect(atoms.fy[j], atoms.py[j], atoms.vy[j], Ekin);
			integrator.advect(atoms.fz[j], atoms.pz[j], atoms.vz[j], Ekin);
		}
	}
	atoms.apply_pbc();
	update_forces(forces, atoms, Epot);
	Etot = Ekin+Epot;

	if (i%n_dump_trajec==0)
	{
		label = std::to_string(i/n_dump_trajec);
		file_name = "conf/conf" + std::string(n_zeros - label.length(), '0') + label + ".gro";
		atoms.dump(file_name, Lx, Ly, Lz);
	}

	if (i%n_dump_energy==0)
	{
		std::cout << "Iteration " << i <<std::endl;
		std::cout << "Ekin = " << Ekin <<std::endl;
		std::cout << "Epot = " << Epot <<std::endl;
		std::cout << "Etot = " << Etot <<std::endl;
		ener.append(Ekin, Epot, Etot);
	}

	Ekin = 0;
	Epot = 0;
	Etot = 0;

}


system("cat conf/*.gro > traj.gro");
system("rm conf/*.gro");
ener.output_xvg("ener.xvg");

return 0;

}

/* ------------------------------------------------------------------- */
/* ----- FUNCTIONS --------------------------------------------------- */
/* ------------------------------------------------------------------- */

void init_ensemble
(Ensemble& ens, EngineWrapper& rng, 
 double kep, double Lx, double Ly, double Lz, double sp)
{
	std::cout << "Initialize in a cubic lattice with spacing sp = " << sp << std::endl;
	// Double-check this one! What is the correct kinetic energy?
	int N = ens.n_particles();
	int N_min = (int)ceil(cbrt(N));
	int n = 0;
	for (int i = 0; i<N_min && n<N; ++i)
	{
		for (int j = 0; j<N_min && n<N; ++j)
		{
			for (int k = 0; k<N_min && n<N; ++k)
			{
				ens.px[n] = i*sp;
				ens.py[n] = j*sp;
				ens.pz[n] = k*sp;
				n++;
			}
		}
	}
	for (int i = 0; i<N; ++i)
	{
		// Uniform spatial distribution
		// Does not work (unless one does a few energy minimization steps before)
		/*
		ens.px[i] = rng.uniform(Lx);
		ens.py[i] = rng.uniform(Ly);
		ens.pz[i] = rng.uniform(Lz);
		*/
		ens.vx[i] = rng.gaussian(kep);
		ens.vy[i] = rng.gaussian(kep);
		ens.vz[i] = rng.gaussian(kep);
	}
}

void update_forces
(LennardJones& lj, Ensemble& ens, double& e)
{
	int np = ens.n_particles();
	std::fill(ens.fx.begin(), ens.fx.end(), 0.0);
	std::fill(ens.fy.begin(), ens.fy.end(), 0.0);
	std::fill(ens.fz.begin(), ens.fz.end(), 0.0);
	double dx, dy, dz, r;
	double fij;
	for (int i=0; i<np; ++i)
	{
		for (int j=i+1; j<np; ++j)
		{
			ens.pbc_dist(i, j, dx, dy, dz, r);
			fij = lj.pair_force(r);
			e += lj.pair_energy(r);
			ens.fx[i] += fij * dx/r;
            		ens.fy[i] += fij * dy/r;
			ens.fz[i] += fij * dz/r;
            		ens.fx[j] -= fij * dx/r;
            		ens.fy[j] -= fij * dy/r;
            		ens.fz[j] -= fij * dz/r;
		}
	}
}
