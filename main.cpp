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
double Lx = ip.get_Lx();
double Ly = ip.get_Ly();
double Lz = ip.get_Lz();
double T0 = ip.get_T0();
int np = ip.get_n_part();
int n_steps_collisions = ip.get_coll_frq();
int n_dump_energy = ip.get_n_dump_energy();
int n_dump_trajec = ip.get_n_dump_trajec();
int N = ip.get_n_steps();
double dt = ip.get_dt();
double m = ip.get_mass();

EngineWrapper rng(ip.get_seed(), np);

LennardJones forces(ip.get_epsilon(), ip.get_sigma());
Ensemble atoms(np, Lx, Ly, Lz);
init_ensemble(atoms, rng, sqrt(T0/m), Lx, Ly, Lz, 1.122462048309373*ip.get_sigma());

// Testing output
std::string file_name("conf/conf00000.gro");
atoms.dump(file_name, Lx, Ly, Lz);

// Testing time marching
LeapFrog integrator(m, dt);
DPD thermostat(ip.get_friction(), ip.get_cutoff(), ip.get_coll_num(), np, m, dt, T0);
int n_zeros = 5;
std::string label;
double Epot=0;
double Ekin=0;
// Difference w.r.t. conserved kinetic energy
double dEkin=0;
double Etot=0;
// Conserved energy
double Econ=0;
Energy ener(dt, n_dump_energy);
double xcom, ycom, zcom;
double vcx,  vcy,  vcz;

system("rm traj.gro");

update_forces(forces, atoms, Epot);
for (int j = 0; j<np; j++)
{
	integrator.half_step_back(atoms.fx[j], atoms.vx[j], Ekin);
	integrator.half_step_back(atoms.fy[j], atoms.vy[j], Ekin);
	integrator.half_step_back(atoms.fz[j], atoms.vz[j], Ekin);
}
atoms.apply_pbc();
Etot = Ekin+Epot;

for (int i = 1; i<N; i++)
{

	update_forces(forces, atoms, Epot);
	
	for (int j = 0; j<np; j++)
	{
		integrator.advect(atoms.fx[j], atoms.px[j], atoms.vx[j], Ekin);
		integrator.advect(atoms.fy[j], atoms.py[j], atoms.vy[j], Ekin);
		integrator.advect(atoms.fz[j], atoms.pz[j], atoms.vz[j], Ekin);
	}

	Etot = Ekin+Epot;

	if (i%n_steps_collisions==0)
		thermostat.dpd_step(atoms, rng, dEkin);

	atoms.apply_pbc();

	// Technically, one need also to re-compute the work because of new POSITIONS
	Econ = Etot-dEkin;

	if (i%n_dump_trajec==0)
	{
		label = std::to_string(i/n_dump_trajec);
		file_name = "conf/conf" + std::string(n_zeros - label.length(), '0') + label + ".gro";
		atoms.dump(file_name, Lx, Ly, Lz);
	}

	if (i%n_dump_energy==0)
	{
		Ekin /= n_dump_energy;
		Epot /= n_dump_energy;
		Etot /= n_dump_energy;
		Econ /= n_dump_energy;
		atoms.drift(vcx, vcy, vcz);
		std::cout << "---------------" <<std::endl;
		std::cout << "Iteration " << i <<std::endl;
		std::cout << "Ekin = " << Ekin <<std::endl;
		std::cout << "Epot = " << Epot <<std::endl;
		std::cout << "Etot = " << Etot <<std::endl;
		std::cout << "Econ = " << Econ <<std::endl;
		std::cout << "v_drift = [" << vcx << "," << vcy << "," << vcz << "]" << std::endl;
		std::cout << "T = " << (2.0/3.0)*(Ekin/np) << std::endl;
		ener.append_energy(Ekin, Epot, Etot);
		atoms.com(xcom, ycom, zcom);
		ener.append_com(xcom, ycom, zcom);
		Ekin = 0.0;
		dEkin = 0.0;
		Epot = 0.0;
		Etot = 0.0;
		Econ = 0.0;
	}

}


system("cat conf/*.gro > traj.gro");
system("rm conf/*.gro");
ener.output_energy("ener.xvg");
ener.output_com("com.xvg");

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
				ens.px[n] = j*sp;
				ens.py[n] = k*sp;
				ens.pz[n] = i*sp;
				n++;
			}
		}
	}
	for (int i = 0; i<N; ++i)
	{
		ens.vx[i] = rng.gaussian(kep);
		ens.vy[i] = rng.gaussian(kep);
		ens.vz[i] = rng.gaussian(kep);
	}
	double xcom, ycom, zcom;
	ens.com(xcom, ycom, zcom);
	double vcx, vcy, vcz;
	ens.drift(vcx, vcy, vcz);
	for (int i = 0; i<N; ++i)
	{
		ens.px[i] += 0.5*Lx-xcom;
		ens.py[i] += 0.5*Ly-ycom;
		ens.pz[i] += 0.5*Lz-zcom;
		ens.vx[i] -= vcx;
		ens.vy[i] -= vcy;
		ens.vz[i] -= vcz;
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
