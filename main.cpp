#include "rng.h"
#include "forces.h"
#include "ensemble.h"
#include "time_marching.h"
#include "input_parameters.h"
#include "energy.h"
#include "thermostat.h"
#include "nonequilibrium.h"
#include "binning.h"

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
 double, double, double, double, double, double);

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
double Tref = ip.get_Tref();
int np = ip.get_n_part();
int n_steps_collisions = ip.get_coll_frq();
int n_dump_energy = ip.get_n_dump_energy();
int n_dump_trajec = ip.get_n_dump_trajec();
int t_binning = ip.get_t_binning();
int N = ip.get_n_steps();
double dt = ip.get_dt();
double m = ip.get_mass();
double vx0 = ip.get_vx0();

EngineWrapper rng(ip.get_seed(), np);

LennardJones forces(ip.get_epsilon(), ip.get_sigma());
CosineForcing forcing(ip.get_amplitude_x(), ip.get_wave_number());
// ConstantForcing forcing(ip.get_amplitude_x(), ip.get_amplitude_y(), ip.get_amplitude_z());
Ensemble atoms(np, Lx, Ly, Lz);

// Initial configuration
std::ifstream input_conf(ip.get_init_conf());
if (input_conf)
{
	input_conf.close();
	std::cout << "Reading an input configuration: " << ip.get_init_conf() << std::endl;
	atoms.read_input_conf(ip.get_init_conf());
}
else
{
	std::cout << "No input file detected" << std::endl; 
	init_ensemble(atoms, rng, sqrt(T0/m), Lx, Ly, Lz, 1.122462048309373*ip.get_sigma(), vx0);
}

// Testing output
std::string file_name("conf/conf00000.gro");
atoms.dump(file_name, Lx, Ly, Lz);

// Testing time marching
LeapFrog integrator(m, dt);
// DPD thermostat(ip.get_friction(), ip.get_cutoff(), ip.get_coll_num(), np, m, dt, Tref);
Andersen thermostat(Tref, m, ip.get_coll_num(), np);
int n_zeros = 5;
std::string label;
double Epot=0;
double Ekin=0;
double Etot=0;
// Conserved energy
double dEkin=0;
double dEpot=0;
double Econ=0;
double T=0;
Energy ener(dt, n_dump_energy);
double xcom, ycom, zcom;
double vcx,  vcy,  vcz;

Binning binning(ip.get_n_bins(), Lz);

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
	forcing.apply(atoms);
	
	for (int j = 0; j<np; j++)
	{
		integrator.advect(atoms.fx[j], atoms.px[j], atoms.vx[j], Ekin);
		integrator.advect(atoms.fy[j], atoms.py[j], atoms.vy[j], Ekin);
		integrator.advect(atoms.fz[j], atoms.pz[j], atoms.vz[j], Ekin);
	}

	Etot = Ekin + Epot;

	if (i%n_steps_collisions==0)
	{
		// thermostat.dpd_step(atoms, rng, dEkin, dEpot);
		thermostat.andersen_step(atoms, rng, dEkin, dEpot);
	}

	atoms.apply_pbc();

	// TO-DO: compute the whole work that the thermostat does on the system, 
	// both the kinetic and the potential contribution (since positions are
	// updated, too)
	Econ = Etot-dEkin+dEpot;

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
		// std::cout << "Econ = " << Econ <<std::endl;
		std::cout << "v_drift = [" << vcx << "," << vcy << "," << vcz << "]" << std::endl;
		T = atoms.temperature(m);
		std::cout << "T = " << T << std::endl;
		ener.append_energy(Ekin, Epot, Etot);
		ener.append_temp(T);
		atoms.com(xcom, ycom, zcom);
		ener.append_com(xcom, ycom, zcom);
		Ekin = 0.0;
		Epot = 0.0;
		Etot = 0.0;
		dEkin = 0.0;
		dEpot = 0.0;
		Econ = 0.0;
	}

	// Perform binning only after 1/2 the total no. of steps
	if ( i%t_binning==0 && i > 0.5*N )
		binning.add_frame(atoms);
		
}

binning.average_over_frames();

system("cat conf/*.gro > traj.gro");
system("rm conf/*.gro");
atoms.dump("confout.gro", Lx, Ly, Lz);
ener.output_energy("ener.xvg");
ener.output_com("com.xvg");
ener.output_temperature("temperature.xvg");
binning.output_density("density.xvg");
binning.output_velocity("velocity.xvg");

return 0;

}

/* ------------------------------------------------------------------- */
/* ----- FUNCTIONS --------------------------------------------------- */
/* ------------------------------------------------------------------- */

// TO-DO: change how the system is initialized in order to accomodate for a 'thin'
// simulation box (or a 'slice'); it does not work with the present configuration
void init_ensemble
(Ensemble& ens, EngineWrapper& rng, 
 double kep, double Lx, double Ly, double Lz, double sp, double vx0)
{
	// Cubic lattice
	std::cout << "Initialize in a cubic lattice with spacing sp = " << sp << std::endl;
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
	// FCC lattice
	/*
	std::cout << "Initialize in a FCC cubic lattice with spacing sp = " << sp << std::endl;
	int N = ens.n_particles();
	int N_min = (int)ceil(cbrt(N));
	int n = 0;
	double r[4][3] = { {.0, .0, .0}, {.5, .5, .0}, {.0, .5, .5}, {.5, .0, .5} }; 
	for (int i = 0; i<N_min && n<N; ++i)
	{
		for (int j = 0; j<N_min && n<N; ++j)
		{
			for (int k = 0; k<N_min && n<N; ++k)
			{
				for (int m = 0; m<4 && n<N; ++m)
				{
					ens.px[n] = (r[m][0]+i)*sp;
					ens.py[n] = (r[m][1]+j)*sp;
					ens.pz[n] = (r[m][2]+k)*sp;
					n++;
				}
			}
		}
	}
	*/
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
	for (int i=0; i<N; ++i)
		ens.vx[i] += vx0;
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
