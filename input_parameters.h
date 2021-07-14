#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

class InputParameters
{

private:
	int seed;
	std::string init_conf="";
	int n_part;
	double Lx;
	double Ly;
	double Lz;
	double T0;
	double Tref;
	double epsilon;
	double sigma;
	int n_steps;
	double dt;
	double mass;
	int n_dump_energy;
	int n_dump_trajec;
	double friction;
	double cutoff;
	int coll_frq;
	int coll_num;
	double ampx=0;
	double ampy=0;
	double ampz=0;
	double wn=0;
	int n_bins_v=1;
	int n_bins_g=1;
	int t_binning=1;
	double vx0=0.0;
	std::string thermostat="";
	std::string integrator="";

public:
	InputParameters(std::string file_name) 
	{
		std::string line;
		std::string label;
  		std::ifstream input(file_name);
  		if (input.is_open())
  		{
    			while ( std::getline (input, line) )
    			{
				std::istringstream iss(line);
				iss >> label;
				if (label=="seed") 			iss >> label >> seed;
				else if (label=="init_conf")		iss >> label >> init_conf;
				else if (label=="n_part") 		iss >> label >> n_part;
				else if (label=="Lx") 			iss >> label >> Lx;
				else if (label=="Ly") 			iss >> label >> Ly;
				else if (label=="Lz") 			iss >> label >> Lz;
				else if (label=="T0") 			iss >> label >> T0;
				else if (label=="Tref")			iss >> label >> Tref;
				else if (label=="epsilon") 		iss >> label >> epsilon;
				else if (label=="sigma") 		iss >> label >> sigma;
				else if (label=="n_steps") 		iss >> label >> n_steps;
				else if (label=="dt") 			iss >> label >> dt;
				else if (label=="mass") 		iss >> label >> mass;
				else if (label=="n_dump_energy") 	iss >> label >> n_dump_energy;
				else if (label=="n_dump_trajec") 	iss >> label >> n_dump_trajec;
				else if (label=="friction")		iss >> label >> friction;
				else if (label=="cutoff")		iss >> label >> cutoff;
				else if (label=="coll_frq")		iss >> label >> coll_frq;
				else if (label=="coll_num")		iss >> label >> coll_num;
				else if (label=="amplitude_x")		iss >> label >> ampx;
				else if (label=="amplitude_y")		iss >> label >> ampy;
				else if (label=="amplitude_z")		iss >> label >> ampz;
				else if (label=="wave_number")		iss >> label >> wn;
				else if (label=="n_bins_v")		iss >> label >> n_bins_v;
				else if (label=="n_bins_g")		iss >> label >> n_bins_g;
				else if (label=="t_binning")		iss >> label >> t_binning;
				else if (label=="drift_velocity")	iss >> label >> vx0;
				else if (label=="thermostat")		iss >> label >> thermostat;
				else if (label=="integrator")		iss >> label >> integrator;
				label = "";
    			}
    			input.close();
  		}
	}

	int get_seed(void) { return seed; }
	std::string get_init_conf(void) { return init_conf; }
	int get_n_part(void) { return n_part; }
	double get_Lx(void) { return Lx; }
	double get_Ly(void) { return Ly; }
	double get_Lz(void) { return Lz; }
	double get_T0(void) { return T0; }
	double get_Tref(void) { return Tref; }
	double get_epsilon(void) { return epsilon; }
	double get_sigma(void) { return sigma; }
	int get_n_steps(void) { return n_steps; }
	double get_dt(void) { return dt; }
	double get_mass(void) { return mass; }
	int get_n_dump_energy(void) { return n_dump_energy; }
	int get_n_dump_trajec(void) { return n_dump_trajec; }
	double get_friction(void) { return friction; }
	double get_cutoff(void) { return cutoff; }
	int get_coll_frq(void) { return coll_frq; }
	int get_coll_num(void) { return coll_num; }
	double get_amplitude_x(void) { return ampx; }
	double get_amplitude_y(void) { return ampy; }
	double get_amplitude_z(void) { return ampz; }
	double get_wave_number(void) { return wn; }
	int get_n_bins_v(void) { return n_bins_v; }
	int get_n_bins_g(void) { return n_bins_g; }
	int get_t_binning(void) { return t_binning; }
	double get_vx0(void) { return vx0; }
	std::string get_thermostat(void) { return thermostat; }
	std::string get_integrator(void) { return integrator; }

	void print_info(void) const
	{
		std::cout << "seed = " << seed << std::endl;
		std::cout << "init. configuration = " << init_conf << std::endl;
		std::cout << "no. part = " << n_part << std::endl;
		std::cout << "Lx = " << Lx << std::endl;
		std::cout << "Ly = " << Ly << std::endl;
		std::cout << "Lz = " << Lz << std::endl;
		std::cout << "init. temperature = " << T0 << std::endl;
		std::cout << "ref. temperature = " << T0 << std::endl;
		std::cout << "epsilon = " << epsilon << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
		std::cout << "no. steps = " << n_steps << std::endl;
		std::cout << "dt = " << dt << std::endl;
		std::cout << "mass = " << mass << std::endl;
		std::cout << "no. steps energy = " << n_dump_energy << std::endl;
		std::cout << "no. steps trajec = " << n_dump_trajec << std::endl;
		std::cout << "friction = " << friction << std::endl;
		std::cout << "dpd cutoff = " << cutoff << std::endl;
		std::cout << "coll. frequency = " << coll_frq << std::endl;
		std::cout << "coll. number = " << coll_num << std::endl;
		std::cout << "force amp. x = " << ampx << std::endl;
		std::cout << "force amp. y = " << ampy << std::endl;
		std::cout << "force amp. z  = " << ampz << std::endl;
		std::cout << "force wave no. = " << wn << std::endl;
		std::cout << "no. bins velocity = " << n_bins_v << std::endl;
		std::cout << "no. bins radial = "  << n_bins_g << std::endl;
		std::cout << "binning frame = " << t_binning << std::endl;
		std::cout << "drift velocity = " << vx0 << std::endl;
		std::cout << "themostat type = " << thermostat << std::endl;
		std::cout << "integrator type = " << integrator << std::endl;
	}

};

#endif /* INPUT_PARAMETERS_H */
