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
	int n_part;
	double Lx;
	double Ly;
	double Lz;
	double T0;
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
				else if (label=="n_part") 		iss >> label >> n_part;
				else if (label=="Lx") 			iss >> label >> Lx;
				else if (label=="Ly") 			iss >> label >> Ly;
				else if (label=="Lz") 			iss >> label >> Lz;
				else if (label=="T0") 			iss >> label >> T0;
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
				label = "";
    			}
    			input.close();
  		}
	}

	int get_seed(void) { return seed; }
	int get_n_part(void) { return n_part; }
	double get_Lx(void) { return Lx; }
	double get_Ly(void) { return Ly; }
	double get_Lz(void) { return Lz; }
	double get_T0(void) { return T0; }
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

	void print_info(void) const
	{
		std::cout << "seed = " << seed << std::endl;
		std::cout << "n_part = " << n_part << std::endl;
		std::cout << "Lx = " << Lx << std::endl;
		std::cout << "Ly = " << Ly << std::endl;
		std::cout << "Lz = " << Lz << std::endl;
		std::cout << "T0 = " << T0 << std::endl;
		std::cout << "epsilon = " << epsilon << std::endl;
		std::cout << "sigma = " << sigma << std::endl;
		std::cout << "n_steps = " << n_steps << std::endl;
		std::cout << "dt = " << dt << std::endl;
		std::cout << "mass = " << mass << std::endl;
		std::cout << "n_dump_energy = " << n_dump_energy << std::endl;
		std::cout << "n_dump_trajec = " << n_dump_trajec << std::endl;
		std::cout << "friction = " << friction << std::endl;
		std::cout << "cutoff = " << cutoff << std::endl;
		std::cout << "coll_frq = " << coll_frq << std::endl;
		std::cout << "coll_num = " << coll_num << std::endl;
		std::cout << "amplitude_x = " << ampx << std::endl;
		std::cout << "amplitude_y = " << ampy << std::endl;
		std::cout << "amplitude_z = " << ampz << std::endl;
		std::cout << "wave_number = " << wn << std::endl;
	}

};

#endif /* INPUT_PARAMETERS_H */
