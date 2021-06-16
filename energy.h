#ifndef ENERGY_H
#define ENERGY_H

#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>

class Energy
{

private:
	double dt;
	std::vector<double> E_kin;
	std::vector<double> E_pot;
	std::vector<double> E_tot;

public:
	Energy(double time_step): dt(time_step) { }
	void append(double ek, double ep, double et)
	{
		E_kin.push_back(ek);
		E_pot.push_back(ep);
		E_tot.push_back(et);
	}
	void output_xvg(std::string file_name)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# t (1) Ekin (1) Epot (1) Etot (1)\n";
			output_file << "@    title \"Energies\"\n";
			output_file << "@    xaxis  label \"time [1]\"\n";
			output_file << "@    yaxis  label \"energy [1]\"\n";
			output_file << "@TYPE xy\n";
 			output_file << "@ s0 legend \"kinetic\"\n";
			output_file << "@ s1 legend \"potential\"\n";
			output_file << "@ s2 legend \"total\"\n";		
			for (int i = 0; i<E_kin.size(); i++)
			{
				sprintf(c, "%.3f %.3f %.3f %.3f\n",  i*dt, E_kin[i], E_pot[i], E_tot[i] );
				output_file << c;
			}		
			output_file.close();
  		}
		// ???	
		// delete[] c;
	}

};

#endif /* ENERGY_H*/
