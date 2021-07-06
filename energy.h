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
	int n_dump;
	std::vector<double> E_kin;
	std::vector<double> E_pot;
	std::vector<double> E_tot;
	std::vector<double> xcom;
	std::vector<double> ycom;
	std::vector<double> zcom;
	std::vector<double> temperature;

public:
	Energy(double time_step, int n): dt(time_step), n_dump(n) { }
	void append_energy(double ek, double ep, double et)
	{
		E_kin.push_back(ek);
		E_pot.push_back(ep);
		E_tot.push_back(et);
	}
	void append_com(double x, double y, double z)
	{
		xcom.push_back(x);
		ycom.push_back(y);
		zcom.push_back(z);
	}
	void append_temp(double T)
	{
		temperature.push_back(T);
	}
	void output_energy(std::string file_name)
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
				sprintf(c, "%.3f %.3f %.3f %.3f\n",  i*n_dump*dt, E_kin[i], E_pot[i], E_tot[i] );
				output_file << c;
			}		
			output_file.close();
  		}
	}
	void output_com(std::string file_name)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# t (1) Xcom (1) Ycom (1) Zcom (1)\n";
			output_file << "@    title \"Center of mass\"\n";
			output_file << "@    xaxis  label \"time [1]\"\n";
			output_file << "@    yaxis  label \"position [1]\"\n";
			output_file << "@TYPE xy\n";
 			output_file << "@ s0 legend \"X\"\n";
			output_file << "@ s1 legend \"Y\"\n";
			output_file << "@ s2 legend \"Z\"\n";		
			for (int i = 0; i<xcom.size(); i++)
			{
				sprintf(c, "%.3f %.3f %.3f %.3f\n",  i*dt, xcom[i], ycom[i], zcom[i] );
				output_file << c;
			}		
			output_file.close();
  		}
	}
	void output_temperature(std::string file_name)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "# t (1) temp\n";
			output_file << "@    title \"Temperature\"\n";
			output_file << "@    xaxis  label \"time [1]\"\n";
			output_file << "@    yaxis  label \"kT [1]\"\n";
			output_file << "@TYPE xy\n";		
			for (int i = 0; i<temperature.size(); i++)
			{
				sprintf(c, "%.3f %.3f\n",  i*n_dump*dt, temperature[i]);
				output_file << c;
			}		
			output_file.close();
  		}
	}

};

#endif /* ENERGY_H*/
