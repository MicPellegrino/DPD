#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <cassert>

// To perform c.o.m. unwrapping
template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

class Ensemble
{

private:
	int np;
	double Lx;
	double Ly;
	double Lz;

public:
	std::vector<double> px;
	std::vector<double> py;
	std::vector<double> pz;
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> vz;
	std::vector<double> fx;
	std::vector<double> fy;
	std::vector<double> fz;
	std::vector<int> period_x;
	std::vector<int> period_y;
	std::vector<int> period_z;

	Ensemble(int n, double lx, double ly, double lz): 
		np(n), Lx(lx), Ly(ly), Lz(lz), 
		px(n, 0.0), py(n, 0.0), pz(n, 0.0), 
		vx(n, 0.0), vy(n, 0.0), vz(n, 0.0),
		fx(n, 0.0), fy(n, 0.0), fz(n, 0.0),
		period_x(n, 0), period_y(n, 0), period_z(n, 0)
		{ }
	
	int n_particles(void) const { return np; }

	void dump(std::string file_name, double Lx, double Ly, double Lz)
	{
		char c[256];
		std::ofstream output_file(file_name);
		if (output_file.is_open())
  		{
    			output_file << "Funny GROMACS header.\n";
    			output_file << np << "\n";		
			for (int i = 0; i<np; i++)
			{
				sprintf(c, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", 
					i+1, "LJS", "LJS", i+1, px[i], py[i], pz[i], vx[i], vy[i], vz[i] );
				output_file << c;
			}
			sprintf(c, "%10.5f%10.5f%10.5f\n", Lx, Ly, Lz);
			output_file << c;			
			output_file.close();
  		}
	}

	void apply_pbc(void)
	{
		for (int j = 0; j<np; ++j)
		{
			period_x[j] += std::floor(px[j]/Lx);
			px[j] = px[j] - Lx*std::floor(px[j]/Lx);
			period_y[j] += std::floor(py[j]/Ly);
			py[j] = py[j] - Ly*std::floor(py[j]/Ly);
			period_z[j] += std::floor(pz[j]/Lz);
			pz[j] = pz[j] - Lz*std::floor(pz[j]/Lz);
		}
	}

	void pbc_dist(int i, int j, double& dx, double& dy, double& dz, double& r)
	{
		dx = px[i] - px[j];
    		dy = py[i] - py[j];
		dz = pz[i] - pz[j];
    		while (dx < -0.5*Lx) { dx += Lx; }
    		while (dx >  0.5*Lx) { dx -= Lx; }
    		while (dy < -0.5*Ly) { dy += Ly; }
    		while (dy >  0.5*Ly) { dy -= Ly; }
		while (dz < -0.5*Lz) { dz += Lz; }
    		while (dz >  0.5*Lz) { dz -= Lz; }
    		r = sqrt(dx*dx + dy*dy + dz*dz);
	}

	void com(double& xc, double& yc, double& zc)
	{
		xc=0;
		yc=0;
		zc=0;
		for (int i = 0; i < np; i++)
		{
			xc += px[i]+Lx*period_x[i];
			yc += py[i]+Ly*period_y[i];
			zc += pz[i]+Lz*period_z[i];
		}
		xc/=np;
		yc/=np;
		zc/=np;
	}

	void drift(double& vcx, double& vcy, double& vcz)
	{
		vcx=0;
		vcy=0;
		vcz=0;
		for (int i = 0; i < np; i++)
		{
			vcx += vx[i];
			vcy += vy[i];
			vcz += vz[i];
		}
		vcx/=np;
		vcy/=np;
		vcz/=np;
	}

	double temperature(double m)
	{
		double vcx=0;
		double vcy=0;
		double vcz=0;
		double T;
		drift(vcx, vcy, vcz);
		for (int i = 0; i < np; i++)
		{
			T += (1.0/(3.0*np)) * m * (vx[i]-vcx)*(vx[i]-vcx);
			T += (1.0/(3.0*np)) * m * (vy[i]-vcy)*(vy[i]-vcy);
			T += (1.0/(3.0*np)) * m * (vz[i]-vcz)*(vz[i]-vcz);
		}
	}

	int n_particles(void) { return np; }

	void read_input_conf(std::string file_name)
	{
		// Reading .gro file
		std::string line;
		std::string label;
  		std::ifstream input(file_name);
		int n = 0;
		int dummy1, dummy4;
		std::string dummy2, dummy3;
		while ( std::getline (input, line) )
    		{
			std::istringstream iss(line);
			// Skip n==0: it's just the header
			if (n==1)
			{
				int n_atoms;
				iss >> n_atoms;
				assert( (n_atoms==np) && "Number of atoms does not match!" );
			}
			else if (n==np+3)
			{
				double lx, ly, lz;
				iss >> lx >> ly >> lz;
				assert( (Lx==lx)&&(Ly==ly)&&(Lz==lz) && "Box dimensions does not match!" );
			}
			else if (n>0)
			{
				iss >> dummy1 >> dummy2 >> dummy3 >> dummy4 
					>> px[n-2] >> py[n-2] >> pz[n-2]
					>> vx[n-2] >> vy[n-2] >> vz[n-2];
			}
			n++;
		}
		input.close();
	}

};

#endif /* ENSEMBLE_H */
