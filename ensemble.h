#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>

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

	Ensemble(int n, double lx, double ly, double lz): 
		np(n), Lx(lx), Ly(ly), Lz(lz), 
		px(n, 0.0), py(n, 0.0), pz(n, 0.0), 
		vx(n, 0.0), vy(n, 0.0), vz(n, 0.0),
		fx(n, 0.0), fy(n, 0.0), fz(n, 0.0)
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
				// ???
				// delete[] c;
			}
			sprintf(c, "%10.5f%10.5f%10.5f\n", Lx, Ly, Lz);
			output_file << c;			
			output_file.close();
  		}
		// ???	
		// delete[] c;
	}

	void apply_pbc(void)
	{
		for (int j = 0; j<np; ++j)
		{
			px[j] = px[j] - Lx*std::floor(px[j]/Lx);
			py[j] = py[j] - Ly*std::floor(py[j]/Ly);
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

	double pbc_dist(int i, int j)
	{
		double dx = px[i] - px[j];
    		double dy = py[i] - py[j];
		double dz = pz[i] - pz[j];
    		while (dx < -0.5*Lx) { dx += Lx; }
    		while (dx >  0.5*Lx) { dx -= Lx; }
    		while (dy < -0.5*Ly) { dy += Ly; }
    		while (dy >  0.5*Ly) { dy -= Ly; }
		while (dz < -0.5*Lz) { dz += Lz; }
    		while (dz >  0.5*Lz) { dz -= Lz; }
    		return sqrt(dx*dx + dy*dy + dz*dz);
	}

	int n_particles(void) { return np; }

};

#endif /* ENSEMBLE_H */
