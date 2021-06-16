#ifndef FORCES_H
#define FORCES_H

#include <cmath>

class LennardJones
{

private:
	
	double epsilon;
    	double sigma;
	double lj_term;

public:

	LennardJones(double e, double s) : epsilon(e), sigma(s) { }
	
	double pair_energy(double r)
	{ 	
		lj_term = std::pow( sigma/r, 6 );
		return 4.0*epsilon*( lj_term*lj_term - lj_term );
	}

	double pair_force(double r)
	{ 	
		lj_term = std::pow( sigma/r, 6 );
		return 24.0*epsilon*( 2.0*lj_term*lj_term - lj_term )/r;
	}

};

#endif /* FORCES_H */
