#ifndef RNG_H
#define RNG_H

#include <random>
#include <math.h>

class EngineWrapper
{

private:
	
	// Mersenne twister engine
	int seed;
	int n_part;
	std::mt19937 eng;
	std::uniform_real_distribution<> unif;
	std::normal_distribution<> norm;
	std::uniform_int_distribution<> idx;

public:
	
	EngineWrapper(int s, int n): seed(s), n_part(n), eng(s), unif(0,1), norm(0,1), idx(0, n-1) {}

	double gaussian(double sigma) {	return sigma*norm(eng); }
	double uniform(double L) { return L*unif(eng); }
	int random_index(void) { return idx(eng); }

};

#endif /* RNG_H */
