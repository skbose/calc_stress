#include "../includes/Optimizer.hpp"
#include <random>

bool Optimizer::searchOptimum()
{
	int dim;

	w = solv->getWeightVector(dim);

	// array of gaussian generators
	std::default_random_engine generator;
	std::vector<std::normal_distribution<double>> dists(n);

	
}