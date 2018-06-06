#include "../includes/Optimizer.hpp"
#include <random>

bool Optimizer::searchOptimum()
{
	std::cout << "\nFinding optimal support regions.\n";

	w = solv->getWeightVector(dim);
	
	// get num samples from param options
	int num_samples = solv->getNumSamplesForOptimizer();

	// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> weights(dim, num_samples);

	// generateRandomVector(weights, num_samples);

	// vector to store the random update vector
	Eigen::Matrix<double, Eigen::Dynamic, 1> dw(dim);

	std::cout << std::endl;
	for (int sample = 0; sample < num_samples; sample++)
	{
		// solv->setWeightVector(w);
		// solv->runImplicitNewmarkDense();

		showWeightVector();
		generateRandomVector(dw);
		
		for (int i = 0; i < dim; i++)
		{

			w[i] += dw(i);
			
			// cap values to appropriate range
			// w[i] = min(1.0, w[i]);
			// w[i] = max(0.0, w[i]);
		}

		// solv->flushImplicitNewmarkDenseData();
	}
}

// wp represents the old value
bool Optimizer::updateWeightVector()
{
	return true;
}

void Optimizer::generateRandomVector(Eigen::Matrix<double, Eigen::Dynamic, 1> & mean)
{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> covar(dim, dim);

	// vector to store the differential weight vector.
	// Eigen::Matrix<double, Eigen::Dynamic, 1> dw(dim);

	covar.setIdentity(dim, dim);
	covar *= 0.01;

	// the "true" is to activate the cholesky factorization
	Eigen::EigenMultivariateNormal<double> gaussian_mult_generator(mean, covar, true);
	
	// update dw, here dw = mean
	mean = gaussian_mult_generator.samples(1);
}

void Optimizer::showWeightVector() const
{
	for (int i = 0; i < dim; i++)
	{
		std::cout << w[i] << " ";
	}
	std::cout << std::endl;
}