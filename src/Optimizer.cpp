#include "../includes/Optimizer.hpp"
#include <random>

bool Optimizer::searchOptimum()
{
	std::cout << "\nFinding optimal support regions.\n";

	w = solv->getWeightVector(dim);
	
	// get num samples from param options
	int num_samples = solv->getNumSamplesForOptimizer();

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> weights(dim, num_samples);

	generateRandomVector(weights, num_samples);

	std::cout << std::endl;
	for (int sample = 0; sample < num_samples; sample++)
	{
		for (int i = 0; i < dim; i ++)
		{

			w[i] += weights(i, sample);
			
			// cap values to appropriate range
			w[i] = min(1.0, w[i]);
			w[i] = max(0.0, w[i]);
		}

		solv->setWeightVector(w);
		solv->runImplicitNewmarkDense();
		solv->flushImplicitNewmarkDenseData();
		showWeightVector();
	}
}

// wp represents the old value
bool Optimizer::updateWeightVector()
{
	return true;
}

void Optimizer::generateRandomVector(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & weights, int num_samples)
{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> covar(dim, dim);

	// vector to store the differential weight vector.
	// Eigen::Matrix<double, Eigen::Dynamic, 1> dw(dim);

	covar.setIdentity(dim, dim);
	covar *= 0.01;

	Eigen::Matrix<double, Eigen::Dynamic, 1> mean(dim);
	mean.fill(0.0);

	// the "true" is to activate the cholesky factorization
	Eigen::EigenMultivariateNormal<double> gaussian_mult_generator(mean, covar, true, time(NULL));
	weights = gaussian_mult_generator.samples(num_samples);
}

void Optimizer::showWeightVector() const
{
	for (int i = 0; i < dim; i++)
	{
		std::cout << w[i] << " ";
	}
	std::cout << std::endl;
}