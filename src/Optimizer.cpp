#include "../includes/Optimizer.hpp"
#include <random>

void Optimizer::fill(double * w, double val)
{
	for (int i = 0; i < dim; i++)
		w[i] = val;
}

bool Optimizer::searchOptimum()
{
	// random number seed
	srand(time(NULL));

	// optimization const
	double const C = 2;
	
	std::cout << "\nFinding optimal support regions.\n";

	w = solv->getWeightVector(dim);
	
	// get num samples from param options
	int num_samples = solv->getNumSamplesForOptimizer();

	// vector to store the random update vector
	Eigen::Matrix<double, Eigen::Dynamic, 1> wp(dim);
	wp.fill(0.0);

	// free later
	double * w_save = (double *) malloc(sizeof(double) * dim);

	double test = 0.0;
	for (int sample = 0; sample < num_samples; sample++)
	{

		generateRandomVector(wp, sample);
		
		for (int i = 0; i < dim; i++)
		{
			w_save[i] = w[i];
			w[i] = wp(i);
			
			// cap values to appropriate range
			w[i] = min(1.0, w[i]);
			w[i] = max(0.0, w[i]);
		}

		solv->setWeightVector(w);
		showWeightVector();
		
		solv->runImplicitNewmarkDense();

		double stress = solv->getStress();

		double rn = ((double) rand() / (double) RAND_MAX);
		double norm = solv->getNormOfWeightVector();

		stress += C * norm;
		double const base = 10;
		double const a = 8, b = 6;

		double cost = 1.0/(1+exp(-(stress/a - b)));
		
		// more the cost, lesser is the prob of getting selected.
		std::cout << "Random number: " << rn << ", Cost: " << cost << std::endl;
		// std::cout << "Stress: " << stress << std::endl;
		
		if (rn >= cost)
		{
			for (int i = 0; i < dim; i++)
				wp(i) = w_save[i];
		}
		else
		{
			for (int i = 0; i < dim; i++)
				wp(i) = w[i];
		}

		test += w[0];

		solv->flushImplicitNewmarkDenseData();
	}

	std::cout << test / num_samples << std::endl;
	// showWeightVector();
	free(w_save);
}

// wp represents the old value
bool Optimizer::updateWeightVector()
{
	return true;
}

void Optimizer::generateRandomVector(Eigen::Matrix<double, Eigen::Dynamic, 1> & mean, int seed)
{
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> covar(dim, dim);

	// vector to store the differential weight vector.
	// Eigen::Matrix<double, Eigen::Dynamic, 1> dw(dim);

	covar.setIdentity(dim, dim);
	covar *= 0.01;

	// the "true" is to activate the cholesky factorization
	Eigen::EigenMultivariateNormal<double> gaussian_mult_generator(mean, covar, true, time(NULL) + seed);
	
	// update dw, here dw = mean
	mean = gaussian_mult_generator.samples(1);
}

void Optimizer::showWeightVector() const
{
	std::cout << "weights: ( ";
	for (int i = 0; i < dim; i++)
	{
		std::cout << w[i] << " ";
	}
	std::cout <<")" << std::endl;
}