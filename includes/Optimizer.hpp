#ifndef __Optimizer__hpp__
#define __Optimizer__hpp__

#include "FEASolve.hpp"
#include "eigenmvn.h"

#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

class Optimizer
{
	private:
		FEASolve * solv;
		double * w;
		int dim;

	public:
		
		Optimizer(FEASolve * solv_) : solv(solv_) {};
		
		bool searchOptimum();

		// generate a random vector by sampling a multivariate gaussian, mean gets updated
		void generateRandomVector(Eigen::Matrix<double, Eigen::Dynamic, 1> & mean);
		
		// update the weight vector
		bool updateWeightVector();

		// show weight vector
		void showWeightVector() const;

};

#endif