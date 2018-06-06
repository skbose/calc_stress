#ifndef __Optimizer__hpp__
#define __Optimizer__hpp__

#include "FEASolve.hpp"

class Optimizer
{
	private:
		FEASolve * solv;
		double * w;

	public:
		Optimizer(FEASolve * solv_) : solv(solv_) {};
		bool searchOptimum();
};

#endif