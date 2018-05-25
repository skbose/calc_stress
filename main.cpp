#include "includes/FEASolve.hpp"

int main(int argc, char * argv[])
{
	SimulatorApp sim;
	if (!sim.parseOptions(argc, argv))
		return -1;
	
	std::cout << sim.optsToString();

	FEASolve solv(&sim);
	
	return 0;
}
