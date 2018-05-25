#include "includes/FEASolve.hpp"

int main(int argc, char * argv[])
{
	SimulatorApp sim;
	if (!sim.parseOptions(argc, argv))
		return -1;
	
	std::cout << sim.optsToString();

	// FEASolve * solv = NULL;
	SimulatorApp::Options const opts_ = sim.options();

	if (opts_.in_slow)
	{
		FEASolve solv(&sim);
		solv.initImplicitBackwardEulerSparse();
		solv.runImplicitBackwardEulerSparse();
		solv.calculatePerVertexEnergy();

		cout << "\nTotal Absolute Strain Energy: " << solv.getAbsoluteStrainEnergy() << std::endl;
	}
	else if (opts_.in_fast)
	{
		cout << "\nTODO\n";
	}
	
	return 0;
}
