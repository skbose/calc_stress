#include "SimulatorApp.hpp"

int main(int argc, char * argv[])
{
	SimulatorApp sim;
	sim.parseOptions(argc, argv);
	std::cout << sim.optsToString();
	
	return 0;
}
