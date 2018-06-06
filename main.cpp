#include "includes/FEASolve.hpp"
#include "includes/Optimizer.hpp"

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
		FEASolve solv(&sim);
		// mass copied from the experimental version for comparison.
		solv.setMass(87.5668);
		solv.initImplicitNewmarkDense();
		
		if (opts_.num_samples_for_optimizer > 0)
		{
			// find the optimal solution
			Optimizer opt(&solv);
			opt.searchOptimum();
		}
		else
			solv.runImplicitNewmarkDense();

		// string o_mesh_file_path = "mesh.deform.obj";
		// solv.applyDeformationAndSaveSurfaceMesh(o_mesh_file_path);

		// string o_deformations_file_path = "deformations.out";
		// solv.saveDeformationsPerVertex(o_deformations_file_path);
	}
	
	return 0;
}
