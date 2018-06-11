#include "includes/FEASolve.hpp"
#include "includes/Optimizer.hpp"
#include "includes/Utils.hpp"

int main(int argc, char * argv[])
{
	SimulatorApp sim;
	if (!sim.parseOptions(argc, argv))
		return -1;
	
	// FEASolve * solv = NULL;
	SimulatorApp::Options const opts_ = sim.options();
	
	if (opts_.to_pts)
	{
		int oneIndexed = opts_.oneIndexed ? 1 : 0;
		Utils::bouToPts(opts_.in_veg_file_path, opts_.support_vertices_files, opts_.in_simulation_mesh_file_path, oneIndexed);

		return 0;
	}

	std::cout << sim.optsToString();

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
		{
			solv.runImplicitNewmarkDense();
			cout << "Stress: " << solv.getStress() << endl;
		}

		string o_mesh_file_path = "mesh.deform.obj";
		solv.applyDeformationAndSaveSurfaceMesh(o_mesh_file_path);

		string o_deformations_file_path = "deformations.out";
		solv.saveDeformationsPerVertex(o_deformations_file_path);

		string o_area_def_feat = "area_deformations.feat";
		solv.m->saveLocalAreaDeformationAboutPointsAsFeature(o_area_def_feat);

		string o_point_displ_norm = "points_displc.feat";
		solv.m->savePointMovementMagnitudeAsFeature(o_point_displ_norm);

		string o_laplac_feat = "points_displc_laplace.feat";
		solv.m->saveLaplacianVectorNormAsFeature(o_laplac_feat);
	}
	
	return 0;
}
