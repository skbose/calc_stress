#ifndef __SimulatorApp__hpp__
#define __SimulatorApp__hpp__


#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include "Thea/FilePath.hpp"
#include "Thea/FileSystem.hpp"

namespace po = boost::program_options;

class SimulatorApp
{	
	public:
		struct Options
		{
			Options();

			std::string in_veg_file_path;
			std::string in_fixed_vertices_file_path;
			std::string in_support_vertices_directory;

			std::vector<std::string> support_vertices_files;

			// file paths for fast deformations
			std::string in_cubic_polynomial_file_path;
			std::string in_simulation_mesh_file_path;
			std::string in_rendering_mesh_file_path;
			std::string in_modal_matrix_file_path;

			int in_steps;
			double in_timestep;

			//filepath to read the output directory path
			std::string output_dir;

			// filepaths for deformations
			std::string out_mean_deformations_file_path;
			std::string out_max_deformations_file_path;
			std::string out_all_deformations_file_path;

			std::string out_deformed_mesh_file_path;
			// generate pts file
			std::string out_point_cloud_file_path;
			// generate the fixed vertices pts file
			std::string out_point_cloud_fixed_vertices_file_path;

			// slow or fast deformation
			bool in_slow;
			bool in_fast;
			bool in_preprocess;

			// misc
			int num_support_regions;
			bool oneIndexed;
		};

		/** Get the set of program options. */
    	Options const & options() const { return opts; }

    	/** Get a textual representation of all the program options. */
    	std::string optsToString() const;
		
		bool parseOptions(std::vector<std::string> const & args);
		bool parseOptions(int argc, char * argv[]);

		// declare the FEA Solver class as a friend of Simulator app so as to access the program options.
		friend class FEASolve;

		private:
	    	Options opts;
};

#endif