#include "../includes/SimulatorApp.hpp"

// options struct constructor
SimulatorApp::Options::Options() {}

std::string
SimulatorApp::optsToString() const
{
	std::ostringstream oss;

	oss << "\nProgram Options"
		<< "\n=============== "
		<< "\n  steps = " << opts.in_steps
		<< "\n  timesteps = " << opts.in_timestep
		<< "\n  fast = " << opts.in_fast
		<< "\n  slow = " << opts.in_slow
		<< "\n  veg-file = " << opts.in_veg_file_path
		<< "\n  fixed-vertices-file = " << opts.in_fixed_vertices_file_path
		<< "\n  cubic-polynomial-file = " << opts.in_cubic_polynomial_file_path
		<< "\n  sim-mesh-file = " << opts.in_simulation_mesh_file_path
		<< "\n  render-mesh-file = " << opts.in_rendering_mesh_file_path
		<< "\n  modal-matrix-file = " << opts.in_modal_matrix_file_path
		<< "\n  output-dir = " << opts.output_dir
		<< "\n  support vertices directory = " << opts.in_support_vertices_directory
		<< "\n  one-indexed = " << opts.oneIndexed
		<< "\n  number of iterations for the optimizer: " << opts.num_samples_for_optimizer
		<< "\n";

		return oss.str();
}

bool
SimulatorApp::parseOptions(int argc, char * argv[])
{
  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i)  // omit the program path
    args.push_back(argv[i]);

  return parseOptions(args);
}

bool
SimulatorApp::parseOptions(std::vector<std::string> const & args)
{
	std::string conf_file;
	std::string usage = "usage: solve --timestep REAL --steps INT [args]";

	po::options_description visible("Allowed options");
	visible.add_options()
		("help,h", "produce help message")
		("version,v", "Print the program version")
		("conf", po::value<std::string>(&conf_file)->default_value("solve.conf"), "Configuration file (overridden by duplicate cmdline options)")
		("steps", po::value<int>(&opts.in_steps), "set number of iterations for the solver.")
		("timestep", po::value<double>(&opts.in_timestep), "Set the timestep per iteration of the solver.")
		("fast,f", "Run the fast reduced solver (get approximate values)")
		("slow,s", "Run the slow deformation solver (get accurate values)")
		("preprocess,p", "Preprocess the volumetric mesh to compute the reduced system.")
		("veg", po::value<std::string>(&opts.in_veg_file_path)->default_value(""), "Vega file (.veg)")
		("fixed", po::value<std::string>(&opts.in_fixed_vertices_file_path), "Fixed Vertices file (.bou)")
		("cub", po::value<std::string>(&opts.in_cubic_polynomial_file_path), "Cubic Polynomial file path (.cub)")
		("sim_mesh", po::value<std::string>(&opts.in_simulation_mesh_file_path), "Simulation mesh file path (.obj)")
		("render_mesh", po::value<std::string>(&opts.in_rendering_mesh_file_path), "Rendering mesh file path (.obj)")
		("modal_matrix", po::value<std::string>(&opts.in_modal_matrix_file_path), "Modal Matrix file path (.URendering.float)")
		("output_dir", po::value<std::string>(&opts.output_dir)->default_value("./new_output_dir"), "Output directory path to write the output files.")
		("s_verts", po::value<std::string>(&opts.in_support_vertices_directory),"path to the directory containing support vertices (turn on the spring penalty based support)")
		("one_indexed, oi", "The vertices are one indexed (default is 0 indexed).")
		("num_samples, ns", po::value<int>(&opts.num_samples_for_optimizer)->default_value(0), "Number of iterations you want the MH algorithm to make to find the optimal region set.")
	;

	po::options_description desc;
	desc.add(visible);

	// get the number of options
	int argc = args.size();
	if (argc < 1)
	{
		std::cout << usage << std::endl;
		std::cerr << visible;  // should be intercepted by out
		return false;
	}

	po::parsed_options cmdline_parsed = po::basic_command_line_parser<char>(args).options(desc).run();
	po::variables_map vm;
	po::store(cmdline_parsed, vm);

	// read the conf file if it is found.
	if (vm.count("conf") > 0 && Thea::FileSystem::fileExists(conf_file))
	{
		THEA_CONSOLE << "Reading options from config file:" << conf_file;

		std::ifstream conf_in(conf_file.c_str());
		po::parsed_options conf_file_parsed = po::parse_config_file(conf_in, desc);
		po::store(conf_file_parsed, vm);
	}

	po::notify(vm);

	bool quit = false;

	if (vm.count("version") > 0)
	{
		std::cout << "Solve version 1.0" << std::endl;
		std::cout << "Sourav Kumar Bose, 2018" << std::endl;
		quit = true;
	}

	if (vm.count("help") > 0)
	{
		if (quit) std::cout << "";
		std::cout << usage << std::endl;
		std::cout << visible;  // should be intercepted by THEA_CONSOLE
		return false;
	}

	if (vm.count("steps") != 1)
	{
		std::cout << usage << std::endl;
		std::cout << visible;  // should be intercepted by THEA_CONSOLE
		return false;	
	}

	if (vm.count("timestep") != 1)
	{
		std::cout << usage << std::endl;
		std::cout << visible;  // should be intercepted by THEA_CONSOLE
		return false;	
	}

	// read boolean
	opts.in_slow = (vm.count("slow") > 0);
	opts.in_fast = (vm.count("fast") > 0);
	opts.in_preprocess = (vm.count("preprocess") > 0);

	bool s_verts = (vm.count("s_verts") > 0);
	opts.oneIndexed = (vm.count("one_indexed") > 0);

	if (opts.in_slow)
	{
		if (s_verts)
		{
			std::cout << "Use support vertices only with fast solver.\n";
			std::cout << usage << std::endl;
			std::cout << visible << std::endl;
		}

		bool is_veg_present = (vm.count("veg") > 0);
		bool is_fixed_vertices_file_present = (vm.count("fixed") > 0);

		if (!is_veg_present || !is_fixed_vertices_file_present)
		{
			std::cout << "Please provide the veg file and fixed vertices file for the volumetric solver.\n";
			std::cout << usage << std::endl;
			std::cout << visible << std::endl;
			// std::cout << visible;
			return false;
		}
		else
		{
			opts.in_veg_file_path = Thea::FileSystem::resolve(opts.in_veg_file_path);
			opts.in_fixed_vertices_file_path = Thea::FileSystem::resolve(opts.in_fixed_vertices_file_path);
		}
	}
	else if (opts.in_fast)
	{
		bool is_cub_present = (vm.count("cub") > 0);
		bool is_sim_mesh_present = (vm.count("sim_mesh") > 0);
		bool is_modal_matrix_present = (vm.count("modal_matrix") > 0);

		if (!is_cub_present || !is_sim_mesh_present || !is_modal_matrix_present)
		{	
			if (!opts.in_preprocess)
			{
				std::cout << "Please provide the cub file, modal matrix and simulation mesh for the reduced system (or turn ON the preprocess flag).\n";
				std::cout << usage << std::endl;
				std::cout << visible << std::endl;
				// std::cout << visible;
				return false;		
			}
		}
		else
		{
			opts.in_cubic_polynomial_file_path = Thea::FileSystem::resolve(opts.in_cubic_polynomial_file_path);
			opts.in_simulation_mesh_file_path = Thea::FileSystem::resolve(opts.in_simulation_mesh_file_path);
			opts.in_rendering_mesh_file_path = opts.in_simulation_mesh_file_path;
			opts.in_modal_matrix_file_path = Thea::FileSystem::resolve(opts.in_modal_matrix_file_path);
		}

		if (s_verts)
		{
			opts.in_support_vertices_directory = Thea::FileSystem::resolve(opts.in_support_vertices_directory);
			int ret = Thea::FileSystem::getDirectoryContents(opts.in_support_vertices_directory, opts.support_vertices_files, -1, "*.bou");
			opts.num_support_regions = ret;
		}
		else
		{
			opts.num_support_regions = 0;
		}
	}
	else
	{
		std::cout << "Which solver do I choose? specify \'v\' or \'f\' as cmd arg.\n";
		std::cout << usage << std::endl;
		std::cout << visible << std::endl;
		return false;
	}

	opts.output_dir = Thea::FileSystem::resolve(opts.output_dir);

	// TODO: later include this to output the reduced rendering mesh
	// opts.in_rendering_mesh_file_path = Thea::FileSystem::resolve(opts.in_rendering_mesh_file_path);

	// if (quit)
	// 	return false;

	return true;
}
