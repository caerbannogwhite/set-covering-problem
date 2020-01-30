
#include "main.hpp"
#include "sc.hpp"

int main(int argc, char **argv)
{

	SCinstance inst;

	main_initialization(&inst);
	if (main_read_params(&inst, argc, argv))
	{
		cout << "Type --help to display available options." << endl;
		return 0;
	}

	SCMILPsolver(&inst);

	return 0;
}

// Initialize all configuration variables in SCinstance for computation
int main_initialization(SCinstance *inst)
{

	inst->nscrows = -1;
	inst->nsccols = -1;
	inst->costs = NULL;

	inst->MIP_nodesel = CPX_NODESEL_BESTBOUND;
	inst->MIP_varsel = CPX_VARSEL_DEFAULT;
	inst->MIP_reduce_prob = 1;

	inst->best_obj_val = 0;
	inst->obj_val = 0;

	inst->time_presolver = 0;
	inst->time_solver = 0;
	inst->time_total = 0;

	return 0;
}

// Scan arguments and set configuration paramentes
int main_read_params(SCinstance *inst, int argc, char **argv)
{
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
		("help", "produce help message")
		("presolver", po::value<string>(&inst->presolver)->default_value("none"), "set presolver, options: 'none' (default), 'cplex', 'dominance'")
		("solver", po::value<string>(&inst->solver)->default_value("cplex"), "set solver, options: 'cplex' (default), 'balasbcrule1', 'balasbcrule1-sparse', 'balasbcrule2', 'maxcol', 'maxcol2', 'maxcol-sparse'")
		("inputFile", po::value<string>(&inst->instance_name), "set input file path")
		("verbosity", po::value<int>(&inst->verbosity)->default_value(2), "set verbosity level")
		("numThreads", po::value<int>(&inst->num_threads)->default_value(1), "set number of threads")
		("timeLimit", po::value<double>(&inst->MIP_time_limit)->default_value(DBL_MAX), "set solver time limit (s)")
		("branchNum", po::value<int>(&inst->SC_BRANCHCB_NBRVARS)->default_value(2), "set number of branching variables")
		("maxBranch", po::value<int>(&inst->SC_BALAS_MAX_BRANCH)->default_value(8), "set maximum number of Balas nodes")
		("maxSingl", po::value<int>(&inst->SC_BALAS_MAX_SINGL)->default_value(2), "set maximum number of singletons in Balas")
		("cutsFactor", po::value<double>(&inst->MIP_cuts_factor)->default_value(-1.0), "set the MIP cuts factor level")
		("seed", po::value<int>(&inst->random_seed)->default_value(0), "set random seed")
		("debug", po::value<bool>(&inst->debug)->default_value(false), "set debug mode");

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			cout << desc << "\n";
			exit(0);
		}

		if (vm.count("presolver"))
		{
			cout << "Presolver set to " << vm["presolver"].as<string>() << ".\n";
		}
		if (vm.count("solver"))
		{
			cout << "Solver set to " << vm["solver"].as<string>() << ".\n";
		}
		if (vm.count("inputFile"))
		{
			cout << "Input file path set to " << vm["inputFile"].as<string>() << ".\n";
		}
		else
		{
			cout << "Input file path not set. Exiting.\n";
			return 1;
		}
		if (vm.count("verbosity"))
		{
			cout << "Verbosity level set to " << vm["verbosity"].as<int>() << ".\n";
		}
		if (vm.count("numThread"))
		{
			cout << "Number of threads set to " << vm["numThread"].as<int>() << ".\n";
		}
		if (vm.count("timeLimit"))
		{
			cout << "MIP time limit set to " << vm["timeLimit"].as<double>() << ".\n";
		}
		if (vm.count("branchNum"))
		{
			cout << "Branch number set to " << vm["branchNum"].as<int>() << ".\n";
		}
		if (vm.count("maxBranch"))
		{
			cout << "Maximum number of branches set to " << vm["maxBranch"].as<int>() << ".\n";
		}
		if (vm.count("maxSingl"))
		{
			cout << "Maximum number of singletons set to " << vm["maxSingl"].as<int>() << ".\n";
		}
		if (vm.count("cutsFactor"))
		{
			cout << "MIP cuts factor set to " << vm["cutsFactor"].as<double>() << ".\n";
		}
		if (vm.count("seed"))
		{
			cout << "Random seed set to " << vm["seed"].as<int>() << ".\n";
		}
		if (vm.count("debug"))
		{
			cout << "Debug value set to " << vm["debug"].as<bool>() << ".\n";
		}
	}
	catch (exception &e)
	{
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch (...)
	{
		cerr << "Exception of unknown type!\n";
		return 1;
	}

	return 0;
}
