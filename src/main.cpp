
#include "main.hpp"
#include "sc.hpp"

int main(int argc, char *argv[])
{

	SCinstance inst;

	main_initialization(inst);
	if (main_read_params(inst, argc, argv))
	{
		cout << "Type --help to display available options." << endl;
		return 0;
	}

	SCMILPsolver(&inst);

	return 0;
}

/** Initialize all configuration variables in SCinstance.
 * 
 * @param inst - SCinstance
 * @return a code
 */
int main_initialization(SCinstance &inst)
{
	inst.mipNodesel = CPX_NODESEL_BESTBOUND;
	inst.mipVarsel = CPX_VARSEL_DEFAULT;
	inst.mipReduceProb = 1;

	inst.bestObjVal = 0;
	inst.objVal = 0;

	inst.timePresolver = 0;
	inst.timeSolver = 0;
	inst.timeTotal = 0;

	return SC_SUCCESFULL;
}

/** Scan arguments and set configuration paramentes.
 * 
 * @param inst - SCinstance
 * @param argc - the length of the argv array
 * @param argv - the input arguments array
 * @return a code
 */
int main_read_params(SCinstance &inst, int argc, char *argv[])
{
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
		("help", "produce help message")
		("presolver", po::value<string>(&inst.presolver)->default_value("none"), "set presolver, options: 'none' (default), 'cplex', 'dominance'")
		("solver", po::value<string>(&inst.solver)->default_value("cplex"), "set solver, options: 'cplex' (default), 'balasbcrule1', 'balasbcrule1-sparse', 'balasbcrule2', 'maxcol', 'maxcol2', 'maxcol-sparse'")
		("inputFile", po::value<string>(&inst.inputFilePath), "set input file path")
		("verbosity", po::value<int>(&inst.verbosityLevel)->default_value(2), "set verbosity level")
		("numThreads", po::value<int>(&inst.numThreads)->default_value(1), "set number of threads")
		("timeLimit", po::value<double>(&inst.mipTimeLimit)->default_value(DBL_MAX), "set solver time limit (s)")
		("branchNum", po::value<int>(&inst.scBranchCallNumVars)->default_value(2), "set number of branching variables")
		("maxBranch", po::value<int>(&inst.scBalasMaxBranch)->default_value(8), "set maximum number of Balas nodes")
		("maxSingl", po::value<int>(&inst.scBalasMaxSingl)->default_value(2), "set maximum number of singletons in Balas")
		("cutsFactor", po::value<double>(&inst.mipCutsFactor)->default_value(-1.0), "set the MIP cuts factor level")
		("seed", po::value<int>(&inst.randomSeed)->default_value(0), "set random seed")
		("debug", po::value<bool>(&inst.debug)->default_value(false), "set debug mode");

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
			return SC_GENERIC_ERROR;
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
		return SC_GENERIC_ERROR;
	}
	catch (...)
	{
		cerr << "Exception of unknown type!\n";
		return SC_GENERIC_ERROR;
	}

	return SC_SUCCESFULL;
}

/**
 * Read a Set Covering Problem instance at inst.inputFile
 * and represent it as a dense matrix.
 *
 * @param inst - SCinstance
 * @returns a code
 */
int main_read_instance_dns(SCinstance &inst)
{
	int m, n, i, j, nz, col;
	std::ifstream fileHandler;
	fileHandler.open(inst.inputFilePath);

	// read number of rows and columns
	fileHandler >> m >> n;

	inst.obj = ublas::vector<double>(n);
	inst.dnsmat = ublas::matrix<double>(m, n);

	// read objective values
	for (auto it = inst.obj.begin(); it != inst.obj.end(); ++it)
	{
		fileHandler >> *it;
	}

	// read matrix ones
	for (i = 0; i < m; ++i)
	{
		fileHandler >> nz;
		for (j = 0; j < nz; ++j)
		{
			fileHandler >> col;
			inst.dnsmat(i, col) = 1.0;
		}
	}
	fileHandler.close();

	// DEBUG: print this instance
	std::cout << "obj = ";
	for (auto it = inst.obj.cbegin(); it != inst.obj.cend(); ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;

	std::cout << "mat = " << std::endl;
	for (auto it = inst.dnsmat.cbegin1(); it != inst.dnsmat.cend1(); ++it)
	{
		for (auto jt = it.cbegin(); jt != it.cend(); ++jt)
		{
			std::cout << *jt << " ";
		}
		std::cout << std::endl;
	}

	return SC_SUCCESFULL;
}

/**
 * Read a Set Covering Problem instance at inst.inputFile
 * and represent it as a sparse matrix.
 *
 * @param &inst - SCinstance
 * @returns code
 */
int main_read_instance_spr(SCinstance &inst)
{

	return SC_SUCCESFULL;
}
