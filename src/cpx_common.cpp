
#include "cpx_common.hpp"

/**
 * Initialize all configuration variables in SCinstance.
 * 
 * @param inst - SCinstance
 * @return a status code
 */
STATUS cpxcomm_initialization(SCinstance &inst)
{
    inst.mipNodesel = CPX_NODESEL_BESTBOUND;
    inst.mipVarsel = CPX_VARSEL_DEFAULT;
    inst.mipReduceProb = 1;

    inst.bestObjVal = 0;
    inst.objVal = 0;

    inst.timePresolver = 0;
    inst.timeSolver = 0;
    inst.timeTotal = 0;

    inst.debug = false;

    return SC_SUCCESFULL;
}

/**
 * Scan arguments and set configuration paramentes.
 * 
 * @param inst - SCinstance
 * @param argc - the length of the argv array
 * @param argv - the input arguments array
 * @return a status code
 */
STATUS cpxcomm_read_params(SCinstance &inst, int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("presolver", po::value<std::string>(&inst.presolver)->default_value("none"), "set presolver, options: 'none' (default), 'cplex', 'dominance'")
        ("solver", po::value<std::string>(&inst.solver)->default_value("cplex"), "set solver, options: 'cplex' (default), 'balasbcrule1', 'balasbcrule1-sparse', 'balasbcrule2', 'maxcol', 'maxcol2', 'maxcol-sparse'")
        ("inputFile", po::value<std::string>(&inst.inputFilePath), "set input file path")
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
            std::cout << desc << std::endl;
            exit(0);
        }

        std::printf("####################          PARAMETERS          ####################\n\n");
        if (vm.count("presolver"))
        {
            std::printf("%40s %s\n", "Presolver set to ", &vm["presolver"].as<std::string>()[0]);
        }
        if (vm.count("solver"))
        {
            std::printf("%40s %s\n", "Solver set to ", &vm["solver"].as<std::string>()[0]);
        }
        if (vm.count("inputFile"))
        {
            std::printf("%40s %s\n", "Input file path set to ", &vm["inputFile"].as<std::string>()[0]);
        }
        else
        {
            std::cerr << "Input file path not set. Exiting.\n";
            return SC_GENERIC_ERROR;
        }
        if (vm.count("verbosity"))
        {
            std::printf("%40s %i\n", "Verbosity level set to ", vm["verbosity"].as<int>());
        }
        if (vm.count("numThread"))
        {
            std::printf("%40s %i\n", "Number of threads set to ", vm["numThread"].as<int>());
        }
        if (vm.count("timeLimit"))
        {
            std::printf("%40s %.4E\n", "MIP time limit set to ", vm["timeLimit"].as<double>());
        }
        if (vm.count("branchNum"))
        {
            std::printf("%40s %i\n", "Branch number set to ", vm["branchNum"].as<int>());
        }
        if (vm.count("maxBranch"))
        {
            std::printf("%40s %i\n", "Maximum number of branches set to ", vm["maxBranch"].as<int>());
        }
        if (vm.count("maxSingl"))
        {
            std::printf("%40s %i\n", "Maximum number of singletons set to ", vm["maxSingl"].as<int>());
        }
        if (vm.count("cutsFactor"))
        {
            std::printf("%40s %lf\n", "MIP cuts factor set to ", vm["cutsFactor"].as<double>());
        }
        if (vm.count("seed"))
        {
            std::printf("%40s %i\n", "Random seed set to ", vm["seed"].as<int>());
        }
        if (vm.count("debug"))
        {
            std::printf("%40s %i\n", "Debug value set to ", vm["debug"].as<bool>());
        }
        std::printf("\n");
    }
    catch (std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
        return SC_GENERIC_ERROR;
    }
    catch (...)
    {
        std::cerr << "Exception of unknown type!\n";
        return SC_GENERIC_ERROR;
    }

    return SC_SUCCESFULL;
}

/**
 * Print a log message based on the level of verbosity
 * set by the user.
 * 
 * @param inst - SCP instance
 * @param level - the verbosity level of the message (0 to 10)
 * @param msg - the message to print
 * @return a status code
 */
STATUS cpxcomm_log(SCinstance &inst, int level, std::string msg)
{
    if (inst.verbosityLevel > level)
    {
        std::cout << msg << std::endl;
    } 
    return SC_SUCCESFULL;
}

/**
 * Read a Set Covering Problem instance at inst.inputFile
 * and represent it as a dense matrix.
 *
 * @param inst - SCinstance
 * @returns a status code
 */
STATUS cpxcomm_read_instance_dns(SCinstance &inst)
{
    int m;
    int n;
    int i;
    int j;
    int nz;
    int col;
    std::ifstream fileHandler;
    fileHandler.open(inst.inputFilePath);

    // Read number of rows and columns
    fileHandler >> m >> n;

    inst.dnsobj = arma::vec(n);
    inst.dnsmat = arma::mat(m, n);

    // Read objective values
    for (auto it = inst.dnsobj.begin(); it != inst.dnsobj.end(); ++it)
    {
        fileHandler >> *it;
    }

    // Read matrix ones
    for (i = 0; i < m; ++i)
    {
        fileHandler >> nz;
        for (j = 0; j < nz; ++j)
        {
            fileHandler >> col;
            inst.dnsmat(i, col - 1) = 1.0;
        }
    }
    fileHandler.close();

    return SC_SUCCESFULL;
}

/**
 * Read a Set Covering Problem instance at inst.inputFile
 * and represent it as a sparse matrix.
 *
 * @param &inst - SCinstance
 * @returns a status code
 */
STATUS cpxcomm_read_instance_spr(SCinstance &inst)
{
    int m;
    int n;
    int i;
    int j;
    int cntNonZero;
    int nz;
    int col;
    double val;

    std::ifstream fileHandler;
    fileHandler.open(inst.inputFilePath);

    // Read number of rows and columns
    fileHandler >> m >> n;

    inst.dnsobj = arma::vec(n);

    // Read objective values
    for (auto it = inst.dnsobj.begin(); it != inst.dnsobj.end(); ++it)
    {
        fileHandler >> *it;
    }

    // Count the number of non-zero elements
    cntNonZero = 0;
    for (i = 0; i < m; ++i)
    {
        fileHandler >> nz;
        cntNonZero += nz;
        for (j = 0; j < nz; ++j)
        {
            fileHandler >> col;
        }
    }
    fileHandler.close();

    // Reopen file to read the matrix
    fileHandler.open(inst.inputFilePath);

    // Read number of rows and columns and skip all objective values
    fileHandler >> m >> n;
    for (auto it = inst.dnsobj.cbegin(); it != inst.dnsobj.cend(); ++it)
    {
        fileHandler >> val;
    }

    // Read the matrix
    std::unique_ptr<arma::umat> locationsPtr(new arma::umat(2, cntNonZero));
    std::unique_ptr<arma::vec> valuesPtr(new arma::vec(cntNonZero));
    cntNonZero = 0;
    for (i = 0; i < m; ++i)
    {
        fileHandler >> nz;
        for (j = 0; j < nz; ++j)
        {
            fileHandler >> col;
            (*locationsPtr)(0, cntNonZero) = i;
            (*locationsPtr)(1, cntNonZero) = col - 1;
            (*valuesPtr)(cntNonZero) = 1.0;
            ++cntNonZero;
        }
    }
    fileHandler.close();

    inst.sprmat = arma::sp_mat(*locationsPtr, *valuesPtr, m, n, true, false);

    return SC_SUCCESFULL;
}
