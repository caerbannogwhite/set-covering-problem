
#ifndef SC_COMMON_H
#define SC_COMMON_H


#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <unordered_set>

#include <armadillo>
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplex.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

#define DEBUG_VERBOSITY 0

#define SC_EPSILON_MED 1e-5
#define SC_EPSILON_SMALL 1e-12
#define BIG_M 1e20

typedef enum STATUS
{
    SC_SUCCESFULL,
    SC_GENERIC_ERROR,
    SC_SOLVER_NOT_FOUND,
    SC_ERR_NOT_COVER,
    SC_ERR_NOT_PRIME_COVER,
    SC_ERR_NOT_DUAL_SOL,
    SC_ERR_PRIMAL,
    SC_ERR_DUAL,
    SC_ERR_DUAL_3,
    SC_ERR_CPXDUAL_0,
    SC_ERR_BALAS_COND_VIOLATED,
    SC_ERR_BALAS_BRANCH_RULE_1,
};

// data structures
typedef struct SCinstance
{
    // input data
    string presolver;
    string solver;
    string inputFilePath;

    arma::vec obj;
    arma::mat dnsmat;

    // CPLEX parameters
    int randomSeed;
    int numThreads;
    int verbosityLevel;
    int mipNodesel;
    int mipVarsel;
    int mipReduceProb;
    double mipTimeLimit;
    double mipCutsFactor;

    // other
    double bestObjVal;
    double objVal;

    // balas branch callback
    int cplexCntNode;
    int balasCntCode;
    int scBranchCallNumVars;
    int scBalasMaxBranch;
    int scBalasMaxSingl;
    int scBalasNodeNumCuts;

    double timePresolver;
    double timeBuild;
    double timeSolver;
    double timeTotal;

    bool debug;
} SCinstance;

typedef struct SCnodedata
{
    int ncols;
    int q;
    int p;
    int seqnum;
    char *qmat;
    int *rqbeg;
    int *rqind;
} SCnodedata;

STATUS comm_initialization(SCinstance &inst);
STATUS comm_read_params(SCinstance &inst, int argc, char *argv[]);
STATUS comm_read_instance_dns(SCinstance &inst);
STATUS comm_read_instance_spr(SCinstance &inst);

#endif //SC_COMMON_H
