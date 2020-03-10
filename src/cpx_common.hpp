#ifndef CPX_COMMON_H
#define CPX_COMMON_H

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

#include "params.hpp"

namespace po = boost::program_options;

// data structures
typedef struct SCinstance
{
    // input data
    std::string presolver;
    std::string solver;
    std::string inputFilePath;

    arma::vec dnsobj;
    arma::mat dnsmat;
    arma::vec sprobj;
    arma::sp_mat sprmat;

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

    double timeStart;
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

STATUS cpxcomm_initialization(SCinstance &inst);
STATUS cpxcomm_read_params(SCinstance &inst, int argc, char *argv[]);
STATUS cpxcomm_log(SCinstance &inst, int level, std::string msg);
STATUS cpxcomm_read_instance_dns(SCinstance &inst);
STATUS cpxcomm_read_instance_spr(SCinstance &inst);

#endif // CPX_COMMON_H
