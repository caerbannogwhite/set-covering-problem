#ifndef BALAS_COMMON_H
#define BALAS_COMMON_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <unordered_set>

#include <armadillo>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "params.hpp"

namespace po = boost::program_options;

// data structures
typedef struct BALSOLEnv
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
} BALSOLEnv;

STATUS balcomm_initialization(BALSOLEnv &inst);
STATUS balcomm_read_params(BALSOLEnv &inst, int argc, char *argv[]);
STATUS balcomm_log(BALSOLEnv &inst, int level, std::string msg);
STATUS balcomm_read_instance_dns(BALSOLEnv &inst);
STATUS balcomm_read_instance_spr(BALSOLEnv &inst);

#endif // BALAS_COMMON_H
