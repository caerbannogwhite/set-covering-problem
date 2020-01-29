
#ifndef SC_COMMON_H
#define SC_COMMON_H

#include <boost/program_options.hpp>
#include <cplex.h>
#include <cpxconst.h>
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

using namespace std;
namespace po = boost::program_options;

#define DEBUG_VERBOSITY				0

#define SC_EPSILON_MED				1e-5
#define SC_EPSILON_SMALL			1e-10
#define BIG_M                   	1e20

#define SC_ERR_NOT_COVER			1
#define SC_ERR_NOT_PRIME_COVER		11
#define SC_ERR_NOT_DUAL_SOL			12
#define SC_ERR_PRIMAL				2
#define SC_ERR_DUAL					3
#define SC_ERR_DUAL_3				31
#define SC_ERR_CPXDUAL_0			32
#define SC_ERR_BALAS_COND_VIOLATED	40
#define SC_ERR_BALAS_BRANCH_RULE_1	41

// data structures
typedef struct SCinstance {

    // input data
    char presolver[100];
    char solver[100];
    char instance_name[100];
    int nscrows;
    int nsccols;
    double *costs;

    // CPLEX parameters
    int random_seed;
    int num_threads;
    int verbosity;
    int MIP_nodesel;
    int MIP_varsel;
    int MIP_reduce_prob;
    double MIP_time_limit;
    double MIP_cuts_factor;

    // other
    double best_obj_val;
    double obj_val;

    // balas branch callback
    int cplexnodecnt;
    int balasnodecnt;
    int SC_BRANCHCB_NBRVARS;
    int SC_BALAS_MAX_BRANCH;
    int SC_BALAS_MAX_SINGL;
    int SC_BALAS_NODE_CUTS_NUM;

    double time_presolver;
    double time_build;
    double time_solver;
    double time_total;

    char debug;

} SCinstance;


typedef struct SCnodedata {
    int ncols;
    int q;
    int p;
    int seqnum;
    char *qmat;
    int *rqbeg;
    int *rqind;
} SCnodedata;


typedef struct SCi2tuple {
    int a;
    int b;
} SCi2tuple;


typedef struct SCi3tuple {
	int a;
	int b;
	int c;
} SCi3tuple;


typedef struct SCidtuple {
	int a;
	double b;
} SCidtuple;


int SCint_cmp(const void *p, const void *q);
int SCi2tuple_cmpa(const void *p, const void *q);
int SCi2tuple_cmparev(const void *p, const void *q);
int SCi3tuple_cmpa(const void *p, const void *q);
int SCi3tuple_cmparev(const void *p, const void *q);
int SCidtuple_cmpb(const void *p, const void *q);

double func1(double c, int k);
double func2(double c, int k);
double func3(double c, int k);
double func4(double c, int k);
double func5(double c, int k);

#endif //SC_COMMON_H
