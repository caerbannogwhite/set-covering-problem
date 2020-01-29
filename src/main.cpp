
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "main.hpp"
#include "common.hpp"
#include "sc.hpp"

int main(int argc, char **argv) {

#if DEBUG_VERBOSITY
	char write = 'w';
	FILE *log = fopen("debug.log", &write);
	fclose(log);
#endif

	SCinstance inst;

    main_initialization(&inst);
    main_read_params(&inst, argc, argv);
    if (inst.verbosity > 0) main_print_info(&inst);

    SCMILPsolver(&inst);

    return 0;
}


// Initialize all configuration variables in SCinstance for computation
int main_initialization(SCinstance *inst) {

    stpcpy(inst->presolver, "none");
    stpcpy(inst->solver, "cplex");
    stpcpy(inst->instance_name, "none");
    inst->nscrows = -1;
    inst->nsccols = -1;
    inst->costs = NULL;

    inst->num_threads = 1;
    inst->random_seed = 0;
    inst->MIP_nodesel = CPX_NODESEL_BESTBOUND;
    inst->MIP_varsel = CPX_VARSEL_DEFAULT;
    inst->MIP_reduce_prob = 1;
    inst->MIP_time_limit = DBL_MAX;
    inst->MIP_cuts_factor = -1; // CPLEX default value: -1.0

    inst->best_obj_val = 0;
    inst->obj_val = 0;

    inst->SC_BRANCHCB_NBRVARS = 2;
    inst->SC_BALAS_MAX_BRANCH = 8;
    inst->SC_BALAS_MAX_SINGL = 2;

    inst->time_presolver = 0;
    inst->time_solver = 0;
    inst->time_total = 0;

    inst->debug = 0;
    return 0;
}


// Scan arguments and set configuration paramentes
int main_read_params(SCinstance *inst, int argc, char **argv) {

    int i = 1;
    char tmp_str[100];

    while (i < argc) {
        sscanf(argv[i], "--%s", tmp_str);

        if (strcmp(tmp_str, "presolver") == 0) { sscanf(argv[++i], "%s", tmp_str); strcpy(inst->presolver, tmp_str); ++i; }
        else if (strcmp(tmp_str, "solver") == 0) { sscanf(argv[++i], "%s", tmp_str); strcpy(inst->solver, tmp_str); ++i; }
        else if (strcmp(tmp_str, "instance") == 0) { sscanf(argv[++i], "%s", tmp_str); strcpy(inst->instance_name, tmp_str); ++i; }
        else if (strcmp(tmp_str, "verbosity") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->verbosity = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "num_threads") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->num_threads = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "random_seed") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->random_seed = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "MIP_time_limit") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->MIP_time_limit = strtod(tmp_str, NULL); ++i; }
        else if (strcmp(tmp_str, "MIP_cuts_factor") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->MIP_cuts_factor = strtod(tmp_str, NULL); ++i; }
        //else if (strcmp(tmp_str, "MIP_nodesel") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->MIP_nodesel = strtol(tmp_str, NULL, 10); ++i; }
        //else if (strcmp(tmp_str, "MIP_varsel") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->MIP_varsel = strtol(tmp_str, NULL, 10); ++i; }
        //else if (strcmp(tmp_str, "MIP_reduce_prob") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->MIP_reduce_prob = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "branchcb_nbrvars") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->SC_BRANCHCB_NBRVARS = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "balas_max_branch") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->SC_BALAS_MAX_BRANCH = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "balas_max_singl") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->SC_BALAS_MAX_SINGL = strtol(tmp_str, NULL, 10); ++i; }
        else if (strcmp(tmp_str, "debug") == 0) { sscanf(argv[++i], "%s", tmp_str); inst->debug = strtol(tmp_str, NULL, 10); ++i; }
        else { ++i; }
    }

    return 0;
}


// Print configuration parameters
int main_print_info(SCinstance *inst) {

    printf("\nSET COVERING - INFO\n");
    printf("presolver:        %s\n", inst->presolver);
    printf("solver:           %s\n", inst->solver);
    printf("SCinstance:       %s\n", inst->instance_name);
    printf("num threads:      %d\n", inst->num_threads);
    printf("random seed:      %d\n", inst->random_seed);
    printf("verbosity:        %d\n", inst->verbosity);
    printf("MIP time limit:   %lf\n", inst->MIP_time_limit);
	printf("MIP cuts factor:  %lf\n", inst->MIP_cuts_factor);
	//printf("MIP node sel:     %d\n", inst->MIP_nodesel);
	//printf("MIP var sel:      %d\n", inst->MIP_varsel);
	printf("MIP reduce prob:  %d\n", inst->MIP_reduce_prob);
    printf("branchcb nbrvars: %d\n", inst->SC_BRANCHCB_NBRVARS);
    printf("balas max branch: %d\n", inst->SC_BALAS_MAX_BRANCH);
	printf("balas max singl:  %d\n", inst->SC_BALAS_MAX_SINGL);
    printf("\n");

    return 0;
}
