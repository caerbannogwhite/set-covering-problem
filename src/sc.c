
#include <cplex.h>
#include <string.h>
#include "sc.h"
#include "aux.h"
#include "callbacks.h"
#include "balas_dense.h"
#include "preprocessing.h"

/**	SCMILPsolver: this function reads the model, starts preprocessing routines
 *	and launches the solver selected by arguments (with the chosen configuration).
 *	The last step prints out some valuable information regarding the computation.
 */
int SCMILPsolver(SCinstance *inst) {
    char problem_name[100];
    char mode = 'w';
    int error;
    double timestamp_s, timestamp_e;

    strcpy(problem_name, "set_covering-");
    strcat(problem_name, inst->instance_name);


    /* ------------------------------- START CPLEX AND IMPORT THE MODEL --------------------------------------------- */
    CPXENVptr env = CPXopenCPLEX(&error);
    CPXLPptr lp = CPXcreateprob(env, &error, problem_name);

	CPXsetintparam(env, CPX_PARAM_THREADS, inst->num_threads);


    /* ------------------------------- PRESOLVE --------------------------------------------------------------------- */
	CPXENVptr pre_env = CPXopenCPLEX(&error);
	CPXLPptr pre_lp = CPXcreateprob(pre_env, &error, "presolver");
	CPXLPptr red_lp = CPXcreateprob(pre_env, &error, "red-lp");
	CPXreadcopyprob(pre_env, pre_lp, inst->instance_name, NULL);

	CPXsetintparam(pre_env, CPX_PARAM_THREADS, inst->num_threads);

	CPXgettime(env, &timestamp_s);

	// CPLEX PRESOLVER: set the maximum number of nodes to 1 and launch the solver
	// Then collect the reduced model and save it in lp file
    if (strcmp(inst->presolver, "cplex") == 0) {
		CPXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);
		CPXmipopt(pre_env, pre_lp);

		//CPXgetredlp(pre_env, pre_lp, &red_lp);

		//char cpx_problem_name[100];
		//strcpy(cpx_problem_name, inst->instance_name);
		//strcat(cpx_problem_name, "_cpx.lp");

		//error = CPXwriteprob(pre_env, red_lp, cpx_problem_name, NULL);
		//if (error) { printf("CPXwriteprob error: %d\n", error); }
		//printf("cols = %d\n", CPXgetnumcols(pre_env, red_lp));
    }

    // DOMINANCE PRESOLVER: scen all the columns to find dominated ones
    else if (strcmp(inst->presolver, "dominance") == 0) {

        CPXreadcopyprob(env, lp, inst->instance_name, NULL);

		SCdominancepresolver(inst, env, lp);

		double ub = 0.0;
		for (int j = CPXgetnumcols(env, lp); j >= 0; j--) {
			CPXgetub(env, lp, &ub, j, j);
			if (ub < 0.5) {
				CPXdelcols(env, lp, j, j);
			}
		}

		char dom_problem_name[100];
		strcpy(dom_problem_name, inst->instance_name);
		strcat(dom_problem_name, "_dom.lp");

		error = CPXwriteprob(env, lp, dom_problem_name, NULL);
		if (error) { printf("CPXwriteprob error: %d\n", error); }
		printf("cols = %d\n", CPXgetnumcols(env, lp));
    }

    // CPLEX+DOMINANCE PRESOLVER: launch cplex and stop the computation at the root node
    // collect the reduced model and launch the dominance presolver on it
    // REMOVE (?): the same result can be obtained by using the cplex presolver and providing
    // the reduced model as the input for the dominance routine
    /*else if (strcmp(inst->presolver, "cplex+dominance") == 0) {
		CPXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);

		CPXmipopt(pre_env, pre_lp);
		CPXgetredlp(pre_env, pre_lp, &red_lp);

		printf("cols = %d\n", CPXgetnumcols(pre_env, red_lp));

		error = CPXwriteprob(pre_env, red_lp, "reduced_model.lp", NULL);
		if (error) { printf("CPXwriteprob error: %d\n", error); }
		error = CPXreadcopyprob(env, lp, "reduced_model.lp", NULL);
		if (error) { printf("CPXreadcopyprob error: %d\n", error); }

		inst->num_selected_cols = CPXgetnumcols(env, lp);
		inst->selected_columns = (char *) malloc(inst->num_selected_cols);
		for (int j = 0; j < inst->num_selected_cols; ++j) inst->selected_columns[j] = 1;

		SCdominancepresolver(inst, env, lp);

		double ub = 0.0;
		for (int j = CPXgetnumcols(env, lp); j >= 0; j--) {
			CPXgetub(env, lp, &ub, j, j);
			if (ub < 0.5) {
				CPXdelcols(env, lp, j, j);
			}
		}

		char dom_problem_name[100];

		strcpy(dom_problem_name, inst->instance_name);
		strcat(dom_problem_name, "_dom.lp");

		error = CPXwriteprob(env, lp, dom_problem_name, NULL);
		if (error) { printf("CPXwriteprob error: %d\n", error); }

		printf("cols = %d\n", CPXgetnumcols(env, lp));
		printf("cols = %d\n", inst->num_selected_cols);
		free(inst->selected_columns);
    }*/

    // NO PRESOLVER
    else {
        if (inst->verbosity > 0) printf("No presolver selected.\n");
        CPXreadcopyprob(env, lp, inst->instance_name, NULL);
    }
	CPXgettime(env, &timestamp_e);

    CPXfreeprob(pre_env, &red_lp);
	CPXfreeprob(pre_env, &pre_lp);
	CPXcloseCPLEX(&pre_env);

    inst->time_presolver = timestamp_e - timestamp_s;


    /* ------------------------------- SET PARAMETERS AND SOLVE ----------------------------------------------------- */
    if (inst->verbosity > 0) {
    	CPXsetlogfilename(env, "log", &mode);
    }

    CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, inst->verbosity);
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->random_seed);
    CPXsetdblparam(env, CPX_PARAM_EPAGAP, 1e-4);
    CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);

	CPXsetintparam(env, CPX_PARAM_THREADS, inst->num_threads);

    CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);				// Disable problem reduction
    CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
    CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
    CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, inst->MIP_cuts_factor);

    CPXsetintparam(env, CPX_PARAM_NODESEL, 0 <= inst->MIP_nodesel && inst->MIP_nodesel <= 3 ? inst->MIP_nodesel : CPX_NODESEL_BESTBOUND);
	CPXsetintparam(env, CPX_PARAM_VARSEL, -1 <= inst->MIP_varsel && inst->MIP_varsel <= 4 ? inst->MIP_varsel : CPX_VARSEL_DEFAULT);

	//CPXsetintparam(env, CPX_PARAM_REDUCE, inst->MIP_reduce_prob ? CPX_ON : CPX_OFF);				// Disable problem reduction

    CPXgettime(env, &timestamp_s);
    if (strcmp(inst->solver, "cplex") == 0) {

    	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);
        CPXmipopt(env, lp);

	} else if (strcmp(inst->solver, "balascuts") == 0) {
		SCsolverbalascuts(inst, env, lp);
    } else if (strcmp(inst->solver, "balasbcrule1") == 0) {
        SCsolverbalasrule1(inst, env, lp);
	} else if (strcmp(inst->solver, "balasbcrule1-test") == 0) {
		SCsolverbalasrule1_test(inst, env, lp);
	} else if (strcmp(inst->solver, "balasbcrule1-sparse") == 0) {
		SCsolverbalasrule1_sparse(inst, env, lp);
	} else if (strcmp(inst->solver, "balasbcrule1maxcol-sparse") == 0) {
		SCsolverbalasrule1maxcol_sparse(inst, env, lp);
	} else if (strcmp(inst->solver, "balasbcrule2") == 0) {
		SCsolverbalasrule2(inst, env, lp);
	} else if (strcmp(inst->solver, "maxcol") == 0) {
		SCsolvermaxcol(inst, env, lp);
	} else if (strcmp(inst->solver, "maxcol2") == 0) {
		SCsolvermaxcol2(inst, env, lp);
	} else if (strcmp(inst->solver, "maxcol-sparse") == 0) {
		SCsolvermaxcol_sparse(inst, env, lp);
	} else if (strcmp(inst->solver, "maxcoldom") == 0) {
		SCsolvermaxcoldom(inst, env, lp);
	} else {
        //perror("solver unknown");
    }
    CPXgettime(env, &timestamp_e);
    inst->time_solver = timestamp_e - timestamp_s;


    /* ------------------------------- GET SOLUTION ----------------------------------------------------------------- */
    //ncols = CPXgetnumcols(env, lp);
    //inst->x_star = (double *) calloc(ncols, sizeof(double));
    //CPXgetx(env, lp, inst->x_star, 0, ncols - 1);
    CPXgetbestobjval(env, lp, &inst->best_obj_val);
    CPXgetobjval(env, lp, &inst->obj_val);


    /* ------------------------------- CLOSE CPLEX ------------------------------------------------------------------ */
    printf("\n");
    //printf("Build time       = %10.4lf\n", inst->time_build);
	printf("Presolver time   = %10.4lf\n", inst->time_presolver);
    printf("Solver time      = %10.4lf\n", inst->time_solver);
    //printf("Total time       = %10.4lf\n", inst->time_total);
    printf("Best obj val     = %10.4lf\n", inst->best_obj_val);
    printf("Obj val          = %10.4lf\n", inst->obj_val);
    printf("Node total       = %10d\n", CPXgetnodecnt(env, lp));
    printf("Node left        = %10d\n\n", CPXgetnodeleftcnt(env, lp));

    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}


int SCsolverbalascuts_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {
    int status;

    inst->nsccols = CPXgetnumcols(env, lp);
    inst->nscrows = CPXgetnumrows(env, lp);

    inst->balasnodecnt = 0;
    inst->cplexnodecnt = 0;

    inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
    status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
    if (status) { perror("Unable to get objective"); goto TERMINATE; }

    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

    inst->SC_BALAS_NODE_CUTS_NUM = 10;

    status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_sparse, inst);
    if (status) { perror("Unable to set user cuts callback"); goto TERMINATE; }

    CPXmipopt(env, lp);

    TERMINATE:
    free(inst->costs);

    return 1;
}


int SCsolverbalascuts(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {
	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	// usato per contare i tagli in questo caso
	inst->balasnodecnt = 0;

	inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	inst->SC_BALAS_NODE_CUTS_NUM = 10;

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_test, inst);
	if (status) { perror("Unable to set user cuts callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	//printf("balas cuts = %d\n", inst->balasnodecnt);

	TERMINATE:
	free(inst->costs);

	return 1;
}


int SCsolverbalasrule1(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

    int status;

    inst->nsccols = CPXgetnumcols(env, lp);
    inst->nscrows = CPXgetnumrows(env, lp);

    inst->balasnodecnt = 0;
    inst->cplexnodecnt = 0;

    inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
    status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
    if (status) { perror("Unable to get objective"); goto TERMINATE; }

    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

    status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1v1, inst);
    if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

    CPXmipopt(env, lp);

    //printf("balas node  = %d\n", inst->balasnodecnt);
    //printf("cplex node  = %d\n", inst->cplexnodecnt);

    TERMINATE:
    free(inst->costs);

    return 1;
}


int SCsolverbalasrule1_test(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	inst->balasnodecnt = 0;
	inst->cplexnodecnt = 0;

	inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_test, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	printf("balas node  = %d\n", inst->balasnodecnt);
	printf("cplex node  = %d\n", inst->cplexnodecnt);

	TERMINATE:
	free(inst->costs);

	return 1;
}


int SCsolverbalasrule1_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	inst->balasnodecnt = 0;
	inst->cplexnodecnt = 0;

	inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	/*CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	CPXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);*/

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_sparse, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst->balasnodecnt);
	//printf("cplex node  = %d\n", inst->cplexnodecnt);

	TERMINATE:
	free(inst->costs);

	return 1;
}


int SCsolverbalasrule1maxcol_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	inst->balasnodecnt = 0;
	inst->cplexnodecnt = 0;

	inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	/*CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	CPXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);*/

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1maxcol_sparse, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst->balasnodecnt);
	//printf("cplex node  = %d\n", inst->cplexnodecnt);

	TERMINATE:
	free(inst->costs);

	return 1;
}


int SCsolverbalasrule2(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);				// Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0);		// Turn off cuts generation
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule2, inst);
	if (status) { perror("Unable to set branch callback"); }

	CPXmipopt(env, lp);

	return 1;
}


int SCsolvermaxcol(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->costs = (double *) malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);				// Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0);		// Turn off cuts generation
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	TERMINATE:
	free(inst->costs);
	return 1;
}


int SCsolvermaxcol2(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->costs = (double *) malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	TERMINATE:
	free(inst->costs);
	return 1;
}


int SCsolvermaxcol_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->costs = (double *) malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol_sparse, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	TERMINATE:
	free(inst->costs);
	return 1;
}


int SCsolvermaxcoldom(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

	int status;

	inst->costs = (double *) malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status) { perror("Unable to get objective"); goto TERMINATE; }

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

	status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcoldom, inst);
	if (status) { perror("Unable to set branch callback"); goto TERMINATE; }

	CPXmipopt(env, lp);

	TERMINATE:
	free(inst->costs);
	return 1;
}
