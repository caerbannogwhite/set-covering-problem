
#include "sc.hpp"
#include "callbacks.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"
#include "preprocessing.hpp"

/**	SCMILPsolver: this function reads the model, starts preprocessing routines
 *	and launches the solver selected by arguments (with the chosen configuration).
 *	The last step prints out some valuable information regarding the computation.
 */

/**
 * In this function preprocessing and optimisation
 * routines are colled.
 * 
 * @param inst - SCP instance
 * @return a code
 */
int SCMILPsolver(SCinstance &inst)
{
	int error;
	double timeStart, timeEnd;
	std::string problemName = "set_covering-";
	std::vector<std::string> tokens;

	boost::split(tokens, inst.inputFilePath, boost::any_of("/."));
	problemName += tokens[tokens.size() - 2];


	/////////////////////////	 START CPLEX AND IMPORT THE MODEL	/////////////////////
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, problemName.c_str);

	CPXXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);


	/////////////////////////	 PRESOLVE	/////////////////////////////////////////////
	CPXENVptr pre_env = CPXopenCPLEX(&error);
	CPXLPptr pre_lp = CPXcreateprob(pre_env, &error, "presolver");
	CPXLPptr red_lp = CPXcreateprob(pre_env, &error, "red-lp");
	CPXreadcopyprob(pre_env, pre_lp, &inst.inputFilePath, NULL);

	CPXXsetintparam(pre_env, CPX_PARAM_THREADS, inst.numThreads);

	CPXgettime(env, &timeStart);

	/*// CPLEX PRESOLVER: set the maximum number of nodes to 1 and launch the solver
	// Then collect the reduced model and save it in lp file
	if (inst.presolver.compare("cplex") == 0)
	{
		CPXXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);
		CPXmipopt(pre_env, pre_lp);

		//CPXgetredlp(pre_env, pre_lp, &red_lp);

		//char cpx_problemName[100];
		//strcpy(cpx_problemName, inst.inputFilePath);
		//strcat(cpx_problemName, "_cpx.lp");

		//error = CPXwriteprob(pre_env, red_lp, cpx_problemName, NULL);
		//if (error) { printf("CPXwriteprob error: %d\n", error); }
		//printf("cols = %d\n", CPXgetnumcols(pre_env, red_lp));
	}

	// DOMINANCE PRESOLVER: scen all the columns to find dominated ones
	else if (inst.presolver.compare("dominance") == 0)
	{

		CPXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);

		SCdominancepresolver(inst, env, lp);

		double ub = 0.0;
		for (int j = CPXgetnumcols(env, lp); j >= 0; j--)
		{
			CPXgetub(env, lp, &ub, j, j);
			if (ub < 0.5)
			{
				CPXdelcols(env, lp, j, j);
			}
		}

		string dom_problemName;
		dom_problemName = inst.inputFilePath + "_dom.lp";

		error = CPXwriteprob(env, lp, &dom_problemName[0], NULL);
		if (error)
		{
			printf("CPXwriteprob error: %d\n", error);
		}
		printf("cols = %d\n", CPXgetnumcols(env, lp));
	}

	// CPLEX+DOMINANCE PRESOLVER: launch cplex and stop the computation at the root node
	// collect the reduced model and launch the dominance presolver on it
	// REMOVE (?): the same result can be obtained by using the cplex presolver and providing
	// the reduced model as the input for the dominance routine
	else if (strcmp(inst.presolver, "cplex+dominance") == 0) {
		CPXXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);

		CPXmipopt(pre_env, pre_lp);
		CPXgetredlp(pre_env, pre_lp, &red_lp);

		printf("cols = %d\n", CPXgetnumcols(pre_env, red_lp));

		error = CPXwriteprob(pre_env, red_lp, "reduced_model.lp", NULL);
		if (error) { printf("CPXwriteprob error: %d\n", error); }
		error = CPXreadcopyprob(env, lp, "reduced_model.lp", NULL);
		if (error) { printf("CPXreadcopyprob error: %d\n", error); }

		inst.num_selected_cols = CPXgetnumcols(env, lp);
		inst.selected_columns = (char *) malloc(inst.num_selected_cols);
		for (int j = 0; j < inst.num_selected_cols; ++j) inst.selected_columns[j] = 1;

		SCdominancepresolver(inst, env, lp);

		double ub = 0.0;
		for (int j = CPXgetnumcols(env, lp); j >= 0; j--) {
			CPXgetub(env, lp, &ub, j, j);
			if (ub < 0.5) {
				CPXdelcols(env, lp, j, j);
			}
		}

		char dom_problemName[100];

		strcpy(dom_problemName, inst.inputFilePath);
		strcat(dom_problemName, "_dom.lp");

		error = CPXwriteprob(env, lp, dom_problemName, NULL);
		if (error) { printf("CPXwriteprob error: %d\n", error); }

		printf("cols = %d\n", CPXgetnumcols(env, lp));
		printf("cols = %d\n", inst.num_selected_cols);
		free(inst.selected_columns);
    }

	// NO PRESOLVER
	else
	{
		if (inst.verbosity > 0)
			printf("No presolver selected.\n");
		CPXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);
	}
	CPXgettime(env, &timeEnd);

	CPXfreeprob(pre_env, &red_lp);
	CPXfreeprob(pre_env, &pre_lp);
	CPXcloseCPLEX(&pre_env);*/

	inst.timePresolver = timeEnd - timeStart;


	/////////////////////////	 SET PARAMETERS AND SOLVE	/////////////////////////////
	if (inst.verbosity > 0)
	{
		CPXsetlogfilename(env, "log", &mode);
	}

	CPXXsetintparam(env, CPX_PARAM_MIPDISPLAY, inst.verbosity);
	CPXXsetintparam(env, CPX_PARAM_RANDOMSEED, inst.random_seed);
	CPXsetdblparam(env, CPX_PARAM_EPAGAP, 1e-4);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);

	CPXXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);

	CPXXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, inst.MIP_cuts_factor);

	CPXXsetintparam(env, CPX_PARAM_NODESEL, 0 <= inst.MIP_nodesel && inst.MIP_nodesel <= 3 ? inst.MIP_nodesel : CPX_NODESEL_BESTBOUND);
	CPXXsetintparam(env, CPX_PARAM_VARSEL, -1 <= inst.MIP_varsel && inst.MIP_varsel <= 4 ? inst.MIP_varsel : CPX_VARSEL_DEFAULT);

	//CPXXsetintparam(env, CPX_PARAM_REDUCE, inst.MIP_reduce_prob ? CPX_ON : CPX_OFF);				// Disable problem reduction

	CPXgettime(env, &timeStart);
	if (inst.solver.compare("cplex") == 0)
	{

		CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);
		CPXmipopt(env, lp);

		/*} else if (strcmp(inst.solver, "balascuts") == 0) {
		SCsolverbalascuts(inst, env, lp);*/
	}
	else if (inst.solver.compare("balasbcrule1") == 0)
	{
		SCsolverbalasrule1(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-test") == 0)
	{
		SCsolverbalasrule1_test(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-sparse") == 0)
	{
		SCsolverbalasrule1_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule2") == 0)
	{
		SCsolverbalasrule2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol") == 0)
	{
		SCsolvermaxcol(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol2") == 0)
	{
		SCsolvermaxcol2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol-sparse") == 0)
	{
		SCsolvermaxcol_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("maxcoldom") == 0)
	{
		SCsolvermaxcoldom(inst, env, lp);
	}
	else
	{
		std::perror("solver unknown");
	}
	CPXgettime(env, &timeEnd);
	inst.timeSolver = timeEnd - timeStart;


	/////////////////////////	 GET SOLUTION 	/////////////////////////////////////////
	//ncols = CPXgetnumcols(env, lp);
	//inst.x_star = (double *) calloc(ncols, sizeof(double));
	//CPXgetx(env, lp, inst.x_star, 0, ncols - 1);
	CPXgetbestobjval(env, lp, &inst.bestObjVal);
	CPXgetobjval(env, lp, &inst.objVal);


	/////////////////////////	CLOSE CPLEX 	/////////////////////////////////////////
	std::printf("\n");
	std::printf("Build time       = %10.4lf\n", inst.timeBuild);
	std::printf("Presolver time   = %10.4lf\n", inst.timePresolver);
	//std::printf("Solver time      = %10.4lf\n", inst.timeSolver);
	std::printf("Total time       = %10.4lf\n", inst.timeTotal);
	std::printf("Best obj val     = %10.4lf\n", inst.bestObjVal);
	//std::printf("Obj val          = %10.4lf\n", inst.objVal);
	std::printf("Node total       = %10d\n", CPXgetnodecnt(env, lp));
	std::printf("Node left        = %10d\n\n", CPXgetnodeleftcnt(env, lp));

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return SC_SUCCESFULL;
}
/*int SCsolverbalascuts_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {
    int status;

    inst.nsccols = CPXgetnumcols(env, lp);
    inst.nscrows = CPXgetnumrows(env, lp);

    inst.balasnodecnt = 0;
    inst.cplexnodecnt = 0;

    inst.costs = (double *) malloc(inst.nsccols * sizeof(double));
    status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
    if (status) { perror("Unable to get objective"); goto TERMINATE; }

    CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

    inst.SC_BALAS_NODE_CUTS_NUM = 10;

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_sparse, inst);
	if (status) { perror("Unable to set user cuts callback"); goto TERMINATE; }

    CPXmipopt(env, lp);

    TERMINATE:
    free(inst.costs);

    return 1;
}

int SCsolverbalascuts(SCinstance &inst, IloEnv env, IloCplex cplex)
{
	int status;

	inst.nsccols = CPXgetnumcols(env, lp);
	inst.nscrows = CPXgetnumrows(env, lp);

	// usato per contare i tagli in questo caso
	inst.balasnodecnt = 0;

	inst.costs = (double *)malloc(inst.nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
	if (status)
	{
		std::cerr << "Unable to get objective" << std::endl;
		return SC_GENERIC_ERROR;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_test, inst);
	if (status)
	{
		std::cerr << "Unable to set user cuts callback" << std::endl;
		return SC_GENERIC_ERROR;
	}

	CPXmipopt(env, lp);

	//printf("balas cuts = %d\n", inst.balasnodecnt);

	return SC_SUCCESFULL;
}

int SCsolverbalasrule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst.costs = (double *)malloc(inst.nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
	if (status)
	{
		std::cerr << "Unable to get objective" << std::endl;
		return SC_GENERIC_ERROR;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1v1, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst.balasnodecnt);
	//printf("cplex node  = %d\n", inst.cplexnodecnt);

	return 1;
}

int SCsolverbalasrule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.nsccols = CPXgetnumcols(env, lp);
	inst.nscrows = CPXgetnumrows(env, lp);

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst.costs = (double *)malloc(inst.nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_test, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst.balasCntCode);
	//printf("cplex node  = %d\n", inst.cplexCntNode);

	return 1;
}

int SCsolverbalasrule1_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst.costs = (double *)malloc(inst.nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	//CPXXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	//CPXXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	//CPXXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_sparse, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst.balasnodecnt);
	//printf("cplex node  = %d\n", inst.cplexnodecnt);

	return 1;
}

int SCsolverbalasrule1maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{

	int status;

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst.costs = (double *)malloc(inst.nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	//CPXXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	//CPXXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	//CPXXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1maxcol_sparse, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst.balasnodecnt);
	//printf("cplex node  = %d\n", inst.cplexnodecnt);

	return 1;
}

int SCsolverbalasrule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	CPXXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0); // Turn off cuts generation
	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule2, inst);
	if (status)
	{
		perror("Unable to set branch callback");
	}

	CPXmipopt(env, lp);

	return 1;
}

int SCsolvermaxcol(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0); // Turn off cuts generation
	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	return 1;
}

int SCsolvermaxcol2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	return 1;
}

int SCsolvermaxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol_sparse, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	return 1;
}

int SCsolvermaxcoldom(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst.costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbranchmaxcoldom, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	return 1;
}*/
