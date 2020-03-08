
#include "sc.hpp"
#include "callbacks.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"
#include "preprocessing.hpp"

/**	SCMILPsolver: this function reads the model, starts preprocessing routines
 *	and launches the solver selected by arguments (with the chosen configuration).
 *	The last step prints out some valuable information regarding the computation.
 */
STATUS SCMILPsolver(SCinstance &inst)
{
	STATUS status = SC_SUCCESFULL; 
	double timeStart, timeEnd;
	std::vector<std::string> tokens;

	boost::split(tokens, inst.inputFilePath, boost::is_any_of("/"));
	string problem_name = "set_covering-" + tokens[tokens.size() - 1];


	/* ---------------------- START CPLEX AND IMPORT THE MODEL ---------------- */
	IloEnv env;
	IloCplex cplex(env);
	IloModel model(env);


	/* ---------------------- PRESOLVE ---------------------------------------- */
	/*CPXENVptr pre_env = CPXopenCPLEX(&error);
	CPXLPptr pre_lp = CPXcreateprob(pre_env, &error, "presolver");
	CPXLPptr red_lp = CPXcreateprob(pre_env, &error, "red-lp");
	CPXreadcopyprob(pre_env, pre_lp, &inst->instance_name[0], NULL);

	cplex.setParam(IloCplex::Param::Threads, inst.numThreads);

	timeStart = cplex.getCplexTime();

	// CPLEX PRESOLVER: set the maximum number of nodes to 1 and launch the solver
	// Then collect the reduced model and save it in lp file
	if (inst.presolver.compare("cplex") == 0)
	{
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
	else if (inst.dnsmatpresolver.compare("dominance") == 0)
	{

		CPXreadcopyprob(env, lp, &inst->instance_name[0], NULL);

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

		string dom_problem_name;
		dom_problem_name = inst->instance_name + "_dom.lp";

		error = CPXwriteprob(env, lp, &dom_problem_name[0], NULL);
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
	else if (strcmp(inst->presolver, "cplex+dominance") == 0) {
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
    }

	// NO PRESOLVER
	else
	{
		if (inst->verbosity > 0)
			printf("No presolver selected.\n");
		CPXreadcopyprob(env, lp, &inst->instance_name[0], NULL);
	}
	CPXgettime(env, &timestamp_e);
	timeEnd = cplex.getCplexTime();

	CPXfreeprob(pre_env, &red_lp);
	CPXfreeprob(pre_env, &pre_lp);
	CPXcloseCPLEX(&pre_env);*/

	inst.timePresolver = timeEnd - timeStart;


	/* ---------------------- SET PARAMETERS AND SOLVE ------------------------ */
	if (inst.verbosityLevel > 0)
	{
		CPXsetlogfilename(env, "log", &mode);
	}

	cplex.setParam(IloCplex::Param::MIP::Display, inst.verbosityLevel);
	cplex.setParam(IloCplex::Param::RandomSeed, inst.randomSeed);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1E-4);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1E-4);

	cplex.setParam(IloCplex::Param::Threads, inst.numThreads);

	//CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	//CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 1);
	cplex.setParam(IloCplex::Param::MIP::Limits::CutsFactor, inst.mipCutsFactor);

	cplex.setParam(IloCplex::Param::MIP::Strategy::NodeSelect, 0 <= inst.mipNodesel && inst.mipNodesel <= 3 ? inst.mipNodesel : CPX_NODESEL_BESTBOUND);
	cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, -1 <= inst.mipVarsel && inst.mipVarsel <= 4 ? inst.mipVarsel : CPX_VARSEL_DEFAULT);

	timeStart = cplex.getCplexTime();
	if (inst.solver.compare("cplex") == 0)
	{
		cplex.setParam(IloCplex::Param::TimeLimit, inst.mipTimeLimit);
		cplex.solve();
	}
	else
	if (inst.solver.compare("balasbcrule1") == 0)
	{
		//SCsolverbalasrule1(inst, env, lp);
	}
	else if (inst->solver.compare("balasbcrule1-test") == 0)
	{
		//SCsolverbalasrule1_test(inst, env, lp);
	}
	else
	if (inst.solver.compare("balasbcrule1-sparse") == 0)
	{
		//SCsolverbalasrule1_sparse(inst, env, lp);
	}
	else
	if (inst.solver.compare("balasbcrule2") == 0)
	{
		//SCsolverbalasrule2(inst, env, lp);
	}
	else
	if (inst.solver.compare("maxcol") == 0)
	{
		//SCsolvermaxcol(inst, env, lp);
	}
	else
	if (inst.solver.compare("maxcol2") == 0)
	{
		//SCsolvermaxcol2(inst, env, lp);
	}
	else
	if (inst.solver.compare("maxcol-sparse") == 0)
	{
		//SCsolvermaxcol_sparse(inst, env, lp);
	}
	else
	if (inst->solver.compare("maxcoldom") == 0)
	{
		//SCsolvermaxcoldom(inst, env, lp);
	}
	else
	{
		std::cerr << "solver unknown." << std::endl;
		status = SC_SOLVER_NOT_FOUND;
	}
	timeEnd = cplex.getCplexTime();
	inst.timeSolver = timeEnd - timeStart;

	/* ---------------------- GET SOLUTION ------------------------------------ */
	CPXgetbestobjval(env, lp, &inst->best_obj_val);
	CPXgetobjval(env, lp, &inst->obj_val);


	/* ---------------------- CLOSE CPLEX ------------------------------------- */
	std::printf("\n");
	//std::printf("Build time       = %10.4lf\n", inst.timeBuild);
	std::printf("Presolver time   = %10.4lf\n", inst.timePresolver);
	std::printf("Solver time      = %10.4lf\n", inst.timeSolver);
	//std::printf("Total time       = %10.4lf\n", inst.timeTotal);
	std::printf("Best obj val     = %10.4lf\n", inst.bestObjVal);
	std::printf("Obj val          = %10.4lf\n", inst.objVal);
	std::printf("Node total       = %10d\n", CPXgetnodecnt(env, lp));
	std::printf("Node left        = %10d\n\n", CPXgetnodeleftcnt(env, lp));

	return status;
}

/*int SCsolverbalascuts_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {
    int status;

    inst->nsccols = CPXgetnumcols(env, lp);
    inst->nscrows = CPXgetnumrows(env, lp);

    inst->balasnodecnt = 0;
    inst->cplexnodecnt = 0;

    inst->costs = (double *) malloc(inst->nsccols * sizeof(double));
    status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
    if (status) { perror("Unable to get objective"); goto TERMINATE; }

    CPXXsetdblparam(env, CPX_PARAM_TILIM, inst->MIP_time_limit);

    inst->SC_BALAS_NODE_CUTS_NUM = 10;

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_sparse, inst);
	if (status) { perror("Unable to set user cuts callback"); goto TERMINATE; }

    CPXmipopt(env, lp);

    TERMINATE:
    free(inst->costs);

    return 1;
}

int SCsolverbalascuts(SCinstance &inst, IloEnv env, IloCplex cplex)
{
	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	// usato per contare i tagli in questo caso
	inst->balasnodecnt = 0;

	inst->costs = (double *)malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
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

	//printf("balas cuts = %d\n", inst->balasnodecnt);

	return SC_SUCCESFULL;
}

int SCsolverbalasrule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst->costs = (double *)malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
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

	//printf("balas node  = %d\n", inst->balasnodecnt);
	//printf("cplex node  = %d\n", inst->cplexnodecnt);

	return 1;
}

int SCsolverbalasrule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	inst->nsccols = CPXgetnumcols(env, lp);
	inst->nscrows = CPXgetnumrows(env, lp);

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst->costs = (double *)malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
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

	inst->costs = (double *)malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	//CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	//CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	//CPXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_sparse, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst->balasnodecnt);
	//printf("cplex node  = %d\n", inst->cplexnodecnt);

	return 1;
}

int SCsolverbalasrule1maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{

	int status;

	inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	inst->costs = (double *)malloc(inst->nsccols * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, inst->nsccols - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	//CPXsetintparam(env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);

	//CPXsetintparam(env, CPX_PARAM_STARTALG, CPX_ALG_DUAL);
	//CPXsetintparam(env, CPX_PARAM_SUBALG, CPX_ALG_DUAL);

	CPXXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1maxcol_sparse, inst);
	if (status)
	{
		perror("Unable to set branch callback");
		goto TERMINATE;
	}

	CPXmipopt(env, lp);

	//printf("balas node  = %d\n", inst->balasnodecnt);
	//printf("cplex node  = %d\n", inst->cplexnodecnt);

	return 1;
}

int SCsolverbalasrule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
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

	inst->costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
	if (status)
	{
		perror("Unable to get objective");
		goto TERMINATE;
	}

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
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

	inst->costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
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

	inst->costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
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

	inst->costs = (double *)malloc(CPXgetnumcols(env, lp) * sizeof(double));
	status = CPXgetobj(env, lp, inst->costs, 0, CPXgetnumcols(env, lp) - 1);
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
