
#include "sc.hpp"
#include "callbacks.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"
#include "preprocessing.hpp"

/**
 * In this function preprocessing and optimisation
 * routines are colled.
 * 
 * @param inst - SCP instance
 * @return a code
 */
int sc_solver(SCinstance &inst)
{
	int error;
	double timeStart, timeEnd;
	std::string thisFuncName = "sc_solver";
	std::string ext;
	std::string problemName = "set_covering-";
	std::vector<std::string> tokens;


	/////////////////////////	READ AND BUILD PROBLEM 		/////////////////////////////
	boost::split(tokens, inst.inputFilePath, boost::is_any_of("/."));
	ext = tokens[tokens.size() - 1];
	problemName += tokens[tokens.size() - 2];

	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXXcreateprob(env, &error, problemName[0]);

	CPXgettime(env, &timeStart);
	if (ext.compare("txt"))	// Read problem in raw text format
	{
		comm_read_instance_dns(inst);
		sc_build_raw2lp(inst, env, lp);
	} else
	if (ext.compare("lp")) // Import model in lp format
	{
		CPXXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);
		sc_build_lp2raw(inst, env, lp);
	} else
	{
		std::cerr << thisFuncName + " - unknown file extention found: " + ext << std::endl;
		std::cerr << thisFuncName + " - unable to read file." << std::endl;
		return SC_GENERIC_ERROR;
	}
	CPXgettime(env, &timeEnd);
	inst.timeBuild = timeEnd - timeStart;


	/////////////////////////	 PRESOLVE	/////////////////////////////////////////////
	/*CPXENVptr pre_env = CPXopenCPLEX(&error);
	CPXLPptr pre_lp = CPXcreateprob(pre_env, &error, "presolver");
	CPXLPptr red_lp = CPXcreateprob(pre_env, &error, "red-lp");
	CPXreadcopyprob(pre_env, pre_lp, &inst.inputFilePath, NULL);

	CPXsetintparam(pre_env, CPX_PARAM_THREADS, inst.numThreads);

	CPXgettime(env, &timeStart);

	// CPLEX PRESOLVER: set the maximum number of nodes to 1 and launch the solver
	// Then collect the reduced model and save it in lp file
	if (inst.presolver.compare("cplex") == 0)
	{
		CPXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);
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
		CPXsetintparam(pre_env, CPX_PARAM_NODELIM, 1);

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
		if (inst.verbosityLevel > 0)
			printf("No presolver selected.\n");
		CPXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);
	}
	CPXgettime(env, &timeEnd);

	CPXfreeprob(pre_env, &red_lp);
	CPXfreeprob(pre_env, &pre_lp);
	CPXcloseCPLEX(&pre_env);

	inst.timePresolver = timeEnd - timeStart;*/


	/////////////////////////	 SET PARAMETERS AND SOLVE	/////////////////////////////

	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, inst.verbosityLevel);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst.randomSeed);
	CPXsetdblparam(env, CPX_PARAM_EPAGAP, 1e-4);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);

	//CPXsetintparam(env, CPX_PARAM_NODESEL, 0 <= inst.mipNodesel && inst.mipNodesel <= 3 ? inst.mipNodesel : CPX_NODESEL_BESTBOUND);
	//CPXsetintparam(env, CPX_PARAM_VARSEL, -1 <= inst.mipVarsel && inst.mipVarsel <= 4 ? inst.mipVarsel : CPX_VARSEL_DEFAULT);

	CPXgettime(env, &timeStart);
	if (inst.solver.compare("cplex") == 0)
	{
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);
		CPXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);

		CPXmipopt(env, lp);
	}
	else if (inst.solver.compare("balasbcrule1") == 0)
	{
		sc_solver_balas_rule1(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-test") == 0)
	{
		sc_solver_balas_rule1_test(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-sparse") == 0)
	{
		sc_solver_balas_rule1_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule2") == 0)
	{
		sc_solver_balas_rule2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol") == 0)
	{
		sc_solver_maxcol(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol2") == 0)
	{
		sc_solver_maxcol2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol-sparse") == 0)
	{
		sc_solver_maxcol_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("maxcoldom") == 0)
	{
		sc_solver_maxcol_dom(inst, env, lp);
	}
	else
	{
		std::cerr << thisFuncName + " - solver unknown: " + inst.solver << std::endl;
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

int sc_solver_balas_cuts_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp)
{
    int status;

    /*inst.nsccols = CPXgetnumcols(env, lp);
    inst.nscrows = CPXgetnumrows(env, lp);

    inst.balasnodecnt = 0;
    inst.cplexnodecnt = 0;

    status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);

    CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

    inst.SC_BALAS_NODE_CUTS_NUM = 10;

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_sparse, inst);

    CPXmipopt(env, lp);*/

    return 1;
}

int sc_solver_balas_cuts(SCinstance &inst, IloEnv env, IloCplex cplex)
{
	int status;

	/*inst.nsccols = CPXgetnumcols(env, lp);
	inst.nscrows = CPXgetnumrows(env, lp);

	// usato per contare i tagli in questo caso
	inst.balasnodecnt = 0;

	status = CPXgetobj(env, lp, inst.costs, 0, inst.nsccols - 1);

	CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);

	status = CPXsetusercutcallbackfunc(env, SCcallbackbalasusercuts_test, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

int sc_solver_balas_rule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1v1, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_balas_rule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_test, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_balas_rule1_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1_sparse, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_balas_rule1_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule1maxcol_sparse, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_balas_rule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*status = CPXsetbranchcallbackfunc(env, SCcallbackbalasbranchrule2, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_maxcol(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;
	size_t i;
	double val;

	for (i = 0; i < CPXgetnumcols(env, lp); ++i)
	{
		status = CPXgetobj(env, lp, &val, i, i);
		inst.obj(i) = val;
	}

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0); // Turn off cuts generation
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit);
	CPXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);

	status = CPXsetbranchcallbackfunc(env, callbacks_branch_maxcol, &inst);

	CPXmipopt(env, lp);

	return SC_SUCCESFULL;
}

int sc_solver_maxcol2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*int status;

	status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

int sc_solver_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcol_sparse, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

int sc_solver_maxcol_dom(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;

	/*status = CPXsetbranchcallbackfunc(env, SCcallbackbranchmaxcoldom, inst);

	CPXmipopt(env, lp);*/

	return 1;
}

STATUS sc_build_lp2raw(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	size_t i;
	size_t j;
	double val;

	// Read objective function
	inst.obj = arma::vec(CPXgetnumcols(env, lp));
	for (j = 0; j < CPXgetnumcols(env, lp); ++j)
	{
		CPXgetobj(env, lp, &val, j, j);
		inst.obj(j) = val;
	}

	// Read matrix
	inst.dnsmat = arma::mat(CPXgetnumrows(env, lp), CPXgetnumcols(env, lp));
	for (i = 0; i < (size_t)CPXgetnumrows(env, lp); ++i)
	{
		for (j = 0; j < (size_t)CPXgetnumcols(env, lp); ++j)
		{
			CPXgetcoef(env, lp, i, j, &val);
			inst.dnsmat = val;
		}
	}

	return SC_SUCCESFULL;
}

STATUS sc_build_raw2lp(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	char binCol = 'B';
	char geSense = 'G';
	size_t i;
	size_t j;
	double rhs = 1.0;
	double val;

	// Set objective function on lp model
	for (j = 0; j < inst.obj.n_elem; ++j)
	{
		val = inst.obj(j);
		int CPXnewcols(env, lp, 1, &val, NULL, NULL, &binCol, NULL);
	}

	// Set matrix on lp model
	for (i = 0; i < inst.dnsmat.n_rows; ++i)
	{
		CPXnewrows(env,lp, 1, &rhs, &geSense, NULL, NULL);
		for (j = 0; j < inst.dnsmat.n_cols; ++j)
		{
			val = inst.obj(j);
			CPXchgcoef(env, lp, i, j, val);
		}
	}

	return SC_SUCCESFULL;
}
