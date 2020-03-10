
#include "cpx_solver.hpp"
#include "cpx_callbacks.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"

/**
 * In this function preprocessing and optimisation
 * routines are colled.
 * 
 * @param inst - SCP instance
 * @return a status code
 */
STATUS cpxsol(SCinstance &inst)
{
	int error;
	double timeCurr;
	double timeStart;
	double timeEnd;
	std::string thisFuncName = "cpxsol";
	std::string ext;
	std::string problemName = "set_covering-";
	std::vector<std::string> tokens;


	/////////////////////////	READ AND BUILD PROBLEM 		/////////////////////////////

	boost::split(tokens, inst.inputFilePath, boost::is_any_of("/."));
	ext = tokens[tokens.size() - 1];
	problemName += tokens[tokens.size() - 2];

	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXXcreateprob(env, &error, &problemName[0]);

	CPXgettime(env, &inst.timeStart); // Start counting time from this moment
	CPXgettime(env, &timeStart);
	if (ext.compare("txt") == 0)	// Read problem in raw text format
	{
		cpxcomm_read_instance_dns(inst);
		cpxsol_build_raw2lp(inst, env, lp);
	} else
	if (ext.compare("lp") == 0) // Import model in lp format
	{
		CPXXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);
		cpxsol_build_lp2raw(inst, env, lp);
	} else
	{
		std::cerr << thisFuncName + " - unknown file extention found: " + ext << std::endl;
		std::cerr << thisFuncName + " - unable to read file." << std::endl;
		return SC_GENERIC_ERROR;
	}
	CPXgettime(env, &timeEnd);
	inst.timeBuild = timeEnd - timeStart;


	/////////////////////////	 PRESOLVE	/////////////////////////////////////////////

	CPXENVptr envPresol = CPXopenCPLEX(&error);
	CPXLPptr lpPresol = CPXcreateprob(envPresol, &error, "presol");
	CPXLPptr lpReduced = CPXcreateprob(envPresol, &error, "reduced");
	CPXreadcopyprob(envPresol, lpPresol, &inst.inputFilePath[0], NULL);

	CPXsetintparam(envPresol, CPX_PARAM_THREADS, inst.numThreads);

	// CPLEX+DOMINANCE PRESOLVER: launch cplex and stop the computation at the root node
	// collect the reduced model and launch the dominance presolver on it
	// REMOVE (?): the same result can be obtained by using the cplex presolver and providing
	// the reduced model as the input for the dominance routine
	CPXgettime(env, &timeStart);
	/*if (strcmp(inst.presolver, "cpxdom") == 0) {
		CPXsetintparam(envPresol, CPX_PARAM_NODELIM, 1);

		CPXmipopt(envPresol, lpPresol);
		CPXgetredlp(envPresol, lpPresol, &lpReduced);

		printf("cols = %d\n", CPXgetnumcols(envPresol, lpReduced));

		error = CPXwriteprob(envPresol, lpReduced, "reduced_model.lp", NULL);
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
	else // No presolver
	{
		comm_log(inst, 0, "No presolver selected.");
		CPXreadcopyprob(env, lp, &inst.inputFilePath[0], NULL);
	}*/
	CPXgettime(env, &timeEnd);

	CPXfreeprob(envPresol, &lpReduced);
	CPXfreeprob(envPresol, &lpPresol);
	CPXcloseCPLEX(&envPresol);

	inst.timePresolver = timeEnd - timeStart;


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
		CPXgettime(env, &timeCurr);
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit - (timeCurr - inst.timeStart));
		CPXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);

		CPXmipopt(env, lp);
	}
	else if (inst.solver.compare("balasbcrule1") == 0)
	{
		cpxsol_balas_rule1(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-test") == 0)
	{
		cpxsol_balas_rule1_test(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule1-sparse") == 0)
	{
		cpxsol_balas_rule1_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("balasbcrule2") == 0)
	{
		cpxsol_balas_rule2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol") == 0)
	{
		cpxsol_maxcol(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol2") == 0)
	{
		cpxsol_maxcol2(inst, env, lp);
	}
	else if (inst.solver.compare("maxcol-sparse") == 0)
	{
		cpxsol_maxcol_sparse(inst, env, lp);
	}
	else if (inst.solver.compare("maxcoldom") == 0)
	{
		cpxsol_maxcol_dom(inst, env, lp);
	}
	else
	{
		std::cerr << thisFuncName + " - solver unknown: " + inst.solver << std::endl;
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);

		return SC_GENERIC_ERROR;
	}
	CPXgettime(env, &timeEnd);
	inst.timeSolver = timeEnd - timeStart;
	inst.timeTotal = timeEnd - inst.timeStart;


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
	std::printf("Solver time      = %10.4lf\n", inst.timeSolver);
	std::printf("Total time       = %10.4lf\n", inst.timeTotal);
	std::printf("Best obj val     = %10.4lf\n", inst.bestObjVal);
	std::printf("Obj val          = %10.4lf\n", inst.objVal);
	std::printf("Node total       = %10d\n", CPXgetnodecnt(env, lp));
	std::printf("Node left        = %10d\n\n", CPXgetnodeleftcnt(env, lp));

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return SC_SUCCESFULL;
}

/** Dominance presolver routine: see report.
 *	The function receives a SCP instance, scanning it to find the last column
 *	with objective equals to 1.0.
 *	The upper bounds of the variables with cost greater than 1.0 are set to 0. The
 *	procedure scans these columns to find the ones that must be included in the
 *	reduced model (resetting the upper bound to 1).
 */
STATUS cpxsol_preproc_dominance(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	char b = 'B', u = 'U';
	int i, j, ncols, nrows, stat, lastcol, removed_cols = 0;
	double coeff, obj_val, obj, one = 1.0, rhs, zero = 0.0;

	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);

	// skip all columns with cost == 1
	for (j = 0; j < ncols; ++j)
	{
		CPXgetobj(env, lp, &coeff, j, j + 1);
		if (coeff > 1.0)
		{
			lastcol = j;
			break;
		}
	}

	// fix all columns with cost > 1 to 0
	zero = 0.0;
	for (j = lastcol; j < ncols; ++j)
	{
		CPXchgbds(env, lp, 1, &j, &b, &zero);
	}

	//printf("Dominance pre-solver receives %d rows and %d columns (selected: %d).\n", nrows, ncols, inst->num_selected_cols);

	// set cplex params
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	//CPXsetdblparam(env, CPX_PARAM_TILIM, 5.0);

	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 2);
	CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);

	// pass every column with cost > 1
	for (j = lastcol; j < ncols; ++j)
	{

		// set rhs of each row to next column values
		for (i = 0; i < nrows; ++i)
		{
			CPXgetcoef(env, lp, i, j, &rhs);
			if (CPXchgrhs(env, lp, 1, &i, &rhs))
				perror("ERROR: error at CPXchgrhs.");
		}

		// optimize subproblem
		CPXmipopt(env, lp);

		CPXgetobjval(env, lp, &obj_val);
		CPXgetobj(env, lp, &obj, j, j);
		stat = CPXgetstat(env, lp);

		// infeasible solution OR
		// optimal solution found but obj val > column cost
		// THEN add this column to the model
		if (stat == CPXMIP_INFEASIBLE || (stat == CPXMIP_OPTIMAL && (obj_val > obj)))
		{

			// reset upper bound of the column to 1
			one = 1.0;
			CPXchgbds(env, lp, 1, &j, &u, &one);
		}
		// column not selected
		else
		{
			removed_cols++;
		}

		//if (j % 100 == 0) printf("Dominance pre-solver analyses %d of %d columns (removed: %d).\n", j, ncols, removed_cols);
	}

	// reset rhs of all rows to 1 (standard set cover)
	one = 1.0;
	for (i = 0; i < nrows; ++i)
	{
		if (CPXchgrhs(env, lp, 1, &i, &one))
			perror("ERROR: error at CPXchgrhs.");
	}

	// count selected columns
	//printf("Dominance pre-solver removed %d columns.\n\n", removed_cols);

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_cuts_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp)
{
    /*inst.balasnodecnt = 0;
    inst.cplexnodecnt = 0;

    inst.SC_BALAS_NODE_CUTS_NUM = 10;

	status = CPXsetusercutcallbackfunc(env, cpxcb_bala_susercuts_sparse, inst);

    CPXmipopt(env, lp);*/

    return 1;
}

STATUS cpxsol_balas_cuts(SCinstance &inst, IloEnv env, IloCplex cplex)
{
	/*// usato per contare i tagli in questo caso
	inst.balasnodecnt = 0;

	status = CPXsetusercutcallbackfunc(env, cpxcb_balas_usercuts_test, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_rule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, cpxcb_balas_branch_rule1v1, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_rule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, cpxcb_balas_branch_rule1_test, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_rule1_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, cpxcb_balas_branch_rule1_sparse, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_rule1_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*inst.balasCntCode = 0;
	inst.cplexCntNode = 0;

	status = CPXsetbranchcallbackfunc(env, cpxcb_balas_branch_rule1_maxcol_sparse, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_balas_rule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*status = CPXsetbranchcallbackfunc(env, cpxcb_balas_branch_rule2, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_maxcol(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	int status;
	double timeCurr;
	std::string thisFuncName = "cpxsol_maxcol";

	CPXgettime(env, &timeCurr);

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // Disable problem reduction
	CPXsetintparam(env, CPX_PARAM_REDUCE, 1);
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	CPXsetdblparam(env, CPXPARAM_MIP_Limits_CutsFactor, 0.0); // Turn off cuts generation
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst.mipTimeLimit - (timeCurr - inst.timeStart));
	CPXsetintparam(env, CPX_PARAM_THREADS, inst.numThreads);

	status = CPXsetbranchcallbackfunc(env, cpxcb_branch_maxcol, &inst);
	if (status)
	{
		std::cerr << thisFuncName + " - CPXsetbranchcallbackfunc return exit status " << status << std::endl;
		return SC_ERR_CALLBACK;
	}

	status = CPXmipopt(env, lp);
	if (status)
	{
		std::cerr << thisFuncName + " - CPXmipopt return exit status " << status << std::endl;
		return SC_ERR_CALLBACK;
	}

	return SC_SUCCESFULL;
}

STATUS cpxsol_maxcol2(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*int status;

	status = CPXsetbranchcallbackfunc(env, cpxcb_branchmaxcol, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*status = CPXsetbranchcallbackfunc(env, cpxcb_branchmaxcol_sparse, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

STATUS cpxsol_maxcol_dom(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	/*status = CPXsetbranchcallbackfunc(env, cpxcb_branchmaxcoldom, inst);

	CPXmipopt(env, lp);*/

	return SC_SUCCESFULL;
}

/**
 * Take the cplex lp model and fill the values on
 * the instance objective and the instance matrix.
 * 
 * @param inst - SCinstance current instance
 * @param env - CPXENVptr cplex environment
 * @param lp - CPXLPptr cplex linear programming model
 * @return a status code
 */
STATUS cpxsol_build_lp2raw(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	size_t i;
	size_t j;
	double val;

	// Read objective function from lp model
	inst.dnsobj = arma::vec(CPXgetnumcols(env, lp));
	for (j = 0; j < CPXgetnumcols(env, lp); ++j)
	{
		CPXgetobj(env, lp, &val, j, j);
		inst.dnsobj(j) = val;
	}

	// Read matrix from lp model
	inst.dnsmat = arma::mat(CPXgetnumrows(env, lp), CPXgetnumcols(env, lp));
	for (i = 0; i < (size_t)CPXgetnumrows(env, lp); ++i)
	{
		for (j = 0; j < (size_t)CPXgetnumcols(env, lp); ++j)
		{
			CPXgetcoef(env, lp, i, j, &val);
			inst.dnsmat(i, j) = val;
		}
	}

	return SC_SUCCESFULL;
}

/**
 * Take the raw model in instance (inst.dnsobj and inst.dnsmat)
 * and fill the cplex lp model.
 * 
 * @param inst - SCinstance current instance
 * @param env - CPXENVptr cplex environment
 * @param lp - CPXLPptr cplex empty linear programming model to build
 * @return a status code
 */
STATUS cpxsol_build_raw2lp(SCinstance &inst, CPXENVptr env, CPXLPptr lp)
{
	char binCol = 'B';
	char geSense = 'G';
	int status;
	size_t i;
	size_t j;
	double rhs = 1.0;
	double val;

	// Set objective function on lp model
	for (j = 0; j < inst.dnsobj.n_elem; ++j)
	{
		val = inst.dnsobj(j);
		status = CPXnewcols(env, lp, 1, &val, NULL, NULL, &binCol, NULL);
	}

	// Set matrix on lp model
	for (i = 0; i < inst.dnsmat.n_rows; ++i)
	{
		CPXnewrows(env, lp, 1, &rhs, &geSense, NULL, NULL);
		for (j = 0; j < inst.dnsmat.n_cols; ++j)
		{
			val = inst.dnsmat(i, j);
			status = CPXchgcoef(env, lp, i, j, val);
		}
	}

	return SC_SUCCESFULL;
}
