
#include "preprocessing.hpp"

/** Dominance presolver routine: see report.
 *	The function receives a set-covering instance, scanning it to find the last column
 *	with objective equals to 1.0.
 *	The upper bounds of the variables with cost greater than 1.0 are set to 0. The
 *	procedure scans these columns to find the ones that must be included in the
 *	reduced model (resetting the upper bound to 1).
 */
int SCdominancepresolver(SCinstance *inst, CPXENVptr env, CPXLPptr lp)
{

    char b = 'B', u = 'U';
    int i, j, ncols, nrows, stat, lastcol, removed_cols = 0;
    double coeff, obj_val, obj, one = 1.0, rhs, zero = 0.0;

    ncols = CPXgetnumcols(env, lp);
    nrows = CPXgetnumrows(env, lp);

    // skip all columns with cost == 1
    for (j = 0; j < ncols; ++j) {
        CPXgetobj(env, lp, &coeff, j, j+1);
        if (coeff > 1.0) {
            lastcol = j;
            break;
        }
    }

    // fix all columns with cost > 1 to 0
    zero = 0.0;
    for (j = lastcol; j < ncols; ++j) {
        CPXchgbds(env, lp, 1, &j, &b, &zero);
    }

    //printf("Dominance pre-solver receives %d rows and %d columns (selected: %d).\n", nrows, ncols, inst->num_selected_cols);

    // set cplex params
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
    //CPXsetdblparam(env, CPX_PARAM_TILIM, 5.0);

    CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 2);
    CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 0);

    // pass every column with cost > 1
    for (j = lastcol; j < ncols; ++j) {

        // set rhs of each row to next column values
        for (i = 0; i < nrows; ++i) {
            CPXgetcoef(env, lp, i, j, &rhs);
            if (CPXchgrhs(env, lp, 1, &i, &rhs)) perror("ERROR: error at CPXchgrhs.");
        }

        // optimize subproblem
        CPXmipopt(env, lp);

        CPXgetobjval(env, lp, &obj_val);
        CPXgetobj(env, lp, &obj, j, j);
        stat = CPXgetstat(env, lp);

        // infeasible solution OR
        // optimal solution found but obj val > column cost
        // THEN add this column to the model
        if (stat == CPXMIP_INFEASIBLE || (stat == CPXMIP_OPTIMAL && (obj_val > obj))) {

            // reset upper bound of the column to 1
            one = 1.0;
            CPXchgbds(env, lp, 1, &j, &u, &one);
        }
            // column not selected
        else {
            removed_cols++;
        }

        //if (j % 100 == 0) printf("Dominance pre-solver analyses %d of %d columns (removed: %d).\n", j, ncols, removed_cols);
    }

    // reset rhs of all rows to 1 (standard set cover)
    one = 1.0;
    for (i = 0; i < nrows; ++i) {
        if (CPXchgrhs(env, lp, 1, &i, &one)) perror("ERROR: error at CPXchgrhs.");
    }

    // count selected columns
    //printf("Dominance pre-solver removed %d columns.\n\n", removed_cols);

    return 0;
}
