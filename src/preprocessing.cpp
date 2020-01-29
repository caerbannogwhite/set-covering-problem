
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

/*
int util_read_input(SCinstance *inst) {
    int line_counter = 0, curr_col_num, col_ind = 0, row_ind = 0, num, i;
    char new_row = 0;

    FILE *file_handler = fopen(inst->input_file_name, "r");

    while (!feof(file_handler)) {

        // first row: number of rows and cols
        if (line_counter == 0) {
            if (!fscanf(file_handler, "%d %d", &inst->num_rows, &inst->num_cols))
                perror("ERROR: fscan error at main_read_input-first row.");

            inst->costs = (double *) calloc(inst->num_cols, sizeof(double));
            inst->table = (double *) calloc(inst->num_rows * inst->num_cols, sizeof(double));

        }

            // costs
        else if (line_counter > 0 && line_counter <= inst->num_cols) {
            if (!fscanf(file_handler, "%d", &num)) perror("ERROR: fscan error at main_read_input-costs.");
            inst->costs[col_ind++] = num;

            if (line_counter == inst->num_cols) {
                col_ind = 0;
                row_ind = 0;
                new_row = 1;
            }
        }

            // new row
        else if (new_row) {
            if (!fscanf(file_handler, "%d", &curr_col_num)) perror("ERROR: fscan error at main_read_input-new row.");
            new_row = 0;
        }

            // entries
        else {

            // num: index of the column whose coefficient must be 1, 1 <= index <= num_cols
            if (!fscanf(file_handler, "%d", &num)) perror("ERROR: fscan error at main_read_input-entries.");
            inst->table[row_ind * inst->num_cols + num - 1] = 1;
            ++col_ind;

            if (curr_col_num == col_ind) {
                col_ind = 0;
                ++row_ind;
                new_row = 1;
            }
        }

        ++line_counter;
    }

    fclose(file_handler);
}

int util_print_solution(SCinstance *inst) {
    int i, cnt = 0;

    FILE *file_handler = fopen("solution", "w");

    printf("\nSOLUTION (non-zero columns):");
    for (i = 0; i < inst->num_cols; ++i) {

        // selected_columns[i] contains the index of the i-th columns in the reduced model
        if (inst->selected_columns[i] > -1 && inst->x_star[inst->selected_columns[i]] > (0 + EPSILON_MED)) {
            fprintf(file_handler, "%d\n", i);
            printf("%d,", i);
            cnt++;
        }
    }
    printf("\n");
    fclose(file_handler);

    if (inst->debug) {
        printf("\nDEBUG:\n");
        for (i = 0; i < inst->num_cols; ++i) {
            printf("i=%d, sel_col[i]=%d, x_star[sel_col[i]]=%.2lf\n", i, inst->selected_columns[i],
                   inst->x_star[inst->selected_columns[i]]);
        }
        printf("\n");
    }

    inst->solution_len = cnt;
    return 0;
}

int util_print_table(SCinstance *inst) {
    int i, j;

    for (i = 0; i < inst->num_rows; ++i) {
        for (j = 0; j < inst->num_cols; ++j) {
            printf("%d ", inst->table[inst->num_cols * i + j]);
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}

int sc_build_model(SCinstance *inst, CPXENVptr env, CPXLPptr lp) {

    char binary = 'B';
    char greater_equal = 'G';
    int i, j, lastrow;
    double obj, one = 1.0, zero = 0.0;

    // COLUMNS
    for (i = 0; i < inst->num_cols; ++i)
    {
        // ccnt: number of new variables being added (always 1 here)
        // obj: objective function coefficient
        // lb: lower bound
        // ub: upper bound, NULL means CPX_INFBOUND
        // xctype: type of the variable (INTEGER, BINARY, CONTINUOUS)
        // colname: name of the variable

        if (inst->use_unit_costs) obj = (double) 1.0;
        else obj = (double) inst->costs[i];

        if (CPXnewcols(env, lp, 1, &obj, &zero, &one, &binary, NULL)) perror("ERROR: error at CPXnewcols.");
    }

    // ROWS
    for (i = 0; i < inst->num_rows; ++i)
    {
        // rcnt: number of new rows being added (always 1 here)
        // rhs: right hand side
        // sense: 'E' for equations
        // rngval: range values for the new constraints
        // rowname: name of the constraint

        lastrow = CPXgetnumrows(env, lp);
        if (CPXnewrows(env, lp, 1, &one, &greater_equal, NULL, NULL)) perror("ERROR: error at CPXnewrows.");
        for (j = 0; j < inst->num_cols; ++j) {

            if (CPXchgcoef(env, lp, lastrow, j, inst->table[i * inst->num_cols + j])) perror("ERROR: error at CPXchgcoef.");
        }
    }

    return 0;
}

int sc_build_dual_model(SCinstance *inst, CPXENVptr env, CPXLPptr lp, CPXLPptr dual_lp) {

    char lower_equal = 'L';
    int i, j, lastrow, nrows, ncols;
    double coeff, rhs, one = 1.0, zero = 0.0;

    nrows = CPXgetnumrows(env, lp);
    ncols = CPXgetnumcols(env, lp);

    // COLUMNS # cols_dual == # rows prim
    for (i = 0; i < nrows; ++i)
    {
        // ccnt: number of new variables being added (always 1 here)
        // obj: objective function coefficient
        // lb: lower bound
        // ub: upper bound, NULL means CPX_INFBOUND
        // xctype: type of the variable (INTEGER, BINARY, CONTINUOUS)
        // colname: name of the variable

        if (CPXnewcols(env, dual_lp, 1, &one, &zero, NULL, &integer, NULL)) perror("ERROR: error at CPXnewcols.");
    }

    // ROWS
    for (j = 0; j < ncols; ++j)
    {
        // rcnt: number of new rows being added (always 1 here)
        // rhs: right hand side
        // sense: 'E' for equations
        // rngval: range values for the new constraints
        // rowname: name of the constraint

        lastrow = CPXgetnumrows(env, dual_lp);
        CPXgetobj(env, lp, &rhs, j, j+1);       // get the cost of col j and put in rhs
        if (CPXnewrows(env, dual_lp, 1, &rhs, &lower_equal, NULL, NULL)) perror("ERROR: error at CPXnewrows.");
        for (i = 0; i < nrows; ++i) {
            CPXgetcoef(env, lp, i, j, &coeff);
            if (CPXchgcoef(env, dual_lp, lastrow, i, coeff)) perror("ERROR: error at CPXchgcoef.");
        }
    }

    // maximize
    CPXchgobjsen(env, dual_lp, CPX_MAX);

    return 0;
}
*/
