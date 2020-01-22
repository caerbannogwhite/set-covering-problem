//
// Created by macs on 19/06/19.
//

//#include <cblas.h>
#include <string.h>
#include <cplex.h>
#include <math.h>
#include "callbacks.h"
#include "aux.h"
#include "balas_dense.h"

int CPXPUBLIC SCcallbackbalasusercuts_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
        int *useraction_p) {

    char sense;
    int nrows, ncols, numnz, nzcnt, status, i, j, cnt, nzcnt_tmp, surplus;
    int nrowsred, ncolsred, numnzred, flag;
    int *rmatbeg, *rmatind, *rmatbegred, *rmatindred, *cutind, *colsred2sc, *colssc2red, *rowsred2sc, *rowssc2red;
    double rhs;
    double *ub, *lb, *rmatval, *x, *cutval, *objred;

    CPXCLPptr lp;

    /*CPXgetcallbacklp(env, cbdata, wherefrom, &lp);

    // Get reduced matrix and objective
    ncols = CPXgetnumcols(env, lp);
    nrows = CPXgetnumrows(env, lp);
    numnz = CPXgetnumnz(env, lp);

    // Get the upper bounds
    ub = (double *) malloc(ncols * sizeof(double));
    status = CPXgetub(env, lp, ub, 0, ncols - 1);
    if (status) { perror("Can't get the ub"); }

    // Get the lower bounds
    lb = (double *) malloc(ncols * sizeof(double));
    status = CPXgetlb(env, lp, lb, 0, ncols - 1);
    if (status) { perror("Can't get the lb"); }

    // Get the matrix of the node lp
    rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
    rmatind = (int *) malloc(numnz * sizeof(int));
    rmatval = (double *) malloc(numnz * sizeof(double));

    rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
    rmatindred = (int *) malloc(numnz * sizeof(int));

    colsred2sc = (int *) malloc(ncols * sizeof(int));
    colssc2red = (int *) malloc(ncols * sizeof(int));
    rowsred2sc = (int *) malloc(nrows * sizeof(int));
    rowssc2red = (int *) malloc(nrows * sizeof(int));

    status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
    if (status) { perror("Can't get rows"); }

    // remove 1 vars rows and cols, remove 0 vars cols
    cnt = 0;
    for (j = 0; j < ncols; ++j) {
        colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
        colsred2sc[cnt] = j;
        cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
    }
    ncolsred = cnt;

    nrowsred = 0;
    numnzred = 0;
    rmatbeg[nrows] = numnz;
    for (i = 0; i < nrows; ++i) {

        flag = 1;
        rmatbegred[nrowsred] = numnzred;
        for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

            // the row is covered: do not include it in reduced matrix
            // clean up the row and set flag to 0
            if (lb[rmatind[j]] > 0.5) {
                flag = 0;
                numnzred = rmatbegred[nrowsred];
                break;
            }

            if (colssc2red[rmatind[j]] > -1) {
                rmatindred[numnzred] = colssc2red[rmatind[j]];
                numnzred++;
            }
        }

        if (flag) {
            rowssc2red[i] = nrowsred;
            rowsred2sc[nrowsred] = i;
            nrowsred++;
        } else {
            rowssc2red[i] = -1;
        }
    }
    rmatbegred[nrowsred] = numnzred;

    free(ub);
    free(lb);

    free(rmatbeg);
    free(rmatind);
    free(rmatval);

    objred = (double *) malloc(ncolsred * sizeof(double));
    for (j = 0; j < ncolsred; ++j) {
        objred[j] = inst->costs[colsred2sc[j]];
    }

    cmatbegred = (int *) malloc((ncolsred + 1) * sizeof(int));
    cmatindred = (int *) malloc(numnzred * sizeof(int));

    cnt = 0;
    for (j = 0; j < ncolsred; ++j) {
        cmatbegred[j] = cnt;
        for (i = 0; i < nrowsred; ++i) {
            item = (int *) bsearch(&j, &rmatindred[rmatbegred[i]], rmatbegred[i+1] - rmatbegred[i], sizeof(int), SCint_cmp);
            if (item != NULL) {
                cmatindred[cnt] = i;
                cnt++;
            }
        }
    }
    cmatbegred[ncolsred] = cnt;


    // Primal
    xindlen = 0;
    xind = (int *) malloc(ncolsred * sizeof(int));
    incumbent = (double *) malloc(ncols * sizeof(double));
    CPXgetcallbackincumbent(env, cbdata, wherefrom, incumbent, 0, ncols - 1);

    for (j = 0; j < ncols; ++j) {
        if (incumbent[j] > 0.5 && colssc2red[j] > -1) {
            xind[xindlen] = colssc2red[j];
            xindlen++;
        }
    }
    free(incumbent);

    zu = SCbalasheurprimal0_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, objred, nrowsred, ncolsred, xind, &xindlen, 3);

    // Dual and reduced costs
    dj = (double *) malloc(ncolsred * sizeof(double));
    pi = (double *) malloc(nrowsred * sizeof(double));

    for (i = 0; i < nrowsred; ++i) {
        CPXgetpi(env, lp, &pi[i], rowsred2sc[i], rowsred2sc[i]);
    }

    // Check if is dual solution
    double acc;
    for (j = 0; j < ncolsred; ++j) {
        acc = 0.0;
        for (i = cmatbegred[j]; i < cmatbegred[j+1]; ++i) {
            acc += pi[cmatindred[i]];
        }

        if ((acc - SC_EPSILON_SMALL) > objred[j]) {
            printf("NOT DUAL SOLUTION. acc=%lf, obj=%lf\n", acc, objred[j]);
            break;
        }
    }

    for (j = 0; j < ncolsred; ++j) {
        dj[j] = objred[j];
    }

    for (i = 0; i < nrowsred; ++i) {
        for (j = rmatbegred[i]; j < rmatbegred[i+1]; ++j) {
            dj[rmatindred[j]] -= pi[i];
        }
    }

    for (j = 0; j < ncolsred; ++j) {
        if (dj[j] < -SC_EPSILON_SMALL) {
            printf("NOT VALID REDUCED COST. dj=%lf, j=%d\n", dj[j], j);
            break;
        }
    }

    //SCbalasheurdual1_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj);
    //SCdual0_sparse(rmatbegred, rmatindred, objred, nrowsred, ncolsred, pi, dj);

    y = 0;
    for (i = 0; i < nrowsred; ++i) {
        y += pi[i];
    }

    if ((y < 0) || (abs(y) == NAN) || (zu <= y)) {
        free(xind);
        free(dj);
        free(pi);

        free(colsred2sc);
        free(colssc2red);
        free(rowsred2sc);
        free(rowssc2red);

        free(rmatbegred);
        free(rmatindred);

        free(cmatbegred);
        free(cmatindred);

        free(objred);

        SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
        *useraction_p = CPX_CALLBACK_SET;
        status = 0;

        goto TERMINATE;
    }

    SCbalasheurdual3_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj, zu);

    y = 0;
    for (i = 0; i < nrowsred; ++i) {
        y += pi[i];
    }

    cutind = (int *) malloc(ncolsred * sizeof(int));
    nzcnt = SCbalascut_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, nrowsred, ncolsred, xind, xindlen, dj,
            zu, y, cutind);

    cutval = (double *) malloc(nzcnt * sizeof(double));
    for (i = 0; i < nzcnt; ++i) {
        cutval[i] = 1.0;
    }

    rhs = 1.0;
    sense = 'G';
    CPXcutcallbackaddlocal(env, cbdata, wherefrom, nzcnt, &rhs, &sense, cutind, cutval);

    free(cutind);
    free(cutval);*/

    TERMINATE:
    return status;
}


int CPXPUBLIC SCcallbackbalasusercuts(CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *useraction_p) {

	unsigned char *mat, *xchr;
	int ncolsred, nrowsred, ncols, nrows, nodedepth, status = 0;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *xsupp;
	unsigned int i, j, xsupplen;
	double objval, zu, y;
	double *dj, *pi, *obj;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// ROOT NODE: do nothing
	if (nodedepth == 0) {
		goto TERMINATE;
	}

	// Get reduced matrix and objective
	nrows = CPXgetnumrows(env, lp);
	ncols = CPXgetnumcols(env, lp);
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	mat = SCgetmatrix(env, lp, rowsred2sc, rowssc2red, colsred2sc, colssc2red, &nrowsred, &ncolsred);
	obj = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		obj[j] = inst->costs[colsred2sc[j]];
	}

	// Primal
	xsupplen = 0;
	xsupp = (int *) malloc(ncolsred * sizeof(int));
	xchr = (unsigned char *) calloc(ncolsred, sizeof(unsigned char));

	zu = SCbalasheurprimal0(mat, inst->costs, nrowsred, ncolsred, xchr, xsupp, &xsupplen, 3);

	if (zu <= 0) {
		free(colsred2sc);
		free(colssc2red);
		free(mat);
		free(obj);
		free(xsupp);
		free(xchr);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;

		goto TERMINATE;
	}

	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (j = 0; j < ncolsred; ++j) {
		dj[j] = obj[j];
	}

	SCbalasheurdual1(mat, nrowsred, ncolsred, xchr, pi, dj);
	//SCdual0(mat, obj, nrowsred, ncolsred, pi, dj);

	status = SCbalasheurdual3(mat, nrowsred, ncolsred, xchr, pi, dj, zu);
	if (status < 0) {
		free(colsred2sc);
		free(colssc2red);
		free(mat);
		free(obj);
		free(xsupp);
		free(xchr);
		free(dj);
		free(pi);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;

		goto TERMINATE;
	}

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	int *cutind = (int *) calloc(ncolsred, sizeof(int));
	int cutnzcnt = SCbalascutprocedure(mat, nrowsred, ncolsred, xchr, dj, zu, y, cutind);

	double *cutval = (double *) malloc(cutnzcnt * sizeof(double));
	for (j = 0; j < cutnzcnt; ++j) {
		cutval[j] = 1.0;
	}

	CPXcutcallbackaddlocal(env, cbdata, wherefrom, cutnzcnt, 1.0, 'G', cutind, cutval);

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);
	free(mat);
	free(obj);
	free(xsupp);
	free(xchr);
	free(dj);
	free(pi);
	free(cutind);
	free(cutval);

	TERMINATE:
	return status;
} /* END SCcallbackbalasusercuts */


int CPXPUBLIC SCcallbackbalasusercuts_test(CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *useraction_p) {

	unsigned char *mat, *xchr;
	int off, ncolsred, nrowsred, ncols, nrows, nodedepth, status = 0;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *xsupp;
	unsigned int i, j, xsupplen;
	double objval, zu, y;
	double *dj, *pi, *obj;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// ROOT NODE: do nothing
	if (nodedepth == 0) {
		goto TERMINATE;
	}

	// Get reduced matrix and objective
	nrows = CPXgetnumrows(env, lp);
	ncols = CPXgetnumcols(env, lp);
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	mat = SCgetmatrix(env, lp, rowsred2sc, rowssc2red, colsred2sc, colssc2red, &nrowsred, &ncolsred);
	obj = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		obj[j] = inst->costs[colsred2sc[j]];
	}

	// Primal
	xchr = (unsigned char *) calloc(ncolsred, sizeof(unsigned char));

	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (i = 0; i < nrowsred; ++i) {
		CPXgetpi(env, lp, &pi[i], rowsred2sc[i], rowsred2sc[i]);
	}

	// Copy obj in dj
	memcpy(dj, obj, ncolsred * sizeof(double));
	//cblas_ccopy(ncolsred, obj, 1, dj, 1);

	for (i = 0; i < nrowsred; ++i) {
		off = i * ncolsred;
		for (j = 0; j < ncolsred; ++j) {
			dj[j] -= mat[off + j] * pi[i];
		}
	}

	//cblas_dgemv(CblasRowMajor, CblasTrans, nrowsred, ncolsred, -1.0, mat, ncolsred, pi, 1, 1.0, dj, 1);

	// Make dual feasible
	int l;
	for (j = 0; j < ncolsred; ++j) {
		for (i = 0; i < nrowsred; ++i) {
			if (dj[j] >= SC_EPSILON_SMALL) {
				break;
			}

			if (mat[i * ncolsred + j]) {
				if (dj[j] < -pi[i]) {
					for (l = 0; l < ncolsred; ++l) {
						dj[l] += pi[i] * mat[i * ncolsred + l];
					}
				} else {
					pi[i] += dj[j];
					for (l = 0; l < ncolsred; ++l) {
						dj[l] -= dj[j] * mat[i * ncolsred + l];
					}
				}
			}
		}
	}

	//zu = SCbalasheurprimal5b(mat, obj, nrowsred, ncolsred, xchr, pi, dj);
	xsupp = (int *) malloc(ncolsred * sizeof(int));
	xsupplen = 0;
	zu = SCbalasheurprimal0(mat, obj, nrowsred, ncolsred, xchr, xsupp, &xsupplen, 3);
	free(xsupp);
	SCbalasheurdual3(mat, nrowsred, ncolsred, xchr, pi, dj, zu);
	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	//printf("zl=%lf, zu=%lf\n", y, zu);

	//if (!SCiscover(mat, nrowsred, ncolsred, xchr)) printf("NOT COVER.\n");
	//if (!SCisdualsolution(mat, obj, nrowsred, ncolsred, pi)) printf("NOT DUAL SOLUTION.\n");

	if ((zu <= y) || (!SCisdualsolution(mat, obj, nrowsred, ncolsred, pi))) {
		free(colsred2sc);
		free(colssc2red);
		free(rowsred2sc);
		free(rowssc2red);
		free(mat);
		free(obj);
		free(xchr);
		free(dj);
		free(pi);

		*useraction_p = CPX_CALLBACK_DEFAULT;
		status = 0;

		goto TERMINATE;
	}

	int *cutind = (int *) calloc(ncolsred, sizeof(int));
	int cutnzcnt = SCbalascutprocedure(mat, nrowsred, ncolsred, xchr, dj, zu, y, cutind);

	double *cutval = (double *) malloc(cutnzcnt * sizeof(double));
	for (j = 0; j < cutnzcnt; ++j) {
		cutval[j] = 1.0;
	}

	CPXcutcallbackaddlocal(env, cbdata, wherefrom, cutnzcnt, 1.0, 'G', cutind, cutval);

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);
	free(mat);
	free(obj);
	free(xchr);
	free(dj);
	free(pi);
	free(cutind);
	free(cutval);

	TERMINATE:
	return status;
} /* END SCcallbackbalasusercuts_test */


int CPXPUBLIC SCcallbackbranchmaxcol(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) cbhandle;
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char flag;
	unsigned char *mat, *matr;
	int ncolsred, nrowsred, ncols, nrows, numnz, nzcnt_tmp, surplus, jt, status = 0, seqnum1, seqnum2;
	int *colsred2sc, *colssc2red, *rmatbeg, *rmatind;
	unsigned int i, j, k, l, cnt, max;
	double objval;
	double *rmatval, *ub, *lb;

	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc(nrows * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }


	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;

	mat = (unsigned char *) calloc(sizeof(unsigned char), (nrows) * (ncolsred));
	nrowsred = 0;
	for (i = 0; i < nrows; ++i) {
		matr = &mat[(nrowsred) * (ncolsred)];
		k = (i == nrows - 1 ? numnz : rmatbeg[i + 1]);
		flag = 1;
		for (j = rmatbeg[i]; j < k; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				for (l = 0; l < ncolsred; ++l) {
					matr[l] = 0;
				}
				flag = 0;
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				matr[colssc2red[rmatind[j]]] = 1;
			}
		}

		if (flag) {
			nrowsred = nrowsred + 1;
		}
	}

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	max = 0;
	jt = -1;
	for (j = ncolsred - 1; j--; ) {
		cnt = 0;
		for (i = nrowsred - 1; i--;  ) {
			cnt += mat[i * ncolsred + j];
		}

		if (cnt >= max) {
			max = cnt;
			jt = j;
		}
	}

	jt = colsred2sc[jt];
	int varind[1] = {jt};
	char varlu[1] = {'U'};
	double varbd[1] = {0.0};
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, varind, varlu, varbd, objval, NULL, &seqnum1);
	if (status) { perror("Can't use CPXbranchcallbackbranchbds (1)"); }

	varlu[0] = 'L';
	varbd[0] = 1.0;
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, varind, varlu, varbd, objval, NULL, &seqnum2);
	if (status) { perror("Can't use CPXbranchcallbackbranchbds (2)"); }

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);
	free(colssc2red);
	free(mat);

	return status;
} /* END SCcallbackbranchmaxcol */


int CPXPUBLIC SCcallbackbranchmaxcol2(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) cbhandle;
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char flag;
	unsigned char *mat, *matr, *acc;
	int ncolsred, nrowsred, ncols, nrows, numnz, nzcnt_tmp, surplus, jt, status = 0, seqnum1, seqnum2;
	int *colsred2sc, *colssc2red, *rmatbeg, *rmatind;
	unsigned int i, j, k, l, cnt, max;
	double objval;
	double *rmatval, *ub, *lb;

	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc(nrows * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }


	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;

	mat = (unsigned char *) calloc(sizeof(unsigned char), (nrows) * (ncolsred));
	nrowsred = 0;
	for (i = 0; i < nrows; ++i) {
		matr = &mat[(nrowsred) * (ncolsred)];
		k = (i == nrows - 1 ? numnz : rmatbeg[i + 1]);
		flag = 1;
		for (j = rmatbeg[i]; j < k; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				for (l = 0; l < ncolsred; ++l) {
					matr[l] = 0;
				}
				flag = 0;
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				matr[colssc2red[rmatind[j]]] = 1;
			}
		}

		if (flag) {
			nrowsred = nrowsred + 1;
		}
	}

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	acc = (unsigned char *) malloc(ncolsred * sizeof(unsigned char));
	max = 0;
	jt = -1;
	for (j = 0; j < ncolsred; ++j) {

		for (k = 0; k < ncolsred; ++k) {
			acc[k] = 0;
		}

		for (i = 0; i < nrowsred; ++i) {
			if (mat[i*ncolsred + j]) {
				matr = &mat[i* ncols];
				for (k = 0; k < ncolsred; ++k) {
					acc[k] |= matr[k];
				}
			}
		}

		cnt = 0;
		for (k = 0; k < ncolsred; ++k) {
			cnt += acc[k];
		}

		if (cnt > max) {
			max = cnt;
			jt = j;
		}
	}

	free(acc);

	jt = colsred2sc[jt];
	int varind[1] = {jt};
	char varlu[1] = {'U'};
	double varbd[1] = {0.0};
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, varind, varlu, varbd, objval, NULL, &seqnum1);
	if (status) { perror("Can't use CPXbranchcallbackbranchbds (1)"); }

	varlu[0] = 'L';
	varbd[0] = 1.0;
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, varind, varlu, varbd, objval, NULL, &seqnum2);
	if (status) { perror("Can't use CPXbranchcallbackbranchbds (2)"); }

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);
	free(colssc2red);
	free(mat);

	return status;
} /* END SCcallbackbranchmaxcol */


int CPXPUBLIC SCcallbackbranchmaxcol_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

// Silence compiler warning(s) about unused parameter(s).
	(void) cbhandle;
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	unsigned char flag;
	int cnt, ncolsred, nrowsred, numnzred, ncols, nrows, numnz, status = 0, seqnum1, seqnum2, nzcnt_tmp, surplus;
	int brcol;
	int *colscnt, *rmatind, *rmatbeg, *rmatindred, *rmatbegred, *colssc2red, *colsred2sc;
	unsigned int i, j, max;
	double objval;
	double *ub, *lb, *rmatval;

	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
	rmatindred = (int *) malloc(numnz * sizeof(int));

	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }

	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;
	colscnt = (int *) calloc(ncolsred, sizeof(int));

	nrowsred = 0;
	numnzred = 0;
	rmatbeg[nrows] = numnz;
	for (i = 0; i < nrows; ++i) {

		flag = 1;
		rmatbegred[nrowsred] = numnzred;
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				flag = 0;
				numnzred = rmatbegred[nrowsred];
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				rmatindred[numnzred] = colssc2red[rmatind[j]];
				numnzred++;
			}
		}

		if (flag) {

			for (j = rmatbeg[i]; j < rmatbeg[i + 1]; ++j) {
				colscnt[colssc2red[rmatind[j]]] += (colssc2red[rmatind[j]] > -1);
			}

			nrowsred++;
		}
	}

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	brcol = -1;
	max = 0;
	for (j = 0; j < ncolsred; ++j) {
		if (max < colscnt[j]) {
			max = colscnt[j];
			brcol = j;
		}
	}

	char nulllu[1] = {'U'};
	int nullind[1] = { colsred2sc[brcol] };
	double nullval[1] = {0.0};

	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, NULL, &seqnum1);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds (1)"); }

	nulllu[0] = 'L';
	nullval[0] = 1.0;

	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, NULL, &seqnum2);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds (2)"); }

	*useraction_p = CPX_CALLBACK_SET;

	free(colscnt);

	free(colsred2sc);
	free(colssc2red);

	free(rmatbegred);
	free(rmatindred);

	return status;
} /* END SCcallbackbranchmaxcol_sparse */


int CPXPUBLIC SCcallbackbranchmaxcoldom(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char integer = 'I', up = 'U', cond;
	char *senses;
	unsigned char flag;
	unsigned char *removedcols;
	int cnt, ncolsred, nrowsred, numnzred, ncols, nrows, numnz, k, status = 0, seqnum1, seqnum2, nzcnt_tmp, surplus, nodedepth;
	int removedcolscnt, ind, onecolscnt, brcol;
	int *colscnt, *rmatind, *rmatbeg, *rmatindred, *rmatbegred, *colssc2red, *colsred2sc, *rowssc2red, *rowsred2sc;
	unsigned int i, j, max;
	double objval, zero = 0.0, val, objvalred;
	double *ub, *lb, *rmatval;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
	rmatindred = (int *) malloc(numnz * sizeof(int));

	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }

	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;
	colscnt = (int *) calloc(ncolsred, sizeof(int));

	nrowsred = 0;
	numnzred = 0;
	rmatbeg[nrows] = numnz;
	for (i = 0; i < nrows; ++i) {

		flag = 1;
		rmatbegred[nrowsred] = numnzred;
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				flag = 0;
				numnzred = rmatbegred[nrowsred];
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				rmatindred[numnzred] = colssc2red[rmatind[j]];
				numnzred++;
			}
		}

		if (flag) {

			for (j = rmatbeg[i]; j < rmatbeg[i + 1]; ++j) {
				colscnt[colssc2red[rmatind[j]]] += (colssc2red[rmatind[j]] > -1);
			}

			rowssc2red[i] = nrowsred;
			rowsred2sc[nrowsred] = i;
			nrowsred++;
		} else {
			rowssc2red[i] = -1;
		}
	}

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);

	removedcolscnt = 0;
	removedcols = (unsigned char *) calloc(ncolsred, sizeof(unsigned char));

	//if (((nrows / 5.0) < nrowsred) && (nrowsred < (nrows / 2.0)) && (ncolsred > (ncols / 2.0)) && (nodedepth < 5)) {
	if (0) {

		CPXENVptr envred = CPXopenCPLEX(&status);
		if (status) { perror("Can't use CPXopenCPLEX"); exit(status); }
		CPXLPptr lpred = CPXcreateprob(envred, &status, "red");
		if (status) { perror("Can't use CPXcreateprob"); exit(status); }

		onecolscnt = 0;
		for (j = 0; j < ncolsred; ++j) {
			if (inst->costs[colsred2sc[j]] <= 1.0) {
				status = CPXnewcols(envred, lpred, 1, &inst->costs[colsred2sc[j]], NULL, NULL, &integer, NULL);
				if (status) { perror("Can't use CPXnewcols"); exit(status); }
				onecolscnt++;
			} else {
				status = CPXnewcols(envred, lpred, 1, &inst->costs[colsred2sc[j]], NULL, &zero, &integer, NULL);
				if (status) { perror("Can't use CPXnewcols"); exit(status); }
			}
		}

		senses = (char *) malloc(nrowsred * sizeof(char));
		for (i = 0; i < nrowsred; ++i) {
			senses[i] = 'G';
		}

		status = CPXaddrows(envred, lpred, 0, nrowsred, numnzred, rmatval, senses, rmatbegred, rmatindred, rmatval, NULL, NULL);
		if (status) { perror("Can't use CPXaddrows"); }
		free(senses);

		CPXsetintparam(envred, CPXPARAM_ParamDisplay, 0);
		CPXsetintparam(envred, CPXPARAM_MIP_Display, 0);
		CPXsetintparam(envred, CPX_PARAM_REPEATPRESOLVE, 0);

		for (j = onecolscnt; j < ncolsred; ++j) {

			// Change rhs of all rows
			for (i = 0; i < nrowsred; ++i) {
				ind = i;
				status = CPXgetcoef(envred, lpred, i, j, &val);
				if (status) { perror("Can't use CPXgetcoef"); }
				status = CPXchgrhs(envred, lpred, 1, &ind, &val);
				if (status) { perror("Can't use CPXchgrhs"); }
			}

			// optimize subproblem
			status = CPXmipopt(envred, lpred);
			if (status) { perror("Can't use CPXmipopt"); }

			status = CPXgetobjval(envred, lpred, &objvalred);
			//if (status) { perror("Can't use CPXgetobjval"); }

			status = CPXgetstat(envred, lpred);
			//printf("stat=%d, obj=%5.1lf, costs=%5.1lf\n", status, objvalred, inst->costs[colsred2sc[j]]);

			// Infeasible solution OR
			// optimal solution found but obj val > column cost
			// THEN add this column to the model
			if (status == CPXMIP_INFEASIBLE || (status == CPXMIP_OPTIMAL && (objvalred > inst->costs[colsred2sc[j]]))) {
				// reset upper bound of the column to infinity
				ind = j;
				CPXchgbds(envred, lpred, 1, &ind, &up, NULL);
			}

			// column not selected
			else {
				removedcolscnt++;
				removedcols[j] = 1;
			}
		}

		//printf("cpxnrows=%4d, cpxncols=%4d, nrowsred=%4d, ncolsred=%4d, depth=%3d, removed cnt=%4d\n", CPXgetnumrows(env, lp), CPXgetnumcols(env, lp), nrowsred, ncolsred, nodedepth, removedcolscnt);

		CPXfreeprob(envred, &lpred);
		CPXcloseCPLEX(&envred);
		status = 0;
	}

	free(rmatval);

	brcol = -1;
	max = 0;
	for (j = 0; j < ncolsred; ++j) {
		if ((removedcols[j] == 0) && (max < colscnt[j])) {
			max = colscnt[j];
			brcol = j;
		}
	}

	char *nulllu;
	int *nullind;
	double *nullval;

	nulllu = (char *) malloc((removedcolscnt + 1) * sizeof(char));
	nullind = (int *) malloc((removedcolscnt + 1) * sizeof(int));
	nullval = (double *) malloc((removedcolscnt + 1) * sizeof(double));

	// UP NODE: upper bounds to 0
	k = 0;
	cond = 1;
	for (j = 0; j < ncolsred; ++j) {
		if ((j >= brcol) && cond) {
			nulllu[k] = 'U';
			nullind[k] = colsred2sc[brcol];
			nullval[k] = 0.0;
			k++;
			cond = 0;
		}

		if (removedcols[j]) {
			nulllu[k] = 'U';
			nullind[k] = colsred2sc[j];
			nullval[k] = 0.0;
			k++;
		}
	}

	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, k, nullind, nulllu, nullval, objval, NULL, &seqnum1);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds 1"); }

	// DOWN NODE: upper bounds to 0 and brcol to 1
	k = 0;
	cond = 1;
	for (j = 0; j < ncolsred; ++j) {
		if ((j >= brcol) && cond) {
			nulllu[k] = 'L';
			nullind[k] = colsred2sc[brcol];
			nullval[k] = 1.0;
			k++;
			cond = 0;
		}

		if (removedcols[j]) {
			nulllu[k] = 'U';
			nullind[k] = colsred2sc[j];
			nullval[k] = 0.0;
			k++;
		}
	}

	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, k, nullind, nulllu, nullval, objval, NULL, &seqnum2);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds 2"); }

	*useraction_p = CPX_CALLBACK_SET;

	free(nulllu);
	free(nullind);
	free(nullval);

	free(removedcols);
	free(colscnt);

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);

	free(rmatbegred);
	free(rmatindred);

	TERMINATE:
	return status;
} /* END SCcallbackbranchmaxcoldom */


int CPXPUBLIC SCcallbackbalasbranchrule1v1(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	unsigned char *mat, *qmat, *qmatred, *xchr;
	int ncolsred, nrowsred, ncols, nrows, p, nodedepth, status = 0, seqnum;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *xsupp;
	unsigned int i, j, xsupplen;
	double objval, zu, y;
	double *dj, *pi, *obj;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	int (*branchfunc)(CPXCENVptr env, void *cbdata, int wherefrom, double objval, unsigned char *qmat,
			int p, int ncols, int q, SCinstance *inst) = &SCmakebalasbranchrule1v1;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// NOT ROOT NODE
	if (nodedepth > 0) {
		SCnodedata *data = (SCnodedata *) malloc(sizeof(SCnodedata));

		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &data);
		if (status) { perror("Can't get node userhandle"); }

		// Get current node seqnum
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &seqnum);
		if (status) { perror("Can't get node depth"); }

		if ((data->p > 0) && (data->q > 0) && (data->seqnum == seqnum)) {
			branchfunc(env, cbdata, wherefrom, objval, data->qmat, data->p, data->ncols, data->q, inst);
			*useraction_p = CPX_CALLBACK_SET;
			free(data);
			goto TERMINATE;
		}
	}

	// Get reduced matrix and objective
	nrows = CPXgetnumrows(env, lp);
	ncols = CPXgetnumcols(env, lp);
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	mat = SCgetmatrix(env, lp, rowsred2sc, rowssc2red, colsred2sc, colssc2red, &nrowsred, &ncolsred);
	obj = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		obj[j] = inst->costs[colsred2sc[j]];
	}

	// Primal
	xsupplen = 0;
	xsupp = (int *) malloc(ncolsred * sizeof(int));
	xchr = (unsigned char *) calloc(ncolsred, sizeof(unsigned char));

	zu = SCbalasheurprimal0(mat, inst->costs, nrowsred, ncolsred, xchr, xsupp, &xsupplen, 3);

	double best_integer;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &best_integer);
	//printf("best int = %lf, zu =%lf\n", best_integer, zu);

	if (zu <= 0 || zu >= best_integer) {
		free(colsred2sc);
		free(colssc2red);
		free(mat);
		free(obj);
		free(xsupp);
		free(xchr);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;

		goto TERMINATE;
	}

	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (j = 0; j < ncolsred; ++j) {
		dj[j] = obj[j];
	}

	SCbalasheurdual1(mat, nrowsred, ncolsred, xchr, pi, dj);

	status = SCbalasheurdual3(mat, nrowsred, ncolsred, xchr, pi, dj, zu);
	if (status < 0) {
		free(colsred2sc);
		free(colssc2red);
		free(mat);
		free(obj);
		free(xsupp);
		free(xchr);
		free(dj);
		free(pi);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;

		goto TERMINATE;
	}

	qmatred = (unsigned char *) calloc(sizeof(unsigned char), (inst->SC_BALAS_MAX_BRANCH) * ncolsred);
	qmat = (unsigned char *) calloc(sizeof(unsigned char), (inst->SC_BALAS_MAX_BRANCH) * ncols);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	p = SCbalasbranchrule1(mat, nrowsred, ncolsred, xchr, dj, zu, y, qmatred, inst->SC_BALAS_MAX_BRANCH, inst->SC_BALAS_MAX_SINGL);

	if (p > 1) {

		for (i = 0; i < p; ++i) {
			for (j = 0; j < ncolsred; ++j) {
				if (qmatred[i*ncolsred + j]) {
					qmat[i*ncols + colsred2sc[j]] = 1;
				}
			}
		}

		branchfunc(env, cbdata, wherefrom, objval, qmat, p, ncols, 0, inst);
	} else {
		free(qmat);
		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
	}

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);
	free(mat);
	free(obj);
	free(xsupp);
	free(xchr);
	free(dj);
	free(pi);
	free(qmatred);

	TERMINATE:
	return status;
} /* END SCcallbackbalasbranchrule1 */


int CPXPUBLIC SCcallbackbalasbranchrule1_test(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	unsigned char *mat, *qmat, *qmatred, *xchr;
	int off, ncolsred, nrowsred, ncols, nrows, p, nodedepth, status = 0, seqnum;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *xsupp;
	unsigned int i, j, k, xsupplen;
	double objval, zu, y;
	double *dj, *pi, *obj;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	int (*branchfunc)(CPXCENVptr env, void *cbdata, int wherefrom, double objval, unsigned char *qmat,
			int p, int ncols, int q, SCinstance *inst) = &SCmakebalasbranchrule1v1;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// NOT ROOT NODE
	if (nodedepth > 0) {
		SCnodedata *data = (SCnodedata *) malloc(sizeof(SCnodedata));

		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &data);
		if (status) { perror("Can't get node userhandle"); }

		// Get current node seqnum
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &seqnum);
		if (status) { perror("Can't get node depth"); }

		if ((data->p > 0) && (data->q > 0) && (data->seqnum == seqnum)) {
			branchfunc(env, cbdata, wherefrom, objval, data->qmat, data->p, data->ncols, data->q, inst);
			*useraction_p = CPX_CALLBACK_SET;
			free(data);
			goto TERMINATE;
		}
		//free(data);
	}

	// Get reduced matrix and objective
	nrows = CPXgetnumrows(env, lp);
	ncols = CPXgetnumcols(env, lp);
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));
	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	mat = SCgetmatrix(env, lp, rowsred2sc, rowssc2red, colsred2sc, colssc2red, &nrowsred, &ncolsred);
	obj = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		obj[j] = inst->costs[colsred2sc[j]];
	}

	// Primal
	xchr = (unsigned char *) calloc(ncolsred, sizeof(unsigned char));

	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (i = 0; i < nrowsred; ++i) {
		CPXgetpi(env, lp, &pi[i], rowsred2sc[i], rowsred2sc[i]);
	}

	// Copy obj in dj
	memcpy(dj, obj, ncolsred * sizeof(double));
	//cblas_ccopy(ncolsred, obj, 1, dj, 1);

	for (i = 0; i < nrowsred; ++i) {
		off = i * ncolsred;
		for (j = 0; j < ncolsred; ++j) {
			dj[j] -= mat[off + j] * pi[i];
		}
	}

	//cblas_dgemv(CblasRowMajor, CblasTrans, nrowsred, ncolsred, -1.0, mat, ncolsred, pi, 1, 1.0, dj, 1);

	// Make dual feasible
	int l;
	for (j = 0; j < ncolsred; ++j) {
		for (i = 0; i < nrowsred; ++i) {
			if (dj[j] >= SC_EPSILON_SMALL) {
				break;
			}

			if (mat[i * ncolsred + j]) {
				if (dj[j] < -pi[i]) {
					for (l = 0; l < ncolsred; ++l) {
						dj[l] += pi[i] * mat[i * ncolsred + l];
					}
				} else {
					pi[i] += dj[j];
					for (l = 0; l < ncolsred; ++l) {
						dj[l] -= dj[j] * mat[i * ncolsred + l];
					}
				}
			}
		}
	}

	//zu = SCbalasheurprimal5b(mat, obj, nrowsred, ncolsred, xchr, pi, dj);
	xsupp = (int *) malloc(ncolsred * sizeof(int));
	xsupplen = 0;
	zu = SCbalasheurprimal0(mat, obj, nrowsred, ncolsred, xchr, xsupp, &xsupplen, 3);
	free(xsupp);

	SCbalasheurdual3(mat, nrowsred, ncolsred, xchr, pi, dj, zu);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	//printf("zl=%lf, zu=%lf\n", y, zu);

	//if (!SCiscover(mat, nrowsred, ncolsred, xchr)) printf("NOT COVER.\n");
	//if (!SCisdualsolution(mat, obj, nrowsred, ncolsred, pi)) printf("NOT DUAL SOLUTION.\n");

	if ((zu <= y) || (!SCisdualsolution(mat, obj, nrowsred, ncolsred, pi))) {
		free(colsred2sc);
		free(colssc2red);
		free(rowsred2sc);
		free(rowssc2red);
		free(mat);
		free(obj);
		free(xchr);
		free(dj);
		free(pi);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;
		status = 0;

		goto TERMINATE;
	}

	qmatred = (unsigned char *) calloc(sizeof(unsigned char), (inst->SC_BALAS_MAX_BRANCH) * ncolsred);
	qmat = (unsigned char *) calloc(sizeof(unsigned char), (inst->SC_BALAS_MAX_BRANCH) * ncols);

	//p = SCbalasbranchrule1(mat, nrowsred, ncolsred, xchr, dj, zu, y, qmatred, inst->SC_BALAS_MAX_BRANCH, inst->SC_BALAS_MAX_SINGL);
	p = SCbalasbranchrule1_test(mat, nrowsred, ncolsred, xchr, dj, zu, y, qmatred, inst->SC_BALAS_MAX_BRANCH, inst->SC_BALAS_MAX_SINGL);

	if (p > 1) {

		for (i = 0; i < p; ++i) {
			for (j = 0; j < ncolsred; ++j) {
				if (qmatred[i*ncolsred + j]) {
					qmat[i*ncols + colsred2sc[j]] = 1;
				}
			}
		}

		//printf("p=%d\n", p);

		branchfunc(env, cbdata, wherefrom, objval, qmat, p, ncols, 0, inst);
	} else {
		free(qmat);
		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		//SCmakebalasbranchrule2(env, cbdata, wherefrom, objval, mat, nrowsred, ncolsred, colsred2sc, inst);
	}

	*useraction_p = CPX_CALLBACK_SET;

	free(colsred2sc);//printf("ok00\n");
	free(colssc2red);//printf("ok01\n");
	free(rowsred2sc);//printf("ok02\n");
	free(rowssc2red);//printf("ok03\n");
	free(mat);//printf("ok04\n");
	free(obj);//printf("ok05\n");
	free(xchr);//printf("ok06\n");
	free(dj);//printf("ok07\n");
	free(pi);//printf("ok08\n");
	free(qmatred);//printf("ok09\n");

	TERMINATE:
	return status;
} /* END SCcallbackbalasbranchrule1_test */


int CPXPUBLIC SCcallbackbalasbranchrule1_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char flag;
	int ncolsred, nrowsred, ncols, nrows, p, nodedepth, status = 0, seqnum, numnz, numnzred, nzcnt_tmp, surplus, cnt;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *item, *xind;
	int *rmatind, *rmatindred, *rmatbeg, *rmatbegred, *cmatbegred, *cmatindred, *rqbeg, *rqind;
	unsigned int i, j, k, xindlen;
	double objval, zu, y;
	double *dj, *pi, *objred, *ub, *lb, *rmatval, *incumbent;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	int (*branchfunc)(CPXCENVptr env, void *cbdata, int wherefrom, double objval, int *rqbeg, int *rqind,
			int p, int ncols, int q, SCinstance *inst) = &SCmakebalasbranchrule1v1_sparse;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// NOT ROOT NODE
	if (nodedepth > 0) {
		SCnodedata *data = (SCnodedata *) malloc(sizeof(SCnodedata));

		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &data);
		if (status) { perror("Can't get node userhandle"); }

		// Get current node seqnum
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &seqnum);
		if (status) { perror("Can't get node depth"); }

		if ((data->p > 0) && (data->q > 0) && (data->seqnum == seqnum)) {
			branchfunc(env, cbdata, wherefrom, objval, data->rqbeg, data->rqind, data->p, data->ncols, data->q, inst);
			*useraction_p = CPX_CALLBACK_SET;
			free(data);
			goto TERMINATE;
		}
		//free(data);
	}


	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
	rmatindred = (int *) malloc(numnz * sizeof(int));

	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }

	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;

	nrowsred = 0;
	numnzred = 0;
	rmatbeg[nrows] = numnz;
	for (i = 0; i < nrows; ++i) {

		flag = 1;
		rmatbegred[nrowsred] = numnzred;
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				flag = 0;
				numnzred = rmatbegred[nrowsred];
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				rmatindred[numnzred] = colssc2red[rmatind[j]];
				numnzred++;
			}
		}

		if (flag) {
			rowssc2red[i] = nrowsred;
			rowsred2sc[nrowsred] = i;
			nrowsred++;
		} else {
			rowssc2red[i] = -1;
		}
	}
	rmatbegred[nrowsred] = numnzred;

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	objred = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		objred[j] = inst->costs[colsred2sc[j]];
	}

	cmatbegred = (int *) malloc((ncolsred + 1) * sizeof(int));
	cmatindred = (int *) malloc(numnzred * sizeof(int));

	cnt = 0;
	for (j = 0; j < ncolsred; ++j) {
		cmatbegred[j] = cnt;
		for (i = 0; i < nrowsred; ++i) {
			item = (int *) bsearch(&j, &rmatindred[rmatbegred[i]], rmatbegred[i+1] - rmatbegred[i], sizeof(int), SCint_cmp);
			if (item != NULL) {
				cmatindred[cnt] = i;
				cnt++;
			}
		}
	}
	cmatbegred[ncolsred] = cnt;


	// Primal
	xindlen = 0;
	xind = (int *) malloc(ncolsred * sizeof(int));
	incumbent = (double *) malloc(ncols * sizeof(double));
	CPXgetcallbackincumbent(env, cbdata, wherefrom, incumbent, 0, ncols - 1);

	for (j = 0; j < ncols; ++j) {
		if (incumbent[j] > 0.5 && colssc2red[j] > -1) {
			xind[xindlen] = colssc2red[j];
			xindlen++;
		}
	}
	free(incumbent);

	zu = SCbalasheurprimal0_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, objred, nrowsred, ncolsred, xind, &xindlen, 3);

	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (i = 0; i < nrowsred; ++i) {
		CPXgetpi(env, lp, &pi[i], rowsred2sc[i], rowsred2sc[i]);
	}

	// Check if is dual solution
	/*double acc;
	for (j = 0; j < ncolsred; ++j) {
		acc = 0.0;
		for (i = cmatbegred[j]; i < cmatbegred[j+1]; ++i) {
			acc += pi[cmatindred[i]];
		}

		if ((acc - SC_EPSILON_SMALL) > objred[j]) {
			printf("NOT DUAL SOLUTION. acc=%lf, obj=%lf\n", acc, objred[j]);
			break;
		}
	}*/

	for (j = 0; j < ncolsred; ++j) {
		dj[j] = objred[j];
	}

	for (i = 0; i < nrowsred; ++i) {
		for (j = rmatbegred[i]; j < rmatbegred[i+1]; ++j) {
			dj[rmatindred[j]] -= pi[i];
		}
	}

	/*for (j = 0; j < ncolsred; ++j) {
		if (dj[j] < -SC_EPSILON_SMALL) {
			printf("NOT VALID REDUCED COST. dj=%lf, j=%d\n", dj[j], j);
			break;
		}
	}*/

	//SCbalasheurdual1_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj);
	//SCdual0_sparse(rmatbegred, rmatindred, objred, nrowsred, ncolsred, pi, dj);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	if ((y < 0) || (abs(y) == NAN) || (zu <= y)) {
		free(xind);
		free(dj);
		free(pi);

		free(colsred2sc);
		free(colssc2red);
		free(rowsred2sc);
		free(rowssc2red);

		free(rmatbegred);
		free(rmatindred);

		free(cmatbegred);
		free(cmatindred);

		free(objred);

		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
		*useraction_p = CPX_CALLBACK_SET;
		status = 0;

		goto TERMINATE;
	}

	SCbalasheurdual3_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj, zu);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	rqbeg = (int *) malloc((inst->SC_BALAS_MAX_BRANCH + 1) * sizeof(int));
	rqind = (int *) malloc(numnzred * sizeof(int));
	p = SCbalasbranchrule1_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, nrowsred, ncolsred, xind, xindlen, dj,
			zu, y, rqbeg, rqind, inst->SC_BALAS_MAX_BRANCH, inst->SC_BALAS_MAX_SINGL);

	if (p > 1) {

		for (i = 0; i < p; ++i) {
			for (j = rqbeg[i]; j < rqbeg[i+1]; ++j) {
				rqind[j] = colsred2sc[rqind[j]];
			}
		}

		printf("p=%d\n", p);
		branchfunc(env, cbdata, wherefrom, objval, rqbeg, rqind, p, ncols, 0, inst);
	} else {
		free(rqbeg);
		free(rqind);
		SCmakecplexbranch(env, cbdata, wherefrom, ncols, inst);
	}

	*useraction_p = CPX_CALLBACK_SET;

	free(xind);
	free(dj);
	free(pi);

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);

	free(rmatbegred);
	free(rmatindred);

	free(cmatbegred);
	free(cmatindred);

	free(objred);

	TERMINATE:
	return status;
} /* END SCcallbackbalasbranchrule1_sparse */


int CPXPUBLIC SCcallbackbalasbranchrule1maxcol_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char flag;
	int ncolsred, nrowsred, ncols, nrows, p, nodedepth, status = 0, seqnum, numnz, numnzred, nzcnt_tmp, surplus, cnt;
	int *rowsred2sc, *rowssc2red, *colsred2sc, *colssc2red, *item, *xind;
	int *rmatind, *rmatindred, *rmatbeg, *rmatbegred, *cmatbegred, *cmatindred, *rqbeg, *rqind, *colscnt;
	unsigned int i, j, k, xindlen;
	double objval, zu, y;
	double *dj, *pi, *objred, *ub, *lb, *rmatval, *incumbent;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	int (*branchfunc)(CPXCENVptr env, void *cbdata, int wherefrom, double objval, int *rqbeg, int *rqind,
					  int p, int ncols, int q, SCinstance *inst) = &SCmakebalasbranchrule1v1_sparse;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get current node depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &nodedepth);
	if (status) { perror("Can't get node depth"); }

	// NOT ROOT NODE
	if (nodedepth > 0) {
		SCnodedata *data = (SCnodedata *) malloc(sizeof(SCnodedata));

		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_USERHANDLE, &data);
		if (status) { perror("Can't get node userhandle"); }

		// Get current node seqnum
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &seqnum);
		if (status) { perror("Can't get node depth"); }

		if ((data->p > 0) && (data->q > 0) && (data->seqnum == seqnum)) {
			branchfunc(env, cbdata, wherefrom, objval, data->rqbeg, data->rqind, data->p, data->ncols, data->q, inst);
			*useraction_p = CPX_CALLBACK_SET;
			free(data);
			goto TERMINATE;
		}
		//free(data);
	}


	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
	rmatindred = (int *) malloc(numnz * sizeof(int));

	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));
	rowsred2sc = (int *) malloc(nrows * sizeof(int));
	rowssc2red = (int *) malloc(nrows * sizeof(int));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }

	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;
	colscnt = (int *) calloc(ncolsred, sizeof(int));

	nrowsred = 0;
	numnzred = 0;
	rmatbeg[nrows] = numnz;
	for (i = 0; i < nrows; ++i) {

		flag = 1;
		rmatbegred[nrowsred] = numnzred;
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				flag = 0;
				numnzred = rmatbegred[nrowsred];
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				rmatindred[numnzred] = colssc2red[rmatind[j]];
				numnzred++;
			}
		}

		if (flag) {

			for (j = rmatbeg[i]; j < rmatbeg[i + 1]; ++j) {
				colscnt[colssc2red[rmatind[j]]] += (colssc2red[rmatind[j]] > -1);
			}

			rowssc2red[i] = nrowsred;
			rowsred2sc[nrowsred] = i;
			nrowsred++;
		} else {
			rowssc2red[i] = -1;
		}
	}
	rmatbegred[nrowsred] = numnzred;

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	objred = (double *) malloc(ncolsred * sizeof(double));
	for (j = 0; j < ncolsred; ++j) {
		objred[j] = inst->costs[colsred2sc[j]];
	}

	cmatbegred = (int *) malloc((ncolsred + 1) * sizeof(int));
	cmatindred = (int *) malloc(numnzred * sizeof(int));

	cnt = 0;
	for (j = 0; j < ncolsred; ++j) {
		cmatbegred[j] = cnt;
		for (i = 0; i < nrowsred; ++i) {
			item = (int *) bsearch(&j, &rmatindred[rmatbegred[i]], rmatbegred[i+1] - rmatbegred[i], sizeof(int), SCint_cmp);
			if (item != NULL) {
				cmatindred[cnt] = i;
				cnt++;
			}
		}
	}
	cmatbegred[ncolsred] = cnt;


	// Primal
	xindlen = 0;
	xind = (int *) malloc(ncolsred * sizeof(int));
	incumbent = (double *) malloc(ncols * sizeof(double));
	CPXgetcallbackincumbent(env, cbdata, wherefrom, incumbent, 0, ncols - 1);

	for (j = 0; j < ncols; ++j) {
		if (incumbent[j] > 0.5 && colssc2red[j] > -1) {
			xind[xindlen] = colssc2red[j];
			xindlen++;
		}
	}
	free(incumbent);

	zu = SCbalasheurprimal0_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, objred, nrowsred, ncolsred, xind, &xindlen, 3);


	// Dual and reduced costs
	dj = (double *) malloc(ncolsred * sizeof(double));
	pi = (double *) malloc(nrowsred * sizeof(double));

	for (j = 0; j < ncolsred; ++j) {
		dj[j] = objred[j];
	}

	SCbalasheurdual1_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	if (zu <= y) {
		free(xind);
		free(dj);
		free(pi);

		free(colsred2sc);
		free(colssc2red);
		free(rowsred2sc);
		free(rowssc2red);

		free(rmatbegred);
		free(rmatindred);

		free(cmatbegred);
		free(cmatindred);

		free(objred);

		int brcol = -1;
		int max = 0;
		for (j = 0; j < ncolsred; ++j) {
			if (max < colscnt[j]) {
				max = colscnt[j];
				brcol = j;
			}
		}

		char nulllu[1] = {'U'};
		int nullind[1] = { colsred2sc[brcol] };
		double nullval[1] = {0.0};

		SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data1->ncols = ncols;
		data1->p = 0;
		data1->q = 0;

		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds (1)"); }

		nulllu[0] = 'L';
		nullval[0] = 1.0;

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = 0;
		data2->q = 0;

		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds (2)"); }

		free(colscnt);

		*useraction_p = CPX_CALLBACK_SET;
		status = 0;

		goto TERMINATE;
	}

	SCbalasheurdual3_sparse(rmatbegred, rmatindred, nrowsred, xind, xindlen, pi, dj, zu);

	y = 0;
	for (i = 0; i < nrowsred; ++i) {
		y += pi[i];
	}

	rqbeg = (int *) malloc((inst->SC_BALAS_MAX_BRANCH + 1) * sizeof(int));
	rqind = (int *) malloc(numnzred * sizeof(int));
	p = SCbalasbranchrule1_sparse(rmatbegred, rmatindred, cmatbegred, cmatindred, nrowsred, ncolsred, xind, xindlen, dj,
			zu, y, rqbeg, rqind, inst->SC_BALAS_MAX_BRANCH, inst->SC_BALAS_MAX_SINGL);

	if (p > 1) {

		for (i = 0; i < p; ++i) {
			for (j = rqbeg[i]; j < rqbeg[i+1]; ++j) {
				rqind[j] = colsred2sc[rqind[j]];
			}
		}

		printf("p=%d\n", p);
		branchfunc(env, cbdata, wherefrom, objval, rqbeg, rqind, p, ncols, 0, inst);
	} else {
		free(rqbeg);
		free(rqind);

		int brcol = -1;
		int max = 0;
		for (j = 0; j < ncolsred; ++j) {
			if (max < colscnt[j]) {
				max = colscnt[j];
				brcol = j;
			}
		}

		char nulllu[1] = {'U'};
		int nullind[1] = { colsred2sc[brcol] };
		double nullval[1] = {0.0};

		SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data1->ncols = ncols;
		data1->p = 0;
		data1->q = 0;

		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds (1)"); }

		nulllu[0] = 'L';
		nullval[0] = 1.0;

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = 0;
		data2->q = 0;

		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, nullind, nulllu, nullval, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds (2)"); }
	}

	*useraction_p = CPX_CALLBACK_SET;

	free(xind);
	free(dj);
	free(pi);

	free(colscnt);

	free(colsred2sc);
	free(colssc2red);
	free(rowsred2sc);
	free(rowssc2red);

	free(rmatbegred);
	free(rmatindred);

	free(cmatbegred);
	free(cmatindred);

	free(objred);

	TERMINATE:
	return status;
} /* END SCcallbackbalasbranchrule1maxcol_sparse */


int CPXPUBLIC SCcallbackbalasbranchrule2(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p) {

	// Silence compiler warning(s) about unused parameter(s).
	(void) brtype;
	(void) sos;
	(void) nodecnt;
	(void) bdcnt;
	(void) nodebeg;
	(void) indices;
	(void) lu;
	(void) bd;
	(void) nodeest;

	char flag;
	char *varlu;
	int ncolsred, nrowsred, ncols, nrows, status = 0, numnz, numnzred, nzcnt_tmp, surplus, cnt;
	int min, max, it, ht, nbrvars, seqnum1, seqnum2;
	int *colsred2sc, *colssc2red;
	int *rmatind, *rmatindred, *rmatbeg, *rmatbegred, *varind;
	unsigned int i, j, k;
	double objval;
	double *ub, *lb, *rmatval, *varbd;

	SCinstance *inst = cbhandle;
	CPXLPptr lp;

	// Initialize useraction to indicate no user action taken
	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Get current node lp and other relative data (ncols, nrows, number of non-zero elements)
	status = CPXgetcallbacknodelp(env, cbdata, wherefrom, &lp);
	if (status) { perror("Can't get node lp"); }

	// Get objval
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	if (status) { perror("Can't get the node objval"); }

	// Get reduced matrix and objective
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc((nrows + 1) * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	rmatbegred = (int *) malloc((nrows + 1) * sizeof(int));
	rmatindred = (int *) malloc(numnz * sizeof(int));

	colsred2sc = (int *) malloc(ncols * sizeof(int));
	colssc2red = (int *) malloc(ncols * sizeof(int));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }

	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	ncolsred = cnt;

	nrowsred = 0;
	numnzred = 0;
	rmatbeg[nrows] = numnz;
	for (i = 0; i < nrows; ++i) {

		flag = 1;
		rmatbegred[nrowsred] = numnzred;
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				flag = 0;
				numnzred = rmatbegred[nrowsred];
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				rmatindred[numnzred] = colssc2red[rmatind[j]];
				numnzred++;
			}
		}

		if (flag) {
			nrowsred++;
		}
	}
	rmatbegred[nrowsred] = numnzred;

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	max = 0;
	it = -1;
	for (i = 0; i < nrowsred; ++i) {
		if ((rmatbegred[i + 1] - rmatbegred[i]) > max) {
			max = rmatbegred[i + 1] - rmatbegred[i];
			it = i;
		}
	}

	min = ncolsred;
	ht = -1;
	for (i = 0; i < nrowsred; ++i) {

		if (i == it) {
			continue;
		}

		cnt = 0;
		j = rmatbegred[it];
		k = rmatbegred[i];
		while (j < rmatbegred[it + 1] && k < rmatbegred[i+1]) {
			if (rmatindred[j] == rmatindred[k]) {
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatindred[j] < rmatindred[k]) {
				j++;
				continue;
			}

			if (rmatindred[j] > rmatindred[k]) {
				k++;
				continue;
			}
		}

		if ((0 < cnt) && (cnt < min)) {
			min = cnt;
			ht = i;
		}

		if (min == 1) {
			break;
		}
	}

	nbrvars = min;
	varind = (int *) malloc(nbrvars * sizeof(int));
	varlu = (char *) malloc(nbrvars * sizeof(char));
	varbd = (double *) malloc(nbrvars * sizeof(double));

	cnt = 0;
	j = rmatbegred[it];
	k = rmatbegred[ht];
	while (j < rmatbegred[it+1] && k < rmatbegred[ht+1]) {
		if (rmatindred[j] == rmatindred[k]) {
			varind[cnt] = colsred2sc[rmatindred[j]];
			varlu[cnt] = 'U';
			varbd[cnt] = 0.0;

			cnt++;
			j++;
			k++;
			continue;
		}

		if (rmatindred[j] < rmatindred[k]) {
			j++;
			continue;
		}

		if (rmatindred[j] > rmatindred[k]) {
			k++;
			continue;
		}
	}

	// Check last row
	cnt = 0;
	j = 0;
	k = rmatbegred[nrowsred - 1];
	while (j < nbrvars && k < rmatbegred[nrowsred]) {
		if (varind[j] == colsred2sc[rmatindred[k]]) {
			cnt++;
			j++;
			k++;
			continue;
		}

		if (varind[j] < colsred2sc[rmatindred[k]]) {
			j++;
			continue;
		}

		if (varind[j] > colsred2sc[rmatindred[k]]) {
			k++;
			continue;
		}
	}

	// Avoid cycles: if vars are the same as las row, remove last var
	if (nbrvars > 1 && (cnt == nbrvars) && (cnt == (rmatbegred[nrowsred] - rmatbegred[nrowsred - 1]))) {
		nbrvars--;
	}

	//printf("\nnrows=%d, ncols=%d\n", nrows, ncols);
	//printf("ind = "); for (i = 0; i < nbrvars; ++i) printf("%d ", varind[i]); printf("\n");
	//printf("last row = "); for (j = 0; j < ncols; ++j) if (mat[(nrows-1)*ncols + j]) printf("%d ", j); printf("\n");

	// Up node
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, objval, NULL, &seqnum1);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds"); }


	// Down node
	for (i = 0; i < nbrvars; ++i) {
		varbd[i] = 1.0;
	}

	double rhs[1] = {1};
	char senses[1] = {'G'};
	int rind[1] = {0};

	status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, nbrvars, rhs, senses, rind, varind, varbd, objval, NULL, &seqnum2);
	if (status) { perror("Can't call CPXbranchcallbackbranchconstraints"); }

	*useraction_p = CPX_CALLBACK_SET;

	free(varind);
	free(varlu);
	free(varbd);

	free(colsred2sc);
	free(colssc2red);

	free(rmatbegred);
	free(rmatindred);

	return status;
} /* END SCcallbackbalasbranchrule2 */


unsigned char *SCgetmatrix(CPXCENVptr env, CPXLPptr lp, int *rowsred2sc,
		int *rowssc2red, int *colsred2sc, int *colssc2red, int *nrowsred,
		int *ncolsred) {

	unsigned char flag;
	unsigned char *mat, *matr;
	int cnt, nrows, ncols, numnz, status, nzcnt_tmp, surplus;
	int *rmatbeg, *rmatind;
	unsigned int i, j, k, l;
	double *ub, *lb, *rmatval;

	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	numnz = CPXgetnumnz(env, lp);

	// Get the upper bounds
	ub = (double *) malloc(ncols * sizeof(double));
	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) { perror("Can't get the ub"); }

	// Get the lower bounds
	lb = (double *) malloc(ncols * sizeof(double));
	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) { perror("Can't get the lb"); }

	// Get the matrix of the node lp
	rmatbeg = (int *) malloc(nrows * sizeof(int));
	rmatind = (int *) malloc(numnz * sizeof(int));
	rmatval = (double *) malloc(numnz * sizeof(double));

	status = CPXgetrows(env, lp, &nzcnt_tmp, rmatbeg, rmatind, rmatval, numnz, &surplus, 0, nrows - 1);
	if (status) { perror("Can't get matrix rows"); }


	// remove 1 vars rows and cols, remove 0 vars cols
	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		colssc2red[j] = (lb[j] < 0.5) && (ub[j] > 0.5) ? cnt: -1;
		colsred2sc[cnt] = j;
		cnt += (lb[j] < 0.5) && (ub[j] > 0.5);
	}
	*ncolsred = cnt;

	mat = (unsigned char *) calloc(sizeof(unsigned char), (nrows) * (*ncolsred));
	*nrowsred = 0;
	for (i = 0; i < nrows; ++i) {
		matr = &mat[(*nrowsred) * (*ncolsred)];
		k = (i == nrows - 1 ? numnz : rmatbeg[i + 1]);
		flag = 1;
		for (j = rmatbeg[i]; j < k; ++j) {

			// the row is covered: do not include it in reduced matrix
			// clean up the row and set flag to 0
			if (lb[rmatind[j]] > 0.5) {
				for (l = 0; l < *ncolsred; ++l) {
					matr[l] = 0;
				}
				flag = 0;
				break;
			}

			if (colssc2red[rmatind[j]] > -1) {
				matr[colssc2red[rmatind[j]]] = 1;
			}
		}

		if (flag) {
			rowssc2red[i] = *nrowsred;
			rowsred2sc[*nrowsred] = i;
			*nrowsred = *nrowsred + 1;
		} else {
			rowssc2red[i] = -1;
		}
	}

	free(ub);
	free(lb);

	free(rmatbeg);
	free(rmatind);
	free(rmatval);

	return mat;
}


int SCmakecplexbranch(CPXCENVptr env, void *cbdata, int wherefrom,
		int ncols, SCinstance *inst) {

	int status = 0;

	SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data1->ncols = ncols;
	data1->p = 0;
	data1->q = 0;

	status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 0, data1, &data1->seqnum);
	if (status) { perror("Can't make CPLEX branch 1"); status = CPXERR_CALLBACK;}
	inst->cplexnodecnt++;

	SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data2->ncols = ncols;
	data2->p = 0;
	data2->q = 0;

	status = CPXbranchcallbackbranchasCPLEX(env, cbdata, wherefrom, 1, data2, &data2->seqnum);
	if (status) { perror("Can't make CPLEX branch 2"); status = CPXERR_CALLBACK; }
	inst->cplexnodecnt++;

	return status;
}


int SCmakebalasbranchrule1v1(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *qmat, int p, int ncols, int q,
		SCinstance *inst) {
	char *varlu;
	unsigned int i, j, nbrvars, off;
	int * varind;
	double *varbd;

	int status = 0;

	varind = (int *) malloc(ncols * sizeof(int));
	varlu = (char *) malloc(ncols * sizeof(char));
	varbd = (double *) malloc(ncols * sizeof(double));

	nbrvars = 0;
	off = q * ncols;
	for (j = 0; j < ncols; ++j) {
		varind[nbrvars] = j * qmat[off + j];
		nbrvars += qmat[off + j];
	}

	// UP NODE: vars upper bound set to 0
	for (i = 0; i < nbrvars; ++i) {
		varlu[i] = 'U';
		varbd[i] = 0.0;
	}

	SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data1->ncols = ncols;
	data1->p = 0;
	data1->q = 0;

	if (q == 0) {
		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;
	} else {
		char *row = (char *) malloc(ncols * sizeof(char));
		int *rind = (int *) malloc(ncols * sizeof(int));
		double *rval = (double *) malloc(ncols * sizeof(double));

		int nzcnt = 0;

		double rhs[1] = {1.0};
		char sense[1] = {'G'};
		int rbeg[1] = {0};

		for (i = 0; i < q; ++i) {
			off = i * ncols;
			for (j = 0; j < ncols; ++j) {
				row[j] = row[j] || qmat[off + j];
			}
		}

		for (j = 0; j < ncols; ++j) {
			rind[nzcnt] = j * row[j];
			nzcnt += row[j];
		}

		status = CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, 1, nzcnt, rhs, sense, rbeg, rind, rval, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;

		free(row);
		free(rind);
		free(rval);
	}

	// DONW NODE: vars >= 1 constraint, varbd used as varval
	q++;
	if (p - 1 == q) {

		char *row = (char *) malloc(ncols * sizeof(char));
		int *rind = (int *) malloc(ncols * sizeof(int));
		double *rval = (double *) malloc(ncols * sizeof(double));

		int nzcnt = 0;

		double rhs[1] = {1.0};
		char sense[1] = {'G'};
		int rbeg[1] = {0};

		for (i = 0; i < q; ++i) {
			off = i * ncols;
			for (j = 0; j < ncols; ++j) {
				row[j] = row[j] || qmat[off + j];
			}
		}

		for (j = 0; j < ncols; ++j) {
			rind[nzcnt] = j * row[j];
			nzcnt += row[j];
		}

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = 0;
		data2->q = 0;

		status = CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, 1, nzcnt, rhs, sense, rbeg, rind, rval, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;

		free(row);
		free(rind);
		free(rval);

		// Only at this point it is possible to free the qmat!!!
		free(qmat);
	} else {
		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = p;
		data2->q = q;
		data2->qmat = qmat;

		status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 0, 0, NULL, NULL, NULL, NULL, NULL, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral"); status = CPXERR_CALLBACK; goto TERMINATE; }
	}


	TERMINATE:
	free(varind);
	free(varlu);
	free(varbd);

	return (status);
}


int SCmakebalasbranchrule1v2(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *qmat, int p, int ncols, int q,
		SCinstance *inst) {
	char *varlu;
	unsigned int i, j, nbrvars, off;
	int * varind;
	double *varbd;

	int status = 0;

	varind = (int *) malloc(ncols * sizeof(int));
	varlu = (char *) malloc(ncols * sizeof(char));
	varbd = (double *) malloc(ncols * sizeof(double));

	nbrvars = 0;
	off = q * ncols;
	for (j = 0; j < ncols; ++j) {
		varind[nbrvars] = j * qmat[off + j];
		nbrvars += qmat[off + j];
	}

	for (i = 0; i < nbrvars; ++i) {
		varlu[i] = 'U';
		varbd[i] = 0.0;
	}

	SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data1->ncols = ncols;
	data1->p = 0;
	data1->q = 0;

	inst->balasnodecnt++;
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, objval, data1, &data1->seqnum);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds"); status = CPXERR_CALLBACK; goto TERMINATE; }

	if (p - 2 == q) {

		nbrvars = 0;
		off = (q+1) * ncols;
		for (j = 0; j < ncols; ++j) {
			varind[nbrvars] = j * qmat[off + j];
			nbrvars += qmat[off + j];
		}

		for (i = 0; i < nbrvars; ++i) {
			varlu[i] = 'U';
			varbd[i] = 0.0;
		}

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = 0;
		data2->q = 0;

		inst->balasnodecnt++;
		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nbrvars, varind, varlu, varbd, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds"); status = CPXERR_CALLBACK; goto TERMINATE; }

		free(qmat);

	} else {

		for (i = 0; i < nbrvars; ++i) {
			varbd[i] = 1.0;
		}

		double rhs[1] = {1};
		char senses[1] = {'G'};
		int rind[1] = {0};

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = p;
		data2->q = q + 1;
		data2->qmat = qmat;

		status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, nbrvars, rhs, senses, rind, varind, varbd, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchconstraints"); status = CPXERR_CALLBACK; goto TERMINATE; }
	}

	TERMINATE:
	free(varind);
	free(varlu);
	free(varbd);

	return (status);
}


int SCmakebalasbranchrule2(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *mat, unsigned int nrows, unsigned int ncols,
		const int *colsred2sc, SCinstance *inst) {

	char cond, flag;
	char *varlu;
	unsigned char *matr, *mati;
	int it, ht, status;
	int *varindsc, *varindred;
	unsigned int cnt, nzrowcnt, i, j, max, min, nbrvars;
	double *varbd;

	max = 0;
	it = -1;
	for (i = 0; i < nrows; ++i) {

		cnt = 0;
		matr = &mat[i*ncols];
		for (j = 0; j < ncols; ++j) {
			cnt += matr[j];
		}

		if (cnt > max) {
			max = cnt;
			it = i;
		}
	}

	min = ncols;
	ht = -1;
	mati = &mat[it*ncols];
	for (i = 0; i < nrows; ++i) {

		cnt = 0;
		matr = &mat[i*ncols];
		for (j = 0; j < ncols; ++j) {
			cnt += matr[j] && mati[j];
		}

		if ((0 < cnt) && (cnt < min) && (i != it)) {
			min = cnt;
			ht = i;
		}

		if (min == 1) {
			break;
		}
	}

	varindsc = (int *) malloc(ncols * sizeof(int));
	varindred = (int *) malloc(ncols * sizeof(int));
	varlu = (char *) malloc(ncols * sizeof(char));
	varbd = (double *) malloc(ncols * sizeof(double));

	nbrvars = 0;
	nzrowcnt = 0;
	matr = &mat[ht*ncols];
	for (j = 0; j < ncols; ++j) {
		cond = mati[j] && matr[j];
		varindsc[nbrvars] = colsred2sc[j * cond];
		varindred[nbrvars] = j * cond;
		nbrvars += cond;
		nzrowcnt += mat[(nrows-1)*ncols + j];
	}

	// Avoid cycles: if vars are the same as las row, remove last var
	flag = 1;
	for (j = 0; j < nbrvars; ++j) {
		flag = flag && mat[(nrows-1)*ncols + varindred[j]];
	}

	if (flag && (nbrvars == nzrowcnt) && (nbrvars > 1)) {
		nbrvars--;
	}


	//printf("\nnrows=%d, ncols=%d\n", nrows, ncols);
	//printf("ind = "); for (i = 0; i < nbrvars; ++i) printf("%d ", varindsc[i]); printf("\n");
	//printf("last row = "); for (j = 0; j < ncols; ++j) if (mat[(nrows-1)*ncols + j]) printf("%d ", j); printf("\n");


	// Up node
	for (i = 0; i < nbrvars; ++i) {
		varlu[i] = 'U';
		varbd[i] = 0.0;
	}

	SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data1->ncols = ncols;
	data1->p = 0;
	data1->q = 0;

	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nbrvars, varindsc, varlu, varbd, objval, data1, &data1->seqnum);
	if (status) { perror("Can't call CPXbranchcallbackbranchbds"); }


	// Down node
	for (i = 0; i < nbrvars; ++i) {
		varbd[i] = 1.0;
	}

	double rhs[1] = {1};
	char senses[1] = {'G'};
	int rind[1] = {0};

	SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data2->ncols = ncols;
	data2->p = 0;
	data2->q = 0;

	status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, nbrvars, rhs, senses, rind, varindsc, varbd, objval, data2, &data2->seqnum);
	if (status) { perror("Can't call CPXbranchcallbackbranchconstraints"); }

	free(varindsc);
	free(varindred);
	free(varlu);
	free(varbd);

	return 0;
}


int SCmakebalasbranchrule1v1_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, int *rqbeg, int *rqind, int p, int ncols, int q,
		SCinstance *inst) {

	char *varlu;
	int nzcnt, status = 0;
	unsigned int i, j;
	double *varbd;

	double rhs[1] = {1.0};
	char sense[1] = {'G'};
	int rbeg[1] = {0};

	char *row = (char *) malloc(ncols * sizeof(char));
	int *rind = (int *) malloc(ncols * sizeof(int));
	double *rval = (double *) malloc(ncols * sizeof(double));

	varlu = (char *) malloc(ncols * sizeof(char));
	varbd = (double *) malloc(ncols * sizeof(double));

	// UP NODE: vars upper bound set to 0
	for (i = 0; i < ncols; ++i) {
		varlu[i] = 'U';
		varbd[i] = 0.0;
	}

	SCnodedata *data1 = (SCnodedata *) malloc(sizeof(SCnodedata));
	data1->ncols = ncols;
	data1->p = 0;
	data1->q = 0;

	if (q == 0) {
		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, rqbeg[q+1] - rqbeg[q], &rqind[q], varlu, varbd, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchbds"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;
	} else {

		for (j = 0; j < ncols; ++j) {
			row[j] = 0;
		}

		for (i = 0; i < q; ++i) {
			for (j = rqbeg[i]; j < rqbeg[i+1]; ++j) {
				row[rqind[j]] = 1;
			}
		}

		nzcnt = 0;
		for (j = 0; j < ncols; ++j) {
			if (row[j]) {
				rind[nzcnt] = j;
				rval[nzcnt] = 1.0;
				nzcnt++;
			}
		}

		status = CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, rqbeg[q+1] - rqbeg[q], &rqind[q], varlu, varbd, 1, nzcnt, rhs, sense, rbeg, rind, rval, objval, data1, &data1->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral (1)"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;
	}

	// DONW NODE: vars >= 1 constraint, varbd used as varval
	q++;
	if (p - 1 == q) {

		for (j = 0; j < ncols; ++j) {
			row[j] = 0;
		}

		for (i = 0; i < q; ++i) {
			for (j = rqbeg[i]; j < rqbeg[i+1]; ++j) {
				row[rqind[j]] = 1;
			}
		}

		nzcnt = 0;
		for (j = 0; j < ncols; ++j) {
			if (row[j]) {
				rind[nzcnt] = j;
				rval[nzcnt] = 1.0;
				nzcnt++;
			}
		}

		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = 0;
		data2->q = 0;

		status = CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, rqbeg[q+1] - rqbeg[q], &rqind[q], varlu, varbd, 1, nzcnt, rhs, sense, rbeg, rind, rval, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral (2)"); status = CPXERR_CALLBACK; goto TERMINATE; }
		inst->balasnodecnt++;

		// Only at this point it is possible to free the q stuff!
		free(rqbeg);
		free(rqind);
	} else {
		SCnodedata *data2 = (SCnodedata *) malloc(sizeof(SCnodedata));
		data2->ncols = ncols;
		data2->p = p;
		data2->q = q;
		data2->rqbeg = rqbeg;
		data2->rqind = rqind;

		status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 0, 0, NULL, NULL, NULL, NULL, NULL, objval, data2, &data2->seqnum);
		if (status) { perror("Can't call CPXbranchcallbackbranchgeneral (3)"); status = CPXERR_CALLBACK; goto TERMINATE; }
	}

	TERMINATE:
	free(varlu);
	free(varbd);

	free(row);
	free(rind);
	free(rval);

	return (status);
}
