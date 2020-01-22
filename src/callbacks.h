//
// Created by macs on 19/06/19.
//

#ifndef SET_COVERING_CALLBACKS_H
#define SET_COVERING_CALLBACKS_H

#include "aux.h"

int CPXPUBLIC SCcallbackbalasusercuts_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
        int *useraction_p);

int CPXPUBLIC SCcallbackbalasusercuts(CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *useraction_p);

int CPXPUBLIC SCcallbackbalasusercuts_test(CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *useraction_p);

int CPXPUBLIC SCcallbackbranchmaxcol(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbranchmaxcol2(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbranchmaxcol_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbranchmaxcoldom(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbalasbranchrule1v1(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbalasbranchrule1_test(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbalasbranchrule1_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbalasbranchrule1maxcol_sparse(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

int CPXPUBLIC SCcallbackbalasbranchrule2(CPXCENVptr env, void *cbdata,
		int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
		int bdcnt, const int *nodebeg, const int *indices, const char *lu,
		const double *bd, const double * nodeest, int *useraction_p);

unsigned char *SCgetmatrix(CPXCENVptr env, CPXLPptr lp, int *rowsred2sc,
		int *rowssc2red, int *colsred2sc, int *colssc2red, int *nrowsred,
		int *ncolsred);

int SCmakecplexbranch(CPXCENVptr env, void *cbdata, int wherefrom,
		int ncols, SCinstance *inst);

int SCmakebalasbranchrule1v1(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *qmat, int p, int ncols, int q,
		SCinstance *inst);

int SCmakebalasbranchrule1v2(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *qmat, int p, int ncols, int q,
		SCinstance *inst);

int SCmakebalasbranchrule2(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, unsigned char *mat, unsigned int nrows, unsigned int ncols,
		const int *colsred2sc, SCinstance *inst);

int SCmakebalasbranchrule1v1_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
		double objval, int *rqbeg, int *rqind, int p, int ncols, int q,
		SCinstance *inst);

#endif //SET_COVERING_CALLBACKS_H
