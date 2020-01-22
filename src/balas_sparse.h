#ifndef SET_COVERING_BALAS_SPARSE_H
#define SET_COVERING_BALAS_SPARSE_H

#include <cpxconst.h>
#include "aux.h"

int SCprimecover_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg,
		const int *cmatind, unsigned int nrows, unsigned int ncols, int *xind,
		unsigned int xindlen);

double SCbalasheurprimal0_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
		const double *obj, unsigned int nrows, unsigned int ncols, int *xind, unsigned int *xindlen,
		int whichfunc);

int SCdual0_sparse(const int *rmatbeg, const int *rmatind, double *obj, unsigned int nrows,
		unsigned int ncols, double *u, double *s);

int SCbalasheurdual1_sparse(const int *rmatbeg, const int *rmatind, unsigned int nrows, const int *xind,
		unsigned int xindlen, double *u, double *s);

int SCbalasheurdual3_sparse(const int *rmatbeg, const int *rmatind, unsigned int nrows,
		const int *xind, unsigned int xindlen, double *u, double *s, double zu);

int SCbalasbranchrule1_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
		unsigned int nrows, unsigned int ncols, const int *xind, unsigned int xindlen, double *dj,
		double zu, double y, int *rqbeg, int *rqind, int max_branch, int max_singl);

#endif //SET_COVERING_BALAS_SPARSE_H
