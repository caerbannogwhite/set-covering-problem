
#ifndef SC_BALAS_SPARSE_H
#define SC_BALAS_SPARSE_H

#include "common.hpp"

int SCprimecover_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg,
						const int *cmatind, int nrows, int ncols, int *xind,
						int xindlen);

double SCbalasheurprimal0_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
								 const double *obj, int nrows, int ncols, int *xind, int *xindlen,
								 int whichfunc);

int SCdual0_sparse(const int *rmatbeg, const int *rmatind, double *obj, int nrows,
				   int ncols, double *u, double *s);

int SCbalasheurdual1_sparse(const int *rmatbeg, const int *rmatind, int nrows, const int *xind,
							int xindlen, double *u, double *s);

int SCbalasheurdual3_sparse(const int *rmatbeg, const int *rmatind, int nrows,
							const int *xind, int xindlen, double *u, double *s, double zu);

int SCbalasbranchrule1_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
							  int nrows, int ncols, const int *xind, int xindlen, double *dj,
							  double zu, double y, int *rqbeg, int *rqind, int max_branch, int max_singl);

#endif //SC_BALAS_SPARSE_H
