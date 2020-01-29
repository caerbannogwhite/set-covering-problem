
#ifndef SC_BALAS_DENSE_H
#define SC_BALAS_DENSE_H

#include "common.hpp"

int SCbalascutprocedure(char *mat, int nrows, int ncols,
		const char *x, double *s, double zu, double y,
		int *cutind);

int SCdual0(const char *mat, double *obj, int nrows,
		int ncols, double *u, double *s);

int SCbalasbranchrule1(char *mat, int nrows, int ncols,
        const char *x, double *s, double zu, double y,
        char *qmat, int max_branch, int max_singl);

int SCbalasbranchrule1_test(char *mat, int nrows, int ncols,
		const char *x, double *s, double zu, double y,
		char *qmat, int max_branch, int max_singl);

int SCprimecover(const char *mat, int nrows,
        int ncols, char *x);

int SCoversatrows(char *mat, int nrows,
		int ncols, char *x, int *xsupp,
		int *xsupplen_p);

int SCiscover(char *mat, int nrows,
		int ncols, const char *x);

double SCbalasheurprimal0(char *mat, const double *costs,
		int nrows, int ncols,
		char *x, int *xsupp, int *xsupplen,
		int whichfunc);

double SCbalasheurprimal12(char *mat, const double *costs,
		int nrows, int ncols,
		char *x, int *xsupp, int *xsupplen_p);

double SCbalasheurprimal5b(char *mat, const double *costs,
		int nrows, int ncols, char *x,
		double *u, double *s);

int SCisdualsolution(const char *mat, const double *obj,
		int nrows, int ncols, double *u);

int SCbalasheurdual1(char *mat,
		int nrows, int ncols,
		const char *x, double *u, double *s);

int SCbalasheurdual3(char *mat,
		int nrows, int ncols,
		const char *x, double *u, double *s, double zu);

#endif //SC_BALAS_DENSE_H
