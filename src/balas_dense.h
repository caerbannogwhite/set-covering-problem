#ifndef SET_COVERING_BALAS_DENSE_H
#define SET_COVERING_BALAS_DENSE_H

#include <cpxconst.h>
#include "aux.h"

int SCbalascutprocedure(unsigned char *mat, int nrows, int ncols,
		unsigned const char *x, double *s, double zu, double y,
		int *cutind);

int SCdual0(unsigned const char *mat, double *obj, unsigned int nrows,
		unsigned int ncols, double *u, double *s);

int SCbalasbranchrule1(unsigned char *mat, int nrows, int ncols,
        unsigned const char *x, double *s, double zu, double y,
        unsigned char *qmat, int max_branch, int max_singl);

int SCbalasbranchrule1_test(unsigned char *mat, int nrows, int ncols,
		unsigned const char *x, double *s, double zu, double y,
		unsigned char *qmat, int max_branch, int max_singl);

int SCprimecover(const unsigned char *mat, unsigned int nrows,
        unsigned int ncols, unsigned char *x);

int SCoversatrows(unsigned char *mat, unsigned int nrows,
		unsigned int ncols, unsigned char *x, int *xsupp,
		unsigned int *xsupplen_p);

int SCiscover(unsigned char *mat, unsigned int nrows,
		unsigned int ncols, unsigned const char *x);

double SCbalasheurprimal0(unsigned char *mat, const double *costs,
		unsigned int nrows, unsigned int ncols,
		unsigned char *x, int *xsupp, unsigned int *xsupplen,
		unsigned int whichfunc);

double SCbalasheurprimal12(unsigned char *mat, const double *costs,
		unsigned int nrows, unsigned int ncols,
		unsigned char *x, int *xsupp, unsigned int *xsupplen_p);

double SCbalasheurprimal5b(unsigned char *mat, const double *costs,
		unsigned int nrows, unsigned int ncols, unsigned char *x,
		double *u, double *s);

int SCisdualsolution(const unsigned char *mat, const double *obj,
		unsigned int nrows, unsigned int ncols, double *u);

int SCbalasheurdual1(unsigned char *mat,
		unsigned int nrows, unsigned int ncols,
		unsigned const char *x, double *u, double *s);

int SCbalasheurdual3(unsigned char *mat,
		unsigned int nrows, unsigned int ncols,
		unsigned const char *x, double *u, double *s, double zu);

#endif //SET_COVERING_BALAS_DENSE_H
