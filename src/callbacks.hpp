
#ifndef SC_CALLBACKS_H
#define SC_CALLBACKS_H

#include "common.hpp"

int CPXPUBLIC callbacks_balas_usercuts_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
											 int *useraction_p);

int CPXPUBLIC callbacks_balassusercuts(CPXCENVptr env, void *cbdata, int wherefrom,
									  void *cbhandle, int *useraction_p);

int CPXPUBLIC callbacks_balas_usercuts_test(CPXCENVptr env, void *cbdata, int wherefrom,
										   void *cbhandle, int *useraction_p);

int CPXPUBLIC callbacks_branch_maxcol(CPXCENVptr env, void *cbdata,
									 int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
									 int bdcnt, const int *nodebeg, const int *indices, const char *lu,
									 const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_branch_maxcol2(CPXCENVptr env, void *cbdata,
									  int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
									  int bdcnt, const int *nodebeg, const int *indices, const char *lu,
									  const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_branch_maxcol_sparse(CPXCENVptr env, void *cbdata,
											int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
											int bdcnt, const int *nodebeg, const int *indices, const char *lu,
											const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_branch_maxcol_dom(CPXCENVptr env, void *cbdata,
										int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
										int bdcnt, const int *nodebeg, const int *indices, const char *lu,
										const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_balas_branch_rule1v1(CPXCENVptr env, void *cbdata,
										   int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
										   int bdcnt, const int *nodebeg, const int *indices, const char *lu,
										   const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_balas_branch_rule1_test(CPXCENVptr env, void *cbdata,
											  int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
											  int bdcnt, const int *nodebeg, const int *indices, const char *lu,
											  const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_balas_branch_rule1_sparse(CPXCENVptr env, void *cbdata,
												int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
												int bdcnt, const int *nodebeg, const int *indices, const char *lu,
												const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_balas_branch_rule1_maxcol_sparse(CPXCENVptr env, void *cbdata,
													  int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
													  int bdcnt, const int *nodebeg, const int *indices, const char *lu,
													  const double *bd, const double *nodeest, int *useraction_p);

int CPXPUBLIC callbacks_balas_branch_rule2(CPXCENVptr env, void *cbdata,
										   int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt,
										   int bdcnt, const int *nodebeg, const int *indices, const char *lu,
										   const double *bd, const double *nodeest, int *useraction_p);

char *callbacks_get_matrix(CPXCENVptr env, CPXLPptr lp, int *rowsred2sc,
						  int *rowssc2red, int *colsred2sc, int *colssc2red, int *nrowsred,
						  int *ncolsred);

int callbacks_make_cplex_branch(CPXCENVptr env, void *cbdata, int wherefrom,
							  int ncols, SCinstance *inst);

int callbacks_make_balas_branch_rule1v1(CPXCENVptr env, void *cbdata, int wherefrom,
									 double objval, char *qmat, int p, int ncols, int q,
									 SCinstance *inst);

int callbacks_make_balas_branch_rule1v2(CPXCENVptr env, void *cbdata, int wherefrom,
									 double objval, char *qmat, int p, int ncols, int q,
									 SCinstance *inst);

int callbacks_make_balas_branch_rule2(CPXCENVptr env, void *cbdata, int wherefrom,
								   double objval, char *mat, int nrows, int ncols,
								   const int *colsred2sc, SCinstance *inst);

int callbacks_make_balas_branch_rule1v1_sparse(CPXCENVptr env, void *cbdata, int wherefrom,
									double objval, int *rqbeg, int *rqind, int p, int ncols, int q,
									SCinstance *inst);

#endif //SC_CALLBACKS_H
