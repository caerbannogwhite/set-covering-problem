
#ifndef SC_H
#define SC_H

#include "common.hpp"

int sc_solver(SCinstance &inst);
int sc_solver_balas_rule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_balas_rule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_balas_rule1_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_balas_rule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_maxcol(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_maxcol2(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
int sc_solver_maxcol_dom(SCinstance &inst, CPXENVptr env, CPXLPptr lp);

STATUS sc_build_lp2raw(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS sc_build_raw2lp(SCinstance &inst, CPXENVptr env, CPXLPptr lp);

#endif //SC_H
