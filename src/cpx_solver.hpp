
#ifndef CPX_SOLVER_H
#define CPX_SOLVER_H

#include "cpx_common.hpp"

STATUS cpxsol(SCinstance &inst);
STATUS cpxsol_preproc_dominance(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_balas_rule1(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_balas_rule1_test(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_balas_rule1_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_balas_rule2(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_maxcol(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_maxcol2(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_maxcol_sparse(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_maxcol_dom(SCinstance &inst, CPXENVptr env, CPXLPptr lp);

STATUS cpxsol_build_lp2raw(SCinstance &inst, CPXENVptr env, CPXLPptr lp);
STATUS cpxsol_build_raw2lp(SCinstance &inst, CPXENVptr env, CPXLPptr lp);

#endif // CPX_SOLVER_H
