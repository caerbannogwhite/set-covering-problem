
#ifndef TESI_SET_COVER_SC_H
#define TESI_SET_COVER_SC_H

#include <cplex.h>
#include "aux.h"

int SCMILPsolver(SCinstance *inst);
int SCsolverbalasrule1(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule1_test(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule1_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule2(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol2(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcoldom(SCinstance *inst, CPXENVptr env, CPXLPptr lp);

#endif //TESI_SET_COVER_SC_H
