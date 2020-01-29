
#ifndef SC_H
#define SC_H

#include "common.hpp"

int SCMILPsolver(SCinstance *inst);
int SCsolverbalasrule1(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule1_test(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule1_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolverbalasrule2(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol2(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcol_sparse(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCsolvermaxcoldom(SCinstance *inst, CPXENVptr env, CPXLPptr lp);

#endif //SC_H
