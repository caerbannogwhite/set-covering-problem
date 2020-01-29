
#ifndef SET_COVERING_PREPROCESSING_H
#define SET_COVERING_PREPROCESSING_H

#include "common.hpp"

//int util_simple_presolver(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
int SCdominancepresolver(SCinstance *inst, CPXENVptr env, CPXLPptr lp);
//int util_read_input(SCinstance *inst);
//int util_print_solution(SCinstance *inst);
//int util_print_table(SCinstance *inst);

#endif //SET_COVERING_PREPROCESSING_H
