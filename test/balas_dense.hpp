
#ifndef SC_BALAS_DENSE_H
#define SC_BALAS_DENSE_H

#include "common.hpp"

int baldns_make_prime_cover(const ublas::matrix<double> &mat, ublas::vector<double> &x);
bool baldns_is_cover(const ublas::matrix<double> &mat, const ublas::vector<double> &x);
double baldns_heur_primal_0(const ublas::matrix<double> &mat, const ublas::vector<double> &obj,
                            ublas::vector<double> &x, std::vector<int> &xSupp,
                            const int whichFunc);

#endif //SC_BALAS_DENSE_H
