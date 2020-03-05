
#ifndef SC_BALAS_DENSE_H
#define SC_BALAS_DENSE_H

#include "common.hpp"

int baldns_make_prime_cover(const arma::mat &mat, arma::vec &x);
bool baldns_is_cover(const arma::mat &mat, arma::vec &x);
double baldns_heur_primal_0(arma::mat &mat, arma::vec &obj,
                            arma::vec &x, std::vector<int> &xSupp,
                            const int whichFunc);

#endif //SC_BALAS_DENSE_H
