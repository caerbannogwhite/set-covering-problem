#ifndef BALAS_DENSE_H
#define BALAS_DENSE_H

#include "balas_common.hpp"

int baldns_cut_proc(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::set<int> &wSet);
int baldns_branch_rule1(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::vector<std::set<int>> branchSet, const int maxBranch, const int maxSingl);
int baldns_branch_rule1_test(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::vector<std::unordered_set<int>> branchSet, const int maxBranch, const int maxSingl);

int baldns_over_sat_rows(arma::mat &mat, arma::vec &x);
int baldns_make_prime_cover(const arma::mat &mat, arma::vec &x);
bool baldns_is_cover(const arma::mat &mat, const arma::vec &x);
double baldns_heur_primal_0(arma::mat &mat, arma::vec &obj, arma::vec &x, const int whichFunc);
double baldns_heur_primal_12(arma::mat &mat, arma::vec &obj, arma::vec &x);
double baldns_heur_primal_5b(arma::mat &mat, arma::vec &obj, arma::vec &x, arma::vec &s, arma::vec &u);

bool baldns_is_dual_sol(arma::mat &mat, arma::vec &obj, arma::vec &u);
int baldns_heur_dual_1(arma::mat &mat, arma::vec &x, arma::vec &u, arma::vec &s);
int baldns_heur_dual_3(arma::mat &mat, arma::vec &x, arma::vec &u, arma::vec &s, const double zUpp);

#endif //BALAS_DENSE_H
