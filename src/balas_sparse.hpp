#ifndef BALAS_SPARSE_H
#define BALAS_SPARSE_H

#include "balas_common.hpp"

int balspr_cut_proc(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::set<int> &wSet);
int balspr_branch_rule1(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::vector<std::set<int>> branchSet, const int maxBranch, const int maxSingl);
int balspr_branch_rule1_test(arma::mat &mat, arma::vec &x, arma::vec &s, double zUpp, double zLow, std::vector<std::unordered_set<int>> branchSet, const int maxBranch, const int maxSingl);

int balspr_over_sat_rows(arma::mat &mat, arma::vec &x);
int balspr_make_prime_cover(const arma::mat &mat, arma::vec &x);
bool balspr_is_cover(const arma::mat &mat, arma::vec &x);

double balspr_heur_primal_0(arma::mat &mat, arma::vec &obj, arma::vec &x, const int whichFunc);
double balspr_heur_primal_12(arma::mat &mat, arma::vec &obj, arma::vec &x);
double balspr_heur_primal_5b(arma::mat &mat, arma::vec &obj, arma::vec &x, arma::vec &s, arma::vec &u);

bool balspr_is_dual_sol(arma::mat &mat, arma::vec &obj, arma::vec &u);
int balspr_heur_dual_1(arma::mat &mat, arma::vec &x, arma::vec &u, arma::vec &s);
int balspr_heur_dual_3(arma::mat &mat, arma::vec &x, arma::vec &u, arma::vec &s, const double zUpp);

#endif // BALAS_SPARSE_H
