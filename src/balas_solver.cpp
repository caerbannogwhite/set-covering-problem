#include "balas_common.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"

STATUS balsol(BALSOLEnv &inst)
{
    double uppVal;
    double lowVal;

    //////////////////////////////      DENSE       //////////////////////////
    balcomm_read_instance_dns(inst);

    arma::vec x = arma::vec(inst.dnsmat.n_cols);
    arma::vec u = arma::vec(inst.dnsmat.n_rows);
    arma::vec s;

    std::set<int> w = std::set<int>();

    std::cout << "BALAS SOLVER - DENSE" << std::endl;

    uppVal = baldns_heur_primal_0(inst.dnsmat, inst.dnsobj, x, 3);
    std::cout << "is prim sol = " << baldns_is_cover(inst.dnsmat, x) << std::endl;
    std::cout << "prim val = " << uppVal << std::endl;

    s = arma::vec(inst.dnsobj);
    lowVal = baldns_heur_dual_1(inst.dnsmat, x, u, s);
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << lowVal << std::endl;

    // generate cut
    baldns_separation_proc(inst.dnsmat, x, s, uppVal, lowVal, w);
    std::cout << "w =\n";
    for (auto it = w.cbegin(); it != w.cend(); ++it)
    {
        std::cout << *it << std::endl;
    }


    //////////////////////////////      SPARSE       //////////////////////////
    balcomm_read_instance_spr(inst);

    arma::sp_mat xSpr = arma::sp_mat(inst.sprmat.n_cols, 1);

    std::cout << "\n\nBALAS SOLVER - SPARSE" << std::endl;

    uppVal = balspr_heur_primal_0(inst.sprmat, inst.dnsobj, xSpr, 3);
    std::cout << "is prim sol = " << balspr_is_cover(inst.sprmat, xSpr) << std::endl;
    std::cout << "prim val = " << uppVal << std::endl;

    s = arma::vec(inst.dnsobj);
    lowVal = balspr_heur_dual_1(inst.sprmat, xSpr, u, s);
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << lowVal << std::endl;
    

    return SC_SUCCESFULL;
}

STATUS balsol_branch_and_bound()
{

    return SC_SUCCESFULL;
}

STATUS balsol_cutting_planes()
{
    return SC_SUCCESFULL;
}