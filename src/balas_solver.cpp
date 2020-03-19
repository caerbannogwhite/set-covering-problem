#include "balas_common.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"

STATUS balsol(BALSOLEnv &inst)
{
    int cutSize;
    double uppVal;
    double lowVal;

    std::unordered_set<int> tSet = std::unordered_set<int>();
    std::unordered_set<int> wSet = std::unordered_set<int>();

    //////////////////////////////      DENSE       //////////////////////////
    balcomm_read_instance_dns(inst);

    arma::vec x = arma::vec(inst.dnsmat.n_cols);
    arma::vec p = arma::vec(inst.dnsmat.n_rows);
    arma::vec u = arma::vec(inst.dnsmat.n_rows);
    arma::vec s;


    std::cout << "BALAS SOLVER - DENSE" << std::endl;

    uppVal = baldns_heur_primal_0(inst.dnsmat, inst.dnsobj, x, 3);
    std::cout << "\n--------        primal 0        --------\n";
    std::cout << "is prim sol = " << baldns_is_cover(inst.dnsmat, x) << std::endl;
    std::cout << "prim val = " << uppVal << std::endl;

    // intialize T(x) set
    p = inst.dnsmat * x;
    for (int i = inst.dnsmat.n_rows; i--;)
    {
        if (fabs(p(i) - 1.0) < SC_EPSILON_SMALL)
        {
            tSet.insert(i);
        }
    }

    s = arma::vec(inst.dnsobj);
    lowVal = baldns_heur_dual_1(inst.dnsmat, x, u, s, tSet);
    std::cout << "\n--------        dual 1        --------\n";
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << lowVal << std::endl;

    lowVal = baldns_heur_dual_3(inst.dnsmat, x, u, s, tSet, uppVal, lowVal);
    std::cout << "\n--------        dual 3        --------\n";
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << lowVal << std::endl;

    // generate cut
    cutSize = baldns_separation_proc(inst.dnsmat, x, s, tSet, wSet, uppVal, lowVal);
    std::cout << "\n--------        separation        --------\n";
    std::cout << "w size = " << cutSize << std::endl;
    //std::cout << "w =\n"; for (auto it = wSet.cbegin(); it != wSet.cend(); ++it) { std::cout << *it << std::endl; }


    //////////////////////////////      SPARSE       //////////////////////////

    tSet.clear();
    wSet.clear();

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