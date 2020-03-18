#include "balas_common.hpp"
#include "balas_dense.hpp"
#include "balas_sparse.hpp"

STATUS balsol(BALSOLEnv &inst)
{
    double val;



    //////////////////////////////      DENSE       //////////////////////////
    balcomm_read_instance_dns(inst);

    arma::vec x = arma::vec(inst.dnsmat.n_cols);
    arma::vec u = arma::vec(inst.dnsmat.n_rows);
    arma::vec s;

    std::cout << "BALAS SOLVER - DENSE" << std::endl;
    //std::cout << "mat =\n"
    //          << inst.dnsmat << std::endl;
    //std::cout << "obj =\n"
    //          << inst.dnsobj << std::endl;

    val = baldns_heur_primal_0(inst.dnsmat, inst.dnsobj, x, 3);
    std::cout << "is prim sol = " << baldns_is_cover(inst.dnsmat, x) << std::endl;
    std::cout << "prim val = " << val << std::endl;
    //std::cout << "x =\n"
    //          << x << std::endl;

    s = arma::vec(inst.dnsobj);
    val = baldns_heur_dual_1(inst.dnsmat, x, u, s);
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << val << std::endl;
    //std::cout << "u =\n"
    //          << u << std::endl;
    //std::cout << "s =\n"
    //          << s << std::endl;



    //////////////////////////////      SPARSE       //////////////////////////
    balcomm_read_instance_spr(inst);

    arma::sp_mat xSpr = arma::sp_mat(inst.sprmat.n_cols, 1);
    //arma::vec u = arma::vec(inst.dnsmat.n_rows);
    //arma::vec s = arma::vec(inst.dnsmat.n_cols);

    std::cout << "\n\nBALAS SOLVER - SPARSE" << std::endl;
    //std::cout << "mat =\n"
    //          << inst.dnsmat << std::endl;
    //std::cout << "obj =\n"
    //          << inst.dnsobj << std::endl;

    val = balspr_heur_primal_0(inst.sprmat, inst.dnsobj, xSpr, 3);
    std::cout << "is prim sol = " << balspr_is_cover(inst.sprmat, xSpr) << std::endl;
    std::cout << "prim val = " << val << std::endl;
    //std::cout << "x =\n"
    //          << x << std::endl;

    s = arma::vec(inst.dnsobj);
    val = balspr_heur_dual_1(inst.sprmat, xSpr, u, s);
    std::cout << "is dual sol = " << baldns_is_dual_sol(inst.dnsmat, inst.dnsobj, u) << std::endl;
    std::cout << "dual val = " << val << std::endl;
    //std::cout << "u =\n"
    //          << u << std::endl;
    //std::cout << "s =\n"
    //          << s << std::endl;

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