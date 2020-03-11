#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "balas_sparse.hpp"

TEST_CASE("BALAS SPARSE - balspr_is_cover", "[BALAS SPARSE]")
{
    arma::umat locs;
    arma::vec vals;
    arma::sp_mat mat;
    arma::sp_mat x;

    // matrix
    //      1 1 0 1
    //      0 1 1 0
    //      1 0 1 0
    //      1 0 0 1
    locs << 0 << 0 << arma::endr
         << 0 << 1 << arma::endr
         << 0 << 3 << arma::endr
         << 1 << 1 << arma::endr
         << 1 << 2 << arma::endr
         << 2 << 0 << arma::endr
         << 2 << 2 << arma::endr
         << 3 << 0 << arma::endr
         << 3 << 3 << arma::endr;

    vals << 1.0 << 1.0 << 1.0 << 1.0 << 1.0 << 1.0 << 1.0 << 1.0 << 1.0;

    mat = arma::sp_mat(locs, vals);

    // COVER
    // x = [1, 1, 0, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr
         << 0 << 1 << arma::endr;
    vals << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [1, 0, 1, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr
         << 0 << 2 << arma::endr;
    vals << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [0, 0, 1, 1]
    locs.clear();
    vals.clear();
    locs << 0 << 2 << arma::endr
         << 0 << 3 << arma::endr;
    vals << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [1, 1, 1, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr
         << 0 << 1 << arma::endr
         << 0 << 2 << arma::endr;
    vals << 1.0 << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [1, 1, 0, 1]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr
         << 0 << 1 << arma::endr
         << 0 << 3 << arma::endr;
    vals << 1.0 << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [0, 1, 1, 1]
    locs.clear();
    vals.clear();
    locs << 0 << 1 << arma::endr
         << 0 << 2 << arma::endr
         << 0 << 3 << arma::endr;
    vals << 1.0 << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // x = [1, 1, 1, 1]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr
         << 0 << 1 << arma::endr
         << 0 << 2 << arma::endr
         << 0 << 3 << arma::endr;
    vals << 1.0 << 1.0 << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(balspr_is_cover(mat, x));

    // NOT COVER
    // x = [1, 0, 0, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 0 << arma::endr;
    vals << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(!balspr_is_cover(mat, x));

    // x = [0, 1, 0, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 1 << arma::endr;
    vals << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(!balspr_is_cover(mat, x));

    // x = [0, 1, 1, 0]
    locs.clear();
    vals.clear();
    locs << 0 << 1 << arma::endr
         << 0 << 2 << arma::endr;
    vals << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(!balspr_is_cover(mat, x));

    // x = [0, 1, 0, 1]
    locs.clear();
    vals.clear();
    locs << 0 << 1 << arma::endr
         << 0 << 3 << arma::endr;
    vals << 1.0 << 1.0 << arma::endr;
    x = arma::sp_mat(locs, vals);
    REQUIRE(!balspr_is_cover(mat, x));
}

/*TEST_CASE("BALAS SPARSE - balspr_make_prime_cover", "[BALAS SPARSE]")
{
    double objVal;
    arma::mat mat = arma::mat(4, 4);
    arma::vec obj = arma::vec(4);
    arma::vec x = arma::vec(4);

    // matrix
    //      1 1 0 1
    //      0 1 1 0
    //      1 0 1 0
    //      1 0 0 1
    mat(0, 0) = 1.0;
    mat(0, 1) = 1.0;
    mat(0, 3) = 1.0;
    mat(1, 1) = 1.0;
    mat(1, 2) = 1.0;
    mat(2, 0) = 1.0;
    mat(2, 2) = 1.0;
    mat(3, 0) = 1.0;
    mat(3, 3) = 1.0;

    obj(0) = 1.0;
    obj(1) = 1.0;
    obj(2) = 1.0;
    obj(3) = 2.0;

    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 1.0;

    balspr_make_prime_cover(mat, x);
    objVal = arma::dot(obj, x);

    REQUIRE(balspr_is_cover(mat, x));
    REQUIRE(fabs(objVal - 2.0) < SC_EPSILON_SMALL);
    REQUIRE(fabs(x(2)) < SC_EPSILON_SMALL);
    REQUIRE(fabs(x(3)) < SC_EPSILON_SMALL);
}

TEST_CASE("BALAS SPARSE - balspr_heur_primal_0", "[BALAS SPARSE]")
{
    double objVal;
    arma::mat mat = arma::mat(4, 4);
    arma::vec obj = arma::vec(4);
    arma::vec x = arma::vec(4);

    // matrix
    //      1 1 0 1
    //      0 1 1 0
    //      1 0 1 0
    //      1 0 0 1
    mat(0, 0) = 1.0;
    mat(0, 1) = 1.0;
    mat(0, 3) = 1.0;
    mat(1, 1) = 1.0;
    mat(1, 2) = 1.0;
    mat(2, 0) = 1.0;
    mat(2, 2) = 1.0;
    mat(3, 0) = 1.0;
    mat(3, 3) = 1.0;

    obj(0) = 1.0;
    obj(1) = 1.0;
    obj(2) = 1.0;
    obj(3) = 2.0;

    x(0) = 0.0;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = .0;

    balspr_heur_primal_0(mat, obj, x, 3);

    std::cout << x << std::endl;
}*/