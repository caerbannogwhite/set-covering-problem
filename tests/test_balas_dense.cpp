#define CATCH_CONFIG_MAIN // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "balas_dense.hpp"

TEST_CASE("BALAS DENSE - baldns_is_cover", "[BALAS DENSE]")
{
    arma::mat mat = arma::mat(4, 4);
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

    // COVER
    // x = [1, 1, 0, 0]
    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 0.0;
    x(3) = 0.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [1, 0, 1, 0]
    x(0) = 1.0;
    x(1) = 0.0;
    x(2) = 1.0;
    x(3) = 0.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [0, 0, 1, 1]
    x(0) = 0.0;
    x(1) = 0.0;
    x(2) = 1.0;
    x(3) = 1.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [1, 1, 1, 0]
    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 0.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [1, 1, 0, 1]
    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 0.0;
    x(3) = 1.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [0, 1, 1, 1]
    x(0) = 0.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 1.0;
    REQUIRE(baldns_is_cover(mat, x));

    // x = [1, 1, 1, 1]
    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 1.0;
    REQUIRE(baldns_is_cover(mat, x));

    // NOT COVER
    // x = [1, 0, 0, 0]
    x(0) = 1.0;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = 0.0;
    REQUIRE_FALSE(baldns_is_cover(mat, x));

    // x = [0, 1, 0, 0]
    x(0) = 0.0;
    x(1) = 1.0;
    x(2) = 0.0;
    x(3) = 0.0;
    REQUIRE_FALSE(baldns_is_cover(mat, x));

    // x = [0, 1, 1, 0]
    x(0) = 0.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 0.0;
    REQUIRE_FALSE(baldns_is_cover(mat, x));

    // x = [0, 1, 0, 1]
    x(0) = 0.0;
    x(1) = 1.0;
    x(2) = 0.0;
    x(3) = 1.0;
    REQUIRE_FALSE(baldns_is_cover(mat, x));
}

TEST_CASE("BALAS DENSE - baldns_make_prime_cover", "[BALAS DENSE]")
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

    objVal = 5.0;
    baldns_make_prime_cover(mat, obj, x, objVal);

    REQUIRE(baldns_is_cover(mat, x));
    REQUIRE(fabs(objVal - 2.0) < SC_EPSILON_SMALL);
    REQUIRE(fabs(x(2)) < SC_EPSILON_SMALL);
    REQUIRE(fabs(x(3)) < SC_EPSILON_SMALL);
}

TEST_CASE("BALAS DENSE - baldns_heur_primal_0", "[BALAS DENSE]")
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
    x(3) = 0.0;

    objVal = baldns_heur_primal_0(mat, obj, x, 3);
    REQUIRE((objVal - 2.0) < SC_EPSILON_SMALL);

    obj(0) = 3.0;
    obj(1) = 3.0;
    obj(2) = 2.0;
    obj(3) = 1.0;

    x(0) = 0.0;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = 0.0;

    objVal = baldns_heur_primal_0(mat, obj, x, 3);
    REQUIRE((objVal - 3.0) < SC_EPSILON_SMALL);

    obj(0) = 3.0;
    obj(1) = 3.0;
    obj(2) = 2.0;
    obj(3) = 1.0;

    x(0) = 1.0;
    x(1) = 0.0;
    x(2) = 0.0;
    x(3) = 0.0;

    objVal = baldns_heur_primal_0(mat, obj, x, 3);
    REQUIRE((objVal - 5.0) < SC_EPSILON_SMALL);

    obj(0) = 3.0;
    obj(1) = 3.0;
    obj(2) = 2.0;
    obj(3) = 1.0;

    x(0) = 1.0;
    x(1) = 1.0;
    x(2) = 1.0;
    x(3) = 1.0;

    objVal = baldns_heur_primal_0(mat, obj, x, 3);
    REQUIRE((objVal - 6.0) < SC_EPSILON_SMALL);
}