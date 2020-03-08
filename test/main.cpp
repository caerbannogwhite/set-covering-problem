
#include "common.hpp"
#include "balas_dense.hpp"

int main()
{
    std::vector<int> xSupp;

    arma::mat mat = arma::mat(4, 4);
    arma::vec obj = arma::vec(4);
    arma::vec x = arma::vec(4);

    mat(0, 0) = 1.0;
    mat(0, 1) = 1.0;
    mat(0, 3) = 1.0;
    mat(1, 1) = 1.0;
    mat(1, 2) = 1.0;
    mat(2, 0) = 1.0;
    mat(2, 2) = 1.0;
    mat(2, 3) = 1.0;
    mat(3, 0) = 1.0;
    mat(3, 3) = 1.0;

    obj(0) = 1.0;
    obj(1) = 1.0;
    obj(2) = 1.0;
    obj(3) = 2.0;

    x(0) = 1.0;
    //x(1) = 1.0;
    x(2) = 1.0;
    //x(3) = 1.0;

    std::cout << "mat = " << std::endl;
    std::cout << mat << std::endl;

    std::cout << "obj = " << std::endl;
    std::cout << obj << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;

    arma::uvec findMat = arma::find(mat > 0.5);
    arma::uvec findX = arma::find(x > 0.5);

    std::cout << "find mat =\n";
    std::cout << findMat << std::endl;

    for (auto it = findX.cbegin(); it != findX.cend(); ++it)
    {
        std::cout << *it << std::endl;
    }

    /*std::cout << "is x cover = " << baldns_is_cover(mat, x) << std::endl;
    std::cout << "removed cols = " << baldns_make_prime_cover(mat, x) << std::endl;

    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;

    x = arma::vec(4);

    std::cout << "score = " << baldns_heur_primal_0(mat, obj, x, xSupp, 3) << std::endl;

    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;

    std::cout << "mat . x = " << std::endl;
    std::cout << mat * x << std::endl;*/

    return 0;
}