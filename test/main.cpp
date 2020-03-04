
#include "common.hpp"
#include "balas_dense.hpp"

int main()
{
    ublas::matrix<double> mat = ublas::matrix<double>(4, 4);
    ublas::vector<double> obj = ublas::vector<double>(4);
    ublas::vector<double> x = ublas::vector<double>(4);

    mat(0, 0) = 1;
    mat(0, 1) = 1;
    mat(0, 3) = 1;
    mat(1, 0) = 1;
    mat(1, 2) = 1;
    mat(2, 3) = 1;
    mat(3, 2) = 1;
    mat(3, 3) = 1;

    obj(0) = 1;
    obj(1) = 1;
    obj(2) = 1;
    obj(3) = 2;

    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
    x(3) = 1;

    std::cout << "x = " << x << std::endl;
    std::cout << baldns_is_cover(mat, x) << std::endl;
    baldns_make_prime_cover(mat, x);
    std::cout << "x = " << x << std::endl;
    std::cout << baldns_is_cover(mat, x) << std::endl;

    return 0;
}