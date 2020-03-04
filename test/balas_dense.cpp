
#include "balas_dense.hpp"

/**
 * Given a SCP matrix and a SCP solution x,
 * return the number of removed columns
 * from x (to make x a prime cover).
 * 
 * @param mat - ublas::matrix<double> a SCP matrix
 * @param x - ublas::vector<double> a SCP solution
 * @return the number of removed columns from x
 */
int baldns_make_prime_cover(const ublas::matrix<double> &mat, ublas::vector<double> &x)
{
    int cntRemoved;
    size_t i, j;

    std::unique_ptr<ublas::vector<double>> matDotX (new ublas::vector<double>(mat.size1()));
    *matDotX = ublas::prod(mat, x);

    cntRemoved = 0;
    for (j = mat.size2(); j--;)
    {
        if (x(j) > SC_EPSILON_SMALL)
        {
            x(j) = 0.0;
            for (i = 0; i < mat.size1(); ++i)
            {
                if ((mat(i, j) > SC_EPSILON_SMALL) && ((*matDotX)(i) < 2.0))
                {
                    x(j) = 1.0;
                    break;
                }
            }

            if (x(j) < SC_EPSILON_SMALL)
            {
                cntRemoved++;
                (*matDotX) -= ublas::column(mat, j);
            }
        }
    }

    return cntRemoved;
}

/**
 * Given a SCP dense matrix mat and a solution vector x,
 * return true if x is a cover for mat, false otherwise.
 * 
 * @param mat - ublas::matrix<double> SCP matrix
 * @param x - ublas::vector<double> a SCP solution vector
 * @return true if x covers mat, false otherwise
 */
bool baldns_is_cover(const ublas::matrix<double> &mat, const ublas::vector<double> &x)
{
    bool flag;
    std::unique_ptr<ublas::vector<double>> res(new ublas::vector<double>(mat.size1()));
    *res = ublas::prod(mat, x);

    flag = true;
    for (auto it = (*res).cbegin(); it != (*res).cend(); ++it)
    {
        if (*it < 1.0)
        {
            flag = false;
        }
    }

    return flag;
}

/**
 * Primal heuristic for the SCP as described in "Set Covering algorithms using cutting
 * planes, heuristics, and subgradient optimization: a computational study"
 * by Egon Balas and Andrew Ho - Carnegie-Mellon University, Pittsburgh, PA, U.S.A.
 * 
 * @param mat - ublas::matrix<double> SCP matrix
 * @param obj - ublas::vector<double> SCP objective values
 * @param x - ublas::vector<double> a solution for the SCP
 * @param xSupp - 
 * @param whichFunc - select function for heuristic
 * @return the value of the best primal solution found
 */
double baldns_heur_primal_0(const ublas::matrix<double> &mat, const ublas::vector<double> &obj,
                            ublas::vector<double> &x, std::vector<int> &xSupp,
                            const int whichFunc)
{
    int col, row, rcnt;
    size_t i, j, coveredRows;
    double val, zUpp, cnt;
    double (*func)(const double, const double);

    std::unique_ptr<std::vector<std::pair<int, int>>> rSet (new std::vector<std::pair<int, int>>(mat.size1()));

    // all the functions defined by Balas and Ho
    // plus the last two defined by Vasko and Wilson
    switch (whichFunc)
    {
    case 1:
        func = [](const double c, const double k) -> double { return c; };
        break;
    case 2:
        func = [](const double c, const double k) -> double { return c / k; };
        break;
    case 3:
        func = [](const double c, const double k) -> double { return k < 2 ? c : c / log2(k); };
        break;
    case 4:
        func = [](const double c, const double k) -> double { return k < 2 ? c / k : c / k * log2(k); };
        break;
    case 5:
        func = [](const double c, const double k) -> double { return k < 3 ? c / k : c / k * log(k); };
        break;
    case 6:
        func = [](const double c, const double k) -> double { return c / (k * k); };
        break;
    case 7:
        func = [](const double c, const double k) -> double { return sqrt(c) / (k * k); };
        break;
    default:
        func = [](const double c, const double k) -> double { return k < 2 ? c : c / log2(k); };
        break;
    }

    // get the starting value of the current solution x
    zUpp = ublas::inner_prod(obj, x);

    // find the number of already covered rowa
    coveredRows = 0;
    for (i = mat.size1(); i--;)
    {
        cnt = ublas::inner_prod(ublas::row(mat, i), x);
        rcnt = round(ublas::sum(ublas::row(mat, i)));

        if (fabs(cnt) < SC_EPSILON_SMALL)
        {
            (*rSet).push_back(std::make_pair(rcnt, i));
        }
        else
        {
            coveredRows++;
        }
    }

    // sort not covered rows by decreasing number of
    // non-zero elements in the row (read rSet from back to front)
    std::sort((*rSet).begin(), (*rSet).end(), [](std::pair<int, int> a, std::pair<int, int> b) -> bool { return a.first > b.first; });

    while (coveredRows < mat.size1())
    {
        // get last element of R
        row = (*rSet).back().second;
        (*rSet).pop_back();

        val = INFINITY;
        col = -1;
        for (j = 0; j < mat.size2(); ++j)
        {
            if (fabs(mat(row, j)) < SC_EPSILON_SMALL)
            {
                continue;
            }

            cnt = ublas::sum(ublas::column(mat, j));

            if (func(obj(j), cnt) < val)
            {
                val = func(obj(j), cnt);
                col = j;
            }
        }

        xSupp.push_back(col);
        x(col) = 1.0;
        zUpp += obj(col);

        // remove covered rows
        for (auto it = (*rSet).begin(); it != (*rSet).end();)
        {
            if (fabs(mat((*it).second, col) - 1.0) < SC_EPSILON_SMALL)
            {
                coveredRows++;
                (*rSet).erase(it);
            } else
            {
                ++it;
            }
        }
    }

    // make the cover found a prime cover
    baldns_make_prime_cover(mat, x);

    for (auto it = xSupp.begin(); it != xSupp.end();)
    {
        if (x[*it] < SC_EPSILON_SMALL)
        {
            zUpp -= obj(*it);
            xSupp.erase(it);
        } else
        {
            ++it;
        }
    }

    delete &rSet;

    return zUpp;
}
