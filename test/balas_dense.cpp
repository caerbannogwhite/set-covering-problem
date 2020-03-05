
#include "balas_dense.hpp"

/**
 * Given a SCP matrix and a SCP solution x,
 * return the number of removed columns
 * from x (to make x a prime cover).
 * 
 * @param mat - arma::mat a SCP matrix
 * @param x - arma::vec a SCP solution
 * @return the number of removed columns from x
 */
int baldns_make_prime_cover(const arma::mat &mat, arma::vec &x)
{
    int cntRemoved;
    size_t i, j;

    std::unique_ptr<arma::vec> matDotX (new arma::vec(mat.n_rows));
    *matDotX = mat * x;

    cntRemoved = 0;
    for (j = mat.n_cols; j--;)
    {
        if (x(j) > SC_EPSILON_SMALL)
        {
            x(j) = 0.0;
            for (i = 0; i < mat.n_rows; ++i)
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
                (*matDotX) -= mat.col(j);
            }
        }
    }

    return cntRemoved;
}

/**
 * Given a SCP dense matrix mat and a solution vector x,
 * return true if x is a cover for mat, false otherwise.
 * 
 * @param mat - arma::mat SCP matrix
 * @param x - arma::vec a SCP solution vector
 * @return true if x covers mat, false otherwise
 */
bool baldns_is_cover(const arma::mat &mat, arma::vec &x)
{
    bool flag;

    std::unique_ptr<arma::vec> matDotX(new arma::vec(mat.n_rows));
    *matDotX = mat * x;

    flag = true;
    for (auto it = (*matDotX).cbegin(); it != (*matDotX).cend(); ++it)
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
 * @param mat - arma::mat SCP matrix
 * @param obj - arma::vec SCP objective values
 * @param x - arma::vec a solution for the SCP
 * @param xSupp - 
 * @param whichFunc - select function for heuristic
 * @return the value of the best primal solution found
 */
double baldns_heur_primal_0(arma::mat &mat, arma::vec &obj,
                            arma::vec &x, std::vector<int> &xSupp,
                            const int whichFunc)
{
    int col, row, rcnt;
    size_t i, j, coveredRows;
    double val, zUpp, cnt;
    double (*func)(const double, const double);

    std::unique_ptr<std::vector<std::pair<int, int>>> rSet (new std::vector<std::pair<int, int>>);

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
    zUpp = arma::dot(obj, x);

    // find the number of already covered rows
    coveredRows = 0;
    for (i = mat.n_rows; i--;)
    {
        cnt = arma::dot(mat.row(i), x);
        rcnt = round(arma::sum(mat.row(i)));

        if (fabs(cnt) < SC_EPSILON_SMALL)
        {
            (*rSet).push_back(std::make_pair(rcnt, i));
        }
        else
        {
            coveredRows++;
        }
    }

    /*std::cout << "cov rows = " << coveredRows << std::endl;
    std::cout << "rows = " << std::endl;
    for (auto it : *rSet)
    {
        std::cout << it.first << " " << it.second << std::endl; 
    }*/

    // sort not covered rows by decreasing number of
    // non-zero elements in the row (read rSet from back to front)
    std::sort((*rSet).begin(), (*rSet).end(), [](std::pair<int, int> a, std::pair<int, int> b) -> bool { return a.first < b.first; });

    std::cout << "\n\nBALAS HEUR 0\n";
    std::cout << "sort rows = " << std::endl;
    for (auto it : *rSet)
    {
        std::cout << it.first << " " << it.second << std::endl;
    }

    while (coveredRows < mat.n_rows)
    {
        // get last element of R
        row = (*rSet).back().second;
        (*rSet).pop_back();

        val = INFINITY;
        col = -1;
        for (j = 0; j < mat.n_cols; ++j)
        {
            // if mat[row, col] is 0 or if x(col) is already taken, skip this column
            if ((fabs(mat(row, j)) < SC_EPSILON_SMALL) || (x(j) > SC_EPSILON_SMALL))
            {
                continue;
            }

            cnt = arma::sum(mat.col(j));

            // use the column j that minimise the value
            // of the selected function
            if (func(obj(j), cnt) < val)
            {
                val = func(obj(j), cnt);
                col = j;
            }
        }

        std::cout << "row = " << row << " col = " << col << " cov rows = " << coveredRows << std::endl;

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

    return zUpp;
}

/*int SCbalasheurdual1(arma::mat &mat, arma::vec &x, arma::vec &u, arma::vec &s)
{

    char *matr;
    int cnt, rcnt, firsttime, itercnt, i, j, row, srtrowssize = 2, srtrowslen = 0;
    SCi2tuple *srtrows;

    for (i = 0; i < nrows; ++i)
    {
        u[i] = 0.0;
    }

    srtrows = (SCi2tuple *)malloc(srtrowssize * sizeof(SCi2tuple));
    for (i = 0; i < nrows; ++i)
    {

        // No need to generate R and T(x) sets: cnt says how
        // much a row is covered and if cnt < 1 (minimum
        // cover) skip the row
        cnt = 0;
        rcnt = 0;
        matr = &mat[i * ncols];
        for (j = 0; j < ncols; ++j)
        {
            cnt += (matr[j] & x[j]);
            rcnt += matr[j];
        }

        if (cnt == 1)
        {
            srtrows[srtrowslen].a = rcnt;
            srtrows[srtrowslen++].b = i;

            if (srtrowslen == srtrowssize)
            {
                srtrowssize <<= 1;
                srtrows = (SCi2tuple *)realloc(srtrows, srtrowssize * sizeof(SCi2tuple));
            }
        }
    }

    qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);

    itercnt = 0;
    firsttime = 1;
    while (itercnt < srtrowslen)
    {

        row = srtrows[itercnt++].b;

        u[row] = INFINITY;
        matr = &mat[row * ncols];
        for (j = 0; j < ncols; ++j)
        {
            u[row] = matr[j] && (s[j] < u[row]) ? s[j] : u[row];
        }

        for (j = 0; j < ncols; ++j)
        {
            s[j] -= u[row] * matr[j];
        }

        if ((itercnt == srtrowslen) && firsttime)
        {

            itercnt = 0;
            firsttime = 0;
            srtrowslen = 0;

            for (i = 0; i < nrows; ++i)
            {

                // No need to generate R and T(x) sets: cnt says how
                // much a row is covered and if cnt < 1 (minimum
                // cover) skip the row
                cnt = 0;
                rcnt = 0;
                matr = &mat[i * ncols];
                for (j = 0; j < ncols; ++j)
                {
                    cnt += (matr[j] & x[j]);
                    rcnt += matr[j];
                }

                if (cnt > 1)
                {
                    srtrows[srtrowslen].a = rcnt;
                    srtrows[srtrowslen++].b = i;

                    if (srtrowslen == srtrowssize)
                    {
                        srtrowssize <<= 1;
                        srtrows = (SCi2tuple *)realloc(srtrows, srtrowssize * sizeof(SCi2tuple));
                    }
                }
            }
            qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);
        }
    }

    return 0;
}

int baldns_heur_dual_3(arma::mat &mat,arma::vec &x, arma::vec &u, arma::vec &s, const double zUpp)
{
    char *matr;
    int cnt, itercnt, i, j, row, srtrowssize = 2, srtrowslen = 0;
    double sdotx, usum, val;
    SCi2tuple *srtrows;

    sdotx = 0.0;
    for (j = 0; j < ncols; ++j)
    {
        sdotx += s[j] * x[j];
    }

    usum = 0.0;
    for (i = 0; i < nrows; ++i)
    {
        usum += u[i];
    }

    if (sdotx >= (zu - usum))
    {
        goto TERMINATE1;
    }

    srtrows = (SCi2tuple *)malloc(srtrowssize * sizeof(SCi2tuple));
    for (i = 0; i < nrows; ++i)
    {

        // No need to generate R and T(x) sets: cnt says how
        // much a row is covered and if cnt < 1 (minimum
        // cover) skip the row
        cnt = 0;
        matr = &mat[i * ncols];
        for (j = 0; j < ncols; ++j)
        {
            cnt += (matr[j] & x[j]);
        }

        if ((u[i] > SC_EPSILON_SMALL) && (cnt > 1))
        {
            srtrows[srtrowslen].a = cnt;
            srtrows[srtrowslen++].b = i;

            if (srtrowslen == srtrowssize)
            {
                srtrowssize <<= (char)1;
                srtrows = (SCi2tuple *)realloc(srtrows, srtrowssize * sizeof(SCi2tuple));
            }
        }
    }

    if (srtrowslen == 0)
    {
        free(srtrows);
        return 0;
    }

    qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);

    itercnt = 0;
    while (sdotx < (zu - usum - SC_EPSILON_SMALL))
    {

        if (itercnt == srtrowslen)
        {
            free(srtrows);
            return -1;
        }

        row = srtrows[itercnt++].b;

        val = u[row];
        matr = &mat[row * ncols];
        for (j = 0; j < ncols; ++j)
        {
            s[j] += val * matr[j];
            sdotx += val * (matr[j] & x[j]);
        }

        usum -= val;
        u[row] = 0.0;
    }

    return 0;
}*/
