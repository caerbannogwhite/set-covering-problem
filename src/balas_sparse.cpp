#include "balas_sparse.hpp"

/**
 * Given a SCP sparse matrix and a SCP solution x,
 * return the number of removed columns
 * from x (to make x a prime cover).
 * 
 * @param mat - arma::sp_mat a SCP sparse matrix
 * @param x - arma::sp_mat a SCP sparse solution COLUMN vector
 * @return the number of removed columns from x
 */
int balspr_make_prime_cover(arma::sp_mat &mat, arma::vec &obj, arma::sp_mat &x, double &zUpp)
{
	bool remove;
	int col;
	int cntRemoved;

	std::unique_ptr<arma::vec> matDotXPtr(new arma::vec(mat.n_rows));
	std::unique_ptr<std::unordered_set<int>> removedColsPtr(new std::unordered_set<int>());
	(*matDotXPtr) = mat * x;

	cntRemoved = 0;
	for (auto jt = --x.end(); ; --jt) // Reverse iteration
	{
		col = jt.row();
		remove = true;
		for (auto it = mat.begin_col(col); it != mat.end_col(col); ++it)
		{
			if ((mat(it.row(), col) > SC_EPSILON_SMALL) && ((*matDotXPtr)(it.row()) < (2.0 - SC_EPSILON_SMALL)))
			{
				remove = false;
				break;
			}
		}

		if (remove)
		{
			zUpp -= obj(col);
			removedColsPtr->insert(col);
			++cntRemoved;
			(*matDotXPtr) -= mat.col(col);
		}


		if (jt == x.begin()) // Check stop condition
		{
			break;
		}
	}

	for (auto jt = removedColsPtr->cbegin(); jt != removedColsPtr->cend(); ++jt)
	{
		x(*jt) = 0.0;
	}
	x.remove_zeros();

	return cntRemoved;
}

/**
 * Given a SCP sparse matrix mat and a solution vector x,
 * return true if x is a cover for mat, false otherwise.
 * 
 * @param mat - arma::sp_mat SCP sparse matrix
 * @param x - arma::sp_mat a SCP sparse solution COLUMN vector
 * @return true if x covers mat, false otherwise
 */
bool balspr_is_cover(arma::sp_mat &mat, arma::sp_mat &x)
{
	std::unique_ptr<arma::vec> matDotXPtr(new arma::vec(mat.n_rows));
	(*matDotXPtr) = mat * x;
	return arma::all((*matDotXPtr) > (1.0 - SC_EPSILON_SMALL));
}

/**
 * Primal heuristic for the SCP as described in "Set Covering algorithms using cutting
 * planes, heuristics, and subgradient optimization: a computational study"
 * by Egon Balas and Andrew Ho - Carnegie-Mellon University, Pittsburgh, PA, U.S.A.
 * 
 * Functions description
 * 
 * 1 - f(c, k) = c
 * 
 * 2 - f(c, k) = c / k
 * 
 * 3 - f(c, k) = c 				 	if k < 2 
 *				 c / log2(k)	 	otherwise

 * 4 - f(c, k) = c / k 			 	if k < 2
 * 				 c / k * log2(k) 	otherwise
 * 
 * 5 - f(c, k) = c / k				if k < 3
 * 				 c / k * ln(k)		otherwise
 * 
 * 6 - f(c, k) = c / k^2
 * 
 * 7 - f(c, k) = sqrt(c) / k^2
 * 
 * @param mat - arma::mat_sp SCP sparse matrix
 * @param obj - arma::vac SCP dense objective values
 * @param x - arma::sp_mat a SCP sparse solution COLUMN vector
 * @param whichFunc - select function for heuristic (from 1 to 7)
 * @return the value of the best primal solution found
 */
double balspr_heur_primal_0(arma::sp_mat &mat, arma::vec &obj, arma::sp_mat &x,
							const int whichFunc)
{
	size_t i;
	size_t j;
	int rowIdx;
	int colIdx;
	int cnt;
	double val;
	double bestVal;
	double zUpp;
	double (*func)(const double, const double);

	std::unique_ptr<std::vector<std::pair<int, int>>> notCoveredRowsPtr(new std::vector<std::pair<int, int>>);

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

	// find rows not already covered and put them in set R
	for (i = mat.n_rows; i--;)
	{
		val = arma::dot(mat.row(i), x.t());
		cnt = mat.row(i).n_nonzero;

		if (fabs(val) < SC_EPSILON_SMALL)
		{
			notCoveredRowsPtr->push_back(std::make_pair(cnt, i));
		}
	}

	// sort not covered rows by decreasing number of
	// non-zero elements in the row (read notCoveredRowsPtr from back to front)
	std::sort(notCoveredRowsPtr->begin(), notCoveredRowsPtr->end(), [](std::pair<int, int> a, std::pair<int, int> b) -> bool { return a.first < b.first; });

	while (notCoveredRowsPtr->size() > 0)
	{
		// get last element of set R
		rowIdx = notCoveredRowsPtr->back().second;
		notCoveredRowsPtr->pop_back();

		bestVal = INFINITY;
		colIdx = -1;
		for (j = 0; j < mat.n_cols; ++j)
		{
			// if mat[row, col] is 0 or if x(col) is already taken, skip this column
			if ((fabs(mat(rowIdx, j)) < SC_EPSILON_SMALL) || (x(j) > SC_EPSILON_SMALL))
			{
				continue;
			}

			val = mat.col(j).n_elem;

			// use the column j that minimise the value
			// of the selected function
			if (func(obj(j), val) < bestVal)
			{
				bestVal = func(obj(j), val);
				colIdx = j;
			}
		}

		// update solution and objective value
		x(colIdx) = 1.0;
		zUpp += obj(colIdx);

		// remove covered rows from set R
		for (auto it = notCoveredRowsPtr->begin(); it != notCoveredRowsPtr->end();)
		{
			if (fabs(mat((*it).second, colIdx) - 1.0) < SC_EPSILON_SMALL)
			{
				notCoveredRowsPtr->erase(it);
			}
			else
			{
				++it;
			}
		}
	}

	// Make the cover a prime cover
	balspr_make_prime_cover(mat, obj, x, zUpp);

	return zUpp;
}

/*int SCbalasheurdual1_sparse(const int *rmatbeg, const int *rmatind, const int nrows, const int *xind,
							const int xindlen, double *u, double *s)
{

	int cnt, firsttime, itercnt, i, j, k, row, srtrowslen = 0;
	SCi2tuple *srtrows;

	for (i = 0; i < nrows; ++i)
	{
		u[i] = 0.0;
	}

	srtrows = (SCi2tuple *)malloc(nrows * sizeof(SCi2tuple));
	for (i = 0; i < nrows; ++i)
	{

		// No need to generate R and T(x) sets: cnt says how
		// much a row is covered and if cnt < 1 (minimum
		// cover) skip the row
		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i + 1] && k < xindlen)
		{
			if (rmatind[j] == xind[k])
			{
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k])
			{
				j++;
				continue;
			}

			if (rmatind[j] > xind[k])
			{
				k++;
				continue;
			}
		}

		if (cnt == 1)
		{
			srtrows[srtrowslen].a = rmatbeg[i + 1] - rmatbeg[i];
			srtrows[srtrowslen++].b = i;
		}
	}

	qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);

	itercnt = 0;
	firsttime = 1;

	while (itercnt < srtrowslen)
	{

		row = srtrows[itercnt++].b;

		u[row] = INFINITY;
		for (j = rmatbeg[row]; j < rmatbeg[row + 1]; ++j)
		{
			u[row] = (s[rmatind[j]] < u[row]) ? s[rmatind[j]] : u[row];
		}

		for (j = rmatbeg[row]; j < rmatbeg[row + 1]; ++j)
		{
			s[rmatind[j]] -= u[row];
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
				j = rmatbeg[i];
				k = 0;
				while (j < rmatbeg[i + 1] && k < xindlen)
				{
					if (rmatind[j] == xind[k])
					{
						cnt++;
						j++;
						k++;
						continue;
					}

					if (rmatind[j] < xind[k])
					{
						j++;
						continue;
					}

					if (rmatind[j] > xind[k])
					{
						k++;
						continue;
					}
				}

				if (cnt > 1)
				{
					srtrows[srtrowslen].a = rmatbeg[i + 1] - rmatbeg[i];
					srtrows[srtrowslen++].b = i;
				}
			}
			qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);
		}
	}

	free(srtrows);
	return 0;
}

int SCbalasheurdual3_sparse(const int *rmatbeg, const int *rmatind, const int nrows,
							const int *xind, int xindlen, double *u, double *s, const double zu)
{

	int cnt, itercnt, i, j, k, row, srtrowslen = 0;
	double sdotx, usum, val;
	SCi2tuple *srtrows;

	sdotx = 0.0;
	for (j = 0; j < xindlen; ++j)
	{
		sdotx += s[xind[j]];
	}

	usum = 0.0;
	for (i = 0; i < nrows; ++i)
	{
		usum += u[i];
	}

	if (sdotx >= (zu - usum))
	{
		goto TERMINATE;
	}

	srtrows = (SCi2tuple *)malloc(nrows * sizeof(SCi2tuple));
	for (i = 0; i < nrows; ++i)
	{

		// No need to generate R and T(x) sets: cnt says how
		// much a row is covered and if cnt < 1 (minimum
		// cover) skip the row
		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i + 1] && k < xindlen)
		{
			if (rmatind[j] == xind[k])
			{
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k])
			{
				j++;
				continue;
			}

			if (rmatind[j] > xind[k])
			{
				k++;
				continue;
			}
		}

		if ((u[i] > SC_EPSILON_SMALL) && (cnt > 1))
		{
			srtrows[srtrowslen].a = cnt;
			srtrows[srtrowslen++].b = i;
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
		j = rmatbeg[row];
		k = 0;
		while (j < rmatbeg[row + 1] && k < xindlen)
		{
			if (rmatind[j] == xind[k])
			{
				s[rmatind[j]] += val;
				sdotx += val;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k])
			{
				s[rmatind[j]] += val;
				j++;
				continue;
			}

			if (rmatind[j] > xind[k])
			{
				k++;
				continue;
			}
		}

		while (j < rmatbeg[row + 1])
		{
			s[rmatind[j]] += val;
			j++;
		}

		usum -= val;
		u[row] = 0.0;
	}

	free(srtrows);

TERMINATE:
	return 0;
}

int SCbalasbranchrule1_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
							  const int nrows, const int ncols, const int *xind, int xindlen, double *dj,
							  const double zu, double y, int *rqbeg, int *rqind, int max_branch, int max_singl)
{

	int cnt, qnzcnt, qnsingl, i, it, j, jt, k, h, p, sindlen, txindlen, jindlen, mjindlen;
	int *sind, *txind, *jind, *mjind, *tmp;
	double min, v, v1, v2;

	jind = (int *)malloc(xindlen * sizeof(int));
	mjind = (int *)malloc(nrows * sizeof(int));
	tmp = (int *)malloc((nrows > xindlen ? nrows : xindlen) * sizeof(int));

	// S set defined with column indices
	sindlen = 0;
	sind = (int *)malloc(xindlen * sizeof(int));
	for (j = 0; j < xindlen; ++j)
	{
		if (dj[xind[j]] > SC_EPSILON_SMALL)
		{
			sind[sindlen] = xind[j];
			sindlen++;
		}
	}

	// Tx set defined with row indices
	txindlen = 0;
	txind = (int *)malloc(nrows * sizeof(int));
	for (i = 0; i < nrows; ++i)
	{

		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i + 1] && k < xindlen)
		{
			if (cnt > 1)
			{
				break;
			}

			if (rmatind[j] == xind[k])
			{
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k])
			{
				j++;
				continue;
			}

			if (rmatind[j] > xind[k])
			{
				k++;
				continue;
			}
		}

		if (cnt == 1)
		{
			txind[txindlen] = i;
			txindlen++;
		}
	}

	p = 0;
	qnzcnt = 0;
	qnsingl = 0;
	while (1)
	{

		v1 = -INFINITY;
		v2 = INFINITY;

		for (j = 0; j < sindlen; ++j)
		{
			v1 = (dj[sind[j]] > v1) ? dj[sind[j]] : v1;
			v2 = (dj[sind[j]] >= zu - y) && (dj[sind[j]] < v2) ? dj[sind[j]] : v2;
		}

		v = v1 > v2 ? v2 : v1;

		// fill J
		jindlen = 0;
		for (j = 0; j < sindlen; ++j)
		{
			if ((v - dj[sind[j]]) < SC_EPSILON_SMALL)
			{
				jind[jindlen] = sind[j];
				jindlen++;
			}
		}

		// fill M_J
		mjindlen = 0;
		for (i = 0; i < nrows; ++i)
		{
			tmp[i] = 0;
		}

		for (j = 0; j < jindlen; ++j)
		{
			for (i = cmatbeg[jind[j]]; i < cmatbeg[jind[j] + 1]; ++i)
			{
				tmp[cmatind[i]] = 1;
			}
		}

		for (i = 0; i < nrows; ++i)
		{
			if (tmp[i])
			{
				mjind[mjindlen] = i;
				mjindlen++;
			}
		}

		//printf("s=["); for (i = 0; i < sindlen; ++i) printf("%d ", sind[i]); printf("]\n");
		//printf("slen=%d, txlen=%d, jlen=%d, mjlen=%d\n", sindlen, txindlen, jindlen, mjindlen);

		// fill Q_p
		rqbeg[p] = qnzcnt;
		for (j = 0; j < ncols; ++j)
		{
			rqind[qnzcnt] = j * (dj[j] >= v);
			qnzcnt += (dj[j] >= v);
		}
		rqbeg[p + 1] = qnzcnt;

		it = -1;
		min = INFINITY;
		i = 0;
		h = 0;
		while (i < txindlen && h < mjindlen)
		{
			if (txind[i] == mjind[h])
			{

				cnt = 0;
				j = rmatbeg[txind[i]];
				k = rqbeg[p];
				while (j < rmatbeg[txind[i] + 1] && k < rqbeg[p + 1])
				{
					if (rmatind[j] == rqind[k])
					{
						j++;
						k++;
						continue;
					}

					if (rmatind[j] < rqind[k])
					{
						cnt++;
						j++;
						continue;
					}

					if (rmatind[j] > rqind[k])
					{
						k++;
						continue;
					}
				}

				while (j < rmatbeg[txind[i] + 1])
				{
					cnt++;
					j++;
				}

				if (cnt < min)
				{
					min = cnt;
					it = txind[i];
				}

				i++;
				h++;
				continue;
			}

			if (txind[i] < mjind[h])
			{
				i++;
				continue;
			}

			if (txind[i] > mjind[h])
			{
				h++;
				continue;
			}
		}

		jt = -1;
		j = 0;
		k = rmatbeg[it];
		while (j < jindlen && k < rmatbeg[it + 1])
		{
			if (jind[j] == rmatind[k])
			{
				jt = jind[j];
				break;
			}

			if (jind[j] < rmatind[k])
			{
				j++;
				continue;
			}

			if (jind[j] > rmatind[k])
			{
				k++;
				continue;
			}
		}

		//printf("j=["); for (i = 0; i < jindlen; ++i) printf("%d ", jind[i]); printf("]\n");
		//printf("r=["); for (i = rmatbeg[it]; i < rmatbeg[it+1]; ++i) printf("%d ", rmatind[i]); printf("]\n");

		y += dj[jt];
		//printf("p=%d, v=%lf, dj=%lf, it=%d, jt=%d, zu=%lf, y=%lf\n", p, v, dj[jt], it, jt, zu, y);

		cnt = 0;
		j = rmatbeg[it];
		k = rqbeg[p];
		while (j < rmatbeg[it + 1] && k < rqbeg[p + 1])
		{
			if (rmatind[j] == rqind[k])
			{
				tmp[cnt] = rmatind[j];

				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < rqind[k])
			{
				j++;
				continue;
			}

			if (rmatind[j] > rqind[k])
			{
				k++;
				continue;
			}
		}

		for (j = 0; j < cnt; ++j)
		{
			rqind[rqbeg[p] + j] = tmp[j];
		}
		qnzcnt = qnzcnt - (rqbeg[p + 1] - rqbeg[p]) + cnt;
		rqbeg[p + 1] = rqbeg[p] + cnt;

		qnsingl += (cnt == 1);

		// check
		if ((y + SC_EPSILON_SMALL < zu) && (p + 1 >= max_branch))
		{
			p = -1;
			//printf("max branch reached\n");
			break;
		}
		if ((y + SC_EPSILON_SMALL < zu) && (qnsingl >= max_singl))
		{
			p = -1;
			//printf("max singl num reached\n");
			break;
		}
		if (y + SC_EPSILON_SMALL >= zu)
		{
			p++;
			//printf("OK!\n");
			break;
		}

		for (j = 0; j < sindlen; ++j)
		{
			if (sind[j] == jt)
			{
				sindlen--;
				memcpy(&sind[j], &sind[j + 1], (sindlen - j) * sizeof(int));
				break;
			}
		}

		p++;
	}

	free(sind);
	free(txind);
	free(jind);
	free(mjind);
	free(tmp);

	if (p < 0)
	{
		return p;
	}
	return (qnzcnt > p * log2(p)) ? p : -1;
}*/