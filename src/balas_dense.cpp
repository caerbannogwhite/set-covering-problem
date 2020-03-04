
#include "balas_dense.hpp"

int SCbalascutprocedure(char *mat, int nrows, int ncols,
						const char *x, double *s, double zu, double y,
						int *cutind)
{

	int cnt, i, it, j, jt, off, p = 0;
	double min, v, v1, v2;

	// init W, J, M_J
	char *wset = (char *)calloc(ncols, sizeof(char));
	char *qset = (char *)malloc(ncols * sizeof(char));
	char *jset = (char *)malloc(ncols * sizeof(char));
	char *mjset = (char *)malloc(nrows * sizeof(char));

	// init set S
	char *sset = (char *)malloc(ncols);
	for (j = ncols; j--;)
	{
		sset[j] = x[j] & (s[j] > 0);
	}

	char *txset = (char *)malloc(nrows);
	for (i = nrows; i--;)
	{
		cnt = 0;
		off = ncols * i;
		for (j = ncols; j--;)
		{
			cnt += mat[off + j] & x[j];
		}

		txset[i] = (cnt == 1);
	}

	char *matr;

	while (1)
	{

		v1 = -INFINITY;
		v2 = INFINITY;
		for (j = ncols; j--;)
		{
			v1 = sset[j] && (s[j] > v1) ? s[j] : v1;
			v2 = sset[j] && (s[j] >= zu - y) && (s[j] < v2) ? s[j] : v2;
		}

		v = v1 > v2 ? v2 : v1;

		// fill J
		for (j = ncols; j--;)
		{
			jset[j] = (sset[j] && (abs(v - s[j]) < SC_EPSILON_SMALL));
		}
		//printf("bcp1 ok 3 v=%lf p=%d y=%lf zu=%lf\n", v, p, y, zu);

		// fill M_J
		for (i = nrows; i--;)
		{
			mjset[i] = 0;
		}

		for (j = ncols; j--;)
		{
			if (jset[j])
			{
				for (i = nrows; i--;)
				{
					mjset[i] |= mat[ncols * i + j];
				}
			}
		}

		// fill Q_p
		for (j = 0; j < ncols; ++j)
		{
			qset[j] = (s[j] >= v);
		}

		it = -1;
		min = INFINITY;
		for (i = nrows; i--;)
		{
			if (txset[i] & mjset[i])
			{
				cnt = 0;
				matr = &mat[i * ncols];
				for (j = ncols; j--;)
				{
					cnt += matr[j] && !(wset[j] | qset[j]);
				}

				if (cnt < min)
				{
					min = cnt;
					it = i;
				}
			}
		}
		//printf("bcp1 ok 4 it=%d nrows=%d\n", it, nrows);

		jt = -1;
		matr = &mat[it * ncols];
		for (j = ncols; j--;)
		{
			if (jset[j] & matr[j])
			{
				jt = j;
				break;
			}
		}
		//printf("bcp1 ok 5 p=%2d it=%3d min=%.1lf jt=%3d\n", p, it, min, jt);

		y += s[jt];

		cnt = 0;
		matr = &mat[it * ncols];
		for (j = ncols; j--;)
		{
			wset[j] = wset[j] | (matr[j] & !(qset[j]));
			s[j] -= s[jt] * (qset[j] & matr[j]);
		}

		// check
		if (y + SC_EPSILON_MED >= zu)
		{
			p++;
			break;
		}

		sset[jt] = 0;
		p++;
	}

	cnt = 0;
	for (j = 0; j < ncols; ++j)
	{
		cutind[cnt] = j;
		cnt += wset[j];
	}

	free(wset);
	free(qset);
	free(jset);
	free(mjset);
	free(sset);
	free(txset);

	return cnt;
}

int SCbalasbranchrule1(char *mat, int nrows, int ncols,
					   const char *x, double *s, double zu, double y,
					   char *qmat, int max_branch, int max_singl)
{

	char cond;
	int cnt, nzcnt, nsingl, i, it, j, jt, off, p = 0;
	double min, v, v1, v2;

	// init J, M_J
	char *jset = (char *)malloc(ncols);
	char *mjset = (char *)malloc(nrows);

	// init set S
	char *sset = (char *)malloc(ncols);
	for (j = ncols; j--;)
	{
		sset[j] = (x[j]) && (s[j] > 0);
	}

	char *txset = (char *)malloc(nrows);
	for (i = nrows; i--;)
	{
		cnt = 0;
		off = ncols * i;
		for (j = ncols; j--;)
		{
			cnt += mat[off + j] && x[j];
		}

		txset[i] = (cnt == 1);
	}

	char *qmatr = qmat;
	char *matr;

	nzcnt = 0;
	nsingl = 0;
	while (1)
	{

		v1 = -INFINITY;
		v2 = INFINITY;
		for (j = ncols; j--;)
		{
			v1 = sset[j] && (s[j] > v1) ? s[j] : v1;
			v2 = sset[j] && (s[j] >= (zu - y)) && (s[j] < v2) ? s[j] : v2;
		}

		v = v1 > v2 ? v2 : v1;

		// fill J
		for (j = ncols; j--;)
		{
			jset[j] = (sset[j] && (abs(v - s[j]) < SC_EPSILON_SMALL));
		}
		//printf("bb1 ok 3 v=%lf p=%d y=%lf zu=%lf\n", v, p, y, zu);

		// fill M_J
		for (i = nrows; i--;)
		{
			mjset[i] = 0;
		}

		for (j = ncols; j--;)
		{
			if (jset[j])
			{
				for (i = nrows; i--;)
				{
					mjset[i] |= mat[ncols * i + j];
				}
			}
		}

		// fill Q_p
		for (j = 0; j < ncols; ++j)
		{
			qmatr[j] = (s[j] >= v);
		}

		it = -1;
		min = INFINITY;
		for (i = nrows; i--;)
		{
			if (txset[i] & mjset[i])
			{
				cnt = 0;
				matr = &mat[i * ncols];
				for (j = ncols; j--;)
				{
					cnt += matr[j] && !qmatr[j];
				}

				if (cnt < min)
				{
					min = cnt;
					it = i;
				}
			}
		}
		//printf("bb1 ok 4 it=%d nrows=%d\n", it, nrows);

		jt = -1;
		matr = &mat[it * ncols];
		for (j = ncols; j--;)
		{
			if (jset[j] & matr[j])
			{
				jt = j;
				break;
			}
		}
		//printf("bb1 ok 5 p=%2d it=%3d min=%.1lf jt=%3d\n", p, it, min, jt);

		y += s[jt];

		cnt = 0;
		matr = &mat[it * ncols];
		for (j = ncols; j--;)
		{
			cond = (qmatr[j] && matr[j]);
			qmatr[j] = cond;
			cnt += cond;

			s[j] -= s[jt] * cond;
		}

		nsingl += (cnt == 1);
		nzcnt += cnt;

		// check
		if ((y + SC_EPSILON_SMALL < zu) && (p + 1 == max_branch))
		{
			return -1;
		}
		if ((y + SC_EPSILON_SMALL < zu) && (nsingl >= max_singl))
		{
			p = -1;
			break;
		}
		if (y + SC_EPSILON_SMALL >= zu)
		{
			p++;
			break;
		}

		sset[jt] = 0;
		p++;
		qmatr = &qmat[p * ncols];
	}

	free(jset);
	free(mjset);
	free(sset);
	free(txset);

	return (nzcnt > p * log2(p)) ? p : -1;
}

int SCbalasbranchrule1_test(char *mat, int nrows, int ncols,
							const char *x, double *s, double zu, double y,
							char *qmat, int max_branch, int max_singl)
{

	#if DEBUG_VERBOSITY
	char append = 'a';
	int i;
	double acc;
	double *s_cp = (double *)malloc(ncols * sizeof(double));
	double *sj = (double *)malloc(nrows * sizeof(double));

	memcpy(s_cp, s, ncols * sizeof(double));
	for (i = 0; i < nrows; ++i)
		sj[i] = 0;

	FILE *log = fopen("debug.log", &append);
	printf("\nSCbalasbranchrule1_test\n");
	if (log != NULL)
		fprintf(log, "\nSCbalasbranchrule1_test\n");
	#endif

	char cond;
	char *qmatr = qmat;
	int jt;
	int j, p = 0, cnt, nzcnt, nsingl;
	double v;

	// init set S
	char *sset = (char *)malloc(ncols);
	for (j = ncols; j--;)
	{
		sset[j] = x[j] && (s[j] > 0);
	}

	nzcnt = 0;
	nsingl = 0;
	while (1)
	{

		v = 0.0;
		jt = -1;
		for (j = ncols; j--;)
		{
			if (sset[j] && (s[j] > v))
			{
				v = s[j];
				jt = j;
			}
		}

	#if DEBUG_VERBOSITY
		if (jt == -1)
		{
			printf("\tjt == -1! p=%d - Closing...\n", p);
			if (log != NULL)
				fprintf(log, "\tjt == -1! p=%d - Closing...\n", p);
			if (log != NULL)
			{
				fprintf(log, "sx=");
				for (i = 0; i < ncols; ++i)
				{
					if (sset[i])
					{
						fprintf(log, "(%d,%lf) ", i, s[i]);
					}
				}
			}
			if (log != NULL)
				fclose(log);
			exit(SC_ERR_BALAS_BRANCH_RULE_1);
		}

		if (fabs(v - 0.0) < SC_EPSILON_SMALL)
		{
			printf("\tv == 0! Closing...\n");
			if (log != NULL)
				fprintf(log, "\tv == 0!\n");
			if (log != NULL)
				fclose(log);
			exit(SC_ERR_BALAS_BRANCH_RULE_1);
		}
	#endif

		// fill Q_p
		cnt = 0;
		for (j = 0; j < ncols; ++j)
		{
			cond = (s[j] >= v) && (!sset[j]);
			qmatr[j] = cond;
			cnt += cond;
		}

		cnt += (!qmatr[jt]);
		qmatr[jt] = 1;

		nsingl += (cnt == 1);
		nzcnt += cnt;

		y += v;

	#if DEBUG_VERBOSITY
		sj[p] = v;
	#endif

		// check
		if ((y + SC_EPSILON_SMALL < zu) && (p + 1 == max_branch))
		{
			p = -1;
			break;
		}
		if ((y + SC_EPSILON_SMALL < zu) && (nsingl >= max_singl))
		{
			p = -1;
			break;
		}
		if (y + SC_EPSILON_SMALL >= zu)
		{
			p++;
			break;
		}

		for (j = ncols; j--;)
		{
			s[j] -= v * qmatr[j];
		}

		//s[jt] = 0;
		sset[jt] = 0;
		p++;
		qmatr = &qmat[p * ncols];
	}

	free(sset);

	#if DEBUG_VERBOSITY
	if (p > -1)
	{
		for (j = 0; j < ncols; ++j)
		{
			acc = 0;
			for (i = 0; i < p; ++i)
			{
				acc += sj[j] * qmat[i * ncols + j];
			}

			if (acc > s_cp[j])
			{
				printf("\tBalas condition (7) violated! Closing...\n");
				if (log != NULL)
					fprintf(log, "\tBalas condition (7) violated!\n");
				if (log != NULL)
					fclose(log);
				exit(SC_ERR_BALAS_COND_VIOLATED);
			}
		}
	}

	free(s_cp);
	free(sj);
	if (log != NULL)
		fclose(log);
	#endif

	return (nzcnt > p * log2(p)) ? p : -1;
}

/**
 * Given a SCP matrix and a SCP solution x,
 * return the number of removed columns
 * from x (to make x a prime cover).
 * 
 * @param mat - ublas::matrix<double> a SCP matrix
 * @param x - ublas::vector<double> a SCP solution
 * @return the number of removed columns from x
 */
int baldns_make_prime_cover(const ublas::matrix<double> &mat, const ublas::vector<double> &x)
{
	int cntRemoved;
	size_t j;

	ublas::vector<double> matDotX = ublas::prod(mat, x);

	cntRemoved = 0;
	for (j = mat.size2; j--;)
	{
		if (x(j) > SC_EPSILON_SMALL)
		{
			x(j) = 0.0;
			for (i = 0; i < nrows; ++i)
			{
				if ((mat(i, j) > SC_EPSILON_SMALL) && (matDotX(i) < 2.0))
				{
					x(j) = 1.0;
					break;
				}
			}

			if (x(j) < SC_EPSILON_SMALL)
			{
				cntRemoved++;
				matDotX -= ublas::matrix_column(mat, j);
			}
		}
	}

	delete matDotX;

	return cntRemoved;
}

int SCoversatrows(char *mat, const int nrows,
				  const int ncols, char *x, int *xsupp,
				  int *xsupplen_p)
{

	char *matr;
	int cnt, i, j, k;

	int *matdotx = (int *)malloc(nrows * sizeof(int));
	for (i = 0; i < nrows; ++i)
	{
		matr = &mat[i * ncols];
		matdotx[i] = 0;
		for (j = 0; j < ncols; ++j)
		{
			matdotx[i] += matr[j] & x[j];
		}
	}

	//printf("osr ok2 *xsupplen_p=%d\n", *xsupplen_p);
	for (k = *xsupplen_p; k--;)
	{
		cnt = 0;
		for (i = 0; i < nrows; ++i)
		{
			cnt += mat[i * ncols + xsupp[k]] & (matdotx[i] > 1);
		}
		//printf("	osr ok25 k=%d xsupp[k]=%d cnt=%d *xsupplen_p=%d\n", k, xsupp[k], cnt, *xsupplen_p);
		if (cnt > 0)
		{
			for (i = 0; i < nrows; ++i)
			{
				matdotx[i] -= mat[i * ncols + xsupp[k]];
			}
			x[xsupp[k]] = 0;
			*xsupplen_p = *xsupplen_p - 1;
			memcpy(&xsupp[k], &xsupp[k + 1], (*xsupplen_p - k) * sizeof(int));
		}
	}
	//printf("osr ok3 *xsupplen_p=%d\n", *xsupplen_p);
	free(matdotx);

	return 0;
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
	ublas::vector<double> res = ublas::prod(mat, x);

	flag = true;
	for (auto it = res.cbegin(); it != it != res.cend(); ++it)
	{
		if (*it < 1.0)
		{
			flag = false;
		}
	}	

	delete res;

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
	int col, row, cnt, rcnt;
	size_t i, j, coveredRows;
	double val, zUpp;
	double (*func)(const double, const int);
	char *matr;

	std::vector<std::pair<int, int>> rSet = std::vector<std::pair<int, int>>(mat.size1);

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
			func = [](const double c, const double k) -> double { return math.sqrt(c) / (k * k); };
			break;
		default:
			func = [](const double c, const double k) -> double { return k < 2 ? c : c / log2(k); };
			break;
	}

	// get the starting value of the current solution x
	zUpp = ublas::prod(obj, x);

	// find the number of already covered rowa
	coveredRows = 0;
	for (i = mat.size1; i--;)
	{
		cnt = math.round(ublas::prod(mat(i), x));
		rcnt = math.round(ublas::vector_sum(mat(i)));

		if (cnt == 0)
		{
			rSet.push_back(std::make_pair(rcnt, i));
		}
		else
		{
			coveredRows++;
		}
	}

	// sort not covered rows by decreasing number of
	// non-zero elements in the row (read rSet from back to front)
	std::sort(rSet.begin(), rSet.end(), [](std::pair<int, int> a, std::pair<int, int> b) -> bool { return a.first > b.first; });

	while (coveredRows < mat.size1)
	{
		// get last element of R
		row = rSet.back().second;
		rSet.pop_back();

		val = INFINITY;
		col = -1;
		for (j = 0; j < mat.size2; ++j)
		{

			if (fabs(mat(row, j)) < SC_EPSILON_SMALL)
			{
				continue;
			}

			cnt = 0;
			for (auto it = rSet.cbegin(); it != rSet.cend(); ++it)
			{
				cnt += mat[ncols * (*it).second + j]);
			}

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
		for (i = 0; i < rsetlen; ++i)
		{
			if ((rset[i].b != -1) && (mat[rset[i].b * ncols + col]))
			{
				rset[i].b = -1;
				coveredRows++;
			}
		}
	}

	// make the cover found a prime cover
	SCprimecover(mat, x);

	for (i = 0; i < *xsupplen; ++i)
	{
		if (!x[xsupp[i]])
		{
			zUpp -= obj[xsupp[i]];
			*xsupplen = *xsupplen - 1;
			memcpy(&xsupp[i], &xsupp[i + 1], *xsupplen * sizeof(int));
			i--;
		}
	}

	delete rSet;

	return zUpp;
}

double SCbalasheurprimal12(char *mat, const double *costs,
						   const int nrows, const int ncols,
						   char *x, int *xsupp, int *xsupplen_p)
{

	char *x1, *x2, *xcp;
	int *xsupp1, *xsupp2, *xsuppcp;
	int xsupplen1, xsupplen2, xsupplencp;
	double zu, zu1, zu2;

	zu = SCbalasheurprimal0(mat, costs, nrows, ncols, x, xsupp, xsupplen_p, 3);

	x1 = (char *)malloc(ncols * sizeof(char));
	x2 = (char *)malloc(ncols * sizeof(char));
	xcp = (char *)malloc(ncols * sizeof(char));
	xsupp1 = (int *)malloc(ncols * sizeof(int));
	xsupp2 = (int *)malloc(ncols * sizeof(int));
	xsuppcp = (int *)malloc(ncols * sizeof(int));

	// First round
	memcpy(xcp, x, ncols * sizeof(char));
	memcpy(xsuppcp, xsupp, *xsupplen_p * sizeof(int));
	xsupplencp = *xsupplen_p;

	SCoversatrows(mat, nrows, ncols, x, xsupp, xsupplen_p);

	memcpy(x1, x, ncols * sizeof(char));
	memcpy(x2, x, ncols * sizeof(char));
	memcpy(xsupp1, xsupp, *xsupplen_p * sizeof(int));
	memcpy(xsupp2, xsupp, *xsupplen_p * sizeof(int));
	xsupplen1 = *xsupplen_p;
	xsupplen2 = *xsupplen_p;
	zu1 = SCbalasheurprimal0(mat, costs, nrows, ncols, x1, xsupp1, &xsupplen1, 1);
	zu2 = SCbalasheurprimal0(mat, costs, nrows, ncols, x2, xsupp2, &xsupplen2, 2);
	//printf("zu=(%lf %lf %lf) xsupplen=(%d %d %d)\n", zu, zu1, zu2, *xsupplen_p, xsupplen1, xsupplen2);
	if ((zu1 < zu2) && (zu1 < zu))
	{
		memcpy(x, x1, ncols * sizeof(char));
		memcpy(xsupp, xsupp1, xsupplen1 * sizeof(int));
		*xsupplen_p = xsupplen1;
		zu = zu1;
	}
	else if ((zu2 < zu1) && (zu2 < zu))
	{
		memcpy(x, x2, ncols * sizeof(char));
		memcpy(xsupp, xsupp2, xsupplen2 * sizeof(int));
		*xsupplen_p = xsupplen2;
		zu = zu2;
	}
	else
	{
		memcpy(x, xcp, ncols * sizeof(char));
		memcpy(xsupp, xsuppcp, xsupplencp * sizeof(int));
		*xsupplen_p = xsupplencp;
	}
	//printf("len=%d [", *xsupplen_p); for (int i = 0; i < *xsupplen_p; ++i) printf("%d ", xsupp[i]); printf("]\n");

	// Second round
	memcpy(xcp, x, ncols * sizeof(char));
	memcpy(xsuppcp, xsupp, *xsupplen_p * sizeof(int));
	xsupplencp = *xsupplen_p;

	SCoversatrows(mat, nrows, ncols, x, xsupp, xsupplen_p);

	memcpy(x1, x, ncols * sizeof(char));
	memcpy(x2, x, ncols * sizeof(char));
	memcpy(xsupp1, xsupp, *xsupplen_p * sizeof(int));
	memcpy(xsupp2, xsupp, *xsupplen_p * sizeof(int));
	xsupplen1 = *xsupplen_p;
	xsupplen2 = *xsupplen_p;
	zu1 = SCbalasheurprimal0(mat, costs, nrows, ncols, x1, xsupp1, &xsupplen1, 4);
	zu2 = SCbalasheurprimal0(mat, costs, nrows, ncols, x2, xsupp2, &xsupplen2, 5);
	//printf("zu=(%lf %lf %lf) xsupplen=(%d %d %d)\n", zu, zu1, zu2, *xsupplen_p, xsupplen1, xsupplen2);
	if ((zu1 < zu2) && (zu1 < zu))
	{
		memcpy(x, x1, ncols * sizeof(char));
		memcpy(xsupp, xsupp1, xsupplen1 * sizeof(int));
		*xsupplen_p = xsupplen1;
		zu = zu1;
	}
	else if ((zu2 < zu1) && (zu2 < zu))
	{
		memcpy(x, x2, ncols * sizeof(char));
		memcpy(xsupp, xsupp2, xsupplen2 * sizeof(int));
		*xsupplen_p = xsupplen2;
		zu = zu2;
	}
	else
	{
		memcpy(x, xcp, ncols * sizeof(char));
		memcpy(xsupp, xsuppcp, xsupplencp * sizeof(int));
		*xsupplen_p = xsupplencp;
	}
	//printf("len=%d [", *xsupplen_p); for (int i = 0; i < *xsupplen_p; ++i) printf("%d ", xsupp[i]); printf("]\n");

	free(x1);
	free(x2);
	free(xcp);
	free(xsupp1);
	free(xsupp2);
	free(xsuppcp);

	return zu;
}

double SCbalasheurprimal5b(char *mat, const double *costs,
						   const int nrows, const int ncols, char *x,
						   double *u, double *s)
{

	char *matr;
	int *rset;
	int cvdrows, cnt, i, it, j, rsetlen = 0, rsetsize = nrows;
	double v, zu;

	for (j = 0; j < ncols; ++j)
	{
		x[j] = s[j] < SC_EPSILON_SMALL ? 1 : 0;
	}

	SCprimecover(mat, nrows, ncols, x);

	rset = (int *)malloc(rsetsize * sizeof(int));
	cvdrows = 0;
	for (i = nrows; i--;)
	{
		cnt = 0;
		matr = &mat[ncols * i];
		for (j = ncols; j--;)
		{
			if (matr[j] & x[j])
			{
				cnt = 1;
				break;
			}
		}

		if (!cnt)
		{
			rset[rsetlen++] = i;
		}
		else
		{
			cvdrows++;
		}
	}
	//printf("bhp5 ok4 cvdrows=%d rsetlen=%d\n", cvdrows, rsetlen);

	while (cvdrows < nrows)
	{

		for (i = 0; i < rsetlen; ++i)
		{
			if (rset[i] > -1)
			{
				it = rset[i];
				break;
			}
		}
		//printf("bhp5 ok5 it=%d cvdrows=%d rsetlen=%d\n", it, cvdrows, rsetlen);
		matr = &mat[it * ncols];
		v = INFINITY;
		for (j = 0; j < ncols; ++j)
		{
			v = matr[j] && (s[j] < v) ? s[j] : v;
		}
		//printf("bhp5 ok6 v=%lf\n", v);
		u[it] += v;
		for (j = 0; j < ncols; ++j)
		{
			s[j] -= v * matr[j];
		}

		for (j = 0; j < ncols; ++j)
		{
			x[j] = s[j] < SC_EPSILON_SMALL ? 1 : 0;
		}

		SCprimecover(mat, nrows, ncols, x);

		rsetlen = 0;
		cvdrows = 0;
		for (i = nrows; i--;)
		{
			cnt = 0;
			matr = &mat[ncols * i];
			for (j = ncols; j--;)
			{
				if (matr[j] & x[j])
				{
					cnt = 1;
					break;
				}
			}

			if (!cnt)
			{
				rset[rsetlen++] = i;
			}
			else
			{
				cvdrows++;
			}
		}

		//printf("bhp5 ok8 cvdrows=%d rsetlen=%d\n", cvdrows, rsetlen);
	}

	cnt = SCprimecover(mat, nrows, ncols, x);
	//printf("bhp5 ok9 remrows=%d\n", cnt);

	zu = 0.0;
	for (j = 0; j < ncols; ++j)
	{
		zu += x[j] * costs[j];
	}

	return zu;
}

int SCisdualsolution(const char *mat, const double *obj,
					 const int nrows, const int ncols, double *u)
{

	#if DEBUG_VERBOSITY
	char append = 'a';
	FILE *log = fopen("debug.log", &append);
	printf("\nSCisdualsolution\n");
	fprintf(log, "\nSCisdualsolution\n");
	#endif

	int i, j;
	double rowsum;

	for (j = 0; j < ncols; ++j)
	{
		rowsum = 0;
		for (i = 0; i < nrows; ++i)
		{
			rowsum += mat[i * ncols + j] * u[i];
		}

		if ((rowsum - SC_EPSILON_MED) > obj[j])
		{

	#if DEBUG_VERBOSITY
			if (DEBUG_VERBOSITY > 4)
			{
				printf("printing mat, obj, u on log file.\n");
				fprintf(log, "\tm=%d\n\tn=%d\n", nrows, ncols);
				fprintf(log, "\tmat=[");
				for (i = 0; i < nrows * ncols - 1; ++i)
					fprintf(log, "%d, ", mat[i]);
				fprintf(log, "%d]\n", mat[nrows * ncols - 1]);
				fprintf(log, "\tobj=[");
				for (i = 0; i < ncols - 1; ++i)
					fprintf(log, "%.1lf, ", obj[i]);
				fprintf(log, "%.1lf]\n", obj[ncols - 1]);
				fprintf(log, "\tu=[");
				for (i = 0; i < nrows - 1; ++i)
					fprintf(log, "%lf, ", u[i]);
				fprintf(log, "%lf]\n", u[nrows - 1]);
			}
	#endif

			return 0;
		}
	}

	#if DEBUG_VERBOSITY
	fclose(log);
	#endif

	return 1;
}

int SCdual0(const char *mat, double *obj, const int nrows,
			const int ncols, double *u, double *s)
{

	char cont = 'C', le = 'L';
	int error, i, j, off;
	double one = 1.0, zero = 0.0;

	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "dual");

	// COLUMNS # cols_dual == # rows prim
	for (i = 0; i < nrows; ++i)
	{
		if (CPXnewcols(env, lp, 1, &one, &zero, NULL, &cont, NULL))
			perror("Can't use CPXnewcols");
	}

	// ROWS
	for (j = 0; j < ncols; ++j)
	{

		if (CPXnewrows(env, lp, 1, &obj[j], &le, NULL, NULL))
			perror("Can't use CPXnewrows");
		for (i = 0; i < nrows; ++i)
		{
			if (mat[i * ncols + j])
			{
				if (CPXchgcoef(env, lp, j, i, 1.0))
					perror("Can't use CPXchgcoef");
			}
		}
	}

	// Maximize
	CPXchgobjsen(env, lp, CPX_MAX);

	CPXmipopt(env, lp);

	// Dual vector
	CPXgetx(env, lp, u, 0, nrows - 1);

	for (j = 0; j < ncols; ++j)
	{
		s[j] = obj[j];
	}

	for (i = 0; i < nrows; ++i)
	{
		off = i * ncols;
		for (j = 0; j < ncols; ++j)
		{
			s[j] -= mat[off + j] * u[i];
		}
	}

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0;
}

int SCbalasheurdual1(ublas::matrix<double> &mat,
					 ublas::vactor<double> &x, ublas::vector<double> &u, ublas::vector<double> &s)
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

int SCbalasheurdual3(ublas::matrix<double> &mat,
					 ublas::vactor<double> &x, ublas::vector<double> &u, ublas::vector<double> &s, const double zUpp)
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
}
