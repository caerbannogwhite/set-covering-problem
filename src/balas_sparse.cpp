
#include "balas_sparse.hpp"

int SCprimecover_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg,
		const int *cmatind, const int nrows, const int ncols, int *xind,
		const int xindlen) {

	char flag;
	int i, j, k, removedcnt;
	int *matdotx = (int *) malloc(nrows * sizeof(int));

	for (i = 0; i < nrows; ++i) {
		matdotx[i] = 0;
		j = rmatbeg[i];
		k = 0;
		while ((j < rmatbeg[i+1]) && (k < xindlen)) {
			if (rmatind[j] == xind[k]) {
				matdotx[i]++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k]) {
				j++;
				continue;
			}

			if (rmatind[j] > xind[k]) {
				k++;
				continue;
			}
		}
	}

	removedcnt = 0;
	for (j = xindlen - 1; j > 0; --j) {
		flag = 1;
		for (k = cmatbeg[xind[j]]; k < cmatbeg[xind[j] + 1]; ++k) {
			if (matdotx[cmatind[k]] < 2) {
				flag = 0;
				break;
			}
		}

		if (flag) {
			for (k = cmatbeg[xind[j]]; k < cmatbeg[xind[j] + 1]; ++k) {
				matdotx[cmatind[k]]--;
			}

			xind[j] = -1;
			removedcnt++;
		}
	}
	free(matdotx);

	return removedcnt;
}


double SCbalasheurprimal0_sparse(const int *rmatbeg, const int *rmatind, const int *cmatbeg, const int *cmatind,
		const double *obj, const int nrows, const int ncols, int *xind, int *xindlen,
		const int whichfunc) {

	char flag;
	int jt, it;
	int *item, *rset, *rsetnew;
	int i, j, k, cnt, cvdrows, rsetlen;
	double v, zu;
	double (*func)(const double, const int);
	SCi2tuple *rsettup;

	switch (whichfunc) {
		case 1: func = &func1;
			break;
		case 2: func = &func2;
			break;
		case 3: func = &func3;
			break;
		case 4: func = &func4;
			break;
		case 5: func = &func5;
			break;
		default: func = &func3;
			break;
	}

	/*cmatbeg = (int *) malloc((ncols + 1) * sizeof(int));
	cmatind = (int *) malloc(rmatbeg[nrows] * sizeof(int));

	cnt = 0;
	for (j = 0; j < ncols; ++j) {
		cmatbeg[j] = cnt;
		for (i = 0; i < nrows; ++i) {
			item = bsearch(&j, &rmatind[rmatbeg[i]], rmatbeg[i+1] - rmatbeg[i], sizeof(int), SCint_cmp);
			if (item != NULL) {
				cmatind[cnt] = i;
				cnt++;
			}
		}
	}
	cmatbeg[ncols] = cnt;*/

	zu = 0.0;
	rset = (int *) malloc(nrows * sizeof(int));
	rsetnew = (int *) malloc(nrows * sizeof(int));
	rsettup = (SCi2tuple *) malloc(nrows * sizeof(SCi2tuple));

	rsetlen = 0;
	cvdrows = 0;

	if ((*xindlen) > 0) {
		for (j = 0; j < (*xindlen); ++j) {
			zu += obj[xind[j]];
		}

		for (i = 0; i < nrows; ++i) {

			flag = 0;
			j = 0;
			k = rmatbeg[i];
			while ((j < (*xindlen)) && (k < rmatbeg[i + 1])) {
				if (xind[j] == rmatind[k]) {
					flag = 1;
					break;
				}

				if (xind[j] < rmatind[k]) {
					j++;
					continue;
				}

				if (xind[j] > rmatind[k]) {
					k++;
					continue;
				}
			}

			if (!flag) {
				rsettup[rsetlen].a = rmatbeg[i + 1] - rmatbeg[i];
				rsettup[rsetlen].b = i;
				rsetlen++;
			} else {
				cvdrows++;
			}
		}
	} else {
		for (i = 0; i < nrows; ++i) {
			rsettup[i].a = rmatbeg[i + 1] - rmatbeg[i];
			rsettup[i].b = i;
		}
		rsetlen = nrows;
	}

	qsort(rsettup, rsetlen, sizeof(SCi2tuple), SCi2tuple_cmpa);
	for (i = 0; i < rsetlen; ++i) {
		rset[i] = rsettup[i].b;
	}

	free(rsettup);

	while (cvdrows < nrows) {

		it = rset[0];

		v = INFINITY;
		jt = -1;
		for (j = 0; j < ncols; ++j) {

			flag = 0;
			cnt = 0;
			for (i = 0; i < rsetlen; ++i) {
				if (rset[i] != -1) {
					item = (int *) bsearch(&rset[i], &cmatind[cmatbeg[j]], cmatbeg[j+1] - cmatbeg[j], sizeof(int), SCint_cmp);
					if (item != NULL) {
						cnt++;
						if (rset[i] == it) {
							flag = 1;
						}
					}
				}
			}

			if (!flag) {
				continue;
			}

			if (func(obj[j], cnt) < v) {
				v = func(obj[j], cnt);
				jt = j;
			}
		}

		/*if (jt == -1) {
			free(rset);
			return -1;
		}*/

		xind[*xindlen] = jt;
		*xindlen += 1;
		zu += obj[jt];

		cnt = 0;
		for (i = 0; i < rsetlen; ++i) {
			item = (int *) bsearch(&rset[i], &cmatind[cmatbeg[jt]], cmatbeg[jt+1] - cmatbeg[jt], sizeof(int), SCint_cmp);
			if (item == NULL) {
				rsetnew[cnt] = rset[i];
				cnt++;
			} else {
				cvdrows++;
			}
		}

		memcpy(rset, rsetnew, cnt * sizeof(int));
		rsetlen = cnt;
	}

	free(rset);
	free(rsetnew);

	// Make the cover a prime cover
	SCprimecover_sparse(rmatbeg, rmatind, cmatbeg, cmatind, nrows, ncols, xind, *xindlen);

	for (i = 0; i < *xindlen; ++i) {
		if (xind[i] == -1) {
			zu -= obj[xind[i]];
			*xindlen = *xindlen - 1;
			memcpy(&xind[i], &xind[i+1], (*xindlen) * sizeof(int));
			i--;
		}
	}

	qsort(xind, *xindlen, sizeof(int), SCint_cmp);

	return zu;
}


int SCdual0_sparse(const int *rmatbeg, const int *rmatind, double *obj, const int nrows,
		const int ncols, double *u, double *s) {

	char le = 'L';
	int error, i, j, nzcnt, status;
	double *cmatval;

	nzcnt = rmatbeg[nrows];
	cmatval = (double *) malloc(nzcnt * sizeof(double));
	for (i = 0; i < nzcnt; ++i) {
		cmatval[i] = 1.0;
	}

	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "dual");

	CPXaddcols(env, lp, (int) nrows, nzcnt, cmatval, rmatbeg, rmatind, cmatval, NULL, NULL, NULL);

	for (j = 0; j < ncols; ++j) {
		CPXnewrows(env, lp, 1, &obj[j], &le, NULL, NULL);
	}

	// Maximize
	CPXchgobjsen(env, lp, CPX_MAX);
	CPXchgprobtype(env, lp, CPXPROB_LP);

	status = CPXlpopt(env, lp);
	if (status) { printf("SCdual0_sparse - CPXmipopt error: %d\n", status); }

	// Dual vector
	status = CPXgetx(env, lp, u, 0, (int) (nrows - 1));
	if (status) { printf("SCdual0_sparse - CPXgetx error: %d\n", status); }

	for (j = 0; j < ncols; ++j) {
		s[j] = obj[j];
	}

	for (i = 0; i < nrows; ++i) {
		for (j = rmatbeg[i]; j < rmatbeg[i+1]; ++j) {
			s[rmatind[j]] -= u[i];
		}
	}

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	free(cmatval);

	/*char *senses;
	int error, i, j, nzcnt, status;
	double *cmatval;

	nzcnt = rmatbeg[nrows];
	cmatval = (double *) malloc(nzcnt * sizeof(double));
	for (i = 0; i < nzcnt; ++i) {
		cmatval[i] = 1.0;
	}

	senses = (char *) malloc(nrows * sizeof(char));
	for (i = 0; i < nrows; ++i) {
		senses[i] = 'G';
	}

	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "dual");

	CPXaddrows(env, lp, (int) ncols, (int) nrows, nzcnt, cmatval, senses, rmatbeg, rmatind, cmatval, NULL, NULL);

	for (j = 0; j < ncols; ++j) {
		CPXchgobj(env, lp, 1, &j, &obj[j]);
	}

	// Maximize
	CPXchgprobtype(env, lp, CPXPROB_LP);

	status = CPXlpopt(env, lp);
	if (status) { printf("SCdual0_sparse - CPXmipopt error: %d\n", status); }

	// Dual vector
	status = CPXgetpi(env, lp, u, 0, (int) (nrows - 1));
	if (status) { printf("SCdual0_sparse - CPXgetx error: %d\n", status); }

	status = CPXgetdj(env, lp, s, 0, (int) (ncols - 1));
	if (status) { printf("SCdual0_sparse - CPXgetdj error: %d\n", status); }


	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	free(cmatval);
	free(senses);*/

	return 0;
}


int SCbalasheurdual1_sparse(const int *rmatbeg, const int *rmatind, const int nrows, const int *xind,
		const int xindlen, double *u, double *s) {

	int cnt, firsttime, itercnt, i, j, k, row, srtrowslen = 0;
	SCi2tuple *srtrows;

	for (i = 0; i < nrows; ++i) {
		u[i] = 0.0;
	}

	srtrows = (SCi2tuple *) malloc(nrows * sizeof(SCi2tuple));
	for (i = 0; i < nrows; ++i) {

		// No need to generate R and T(x) sets: cnt says how
		// much a row is covered and if cnt < 1 (minimum
		// cover) skip the row
		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i+1] && k < xindlen) {
			if (rmatind[j] == xind[k]) {
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k]) {
				j++;
				continue;
			}

			if (rmatind[j] > xind[k]) {
				k++;
				continue;
			}
		}

		if (cnt == 1) {
			srtrows[srtrowslen].a = rmatbeg[i+1] - rmatbeg[i];
			srtrows[srtrowslen++].b = i;
		}
	}

	qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);

	itercnt = 0;
	firsttime = 1;

	while (itercnt < srtrowslen) {

		row = srtrows[itercnt++].b;

		u[row] = INFINITY;
		for (j = rmatbeg[row]; j < rmatbeg[row + 1]; ++j) {
			u[row] = (s[rmatind[j]] < u[row]) ? s[rmatind[j]] : u[row];
		}

		for (j = rmatbeg[row]; j < rmatbeg[row + 1]; ++j) {
			s[rmatind[j]] -= u[row];
		}

		if ((itercnt == srtrowslen) && firsttime) {

			itercnt = 0;
			firsttime = 0;
			srtrowslen = 0;

			for (i = 0; i < nrows; ++i) {

				// No need to generate R and T(x) sets: cnt says how
				// much a row is covered and if cnt < 1 (minimum
				// cover) skip the row
				cnt = 0;
				j = rmatbeg[i];
				k = 0;
				while (j < rmatbeg[i+1] && k < xindlen) {
					if (rmatind[j] == xind[k]) {
						cnt++;
						j++;
						k++;
						continue;
					}

					if (rmatind[j] < xind[k]) {
						j++;
						continue;
					}

					if (rmatind[j] > xind[k]) {
						k++;
						continue;
					}
				}

				if (cnt > 1) {
					srtrows[srtrowslen].a = rmatbeg[i+1] - rmatbeg[i];
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
		const int *xind, int xindlen, double *u, double *s, const double zu) {

	int cnt, itercnt, i, j, k, row, srtrowslen = 0;
	double sdotx, usum, val;
	SCi2tuple *srtrows;

	sdotx = 0.0;
	for (j = 0; j < xindlen; ++j) {
		sdotx += s[xind[j]];
	}

	usum = 0.0;
	for (i = 0; i < nrows; ++i) {
		usum += u[i];
	}

	if (sdotx >= (zu - usum)) {
		goto TERMINATE;
	}

	srtrows = (SCi2tuple *) malloc(nrows * sizeof(SCi2tuple));
	for (i = 0; i < nrows; ++i) {

		// No need to generate R and T(x) sets: cnt says how
		// much a row is covered and if cnt < 1 (minimum
		// cover) skip the row
		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i+1] && k < xindlen) {
			if (rmatind[j] == xind[k]) {
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k]) {
				j++;
				continue;
			}

			if (rmatind[j] > xind[k]) {
				k++;
				continue;
			}
		}

		if ((u[i] > SC_EPSILON_SMALL) && (cnt > 1)) {
			srtrows[srtrowslen].a = cnt;
			srtrows[srtrowslen++].b = i;
		}
	}

	if (srtrowslen == 0) {
		free(srtrows);
		return 0;
	}

	qsort(srtrows, srtrowslen, sizeof(SCi2tuple), SCi2tuple_cmpa);

	itercnt = 0;
	while (sdotx < (zu - usum - SC_EPSILON_SMALL)) {

		if (itercnt == srtrowslen) {
			free(srtrows);
			return -1;
		}

		row = srtrows[itercnt++].b;

		val = u[row];
		j = rmatbeg[row];
		k = 0;
		while (j < rmatbeg[row+1] && k < xindlen) {
			if (rmatind[j] == xind[k]) {
				s[rmatind[j]] += val;
				sdotx += val;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k]) {
				s[rmatind[j]] += val;
				j++;
				continue;
			}

			if (rmatind[j] > xind[k]) {
				k++;
				continue;
			}
		}

		while (j < rmatbeg[row+1]) {
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
		const double zu, double y, int *rqbeg, int *rqind, int max_branch, int max_singl) {

	int cnt, qnzcnt, qnsingl, i, it, j, jt, k, h, p, sindlen, txindlen, jindlen, mjindlen;
	int *sind, *txind, *jind, *mjind, *tmp;
	double min, v, v1, v2;

	jind = (int *) malloc(xindlen * sizeof(int));
	mjind = (int *) malloc(nrows * sizeof(int));
	tmp = (int *) malloc((nrows > xindlen ? nrows : xindlen) * sizeof(int));

	// S set defined with column indices
	sindlen = 0;
	sind = (int *) malloc(xindlen * sizeof(int));
	for (j = 0; j < xindlen; ++j) {
		if (dj[xind[j]] > SC_EPSILON_SMALL) {
			sind[sindlen] = xind[j];
			sindlen++;
		}
	}

	// Tx set defined with row indices
	txindlen = 0;
	txind = (int *) malloc(nrows * sizeof(int));
	for (i = 0; i < nrows; ++i) {

		cnt = 0;
		j = rmatbeg[i];
		k = 0;
		while (j < rmatbeg[i+1] && k < xindlen) {
			if (cnt > 1) {
				break;
			}

			if (rmatind[j] == xind[k]) {
				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < xind[k]) {
				j++;
				continue;
			}

			if (rmatind[j] > xind[k]) {
				k++;
				continue;
			}
		}

		if (cnt == 1) {
			txind[txindlen] = i;
			txindlen++;
		}
	}

	p = 0;
	qnzcnt = 0;
	qnsingl = 0;
	while (1) {

		v1 = -INFINITY;
		v2 = INFINITY;

		for (j = 0; j < sindlen; ++j) {
			v1 = (dj[sind[j]] > v1) ? dj[sind[j]] : v1;
			v2 = (dj[sind[j]] >= zu - y) && (dj[sind[j]] < v2) ? dj[sind[j]] : v2;
		}

		v = v1 > v2 ? v2 : v1;

		// fill J
		jindlen = 0;
		for (j = 0; j < sindlen; ++j) {
			if ((v - dj[sind[j]]) < SC_EPSILON_SMALL) {
				jind[jindlen] = sind[j];
				jindlen++;
			}
		}

		// fill M_J
		mjindlen = 0;
		for (i = 0; i < nrows; ++i) {
			tmp[i] = 0;
		}

		for (j = 0; j < jindlen; ++j) {
			for (i = cmatbeg[jind[j]]; i < cmatbeg[jind[j] + 1]; ++i) {
				tmp[cmatind[i]] = 1;
			}
		}

		for (i = 0; i < nrows; ++i) {
			if(tmp[i]) {
				mjind[mjindlen] = i;
				mjindlen++;
			}
		}

		//printf("s=["); for (i = 0; i < sindlen; ++i) printf("%d ", sind[i]); printf("]\n");
		//printf("slen=%d, txlen=%d, jlen=%d, mjlen=%d\n", sindlen, txindlen, jindlen, mjindlen);

		// fill Q_p
		rqbeg[p] = qnzcnt;
		for (j = 0; j < ncols; ++j) {
			rqind[qnzcnt] = j * (dj[j] >= v);
			qnzcnt += (dj[j] >= v);
		}
		rqbeg[p+1] = qnzcnt;

		it = -1;
		min = INFINITY;
		i = 0;
		h = 0;
		while (i < txindlen && h < mjindlen) {
			if (txind[i] == mjind[h]) {

				cnt = 0;
				j = rmatbeg[txind[i]];
				k = rqbeg[p];
				while (j < rmatbeg[txind[i] + 1] && k < rqbeg[p+1]) {
					if (rmatind[j] == rqind[k]) {
						j++;
						k++;
						continue;
					}

					if (rmatind[j] < rqind[k]) {
						cnt++;
						j++;
						continue;
					}

					if (rmatind[j] > rqind[k]) {
						k++;
						continue;
					}
				}

				while (j < rmatbeg[txind[i] + 1]) {
					cnt++;
					j++;
				}

				if (cnt < min) {
					min = cnt;
					it = txind[i];
				}

				i++;
				h++;
				continue;
			}

			if (txind[i] < mjind[h]) {
				i++;
				continue;
			}

			if (txind[i] > mjind[h]) {
				h++;
				continue;
			}
		}

		jt = -1;
		j = 0;
		k = rmatbeg[it];
		while (j < jindlen && k < rmatbeg[it+1]) {
			if (jind[j] == rmatind[k]) {
				jt = jind[j];
				break;
			}

			if (jind[j] < rmatind[k]) {
				j++;
				continue;
			}

			if (jind[j] > rmatind[k]) {
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
		while (j < rmatbeg[it+1] && k < rqbeg[p+1]) {
			if (rmatind[j] == rqind[k]) {
				tmp[cnt] = rmatind[j];

				cnt++;
				j++;
				k++;
				continue;
			}

			if (rmatind[j] < rqind[k]) {
				j++;
				continue;
			}

			if (rmatind[j] > rqind[k]) {
				k++;
				continue;
			}
		}

		for (j = 0; j < cnt; ++j) {
			rqind[rqbeg[p] + j] = tmp[j];
		}
		qnzcnt = qnzcnt - (rqbeg[p+1] - rqbeg[p]) + cnt;
		rqbeg[p+1] = rqbeg[p] + cnt;

		qnsingl += (cnt == 1);

		// check
		if ((y + SC_EPSILON_SMALL < zu) && (p+1 >= max_branch)) {
			p = -1;
			//printf("max branch reached\n");
			break;
		}
		if ((y + SC_EPSILON_SMALL < zu) && (qnsingl >= max_singl)) {
			p = -1;
			//printf("max singl num reached\n");
			break;
		}
		if (y + SC_EPSILON_SMALL >= zu) {
			p++;
			//printf("OK!\n");
			break;
		}

		for (j = 0; j < sindlen; ++j) {
			if (sind[j] == jt) {
				sindlen--;
				memcpy(&sind[j], &sind[j+1], (sindlen - j) * sizeof(int));
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

	if (p < 0) {
		return p;
	}
	return (qnzcnt > p * log2(p)) ? p : -1;
}