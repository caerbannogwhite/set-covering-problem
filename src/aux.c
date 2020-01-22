
#include "aux.h"

int SCint_cmp(const void *p, const void *q) {
	return *((int *) p) - *((int *) q);
}


int SCi2tuple_cmpa(const void *p, const void *q) {
	return ((SCi2tuple *) p)->a - ((SCi2tuple *) q)->a;
}


int SCi2tuple_cmparev(const void *p, const void *q) {
	return ((SCi2tuple *) q)->a - ((SCi2tuple *) p)->a;
}

int SCi3tuple_cmpa(const void *p, const void *q) {
	return ((SCi3tuple *) p)->a - ((SCi3tuple *) q)->a;
}

int SCi3tuple_cmparev(const void *p, const void *q) {
	return ((SCi3tuple *) q)->a - ((SCi3tuple *) p)->a;
}

int SCidtuple_cmpb(const void *p, const void *q) {
	double cmp = ((SCidtuple *) p)->b - ((SCidtuple *) q)->b;
	if (cmp > SC_EPSILON_SMALL) {
		return 1;
	} else if (cmp < -SC_EPSILON_SMALL) {
		return -1;
	} else {
		return 0;
	}
}

double func1(const double c, const unsigned int k) { return c; }
double func2(const double c, const unsigned int k) { return c / k; }
double func3(const double c, const unsigned int k) { return k < 2 ? c : c / log2(k); }
double func4(const double c, const unsigned int k) { return k < 2 ? c / k : c / k * log2(k); }
double func5(const double c, const unsigned int k) { return k < 3 ? c / k : c / k * log(k); }
