# Set Covering Problem
### A computational study of different approaches with CPLEX

The present work deals with the implementation of several methods for the preprocessing and the exact resolution of the **Set Covering Problem** (SCP) with techniques of **Mixed-Integer Linear Programming** (MILP).

The main ideas used in this project were provided by the works of E. Balas and A. Ho ([1](https://link.springer.com/chapter/10.1007/BFb0120886), [2](https://link.springer.com/chapter/10.1007/BFb0120885)).
The implemented approaches have been tested on the `scpnre1-scpnrf5` instances available on the [OR-Library](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html).
The results obtained are described in this [report](https://github.com/caerbannogwhite/set-covering-problem/tree/master/report/report.pdf).

The software requires the following software:
- [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)
- [Boost](https://www.boost.org/)
- [Armadillo](http://arma.sourceforge.net/)
- [OpenBLAS](http://www.openblas.net/)
- [ARPACK](https://www.caam.rice.edu/software/ARPACK/)
- [LAPACK](https://github.com/Reference-LAPACK/lapack)
- [SuperLU](https://github.com/xiaoyeli/superlu)

On Ubuntu or Debian, you should run

sudo apt install cmake libopenblas-dev liblapack-dev libarpack2-dev libsuperlu-dev

Once you have installed them in your computer, you should change the `CPLEX_HOME` and the `BOOST_HOME` variables in `Makefile`.

### Extensions

Future extensions of the project include performance improvements and the implementation of a solver (both heuristic and with optimality test) independent of CPLEX.