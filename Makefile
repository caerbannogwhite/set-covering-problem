
###################          CHANGE THESE PATHS            ####################

SHELL=/bin/sh
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

###############################################################################

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive -DIL_STD

LDLIBS=-L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
LDLIBSCPX=-L$(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_HOME)/concert/lib/x86-64_linux/static_pic -lcplex $(LDLIBS) 
INC=-I$(BOOST_HOME)
INCCPX=-I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ $(INC)


########################          CPXSOL         ##############################

cpxsol : balas_dense balas_sparse cpx_callbacks cpx_common cpx_main cpx_solver
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/cpx_callbacks.o lib/cpx_common.o lib/cpx_main.o lib/cpx_solver.o -o $@ $(LDLIBSCPX)
	mv $@ lib
	rule cpxsol COMPLETE

cpx_callbacks : src/cpx_callbacks.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	mv $@.o lib

cpx_common : src/cpx_common.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	mv $@.o lib

cpx_main : src/cpx_main.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	mv $@.o lib

cpx_solver : src/cpx_solver.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	mv $@.o lib

########################          END CPXSOL         ##########################

balas_dense : src/balas_dense.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

balas_sparse : src/balas_sparse.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib


########################          BALSOL         ##############################

balsol : balas_dense balas_sparse balas_common balas_main balas_solver
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/balas_common.o lib/balas_main.o lib/balas_solver.o -o $@ $(LDLIBS)
	mv $@ lib
	rule balsol COMPLETE

balas_common : src/balas_common.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

balas_main : src/balas_main.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

balas_solver : src/balas_solver.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

clean:
	rm -f lib/balsol lib/cpxsol lib/*.o

.PHONY: clean