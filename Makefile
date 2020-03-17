
###################          CHANGE THESE PATHS            ####################

SHELL=/bin/sh

BOOST_HOME=/usr/local/boost_1_71_0
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
PROJ_HOME=/home/macs/coding/optimization/set-covering-problem

###############################################################################

CC=g++
CFLAGS=-g -Wall# -O2 -fpermissive

LDLIBS=-L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
LDLIBSCPX=-L$(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic -lcplex $(LDLIBS) 
INC=-I$(BOOST_HOME) -I$(PROJ_HOME)/src
INCCPX=-I$(CPLEX_HOME)/cplex/include/ $(INC)

####    COLORS AND FORMATS
BOLD=\e[1m
BOLDEND=\e[21m
COLEND=\e[39m
COLGREEN=\e[32m
COLRED=\e[31m
COLYELLOW=\e[33m


########################          CPXSOL         ##############################

cpxsol : balas_dense balas_sparse cpx_callbacks cpx_common cpx_main cpx_solver
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/cpx_callbacks.o lib/cpx_common.o lib/cpx_main.o lib/cpx_solver.o -o $@ $(LDLIBSCPX)
	@mv $@ lib
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule cpxsol$(BOLDEND)"

cpx_callbacks : src/cpx_callbacks.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	@mv $@.o lib

cpx_common : src/cpx_common.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	@mv $@.o lib

cpx_main : src/cpx_main.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	@mv $@.o lib

cpx_solver : src/cpx_solver.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INCCPX)
	@mv $@.o lib


########################          BALAS LIB         ###########################

balas_dense : src/balas_dense.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	@mv $@.o lib

balas_sparse : src/balas_sparse.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	@mv $@.o lib


########################          BALSOL         ##############################

balsol : balas_dense balas_sparse balas_common balas_main balas_solver
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/balas_common.o lib/balas_main.o lib/balas_solver.o -o $@ $(LDLIBS)
	@mv $@ lib
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule balsol$(BOLDEND)"

balas_common : src/balas_common.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	@mv $@.o lib

balas_main : src/balas_main.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	@mv $@.o lib

balas_solver : src/balas_solver.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	@mv $@.o lib


#########################         UNIT TESTS       ############################

test_balas_dense : balas_dense
	$(CC) $(CFLAGS) -c tests/$@.cpp $(INC)
	$(CC) $(CFLAGS) $@.o lib/balas_dense.o $(LDLIBS) -o $@
	@rm $@.o
	@mv $@ tests
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule test_balas_dense$(BOLDEND)"

test_balas_sparse : balas_sparse
	$(CC) $(CFLAGS) -c tests/$@.cpp $(INC)
	$(CC) $(CFLAGS) $@.o lib/balas_sparse.o $(LDLIBS) -o $@
	@rm $@.o
	@mv $@ tests
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule test_balas_sparse$(BOLDEND)"

unittests-all: test_balas_dense test_balas_sparse
	@cd tests && ./test_balas_dense
	@cd tests && ./test_balas_sparse

unittests-balas_dense: test_balas_dense
	@cd tests && ./test_balas_dense

unittests-balas_sparse: test_balas_sparse
	@cd tests && ./test_balas_sparse

clean:
	@rm -f lib/balsol lib/cpxsol lib/*.o
	@cd tests && rm test_balas_dense test_balas_sparse

.PHONY: clean