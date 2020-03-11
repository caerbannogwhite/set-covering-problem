
SHELL=/bin/sh

PROJ_HOME=/home/macs/coding/optimization/set-covering-problem
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive

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

test_balas_dense : test_balas_dense.cpp $(PROJ_HOME)/lib/balas_dense.o
	$(CC) $(CFLAGS) -c $@.cpp $(INC)
	$(CC) $(CFLAGS) $@.o $(PROJ_HOME)/lib/balas_dense.o $(LDLIBS) -o $@
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule test_balas_dense$(BOLDEND)"

test_balas_sparse : test_balas_sparse.cpp $(PROJ_HOME)/lib/balas_sparse.o
	$(CC) $(CFLAGS) -c $@.cpp $(INC)
	$(CC) $(CFLAGS) $@.o $(PROJ_HOME)/lib/balas_sparse.o $(LDLIBS) -o $@
	@echo "$(BOLD)[$(COLGREEN)COMPLETE$(COLEND)] Rule test_balas_sparse$(BOLDEND)"

clean:
	rm -f *.o test_balas_dense test_balas_sparse

.PHONY: clean