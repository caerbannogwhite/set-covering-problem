
SHELL=/bin/sh

PROJ_HOME=/home/macs/coding/optimization/set-covering-problem
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive -DIL_STD

LDLIBS=-L$(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_HOME)/concert/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
INC=-I$(PROJ_HOME)/src/ -I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)


test_balas_dense : test_balas_dense.cpp $(PROJ_HOME)/lib/balas_dense.o
	$(CC) $(CFLAGS) -c $@.cpp $(INC)
	$(CC) $(CFLAGS) $@.o $(PROJ_HOME)/lib/balas_dense.o $(LDLIBS) -o $@

clean:
	rm -f *.o test_balas_dense

.PHONY: clean