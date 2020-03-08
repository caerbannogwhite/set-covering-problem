
SHELL=/bin/sh

HOME=/home/macs/coding/optimization/set-covering-problem
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive

LDLIBS=-L$(CPLEX_HOME)/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
INC=-I$(HOME) -I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)


test_balas_dense : test_balas_dense.cpp
	$(CC) $(CFLAGS) -c $*.cpp $(INC)
	$(CC) $(CFLAGS) $*.o $(LDLIBS) -o $*

clean:
	rm -f *.o

.PHONY: clean