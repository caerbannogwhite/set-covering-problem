P=test
SHELL=/bin/sh

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive

LDLIBS=-L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
INC=-I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)

$(P) : balas_dense.o main.o
	$(CC) $(CFLAGS) balas_dense.o main.o -o $@ $(LDLIBS)

balas_dense.o : balas_dense.cpp balas_dense.hpp
	$(CC) $(CFLAGS) -c $*.cpp $(INC)

main.o : main.cpp
	$(CC) $(CFLAGS) -c $*.cpp $(INC)

clean:
	rm -f $(P) *.o

.PHONY: clean