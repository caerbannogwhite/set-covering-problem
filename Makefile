P=scsolver
SHELL=/bin/sh

###################          CHANGE THESE PATHS            ####################

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

###############################################################################

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive -DIL_STD

LDLIBS=-L$(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_HOME)/concert/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
INC=-I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)

$(P) : balas_dense balas_sparse callbacks common main preprocessing sc
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/callbacks.o lib/common.o lib/main.o lib/preprocessing.o lib/sc.o -o $@ $(LDLIBS)
	mv $@ lib

balas_dense : src/balas_dense.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

balas_sparse : src/balas_sparse.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

callbacks : src/callbacks.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

common : src/common.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

main : src/main.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

preprocessing : src/preprocessing.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

sc : src/sc.cpp
	$(CC) $(CFLAGS) -c src/$@.cpp $(INC)
	mv $@.o lib

clean:
	rm -f lib/$(P) lib/*.o

.PHONY: clean