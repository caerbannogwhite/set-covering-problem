P=scsolver
SHELL=/bin/sh

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive -DIL_STD

LDLIBS=-L$(CPLEX_HOME)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_HOME)/concert/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl -larmadillo
INC=-I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)

$(P) : balas_dense.o balas_sparse.o callbacks.o common.o main.o preprocessing.o sc.o
	$(CC) $(CFLAGS) lib/balas_dense.o lib/balas_sparse.o lib/callbacks.o lib/common.o lib/main.o lib/preprocessing.o lib/sc.o -o $@ $(LDLIBS)
	mv $@ lib

balas_dense.o : src/balas_dense.cpp src/balas_dense.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

balas_sparse.o : src/balas_sparse.cpp src/balas_sparse.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

callbacks.o : src/callbacks.cpp src/callbacks.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

common.o : src/common.cpp src/common.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

main.o : src/main.cpp src/main.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

preprocessing.o : src/preprocessing.cpp src/preprocessing.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

sc.o : src/sc.cpp src/sc.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ lib

clean:
	rm -f lib/$(P) lib/*.o

.PHONY: clean