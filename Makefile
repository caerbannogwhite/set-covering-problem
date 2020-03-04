P=scsolver
SHELL=/bin/sh

CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio128
BOOST_HOME=/usr/local/boost_1_71_0

CC=g++
CFLAGS=-g -Wall -O2 -fpermissive

LDLIBS=-L$(CPLEX_HOME)/lib/x86-64_linux/static_pic -lcplex -L$(BOOST_HOME)/stage/lib -lboost_program_options -lm -lpthread -ldl
INC=-I$(CPLEX_HOME)/cplex/include/ -I$(CPLEX_HOME)/concert/include/ -I$(BOOST_HOME)

$(P) : balas_dense.o balas_sparse.o callbacks.o common.o main.o preprocessing.o sc.o
	$(CC) $(CFLAGS) bin/balas_dense.o bin/balas_sparse.o bin/callbacks.o bin/common.o bin/main.o bin/preprocessing.o bin/sc.o -o $@ $(LDLIBS)
	mv $@ bin

balas_dense.o : src/balas_dense.cpp src/balas_dense.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

balas_sparse.o : src/balas_sparse.cpp src/balas_sparse.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

callbacks.o : src/callbacks.cpp src/callbacks.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

common.o : src/common.cpp src/common.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

main.o : src/main.cpp src/main.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

preprocessing.o : src/preprocessing.cpp src/preprocessing.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

sc.o : src/sc.cpp src/sc.hpp
	$(CC) $(CFLAGS) -c src/$*.cpp $(INC)
	mv $@ bin

clean:
	rm -f bin/$(P) bin/*.o

.PHONY: clean