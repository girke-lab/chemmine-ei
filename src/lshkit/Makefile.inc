AR = ar
ARFLAGS = cru
PG= #-pg
OPENMP = #-fopenmp
#CC = gcc
#CXX = g++ 
AM_CXXFLAGS = $(PG) $(OPENMP)  -O3 -mtune=nocona
AM_LDFLAGS = $(PG) $(OPENMP)
#LDLIBS = -lboost_program_options-gcc41-mt -lboost_thread-gcc41-mt -lm -ldb  -lgsl -lgslcblas 
LIBS = -lboost_program_options -lboost_thread -lm -ldb-4.8  -lgsl -lgslcblas 

#.SUFFIXES:    .cpp .h .o .d .a
