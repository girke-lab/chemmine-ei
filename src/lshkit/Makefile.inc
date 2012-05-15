AR = ar
ARFLAGS = cru
PG= #-pg
OPENMP = #-fopenmp
CC = gcc
CXX = g++ 
CXXFLAGS += $(PG) $(OPENMP)  -O3 -mtune=nocona
LDFLAGS += $(PG) $(OPENMP)
#LDLIBS += -lboost_program_options-gcc41-mt -lboost_thread-gcc41-mt -lm -ldb  -lgsl -lgslcblas 
LDLIBS += -lboost_program_options -lboost_thread -lm -ldb-4.8  -lgsl -lgslcblas 

.SUFFIXES:    .cpp .h .o .d .a
