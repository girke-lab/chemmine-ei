include ../Makefile.inc


AM_CXXFLAGS+= -I../include 

SRC = btune.cpp  cache.cpp  file-glibc.cpp  mplsh.cpp  mplsh-file.cpp  mplsh-model.cpp


if BUILD_LSHKIT
noinst_LIBRARIES =  liblshkit.a
liblshkit_a_SOURCES = $(SRC)
endif

