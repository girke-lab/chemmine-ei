include ../Makefile.inc

# change this to help me locate header files for boost and boost-threadpool
#AM_CXXFLAGS+= -I../include -I/home/ycao/opt/src/boost -I/home/ycao/opt/src/threadpool/boost
AM_CXXFLAGS+= -I../include -I/usr/include/boost -I/home/khoran/programs/threadpool-0_2_5-src/threadpool/boost -I/home/khoran/programs/libdb/usr/include/

#LIB = liblshkit.a
SRC = btune.cpp  cache.cpp  file-glibc.cpp  mplsh.cpp  mplsh-file.cpp  mplsh-model.cpp


noinst_LIBRARIES =  liblshkit.a
liblshkit_a_SOURCES = $(SRC)

