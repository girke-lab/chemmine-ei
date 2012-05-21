
SRC = solver.cc routines.c ../atompair/profiling.cc

#AM_CFLAGS = -g -D_TIMER_
#AM_CXXFLAGS = -g -D_TIMER_
AM_CFLAGS = -mtune=nocona  -O3  -fomit-frame-pointer -funroll-loops -mfpmath=sse  -D_TIMER_ -I ../atompair/
AM_CXXFLAGS = -mtune=nocona -O3  -fomit-frame-pointer -funroll-loops -mfpmath=sse  -D_TIMER_ -I ../atompair/
#AM_LDFLAGS = -static
LIBS = -lm -lf2c  

bin_PROGRAMS = coord icoord coord_server

coord_SOURCES = coord.cc $(SRC)
icoord_SOURCES = icoord.cc $(SRC)
icoord_LDADD = -lreadline
coord_server_SOURCES = coord_server.cc $(SRC)
