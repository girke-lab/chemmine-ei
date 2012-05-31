include ../Makefile.inc

program_transform_name = s&^&ei-&

AM_CXXFLAGS+= -I../include
AM_LDFLAGS = -L../src
LIBS += -llshkit 
SRC = flat-index.cpp   mkmplsh.cpp  runmplsh.cpp  test-probe.cpp


if BUILD_LSHKIT
bin_PROGRAMS =  runmplsh  test-probe 
#mkmplsh_SOURCES=mkmplsh.cpp
runmplsh_SOURCES=runmplsh.cpp
test_probe_SOURCES=test-probe.cpp
endif
