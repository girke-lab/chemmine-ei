include ../Makefile.inc


AM_CXXFLAGS+= -I../include
AM_LDFLAGS = -L../src  -L ../../atompair 
LIBS += -llshkit -ldescriptors
SRC = fitdata.cpp  flat-index.cpp  isearch.cpp  mkmplsh.cpp  mplsh.cpp  mplsh-predict.cpp  mplsh-tune.cpp  runmplsh.cpp  scan.cpp  search.cpp  searchp.cpp  test-probe.cpp


if BUILD_LSHKIT
bin_PROGRAMS = fitdata  flat-index  isearch  mkmplsh  mplsh  mplsh-predict  mplsh-tune  runmplsh  scan  search  searchp  test-probe
fitdata_SOURCES= fitdata.cpp
flat_index_SOURCES=flat-index.cpp
isearch_SOURCES=isearch.cpp
mkmplsh_SOURCES=mkmplsh.cpp
mplsh_SOURCES=mplsh.cpp
mplsh_predict_SOURCES= mplsh-predict.cpp
mplsh_tune_SOURCES=mplsh-tune.cpp
runmplsh_SOURCES=runmplsh.cpp
scan_SOURCES=scan.cpp
search_SOURCES=search.cpp
searchp_SOURCES=searchp.cpp
test_probe_SOURCES=test-probe.cpp
endif
