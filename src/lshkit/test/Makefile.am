include ../Makefile.inc

CC=$(CXX)

# change this to help me locate header files for boost and boost-threadpool
#AM_CXXFLAGS += -I../include -I ../../atompair/ -I/home/ycao/opt/src/boost -I/home/ycao/opt/src/threadpool/boost
AM_CXXFLAGS+= -I../include -I/usr/include/boost -I/home/khoran/programs/threadpool-0_2_5-src/threadpool/boost -I/home/khoran/programs/libdb/usr/include/

# change this to help me locate boost library files
#LDFLAGS += -L../src  -L ../../atompair -L/home/ycao/opt64/lib/ -L /home/khoran/programs/libdb/usr/lib/
AM_LDFLAGS += -L../src  -L ../../atompair -L /usr/lib/  # -L /home/khoran/programs/libdb/usr/lib/

LIBS += -llshkit -ldescriptors

SRC = fitdata.cpp  flat-index.cpp  isearch.cpp  mkmplsh.cpp  mplsh.cpp  mplsh-predict.cpp  mplsh-tune.cpp  runmplsh.cpp  scan.cpp  search.cpp  searchp.cpp  test-probe.cpp


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
