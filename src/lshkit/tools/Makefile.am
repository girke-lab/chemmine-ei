include ../Makefile.inc

program_transform_name = s&^&ei-&

AM_CXXFLAGS+= -I../include -fPIC
AM_LDFLAGS = -L../src  -L../../atompair 
LIBS += -llshkit -ldescriptors 
#SRC = fitdata.cpp isearch.cpp mplsh.cpp mplsh-predict.cpp mplsh-tune.cpp scan.cpp search.cpp searchp.cpp search_server.cpp


if BUILD_LSHKIT
bin_PROGRAMS = fitdata mplsh-predict mplsh-tune scan search_server single_search
fitdata_SOURCES= fitdata.cpp
#isearch_SOURCES=isearch.cpp
mplsh_predict_SOURCES= mplsh-predict.cpp
mplsh_tune_SOURCES=mplsh-tune.cpp
scan_SOURCES=scan.cpp
#search_SOURCES=search.cpp
#searchp_SOURCES=searchp.cpp
search_server_SOURCES=search_server.cpp
#txt2bin_SOURCES=txt2bin.cpp
#mplsh_run_SOURCES=mplsh-run.cpp
single_search_SOURCES=single_search.cpp
endif
