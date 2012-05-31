include ../Makefile.inc

program_transform_name = s&^&ei-&

AM_CXXFLAGS+= -I../include
AM_LDFLAGS = -L../src  -L../../atompair 
LIBS += -llshkit -ldescriptors
SRC = fitdata.cpp isearch.cpp mplsh.cpp mplsh-predict.cpp mplsh-tune.cpp scan.cpp search.cpp searchp.cpp search_server.cpp


if BUILD_LSHKIT
bin_PROGRAMS = fitdata mplsh-predict mplsh-tune scan search_server
fitdata_SOURCES= fitdata.cpp
#isearch_SOURCES=isearch.cpp
mplsh_predict_SOURCES= mplsh-predict.cpp
mplsh_tune_SOURCES=mplsh-tune.cpp
scan_SOURCES=scan.cpp
#search_SOURCES=search.cpp
#searchp_SOURCES=searchp.cpp
search_server_SOURCES=search_server.cpp
endif
