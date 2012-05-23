include ../Makefile.inc


INCLUDEDIR = ../include
AM_CXXFLAGS+= -I$(INCLUDEDIR)

SRC = btune.cpp  cache.cpp  file-glibc.cpp  mplsh.cpp  mplsh-file.cpp  mplsh-model.cpp


if BUILD_LSHKIT
noinst_HEADERS = $(INCLUDEDIR)/btune.h $(INCLUDEDIR)/lshkit.h $(INCLUDEDIR)/lshkit/adaptor.h $(INCLUDEDIR)/lshkit/common.h $(INCLUDEDIR)/lshkit/composite.h $(INCLUDEDIR)/lshkit/concept.h $(INCLUDEDIR)/lshkit/eval.h $(INCLUDEDIR)/lshkit/flat-index.h $(INCLUDEDIR)/lshkit/forest-index.h $(INCLUDEDIR)/lshkit/histogram.h $(INCLUDEDIR)/lshkit/lsh-file.h $(INCLUDEDIR)/lshkit/lsh.h $(INCLUDEDIR)/lshkit/matrix.h $(INCLUDEDIR)/lshkit/matrix-io.h $(INCLUDEDIR)/lshkit/metric.h $(INCLUDEDIR)/lshkit/mplsh-file.h $(INCLUDEDIR)/lshkit/mplsh.h $(INCLUDEDIR)/lshkit/mplsh-model.h $(INCLUDEDIR)/lshkit/param-fitting.h $(INCLUDEDIR)/lshkit/sketch.h $(INCLUDEDIR)/lshkit/storage.h $(INCLUDEDIR)/lshkit/topk.h
noinst_LIBRARIES =  liblshkit.a
liblshkit_a_SOURCES = $(SRC)
endif

