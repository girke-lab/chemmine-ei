include ../Makefile.inc


INCLUDEDIR = ../include
AM_CXXFLAGS+= -I$(INCLUDEDIR)

SRC = apost.cpp char_bit_cnt.cpp  mplsh.cpp   mplsh-model.cpp vq.cpp


if BUILD_LSHKIT
noinst_HEADERS = $(INCLUDEDIR)/lshkit.h $(INCLUDEDIR)/lshkit/apost.h $(INCLUDEDIR)/lshkit/archive.h $(INCLUDEDIR)/lshkit/common.h $(INCLUDEDIR)/lshkit/composite.h $(INCLUDEDIR)/lshkit/concept.h $(INCLUDEDIR)/lshkit/eval.h $(INCLUDEDIR)/lshkit/forest.h $(INCLUDEDIR)/lshkit/histogram.h $(INCLUDEDIR)/lshkit/kernel.h $(INCLUDEDIR)/lshkit/lsh-index.h $(INCLUDEDIR)/lshkit/lsh.h $(INCLUDEDIR)/lshkit/matrix.h $(INCLUDEDIR)/lshkit/matrix-io.h $(INCLUDEDIR)/lshkit/metric.h $(INCLUDEDIR)/lshkit/mplsh-model.h $(INCLUDEDIR)/lshkit/mplsh.h   $(INCLUDEDIR)/lshkit/sketch.h $(INCLUDEDIR)/lshkit/spectral-hash.h $(INCLUDEDIR)/lshkit/topk.h $(INCLUDEDIR)/lshkit/tune.h $(INCLUDEDIR)/lshkit/vq.h $(INCLUDEDIR)/lshkit/adaptor.h $(INCLUDEDIR)/lshkit/lsh-file.h $(INCLUDEDIR)/lshkit/mplsh-file.h $(INCLUDEDIR)/lshkit/storage.h
noinst_LIBRARIES =  liblshkit.a
liblshkit_a_SOURCES = $(SRC)
endif

