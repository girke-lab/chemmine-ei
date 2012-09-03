include ../Makefile.inc


INCLUDEDIR = ../include 
AM_CXXFLAGS+= -I../include -fPIC

#SRC = apost.cpp char_bit_cnt.cpp  mplsh.cpp   mplsh-model.cpp vq.cpp
SRC = apost.cpp char_bit_cnt.cpp  mplsh.cpp   mplsh-model.cpp


if BUILD_LSHKIT
#noinst_HEADERS = kdtree.h ../include/lshkit.h ../include/lshkit/apost.h ../include/lshkit/archive.h ../include/lshkit/common.h ../include/lshkit/composite.h ../include/lshkit/concept.h ../include/lshkit/eval.h ../include/lshkit/forest.h ../include/lshkit/histogram.h ../include/lshkit/kernel.h ../include/lshkit/lsh-index.h ../include/lshkit/lsh.h ../include/lshkit/matrix.h ../include/lshkit/matrix-io.h ../include/lshkit/metric.h ../include/lshkit/mplsh-model.h ../include/lshkit/mplsh.h   ../include/lshkit/sketch.h ../include/lshkit/spectral-hash.h ../include/lshkit/topk.h ../include/lshkit/tune.h ../include/lshkit/vq.h ../include/lshkit/adaptor.h ../include/lshkit/lsh-file.h ../include/lshkit/mplsh-file.h ../include/lshkit/storage.h
noinst_HEADERS = ../include/lshkit.h ../include/lshkit/apost.h ../include/lshkit/archive.h ../include/lshkit/common.h ../include/lshkit/composite.h ../include/lshkit/concept.h ../include/lshkit/eval.h ../include/lshkit/forest.h ../include/lshkit/histogram.h ../include/lshkit/kernel.h ../include/lshkit/lsh-index.h ../include/lshkit/lsh.h ../include/lshkit/matrix.h ../include/lshkit/matrix-io.h ../include/lshkit/metric.h ../include/lshkit/mplsh-model.h ../include/lshkit/mplsh.h   ../include/lshkit/sketch.h ../include/lshkit/spectral-hash.h ../include/lshkit/topk.h ../include/lshkit/tune.h ../include/lshkit/vq.h ../include/lshkit/adaptor.h ../include/lshkit/lsh-file.h ../include/lshkit/mplsh-file.h ../include/lshkit/storage.h
noinst_LIBRARIES =  liblshkit.a
liblshkit_a_SOURCES = $(SRC)
endif

