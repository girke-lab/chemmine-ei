
SRC = recordfile.cc search.cc ../atompair/desc.cc ../atompair/formats.cc ../atompair/molecule.cc ../atompair/similarity.cc ../atompair/simpledb.cc  ../atompair/profiling.cc eucsearch.h search.h
AM_CXXFLAGS = -O3 -mtune=nocona -D_REFINE_ -I ../atompair #-g
AM_LDFLAGS = -static

bin_PROGRAMS = rnn knn knn-self scan self split
rnn_SOURCES = $(SRC) rnn.cc
knn_SOURCES = $(SRC) knn.cc
knn_self_SOURCES = $(SRC) knn-self.cc
scan_SOURCES = $(SRC) scan.cc
self_SOURCES = $(SRC) self.cc
split_SOURCES = $(SRC) split.cc
