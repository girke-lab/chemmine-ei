
program_transform_name = s&^&ei-&
SUBDIRS = python R
SRC = molecule.cc formats.cc desc.cc similarity.cc simpledb.cc search.cc profiling.cc db_build.h debug.h desc.h molecule.h profiling.h script.h search.h simpledb.h helpers.cc

AM_CXXFLAGS = -g -Wall -fPIC -D_BUILD_DB_BUILD_ -DPRELOAD_DATABASE -O3 -mtune=nocona #-DDEBUG
LIBS =    -lreadline -static

if HAVE_OPENBABEL
AM_CXXFLAGS += -DHAS_OPENBABEL $(OPENBABEL_CFLAGS)
LIBS += -lopenbabel
endif

bin_PROGRAMS = db_builder db2db_distance db_subset 

noinst_LIBRARIES = libdescriptors.a
libdescriptors_a_SOURCES =  $(SRC)

db_builder_SOURCES = db_build.cc $(SRC)
db2db_distance_SOURCES = db2db_distance.cc  $(SRC)
db_subset_SOURCES = db_subset.cc   $(SRC)

#fp_db_compare_server_SOURCES = db_compare_server.cc $(SRC)
#db_search_SOURCES = db_search.cc  $(SRC)
#db_isearch_SOURCES = db_isearch.cc $(SRC)
#compare_SOURCES = compare.cc $(SRC)
#db_search_self_SOURCES = db_search_self.cc  $(SRC)
#db_search_range_SOURCES = db_search_range.cc $(SRC)
#fp_search_server_SOURCES = fp_search_server.cc $(SRC)

