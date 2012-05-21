
SRC = molecule.cc formats.cc desc.cc similarity.cc simpledb.cc search.cc profiling.cc 

AM_CXXFLAGS = -g -Wall -fPIC -D_BUILD_DB_BUILD_ -DPRELOAD_DATABASE -O3 -mtune=nocona #-DDEBUG
LIBS =    -lreadline -static
bin_PROGRAMS = compare db_build db_search db_isearch db_search_self db_search_range db2db_distance db_subset

noinst_LIBRARIES = libdescriptors.a
libdescriptors_a_SOURCES =  $(SRC)

compare_SOURCES = compare.cc $(SRC)
db_build_SOURCES = db_build.cc $(SRC)
db_search_SOURCES = db_search.cc  $(SRC)
db_isearch_SOURCES = db_isearch.cc $(SRC)
db_search_self_SOURCES = db_search_self.cc  $(SRC)
db_search_range_SOURCES = db_search_range.cc $(SRC)
db2db_distance_SOURCES = db2db_distance.cc  $(SRC)
db_subset_SOURCES = db_subset.cc   $(SRC)

