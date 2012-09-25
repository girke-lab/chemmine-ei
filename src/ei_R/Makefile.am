PACKAGEFOLDER = eiR
SRCDIR = $(PACKAGEFOLDER)/src
#CXXFLAGS +=  -fPIC -DNO_MAIN -I../../../lshkit/include
#LIBS = -lm -lf2c   -L../../../lshkit/src -Wl,--whole-archive  -llshkit  -Wl,--no-whole-archive -lgsl -lgslcblas -O1
 #$(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB)  -lm -ldb  -lgsl -lgslcblas
MODULENAME = eiR
CCSOURCES  = \
	../atompair/molecule.cc \
	../atompair/formats.cc \
	../atompair/desc.cc \
	../atompair/similarity.cc \
	../atompair/simpledb.cc \
	../atompair/search.cc \
	../atompair/profiling.cc \
	../atompair/db_build.h \
	../atompair/debug.h \
	../atompair/desc.h \
	../atompair/molecule.h \
	../atompair/profiling.h \
	../atompair/script.h \
	../atompair/search.h \
	../atompair/simpledb.h \
	../atompair/db_build.cc \
	../atompair/db2db_distance.cc \
	../atompair/db2db_distance_r.cc \
	../atompair/helpers.h \
	../atompair/helpers.cc \
	../atompair/db_subset.cc \
	../matrix/format.cc \
	../embedder/embedCoords.cc \
	../embedder/solver.cc \
	../embedder/solver.h \
	../embedder/routines.c \
	../embedder/config.h \
	../lshkit/tools/search_r.cpp \
	../euclidean-search/scan.cc \
	../euclidean-search/eucsearch.h \
	../euclidean-search/eucsearch.cc \
	../euclidean-search/recordfile.cc \
	../lshkit/include \
	../lshkit/src/liblshkit.a

EXTRA_DIST = $(MODULENAME)/DESCRIPTION $(MODULENAME)/NAMESPACE $(MODULENAME)/README  $(MODULENAME)/src/Makevars

LOCAL_INST=~/testing

#empty=
#space=$(empty) $(empty)
#ESC_CXXFLAGS=$(subst $(space),\$(space),$(CXXFLAGS))
#ESC_LIBS=$(subst $(space),\$(space),$(LIBS))

if BUILD_EIR
all : $(SRCDIR)/r_wrap.cpp ccsource

$(SRCDIR)/r_wrap.cpp : swig.i
	mkdir -p $(MODULENAME)/R $(MODULENAME)/src
	$(SWIG) -c++  -module $(MODULENAME) -o $(SRCDIR)/r_wrap.cpp -r swig.i
	mv $(SRCDIR)/$(MODULENAME).R $(PACKAGEFOLDER)/R/

ccsource : $(CCSOURCES)
	cp -a $(CCSOURCES) $(SRCDIR)

#MAKEFLAGS='CXXFLAGS=$(ESC_CXXFLAGS) LDFLAGS=$(ESC_LIBS)' 
check: all
	$(R_PATH) CMD check $(PACKAGEFOLDER)
test: all local
	export RCMDCHECK=FALSE; \
	export R_LIBS_USER=$(LOCAL_INST); \
	$(R_PATH) --vanilla --slave < eiR/tests/test.R

local: all
	 $(R_PATH) CMD INSTALL -l $(LOCAL_INST) $(PACKAGEFOLDER)

install : all
	 $(R_PATH) CMD INSTALL $(PACKAGEFOLDER)
endif

dist : all
	rm -f  $(PACKAGEFOLDER)/src/*.o $(PACKAGEFOLDER)/src/*.so
	tar czf $(PACKAGEFOLDER).tar.gz --exclude=*.swp $(PACKAGEFOLDER)


clean-generic :
	rm -f $(SRCDIR)/r_wrap.cpp
	rm -f $(PACKAGEFOLDER)/R/$(MODULENAME).R
	(cd $(SRCDIR);\
	rm -f *.so *.o;\
	rm -rf include;\
	for i in $(CCSOURCES); do \
		rm -f `basename $$i`;\
	done)
