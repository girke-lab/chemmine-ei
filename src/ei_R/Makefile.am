PACKAGEFOLDER = eiR
SRCDIR = $(PACKAGEFOLDER)/src
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

EXTRA_DIST = $(MODULENAME)/DESCRIPTION \
				 $(MODULENAME)/NAMESPACE \
				 $(MODULENAME)/README  \
				 $(MODULENAME)/src/Makevars \
				 $(MODULENAME)/man/eiInit.Rd \
				 $(MODULENAME)/man/eiMakeDb.Rd \
				 $(MODULENAME)/man/eiPerformanceTest.Rd \
				 $(MODULENAME)/man/eiQuery.Rd \
				 $(MODULENAME)/R/core.R  \
				 $(MODULENAME)/tests/test.R  \
				 $(MODULENAME)/inst/runit/runit.eiInit.r  \
				 $(MODULENAME)/data/example_compounds.RData  

LOCAL_INST=~/testing


if BUILD_EIR
all : $(SRCDIR)/r_wrap.cpp ccsource

$(SRCDIR)/r_wrap.cpp : swig.i
	mkdir -p $(MODULENAME)/R $(MODULENAME)/src
	$(SWIG) -c++  -module $(MODULENAME) -o $(SRCDIR)/r_wrap.cpp -r swig.i
	mv $(SRCDIR)/$(MODULENAME).R $(PACKAGEFOLDER)/R/

ccsource : $(CCSOURCES)
	cp -a $(CCSOURCES) $(SRCDIR)

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

dist : all
	 $(R_PATH) CMD build $(PACKAGEFOLDER)

#rm -f  $(PACKAGEFOLDER)/src/*.o $(PACKAGEFOLDER)/src/*.so
#tar czf $(PACKAGEFOLDER).tar.gz --exclude=*.swp $(PACKAGEFOLDER)

endif

clean-generic :
	rm -f $(SRCDIR)/r_wrap.cpp
	rm -f $(PACKAGEFOLDER)/R/$(MODULENAME).R
	(cd $(SRCDIR);\
	rm -f *.so *.o;\
	rm -rf include;\
	for i in $(CCSOURCES); do \
		rm -f `basename $$i`;\
	done)
