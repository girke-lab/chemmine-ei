PACKAGEFOLDER = eiR
SRCDIR = $(PACKAGEFOLDER)/src
MODULENAME = eiR
CCSOURCES  = \
	../embedder/embedCoords.cc \
	../embedder/solver.cc \
	../embedder/solver.h \
	../embedder/routines.c \
	../embedder/config.h \
	../matrix/format.cc \
	../atompair/profiling.cc \
	../atompair/profiling.h \
	../euclidean-search/eucsearch.h \
	../euclidean-search/eucsearch.cc \
	../euclidean-search/scan.cc \
	../euclidean-search/recordfile.cc \
	../lshkit/tools/search_r.cpp \
	../lshkit/include \
	../lshkit/src/apost.cpp \
	../lshkit/src/char_bit_cnt.cpp  \
	../lshkit/src/mplsh.cpp   \
	../lshkit/src/mplsh-model.cpp

EXTRA_DIST = $(MODULENAME)/DESCRIPTION \
				 $(MODULENAME)/NAMESPACE \
				 $(MODULENAME)/README  \
				 $(MODULENAME)/NEWS\
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
	patch $(PACKAGEFOLDER)/R/$(MODULENAME).R $(MODULENAME).R.patch

ccsource : $(CCSOURCES)
	cp -a $(CCSOURCES) $(SRCDIR)

check: all
	$(R_PATH) CMD check $(PACKAGEFOLDER)
test: all local
	export RCMDCHECK=FALSE; \
	export R_LIBS_USER=$(LOCAL_INST); \
	$(R_PATH) --vanilla --slave < eiR/tests/runTests.R

local: all
	 $(R_PATH) CMD INSTALL -l $(LOCAL_INST) $(PACKAGEFOLDER)

install : all
	 $(R_PATH) CMD INSTALL $(PACKAGEFOLDER)

dist : all
	 $(R_PATH) CMD build $(PACKAGEFOLDER)

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
	(cd $(SRCDIR)/f2c; make clean)
