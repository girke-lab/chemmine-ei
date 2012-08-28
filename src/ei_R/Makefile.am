PACKAGEFOLDER = eiR
SRCDIR = $(PACKAGEFOLDER)/src
CXXFLAGS +=  -fPIC -DNO_MAIN
LIBS = -lm -lf2c  

MODULENAME = eiR
CCSOURCEDIR = ../atompair
CCSOURCES = $(CCSOURCEDIR)/molecule.cc $(CCSOURCEDIR)/formats.cc $(CCSOURCEDIR)/desc.cc $(CCSOURCEDIR)/similarity.cc $(CCSOURCEDIR)/simpledb.cc $(CCSOURCEDIR)/search.cc $(CCSOURCEDIR)/profiling.cc $(CCSOURCEDIR)/db_build.h $(CCSOURCEDIR)/debug.h $(CCSOURCEDIR)/desc.h $(CCSOURCEDIR)/molecule.h $(CCSOURCEDIR)/profiling.h $(CCSOURCEDIR)/script.h $(CCSOURCEDIR)/search.h $(CCSOURCEDIR)/simpledb.h $(CCSOURCEDIR)/db_build.cc $(CCSOURCEDIR)/db2db_distance.cc ../matrix/format.cc ../embedder/embedCoords.cc ../embedder/solver.cc ../embedder/solver.h ../embedder/routines.c ../embedder/config.h ../atompair/db2db_distance_r.cc ../atompair/helpers.h ../atompair/helpers.cc ../atompair/db_subset.cc
EXTRA_DIST = $(MODULENAME)/DESCRIPTION $(MODULENAME)/NAMESPACE $(MODULENAME)/README  $(MODULENAME)/src/Makevars

LOCAL_INST=~/testing

empty=
space=$(empty) $(empty)
ESC_CXXFLAGS=$(subst $(space),\$(space),$(CXXFLAGS))
ESC_LIBS=$(subst $(space),\$(space),$(LIBS))

if BUILD_EIR
all : $(SRCDIR)/r_wrap.cpp ccsource

$(SRCDIR)/r_wrap.cpp : swig.i
	mkdir -p $(MODULENAME)/R $(MODULENAME)/src
	$(SWIG) -c++ -I$(CCSOURCEDIR) -module $(MODULENAME) -o $(SRCDIR)/r_wrap.cpp -r swig.i
	mv $(SRCDIR)/$(MODULENAME).R $(PACKAGEFOLDER)/R/

ccsource : $(CCSOURCES)
	cp -a $(CCSOURCES) $(SRCDIR)

check: all
	MAKEFLAGS='CXXFLAGS=$(ESC_CXXFLAGS) LDFLAGS=$(ESC_LIBS)' $(R_PATH) CMD check $(PACKAGEFOLDER)
test: all local
	export RCMDCHECK=FALSE; \
	export R_LIBS_USER=$(LOCAL_INST); \
	MAKEFLAGS='CXXFLAGS=$(ESC_CXXFLAGS) LDFLAGS=$(ESC_LIBS)' $(R_PATH) --vanilla --slave < eiR/tests/test.R

local: all
	MAKEFLAGS='CXXFLAGS=$(ESC_CXXFLAGS) LDFLAGS=$(ESC_LIBS)' $(R_PATH) CMD INSTALL -l $(LOCAL_INST) $(PACKAGEFOLDER)

install : all
	MAKEFLAGS='CXXFLAGS=$(ESC_CXXFLAGS) LDFLAGS=$(ESC_LIBS)' $(R_PATH) CMD INSTALL $(PACKAGEFOLDER)
endif


clean-generic :
	rm -f $(SRCDIR)/r_wrap.cpp
	rm -f $(PACKAGEFOLDER)/R/$(MODULENAME).R
	(cd $(SRCDIR);\
	rm -f *.so *.o;\
	for i in $(CCSOURCES); do \
		rm -f `basename $$i`;\
	done)
