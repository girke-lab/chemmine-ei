PACKAGEFOLDER = eiR
SRCDIR = $(PACKAGEFOLDER)/src
CXXFLAGS +=  -fPIC
MODULENAME = eiR
CCSOURCEDIR = ../atompair
CCSOURCES = $(CCSOURCEDIR)/molecule.cc $(CCSOURCEDIR)/formats.cc $(CCSOURCEDIR)/desc.cc $(CCSOURCEDIR)/similarity.cc $(CCSOURCEDIR)/simpledb.cc $(CCSOURCEDIR)/search.cc $(CCSOURCEDIR)/profiling.cc $(CCSOURCEDIR)/db_build.h $(CCSOURCEDIR)/debug.h $(CCSOURCEDIR)/desc.h $(CCSOURCEDIR)/molecule.h $(CCSOURCEDIR)/profiling.h $(CCSOURCEDIR)/script.h $(CCSOURCEDIR)/search.h $(CCSOURCEDIR)/simpledb.h $(CCSOURCEDIR)/db_build.cc
EXTRA_DIST = $(MODULENAME)/DESCRIPTION $(MODULENAME)/NAMESPACE $(MODULENAME)/README  $(MODULENAME)/src/Makevars

LOCAL_INST=~/testing

if BUILD_EIR
all : $(SRCDIR)/r_wrap.cpp ccsource

$(SRCDIR)/r_wrap.cpp : swig.i
	mkdir -p $(MODULENAME)/R $(MODULENAME)/src
	$(SWIG) -c++ -I$(CCSOURCEDIR) -module $(MODULENAME) -o $(SRCDIR)/r_wrap.cpp -r swig.i
	mv $(SRCDIR)/$(MODULENAME).R $(PACKAGEFOLDER)/R/

ccsource : $(CCSOURCES)
	cp -a $(CCSOURCES) $(SRCDIR)

check: all
	$(R_PATH) CMD check $(PACKAGEFOLDER)
test: all local
	export RCMDCHECK=FALSE; \
	export R_LIBS_USER=$(LOCAL_INST); \
	${R_PATH} --vanilla --slave < eiR/tests/test.R


local: all
	$(R_PATH) CMD INSTALL -l $(LOCAL_INST) $(PACKAGEFOLDER)

install : all
	$(R_PATH) CMD INSTALL $(PACKAGEFOLDER)
endif


clean-generic :
	rm -f $(SRCDIR)/r_wrap.cpp
	rm -f $(PACKAGEFOLDER)/R/$(MODULENAME).R
	(cd $(SRCDIR);\
	rm -f *.so *.o;\
	for i in $(CCSOURCES); do \
		rm -f `basename $$i`;\
	done)
