
PACKAGEFOLDER = ChemmineRpp
SRCDIR = $(PACKAGEFOLDER)/src
CXXFLAGS +=  -fPIC
MODULENAME = ChemmineRpp
CCSOURCEDIR = ..
CCSOURCES = $(CCSOURCEDIR)/desc.h $(CCSOURCEDIR)/debug.h $(CCSOURCEDIR)/molecule.h $(CCSOURCEDIR)/script.h $(CCSOURCEDIR)/simpledb.h $(CCSOURCEDIR)/desc.cc $(CCSOURCEDIR)/molecule.cc $(CCSOURCEDIR)/script.cc $(CCSOURCEDIR)/simpledb.cc $(CCSOURCEDIR)/formats.cc $(CCSOURCEDIR)/similarity.cc $(CCSOURCEDIR)/db_build.h $(CCSOURCEDIR)/db_build.cc
EXTRA_DIST = $(MODULENAME)/DESCRIPTION $(MODULENAME)/NAMESPACE $(MODULENAME)/README  $(MODULENAME)/src/Makevars


if BUILD_CHEMMINER
all : $(SRCDIR)/r_wrap.cc ccsource

$(SRCDIR)/r_wrap.cc : ../swig.i
	mkdir -p $(MODULENAME)/R
	$(SWIG) -c++ -I$(CCSOURCEDIR) -module $(MODULENAME) -o $(SRCDIR)/r_wrap.cc -r ../swig.i
	mv $(SRCDIR)/$(MODULENAME).R $(PACKAGEFOLDER)/R/

ccsource : $(CCSOURCES)
	cp $(CCSOURCES) $(SRCDIR)

install : all
	$(R_PATH) CMD INSTALL $(PACKAGEFOLDER)
endif


clean-generic :
	rm -f $(SRCDIR)/r_wrap.cc
	rm -f $(PACKAGEFOLDER)/R/$(MODULENAME).R
	(cd $(SRCDIR);\
	for i in $(CCSOURCES); do \
	rm -f `basename $$i`;\
	done)
