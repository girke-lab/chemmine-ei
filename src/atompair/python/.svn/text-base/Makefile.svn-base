PACKAGEFOLDER = descriptors
SRCDIR = src
PYINC = -I/usr/include/python2.4
CXXFLAGS += $(PYINC) -I.. -fPIC
MODULENAME = CDescriptors
LINKER = g++ 
CCSOURCEDIR = ..
CCSOURCES = $(CCSOURCEDIR)/desc.h $(CCSOURCEDIR)/debug.h $(CCSOURCEDIR)/molecule.h $(CCSOURCEDIR)/script.h $(CCSOURCEDIR)/simpledb.h $(CCSOURCEDIR)/desc.cc $(CCSOURCEDIR)/molecule.cc $(CCSOURCEDIR)/script.cc $(CCSOURCEDIR)/simpledb.cc $(CCSOURCEDIR)/formats.cc $(CCSOURCEDIR)/similarity.cc $(CCSOURCEDIR)/db_build.h $(CCSOURCEDIR)/db_build.cc

all : $(SRCDIR)/py_wrap.cc ccsource

$(SRCDIR)/py_wrap.cc : ../swig.i
	swig -c++ -module $(MODULENAME) -o $@ -python $^
	mv $(SRCDIR)/$(MODULENAME).py $(PACKAGEFOLDER)/

ccsource : $(CCSOURCES)
	cp $(CCSOURCES) $(SRCDIR)

dist : all
	rm -rf dist
	mkdir dist
	cp -r descriptors setup.py src README dist
	tar -czf py-descriptor.tar.gz dist
	rm -rf dist

clean :
	rm -f $(SRCDIR)/py_wrap.cc
	rm -f $(PACKAGEFOLDER)/$(MODULENAME).py $(PACKAGEFOLDER)/*.pyc
	rm -rf dist
	rm -rf build
	(cd $(SRCDIR);\
	for i in $(CCSOURCES); do \
	rm -f `basename $$i`;\
	done)
