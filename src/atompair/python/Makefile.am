

PACKAGENAME = descriptors
#SRCDIR = src
SRCDIR = .
MODULENAME = CDescriptors
BUILT_SOURCES = $(SRCDIR)/py_wrap.cc
CCSOURCEDIR = ..
SWIG_SOURCES = $(CCSOURCEDIR)/swig.i
CCSOURCES = $(CCSOURCEDIR)/desc.h $(CCSOURCEDIR)/debug.h $(CCSOURCEDIR)/molecule.h $(CCSOURCEDIR)/script.h $(CCSOURCEDIR)/simpledb.h $(CCSOURCEDIR)/desc.cc $(CCSOURCEDIR)/molecule.cc $(CCSOURCEDIR)/script.cc $(CCSOURCEDIR)/simpledb.cc $(CCSOURCEDIR)/formats.cc $(CCSOURCEDIR)/similarity.cc $(CCSOURCEDIR)/db_build.h $(CCSOURCEDIR)/db_build.cc
LIBS=

pkgpyexecdir = $(pyexecdir)/$(PACKAGENAME)
pkgpythondir = $(pythondir)/$(PACKAGENAME)
pkgpython_PYTHON =  $(MODULENAME).py $(PACKAGENAME)/__init__.py
pkgpyexec_LTLIBRARIES  = _CDescriptors.la
_CDescriptors_la_SOURCES = $(SRCDIR)/py_wrap.cc  $(SWIG_SOURCES) $(CCSOURCES)
_CDescriptors_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I/usr/include/python2.6 -I$(CCSOURCEDIR) 
_CDescriptors_la_LDFLAGS = -module



$(SRCDIR)/py_wrap.cc : $(SWIG_SOURCES)
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(CCSOURCEDIR) -module $(MODULENAME) -o $(SRCDIR)/py_wrap.cc -python $(SWIG_SOURCES)

clean-generic:
	rm -f $(MODULENAME).py py_wrap.cc


#mv $(SRCDIR)/$(MODULENAME).py $(PACKAGEFOLDER)/
#mkdir $(SRCDIR)
