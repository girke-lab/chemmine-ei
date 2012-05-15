#!/usr/bin/env python
from distutils.core import setup, Extension
import os
##############################################################
# uncomment this if you encounter error about `cc1plus'
os.environ['CC'] = 'g++'

##############################################################
# change the following to True to use OpenBabel
use_openbabel = False
# OpenBabel library location if you want to manually set it
openbabel_libdir = '/home/ycao/opt64/lib'
# OpenBabel header location if you want to manually set it
openbabel_incdir = '/home/ycao/opt/src/openbabel-2.1.1/include'

##############################################################
library_dirs=[]
define_macros=[]
include_dirs=[]
libraries = []
if use_openbabel:
	if openbabel_libdir:
		library_dirs.append(openbabel_libdir)
	if openbabel_incdir:
		include_dirs.append(openbabel_incdir)
	libraries.append('openbabel')
	define_macros.append(('HAS_OPENBABEL', None))


setup(name='descriptors',
      version='1.0',
      packages=['descriptors'],
	  ext_modules = [Extension('descriptors._CDescriptors', [
	  	'src/py_wrap.cc',
	  	'src/desc.cc',
	  	'src/formats.cc',
	  	'src/molecule.cc',
	  	'src/script.cc',
	  	'src/similarity.cc',
	  	'src/simpledb.cc',
	  	'src/db_build.cc',
	  ],
		  define_macros=define_macros,
		  library_dirs=library_dirs,
		  libraries=libraries,
		  include_dirs=include_dirs
	  )]
      )
