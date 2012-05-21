
AM_CXXFLAGS = -O3 -mtune=nocona
LIBS = $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB)  -lm -ldb  -lgsl -lgslcblas 

# this could work, but includes stuff we don't need, including -lf2c which must be statically linked
#LIBS += $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB)

