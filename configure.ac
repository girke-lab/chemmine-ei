#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([ei], [1.0])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CC
AM_PATH_PYTHON([2.5])
#AX_PKG_SWIG([1.3.17])
#AX_SWIG_PYTHON

# Checks for libraries.
AC_CHECK_LIB([f2c], [pow])
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([readline], [rl_complete])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor ftruncate gettimeofday memset pow sqrt])

AC_CONFIG_FILES([Makefile
					  src/Makefile
                 src/atompair/Makefile
                 src/euclidean-search/Makefile
                 src/lshkit/Makefile
                 src/lshkit/src/Makefile
                 src/lshkit/test/Makefile
                 src/matrix/Makefile])
                 #src/atompair/R/Makefile
                 #src/atompair/python/Makefile
                 #src/embedder/Makefile
                 #src/euclidean-search/Makefile
                 #src/lshkit/Makefile
                 #src/lshkit/src/Makefile
                 #src/lshkit/test/Makefile
                 #src/matrix/Makefile])
AC_OUTPUT
