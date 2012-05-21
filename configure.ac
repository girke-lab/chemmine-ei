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
#AX_PKG_SWIG
#AX_SWIG_PYTHON

# Checks for libraries.
AC_CHECK_LIB([m], [pow],[],[AC_MSG_FAILURE])
AC_CHECK_LIB([readline], [rl_complete],[],[AC_MSG_FAILURE])
AC_CHECK_LIB([db], [db_create],[],[AC_MSG_FAILURE])

AC_ARG_WITH([lshkit],AS_HELP_STRING([--without-lshkit],[Build without lshkit]))
AS_IF([test "x$with_lshkit" != "xno"],[
	AX_BOOST_BASE
	AX_BOOST_PROGRAM_OPTIONS
	AX_BOOST_THREAD

	AC_LANG_PUSH([C++])
	AC_CHECK_HEADERS([boost/threadpool.hpp],[],[AC_MSG_FAILURE])
	AC_CHECK_HEADERS([gsl/gsl_integration.h],[],[AC_MSG_FAILURE])
	AC_CHECK_HEADERS([gsl/gsl_multifit.h],[],[AC_MSG_FAILURE])
	AC_LANG_POP([C++])

	AC_CHECK_LIB([gslcblas],[cblas_ctpmv],[],[AC_MSG_FAILURE])
	AC_CHECK_LIB([gsl],[gsl_matrix_set],[],[AC_MSG_FAILURE],[-lgslcblas])
])
AM_CONDITIONAL([BUILD_LSHKIT], [test "x$with_lshkit" != "xno"])

#needs to be last library always
AC_CHECK_LIB([f2c], [main],[],[AC_MSG_FAILURE],[-static])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdint.h stdlib.h string.h sys/time.h unistd.h])
AC_CHECK_HEADERS([db.h],[],[AC_MSG_FAILURE])


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
                 src/embedder/Makefile
                 src/euclidean-search/Makefile
                 src/lshkit/Makefile
                 src/lshkit/src/Makefile
                 src/lshkit/test/Makefile
                 src/matrix/Makefile])
                 #src/atompair/R/Makefile
                 #src/atompair/python/Makefile

AC_OUTPUT
