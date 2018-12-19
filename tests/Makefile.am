
bin_PROGRAMS = timetest

if ENABLE_PAPI
bin_PROGRAMS += perftest
endif 

if ENABLE_MPFR
bin_PROGRAMS += acctest
endif

timetest_SOURCES = timetest.cpp timing.cpp timing.h
timetest_CXXFLAGS = -O2 -I../src/
timetest_LDADD = ../src/libdot.la -lqd

if ENABLE_MPFR
acctest_SOURCES = acctest.cpp
acctest_CXXFLAGS = -O2 -I../src/
acctest_LDADD = ../src/libdot.la libtools.la -lgmp -lmpfr
endif

perftest_SOURCES = perftest.cpp
perftest_CXXFLAGS = -O2 -I../src/
perftest_LDADD = ../src/libdot.la -lqd -lpapi

# indivtest_SOURCES = indivtest.cpp
# indivtest_CXXFLAGS = $(OPT1_CXXFLAGS) -I../src/
# indivtest_LDADD = ../src/libdot.la libtools.la -lpapi

if ENABLE_MPFR
noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = mpfr_tools.cpp mpfr_tools.h gendot.cpp gendot.h
libtools_la_CXXFLAGS = $(AM_CXXFLAGS) -I../src/
endif
