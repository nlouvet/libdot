lib_LTLIBRARIES = libdot.la
libdot_la_SOURCES = 
nodist_EXTRA_libdot_la_SOURCES = dummy.cpp
libdot_la_LIBADD = libdot1.la libdot2.la

noinst_LTLIBRARIES = libdot1.la libdot2.la

#############################################

libdot1_la_SOURCES = libdot.h \
                     eft.h \
                     vec.h \
                     libdot_align.cpp \
                     dotprod_vec.cpp \
                     dotprod2.cpp \
                     dotprod2_vec.cpp \
                     dotprod2_dd.cpp \
                     dotprod_dd_vec.cpp \
                     dotprod_dd.cpp \
                     dotprod_comp2.cpp \
                     dotprod_comp2_vec.cpp \
                     dotprod_extr2.cpp \
                     dotprod_extr2_vec.cpp \
                     dotprod_qd.cpp \
                     dotprod_comp4.cpp \
                     dotprod_extr4.cpp \
                     dotprod_comp4_vec.cpp \
                     dotprod_extr4_vec.cpp
libdot1_la_CPPFLAGS =
libdot1_la_CXXFLAGS = $(CXXFLAGS1)

##############################################

libdot2_la_SOURCES = dotprod.cpp
libdot2_la_CPPFLAGS =
libdot2_la_CXXFLAGS = $(CXXFLAGS2)
