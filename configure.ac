AC_INIT([libdot], [1.0])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# If you want to overwrite/redefine FLAGS_X, it is
# important to do so before invocating AC_PROG_X
: ${CXXFLAGS=""}

AC_LANG_CPLUSPLUS

AC_PROG_CXX

AM_PROG_AR

LT_INIT

#Recompile pour toute modif dans la chaine autotool
AM_MAINTAINER_MODE(enable)

AC_CONFIG_HEADERS([config.h])

dnl *************************************************************
dnl Selecting compiler options, with respect to CXX and ARCH
dnl *************************************************************

dnl At the end,
dnl - CXXFLAGS1 will be used to compile most of the code, including "fragile" code,
dnl - CXXFLAGs will be set to CXXFLAGS1
dnl - CXXFLAGS2 will be used for "not too fragile" code.

AC_ARG_VAR(ARCH, [specify the target processor for the library])

dnl *********************************************
dnl Compiler options controling code optimization
dnl *********************************************

if $CXX --version | grep g++; then
  COMPXX=g++
elif $CXX --version | grep clang; then
  COMPXX=clang++
elif $CXX --version | grep icpc; then
  COMPXX=icpc
else
  COMPXX=unknown
fi

case "$COMPXX" in
  dnl GCC compiler options
  g++*) OPT1_CXXFLAGS="-O3 -fno-associative-math -ffp-contract=off -fprefetch-loop-arrays -funroll-all-loops"
        OPT2_CXXFLAGS="-O3 -fno-signed-zeros -fno-trapping-math -fassociative-math -ftree-vectorize -fprefetch-loop-arrays -funroll-all-loops"
	COMP="g++"
  ;;
  dnl Clang compiler options
  clang++*) OPT1_CXXFLAGS="-O2"
             OPT2_CXXFLAGS="-O2"
	     COMP="clang++"
  ;;
  dnl Intel c++ compiler options
  icpc*) OPT1_CXXFLAGS="-O3 -fp-model precise"
         OPT2_CXXFLAGS="-O3 -fp-model fast=2"
	 COMP="icpc"
  ;;
  dnl Another compiler
  *) OPT1_CXXFLAGS="-O2"
     OPT2_CXXFLAGS="-O2"
     COMP="def"
  ;;
esac

dnl Adjusting compiler's options to the target architecture

if test "$COMP" = "g++"; then
  case "$ARCH" in
    "sse4")      ARCH_CXXFLAGS="-msse2 -msse4"
    ;;
    "avx2")      ARCH_CXXFLAGS="-mavx2 -mfma"
    ;;
    "avx512")    ARCH_CXXFLAGS="-mavx512 -mfma"
    ;;
    *)           ARCH_CXXFLAGS="-msse2 -msse4"
  esac
fi

if test "$COMP" = "clang++"; then
  case "$ARCH" in
    "sse4")      ARCH_CXXFLAGS="-msse2 -msse4"
    ;;
    "avx2")      ARCH_CXXFLAGS="-mavx2 -mfma"
    ;;
    "avx512")    ARCH_CXXFLAGS="-mavx512 -mfma"
    ;;
    *)           ARCH_CXXFLAGS="-msse2 -msse4"
  esac
fi

if test "$COMP" = "icpc"; then
  case "$ARCH" in
    "sse4")      ARCH_CXXFLAGS="-msse4"
    ;;
    "avx2")      ARCH_CXXFLAGS="-axCORE_AVX2 -fma -D__AVX2__ -D__FP_FAST_FMA"
    ;;
    "avx512")    ARCH_CXXFLAGS="-axCORE_AVX512 -fma -D__AVX2__ -D__FP_FAST_FMA"
    ;;
    *)           ARCH_CXXFLAGS="-msse4"
  esac
fi

CXXFLAGS1="$ARCH_CXXFLAGS $OPT1_CXXFLAGS"
CXXFLAGS2="$ARCH_CXXFLAGS $OPT2_CXXFLAGS"

dnl ***************************************************************
dnl Trying to check that CXX_FLAGS1 indeed preserves "fragile code"
dnl ***************************************************************

AC_MSG_NOTICE([Checking whether the compiler options used indeed preserve floating-point codes])

cat > fcts.cpp <<EOF
double fct0(double a, double b) {
  return a+b;
}

double fct1(double a, double b) {
  return (a+b)-a;
}

double fct2(double a, double b, double c) {
  return (a+b)-c;
}
EOF

cat > test.cpp <<EOF
#include <cstdio>
#include <cmath>

double fct0(double a, double b);
double fct1(double a, double b);
double fct2(double a, double b, double c);

int main(void) {
  double a, b, c, x;
  int numerr = 0;

  printf("Testing if rounding is done to the nearest, with ties to even : ");
  a = 1.0; b = ldexp(1.0, -53);
  x = fct0(a, b);
  if(x != 1.0) {
    printf("no\n");
    numerr++;
  }
  else printf("yes\n");

  printf("Testing whether parentheses are respected : ");
  a = 1.0; b = ldexp(1.0, -53);
  x = fct1(a, b);
  if(x != 0.0){
    printf("no\n");
    numerr++;
  }
  else printf("yes\n");

  printf("Testing if double rounding is avoided : ");
  a = 1.0; b = ldexp(1.0, -53) + ldexp(1.0, -105); c = ldexp(1.0, -52);
  x = fct2(a, b, c);
  if(x != 1.0) {
    printf("no\n");
    numerr++;
  }
  else printf("yes\n");

  
  if(numerr > 0) return 1;
  else return 0;
}
EOF

if ! $CXX $CXXFLAGS1 -c fcts.cpp; then AC_MSG_FAILURE([compilation failed]) fi
if ! $CXX $CXXFLAGS1 test.cpp -o test fcts.o -lm; then AC_MSG_FAILURE([compilation failed]) fi
if ./test; then AC_MSG_NOTICE([test passed]) else AC_MSG_FAILURE([test failed]) fi
rm -f fcts.{o,cpp}

dnl Checking headers needed for using intrinsics

AC_CHECK_HEADER(
 [immintrin.h],
 AC_MSG_NOTICE([immintrin.h usable]),
 AC_MSG_FAILURE([no usable immintrin.h found])
)

dnl *******
dnl Summary
dnl *******

AC_SUBST(CXXFLAGS1)
AC_SUBST(CXXFLAGS2)

AC_MSG_NOTICE([Compiler: $CXX])
AC_MSG_NOTICE([Default options: $CXXFLAGS1])
AC_MSG_NOTICE([Vectorizing options: $CXXFLAGS2])

dnl *************************************************************
dnl qd/dd_real.h and qd/qd_real.h must be available and usable...
dnl *************************************************************

AC_ARG_VAR(QDPATH, [specify the path of the qd library])

if test -n "$QDPATH"; then
  AC_MSG_NOTICE([qd library path provided: $QDPATH])
  CPPFLAGS="$CPPFLAGS -I$QDPATH/include"
  CXXFLAGS="$CXXFLAGS -L$QDPATH/lib"
else
  AC_MSG_NOTICE([no qd library path provided])
fi

dnl An adequate test for the qd library files should be inserted here...

AC_CHECK_HEADER([qd/dd_real.h], , AC_MSG_ERROR([no usable qd/dd_real.h header found]))

AC_CHECK_HEADER([qd/qd_real.h], , AC_MSG_ERROR([no usable qd/qd_real.h header found]))

AC_MSG_CHECKING(dd_real type and linking against the qd library)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <qd/dd_real.h>], [dd_real x = 1.0;])],
  AC_MSG_RESULT([yes]),
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([unable to use the dd_real type or to link against the qd library])
)

AC_MSG_CHECKING(qd_real type and linking against the qd library)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <qd/qd_real.h>], [qd_real x = 1.0;])],
  AC_MSG_RESULT([yes]),
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([unable to use the qd_real type or to link against the qd library])
)

dnl *************************************************************

AC_CHECK_HEADER(papi.h,
      have_papi="yes",
      have_papi="no"
    )

AM_CONDITIONAL([ENABLE_PAPI], [test "$have_papi" = "yes"])

dnl *************************************************************

AC_ARG_ENABLE([mpfr], [--enable-mpfr],
 [case "${enableval}" in
    yes) enable_mpfr="yes" ;;
    no)  enable_mpfr="no" ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-mpfr]) ;;
  esac
 ],
 [enable_mpfr=no]
)

AM_CONDITIONAL([ENABLE_MPFR], [test "$enable_mpfr" = "yes"])

dnl *************************************************************

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 tests/Makefile
])

AC_OUTPUT

dnl AC_PROG_CXX have been used, so the variable CXX is set here to a fonctionning C++ compiler.
dnl Note that AC_PROG_CXX first whether the environnement variable was already set.


