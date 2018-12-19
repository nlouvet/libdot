/*

Created on Wed, 19 Dec 2018 12:12:38 +0000

Copyright (C) 2018 Nicolas Louvet

This file is part of the libdot library 

The libdot library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The libdot library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the hplll library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <mpfr.h>
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#include "mpfr_tools.h"
#include "timing.h"
#include "dotprod.h"
#include "vec.h"
#include "eft.h"

const int nb_tests = 10000;

void testfct_d(double (* const ptrfct)(const double *, const double *, int), const char fctname[], int minlen, int maxlen, int step);
void testfct_dd(dd_real (* const ptrfct)(const dd_real *, const dd_real *, int), const char fctname[], int minlen, int maxlen, int step);
void testfct_qd(qd_real (* const ptrfct)(const qd_real *, const qd_real *, int), const char fctname[], int minlen, int maxlen, int step);

int numfcts_d = 6;
const char *namefct_d[] = {"dotprod", "dotprod_pw", "dotprod_vec", "dotprod2", "dotprod2_vec", "dotprod2_dd"};
double (* const ptrfct_d[])(const double *, const double *, int) = {dotprod, dotprod_pw, dotprod_vec, dotprod2, dotprod2_vec, dotprod2_dd};
  
int numfcts_dd = 6;
const char *namefct_dd[] = {"dotprod_dd", "dotprod_comp2", "dotprod_extr2", "dotprod_comp2_vec", "dotprod_dd_vec", "dotprod_extr2_vec"};
dd_real (* const ptrfct_dd[])(const dd_real *, const dd_real *, int) = {dotprod_dd, dotprod_comp2, dotprod_extr2, dotprod_comp2_vec, dotprod_dd_vec, dotprod_extr2_vec};

int numfcts_qd = 5;
const char *namefct_qd[] = {"dotprod_qd", "dotprod_comp4", "dotprod_extr4", "dotprod_comp4_vec", "dotprod_extr4_vec"};
qd_real (* const ptrfct_qd[])(const qd_real *, const qd_real *, int) = {dotprod_qd, dotprod_comp4, dotprod_extr4, dotprod_comp4_vec, dotprod_extr4_vec};

//
// main function :D
//
int main(int argc, const char* argv[]) {
  int minlen = 500, maxlen = 500, step = 1;
  
  setaffinity(0);
  srand(time(NULL));
  
  if(argc > 1) {
    int i;
    
    for(i = 0; i < numfcts_d; i++) if(strcmp(argv[1], namefct_d[i]) == 0) break;
    if(i < numfcts_d) {
      testfct_d(ptrfct_d[i], namefct_d[i], minlen, maxlen, step);
      return(0);
    }

    for(i = 0; i < numfcts_dd; i++) if(strcmp(argv[1], namefct_dd[i]) == 0) break;
    if(i < numfcts_dd) {
      testfct_dd(ptrfct_dd[i], namefct_dd[i], minlen, maxlen, step);
      return(0);
    }

    for(i = 0; i < numfcts_qd; i++) if(strcmp(argv[1], namefct_qd[i]) == 0) break;
    if(i < numfcts_qd) {
      testfct_qd(ptrfct_qd[i], namefct_qd[i], minlen, maxlen, step);
      return(0);
    }
  }

  printf("Functions available for timing: ");
  for(int i = 0; i < numfcts_d; i++) printf("%s\n", namefct_d[i]);
  for(int i = 0; i < numfcts_dd; i++) printf("%s\n", namefct_dd[i]);
  for(int i = 0; i < numfcts_qd; i++) printf("%s\n", namefct_qd[i]);

  return(0);
}

//
// Testing functions
//

static inline double myrand_d() {
  return 1.0 - 2.0 * ((double)rand())/((double)RAND_MAX);
}

static inline dd_real myrand_dd() {
  double x, y, u, v;
  x = 1.0 - 2.0 * rand()/(double)RAND_MAX;
  y = 1e-16 * (1.0 - 2.0 * rand()/(double)RAND_MAX);
  TwoSum(u, v, x, y);  
  return dd_real(u, v);
}

static inline qd_real myrand_qd() {
  double u;
  qd_real x(0.0);

  u = myrand_d(); x += u;
  u = myrand_d(); x += u * 1e-16;
  u = myrand_d(); x += u * 1e-32;
  u = myrand_d(); x += u * 1e-48; 
  return x;
}

static inline int myrand_int(int n) {
  return rand() % n;
}

//
// Function for timing dot products with double precision inputs and output
//
void testfct_d(double (* const ptrfct)(const double *, const double *, int), const char fctname[], int minlen, int maxlen, int step) {
  FILE *fd;
  uint64_t t, offset, *measure;
  int len, i, j;    
  double *x_d = NULL, *y_d = NULL, z_d;


  if( !(fd=fopen("/dev/null", "a")) ) fprintf(stderr, "Error with fopen /dev/null.\n");
  if( !(measure=(uint64_t *)malloc(sizeof(uint64_t)*nb_tests)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&x_d, VEC_ALIGN_CST, maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&y_d, VEC_ALIGN_CST, maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
  
  // 'Calibrating' the time stamp counter...
  for(i=0; i<nb_tests; i++) {
    t = rdtscA();
    measure[i] = rdtscB() - t;
  }
  qsort(measure, nb_tests, sizeof (uint64_t), reverse_compare_uint64);
  offset = measure[0];
  printf("# offset = %" PRId64 "\n", offset);
  fprintf(fd, "%" PRId64 "\n", offset);
  
  // Taking the measures.
  printf("# %6s %8s\n", "len", fctname);
  for(len = minlen; len <= maxlen; len += step) {
    printf("%8d ", len);
    
    for(i=0; i<nb_tests; i++) {
      for(j=0; j<len; j++) {
	x_d[j] = myrand_d();
	y_d[j] = myrand_d();
      }
      t = rdtscA();
      z_d = ptrfct(x_d, y_d, len);
      measure[i] = rdtscB() - t - offset;
      fprintf(fd, "%f\n", z_d);
    }
    print_measure_tsc(measure, nb_tests, 0);
    putchar('\n');
  }
  
  free(x_d);
  free(y_d);
  free(measure);
  fclose(fd);
}

//
// Function for timing dot products with double-double inputs and output, and double-double computations.
//
void testfct_dd(dd_real (* const ptrfct)(const dd_real *, const dd_real *, int), const char fctname[], int minlen, int maxlen, int step) {
  FILE *fd;
  uint64_t t, offset, *measure;
  int len, i;
  dd_real *x_dd = NULL, *y_dd = NULL, z_dd;
  
  if( !(fd=fopen("/dev/null", "a")) ) fprintf(stderr, "Error with fopen /dev/null.\n");
  if( !(measure=(uint64_t *)malloc(sizeof(uint64_t)*nb_tests)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&x_dd, VEC_ALIGN_CST, maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&y_dd, VEC_ALIGN_CST, maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
 
  // 'Calibrating' the time stamp counter...
  for(i=0; i<nb_tests; i++) {
    t = rdtscA();
    measure[i] = rdtscB() - t;
  }
  qsort(measure, nb_tests, sizeof (uint64_t), reverse_compare_uint64);
  offset = measure[0];
  printf("# offset = %" PRId64 "\n", offset);
  fprintf(fd, "%" PRId64 "\n", offset);

  // Taking the measures.
  printf("# %6s %s\n", "len", fctname);
  for(len = minlen; len <= maxlen; len += step) {
    printf("%8d ", len);
    for(i=0; i<nb_tests; i++) {
      x_dd[myrand_int(len)] = myrand_d();
      y_dd[myrand_int(len)] = myrand_d();
      t = rdtscA();
      z_dd = ptrfct(x_dd, y_dd, len);
      measure[i] = rdtscB() - t - offset;
    }
    print_measure_tsc(measure, nb_tests, 0);
    putchar('\n');
  }
  
  free(x_dd);
  free(y_dd);
  free(measure);
  fclose(fd);
}

//
// Function for timing dot products with quad-double inputs and output, and quad-double computations.
//
void testfct_qd(qd_real (* const ptrfct)(const qd_real *, const qd_real *, int), const char fctname[], int minlen, int maxlen, int step) {
  FILE *fd;
  uint64_t t, offset, *measure;
  int len, i;
  qd_real *x_qd = NULL, *y_qd = NULL, z_qd;
    
  if( !(fd=fopen("/dev/null", "a")) ) fprintf(stderr, "Error with fopen /dev/null.\n");
  if( !(measure=(uint64_t *)malloc(sizeof(uint64_t)*nb_tests)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&x_qd, VEC_ALIGN_CST, maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
  if( posix_memalign((void **)&y_qd, VEC_ALIGN_CST, maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
 
  // 'Calibrating' the time stamp counter...
  for(i=0; i<nb_tests; i++) {
    t = rdtscA();
    measure[i] = rdtscB() - t;
  }
  qsort(measure, nb_tests, sizeof (uint64_t), reverse_compare_uint64);
  offset = measure[0];
  printf("# offset = %" PRId64 "\n", offset);
  fprintf(fd, "%" PRId64 "\n", offset);

  // Taking the measures.

  printf("# %6s %s\n", "len", fctname);
  for(len = minlen; len <= maxlen; len += step) {
    printf("%8d ", len);
    for(i=0; i<nb_tests; i++) {
      x_qd[myrand_int(len)] = myrand_d();
      y_qd[myrand_int(len)] = myrand_d();
      t = rdtscA();
      z_qd = ptrfct(x_qd, y_qd, len);
      measure[i] = rdtscB() - t - offset;
    }
    print_measure_tsc(measure, nb_tests, 0);
    putchar('\n');
  }
  
  free(x_qd);
  free(y_qd);
  free(measure);
  fclose(fd);
}
