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
#include <qd/dd_real.h>
#include <qd/qd_real.h>

#include <cstdio>
#include <cstring>

#include "timing.h"
#include "libdot.h"

static inline double myrand_d() {
  return 1.0 - 2.0 * ((double)rand())/((double)RAND_MAX);
}

int main(int argc, const char* argv[]) {
  FILE *fd = NULL, *fdnull = NULL;
  char *file = NULL;

  if(argc < 2) {
    fd = stdout;
    printf("#output: stdout\n");
  }
  else {
    file = static_cast<char *>(malloc(sizeof(char) * (strlen(argv[1])+1+32)));
    sprintf(file, "perf-%s.dat", argv[1]);
    printf("#output: %s\n", file);
    if(!(fd = fopen(file, "w"))) {
      perror("fopen the output file");
      return(-1);
    }
  }

  if(!(fdnull = fopen("/dev/null", "a"))) {
    perror("fopen /dev/null");
    return(-1);
  }
  
  int numfcts_d = 5;
  const char *namefct_d[] = {"dotprod", "dotprod_vec", "dotprod2", "dotprod2_vec", "dotprod2_dd"};
  double (* const ptrfct_d[])(const double *, const double *, int) = {dotprod, dotprod_vec, dotprod2, dotprod2_vec, dotprod2_dd};

  int numfcts_dd = 6;
  const char *namefct_dd[] = {"dotprod_dd", "dotprod_dd_vec", "dotprod_comp2", "dotprod_extr2", "dotprod_comp2_vec", "dotprod_extr2_vec"};
  dd_real (* const ptrfct_dd[])(const dd_real *, const dd_real *, int) = {dotprod_dd, dotprod_dd_vec, dotprod_comp2, dotprod_extr2, dotprod_comp2_vec, dotprod_extr2_vec};

  int numfcts_qd = 5;
  const char *namefct_qd[] = {"dotprod_qd", "dotprod_comp4", "dotprod_extr4", "dotprod_comp4_vec", "dotprod_extr4_vec"};
  qd_real (* const ptrfct_qd[])(const qd_real *, const qd_real *, int) = {dotprod_qd, dotprod_comp4, dotprod_extr4, dotprod_comp4_vec, dotprod_extr4_vec};

  srand(time(NULL));
  
  {
    uint64_t t, t1, t2, offset, *measure;
    int nb_tests = 200, maxlen = 100, len, i, j, k;

    double *x_d, *y_d, z;
    x_d = y_d = NULL;

    dd_real *x_dd, *y_dd, z_dd;
    x_dd = y_dd = NULL;

    qd_real *x_qd, *y_qd, z_qd;
    x_qd = y_qd = NULL;

    fprintf(fd, "# Timings in clock cycles with rdtsc\n");
    if(file != NULL) fprintf(fd, "# filename: %s\n", file);
    else fprintf(fd, "# filename: stdout\n");
    fprintf(fd, "# nb_tests: %d\n", nb_tests);
    fprintf(fd, "# col  1: %17s\n", "vector len");
    for(k=0; k < numfcts_d; k++)  fprintf(fd, "# col %2d: %17s \n", k + 2, namefct_d[k]);
    for(k=0; k < numfcts_dd; k++) fprintf(fd, "# col %2d: %17s \n", numfcts_d + k + 2, namefct_dd[k]);
    for(k=0; k < numfcts_qd; k++) fprintf(fd, "# col %2d: %17s \n", numfcts_d + numfcts_dd + k + 2, namefct_qd[k]);
    
    if( !(measure=(uint64_t *)malloc(sizeof(uint64_t)*nb_tests)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&x_d, libdot_align(), maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_d, libdot_align(), maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&x_dd, libdot_align(), maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_dd, libdot_align(), maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&x_qd, libdot_align(), maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_qd, libdot_align(), maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");

    for(i=0; i<nb_tests; i++) {
      rdtsc(t1);
      rdtscp(t2);
      measure[i] = t2-t1;
    }
    qsort(measure, nb_tests, sizeof (uint64_t), reverse_compare_uint64);
    offset = measure[0];
    fprintf(fd, "# offset = %" PRId64 "\n", offset);

    fprintf(fd, "#   col 1");
    for(k=2; k<=1+numfcts_d + numfcts_dd + numfcts_qd; k++) fprintf(fd, "   col %2d", k);
    fprintf(fd, "\n");

    for(len = 10; len <= maxlen; len += 5) {
      fprintf(fd, "%8d ", len);

      for(k = 0; k < numfcts_d; k++) {
        for(i = 0; i < nb_tests; i++) {
          for(j = 0; j < len; j++) {
            x_d[j] = myrand_d();
            y_d[j] = myrand_d();
          }
	  rdtsc(t1);
          z = ptrfct_d[k](x_d, y_d, len);
	  rdtscp(t2);
	  measure[i] = t2-t1;
          fprintf(fdnull, "%f\n", z);
        }
        fprint_measure_tsc(fd, measure, nb_tests, 0);
      }
    
      for(k = 0; k < numfcts_dd; k++) {
        for(i = 0; i < nb_tests; i++) {
          for(j = 0; j < len; j++) {
            x_dd[j] = ddrand();
            y_dd[j] = ddrand();
          }
	  rdtsc(t1);
          z_dd = ptrfct_dd[k](x_dd, y_dd, len);
	  rdtscp(t2);
	  measure[i] = t2-t1;
          fprintf(fdnull, "%f %f\n", z_dd._hi(), z_dd._lo());
        }
        fprint_measure_tsc(fd, measure, nb_tests, offset);
      }

      for(k = 0; k < numfcts_qd; k++) {
        for(i = 0; i < nb_tests; i++) {
          for(j = 0; j < len; j++) {
            x_qd[j] = qdrand();
            y_qd[j] = qdrand();
          }
	  rdtsc(t1);
          z_qd = ptrfct_qd[k](x_qd, y_qd, len);
	  rdtscp(t2);
	  measure[i] = t2-t1;
          fprintf(fdnull, "%f %f %f %f\n", z_qd[0], z_qd[1], z_qd[2], z_qd[3]);
        }
        fprint_measure_tsc(fd, measure, nb_tests, offset);
      }

      fprintf(fd, "\n");
    }
    fprintf(fd, "# end\n");
    
    free(x_d); free(y_d);
    free(x_dd); free(y_dd);
    free(x_qd); free(y_qd);
    free(measure);
  }

  fclose(fdnull);
  if(fd != stdout) fclose(fd);
  if(file != NULL) free(file);
  
  return(0);
}
