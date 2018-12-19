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
along with the libdot library; see the file COPYING.LESSER.  If not, see
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
#include "libdot.h"
#include <papi.h>

static inline double myrand_d() {
  return 1.0 - 2.0 * ((double)rand())/((double)RAND_MAX);
}

void print_values(const char* namefct, long long values[], int num_events, int nb_tests) {
  printf("%20s ", namefct);
  for(int i = 0; i < num_events; i++) {
    printf("& %12.0f ", (double)values[i]/(double)nb_tests);
  }
  printf("& %12.2f ", (double)values[0] / (double)values[1]);
  printf("\\\\\n");
}

void print_evnames(const char *evname[], int num_events) {
  printf("%20s ", "function");
  for(int i = 0; i < num_events; i++)
    printf("  %12s ", evname[i]);
  printf("  %12s ", "IPC");
  printf("\n");
}

int main() {
  FILE *fdnull;

  if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
    fprintf(stderr, "WARNING: papi init error!\n");

  int numfcts_d = 1;
  const char *namefct_d[] = {"dotprod"};
  double (* const ptrfct_d[])(const double *, const double *, int) = {dotprod};

  int numfcts_dd = 3;
  const char *namefct_dd[] = {"dotprod_dd", "dotprod_comp2", "dotprod_extr2"};
  dd_real (* const ptrfct_dd[])(const dd_real *, const dd_real *, int) = {dotprod_dd, dotprod_comp2, dotprod_extr2};

  int numfcts_qd = 3;
  const char *namefct_qd[] = {"dotprod_qd", "dotprod_comp4", "dotprod_extr4"};
  qd_real (* const ptrfct_qd[])(const qd_real *, const qd_real *, int) = {dotprod_qd, dotprod_comp4, dotprod_extr4};

  int num_events = 2;
  int events[] = {PAPI_TOT_INS, PAPI_TOT_CYC};
  const char *evnames[] = {"TOT_INS", "TOT_CYC"};
  long long values[num_events];

  srand(time(NULL));

  if(!(fdnull = fopen("/dev/null", "a"))) {
    perror("fopen /dev/null");
    return(-1);
  }

  {
    int nb_tests = 100000, maxlen = 500, len, i, j, k;

    double *x_d, *y_d, *t_d, z;
    x_d = y_d = t_d = NULL;

    dd_real *x_dd, *y_dd, *t_dd, z_dd;
    x_dd = y_dd = t_dd = NULL;

    qd_real *x_qd, *y_qd, *t_qd, z_qd;
    x_qd = y_qd = t_qd = NULL;

    if( posix_memalign((void **)&x_d, libdot_align(), maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_d, libdot_align(), maxlen*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&t_d, libdot_align(), 7*sizeof(double)) ) fprintf(stderr, "Memory allocation problem.\n");

    if( posix_memalign((void **)&x_dd, libdot_align(), maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_dd, libdot_align(), maxlen*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&t_dd, libdot_align(), 7*sizeof(dd_real)) ) fprintf(stderr, "Memory allocation problem.\n");

    if( posix_memalign((void **)&x_qd, libdot_align(), maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&y_qd, libdot_align(), maxlen*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");
    if( posix_memalign((void **)&t_qd, libdot_align(), 7*sizeof(qd_real)) ) fprintf(stderr, "Memory allocation problem.\n");

    len = maxlen;
    {
      printf("# len = %3d\n", len);

      print_evnames(evnames, num_events);

      for(k = 0; k < numfcts_d; k++) {

	for(i = 0; i < len; i++) { x_d[i] = myrand_d(); y_d[i] = myrand_d(); }
	for(i = 0; i < 7; i++) { t_d[i] = myrand_d(); }
	if(PAPI_start_counters(events, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
        for(j = 0; j < nb_tests; j++) {
          x_d[j % len] = t_d[(j)   % 7];
	  y_d[j % len] = t_d[(j+1) % 7];
          z = ptrfct_d[k](x_d, y_d, len);
        }
	if(PAPI_stop_counters(values, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
	fprintf(fdnull, "%f\n", z);

	print_values(namefct_d[k], values, num_events, nb_tests);

      }

      for(k = 0; k < numfcts_dd; k++) {

	for(i = 0; i < len; i++) { x_dd[i] = ddrand(); y_dd[i] = ddrand(); }
	for(i = 0; i < 7; i++) { t_dd[i] = ddrand(); }
	if(PAPI_start_counters(events, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
        for(j = 0; j < nb_tests; j++) {
	  x_dd[j % len] = t_dd[(j)   % 7];
	  y_dd[j % len] = t_dd[(j+1) % 7];
          z_dd = ptrfct_dd[k](x_dd, y_dd, len);
        }
	if(PAPI_stop_counters(values, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
	fprintf(fdnull, "%f %f\n", z_dd._hi(), z_dd._lo());

	print_values(namefct_dd[k], values, num_events, nb_tests);

      }

      for(k = 0; k < numfcts_qd; k++) {

	for(i = 0; i < len; i++) { x_qd[i] = qdrand(); y_qd[i] = qdrand(); }
	for(i = 0; i < 7; i++) { t_qd[i] = qdrand(); }
	if(PAPI_start_counters(events, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
        for(j = 0; j < nb_tests; j++) {
          x_qd[j % len] = t_qd[(j)   % 7];
	  y_qd[j % len] = t_qd[(j+1) % 7];
          z_qd = ptrfct_qd[k](x_qd, y_qd, len);
        }
	if(PAPI_stop_counters(values, num_events) != PAPI_OK)
	  fprintf(stderr, "WARNING: problem while starting papi counters.");
	fprintf(fdnull, "%f %f %f %f\n", z_qd[0], z_qd[1], z_qd[2], z_qd[3]);

	print_values(namefct_qd[k], values, num_events, nb_tests);

      }

    }

    free(x_d); free(y_d); free(t_d);
    free(x_dd); free(y_dd); free(t_dd);
    free(x_qd); free(y_qd); free(t_qd);
  }

  fclose(fdnull);

  return(0);
}

