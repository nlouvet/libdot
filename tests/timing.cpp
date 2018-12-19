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

#include "timing.h"


void setaffinity(int cpu) {
#ifdef _GNU_SOURCE
  cpu_set_t mask1, mask2;
  printf("# Trying to set the process's CPU affinity mask to CPU%d: ", cpu);
  CPU_ZERO(&mask1);
  CPU_SET(cpu, &mask1);
  if(sched_setaffinity(0, sizeof(mask1), &mask1)) perror("Problem with sched_setaffinity");
  if(sched_getaffinity(0, sizeof(mask2), &mask2)) perror("Problem with sched_getaffinity");
  if(CPU_EQUAL(&mask1, &mask2)) printf("Success.\n");
  else printf("Failure.\n");
#else
  printf("# No way to set the process's CPU affinity mask to CPU%d...\n", cpu);
#endif
}

void tic(tv_t *tv1) {
  gettimeofday(tv1, NULL);
}

double toc(tv_t tv1) {
  tv_t tv2;
  gettimeofday(&tv2, NULL);
  return((double)(tv2.tv_sec-tv1.tv_sec)*1000.0 + (double)(tv2.tv_usec-tv1.tv_usec)/1000.0);
}

int reverse_compare_double(const void *a, const void *b) {
  const double *aa = (const double *)a;
  const double *bb = (const double *)b;
  return (*aa < *bb) ? -1 : (*aa == *bb) ? 0 : +1;
}

void print_measure_cw(double *measure, int nb_tests) {
  double avg;
  int i;

  qsort(measure, nb_tests, sizeof (double), reverse_compare_double);
  for(i=0, avg=0.0; i<nb_tests; i++) avg += measure[i];
  avg /= nb_tests;
  printf("%7.2f %7.2f %7.2f", avg, measure[0], measure[nb_tests-1]);
}

int reverse_compare_uint64 (const void *a, const void *b)
{
  const uint64_t * aa = (const uint64_t *)a;
  const uint64_t * bb = (const uint64_t *)b;
  return (*aa < *bb) ? -1 : (*aa == *bb) ? 0 : +1;
}

void fprint_measure_tsc(FILE *fd, uint64_t *measure, int nb_tests, uint64_t offset) {
  int nb_kept = (nb_tests*9)/10;
  double avg;
  int i;

  for(i=0; i<nb_tests; i++) measure[i] -= offset;
  qsort(measure, nb_tests, sizeof (uint64_t), reverse_compare_uint64);
  for(i=0, avg=0.0; i<nb_kept; i++) avg += (double)(measure[i]);
  avg = ceil(avg/(double)(nb_kept));
  //printf("%8.0f %8" PRId64 " %8" PRId64 " ", avg, measure[0], measure[nb_kept-1]);
  fprintf(fd, "%8.0f ", avg);
}

void print_measure_tsc(uint64_t *measure, int nb_tests, uint64_t offset) {
  fprint_measure_tsc(stdout, measure, nb_tests, offset);
}

