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
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#ifndef __TIMING__
#define __TIMING__

#ifdef _GNU_SOURCE
#include <sched.h>
#endif
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>

typedef struct timeval tv_t;

#define rdtsc(__time__)                                                                            \
{                                                                                                  \
  uint32_t lo, hi;                                                                                 \
  __asm__ __volatile__ ("cpuid\n"                                                                  \
                        "rdtsc\n"                                                                  \
                        "mov %%edx, %0\n"                                                          \
                        "mov %%eax, %1\n" : "=r"(hi), "=r"(lo) :: "%rax", "%rbx", "%rcx", "%rdx"); \
  __time__ = (uint64_t)hi << 32 | (uint64_t)lo;	                                                   \
}

#define rdtscp(__time__)                                                                   \
{                                                                                          \
  uint32_t lo, hi;                                                                         \
  __asm__ __volatile__ ("rdtscp\n"                                                         \
                        "mov %%edx, %0\n"                                                  \
                        "mov %%eax, %1\n"                                                  \
                        "cpuid\n" : "=r"(hi), "=r"(lo) :: "%rax", "%rbx", "%rcx", "%rdx"); \
  __time__ = (uint64_t)hi << 32 | (uint64_t)lo;	                                           \
}


// CPUID
// CPUID can be executed at any privilege level to serialize instruction
// execution. Serializing instruction execution guarantees that any
// modifications to flags, registers, and memory for previous instructions
// are completed before the next instruction is fetched and executed.

// RDTSC
// The RDTSC instruction is not a serializing instruction.
// It does not necessarily wait until all previous instructions have
// been executed before reading the counter. Similarly, subsequent
// instructions may begin execution before the read operation is performed.
// If software requires RDTSC to be executed only after all previous
// instructions have completed locally, it can either use RDTSCP
// (if the processor supports that instruction) or execute the sequence
// LFENCE;RDTSC. This instruction was introduced by the Pentium processor.

// RDTSCP
// The RDTSCP instruction waits until all previous instructions have
// been executed before reading the counter. However, subsequent
// instructions may begin execution before the read operation is
// performed.

void setaffinity(int cpu);
void tic(tv_t *tv1);
double toc(tv_t tv1);
int reverse_compare_double(const void *a, const void *b);
void print_measure_cw(double *measure, int nb_tests);
int reverse_compare_uint64 (const void *a, const void *b);
void fprint_measure_tsc(FILE *fd, uint64_t *measure, int nb_tests, uint64_t offset);
void print_measure_tsc(uint64_t *measure, int nb_tests, uint64_t offset);

#endif

