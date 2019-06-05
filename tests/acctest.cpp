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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <mpfr.h>
#include <time.h>
#include <unistd.h>
#include "gendot.h"
#include "libdot.h"
#include "mpfr_tools.h"

#define N 997
#define K 1
#define P 128

#define RUN_D  1
#define RUN_D2 2
#define RUN_DD 3
#define RUN_QD 4
#define RUN_EX 5

void print_usage(char *str) {
  fprintf(stderr, "Usage: %s [-r d2|dd|qd] [-n n] [-k k]\n", str);
  fprintf(stderr, "  -r experiments to be run (d2|dd|qd|ex)\n");
  fprintf(stderr, "  -n length of the vectors used\n");
  fprintf(stderr, "  -k number of dot product generated for each conditionning\n");
}

int main(int argc, char **argv) {
  gmp_randstate_t state;
  int run = RUN_D2, n = N, k = K, p = P;
  
  extern char *optarg; 
  extern int optind, opterr; 
  int opt, errflag = 0;

  while((opt = getopt(argc, argv, "hr:n:k:p:")) != -1) {
    switch(opt) {
    case 'h':
      print_usage(argv[0]);
      exit(EXIT_SUCCESS);
    case 'r':
      if(strcmp(optarg, "d") == 0) run = RUN_D;
      else if(strcmp(optarg, "d2") == 0) run = RUN_D2;
      else if(strcmp(optarg, "dd") == 0) run = RUN_DD;
      else if(strcmp(optarg, "qd") == 0) run = RUN_QD;
      else if(strcmp(optarg, "ex") == 0) run = RUN_EX;
      else {
	fprintf(stderr, "wrong argument to -%c\n", opt);
	errflag = 1;
      }
      break;
    case 'n':
      if(sscanf(optarg, "%d", &n) != 1) {
	fprintf(stderr, "wrong argument to -%c\n", opt);
	errflag = 1;
      }
      break;
    case 'k':
      if(sscanf(optarg, "%d", &k) != 1) {
	fprintf(stderr, "wrong argument to -%c\n", opt);
	errflag = 1;
      }
      break;
    case 'p':
      if(optind >= argc) {
	fprintf(stderr, "expected argument after -%c\n", opt);
	errflag = 1;
	break;
      }
      if(sscanf(optarg, "%d", &k) != 1) {
	fprintf(stderr, "wrong argument to -%c\n", opt);
	errflag = 1;
      }
      break;
    case '?':  
      fprintf(stderr, "Unknown option: %c\n", optopt); 
    break;  
    }  
  }  

  if(errflag) {
    fprintf(stderr, "Wrong arguments.");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  //srand(time(NULL));
  srand(0);
  gmp_randinit_default(state);
  gmp_randseed_ui(state, rand());

  if(run == RUN_D) {
    
    int numfcts = 3, i, j, k;
    const char *namefct[] = {"d", "d_vec", "d_pw"};
    double (* const ptrfct[])(const double *, const double *, int) = {dotprod, dotprod_vec, dotprod_pw};
    double *x = NULL, *y = NULL, zh, c, e;
    
    if(posix_memalign((void **)&x, libdot_align(), N*sizeof(double))) fprintf(stderr, "Memory allocation problem.\n");
    if(posix_memalign((void **)&y, libdot_align(), N*sizeof(double))) fprintf(stderr, "Memory allocation problem.\n");
    
    printf("# \"Working precision...\"\n");
    printf("# %10s ", "cond");
    for(j=0; j<numfcts; j++) printf("%12s ", namefct[j]);
    printf("\n");
    
    for(i=0; i<20; i++) {
      for(k = 0; k < 2; k++) {
        c = GenDot2_double(x, y, N, pow(10, i), state);
        printf("%12.5e ", c);
        for(j=0; j<numfcts; j++) {
          zh = ptrfct[j](x, y, N);
          e = DotErr2_double(x, y, N, zh);
          printf("%12.5e ", e);
        }
        printf("\n");
      }
    }
    
    free(x);
    free(y);
    
  }
  
  else if(run == RUN_D2) {
    
    int numfcts = 5, i, j, k;
    const char *namefct[] = {"d", "d_vec", "d2", "d2_vec", "d2_dd"};
    double (* const ptrfct[])(const double *, const double *, int) = {dotprod, dotprod_vec, dotprod2, dotprod2_vec, dotprod2_dd};
    double *x = NULL, *y = NULL, zh, c, e;
    
    if(posix_memalign((void **)&x, libdot_align(), N*sizeof(double))) fprintf(stderr, "Memory allocation problem.\n");
    if(posix_memalign((void **)&y, libdot_align(), N*sizeof(double))) fprintf(stderr, "Memory allocation problem.\n");

    printf("# \"twice the working precision...\"\n");
    printf("# %10s ", "cond");
    for(j=0; j<numfcts; j++)
      printf("%12s ", namefct[j]);
    printf("\n");
    
    for(i=0; i<32; i++) {
      for(k = 0; k < 2; k++) {
        c = GenDot2_double(x, y, N, pow(10, i), state);
        printf("%12.5e ", c);
        for(j=0; j<numfcts; j++) {
          zh = ptrfct[j](x, y, N);
          e = DotErr2_double(x, y, N, zh);
          printf("%12.5e ", e);
        }
        printf("\n");
      }
    }
    
    free(x);
    free(y);
    
  }
  else if(run == RUN_DD) {
    
    mpfr_t *mpfx, *mpfy, mpfz;
    dd_real *x, *y, xh;
    int i, j;
  
    x = y = NULL;
    if(posix_memalign((void **)&x, libdot_align(), n*sizeof(dd_real))) fprintf(stderr, "Memory allocation problem.\n");
    if(posix_memalign((void **)&y, libdot_align(), n*sizeof(dd_real))) fprintf(stderr, "Memory allocation problem.\n");
    mpfx = mpfr_array_new(n, 106);
    mpfy = mpfr_array_new(n, 106);
    mpfr_init2(mpfz, 106);

    printf("# K = 2\n");
    printf("# n = %d\n", n);
    printf("# col 1: cond\n");
    printf("# col 2: mpfr_dot\n");
    printf("# col 3: dd\n");
    printf("# col 4: comp2\n");
    printf("# col 5: extr2\n");
    printf("# col 6: dd_vec\n");
    printf("# col 7: comp2_vec\n");
    printf("# col 8: extr2_vec\n");
   
    for(i = 0; i < 40; i++) {
      for(j = 0; j < k; j++) {
        printf("%12.5e ", GenDot2_dd_real(x, y, n, rand01() * pow(10, i), state));
        mpfr_array_set_dd(mpfx, x, n);
        mpfr_array_set_dd(mpfy, y, n);    

        mpfr_dot(mpfz, mpfx, mpfy, n);
        printf("%12.5e ", DotErr2_mpfr(mpfx, mpfy, n, mpfz));
	
        xh = dotprod_dd(x, y, n);
    	printf("%12.5e ", DotErr2_dd_real(x, y, n, xh));

        xh = dotprod_comp2(x, y, n);
        printf("%12.5e ", DotErr2_dd_real(x, y, n, xh));

        xh = dotprod_extr2(x, y, n);
        printf("%12.5e ", DotErr2_dd_real(x, y, n, xh));
	
        xh = dotprod_dd_vec(x, y, n);
        printf("%12.5e ", DotErr2_dd_real(x, y, n, xh));
	
        xh = dotprod_comp2_vec(x, y, n);
        printf("%12.5e ", DotErr2_dd_real(x, y, n, xh));
	
        xh = dotprod_extr2_vec(x, y, n);
        printf("%12.5e\n", DotErr2_dd_real(x, y, n, xh));
      }
    }
  
    free(x);
    free(y);
    mpfr_array_clear(mpfx, n);
    mpfr_array_clear(mpfy, n);
    mpfr_clear(mpfz);

  }
  else if(run == RUN_QD) {
    
    mpfr_t *mpfx, *mpfy, mpfz;
    qd_real *x, *y, xh;
    int i, j;
  
    x = y = NULL;
    if(posix_memalign((void **)&x, libdot_align(), n*sizeof(qd_real))) fprintf(stderr, "Memory allocation problem.\n");
    if(posix_memalign((void **)&y, libdot_align(), n*sizeof(qd_real))) fprintf(stderr, "Memory allocation problem.\n");
    mpfx = mpfr_array_new(n, 212);
    mpfy = mpfr_array_new(n, 212);
    mpfr_init2(mpfz, 212);

    printf("# K = 4\n");
    printf("# n = %d\n", n);
    printf("# col 1: cond\n");
    printf("# col 2: qd\n");
    printf("# col 3: comp4\n");
    printf("# col 4: extr4\n");
    printf("# col 5: comp4_vec\n");
    printf("# col 6: extr4_vec\n");
   
    for(i = 2; i < 72; i++) {
      for(j = 0; j < k; j++) {
        printf("%12.5e ", GenDot2_qd_real(x, y, n, pow(10, i), state));
        mpfr_array_set_qd(mpfx, x, n);
        mpfr_array_set_qd(mpfy, y, n);    

        xh = dotprod_qd(x, y, n);
        printf("%12.5e ", DotErr2_qd_real(x, y, n, xh));

        xh = dotprod_comp4(x, y, n);
        printf("%12.5e ", DotErr2_qd_real(x, y, n, xh));

        xh = dotprod_extr4(x, y, n);
        printf("%12.5e ", DotErr2_qd_real(x, y, n, xh));
	
        xh = dotprod_comp4_vec(x, y, n);
        printf("%12.5e ", DotErr2_qd_real(x, y, n, xh));

        xh = dotprod_extr4_vec(x, y, n);
        printf("%12.5e\n", DotErr2_qd_real(x, y, n, xh));
      }
    }
  
    free(x);
    free(y);
    mpfr_array_clear(mpfx, n);
    mpfr_array_clear(mpfy, n);
    mpfr_clear(mpfz);

  }
  else if(run = RUN_EX) {
    
    mpfr_t *mpfx, *mpfy, mpfz, *mpfx512, *mpfy512;
    int i, j;

    mpfx = mpfr_array_new(n, p);
    mpfy = mpfr_array_new(n, p);
    mpfr_init2(mpfz, p);
    mpfx512 = mpfr_array_new(n, 512);
    mpfy512 = mpfr_array_new(n, 512);

    for(i = 0; i < 80; i++) {
      for(j = 0; j < k; j++) {
        printf("%12.5e ", GenDot2_mpfr(mpfx, mpfy, n, rand01() * pow(10, i), state));
        mpfr_array_set(mpfx512, mpfx, n);
        mpfr_array_set(mpfy512, mpfy, n);    

	mpfr_dot(mpfz, mpfx, mpfy, n);
        printf("%12.5e\n", DotErr2_mpfr(mpfx, mpfy, n, mpfz));
      }
    }
  
    mpfr_array_clear(mpfx, n);
    mpfr_array_clear(mpfy, n);
    mpfr_array_clear(mpfx512, n);
    mpfr_array_clear(mpfy512, n);
    mpfr_clear(mpfz);

  }
  else 
    fprintf(stderr, "This should not have happened...\n");
  
  gmp_randclear(state);
        
  return 0;
}

