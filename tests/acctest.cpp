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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <mpfr.h>
#include <time.h>
#include "gendot.h"
#include "libdot.h"
#include "mpfr_tools.h"

#define N 997
#define K 1
#define P 128

int main(int argc, char **argv) {
  gmp_randstate_t state;
  int n = N, k = K, p = P;
  
  //srand(time(NULL));
  srand(0);
  gmp_randinit_default(state);
  gmp_randseed_ui(state, rand());

  if( (argc == 1) || (strcmp(argv[1], "d2") == 0) ) {
    if( (argc == 3) && (sscanf(argv[2], "%d", &n) != 1)) n = N;

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
  else if( (argc >= 2) && (strcmp(argv[1], "dd") == 0)) {
    if( (argc >= 3) && (sscanf(argv[2], "%d", &n) != 1)) n = N;
    if( (argc >= 4) && (sscanf(argv[3], "%d", &k) != 1)) k = K;

    
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
  else if( (argc > 1) && (strcmp(argv[1], "qd") == 0) ) {
    if( (argc >= 3) && (sscanf(argv[2], "%d", &n) != 1)) n = N;
    if( (argc >= 4) && (sscanf(argv[3], "%d", &k) != 1)) k = K;
    
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
  else if( (argc >= 2) && (strcmp(argv[1], "exp") == 0)) {
    if( (argc >= 3) && (sscanf(argv[2], "%d", &n) != 1)) n = N;
    if( (argc >= 4) && (sscanf(argv[3], "%d", &k) != 1)) k = K;
    if( (argc >= 5) && (sscanf(argv[4], "%d", &p) != 1)) p = P;
    
    
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
    printf("Wrong arguments...\n");
  
  gmp_randclear(state);
        
  return 0;
}

