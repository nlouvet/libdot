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

#include "mpfr_tools.h"

/* ************************************************ */

double rand01(void) {
  return (double)rand() / (double)RAND_MAX;
}

long int rand_exp(double b) {
  return (long int)rint(b * rand01());
}

long int rand_idx(long int n) {
  return (long int)rint(n * (double)rand() / (double)(RAND_MAX+1.0));
}

void rand_unif01(mpfr_t x, gmp_randstate_t state) {
  mpfr_urandom(x, state, MPFR_RNDN);
}

void rand_pos(mpfr_t x, long int e, gmp_randstate_t state) {
  rand_unif01(x, state);
  mpfr_mul_2si(x, x, e, MPFR_RNDN);
}

void rand_sgn(mpfr_t x, gmp_randstate_t state) {
  rand_unif01(x, state);
  if(rand01() < 0.5) mpfr_neg(x, x, MPFR_RNDN);
}

void rand_sgn(mpfr_t x, long int e, gmp_randstate_t state) {
  rand_sgn(x, state);
  mpfr_mul_2si(x, x, e, MPFR_RNDN);
}


/* ************************************************ */

mpfr_t *mpfr_array_new(int n, mpfr_prec_t prec) {
  mpfr_t *a;

  a = (mpfr_t *)malloc(n * sizeof(mpfr_t));
  for(int i = 0; i < n; i++) mpfr_init2(a[i], prec);
  return(a);
}

void mpfr_array_clear(mpfr_t *a, int n) {
  for(int i=0; i<n; i++) mpfr_clear(a[i]);
  free(a);
}

void mpfr_array_abs(mpfr_t *r, mpfr_t *a, int n) {
  for(int i=0; i<n; i++) mpfr_abs(r[i], a[i], MPFR_RNDN);
}

void mpfr_dot(mpfr_t r, mpfr_t *x, mpfr_t *y, int n) {
  int i, p = mpfr_get_prec(x[0]);
  mpfr_t t;
  
  mpfr_init2(t, p);
  mpfr_mul(r, x[0], y[0], MPFR_RNDN);
  for(i = 1; i < n; i++) {
    mpfr_mul(t, x[i], y[i], MPFR_RNDN);
    mpfr_add(r, r, t, MPFR_RNDN);
  }
  
  mpfr_clear(t);
}

void mpfr_array_set(mpfr_t *d, mpfr_t *s, int n) {
  for(int i = 0; i < n; i++) mpfr_set(d[i], s[i], MPFR_RNDN);
}

/* double */

void mpfr_array_set_d(mpfr_t *a, double *v, int n) {
  for(int i = 0; i < n; i++) mpfr_set_d(a[i], v[i], MPFR_RNDN);
}

mpfr_t *mpfr_array_new_set_d(double *v, int n) {
  mpfr_t *a = mpfr_array_new(n, 53);
  mpfr_array_set_d(a, v, n);
  return a;
}

void mpfr_array_get_d(double *v, mpfr_t *a,  int n) {
  for(int i = 0; i < n; i++) v[i] = mpfr_get_d(a[i], MPFR_RNDN);
}

/* ******************* dd_real ******************** */

void mpfr_set_dd(mpfr_t r, dd_real x) {
  mpfr_set_d(r, x._hi(), MPFR_RNDN);
  mpfr_add_d(r, r, x._lo(), MPFR_RNDN);
}

dd_real mpfr_get_dd(mpfr_t r) {
  mpfr_t t;
  mpfr_init2(t, mpfr_get_prec(r));
  mpfr_set(t, r, MPFR_RNDN);
  double u = mpfr_get_d(t, MPFR_RNDN);
  mpfr_sub_d(t, t, u, MPFR_RNDN);
  double v = mpfr_get_d(t, MPFR_RNDN);
  mpfr_clear(t);
  
  return dd_real(u, v);
}

void mpfr_array_set_dd(mpfr_t *a, dd_real *v, int n) {
  for(int i = 0; i < n; i++) mpfr_set_dd(a[i], v[i]);
}

mpfr_t *mpfr_array_new_set_dd(dd_real *x, int n) {
  mpfr_t *a = mpfr_array_new(n, 106);
  mpfr_array_set_dd(a, x, n);
  return a;
}

void mpfr_array_get_dd(dd_real *v, mpfr_t *a,  int n) {
  for(int i = 0; i < n; i++) v[i] = mpfr_get_dd(a[i]);
}

/* ******************* qd_real ******************** */

void mpfr_set_qd(mpfr_t r, qd_real x) {
  mpfr_set_d(r, x[0], MPFR_RNDN);
  for(int i = 1; i < 4; i++) 
    mpfr_add_d(r, r, x[i], MPFR_RNDN);
}

qd_real mpfr_get_qd(mpfr_t x) {
  qd_real r = 0.0;
  double u;
  mpfr_t t;

  mpfr_init2(t, mpfr_get_prec(x));
  mpfr_set(t, x, MPFR_RNDN);
  for(int i = 0; i < 4; i++) {
    u = mpfr_get_d(t, MPFR_RNDN);
    mpfr_sub_d(t, t, u, MPFR_RNDN);
    r[i] = u;
  }
  mpfr_clear(t);
  
  return r;
}

void mpfr_array_set_qd(mpfr_t *a, qd_real *v, int n) {
  for(int i = 0; i < n; i++) mpfr_set_qd(a[i], v[i]);
}

mpfr_t *mpfr_array_new_set_qd(qd_real *x, int n) {
  mpfr_t *a = mpfr_array_new(n, 212);
  mpfr_array_set_qd(a, x, n);
  return a;
}

void mpfr_array_get_qd(qd_real *v, mpfr_t *a,  int n) {
  for(int i = 0; i < n; i++) v[i] = mpfr_get_qd(a[i]);
}
