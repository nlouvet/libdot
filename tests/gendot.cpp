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

#include "gendot.h"
#include "mpfr_tools.h"

// *********************************************************************

void Dot(mpfr_t r, mpfr_t *x, mpfr_t *y, int n) {
  int i, p = mpfr_get_prec(x[0]);
  
  mpfr_ptr *t = (mpfr_ptr *)malloc(n*sizeof(mpfr_ptr));
  
  for(i=0; i<n; i++) {
    t[i] = (mpfr_ptr)malloc(sizeof(mpfr_t));
    mpfr_init2(t[i], 2*p);
    mpfr_mul(t[i], x[i], y[i], MPFR_RNDN);
  }
  mpfr_sum(r, t, n, MPFR_RNDN);
  for(i=0; i<n; i++) {
    mpfr_clear(t[i]);
    free(t[i]);
  }
  free(t);
}

double DotErr2_mpfr(mpfr_t *u, mpfr_t *v, int n, mpfr_t sh) {
  mpfr_prec_t p = mpfr_get_prec(sh);
  mpfr_t s, t;
  double e;

  mpfr_init2(s, 4*p);
  mpfr_init2(t, 4*p);

  Dot(s, u, v, n);
  mpfr_sub(t, s, sh, MPFR_RNDN);
  e = fabs(mpfr_get_d(t, MPFR_RNDN) / mpfr_get_d(s, MPFR_RNDN));

  //if( e >= 2.0 ) e = 2.0;

  mpfr_clear(s);
  mpfr_clear(t);
  
  return(e);
}

double CondDot2_mpfr(mpfr_t *u, mpfr_t *v, int n) {
  mpfr_prec_t p = mpfr_get_prec(u[0]);
  mpfr_t s, a;
  double c;

  mpfr_init2(s, 4*p);
  mpfr_init2(a, p);

  Dot(s, u, v, n);

  mpfr_array_abs(u, u, n);
  mpfr_array_abs(v, v, n);
  Dot(a, u, v, n);
  
  c = fabs(mpfr_get_d(a, MPFR_RNDN) / mpfr_get_d(s, MPFR_RNDN));
    
  return c;
}

double GenDot2_mpfr(mpfr_t *x, mpfr_t *y, int n, double c, gmp_randstate_t state) {
  int i, j, k, e, n2 = n/2;
  double b;
  mpfr_prec_t p = mpfr_get_prec(x[0]);
  mpfr_t s, t;

  mpfr_init2(s, 4*p);
  mpfr_init2(t, p);

  mpfr_set_ui(s, 0, MPFR_RNDN);
 
  b = log2(c);
  
  rand_sgn(x[0], state);
  rand_sgn(y[0], state);
  mpfr_fma(s, x[0], y[0], s, MPFR_RNDN);


  rand_sgn(x[1], round(b / 2) + 1, state);
  rand_sgn(y[1], round(b / 2) + 1, state);
  mpfr_fma(s, x[1], y[1], s, MPFR_RNDN);

  for(i = 2; i < n2; i++) {
    rand_sgn(x[i], round(rand01() * b / 2), state);
    rand_sgn(y[i], round(rand01() * b / 2), state);
    mpfr_fma(s, x[i], y[i], s, MPFR_RNDN);
  }
  
  for(i = n2; i < n; i++) {
    e = round(b - (i-n2) * b / (n-n2-1));
    rand_sgn(x[i], round(rand01() * e), state);
    rand_sgn(y[i], round(rand01() * e), state);
    mpfr_set(t, s, MPFR_RNDN);
    mpfr_sub(t, y[i], t, MPFR_RNDN);
    mpfr_div(y[i], t, x[i], MPFR_RNDN);
    mpfr_fma(s, x[i], y[i], s, MPFR_RNDN);
  }

  for(i = 1; i < n/2; i++) {
    j = rand_idx(n-1);
    while((k = rand_idx(n-1)) == j);
    mpfr_swap (x[j], x[k]);
    mpfr_swap (y[j], y[k]);
  }

  Dot(t, x, y, n);
  for(i=0, c=0.0; i<n; i++)
    c += fabs(mpfr_get_d(x[i], MPFR_RNDN) * mpfr_get_d(y[i], MPFR_RNDN));
  c /= 2 * fabs(mpfr_get_d(t, MPFR_RNDN));
  
  mpfr_clear(s);
  mpfr_clear(t);
  
  return c;
}

// *********************************************************************

double DotErr2_double(double *x, double *y, int n, double sh) {
  MPFR_DECL_INIT(t, 53);
  mpfr_set_d(t, sh, MPFR_RNDN);
  mpfr_t *u, *v;
  double e;

  u = mpfr_array_new_set_d(x, n);
  v = mpfr_array_new_set_d(y, n);
  e = DotErr2_mpfr(u, v, n, t);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  return(e);
}

double CondDot2_double(double *x, double *y, int n) {
  mpfr_t *u, *v;
  double c;

  u = mpfr_array_new_set_d(x, n);
  v = mpfr_array_new_set_d(y, n);
  c = CondDot2_mpfr(u, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  return c;
}

double GenDot2_double(double *x, double *y, int n, double c, gmp_randstate_t state) {
  mpfr_t *u, *v;
  
  u = mpfr_array_new(n, 53);
  v = mpfr_array_new(n, 53);
  c = GenDot2_mpfr(u, v, n, c, state);
  mpfr_array_get_d(x, u, n);
  mpfr_array_get_d(y, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  return CondDot2_double(x, y, n);
}

// *********************************************************************

double DotErr2_dd_real(dd_real *x, dd_real *y, int n, dd_real sh) {
  MPFR_DECL_INIT(t, 106);
  mpfr_set_dd(t, sh);
  mpfr_t *u, *v;
  double e;

  u = mpfr_array_new_set_dd(x, n);
  v = mpfr_array_new_set_dd(y, n);
  e = DotErr2_mpfr(u, v, n, t);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  return e;
}

double CondDot2_dd_real(dd_real *x, dd_real *y, int n) {
  mpfr_t *u, *v;
  double c;

  u = mpfr_array_new_set_dd(x, n);
  v = mpfr_array_new_set_dd(y, n);
  c = CondDot2_mpfr(u, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  return c;
}

double GenDot2_dd_real(dd_real *x, dd_real *y, int n, double c, gmp_randstate_t state) {
  mpfr_t *u, *v;
  
  u = mpfr_array_new(n, 106);
  v = mpfr_array_new(n, 106);
  c = GenDot2_mpfr(u, v, n, c, state);
  mpfr_array_get_dd(x, u, n);
  mpfr_array_get_dd(y, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  
  return CondDot2_dd_real(x, y, n);
}

// *********************************************************************

double DotErr2_qd_real(qd_real *x, qd_real *y, int n, qd_real sh) {
  MPFR_DECL_INIT(t, 212);
  mpfr_set_qd(t, sh);
  mpfr_t *u, *v;
  double e;

  u = mpfr_array_new_set_qd(x, n);
  v = mpfr_array_new_set_qd(y, n);
  e = DotErr2_mpfr(u, v, n, t);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  
  return e;
}

double CondDot2_qd_real(qd_real *x, qd_real *y, int n) {
  mpfr_t *u, *v;
  double c;

  u = mpfr_array_new_set_qd(x, n);
  v = mpfr_array_new_set_qd(y, n);
  c = CondDot2_mpfr(u, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  
  return c;
}

double GenDot2_qd_real(qd_real *x, qd_real *y, int n, double c, gmp_randstate_t state) {
  mpfr_t *u, *v;
  
  u = mpfr_array_new(n, 212);
  v = mpfr_array_new(n, 212);
  c = GenDot2_mpfr(u, v, n, c, state);
  mpfr_array_get_qd(x, u, n);
  mpfr_array_get_qd(y, v, n);
  mpfr_array_clear(u, n);
  mpfr_array_clear(v, n);
  
  return CondDot2_qd_real(x, y, n);
}
