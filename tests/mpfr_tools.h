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

#ifndef __MPFR_TOOLS__
#define __MPFR_TOOLS__

#include <stdlib.h>
#include <mpfr.h>
#include <qd/dd_real.h>
#include <qd/qd_real.h>

double rand01(void);
long int rand_exp(double b);
long int rand_idx(long int n);

void rand_unif01(mpfr_t x, gmp_randstate_t state);
void rand_pos(mpfr_t x, long int e, gmp_randstate_t state);
void rand_sgn(mpfr_t x, gmp_randstate_t state);
void rand_sgn(mpfr_t x, long int e, gmp_randstate_t state);

mpfr_t *mpfr_array_new(int n, mpfr_prec_t prec);
void mpfr_array_clear(mpfr_t *a, int n);
void mpfr_array_abs(mpfr_t *r, mpfr_t *a, int n);
void mpfr_array_set(mpfr_t *d, mpfr_t *s, int n);
void mpfr_dot(mpfr_t r, mpfr_t *x, mpfr_t *y, int n);

/* with double */

void mpfr_array_set_d(mpfr_t *a, double *v, int n);
mpfr_t *mpfr_array_new_set_d(double *v, int n);
void mpfr_array_get_d(double *v, mpfr_t *a,  int n);

/* with dd_real */

void mpfr_set_dd(mpfr_t r, dd_real x);
dd_real mpfr_get_dd(mpfr_t r);
void mpfr_array_set_dd(mpfr_t *a, dd_real *v, int n);
mpfr_t *mpfr_array_new_set_dd(dd_real *x, int n);
void mpfr_array_get_dd(dd_real *v, mpfr_t *a,  int n);

/* with dd_real */

void mpfr_set_qd(mpfr_t r, qd_real x);
qd_real mpfr_get_qd(mpfr_t r);
void mpfr_array_set_qd(mpfr_t *a, qd_real *v, int n);
mpfr_t *mpfr_array_new_set_qd(qd_real *x, int n);
void mpfr_array_get_qd(qd_real *v, mpfr_t *a,  int n);

#endif /* __MPFR_TOOLS__ */

