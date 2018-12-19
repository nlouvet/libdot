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

#ifndef __GEN2__
#define __GEN2__

#include <stdlib.h>
#include <math.h>
#include "mpfr.h"
#include "gmp.h"
#include <qd/dd_real.h>
#include <qd/qd_real.h>

double GenDot2_mpfr(mpfr_t *x, mpfr_t *y, int n, double c, gmp_randstate_t state);
double DotErr2_mpfr(mpfr_t *u, mpfr_t *v, int n, mpfr_t sh);

double GenDot2_double(double *x, double *y, int n, double c, gmp_randstate_t state);
double DotErr2_double(double *x, double *y, int n, double sh);

double GenDot2_dd_real(dd_real *x, dd_real *y, int n, double c, gmp_randstate_t state);
double DotErr2_dd_real(dd_real *x, dd_real *y, int n, dd_real sh);

double GenDot2_qd_real(qd_real *x, qd_real *y, int n, double c, gmp_randstate_t state);
double DotErr2_qd_real(qd_real *x, qd_real *y, int n, qd_real sh);

void Dot(mpfr_t r, mpfr_t *x, mpfr_t *y, int n);

#endif

