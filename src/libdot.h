/*

Created on Wed, 19 Dec 2018 12:08:02 +0000

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

#if !defined(__LIBDOT__)
#define __LIBDOT__

#include <cstdlib>
#include <qd/dd_real.h>
#include <qd/qd_real.h>

size_t libdot_align(void);

double dotprod(const double *x, const double *y, int len);

double dotprod_vec(const double *x, const double *y, int len);
double dotprod2(const double *x, const double *y, int len);
double dotprod2_vec(const double *x, const double *y, int len);
double dotprod2_dd(const double *x, const double *y, int len);

dd_real dotprod_dd(const dd_real *x, const dd_real *y, int n);
dd_real dotprod_dd_vec(const dd_real *x, const dd_real *y, int len);

dd_real dotprod_comp2(const dd_real *x, const dd_real *y, int len);
dd_real dotprod_extr2(const dd_real *x, const dd_real *y, int len);

dd_real dotprod_comp2_vec(const dd_real *x, const dd_real *y, int len);
dd_real dotprod_extr2_vec(const dd_real *x, const dd_real *y, int n);

qd_real dotprod_qd(const qd_real *x, const qd_real *y, int n);

qd_real dotprod_comp4(const qd_real *x, const qd_real *y, int len);
qd_real dotprod_extr4(const qd_real *x, const qd_real *y, int len);

qd_real dotprod_comp4_vec(const qd_real *x, const qd_real *y, int n);
qd_real dotprod_extr4_vec(const qd_real *x, const qd_real *y, int n);

#endif
