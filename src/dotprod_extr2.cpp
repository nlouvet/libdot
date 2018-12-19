/*

Created on Wed, 19 Dec 2018 12:08:01 +0000

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

#include <cstdio>
#include <qd/dd_real.h>
#include "vec.h"
#include "eft.h"

dd_real dotprod_extr2(const dd_real *x, const dd_real *y, int n)
{
  double p, pi, sh, sl, T0, sig0;
  const double eps = 0x1.p-53;
  int i;

  T0 = 0.0;
  for(i = n-1; i >= 0; i--)
    T0 += fabs(x[i]._hi() * y[i]._hi());
  T0 = T0 / (1 - n * eps);
  sig0 = (2.0 * T0) / (1 - (3*n+1) * eps);

  sh = sig0;
  sl = 0.0;
  for(i = 0; i < n; i++) {
    TwoProd(p, pi, x[i]._hi(), y[i]._hi());
    FastTwoSumIn(sh, p, p);
    sl += p + pi + x[i]._hi() * y[i]._lo() + x[i]._lo() * y[i]._hi();
  }
  sh -= sig0;
  
  TwoSumIn(sh, sl, sl);

  return dd_real(sh, sl);
}
