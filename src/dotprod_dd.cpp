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
along with the libdot library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include <qd/dd_real.h>
#include "eft.h"

dd_real dotprod_dd(const dd_real *x, const dd_real *y, int n) {
  double xh, xl, yh, yl, sh, sl, ph, pl, e;

  sh = sl = 0.0;
  for(int i=0; i<n; i++) {
    xh = x[i]._hi(); xl = x[i]._lo();
    yh = y[i]._hi(); yl = y[i]._lo();
    // (ph+pl) <- (xh+xl) * (yh+yl)
    TwoProd(ph, pl, xh, yh);
    pl += xh*yl + xl*yh;
    FastTwoSumIn(ph, pl, pl);
    // (sh+sl) <- (sh+sl) + (ph+pl)
    TwoSumIn(sh, e, ph);
    e += (sl + pl);
    FastTwoSumIn(sh, sl, e);
  }
  return dd_real(sh, sl);
}

/*
dd_real dotprod_dd(const dd_real *x, const dd_real *y, int n) {
  dd_real r(0.0);

  for(int i=0; i<n; i++)
    r += x[i] * y[i];
  return r;
}
*/
