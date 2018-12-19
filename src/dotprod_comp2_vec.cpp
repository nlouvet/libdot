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

#include <cstdio>
#include <qd/dd_real.h>
#include "vec.h"
#include "eft.h"

dd_real dotprod_comp2_vec(const dd_real *x, const dd_real *y, int n) {
  double *xx = (double *)x;
  double *yy = (double *)y;
  vec<double> x1, x2, y1, y2, xh, yh, xl, yl, p, sh, sl, e1, e2, t1, t2, t3;
  int i, j;

  sh.set1(0.0);
  sl.set1(0.0);

  for(i = 0; i <= 2*(n-vec<double>::len); i += 2*vec<double>::len) {
    x1.load(xx + i);
    x2.load(xx + i + vec<double>::len);
    xh.unpacklo(x1, x2);
    xl.unpackhi(x1, x2);

    y1.load(yy + i);
    y2.load(yy + i + vec<double>::len);
    yh.unpacklo(y1, y2);
    yl.unpackhi(y1, y2);

    p.twoprod(e1, xh, yh);
    sh.twosum_in(e2, p);
    t1.muladd(xh, yl, e1);
    t2.muladd(xl, yh, e2);
    t3.add(t1, t2);
    sl.add(sl, t3);
    
  }
    
  vec_align double tsh[vec<double>::len];
  vec_align double tsl[vec<double>::len];
  double dsh(0.0), dsl(0.0), dp1, de1, de;

  sh.store(tsh);
  sl.store(tsl);

  for(j=i/2; j<n; j++) {
    TwoProd(dp1, de1, x[j]._hi(), y[j]._hi());
    TwoSumIn(dsh, de, dp1);
    dsl += ((x[j]._hi())*(y[j]._lo()) + (x[j]._lo())*(y[j]._hi()) + de1 + de);
  }

  for(i = 0; i < vec<double>::len; i++) {
    TwoSumIn(dsh, de, tsh[i]);
    dsl += (tsl[i] + de);
  }

  TwoSumIn(dsh, dsl, dsl);

  return dd_real(dsh, dsl);
}

