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
along with the hplll library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include <cstdio>
#include <qd/dd_real.h>
#include "vec.h"
#include "eft.h"

// xx and yy allow us to look at x and y as arrays of doubles:
// x = [{x0h, x0l}, {x1h, x1l}, {x2h, x2l}, {x3h, x3l}, ...]
// xx = [x0h, x0l, x1h, x1l, x2h, x2l, x3h, x3l, ...]
// At each iteration of the first for loop, we load 2*vec<double>::len
// elements both from xx and from yy, to form 4 simd registers xh, xl, yh, yl
// of size vec<double>::len.

dd_real dotprod_dd_vec(const dd_real *x, const dd_real *y, int n) {
  double *xx = (double *)x;
  double *yy = (double *)y;
  vec<double> x1, x2, y1, y2, xh, yh, xl, yl, ph, pl, sh, sl, t1, t2;
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
    /*
    ddmul(ph, pl, xh, xl, yh, yl);
    ddadd_in(sh, sl, ph, pl);
    */
    
    //ddmul(ph, pl, xh, xl, yh, yl);
    ph.twoprod(pl, xh, yh);
    t1.mul(xh, yl);
    pl.add(pl, t1);
    t2.mul(xl, yh);
    pl.add(pl, t2);
    ph.fasttwosum_in(pl, pl);
    //ddadd_in(sh, sl, ph, pl);
    sh.twosum(t1, sh, ph);
    t1.add(t1, sl);
    t1.add(t1, pl);
    sh.fasttwosum_in(sl, t1);
  }
  vec_align double tsh[vec<double>::len];
  vec_align double tsl[vec<double>::len];
  dd_real s(0.0);

  sh.store(tsh);
  sl.store(tsl);
  
  for(j=i/2; j<n; j++) s += x[j] * y[j];
  
  for(i=0; i<vec<double>::len; i++) s += dd_real(tsh[i], tsl[i]);

  return dd_real(s);
}

