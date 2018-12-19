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
#include <qd/qd_real.h>
#include "libdot.h"
#include "vec.h"
#include "eft.h"

qd_real dotprod_comp4_vec(const qd_real *x, const qd_real *y, int n) {
  double *xx = (double *)x;
  double *yy = (double *)y;
  vec<double> x0, x1, x2, x3, y0, y1, y2, y3, p, e, s0, s1, s2, s3, t;
  int i, j;

  // At each iteration, we treat vec<double>::len consecutive qd_reals
  // (composed of 4 consecutive doubles) in each input vectors: at each
  // iteration the offset i to the pointers xx and yy is incremented by
  // 4*vec<double>::len.
  s0.set1(0.0);
  s1.set1(0.0);
  s2.set1(0.0);
  s3.set1(0.0);
  for(i = 0; i <= 4*(n-vec<double>::len); i += 4*vec<double>::len) {
    load_qd_real(x0, x1, x2, x3, xx+i);
    load_qd_real(y0, y1, y2, y3, yy+i);

    p.twoprod(e, x0, y0);
    s0.twosum_in(p, p); s1.twosum_in(p, p); s2.twosum_in(p, p);
    s1.twosum_in(e, e); s2.twosum_in(e, e);
    t.add(p, e); s3.add_in(t);
    
    p.twoprod(e, x0, y1);
    s1.twosum_in(p, p); s2.twosum_in(p, p);
    s2.twosum_in(e, e);
    t.add(p, e); s3.add_in(t);

    p.twoprod(e, x1, y0);
    s1.twosum_in(p, p); s2.twosum_in(p, p);
    s2.twosum_in(e, e);
    t.add(p, e); s3.add_in(t);

    p.twoprod(e, x0, y2);
    s2.twosum_in(p, p);
    t.add(p, e); s3.add_in(t);

    p.twoprod(e, x1, y1);
    s2.twosum_in(p, p);
    t.add(p, e); s3.add_in(t);
    
    p.twoprod(e, x2, y0);
    s2.twosum_in(p, p);
    t.add(p, e); s3.add_in(t);

    p.mul(x0, y3); s3.add_in(p);
    p.mul(x1, y2); s3.add_in(p);
    p.mul(x2, y1); s3.add_in(p);
    p.mul(x3, y0); s3.add_in(p);
  }

  // So far, we have treated i/4 elements in each input
  // vector, so we are left with n-i/4 elements to deal with.
  qd_real r = dotprod_comp4(x+i/4, y+i/4, n-i/4);

  // Final reduction.
  vec_align double tval[4][vec<double>::len], tv[4];
  
  s0.store(tval[0]);
  s1.store(tval[1]);
  s2.store(tval[2]);
  s3.store(tval[3]);
  for(i = 0; i < vec<double>::len; i++) {
    for(j = 0; j < 4; j++) tv[j] = tval[j][i];
    RenormQd(tv);
    r += qd_real(tv);
  }
  return r;
}
