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

#include <cstdio>
#include <qd/dd_real.h>
#include "libdot.h"
#include "vec.h"
#include "eft.h"

dd_real dotprod_extr2_vec(const dd_real *x, const dd_real *y, int n) {
  const double eps = 0x1.p-53;
  double *xx = (double *)x;
  double *yy = (double *)y;
  vec<double> x1, x2, y1, y2, xh, yh, xl, yl;
  vec<double> vT0, vsig0, vsh, vsl, vp, ve, vpi, vp1, vp2, vt;
  int i, m = n / vec<double>::len;
  
  // first pass, to determine T0 and sig0
  vT0.set1(0.0);
  for(i = 0; i <= 2*(n-vec<double>::len); i += 2*vec<double>::len) {
    x1.load(xx + i);
    x2.load(xx + i + vec<double>::len);
    xh.unpacklo(x1, x2);
  
    y1.load(yy + i);
    y2.load(yy + i + vec<double>::len);
    yh.unpacklo(y1, y2);

    vp.mul(xh, yh);
    vp.abs(vp);
    vT0.add(vT0, vp);
  }
  
  vt.set1(1 - m * eps);
  vT0.div(vT0, vt);
  
  vsig0.set1(2.0);
  vsig0.mul(vsig0, vT0);
  vt.set1(1 - (3*m+1) * eps);
  vsig0.div(vsig0, vt);
  
  // second pass through the vector...
  vsh = vsig0;
  //vsh.set1(0.0);
  vsl.set1(0.0);
  for(i = 0; i <= 2*(n-vec<double>::len); i += 2*vec<double>::len) {
    x1.load(xx + i);
    x2.load(xx + i + vec<double>::len);
    xh.unpacklo(x1, x2);
    xl.unpackhi(x1, x2);

    y1.load(yy + i);
    y2.load(yy + i + vec<double>::len);
    yh.unpacklo(y1, y2);
    yl.unpackhi(y1, y2);
    
    vp.twoprod(vpi, xh, yh);
    vsh.fasttwosum_in(ve, vp);

    vp1.muladd(xh, yl, ve);
    vp2.muladd(xl, yh, vpi);
    vsl.add(vsl, vp1);
    vsl.add(vsl, vp2);
  }
  vsh.sub(vsh, vsig0);

  // final reduction...
  double dsh(0.0), dsl(0.0), dp, de, df;
  vec_align double tsh[vec<double>::len];
  vec_align double tsl[vec<double>::len];
  vsh.store(tsh);
  vsl.store(tsl);

  for(int j = i/2; j < n; j++) {
    TwoProd(dp, de, x[j]._hi(), y[j]._hi());
    TwoSumIn(dsh, df, dp);
    dsl += (x[j]._hi())*(y[j]._lo()) + (x[j]._lo())*(y[j]._hi()) + de + df;
  }

  for(i = 0; i < vec<double>::len; i++) {
    TwoSumIn(dsh, de, tsh[i]);
    dsl += tsl[i] + de;
  }
  
  TwoSumIn(dsh, dsl, dsl);

  return dd_real(dsh, dsl);
}
