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
#include <qd/qd_real.h>
#include "libdot.h"
#include "vec.h"
#include "eft.h"

inline void acc_extr4_vec(vec<double> vsig[4], vec<double> v, int k) {
  for(int j = k; j < 3; j++)
    vsig[j].fasttwosum_in(v, v);
  vsig[3].add(vsig[3], v);
}

qd_real dotprod_extr4_vec(const qd_real *x, const qd_real *y, int n) {
  // We will only use the new algorithm to deal with the first
  // m = n / vec<double>::len qd_reals in x and y. The rest of
  // the dot-product, and the reduction, will be computed using
  // non-vectorized algorithms.
  const double eps = 0x1.p-53;
  double *xx = (double *)x;
  double *yy = (double *)y;
  vec<double> x0, x1, x2, x3, y0, y1, y2, y3;
  vec<double> vT, vT0, vt1, vt2;
  vec<double> vp00, vp01, vp10, vp02, vp11, vp20;
  vec<double> vpi00, vpi01, vpi10, vpi02, vpi11, vpi20;  
  vec<double> vsig0[4], vsig[4];
  int i, j, m = n / vec<double>::len;

  // First pass through the data.
  // At each iteration, we access vec<double>::len consecutive qd_reals
  // (composed of 4 consecutive doubles): at each iteration the offset
  // i to the pointers xx and yy is incremented by 4*vec<double>::len.
  // In other words, the loop deals with the first n / vec<double>::len
  // qd_reals of the input vectors.
  vT0.set1(0.0);
  for(i = 0; i <= 4*(n-vec<double>::len); i += 4*vec<double>::len) {
    loadhi_qd_real(x0, xx+i);
    loadhi_qd_real(y0, yy+i);
    vt1.mul(x0, y0);
    vt1.abs(vt1);
    vT0.add(vT0, vt1);
  }

  // vT0 = t / (1 - m*eps);
  vt1.set1(1 - m * eps);
  vT0.div(vT0, vt1);
  // vsig[0] = (2*vT0) / (1 - (3*m+1) * eps);
  vt1.set1(2.0);
  vt1.mul(vt1, vT0);
  vt2.set1(1 - (3*m+1) * eps);
  vsig[0].div(vt1, vt2);

  // vT = (3/2 * m * vsig[0] + 3 * vT0) * eps;
  vt1.set1(1.5 * m);
  vt1.mul(vt1, vsig[0]);
  vt2.set1(3);
  vt2.mul(vt2, vT0);
  vt1.add(vt1, vt2);
  vt2.set1(eps);
  vT.mul(vt1, vt2);
  // vsig[1] = (2*vT) / (1 - (12*m+1) * eps);
  vt1.set1(2.0);
  vt1.mul(vt1, vT);
  vt2.set1(1 - (12*m+1) * eps);
  vsig[1].div(vt1, vt2);

  // vT = (6 * m * vsig[1] + 5 * vT0 * eps) * eps;
  vt1.set1(6.0 * m);
  vt1.mul(vt1, vsig[1]);
  vt2.set1(5.0 * eps);
  vt2.mul(vt2, vT0);
  vt1.add(vt1, vt2);
  vt2.set1(eps);
  vT.mul(vt1, vt2);
  // vsig[2] = (2*vT) / (1 - (27*m+1) * eps);
  vt1.set1(2.0);
  vt1.mul(vt1, vT);
  vt2.set1(1 - (27*m+1) * eps);
  vsig[2].div(vt1, vt2);

  vsig[3].set1(0.0);

  for(i = 0; i < 4; i++) vsig0[i] = vsig[i];
  
  // First pass through the data
  // Again, at each iteration the offset i to the pointers xx and yy
  // is incremented by 4*vec<double>::len.
  for(i = 0; i <= 4*(n-vec<double>::len); i += 4*vec<double>::len) {
    load_qd_real(x0, x1, x2, x3, xx+i);
    load_qd_real(y0, y1, y2, y3, yy+i);
    
    vp00.twoprod(vpi00, x0, y0);
    acc_extr4_vec(vsig, vp00, 0);
    acc_extr4_vec(vsig, vpi00, 1);

    vp01.twoprod(vpi01, x0, y1);
    acc_extr4_vec(vsig, vp01, 1);
    acc_extr4_vec(vsig, vpi01, 2);
    vp10.twoprod(vpi10, x1, y0);
    acc_extr4_vec(vsig, vp10, 1);
    acc_extr4_vec(vsig, vpi10, 2);
    
    vp02.twoprod(vpi02, x0, y2);
    acc_extr4_vec(vsig, vp02, 2);
    acc_extr4_vec(vsig, vpi02, 3);    
    vp11.twoprod(vpi11, x1, y1);
    acc_extr4_vec(vsig, vp11, 2);
    acc_extr4_vec(vsig, vpi11, 3);    
    vp20.twoprod(vpi20, x2, y0);
    acc_extr4_vec(vsig, vp20, 2);
    acc_extr4_vec(vsig, vpi20, 3);

    vt1.mul(x0, y3);
    vt2.mul(x1, y2);
    vt1.add(vt1, vt2);
    vt2.mul(x2, y1);
    vt1.add(vt1, vt2);
    vt2.mul(x3, y0);
    vt1.add(vt1, vt2);
    vsig[3].add(vsig[3], vt1);
  }

  // We have treated i/4 elements in each input vector,
  // so we are left with n-i/4 elements to deal with.
  qd_real r = dotprod_comp4(x+i/4, y+i/4, n-i/4);

  // Final reduction.
  double tsig[4][vec<double>::len], ts[4];
  
  for(i = 0; i < 4; i++) {
    vsig[i].sub(vsig[i], vsig0[i]);
    vsig[i].store(&tsig[i][0]);
  }

  for(j = 0; j < vec<double>::len; j++) {
    for(i = 0; i < 4; i++) ts[i] = tsig[i][j];
    RenormQd(ts);
    r += qd_real(ts);
  }

  return r;
}
