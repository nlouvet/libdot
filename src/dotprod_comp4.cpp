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

#include "libdot.h"
#include "vec.h"
#include "eft.h"

qd_real dotprod_comp4(const qd_real *x, const qd_real *y, int n) {
  double p_1, e_1, p_2, e_2, p_3, e_3, s[4];
  double x0, x1, x2, x3;
  double y0, y1, y2, y3;
  double s0, s1, s2, s3;
  
  s0 = s1 = s2 = s3 = 0.0;

  for(int k = 0; k < n; k++) {
    // products x[k][i] * y[k][j], for i, j = 0..3 and i + j <= 3
    // i + j = 0
    x0 = x[k][0]; y0 = y[k][0];
    TwoProd(p_1, e_1, x0, y0);
    TwoSumIn(s0, p_1, p_1); TwoSumIn(s1, p_1, p_1); TwoSumIn(s2, p_1, p_1);
    TwoSumIn(s1, e_1, e_1); TwoSumIn(s2, e_1, e_1);
    s3 += p_1 + e_1;
    
    y1 = y[k][1]; x1 = x[k][1];
    TwoProd(p_1, e_1, x0, y1);
    TwoProd(p_2, e_2, x1, y0);
    TwoSumIn(s1, p_1, p_1); TwoSumIn(s2, p_1, p_1);
    TwoSumIn(s1, p_2, p_2); TwoSumIn(s2, p_2, p_2);
    TwoSumIn(s2, e_1, e_1);
    TwoSumIn(s2, e_2, e_2);
    s3 += e_1 + p_1;
    s3 += e_2 + p_2;
    
    // i + j = 2
    y2 = y[k][2]; x2 = x[k][2];
    TwoProd(p_1, e_1, x0, y2); s3 += e_1;
    TwoProd(p_2, e_2, x1, y1); s3 += e_2;
    TwoProd(p_3, e_3, x2, y0); s3 += e_3;
    TwoSumIn(s2, p_1, p_1); s3 += p_1;
    TwoSumIn(s2, p_2, p_2); s3 += p_2;
    TwoSumIn(s2, p_3, p_3); s3 += p_3;
    
    // i + j = 3
    x3 = x[k][3]; y3 = y[k][3];
    s3 += x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0;
  }

  s[0] = s0; s[1] = s1; s[2] = s2; s[3] = s3;

  RenormQd(s);

  return qd_real(s);
}

/*
qd_real dotprod_comp4(const qd_real *x, const qd_real *y, int n) {
  double p, e, s[4];
    
  for(int j = 0; j < 4; j++) s[j] = 0.0;
  for(int k = 0; k < n; k++) {
    // We consider the products x[k][i] * y[k][j], for i, j = 0..3 and i + j <= 3
    // i + j = 0
    TwoProd(p, e, x[k][0], y[k][0]);
    TwoSumIn(s[0], p, p); TwoSumIn(s[1], p, p); TwoSumIn(s[2], p, p);
    TwoSumIn(s[1], e, e); TwoSumIn(s[2], e, e);
    s[3] += (p + e);
    
    // i + j = 1
    TwoProd(p, e, x[k][0], y[k][1]);
    TwoSumIn(s[1], p, p); TwoSumIn(s[2], p, p);
    TwoSumIn(s[2], e, e);
    s[3] += (p + e);

    TwoProd(p, e, x[k][1], y[k][0]);
    TwoSumIn(s[1], p, p); TwoSumIn(s[2], p, p);
    TwoSumIn(s[2], e, e);
    s[3] += (p + e);
    
    // i + j = 2
    TwoProd(p, e, x[k][0], y[k][2]);
    TwoSumIn(s[2], p, p);
    s[3] += (p + e);

    TwoProd(p, e, x[k][1], y[k][1]);
    TwoSumIn(s[2], p, p);
    s[3] += (p + e);

    TwoProd(p, e, x[k][2], y[k][0]);
    TwoSumIn(s[2], p, p);
    s[3] += (p + e);
    
    // i + j = 3
    s[3] += (x[k][0] * y[k][3] + x[k][1] * y[k][2] + x[k][2] * y[k][1] + x[k][3] * y[k][0]);
  }

  RenormQd(s);

  return qd_real(s);
}
*/
