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

#include "libdot.h"
#include "vec.h"
#include "eft.h"

qd_real dotprod_extr4(const qd_real *x, const qd_real *y, int n) {
  const double eps = 0x1.p-53;
  double p_1, e_1, p_2, e_2, p_3, e_3, sigma[4];
  double x0, x1, x2, x3, y0, y1, y2, y3;
  double s0, s1, s2, s3;  
  double T0, T, t;
  int k;

  for(t = 0.0, k = n-1; k >= 0; k--) t += fabs(x[k][0]) * fabs(y[k][0]);

  T0 = t / (1 - n*eps);
  sigma[0] = (2*T0) / (1 - (3*n+1) * eps);

  T = (3/2 * n * sigma[0] + 3 * T0) * eps;
  sigma[1] = (2*T) / (1 - (12*n+1) * eps);

  T = (6 * n * sigma[1] + 5 * T0 * eps) * eps;
  sigma[2] = (2*T) / (1 - (27*n+1) * eps);
  
  s0 = sigma[0];
  s1 = sigma[1];
  s2 = sigma[2];
  s3 = 0.0;

  for(k = 0; k < n; k++) {
    x0 = x[k][0]; y0 = y[k][0];
    TwoProd(p_1, e_1, x0, y0);
    FastTwoSumIn(s0, p_1, p_1); FastTwoSumIn(s1, p_1, p_1); FastTwoSumIn(s2, p_1, p_1);
    FastTwoSumIn(s1, e_1, e_1); FastTwoSumIn(s2, e_1, e_1);
    s3 += p_1 + e_1;

    y1 = y[k][1]; x1 = x[k][1];
    TwoProd(p_1, e_1, x0, y1);
    TwoProd(p_2, e_2, x1, y0);
    FastTwoSumIn(s1, p_1, p_1); FastTwoSumIn(s2, p_1, p_1);
    FastTwoSumIn(s1, p_2, p_2); FastTwoSumIn(s2, p_2, p_2);
    FastTwoSumIn(s2, e_1, e_1);
    FastTwoSumIn(s2, e_2, e_2);
    s3 += e_1 + p_1;
    s3 += e_2 + p_2;
    
    y2 = y[k][2]; x2 = x[k][2];
    TwoProd(p_1, e_1, x0, y2); s3 += e_1;
    TwoProd(p_2, e_2, x1, y1); s3 += e_2;
    TwoProd(p_3, e_3, x2, y0); s3 += e_3;
    FastTwoSumIn(s2, p_1, p_1); s3 += p_1;
    FastTwoSumIn(s2, p_2, p_2); s3 += p_2;
    FastTwoSumIn(s2, p_3, p_3); s3 += p_3;
    
    x3 = x[k][3]; y3 = y[k][3];
    s3 += x0 * y3 + x1 * y2 + x2 * y1 + x3 * y0;
  }

  sigma[0] = s0 - sigma[0];
  sigma[1] = s1 - sigma[1];
  sigma[2] = s2 - sigma[2];
  sigma[3] = s3;
  
  RenormQd(sigma);

  return qd_real(sigma);
}

/*
inline void acc_extr4(double sigma[4], double v, int k) {
  for(int j = k; j < 3; j++)
    FastTwoSumIn(sigma[j], v, v);
  sigma[3] += v;
}

qd_real dotprod_extr4(const qd_real *x, const qd_real *y, int n) {
  const double eps = 0x1.p-53;
  double sigma0[4], sigma[4];
  double p00, p01, p10, p02, p11, p20;
  double pi00, pi01, pi10, pi02, pi11, pi20;
  double T0, T, t;
  int i;

  for(t = 0.0, i = n-1; i >= 0; i--) t += fabs(x[i][0]) * fabs(y[i][0]);

  T0 = t / (1 - n*eps);
  sigma[0] = (2*T0) / (1 - (3*n+1) * eps);

  T = (3/2 * n * sigma[0] + 3 * T0) * eps;
  sigma[1] = (2*T) / (1 - (12*n+1) * eps);

  T = (6 * n * sigma[1] + 5 * T0 * eps) * eps;
  sigma[2] = (2*T) / (1 - (27*n+1) * eps);
  
  sigma[3] = 0.0;

  for(i = 0; i < 4; i++) sigma0[i] = sigma[i];

  for(i = 0; i < n; i++) {
    TwoProd(p00, pi00, x[i][0], y[i][0]);
    acc_qdcmp2(sigma, p00, 0);
    acc_qdcmp2(sigma, pi00, 1);
    
    TwoProd(p01, pi01, x[i][0], y[i][1]);
    acc_qdcmp2(sigma, p01, 1);
    acc_qdcmp2(sigma, pi01, 2);
    TwoProd(p10, pi10, x[i][1], y[i][0]);
    acc_qdcmp2(sigma, p10, 1);
    acc_qdcmp2(sigma, pi10, 2);

    TwoProd(p02, pi02, x[i][0], y[i][2]);
    acc_qdcmp2(sigma, p02, 2);
    acc_qdcmp2(sigma, pi02, 3);
    TwoProd(p11, pi11, x[i][1], y[i][1]);
    acc_qdcmp2(sigma, p11, 2);
    acc_qdcmp2(sigma, pi11, 3);
    TwoProd(p20, pi20, x[i][2], y[i][0]);
    acc_qdcmp2(sigma, p20, 2);
    acc_qdcmp2(sigma, pi20, 3);
    
    sigma[3] += x[i][0] * y[i][3] + x[i][1] * y[i][2] + x[i][2] * y[i][1] + x[i][3] * y[i][0];
  }
  
  for(i = 0; i < 4; i++) sigma[i] -= sigma0[i];
  RenormQd(sigma);

  return qd_real(sigma);
}
*/
