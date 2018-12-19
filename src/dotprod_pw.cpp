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
along with the libdot library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include <math.h>
#include "vec.h"

double dotprod_pw(const double *x, const double *y, int n) {
  unsigned int len = n, m = n & ~(2 * vec<double>::len - 1), i;
  vec_align vec<double> stack[64];
  vec<double> u1, v1, u2, v2, w;
  int p = 0;
  
  for(i = 0; i < m; i += 2*vec<double>::len) {
    u1.load(x+i);
    v1.load(y+i);
    u1.mul(u1, v1);
        
    u2.load(x+i+vec<double>::len);
    v2.load(y+i+vec<double>::len);
    u2.mul(u2, v2);
    
    w.add(u1, u2);
    for(int b = 2*vec<double>::len; i & b; b <<= 1, --p) w.add(w, stack[p-1]);
    stack[p++] = w;
  }
  
  vec<double> vsum;
  vsum.set1(0.0);
  while(p) vsum.add(vsum, stack[--p]);
  double sum = vsum.add_red();
  for(i = m; i < len; i++) sum += x[i] * y[i];
  
  return sum;
}

double dotprod_pw0(double *x, double *y, int n) {
  double stack[64], v;
  int p = 0;
  
  for(int i = 0; i < n; ++i) {
    v = x[i] * y[i];
    for(int b = 1; i & b; b <<= 1, --p)
      v += stack[p-1];
    stack[p++] = v;
  }
  double s = 0.0;
  while(p) s += stack[--p];
  return s;
}

/*
double dotprod_pw(double *x, double *y, int n) {
  double stack[64], v, w;
  int p = 0;
  
  for(int i = 0; i < n; i += 4) {
    v = x[i+0] * y[i+0] + x[i+1] * y[i+1];
    w = x[i+2] * y[i+2] + x[i+3] * y[i+3]; 
    v+= w;
    for(int b = 4; i & b; b <<= 1, --p)
      v += stack[p-1];
    stack[p++] = v;
  }
  double s = 0.0;
  int epilogue = n & ~3ul;
  for(int i = epilogue; i < n; ++i) s += x[i] * y[i];
  while(p) s += stack[--p];
  return s;
}
*/

/*
double dotprod_pw(double *x, double *y, int n) {
  unsigned int m = n & ~(2 * vec<double>::len - 1);
  vec_align vec<double> stack[16];
  vec<double> u1, v1, u2, v2, w;
  int i, p = 0;
  
  for(i = 0; i < m; i += 2*vec<double>::len) {
    u1.load(x+i);
    v1.load(y+i);
    u1.mul(u1, v1);
        
    u2.load(x+i+vec<double>::len);
    v2.load(y+i+vec<double>::len);
    u2.mul(u2, v2);
    
    w.add(u1, u2);
    for(int b = 2*vec<double>::len; i & b; b <<= 1, --p) w.add(w, stack[p-1]);
    stack[p++] = w;
  }
  
  vec<double> vsum;
  vsum.set1(0.0);
  while(p) vsum.add(vsum, stack[--p]);
  
  vec_align double t[vec<double>::len];
  vsum.store(t);
  
  double sum = 0.0;
  for(i = m; i < n; i++) sum += x[i] * y[i];
  for(i = 0; i < vec<double>::len; i++) sum += t[i];
  
  return sum;
}
*/






