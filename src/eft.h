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
along with the hplll Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

*/

///////////////
#ifndef __EFT__
#define __EFT__
///////////////

/* (s, e) = TwoSum(a, b) */
inline void TwoSum(double &s, double &e, double a, double b) {
  s = a + b;
  double t0 = s - a;
  double t1 = s - t0;
  double t2 = b - t0;
  double t3 = a - t1;
  e = t2 + t3;
}

/* (s, e) = TwoSum(s, b) */
inline void TwoSumIn(double &s, double &e, double b) {
  double ts = s + b;
  double t0 = ts - s;
  double t1 = ts - t0;
  double t2 = b - t0;
  double t3 = s - t1;
  e = t2 + t3;
  s = ts;
}

/* (s, e) = TwoSum(a, b), assuming |a| >= |b| */
inline void FastTwoSum(double &s, double &e, double a, double b) {
  double t;
  s = a + b;
  t = a - s;
  e = t + b;
}

/* (s, e) = FastTwoSum(s, b), assuming |a| >= |b| */
inline void FastTwoSumIn(double &s, double &e, double b) {
  double r = s + b;
  double t = s - r;
  e = t + b ;
  s = r;
}

//////////////////////////
#if defined(__FP_FAST_FMA) // See "Common Predefined Macros" in GNU CPP documentation
//////////////////////////

/* (p, e) = TwoProd(a, b) */
inline void TwoProd(double &p, double &e, double a, double b) {
  p = a * b;
  e = fma(a, b, -p);
}

/////
#else
/////

/* (p, e) = TwoProd(a, b) */
inline void TwoProd(double &p, double &e, double a, double b) {
  const static double splitter = 134217729.0; /*2^27+1*/
  double ashift = a * splitter; double bshift = b * splitter;
  double t1 = ashift - a;       double t2 = bshift - b;
  double ahi = ashift - t1;     double bhi = bshift - t2;
  double alo = a - ahi;         double blo = b - bhi;
  p = a * b;
  e = (((ahi*bhi - p) + ahi*blo) + alo*bhi) + alo*blo;
}

//////
#endif
//////

// Renormalization of a quad-double expansion. Algorithm 25 in V. Popescu's
// PhD thesis (Towards fast and certified multiple-precision libraries, 2017).
inline void RenormQd(double x[4]) {
  double e, t;
  int i, j;

  for(i = 1; i < 4; i++) {
    t = x[i];
    for(j = i-1; j >= 0; j--)
      TwoSumIn(t, x[j+1], x[j]);
    x[0] = t;
  }
  e = x[0];
  j = 0;
  for(i = 1; i < 4; i++) {
    TwoSum(t, e, x[i], e);
    if(e == 0.0) {
      e = t;
    }
    else {
      x[j] = t;
      j++;
    }
  }
  if(e != 0.0) x[j] = e;
}

/////////////////
#endif // __EFT__
/////////////////
