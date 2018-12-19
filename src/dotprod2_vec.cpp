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
#include "vec.h"
#include "eft.h"

double dotprod2_vec(const double *x, const double *y, int n) {
  vec<double> vs, ves, vp, vep, vc, v1, v2;
  int i;
  
  vs.set1(0.0);
  vc.set1(0.0);

  for(i=0; i<=n-vec<double>::len; i+=vec<double>::len) {
    v1.load(x+i);
    v2.load(y+i);
    vp.twoprod(vep, v1, v2);
    vc.add(vc, vep);
    vs.twosum_in(ves, vp);
    vc.add(vc, ves);
  }

  vec_align double ts[vec<double>::len];
  vec_align double tc[vec<double>::len];
  double p, ep, s, es, c;
  
  vs.store(ts);
  vc.store(tc);
  
  for(s=0.0, c=0.0; i<n; i++) {
    TwoProd(p, ep, x[i], y[i]);
    TwoSumIn(s, es, p);
    c += (es + ep);
  }
  
  for(i=0; i<vec<double>::len; i++) {
    TwoSumIn(s, es, ts[i]);
    c += (tc[i] + es);
  }

  return s+c;
}

