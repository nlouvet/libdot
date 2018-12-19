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

#include <cstdio>
#include "vec.h"

double dotprod_vec(const double *x, const double *y, int n) {
  vec<double> vmm0, vmm1, vmm2;
  int i;
  
  vmm0.set1(0.0);
  
  for(i=0; i<n-vec<double>::len; i+=vec<double>::len) {
    vmm1.load(x+i);
    vmm2.load(y+i);
    vmm0.muladd(vmm1, vmm2, vmm0);
  }

  vec_align double t[vec<double>::len];
  double s;

  vmm0.store(t);
  
  for(s=0; i<n; i++) s += x[i] * y[i];
  
  for(i=0; i<vec<double>::len; i++) s += t[i];
  
  return s;
}
