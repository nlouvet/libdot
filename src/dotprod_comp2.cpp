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
#include "eft.h"

dd_real dotprod_comp2(const dd_real *x, const dd_real *y, int n) {
  double xh, xl, yh, yl, s, c, p1, e1, e2;
  int i;
  
  s = c = 0.0; 
  for(i=0; i<n; i++) {
    xh = x[i]._hi();
    xl = x[i]._lo();
    
    yh = y[i]._hi();
    yl = y[i]._lo();
    
    TwoProd(p1, e1, xh, yh);
    TwoSumIn(s, e2, p1);
    c += (xh*yl + xl*yh + e1 + e2);
  }
  
  TwoSumIn(s, c, c);
  
  return dd_real(s, c);
}

