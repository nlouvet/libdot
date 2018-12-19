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

// 2^(-53) = 1/9007199254740992 = 0x1.0p-53

void print_result(const char *msg, int res) {
  printf("- %s: ", msg);
  if(res == 0) printf("failed\n");
  else printf("passed\n");
}
// ------------------------------------------------------------
int test_rtn1(void) {
  double ta[] = {1.0,       1.0      };
  double tb[] = {0x1.0p-54, 0x1.0p-53};
  double tr[] = {1.0,       1.0      };
  double a, b;
  int i, cnt = 0;

  for(i=0; i<2; i++) {
    a = ta[i]; b = tb[i];
    if(a + b != tr[i]) cnt++;
  }

  // if cnt > 0, then the test has failed, and 0 is returned
  if(cnt > 0) return(0);
  else return(1);
}
// ------------------------------------------------------------
int test_reasso1() {
  double a = 1.0, b = 0x1.0p-53;
  double x = (a + b) - a;
  // the expected result is 0.0
  if(x == 0.0) return(1);
  else return(0);
}
// ------------------------------------------------------------
double test_reasso2_aux(double a, double b) {
  return ((a + b) - a);
}

int test_reasso2() {
  double a = 1.0, b = 0x1.0p-53;
  // the expected result is 0.0
  if(test_reasso2_aux(a, b) == 0.0) return(1);
  else return(0);
}
// ------------------------------------------------------------
int test_reasso3_aux(double a[], double b[], double y[], int n) {
  double x;
  int i;

  for(i=0; i<n; i++) {
    x = (a[i] + b[i]) - a[i];
    if(x != y[i]) return(0);
  }
  return (1);
}

int test_reasso3() {
  double a[] = {1.0};
  double b[] = {0x1.0p-53};
  double y[] = {0.0};
  return test_reasso3_aux(a, b, y, 1);

}
// ------------------------------------------------------------
int main(void) {
  int res, cnt = 0;

  printf("* checking whether 'rounding to the nearest with ties to even' is the default rounding mode:\n");
  res = test_rtn1(); print_result("test 1", res); cnt += res;
  printf("* checking whether re-association of operands in floating-point expressions is disabled:\n");
  res = test_reasso1(); print_result("test 1", res); cnt += res;
  res = test_reasso2(); print_result("test 2", res); cnt += res;
  res = test_reasso3(); print_result("test 3", res); cnt += res;

  if(cnt == 4) return(0);
  else return(1);
}

