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
#ifndef __VEC__
#define __VEC__
///////////////

#include <math.h>
#include <immintrin.h>

template<class T>
class vec {
};

#define _LDST(CT, BL, IT) \
inline void load(const CT *p) { r = _mm ## BL ## _loadu_ ## IT(p); }   \
inline void store(CT *p) { _mm ## BL ## _storeu_ ## IT (p, r); }

#define _UNPACK(CT, BL, IT) \
inline void unpackhi(vec<CT> v1, vec<CT> v2) {   \
  r = _mm ## BL ## _unpackhi_ ## IT(v1.r, v2.r); \
}                                                \
                                                 \
inline void unpacklo(vec<CT> v1, vec<CT> v2) {   \
  r = _mm ## BL ## _unpacklo_ ## IT(v1.r, v2.r); \
}

#define _SET1(CT, BL, IT) \
inline void set1(CT x) { r = _mm ## BL ## _set1_ ## IT(x); }

#define _ARITH(CT, BL, IT) \
inline void add(vec<CT> v1, vec<CT> v2) { r = _mm ## BL ## _add_ ## IT(v1.r, v2.r); } \
inline void add_in(vec<CT> v1)          { r = _mm ## BL ## _add_ ## IT(r, v1.r); }    \
inline void sub(vec<CT> v1, vec<CT> v2) { r = _mm ## BL ## _sub_ ## IT(v1.r, v2.r); } \
inline void mul(vec<CT> v1, vec<CT> v2) { r = _mm ## BL ## _mul_ ## IT(v1.r, v2.r); } \
inline void div(vec<CT> v1, vec<CT> v2) { r = _mm ## BL ## _div_ ## IT(v1.r, v2.r); }

#define _FMA(CT, BL, IT) \
inline void muladd(vec<CT> v1, vec<CT> v2, vec<CT> v3) { r = _mm ## BL ## _fmadd_ ## IT(v1.r, v2.r, v3.r); } \
inline void mulsub(vec<CT> v1, vec<CT> v2, vec<CT> v3) { r = _mm ## BL ## _fmsub_ ## IT(v1.r, v2.r, v3.r); }

#define _ABS(CT, BL, IT, MT)                                     \
inline void abs(vec<CT> &v) {                                    \
    static const MT sgnmask  = _mm ## BL ## _set1_ ## IT(-0.0f); \
    r = _mm ## BL ## _andnot_ ## IT(sgnmask, v.r);               \
}

#define _TWOPROD_DEKKER(CT, BL, IT, MT) \
inline void twoprod(vec<CT> &e, vec<CT> a, vec<CT> b) {         \
  MT splitter;                                                  \
  MT ashift, bshift;                                            \
  MT ah, al;                                                    \
  MT bh, bl;                                                    \
  MT t1, t2;                                                    \
  splitter = _mm ## BL ## _set1_ ## IT(134217729.0);            \
  ashift = _mm ## BL ## _mul_ ## IT(a.r, splitter);  bshift = _mm ## BL ## _mul_ ## IT(b.r, splitter); \
  t1 = _mm ## BL ## _sub_ ## IT(ashift, a.r);        t2 = _mm ## BL ## _sub_ ## IT(bshift, b.r);       \
  ah = _mm ## BL ## _sub_ ## IT(ashift, t1);         bh = _mm ## BL ## _sub_ ## IT(bshift, t2);        \
  al = _mm ## BL ## _sub_ ## IT(a.r, ah);            bl = _mm ## BL ## _sub_ ## IT(b.r, bh);           \
  r =  _mm ## BL ## _mul_ ## IT(a.r, b.r);                                                             \
  e.r = _mm ## BL ## _add_ ## IT(_mm ## BL ## _add_ ## IT(_mm ## BL ## _add_ ## IT(_mm ## BL ## _sub_ ## IT(_mm ## BL ## _mul_ ## IT(ah, bh), r), _mm ## BL ## _mul_ ## IT(ah, bl)), _mm ## BL ## _mul_ ## IT(al, bh)), _mm ## BL ## _mul_ ## IT(al, bl)); \
}

#define _TWOPROD_FMA(CT, BL, IT, MT)                                                      \
inline void twoprod(vec<CT> &ve, vec<CT> v1, vec<CT> v2) {                                \
  r = _mm ## BL ## _mul_ ## IT(v1.r, v2.r);         /* r = fl(v1.r * v2.r) */             \
  ve.r = _mm ## BL ## _fmsub_ ## IT(v1.r, v2.r, r); /* r + ve.r = v1.r * v2.r exactely */ \
}

#define _TWOSUM_KNUTH(CT, BL, IT, MT)                     \
inline void twosum(vec<CT> &ve, vec<CT> v1, vec<CT> v2) { \
  MT t0, t1, t2, t3;                                      \
  r    = _mm ## BL ## _add_ ## IT(v1.r, v2.r);            \
  t0   = _mm ## BL ## _sub_ ## IT(r, v1.r);               \
  t1   = _mm ## BL ## _sub_ ## IT(r, t0);                 \
  t2   = _mm ## BL ## _sub_ ## IT(v2.r, t0);              \
  t3   = _mm ## BL ## _sub_ ## IT(v1.r, t1);              \
  ve.r = _mm ## BL ## _add_ ## IT(t2, t3);                \
}                                                         \
                                                          \
inline void twosum_in(vec<CT> &ve, vec<CT> v1) {          \
  MT s, t0, t1, t2, t3;                                   \
  s    = _mm ## BL ## _add_ ## IT(r, v1.r);               \
  t0   = _mm ## BL ## _sub_ ## IT(s, r);                  \
  t1   = _mm ## BL ## _sub_ ## IT(s, t0);                 \
  t2   = _mm ## BL ## _sub_ ## IT(v1.r, t0);              \
  t3   = _mm ## BL ## _sub_ ## IT(r, t1);                 \
  ve.r = _mm ## BL ## _add_ ## IT(t2, t3);                \
  r = s;                                                  \
}

#define _FASTTWOSUM(CT, BL, IT, MT)                                          \
/* Assume |v1.r| >= |v2.r| */                                                \
inline void fasttwosum(vec<CT> &ve, vec<CT> v1, vec<CT> v2) {                \
  r = _mm ## BL ## _add_ ## IT(v1.r, v2.r);                                  \
  ve.r = _mm ## BL ## _add_ ## IT(_mm ## BL ## _sub_ ## IT (v1.r, r), v2.r); \
}                                                                            \
                                                                             \
/* Assume |r| >= |v2.r| */                                                   \
inline void fasttwosum_in(vec<CT> &ve, vec<CT> v2) {                         \
  MT s;                                                                      \
  s = _mm ## BL ## _add_ ## IT(r, v2.r);                                     \
  ve.r = _mm ## BL ## _add_ ## IT(_mm ## BL ## _sub_ ## IT (r, s), v2.r);    \
  r = s;                                                                     \
}

////////////////////////
#if defined(__AVX512F__)
////////////////////////
#define VEC_ALIGN_CST 64

template<> class vec<double> {
  public:
  __m512d r;
  static const int len = 8;
  _LDST(double, 512, pd)
  _UNPACK(double, 512, pd)
  _SET1(double, 512, pd)
  _ARITH(double, 512, pd)
  _FMA(double, 512, pd)
  _TWOPROD_FMA(double, 512, pd, __m512d)
  _FASTTWOSUM(double, 512, pd, __m512d)

  //
  inline void abs() {
    
  }

  inline void abs(vec<double> &v) {
    static const __m512i sgnmask = _mm512_set1_epi64(0x7fffffffffffffffLL);
    r = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v.r), sgnmask));
  }
  
  // (this, ve) <- TwoProd(v1, v2) componentwise.
  // We must have &v1 != &this (but we may have &v2 == &ve).
  // Note that the castings (__m512i) and (__m512d) may be replaced
  // by _mm512_castpd_si512 and _mm512_castsi512_pd respectively; but
  // I prefer keeping the casting since it seems to have a slight
  // positive impact on the running time...
  inline void twosum(vec<double> &ve, vec<double> v1, vec<double> v2) {
    r = _mm512_add_pd(v1.r, v2.r);
    static const __m512i sgnmask = _mm512_set1_epi64(0x7fffffffffffffffLL);
    __m512d abs1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v1.r), sgnmask));
    __m512d abs2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v2.r), sgnmask));
    __mmask8 mask = _mm512_cmp_pd_mask(abs1, abs2, _CMP_GE_OS);
    __m512d maxabs = _mm512_mask_blend_pd(mask, v2.r, v1.r);
    __m512d minabs = _mm512_mask_blend_pd(mask, v1.r, v2.r);
    ve.r = _mm512_add_pd(_mm512_sub_pd(maxabs, r), minabs);
  }
  
  inline void twosum_in(vec<double> &ve, vec<double> v1) {
    __m512d s = _mm512_add_pd(r, v1.r);
    static const __m512i sgnmask = _mm512_set1_epi64(0x7fffffffffffffffLL);
    __m512d abs1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(r), sgnmask));
    __m512d abs2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(v1.r), sgnmask));
    __mmask8 mask = _mm512_cmp_pd_mask(abs1, abs2, _CMP_GE_OS);
    __m512d maxabs = _mm512_mask_blend_pd(mask, v1.r, r);
    __m512d minabs = _mm512_mask_blend_pd(mask, r, v1.r);
    r = s;
    ve.r = _mm512_add_pd(_mm512_sub_pd(maxabs, s), minabs);
  }
  
  inline double add_red(void) {
    // [r0, r1, r2, r3, r4, r5, r6, r7] <- r
    // [r0, r1, r2, r3] <- _mm512_castpd512_pd256(r)
    // [r4, r5, r6, r7] <- _mm512_extractf64x4_pd(r, 1)
    // [r0+r4, r1+r5, r2+r6, r3+r7] <- u
    // [r0+r1+r4+r5, r0+r1+r4+r5, r2+r3+r6+r7, r2+r3+r6+r7] <- v
    // [r2+r3+r6+r7, r2+r3+r6+r7, r0+r1+r4+r5, r0+r1+r4+r5] <- _mm256_permute2f128_pd(t, t, 0x01)
    // [r0+ ... +r7, r0+ ... +r7, r0+ ... +r7, r0+ ... +r7] <- w
    __m256d u = _mm256_add_pd(_mm512_castpd512_pd256(r), _mm512_extractf64x4_pd(r, 1));
    __m256d v = _mm256_hadd_pd(u, u);
    __m256d w = _mm256_add_pd(v, _mm256_permute2f128_pd(v, v, 0x01));
    return _mm_cvtsd_f64(_mm256_castpd256_pd128(w));
  }
};

inline void load_qd_real(vec<double> &v0, vec<double> &v1, vec<double> &v2, vec<double> &v3, double const *addr) {
  static const __m256i vidx0 = _mm256_set_epi32( 0, 32, 64,  96, 128, 160, 192, 224);
  static const __m256i vidx1 = _mm256_set_epi32( 8, 40, 72, 104, 136, 168, 200, 232);
  static const __m256i vidx2 = _mm256_set_epi32(16, 48, 80, 112, 144, 176, 208, 240);
  static const __m256i vidx3 = _mm256_set_epi32(24, 56, 88, 120, 152, 184, 216, 248);
  v0.r = _mm512_i32gather_pd(vidx0, addr, 1);
  v1.r = _mm512_i32gather_pd(vidx1, addr, 1);
  v2.r = _mm512_i32gather_pd(vidx2, addr, 1);
  v3.r = _mm512_i32gather_pd(vidx3, addr, 1);
}

inline void loadhi_qd_real(vec<double> &v0, double const *addr) {
  static const __m256i vidx0 = _mm256_set_epi32( 0, 32, 64,  96, 128, 160, 192, 224);
  v0.r = _mm512_i32gather_pd(vidx0, addr, 1);
}

///////////////////////
#elif defined(__AVX2__)
///////////////////////
#define VEC_ALIGN_CST 32

template<> class vec<double> {
  public:
  __m256d r;
  static const int len = 4;
  _LDST(double, 256, pd)
  _UNPACK(double, 256, pd)
  _SET1(double, 256, pd)
  _ARITH(double, 256, pd)
  _ABS(double, 256, pd, __m256d)  
  _FMA(double, 256, pd)
  _TWOPROD_FMA(double, 256, pd, __m256d)
  _FASTTWOSUM(double, 256, pd, __m256d)

  // (this, ve) <- TwoProd(v1, v2) componentwise.
  // We must have &v1 != &this (but we may have &v2 == &ve).
  inline void twosum(vec<double> &ve, vec<double> v1, vec<double> v2) {
    r = _mm256_add_pd(v1.r, v2.r);
    __m256d vmin = _mm256_min_pd(v1.r, v2.r);
    __m256d vmax = _mm256_max_pd(v1.r, v2.r);
    __m256d vminmag = _mm256_blendv_pd(vmin, vmax, r);
    __m256d vmaxmag = _mm256_blendv_pd(vmax, vmin, r);
    __m256d vt = _mm256_sub_pd(r, vmaxmag);
    ve.r = _mm256_sub_pd(vminmag, vt);
  }

  // (this, ve) <- TwoProd(this, v1) componentwise.
  inline void twosum_in(vec<double> &ve, vec<double> v1) {
    __m256d vmin = _mm256_min_pd(r, v1.r);
    __m256d vmax = _mm256_max_pd(r, v1.r);
    r = _mm256_add_pd(r, v1.r);
    __m256d vminmag = _mm256_blendv_pd(vmin, vmax, r);
    __m256d vmaxmag = _mm256_blendv_pd(vmax, vmin, r);
    __m256d vt = _mm256_sub_pd(r, vmaxmag);
    ve.r = _mm256_sub_pd(vminmag, vt);
  }
  
  inline double add_red(void) {
    // Assume r = [r0, r1, r2, r3]:
    __m256d u = _mm256_hadd_pd(r, r);                // u = [r0+r1, r0+r1, r2+r3, r2+r3]
    __m256d v = _mm256_permute2f128_pd(u, u, 0x01);  // v = [r2+r3, r2+r3, r0+r1, r0+r1]
    __m256d w = _mm256_add_pd(u, v);                 // w = [r0+r1+r2+r3, ... ]
    return _mm_cvtsd_f64(_mm256_castpd256_pd128(w)); // Convert w0 in w = [w0, w1, w2, w3] to a double
  }

  /* 
  // Previous version of the code for twoprod and twoprod_in:
  // this code works, but the current version is a bit faster...

  inline void twosum(vec<double> &ve, vec<double> v1, vec<double> v2) {
    const __m256d sgnmask  = _mm256_set1_pd(-0.0f);
    __m256d mask, maxabs, minabs;
    
    r = _mm256_add_pd(v1.r, v2.r);
    mask = _mm256_cmp_pd(_mm256_andnot_pd(sgnmask, v1.r), _mm256_andnot_pd(sgnmask, v2.r), _CMP_GE_OS);
    maxabs = _mm256_blendv_pd(v2.r, v1.r, mask);
    minabs = _mm256_blendv_pd(v1.r, v2.r, mask);
    ve.r = _mm256_add_pd(_mm256_sub_pd(maxabs, r), minabs);
  }

  inline void twosum_in(vec<double> &ve, vec<double> v1) {
    const __m256d sgnmask  = _mm256_set1_pd(-0.0f);
    __m256d s, mask, maxabs, minabs;
    
    s = _mm256_add_pd(r, v1.r);
    mask = _mm256_cmp_pd(_mm256_andnot_pd(sgnmask, r), _mm256_andnot_pd(sgnmask, v1.r), _CMP_GE_OS);
    maxabs = _mm256_blendv_pd(v1.r, r, mask);
    minabs = _mm256_blendv_pd(r, v1.r, mask);
    r = s;
    ve.r = _mm256_add_pd(_mm256_sub_pd(maxabs, s), minabs);
  }
  */
};

// Load 4 consecutive qd_real starting from address addr into 4 vec<double> v0, v1, v2, v3.
// Memory: x[0], ..., x[3], y[0], ..., y[3], z[0], ..., z[3], t[0], ..., t[t3]
// v0 <- [x[0], y[0], z[0], t[0]]
// ...
// v3 <- [x[3], y[3], z[3], t[3]]
inline void load_qd_real(vec<double> &v0, vec<double> &v1, vec<double> &v2, vec<double> &v3, double const *addr) {
  static const __m128i vidx0 = _mm_set_epi32( 0, 32, 64,  96);
  static const __m128i vidx1 = _mm_set_epi32( 8, 40, 72, 104);
  static const __m128i vidx2 = _mm_set_epi32(16, 48, 80, 112);
  static const __m128i vidx3 = _mm_set_epi32(24, 56, 88, 120);
  v0.r = _mm256_i32gather_pd(addr, vidx0, 1);
  v1.r = _mm256_i32gather_pd(addr, vidx1, 1);
  v2.r = _mm256_i32gather_pd(addr, vidx2, 1);
  v3.r = _mm256_i32gather_pd(addr, vidx3, 1);
}

// Load the high order part of 4 consecutive qd_real starting from address addr into vec<double> v0.
// Memory: x[0], ..., x[3], y[0], ..., y[3], z[0], ..., z[3], t[0], ..., t[t3]
// v0 <- [x[0], y[0], z[0], t[0]]
inline void loadhi_qd_real(vec<double> &v0, double const *addr) {
  static const __m128i vidx0 = _mm_set_epi32( 0, 32, 64,  96);
  v0.r = _mm256_i32gather_pd(addr, vidx0, 1);
}


/* Seems that this kind of loading procedure does not improve the speed..
inline void load_qd_real_(vec<double> &v0, vec<double> &v1, vec<double> &v2, vec<double> &v3, double const *addr) {
  __m256d x0 = _mm256_load_pd(addr);
  __m256d x1 = _mm256_load_pd(addr+4);
  __m256d x2 = _mm256_load_pd(addr+8);
  __m256d x3 = _mm256_load_pd(addr+12);

  __m256d y0 = _mm256_permute2f128_pd(x0, x1, 0x20);
  __m256d y1 = _mm256_permute2f128_pd(x2, x3, 0x20);
  __m256d y2 = _mm256_permute2f128_pd(x0, x1, 0x31);
  __m256d y3 = _mm256_permute2f128_pd(x2, x3, 0x31);
  
  v0.r = _mm256_unpacklo_pd(y0, y1);
  v1.r = _mm256_unpackhi_pd(y0, y1);
  v2.r = _mm256_unpacklo_pd(y2, y3);
  v3.r = _mm256_unpackhi_pd(y2, y3);
}
*/

/////////////////////////////////////////////////////////////////////
#elif defined(__SSE2__) && defined(__SSE4_1__) && defined(__SSE4_2__)
/////////////////////////////////////////////////////////////////////

#define VEC_ALIGN_CST 16

template<> class vec<double> {
  public:
  __m128d r;
  static const int len = 2;
  _LDST(double, , pd)
  _UNPACK(double, , pd)
  _SET1(double, , pd)
  _ARITH(double, , pd)
  _ABS(double, , pd, __m128d)  
  _TWOPROD_DEKKER(double, , pd, __m128d)  
  _FASTTWOSUM(double, , pd, __m128d)

  // (this, ve) <- TwoProd(v1, v2) componentwise.
  // We must have &v1 != &this (but we may have &v2 == &ve).
  inline void twosum(vec<double> &ve, vec<double> v1, vec<double> v2) {
    static const __m128d sign_mask  = _mm_set1_pd(-0.0f);
    __m128d mask, maxabs, minabs;
    r = _mm_add_pd(v1.r, v2.r);
    mask = _mm_cmpge_pd(_mm_andnot_pd(sign_mask, v1.r), _mm_andnot_pd(sign_mask,v2.r));
    maxabs = _mm_blendv_pd(v2.r, v1.r, mask);
    minabs = _mm_blendv_pd(v1.r, v2.r, mask);
    ve.r = _mm_add_pd(_mm_sub_pd(maxabs, r), minabs);
  }

  // (this, ve) <- TwoProd(this, v1) componentwise.
  inline void twosum_in(vec<double> &ve, vec<double> v1) {
    static const __m128d sign_mask  = _mm_set1_pd(-0.0f);
    __m128d mask, maxabs, minabs;
    mask = _mm_cmpge_pd(_mm_andnot_pd(sign_mask, r), _mm_andnot_pd(sign_mask, v1.r));
    maxabs = _mm_blendv_pd(v1.r, r, mask); // Available with sse4
    minabs = _mm_blendv_pd(r, v1.r, mask);
    r = _mm_add_pd(r, v1.r);
    ve.r = _mm_add_pd(_mm_sub_pd(maxabs, r), minabs);
  }

  inline double add_red(void) {
    // Assume r = [r0, r1]:
    __m128d u = _mm_hadd_pd(r, r); // u = [r0+r1, r0+r1]
    return _mm_cvtsd_f64(u);       // Convert u0 in u = [u0, u1] to a double
  }

  // no FMA available, so we replace it by an unfused MA...
  inline void muladd(vec<double> v1, vec<double> v2, vec<double> v3) {
    r = _mm_mul_pd(v1.r, v2.r);
    r = _mm_add_pd(r, v3.r);
  }

  inline void mulsub(vec<double> v1, vec<double> v2, vec<double> v3) {
    r = _mm_mul_pd(v1.r, v2.r);
    r = _mm_sub_pd(r, v3.r);
  }
  
  /*
  // Previous version of the code for twoprod and twoprod_in:
  // since the instruction blendvpd is available in SSE4.1,
  // it seems reasonnable to use an implementation base on
  // this instruction (both implementations seem to run
  // equally fast).

  #define _mm_select_pd(t, x, y) (_mm_or_pd(_mm_and_pd(t, x), _mm_andnot_pd(t, y)))

  // (this, ve) <- TwoProd(v1, v2) componentwise
  // We must have &v1 != &this (but we may have &v2 == &ve).
  inline void twosum(vec<double> &ve, vec<double> v1, vec<double> v2) {
    static const __m128d sign_mask  = _mm_set1_pd(-0.0f);
    __m128d tau, maxabs, minabs;
    
    r = _mm_add_pd(v1.r, v2.r);
    tau = _mm_cmpge_pd(_mm_andnot_pd(sign_mask, v1.r), _mm_andnot_pd(sign_mask,v2.r));
    maxabs = _mm_select_pd(tau, v1.r, v2.r);
    minabs = _mm_select_pd(tau, v2.r, v1.r);
    ve.r = _mm_add_pd(_mm_sub_pd(maxabs, r), minabs);
  }

  // (this, ve) <- TwoProd(v1, v2) componentwise
  // We must have &v1 != &this (but we may have &v2 == &ve).
  inline void twosum_in(vec<double> &ve, vec<double> v1) {
    static const __m128d sign_mask  = _mm_set1_pd(-0.0f);
    __m128d tau, maxabs, minabs;
    tau = _mm_cmpge_pd(_mm_andnot_pd(sign_mask, r), _mm_andnot_pd(sign_mask, v1.r));
    maxabs = _mm_select_pd(tau, r, v1.r);
    minabs = _mm_select_pd(tau, v1.r, r);
    r = _mm_add_pd(r, v1.r);
    ve.r = _mm_add_pd(_mm_sub_pd(maxabs, r), minabs);
  }
  
  // This code may be used as an alternative to Dekker's twoprod,
  // but it seems that it does not improves the running time...
  inline void twoprod(vec<double> &e, vec<double> a, vec<double> b) {
    static const __m128d mask = _mm_castsi128_pd(_mm_set1_epi64x(0xFFFFFFFFF8000000));
    r = _mm_mul_pd(a.r, b.r);
    __m128d ah = _mm_and_pd(a.r, mask);
    __m128d bh = _mm_and_pd(b.r, mask);
    __m128d al = _mm_sub_pd(a.r, ah);
    __m128d bl = _mm_sub_pd(b.r, bh);
    e.r = _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_sub_pd(_mm_mul_pd(ah, bh), r), _mm_mul_pd(ah, bl)), _mm_mul_pd(al, bh)), _mm_mul_pd(al, bl));
  }
  */
};

// Load 2 consecutive qd_real starting from address addr into 4 vec<double> v0, v1, v2, v3.
// Memory: x[0], ..., x[3], y[0], ..., y[3]
// v0 <- [x[0], y[0]]
// ...
// v3 <- [x[3], y[3]]
inline void load_qd_real(vec<double> &v0, vec<double> &v1, vec<double> &v2, vec<double> &v3, double const *addr) {
  __m128d x0 = _mm_load_pd(addr+0);
  __m128d x1 = _mm_load_pd(addr+2);
  __m128d x2 = _mm_load_pd(addr+4);
  __m128d x3 = _mm_load_pd(addr+6);
  v0.r = _mm_unpacklo_pd(x0, x2);
  v1.r = _mm_unpackhi_pd(x0, x2);
  v2.r = _mm_unpacklo_pd(x1, x3);
  v3.r = _mm_unpackhi_pd(x1, x3);
}

// Load the high order part of 2 consecutive qd_real starting from address addr into vec<double> v0.
// Memory: x[0], ..., x[3], y[0], ..., y[3]
// v0 <- [x[0], y[0]]
inline void loadhi_qd_real(vec<double> &v0, double const *addr) {
  __m128d x0 = _mm_load_pd(addr+0);
  __m128d x2 = _mm_load_pd(addr+4);
  v0.r = _mm_unpacklo_pd(x0, x2);
}

/////
#else
/////
#define VEC_ALIGN_CST 1

#error "All of this won't work with a non-simd isa..."

//////
#endif
//////

#define vec_align __attribute__((aligned(VEC_ALIGN_CST)))

// double-double = double-double * double-double
// (ph, pl) = (ah, al) * (bh, bl)
inline void ddmul(vec<double> &ph, vec<double> &pl, vec<double> ah, vec<double> al, vec<double> bh, vec<double> bl) {
  vec<double> p1, p2;

  ph.twoprod(pl, ah, bh);
  p1.mul(ah, bl);
  pl.add(pl, p1);
  p2.mul(al, bh);
  pl.add(pl, p2);
  ph.fasttwosum_in(pl, pl);
}

// double-double = double-double + double-double
// (ph, pl) = (ah, al) + (bh, bl)
// Note this is the "sloppy" version from the qd/dd library.
inline void ddadd(vec<double> &sh, vec<double> &sl, vec<double> ah, vec<double> al, vec<double> bh, vec<double> bl) {
  vec<double> t;
  
  sh.twosum(sl, ah, bh);
  t.add(al, bl);
  sl.add(sl, t);
  sh.fasttwosum_in(sl, sl);
}

// double-double += double-double
// (ph, pl) += (ah, al)
// Note this is the "sloppy" version from the qd/dd library. Also note that the difference
// with ddadd in really subtle...
inline void ddadd_in(vec<double> &sh, vec<double> &sl, vec<double> ah, vec<double> al) {
  vec<double> t;
  
  sh.twosum(t, sh, ah);
  t.add(t, sl);
  t.add(t, al);
  sh.fasttwosum_in(sl, t);
}

////////////////////////
#endif //#ifndef __VEC__
////////////////////////
