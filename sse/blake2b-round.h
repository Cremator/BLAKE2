/*
   BLAKE2 reference source code package - optimized C implementations

   Copyright 2012, Samuel Neves <sneves@dei.uc.pt>.  You may use this under the
   terms of the CC0, the OpenSSL Licence, or the Apache Public License 2.0, at
   your option.  The terms of these licenses can be found at:

   - CC0 1.0 Universal : http://creativecommons.org/publicdomain/zero/1.0
   - OpenSSL license   : https://www.openssl.org/source/license.html
   - Apache 2.0        : http://www.apache.org/licenses/LICENSE-2.0

   More information about the BLAKE2 hash function can be found at
   https://blake2.net.
*/
#ifndef BLAKE2B_ROUND_H
#define BLAKE2B_ROUND_H

#define LOADU(p)  vec_loadu1q( (const __m128i *)(p) )
#define STOREU(p,r) vec_storeu1q((__m128i *)(p), r)

#define TOF(reg) vec_cast4spto1q((reg))
#define TOI(reg) vec_cast4spto1q((reg))

#define LIKELY(x) __builtin_expect((x),1)


/* Microarchitecture-specific macros */
#ifndef HAVE_XOP
#ifdef HAVE_SSSE3
#define _mm_roti_epi64(x, c) \
    (-(c) == 32) ? _mm_shuffle_epi32((x), _MM_SHUFFLE(2,3,0,1))  \
    : (-(c) == 24) ? _mm_shuffle_epi8((x), r24) \
    : (-(c) == 16) ? _mm_shuffle_epi8((x), r16) \
    : (-(c) == 63) ? vec_bitxor1q(_mm_srli_epi64((x), -(c)), vec_add2sd((x), (x)))  \
    : vec_bitxor1q(_mm_srli_epi64((x), -(c)), _mm_slli_epi64((x), 64-(-(c))))
#else
#define _mm_roti_epi64(r, c) vec_bitxor1q(_mm_srli_epi64( (r), -(c) ),_mm_slli_epi64( (r), 64-(-(c)) ))
#endif
#else
/* ... */
#endif
#ifdef __VSX__
#define _mm_roti_epi64(r, c) vec_bitxor1q(vec_shiftrightimmediate2sd( (r), -(c) ),vec_shiftleftimmediate2sd( (r), 64-(-(c)) ))
#endif


#define G1(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1) \
  row1l = vec_add2sd(vec_add2sd(row1l, b0), row2l); \
  row1h = vec_add2sd(vec_add2sd(row1h, b1), row2h); \
  \
  row4l = vec_bitxor1q(row4l, row1l); \
  row4h = vec_bitxor1q(row4h, row1h); \
  \
  row4l = _mm_roti_epi64(row4l, -32); \
  row4h = _mm_roti_epi64(row4h, -32); \
  \
  row3l = vec_add2sd(row3l, row4l); \
  row3h = vec_add2sd(row3h, row4h); \
  \
  row2l = vec_bitxor1q(row2l, row3l); \
  row2h = vec_bitxor1q(row2h, row3h); \
  \
  row2l = _mm_roti_epi64(row2l, -24); \
  row2h = _mm_roti_epi64(row2h, -24); \

#define G2(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1) \
  row1l = vec_add2sd(vec_add2sd(row1l, b0), row2l); \
  row1h = vec_add2sd(vec_add2sd(row1h, b1), row2h); \
  \
  row4l = vec_bitxor1q(row4l, row1l); \
  row4h = vec_bitxor1q(row4h, row1h); \
  \
  row4l = _mm_roti_epi64(row4l, -16); \
  row4h = _mm_roti_epi64(row4h, -16); \
  \
  row3l = vec_add2sd(row3l, row4l); \
  row3h = vec_add2sd(row3h, row4h); \
  \
  row2l = vec_bitxor1q(row2l, row3l); \
  row2h = vec_bitxor1q(row2h, row3h); \
  \
  row2l = _mm_roti_epi64(row2l, -63); \
  row2h = _mm_roti_epi64(row2h, -63); \

#if defined(HAVE_SSSE3)
#define DIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h) \
  t0 = _mm_alignr_epi8(row2h, row2l, 8); \
  t1 = _mm_alignr_epi8(row2l, row2h, 8); \
  row2l = t0; \
  row2h = t1; \
  \
  t0 = row3l; \
  row3l = row3h; \
  row3h = t0;    \
  \
  t0 = _mm_alignr_epi8(row4h, row4l, 8); \
  t1 = _mm_alignr_epi8(row4l, row4h, 8); \
  row4l = t1; \
  row4h = t0;

#define UNDIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h) \
  t0 = _mm_alignr_epi8(row2l, row2h, 8); \
  t1 = _mm_alignr_epi8(row2h, row2l, 8); \
  row2l = t0; \
  row2h = t1; \
  \
  t0 = row3l; \
  row3l = row3h; \
  row3h = t0; \
  \
  t0 = _mm_alignr_epi8(row4l, row4h, 8); \
  t1 = _mm_alignr_epi8(row4h, row4l, 8); \
  row4l = t1; \
  row4h = t0;
#else

#define DIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h) \
  t0 = row4l;\
  t1 = row2l;\
  row4l = row3l;\
  row3l = row3h;\
  row3h = row4l;\
  row4l = vec_unpackhigh1sd(row4h, vec_unpacklow1sd(t0, t0)); \
  row4h = vec_unpackhigh1sd(t0, vec_unpacklow1sd(row4h, row4h)); \
  row2l = vec_unpackhigh1sd(row2l, vec_unpacklow1sd(row2h, row2h)); \
  row2h = vec_unpackhigh1sd(row2h, vec_unpacklow1sd(t1, t1))

#define UNDIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h) \
  t0 = row3l;\
  row3l = row3h;\
  row3h = t0;\
  t0 = row2l;\
  t1 = row4l;\
  row2l = vec_unpackhigh1sd(row2h, vec_unpacklow1sd(row2l, row2l)); \
  row2h = vec_unpackhigh1sd(t0, vec_unpacklow1sd(row2h, row2h)); \
  row4l = vec_unpackhigh1sd(row4l, vec_unpacklow1sd(row4h, row4h)); \
  row4h = vec_unpackhigh1sd(row4h, vec_unpacklow1sd(t1, t1))

#endif

#if defined(HAVE_SSE41)
#include "blake2b-load-sse41.h"
#elseif defined(__VSX__)
#include "blake2b-load-vsx.h"
#else
#include "blake2b-load-sse2.h"
#endif

#define ROUND(r) \
  LOAD_MSG_ ##r ##_1(b0, b1); \
  G1(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1); \
  LOAD_MSG_ ##r ##_2(b0, b1); \
  G2(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1); \
  DIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h); \
  LOAD_MSG_ ##r ##_3(b0, b1); \
  G1(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1); \
  LOAD_MSG_ ##r ##_4(b0, b1); \
  G2(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h,b0,b1); \
  UNDIAGONALIZE(row1l,row2l,row3l,row4l,row1h,row2h,row3h,row4h);

#endif
