/************************************************/
#pragma once
/************************************************/
#include <math.h>
#include <stddef.h> // NULL
/************************************************/
#ifdef __SSE__
# include <xmmintrin.h>
#endif
/************************************************/

typedef struct
#ifdef __SSE__
__attribute__((aligned(16)))
#endif
{
	float f32[4];
} Vec4f_t;

/************************************************/

#define VEC4F_EMPTY (Vec4f_t){{0,0,0,0}}

/************************************************/

static inline Vec4f_t Vec4f_Add(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_add_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] + b->f32[0];
	y.f32[1] = a->f32[1] + b->f32[1];
	y.f32[2] = a->f32[2] + b->f32[2];
	y.f32[3] = a->f32[3] + b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_Addi(const Vec4f_t *a, float b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_add_ps(_mm_load_ps(a->f32), _mm_set1_ps(b));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] + b;
	y.f32[1] = a->f32[1] + b;
	y.f32[2] = a->f32[2] + b;
	y.f32[3] = a->f32[3] + b;
#endif
	return y;
}

/************************************************/

static inline Vec4f_t Vec4f_Sub(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_sub_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] - b->f32[0];
	y.f32[1] = a->f32[1] - b->f32[1];
	y.f32[2] = a->f32[2] - b->f32[2];
	y.f32[3] = a->f32[3] - b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_Subi(const Vec4f_t *a, float b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_sub_ps(_mm_load_ps(a->f32), _mm_set1_ps(b));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] - b;
	y.f32[1] = a->f32[1] - b;
	y.f32[2] = a->f32[2] - b;
	y.f32[3] = a->f32[3] - b;
#endif
	return y;
}

/************************************************/

static inline Vec4f_t Vec4f_Mul(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_mul_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] * b->f32[0];
	y.f32[1] = a->f32[1] * b->f32[1];
	y.f32[2] = a->f32[2] * b->f32[2];
	y.f32[3] = a->f32[3] * b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_Muli(const Vec4f_t *a, float b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_mul_ps(_mm_load_ps(a->f32), _mm_set1_ps(b));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] * b;
	y.f32[1] = a->f32[1] * b;
	y.f32[2] = a->f32[2] * b;
	y.f32[3] = a->f32[3] * b;
#endif
	return y;
}

/************************************************/

static inline Vec4f_t Vec4f_Div(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_div_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] / b->f32[0];
	y.f32[1] = a->f32[1] / b->f32[1];
	y.f32[2] = a->f32[2] / b->f32[2];
	y.f32[3] = a->f32[3] / b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_DivSafe(const Vec4f_t *a, const Vec4f_t *b, const Vec4f_t *DivByZeroValue) {
	static const Vec4f_t Zero = VEC4F_EMPTY;
	if(!DivByZeroValue) DivByZeroValue = &Zero;

	Vec4f_t y;
#ifdef __SSE__
	__m128 z  = _mm_load_ps(DivByZeroValue->f32);
	__m128 c  = _mm_cmpeq_ps(_mm_load_ps(b->f32), _mm_setzero_ps());
	__m128 yv = _mm_div_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	       z  = _mm_and_ps(c, z);
	       yv = _mm_andnot_ps(c, yv);
	       yv = _mm_or_ps(yv, z);
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = (b->f32[0] == 0.0f) ? DivByZeroValue->f32[0] : (a->f32[0] / b->f32[0]);
	y.f32[1] = (b->f32[1] == 0.0f) ? DivByZeroValue->f32[1] : (a->f32[1] / b->f32[1]);
	y.f32[2] = (b->f32[2] == 0.0f) ? DivByZeroValue->f32[2] : (a->f32[2] / b->f32[2]);
	y.f32[3] = (b->f32[3] == 0.0f) ? DivByZeroValue->f32[3] : (a->f32[3] / b->f32[3]);
#endif
	return y;
}

static inline Vec4f_t Vec4f_Divi(const Vec4f_t *a, float b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_div_ps(_mm_load_ps(a->f32), _mm_set1_ps(b));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = a->f32[0] / b;
	y.f32[1] = a->f32[1] / b;
	y.f32[2] = a->f32[2] / b;
	y.f32[3] = a->f32[3] / b;
#endif
	return y;
}

static inline Vec4f_t Vec4f_InverseDivi(const Vec4f_t *a, float b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_div_ps(_mm_set1_ps(b), _mm_load_ps(a->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = b / a->f32[0];
	y.f32[1] = b / a->f32[1];
	y.f32[2] = b / a->f32[2];
	y.f32[3] = b / a->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_InverseDiviSafe(const Vec4f_t *a, float b, const Vec4f_t *DivByZeroValue) {
	static const Vec4f_t Zero = VEC4F_EMPTY;
	if(!DivByZeroValue) DivByZeroValue = &Zero;

	Vec4f_t y;
#ifdef __SSE__
	__m128 z  = _mm_load_ps(DivByZeroValue->f32);
	__m128 c  = _mm_cmpeq_ps(_mm_load_ps(a->f32), _mm_setzero_ps());
	__m128 yv = _mm_div_ps(_mm_set1_ps(b), _mm_load_ps(a->f32));
	       z  = _mm_and_ps(c, z);
	       yv = _mm_andnot_ps(c, yv);
	       yv = _mm_or_ps(yv, z);
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = (a->f32[0] == 0.0f) ? DivByZeroValue->f32[0] : (b / a->f32[0]);
	y.f32[1] = (a->f32[1] == 0.0f) ? DivByZeroValue->f32[1] : (b / a->f32[1]);
	y.f32[2] = (a->f32[2] == 0.0f) ? DivByZeroValue->f32[2] : (b / a->f32[2]);
	y.f32[3] = (a->f32[3] == 0.0f) ? DivByZeroValue->f32[3] : (b / a->f32[3]);
#endif
	return y;
}

/************************************************/

static inline Vec4f_t Vec4f_Abs(const Vec4f_t *x) {
	Vec4f_t y;
#ifdef __SSE__
	//! Assumes no infinities/NaN
	__m128 yv = _mm_andnot_ps(_mm_set1_ps(-0.0f), _mm_load_ps(x->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = fabsf(x->f32[0]);
	y.f32[1] = fabsf(x->f32[1]);
	y.f32[2] = fabsf(x->f32[2]);
	y.f32[3] = fabsf(x->f32[3]);
#endif
	return y;
}

static inline Vec4f_t Vec4f_Sqrt(const Vec4f_t *x) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_sqrt_ps(_mm_load_ps(x->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = sqrtf(x->f32[0]);
	y.f32[1] = sqrtf(x->f32[1]);
	y.f32[2] = sqrtf(x->f32[2]);
	y.f32[3] = sqrtf(x->f32[3]);
#endif
	return y;
}

static inline float Vec4f_SumOf(const Vec4f_t *x) {
	return x->f32[0] +
	       x->f32[1] +
	       x->f32[2] +
	       x->f32[3] ;
}

static inline float Vec4f_Dot(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y = Vec4f_Mul(a, b);
	return Vec4f_SumOf(&y);
}

static inline float Vec4f_Length2(const Vec4f_t *x) {
	return Vec4f_Dot(x, x);
}

static inline float Vec4f_Length(const Vec4f_t *x) {
	return sqrtf(Vec4f_Length2(x));
}

static inline float Vec4f_Dist2(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t x = Vec4f_Sub(a, b);
	return Vec4f_Length2(&x);
}

static inline float Vec4f_Dist(const Vec4f_t *a, const Vec4f_t *b) {
	return sqrtf(Vec4f_Dist2(a, b));
}

static inline float Vec4f_DistL1(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t x = Vec4f_Sub(a, b);
	        x = Vec4f_Abs(&x);
	return Vec4f_SumOf(&x);
}

/************************************************/

static inline Vec4f_t Vec4f_Broadcast(float x) {
	Vec4f_t y;
	y.f32[0] = y.f32[1] = y.f32[2] = y.f32[3] = x;
	return y;
}

static inline float Vec4f_MinOf(const Vec4f_t *x) {
	float y = x->f32[0];
	if(x->f32[1] < y) y = x->f32[1];
	if(x->f32[2] < y) y = x->f32[2];
	if(x->f32[3] < y) y = x->f32[3];
	return y;
}

static inline float Vec4f_MaxOf(const Vec4f_t *x) {
	float y = x->f32[0];
	if(x->f32[1] > y) y = x->f32[1];
	if(x->f32[2] > y) y = x->f32[2];
	if(x->f32[3] > y) y = x->f32[3];
	return y;
}

static inline Vec4f_t Vec4f_Min(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_min_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = (a->f32[0] < b->f32[0]) ? a->f32[0] : b->f32[0];
	y.f32[1] = (a->f32[1] < b->f32[1]) ? a->f32[1] : b->f32[1];
	y.f32[2] = (a->f32[2] < b->f32[2]) ? a->f32[2] : b->f32[2];
	y.f32[3] = (a->f32[3] < b->f32[3]) ? a->f32[3] : b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_Max(const Vec4f_t *a, const Vec4f_t *b) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_max_ps(_mm_load_ps(a->f32), _mm_load_ps(b->f32));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = (a->f32[0] > b->f32[0]) ? a->f32[0] : b->f32[0];
	y.f32[1] = (a->f32[1] > b->f32[1]) ? a->f32[1] : b->f32[1];
	y.f32[2] = (a->f32[2] > b->f32[2]) ? a->f32[2] : b->f32[2];
	y.f32[3] = (a->f32[3] > b->f32[3]) ? a->f32[3] : b->f32[3];
#endif
	return y;
}

static inline Vec4f_t Vec4f_Round(const Vec4f_t *x) {
	Vec4f_t y;
	y.f32[0] = roundf(x->f32[0]);
	y.f32[1] = roundf(x->f32[1]);
	y.f32[2] = roundf(x->f32[2]);
	y.f32[3] = roundf(x->f32[3]);
	return y;
}

static inline Vec4f_t Vec4f_Quantize(const Vec4f_t *x, const Vec4f_t *Depth) {
	Vec4f_t y = Vec4f_Mul(x, Depth);
	        y = Vec4f_Round(&y);
	return Vec4f_DivSafe(&y, Depth, NULL);
}

static inline Vec4f_t Vec4f_Clamp(const Vec4f_t *x, float Min, float Max) {
	Vec4f_t y;
#ifdef __SSE__
	__m128 yv = _mm_load_ps(x->f32);
	yv = _mm_max_ps(yv, _mm_set1_ps(Min));
	yv = _mm_min_ps(yv, _mm_set1_ps(Max));
	_mm_store_ps(y.f32, yv);
#else
	y.f32[0] = (x->f32[0] < Min) ? Min : (x->f32[0] > Max) ? Max : x->f32[0];
	y.f32[1] = (x->f32[1] < Min) ? Min : (x->f32[1] > Max) ? Max : x->f32[1];
	y.f32[2] = (x->f32[2] < Min) ? Min : (x->f32[2] > Max) ? Max : x->f32[2];
	y.f32[3] = (x->f32[3] < Min) ? Min : (x->f32[3] > Max) ? Max : x->f32[3];
#endif
	return y;
}

/************************************************/
//! EOF
/************************************************/
