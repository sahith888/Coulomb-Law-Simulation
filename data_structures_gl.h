#pragma once
#include<immintrin.h>
#include<stdio.h>

typedef struct mat2 {
	float m[4];

	
};

struct mat2 e = { 1.0, 1.0, 1.0 ,1.0 };

typedef struct vec3 {
	float *x;
	float *y;
	float *z;
	void initvec3(int width, int height) {
		this->x = (float*)malloc(sizeof(float) * width * height + 32);
		this->y = (float*)malloc(sizeof(float) * width * height + 32);
		this->z = (float*)malloc(sizeof(float) * width * height + 32);
		for (int i = 0; i < width * height; i++) { // 0 -> 1
			this->x[i] = 0.0f;
			this->y[i] = 0.0f;
			this->z[i] = 0.0f;
		}
	}
	void initfvec3() {
		this->x = (float*)malloc(sizeof(float));
		this->y = (float*)malloc(sizeof(float));
		this->z = (float*)malloc(sizeof(float));
	}
	void set(float r, float g, float b) {
		this->x[0] = r;
		this->y[0] = g;
		this->z[0] = b;
	}
	/*
	void ADD(float addend) {
		for (int i = 0; i < w * h; i++) { // 0 -> 1
			this->x[i] += addend;
			this->y[i] += addend;
		}
	}
	void MUL(float mul) {
		for (int i = 0; i < w * h; i++) { // 0 -> 1
			this->x[i] *= mul;
			this->y[i] *= mul;
		}
	}
	*/

} vec3;

#define float32 __m256

typedef struct _declspec(align(32)) float2_vector_SIMD {
	__m256 x;
	__m256 y;

} fvec2;

typedef struct _declspec(align(32)) float3_vector_SIMD {
	__m256 r;
	__m256 g;
	__m256 b;
} fvec3; // temp name

fvec2 add(fvec2 a, fvec2 b) {
	a.x = _mm256_add_ps(a.x, b.x);
	a.y = _mm256_add_ps(a.y, b.y);
	return a;
}
fvec3 add(fvec3 a, fvec3 b) {
	a.r = _mm256_add_ps(a.r, b.r);
	a.g = _mm256_add_ps(a.g, b.g);
	a.b = _mm256_add_ps(a.b, b.b);
	return a;
}

fvec2 mul(fvec2 a, fvec2 b) {
	a.x = _mm256_mul_ps(a.x, b.x);
	a.y = _mm256_mul_ps(a.y, b.y);
	return a;
}
fvec3 mul(fvec3 a, fvec3 b) {
	a.r = _mm256_mul_ps(a.r, b.r);
	a.g = _mm256_mul_ps(a.g, b.g);
	a.b = _mm256_mul_ps(a.b, b.b);
	return a;
}

fvec2 set(__m256 s, __m256 t) {
	fvec2 st = { s, t };
	return st;
}
fvec3 set(__m256 s, __m256 t, __m256 u) {
	fvec3 stu = { s, t, u };
	return stu;
}

fvec2 clamp(fvec2 st, __m256 min, __m256 max) {

	__m256 s = _mm256_min_ps(max, st.x);
	s = _mm256_max_ps(min, s);
	__m256 t = _mm256_min_ps(max, st.y);
	t = _mm256_max_ps(min, t);
	fvec2 xy = { s, t };
	return xy;
}
fvec3 clamp(fvec3 stu, __m256 min, __m256 max) {

	stu.r = _mm256_min_ps(max, stu.r);
	stu.r = _mm256_max_ps(min, stu.r);
	stu.g = _mm256_min_ps(max, stu.g);
	stu.g = _mm256_max_ps(min, stu.g);
	stu.b = _mm256_min_ps(max, stu.b);
	stu.b = _mm256_max_ps(min, stu.b);
	return stu;
}

fvec2 sub(fvec2 a, fvec2 b) {
	a.x = _mm256_sub_ps(a.x, b.x);
	a.y = _mm256_sub_ps(a.y, b.y);
	return a;
}
fvec3 sub(fvec3 a, fvec3 b) {
	a.r = _mm256_sub_ps(a.r, b.r);
	a.g = _mm256_sub_ps(a.g, b.g);
	a.b = _mm256_sub_ps(a.b, b.b);
	return a;
}

float32 length(fvec2 st) {
	st.x = _mm256_fmadd_ps(st.x, st.x, _mm256_mul_ps(st.y, st.y));
	return _mm256_sqrt_ps(st.x);
}
float32 length(fvec3 rgb) {
	rgb.r = _mm256_fmadd_ps(rgb.r, rgb.r, _mm256_fmadd_ps(rgb.g, rgb.g, _mm256_mul_ps(rgb.b, rgb.b))); //x^2+y^2+z^2
	return _mm256_sqrt_ps(rgb.r);
}

float32 distance(fvec2 p0, fvec2 p1) {
	return length(sub(p0, p1));
}
float32 distance(fvec3 p0, fvec3 p1) {
	return length(sub(p0, p1));
}

fvec2 normalize(fvec2 st) {
	float32 magnitude = _mm256_fmadd_ps(st.x, st.x, _mm256_mul_ps(st.y, st.y));
	magnitude = _mm256_rsqrt_ps(magnitude);
	st.x = _mm256_mul_ps(magnitude, st.x);
	st.y = _mm256_mul_ps(magnitude, st.y);
	return st;
}
fvec3 normalize(fvec3 rgb) {
	float32 magnitude = _mm256_fmadd_ps(rgb.r, rgb.r, _mm256_fmadd_ps(rgb.g, rgb.g, _mm256_mul_ps(rgb.b, rgb.b)));
	magnitude = _mm256_rsqrt_ps(magnitude);
	rgb.r = _mm256_mul_ps(magnitude, rgb.r);
	rgb.g = _mm256_mul_ps(magnitude, rgb.g);
	rgb.b = _mm256_mul_ps(magnitude, rgb.b);
	return rgb;
}

float32 dot(fvec2 a, fvec2 b) {
	return _mm256_fmadd_ps(a.x, b.x, _mm256_mul_ps(a.y, b.y));
}
float32 dot(fvec3 a, fvec3 b) {
	return _mm256_fmadd_ps(a.r, b.r, _mm256_fmadd_ps(a.g, b.g, _mm256_mul_ps(a.b, b.b)));
}

fvec3 cross(fvec3 x, fvec3 y) {
	fvec3 product;
	product.r = _mm256_fmsub_ps(x.g, y.b, _mm256_mul_ps(y.g, x.b));
	product.g = _mm256_fmsub_ps(x.b, y.r, _mm256_mul_ps(y.b, x.r));
	product.b = _mm256_fmsub_ps(x.r, y.g, _mm256_mul_ps(y.r, x.g));
	return product;
}