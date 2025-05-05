#pragma once
#include<stdlib.h>
#include<math.h>
#include"data_structures_gl.h"
#define f(x) _mm256_set1_ps(x) // avx2 quick const definition
#define i(x) _mm256_set1_epi32(x) // avx2 quick const definition


struct texture {
	float *u;
	float *v;
	float *c;
	int w;
	int h;

	inline void TEXTURE_INITIALIZATION(int w, int h, float* color) {
		this->w = w;
		this->h = h;
		this->u = (float*)malloc(sizeof(float) * w * h + 32);
		this->v = (float*)malloc(sizeof(float) * w * h + 32);

		this->c = color;

		for (int i = 0; i < w * h; i++) {
			this->u[i] = (float)(i % w) / w;
			this->v[i] = (float)(i / w) / h;
		}
	}

	void FLIP_TEXTURE_VERT() {
		int width = this->w;
		int height = this->h;
		float* copy = this->c;
		float* newtex = (float*)malloc(sizeof(float) * w * h + 32);

		for (int i = 0; i < w * h; i++) {
			int v = i / width;
			int dv = height - v;
			int u = i % width;
			int du = width - u;
			int vert_i = u + (v * width); // vertical index ( pixel column index)
			int fvert_i = (dv * width) - du; // opposite vertical index
			newtex[vert_i] = copy[fvert_i];

		}
		this->c = newtex;
		//free(copy);
		copy = NULL;
	}
};

float  *TEXTURE1_DATA = new float[10 * 10] {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, // 10x10
						 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.2,
						 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.3,
						 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4,
						 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5,
						 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6,
						 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7,
						 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8,
						 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.9,
						 1.0, 0.2, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
};
/*
void MAP_TEXTURE(texture* tex, vec3 col, int w, int h, float * u, float * v) {
	float* s = tex->u; // texCoord S x-axis
	float* t = tex->v; // texCoord T y-axis
	float* texCol = tex->c;
	int texWidth = tex->w;
	int texHeight = tex->h;
	int j = 0;
	for (int i = 0; i < w * h * 3; i += 3) {
		//int iX = (int)trunc(u[i] * texWidth); // get the indexes
		//int iY = (int)trunc(v[i] * texHeight); // additionally it must be truncated rather than rounded as it wont map correctly
		int iX = (int)(u[j] * texWidth); // get the indexes
		int iY = (int)(v[j] * texHeight); // additionally it must be truncated rather than rounded as it wont map correctly

		float fragColor = texCol[(iX + (texWidth * iY)) * 3 + 2];

		col.x[j] += fragColor; // index into texture color vector and put that into the frame buffer
		col.y[j] += texCol[(iX + (texWidth * iY)) * 3 +1];
		col.z[j] += texCol[(iX + (texWidth * iY)) * 3];

		//frame_buffer[i] = texCol[iX + (texWidth * iY)]; // index into texture color vector and put that into the frame buffer
		j++;
	}
}
*/
struct sampler2D {
	int width;
	int height;
	float w;
	float h;
	int subpixelcount;
	float* data;

	void init(int wid, int hei, float* col) {
		this->width = wid;
		this->height = hei;
		this->w = (float)wid;
		this->h = (float)hei;

		this->data = col;
		this->subpixelcount = wid * hei * 3;
	}
};

struct vec3 texture2d(sampler2D tex, float u, float v) {
	struct vec3 col;
	float rgb[3];
	// Point to the same pre-allocated array
	col.x = &rgb[0];
	col.y = &rgb[1];
	col.z = &rgb[2];

	int x = (int)(tex.w * min(max(u, 0.0), 0.99)); // temp clamp of coords, it needs to be changed for SIMD
	int y = (int)(tex.h * min(max(v, 0.0), 0.99));
	int index = 3 * (x + tex.w * y);

	col.x[0] = tex.data[index];      // red channel
	col.y[0] = tex.data[index + 1];  // green channel
	col.z[0] = tex.data[index + 2];  // blue channel

	return col;
}

fvec3 texture(sampler2D tex, fvec2 uv) {
	float w = tex.w;
	float h = tex.h;
	uv = clamp(uv, f(0.0), f(0.99));
	__m256 u = _mm256_mul_ps(uv.x, f(w));
	__m256 v = _mm256_mul_ps(uv.y, f(h));

	__m256i x = _mm256_cvtps_epi32(u);
	__m256i y = _mm256_cvtps_epi32(v);

	y = _mm256_mullo_epi32(i(tex.width), y);
	y = _mm256_add_epi32(x, y);
	y = _mm256_mullo_epi32(i(3), y);
	y = _mm256_min_epi32(i(tex.subpixelcount), y);

	__m256 col_r = _mm256_i32gather_ps(tex.data, y, 4);
	__m256 col_g = _mm256_i32gather_ps(tex.data, _mm256_add_epi32(i(1), y), 4);
	__m256 col_b = _mm256_i32gather_ps(tex.data, _mm256_add_epi32(i(2), y), 4);
	return fvec3{ col_r, col_g, col_b };

}

fvec3 gather_texel(const sampler2D* tex, __m256i idx) {
	__m256 r = _mm256_i32gather_ps(tex->data, idx, 4);
	__m256 g = _mm256_i32gather_ps(tex->data, _mm256_add_epi32(idx, _mm256_set1_epi32(1)), 4);
	__m256 b = _mm256_i32gather_ps(tex->data, _mm256_add_epi32(idx, _mm256_set1_epi32(2)), 4);
	return fvec3{ r, g, b };
}

fvec3 bilinear(const sampler2D tex, fvec2 uv) {
	float w = (float)tex.width;
	float h = (float)tex.height;
	uv = clamp(uv, f(0.0f), f(0.99f));

	// Compute floating point pixel positions
	__m256 u = _mm256_mul_ps(uv.x, _mm256_set1_ps(w - 1));
	__m256 v = _mm256_mul_ps(uv.y, _mm256_set1_ps(h - 1));

	// Calculate integer and fractional parts
	__m256 x = _mm256_floor_ps(u);
	__m256 y = _mm256_floor_ps(v);
	__m256i x0 = _mm256_cvtps_epi32(u);
	__m256i y0 = _mm256_cvtps_epi32(v);

	__m256 u_ratio = _mm256_sub_ps(u, _mm256_cvtepi32_ps(x0));
	__m256 v_ratio = _mm256_sub_ps(v, _mm256_cvtepi32_ps(y0));

	__m256 u_opposite = _mm256_sub_ps(_mm256_set1_ps(1.0f), u_ratio);
	__m256 v_opposite = _mm256_sub_ps(_mm256_set1_ps(1.0f), v_ratio);

	// Calculate indices of the four surrounding texels
	__m256i texWidth = _mm256_set1_epi32(tex.width);
	__m256i texHeight = _mm256_set1_epi32(tex.height);

	__m256i x1 = _mm256_min_epi32(texWidth, _mm256_add_epi32(x0, _mm256_set1_epi32(1)));
	__m256i y1 = _mm256_min_epi32(texHeight, _mm256_add_epi32(y0, _mm256_set1_epi32(1)));

	x0 = _mm256_max_epi32(_mm256_set1_epi32(0), x0);
	y0 = _mm256_max_epi32(_mm256_set1_epi32(0), y0);
	x1 = _mm256_max_epi32(_mm256_set1_epi32(0), x1);
	y1 = _mm256_max_epi32(_mm256_set1_epi32(0), y1);

	// Compute linear indices of the four surrounding texels
	__m256i idx00 = _mm256_mullo_epi32(y0, texWidth);
	idx00 = _mm256_add_epi32(idx00, x0);
	idx00 = _mm256_mullo_epi32(idx00, _mm256_set1_epi32(3));

	__m256i idx10 = _mm256_mullo_epi32(y0, texWidth);
	idx10 = _mm256_add_epi32(idx10, x1);
	idx10 = _mm256_mullo_epi32(idx10, _mm256_set1_epi32(3));

	__m256i idx01 = _mm256_mullo_epi32(y1, texWidth);
	idx01 = _mm256_add_epi32(idx01, x0);
	idx01 = _mm256_mullo_epi32(idx01, _mm256_set1_epi32(3));

	__m256i idx11 = _mm256_mullo_epi32(y1, texWidth);
	idx11 = _mm256_add_epi32(idx11, x1);
	idx11 = _mm256_mullo_epi32(idx11, _mm256_set1_epi32(3));

	// Gather texels
	fvec3 texel00 = gather_texel(&tex, idx00);
	fvec3 texel10 = gather_texel(&tex, idx10);
	fvec3 texel01 = gather_texel(&tex, idx01);
	fvec3 texel11 = gather_texel(&tex, idx11);

	// Perform bilinear interpolation
	__m256 r0 = _mm256_add_ps(_mm256_mul_ps(texel00.r, u_opposite), _mm256_mul_ps(texel10.r, u_ratio));
	__m256 r1 = _mm256_add_ps(_mm256_mul_ps(texel01.r, u_opposite), _mm256_mul_ps(texel11.r, u_ratio));
	__m256 r = _mm256_add_ps(_mm256_mul_ps(r0, v_opposite), _mm256_mul_ps(r1, v_ratio));

	__m256 g0 = _mm256_add_ps(_mm256_mul_ps(texel00.g, u_opposite), _mm256_mul_ps(texel10.g, u_ratio));
	__m256 g1 = _mm256_add_ps(_mm256_mul_ps(texel01.g, u_opposite), _mm256_mul_ps(texel11.g, u_ratio));
	__m256 g = _mm256_add_ps(_mm256_mul_ps(g0, v_opposite), _mm256_mul_ps(g1, v_ratio));

	__m256 b0 = _mm256_add_ps(_mm256_mul_ps(texel00.b, u_opposite), _mm256_mul_ps(texel10.b, u_ratio));
	__m256 b1 = _mm256_add_ps(_mm256_mul_ps(texel01.b, u_opposite), _mm256_mul_ps(texel11.b, u_ratio));
	__m256 b = _mm256_add_ps(_mm256_mul_ps(b0, v_opposite), _mm256_mul_ps(b1, v_ratio));

	return fvec3{ r, g, b };
}

vec3 bilinear2d(sampler2D tex, float u, float v) {
	vec3 p1;
	vec3 p2;
	vec3 p3;
	vec3 p4;
	vec3 col;
	float rgb[15];
	p1.x = &rgb[0];
	p1.y = &rgb[1];
	p1.z = &rgb[2];
	p2.x = &rgb[3];
	p2.y = &rgb[4];
	p2.z = &rgb[5];
	p3.x = &rgb[6];
	p3.y = &rgb[7];
	p3.z = &rgb[8];
	p4.x = &rgb[9];
	p4.y = &rgb[10];
	p4.z = &rgb[11];
	col.x = &rgb[12];
	col.y = &rgb[13];
	col.z = &rgb[14];
	// Clamp coordinates to the [0, 1) range
	u = min(max(u, 0.0), 0.99);
	v = min(max(v, 0.0), 0.99);
	// Calculate the exact position within the texture
	float x = u * (tex.w - 1);
	float y = v * (tex.h - 1);

	// Calculate the integer part of the position
	int x0 = (int)x;
	int y0 = (int)y;

	// Calculate the fractional part of the position
	float u_ratio = x - x0;
	float v_ratio = y - y0;

	// Calculate the coordinates of the 4 surrounding texels
	int x1 = x0 + 1;
	int y1 = y0 + 1;

	// Ensure the coordinates don't exceed the texture dimensions
	if (x1 >= tex.w) x1 = tex.w - 1;
	if (y1 >= tex.h) y1 = tex.h - 1;

	// Fetch the colors of the 4 surrounding texels
	p1.set(tex.data[3 * (int)(x0 + tex.w * y0)], tex.data[3 * (int)(x0 + tex.w * y0) + 1], tex.data[3 * (int)(x0 + tex.w * y0) + 2]);
	p2.set(tex.data[3 * (int)(x1 + tex.w * y0)], tex.data[3 * (int)(x1 + tex.w * y0) + 1], tex.data[3 * (int)(x1 + tex.w * y0) + 2]);
	p3.set(tex.data[3 * (int)(x0 + tex.w * y1)], tex.data[3 * (int)(x0 + tex.w * y1) + 1], tex.data[3 * (int)(x0 + tex.w * y1) + 2]);
	p4.set(tex.data[3 * (int)(x1 + tex.w * y1)], tex.data[3 * (int)(x1 + tex.w * y1) + 1], tex.data[3 * (int)(x1 + tex.w * y1) + 2]);

	// Interpolate between the texels
	col.x[0] = (1 - u_ratio) * (1 - v_ratio) * p1.x[0] +
		u_ratio * (1 - v_ratio) * p2.x[0] +
		(1 - u_ratio) * v_ratio * p3.x[0] +
		u_ratio * v_ratio * p4.x[0];
	col.y[0] = (1 - u_ratio) * (1 - v_ratio) * p1.y[0] +
		u_ratio * (1 - v_ratio) * p2.y[0] +
		(1 - u_ratio) * v_ratio * p3.y[0] +
		u_ratio * v_ratio * p4.y[0];
	col.z[0] = (1 - u_ratio) * (1 - v_ratio) * p1.z[0] +
		u_ratio * (1 - v_ratio) * p2.z[0] +
		(1 - u_ratio) * v_ratio * p3.z[0] +
		u_ratio * v_ratio * p4.z[0];

	return col;
}