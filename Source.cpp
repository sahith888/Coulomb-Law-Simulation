#include<iostream>
#include<utility>

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<Windows.h>
#include<tchar.h>
#include<wchar.h>
#include"intcast.h"
#include"data_structures_gl.h"

#include<immintrin.h>
#include"texture.h"

#include"file_loader.h"

#define w 300
#define h 300
constexpr auto PI = 3.14159265358979323846  /* pi */;
constexpr auto TAU = 2.0 * 3.14159265358979323846  /* pi */;
#define WHITE 16777215

#define UNICODE
#define _UNICODE
#include <windows.h>
#include <stdbool.h>
#include <stdint.h>

#define UWIDTH 1280
#define UHEIGHT 720
float zoom = 1.0;
static bool quit = false;
struct {
	int width;
	int height;
	uint32_t* pixels;
} frame = { 0 };

LRESULT CALLBACK WindowProcessMessage(HWND, UINT, WPARAM, LPARAM);
#if RAND_MAX == 32767
#define Rand32() ((rand() << 16) + (rand() << 1) + (rand() & 1))
#else
#define Rand32() rand()
#endif

static BITMAPINFO frame_bitmap_info;
static HBITMAP frame_bitmap = 0;
static HDC frame_device_context = 0;

#define f(x) _mm256_set1_ps(x) // avx2 quick const definition
//deprecated shader functions
__m256 abs2d(__m256 x) {
	return _mm256_andnot_ps(f(-0.0f), x);
}
__m256 clamp2d(__m256 x, __m256 min, __m256 max) {
	__m256 x1 = _mm256_min_ps(max, x);
	return _mm256_max_ps(min, x1);
}

__m256 smoothstep2d(__m256 edge0, __m256 edge1, __m256 x) {
	__m256 x_e0 = _mm256_sub_ps(x, edge0);
	__m256 e1_e0 = _mm256_sub_ps(edge1, edge0);
	__m256 div = _mm256_div_ps(x_e0, e1_e0);
	__m256 t = clamp2d(div, f(0.0), f(1.0));

	__m256 mintwo_t = _mm256_mul_ps(f(-2.0), t);
	mintwo_t = _mm256_add_ps(f(3.0), mintwo_t);
	t = _mm256_mul_ps(t, t);
	return _mm256_mul_ps(t, mintwo_t);
}

__m256 step2d(__m256 edge, __m256 x) {
	__m256 comparison_mask = _mm256_cmp_ps(x, edge, _CMP_GE_OQ);
	return _mm256_and_ps(comparison_mask, f(1.0));
}

float sigmoid_bonding(float x) {
	return 2.0 / (1.0 + exp(-1.0 * (x - 0.05))) - 1; // the "Locking Function" so it is at a stable bond energy at 0.05 distance where there will be no forces if bond length is shorter the particles get repelled
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PSTR pCmdLine, int nCmdShow) {
	const wchar_t window_class_name[] = L"My Window Class";
	static WNDCLASS window_class = { 0 };
	window_class.lpfnWndProc = WindowProcessMessage;
	window_class.hInstance = hInstance;
	window_class.lpszClassName = window_class_name;
	RegisterClass(&window_class);

	frame_bitmap_info.bmiHeader.biSize = sizeof(frame_bitmap_info.bmiHeader);
	frame_bitmap_info.bmiHeader.biPlanes = 1;
	frame_bitmap_info.bmiHeader.biBitCount = 32;
	frame_bitmap_info.bmiHeader.biCompression = BI_RGB;
	frame_device_context = CreateCompatibleDC(0);

	static HWND window_handle;
	window_handle = CreateWindow(window_class_name, L"Drawing Pixels", WS_OVERLAPPEDWINDOW | WS_VISIBLE,
		640, 480, UWIDTH, UHEIGHT, NULL, NULL, hInstance, NULL);
	if (window_handle == NULL) { return -1; }

	/*for (int i = 0; i < frame.width * frame.height; i++) {   // debug shader
		if (frame.pixels[i] > 254) {
			frame.pixels[i] = WHITE;
		}
	}*/
	float* u = (float*)malloc(sizeof(float) * frame.width * frame.height + 32);
	float* v = (float*)malloc(sizeof(float) * frame.height * frame.width + 32);
	for (int i = 0; i < frame.width * frame.height; i++) { u[i] = (float)(i % frame.width) / frame.width; } // generation of normalized
	for (int j = 0; j < frame.height * frame.width; j++) { v[j] = (float)(j / frame.width) / frame.height; } // generation fo normalized 
	//for (int j = 0; j < frame.width * frame.height; j++) { u[j] = u[j] * (float)frame.width / frame.height; }
	float time = 0.0f;
	__m256 mask = _mm256_set1_ps(-0.0f); // mask used for absolute value
	__m256 const_vec = _mm256_set1_ps(0.5f);

	struct rotation_matrix { // 2x2
		__m256 r;
		__m256 g;
		__m256 b;
		__m256 a;

		void update(float a) {
			float cosine = cos(a);
			float sine = sin(a);
			this->r = _mm256_set1_ps(cosine);
			this->g = _mm256_set1_ps(-sine);
			this->b = _mm256_set1_ps(sine);
			this->a = _mm256_set1_ps(cosine);
		}

	};
	struct rotation_matrix romat = { 0.0 };
	romat.update(1.0);

	struct texture tex1 = { 0 };
	tex1.TEXTURE_INITIALIZATION(10, 10, TEXTURE1_DATA);

	static unsigned char* opened = load_raw_filedata("cube.bmp"); // first BGR color index starts at 18
	unsigned int* ivec = (unsigned int*)opened;
	

	struct texture tex2 = {0};
	float* positionliteral = new float[1920 * 1080 * 3] { 0.0 };
	for (int i = 0; i < 1920 * 1080 * 3; i += 3) {
		positionliteral[i] = (float)(opened[(i) + 18 + 2]) / 255.0; //R texture is in BGR, but RGB is needed; 
		positionliteral[i+1] = (float)(opened[(i) + 18 + 1]) / 255.0; //G
		positionliteral[i+2] = (float)(opened[(i) + 18]) / 255.0; //B
	}
	struct vec3 rgb = {};
	rgb.initvec3(1920, 1080);
	for (int i = 0; i < 1920 * 1080 * 3; i += 3) {
		rgb.z[i/3] = positionliteral[i];
		rgb.y[i/3] = positionliteral[i + 1];
		rgb.x[i/3] = positionliteral[i + 2];
	}
	struct sampler2D channel0 = {};
	channel0.init(1920, 1080, positionliteral);
	//rgb = texture2d(channel0, 0.1, 0.1);


	
	tex2.TEXTURE_INITIALIZATION(1920, 1080, positionliteral);
	//tex2.FLIP_TEXTURE_VERT(); */

	SetThreadExecutionState(ES_DISPLAY_REQUIRED | ES_CONTINUOUS); // sets the state so that the display is kept on and that this state is continous
	struct vec3 color = {};
	color.initvec3(frame.width, frame.height);
	float texturebuffer[3];
	fvec3* buffer_a = new fvec3[((UWIDTH * UHEIGHT) / 8) + 2]{ f(0.0), f(0.0), f(0.0) };
	std::pair<float, float> ball = { 0.2, 0.0 }; // coords
	std::pair<float, float> velocity = { 200.047, -0.012 }; //xy velocity
	std::pair<float, float> acceleration = { 0.0, 0.0 };

	std::pair<float, float> ball2 = { 0.25, 0.0 }; // coords
	std::pair<float, float> velocity2 = { 0.0241, 0.012 }; //xy velocity
	std::pair<float, float> acceleration2 = { 0.0, 0.0 };
	float dt = 0.01f;
	while (!quit) {
		static MSG message = { 0 };
		while (PeekMessage(&message, NULL, 0, 0, PM_REMOVE)) { DispatchMessage(&message); }
		//romat.update(time * 0.0025);
		POINT cursor;
		GetCursorPos(&cursor);
		float cursor_x = (float)cursor.x / UWIDTH - 0.5f;
		float cursor_y = (float)(1.0 - cursor.y) / UHEIGHT + 0.5f;
		float bd = 0.0;
		if (GetKeyState(VK_LMENU) & 0x8000) {
			zoom = 0.98f * zoom;
			bd = 1.0;
		}
		if (GetKeyState(VK_RMENU) & 0x8000) {
			zoom = 1.02f * zoom;
		}
		for (int i = 0; i < frame.width * frame.height; i += 8) { // initial shader
			__m256 vec_time = _mm256_set1_ps(time);
			__m256 vec_u = _mm256_load_ps((float*)u + i); // uv.x X component of uv
			__m256 vec_v = _mm256_load_ps((float*)v + i);
			fvec2 uv = { vec_u, vec_v };
			uv = add(uv, fvec2{ f(-0.5), f(-0.5) });
			uv.x = _mm256_mul_ps(uv.x, f((float)UWIDTH / UHEIGHT));
			__m256 vec_ut = _mm256_sub_ps(vec_u, const_vec);
			__m256 vec_vt = _mm256_sub_ps(vec_v, const_vec);
			vec_ut = _mm256_mul_ps(vec_ut, f(25.0));
			vec_vt = _mm256_mul_ps(vec_vt, f(25.0));
			fvec3 col = {f(0.0), f(0.0), f(0.0)};
			fvec2 obj = { f(ball.first), f(ball.second)};
			__m256 dist = distance(uv, obj);
			dist = clamp2d(dist, f(0.0), f(1.0));
			dist = _mm256_sub_ps(f(1.0), dist);
			dist = step2d(f(0.97), dist);
			fvec2 obj2 = { f(ball2.first), f(ball2.second) };
			__m256 dist2 = distance(uv, obj2);
			dist2 = clamp2d(dist2, f(0.0), f(1.0));
			dist2 = _mm256_sub_ps(f(1.0), dist2);
			dist2 = step2d(f(0.97), dist2);

			//dist = _mm256_mul_ps(dist, f(sqrt((velocity.first * velocity.first) + velocity.second * velocity.second)));
			//dist2 = _mm256_mul_ps(dist2, f(sqrt((velocity2.first * velocity2.first) + velocity2.second * velocity2.second)));
			col = fvec3{ dist,f(0.0),dist2};

			//fvec2 mouse = { f(cursor_x), f(cursor_y) };

			//fvec3 ro = fvec3{ f(0.0), f(0.0), f(-1.0) }; // ray origin
			//fvec3 rd = fvec3{ uv.x, uv.y, f(1.0) }; // ray direction
			//fvec3 p = fvec3{ f(sin(time)), f(0.0), f(cos(time) * 2.0 + 2.0)};

			//float32 d = _mm256_div_ps(length(cross(sub(p, ro), rd)), length(rd));
			//d = smoothstep2d(f(0.1), f(0.05), d);
			//col = { d,d,d };
			//col = bilinear(channel0, fvec2{_mm256_add_ps(vec_u, _mm256_mul_ps(_mm256_sin_ps(_mm256_add_ps(_mm256_mul_ps(vec_v, f(20.0)), f(time))),f(0.1))), vec_v});

			//col = clamp(col, f(0.0), f(1.0));
			//buffer_a[i / 8] = col;
			//col = mul(col, texture(channel0, fvec2{vec_u, vec_v}));
			_mm256_store_ps((float*)color.x + i, col.r);
			_mm256_store_ps((float*)color.y + i, col.g);
			_mm256_store_ps((float*)color.z + i, col.b);
			//col[i] = u[i%frame.width];
		}

		
		float dx = ball2.first - ball.first;
		float dy = ball2.second - ball.second;

		float distSquared = dx * dx + dy * dy;
		float distance = sqrt(distSquared);
		// Prevent division by zero
		if (distance != 0.0f) {
			// Normalize direction vector
			float nx = dx / distance;
			float ny = dy / distance;

			// Choose acceleration magnitude (e.g., inverse of distance or inverse-square)
			float accelMag = (1.0f / distSquared) * sigmoid_bonding(distance); // or 1.0f / distance;
			//printf("%f\n", accelMag);
			// Set acceleration vectors in opposite directions
			acceleration.first = nx * accelMag;
			acceleration.second = ny * accelMag;

			acceleration2.first = -acceleration.first;
			acceleration2.second = -acceleration.second;
		}
		//update accel		
		velocity.first += acceleration.first * dt;
		velocity.second += acceleration.second * dt;
		ball.first += velocity.first * dt;
		ball.second += velocity.second * dt;

		float loss = 0.5;
		if (ball.first < -0.8888f) {
			ball.first = -0.8888f;
			velocity.first = -velocity.first * loss;
		}
		if (ball.first > 0.8888f) {
			ball.first = 0.8888f;
			velocity.first = -velocity.first * loss;
		}
		if (ball.second < -0.5f) {
			ball.second = -0.5f;
			velocity.second = -velocity.second * loss;
		}
		if (ball.second > 0.5f) {
			ball.second = 0.5f;
			velocity.second = -velocity.second * loss;
		}
		velocity2.first += acceleration2.first * dt;
		velocity2.second += acceleration2.second * dt;
		ball2.first += velocity2.first * dt;
		ball2.second += velocity2.second * dt;

		if (ball2.first < -0.8888f) {
			ball2.first = -0.8888f;
			velocity2.first = -velocity2.first * loss;
		}
		if (ball2.first > 0.8888f) {
			ball2.first = 0.8888f;
			velocity2.first = -velocity2.first * loss;
		}
		if (ball2.second < -0.5f) {
			ball2.second = -0.5f;
			velocity2.second = -velocity2.second * loss;
		}
		if (ball2.second > 0.5f) {
			ball2.second = 0.5f;
			velocity2.second = -velocity2.second * loss;
		}
		
		//MAP_TEXTURE(&tex1, color, frame.width, frame.height, u, v);
		//MAP_TEXTURE(&tex2, color, frame.width, frame.height, u, v);
		for (int i = 0; i < frame.width * frame.height; i += 8) {//denormalization
			//int out = (int)(col[i] * 255.0);
			//frame.pixels[i] = out << 16 | out << 8 | out; 
			__m256 vec_colr = _mm256_load_ps(color.x + i);
			__m256 vec_colg = _mm256_load_ps(color.y + i);
			__m256 vec_colb = _mm256_load_ps(color.z + i);

			vec_colr = clamp2d(vec_colr, f(0.0), f(1.0)); vec_colg = clamp2d(vec_colg, f(0.0), f(1.0)); vec_colb = clamp2d(vec_colb, f(0.0), f(1.0));
			__m256 mul = f(255.0f);

			vec_colr = _mm256_mul_ps(vec_colr, mul);
			vec_colg = _mm256_mul_ps(vec_colg, mul);
			vec_colb = _mm256_mul_ps(vec_colb, mul);

			__m256i vec_b = _mm256_cvtps_epi32(vec_colb);
			__m256i vec_g = _mm256_cvtps_epi32(vec_colg);
			__m256i vec_r = _mm256_cvtps_epi32(vec_colr);

			vec_g = _mm256_slli_epi32(vec_g, 8); // shift left 8
			vec_r = _mm256_slli_epi32(vec_r, 16); // shift left 16
			vec_b = _mm256_add_epi32(vec_b, vec_g);
			vec_b = _mm256_add_epi32(vec_b, vec_r);

			_mm256_storeu_epi32(frame.pixels + i, vec_b);


		}// colorspace
		static unsigned int p = 0;
		p++;
		//frame.pixels[(p++) % (frame.width * frame.height)] = 255;
		//frame.pixels[Rand32() % (frame.width * frame.height)] = 0;
		// delete em
		if (p % 100 == 0) { printf("\r<%i>", p); } // \r return carriage moves the cursor back to 0
		InvalidateRect(window_handle, NULL, FALSE);
		UpdateWindow(window_handle);
		time += dt;
	}
	return 0;
}

LRESULT CALLBACK WindowProcessMessage(HWND window_handle, UINT message, WPARAM wParam, LPARAM lParam) {
	switch (message) {
	case WM_QUIT:
	case WM_DESTROY: {
		quit = true;
	} break;
	case WM_PAINT: {
		static PAINTSTRUCT paint;
		static HDC device_context;
		device_context = BeginPaint(window_handle, &paint);
		BitBlt(device_context,
			paint.rcPaint.left, paint.rcPaint.top,
			paint.rcPaint.right - paint.rcPaint.left, paint.rcPaint.bottom - paint.rcPaint.top,
			frame_device_context,
			paint.rcPaint.left, paint.rcPaint.top,
			SRCCOPY);
		EndPaint(window_handle, &paint);
	} break;

	case WM_SIZE: {
		frame_bitmap_info.bmiHeader.biWidth = LOWORD(lParam);
		frame_bitmap_info.bmiHeader.biHeight = HIWORD(lParam);

		if (frame_bitmap) DeleteObject(frame_bitmap);
		frame_bitmap = CreateDIBSection(NULL, &frame_bitmap_info, DIB_RGB_COLORS, (void**)&frame.pixels, 0, 0);
		SelectObject(frame_device_context, frame_bitmap);

		frame.width = LOWORD(lParam);
		frame.height = HIWORD(lParam);
	} break;

	default: {
		return DefWindowProc(window_handle, message, wParam, lParam);
	}
	}
	return 0;
}

























































































void character_map_printf(int index, bool * fragColor, char * buffer) {
	char out[] = { ' ',0b11011111,0b11011100, 0b11011011}; // none, top half, bottom half, full
	for (int i = 0; i < w; i++) {
		int indexb = (fragColor[index + i] << 1) | (fragColor[index + i + w]);
		buffer[i] = out[indexb];
	}
	//printf("%c%c%c%c%c", out[fragColor[index]], out[fragColor[index+1]], out[fragColor[index + 2]], out[fragColor[index + 3]], out[fragColor[index + 4]]);
	printf(buffer);
	return;
}

void FILL_FRAME(bool* data, char * buffer) {
	int width = w;
	int height = h;
	for (int i = height-1; i >= 0; i -= 2) { // 50 indexes 49-0 inclusive 2 lines rendered per map printf call
		for (int j = 0; j < width; j += w) { // 100 indexes 0-99 inclusive (printing 5 chars at a time so increment by 5)
			//j + (i * width) is the index
			character_map_printf(j + (i * width), data, buffer);
		}
		printf("\n");
	}
}

 void plot_line(int x0, int y0, int x1, int y1, bool * fragColor)
{
	int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
	int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
	int err = dx + dy, e2; /* error value e_xy */

	for (;;) {  /* loop */
		fragColor[x0 + (y0 * w)] = 1;
		if (x0 == x1 && y0 == y1) break;
		e2 = 2 * err;
		if (e2 >= dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
		if (e2 <= dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
	}
}

 void plot_triangle(int x0, int y0, int x1, int y1, int x2, int y2, bool * fragColor) {
	 plot_line(x0, y0, x1, y1, fragColor);
	 plot_line(x1, y1, x2, y2, fragColor);
	 plot_line(x0, y0, x2, y2, fragColor);
}

 int rotation_x(int x0, int y0, float t) {
	 int width = w;
	 int height = h;
	 float m[] = { cos(t),-sin(t),
						  sin(t), cos(t) }; // rotation matrix
	 float x = (float)x0 / width; // normalization of coordinates
	 float y = (float)y0 / height;
	 int ex0 = (int)(((m[0] * x) + (m[1] * y)) * width);
	 return ex0 % w;
}	

 int rotation_y(int x0, int y0, float t) {
	 int width = w;
	 int height = h;
	 float m[] = { cos(t),-sin(t),
						  sin(t), cos(t) }; // rotation matrix
	 float x = (float)x0 / width; // normalization of coordinates
	 float y = (float)y0 / height;
	 int ey0 = (int)(((m[2] * x) + (m[3] * y)) * height);
	 return ey0 % h;
 }

 void CLEAR_FRAME(bool* fragColor) {
	 for (int i = 0; i < w*h; i++) {
		 fragColor[i] = 0;
	 }
 }

 void
	 raster_circle(int x0, int y0, int radius, bool* fragColor)
 {
	 int f = 1 - radius;
	 int ddF_x = 1;
	 int ddF_y = -2 * radius;
	 int x = 0;
	 int y = radius;

	 fragColor[x0+ ((y0 + radius) * w)] = 1;
	 fragColor[x0+ ((y0 - radius) * w)] = 1;
	 fragColor[x0 + radius+ (y0 * w)] = 1;
	 fragColor[x0 - radius+ (y0 * w)] = 1;
	 while (x < y)
	 {
		 // ddF_x == 2 * x + 1;
		 // ddF_y == -2 * y;
		 // f == x*x + y*y - radius*radius + 2*x - y + 1;
		 if (f >= 0)
		 {
			 y--;
			 ddF_y += 2;
			 f += ddF_y;
		 }
		 x++;
		 ddF_x += 2;
		 f += ddF_x;
		 fragColor[x0 + x + ((y0 + y) * w)] = 1;
		 fragColor[x0 - x + ((y0 + y) * w)] = 1;
		 fragColor[x0 + x + ((y0 - y) * w)] = 1;
		 fragColor[x0 - x + ((y0 - y) * w)] = 1;
		 fragColor[x0 + y + ((y0 + x) * w)] = 1;
		 fragColor[x0 - y + ((y0 + x) * w)] = 1;
		 fragColor[x0 + y + ((y0 - x) * w)] = 1;
		 fragColor[x0 - y + ((y0 - x) * w)] = 1;
	 }
 }

/* EOF */

int main() {
	bool fragColor[w*h] = {0}; // 5000 element vec1 (0 - 4999)
	//bool* fragColor = (bool*)malloc(w * h);

	WinMain(NULL, NULL, NULL, NULL);

	for (int i = 0; i < w*h; i++) {
		fragColor[i] = 0;
	}

	float t = 0.0; // theta/angle inputted to rotation matrix


	fragColor[99] = 1;
	fragColor[4998] = 1;

	int x0 = 100;
	int y0 = 100;
	int x1 = 100;
	int y1 = 200;
	int x2 = 250;
	int y2 = 100;

	int x3 = 299;
	int y3 = 200;
	
	int cx = (x0+x3)/2; // the centering of uv
	int cy = (y0+y3)/2;

	int x4 = 200;
	int y4 = 200;
	char* buffer = (char*)malloc(sizeof(char) * w);
	/*
	for (;;) {
		int ex0 = rotation_x(x0 - cx, y0 - cy, t); // minus to center coords before tranformation
		int ey0 = rotation_y(x0 - cx, y0 - cy, t);
		int ex1 = rotation_x(x1 - cx, y1 - cy, t);
		int ey1 = rotation_y(x1 - cx, y1 - cy, t);
		int ex2 = rotation_x(x2 - cx, y2 - cy, t);
		int ey2 = rotation_y(x2 - cx, y2 - cy, t);
		int ex3 = rotation_x(x3 - cx, y3 - cy, t);
		int ey3 = rotation_y(x3 - cx, y3 - cy, t);
		plot_triangle(ex0 + cx, ey0 + cy, ex1 + cx, ey1 + cy, ex2 + cx, ey2 + cy, fragColor); // adding back to coords after tranformation
		plot_triangle(ex1 + cx, ey1 + cy, ex2 + cx, ey2 + cy, ex3 + cx, ey3 + cy, fragColor);

		raster_circle(x4, y4, (int)(abs(sin(t)) * 50), fragColor);
		system("cls");
		FILL_FRAME(fragColor, buffer);
		CLEAR_FRAME(fragColor);
		t += 0.03;
		Sleep(1);
	}
	*/
	
	
}
