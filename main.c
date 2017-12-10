#include "b_cos.h"

// Current Ubuntu LTS provides both libpng 1.2 and 1.6, but many of the other development libraries still require 1.2, and
// libpng16-dev can't be installed in parallel with libpng12-dev. This code uses 1.2.
#include "png.h"

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdnoreturn.h>
#include <string.h>

// #define BENCHMARK

#ifdef BENCHMARK
#	include <sys/time.h>
#endif

#define arraysize(a) (sizeof(a) / sizeof(*(a)))

static const float *const b_cos_ifacs[] =
{
	b_cos_ifac0,
	b_cos_ifac1,
	b_cos_ifac2,
	b_cos_ifac3,
	b_cos_ifac4,
	b_cos_ifac5,
	b_cos_ifac6,
	b_cos_ifac7,
	b_cos_ifac8,
	b_cos_ifac9,
};

typedef uint8_t b_cos_2d_src_px; // Too many underscores?
typedef uint8_t b_cos_2d_dest_px;

struct b_cos_2d
{
	png_uint_32 src_width;
	png_uint_32 src_height;
	const b_cos_2d_src_px *src;
	png_uint_32 dest_width;
	png_uint_32 dest_height;
	unsigned radius;

	float scale_fac;
	float *xfac;
	float *x_buffer;
	unsigned x_buffer_start, x_buffer_offset;

	float *yfac0;
	float *y_buf_expand;
	float *y_sum_fac;
	float *y_sum;
	unsigned y;

	unsigned src_y0, y_src_offset, y_src_step, y_src_jump;
};

static inline b_cos_2d_dest_px _clip(float x)
{
	x += 0.5f;
	if(x < 0)
		return 0;
	if(x > 255)
		return 255;
	return (b_cos_2d_dest_px)x;
}

static void _edge_x(
	const struct b_cos_2d *self,
	const float *xfac,
	const b_cos_2d_src_px *src_row,
	int src_x,
	float *sum_fac)
{
	for(unsigned ch = 0; ch != 4; ++ch)
		sum_fac[ch] = 0;

	for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
	{
		int lo = src_x + i - self->radius;
		if(lo < 0)
			lo = 0;
		if((unsigned)lo >= self->src_width)
			lo = self->src_width - 1;
		assert(lo >= 0);
		assert((unsigned)lo < self->src_width);
		for(unsigned ch = 0; ch != 4; ++ch)
			sum_fac[ch] += src_row[lo * 4 + ch] * xfac[i];
	}
}

static void _expand_x(
	const struct b_cos_2d *self,
	const b_cos_2d_src_px *src_row,
	float *dest_row)
{
	// printf("%zd\n", (src_row - self->src) / (self->src_width * 4));

	const float *xfac0 = self->xfac;

	float sum_fac[4];
	_edge_x(self, xfac0, src_row, 0, sum_fac);

	unsigned
		src_x0 = 0,
		src_offset = 0,
		src_step = self->src_width % self->dest_width,
		src_jump = self->src_width / self->dest_width;

	// Not at all friendly to SIMD, this bit right here.
	for(unsigned x = 0; x != self->dest_width; ++x)
	{
		unsigned src_x1 = src_x0 + src_jump;
		src_offset += src_step;
		if(src_offset >= self->dest_width)
		{
			src_offset -= self->dest_width;
			++src_x1;
		}

		assert(src_x1 == (x + 1) * self->src_width / self->dest_width);
		assert(src_x0 < self->src_width);
		assert(src_x1 <= self->src_width);

		int isum[4] = {0};

		for(int i = (int)src_x0 - (int)self->radius; i < (int)src_x1 - (int)self->radius; ++i)
		{
			int src_x = i;
			if(src_x < 0)
				src_x = 0;
			assert(src_x >= 0);
			assert((unsigned)src_x < self->src_width);
			for(unsigned ch = 0; ch != 4; ++ch)
				isum[ch] += src_row[src_x * 4 + ch];
		}

		float fsum[4];
		for(unsigned ch = 0; ch != 4; ++ch)
			fsum[ch] = isum[ch] - sum_fac[ch];

		xfac0 += B_COS_EDGE(self->radius);
		_edge_x(self, xfac0, src_row, src_x1, sum_fac);

		for(unsigned ch = 0; ch != 4; ++ch)
			dest_row[x * 4 + ch] = fsum[ch] + sum_fac[ch];
		// dest_row[x * 4 + 3] = src_row[src_x0 * 4 + 3];

		src_x0 = src_x1;
	}
}

void _clear_x(float *row, unsigned dest_width)
{
	for(unsigned i = 0; i != dest_width * 4; ++i)
		row[i] = 0;
}

static inline const b_cos_2d_src_px *_y_seek(const struct b_cos_2d *self, unsigned y)
{
	return self->src + self->src_width * 4 * y;
}

static float *_x_buffer_cache(const struct b_cos_2d *self, unsigned y)
{
	// TODO Replace these.
	assert(y >= self->x_buffer_start);
	assert(y < self->x_buffer_start + B_COS_EDGE(self->radius));
	return
		self->x_buffer +
		self->dest_width * 4 * ((self->x_buffer_offset + y - self->x_buffer_start) % B_COS_EDGE(self->radius));
}

static void _edge_y(const struct b_cos_2d *self, int src_y, int dest_y)
{
	b_cos_edge(self->yfac0, self->radius, b_cos_ifacs[self->radius], self->src_height, self->dest_height, dest_y);

	for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
		self->yfac0[i] -= 0.5f;

	for(unsigned x = 0; x != self->dest_width * 4; ++x)
		self->y_sum_fac[x] = 0;

	for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
	{
		int lo = src_y + i - self->radius;
		if(lo < 0)
			lo = 0;
		if((unsigned)lo >= self->src_height)
			lo = self->src_height - 1;
		assert(lo >= 0);
		assert((unsigned)lo < self->src_height);

		// _expand_x(self, _y_seek(self, lo), self->y_buf_expand);
		// const float *x_src = self->y_buf_expand;

		const float *x_src = _x_buffer_cache(self, lo);

		for(unsigned x = 0; x != self->dest_width * 4; ++x)
			self->y_sum_fac[x] += x_src[x] * self->yfac0[i];
	}
}

void _mul_add_x(unsigned dest_width, float *src, float *dest, float fac)
{
	for(unsigned i = 0; i != dest_width * 4; ++i)
		dest[i] += src[i] * fac;
}

void b_cos_2d_free(struct b_cos_2d *self)
{
	free(self->xfac);
	free(self->x_buffer);
	free(self->yfac0);
	free(self->y_buf_expand);
	free(self->y_sum_fac);
	free(self->y_sum);
}

bool b_cos_2d_init(struct b_cos_2d *self)
{
	assert(self->src_width);
	assert(self->src_height);

	// TODO: Take repeat into account.
	self->xfac = malloc(B_COS_EDGE(self->radius) * (self->dest_width + 1) * sizeof(float));
	self->x_buffer = malloc(self->dest_width * 4 * sizeof(float) * B_COS_EDGE(self->radius));
	self->yfac0 = malloc(B_COS_EDGE(self->radius) * sizeof(float));
	self->y_buf_expand = malloc(self->dest_width * 4 * sizeof(float));
	self->y_sum_fac = malloc(self->dest_width * 4 * sizeof(float));
	self->y_sum = malloc(self->dest_width * 4 * sizeof(float));

	if(!self->xfac || !self->x_buffer || !self->yfac0 || !self->y_buf_expand || !self->y_sum_fac || !self->y_sum)
	{
		b_cos_2d_free(self);
		return false;
	}

	for(unsigned x = 0; x != self->dest_width + 1; ++x)
	{
		float *xfac0 = self->xfac + x * B_COS_EDGE(self->radius);
		b_cos_edge(xfac0, self->radius, b_cos_ifacs[self->radius], self->src_width, self->dest_width, x);
		for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
			xfac0[i] += 0.5;
	}

	self->scale_fac = /* (65535 / 255) * */ (float)self->dest_width * self->dest_height / (self->src_width * self->src_height);

	self->x_buffer_start = 0;
	self->x_buffer_offset = 0;
	{
		const b_cos_2d_src_px *src_row = self->src;
		float *dest_row = self->x_buffer;
		for(unsigned y = 0; y != B_COS_EDGE(self->radius); ++y)
		{
			if(y < self->dest_height)
			{
				_expand_x(self, src_row, dest_row);
				src_row += self->src_width * 4;
			}
			else
			{
				memcpy(dest_row, dest_row - self->dest_width * 4, sizeof(float) * 4 * self->dest_width);
			}

			dest_row += self->dest_width * 4;
		}
	}

	_edge_y(self, 0, 0);

	self->src_y0 = 0;
	self->y_src_offset = 0;
	self->y_src_step = self->src_height % self->dest_height;
	self->y_src_jump = self->src_height / self->dest_height;
	self->y = 0;

	return true;
}

void b_cos_2d_row(struct b_cos_2d *self, b_cos_2d_dest_px *dest_row)
{
	unsigned src_y1 = self->src_y0 + self->y_src_jump;
	self->y_src_offset += self->y_src_step;
	if(self->y_src_offset >= self->dest_height)
	{
		self->y_src_offset -= self->dest_height;
		++src_y1;
	}

	assert(src_y1 == (self->y + 1) * self->src_height / self->dest_height);
	assert(self->src_y0 < self->src_height);
	assert(src_y1 <= self->src_height);

	for(unsigned x = 0; x != self->dest_width * 4; ++x)
		self->y_sum[x] = 0;

	// Slide the row cache forwards.
	// This could be refactored with what's in b_cos_2d_init...?
	// TODO: Clean up the logic here.
	{
		int new_x_buffer_start = (int)src_y1 - (int)self->radius;
		if(new_x_buffer_start < 0)
			new_x_buffer_start = 0;
		if((unsigned)new_x_buffer_start >= self->src_height)
			new_x_buffer_start = self->src_height - 1;
		assert((unsigned)new_x_buffer_start >= self->x_buffer_start);
		unsigned advance = new_x_buffer_start - self->x_buffer_start;
		if(advance >= B_COS_EDGE(self->radius))
		{
			advance = B_COS_EDGE(self->radius);
			self->x_buffer_offset = 0;
		}
		else
		{
			self->x_buffer_offset = (self->x_buffer_offset + advance) % B_COS_EDGE(self->radius);
		}

		self->x_buffer_start = new_x_buffer_start;
		unsigned new_end = new_x_buffer_start + B_COS_EDGE(self->radius);

		for(unsigned y = new_end - advance; y != new_end; ++y)
			_expand_x(self, _y_seek(self, y < self->src_height ? y : self->src_height - 1), _x_buffer_cache(self, y));
	}

	// TODO: Merge this with _edge_y.
	for(int i = (int)self->src_y0 + (int)self->radius + 1; i < (int)src_y1 + (int)self->radius + 1; ++i)
	{
		unsigned src_y = i;
		if(src_y >= self->src_height)
			src_y = self->src_height - 1;
		assert(src_y >= 0);
		assert((unsigned)src_y < self->src_height);

		const float *row;
		if(src_y >= self->x_buffer_start)
		{
			row = _x_buffer_cache(self, src_y);
		}
		else
		{
			_expand_x(self, _y_seek(self, src_y), self->y_buf_expand);
			row = self->y_buf_expand;
		}

		for(unsigned x = 0; x != self->dest_width * 4; ++x)
			self->y_sum[x] += row[x];
	}

	for(unsigned x = 0; x != self->dest_width * 4; ++x)
		self->y_sum[x] -= self->y_sum_fac[x];

	_edge_y(self, src_y1, self->y + 1);

	for(unsigned x = 0; x != self->dest_width * 4; ++x)
		dest_row[x] = _clip(self->scale_fac * (self->y_sum[x] + self->y_sum_fac[x]));

//	for(unsigned x = 0; x != self->dest_width; ++x)
//		dest_row[x * 4 + 3] = 0xffff;

	self->src_y0 = src_y1;
	++self->y;
}

unsigned _strtou(const char *nptr, const char **endptr)
{
	unsigned long result = strtoul(nptr, (char **)endptr, 10);
	if(result > UINT_MAX)
		*endptr = nptr;
	return result;
}

#ifdef BENCHMARK
double _double_time()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}
#endif

static int _error(const char *prefix, const char *message)
{
	fprintf(stderr, "%s: %s\n", prefix, message);
	return EXIT_FAILURE;
}

static void _error_errno(const char *prefix, int error)
{
	_error(prefix, strerror(error));
}

static void _error_nomem(char **argv)
{
	_error_errno(*argv, ENOMEM);
}

noreturn void _png_raise(png_structp png)
{
	longjmp(png_jmpbuf(png), 1);
}

noreturn void _png_raise_nomem(png_structp png, char **argv)
{
	_error_nomem(argv);
	_png_raise(png);
}

noreturn void _png_raise_errno(png_structp png, const char *prefix)
{
	_error_errno(prefix, errno);
	_png_raise(png);
}

int main(int argc, char **argv)
{
	if(argc != 6)
	{
		printf("Usage: %s src.png dest.png dest_width dest_height radius\n", *argv);
		return EXIT_FAILURE;
	}

	const char *src_png_name = argv[1];
	const char *dest_png_name = argv[2];
	const char *str_end;

	struct b_cos_2d b_cos;
	b_cos.dest_width = _strtou(argv[3], &str_end);
	if(*str_end || !b_cos.dest_width)
		return _error(*argv, "dest_width must an integer greater than 0.");

	b_cos.dest_height = _strtou(argv[4], &str_end);
	if(*str_end || !b_cos.dest_height)
		return _error(*argv, "dest_height must an integer greater than 0.");

	b_cos.radius = _strtou(argv[5], &str_end);
	if(*str_end || b_cos.radius >= arraysize(b_cos_ifacs))
	{
		fprintf(stderr, "%s: radius must an integer between 0 and %zd.\n", *argv, arraysize(b_cos_ifacs) - 1);
		return EXIT_FAILURE;
	}

	// Rescaling an image can work from a moving window of several adjacent rows to output rows in the resampled image. If PNG
	// always supported row-at-a-time loading, then there wouldn't be any need to load the entire image into memory. But PNG
	// supports interlacing, so to handle that, let's just use png_read_image. png_image is another possibility, but current
	// Ubuntu LTS doesn't have that.

	int result = EXIT_FAILURE;

	png_structp src_png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!src_png)
	{
		_error_nomem(argv);
	}
	else
	{
		FILE *src_fp = NULL;
		png_byte **row_pointers = NULL;
		png_infop src_info = NULL;

		if(!setjmp(png_jmpbuf(src_png)))
		{
			src_info = png_create_info_struct(src_png);
			if(!src_info)
				_png_raise_nomem(src_png, argv);

			src_fp = fopen(src_png_name, "rb");
			if(!src_fp)
				_png_raise_errno(src_png, src_png_name);
			png_init_io(src_png, src_fp);

			// Using the low-level libpng APIs here because the high-level APIs free the image memory on shut down.
			png_read_info(src_png, src_info);
			int bit_depth, color_type;
			png_get_IHDR(
				src_png,
				src_info,
				&b_cos.src_width,
				&b_cos.src_height,
				&bit_depth,
				&color_type,
				NULL,
				NULL,
				NULL);

			// libpng 1.2 doesn't have png_set_expand_16().
			png_set_expand(src_png);
			if(bit_depth == 16)
				png_set_strip_16(src_png);
			png_color_8 *sig_bit;
			if(png_get_sBIT(src_png, src_info, &sig_bit))
				png_set_shift(src_png, sig_bit);
			if(bit_depth < 8)
				png_set_packing(src_png);
			if(!(color_type & PNG_COLOR_MASK_ALPHA))
				png_set_filler(src_png, 0xff, PNG_FILLER_AFTER);
			if((color_type & ~(int)PNG_COLOR_MASK_ALPHA) == PNG_COLOR_TYPE_GRAY)
				png_set_gray_to_rgb(src_png);

			b_cos.src = malloc(b_cos.src_width * b_cos.src_height * 4);
			if(!b_cos.src)
				_png_raise_nomem(src_png, argv);
			row_pointers = malloc(sizeof(png_byte *) * b_cos.src_height);
			if(!row_pointers)
				_png_raise_nomem(src_png, argv);

			for(unsigned y = 0; y != b_cos.src_height; ++y)
				row_pointers[y] = (uint8_t *)b_cos.src + b_cos.src_width * 4 * y;
			png_read_image(src_png, row_pointers);

			result = EXIT_SUCCESS;

			/*
			png_read_png(
				src_png,
				src_info,
				PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_EXPAND | PNG_TRANSFORM_SHIFT | PNG_TRANSFORM_GRAY_TO_RGB,
				NULL);

			png_get_rows(src_png, src_info);
			*/
		}

		png_destroy_read_struct(&src_png, &src_info, (png_infopp)NULL);
		free(row_pointers);
		fclose(src_fp);
	}

	if(result)
	{
	}
	else if(!b_cos_2d_init(&b_cos))
	{
		_error_nomem(argv);
	}
	else
	{
		result = EXIT_FAILURE;

		const unsigned channels = 4;

		png_structp dest_png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if(!dest_png)
		{
			_error_nomem(argv);
		}
		else
		{
			FILE *dest_fp = NULL;
			b_cos_2d_dest_px *dest_row = NULL;
			png_infop dest_info = NULL;

			if(!setjmp(png_jmpbuf(dest_png)))
			{
				dest_info = png_create_info_struct(dest_png);
				if(!dest_info)
					_png_raise_nomem(dest_png, argv);

#ifdef BENCHMARK
				dest_png_name = "/dev/null";
#endif
				dest_fp = fopen(dest_png_name, "wb");
				if(!dest_fp)
					_png_raise_errno(dest_png, dest_png_name);
				png_init_io(dest_png, dest_fp);

				dest_row = malloc(channels * b_cos.dest_width * sizeof(*dest_row));
				if(!dest_row)
					_png_raise_nomem(dest_png, argv);

				int color_type = PNG_COLOR_TYPE_RGB_ALPHA;
				// int color_type = src_image.format & PNG_FORMAT_FLAG_ALPHA ? PNG_COLOR_MASK_ALPHA : 0;
				// if(src_image.format & PNG_FORMAT_FLAG_COLOR)
				//	color_type |= PNG_COLOR_MASK_COLOR;
				png_set_IHDR(
					dest_png,
					dest_info,
					b_cos.dest_width,
					b_cos.dest_height,
					sizeof(*dest_row) * CHAR_BIT,
					color_type,
					PNG_INTERLACE_NONE,
					PNG_COMPRESSION_TYPE_DEFAULT,
					PNG_FILTER_TYPE_DEFAULT);

				// TODO: Inherit gamma/chromaticity/sRGB from source image.
				// png_set_sRGB_gAMA_and_cHRM(dest_png, dest_info, PNG_sRGB_INTENT_RELATIVE);
				// png_set_gAMA(dest_png, dest_info, 1.0);
				// png_set_cHRM(dest_png, dest_info, 0.3127, 0.3290, 0.6400, 0.3300, 0.3000, 0.6000, 0.1500, 0.0600);

				png_write_info(dest_png, dest_info);

#ifndef __BYTE_ORDER__
#	error #define __BYTE_ORDER__ to 1234 (little endian) or 4321 (big endian).
#elif __BYTE_ORDER__ == 1234 || __BYTE_ORDER__ == 3412
				png_set_swap(dest_png);
#endif

#ifdef BENCHMARK
				double then = _double_time();
#endif

				png_write_info(dest_png, dest_info);

				for(unsigned y = 0; y != b_cos.dest_height; ++y)
				{
					b_cos_2d_row(&b_cos, dest_row);

#ifndef BENCHMARK
					// PNG compression takes up most of the time here.
					png_write_rows(dest_png, (png_bytepp)&dest_row, 1);
#endif
				}

#ifdef BENCHMARK
				double elapsed = _double_time() - then;
				printf(
					"Elapsed time: %g seconds, %g Mpixels/second\n",
					elapsed,
					1e-6 * ((b_cos.src_width * b_cos.src_height) + (b_cos.dest_width * b_cos.dest_height)) / elapsed);
#endif
				png_write_end(dest_png, NULL);
				result = EXIT_SUCCESS;
			}

			png_destroy_write_struct(&dest_png, &dest_info);
			free(dest_row);
			fclose(dest_fp);

			if(dest_fp && result)
				remove(dest_png_name);
		}

		b_cos_2d_free(&b_cos);
	}

	free((void *)b_cos.src);
	return result;
}
