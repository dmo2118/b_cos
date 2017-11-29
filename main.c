// TODO: Test on macOS.

#include "b_cos.h"

#include "png.h"

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define arraysize(a) (sizeof(a) / sizeof(*(a)))

static int _error(const char *prefix, const char *message)
{
	fprintf(stderr, "%s: %s\n", prefix, message);
	return EXIT_FAILURE;
}

static int _png_error(const char *png_name, png_image *image)
{
	return _error(png_name, image->message);
}

static int _error_errno(const char *prefix, int error)
{
	return _error(prefix, strerror(error));
}

static inline png_uint_16 _clip(double x)
{
	if(x < 0)
		return 0;
	if(x > 65535)
		return 65535;
	return (png_uint_16)x;
}

static const float *const b_cos_ifacs[] =
{
	b_cos_ifac0,
	b_cos_ifac1,
	b_cos_ifac2,
	b_cos_ifac3,
	b_cos_ifac4
};

struct b_cos_2d
{
	unsigned src_width;
	unsigned src_height;
	const png_uint_16 *src;
	unsigned dest_width;
	unsigned dest_height;
	unsigned radius;

	float scale_fac;
	float *xfac;
	float *yfac0;
	float *y_buf_expand;
	float *y_sum_fac;
	float *y_sum;
	unsigned y;

	unsigned src_y0, y_src_offset, y_src_step, y_src_jump;
};

static void _edge_x(
	const struct b_cos_2d *self,
	const float *xfac,
	const png_uint_16 *src_row,
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
	const png_uint_16 *src_row,
	float *dest_row)
{
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

inline const png_uint_16 *_y_seek(const struct b_cos_2d *self, unsigned y)
{
	return self->src + self->src_width * 4 * y;
}

static void _edge_y(const struct b_cos_2d *self, int src_y, int dest_y)
{
	b_cos_edge(self->yfac0, self->radius, b_cos_ifacs[self->radius], self->src_height, self->dest_height, dest_y);

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
		_expand_x(self, _y_seek(self, lo), self->y_buf_expand);
		for(unsigned x = 0; x != self->dest_width * 4; ++x)
			self->y_sum_fac[x] += self->y_buf_expand[x] * (self->yfac0[i] + 0.5f);
	}
}

void _mul_add_x(unsigned dest_width, float *src, float *dest, float fac)
{
	for(unsigned i = 0; i != dest_width * 4; ++i)
		dest[i] += src[i] * fac;
}

bool b_cos_2d_init(struct b_cos_2d *self)
{
	// TODO: Take repeat into account.
	self->xfac = malloc(B_COS_EDGE(self->radius) * (self->dest_width + 1) * sizeof(float));

	for(unsigned x = 0; x != self->dest_width + 1; ++x)
	{
		float *xfac0 = self->xfac + x * B_COS_EDGE(self->radius);
		b_cos_edge(xfac0, self->radius, b_cos_ifacs[self->radius], self->src_width, self->dest_width, x);
		for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
			xfac0[i] += 0.5;
	}

	self->scale_fac = (float)self->dest_width * self->dest_height / (self->src_width * self->src_height);

	self->yfac0 = malloc(B_COS_EDGE(self->radius) * sizeof(float));

	self->y_buf_expand = malloc(self->dest_width * 4 * sizeof(float));
	self->y_sum_fac = malloc(self->dest_width * 4 * sizeof(float));
	self->y_sum = malloc(self->dest_width * 4 * sizeof(float));

	_edge_y(self, 0, 0);

	self->src_y0 = 0;
	self->y_src_offset = 0;
	self->y_src_step = self->src_height % self->dest_height;
	self->y_src_jump = self->src_height / self->dest_height;
	self->y = 0;

	return true;
}

#if 0
void b_cos_2d_row(struct b_cos_2d *self, png_uint_16 *dest_row)
{
#if 0
	if(y < src_image.height)
	{
		memcpy(dest_row, src_row, (src_image.width < self->dest_width ? src_image.width : self->dest_width) * 8);
		_expand_x(xfac, self->radius, src_image.width, src_row, self->dest_width, dest_row);
		src_row += src_image.width * PNG_IMAGE_PIXEL_CHANNELS(src_image.format);
	}
#endif

	unsigned src_y0 = self->y * self->src_height / self->dest_height; // TODO: Division's slow.
	b_cos_edge(self->yfac1, self->radius, b_cos_ifacs[self->radius], self->src_height, self->dest_height, self->y + 1);

	// unsigned ymags_len = gap ? result = EXIT_SUCCESS;B_COS_EDGE(self->radius) * 2 : src_y1 + self->radius + 1 - (src_y0 - self->radius);

	/*
	// TODO Trash, right?
	for(unsigned i = 0; i != ymags_len; ++i)
		ymags[i] = 1;
	for(unsigned i = 0; i != B_COS_EDGE(self->radius); ++i)
	{
		ymags[i] += yfac0[i] - 0.5f;
		ymags[i + ymags_len - self->radius - 1] -= yfac1[i] + 0.5f;
	}
	*/

	// TODO: Probably better to optimize expand_x first.
	// Next up: add_mul existing rows.

	const png_uint_16 *src_row = self->src + self->src_width * 4 * src_y0;

	ptrdiff_t stride = self->dest_width * 4;
#if 1
	_expand_x(self, src_row, self->y_buffer);
	for(unsigned i = 0; i != stride; ++i)
		self->y_buffer[i] *= (float)self->src_height / self->dest_height;
#endif

#if 0
	unsigned src_y1 = (self->y + 1) * self->src_height / self->dest_height;
	const png_uint_16 *src_row1 = src_buffer + src_image.width * 4 * src_y1;
	// bool gap = (int)src_y0 + (int)self->radius < (int)src_y1 - (int)self->radius;

	_expand_x(xfac, self->radius, src_image.width, src_row, self->dest_width, y_buffer + stride * 2);

	_expand_x(xfac, self->radius, src_image.width, src_row1, self->dest_width, y_buffer + stride * 3);

	_clear_x(y_buffer, self->dest_width);
	_mul_add_x(self->dest_width, y_buffer + stride * 2, y_buffer, -yfac0[0]);
	_mul_add_x(self->dest_width, y_buffer + stride * 3, y_buffer, yfac1[0]);
#endif

#if 0
	unsigned y_extent = gap ?
		while(y_buffer_end <)
#endif

	for(unsigned i = 0; i != self->dest_width; ++i)
	{
		dest_row[i * 4    ] = _clip(self->scale_fac * self->y_buffer[i * 4]);
		dest_row[i * 4 + 1] = _clip(self->scale_fac * self->y_buffer[i * 4 + 1]);
		dest_row[i * 4 + 2] = _clip(self->scale_fac * self->y_buffer[i * 4 + 2]);
		dest_row[i * 4 + 3] = _clip(self->scale_fac * self->y_buffer[i * 4 + 3]);
	}

	float *temp = self->yfac0;
	self->yfac0 = self->yfac1;
	self->yfac1 = temp;

	++self->y;
}
#endif

#if 1
void b_cos_2d_row(struct b_cos_2d *self, png_uint_16 *dest_row)
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

	for(int i = (int)self->src_y0 - (int)self->radius; i < (int)src_y1 - (int)self->radius; ++i)
	{
		int src_y = i;
		if(src_y < 0)
			src_y = 0;
		assert(src_y >= 0);
		assert((unsigned)src_y < self->src_height);
		_expand_x(self, _y_seek(self, src_y), self->y_buf_expand);
		for(unsigned x = 0; x != self->dest_width * 4; ++x)
			self->y_sum[x] += self->y_buf_expand[x];
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
#endif

int main(int argc, char **argv)
{
	if(argc != 6)
	{
		printf("Usage: %s src.png dest.png dest_width dest_height radius\n", *argv);
		return EXIT_FAILURE;
	}

	const char *src_png_name = argv[1];
	const char *dest_png_name = argv[2];
	char *str_end;

	struct b_cos_2d b_cos;
	b_cos.dest_width = strtoul(argv[3], &str_end, 10);
	if(*str_end || !b_cos.dest_width)
		return _error(*argv, "dest_width must an integer greater than 0\n");

	b_cos.dest_height = strtoul(argv[4], &str_end, 10);
	if(*str_end || !b_cos.dest_height)
		return _error(*argv, "dest_height must an integer greater than 0\n");

	b_cos.radius = strtoul(argv[5], &str_end, 10);
	if(*str_end || b_cos.radius >= arraysize(b_cos_ifacs))
	{
		fprintf(stderr, "%s: radius must an integer between 0 and %zd\n", *argv, arraysize(b_cos_ifacs) - 1);
		return EXIT_FAILURE;
	}

	// Rescaling an image can work from a moving window of several adjacent rows to output rows in the resampled image. If PNG
	// always supported row-at-a-time loading, then there wouldn't be any need to load the entire image into memory. But PNG
	// supports interlacing, so to handle that, let's just use png_image.

	png_image src_image = {.version = PNG_IMAGE_VERSION};
	if(!png_image_begin_read_from_file(&src_image, src_png_name))
		return _png_error(src_png_name, &src_image);

	// src_image.format = src_image.format & (PNG_FORMAT_FLAG_ALPHA | PNG_FORMAT_FLAG_COLOR) | PNG_FORMAT_FLAG_LINEAR;
	src_image.format = PNG_FORMAT_FLAG_ALPHA | PNG_FORMAT_FLAG_COLOR | PNG_FORMAT_FLAG_LINEAR;
	// src_image.flags |= PNG_IMAGE_FLAG_16BIT_sRGB;

	png_uint_16 *src_buffer = malloc(PNG_IMAGE_SIZE(src_image));
	if(!src_buffer)
	{
		png_image_free(&src_image);
		return _error_errno(*argv, ENOMEM);
	}

	if(!png_image_finish_read(&src_image, 0, src_buffer, PNG_IMAGE_ROW_STRIDE(src_image), NULL))
	{
		free(src_buffer);
		png_image_free(&src_image);
		return _png_error(src_png_name, &src_image);
	}

	b_cos.src = src_buffer;
	b_cos.src_width = src_image.width;
	b_cos.src_height = src_image.height;

	int result = EXIT_FAILURE;

	// TODO: Unchecked unfreed malloc().
	uint16_t *dest_row = malloc(PNG_IMAGE_PIXEL_CHANNELS(src_image.format) * b_cos.dest_width * sizeof(*dest_row));

	FILE *dest_fp = fopen(dest_png_name, "wb");
	if(!dest_fp)
	{
		_error(dest_png_name, strerror(errno));
	}
	else
	{
		// Some apps can't handle 16-bit PNG w/ gamma = 1.0.
		static_assert(CHAR_BIT == 8, "8 bits per byte. No exceptions.");

		png_structp dest_png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		png_infop dest_info = NULL;
		if(dest_png)
		{
			if(!setjmp(png_jmpbuf(dest_png)))
			{
				dest_info = png_create_info_struct(dest_png);
				if(!dest_info)
				{
					_error_errno(*argv, ENOMEM);
					longjmp(png_jmpbuf(dest_png), 1);
				}

				png_init_io(dest_png, dest_fp);

				int color_type = src_image.format & PNG_FORMAT_FLAG_ALPHA ? PNG_COLOR_MASK_ALPHA : 0;
				if(src_image.format & PNG_FORMAT_FLAG_COLOR)
					color_type |= PNG_COLOR_MASK_COLOR;
				png_set_IHDR(
					dest_png,
					dest_info,
					b_cos.dest_width,
					b_cos.dest_height,
					sizeof(*dest_row) * 8,
					color_type,
					PNG_INTERLACE_NONE,
					PNG_COMPRESSION_TYPE_DEFAULT,
					PNG_FILTER_TYPE_DEFAULT);

				// png_set_sRGB_gAMA_and_cHRM(dest_png, dest_info, PNG_sRGB_INTENT_RELATIVE);
				png_set_gAMA(dest_png, dest_info, 1.0);
				// png_set_cHRM(dest_png, dest_info, 0.3127, 0.3290, 0.6400, 0.3300, 0.3000, 0.6000, 0.1500, 0.0600);

				png_write_info(dest_png, dest_info);

#ifndef __BYTE_ORDER__
#	error #define __BYTE_ORDER__ to 1234 (little endian) or 4321 (big endian).
#elif __BYTE_ORDER__ == 1234 || __BYTE_ORDER__ == 3412
				png_set_swap(dest_png);
#endif

				if(b_cos_2d_init(&b_cos))
				{
					for(unsigned y = 0; y != b_cos.dest_height; ++y)
					{
						b_cos_2d_row(&b_cos, dest_row);
						png_write_rows(dest_png, (png_bytepp)&dest_row, 1);
					}

					png_write_info(dest_png, dest_info);
					result = EXIT_SUCCESS;
				}
			}

			png_destroy_write_struct(&dest_png, &dest_info);
		}

		fclose(dest_fp);
	}

	free(src_buffer);
	return result;
}
