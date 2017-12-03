b_cos
=====
An image resampling algorithm.

## Goals
1. The math ought to be the same for both zoomed-in views and zoomed-out views,
   and at 1:1 the resampled image should exactly match the original image.
2. The resampled image, when zoomed in, should not show discontinuities
   between pixels. (i.e. it has to look pretty.)
3. The resampled image, when zoomed out, should average pixels together in such
   a way as to minimize aliasing artifacts. (Again, "pretty.")

## An assumption
The color at each pixel does not represent a point sample from a continuous 2D
surface, but instead represents the average of a small square area* ±0.5
pixel-sized units from the pixel origin in the x and y direction, though still
from a 2D surface. This is held to be true for both the input and output image.

*Counterpoint: Smith, Alvy Ray (July 17, 1995). [A Pixel is Not a Little
Square, a Pixel is Not a Little Square, a Pixel is Not a Little Square! (And a
Voxel is Not a Little Cube)](http://alvyray.com/Memos/CG/Microsoft/6_pixel.pdf).

## Design
Ultimately, this calls for a reconstruction kernel — let's call it B(x) — such
that ∫B(x) dx = 1.0 in the range [-0.5, 0.5], 0.0 for the ranges [-1.5, -0.5],
(0.5, 1.5), extending outward from the center ([1.5, 2.5] is 0, [2.5, 3.5] is 0,
and so on) to some predefined maximum radius.

The kernel implemented here is built from a sum of cosine waves of varying
frequencies and magnitudes. The result looks a bit like a
[Lanczos kernel](https://en.wikipedia.org/wiki/Lanczos_resampling), though the
b_cos kernel stops at half-integers rather than integers. The b_cos kernel has a
greater magnitude as well.

## Usage

```
./b_cos src.png dest.png dest_width dest_height radius
```

`dest_width` and `dest_height` must both be integers greater than 0. `radius` is
the number of neighboring pixels included in the kernel (try `2`) — larger is
smoother, but slower to run, and with more ringing artifacts.

<!-- TODO: Needs pictures, equations, the whole 9 yards. -->
