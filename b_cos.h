#ifndef B_COS_H
#define B_COS_H

// This is not an application for which double-precision floating point is typically necessary.

extern const float b_cos_ifac0[];
extern const float b_cos_ifac1[];
extern const float b_cos_ifac2[];
extern const float b_cos_ifac3[];
extern const float b_cos_ifac4[];
extern const float b_cos_ifac5[];
extern const float b_cos_ifac6[];
extern const float b_cos_ifac7[];
extern const float b_cos_ifac8[];
extern const float b_cos_ifac9[];

void b_cos_edge(float *edge_fac, unsigned r, const float *icos_fac, unsigned src_size, unsigned dst_size, unsigned dest_x);
unsigned b_cos_repeat(unsigned src, unsigned dst);

#define B_COS_EDGE(r) ((r) * 2 + 1)

#endif
