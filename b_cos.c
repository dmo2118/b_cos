#include "b_cos.h"

#include <assert.h>
#include <math.h>

/* See ifac(r) in Notes.wxm. */
const float b_cos_ifac0[] = {0.1591549430919f};
const float b_cos_ifac1[] = {0.30965176733713f, 0.075248412122618f};
const float b_cos_ifac2[] = {0.34026032334082f, 0.16239654295838f, 0.047895901889277f};
const float b_cos_ifac3[] = {0.32925212442321f, 0.1827211439557f,  0.11163843263573f, 0.034892547831776f};
const float b_cos_ifac4[] = {0.32486715557368f, 0.17285820298449f, 0.12830005981992f, 0.085446875718292f, 0.027378303380124f};

const float b_cos_ifac5[] =
{
	0.32267868480766f, 0.16815062417398f, 0.12028996640424f, 0.099940515915626f, 0.069339386290332f, 0.022504543394944f
};

const float b_cos_ifac6[] =
{
	0.32142934375849f, 0.16552451341039f, 0.11600127965337f, 0.093468511487219f, 0.082269241336546f, 0.058393672154302f,
	0.019094379912289f
};

const float b_cos_ifac7[] =
{
	0.32064895631628f,  0.16390622237162f, 0.11342010778027f, 0.089708848640425f, 0.07698003589195f,  0.070097481615884f,
	0.050456619693944f, 0.01657726561029f
};

const float b_cos_ifac8[] =
{
	0.32012889127952f,  0.16283705595935f, 0.11173999911071f, 0.087314694649599f, 0.073712104702606f, 0.065712624701937f,
	0.061158204534577f, 0.044431357595242f, 0.014644177635092f
};

const float b_cos_ifac9[] =
{
	0.31976493794602f,  0.16209320907201f, 0.1105827654849f,  0.085689497888924f, 0.071537132589875f, 0.062868712871305f,
	0.057472277708297f, 0.054292928108874f, 0.03969841048389f, 0.013113548077276f
};

#if 1
static float _ib_cos(unsigned r, const float *icos_fac, float x)
{
	/*
	if(x < -(r + 0.5))
		return -0.5;
	if(x > r + 0.5)
		return 0.5;
	*/

	assert(x >= -(r + 0.5));
	assert(x <=   r + 0.5);

	float result = x / (r * 2 + 1);
	for(unsigned i = 0; i <= r; ++i)
		result += icos_fac[i] * sinf(x * (2 * i + 2) * (float)M_PI / (r * 2 + 1));
	return result;
}
#endif

// Nearest-neighbor
#if 0
static float _ib_cos(unsigned r, const float *icos_fac, float x)
{
	if(x < -0.5)
		return -0.5;
	if(x > 0.5)
		return 0.5;
	return x;
}
#endif

// Linear
#if 0
static float _ib_cos(unsigned r, const float *icos_fac, float x)
{
	if(x < 0)
	{
		if(x > -1)
			return (x + 2) * x / 2;
		return -0.5;
	}

	if(x < 1)
		return -(x - 2) * x / 2;
	return 0.5;
}
#endif

// Cubic
#if 0
static float _ib_cos(unsigned r, const float *icos_fac, float x)
{
	float x2 = x * x, x3 = (5.f/6.f) * x2*x, x4 = x2*x2 / 8;

	if(x < 0)
	{
		if(x > -1)
			return -3 * x4 - x3 + x;
		if(x > -2)
			return x4 + x3 + 2 * x2 + 2 * x + (1.f/6.f);
		return -0.5;
	}

	if(x < 1)
		return 3 * x4 - x3 + x;
	if(x < 2)
		return -x4 + x3 - 2 * x2 + 2 * x - (1.f/6.f);
	return 0.5;
}
#endif

void b_cos_edge(float *fac, unsigned r, const float *icos_fac, unsigned src, unsigned dst, unsigned dest_x)
{
	// TODO slow
	float src_x0d = (float)(dest_x * src) / dst;
	unsigned src_x0 = dest_x * src / dst;

	for(unsigned i = 0; i != B_COS_DIAM(r); ++i)
		fac[i] = _ib_cos(r, icos_fac, (src_x0d + (int)(r - i - src_x0)) - 0.5f);
}

// Binary GCD/Stein's algorithm, adapted from Wikipedia. (CC BY-SA)
// https://en.wikipedia.org/wiki/Binary_GCD_algorithm#Iterative_version_in_C
static inline unsigned _gcd(unsigned u, unsigned v)
{
	if(u == 0)
		return v;
	if(v == 0)
		return u;

	// Let shift := lg K, where K is the greatest power of 2
	// dividing both u and v.
	int shift = 0;
	while(!((u | v) & 1))
	{
		u >>= 1;
		v >>= 1;
		++shift;
	}

	while(!(u & 1))
		u >>= 1;

	// From here on, u is always odd.
	do
	{
		// remove all factors of 2 in v -- they are not common
		//   note: v is not zero, so while will terminate
		while(!(v & 1))  // Loop X
			v >>= 1;

		// Now u and v are both odd. Swap if necessary so u <= v,
		// then set v = v - u (which is even). For bignums, the
		// swapping is just pointer movement, and the subtraction
		// can be done in-place.
		if(u > v)
		{
			unsigned t = v;
			v = u;
			u = t;
		}  // Swap u and v.

		v = v - u;                       // Here v >= u.
	} while(v);

	// restore common factors of 2
	return u << shift;
}

unsigned b_cos_repeat(unsigned src, unsigned dst)
{
	return dst / _gcd(src % dst, dst);
}
