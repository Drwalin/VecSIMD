// Tests for the watertight ray/triangle intersection implemented in
// RayTriangleTest.cpp.
//
// The suite checks four properties:
//   1. The SIMD test with its double-precision fallback returns exactly the
//      same hit/miss decision as the double-precision reference when both are
//      given identical single-precision inputs.
//   2. A plane tiled by adjacent triangles has no cracks: every ray that
//      crosses the tiled region hits at least one triangle, including rays
//      that pass exactly through shared edges and grid lines.
//   3. Two triangles sharing an edge report a consistent result for rays that
//      pass through a point on that shared edge.
//   4. Special cases behave as documented: a ray parallel to the triangle
//      plane, a zero-length segment, a segment that stops before the triangle,
//      reversed direction, and the exact intersection point and normal.
//
// Build and run with CMake from the repository root:
//     cmake -S . -B build
//     cmake --build build
//     ctest --test-dir build
//
// The process exits with status 0 only when every check passes.

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

#include "../RayTriangleTest.cpp"

using namespace simd;

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(cond, fmt, ...)                                                  \
	do {                                                                       \
		++g_checks;                                                            \
		if (!(cond)) {                                                         \
			++g_failures;                                                      \
			std::printf("FAIL %s:%d: " fmt "\n", __FILE__, __LINE__,           \
						##__VA_ARGS__);                                        \
		}                                                                      \
	} while (0)

// A trivial double-precision vector used to build test geometry.
struct D3 {
	double x, y, z;
};

static D3 sub(const D3 &a, const D3 &b)
{
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// Round to the nearest single-precision value but keep it in a double so the
// reference test can be fed the exact same inputs as the SIMD test.
static double asFloat(double v) { return (double)(float)v; }
static D3 asFloat(const D3 &v)
{
	return {asFloat(v.x), asFloat(v.y), asFloat(v.z)};
}

static vec3f mk4(const D3 *v)
{
	float xs[4], ys[4], zs[4];
	for (int i = 0; i < 4; ++i) {
		xs[i] = (float)v[i].x;
		ys[i] = (float)v[i].y;
		zs[i] = (float)v[i].z;
	}
	return vec3f::load(xs, ys, zs);
}

// Single-lane call: replicate the case into all four lanes. Lane 0 of the
// returned intersection and normal is written back through I and N.
static bool castFull(const D3 &start, const D3 &end, const D3 &v0, const D3 &v1,
					 const D3 &v2, vec3f &I, vec3f &N)
{
	const D3 d1 = sub(v1, v2);
	const D3 d2 = sub(v2, v0);
	const D3 S[4] = {start, start, start, start};
	const D3 E[4] = {end, end, end, end};
	const D3 V0[4] = {v0, v0, v0, v0};
	const D3 D1[4] = {d1, d1, d1, d1};
	const D3 D2[4] = {d2, d2, d2, d2};
	uint4 mask = RayTriangleTest(mk4(S), mk4(E), mk4(V0), mk4(D1), mk4(D2), I, N);
	uint32_t h[4];
	mask.store(h);
	return h[0] != 0u;
}

static bool castOne(const D3 &start, const D3 &end, const D3 &v0, const D3 &v1,
					const D3 &v2)
{
	vec3f I, N;
	return castFull(start, end, v0, v1, v2, I, N);
}

// Double-precision reference on the given inputs (no rounding).
static bool castOneDouble(const D3 &start, const D3 &end, const D3 &v0,
						  const D3 &v1, const D3 &v2)
{
	const D3 d1 = sub(v1, v2);
	const D3 d2 = sub(v2, v0);
	vec3d I, N;
	return RayTriangleTestDouble({start.x, start.y, start.z},
	                             {end.x, end.y, end.z}, {v0.x, v0.y, v0.z},
	                             {d1.x, d1.y, d1.z}, {d2.x, d2.y, d2.z}, I, N);
}

// 1. Random fuzzing: the SIMD test must match the double reference on
// identical single-precision inputs. Four independent cases are packed into
// every call so the batched code path and its per-lane fallback are exercised.
static void test_fuzz_simd_matches_double()
{
	std::mt19937_64 rng(0xC0FFEEu);
	std::uniform_real_distribution<double> U(-3.0, 3.0);
	const long long iterations = 10000000; // 10,000,000 iterations x 4 lanes
	long long mismatches = 0;

	for (long long it = 0; it < iterations; ++it) {
		D3 s[4], e[4], v0[4], v1[4], v2[4], d1[4], d2[4];
		for (int i = 0; i < 4; ++i) {
			s[i] = asFloat({U(rng), U(rng), U(rng)});
			e[i] = asFloat({s[i].x + U(rng) * 0.5, s[i].y + U(rng) * 0.5,
							s[i].z + U(rng) * 0.5});
			v0[i] = asFloat({U(rng), U(rng), U(rng)});
			v1[i] = asFloat({U(rng), U(rng), U(rng)});
			v2[i] = asFloat({U(rng), U(rng), U(rng)});
			d1[i] = sub(v1[i], v2[i]);
			d2[i] = sub(v2[i], v0[i]);
		}

		vec3f I, N;
		uint4 mask = RayTriangleTest(mk4(s), mk4(e), mk4(v0), mk4(d1), mk4(d2),
									 I, N);
		uint32_t h[4];
		mask.store(h);

		for (int i = 0; i < 4; ++i) {
			const bool got = h[i] != 0u;
			const bool want = castOneDouble(s[i], e[i], v0[i], v1[i], v2[i]);
			if (got != want) {
				++mismatches;
				if (mismatches <= 10) {
					std::printf("  mismatch it=%lld lane=%d got=%d want=%d\n",
								it, i, got, want);
				}
			}
		}
	}

	CHECK(mismatches == 0, "SIMD vs double over %lld cases: %lld mismatches",
		  iterations * 4, mismatches);
	std::printf("  fuzz: %lld cases, %lld mismatches\n", iterations * 4,
				mismatches);
}

struct Tri {
	D3 a, b, c;
};

// 2. A plane tiled by triangles must have no cracks. Rays are cast straight
// down through exactly representable sample coordinates, including the shared
// edges and grid lines, so the boundary handling is stressed directly.
static void test_watertight_grid_no_cracks()
{
	const int grid = 4;
	std::vector<Tri> tris;
	for (int i = 0; i < grid; ++i) {
		for (int j = 0; j < grid; ++j) {
			D3 p00 = {1.0 * i, 1.0 * j, 0.0};
			D3 p10 = {1.0 * (i + 1), 1.0 * j, 0.0};
			D3 p11 = {1.0 * (i + 1), 1.0 * (j + 1), 0.0};
			D3 p01 = {1.0 * i, 1.0 * (j + 1), 0.0};
			tris.push_back({p00, p10, p11});
			tris.push_back({p00, p11, p01});
		}
	}

	long long tested = 0, cracks = 0;
	const int steps = grid * 8;
	for (int xi = 0; xi <= steps; ++xi) {
		for (int yi = 0; yi <= steps; ++yi) {
			const double x = xi / 8.0;
			const double y = yi / 8.0;
			const D3 s = {x, y, 5.0};
			const D3 e = {x, y, -5.0};
			int hits = 0;
			for (const Tri &t : tris) {
				if (castOne(s, e, t.a, t.b, t.c))
					++hits;
			}
			++tested;
			if (hits == 0) {
				++cracks;
				if (cracks <= 5)
					std::printf("  crack at (%.4f, %.4f)\n", x, y);
			}
		}
	}

	CHECK(cracks == 0, "watertight grid: %lld cracks over %lld rays", cracks,
		  tested);
	std::printf("  grid: %lld rays, %lld cracks\n", tested, cracks);
}

// 3. Two triangles sharing an edge must agree for rays through that edge: no
// crack (both miss) is allowed. The shared edge is the x axis so points on it
// are exactly representable.
static void test_shared_edge_consistency()
{
	const D3 P0 = {0, 0, 0};
	const D3 P1 = {4, 0, 0};
	const D3 Q = {2, -3, 0};
	const D3 R = {2, 3, 0};

	long long both = 0, exactlyOne = 0, none = 0;
	const int samples = 65;
	for (int k = 0; k < samples; ++k) {
		const double t = k / (double)(samples - 1);
		const D3 Pt = {4.0 * t, 0, 0};
		const D3 s = {Pt.x, 0, -5};
		const D3 e = {Pt.x, 0, 5};
		const bool h1 = castOne(s, e, P0, P1, Q);
		const bool h2 = castOne(s, e, P0, P1, R);
		if (h1 && h2)
			++both;
		else if (h1 || h2)
			++exactlyOne;
		else
			++none;
	}

	CHECK(none == 0, "shared edge: %lld rays missed both triangles", none);
	std::printf("  shared edge: both=%lld one=%lld none=%lld\n", both, exactlyOne,
				none);
}

// 4. Documented special cases.
static void test_special_cases()
{
	const D3 A = {0, 0, 0};
	const D3 B = {2, 0, 0};
	const D3 C = {0, 2, 0};

	// Ray lying in the triangle plane: no valid perpendicular hit.
	CHECK(!castOne({-1, 0.5, 0}, {2, 0.5, 0}, A, B, C),
		  "ray parallel to the triangle plane must miss");

	// Zero-length segment.
	CHECK(!castOne({0.5, 0.5, 0}, {0.5, 0.5, 0}, A, B, C),
		  "zero-length segment must miss");

	// Segment that stops before reaching the triangle.
	CHECK(!castOne({0.5, 0.5, 5}, {0.5, 0.5, 1}, A, B, C),
		  "segment ending before the triangle must miss");

	// Straight hit from above, dominant axis +z.
	{
		vec3f I, N;
		const bool hit = castFull({0.5, 0.5, 5}, {0.5, 0.5, -5}, A, B, C, I, N);
		float ix[4], iy[4], iz[4], nx[4], ny[4], nz[4];
		I.store(ix, iy, iz);
		N.store(nx, ny, nz);
		CHECK(hit, "straight hit from above must be reported");
		CHECK(std::fabs(ix[0] - 0.5f) < 1e-6f && std::fabs(iy[0] - 0.5f) < 1e-6f &&
				  std::fabs(iz[0]) < 1e-6f,
			  "intersection must be (0.5, 0.5, 0), got (%g, %g, %g)", ix[0],
			  iy[0], iz[0]);
		CHECK(std::fabs(std::fabs(nz[0]) - 4.0f) < 1e-6f &&
				  std::fabs(nx[0]) < 1e-6f && std::fabs(ny[0]) < 1e-6f,
			  "normal must be +/-(0, 0, 4), got (%g, %g, %g)", nx[0], ny[0],
			  nz[0]);
		std::printf("  straight hit: I=(%g,%g,%g) N=(%g,%g,%g)\n", ix[0], iy[0],
					iz[0], nx[0], ny[0], nz[0]);
	}

	// Reversed direction (negative dominant axis) must still hit.
	CHECK(castOne({0.5, 0.5, -5}, {0.5, 0.5, 5}, A, B, C),
		  "reversed-direction hit must be reported");
}

int main()
{
	std::printf("ray/triangle watertight test suite\n");

	test_fuzz_simd_matches_double();
	test_watertight_grid_no_cracks();
	test_shared_edge_consistency();
	test_special_cases();

	std::printf("%d checks, %d failures\n", g_checks, g_failures);
	return g_failures == 0 ? 0 : 1;
}
// End of test suite.
