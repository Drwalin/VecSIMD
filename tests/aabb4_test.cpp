// Tests for simd::aabb4 (Aabb.hpp).
//
// The suite checks:
//   1. Construction and IsValid for valid, degenerate and inverted boxes.
//   2. Accessors: GetMin/GetMax/GetSizes/GetCenter/GetVolume/GetSurface.
//   3. Expanded.
//   4. HasIntersection (plain, float eps, float4 eps) including touching edges.
//   5. IsIn (plain, float eps, float4 eps) including boundary points.
//   6. Intersection and Sum (box-box and box-point).
//   7. ContainsAll (plain, float eps, float4 eps).
//   8. Operators ==, !=, &&, * and +.
//   9. GetChildRelativeAabb.
//  10. SlowRayTest2/FastRayTest2: concrete cases and a randomized comparison
//      against an independent double-precision slab reference, including the
//      near/far parameter range.
//
// Every check runs on all four SIMD lanes so per-lane independence is tested.
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

#include "../Aabb.hpp"

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

struct F3 {
	float x, y, z;
};
struct FBox {
	F3 mn, mx;
};

static float lf(const float4 &f, int i)
{
	float t[4];
	f.store(t);
	return t[i];
}
static uint32_t lu(const uint4 &m, int i)
{
	uint32_t h[4];
	m.store(h);
	return h[i];
}
static bool laneT(const uint4 &m, int i) { return lu(m, i) != 0u; }

static bool feq(float a, float b) { return std::fabs(a - b) <= 1e-5f; }

static vec3f rep3(const F3 &v)
{
	float x[4], y[4], z[4];
	for (int i = 0; i < 4; ++i) {
		x[i] = v.x;
		y[i] = v.y;
		z[i] = v.z;
	}
	return vec3f::load(x, y, z);
}

static aabb4 makeBox(const FBox b[4])
{
	float mnx[4], mny[4], mnz[4], mxx[4], mxy[4], mxz[4];
	for (int i = 0; i < 4; ++i) {
		mnx[i] = b[i].mn.x;
		mny[i] = b[i].mn.y;
		mnz[i] = b[i].mn.z;
		mxx[i] = b[i].mx.x;
		mxy[i] = b[i].mx.y;
		mxz[i] = b[i].mx.z;
	}
	return {vec3f::load(mnx, mny, mnz), vec3f::load(mxx, mxy, mxz)};
}
static aabb4 make1(const FBox &b)
{
	FBox a[4] = {b, b, b, b};
	return makeBox(a);
}

static bool boxEqLane(const aabb4 &r, const FBox &e, int i)
{
	return feq(lf(r.min.x, i), e.mn.x) && feq(lf(r.min.y, i), e.mn.y) &&
		   feq(lf(r.min.z, i), e.mn.z) && feq(lf(r.max.x, i), e.mx.x) &&
		   feq(lf(r.max.y, i), e.mx.y) && feq(lf(r.max.z, i), e.mx.z);
}

// --- scalar references used by the randomized tests ------------------------

static bool refValid(const FBox &b)
{
	return b.mn.x <= b.mx.x && b.mn.y <= b.mx.y && b.mn.z <= b.mx.z;
}
static bool refHit(const FBox &a, const FBox &b)
{
	return a.mx.x >= b.mn.x && b.mx.x >= a.mn.x && a.mx.y >= b.mn.y &&
		   b.mx.y >= a.mn.y && a.mx.z >= b.mn.z && b.mx.z >= a.mn.z;
}
static bool refContains(const FBox &a, const FBox &b)
{
	return a.mn.x <= b.mn.x && a.mn.y <= b.mn.y && a.mn.z <= b.mn.z &&
		   a.mx.x >= b.mx.x && a.mx.y >= b.mx.y && a.mx.z >= b.mx.z;
}
static bool refIn(const FBox &a, const F3 &p)
{
	return p.x >= a.mn.x && p.x <= a.mx.x && p.y >= a.mn.y && p.y <= a.mx.y &&
		   p.z >= a.mn.z && p.z <= a.mx.z;
}

// --- 1. IsValid -------------------------------------------------------------

static void test_is_valid()
{
	FBox b[4] = {
		{{0, 0, 0}, {1, 1, 1}}, // valid
		{{0, 0, 0}, {0, 0, 0}}, // degenerate, still valid
		{{-1, -1, -1}, {1, 1, 1}},
		{{2, 2, 2}, {1, 1, 1}}, // inverted -> invalid
	};
	uint4 v = makeBox(b).IsValid();
	CHECK(laneT(v, 0), "plain box must be valid");
	CHECK(laneT(v, 1), "degenerate box (min==max) must be valid");
	CHECK(laneT(v, 2), "plain box must be valid");
	CHECK(!laneT(v, 3), "inverted box (min>max) must be invalid");

	// A single inverted axis invalidates the box.
	CHECK(!laneT(make1({{2, 0, 0}, {1, 1, 1}}).IsValid(), 0),
		  "inverted x axis must be invalid");
	CHECK(!laneT(make1({{0, 2, 0}, {1, 1, 1}}).IsValid(), 0),
		  "inverted y axis must be invalid");
	CHECK(!laneT(make1({{0, 0, 2}, {1, 1, 1}}).IsValid(), 0),
		  "inverted z axis must be invalid");
}

// --- 2. accessors -----------------------------------------------------------

static void test_accessors()
{
	FBox b[4] = {
		{{1, 2, 3}, {5, 6, 9}}, // sizes (4,4,6) center (3,4,6) vol 96 surf 128
		{{0, 0, 0}, {1, 1, 1}}, // sizes (1,1,1) center (.5) vol 1 surf 6
		{{-1, -1, -1}, {1, 1, 1}}, // sizes (2,2,2) vol 8 surf 24
		{{2, 3, 4}, {2, 3, 4}}, // sizes 0 vol 0 surf 0
	};
	aabb4 box = makeBox(b);

	vec3f mn = box.GetMin();
	vec3f mx = box.GetMax();
	vec3f sz = box.GetSizes();
	vec3f c = box.GetCenter();
	float4 vol = box.GetVolume();
	float4 surf = box.GetSurface();

	CHECK(feq(lf(mn.x, 0), 1) && feq(lf(mn.y, 0), 2) && feq(lf(mn.z, 0), 3),
		  "GetMin lane0");
	CHECK(feq(lf(mx.x, 0), 5) && feq(lf(mx.y, 0), 6) && feq(lf(mx.z, 0), 9),
		  "GetMax lane0");
	CHECK(feq(lf(sz.x, 0), 4) && feq(lf(sz.y, 0), 4) && feq(lf(sz.z, 0), 6),
		  "GetSizes lane0");
	CHECK(feq(lf(c.x, 0), 3) && feq(lf(c.y, 0), 4) && feq(lf(c.z, 0), 6),
		  "GetCenter lane0");
	CHECK(feq(lf(vol, 0), 96.f), "GetVolume lane0 = 96, got %g", lf(vol, 0));
	CHECK(feq(lf(surf, 0), 128.f), "GetSurface lane0 = 128, got %g",
		  lf(surf, 0));

	CHECK(feq(lf(sz.x, 1), 1) && feq(lf(sz.y, 1), 1) && feq(lf(sz.z, 1), 1),
		  "GetSizes lane1");
	CHECK(feq(lf(c.x, 1), 0.5f), "GetCenter lane1");
	CHECK(feq(lf(vol, 1), 1.f) && feq(lf(surf, 1), 6.f),
		  "GetVolume/GetSurface lane1");

	CHECK(feq(lf(sz.x, 2), 2) && feq(lf(vol, 2), 8.f) && feq(lf(surf, 2), 24.f) &&
			  feq(lf(c.x, 2), 0.f),
		  "lane2 accessors");

	CHECK(feq(lf(sz.x, 3), 0) && feq(lf(sz.y, 3), 0) && feq(lf(sz.z, 3), 0),
		  "degenerate sizes lane3");
	CHECK(feq(lf(vol, 3), 0.f) && feq(lf(surf, 3), 0.f),
		  "degenerate volume/surface lane3");
	CHECK(feq(lf(c.x, 3), 2) && feq(lf(c.y, 3), 3) && feq(lf(c.z, 3), 4),
		  "degenerate center lane3");
}

// --- 3. Expanded ------------------------------------------------------------

static void test_expanded()
{
	aabb4 box = make1({{0, 1, 2}, {3, 4, 5}});

	aabb4 grow = box.Expanded(2.f);
	CHECK(boxEqLane(grow, {{-2, -1, 0}, {5, 6, 7}}, 0),
		  "Expanded(+2) must grow min by -2 and max by +2");

	aabb4 shrink = box.Expanded(-1.f);
	CHECK(boxEqLane(shrink, {{1, 2, 3}, {2, 3, 4}}, 0),
		  "Expanded(-1) must shrink the box by 1 on every side");

	CHECK(boxEqLane(box.Expanded(0.f), {{0, 1, 2}, {3, 4, 5}}, 0),
		  "Expanded(0) must be identity");
}

// --- 4. HasIntersection -----------------------------------------------------

static void test_has_intersection()
{
	aabb4 a = make1({{0, 0, 0}, {2, 2, 2}});
	FBox r[4] = {
		{{1, 1, 1}, {3, 3, 3}}, // overlap
		{{2, 2, 2}, {4, 4, 4}}, // touch at the corner (2,2,2)
		{{3, 3, 3}, {4, 4, 4}}, // disjoint
		{{0.5f, 0.5f, 0.5f}, {1.5f, 1.5f, 1.5f}}, // contained
	};
	uint4 m = a.HasIntersection(makeBox(r));
	CHECK(laneT(m, 0), "overlapping boxes must intersect");
	CHECK(laneT(m, 1), "boxes touching at a face corner must intersect");
	CHECK(!laneT(m, 2), "disjoint boxes must not intersect");
	CHECK(laneT(m, 3), "contained box must intersect");

	// eps widens the test: 1 + eps >= 1.5 only for eps >= 0.5.
	aabb4 ea = make1({{0, 0, 0}, {1, 1, 1}});
	aabb4 eb = make1({{1.5f, 0, 0}, {2.f, 1, 1}});
	CHECK(!laneT(ea.HasIntersection(eb), 0),
		  "gap of 0.5 must not intersect without eps");
	CHECK(laneT(ea.HasIntersection(eb, 0.5f), 0),
		  "gap of 0.5 must intersect with eps = 0.5");
	CHECK(!laneT(ea.HasIntersection(eb, 0.49f), 0),
		  "gap of 0.5 must not intersect with eps = 0.49");

	// Per-lane float4 eps.
	float eps[4] = {0.5f, 0.49f, 1.0f, 0.0f};
	uint4 me = ea.HasIntersection(eb, float4::load(eps));
	CHECK(laneT(me, 0) && !laneT(me, 1) && laneT(me, 2) && !laneT(me, 3),
		  "per-lane eps must be applied independently");
}

// --- 5. IsIn ----------------------------------------------------------------

static void test_is_in()
{
	aabb4 a = make1({{0, 0, 0}, {2, 2, 2}});
	float px[4] = {1.f, 0.f, 2.f, 2.5f};
	float py[4] = {1.f, 0.f, 2.f, 1.f};
	float pz[4] = {1.f, 0.f, 2.f, 1.f};
	uint4 m = a.IsIn(vec3f::load(px, py, pz));
	CHECK(laneT(m, 0), "interior point must be inside");
	CHECK(laneT(m, 1), "min corner must be inside (inclusive)");
	CHECK(laneT(m, 2), "max corner must be inside (inclusive)");
	CHECK(!laneT(m, 3), "point outside must not be inside");

	// Point just outside the +x face widened by eps.
	F3 p = {2.2f, 1.f, 1.f};
	CHECK(!laneT(a.IsIn(rep3(p), 0.1f), 0),
		  "point at 0.2 outside must not be inside with eps = 0.1");
	CHECK(laneT(a.IsIn(rep3(p), 0.3f), 0),
		  "point at 0.2 outside must be inside with eps = 0.3");

	// Per-lane eps.
	float eps[4] = {0.1f, 0.3f, 0.2f, 0.0f};
	uint4 me = a.IsIn(rep3(p), float4::load(eps));
	CHECK(!laneT(me, 0) && laneT(me, 1) && laneT(me, 2) && !laneT(me, 3),
		  "per-lane IsIn eps must be applied independently");
}

// --- 6. Intersection and Sum ------------------------------------------------

static void test_intersection_and_sum()
{
	aabb4 a = make1({{0, 0, 0}, {2, 2, 2}});
	aabb4 b = make1({{1, 1, 1}, {3, 3, 3}});
	CHECK(boxEqLane(a.Intersection(b), {{1, 1, 1}, {2, 2, 2}}, 0),
		  "Intersection of overlapping boxes");

	// Disjoint boxes: min > max, so the result is invalid.
	aabb4 da = make1({{0, 0, 0}, {1, 1, 1}});
	aabb4 db = make1({{2, 2, 2}, {3, 3, 3}});
	aabb4 di = da.Intersection(db);
	CHECK(boxEqLane(di, {{2, 2, 2}, {1, 1, 1}}, 0),
		  "Intersection of disjoint boxes must be an invalid box");
	CHECK(!laneT(di.IsValid(), 0), "disjoint intersection must be invalid");

	aabb4 s = da.Sum(db);
	CHECK(boxEqLane(s, {{0, 0, 0}, {3, 3, 3}}, 0),
		  "Sum of disjoint boxes must enclose both");

	// Sum with a point outside grows, with a point inside is identity.
	CHECK(boxEqLane(da.Sum(rep3({5, 5, 5})), {{0, 0, 0}, {5, 5, 5}}, 0),
		  "Sum with outside point must grow the box");
	CHECK(boxEqLane(da.Sum(rep3({0.5f, 0.5f, 0.5f})), {{0, 0, 0}, {1, 1, 1}}, 0),
		  "Sum with inside point must be identity");
	CHECK(boxEqLane(da.Sum(rep3({-1, -1, -1})), {{-1, -1, -1}, {1, 1, 1}}, 0),
		  "Sum with outside point must grow the box on the min side");
}

// --- 7. ContainsAll ---------------------------------------------------------

static void test_contains_all()
{
	aabb4 a = make1({{0, 0, 0}, {4, 4, 4}});
	FBox r[4] = {
		{{1, 1, 1}, {2, 2, 2}}, // fully inside
		{{3, 3, 3}, {5, 5, 5}}, // sticks out on max side
		{{0, 0, 0}, {4, 4, 4}}, // identical
		{{-1, -1, -1}, {1, 1, 1}}, // sticks out on min side
	};
	uint4 m = a.ContainsAll(makeBox(r));
	CHECK(laneT(m, 0), "inner box must be contained");
	CHECK(!laneT(m, 1), "box with larger max must not be contained");
	CHECK(laneT(m, 2), "identical box must be contained");
	CHECK(!laneT(m, 3), "box with smaller min must not be contained");

	// eps allows a small overshoot: min 0 <= -0.1 + eps only for eps >= 0.1.
	aabb4 inner = make1({{-0.1f, -0.1f, -0.1f}, {2, 2, 2}});
	CHECK(!laneT(a.ContainsAll(inner, 0.05f), 0),
		  "overshoot of 0.1 must not be contained with eps = 0.05");
	CHECK(laneT(a.ContainsAll(inner, 0.2f), 0),
		  "overshoot of 0.1 must be contained with eps = 0.2");

	float eps[4] = {0.05f, 0.2f, 0.1f, 0.0f};
	uint4 me = a.ContainsAll(inner, float4::load(eps));
	CHECK(!laneT(me, 0) && laneT(me, 1) && laneT(me, 2) && !laneT(me, 3),
		  "per-lane ContainsAll eps must be applied independently");
}

// --- 8. operators -----------------------------------------------------------

static void test_operators()
{
	aabb4 a = make1({{0, 0, 0}, {2, 2, 2}});
	aabb4 same = make1({{0, 0, 0}, {2, 2, 2}});
	aabb4 over = make1({{1, 1, 1}, {3, 3, 3}});

	uint4 eq = (a == same);
	CHECK(laneT(eq, 0), "a == a must be true on every lane");
	uint4 ne = (a != same);
	CHECK(!laneT(ne, 0), "a != a must be false");
	uint4 neq = (a != over);
	CHECK(laneT(neq, 0), "a != b must be true when any component differs");
	uint4 eqo = (a == over);
	CHECK(!laneT(eqo, 0), "a == b must be false when any component differs");

	// Operator results must match the named methods on every lane.
	uint4 m1 = (a && over);
	uint4 m2 = a.HasIntersection(over);
	uint4 m3 = (a && rep3({1, 1, 1}));
	uint4 m4 = a.IsIn(rep3({1, 1, 1}));
	for (int i = 0; i < 4; ++i) {
		CHECK(lu(m1, i) == lu(m2, i), "operator&& lane %d must equal "
									  "HasIntersection",
			  i);
		CHECK(lu(m3, i) == lu(m4, i), "operator&&(vec3f) lane %d must equal "
									  "IsIn",
			  i);
	}
	CHECK(boxEqLane(a * over, {{1, 1, 1}, {2, 2, 2}}, 0),
		  "operator* must equal Intersection");
	CHECK(boxEqLane(a + over, {{0, 0, 0}, {3, 3, 3}}, 0),
		  "operator+ must equal Sum");
	CHECK(boxEqLane(a + rep3({5, 5, 5}), {{0, 0, 0}, {5, 5, 5}}, 0),
		  "operator+(vec3f) must equal Sum");
}

// --- 9. GetChildRelativeAabb ------------------------------------------------

static void test_child_relative()
{
	aabb4 parent = make1({{0, 0, 0}, {4, 4, 4}});
	aabb4 child = make1({{0.25f, 0.25f, 0.25f}, {0.5f, 0.5f, 0.5f}});
	aabb4 r = parent.GetChildRelativeAabb(child);
	// new_min = 0 + 4 * 0.25 = 1 ; new_max = 1 + (4-1)*0.5 = 2.5
	CHECK(boxEqLane(r, {{1, 1, 1}, {2.5f, 2.5f, 2.5f}}, 0),
		  "GetChildRelativeAabb must map a unit-space child into the parent");

	// A child covering the whole unit range maps onto the parent.
	CHECK(boxEqLane(parent.GetChildRelativeAabb(make1({{0, 0, 0}, {1, 1, 1}})),
					{{0, 0, 0}, {4, 4, 4}}, 0),
		  "GetChildRelativeAabb of the full unit box must equal the parent");
}

// --- 10. ray tests ----------------------------------------------------------

static void test_slow_ray_basic()
{
	aabb4 box = make1({{0, 0, 0}, {2, 2, 2}});
	float4 near, far;

	// Crossing the box along +x: t = (0 - (-1)) / 4 = 0.25 .. 0.75.
	uint4 m = box.SlowRayTest2(rep3({-1, 1, 1}), rep3({3, 1, 1}), near, far);
	CHECK(laneT(m, 0), "segment crossing the box must hit");
	CHECK(feq(lf(near, 0), 0.25f) && feq(lf(far, 0), 0.75f),
		  "near/far must be 0.25/0.75, got %g/%g", lf(near, 0), lf(far, 0));

	// Segment fully inside the box: clipped to the whole segment.
	m = box.SlowRayTest2(rep3({1, 1, 1}), rep3({1.5f, 1.5f, 1.5f}), near, far);
	CHECK(laneT(m, 0), "segment inside the box must hit");
	CHECK(feq(lf(near, 0), 0.f) && feq(lf(far, 0), 1.f),
		  "inside near/far must be 0/1, got %g/%g", lf(near, 0), lf(far, 0));

	// Segment entirely beyond +x: misses.
	m = box.SlowRayTest2(rep3({3, 1, 1}), rep3({5, 1, 1}), near, far);
	CHECK(!laneT(m, 0), "segment past the box must miss");

	// Segment before the box: misses.
	m = box.SlowRayTest2(rep3({-5, 1, 1}), rep3({-3, 1, 1}), near, far);
	CHECK(!laneT(m, 0), "segment before the box must miss");

	// Parallel to y but offset outside the y slab: misses.
	m = box.SlowRayTest2(rep3({-1, 5, 1}), rep3({3, 5, 1}), near, far);
	CHECK(!laneT(m, 0), "ray parallel to an axis outside the slab must miss");

	// Zero-length segment starting inside hits, starting outside misses.
	m = box.SlowRayTest2(rep3({1, 1, 1}), rep3({1, 1, 1}), near, far);
	CHECK(laneT(m, 0), "zero-length segment inside the box must hit");
	m = box.SlowRayTest2(rep3({3, 1, 1}), rep3({3, 1, 1}), near, far);
	CHECK(!laneT(m, 0), "zero-length segment outside the box must miss");
}

// Boundary-start and boundary-end ray cases on the box [0,2]^3.
static void test_ray_boundary_cases()
{
	aabb4 box = make1({{0, 0, 0}, {2, 2, 2}});

	// wantT = compare near/far against (wn, wf) as well as the hit flag.
	auto ray = [&](F3 s, F3 e, bool wantHit, bool wantT, float wn, float wf,
				   const char *what) {
		float4 n, f;
		uint4 m = box.SlowRayTest2(rep3(s), rep3(e), n, f);
		CHECK(laneT(m, 0) == wantHit, "hit mismatch (%s): got %d want %d",
			  what, (int)laneT(m, 0), (int)wantHit);
		if (wantHit && wantT)
			CHECK(feq(lf(n, 0), wn) && feq(lf(f, 0), wf),
				  "t-range mismatch (%s): got [%g,%g] want [%g,%g]", what,
				  lf(n, 0), lf(f, 0), wn, wf);
	};

	// 1. Start strictly inside and exit through the +x face: near = 0,
	//    far = (2 - 1) / 4 = 0.25, with 0 <= far <= 1.
	ray({1, 1, 1}, {5, 1, 1}, true, true, 0.f, 0.25f, "start inside, exits");

	// 2. Start strictly outside and point further away: no collision.
	ray({3, 1, 1}, {4, 1, 1}, false, false, 0, 0, "outside, away on +x");
	ray({1, 1, 5}, {1, 1, 6}, false, false, 0, 0, "outside, away on +z");
	ray({1, 1, -3}, {1, 1, -5}, false, false, 0, 0, "outside, away on -z");

	// 3a. Start on the min vertex and head inside: near = 0, far = 1.
	ray({0, 0, 0}, {1, 1, 1}, true, true, 0.f, 1.f, "vertex start, inward");
	// 3b. Start on the min vertex and head directly away: the ray only
	//     touches the corner at t = 0, so it is a degenerate hit (near=far=0).
	ray({0, 0, 0}, {-1, -1, -1}, true, true, 0.f, 0.f, "vertex start, outward");
	// 3c. Start on the min-x face and head inside / away. The away ray
	//     grazes the face plane at t = 0 only -> degenerate hit.
	ray({0, 1, 1}, {2, 1, 1}, true, true, 0.f, 1.f, "face start, inward");
	ray({0, 1, 1}, {-2, 1, 1}, true, true, 0.f, 0.f, "face start, outward");
	// 3d. Start on an edge (x=0, y=0) and head inside.
	ray({0, 0, 1}, {2, 2, 1}, true, true, 0.f, 1.f, "edge start, inward");
	// 3e. Start outside and travel alongside the y=0 face (parallel, on the
	//     boundary): the face plane is inclusive, so this is a hit.
	ray({-1, 0, 1}, {3, 0, 1}, true, true, 0.25f, 0.75f, "alongside y=0 face");

	// 4a. End exactly on the min-x face: the hit is at t = 1.
	ray({-1, 1, 1}, {0, 1, 1}, true, true, 1.f, 1.f, "end on face");
	// 4b. End exactly on an edge (x=0, y=0).
	ray({-1, -1, 1}, {0, 0, 1}, true, true, 1.f, 1.f, "end on edge");
	// 4c. End exactly on the min vertex.
	ray({-1, -1, -1}, {0, 0, 0}, true, true, 1.f, 1.f, "end on vertex");

	// --- Grazing: ray parallel to a face, lying exactly on that face
	//     (y = 0). The intersection is the segment x in [0, 2] at y = 0,
	//     z = 1, so the near/far range depends on where the ends fall.

	// start outside, end outside (segment strictly inside the ends).
	ray({-1, 0, 1}, {3, 0, 1}, true, true, 0.25f, 0.75f, "on-face out->out");
	// start outside, end on the face (inside the face region).
	ray({-2, 0, 1}, {1, 0, 1}, true, true, 2.f / 3.f, 1.f, "on-face out->face");
	// start on the face, end outside.
	ray({1, 0, 1}, {4, 0, 1}, true, true, 0.f, 1.f / 3.f, "on-face face->out");
	// start and end both on the face.
	ray({0.5f, 0, 1}, {1.5f, 0, 1}, true, true, 0.f, 1.f, "on-face face->face");
	// along the face edge (y = 0 and z = 0), starting on a vertex.
	ray({0, 0, 0}, {3, 0, 0}, true, true, 0.f, 2.f / 3.f, "on-edge vertex->out");
	// along the face edge, ending on a vertex.
	ray({-1, 0, 0}, {0, 0, 0}, true, true, 1.f, 1.f, "on-edge out->vertex");

	// --- Grazing: single intersection point at the middle of an edge.
	//     Edge y = 0, z = 0, x = 1 (interior). The ray leaves the box
	//     immediately, so the only contact is at t = 0.
	ray({1, 0, 0}, {1, -1, -1}, true, true, 0.f, 0.f, "edge-middle graze t=0");
	// Same single point reached at the far end of the segment (t = 1).
	ray({1, -1, -1}, {1, 0, 0}, true, true, 1.f, 1.f, "edge-middle graze t=1");
}

// Independent double-precision slab reference for the segment ray test.
// The result is clipped to t in [0, 1], matching SlowRayTest2.
struct RefSeg {
	bool hit;
	double lo, hi;   // clipped t range
	bool ambiguous;  // true near a decision boundary: skip exact comparison
};

static RefSeg refSeg(const FBox &b, const F3 &s, const F3 &e)
{
	double d[3] = {e.x - s.x, e.y - s.y, e.z - s.z};
	double o[3] = {s.x, s.y, s.z};
	double mn[3] = {b.mn.x, b.mn.y, b.mn.z};
	double mx[3] = {b.mx.x, b.mx.y, b.mx.z};
	double lo = -1e30, hi = 1e30;
	for (int a = 0; a < 3; ++a) {
		if (d[a] == 0.0) {
			if (o[a] < mn[a] || o[a] > mx[a])
				return {false, 0, 0, false};
		} else {
			double inv = 1.0 / d[a];
			double ta = (mn[a] - o[a]) * inv;
			double tb = (mx[a] - o[a]) * inv;
			if (ta > tb) {
				double t = ta;
				ta = tb;
				tb = t;
			}
			if (ta > lo)
				lo = ta;
			if (tb < hi)
				hi = tb;
		}
	}
	bool hit = lo <= hi && hi >= 0.0 && lo <= 1.0;
	double lo2 = lo < 0.0 ? 0.0 : lo;
	double hi2 = hi > 1.0 ? 1.0 : hi;
	double m1 = fabs(hi - lo), m2 = fabs(hi), m3 = fabs(1.0 - lo);
	double margin = m1 < m2 ? m1 : m2;
	if (m3 < margin)
		margin = m3;
	return {hit, lo2, hi2, margin < 1e-4};
}

static void test_ray_fuzz()
{
	std::mt19937_64 rng(0xABCDEF12u);
	std::uniform_real_distribution<float> C(-2.f, 2.f);
	std::uniform_real_distribution<float> H(0.f, 2.f);
	std::uniform_real_distribution<float> Z(0.f, 1.f);

	const long long iterations = 100000;
	long long compared = 0, mismatches = 0;

	for (long long it = 0; it < iterations; ++it) {
		FBox bx[4];
		F3 ro[4], en[4];
		for (int i = 0; i < 4; ++i) {
			F3 c = {C(rng), C(rng), C(rng)};
			F3 h = {H(rng), H(rng), H(rng)};
			bx[i] = {{c.x - h.x, c.y - h.y, c.z - h.z},
					 {c.x + h.x, c.y + h.y, c.z + h.z}};
			ro[i] = {C(rng), C(rng), C(rng)};
			if (Z(rng) < (1.f / 16.f))
				en[i] = ro[i]; // occasionally a zero-length segment
			else
				en[i] = {C(rng), C(rng), C(rng)};
		}

		float rx[4], ry[4], rz[4], ex[4], ey[4], ez[4];
		for (int i = 0; i < 4; ++i) {
			rx[i] = ro[i].x;
			ry[i] = ro[i].y;
			rz[i] = ro[i].z;
			ex[i] = en[i].x;
			ey[i] = en[i].y;
			ez[i] = en[i].z;
		}
		aabb4 box = makeBox(bx);
		vec3f roV = vec3f::load(rx, ry, rz);
		vec3f enV = vec3f::load(ex, ey, ez);
		vec3f dir = enV - roV;
		float4 zero = float4::load(0.0f);
		float4 one = float4::load(1.0f);
		vec3f inv = {one / dir.x, one / dir.y, one / dir.z};
		uint4 par[3] = {dir.x == zero, dir.y == zero, dir.z == zero};
		int4 sign[3] = {int4::reinterpret(inv.x < zero),
						int4::reinterpret(inv.y < zero),
						int4::reinterpret(inv.z < zero)};

		float4 nearA, farA, nearB, farB, nearC, farC;
		uint4 mA = box.SlowRayTest2(roV, enV, nearA, farA);
		uint4 mB = box.FastRayTest2(roV, inv, par, nearB, farB);
		uint4 mC = box.FastRayTest2(roV, inv, sign, par, nearC, farC);

		for (int i = 0; i < 4; ++i) {
			// The three entry points implement the same test: identical output.
			CHECK(lu(mA, i) == lu(mB, i) && lu(mB, i) == lu(mC, i),
				  "ray entry points must agree (it=%lld lane=%d)", it, i);
			CHECK(feq(lf(nearA, i), lf(nearB, i)) &&
					  feq(lf(farA, i), lf(farB, i)),
				  "SlowRayTest2 and FastRayTest2 near/far must agree");

			RefSeg rf = refSeg(bx[i], ro[i], en[i]);
			if (rf.ambiguous)
				continue;
			++compared;
			bool got = laneT(mA, i);
			if (got != rf.hit) {
				++mismatches;
				if (mismatches <= 10)
					std::printf("  ray mismatch it=%lld lane=%d got=%d want=%d\n",
								it, i, got, rf.hit);
				continue;
			}
			if (rf.hit) {
				if (!feq(lf(nearA, i), (float)rf.lo) ||
					!feq(lf(farA, i), (float)rf.hi)) {
					++mismatches;
					if (mismatches <= 10)
						std::printf("  t-range mismatch it=%lld lane=%d "
									"got=[%g,%g] want=[%g,%g]\n",
									it, i, lf(nearA, i), lf(farA, i),
									(float)rf.lo, (float)rf.hi);
				}
			}
		}
	}

	CHECK(mismatches == 0,
		  "ray fuzz vs double reference: %lld mismatches over %lld cases",
		  mismatches, compared);
	std::printf("  ray fuzz: %lld cases compared, %lld mismatches\n", compared,
				mismatches);
}

// Randomized comparison of the boolean AABB queries against scalar references.
static void test_boolean_fuzz()
{
	std::mt19937_64 rng(0x11223344u);
	std::uniform_real_distribution<float> C(-2.f, 2.f);
	std::uniform_real_distribution<float> H(0.f, 2.f);

	const long long iterations = 200000;
	long long mismatches = 0;

	for (long long it = 0; it < iterations; ++it) {
		FBox a[4], b[4];
		F3 p[4];
		for (int i = 0; i < 4; ++i) {
			F3 c1 = {C(rng), C(rng), C(rng)};
			F3 h1 = {H(rng), H(rng), H(rng)};
			F3 c2 = {C(rng), C(rng), C(rng)};
			F3 h2 = {H(rng), H(rng), H(rng)};
			a[i] = {{c1.x - h1.x, c1.y - h1.y, c1.z - h1.z},
					{c1.x + h1.x, c1.y + h1.y, c1.z + h1.z}};
			b[i] = {{c2.x - h2.x, c2.y - h2.y, c2.z - h2.z},
					{c2.x + h2.x, c2.y + h2.y, c2.z + h2.z}};
			p[i] = {C(rng), C(rng), C(rng)};
		}

		aabb4 A = makeBox(a);
		aabb4 B = makeBox(b);
		float px[4], py[4], pz[4];
		for (int i = 0; i < 4; ++i) {
			px[i] = p[i].x;
			py[i] = p[i].y;
			pz[i] = p[i].z;
		}
		vec3f P = vec3f::load(px, py, pz);

		uint4 v = A.IsValid();
		uint4 h = A.HasIntersection(B);
		uint4 c = A.ContainsAll(B);
		uint4 in = A.IsIn(P);

		for (int i = 0; i < 4; ++i) {
			bool ok = true;
			if (laneT(v, i) != refValid(a[i]))
				ok = false;
			if (laneT(h, i) != refHit(a[i], b[i]))
				ok = false;
			if (laneT(c, i) != refContains(a[i], b[i]))
				ok = false;
			if (laneT(in, i) != refIn(a[i], p[i]))
				ok = false;
			if (!ok) {
				++mismatches;
				if (mismatches <= 10)
					std::printf("  bool mismatch it=%lld lane=%d\n", it, i);
			}
		}
	}

	CHECK(mismatches == 0,
		  "boolean fuzz: %lld mismatches over %lld lane cases", mismatches,
		  iterations * 4);
	std::printf("  boolean fuzz: %lld lane cases, %lld mismatches\n",
				iterations * 4, mismatches);
}

int main()
{
	std::printf("aabb4 test suite\n");

	test_is_valid();
	test_accessors();
	test_expanded();
	test_has_intersection();
	test_is_in();
	test_intersection_and_sum();
	test_contains_all();
	test_operators();
	test_child_relative();
	test_slow_ray_basic();
	test_ray_boundary_cases();
	test_ray_fuzz();
	test_boolean_fuzz();

	std::printf("%d checks, %d failures\n", g_checks, g_failures);
	return g_failures == 0 ? 0 : 1;
}
