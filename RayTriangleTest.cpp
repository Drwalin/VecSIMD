// This file is part of VecSIMD
// Copyright (c) 2026 Marek Zalewski aka Drwalin
// Distributed under the MIT License.

#include <cmath>
#include <cstdint>

#include "Vec3.hpp"

namespace simd
{

// Scalar double-precision reference vector for the fallback test.
struct vec3d {
	double x, y, z;
};

// Double-precision watertight ray/triangle intersection (Woop, Benthin,
// Wald, JCGT 2(1), 2013).
//
// Triangle: v0 with edge deltas d1 = v1 - v2 and d2 = v2 - v0, so the edges
// from v0 are d2 and d1 + d2.  Returns true on a hit on [start, end].  On a
// miss intersection and normal are still assigned, with unspecified values.
// The normal is the unnormalized geometric normal.
static bool RayTriangleTestDouble(const vec3d &start, const vec3d &end,
								  const vec3d &v0, const vec3d &d1,
								  const vec3d &d2, vec3d &intersection,
								  vec3d &normal)
{
	auto comp = [](const vec3d &v, int axis) -> double {
		return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
	};

	const vec3d ray = {end.x - start.x, end.y - start.y, end.z - start.z};

	// Pick the dominant ray axis, move it to z.
	const double ax = std::fabs(ray.x);
	const double ay = std::fabs(ray.y);
	const double az = std::fabs(ray.z);
	int kz = (ax >= ay && ax >= az) ? 0 : (ay >= az ? 1 : 2);
	int kx = (kz + 1) % 3;
	int ky = (kx + 1) % 3;

	// Swap kx/ky to preserve winding when the dominant direction is negative.
	if (comp(ray, kz) < 0.0) {
		const int tmp = kx;
		kx = ky;
		ky = tmp;
	}

	const double invDirZ = 1.0 / comp(ray, kz);
	const double Sx = comp(ray, kx) * invDirZ;
	const double Sy = comp(ray, ky) * invDirZ;
	const double Sz = invDirZ;

	const vec3d e1 = {d1.x + d2.x, d1.y + d2.y, d1.z + d2.z};
	const vec3d e2 = d2;
	const vec3d A = {v0.x - start.x, v0.y - start.y, v0.z - start.z};
	const vec3d B = {A.x + e1.x, A.y + e1.y, A.z + e1.z};
	const vec3d C = {A.x + e2.x, A.y + e2.y, A.z + e2.z};

	// Shear and scale so the ray starts at the origin, pointing along +Z.
	const double Ax = comp(A, kx) - Sx * comp(A, kz);
	const double Ay = comp(A, ky) - Sy * comp(A, kz);
	const double Bx = comp(B, kx) - Sx * comp(B, kz);
	const double By = comp(B, ky) - Sy * comp(B, kz);
	const double Cx = comp(C, kx) - Sx * comp(C, kz);
	const double Cy = comp(C, ky) - Sy * comp(C, kz);

	// Scaled barycentric coordinates (2D edge functions at the origin).
	const double U = Bx * Cy - By * Cx;
	const double V = Cx * Ay - Cy * Ax;
	const double W = Ax * By - Ay * Bx;
	const double det = U + V + W;

	const bool mixedSigns = ((U < 0.0) || (V < 0.0) || (W < 0.0)) &&
							((U > 0.0) || (V > 0.0) || (W > 0.0));
	const double T = U * (comp(A, kz) * Sz) + V * (comp(B, kz) * Sz) +
					 W * (comp(C, kz) * Sz);

	const bool positiveDet = det > 0.0;
	const bool validDet = (det != 0.0) && !mixedSigns;
	const bool validT =
		positiveDet ? (T >= 0.0 && T <= det) : (T <= 0.0 && T >= det);
	const bool hit = validDet && validT;

	const double t = T / det;
	intersection = {start.x + ray.x * t, start.y + ray.y * t,
					start.z + ray.z * t};
	normal = {e1.y * e2.z - e1.z * e2.y, e1.z * e2.x - e1.x * e2.z,
			  e1.x * e2.y - e1.y * e2.x};
	return hit;
}

// Lanes whose single-precision result is numerically ambiguous (degenerate
// or near-degenerate determinant, or a barycentric coordinate within the
// rounding error of a triangle edge) are recomputed by RayTriangleTestDouble.
uint4 RayTriangleTest(vec3f start, vec3f end, vec3f v0, vec3f d1, vec3f d2,
					  vec3f &intersection, vec3f &normal)
{
	const float4 zero = float4::load(0.0f);
	const vec3f ray = end - start;

	// Select the dominant ray axis per lane (first step of the shear).
	const float4 ax = max(ray.x, zero - ray.x);
	const float4 ay = max(ray.y, zero - ray.y);
	const float4 az = max(ray.z, zero - ray.z);
	const uint4 useX = (ax >= ay) & (ax >= az);
	const uint4 useY = (~useX) & (ay >= az);

	float4 kx = select(useX, select(useY, ray.x, ray.z), ray.y);
	float4 ky = select(useX, select(useY, ray.y, ray.x), ray.z);
	float4 kz = select(useX, select(useY, ray.z, ray.y), ray.x);

	// Keep the projected orientation consistent for negative dominant
	// direction, per lane.
	const uint4 reverse = kz < zero;
	const float4 oldKx = kx;
	kx = select(reverse, kx, ky);
	ky = select(reverse, ky, oldKx);

	const float4 sx = zero - kx / kz;
	const float4 sy = zero - ky / kz;
	const float4 sz = float4::load(1.0f) / kz;

	const vec3f e1 = d1 + d2;
	const vec3f e2 = d2;
	const vec3f p0 = v0 - start;
	const vec3f p1 = p0 + e1;
	const vec3f p2 = p0 + e2;

	// Project and shear the triangle into the dominant-axis plane.
	const float4 p0x = select(useX, select(useY, p0.x, p0.z), p0.y);
	const float4 p0y = select(useX, select(useY, p0.y, p0.x), p0.z);
	const float4 p0z = select(useX, select(useY, p0.z, p0.y), p0.x);
	const float4 p1x = select(useX, select(useY, p1.x, p1.z), p1.y);
	const float4 p1y = select(useX, select(useY, p1.y, p1.x), p1.z);
	const float4 p1z = select(useX, select(useY, p1.z, p1.y), p1.x);
	const float4 p2x = select(useX, select(useY, p2.x, p2.z), p2.y);
	const float4 p2y = select(useX, select(useY, p2.y, p2.x), p2.z);
	const float4 p2z = select(useX, select(useY, p2.z, p2.y), p2.x);

	// For negative dominant direction kx/ky are swapped to preserve the
	// winding; sx/sy already came from the swapped direction, so only the
	// projected x/y swap here.
	const float4 q0x = select(reverse, p0x + sx * p0z, p0y + sx * p0z);
	const float4 q0y = select(reverse, p0y + sy * p0z, p0x + sy * p0z);
	const float4 q1x = select(reverse, p1x + sx * p1z, p1y + sx * p1z);
	const float4 q1y = select(reverse, p1y + sy * p1z, p1x + sy * p1z);
	const float4 q2x = select(reverse, p2x + sx * p2z, p2y + sx * p2z);
	const float4 q2y = select(reverse, p2y + sy * p2z, p2x + sy * p2z);

	// Oriented edge functions.  Mixed-sign rejection is the watertight
	// inside test; zero is accepted.
	const float4 U = q1x * q2y - q1y * q2x;
	const float4 V = q2x * q0y - q2y * q0x;
	const float4 W = q0x * q1y - q0y * q1x;
	const float4 det = U + V + W;
	const uint4 mixedSigns = ((U < zero) | (V < zero) | (W < zero)) &
							 ((U > zero) | (V > zero) | (W > zero));

	const float4 T = U * (p0z * sz) + V * (p1z * sz) + W * (p2z * sz);
	const uint4 positiveDet = det > zero;
	const uint4 validDet = (det != zero) & ~mixedSigns;
	const uint4 validT = (positiveDet & (T >= zero) & (T <= det)) |
						 ((~positiveDet) & (T <= zero) & (T >= det));
	const uint4 hit = validDet & validT;

	const float4 t = T / det;
	intersection = start + ray * t;
	normal = e1.cross(e2);

	// Flag lanes whose single-precision decision is within kFallbackEps of a
	// sign flip; those are recomputed exactly.
	const float4 kFallbackEps = float4::load(1.0e-6f);
	auto absF = [&](float4 v) { return max(v, zero - v); };

	const float4 uErr = kFallbackEps * (absF(q1x * q2y) + absF(q1y * q2x));
	const float4 vErr = kFallbackEps * (absF(q2x * q0y) + absF(q2y * q0x));
	const float4 wErr = kFallbackEps * (absF(q0x * q1y) + absF(q0y * q1x));
	const float4 detErr = uErr + vErr + wErr + kFallbackEps * absF(det);

	const float4 T0 = U * (p0z * sz);
	const float4 T1 = V * (p1z * sz);
	const float4 T2 = W * (p2z * sz);
	const float4 tErr = kFallbackEps * (absF(T0) + absF(T1) + absF(T2));

	const uint4 edgeUncertain =
		(absF(U) <= uErr) | (absF(V) <= vErr) | (absF(W) <= wErr);
	const uint4 detUncertain = (det == zero) | (absF(det) <= detErr);
	const uint4 tUncertain =
		(absF(T) <= tErr) | (absF(T - det) <= (tErr + detErr));
	const uint4 uncertain = edgeUncertain | detUncertain | tUncertain;

	uint32_t uncertainLanes[4];
	uncertain.store(uncertainLanes);
	const bool anyUncertain = (uncertainLanes[0] | uncertainLanes[1] |
							   uncertainLanes[2] | uncertainLanes[3]);

	uint4 result = hit;
	if (anyUncertain) {
		float sxa[4], sya[4], sza[4];
		float exa[4], eya[4], eza[4];
		float v0xa[4], v0ya[4], v0za[4];
		float d1xa[4], d1ya[4], d1za[4];
		float d2xa[4], d2ya[4], d2za[4];
		start.store(sxa, sya, sza);
		end.store(exa, eya, eza);
		v0.store(v0xa, v0ya, v0za);
		d1.store(d1xa, d1ya, d1za);
		d2.store(d2xa, d2ya, d2za);

		float ixa[4], iya[4], iza[4];
		float nxa[4], nya[4], nza[4];
		intersection.store(ixa, iya, iza);
		normal.store(nxa, nya, nza);

		uint32_t hitLanes[4];
		hit.store(hitLanes);

		for (int i = 0; i < 4; ++i) {
			if (uncertainLanes[i] == 0u)
				continue;

			const vec3d s = {sxa[i], sya[i], sza[i]};
			const vec3d e = {exa[i], eya[i], eza[i]};
			const vec3d V0 = {v0xa[i], v0ya[i], v0za[i]};
			const vec3d D1 = {d1xa[i], d1ya[i], d1za[i]};
			const vec3d D2 = {d2xa[i], d2ya[i], d2za[i]};

			vec3d I, N;
			if (RayTriangleTestDouble(s, e, V0, D1, D2, I, N)) {
				hitLanes[i] = 0xFFFFFFFFu;
				ixa[i] = (float)I.x;
				iya[i] = (float)I.y;
				iza[i] = (float)I.z;
				nxa[i] = (float)N.x;
				nya[i] = (float)N.y;
				nza[i] = (float)N.z;
			} else {
				hitLanes[i] = 0u;
			}
		}

		intersection = vec3f::load(ixa, iya, iza);
		normal = vec3f::load(nxa, nya, nza);
		result = uint4::load(hitLanes);
	}

	return result;
}

} // namespace simd
