// This file is part of VecSIMD
// Copyright (c) 2026 Marek Zalewski aka Drwalin
// Distributed under the MIT License.

#include "Vec3.hpp"

namespace simd
{

// Tests four segment rays against four triangles at once.
//
// d1 is v1 - v2 and d2 is v2 - v0, hence the two edges from v0 are
// d2 and d1 + d2.  The returned mask has one all-bits-set lane per hit.
// On a miss, intersection and normal are still assigned, but their values
// are unspecified.  The normal is the (unnormalized) geometric normal.
uint4 RayTriangleTest(vec3f start, vec3f end, vec3f v0, vec3f d1,
					  vec3f d2, vec3f &intersection, vec3f &normal)
{
	const float4 zero = float4::load(0.0f);
	const vec3f ray = end - start;

	// Select the dominant ray axis independently in every SIMD lane.  The
	// axis permutation is the first step of the watertight shear transform.
	const float4 ax = max(ray.x, zero - ray.x);
	const float4 ay = max(ray.y, zero - ray.y);
	const float4 az = max(ray.z, zero - ray.z);
	const uint4 useX = (ax >= ay) & (ax >= az);
	const uint4 useY = (~useX) & (ay >= az);

	float4 kx = select(useX, select(useY, ray.x, ray.z), ray.y);
	float4 ky = select(useX, select(useY, ray.y, ray.x), ray.z);
	float4 kz = select(useX, select(useY, ray.z, ray.y), ray.x);

	// Keep the projected triangle orientation consistent when the dominant
	// direction is negative.  This is also done independently per lane.
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

	// Project and shear the triangle into the dominant-axis coordinate plane.
	const float4 p0x = select(useX, select(useY, p0.x, p0.z), p0.y);
	const float4 p0y = select(useX, select(useY, p0.y, p0.x), p0.z);
	const float4 p0z = select(useX, select(useY, p0.z, p0.y), p0.x);
	const float4 p1x = select(useX, select(useY, p1.x, p1.z), p1.y);
	const float4 p1y = select(useX, select(useY, p1.y, p1.x), p1.z);
	const float4 p1z = select(useX, select(useY, p1.z, p1.y), p1.x);
	const float4 p2x = select(useX, select(useY, p2.x, p2.z), p2.y);
	const float4 p2y = select(useX, select(useY, p2.y, p2.x), p2.z);
	const float4 p2z = select(useX, select(useY, p2.z, p2.y), p2.x);

	const float4 q0x = select(reverse, p0x + sx * p0z, p0y + sy * p0z);
	const float4 q0y = select(reverse, p0y + sy * p0z, p0x + sx * p0z);
	const float4 q1x = select(reverse, p1x + sx * p1z, p1y + sy * p1z);
	const float4 q1y = select(reverse, p1y + sy * p1z, p1x + sx * p1z);
	const float4 q2x = select(reverse, p2x + sx * p2z, p2y + sy * p2z);
	const float4 q2y = select(reverse, p2y + sy * p2z, p2x + sx * p2z);

	// Three oriented edge functions.  The mixed-sign rejection is the
	// watertight inside-triangle test; zero is deliberately accepted.
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
	return hit;
}

} // namespace simd
