// Copyright (c) 2026 Marek Zalewski aka Drwalin
// Distributed under the MIT License.

#pragma once

#include "Vec3.hpp"

namespace simd
{
// Tests four segment rays against four triangles at once.
// d1 = v1 - v2 and d2 = v2 - v0, so the edges from v0 are d2 and d1 + d2.
// The returned mask has one all-bits-set lane per hit.  On a miss
// intersection and normal are still assigned, with unspecified values.  The
// normal is the unnormalized geometric normal.
uint4 RayTriangleTest(vec3f start, vec3f end, vec3f v0, vec3f d1, vec3f d2,
					  vec3f &intersection, vec3f &normal);
}
