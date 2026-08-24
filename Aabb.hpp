// This file is part of VecSIMD
// Copyright (c) 2026 Marek Zalewski aka Drwalin
// You should have received a copy of the MIT License along with this program.

#pragma once

#include "Vec3.hpp"

namespace simd
{
inline const float EPSILON = 0.000001f;
inline const float BIG_EPSILON = 0.02f;

struct aabb4 {
	simd::vec3f min;
	simd::vec3f max;

public:
	uint4 IsValid() const;

	float4 GetVolume() const;
	float4 GetSurface() const;

	simd::vec3f GetCenter() const;
	// returns max-min
	simd::vec3f GetSizes() const;
	simd::vec3f GetMin() const;
	simd::vec3f GetMax() const;

	aabb4 Expanded(float by) const;

	uint4 HasIntersection(const aabb4 &r) const;
	uint4 HasIntersection(const aabb4 &r, float eps) const;
	uint4 HasIntersection(const aabb4 &r, float4 eps) const;

	uint4 IsIn(const simd::vec3f &r) const;
	uint4 IsIn(const simd::vec3f &r, float eps) const;
	uint4 IsIn(const simd::vec3f &r, float4 eps) const;

	aabb4 Intersection(const aabb4 &r) const;
	aabb4 Sum(const aabb4 &r) const;
	aabb4 Sum(const simd::vec3f &r) const;

	uint4 ContainsAll(const aabb4 &r) const;
	uint4 ContainsAll(const aabb4 &r, float eps) const;
	uint4 ContainsAll(const aabb4 &r, float4) const;

	uint4 FastRayTest2(const simd::vec3f &ro, const simd::vec3f &invDir,
					   const int4 raySign[3], const uint4 parallel[3],
					   float4 &near, float4 &far) const;
	uint4 FastRayTest2(const simd::vec3f &ro, const simd::vec3f &invDir,
					   const uint4 parallel[3],
					   float4 &near, float4 &far) const;

	uint4 SlowRayTest2(const simd::vec3f start, const simd::vec3f end,
					   float4 &near, float4 &far) const;

	// returns:
	//  new_min = parent_min + (parent_extent * (current_min * scale_inversion))
	//  new_max = new_min + ((parent_max - new_min) * (current_max *
	//                scale_inversion)
	aabb4 GetChildRelativeAabb(aabb4 scaled_child) const;

public:
	// returns ::HasIntersection
	uint4 operator&&(const aabb4 &r) const;
	uint4 operator&&(const simd::vec3f &r) const;
	// returns ::Intersection
	aabb4 operator*(const aabb4 &r) const;
	// returns ::Sum
	aabb4 operator+(const aabb4 &r) const;
	aabb4 operator+(const simd::vec3f &r) const;

	uint4 operator==(const aabb4 &r) const;
	uint4 operator!=(const aabb4 &r) const;
};
} // namespace simd

#include "Aabb.impl.hpp" // IWYU pragma: export
#include "Aabb_ray.impl.hpp" // IWYU pragma: export
