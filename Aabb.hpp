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
	simd::vec3f GetSizes() const;
	simd::vec3f GetMin() const;
	simd::vec3f GetMax() const;

	aabb4 Expanded(float by) const;

	uint4 HasIntersection(const aabb4 &r, float eps = 0.0f) const;

	uint4 IsIn(const simd::vec3f &r, float eps = 0.0f) const;
	aabb4 Intersection(const aabb4 &r) const;
	aabb4 Sum(const aabb4 &r) const;
	aabb4 Sum(const simd::vec3f &r) const;

	uint4 ContainsAll(const aabb4 &r, float eps = 0.0f) const;

	uint4 FastRayTestCenter(const simd::vec3f &ro, const simd::vec3f &rd,
							const simd::vec3f &invDir, float4 length,
							float4 &near, float4 &far) const;
	uint4 SlowRayTestCenter(const simd::vec3f &start, const simd::vec3f &end,
							float4 &near, float4 &far) const;

	uint4 FastRayTest2(const simd::vec3f &ro, const simd::vec3f &invDir,
					   const int4 raySign[3], float4 &near, float4 &far) const;
	uint4 FastRayTest2(const simd::vec3f ro, const simd::vec3f invDir,
					   float4 &near, float4 &far) const;

	uint4 SlowRayTest2(const simd::vec3f start, const simd::vec3f end,
					   float4 &near, float4 &far) const;

	// returns:
	//  new_min = parent_min + (parent_extent * (current_min * scale_inversion))
	//  new_max = new_min + ((parent_max - new_min) * (current_max *
	//                scale_inversion)
	aabb4 GetChildRelativeAabb(aabb4 scaled_child) const;

public:
	// returns wheather both aabb/vec3 have any common part
	uint4 operator&&(const aabb4 &r) const;
	uint4 operator&&(const simd::vec3f &r) const;
	// returns aabb common with passed aabb
	aabb4 operator*(const aabb4 &r) const;
	// returns sum (aabb that contains both this and passed in aabb)
	aabb4 operator+(const aabb4 &r) const;
	// returns sum (aabb that contains both this and passed in vec3)
	aabb4 operator+(const simd::vec3f &r) const;

	uint4 operator==(const aabb4 &r) const;
	uint4 operator!=(const aabb4 &r) const;
};
} // namespace simd
