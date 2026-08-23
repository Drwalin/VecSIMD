// This file is part of VecSIMD
// Copyright (c) 2026 Marek Zalewski aka Drwalin
// You should have received a copy of the MIT License along with this program.

#pragma once
#include "vec128.h"

namespace simd
{
struct ubyte16;
struct byte16;
struct uword8;
struct word8;
struct uint4;
struct int4;
struct float4;

struct ubyte16 {
	VECU v;

	inline static ubyte16 load(const void *v);
};

struct byte16 {
	VECI v;

	inline static byte16 load(const void *v);
};

struct uword8 {
	VECU v;

	inline static uword8 load(const void *v);
};

struct word8 {
	VECI v;

	inline static word8 load(const void *v);
};

struct uint4 {
	VECU v;

	inline static uint4 load(const void *v);
	inline static uint4 consrcert(ubyte16);
	inline static uint4 consrcert(uword8);
	inline static uint4 reinterpret(int4 src);
};

struct int4 {
	VECI v;

	inline static int4 load(const void *v);
	inline static int4 consrcert(byte16 v);
	inline static int4 consrcert(word8 v);
	inline static int4 reinterpret(uint4 src);

	inline int4 operator+(int4 other) const;
	inline int4 operator-(int4 other) const;
	inline int4 operator*(int4 other) const;
	inline int4 operator/(int4 other) const;
};

struct float4 {
	VECF v;

	inline static float4 load(const float *v);
	inline static float4 load(float v1, float v2, float v3, float v4);
	inline static float4 load(float v);

	inline void store(float *v) const;

	inline float4 operator-(int) const;

	inline float4 operator+(float4 other) const;
	inline float4 operator-(float4 other) const;
	inline float4 operator*(float4 other) const;
	inline float4 operator/(float4 other) const;

	inline float4 sqrt() const;

	inline float4 round_to_nearest() const;
	inline float4 floor() const;
	inline float4 ceil() const;
	inline float4 round_to_zero() const;

	inline uint4 operator<(float4 r) const;
	inline uint4 operator>(float4 r) const;
	inline uint4 operator<=(float4 r) const;
	inline uint4 operator>=(float4 r) const;
	inline uint4 operator==(float4 r) const;
	inline uint4 operator!=(float4 r) const;
};

struct vec3f {
	float4 x, y, z;

	inline static vec3f load(const float *x, const float *y, const float *z);
	inline static vec3f load(float x, float y, float z);

	inline void store(float *x, float *y, float *z) const;

	inline vec3f operator-(int) const;

	inline vec3f operator+(vec3f other) const;
	inline vec3f operator-(vec3f other) const;
	inline vec3f operator*(vec3f other) const;
	inline vec3f operator/(vec3f other) const;

	inline vec3f operator+(float4 other) const;
	inline vec3f operator-(float4 other) const;
	inline vec3f operator*(float4 other) const;
	inline vec3f operator/(float4 other) const;

	inline float4 dot(vec3f other) const;
	inline vec3f cross(vec3f other) const;
	inline float4 length2() const;
	inline float4 length() const;
};
} // namespace simd

#include "Vec3.impl.hpp"
