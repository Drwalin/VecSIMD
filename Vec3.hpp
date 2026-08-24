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
	inline static ubyte16 load4(const void *v);
};

struct byte16 {
	VECI v;

	inline static byte16 load(const void *v);
	inline static byte16 load4(const void *v);
};

struct uword8 {
	VECU v;

	inline static uword8 load(const void *v);
	inline static uword8 load4(const void *v);
};

struct word8 {
	VECI v;

	inline static word8 load(const void *v);
	inline static word8 load4(const void *v);
};

struct uint4 {
	VECU v;

	inline static uint4 load(const void *v);
	inline static uint4 convert_from(ubyte16 src);
	inline static uint4 convert_from(uword8 src);
	inline static uint4 reinterpret(int4 src);

	inline uint4 operator&(uint4 other) const;
	inline uint4 operator|(uint4 other) const;
	inline uint4 operator^(uint4 other) const;
	inline uint4 operator~() const;

	inline uint4 operator<<(int bits) const;
	inline uint4 operator>>(int bits) const;

	template <int bits> inline uint4 shl() const;
	template <int bits> inline uint4 shr() const;
};

struct int4 {
	VECI v;

	inline static int4 load(const void *v);
	inline static int4 convert_from(byte16 src);
	inline static int4 convert_from(word8 src);
	inline static int4 reinterpret(uint4 src);

	inline int4 operator+(int4 other) const;
	inline int4 operator-(int4 other) const;
	inline int4 operator*(int4 other) const;
};

struct float4 {
	VECF v;

	inline static float4 load(const float *v);
	inline static float4 load(float v1, float v2, float v3, float v4);
	inline static float4 load(float v);
	inline static float4 convert_from(int4 src);

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

inline uint4 min(uint4 l, uint4 r);
inline uint4 max(uint4 l, uint4 r);
inline int4 min(int4 l, int4 r);
inline int4 max(int4 l, int4 r);
inline float4 min(float4 l, float4 r);
inline float4 max(float4 l, float4 r);
inline vec3f min(vec3f l, vec3f r);
inline vec3f max(vec3f l, vec3f r);
} // namespace simd

#include "Vec3.impl.hpp" // IWYU pragma: export
