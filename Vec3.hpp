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

	static ubyte16 load(const void *v);
	static ubyte16 load4(const void *v);

	void store(void *v) const;

	ubyte16 operator+(ubyte16 other) const;
	ubyte16 operator-(ubyte16 other) const;
	ubyte16 operator*(ubyte16 other) const;

	ubyte16 operator&(ubyte16 other) const;
	ubyte16 operator|(ubyte16 other) const;
	ubyte16 operator^(ubyte16 other) const;
	ubyte16 operator~() const;
};

struct byte16 {
	VECI v;

	static byte16 load(const void *v);
	static byte16 load4(const void *v);

	void store(void *v) const;

	byte16 operator+(byte16 other) const;
	byte16 operator-(byte16 other) const;
	byte16 operator*(byte16 other) const;

	byte16 operator&(byte16 other) const;
	byte16 operator|(byte16 other) const;
	byte16 operator^(byte16 other) const;
	byte16 operator~() const;
};

struct uword8 {
	VECU v;

	static uword8 load(const void *v);
	static uword8 load4(const void *v);

	void store(void *v) const;

	uword8 operator+(uword8 other) const;
	uword8 operator-(uword8 other) const;
	uword8 operator*(uword8 other) const;

	uword8 operator&(uword8 other) const;
	uword8 operator|(uword8 other) const;
	uword8 operator^(uword8 other) const;
	uword8 operator~() const;
};

struct word8 {
	VECI v;

	static word8 load(const void *v);
	static word8 load4(const void *v);

	void store(void *v) const;

	word8 operator+(word8 other) const;
	word8 operator-(word8 other) const;
	word8 operator*(word8 other) const;

	word8 operator&(word8 other) const;
	word8 operator|(word8 other) const;
	word8 operator^(word8 other) const;
	word8 operator~() const;
};

struct uint4 {
	VECU v;

	static uint4 load(const void *v);
	static uint4 convert_from(ubyte16 src);
	static uint4 convert_from(uword8 src);
	static uint4 reinterpret(int4 src);

	void store(void *v) const;

	uint4 operator+(uint4 other) const;
	uint4 operator-(uint4 other) const;
	uint4 operator*(uint4 other) const;

	uint4 operator&(uint4 other) const;
	uint4 operator|(uint4 other) const;
	uint4 operator^(uint4 other) const;
	uint4 operator~() const;

	uint4 operator<<(int bits) const;
	uint4 operator>>(int bits) const;

	template <int bits> uint4 shl() const;
	template <int bits> uint4 shr() const;
};

struct int4 {
	VECI v;

	static int4 load(const void *v);
	static int4 convert_from(byte16 src);
	static int4 convert_from(word8 src);
	static int4 reinterpret(uint4 src);

	void store(void *v) const;

	int4 operator+(int4 other) const;
	int4 operator-(int4 other) const;
	int4 operator*(int4 other) const;

	int4 operator&(int4 other) const;
	int4 operator|(int4 other) const;
	int4 operator^(int4 other) const;
	int4 operator~() const;
};

struct float4 {
	VECF v;

	static float4 load(const float *v);
	static float4 load(float v1, float v2, float v3, float v4);
	static float4 load(float v);
	static float4 convert_from(int4 src);

	void store(float *v) const;

	float4 operator-(int) const;

	float4 operator+(float4 other) const;
	float4 operator-(float4 other) const;
	float4 operator*(float4 other) const;
	float4 operator/(float4 other) const;

	float4 sqrt() const;

	float4 round_to_nearest() const;
	float4 floor() const;
	float4 ceil() const;
	float4 round_to_zero() const;

	uint4 operator<(float4 r) const;
	uint4 operator>(float4 r) const;
	uint4 operator<=(float4 r) const;
	uint4 operator>=(float4 r) const;
	uint4 operator==(float4 r) const;
	uint4 operator!=(float4 r) const;
};

struct vec3f {
	float4 x, y, z;

	static vec3f load(const float *x, const float *y, const float *z);
	static vec3f load(float x, float y, float z);

	void store(float *x, float *y, float *z) const;

	vec3f operator-(int) const;

	vec3f operator+(vec3f other) const;
	vec3f operator-(vec3f other) const;
	vec3f operator*(vec3f other) const;
	vec3f operator/(vec3f other) const;

	vec3f operator+(float4 other) const;
	vec3f operator-(float4 other) const;
	vec3f operator*(float4 other) const;
	vec3f operator/(float4 other) const;

	float4 dot(vec3f other) const;
	vec3f cross(vec3f other) const;
	float4 length2() const;
	float4 length() const;
};

uint4 min(uint4 l, uint4 r);
uint4 max(uint4 l, uint4 r);
int4 min(int4 l, int4 r);
int4 max(int4 l, int4 r);
float4 min(float4 l, float4 r);
float4 max(float4 l, float4 r);
float4 select(uint4 mask, float4 falseValue, float4 trueValue);
vec3f min(vec3f l, vec3f r);
vec3f max(vec3f l, vec3f r);
} // namespace simd
