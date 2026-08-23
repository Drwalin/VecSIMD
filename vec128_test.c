#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vec128.h"

static int fails = 0;

#define CHECK(cond)                                                            \
	do {                                                                       \
		if (!(cond)) {                                                         \
			printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);             \
			fails++;                                                           \
		}                                                                      \
	} while (0)

static void test_convert(void)
{
	const int8_t b[4] = {5, -2, -128, 127};
	const uint8_t ub[4] = {5, 250, 128, 255};
	int32_t r[4];
	uint32_t ru[4];
	VEC_STORE_I(r, VEC_BI_TO_I(b));
	VEC_STORE_U(ru, VEC_BU_TO_U(ub));
	for (int i = 0; i < 4; i++) {
		CHECK(r[i] == b[i]);
		CHECK(ru[i] == ub[i]);
	}

	const int32_t si[4] = {0, 1, -1, 2147483647};
	const uint32_t su[4] = {0, 1, 2147483648u, 4294967295u};
	float rf[4];
	VEC_STORE_F(rf, VEC_ITOF(VEC_LOAD_I(si)));
	CHECK(rf[0] == 0.0f && rf[1] == 1.0f && rf[2] == -1.0f &&
		  rf[3] == 2147483648.0f);
	VEC_STORE_F(rf, VEC_UTOF(VEC_LOAD_U(su)));
	CHECK(rf[0] == 0.0f && rf[1] == 1.0f && rf[2] == 2147483648.0f &&
		  rf[3] == 4294967296.0f);

	const int16_t s16[4] = {5, -2, -32768, 32767};
	const uint16_t us16[4] = {5, 250, 32768, 65535};
	VEC_STORE_I(r, VEC_BI16_TO_I(s16));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == s16[i]);
	VEC_STORE_U(ru, VEC_BU16_TO_U(us16));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == us16[i]);
}

static void test_arith(void)
{
	const float fa[4] = {1.5f, -2.0f, 3.25f, 100.0f};
	const float fb[4] = {2.0f, 0.5f, -4.0f, 10.0f};
	const int32_t ia[4] = {10, -20, 300, 7};
	const int32_t ib[4] = {3, 5, -10, 2};
	const uint32_t ua[4] = {10, 20, 300, 4294967295u};
	const uint32_t ub[4] = {3, 5, 7, 1};
	float rf[4];
	int32_t ri[4];
	uint32_t ru[4];

	VEC_STORE_F(rf, VEC_ADD_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	for (int i = 0; i < 4; i++)
		CHECK(rf[i] == fa[i] + fb[i]);
	VEC_STORE_F(rf, VEC_SUB_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	for (int i = 0; i < 4; i++)
		CHECK(rf[i] == fa[i] - fb[i]);
	VEC_STORE_F(rf, VEC_MUL_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	for (int i = 0; i < 4; i++)
		CHECK(rf[i] == fa[i] * fb[i]);
	VEC_STORE_F(rf, VEC_DIV_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	for (int i = 0; i < 4; i++)
		CHECK(rf[i] == fa[i] / fb[i]);

	VEC_STORE_I(ri, VEC_ADD_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	for (int i = 0; i < 4; i++)
		CHECK(ri[i] == ia[i] + ib[i]);
	VEC_STORE_I(ri, VEC_SUB_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	for (int i = 0; i < 4; i++)
		CHECK(ri[i] == ia[i] - ib[i]);
	VEC_STORE_I(ri, VEC_MUL_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	for (int i = 0; i < 4; i++)
		CHECK(ri[i] == ia[i] * ib[i]);

	VEC_STORE_U(ru, VEC_ADD_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ua[i] + ub[i]);
	VEC_STORE_U(ru, VEC_SUB_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ua[i] - ub[i]);
	VEC_STORE_U(ru, VEC_MUL_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ua[i] * ub[i]);
}

static void test_sign(void)
{
	const float fa[4] = {1.5f, -2.0f, 0.0f, -0.0f};
	const int32_t ia[4] = {7, -3, 0, -1};
	const uint32_t ua[4] = {5, 0, 300, 1};
	float rf[4];
	int32_t ri[4];
	uint32_t ru[4];

	VEC_STORE_F(rf, VEC_SIGN_F(VEC_LOAD_F(fa)));
	CHECK(rf[0] == 1.0f && rf[1] == -1.0f && rf[2] == 0.0f && rf[3] == 0.0f);
	VEC_STORE_F(rf, VEC_SIGNBIT_F(VEC_LOAD_F(fa)));
	CHECK(!signbit(rf[0]) && signbit(rf[1]) && !signbit(rf[2]) &&
		  signbit(rf[3]));

	VEC_STORE_I(ri, VEC_SIGN_I(VEC_LOAD_I(ia)));
	CHECK(ri[0] == 1 && ri[1] == -1 && ri[2] == 0 && ri[3] == -1);
	VEC_STORE_I(ri, VEC_SIGNBIT_I(VEC_LOAD_I(ia)));
	CHECK(ri[0] == 0 && ri[1] == (int32_t)0x80000000 && ri[2] == 0 &&
		  ri[3] == (int32_t)0x80000000);

	VEC_STORE_U(ru, VEC_SIGN_U(VEC_LOAD_U(ua)));
	CHECK(ru[0] == 1 && ru[1] == 0 && ru[2] == 1 && ru[3] == 1);
	VEC_STORE_U(ru, VEC_SIGNBIT_U(VEC_LOAD_U(ua)));
	CHECK(ru[0] == 0 && ru[1] == 0 && ru[2] == 0 && ru[3] == 0);
	{
		const uint32_t u2[4] = {0x80000000u, 5, 0xFFFFFFFFu, 0};
		VEC_STORE_U(ru, VEC_SIGNBIT_U(VEC_LOAD_U(u2)));
		CHECK(ru[0] == 0x80000000u && ru[1] == 0 && ru[2] == 0x80000000u &&
			  ru[3] == 0);
	}
}

static void test_cmp(void)
{
	const float fa[4] = {1.0f, 5.0f, 3.0f, 9.0f};
	const float fb[4] = {2.0f, 4.0f, 6.0f, 8.0f};
	const int32_t ia[4] = {1, -5, 3, 9};
	const int32_t ib[4] = {2, -4, 6, 8};
	const uint32_t ua[4] = {1, 5, 0x80000000u, 9};
	const uint32_t ub[4] = {2, 4, 7, 8};
	int32_t r[4];
	uint32_t ru[4];

	VEC_STORE_I(r, VEC_F2I(VEC_CMP_GT_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb))));
	CHECK(r[0] == 0 && r[1] == -1 && r[2] == 0 && r[3] == -1);
	VEC_STORE_I(r, VEC_F2I(VEC_CMP_LT_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb))));
	CHECK(r[0] == -1 && r[1] == 0 && r[2] == -1 && r[3] == 0);
	VEC_STORE_I(r, VEC_F2I(VEC_CMP_GE_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb))));
	CHECK(r[0] == 0 && r[1] == -1 && r[2] == 0 && r[3] == -1);
	VEC_STORE_I(r, VEC_F2I(VEC_CMP_LE_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb))));
	CHECK(r[0] == -1 && r[1] == 0 && r[2] == -1 && r[3] == 0);
	VEC_STORE_I(r, VEC_F2I(VEC_CMP_EQ_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb))));
	CHECK(r[0] == 0 && r[1] == 0 && r[2] == 0 && r[3] == 0);
	{
		const float fc[4] = {1.0f, 5.0f, 6.0f, 8.0f};
		VEC_STORE_I(r, VEC_F2I(VEC_CMP_NE_F(VEC_LOAD_F(fa), VEC_LOAD_F(fc))));
		CHECK(r[0] == 0 && r[1] == 0 && r[2] == -1 && r[3] == -1);
	}

	VEC_STORE_I(r, VEC_CMP_GT_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == 0 && r[1] == 0 && r[2] == 0 && r[3] == -1);
	VEC_STORE_I(r, VEC_CMP_LT_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == -1 && r[1] == -1 && r[2] == -1 && r[3] == 0);
	VEC_STORE_I(r, VEC_CMP_GE_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == 0 && r[1] == 0 && r[2] == 0 && r[3] == -1);
	VEC_STORE_I(r, VEC_CMP_LE_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == -1 && r[1] == -1 && r[2] == -1 && r[3] == 0);
	VEC_STORE_I(r, VEC_CMP_EQ_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == 0 && r[1] == 0 && r[2] == 0 && r[3] == 0);
	VEC_STORE_I(r, VEC_CMP_NE_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)));
	CHECK(r[0] == -1 && r[1] == -1 && r[2] == -1 && r[3] == -1);

	VEC_STORE_U(ru, VEC_CMP_GT_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	CHECK(ru[0] == 0 && ru[1] == -1u && ru[2] == -1u && ru[3] == -1u);
	VEC_STORE_U(ru, VEC_CMP_LT_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	CHECK(ru[0] == -1u && ru[1] == 0 && ru[2] == 0 && ru[3] == 0);
	VEC_STORE_U(ru, VEC_CMP_GE_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	CHECK(ru[0] == 0 && ru[1] == -1u && ru[2] == -1u && ru[3] == -1u);
	VEC_STORE_U(ru, VEC_CMP_LE_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	CHECK(ru[0] == -1u && ru[1] == 0 && ru[2] == 0 && ru[3] == 0);
	VEC_STORE_U(ru, VEC_CMP_EQ_U(VEC_LOAD_U(ua), VEC_LOAD_U(ub)));
	CHECK(ru[0] == 0 && ru[1] == 0 && ru[2] == 0 && ru[3] == 0);
	{
		const uint32_t uc[4] = {1, 5, 7, 8};
		VEC_STORE_U(ru, VEC_CMP_NE_U(VEC_LOAD_U(ua), VEC_LOAD_U(uc)));
		CHECK(ru[0] == 0 && ru[1] == 0 && ru[2] == -1u && ru[3] == -1u);
	}
}

static void test_logic(void)
{
	const int32_t a[4] = {0x0F0F0F0F, 0x12345678, -1, (int32_t)0x80000000};
	const int32_t b[4] = {(int32_t)0xF0F0F0F0, (int32_t)0x0F0F0F0F,
						  (int32_t)0x00FF00FF, 1};
	int32_t r[4];
	uint32_t ru[4];
	float fa[4], fb[4], rf[4];

	memcpy(fa, a, sizeof(fa));
	memcpy(fb, b, sizeof(fb));

	VEC_STORE_I(r, VEC_AND_I(VEC_LOAD_I(a), VEC_LOAD_I(b)));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] & b[i]));
	VEC_STORE_I(r, VEC_OR_I(VEC_LOAD_I(a), VEC_LOAD_I(b)));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] | b[i]));
	VEC_STORE_I(r, VEC_XOR_I(VEC_LOAD_I(a), VEC_LOAD_I(b)));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] ^ b[i]));
	VEC_STORE_I(r, VEC_NOT_I(VEC_LOAD_I(a)));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == ~a[i]);

	VEC_STORE_U(ru, VEC_AND_U(VEC_LOAD_U((const uint32_t *)a),
							  VEC_LOAD_U((const uint32_t *)b)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ((uint32_t)a[i] & (uint32_t)b[i]));
	VEC_STORE_U(ru, VEC_OR_U(VEC_LOAD_U((const uint32_t *)a),
							 VEC_LOAD_U((const uint32_t *)b)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ((uint32_t)a[i] | (uint32_t)b[i]));
	VEC_STORE_U(ru, VEC_XOR_U(VEC_LOAD_U((const uint32_t *)a),
							  VEC_LOAD_U((const uint32_t *)b)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ((uint32_t)a[i] ^ (uint32_t)b[i]));
	VEC_STORE_U(ru, VEC_NOT_U(VEC_LOAD_U((const uint32_t *)a)));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == ~(uint32_t)a[i]);

	VEC_STORE_F(rf, VEC_AND_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	memcpy(r, rf, sizeof(r));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] & b[i]));
	VEC_STORE_F(rf, VEC_OR_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	memcpy(r, rf, sizeof(r));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] | b[i]));
	VEC_STORE_F(rf, VEC_XOR_F(VEC_LOAD_F(fa), VEC_LOAD_F(fb)));
	memcpy(r, rf, sizeof(r));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] ^ b[i]));
	VEC_STORE_F(rf, VEC_NOT_F(VEC_LOAD_F(fa)));
	memcpy(r, rf, sizeof(r));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == ~a[i]);
}

static void test_shift(void)
{
	const int32_t a[4] = {1, -1, 0x12345678, 1000};
	const uint32_t u[4] = {1, 0x80000000u, 0xFFFFFFFFu, 0x12345678u};
	int32_t r[4];
	uint32_t ru[4];

	VEC_STORE_I(r, VEC_SHL_I(VEC_LOAD_I(a), 4));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] << 4));
	VEC_STORE_I(r, VEC_SHR_I(VEC_LOAD_I(a), 4));
	for (int i = 0; i < 4; i++)
		CHECK(r[i] == (a[i] >> 4));

	VEC_STORE_U(ru, VEC_SHL_U(VEC_LOAD_U(u), 4));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == (u[i] << 4));
	VEC_STORE_U(ru, VEC_SHR_U(VEC_LOAD_U(u), 4));
	for (int i = 0; i < 4; i++)
		CHECK(ru[i] == (u[i] >> 4));
}

static void test_shuffle(void)
{
	const int32_t a[4] = {10, 20, 30, 40};
	const float fa[4] = {1.0f, 2.0f, 3.0f, 4.0f};
	const uint32_t ua[4] = {10, 20, 30, 40};
	int32_t r[4];
	float rf[4];
	uint32_t ru[4];

	VEC_STORE_I(r, VEC_SHUFFLE_I(VEC_LOAD_I(a), 3, 2, 1, 0));
	CHECK(r[0] == 40 && r[1] == 30 && r[2] == 20 && r[3] == 10);
	VEC_STORE_I(r, VEC_SHUFFLE_I(VEC_LOAD_I(a), 0, 0, 3, 1));
	CHECK(r[0] == 10 && r[1] == 10 && r[2] == 40 && r[3] == 20);

	VEC_STORE_F(rf, VEC_SHUFFLE_F(VEC_LOAD_F(fa), 1, 3, 2, 0));
	CHECK(rf[0] == 2.0f && rf[1] == 4.0f && rf[2] == 3.0f && rf[3] == 1.0f);

	VEC_STORE_U(ru, VEC_SHUFFLE_U(VEC_LOAD_U(ua), 2, 1, 1, 0));
	CHECK(ru[0] == 30 && ru[1] == 20 && ru[2] == 20 && ru[3] == 10);
}

static void test_masks(void)
{
	const float a[4] = {1.0f, 5.0f, 3.0f, 9.0f};
	const float b[4] = {2.0f, 4.0f, 6.0f, 8.0f};
	const float c[4] = {10.0f, 20.0f, 30.0f, 40.0f};
	float r[4];

	VEC_STORE_F(r, VEC_AND_F(VEC_CMP_GT_F(VEC_LOAD_F(a), VEC_LOAD_F(b)),
							 VEC_LOAD_F(c)));
	CHECK(r[0] == 0.0f && r[1] == 20.0f && r[2] == 0.0f && r[3] == 40.0f);
	VEC_STORE_F(r,
				VEC_ADD_F(VEC_LOAD_F(a),
						  VEC_AND_F(VEC_CMP_GT_F(VEC_LOAD_F(a), VEC_LOAD_F(b)),
									VEC_LOAD_F(c))));
	CHECK(r[0] == 1.0f && r[1] == 25.0f && r[2] == 3.0f && r[3] == 49.0f);

	{
		const int32_t ia[4] = {1, 5, 3, 9};
		const int32_t ib[4] = {2, 4, 6, 8};
		const int32_t ic[4] = {10, 20, 30, 40};
		int32_t ri[4];
		VEC_STORE_I(ri, VEC_AND_I(VEC_CMP_GT_I(VEC_LOAD_I(ia), VEC_LOAD_I(ib)),
								  VEC_LOAD_I(ic)));
		CHECK(ri[0] == 0 && ri[1] == 20 && ri[2] == 0 && ri[3] == 40);
		VEC_STORE_I(ri,
					VEC_AND_I(VEC_CMP_GT_U(VEC_LOAD_U((const uint32_t *)ia),
										   VEC_LOAD_U((const uint32_t *)ib)),
							  VEC_LOAD_I(ic)));
		CHECK(ri[0] == 0 && ri[1] == 20 && ri[2] == 0 && ri[3] == 40);
	}
}

static void test_loadu_storeu(void)
{
	float f[8];
	int32_t i[8];
	uint32_t u[8];
	float rf[4];
	int32_t ri[4];
	uint32_t ru[4];
	for (int k = 0; k < 8; k++) {
		f[k] = (float)k + 0.25f;
		i[k] = k - 3;
		u[k] = (uint32_t)(k * 1000u);
	}
	for (int off = 0; off <= 4; off++) {
		VEC_STORE_F(rf, VEC_LOADU_F(f + off));
		for (int k = 0; k < 4; k++)
			CHECK(rf[k] == f[off + k]);
		VEC_STOREU_F(f + off, VEC_LOADU_F(f + off));
		VEC_STORE_I(ri, VEC_LOADU_I(i + off));
		for (int k = 0; k < 4; k++)
			CHECK(ri[k] == i[off + k]);
		VEC_STOREU_I(i + off, VEC_LOADU_I(i + off));
		VEC_STORE_U(ru, VEC_LOADU_U(u + off));
		for (int k = 0; k < 4; k++)
			CHECK(ru[k] == u[off + k]);
		VEC_STOREU_U(u + off, VEC_LOADU_U(u + off));
	}
}

int main(void)
{
	test_convert();
	test_arith();
	test_sign();
	test_cmp();
	test_logic();
	test_shift();
	test_shuffle();
	test_masks();
	test_loadu_storeu();
	if (fails == 0) {
		printf("PASS all\n");
		return 0;
	}
	printf("%d FAILURES\n", fails);
	return 1;
}
