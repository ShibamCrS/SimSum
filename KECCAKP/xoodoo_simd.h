#include <immintrin.h>
#include "matrix.h"

typedef unsigned long long int ull;
typedef uint32_t U32;
typedef unsigned int ui;
#define STATE_SIZE 12
#define PLANE_SIZE 4
#define SHEET_SIZE 3

typedef uint32_t Lane;
#define LANE 32
#define HALFLANE 16
#define MASK 0xFFFFFFFF
#define HALFMASK 0xFFFF
   
////////////////////     UTILITIES     ////////////////////

void rotate_lanes(ui *plane, int rotate_by)
{
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(plane + i) = ((*(plane + i) << rotate_by) ^ (*(plane + i) >> (32 - rotate_by)));
	}
}

void rotate_sheets(ui *plane, int rotate_by)
{
	ui temp[PLANE_SIZE];
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(temp + i) = *(plane + i);
	}
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(plane + i) = *(temp + ((i + PLANE_SIZE - rotate_by) % PLANE_SIZE));
	}
}

void display_state(ui *A0, ui *A1, ui *A2)
{
	printf("%x %x %x %x\n", A0[0], A0[1], A0[2], A0[3]);
	printf("%x %x %x %x\n", A1[0], A1[1], A1[2], A1[3]);
	printf("%x %x %x %x\n", A2[0], A2[1], A2[2], A2[3]);
	printf("===================================\n");
}

ui reverse(ui a)
{
        ui c = 0;
        for (int l = 0; l < 32; l++)
        {
                ui bit = a >> l & 1;
                c = c ^ (bit << (31 - l));
        }
        return c;
}

void disp(ui *A0, ui *A1, ui *A2)
{
	printf("%x %x %x %x\n", reverse(A2[0]), reverse(A2[1]), reverse(A2[2]), reverse(A2[3]));
	printf("%x %x %x %x\n", reverse(A1[0]), reverse(A1[1]), reverse(A1[2]), reverse(A1[3]));
	printf("%x %x %x %x\n", reverse(A0[0]), reverse(A0[1]), reverse(A0[2]), reverse(A0[3]));
	printf("===================================\n");
}
void display_lanes(const char *text, void *data)
{
    U32 i;
    U32 *state = data;
    printf("%s:\n", text);
    for(U32 i=0; i<12; i++){
        printf("%08X ", (unsigned int)(state[i]));

        if ((i%4) == 3){
            printf("\n");
        }
        else{
            printf(" ");
        }
    }
    printf("\n");
}
//_______________________________________________________//

////////////////////     Theta     ////////////////////

void theta(ui *A0, ui *A1, ui *A2)
{
	ui P1[PLANE_SIZE], P2[PLANE_SIZE], E[PLANE_SIZE];
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(P1 + i) = *(A0 + i) ^ *(A1 + i) ^ *(A2 + i);
		*(P2 + i) = *(A0 + i) ^ *(A1 + i) ^ *(A2 + i);
	}
	rotate_sheets(P1, 1);
	rotate_lanes(P1, 5);
	rotate_sheets(P2, 1);
	rotate_lanes(P2, 14);
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(E + i) = *(P1 + i) ^ *(P2 + i);
	}
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(A0 + i) = *(A0 + i) ^ *(E + i);
		*(A1 + i) = *(A1 + i) ^ *(E + i);
		*(A2 + i) = *(A2 + i) ^ *(E + i);
	}
}

//___________________________________________________//

////////////////////     ThetaInverseNormal      ////////////////////

int countmult(ui a, ui b)
{
        int cnt = 0;
        ui c = reverse(a);
        ui d = c & b;
	cnt = __builtin_popcount(d);
        return cnt;
}

void gf2mult(ui *A0, ui *A1, ui *A2, const ui inv[384][12])
{
	ui output[12];
	for (int i = 0; i < 12; i++)
	{
		ui inter = 0;
		for (int j = 0; j < 32; j++)
		{
			int s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,sum;
			int num = i * 32 + j;
			s0 = countmult(A0[0], inv[num][0]);
			s1 = countmult(A0[1], inv[num][1]);
			s2 = countmult(A0[2], inv[num][2]);
			s3 = countmult(A0[3], inv[num][3]);
			s4 = countmult(A1[0], inv[num][4]);
			s5 = countmult(A1[1], inv[num][5]);
			s6 = countmult(A1[2], inv[num][6]);
			s7 = countmult(A1[3], inv[num][7]);
			s8 = countmult(A2[0], inv[num][8]);
			s9 = countmult(A2[1], inv[num][9]);
			s10 = countmult(A2[2], inv[num][10]);
			s11 = countmult(A2[3], inv[num][11]);
			sum = s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11;
			sum = sum % 2;
			inter = inter ^ (sum << j);
		}
		output[i] = inter;
	}
	A0[0] = output[0];
	A0[1] = output[1];
	A0[2] = output[2];
	A0[3] = output[3];
	A1[0] = output[4];
        A1[1] = output[5];
        A1[2] = output[6];
        A1[3] = output[7];
	A2[0] = output[8];
        A2[1] = output[9];
        A2[2] = output[10];
        A2[3] = output[11];
}

void theta_inv(ui *A0, ui *A1, ui *A2)
{
	gf2mult(A0, A1, A2, theta_inv_matrix);
}

//_________________________________________________________________//

///////////////////     ThetaInverseSIMD     ////////////////////

void gf2mult_simd(ui *A0, ui *A1, ui *A2, const ui inv[384][12])
{
	ui output[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

	__m128i state0 = _mm_setr_epi32(reverse(A0[0]), reverse(A0[1]), reverse(A0[2]), reverse(A0[3]));
	__m128i state1 = _mm_setr_epi32(reverse(A1[0]), reverse(A1[1]), reverse(A1[2]), reverse(A1[3]));
	__m128i state2 = _mm_setr_epi32(reverse(A2[0]), reverse(A2[1]), reverse(A2[2]), reverse(A2[3]));

	__m128i matrix0, matrix1, matrix2;
	__m128i result0, result1, result2;
	int *and0, *and1, *and2;

	for (int row = 0; row < 384; row++)
	{
		matrix0 = _mm_setr_epi32(inv[row][0], inv[row][1], inv[row][2], inv[row][3]);
		matrix1 = _mm_setr_epi32(inv[row][4], inv[row][5], inv[row][6], inv[row][7]);
		matrix2 = _mm_setr_epi32(inv[row][8], inv[row][9], inv[row][10], inv[row][11]);

		result0 = _mm_and_si128(state0, matrix0);
		result1 = _mm_and_si128(state1, matrix1);
		result2 = _mm_and_si128(state2, matrix2);

		ull result00 = _mm_extract_epi64(result0, 0);
		ull result01 = _mm_extract_epi64(result0, 1);
		ull result10 = _mm_extract_epi64(result1, 0);
		ull result11 = _mm_extract_epi64(result1, 1);
		ull result20 = _mm_extract_epi64(result2, 0);
		ull result21 = _mm_extract_epi64(result2, 1);

		int sum = 0;
		sum += __builtin_popcountll(result00) + __builtin_popcountll(result01);
		sum += __builtin_popcountll(result10) + __builtin_popcountll(result11);
		sum += __builtin_popcountll(result20) + __builtin_popcountll(result21);

		int bit = sum % 2;
		int word = row / 32;
		int pos = row % 32;

		output[word] ^= bit << pos;
	}
	A0[0] = output[0];
	A0[1] = output[1];
	A0[2] = output[2];
	A0[3] = output[3];
	A1[0] = output[4];
	A1[1] = output[5];
	A1[2] = output[6];
	A1[3] = output[7];
	A2[0] = output[8];
	A2[1] = output[9];
	A2[2] = output[10];
	A2[3] = output[11];
}


void theta_inv_simd(ui *A0, ui *A1, ui *A2)
{
	gf2mult_simd(A0, A1, A2, theta_inv_matrix);
}

//_____________________________________________________________//

////////////////////     WRhoRhoInverse     ////////////////////

void rho_west(ui *A0, ui *A1, ui *A2)
{
	rotate_sheets(A1, 1);
	rotate_lanes(A2, 11);
}

void rho_west_inv(ui *A0, ui *A1, ui *A2)
{
	rotate_lanes(A2, -11);
	rotate_sheets(A1, -1);
}

//____________________________________________________________//

////////////////////     IotaIotaInverse     ////////////////////

void iota(ui *A0, ui *A1, ui *A2, int round)
{
	ui rcons[12] = {0x58, 0x38, 0x3C0, 0xD0, 0x120, 0x14, 0x60, 0x2C, 0x380, 0xF0, 0x1A0, 0x12};
	*(A0) = *(A0) ^ *(rcons + 11 + round);
}

void iota_inv(ui *A0, ui *A1, ui *A2, int round)
{
	ui rcons[12] =  {0x58, 0x38, 0x3C0, 0xD0, 0x120, 0x14, 0x60, 0x2C, 0x380, 0xF0, 0x1A0, 0x12};
	*(A0) = *(A0) ^ *(rcons + 11 + round);
}

//_____________________________________________________________//

////////////////////     ChiChiInverse     ////////////////////

void chi(ui *A0, ui *A1, ui *A2)
{
	ui B0[PLANE_SIZE], B1[PLANE_SIZE], B2[PLANE_SIZE];
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(B0 + i) = ~(*(A1 + i)) & *(A2 + i);
		*(B1 + i) = ~(*(A2 + i)) & *(A0 + i);
		*(B2 + i) = ~(*(A0 + i)) & *(A1 + i);
	}
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(A0 + i) = *(A0 + i) ^ *(B0 + i);
		*(A1 + i) = *(A1 + i) ^ *(B1 + i);
		*(A2 + i) = *(A2 + i) ^ *(B2 + i);
	}
}

void chi_inv(ui *A0, ui *A1, ui *A2)
{
	ui B0[PLANE_SIZE], B1[PLANE_SIZE], B2[PLANE_SIZE];
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(B0 + i) = ~(*(A1 + i)) & *(A2 + i);
		*(B1 + i) = ~(*(A2 + i)) & *(A0 + i);
		*(B2 + i) = ~(*(A0 + i)) & *(A1 + i);
	}
	for (int i = 0; i < PLANE_SIZE; i++)
	{
		*(A0 + i) = *(A0 + i) ^ *(B0 + i);
		*(A1 + i) = *(A1 + i) ^ *(B1 + i);
		*(A2 + i) = *(A2 + i) ^ *(B2 + i);
	}
}

//___________________________________________________________//

////////////////////     ChiChiInverseSIMD     ////////////////////

void chi_simd(ui *A0, ui *A1, ui *A2)
{
	__m128i state0 = _mm_setr_epi32(A0[0], A0[1], A0[2], A0[3]);
	__m128i state1 = _mm_setr_epi32(A1[0], A1[1], A1[2], A1[3]);
	__m128i state2 = _mm_setr_epi32(A2[0], A2[1], A2[2], A2[3]);

	__m128i andnot0 = _mm_andnot_si128(state1, state2);
	__m128i andnot1 = _mm_andnot_si128(state2, state0);
	__m128i andnot2 = _mm_andnot_si128(state0, state1);

	__m128i result0 = _mm_xor_si128(state0, andnot0);
	__m128i result1 = _mm_xor_si128(state1, andnot1);
	__m128i result2 = _mm_xor_si128(state2, andnot2);

	A0[0] = _mm_extract_epi32(result0, 0);
	A0[1] = _mm_extract_epi32(result0, 1);
	A0[2] = _mm_extract_epi32(result0, 2);
	A0[3] = _mm_extract_epi32(result0, 3);
	A1[0] = _mm_extract_epi32(result1, 0);
	A1[1] = _mm_extract_epi32(result1, 1);
	A1[2] = _mm_extract_epi32(result1, 2);
	A1[3] = _mm_extract_epi32(result1, 3);
	A2[0] = _mm_extract_epi32(result2, 0);
	A2[1] = _mm_extract_epi32(result2, 1);
	A2[2] = _mm_extract_epi32(result2, 2);
	A2[3] = _mm_extract_epi32(result2, 3);
}

void chi_inv_simd(ui *A0, ui *A1, ui *A2)
{
	__m128i state0 = _mm_setr_epi32(A0[0], A0[1], A0[2], A0[3]);
	__m128i state1 = _mm_setr_epi32(A1[0], A1[1], A1[2], A1[3]);
	__m128i state2 = _mm_setr_epi32(A2[0], A2[1], A2[2], A2[3]);

	__m128i andnot0 = _mm_andnot_si128(state1, state2);
	__m128i andnot1 = _mm_andnot_si128(state2, state0);
	__m128i andnot2 = _mm_andnot_si128(state0, state1);

	__m128i result0 = _mm_xor_si128(state0, andnot0);
	__m128i result1 = _mm_xor_si128(state1, andnot1);
	__m128i result2 = _mm_xor_si128(state2, andnot2);

	A0[0] = _mm_extract_epi32(result0, 0);
	A0[1] = _mm_extract_epi32(result0, 1);
	A0[2] = _mm_extract_epi32(result0, 2);
	A0[3] = _mm_extract_epi32(result0, 3);
	A1[0] = _mm_extract_epi32(result1, 0);
	A1[1] = _mm_extract_epi32(result1, 1);
	A1[2] = _mm_extract_epi32(result1, 2);
	A1[3] = _mm_extract_epi32(result1, 3);
	A2[0] = _mm_extract_epi32(result2, 0);
	A2[1] = _mm_extract_epi32(result2, 1);
	A2[2] = _mm_extract_epi32(result2, 2);
	A2[3] = _mm_extract_epi32(result2, 3);
}

//_______________________________________________________________//

////////////////////     ERhoRhoInverse     ////////////////////

void rho_east(ui *A0, ui *A1, ui *A2)
{
	rotate_lanes(A1, 1);
	rotate_sheets(A2, 2);
	rotate_lanes(A2, 8);
}

void rho_east_inv(ui *A0, ui *A1, ui *A2)
{
	rotate_lanes(A2, -8);
	rotate_sheets(A2, -2);
	rotate_lanes(A1, -1);
}

//____________________________________________________________//

////////////////////     XoodooNormal     ////////////////////

void xoodoo(ui *A0, ui *A1, ui *A2, int round_initial, int round_final)
{
	for (int round = round_initial; round <= round_final; round++)
	{
		theta(A0, A1, A2);
		//printf("Theta\n");
		//display_state(A0, A1, A2);

		rho_west(A0, A1, A2);
		//printf("Rho West\n");
		//display_state(A0, A1, A2);

		iota(A0, A1, A2, round);
		//printf("Iota\n");
		//display_state(A0, A1, A2);

		chi(A0, A1, A2);
		//printf("Chi\n");
		//display_state(A0, A1, A2);

		rho_east(A0, A1, A2);
		//printf("Rho East\n");
		//display_state(A0, A1, A2);
	}
}

//__________________________________________________________//

////////////////////     XoodooSIMD     ////////////////////

void xoodoo_simd(ui *A0, ui *A1, ui *A2, int round_initial, int round_final)
{
	for (int round = round_initial; round <= round_final; round++)
	{
		theta(A0, A1, A2);

		rho_west(A0, A1, A2);

		iota(A0, A1, A2, round);

		chi_simd(A0, A1, A2);

		rho_east(A0, A1, A2);
	}
}

//________________________________________________________//

////////////////////     XoodooInverseNormal     ////////////////////

void xoodoo_inv(ui *A0, ui *A1, ui *A2, int round_initial, int round_final)
{
	for (int round = round_initial; round >= round_final; round--)
	{
		rho_east_inv(A0, A1, A2);
		//printf("Inverse Rho East\n");
		//display_state(A0, A1, A2);

		chi_inv(A0, A1, A2);
		//printf("Inverse Chi\n");
		//display_state(A0, A1, A2);

		iota_inv(A0, A1, A2, round);
		//printf("Inverse Iota\n");
		//display_state(A0, A1, A2);

		rho_west_inv(A0, A1, A2);
		//printf("Inverse Rho West\n");
		//display_state(A0, A1, A2);

		theta_inv(A0, A1, A2);
		//printf("Inverse Theta\n");
		//display_state(A0, A1, A2);
	}
}

//_________________________________________________________________//

////////////////////     XoodooInverseSIMD     ////////////////////

void xoodoo_inv_simd(ui *A0, ui *A1, ui *A2, int round_initial, int round_final)
{
	for (int round = round_initial; round >= round_final; round--)
	{
		rho_east_inv(A0, A1, A2);

		chi_inv_simd(A0, A1, A2);

		iota_inv(A0, A1, A2, round);

		rho_west_inv(A0, A1, A2);

		theta_inv_simd(A0, A1, A2);
	}
}

//_______________________________________________________________//

void xoodoo_on_lane(Lane *state, int rounds){
    int ri = 1-rounds;
    int rf = 0;
    Lane A0[4];
    Lane A1[4];
    Lane A2[4];
    
    memcpy(A0, state, 4*sizeof(Lane));
    memcpy(A1, state+4, 4*sizeof(Lane));
    memcpy(A2, state+8, 4*sizeof(Lane));
    
    /* display_states(A0, A1, A2); */
    /* display_lanes("checking conversion ", state); */
    xoodoo_simd(A0, A1, A2, ri, rf);
    
    memcpy(state,   A0, 4*sizeof(Lane));
    memcpy(state+4, A1, 4*sizeof(Lane));
    memcpy(state+8, A2, 4*sizeof(Lane));

}
void xoodoo_inv_on_lane(Lane *state, int rounds){
    int ri = 0 - rounds;
    int rf = 1 - (2*rounds);
    /* int ri = 0; */
    /* int rf = 1 - rounds; */

    Lane A0[4];
    Lane A1[4];
    Lane A2[4];

    memcpy(A0, state, 4*sizeof(Lane));
    memcpy(A1, state+4, 4*sizeof(Lane));
    memcpy(A2, state+8, 4*sizeof(Lane));

    /* display_states(A0, A1, A2); */
    /* display_lanes("checking conversion ", state); */
    xoodoo_inv_simd(A0, A1, A2, ri, rf);

    memcpy(state,   A0, 4*sizeof(Lane));
    memcpy(state+4, A1, 4*sizeof(Lane));
    memcpy(state+8, A2, 4*sizeof(Lane));

}



//int main(void)
//{
//	ui a0[4] = {0x13,0x0,0x0,0x10000080};
//	ui a1[4] = {0x24,0x0,0x0,0x20000000};
//	ui a2[4] = {0x0,0x8000,0x0,0x0};
//
//	int max = 1 << 16;
//	max = 100000;
//
//	clock_t start = clock();
//
//	for (int i = 0; i < max; i++)
//	{
////		display_state(a0,a1,a2);
//		xoodoo_simd(a0, a1, a2, -4, 0);
////		display_state(a0,a1,a2);
//		xoodoo_inv_simd(a0, a1, a2, 0, -4);
////		display_state(a0,a1,a2);
//	}
//
//	clock_t end = clock();
//	double time = (double)(end - start) / CLOCKS_PER_SEC;
//	printf("Time  = %f\n", time);
//}
