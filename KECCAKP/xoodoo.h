#include "matrix.h"
typedef uint32_t U32;
#define ui unsigned int
#define STATE_SIZE 12
#define PLANE_SIZE 4
#define SHEET_SIZE 3

typedef uint32_t Lane;
#define LANE 32
#define HALFLANE 16
#define MASK 0xFFFFFFFF
#define HALFMASK 0xFFFF

/* extern ui theta_inv_matrix[384][12]; */

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

void display_states(ui *A0, ui *A1, ui *A2)
{
	printf("%x %x %x %x\n", A2[0], A2[1], A2[2], A2[3]);
	printf("%x %x %x %x\n", A1[0], A1[1], A1[2], A1[3]);
	printf("%x %x %x %x\n", A0[0], A0[1], A0[2], A0[3]);
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

int countmult(ui a, ui b)
{
        int cnt = 0;
        ui c = reverse(a);
        ui d = c & b;
        cnt = __builtin_popcount(d);
        /* for (int k = 0; k < 32; k++) */
        /* { */
        /*         if (d >> k & 1) */
        /*         { */
        /*                 cnt++; */
        /*         } */
        /* } */
        return cnt;
}

void gf2mult(ui *A0, ui *A1, ui *A2, ui inv[384][12])
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

void theta_inv(ui *A0, ui *A1, ui *A2)
{
	gf2mult(A0, A1, A2, theta_inv_matrix);
}

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
    xoodoo(A0, A1, A2, ri, rf);
    
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
    xoodoo_inv(A0, A1, A2, ri, rf);

    memcpy(state,   A0, 4*sizeof(Lane));
    memcpy(state+4, A1, 4*sizeof(Lane));
    memcpy(state+8, A2, 4*sizeof(Lane));

}
