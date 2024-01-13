#if WIDTH == 1600
    typedef uint64_t Lane;
    #define LANE 64
    #define HALFLANE 32
    #define MASK 0xFFFFFFFFFFFFFFFF
    #define HALFMASK 0xFFFFFFFF
    #define ROL(a, offset) ((offset != 0) ? ((((Lane)a) << offset) ^ (((Lane) a) >> (sizeof(Lane)*8-offset))) : a)
#elif WIDTH == 800
    typedef uint32_t Lane;
    #define LANE 32
    #define HALFLANE 16
    #define MASK 0xFFFFFFFF
    #define HALFMASK 0xFFFF
    #define ROL(a, offset) ((offset != 0) ? ((((Lane)a) << offset) ^ (((Lane)a) >> (sizeof(Lane)*8-offset))) : a)
#elif WIDTH == 400
    typedef uint16_t Lane;
    #define LANE 16
    #define HALFLANE 8
    #define MASK 0xFFFF
    #define HALFMASK 0xFF
    #define ROL(a, offset) ((offset != 0) ? ((((Lane)a) << offset) ^ (((Lane)a) >> (sizeof(Lane)*8-offset))) : a)
#elif WIDTH == 200
    typedef uint8_t Lane;
    #define LANE 8
    #define HALFLANE 4
    #define MASK 0xFF
    #define HALFMASK 0xF
    #define ROL(a, offset) ((offset != 0) ? ((((Lane)a) << offset) ^ (((Lane)a) >> (sizeof(Lane)*8-offset))) : a)
#endif

#include "parameters.h"
#include "display.h"
//add parity
static void theta(Lane *A){

    U32 x, y;
    Lane C[5];
    Lane D[5];
    
    //for each row compute column parity
    //Note that we are working on the full lane
    for(x=0; x<5; x++){
        C[x] = 0;
        for(y=0; y<5; y++){
            C[x] ^= A[index(x, y)];

        }
    }
    //C[(x+1) mod 5, (z – 1) mod w] + C[(x-1) mod 5, z]
    for(x=0; x<5; x++){
        D[x] = ROL(C[(x+1)%5], 1) ^ C[(x+4)%5];
    }
    
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
            A[index(x, y)] ^= D[x];
        }
    }
}

//shift inside lane
static void rho(Lane *A){
    U32 x, y;
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
           A[index(x, y)] = ROL(A[index(x, y)], RO[index(x, y)]);
        }
    }
    
}

//change positions of the lanes
static void pi(Lane *A){
    U32 x, y;
    Lane A_prime[25];
    
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
           A_prime[index(x, y)] = A[index(x, y)];
        }
    }
    
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
            A[index(0*x + 1*y, 2*x + 3*y)] = A_prime[index(x, y)];
        }
    }
}

//sbox
static void chi(Lane *A){
    U32 x, y;
    Lane S[5];
    for(y=0; y<5; y++){
        //Apply over each row
        for(x=0; x<5; x++){
            S[x] = A[index(x, y)]^((~A[index(x+1, y)])&A[index(x+2, y)]);
        }
        for(x=0; x<5; x++){
            A[index(x, y)] = S[x];
        }
    }
}

//add round constants
static void iota(Lane *A, U32 indexRound)
{
    //Only on the first lane
    A[index(0, 0)] ^= RC[indexRound];
}

void KeccakP(void *State, U32 rounds){
    U32 index;
    Lane *state = (Lane *)State;
    
    #ifdef DISPLAY
        display_lanes("Starting: ", state, WIDTH);
    #endif
    
    //Main Permutation
    //for(index=(maxNrRounds - rounds); index<maxNrRounds; index++){
    for(index=0; index<rounds; index++){
    #ifdef DISPLAY
        display_round_number(index);
    #endif
        theta(state);
    #ifdef DISPLAY
        display_lanes("After theta: ", state, WIDTH);
    #endif
        rho  (state);
    #ifdef DISPLAY
        display_lanes("After rho: ", state, WIDTH);
    #endif
        pi   (state);
    #ifdef DISPLAY
        display_lanes("After pi: ", state, WIDTH);
    #endif
        chi  (state);
    #ifdef DISPLAY
        display_lanes("After chi: ", state, WIDTH);
    #endif
        /* if(index !=0 ) */
        iota (state, index);
    #ifdef DISPLAY
        display_lanes("After iota: ", state, WIDTH);
    #endif
    }
}
