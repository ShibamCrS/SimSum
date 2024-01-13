//Note that U64 is Lane

//Set all bytes to zero
void state_initialize(void *state){
    memset(state, 0, 1600/8);
}

//add value to specific bytes
void xor_byte(void *state, U8 byte, U32 offset){
    ( (U8 *)state )[offset] ^= byte;
}

//add "length" number of bytes from byte number "offset"
void xor_bytes(void *state, U8 *bytes, U32 offset, U32 length){
    U32 i;
    for (i=0; i<length; i++){
        ( (U8 *)state )[offset] ^= bytes[i];
    }
}

//replace "length" number of bytes from byte number "offset"
void overwrite_bytes(void *state, U8 *bytes, U32 offset, U32 length){
    U32 i;
    memcpy( (U8*)state+offset, bytes, length );
}

void overwrite_with_zeros(void *state, U32 length){
    memset(state, 0, length);
}

//add parity
static void theta(U64 *A){
    U32 x, y;
    U64 C[5];
    U64 D[5];
    
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
        D[x] = ROL64(C[(x+1)%5], 1) ^ C[(x+4)%5];
    }
    
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
            A[index(x, y)] ^= D[x];
        }
    }
}

//shift inside lane
static void rho(U64 *A){
    U32 x, y;
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
           A[index(x, y)] = ROL64(A[index(x, y)], RO[index(x, y)]);
        }
    }
    
}

//change positions of the lanes
static void pi(U64 *A){
    U32 x, y;
    U64 A_prime[25];
    
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
static void chi(U64 *A){
    U32 x, y;
    U64 S[5];
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
static void iota(U64 *A, U32 indexRound)
{
    //Only on the first lane
    A[index(0, 0)] ^= RC[indexRound];
}

void KeccakP1600(void *State, U32 rounds){
    U32 index;
    U64 *state = (U64*)State;
    
    #ifdef DISPLAY
        display_lanes("Starting: ", state, 1600);
    #endif
    
    //Main Permutation
    //for(index=(maxNrRounds - rounds); index<maxNrRounds; index++){
    for(index=0; index<rounds; index++){
    #ifdef DISPLAY
        display_round_number(index);
    #endif
        theta(state);
    #ifdef DISPLAY
        display_lanes("After theta: ", state, 1600);
    #endif
        rho  (state);
    #ifdef DISPLAY
        display_lanes("After rho: ", state, 1600);
    #endif
        pi   (state);
    #ifdef DISPLAY
        display_lanes("After pi: ", state, 1600);
    #endif
        chi  (state);
    #ifdef DISPLAY
        display_lanes("After chi: ", state, 1600);
    #endif
        iota (state, index);
    #ifdef DISPLAY
        display_lanes("After iota: ", state, 1600);
    #endif
    }
}
