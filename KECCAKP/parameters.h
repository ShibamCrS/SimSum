typedef uint64_t U64;
typedef uint32_t U32;
typedef uint16_t U16;
typedef uint8_t U8;

#define maxNrRounds 24
#define nrLanes 25
#define index(x, y) (((x)%5)+5*((y)%5))
#define ROL64(a, offset) ((offset != 0) ? ((((U64) a) << offset) ^ (((U64) a) >> (64-offset))) : a)


//#define KeccakReference

#ifdef KeccakReference
static U64 RC[maxNrRounds];
static U32  RO[nrLanes];

static int LFSR86540(U8 *LFSR){
    int result = ((*LFSR) & 0X01) != 0 ;
    if ( ((*LFSR) & 0x80) != 0 ){
        /* Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1 */
        (*LFSR) = ( (*LFSR) << 1 ) ^ 0x71;
    }
    else{
        (*LFSR) <<= 1;
    }
    return result;
}
void KeccakP1600_rc(void){
    U8 LFSRstate = 0x01;
    U32 i, j, bitPos;
    
    for(i=0; i<maxNrRounds; i++){
        RC[i] = 0;
        for( j=0; j<7; j++){
            bitPos = (1 << j) - 1;  // 0, 1, 3, 7, 15, 31, 63
            if (LFSR86540(&LFSRstate)){
                RC[i] ^= (U64) (1UL << bitPos);
            }
        }
    }
}

void KeccakP1600_ro(void){
    U32 x, y, t, newX, newY;
    RO[index(0, 0)] = 0;
    x = 1;
    y = 0;
    for(t=0; t<24; t++){
        RO[index(x, y)] = ( (t+1)*(t+2)/2 ) % 64;
        newX = (0*x + 1*y) % 5;
        newY = (2*x + 3*y) % 5;
        x = newX;
        y = newY;
    }
}

void KeccakP1600_StaticInitialize(void){
    KeccakP1600_rc();
    KeccakP1600_ro();
}

#else
    #if WIDTH == 1600
    static const U64 RC[maxNrRounds] =
    {
        0x0000000000000001,
        0x0000000000008082,
        0x800000000000808a,
        0x8000000080008000,
        0x000000000000808b,
        0x0000000080000001,
        0x8000000080008081,
        0x8000000000008009,
        0x000000000000008a,
        0x0000000000000088,
        0x0000000080008009,
        0x000000008000000a,
        0x000000008000808b,
        0x800000000000008b,
        0x8000000000008089,
        0x8000000000008003,
        0x8000000000008002,
        0x8000000000000080,
        0x000000000000800a,
        0x800000008000000a,
        0x8000000080008081,
        0x8000000000008080,
        0x0000000080000001,
        0x8000000080008008,
    };

    static const unsigned int RO[nrLanes] =
    {
         0,  1, 62, 28, 27, 36, 44,  6, 55, 20,  3, 10, 43, 25, 39, 41, 45, 15, 21,  8, 18,  2, 61, 56, 14
    };
    #elif WIDTH == 800
    static const Lane RC[maxNrRounds] =
    {
        0x00000001,
        0x00008082,
        0x0000808a,
        0x80008000,
        0x0000808b,
        0x80000001,
        0x80008081,
        0x00008009,
        0x0000008a,
        0x00000088,
        0x80008009,
        0x8000000a,
        0x8000808b,
        0x0000008b,
        0x00008089,
        0x00008003,
        0x00008002,
        0x00000080,
        0x0000800a,
        0x8000000a,
        0x80008081,
        0x00008080,
    };

    static const unsigned int RO[nrLanes] =
    {
         0,  1, 30, 28, 27,  4, 12,  6, 23, 20,  3, 10, 11, 25,  7,  9, 13, 15, 21,  8, 18,  2, 29, 24, 14
    };
    #elif WIDTH == 400
    static const Lane RC[maxNrRounds] =
    {
        0x0001,
        0x8082,
        0x808a,
        0x8000,
        0x808b,
        0x0001,
        0x8081,
        0x8009,
        0x008a,
        0x0088,
        0x8009,
        0x000a,
        0x808b,
        0x008b,
        0x8089,
        0x8003,
        0x8002,
        0x0080,
        0x800a,
        0x000a,
    };

    static const unsigned int RO[nrLanes] =
    {
         0,  1, 14, 12, 11,  4, 12,  6,  7,  4,  3, 10, 11,  9,  7,  9, 13, 15,  5,  8,  2,  2, 13,  8, 14
    };
    #elif WIDTH == 200
    static const Lane RC[maxNrRounds] =
    {
        0x01,
        0x82,
        0x8a,
        0x00,
        0x8b,
        0x01,
        0x81,
        0x09,
        0x8a,
        0x88,
        0x09,
        0x0a,
        0x8b,
        0x8b,
        0x89,
        0x03,
        0x02,
        0x80,
    };

    static const unsigned int RO[nrLanes] =
    {
         0, 1, 6, 4, 3, 4, 4, 6, 7, 4, 3, 2, 3, 1, 7, 1, 5, 7, 5, 0, 2, 2, 5, 0, 6
    };
    #endif
#endif
