#include <stdint.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define WIDTH 1600
#define ROUNDS 4
#define DIM 8
#define NROF_THREADS 32
#include "../KECCAKP/utility.h"
#include "../KECCAKP/KeccakP.h"
/*
SHA-512 => RATE 9
SHA-384 => RATE 13
SHA-256 => RATE 17
SHA-224 => RATE 18
permutation => RATE 25
*/
#define RATE 9 //Number of lane that can be accessed

struct SumArgs{
    Lane *state;
    uint64_t start;
    uint64_t end;
    Lane result[25];
};
typedef struct SumArgs SumArgs;

Lane generate_random_symmetric_lane(){
    uint64_t temp = xx_next(SEED);
    temp = temp & HALFMASK;
    Lane lane =  (temp << HALFLANE) | temp;
    return lane;
}

void generateRandomSymmetricState(Lane *state){
    memset(state, 0x00, 25*sizeof(Lane));
    int nrof_accessed_lanes = RATE;
    /* printf("nrof_accessed_lanes %d RATE %d SIZELANE %d \n", nrof_accessed_lanes, RATE, LANE); */
    for(int i=0; i<nrof_accessed_lanes; i++){
        state[i] = generate_random_symmetric_lane();
    }  
    /* state[index(0,4)] = state[index(0,3)]^state[index(0,2)]; */ 
     
}

uint64_t get_property(Lane *sum){
    int asym_bits = 0UL;
    for(int i=0; i<25; i++){
        uint32_t up = (sum[i] >> HALFLANE) & HALFMASK;
        uint32_t low = sum[i] & HALFMASK;
        uint32_t temp = up^low;
        asym_bits = asym_bits + __builtin_popcount(temp);
    }
    if (asym_bits == 0){
        printf("Yah! Symmetric Property\n");
    }
    else{
        printf("# of asymmetric bits = %d\n",asym_bits);
    }
    
    int nzero_counter = 0UL;
    for(int i=0; i<25; i++){
        nzero_counter = nzero_counter + __builtin_popcountll(sum[i]);
    }
    if (nzero_counter == 0){
        printf("Yah! Zerosum Property\n");
    }
    else{
        printf("# of nonzero bits = %d\n",nzero_counter);
    }
}
void check_equations(Lane *state){
    uint64_t eq0 = state[index(0,0)];
    eq0 = eq0 ^ state[index(0,1)];
    eq0 = eq0 ^ state[index(0,2)];
    eq0 = eq0 ^ state[index(0,3)];
    eq0 = eq0 ^ state[index(0,4)];
    printf("Eqs = %ld \n",eq0);
    get_property(state);
}
void SUM(Lane *sum, Lane *state){
    for(int i=0; i<25; i++){
        sum[i] = sum[i] ^ state[i];
    }
}
void constraint_1linear_sym(Lane *state, uint64_t mask){
    //create symmetric mask
    uint64_t mask_up = (mask >> (DIM/2));
    uint64_t mask_low = mask ^ (mask_up << (DIM/2));
    printf("Mask %016lX Mask_up %016lX Mask_low %016lX\n",mask, mask_up, mask_low);
     
    state[index(0,0)] = state[index(0,0)] ^ mask_low ^ (mask_up << HALFLANE);
    state[index(0,1)] = state[index(0,0)];
}
void* sum_over_range(void *args){
    SumArgs *sumargs = (SumArgs*) args;
    Lane temp[25];
    for(uint64_t i=sumargs->start; i<sumargs->end; i++){
        /* printf("**************************Input no %ld*************************\n",i); */
        memcpy(temp, sumargs->state, 25*sizeof(Lane));
        uint64_t mask = i ;
        constraint_1linear_sym(temp, mask);
        KeccakP(temp, ROUNDS);
        SUM(sumargs->result, temp);
    }
    pthread_exit(NULL);
}
void compute_sum_threaded(Lane *state, Lane *sum){
    pthread_t thread_ids[NROF_THREADS]; 
    SumArgs thread_args[NROF_THREADS];
    
    uint64_t cube_size = (1ULL << DIM);
    for(int i=0; i < NROF_THREADS; i++){
        memset(thread_args[i].result, 0x00, 25*sizeof(Lane));
        thread_args[i].state = state;
        thread_args[i].start = i * (cube_size / NROF_THREADS);
        thread_args[i].end = (i+1) * (cube_size / NROF_THREADS);
    }

    for(int i=0; i < NROF_THREADS; i++){
        pthread_create(thread_ids + i, NULL, sum_over_range, (void*) (thread_args + i));
    }
    for(int i=0; i < NROF_THREADS; i++){
        pthread_join(thread_ids[i], NULL);
        SUM(sum, thread_args[i].result);
    }
}

void compute_sum(Lane *state, Lane *sum){
    uint64_t cube = (1UL << DIM);
    Lane temp[25];
    Lane input_sum[25];
    memset(input_sum, 0UL, 25*sizeof(Lane));

    for(uint64_t i=0UL; i<cube; i++){
        /* printf("**********************Input no %ld*************************\n",i); */
        memcpy(temp, state, 25*sizeof(Lane));
        uint64_t mask = i ;
        constraint_1linear_sym(temp, mask);
        /* display_lanes("INPUT", temp, WIDTH); */
        /* check_equations(temp); */
        KeccakP(temp, ROUNDS);
        SUM(sum, temp);
        /* if (i == (cube-3)){ */
        /*     display_lanes("SUM", sum, WIDTH); */
        /*     get_property(sum); */
        /* } */
    }
    printf("*******************************************************************************\n");
}
void test_sysmsum_w_linear(){
    Lane state[25];
    generateRandomSymmetricState(state);
    display_lanes("Starting: ", state, WIDTH);

    Lane sum[25];
    memset(sum, 0x00, 25*sizeof(Lane));
    display_lanes("SUM", sum, WIDTH);
    compute_sum(state, sum);
    /* compute_sum_threaded(state, sum); */
    display_lanes("SUM", sum, WIDTH);
    get_property(sum);
}
int main(){
    xx_initialize(SEED);
    printf("SEED: ");
    printreg(SEED, 32);

    FILE *fp = stdout;
    set_file(fp);   
    test_sysmsum_w_linear();
}
