#include <stdint.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define WIDTH 1600
#define ROUNDS 6
#include "../../KECCAKP/utility.h"
#include "../../KECCAKP/KeccakP.h"
#include "set_equations.h"  //generate constraints states that maintain 2 round linearization conditions

void generateLinearConstraintState(Lane *state){
    memset(state, 0x00, 25*sizeof(Lane));
    int nrof_accessed_lanes = 25;
    Lane ZERO = 0UL;
    Lane ONE = 0xFFFFFFFFFFFFFFFFUL;
    int one_lane[3] = {index(1,1),index(1,2),index(1,3)};
    for(int i=0; i<nrof_accessed_lanes; i++){
        state[i] = ZERO;
    }
    for(int i=0; i<3;i++){
        state[one_lane[i]] = ONE;
    }
    //generate random symmetric 54 bit value
    uint64_t rr = rand();
    uint64_t mask27 = rr & 0x7FFFFFFUL;
    printf("mask27 = 0x%016lX \n",mask27);
    state[index(0,0)] ^= (mask27 <<5) | (mask27 << 37); 
}
void SUM(Lane *sum, Lane *state){
    for(int i=0; i<25; i++){
        sum[i] = sum[i] ^ state[i];
    }
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
/*
A[0, 0] ⊕ A[0, 1] ⊕ A[0, 2] ⊕ A[0, 3] = 0xff· · · f
A[2, 0] ⊕ A[2, 1] ⊕ A[2, 2] ⊕ A[2, 3] = 0xff· · · f
A[2, 0]≪ 62 = A[0, 0] ⊕ A[2, 2]≪ 43
A[2, 1]≪ 6 = A[0, 1]≪ 36
A[2, 2]≪ 43 = A[0, 2]≪ 3
*/

void check_equations(Lane *state){
    uint64_t eq0 = state[index(0,0)] ^ state[index(0,1)] ^ state[index(0,2)];
    uint64_t eq1 = state[index(2,0)] ^ state[index(2,1)] ^ state[index(2,2)];
    uint64_t eq2 = state[index(0,0)] ^ ROL64(state[index(2,0)],62) ^ ROL64(state[index(2,2)],43);
    uint64_t eq3 = ROL64(state[index(2,1)],6) ^ ROL64(state[index(0,1)],36);
    uint64_t eq4 = ROL64(state[index(2,2)],43) ^ ROL64(state[index(0,2)],3);
    printf("Eqs = %ld %ld %ld %ld %ld\n",eq0,eq1,eq2,eq3,eq4);
}

/* FOR 7 rounds we implement the same functuin in threads */
struct SumArgs{
    Lane state[25];
    uint64_t start;
    uint64_t end;
    Lane result[25];
};
typedef struct SumArgs SumArgs;
void* sum_over_range(void *args){
    /* int hdim   = 8; */
    /* uint64_t CMASK = 0xFFUL; */
    
    int hdim   = 16;
    uint64_t CMASK = 0xFFFFUL;

    SumArgs *sumargs = (SumArgs*) args;
    Lane temp[25];
    for(uint64_t i=sumargs->start; i<sumargs->end; i++){
        memcpy(temp, sumargs->state, 25*sizeof(Lane));
        uint64_t mask1 = (i & CMASK);
        uint64_t mask2 = ((i>>hdim) & CMASK);
        temp[index(0,0)] =temp[index(0,0)] ^ (mask1 << 8) ^ (mask2 << 40);
        set_sym_eqns(temp);
        KeccakP(temp, ROUNDS);
        SUM(sumargs->result, temp);
    }
    pthread_exit(NULL);
}
void compute_sum_threaded(Lane *state, Lane *sum){
    int THREAD = 32;
    int dim    = 32;
    pthread_t thread_ids[THREAD]; 
    SumArgs thread_args[THREAD];
    
    uint64_t size = (1ULL << dim);

    uint64_t data_each_thread = size / THREAD;
    uint64_t data_last_thread = size % THREAD;
    printf("%f \n", log(data_each_thread)/log(2));
    printf("%f \n", log(data_last_thread)/log(2));

    for(int i=0; i < THREAD; i++){
        memset(thread_args[i].result, 0x00, 25*sizeof(Lane));
        memcpy(thread_args[i].state,state, 25*sizeof(Lane));

        if( i == (THREAD-1) ){
            thread_args[i].start = i * data_each_thread;
            thread_args[i].end = (i+1) * data_each_thread + data_last_thread;
        }
        else{
            thread_args[i].start = i * data_each_thread;
            thread_args[i].end = (i+1) * data_each_thread;
        }
    }

    for(int i=0; i < THREAD; i++){
        pthread_create(thread_ids + i, NULL, sum_over_range, (void*) (thread_args + i));
    }
    for(int i=0; i < THREAD; i++){
        pthread_join(thread_ids[i], NULL);
        SUM(sum, thread_args[i].result);
    }
    printf("ROUNDS = %d DIM = %d \n", ROUNDS, dim);
    printf("*******************************************************************************\n");
}
void compute_sum(Lane *state, Lane *sum){
    int lane_no = controll_lane[5];
    int dim, hdim;
    uint64_t CMASK;
    if (ROUNDS == 5){
        dim = 8;
        hdim = 4;
        CMASK = 0x0FUL;
    }
    else if (ROUNDS == 6){
        dim = 16;
        hdim = 8;
        CMASK = 0xFF;
    }
    else{
        printf("Invalid Rounds\n");
        exit(1);
    }
    uint64_t cube = (1UL << dim);
    Lane temp[25];
    Lane input_sum[25];
    memset(input_sum, 0UL, 25*sizeof(Lane));


    
    for(uint64_t i=0UL; i<cube; i++){
        memcpy(temp, state, 25*sizeof(Lane));
        uint64_t mask1 = (i & CMASK);
        uint64_t mask2 = ((i>>hdim) & CMASK);
        temp[index(0,0)] =temp[index(0,0)] ^ (mask1 << 8) ^ (mask2 << 40);
        set_sym_eqns(temp);
        /* display_lanes("state", temp, WIDTH); */
        /* check_equations(temp); */
        KeccakP(temp, ROUNDS);
        SUM(sum, temp);
    }
    printf("ROUNDS = %d DIM = %d \n", ROUNDS, dim);
    printf("*******************************************************************************\n");
}
void test_sysmsum_w_linear(){
    Lane state[25];
    generateLinearConstraintState(state);
    display_lanes("Valid state for two rounds linear(independent bits taking random symmetric values)", state, WIDTH);

    Lane sum[25];
    memset(sum, 0x00, WIDTH/8);
    display_lanes("SUM", sum, WIDTH);
    if (ROUNDS == 7){
        compute_sum_threaded(state, sum);
    }
    else{
        compute_sum(state, sum);
    }
    display_lanes("SUM", sum, WIDTH);
    get_property(sum);
}
int main(){
    srand(time(NULL));

    FILE *fp = stdout;
    set_file(fp);
    test_sysmsum_w_linear();
}
