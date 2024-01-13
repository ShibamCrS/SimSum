#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define ROUNDS 5
#define DIM 16
#define NROF_THREADS 32
#include "../KECCAKP/xoodoo_simd.h"
#include "../KECCAKP/utility.h"

struct SumArgs{
    Lane *state;
    uint64_t start;
    uint64_t end;
    Lane result[12];
    Lane result_inverse[12];
};
typedef struct SumArgs SumArgs;

void get_property(Lane *sum){
    int asym_bits = 0UL;
    for(int i=0; i<12; i++){
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
    for(int i=0; i<12; i++){
        nzero_counter = nzero_counter + __builtin_popcount(sum[i]);
    }
    if (nzero_counter == 0){
        printf("Yah! Zerosum Property\n");
    }
    else{
        printf("# of nonzero bits = %d\n",nzero_counter);
    }
}
void check_equations(Lane *state){
    uint64_t eq0 = state[0] ^ state[4] ^ state[8];
    printf("Eqs = %ld \n",eq0);
    get_property(state);
}
Lane generate_random_lane(){
    uint64_t temp = xx_next(SEED);
    Lane lane = temp & MASK;
    return lane;
}
Lane generate_random_symmetric_lane(){
    uint64_t temp = xx_next(SEED);
    temp = temp & HALFMASK;
    Lane lane =  (temp << HALFLANE) | temp;
    return lane;
}
void generateLinearConstraintSymState(Lane *state){
    memset(state, 0x00, 12*sizeof(Lane));
    int nrof_accessed_lanes = 12;
    for(int i=0; i<nrof_accessed_lanes; i++){
        state[i] = generate_random_symmetric_lane();
    }
    /* state[controll_lane[0]] = 0; */
}
void constraint_1linear_sym(Lane *state, uint64_t mask){
    //create symmetric mask
    uint64_t mask_up = (mask >> (DIM/2));
    uint64_t mask_low = mask ^ (mask_up << (DIM/2));
    /* printf("Mask %016lX Mask_up %016lX Mask_low %016lX\n",mask, mask_up, mask_low); */
     
    state[5] = state[5] ^ mask_low ^ (mask_up << HALFLANE);
    state[9] = state[5];
}
void SUM(Lane *sum, Lane *state){
    for(int i=0; i<12; i++){
        sum[i] = sum[i] ^ state[i];
    }
}
void compute_sum(Lane *state, Lane *sum, Lane *sum_inv){
    uint64_t cube = (1UL << DIM);
    Lane temp[12];
    Lane temp_inv[12];

    for(uint64_t i=0UL; i<cube; i++){
        /* printf("**************************Input no %ld*************************\n",i); */
        memcpy(temp, state, 12*sizeof(Lane));
        uint64_t mask = i ;
        constraint_1linear_sym(temp, mask);
        /* constraint_1linear(temp, mask); */
        memcpy(temp_inv, temp, 12*sizeof(Lane));
        /* display_lanes("INPUT", temp); */
        /* check_equations(temp); */
        xoodoo_on_lane(temp, ROUNDS);
        xoodoo_inv_on_lane(temp_inv, ROUNDS);
        /* display_lanes("OUTPUT", temp); */
        /* display_lanes("OUTPUT", temp_inv); */
        SUM(sum, temp);
        SUM(sum_inv, temp_inv);
        /* display_lanes("SUM", sum, WIDTH); */
    }
    /* display_lanes("input_sum",input_sum,WIDTH); */
    printf("*******************************************************************************\n");
}

void* sum_over_range(void *args){
    SumArgs *sumargs = (SumArgs*) args;
    
    Lane temp[12];
    Lane temp_inv[12];
    Lane sum[12];
    memset(sum, 0x00, 12*sizeof(Lane));
    Lane sum_inv[12];
    memset(sum_inv, 0x00, 12*sizeof(Lane));

    for(uint64_t i=sumargs->start; i<sumargs->end; i++){
        /* printf("**************************Input no %ld*************************\n",i); */
        memcpy(temp, sumargs->state, 12*sizeof(Lane));
        uint64_t mask = i ;
        constraint_1linear_sym(temp, mask);
        /* constraint_1linear(temp, mask); */
        memcpy(temp_inv, temp, 12*sizeof(Lane));
        xoodoo_on_lane(temp, ROUNDS);
        xoodoo_inv_on_lane(temp_inv, ROUNDS);
        SUM(sumargs->result, temp);
        SUM(sumargs->result_inverse, temp_inv);
    }
    pthread_exit(NULL);
}
void compute_sum_threaded(Lane *state, Lane *sum, Lane *sum_inv){
    pthread_t thread_ids[NROF_THREADS]; 
    SumArgs thread_args[NROF_THREADS];
    
    uint64_t cube_size = (1ULL << DIM);
    for(int i=0; i < NROF_THREADS; i++){
        memset(thread_args[i].result, 0x00, 12*sizeof(Lane));
        memset(thread_args[i].result_inverse, 0x00, 12*sizeof(Lane));
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
        SUM(sum_inv, thread_args[i].result_inverse);
    }
}
void test_sysmsum_w_linear(){
    Lane state[12];
    generateLinearConstraintSymState(state);
    display_lanes("Starting: ", state);

    Lane sum[12];
    memset(sum, 0x00, 12*sizeof(Lane));
    display_lanes("SUM", sum);
    Lane sum_inv[12];
    memset(sum_inv, 0x00, 12*sizeof(Lane));
    display_lanes("SUM_inv", sum_inv);

    compute_sum_threaded(state, sum, sum_inv);
    display_lanes("SUM", sum);
    get_property(sum);
    display_lanes("SUM_inv", sum_inv);
    get_property(sum_inv);
}
int main(){
    xx_initialize(SEED);
    printf("SEED: ");
    printreg(SEED, 32);
    test_sysmsum_w_linear();
}

