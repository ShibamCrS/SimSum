#include "../dep_matrix_54_s.h"
#define NR_CLANE 6
int controll_lane[NR_CLANE] = {index(2,2),index(2,1),index(2,0),index(0,2),index(0,1),index(0,0)};

void get_dependent_value(uint64_t *dep_val, uint64_t free_val, int dep_var){
    for(int i=0; i<dep_var; i++){
        uint64_t temp = dep_mat[i] & free_val;
        dep_val[i] =(uint64_t) ((__builtin_popcountll(temp) ) % 2);
    }
}
//The dependent vars of the first lane are
//0,1,2,3,4,..................,32,33,34,35,35,..................
void set_sym_eqns(Lane *state){
    int dep_var = 330;
    uint64_t dep_val[dep_var];
    memset(dep_val, 0x00, dep_var*sizeof(uint64_t));

    //prepare mask
    uint64_t mask1 = (state[controll_lane[5]] >> 5) & 0x7FFFFFFUL;
    uint64_t mask2 = (state[controll_lane[5]] >> 37) & 0x7FFFFFFUL;
    uint64_t free_val = mask1 | (mask2 << 27)  | (1UL << (54));
    get_dependent_value(dep_val, free_val, dep_var);

    //Construct Lane from all_value
    for(int laneno=0; laneno<(NR_CLANE-1); laneno++){
        Lane val = 0UL;
        for(int i=0; i<LANE; i++){
            val = val | (dep_val[(LANE*laneno) + i] << i);
        }
        state[controll_lane[laneno]] = val;
    }

    //Now add .........36,35,34,33,32,............. 4,3,2,1,0
    state[controll_lane[5]] ^= dep_val[(LANE*5)];               //0
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 1] << 1);    //1
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 2] << 2);   //2 
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 3] << 3);   //3 
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 4] << 4);   //4 

    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 5] << 32);  //32
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 6] << 33);  //33
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 7] << 34);  //34 
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 8] << 35);  //35 
    state[controll_lane[5]] ^= (dep_val[(LANE*5) + 9] << 36);  //36 
}

