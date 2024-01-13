#include <stdint.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dep_matrix_w_extra_w_sym.h"

void generate_smart_matrix_54(){
    FILE *fp = fopen("dep_matrix_54_s.h","w");
    fprintf(fp, "const uint64_t dep_mat[330] = {\n");
    for(int i=0; i<330; i++){
        uint64_t row = 0UL;
        for(int j=0; j<55; j++){
            uint64_t temp = (uint64_t) dep_mat[i][j];
            row = row |(temp << j); 
        }
        if(i != 329)
            fprintf(fp, "0x%016lX,\n", row);
        else
            fprintf(fp, "0x%016lX};\n", row);
    }
    fclose(fp);
}
int main(){
    generate_smart_matrix_54();
}
