#include "xoshiro256plusplus.h"
#define NROF_THREADS_FHT  8
#define NROF_THREADS_DATA 8
#define NROF_THREADS_SUM 8        

uint64_t SEED[4];
//Generating Randoms
void generate_random(uint8_t *random_data, int nrof_byte){
    int required_call = nrof_byte/8 + 1;
    uint64_t *random64 = (uint64_t *)malloc(required_call*sizeof(uint64_t));
    for(int i=0; i<required_call; i++){
        random64[i] = xx_next(SEED);
    }
    memcpy(random_data, random64, nrof_byte);
    free(random64);
}
void generate_random_nibble(uint8_t *random_data){
    uint64_t random64 = xx_next(SEED);
    for(int i=0; i<16; i++){
        random_data[i] = (random64 >> (4*i)) & 0x0F; 
    }
}
void generate_random_round_keys(uint8_t **rk, int rounds){
    for(int i=0; i<rounds+1; i++){
        rk[i] = (uint8_t *)malloc(16*sizeof(uint8_t));
        generate_random_nibble(rk[i]);
    }
}
void value_to_array(uint8_t *arr, uint64_t v, int arr_size){
    for(int i=0; i<arr_size; i++){
	    arr[i] = (v >> (8*i)) & 0xFF;
	}

}
void value_to_nibble_array(uint8_t *arr, uint32_t v, int arr_size){
    for(int i=0; i<arr_size; i++){
	    arr[i] = (v >> (4*i)) & 0x0F;
	}

}

//Saving, reading, Checking Staffs
void save_array_double(FILE *fp, double *data, uint64_t size){
    uint64_t written = 0UL;
    written = fwrite(data, sizeof(double), size, fp);
    if (written < size){
        printf("Error during writing to file !");
    }
    fclose(fp);
}

void read_array_double(FILE *fp, double *data, uint64_t size){
    uint64_t read = 0;
    read = fread(data, sizeof(double), size, fp);
    
    if (read < size){
        printf("Read %ld ", (unsigned long int)read);
        printf("Error during reading from file !\n");
    }
    fclose(fp);
}
void save_array_float(FILE *fp, float *data, uint64_t size){
    uint64_t written = 0UL;
    written = fwrite(data, sizeof(float), size, fp);
    if (written < size){
        printf("Error during writing to file !");
    }
    fclose(fp);
}

void read_array_float(FILE *fp, float *data, uint64_t size){
    uint64_t read = 0;
    read = fread(data, sizeof(float), size, fp);
    
    if (read < size){
        printf("Read %ld ", (unsigned long int)read);
        printf("Error during reading from file !\n");
    }
    fclose(fp);
}
void save_array_int64(FILE *fp, int64_t *data, uint64_t size){
    uint64_t written = 0UL;
    written = fwrite(data, sizeof(int64_t), size, fp);
    if (written < size){
        printf("Error during writing to file !");
    }
    fclose(fp);
}

void read_array_int64(FILE *fp, int64_t *data, uint64_t size){
    uint64_t read = 0;
    read = fread(data, sizeof(int64_t), size, fp);
    
    if (read < size){
        printf("Read %ld ", (unsigned long int)read);
        printf("Error during reading from file !\n");
    }
    fclose(fp);
}

void save_array_uint64(FILE *fp, uint64_t *data, uint64_t size){
    uint64_t written = 0UL;
    written = fwrite(data, sizeof(uint64_t), size, fp);
    if (written < size){
        printf("Error during writing to file !");
    }
    fclose(fp);
}

void read_array_uint64(FILE *fp, uint64_t *data, uint64_t size){
    uint64_t read = 0;
    read = fread(data, sizeof(uint64_t), size, fp);
    
    if (read < size){
        printf("Read %ld ", (unsigned long int)read);
        printf("Error during reading from file !\n");
    }
    fclose(fp);
}




int mem_cmp(const void *a, const void *b, uint32_t nrof_byte){
    if(memcmp(a, b, nrof_byte) == 0){
        return 1;
    }
    else{
        return 0;
    }
}

//Printing Staffs
void printreg_to_file(const void *a, int nrof_byte, FILE *fp){
    int i;
    unsigned char *f = (unsigned char *)a;
    for(i=0; i < nrof_byte; i++){
        fprintf(fp, "%02X ",(unsigned char) f[i]); //uint8_t c[4+8];
    }
    fprintf(fp, "\n");
}
void printreg(const void *a, int nrof_byte){
    printreg_to_file(a, nrof_byte, stdout);
}
void print_arr(int *arr, int size){
    for (int i=0; i<size; i++){
        printf("%d ",arr[i]);
    }
    printf("\n");
}
void print_arr_uint8(uint8_t *arr, int size){
    for (int i=0; i<size; i++){
        printf("0x%02X,",(unsigned char)arr[i]);
    }
    printf("\n");
}

void print_arr_int64(int64_t *arr, int size){
    for (int i=0; i<size; i++){
        printf("%ld ",(long int)arr[i]);
    }
    printf("\n");
}
void print_arr_float(float *arr, int size){
    for (int i=0; i<size; i++){
        printf("%f ",arr[i]);
    }
    printf("\n");
}
void print_arr_double(double *arr, int size){
    for (int i=0; i<size; i++){
        printf("%lf ",arr[i]);
    }
    printf("\n");
}
