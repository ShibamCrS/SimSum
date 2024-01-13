void index_map(){
    int x, y;
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
            printf("%d %d %d\n", x, y, index(x, y));
        }
    }
    printf("\n");
    for(y=0; y<5; y++){
        int temp[15];
        for(x=0; x<5; x++){
            temp[x] = x;
            temp[x+5] = y;
            temp[x+10] = index(x, y);
        }
        for(x=0; x<5; x++){
            printf("%d, %d | ", temp[x], temp[x+5]);
        }
        printf("   ====>>>   ");
        for(x=0; x<5; x++){
            printf("%02d | ", temp[x+10]);
        }
        printf("\n");
    }
    printf("\n");
}

/* void test_rol64(){ */
/*     U64 one = (1UL << 63); */
/*     U64 rotate; */
/*     for(U32 i=0; i<64; i++){ */
/*         print_bits(one, 64); */
/*         printf(" => "); */
/*         rotate = ROL64(one, i); */
/*         //printf("%d => %lu => ",i, (unsigned long)rotate); */
/*         //print_bits(rotate, 64); */
/*         printf("\n"); */
/*     } */
    
/* } */

void pi_position_map(){
    int x,y, newX, newY;
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
            newX = 0*x + 1*y;
            newY = 2*x + 3*y;
            printf("(%02d %02d) %02d => (%02d %02d) %02d\n", x, y, index(x, y), newX, newY, index(newX, newY));    
        }
    }
    
    printf("\n");
    int old_stateX[5][5];
    int old_stateY[5][5];
    int new_stateX[5][5];
    int new_stateY[5][5];
    for(x=0; x<5; x++){
        for(y=0; y<5; y++){
           old_stateX[x][y] = x;
           old_stateY[x][y] = y; 
           newX = 0*x + 1*y;
           newY = (2*x + 3*y) % 5;
           new_stateX[newX][newY] = x;
           new_stateY[newX][newY] = y;
        }
    }
    
    for(y=0; y<5; y++){
        for(x=0; x<5; x++){
            printf("%d, %d | ", old_stateX[x][y], old_stateY[x][y]);
        }
        printf("   ====>>>   ");
        for(x=0; x<5; x++){
            printf("%d, %d | ", new_stateX[x][y], new_stateY[x][y]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_constants_rho_offsets(){
    for (int i=0; i<maxNrRounds; i++){
        printreg(RC+i, 8);
    }
    for (int i=0; i<5; i++){
        for (int j=0; j<5; j++){
            printf("%02d |",(unsigned int)RO[5*i+j]);
        }
        printf("\n");
    }
    printf("\n");
}

FILE *data_file = 0;

void set_file(FILE *f){
    data_file = f;
}

void display_bytes(const char *text, const U8 *bytes, const U32 size){
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for(U32 i=0; i<size; i++){
            fprintf(data_file, "%02X ", bytes[i]);
        }
        fprintf(data_file, "\n");
    }
}

void display_state_bytes(const char *text, const unsigned char *state, unsigned int width)
{
    display_bytes(text, state, width/8);
}

void display_state_64bit_lanes(const char *text, void *data)
{
    U32 i;
    U64 *state = data;
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for(U32 i=0; i<25; i++){
            fprintf(data_file, "%08X", (unsigned int)(state[i] >> 32));
            fprintf(data_file, "%08X ", ( (unsigned int)(state[i] & 0XFFFFFFFF) ) );
            
            if ((i%5) == 4){
                fprintf(data_file, "\n");
            }
            else{
                fprintf(data_file, " ");
            }
        }
        fprintf(data_file, "\n");
    }
}

void display_state_32bit_lanes(const char *text, void *data)
{
    U32 i;
    U32 *state = data;
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for(U32 i=0; i<25; i++){
            fprintf(data_file, "%08X ", (unsigned int)(state[i]));
            
            if ((i%5) == 4){
                fprintf(data_file, "\n");
            }
            else{
                fprintf(data_file, " ");
            }
        }
        fprintf(data_file, "\n");
    }
}

void display_state_16bit_lanes(const char *text, void *data)
{
    U32 i;
    U16 *state = data;
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for(U32 i=0; i<25; i++){
            fprintf(data_file, "%04X ", (unsigned int)(state[i]));
           
            if ((i%5) == 4){
                fprintf(data_file, "\n");
            }
            else{
                fprintf(data_file, " ");
            }
        }
        fprintf(data_file, "\n");
    }
}

void display_state_8bit_lanes(const char *text, void *data)
{
    U32 i;
    U8 *state = data;
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for(U32 i=0; i<25; i++){
            fprintf(data_file, "%02X ", (unsigned int)(state[i]));
            
            if ((i%5) == 4){
                fprintf(data_file, "\n");
            }
            else{
                fprintf(data_file, " ");
            }
        }
        fprintf(data_file, "\n");
    }
}

void display_lanes(const char *text, void *data, unsigned int width){
    if (width == 1600)
        display_state_64bit_lanes(text, data);
    if (width == 800)
        display_state_32bit_lanes(text, data);
    if (width == 400)
        display_state_16bit_lanes(text, data);
    if (width == 200)
        display_state_8bit_lanes(text, data);
}

//Disply bits: LSB FIRST
void display_bits(const char *text, const U8 *data, const U32 size){
    U32 i, byte_index, bit_index;
    
    if(data_file){
        fprintf(data_file, "%s:\n", text);
        for( i=0; i<size; i++){
            byte_index = i/8;
            bit_index = i%8;
            fprintf(data_file, "%d", ( (data[byte_index] >> bit_index) & 0x01 ) );
        }
        fprintf(data_file, "\n");
    }

}


void display_round_number(U32 i)
{
    if(data_file){
        fprintf(data_file, "\n");
        fprintf(data_file, "--- Round %d ---\n", (unsigned int)i);
        fprintf(data_file, "\n");
    }
}
