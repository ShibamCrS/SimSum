This repository contains the supporting material of the paper. This contains three directories. 
The description is given below.

1) Directory: KECCAKP
    Description: This directory contains the reference implementation of the Keccak algorithm and xoodoo. 

2) Directory: one_round_linearization
    Description: This directory contains an experiment of SymSumsim with the Keccak and Xoodoo using one-round linearization. 
    To customize the number of rounds, simply modify the "#define ROUNDS" parameter in the keccak.c/xoodoo.c file. 

3) Directory: two_round_linearization
    Description: This directory contains an experiment of SymSumsim with Keccak using two-round linearization.
   
    In this directory, you will find the program "find_free_variables.py" which is designed to identify the free variables and 
    dependency equations needed for the two-round linearization process. Additionally, the directory includes "extra_eq.py", where 
    you can find the supplementary equations required to neutralize the impact of round constants. The main verification codes are stored 
    in "keccak.c". 
    
    To customize the number of rounds, simply modify the "#define ROUNDS" parameter in the keccak.c file. 
