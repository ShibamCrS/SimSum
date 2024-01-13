from math import *
from extra_eq import *
from sage.all import *
from termcolor import colored
###### UTILITY FUNCTIONS #########################
Lane =  64
F = GF(2)
def rotate_list_left(lane, k):
    k = k % Lane  # Handle rotation larger than list length
    return lane[k:] + lane[:k]
def rotate_list_right(lane, k):
    k = k % Lane  # Handle rotation larger than list length
    return lane[-k:] + lane[:-k]

Var =  {"22":["x22 "+str(i) for i in range(Lane)],
        "21":["x21 "+str(i) for i in range(Lane)],
        "20":["x20 "+str(i) for i in range(Lane)],
        "02":["x02 "+str(i) for i in range(Lane)],
        "01":["x01 "+str(i) for i in range(Lane)],
        "00":["x00 "+str(i) for i in range(Lane)],
        "O" :["F"]
        }
all_var = []
for key in Var.keys():
    all_var += Var[key]
print(all_var)
nVar = len(all_var)
def get_eq_3_var(a, b, c):
    equations = []
    for i in range(Lane):
        equations.append([a[i],b[i],c[i]])
    return equations
def get_eq_3_var_one(a, b, c):
    equations = []
    for i in range(Lane):
        equations.append([a[i],b[i],c[i],Var["O"][0]])
    return equations

def get_eq_2_var(a, b):
    equations = []
    for i in range(len(a)):
        equations.append([a[i],b[i]])
    return equations

def eq_2_row(eq):
    var_len = len(all_var)
    row = [0 for i in range(var_len)]
    for v in all_var:
        if v in eq:
            row[all_var.index(v)] = 1
    # print_eq_vs_row(eq, row)
    return row

def print_bit(value):
    if value:
        return colored("1", "red")
    else:
        return "."


def print_row_color(row):
    s = ""
    for i, bit in enumerate(row):
        # if (i == 32):
        #     s += "  "
        bit = print_bit(bit)
        s += f"{bit}"
    print(s)
    # print(row.count(1))
def print_matrix_color(M):
    for row in M:
        print_row_color(row)

def print_eqns(eqs):
    for eq in eqs:
        s = " + ".join(eq)
        print(s)
    print("************************************************************")
def print_comapct_matrix(L):
    for l in L:
        m = [str(i) for i in l]
        s = "".join(m)
        print(s)

def print_eq_vs_row(eq, row):
    print(len(row))
    s = "+".join(eq)
    print(s,"=>",end="")
    print_row_color(row)
    non_zero = []
    for i in range(len(row)):
        r = row[i]
        if r==1:
            non_zero.append(i)
    print(non_zero)

def find_first_nonzero_row(row):
    for i in range(len(row)):
        if row[i]:
            return i
def find_pivot(M):
    pivot =[]
    for row in M:
        p = find_first_nonzero_row(row)
        pivot.append(p)
    return pivot
def print_c_formated_matrix(M,f):
    s = "const uint8_t dep_mat[{0}][{1}] = ".format(len(M),len(M[0]))
    s += "{\n"
    t = []
    for row in M:
        row1 = [str(i) for i in row]
        row1 = "{" + ",".join(row1) + "}"
        t.append(row1)
    s += ",\n".join(t) + "\n};"
    f.write(s)

def check_eq_symm(row):
    for i in range(len(row)//2):
        r = row[i]
        r1 = row[i+32]
        if ( (r == 1) and (r1 == 1) ) or ( (r == 0) and (r1 == 0) ): 
            print(0, end="")
        else:
            print(1, end="")
    print("")

def eqs_2_matrix(equations):
    rows = []
    for eqs in equations:
        for eq in eqs:
            row = eq_2_row(eq)
            # row.append(0)
            rows.append(row)
    M = Matrix(F, rows) #, sparse = True)
    return M

def swap_list(LL, swap):
    L = LL[:]
    for s in swap:
        temp = L[s[0]]
        L[s[0]] = L[s[1]]
        L[s[1]] = temp
    return L

def find_free(N, swap =[]):
    print("*************************************************")
    print(all_var)
    all_var_1 = swap_list(all_var, swap)
    print(all_var_1)
    free_var = []
    free_var_idx = []
    pivot = find_pivot(N)
    for i in range(nVar):
        if i not in pivot:
            free_var.append(all_var_1[i])
            free_var_idx.append(i)

    return free_var, free_var_idx
def find_dep_eqns(N, free_var, free_var_idx, swap = []):
    all_var_1 = swap_list(all_var, swap)
    nRows = N.nrows()
    nCols = N.ncols()
    dep_eqns = []
    #describe dependent variables
    for row in N.rows():
        eqn = []
        for c in range(nCols):
            if row[c]:
                eqn.append(all_var_1[c])
        if len(eqn) !=0:
            dep_eqns.append(eqn)
    print_eqns([dep_eqns[0]])
    eqns_mat = {}
    for row in N.rows():
        eqn = []
        first_nonzero = False
        for c in range(nCols):
            if row[c]:
                if first_nonzero == False:
                    first_nonzero = True
                    key = all_var_1[c]
                else:
                    if c not in free_var_idx:
                        print(key, c,"ERROR: There are other free vars")
            # if first_nonzero == True:
            if c in free_var_idx:
                # print(c, row)
                eqn.append(row[c])
        if len(eqn) !=0: 
            eqns_mat[key] = eqn
    # print(eqns_mat)
    M = []
    zero_var = []
    for key in eqns_mat.keys():
        print(key,"=",end="")
        row = eqns_mat[key]
        #get 5,6,7,8,9 at first 5 positions
        row = row[22:27] + row[0:22] + row[27:]
        # row.reverse()
        M.append(row)
        print_row_color(row)
        if len(row) < len(free_var):
            print(len(row))
        all_zero = all(r == 0 for r in row)
        if all_zero:
            zero_var.append(key)
    print(len(eqns_mat))
    print(zero_var)

    return eqns_mat, M


"""
variable conventions
first 3*64 variables are from 0,0, 0,1, 0,2, then 2,0, 2,1, 2,2
first Lane is 0...63 which is the first set of variables
so variables are x00,0....x00,63
"""
def check_sym(row):
    hLane = Lane//2
    x = []
    for i in range(hLane):
        v = int(row[i]) ^ int(row[i+hLane])
        x.append(v)
    print(x)
#get equations
"""
A[0, 0] ⊕ A[0, 1] ⊕ A[0, 2] = F
"""
a = Var["00"]
b = Var["01"]
c = Var["02"]
equations0 = get_eq_3_var_one(a,b,c)
print_eqns(equations0)
"""
A[2, 0] ⊕ A[2, 1] ⊕ A[2, 2] = 0
"""
a = Var["20"]
b = Var["21"]
c = Var["22"]
equations1 = get_eq_3_var(a,b,c)
print_eqns(equations1)
"""
A[2, 0]≪ 62 = A[0, 0] ⊕ A[2, 2]≪ 43 
"""
a = Var["00"]
b = rotate_list_right(Var["20"], 62)
c = rotate_list_right(Var["22"], 43)
equations2 = get_eq_3_var(a,b,c)
print_eqns(equations2)
"""
A[2, 1]≪ 6 = A[0, 1]≪ 36
"""
a = rotate_list_right(Var["21"], 6)
b = rotate_list_right(Var["01"], 36)
equations3 = get_eq_2_var(a,b)
print_eqns(equations3)
"""
A[2, 2]≪ 43 = A[0, 2]≪ 3
"""
a = rotate_list_right(Var["22"], 43)
b = rotate_list_right(Var["02"], 3)
equations4 = get_eq_2_var(a,b)
print_eqns(equations4)


"""
Adding extra equations
"""
equations5 = []
for e in range(len(EQ)):
    if e == 0:
        continue
    l = EQ[e]
    eq = []
    sym_eq = []
    for i in range(Lane):
        if i in l:
            eq.append(Var["00"][i])
            sym_pos = (i+32)%Lane
            print("sym pos ",sym_pos)
            sym_eq.append(Var["00"][sym_pos])
    equations5.append(eq)
    equations5.append(sym_eq)
print_eqns(equations5)

equations = [equations0,equations1,equations2,equations3,equations4,equations5]

M = eqs_2_matrix(equations)
# S = M.right_kernel()
# print(S[1])
# print(S[513])
# check_sym(S[513])
# print(math.log(len(S),2) )
# print(str(len(S)))
# print(S.basis())


print_matrix_color(M)
print(M.rank()) 
swap = []
for i in range(5,10):
    s=[]
    s.append(320+i)
    s.append(320+i-5+32)
    swap.append(s)
for s in swap:
    M.swap_columns(s[0], s[1])
 
N = M.echelon_form()
print_matrix_color(N)
free_var, free_var_idx = find_free(N, swap)
eqns_mat, M = find_dep_eqns(N, free_var, free_var_idx, swap) 

print(f"{free_var_idx=}", f"{len(free_var_idx)=}")  
print(f"{free_var=}", f"{len(free_var)=}")  
print(free_var, len(free_var))
print(nVar)

# print(M)
f = open("mat54.py","w");
s = "M="+str(M)
f.write(s)
f.close()

f = open("dep_matrix_w_extra_w_sym.h","w")
print_c_formated_matrix(M,f)
f.close()

print(free_var_idx, len(free_var_idx))
print(free_var, len(free_var))
print(nVar)
#find dep var 0,0 lane
dep_var=[]
for i in Var["00"]:
    if i not in free_var:
        dep_var.append(i)
print(f"{dep_var=}")
print(f"{swap=}")

M1 = eqs_2_matrix([equations5])
print(M1.rank())
