# Problem V.4 from poly-flag (version March 2024)

field = GF(7)

def sequence(i):
    # input: integer index i
    # output: sequence term u_i, in GF(7)
    return field(2**i + i*i - 1)

# build first 20 terms
terms20 = [sequence(i) for i in range(20)]

print("######################################")
print("#  solution based on linear algebra  #")
print("######################################")

print(f"\nlinear system with 5 equations and 5 unknowns,")

# build "special" matrix from the sequence terms
# (corresponds to equations giving recurrence of order 4)
mat = matrix.toeplitz(terms20[4:9], list(reversed(terms20[:4])))
print(f"matrix of the system:\n{mat}")
# sage: mat
# [3 2 0 2 0]
# [0 3 2 0 2]
# [1 0 3 2 0]
# [1 1 0 3 2]
# [4 1 1 0 3]

# solve system == consider right kernel
v = mat.right_kernel_matrix()
# v is [1 2 2 0 2], gives a recurrence
print(f"recurrence given by coefficients of solution: {v}")


# observe that one could use more equations:
# we always find a single solution (and the same one)
# this is summarized in the fact that the matrix of the system has rank 4
# whenever we use 5 or more equations
for k in range(10, 15):
    matk = matrix.toeplitz(terms20[4:k], list(reversed(terms20[:4])))
    print(f"same solution if using {k-4} equations? --> {matk.right_kernel_matrix() == mat.right_kernel_matrix()}")


print()
print("##############################################################")
print("#  solution based on rational reconstruction / extended GCD  #")
print("##############################################################")

ring.<x> = field[]
A = x**10
B = 2*x^8 + 2*x^6 + 3*x^5 + x^3 + x^2 + 4*x + 4
# B is defined from the first 10 terms, reversed
load("euclidean_algorithm.sage")
ExtendedEuclidAlgorithm(A,B,verbose=True)
print("We see that at step 4, deg(r0) < deg(v) (we could have stopped the computation here")
print("--> the corresponding polynomial v = 4*x^4 + x^3 + x^2 + 1 gives the recurrence of order 4")
print("    (note that 2*v = x^4 + 2*x^3 + 2*x^2 + 2, this is the same recurrence as with the linear algebra solution)")
