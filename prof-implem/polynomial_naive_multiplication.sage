################################################################
#  this file run tests and benchmarks for the naive algorithm  #
#  for multiplying two univariate polynomials                  #
################################################################

def mul_naive(f,g):
    # Input: f and g are two univariate polynomials
    # over some arbitrary ring (supported by SageMath)
    # Note: in this toy implementation we are supposed not to use
    # f*g anywhere, would be cheating !!
    ring = f.parent()

    if f.is_zero():
        return f
    if g.is_zero():
        return g

    # retrieve degrees and vectors of coefficients
    df = f.degree()
    dg = g.degree()

    coeff_f = f.list()
    coeff_g = g.list()

    # compute
    coeff_fg = [ring(0)]*(df+dg+1)
    for k in range(df+dg+1):
        for i in range(k+1):
            if i <= df and k-i <= dg:
                coeff_fg[k] += coeff_f[i]*coeff_g[k-i]

    # convert list to poly
    fg = ring(coeff_fg)
    return fg

# EXAMPLES
field = FiniteField(97)
ring.<x> = PolynomialRing(field)

# will not try our implementation beyond this degree
# (because this becomes very slow)
deg_max = 5000

import time

print(f"d\ttest\tnaive\tsage")

for e in range(1,20):
    d = 2**e
    print(f"{d}\t", end="")
    f = ring.random_element(degree=d)
    g = ring.random_element(degree=d)

    # our naive implementation
    if d < deg_max :
        tnaive = time.perf_counter()
        fg_naive = mul_naive(f,g)
        tnaive = time.perf_counter() - tnaive

    # SageMath's native code
    tsage = time.perf_counter()
    fg_sage = f*g
    tsage = time.perf_counter() - tsage

    # Test
    if d < deg_max and fg_sage != fg_naive:
        print("error\t",end="")
    elif d < deg_max:
        print("ok\t",end="")
    else :
        print("notest\t",end="")

    if d < deg_max :
        print(f"{tnaive:.4f}\t",end="")
    else :
        print("inf\t",end="")
    print(f"{tsage:.4f}")
