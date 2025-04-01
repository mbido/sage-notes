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

def mul_karatsuba(f,g,threshold=0):
    # Input: f and g are two univariate polynomials
    # over some arbitrary ring (supported by SageMath)
    # Note: in this toy implementation we are supposed not to use
    # f*g anywhere, would be cheating !!
    if f.is_zero():
        return f
    if g.is_zero():
        return g
    if f.degree() == 0:
        c = f.constant_coefficient()
        return c * g
    if g.degree() == 0:
        c = g.constant_coefficient()
        return c * f

    deg = max(f.degree(), g.degree())
    if deg <= threshold:
        return mul_naive(f, g)

    deg_half = ceil(deg/2)

    # f = flow + x**(deg/2) fhi
    flow = f.truncate(deg_half)
    fhi = f.shift(-deg_half)
    # g = glow + x**(deg/2) ghi
    glow = g.truncate(deg_half)
    ghi = g.shift(-deg_half)

    # apply the formulae
    fg0 = mul_karatsuba(flow, glow, threshold)
    fg2 = mul_karatsuba(fhi, ghi, threshold)
    fg1 = mul_karatsuba(flow+fhi, glow+ghi, threshold) - fg2 - fg0

    # gather and return
    fg = fg2.shift(2 * deg_half) + fg1.shift(deg_half) + fg0
    return fg


# EXAMPLES
field = FiniteField(97)
ring.<x> = PolynomialRing(field)

import time

thresholds = [0,1,2,5,10,20,40,60,100]
print(f"d\ttest\tnaive\t",end="")
for thres in thresholds:
    print(f"kar-{thres}\t",end="")
print(f"sage")

for e in range(1,20):
    d = 2**e
    print(f"{d}\t", end="")
    f = ring.random_element(degree=d)
    g = ring.random_element(degree=d)

    # own, naive
    if d < 3000:
        tnaive = time.perf_counter()
        fg_naive = mul_naive(f,g)
        tnaive = time.perf_counter() - tnaive

    # own, Karatsuba, varying thresholds
    if d < 9000:
        tkara = []
        fgkara = []
        for thres in thresholds:
            t = time.perf_counter()
            fg_kara = mul_karatsuba(f,g,thres)
            t = time.perf_counter() - t
            tkara.append(t)
            fgkara.append(fg_kara)

    # naive Sage
    tsage = time.perf_counter()
    fg_sage = f*g
    tsage = time.perf_counter() - tsage

    # Test
    if d < 3000 and fg_sage != fg_naive:
        print("error\t",end="")
    elif d < 9000 and any([fg_sage != fgk for fgk in fgkara]):
        print("error\t",end="")
    else :
        print("ok\t",end="")

    if d < 3000:
        print(f"{tnaive:.4f}\t",end="")
    else :
        print("inf\t",end="")
    if d < 9000:
        for tk in tkara:
            print(f"{tk:.4f}\t",end="")
    else :
        for t in thresholds:
            print(f"inf\t",end="")
    print(f"{tsage:.4f}")
