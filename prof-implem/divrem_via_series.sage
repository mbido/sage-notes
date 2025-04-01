# Reversal / mirror:
# in sagemath, see the documentation of reverse() for univariate polynomials
# --> p.reverse() reverses p in the usual definition ( x**deg(p) * p(1/x) )
# --> p.reverse(d) reverses p w.r.t degree d:  if d >= deg(p) this is x**d * p(1/x)

# Multiplying by x**d:
# it is faster to use p.shift(d) than to explicitly multiply by x**d

# Polynomial truncation:
# one can use p.truncate(n) to truncate a polynomial modulo x**n
# (in particular, this is useful to work with polynomials as power series
# without having to explicitly make conversions to power series rings)

def divrem_iterative(a,b):
    """division with remainder, iterative algorithm

    Input: univariate polynomials a,b with b != 0
    Output: univariate polynomials (q,r) which are the quotient
    and remainder in the division of a by b
    """
    # raise exception if b is zero
    if b.is_zero():
        raise ValueError("cannot divide by zero")

    # initially, q is 0 and r is a
    q = a.parent().zero()
    r = a

    # while deg(r) >= deg(b), we reduce the degree of r
    m = r.degree()
    n = b.degree()
    # inversion is expensive: we precompute the inverse of lc(b)
    # since we will use it several times
    inv_lc_b = 1 / b.leading_coefficient()
    while m >= n:
        cst = r.leading_coefficient() * inv_lc_b
        r -= cst * b.shift(m-n)
        q += cst * x**(m-n)
        m = r.degree()
    return (q,r)

def divrem_iterative_bis(a,b):
    """division with remainder, iterative algorithm

    Input: univariate polynomials a,b with b != 0
    Output: univariate polynomials (q,r) which are the quotient
    and remainder in the division of a by b

    Difference with divrem_iterative:
    q is represented via a list of coefficients, which improves timings because
    SageMath is not very efficient on the operation "adding a monomial"
    (used in the line q += cst * x**(m-n))
    """
    # raise exception if b is zero
    if b.is_zero():
        raise ValueError("cannot divide by zero")

    # initially, q is 0 and r is a
    q = [0]*(a.degree()-b.degree()+1)
    r = a

    # while deg(r) >= deg(b), we reduce the degree of r
    m = r.degree()
    n = b.degree()
    # inversion is expensive: we precompute the inverse of lc(b)
    # since we will use it several times
    inv_lc_b = 1 / b.leading_coefficient()
    while m >= n:
        cst = r.leading_coefficient() * inv_lc_b
        r -= cst * b.shift(m-n)
        q[m-n] = cst
        m = r.degree()
    return (a.parent()(q),r)



def divrem_via_series(a,b):
    """division with remainder via power series inversion

    Input: univariate polynomials a,b with b != 0
    Output: univariate polynomials (q,r) which are the quotient
    and remainder in the division of a by b
    """
    # raise exception if b is zero
    if b.is_zero():
        raise ValueError("cannot divide by zero")

    # get degrees and solve basic case deg(a) < deg(b)
    m = a.degree()
    n = b.degree()
    if m < n:
        return (0,a)

    # reverse a and b w.r.t their degrees
    arev = a.reverse()
    brev = b.reverse()

    # compute inverse + multiplication at precision m-n+1
    prec = m-n+1
    qrev = arev * brev.inverse_series_trunc(prec)
    qrev = qrev.truncate(m-n+1)

    # obtain q by reversing qrev
    # warning! the formulas in the algorithm show that qrev has to be reversed
    # w.r.t degree m-n, and possibly deg(qrev) < m-n, so we specify this explicitly
    q = qrev.reverse(m-n)

    # obtain r by the formula; observe that deg(r) < n by definition
    # so we can work modulo x**n, which may speed up computations if n << m
    r = (a.truncate(n) - b * q.truncate(n)).truncate(n)

    return (q,r)






test_bench = True # set to False if tests/benchs not wanted
field = GF(997)
ring.<x> = field[]

# lists of degrees for polynomial a
degs = sorted([2**e for e in range(10,20)] + [2**e + 2**(e+1) for e in range(10,20)])
# for polynomial b, we will use deg(a)+1, deg(a), deg(a) - 1, deg(a) - 5, deg(a)/2, deg(a)/5, deg(a)/20

limit_iter = 30000 # 30000 seems good for small prime fields

if test_bench:
    print("############################################################################")
    print("#  benchmarking and testing univariate polynomial division with remainder  #")
    print("#  - input polynomials a, b of respective degrees da, db                   #")
    print("#  - compare results with SageMath's result of quo_rem                     #")
    print("#  - compared time with:                                                   #")
    print("#         -- iterative division algorithm                                  #")
    print("#         -- SageMath's native quo_rem method                              #")
    print("############################################################################")
    import time
    print(f"da\tdb\ttest\titer\tfast\tSageMath")

    for da in degs:
        #a = ring.random_element(degree=da, num_bound=1000000) # num_bound only over QQ
        a = ring.random_element(degree=da)
        for db in [da+1, da, da-1, da - 5, da//2, da//5, da//20]:
            db = max(0, db) # we do not want negative db
            b = ring.random_element(degree=db)
            print(f"{da}\t{db}\t", end="")

            # own, iterative
            if da < limit_iter :
                tIter = time.perf_counter()
                qrIter = divrem_iterative(a,b)
                tIter = time.perf_counter() - tIter

            # own, iterative
            if da < limit_iter :
                tIterBis = time.perf_counter()
                qrIterBis = divrem_iterative_bis(a,b)
                tIterBis = time.perf_counter() - tIterBis

            # own, Fast
            tFast = time.perf_counter()
            qrFast = divrem_via_series(a,b)
            tFast = time.perf_counter() - tFast

            # native Sage
            tSage = time.perf_counter()
            qrSage = a.quo_rem(b)
            tSage = time.perf_counter() - tSage

            # Test
            if da < limit_iter and qrSage != qrIter:
                print("\n\n\n\nerror Iter: da = {da}, db = {db}\n\n\n")
            if da < limit_iter and qrSage != qrIterBis:
                print("\n\n\n\nerror IterBis: da = {da}, db = {db}\n\n\n")
            elif qrSage != qrFast:
                print("\n\n\n\nerror Fast: da = {da}, db = {db}\n\n\n")
            else :
                print("ok\t",end="")

            if da < limit_iter:
                print(f"{tIter:.4f}\t{tIterBis:.4f}\t",end="")
            else :
                print("inf\tinf\t",end="")
            print(f"{tFast:.4f}\t",end="")
            print(f"{tSage:.4f}")




# # here, field was GF(997)
# sage: %runfile divrem_via_series.sage
# ############################################################################
# #  benchmarking and testing univariate polynomial division with remainder  #
# #  - input polynomials a, b of respective degrees da, db                   #
# #  - compare results with SageMath's result of quo_rem                     #
# #  - compared time with:                                                   #
# #         -- iterative division algorithm                                  #
# #         -- SageMath's native quo_rem method                              #
# ############################################################################
# da      db      test    iter    fast    SageMath
# 1024    1025    ok      0.0000  0.0000  0.0000  0.0070
# 1024    1024    ok      0.0000  0.0000  0.0001  0.0000
# 1024    1023    ok      0.0001  0.0000  0.0001  0.0000
# 1024    1019    ok      0.0001  0.0001  0.0001  0.0000
# 1024    512     ok      0.0191  0.0150  0.0002  0.0001
# 1024    204     ok      0.0321  0.0227  0.0002  0.0001
# 1024    51      ok      0.0325  0.0205  0.0002  0.0001
# 2048    2049    ok      0.0000  0.0000  0.0000  0.0000
# 2048    2048    ok      0.0000  0.0000  0.0000  0.0000
# 2048    2047    ok      0.0001  0.0001  0.0001  0.0000
# 2048    2043    ok      0.0002  0.0001  0.0001  0.0000
# 2048    1024    ok      0.0629  0.0434  0.0003  0.0002
# 2048    409     ok      0.1111  0.0714  0.0003  0.0002
# 2048    102     ok      0.1311  0.0785  0.0003  0.0001
# 3072    3073    ok      0.0000  0.0000  0.0000  0.0000
# 3072    3072    ok      0.0000  0.0000  0.0000  0.0000
# 3072    3071    ok      0.0000  0.0000  0.0001  0.0000
# 3072    3067    ok      0.0002  0.0002  0.0001  0.0000
# 3072    1536    ok      0.1370  0.1102  0.0004  0.0003
# 3072    614     ok      0.2622  0.1734  0.0005  0.0003
# 3072    153     ok      0.3069  0.1852  0.0004  0.0002
# 4096    4097    ok      0.0000  0.0000  0.0000  0.0000
# 4096    4096    ok      0.0000  0.0000  0.0001  0.0000
# 4096    4095    ok      0.0001  0.0001  0.0001  0.0000
# 4096    4091    ok      0.0004  0.0003  0.0001  0.0000
# 4096    2048    ok      0.2656  0.1945  0.0006  0.0004
# 4096    819     ok      0.4840  0.3391  0.0006  0.0004
# 4096    204     ok      0.5736  0.3619  0.0006  0.0004
# 6144    6145    ok      0.0000  0.0000  0.0000  0.0000
# 6144    6144    ok      0.0000  0.0000  0.0001  0.0000
# 6144    6143    ok      0.0001  0.0001  0.0001  0.0000
# 6144    6139    ok      0.0003  0.0003  0.0001  0.0001
# 6144    3072    ok      0.6247  0.4900  0.0010  0.0006
# 6144    1228    ok      1.2099  0.8164  0.0011  0.0007
# 6144    307     ok      1.5028  0.8925  0.0010  0.0006
# 8192    8193    ok      0.0000  0.0000  0.0000  0.0000
# 8192    8192    ok      0.0000  0.0000  0.0001  0.0000
# 8192    8191    ok      0.0001  0.0001  0.0001  0.0000
# 8192    8187    ok      0.0004  0.0004  0.0002  0.0001
# 8192    4096    ok      1.1993  0.9097  0.0015  0.0011
# 8192    1638    ok      2.3048  1.5587  0.0016  0.0009
# 8192    409     ok      2.8600  1.7685  0.0015  0.0009
# 12288   12289   ok      0.0000  0.0000  0.0000  0.0000
# 12288   12288   ok      0.0001  0.0001  0.0001  0.0000
# 12288   12287   ok      0.0001  0.0001  0.0002  0.0000
# 12288   12283   ok      0.0005  0.0005  0.0002  0.0001
# 12288   6144    ok      3.1018  2.4672  0.0026  0.0017
# 12288   2457    ok      6.0322  4.1278  0.0028  0.0018
# 12288   614     ok      7.4559  4.6161  0.0025  0.0016
# 16384   16385   ok      0.0000  0.0000  0.0000  0.0000
# 16384   16384   ok      0.0001  0.0001  0.0002  0.0000
# 16384   16383   ok      0.0003  0.0002  0.0003  0.0001
# 16384   16379   ok      0.0007  0.0007  0.0003  0.0001
# 16384   8192    ok      6.3559  5.1590  0.0039  0.0028
# 16384   3276    ok      12.3870 8.4295  0.0042  0.0026
# 16384   819     ok      17.8804 11.6168 0.0047  0.0029
# 24576   24577   ok      0.0000  0.0000  0.0000  0.0000
# 24576   24576   ok      0.0001  0.0001  0.0004  0.0001
# 24576   24575   ok      0.0003  0.0003  0.0005  0.0001
# 24576   24571   ok      0.0014  0.0014  0.0007  0.0003
# 24576   12288   ok      25.2772 16.3583 0.0082  0.0053
# 24576   4915    ok      37.6846 25.0620 0.0090  0.0058
# 24576   1228    ok      41.6272 25.8723 0.0068  0.0043
# 32768   32769   ok      inf     inf     0.0000  0.0000
# 32768   32768   ok      inf     inf     0.0003  0.0001
# 32768   32767   ok      inf     inf     0.0006  0.0001
# 32768   32763   ok      inf     inf     0.0008  0.0003
# 32768   16384   ok      inf     inf     0.0103  0.0073
# 32768   6553    ok      inf     inf     0.0122  0.0072
# 32768   1638    ok      inf     inf     0.0095  0.0073
# 49152   49153   ok      inf     inf     0.0000  0.0001
# 49152   49152   ok      inf     inf     0.0007  0.0001
# 49152   49151   ok      inf     inf     0.0009  0.0002
# 49152   49147   ok      inf     inf     0.0011  0.0004
# 49152   24576   ok      inf     inf     0.0177  0.0119
# 49152   9830    ok      inf     inf     0.0190  0.0129
# 49152   2457    ok      inf     inf     0.0187  0.0108
# 65536   65537   ok      inf     inf     0.0000  0.0001
# 65536   65536   ok      inf     inf     0.0006  0.0002
# 65536   65535   ok      inf     inf     0.0011  0.0003
# 65536   65531   ok      inf     inf     0.0014  0.0008
# 65536   32768   ok      inf     inf     0.0239  0.0170
# 65536   13107   ok      inf     inf     0.0287  0.0177
# 65536   3276    ok      inf     inf     0.0251  0.0175
# 98304   98305   ok      inf     inf     0.0000  0.0001
# 98304   98304   ok      inf     inf     0.0012  0.0002
# 98304   98303   ok      inf     inf     0.0016  0.0004
# 98304   98299   ok      inf     inf     0.0018  0.0009
# 98304   49152   ok      inf     inf     0.0435  0.0269
# 98304   19660   ok      inf     inf     0.0437  0.0288
# 98304   4915    ok      inf     inf     0.0386  0.0254
# 131072  131073  ok      inf     inf     0.0000  0.0001
# 131072  131072  ok      inf     inf     0.0012  0.0003
# 131072  131071  ok      inf     inf     0.0029  0.0006
# 131072  131067  ok      inf     inf     0.0028  0.0014
# 131072  65536   ok      inf     inf     0.0586  0.0395
# 131072  26214   ok      inf     inf     0.0640  0.0417
# 131072  6553    ok      inf     inf     0.0634  0.0363
# 196608  196609  ok      inf     inf     0.0000  0.0002
# 196608  196608  ok      inf     inf     0.0021  0.0007
# 196608  196607  ok      inf     inf     0.0035  0.0009
# 196608  196603  ok      inf     inf     0.0043  0.0018
# 196608  98304   ok      inf     inf     0.0832  0.0568
# 196608  39321   ok      inf     inf     0.0951  0.0623
# 196608  9830    ok      inf     inf     0.0850  0.0549
# 262144  262145  ok      inf     inf     0.0000  0.0003
# 262144  262144  ok      inf     inf     0.0028  0.0006
# 262144  262143  ok      inf     inf     0.0041  0.0011
# 262144  262139  ok      inf     inf     0.0048  0.0023
# 262144  131072  ok      inf     inf     0.1290  0.0861
# 262144  52428   ok      inf     inf     0.1394  0.0898
# 262144  13107   ok      inf     inf     0.1274  0.0846
# 393216  393217  ok      inf     inf     0.0000  0.0005
# 393216  393216  ok      inf     inf     0.0034  0.0014
# 393216  393215  ok      inf     inf     0.0074  0.0019
# 393216  393211  ok      inf     inf     0.0085  0.0038
# 393216  196608  ok      inf     inf     0.1816  0.1317
# 393216  78643   ok      inf     inf     0.2124  0.1455
# 393216  19660   ok      inf     inf     0.1822  0.1227
# 524288  524289  ok      inf     inf     0.0021  0.0014
# 524288  524288  ok      inf     inf     0.0063  0.0016
# 524288  524287  ok      inf     inf     0.0102  0.0026
# 524288  524283  ok      inf     inf     0.0116  0.0055
# 524288  262144  ok      inf     inf     0.2393  0.1699
# 524288  104857  ok      inf     inf     0.3023  0.1935
# 524288  26214   ok      inf     inf     0.2683  0.1829
# 786432  786433  ok      inf     inf     0.0020  0.0022
# 786432  786432  ok      inf     inf     0.0095  0.0029
# 786432  786431  ok      inf     inf     0.0203  0.0048
# 786432  786427  ok      inf     inf     0.0181  0.0101
# 786432  393216  ok      inf     inf     0.4346  0.2961
# 786432  157286  ok      inf     inf     0.4604  0.3323
# 786432  39321   ok      inf     inf     0.4100  0.3079
# 1572864 1572865 ok      inf     inf     0.0000  0.0047
# 1572864 1572864 ok      inf     inf     0.0325  0.0043
# 1572864 1572863 ok      inf     inf     0.0456  0.0069
# 1572864 1572859 ok      inf     inf     0.0376  0.0206
# 1572864 786432  ok      inf     inf     1.0543  0.7199
# 1572864 314572  ok      inf     inf     1.0801  0.7625
# 1572864 78643   ok      inf     inf     1.1942  0.7332
