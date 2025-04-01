import time

field = GF(97)
ring.<x> = field[[]]

def power_series_inversion_iterative(s):
    """
    Input: a power series `s` with coefficients in a field, at precision `p`,
    with nonzero constant coefficient
    Output: the inverse of `s` at precision `p`
    """
    if s[0] == 0:
        raise ValueError("input power series must be invertible")

    p = s.prec()

    # -> we have access to the coefficients k = 0, 1, 2, ..., p-1 
    # of s using s[k]
    # for more efficiency we retrieve a copy of the list of coefficients
    # via s.padded_list(p)  (padded -> length is p even if zeros in higher
    # degree coefficients)
    s_coeffs = s.padded_list(p)

    # we will build the list of coefficients of the inverse
    s_inv_coeffs = [field.zero() for i in range(p)]

    # constant coefficient of inverse is inverse of constant coefficient
    cst_inv = s[0].inverse_of_unit()
    s_inv_coeffs[0] = cst_inv

    # then, others are deduced iteratively, essentially solving a
    # linear system of equations but with triangular matrix
    for i in range(1,p):
        s_inv_coeffs[i] = - cst_inv * sum(s_inv_coeffs[k] * s_coeffs[i-k] for k in range(i))

    # finally convert s_inv_coeffs into a power series s_inv
    s_inv = ring(s_inv_coeffs) + O(x**p)

    return s_inv

# now let's use SageMath's linear system solving; even though SageMath is
# calling very efficient low-level software for linear system solving, we will
# see that in the present case it is not enough to beat the above iterative
# algorithm, although it is essentially non-optimized python code. The reason
# is that this is system solving with a p x p matrix which is triangular. This has
# complexity O(p^2). The fast solver would not have better complexity, but
# could have better practical efficiency.

def power_series_inversion_via_system_solving(s):
    """
    Input: a power series `s` with coefficients in a field, at precision `p`,
    with nonzero constant coefficient
    Output: the inverse of `s` at precision `p`
    """
    if s[0] == 0:
        raise ValueError("input power series must be invertible")

    p = s.prec()

    # write the matrix of the linear system
    mat = Matrix(field, p, p)
    # retrieve coefficients of s in a list, for efficiency
    s_coeffs = s.padded_list(p) # length p even if there are some zeros in the rightmost coefficients
    # row i = first i coefficients of s, in reversed order
    for i in range(p):
        for j in range(i+1):
            mat[i,j] = s_coeffs[i-j]

    # solve linear system, be careful that depending on how the matrix was
    # written one should use solve_left or solve_right, and the vector of may
    # have to be either [1,0,...,0] or [0,...,0,1]
    # -> gives the list of coefficients of the inverse, as a vector of length p
    #vec = vector(field, [0 for i in range(p-1)] + [1])
    vec = vector(field, [1] + [0 for i in range(p-1)])
    s_inv_coeffs = mat.solve_right(vec)

    # finally convert s_inv_coeffs into a power series s_inv
    # -> note that due to how we expressed the system, we have to reverse
    # the list of coefficients
    #s_inv = ring(list(reversed(s_inv_coeffs))) + O(x**p)
    s_inv = ring(list(s_inv_coeffs)) + O(x**p)

    return s_inv

def power_series_inversion_newton(s, threshold=32):
    """
    Input: a power series `s` with coefficients in a field, at precision `p`,
    with nonzero constant coefficient
    Output: the inverse of `s` at precision `p`
    """
    if s[0] == 0:
        raise ValueError("input power series must be invertible")

    p = s.prec()

    # base case: invert constant coefficient
    if p == 1:
        return s[0].inverse_of_unit() + O(x)

    # second base case: under threshold, iterative algo
    if p <= threshold:
        return power_series_inversion_iterative(s)

    # recursion at half precision
    prec = ceil(p/2)

    # truncate s to obtain it at precision p/2
    srec = s.truncate_powerseries(prec)

    # compute recursively the inverse at precision p/2
    inv_rec = power_series_inversion_newton(srec)

    # we lift inv_rec to precision p (coefficients p/2 ... p-1 are zero)
    inv_rec = inv_rec.lift_to_precision(p)

    # we use the formula of Newton iteration and return
    s_inv = inv_rec + (1 - inv_rec * s) * inv_rec

    return s_inv

#############
#  testing  #
#############

# iterate on "i" random examples for a few precisions "p"
for p in [1,2,5,7,10,15,20]:
    for i in range(5):
        s = ring.random_element(prec=p)
        while not s.is_unit():
            s = ring.random_element(prec=p)
        inv1 = s.inverse_of_unit()
        inv2 = power_series_inversion_newton(s)
        inv3 = power_series_inversion_iterative(s)
        inv4 = power_series_inversion_via_system_solving(s)
        if inv1 != inv2:
            raise ValueError(f"Newton returned wrong output, input {s}")
        if inv1 != inv3:
            raise ValueError(f"Iterative returned wrong output, input {s}")
        if inv1 != inv4:
            raise ValueError(f"System returned wrong output, input {s}")
print("All tests passed.\n")


##################
#  benchmarking  #
##################

# some basic benchmarking
print("prec\tsystem\t\titer\t\tnewton\t\tsage")
for p in [2**k for k in range(2,10)]:
    s = ring.random_element(prec=p)
    if s[0] == 0: # make sure invertible
        s = s + 1
    t_iter = timeit('power_series_inversion_iterative(s)',seconds=True)
    t_newton = timeit('power_series_inversion_newton(s)',seconds=True)
    t_sage = timeit('s.inverse_of_unit()',seconds=True)
    t_system = timeit('power_series_inversion_via_system_solving(s)',seconds=True)
    print(f"{p}\t{t_system:.9f}\t{t_iter:.9f}\t{t_newton:.9f}\t{t_sage:.9f}")

# second loop with fewer repeats to make it finish in reasonable time
for p in [2**k for k in range(10,13)]:
    s = ring.random_element(prec=p)
    if s[0] == 0: # make sure invertible
        s = s + 1
    t_iter = timeit('power_series_inversion_iterative(s)',seconds=True,number=1,repeat=1)
    t_newton = timeit('power_series_inversion_newton(s)',seconds=True,number=1,repeat=1)
    t_sage = timeit('s.inverse_of_unit()',seconds=True,number=2,repeat=2)
    t_system = timeit('power_series_inversion_via_system_solving(s)',seconds=True,number=2,repeat=2)
    print(f"{p}\t{t_system:.9f}\t{t_iter:.9f}\t{t_newton:.9f}\t{t_sage:.9f}")

# second loop with fewer repeats to make it finish in reasonable time
for p in [2**k for k in range(13,20)]:
    s = ring.random_element(prec=p)
    if s[0] == 0: # make sure invertible
        s = s + 1
    t_iter = "inf"
    t_system = "inf"
    t_sage = timeit('s.inverse_of_unit()',seconds=True,number=2,repeat=2)
    t_newton = timeit('power_series_inversion_newton(s)',seconds=True,number=2,repeat=2)
    print(f"{p}\t{t_system}\t{t_iter}\t{t_newton:.9f}\t{t_sage:.9f}")

# Output on a laptop, Sage 9.8i.rc1 2023-02-05:
## NEWTON THRESHOLD == 1
# sage: %runfile power_series_inversion.sage
# prec    system          iter            newton          sage
# 4       0.000191775     0.000040355     0.000130129     0.000016534
# 8       0.000218255     0.000059316     0.000189696     0.000017703
# 16      0.000319391     0.000077964     0.000255348     0.000019067
# 32      0.000523280     0.000147982     0.000341041     0.000024542
# 64      0.001106322     0.000366664     0.000463297     0.000022959
# 128     0.002875087     0.001114878     0.000656464     0.000027114
# 256     0.009411980     0.003878329     0.001001829     0.000036272
# 512     0.035600526     0.015372913     0.001638065     0.000056994
# 1024    0.155295716     0.065962422     0.002992279     0.000117508
# 2048    0.686463910     0.272764295     0.005582228     0.000249340
# 4096    3.399530489     1.108522120     0.010759986     0.000604690
# 8192    inf             inf             0.020985775     0.001552539
# 16384   inf             inf             0.041717443     0.004084622
# 32768   inf             inf             0.084543717     0.010619092
# 65536   inf             inf             0.173105048     0.020976916
# 131072  inf             inf             0.354077917     0.044254696
# 262144  inf             inf             0.714133797     0.106629628
# 524288  inf             inf             1.461050352     0.257780140

# prec    system          iter            newton          sage
# 4       0.000191202     0.000040514     0.000042241     0.000016537
# 8       0.000215976     0.000051615     0.000053037     0.000017531
# 16      0.000318153     0.000077267     0.000078790     0.000019134
# 32      0.000523374     0.000146734     0.000148292     0.000020757
# 64      0.001107259     0.000362081     0.000277887     0.000022919
# 128     0.002889170     0.001106462     0.000474640     0.000026876
# 256     0.009393220     0.003835651     0.000819883     0.000036195
# 512     0.035677191     0.015277246     0.001456396     0.000056957
# 1024    0.155088359     0.065674293     0.002904873     0.000112585
# 2048    0.684066080     0.272640355     0.005304932     0.000249277
# 4096    3.422337360     1.112278218     0.010939700     0.000606120
# 8192    inf             inf             0.020694443     0.001542621
# 16384   inf             inf             0.041565020     0.004059115
# 32768   inf             inf             0.083286922     0.010556829
# 65536   inf             inf             0.171705221     0.020951680
# 131072  inf             inf             0.350583752     0.044680194
# 262144  inf             inf             0.711033386     0.107257531
# 524288  inf             inf             1.458671547     0.257619065

# --> without surprises, setting a threshold in Newton is useful for low-precision
# cases, but not at all for larger ones: there is only one recursive call, so
# unlike e.g. Karatsuba, the time for the base case is negligible in the total time

# Note: the comparison is striking between iter/Newton, but is also a bit
# unfair: the Newton version is not only a better algorithm, it also exploits
# extremely optimized polynomial multiplication. Yet:
# --> forgetting how the two methods compare to each other, we clearly see that
# iter has a quadratic behaviour whereas Newton has an almost linear behaviour
# --> with naive multiplication, both algorithms would be O(n^2) (and probably
# with comparable timings if implemented naively in Sage): this still
# illustrates that the strategy or reducing the bulk of the work to fundamental
# building blocks (here univariate polynomial multiplication, in some other
# algos matrix multiplication) is indeed very efficient!

# finally, we notice that Sage's native version of Newton inversion (which must
# be based on Flint or NTL) is significantly better than our Sage-written Newton
