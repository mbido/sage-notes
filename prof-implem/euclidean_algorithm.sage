def EuclidAlgorithm(a,b):
    r0 = a
    r1 = b
    i = 1
    while r1 != 0:
        # Note: the next three lines could be done shortly, hiding the temporary,
        # using:   (q, r1), r0 = r0.quo_rem(r1), r1
        tmp = r0
        r0 = r1
        (q,r1) = tmp.quo_rem(r1)
        i += 1
    return r0

def ExtendedEuclidAlgorithm(a,b,verbose=True):
    if verbose:
        print("Running Extended Euclidean algorithm...")
    r0, u0, v0 = a, 1, 0
    r1, u1, v1 = b, 0, 1
    i=1
    while r1 != 0:
        # Note: the next three lines could be done shortly, hiding the temporary,
        # using:   (q, r1), r0 = r0.quo_rem(r1), r1
        tmp = r0
        r0 = r1
        (q,r1) = tmp.quo_rem(r1)
        # Note: the next three lines could be done shortly, hiding the temporary,
        # using:    u1, u0 = u0 - q*u1, u1
        tmp = u1
        u1 = u0 - q*u1
        u0 = tmp
        # Note: the next three lines could be done shortly, hiding the temporary,
        # using:    v1, v0 = v0 - q*v1, v1
        tmp = v1
        v1 = v0 - q*v1
        v0 = tmp
        if verbose:
            print(f"Step {i}:\n\tu = {u0},\n\tv = {v0}")
        i += 1
    if verbose:
        print("---------------------")
    return (r0, u0, v0)

## comment/uncomment here to have the field you want
field = QQ
#field = GF(7)
pring.<x> = field[]
p=x^3+2*x^2-x-2
q=x^2+1

print(f"polynomials in input:\n\tp = {p}\n\tq = {q}")
print(f"base field: {field}")

print("---------------------")

g = EuclidAlgorithm(p,q)
print(f"**their GCD via the Euclidean algorithm:\n\t{g}")
print("---------------------")

(gext,u,v) = ExtendedEuclidAlgorithm(p,q)
print("**their GCD g and cofactors (u,v) such that u p + v q = g via the extended Euclidean algorithm:")
print(f"\tu = {u}\n\tv = {v}\n\tg = {gext}")
print("---------------------")

print(f"**GCD returned by SageMath is ``normalized'': gcd(p,q) = {gcd(p,q)}")
print("---------------------")

(gsage,usage,vsage) = xgcd(p,q)
print("furthermore, cofactors returned by SageMath are:")
print(f"\tu = {usage}\n\tv = {vsage}")


# OUTPUT (FIELD = RATIONAL NUMBERS)
#polynomials in input:
#        p = x^3 + 2*x^2 - x - 2
#        q = x^2 + 1
#base field: Rational Field
#---------------------
#**their GCD via the Euclidean algorithm:
#        5
#---------------------
#Running Extended Euclidean algorithm...
#Step 1:
#        u = 0,
#        v = 1
#Step 2:
#        u = 1,
#        v = -x - 2
#Step 3:
#        u = 1/2*x - 1,
#        v = -1/2*x^2 + 3
#---------------------
#**their GCD g and cofactors (u,v) such that u p + v q = g via the extended Euclidean algorithm:
#        u = 1/2*x - 1
#        v = -1/2*x^2 + 3
#        g = 5
#---------------------
#**GCD returned by SageMath is ``normalized'': gcd(p,q) = 1
#---------------------
#furthermore, cofactors returned by SageMath are:
#        u = 1/10*x - 1/5
#        v = -1/10*x^2 + 3/5


# OUTPUT (FINITE FIELD Z/7Z)
#polynomials in input:
#        p = x^3 + 2*x^2 + 6*x + 5
#        q = x^2 + 1
#base field: Finite Field of size 7
#---------------------
#**their GCD via the Euclidean algorithm:
#        5
#---------------------
#Running Extended Euclidean algorithm...
#Step 1:
#        u = 0,
#        v = 1
#Step 2:
#        u = 1,
#        v = 6*x + 5
#Step 3:
#        u = 4*x + 6,
#        v = 3*x^2 + 3
#---------------------
#**their GCD g and cofactors (u,v) such that u p + v q = g via the extended Euclidean algorithm:
#        u = 4*x + 6
#        v = 3*x^2 + 3
#        g = 5
#---------------------
#**GCD returned by SageMath is ``normalized'': gcd(p,q) = 1
#---------------------
#furthermore, cofactors returned by SageMath are:
#        u = 5*x + 4
#        v = 2*x^2 + 2

