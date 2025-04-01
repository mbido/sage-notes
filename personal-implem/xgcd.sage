
def xgcd(a, b):
    # a and b a two polynomials
    # xgcd(a, b) returns (gcd, u, v) such as u*a + v*b = gcd
    
    ring = a.parent() # in most cases it should be a field

    if b == 0:
        return a, ring.one(), ring.zero()
    s0, s1 = ring.one(), ring.zero()
    t0, t1 = ring.zero(), ring.one()
    while b:
        q, r = a.quo_rem(b)
        a, b = b, r
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    return a, s0, t0
