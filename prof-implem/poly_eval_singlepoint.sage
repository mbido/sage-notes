# Evaluating a univariate polynomial at a single point

def eval_verynaive(pol, point):
    ev = 0
    for k in range(len(pol)):
        ev = ev + pol[k] * point**k
    return ev

def eval_naive(pol, point):
    ev = 0
    ptpow = 1 # at iteration k, ptpow = point**k
    for k in range(len(pol)):
        ev = ev + pol[k] * ptpow
        ptpow = point * ptpow 
    return ev

def eval_horner(pol, point):
    ev = pol[len(pol)-1]
    for k in range(len(pol)-2,-1,-1):
        ev = ev * point + pol[k]
    return ev

def eval_horner_memfriendly(pol, point):
    ev = 0
    for coeff in reversed(pol):
        ev = ev * point + coeff
    return ev

# you can try several fields; conclusions are not always similar:
#field = GF(97) # -> very small prime
field = GF(1125899906842679) # -> 51 bit prime
#field = GF(1267650600228229401496703205653) # -> 101 bit prime
#field = GF(3273390607896141870013189696827599152216642046043064789483291368096133796404674554883270092325904157150886684127560071009217256545885393053328527589431) # -> 501 bit prime
ring.<x> = field[]

# some basic benchmarking
print("prec\tverynaive\tnaive\t\thorner\t\thorner_mem\tsage")
for deg in (2**k for k in range(2,20)):
    pol = ring.random_element(degree=deg)
    lpol = pol.list()
    point = field.random_element()
    t_verynaive = timeit('eval_verynaive(lpol,point)',seconds=True,number=3,repeat=3)
    t_naive = timeit('eval_naive(lpol,point)',seconds=True,number=3,repeat=3)
    t_horner = timeit('eval_horner(lpol,point)',seconds=True,number=3,repeat=3)
    t_horner_memfriendly = timeit('eval_horner_memfriendly(lpol,point)',seconds=True,number=3,repeat=3)
    t_sage = timeit('pol(point)',seconds=True,number=3,repeat=3)
    print(f"{deg}\t{t_verynaive:.9f}\t{t_naive:.9f}\t{t_horner:.9f}\t{t_horner_memfriendly:.9f}\t{t_sage:.9f}")

