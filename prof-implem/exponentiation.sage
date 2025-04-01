############################
#  functions for powering  #
############################

# they all have the following specification:
""" computes the power a**n
:a: element of a field
:n: nonnegative integer
:returns: a**n
"""


def powering_naive(a, n):
    # idea: a**n = a * a * ... * a, with a appearing n times
    b = 1
    for i in range(n):
        b = a * b
    return b

def powering_naive_rec(a, n):
    # idea: a**n = a**(n-1) * a, used recursively
    if n == 0:
        return 1
    b = powering_naive_rec(a, n-1)
    return a * b

def powering_fast(a, n):
    # idea: if n = 2*nn, a**n == (a**nn)**2
    #       if n = 2*nn + 1, a**n == a * (a**nn)**2
    if n == 0:
        return 1
    nn = n // 2
    b = powering_fast(a, nn)
    if is_even(n):
        return b * b
    else:
        return a * b * b

    


##############
#  examples  #
##############

import time 

# naive_rec will fail beyond ~2**10 due to maximum recursion depth exceeded
# set to False
include_naive_rec = False

field = FiniteField(17)
a = field.random_element()
max_k = 30 if not include_naive_rec else 11
exps = [2**k for k in range(max_k)]

times_naive = []
times_naive_rec = []
times_fast = []

for n in exps:
    t_start = time.time()
    b = powering_naive(a, n)
    t_end = time.time()
    times_naive.append(t_end - t_start)

    if include_naive_rec:
        t_start = time.time()
        b = powering_naive_rec(a, n)
        t_end = time.time()
        times_naive_rec.append(t_end - t_start)
    else:
        times_naive_rec.append(-1.)

    t_start = time.time()
    b = powering_fast(a, n)
    t_end = time.time()
    times_fast.append(t_end - t_start)

k = 0
print("k\tn == 2**k\tnaive\t\tnaive_rec\tfast")
for n in exps:
    print(f"{k}\t{n:9n}\t{times_naive[k]:.2e}\t{times_naive_rec[k]:.2e}\t{times_fast[k]:.2e}")
    k += 1

## sage: %runfile exponentiation.sage
## k       n == 2**k       naive           naive_rec       fast
## 0               1       1.07e-05        -1.00e+0        2.72e-05
## 1               2       3.34e-06        -1.00e+0        1.07e-05
## 2               4       2.62e-06        -1.00e+0        1.05e-05
## 3               8       3.58e-06        -1.00e+0        1.31e-05
## 4              16       3.81e-06        -1.00e+0        1.43e-05
## 5              32       2.38e-06        -1.00e+0        2.74e-05
## 6              64       3.58e-06        -1.00e+0        1.91e-05
## 7             128       7.15e-06        -1.00e+0        2.48e-05
## 8             256       9.54e-06        -1.00e+0        2.50e-05
## 9             512       2.17e-05        -1.00e+0        2.60e-05
## 10           1024       5.48e-05        -1.00e+0        3.10e-05
## 11           2048       9.97e-05        -1.00e+0        3.34e-05
## 12           4096       2.03e-04        -1.00e+0        3.29e-05
## 13           8192       3.71e-04        -1.00e+0        3.31e-05
## 14          16384       6.54e-04        -1.00e+0        3.19e-05
## 15          32768       1.22e-03        -1.00e+0        4.15e-05
## 16          65536       2.01e-03        -1.00e+0        2.96e-05
## 17         131072       4.14e-03        -1.00e+0        3.31e-05
## 18         262144       8.08e-03        -1.00e+0        3.22e-05
## 19         524288       1.67e-02        -1.00e+0        6.08e-05
## 20        1048576       3.29e-02        -1.00e+0        4.39e-05
## 21        2097152       6.61e-02        -1.00e+0        4.98e-05
## 22        4194304       1.31e-01        -1.00e+0        5.36e-05
## 23        8388608       2.62e-01        -1.00e+0        7.84e-05
## 24       16777216       5.37e-01        -1.00e+0        1.54e-04
## 25       33554432       1.03e+00        -1.00e+0        9.78e-05
## 26       67108864       2.11e+00        -1.00e+0        1.06e-04
## 27      134217728       4.19e+00        -1.00e+0        1.51e-04
## 28      268435456       7.70e+00        -1.00e+0        1.52e-04
## 29      536870912       1.69e+01        -1.00e+0        1.07e-04

