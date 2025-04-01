load("subproduct_tree.sage")

# Multipoint evaluation, recursion
# assumes the subproduct tree has been correctly computed from the points
# does not assume, but is better if pol already reduced modulo the root of tree
#  (otherwise wastes time with two reductions if degree of pol is larger than that of this root)
def multipoint_evaluation_rec(pol, tree):
    # base case: constant evaluates to itself
    if pol.degree() <= 0:
        return [pol.constant_coefficient()]
    # compute remainders w.r.t left and right sub-trees
    # (note these cannot be "None" since otherwise pol would have degree <= 0)
    rem_left = pol.quo_rem(tree.left_tree.data)[1]
    rem_right = pol.quo_rem(tree.right_tree.data)[1]
    # perform recursive calls
    mpe_remleft = multipoint_evaluation_rec(rem_left, tree.left_tree)
    mpe_remright = multipoint_evaluation_rec(rem_right, tree.right_tree)
    return mpe_remleft + mpe_remright

# Multipoint evaluation
# assumes the subproduct tree has been correctly computed from the points
# no constraint on the degree of pol
def multipoint_evaluation(pol, points):
    # build tree
    tree = subproduct_tree(points)
    # first reduce by root if necessary
    rem = pol.quo_rem(tree.data)[1]
    # perform evaluation and return only wanted evaluations
    return multipoint_evaluation_rec(pol, tree)

## Example:
field = GF(97)
pring.<x> = field[]

# create subproduct tree
n = 7
points = [field.random_element() for i in range(n)]
print(f"points: {points}\n")
tree = subproduct_tree(points)
print(f"Subproduct tree:\n{tree}\n")

# take random polynomial (note degree does not have to be exactly n)
pol = pring.random_element(degree=n+3)
print(f"polynomial: {pol}\n")

# let's go, and test
evals = multipoint_evaluation(pol, points)
print(f"evaluations: {evals}\n")

print(f"testing...", end="")
check = "correct"
for i in range(n):
    if evals[i] != pol(points[i]):
        check = "wrong"
print(f" {check}\n")


# basic benchmarking
print("nb pts\tbuild tree\ttotal mpe\tsage repeated evals")
lengths = sorted([1<<k for k in range(3,11)] + [(1<<k) + (1<<(k-1)) for k in range(3,11)])
for n in lengths:
    points = [field.random_element() for i in range(n)]
    pol = pring.random_element(degree=n-1)
    t_tree = timeit('tree = subproduct_tree(points)', seconds=True, number=10, repeat=10)
    t_mpe = timeit('evals = multipoint_evaluation(pol, points)', seconds=True, number=10, repeat=10)
    t_sage = timeit('evals = [pol(pt) for pt in points]', seconds=True, number=10, repeat=10)
    print(f"{n}\t{t_tree:.9f}\t{t_mpe:.9f}\t{t_sage:.9f}")

#lengths = sorted([1<<k for k in range(11,18)] + [(1<<k) + (1<<(k-1)) for k in range(11,18)])
#for n in lengths:
#    points = [field.random_element() for i in range(n)]
#    pol = pring.random_element(degree=n-1)
#    t_tree = timeit('tree = subproduct_tree(points)', seconds=True, number=1, repeat=1)
#    t_mpe = timeit('evals = multipoint_evaluation(pol, points)', seconds=True, number=1, repeat=1)
#    t_sage = timeit('evals = [pol(pt) for pt in points]', seconds=True, number=1, repeat=1)
#    print(f"{n}\t{t_tree:.9f}\t{t_mpe:.9f}\t{t_sage:.9f}")


#### OUTPUT ####

# points: [27, 86, 74, 13, 29, 61, 73]
# 
# Subproduct tree:
# BTree with root x^7 + 25*x^6 + 25*x^5 + 94*x^4 + 45*x^3 + 5*x^2 + 47*x + 15
# --+--left: BTree with root x^4 + 91*x^3 + 20*x^2 + 68*x + 48
# --+----+--left: BTree with root x^2 + 81*x + 91
# --+----+----+--left: BTree with root x + 70
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+----+----+--right: BTree with root x + 11
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+----+--right: BTree with root x^2 + 10*x + 89
# --+----+----+--left: BTree with root x + 23
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+----+----+--right: BTree with root x + 84
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+--right: BTree with root x^3 + 31*x^2 + 94*x + 67
# --+----+--left: BTree with root x^2 + 7*x + 23
# --+----+----+--left: BTree with root x + 68
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+----+----+--right: BTree with root x + 36
# --+----+----+----+--left: None
# --+----+----+----+--right: None
# --+----+--right: BTree with root x + 24
# --+----+----+--left: None
# --+----+----+--right: None
# 
# polynomial: 77*x^10 + 75*x^9 + 86*x^8 + 8*x^7 + 38*x^6 + 80*x^5 + 74*x^4 + 38*x^3 + 79*x^2 + 25*x + 8
# 
# evaluations: [8, 4, 74, 42, 23, 73, 8]
# 
# testing... correct
# 
# nb pts  build tree      total mpe       sage repeated evals
# 8       0.000018410     0.000033358     0.000003354
# 12      0.000025303     0.000046278     0.000005271
# 16      0.000031391     0.000060309     0.000007479
# 24      0.000044132     0.000086685     0.000012087
# 32      0.000056789     0.000111910     0.000017736
# 48      0.000082789     0.000165851     0.000032595
# 64      0.000108147     0.000222995     0.000049602
# 96      0.000165096     0.000343822     0.000094379
# 128     0.000223064     0.000459007     0.000151819
# 192     0.000335962     0.000696136     0.000311365
# 256     0.000444509     0.000951723     0.000521299
# 384     0.000661925     0.001438777     0.001096490
# 512     0.000881019     0.001992610     0.001872929
# 768     0.001329218     0.003081362     0.004083765
# 1024    0.001767169     0.004125975     0.007102948
# 1536    0.002682416     0.006341483     0.015757081
# 2048    0.003471554     0.008281988     0.028230878
# 3072    0.004911169     0.013287760     0.062925702
# 4096    0.006601597     0.017985363     0.110103917
# 6144    0.009942774     0.028404241     0.246030954
# 8192    0.013437557     0.038805128     0.433214495
# 12288   0.019753844     0.060389129     0.974039341
# 16384   0.028109526     0.083973146     1.724792804
# 24576   0.043137406     0.132242080     3.873701516
# 32768   0.057423474     0.180579382     6.886955834
# 49152   0.087884333     0.285646693     15.512226312
# 65536   0.118288267     0.391896166     27.570235391
# 98304   0.183853452     0.616215112     61.841154749
# 131072  0.244853883     0.836325914     109.980323771
# 196608  0.377587017     1.321852976     248.796647431

