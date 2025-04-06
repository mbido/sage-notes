import time

load("subproduct_tree.sage")

# returns n distincts elements in field
# requires field to be large enough
def distincts_elements(n, field):
    l = []
    while len(l) < n :
        x = field.random_element()
        if x not in l:
            l.append(x)
    return l

# returns n nonzero random elements in field
# requires field to be large enough
def nonzero_elements(n, field):
    l = []
    while len(l) < n :
        x = field.random_element()
        if x != 0:
            l.append(x)
    return l

def is_leaf(spt):
    return not(spt.left_tree) and not(spt.right_tree)

# leaves of the remainder tree constructed from the root of spt and its derivative
# Note: we don't care about building the tree, we just want to retrieve the leaves
def leaves_rmt(dpol, spt):
    root = dpol.quo_rem(spt.data)[1]
    if (root.degree() <= 0) and (is_leaf(spt)):
        # base case: leaf
        return [root.constant_coefficient()]

    left_rem  = leaves_rmt(root, spt.left_tree)
    right_rem = leaves_rmt(root, spt.right_tree)

    return left_rem + right_rem

# construct partial sum tree from list of leaves
def construct_pst(trees):
    n = len(trees)
    # base cases, n is 0 or 1
    if n == 0:
        raise ValueError("list of trees must have positive length")
    if n == 1:
        return trees[0]

    # recursion, n >= 2
    # construct list of new trees, one level up
    newtrees = []
    for i in range(0,n-1,2):
        newroot = trees[i].data + trees[i+1].data
        newtrees.append(BTree(newroot, trees[i], trees[i+1]))
    # if n was odd, add dangling tree
    if (n % 2) == 1:
        newtrees.append(trees[n-1])

    # call recursively
    return construct_pst(newtrees)

# construct partial sum tree from x_pts, y_pts, leaves of remainder tree
def partialsum_tree(x_pts, y_pts, rmt_leaves, pring=None):
    field = x_pts[0].parent()
    if pring == None:
        pring = PolynomialRing(field, 'x')
    x = pring.gen()

    # build leaves of partial sum tree
    pst_leaves = [BTree((y/l) / (x-pt)) for pt, y, l in zip(x_pts, y_pts, rmt_leaves)]

    # construct and return whole subproduct tree
    return construct_pst(pst_leaves)


def interpolation(x_pts, y_pts):
    # step 1: compute spt
    spt = subproduct_tree(x_pts)
    #print(spt)

    # step 2: compute remainder tree
    da = derivative(spt.data)
    leaves = leaves_rmt(da, spt)

    # step 3: compute partial sums tree
    pst = partialsum_tree(x_pts, y_pts, leaves)

    # return numerator of the root of pst
    return pst.data.numerator()


########### TESTS COURS : EX VI.9 ############

if False:
    field = GF(7)
    pring.<x> = field[]

    x_pts = [field(e) for e in [2,5,3,4]]
    y_pts = [field(e) for e in [1,2,3,6]]
    points = [(x,y) for x,y in zip(x_pts,y_pts)]

    ts = time.perf_counter()
    P_sage = pring.lagrange_polynomial(points)
    t_sage = time.perf_counter() - ts
    print(f'sage={t_sage:.10f}; {P_sage}')

    ts = time.perf_counter()
    P_fast = interpolation(x_pts,y_pts)
    t_fast = time.perf_counter() - ts
    print(f'fast={t_fast:.10f}; {P_fast}')

############# BENCHMARK ##############

if True: 
    print(f'nb\tsage\t\tfast\t\tratio:s-f')
    for k in range(2, 14):
        n = 1 << k
        
        field = GF(next_prime(n))
        pring.<x> = field[]

        x_pts = distincts_elements(n, field)
        y_pts = nonzero_elements(n, field) #[field.random_element() for _ in range(n)] 
        points =  [(x,y) for x,y in zip(x_pts,y_pts)]

        ts = time.perf_counter()
        sage = pring.lagrange_polynomial(points)
        t_sage = time.perf_counter() - ts

        ts = time.perf_counter()
        fast = interpolation(x_pts, y_pts)
        t_fast = time.perf_counter() - ts

        if sage != fast:
            # print parameters for manual exec
            print("Something went wrong")
            print(x_pts, y_pts)
            print(sage)
            print(fast)
            break
        
        print(f'{n}\t{t_sage:.10f}\t{t_fast:.10f}\t{t_sage/t_fast:.3f}')

    reset()