import textwrap

# basic class for binary tree
class BTree:
    def __init__(self, data, left=None, right=None):
        self.data = data        # root node data
        self.left_tree = left   # left sub-tree
        self.right_tree = right # right sub-tree
    def __repr__(self):
        indent_str = "--+--"
        data = self.data.__str__()
        leftstr = textwrap.indent("left: "+str(self.left_tree), indent_str)
        rightstr = textwrap.indent("right: "+str(self.right_tree), indent_str)
        return "BTree with root " + data + "\n" + leftstr + "\n" + rightstr
    def __str__(self):
        indent_str = "--+--"
        data = self.data.__str__()
        leftstr = textwrap.indent("left: "+str(self.left_tree), indent_str)
        rightstr = textwrap.indent("right: "+str(self.right_tree), indent_str)
        return "BTree with root " + data + "\n" + leftstr + "\n" + rightstr

## Example:
field = GF(97)
pring.<x> = field[]

# let's build a subproduct tree with n=4 points
n = 4
points = [field.random_element() for i in range(n)]

# first some examples to show what it looks like:
print(f"Tree which is a leaf (no left/right trees):\n{BTree(x-points[0])}\n")
print(f"Another tree which is also a leaf (no left/right trees):\n{BTree((x-points[1])*(x-points[0]))}\n")
print(f"Tree which is not a leaf, height 1:\n{BTree((x-points[1])*(x-points[0]), BTree(x-points[1]), BTree(x-points[0]))}\n")

# now let's go
leaves = [BTree(x-pt) for pt in points]
level1 = []
level1.append(BTree(leaves[0].data * leaves[1].data, leaves[0], leaves[1]))
level1.append(BTree(leaves[2].data * leaves[3].data, leaves[2], leaves[3]))
root = BTree(level1[0].data*level1[1].data, level1[0], level1[1])

print(f"Subproduct tree of 4 points:\n{root}\n")


############
#  OUTPUT  #
############

# Tree which is a leaf (no left/right trees):
# BTree with root x + 16
# --+--left: None
# --+--right: None
# 
# Another tree which is a leaf (no left/right trees):
# BTree with root x^2 + 18*x + 32
# --+--left: None
# --+--right: None
# 
# Tree which is not a leaf, height 1:
# BTree with root x^2 + 18*x + 32
# --+--left: BTree with root x + 2
# --+----+--left: None
# --+----+--right: None
# --+--right: BTree with root x + 16
# --+----+--left: None
# --+----+--right: None
# 
# Subproduct tree of 4 points:
# BTree with root x^4 + 43*x^3 + 41*x^2 + 40*x + 50
# --+--left: BTree with root x^2 + 18*x + 32
# --+----+--left: BTree with root x + 16
# --+----+----+--left: None
# --+----+----+--right: None
# --+----+--right: BTree with root x + 2
# --+----+----+--left: None
# --+----+----+--right: None
# --+--right: BTree with root x^2 + 25*x + 44
# --+----+--left: BTree with root x + 92
# --+----+----+--left: None
# --+----+----+--right: None
# --+----+--right: BTree with root x + 30
# --+----+----+--left: None
# --+----+----+--right: None


