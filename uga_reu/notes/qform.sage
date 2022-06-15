import numpy as np

# Genus
g = 2

# Intersection matrices of the form
# [[<a_1, b_1>, ..., <a_1, b_g>],
#  [     ⁞    ,  ⋱ ,     ⁞     ],
#  [<a_g, b_1>, ..., <a_g, b_g>]]
mats = [
    [[1,2],[2,1]],   # a,b
    [[3,4],[-1,-2]], # a,c
    [[-1,-2],[3,4]], # a,d
    [[4,5],[6,7]],   # b,c
    [[0,1],[1,0]],   # b,d
    [[1,5],[-2,-3]], # c,d
]

#mats = [
#    [[1,2,3],[4,5,6],[2,1,3]],       # a,b
#    [[3,4,5],[6,7,8],[-1,-2,-3]],    # a,c
#    [[-1,-2,-3],[-4,-5,-6],[3,4,5]], # a,d
#    [[4,5,6],[7,8,9],[6,7,8]],       # b,c
#    [[0,1,2],[3,4,5],[1,0,-1]],      # b,d
#    [[1,5,9],[2,3,4],[-2,-3,-4]],    # c,d
#]

# mats = [
#     [[1, 0], [0, 1]], # a, b
#     [[1, 0], [0, 1]], # a, c
#     [[0, 1], [1, 0]], # b, c
# ]

def ind(y, x):
    return len(mats) - (n - y - 1) * (n - y) / 2 + x - y - 1

def uncount(l):
    # Turns out this is the right expression o.O
    return (-1 + sqrt(1 + 8 * float(l))) / 2

mats = list(map(matrix, mats))

# Standard assertions
assert g >= 0, "Invalid genus"
for i,mat in enumerate(mats):
    assert mat.ncols() == mat.nrows() == g, \
            f"The block\n{mat}\nshould be {g}x{g}, but is actually " \
            f"{mat.nrows()}x{mat.ncols()}"
assert uncount(len(mats)).is_integer(), \
        "Intersection blocks don't fit into skew-symmetric matrix"

# Generate our base intersection matrix
blocks = []
n = int(uncount(len(mats))) + 1
for i in range(n):
    blocks.append([])
    for j in range(n):
        if i == j:
            blocks[i].append(matrix(g))
        elif i > j:
            blocks[i].append(-mats[ind(j, i)].transpose())
        else:
            blocks[i].append(mats[ind(i, j)])

base_blocks = []
for i in range(g):
    base_blocks.append([])
    for j in range(g):
        base_blocks[i].append(blocks[i][j])

I = block_matrix(blocks)
A = block_matrix(base_blocks)

print(f'Intersection matrix:\n{I}\n')

assert I.is_skew_symmetric(), "Intersection matrix is not skew-symmetric"

print(f'{A}\n')

# Build the kernel of our map
ker = [[] for _ in range(g * (n - g))]

for i in range(g, n):
    for j in range(g):
        coeffs = np.linalg.solve(
            A,
            [blocks[a][i][k][j] for a in range(g) for k in range(g)]
        ).tolist()
        e = g * (i - g) + j
        ker[e].extend(coeffs)
        ker[e].extend(np.zeros(g * (n - g)).tolist())
        ker[e][g*g+j] = -1
        ker[e] = list(map(lambda x: int(x) if x.is_integer() else x, ker[e]))

print(f'Kernel generators:\n{matrix(ker)}\n')

# Compute the intersection form
Q = []

for a in range(g):
    Q.append([])
    for b in range(g):
        s = 0
        for i in range(len(ker[a])):
            for j in range(i, len(ker[b])):
                s += (ker[a][i] * ker[b][j]) * I[i][j]
        Q[a].append(s)

print(f'Q =\n{matrix(Q)}')
