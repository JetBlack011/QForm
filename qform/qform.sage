from functools import reduce
import numpy as np
import random

# For memoizing the elementary symmetric matrices
ESM_dict = {}

def _J(g):
    return block_matrix([[0, matrix.identity(g)],[-matrix.identity(g), 0]])

def _Sp(n):
    return Sp(n, ZZ, invariant_form=_J(ZZ(n / 2)))

def cuts_to_intersection(g, cuts):
    """
    Given representatives of classes in H_1(S_g) expressed in terms of the
    standard generators (thought of as simple closed curves on S_g), compute
    their pairwise intersections.
    :param int g: The genus of this multisection
    :param list[matrix]: A list of the cut systems in terms of the standard
           generators of H_1(S_g)
    :return list[matrix]: The pairwise intersection matrices that fit into
           the whole intersection matrix as its blocks. That is, if we have
           cut systems a, b, and c, then the return value is [(a,b), (a,c), (b,c)],
           where (a,b) is the matrix of intersection values for the a and b
           cut systems.
    """
    J = _J(g)
    I = matrix(len(cuts) * 2)

    for i in range(I.nrows()):
        for j in range(I.ncols()):
            I[i,j] = (matrix(cuts[i//g][i%g]) * J * cuts[j//g][j%g])[0]
    
    mats = [I[i:i+g,j:j+g]
            for i in range(0, I.nrows() - g, g)
            for j in range(i + g, I.ncols(), g)]
    
    return mats

def _ind(y, x, l, n):
    # Translate from an index into a skew-symmetric matrix to the block
    # containing the respective element.
    return l - (n - y - 1) * (n - y) / 2 + x - y - 1

def _count(n):
    # Go from the dimension of a skew-symmetric matrix to the number of
    # its upper-triangular entries.
    return n * (n + 1) / 2

def _uncount(l):
    # Go from the number of upper-triangular entries of a skew-symmetric
    # matrix to its dimension.
    # Turns out this is the right expression o.O
    return (-1 + sqrt(1 + 8 * float(l))) / 2

def _build_intersection_matrix(g, mats):
    """
    Go from len(mats) blocks to the n = _uncount(len(mats) + 1 dimensional
    intersection matrix. This just takes the blocks of the intersection matrix
    and adds in the skew-symmetric side.
    :param int g: The genus of this multisection
    :param list[list[int]]: The intersection data of each cut system,
           as in intersection_form()
    :return (list[list[int]], list[int], int, list[int]): In order:
           The blocks of the intersection matrix, for ease of computation.
           The intersection matrix itself --- curve i and curve j intersect exactly
           I_{i,j} times (signed) for 0 <= i,j < n.
           The dimension of the intersection matrix.
           The basis of H_1(S_g) in terms of the standard alpha, beta cut systems
           as a 2x2 matrix.
    """
    blocks = []
    n = int(_uncount(len(mats))) + 1

    for i in range(n):
        blocks.append([])
        for j in range(n):
            if i == j:
                blocks[i].append(matrix(g))
            elif i > j:
                blocks[i].append(-mats[_ind(j, i, len(mats), n)].transpose())
            else:
                blocks[i].append(mats[_ind(i, j, len(mats), n)])

    base_blocks = []

    for i in range(2):
        base_blocks.append([])
        for j in range(2):
            base_blocks[i].append(blocks[i][j])

    I = block_matrix(blocks)
    P = block_matrix(base_blocks)

    assert I.is_skew_symmetric(), "Intersection matrix is not skew-symmetric"

    return (blocks,I,n,P)

def _all_integers(v):
    return all(list(map(lambda x: x.is_integer(), v)))

def _is_unimodular(M):
    return M.determinant() == 1 or M.determinant() == -1

def _build_kernel_generators(g, blocks, P, n):
    """
    Compute the generators of the kernel of 
    """
    ker = np.zeros((g * (n - 2), g * n))

    for i in range(2, n):
        for j in range(g):
            coeffs = P.inverse() * vector(QQ, [blocks[i][k][j][l]
                for k in range(2) for l in range(g)])
            a = 1 if _all_integers(coeffs) else 1 / reduce(gcd, coeffs) 
            coeffs = a * coeffs

            assert _all_integers(coeffs), \
                    "Kernel generators cannot be given in integer coefficients"

            e = g * (i - 2) + j
            ker[e,:g*2] = coeffs
            ker[e,e+g*2] = a
            ker[e] = list(map(int, ker[e]))

    return ker

def _compute_q_form(I, ker):
    """
    Given a basis for H_2(M) where M is the 4-manifold got from this multisection
    and the intersection data, compute the intersection pairing.
    :return list[int]: The matrix representation of the intersection form with respect
            to the given basis.
    """
    Q = matrix(len(ker))

    for i,e1 in enumerate(ker):
        for j,e2 in enumerate(ker):
            s = 0
            for k in range(len(e1)):
                for l in range(k + 1, len(e2)):
                    s += (e1[k] * e2[l]) * I[k][l]
            Q[i,j] = s

    assert Q.is_symmetric(), "Intersection form must be symmetric"

    return Q

def intersection_form(g, mats):
    """
    Solve for the intersection form given a multisection
    :param int g: Genus of the multisection diagram
    :param list[list[int]] mats: Blocks of the intersection matrix, in the
            form e.g. (a, b), (a, c), (b, c) where (a, b) looks like
             [[<a_1, b_1>, ..., <a_1, b_g>],
              [     ⁞    ,  ⋱ ,     ⁞     ],
              [<a_g, b_1>, ..., <a_g, b_g>]]
            That is, these are the pairwise intersection numbers between each
            cut system. (a, b) is the intersection matrix of the alpha and
            beta cut systems, and these blocks are joined together skew-symmetrically
            in the whole intersection matrix. This way, data entry doesn't have
            repeated (although skew-symmetric) values.
    :return: Returns the intersection form
    """
    mats = list(map(matrix, mats))

    # Standard assertions
    assert g >= 0, "Genus must be >= 0"
    for i,mat in enumerate(mats):
        assert mat.ncols() == mat.nrows() == g, \
                f"The block\n{mat}\nshould be {g}x{g}, but is actually " \
                f"{mat.nrows()}x{mat.ncols()}"
    assert _uncount(len(mats)).is_integer(), \
            "Intersection blocks don't fit into skew-symmetric matrix"

    blocks,I,n,P = _build_intersection_matrix(g, mats)
    #print(f'Intersection matrix:\n{I}\n')
    ker = _build_kernel_generators(g, blocks, P, n)
    #print(f'Kernel generators:\n{matrix(ZZ, ker)}\n')

    return _compute_q_form(I, ker)

def signature(B):
    """
    Find the signature of a bilinear form given its matrix representation.
    :param matrix: Matrix representation of a bilinear form B
    :return (int, int, int): The signature of B
    """
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.signature_vector()

def definiteness(B):
    """
    Find the definiteness of a bilinear form given its matrix representation.
    :param matrix: Matrix representation of a bilinear form B
    :return str: The definiteness of B, either 'pos_def', 'zero', 'degenerate',
           'indefinite', or 'neg_def'
    """
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.compute_definiteness_string_by_determinants()

def parity(B):
    """
    Find the pairty of a bilinear form given its matrix representation.
    :param matrix: Matrix representation of a bilinear form B
    :return str: Either 'even' or 'odd'
    """
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.parity()

def _ESMs(g):
    """
    Generate the elementary symmetric matrices (ESMs) which form the
    generators for Sp(2g, Z).
    :param int g: The genus of this multisection
    :return list(matrix): The ESMs for Sp(2g, Z)
    """
    global ESM_dict

    # Memoize the generators of Sp(2g, Z)
    if g in ESM_dict:
        return ESM_dict[g]
    
    ESMs = []
    I2g = matrix.identity(2 * g)
    X = matrix.identity(2 * g)
    for i in range(1, 2 * g - 1, 2):
        X.swap_rows(i, i + 1)

    def e(i, j):
        E = matrix(2 * g)
        E[i,j] = 1
        return E

    # Element of S_2g transposing 2i and 2i - 1 for all 1 <= i <= g
    def sigma(i):
        return i + 1 if i % 2 == 0 else i - 1

    for i in range(2 * g):
        for j in range(2 * g):
            if i == j:
                continue
            if i == sigma(j):
                ESM = I2g + e(i, j)
            else:
                ESM = I2g + e(i, j) - (-1)^(i + j) * e(sigma(j), sigma(i))

            ESMs.append(X.transpose() * ESM * X)
            #assert ESMs[-1] in _Sp(2 * g), f"ESM not in Sp({2 * g}, Z)"

    ESM_dict[g] = ESMs
    return ESMs

def random_symplectic_matrix(g, p):
    """
    Generate a random element of Sp(2g, Z). THIS IS THE OPPOSITE OF UNIFORMLY
    DISTRIBUTED! All we do here is multiply some random ESMs together to
    get some ""random"" element of Sp(2g, Z). These act on the basis of
    H_1(S_g) while respecting a symplectic form thereon, say an
    intersection form.
    :param int g: The genus of this multisection
    :param int p: The limit of how many ESMs to multiply together
           (p for "product")
    """
    # transvection = matrix([[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,0,0,1]])
    # rotation     = matrix([[0,0,-1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]])
    # mix          = matrix([[1,0,0,0],[0,1,0,0],[0,-1,1,0],[-1,0,0,1]])
    # swap         = matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    # K = matrix([[1,0,0,0],[1,-1,0,0],[0,0,1,1],[0,0,0,-1]])
    # L = matrix([[0,0,-1,0],[0,0,0,-1],[1,0,1,0],[0,1,0,0]])

    ESMs = _ESMs(g)
    P = [random.choice(ESMs) for _ in range(random.randint(1, p))]
    M = reduce(lambda a, x: a * x, P)

    assert M in _Sp(2 * g), f"Generated matrix is not in Sp({2 * g}, Z)"

    return (M, P)

def random_unimodular_matrix(n):
    """
    Generate a random matrix with determinant +1 or -1.
    :param int n: The dimension of the desired matrix
    """
    matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, n)
    return sage.matrix.constructor.random_unimodular_matrix(matrix_space,
            upper_bound=10)

def random_diagram(g, n, limit=100):
    """
    Generate a random genus 2 multisection diagram with n cut systems,
    the first two being the standard alpha and beta cut systems.
    :param int g: The genus of the cut system
    :param int n: The number of cut systems
    :return list[list[vector]: The cut systems expressed in the standard
           basis. The first two cut systems will always be the identity.
    """
    assert n >= 2, "n should be >= 2"

    x = 10
    B = matrix.identity(2 * g) # The standard basis
    
    """
    What we do here is a bit awkward. The idea is to generate random symplectic
    maps which land us in new valid cut systems when we act on the previous one.
    However, it can be harder to generate valid maps for larger values (one
    potential improvement to this design would be to check if certain cuts
    can be reduced to known ones, but this is very hard, as it probably
    involves understanding the resulting 4-manifold). In this case, it might take
    a long while, or even forever, for this process to finish.
    
    So, we opt to use a limiting variable. If we get stuck on too many iterations,
    we start over. Hence the loop. This has been tuned through many run-throughs
    with reasonable values, but will not work for every application, especially
    larger values of n and g.
    """
    while True:
        limit_flag = True
        cuts = [[B.column(i + j) for j in range(g)] for i in range(0, B.ncols(), g)]

        for i in range(len(cuts)):
            cuts[i] = list(map(vector, cuts[i]))

        while len(cuts) < n - 1: # Keep going until we have enough cuts
            l = 0
            while l < limit:
                l += 1
                symp_map = random_symplectic_matrix(g, x)[0]
                image = [symp_map * cuts[-1][i] for i in range(g)]
                # det2 = matrix(image + cuts[0]).determinant()
                # #det = cuts_to_intersection(g, [cuts[-1], image])[0].determinant()
                Q = intersection_form(g, cuts_to_intersection(g, cuts + [image]))
                diagonal = Q.numpy().diagonal()
                b = 2 * (len(cuts) - 1)
                block = Q[:b,:b]

                if _is_unimodular(matrix(cuts[-1] + image)) and all(map(is_even, diagonal)) and not _is_unimodular(block):
                    # We found a valid one! Time to move on.
                    break
            limit_flag = limit_flag and l < limit
            cuts.append(image)

        # The condition that neighboring cut systems form Lagrangian subspaces with respect
        # to the intersection form wraps around. This second loops just checks that the last
        # cut system we generate agrees with the first.
        l = 0
        while l < limit:
            l += 1
            symp_map = random_symplectic_matrix(g, x)[0]
            image = [symp_map * cuts[-1][i] for i in range(g)]
            #det1 = cuts_to_intersection(g, [cuts[-1], image])[0].determinant()
            #det2 = cuts_to_intersection(g, [image, cuts[0]])[0].determinant()
            Q = intersection_form(g, cuts_to_intersection(g, cuts + [image]))
            diagonal = Q.numpy().diagonal()
            # d = definiteness(Q)
            #is_connect_sum = False
            blocks = [Q[:i,:i] for i in range(2, Q.nrows() - 2, 2)]

            #for A in blocks:
            #    if A.determinant() == 1 or A.determinant() == -1:
            #        is_connect_sum = True

            if _is_unimodular(matrix(cuts[-1] + image)) and _is_unimodular(matrix(image + cuts[0])) and all(map(is_even, diagonal)):
                break
        cuts.append(image)

        limit_flag = limit_flag and l < limit
        if limit_flag:
            break

    return cuts

if __name__ == "__main__":
    mats = [
        [[1, 0], [0, 1]],
        [[1, 0], [0, 1]],
        [[-1, 0], [0, -1]],
    ]

    #intersection_form(2, mats)
    _compute_q_form([], [])
