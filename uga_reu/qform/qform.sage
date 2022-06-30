from functools import reduce
import numpy as np
import random

I4 = matrix.identity(4)

def _J(g):
    return block_matrix([[0, matrix.identity(g)],[-matrix.identity(g), 0]])

def _Sp(n):
    return Sp(n, ZZ, invariant_form=_J(ZZ(n / 2)))

def _ind(y, x, l, n):
    return l - (n - y - 1) * (n - y) / 2 + x - y - 1

def _count(n):
    return n * (n + 1) / 2

def _uncount(l):
    # Turns out this is the right expression o.O
    return (-1 + sqrt(1 + 8 * float(l))) / 2

def _build_intersection_matrix(mats):
    # Generate our base intersection matrix
    blocks = []
    n = int(_uncount(len(mats))) + 1

    for i in range(n):
        blocks.append([])
        for j in range(n):
            if i == j:
                blocks[i].append(matrix(_g))
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
    B = block_matrix(base_blocks)

    #print(f'Intersection matrix:\n{I}\n')

    assert I.is_skew_symmetric(), "Intersection matrix is not skew-symmetric"

    #print(f'{A}\n')

    return (blocks,I,n,B)

def _all_integers(v):
    return all(list(map(lambda x: x.is_integer(), v)))

def _build_kernel_generators(blocks, B, n):
    # Build the kernel of our map
    ker = np.zeros((_g * (n - 2), _g * n))

    for i in range(2, n):
        for j in range(_g):
            coeffs = B.inverse() * vector(QQ, [blocks[i][k][j][l]
                for k in range(2) for l in range(_g)])
            a = 1 if _all_integers(coeffs) else 1 / reduce(gcd, coeffs) 
            coeffs = a * coeffs

            assert _all_integers(coeffs), \
                    "Kernel generators cannot be given in integer coefficients"

            e = _g * (i - 2) + j
            ker[e,:_g * 2] = coeffs
            ker[e,e + _g * 2] = a
            ker[e] = list(map(int, ker[e]))

    return ker

def _compute_q_form(I, ker):
    # Compute the intersection form
    Q = matrix(len(ker))

    for i,e1 in enumerate(ker):
        for j,e2 in enumerate(ker):
            s = 0
            for k in range(len(e1)):
                for l in range(k + 1, len(e2)):
                    s += (e1[k] * e2[l]) * I[k][l]
            Q[i,j] = s

    assert Q.is_symmetric(), \
            "Intersection form is not symmetric, this likely means the kernel "\
            "was generated incorrectly"

    return Q

def intersection_data(g, cuts):
    """
    Given some elements of (?), find the intersection data using
    basis B.
    :param T: Basis of Lagrangian decomposition in the form of a symplectic
    matrix
    :param B: Initial Lagrangian decomposition, by default this is the 4x4
    identity matrix
    :param M: The (symplectic) intersection pairing for the given basis, by
    default this is J
    :return: A block matrix of intersection data
    """
    B = matrix.identity(2 * g)
    J = _J(g)

    dim = len(cuts) - 1
    I = matrix(len(cuts) * 2)

    for i in range(I.nrows()):
        for j in range(I.ncols()):
            I[i,j] = (matrix(cuts[i//g][i%g]) * J * cuts[j//g][j%g])[0]
    
    mats = [I[i:i+g,j:j+g]
            for i in range(0, I.nrows() - g, g)
            for j in range(i + g, I.ncols(), g)]
    
    return mats

def intersection_form(g, mats):
    """
    Solve for the intersection form given a multisection
    :param int g: Genus of the multisection diagram
    :param list[list[int]] mats: Blocks of the intersection matrix, in the \
            form e.g. (a, b), (a, c), (b, c) where (a, b) looks like \
             [[<a_1, b_1>, ..., <a_1, b_g>],
              [     ⁞    ,  ⋱ ,     ⁞     ],
              [<a_g, b_1>, ..., <a_g, b_g>]]
    :return: Returns the intersection form
    """
    global _g
    _g = g
    mats = list(map(matrix, mats))

    # Standard assertions
    assert _g >= 0, "Genus must be >= 0"
    for i,mat in enumerate(mats):
        assert mat.ncols() == mat.nrows() == _g, \
                f"The block\n{mat}\nshould be {_g}x{_g}, but is actually " \
                f"{mat.nrows()}x{mat.ncols()}"
    assert _uncount(len(mats)).is_integer(), \
            "Intersection blocks don't fit into skew-symmetric matrix"

    blocks,I,n,B = _build_intersection_matrix(mats)
    #print(f'Intersection matrix:\n{I}\n')
    ker = _build_kernel_generators(blocks, B, n)
    #print(f'Kernel generators:\n{matrix(ZZ, ker)}\n')

    return _compute_q_form(I, ker)

def signature(B):
    """
    Find the signature of a bilinear form given its matrix representation.
    :param matrix: Matrix representation of a bilinear form B
    :return: Returns the signature of B
    """
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.signature_vector()

def definiteness(B):
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.compute_definiteness_string_by_determinants()

def parity(B):
    Q = QuadraticForm(ZZ, B + B.transpose())
    return Q.parity()

def random_symplectic_matrix(n):
    # transvection = matrix([[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,0,0,1]])
    # rotation     = matrix([[0,0,-1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]])
    # mix          = matrix([[1,0,0,0],[0,1,0,0],[0,-1,1,0],[-1,0,0,1]])
    # swap         = matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    K = matrix([[1,0,0,0],[1,-1,0,0],[0,0,1,1],[0,0,0,-1]])
    L = matrix([[0,0,-1,0],[0,0,0,-1],[1,0,1,0],[0,1,0,0]])

    P = [random.choice([K, L]) for _ in range(random.randint(1,100))]
    M = reduce(lambda a, x: a * x, P)

    assert M in _Sp(4), "Generated matrix is not in Sp(Z, 4)"

    return (M, P)

def random_unimodular_matrix(n):
    matrix_space = sage.matrix.matrix_space.MatrixSpace(ZZ, n)
    return sage.matrix.constructor.random_unimodular_matrix(matrix_space,
            upper_bound=3)

def random_diagram(n, g):
    """
    Generate a random genus 2 multisection diagram with n + 2 cut systems,
    the first two being the standard alpha and beta cut systems.
    """
    global _g
    _g = g
    B = matrix.identity(2 * _g)

    limit = 100

    while True:
        limit_flag = True
        cuts = [[B.column(i + j) for j in range(_g)] for i in range(0, B.ncols(), _g)]

        for _ in range(n - 1):
            l = 0
            while True and l < limit:
                l += 1
                symp_map = random_symplectic_matrix(2 * _g)[0]
                image = [symp_map * cuts[-1][i] for i in range(_g)]
                det = matrix(cuts[-1] + image).determinant()
                if det == 1 or det == -1:
                    break
            limit_flag = limit_flag and l < limit
            cuts.append(image)

        l = 0
        while True and l < limit:
            l += 1
            symp_map = random_symplectic_matrix(2 * _g)[0]
            image = [symp_map * cuts[-1][i] for i in range(_g)]
            det1 = matrix(cuts[-1] + image).determinant()
            det2 = matrix(image + cuts[0]).determinant()

            if (det1 == 1 or det1 == -1) and (det2 == 1 or det2 == -1):
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
