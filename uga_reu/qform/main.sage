import qform
import os

# Genus
g = 2

mats = [
    [[1,0],[0,1]],
    [[0,-1],[1,1]],
    [[1,0],[1,1]],
    [[-1,0],[-2,-1]],
    [[-1,0],[-1,-1]],
    [[2,1],[1,0]]
]

# mats = [
#     [[1,0],[0,1]],
#     [[-1,0],[0,-1]],
#     [[-1,-1],[0,-1]],
#     [[1,0],[0,-1]],
#     [[0,1],[1,0]],
#     [[1,1],[0,1]]
# ]

# mats = [
#     [[1,0],[0,1]],
#     [[1,0],[0,1]],
#     [[1,-1],[0,-1]],
#     [[-1,0],[0,-1]],
#     [[0,-1],[1,0]],
#     [[-1,0],[1,1]],
# ]

# mats = [
#     [[4,3],[1,1]],
#     [[4,3],[1,1]],
#     [[1,1],[1,0]],
# ]

# mats = [
#     [[1, 1], [2, 1]], # a, b
#     [[1, 1], [2, 1]], # a, c
#     [[1, 1], [2, 1]], # b, c
# ]

# mats = [
#     [[1,2],[2,1]],   # a,b
#     [[3,4],[-1,-2]], # a,c
#     [[-1,-2],[3,4]], # a,d
#     [[4,5],[6,7]],   # b,c
#     [[0,1],[1,0]],   # b,d
#     [[1,5],[-2,-3]], # c,d
# ]

#mats = [
#    [[1,2,3],[4,5,6],[2,1,3]],       # a,b
#    [[3,4,5],[6,7,8],[-1,-2,-3]],    # a,c
#    [[-1,-2,-3],[-4,-5,-6],[3,4,5]], # a,d
#    [[4,5,6],[7,8,9],[6,7,8]],       # b,c
#    [[0,1,2],[3,4,5],[1,0,-1]],      # b,d
#    [[1,5,9],[2,3,4],[-2,-3,-4]],    # c,d
#]

# Q = qform.intersection_form(g, mats)
#print(f'Signature: {qform.signature(Q)}')
#print(f'Definiteness: {qform.definiteness(Q)}')
#print(f'Parity: {qform.parity(Q)}\n')

n = 4
g = 2

while True:
    cuts = qform.random_diagram(n, g)
    random_Q = qform.intersection_form(g, qform.intersection_data(g, cuts))
    # determinant = random_Q.determinant()
    definiteness = qform.definiteness(random_Q)
    parity = qform.parity(random_Q)
    signature = qform.signature(random_Q)
    sign = signature[0] - signature[1]

    if parity == 'odd':
        print(f'Cuts = {cuts}\nQ =\n{random_Q}')
        print(f'Q = {signature[0]}<1> + {signature[1]}<-1>')
        print()
    if parity =='even':
        print(f'Definiteness = {definiteness}')
        print(f'Parity = {parity}')
        print(f'Signature = {signature}, {sign}')
        print(f'Q = {ZZ(sign / 8) if sign % 8 == 0 else 0}*E_8 + {n - sign}*[[0,1],[1,0]]\n')
        if sign == 8:
            print('*********')
            print(f'Cuts = {cuts}\nQ =\n{random_Q}')
            break
        print(f'Cuts = {cuts}\nQ =\n{random_Q}')
        print()
        #break

    if definiteness != 'indefinite':
        break

    det1 = random_Q[0:2,0:2].determinant()
    det2 = random_Q[2:4,2:4].determinant()
    if (det1 == 1 or det1 == -1) and (det2 == 1 or det2 == -1) and definiteness != 'indefinite':
        break
