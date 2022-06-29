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

while True:
    cuts = qform.random_diagram(n)
    # print(f'{cuts}\n')
    random_Q = qform.intersection_form(2, qform.intersection_data(cuts))
    #print(f'Q =\n{random_Q}')
    definiteness = qform.definiteness(random_Q)
    if definiteness != 'degenerate':
        parity = qform.parity(random_Q)
        try:
            signature = qform.signature(random_Q)
            sign = signature[0] - signature[1]
        except Exception:
            sign = 0
            print(f'Signature = ?')
        #if parity == 'odd':
            #print(f'Q = {signature[0]}<1> + {signature[1]}<-1>')
        if parity =='even':
            print(f'Definiteness = {definiteness}')
            print(f'Parity = {parity}')
            print(f'Signature = {signature}, {sign}')
            print(f'Q = {ZZ(sign / 8) if sign % 8 == 0 else 0}*E_8 + {n - sign}*[[0,1],[1,0]]')
            if sign == 8:
                print('*********')
                print(f'cuts = {cuts}\nQ =\n{random_Q}')
                break
