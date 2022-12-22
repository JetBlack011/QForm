import qform
import os
import sys
import time
import multiprocessing

# Genus
g = 2

# How many processes to run this computation on
process_count = 1

# The intersection blocks
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

def f(run):
    n = 2
    g = 2

    while True:
        cuts = qform.random_diagram(g, n, limit=100)
        print(cuts)
        return
        Q = qform.intersection_form(g, qform.cuts_to_intersection(g, cuts))
        # determinant = Q.determinant()
        definiteness = qform.definiteness(Q)
        parity = qform.parity(Q)
        signature = qform.signature(Q)
        sign = signature[0] - signature[1]

        if not run.is_set():
            break

        print(f'Definiteness = {definiteness}')
        print(f'Parity = {parity}')
        print(f'Signature = {signature}, {sign}')
        print(f'Cuts = {cuts}\nQ =\n{Q}')
        if parity == 'even' and sign != 0:
            print("*******")
            run.clear()
            return

def g(run):
    E8 = matrix([[2, -1, 0, 0, 0, 0, 0, 0],
        [-1,2,-1,0,0,0,0,0],
        [0,-1,2,-1,0,0,0,0],
        [0,0,-1,2,-1,0,0,0],
        [0,0,0,-1,2,-1,0,-1],
        [0,0,0,0,-1,2,-1,0],
        [0,0,0,0,0,-1,2,0],
        [0,0,0,0,-1,0,0,2]])

    Mp = matrix(2)

    # Q8 = QuadraticForm(ZZ, E8 + E8.transpose())
    # signature8 = Q8.signature_vector()
    # sign8 = signature8[0] - signature8[1]
    # 
    # print(f'E8 Definiteness: {Q8.compute_definiteness_string_by_determinants()}')
    # print(f'E8 Signature: {signature8}, {sign8}')

    while run.is_set():
        M = qform.random_unimodular_matrix(8)
        #print(M)
        Mp = M.transpose() * E8 * M
        #print(Mp)

        Q = QuadraticForm(ZZ, Mp + Mp.transpose())
        signature = Q.signature_vector()
        sign = signature[0] - signature[1]

        # print(f'Mp = {Mp}')
        # print(f'Definiteness: {Q.compute_definiteness_string_by_determinants()}')
        # print(f'Signature: {signature}, {sign}')
        det = Mp[:4,:4].determinant()
        Qsmall = QuadraticForm(ZZ, Mp[:4,:4] + Mp[:4,:4].transpose())
        small_sig = Qsmall.signature_vector()
        small_sign = small_sig[0] - small_sig[1]

        print(f'Mp = {Mp[:2,:2]}')
        print(f'Definiteness: {Qsmall.compute_definiteness_string_by_determinants()}')
        print(f'Signature: {small_sig}, {small_sign}')

        if det == 1 or det == -1:
            print("*********")
            run.clear()
            return

if __name__ == '__main__':
    processes = []
    manager = multiprocessing.Manager()
    return_code = manager.dict()
    run = manager.Event()
    run.set()
    for i in range(process_count):
        process = multiprocessing.Process(target=f, args=[run])
        processes.append(process)
        process.start()

    for process in processes:
        process.join()
