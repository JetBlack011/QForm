import qform
import os
import sys
import time
import multiprocessing

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

def f(run):
    n = 6
    g = 2

    while run.is_set():
        cuts = qform.random_diagram(g, n)
        random_Q = qform.intersection_form(g, qform.intersection_data(g, cuts))
        # determinant = random_Q.determinant()
        definiteness = qform.definiteness(random_Q)
        parity = qform.parity(random_Q)
        signature = qform.signature(random_Q)
        sign = signature[0] - signature[1]

        print(f'Definiteness = {definiteness}')
        print(f'Parity = {parity}')
        print(f'Signature = {signature}, {sign}')
        print(f'Cuts = {cuts}\nQ =\n{random_Q}')
        if parity == 'odd':
            print(f'Q = {signature[0]}<1> + {signature[1]}<-1>')
            print()
        if parity =='even':
            print(f'Q = {ZZ(sign / 8) if sign % 8 == 0 else 0}*E_8 + {n - sign}*[[0,1],[1,0]]\n')
            if sign != 0 and sign % 8 == 0:
                print('*************************')
                run.clear()
                return (random_Q, definiteness, signature, sign, cuts)
            print()

        #if sign == -4:
        #    print("***********")
        #    run.clear()
        #    return

        # if sign != 0 and sign % 8 == 0:
        #     run.clear()
        #     print("***************")
        #     return (random_Q, definiteness, signature, sign, cuts)

        # if definiteness != 'indefinite':
        #     run.clear()
        #     print("***************")
        #     return (random_Q, definiteness, signature, sign, cuts)

if __name__ == '__main__':
    processes = []
    manager = multiprocessing.Manager()
    return_code = manager.dict()
    run = manager.Event()
    run.set()
    for i in range(16):
        process = multiprocessing.Process(target=f, args=[run])
        processes.append(process)
        process.start()

    for process in processes:
        process.join()
