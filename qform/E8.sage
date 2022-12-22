import qform

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

while not (Mp[:2,:2].determinant() == 1 or Mp[:2,:2].determinant() == -1):
    M = qform.random_unimodular_matrix(8)
    print(M)
    Mp = M.transpose() * E8 * M
    print(Mp)

    Q = QuadraticForm(ZZ, Mp + Mp.transpose())
    signature = Q.signature_vector()
    sign = signature[0] - signature[1]

    print(f'Definiteness: {Q.compute_definiteness_string_by_determinants()}')
    print(f'Signature: {signature}, {sign}')

#print(qform.random_symplectic_matrix(2, 10)[0])

# g = 2
# n = 2
# 
# while True:
#     cuts = qform.random_diagram(g, n)
#     Q = qform.intersection_form(g, qform.intersection_data(g, cuts))
#     signature = qform.signature(Q)
#     sign = signature[0] - signature[1]
# 
#     if sign == 4:
#         print(f'Q =\n{Q}\nCuts = {cuts}')
#         print(f'Signature = {signature}, {sign}')
#         print(f'Parity = {qform.parity(Q)}')
#         break

# g = 2
# n = 3
# 
# while True:
#     B = matrix.identity(2 * g)
#     limit = 100
# 
#     while True:
#         limit_flag = True
#         cuts = [[(1, 0, 0, 0), (0,1,0,0)], [(0,0,1,0),(0,0,0,1)], [(1,1,0,0),(-2,1,1,-1)]]
# 
#         for i in range(len(cuts)):
#             cuts[i] = list(map(lambda x: vector((x)), cuts[i]))
# 
#         for _ in range(n - 1):
#             l = 0
#             while True and l < limit:
#                 l += 1
#                 symp_map = qform.random_symplectic_matrix(g, 50)[0]
#                 image = [symp_map * cuts[-1][i] for i in range(g)]
#                 det = matrix(cuts[-1] + image).determinant()
#                 if det == 1 or det == -1:
#                     break
#             limit_flag = limit_flag and l < limit
#             cuts.append(image)
# 
#         l = 0
#         while True and l < limit:
#             l += 1
#             symp_map = qform.random_symplectic_matrix(g, 50)[0]
#             image = [symp_map * cuts[-1][i] for i in range(g)]
#             det1 = matrix(cuts[-1] + image).determinant()
#             det2 = matrix(image + cuts[0]).determinant()
# 
#             if (det1 == 1 or det1 == -1) and (det2 == 1 or det2 == -1):
#                 break
#         cuts.append(image)
# 
#         limit_flag = limit_flag and l < limit
#         if limit_flag:
#             break
# 
#     Q = qform.intersection_form(g, qform.intersection_data(g, cuts))
#     signature = qform.signature(Q)
#     sign = signature[0] - signature[1]
#     parity = qform.parity(Q)
#     print(f'Q =\n{Q},\nCuts = {cuts}\nSignature = {signature}, {sign}\nParity = {qform.parity(Q)}')
#     if parity == 'even':
#         break
