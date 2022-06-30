E8 = matrix([[2, -1, 0, 0, 0, 0, 0, 0],
    [-1,2,-1,0,0,0,0,0],
    [0,-1,2,-1,0,0,0,0],
    [0,0,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,-1],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,0],
    [0,0,0,0,-1,0,0,2]])

B = matrix([[

Q = QuadraticForm(ZZ, B + B.transpose())
Q.compute_definiteness()
Q8 = QuadraticForm(ZZ, E8 + E8.transpose())
Q8.compute_definiteness()

signature = Q.signature_vector()
sign = signature[0] - signature[1]

# print(f'E8 Definiteness: {Q8.compute_definiteness_string_by_determinants()}')
# print(f'E8 Signature: {Q8.signature_vector()}')
#print(f'E8 Parity: {Q8.parity()}')
print(f'Definiteness: {Q.compute_definiteness_string_by_determinants()}')
print(f'Signature: {signature}, {sign}')
print(f'Parity: {Q.parity()}')
