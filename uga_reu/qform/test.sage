B = matrix([[  0, -1,  0,  0,  9,  0, -2,  0],
[ -1,  0,  1,  0,-10,  2,  5, -1],
[  0,  1, -2,  0, 20, -4,-10,  2],
[  0,  0,  0,  0,  9,  0, -2,  0],
[  9,-10, 20,  9, 32,-10,-21,  5],
[  0,  2, -4,  0,-10,  2,  5, -1],
[ -2,  5,-10, -2,-21,  5, 68,-10],
[  0, -1,  2,  0,  5, -1,-10,  2]])

E8 = matrix([[2, -1, 0, 0, 0, 0, 0, 0],
    [-1,2,-1,0,0,0,0,0],
    [0,-1,2,-1,0,0,0,0],
    [0,0,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,-1],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,0],
    [0,0,0,0,-1,0,0,2]])

Q = QuadraticForm(ZZ, B + B.transpose())
Q.compute_definiteness()
Q8 = QuadraticForm(ZZ, E8 + E8.transpose())
Q8.compute_definiteness()

print(f'E8 Definiteness: {Q8.compute_definiteness_string_by_determinants()}')
print(f'E8 Signature: {Q8.signature()}')
#print(f'E8 Parity: {Q8.parity()}')
print(f'Definiteness: {Q.compute_definiteness_string_by_determinants()}')
print(f'Signature: {Q.signature()}')
#print(f'Parity: {Q.parity()}')
