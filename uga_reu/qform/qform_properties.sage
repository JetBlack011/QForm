import qform

Qs = [
    [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]],
    [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]],
    [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,-1]],
    [[0,1,0,0],[1,0,0,0],[0,0,-1,0],[0,0,0,-1]],
    [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
    [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]],
    [[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]],
    [[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]],
    [[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]],
    [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]
]

for Q in Qs:
    Q = matrix(Q)
    signature = qform.signature(Q)
    sign = signature[0] - signature[1]
    parity = qform.parity(Q)
    print(f'Q =\n{Q}\nSignature = {signature},{sign}\nParity = {parity}')
