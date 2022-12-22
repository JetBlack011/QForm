D = matrix(QQ, [[1,2,3],[2,1,0],[3,0,2]])
D = matrix(QQ, [[2,1,1,0,1,0,0,0,0],[1,2,0,0,0,0,0,0,0],[1,0,2,1,0,0,0,0,0],[0,0,2,1,0,0,0,0,0],[
A = copy(D)
T = matrix.identity(QQ, 3)

# for i in reversed(range(1, D.nrows())):
#     D[i-1,:] += D[i,:]
#     D[:,i-1] += D[:,i]
#     T[:,i-1] += T[:,i]

for i in range(D.nrows()):
    for j in range(i + 1, D.ncols()):
        a = D[j][i] / D[i][i]
        D[j,:] -= a * D[i,:]
        #T[j,:] -= a * T[i,:]
        D[:,j] -= a * D[:,i]
        T[:,j] -= a * T[:,i]

print(f'D =\n{D}\n')
print(f'T =\n{T}\n')

print(f'T^t * A * T\n{T.transpose() * A * T}')
