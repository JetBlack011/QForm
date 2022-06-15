D = matrix(QQ, [[-1,-1], [-1,0]])
A = copy(D)
T = matrix.identity(QQ, 2)

for i in reversed(range(1, D.nrows())):
    D[i-1,:] += D[i,:]
    T[i-1,:] += T[i,:]

for i in range(D.nrows()):
    for j in range(i + 1, D.ncols()):
        a = D[j][i] / D[i][i]
        D[j,:] -= a * D[i,:]
        T[j,:] -= a * T[i,:]

        b = D[i][j] / D[i][i]
        D[:,j] -= b * D[:,i]
        T[:,j] -= b * T[:,i]

T = T.transpose()

print(D)
print(T)

print(T.transpose() * A * T)
