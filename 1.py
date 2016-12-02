from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math

l = 15.0
n = 100
# m = 1000
tau = 5.0/n
h = l/n
a = 1.0

U = lambda x, t: np.sin(x)*t
f = lambda x, t: np.sin(x)*t + np.sin(x)
phi = lambda x: U(x, 0)
psi0 = lambda t: U(0, t)
psil = lambda t: U(l, t)

X1, T1 = np.mgrid[0:15:100j, 0:5:100j]
Z_sol = U(X1, T1)
F_res = f(X1, T1)
# print(Z_sol[10,:])

def tridiagonal_alg(A, C, B, F):
    n = len(F)
    alf = [0 for i in range(n)]
    bet = [0 for i in range(n)]
    X = [0 for i in range(n)]
    alf[1], bet[1] = -B[0]/C[0], F[0]/C[0]
    for i in range(2, n):
        alf[i] = -B[i-1]/(A[i-1]*alf[i-1] + C[i-1])
        bet[i] = (F[i-1] - A[i-1]*bet[i-1])/(A[i-1]*alf[i-1] + C[i-1])
    X[n-1] = (F[n-1] - A[n-1]*bet[n-1])/(A[n-1]*alf[n-1] + C[n-1])
    for i in reversed(range(n-1)):
        # print i
        X[i] = alf[i+1]*X[i+1] + bet[i+1]
    # print alf
    # print bet
    return X

ar = [1.0 for i in range(4)]
ar2 = [3.0 for i in range(4)]
f2 = [float(i) for i in range(5, 9)]
print(tridiagonal_alg(ar, ar2, ar, f2))
print(f2)

X = [i*h for i in range(n)]
T = [i*tau for i in range(n)]
U_res = [[0 for i in range(n)] for j in range(n)]
U_res2 = [[0 for i in range(n)] for j in range(n)]
# F_res = [[f(X[j], T[i]) for j in range(n)] for i in range(n)]

for i in range(n):
    U_res[i][0] = phi(X[i])
    U_res[0][i] = psi0(T[i])
    U_res[n-1][i] = psil(T[i])

for i in range(n):
    for j in range(n):
        U_res2[i][j] = U(X[i], T[j])

A = [1 for i in range(n-2)]
B = [-(2 + h**2/tau/a**2) for i in range(n-2)]
C = [1 for i in range(n-2)]

for i in range(1, n):
    F = [-(F_res[j][i] + U_res[j][i-1]/tau)*(h**2)/(a**2) for j in range(1, n-1)]
    F[0] -= U_res[0][i]
    F[n-3] -= U_res[n-1][i]
    temp = tridiagonal_alg(A, B, C, F)
    for j in range(1, n-1):
        U_res[j][i] = temp[j-1]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X1, T1, U_res, rstride=8, cstride=8, alpha=0.3, color='green')
ax.plot_surface(X1, T1, Z_sol, rstride=8, cstride=8, alpha=0.3, color='red')

ax.set_xlabel('X')
ax.set_xlim(0*h, 15)
ax.set_ylabel('T')
ax.set_ylim(0, 5)
ax.set_zlabel('U')
ax.set_zlim(-7, 7)

plt.show()
