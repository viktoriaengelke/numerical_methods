import matplotlib.pyplot as plt
import numpy as np
import math

def f_true(x):
    return math.sin(x)-0.001*x**2

def polinom(X, n, m):
    t_prev = [1 for xi in X]
    polys = [t_prev] 
    a = [None]*(m+1)
    b = [None]*(m+1)
    c = [None]*(m+1)
    A = [[0.0] * (m + 1) for _ in range(m + 1)]
    num = sum(y * tp for y, tp in zip(Y, t_prev))
    A[0][0] = sum(tp**2 for tp in t_prev)
    c[0] = num / A[0][0]
    a[1] = sum(X)/n
    t_curr = [xi - a[1] for xi in X]
    polys.append(t_curr)
    num = sum(y * tp for y, tp in zip(Y, t_curr))
    A[1][1] = sum(tp**2 for tp in t_curr)
    c[1] = num / A[1][1]
    for k in range(1, m):
        a[k+1] = sum(X[i]*t_curr[i]**2 for i in range(n))/sum(t_curr[i]**2 for i in range(n))
        b[k+1] = sum(t_curr[i]**2 for i in range(n))/sum(t_prev[i]**2 for i in range(n))
        t_next = [(X[i] - a[k+1])*t_curr[i] - b[k+1]*t_prev[i] for i in range(n)]
        polys.append(t_next)
        num = sum(y*t for y, t in zip(Y, t_next))
        A[k+1][k+1] = sum(t**2 for t in t_next)
        c[k+1] = num/A[k+1][k+1]
        t_prev, t_curr = t_curr, t_next
    return a, b, c, polys, A

def fast_polinom(x, X, n, m, a, b, c):
    t_prev = 1
    t_curr = x - sum(X)/n
    result = c[0]*t_prev + c[1]*t_curr
    for k in range(1, m):
        t_next = (x - a[k+1])*t_curr - b[k+1]*t_prev
        result += c[k+1]*t_next
        t_prev, t_curr = t_curr, t_next
    return result


# def gram_matrix(polys, m):
#     A = [[0.0] * (m + 1) for _ in range(m + 1)]
#     for i in range(m + 1):
#         for k in range(i, m + 1):             
#             A[i][k] = sum(pi * pk for pi, pk in zip(polys[i], polys[k]))
#             A[k][i] = A[i][k]                     
#     return A

def print_gram(A, m):
    print("\n Матрица Грама\n")
    for i in range(m):
        for j in range(m):
            print(f"{A[i][j]:<12.3f}", end="  ")
        print('\n')
    return 

x_point = 3.053
X = [0.000, 1.500, 3.000, 4.500, 6.000]
Y = [0.000, 0.995, 0.132, -0.998, -0.315]
n = len(X)
m = int(input("Введите степень многочлена: "))

a, b, c, polys, A = polinom(X, n, m)
y_point = fast_polinom(x_point, X, n, m, a, b, c)
print(f"\nf({x_point}) = {y_point:.6f}")

# A = gram_matrix(polys, m)
print_gram(A, m)

t = np.linspace(min(X) - 0.5, max(X) + 0.5, 300)
phi_vals = [fast_polinom(x, X, n, m, a, b, c) for x in t]

plt.figure(figsize=(10, 6))
plt.plot(t, phi_vals, 'b-', linewidth=1, label=f'f(x), m={m}')
plt.plot(t, [f_true(x) for x in t], 'g--', linewidth=1, label='f_true(x)')
plt.plot(X, Y, 'ro', markersize=8, label='Данные')
plt.plot(x_point, y_point, 'm*', markersize=10, label=f'x={x_point}')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()