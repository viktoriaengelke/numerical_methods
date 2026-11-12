import matplotlib.pyplot as plt
import numpy as np
import math
import time
import tracemalloc

def round_(x, x_):
    if x == 0:
        return 0.0
    k = 0
    m = -math.floor(math.log10(abs(x)))
    while x_ <= 0.5*10**k: k -= 1
    k += 1
    n = m - k + 1
    return round(x/10**m, n-1)*10**m

def newton(x, X, Y):
    n = len(X)
    N = [[0.0]*(n+1) for i in range(n)] алмр
    for i in range(n):
        N[i][0] = X[i]
        N[i][1] = Y[i]
    for p in range(2, n + 1):
        for k in range(0, n - p + 1):
            N[k][p] = (N[k+1][p-1]-N[k][p-1])/(N[k+p-1][0] - N[k][0])
    return N

def newton_value(x, X, Y):
    N = newton(x, X, Y)
    n = len(X)
    S = 0
    P = 1
    for i in range(1, n+1):
        S += N[0][i]*P
        P *= x - N[i-1][0]
    return S

def print_newton_table(N, n):
    for i in range(n):
        for j in range(n+1):
            print(f"{N[i][j]:<12.6f}", end="  ")
        print('\n')
    return 

def newton_Al(x, X, a):
    n = len(X)
    A = 0.0
    for i in range(n):
        Y_test = [0.0]*n
        Y_test[i] = 1.0
        Ni = newton_value(x, X, Y_test)
        A += abs(Ni)
    return A*a 

def f_true(x):
    return math.sin(x)-0.001*x**2

x_point = 3.053
X = [0.000, 1.500, 3.000, 5.000, 6.000]
Y = [0.000, 0.995, 0.132, -0.983, -0.315]
n = len(X)
a = 0.5 * 10**(-3)

tracemalloc.start()
start = time.perf_counter()
N = newton(x_point, X, Y)
end = time.perf_counter()
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

N_value = newton_value(x_point, X, Y)

print_newton_table(N, n)
AN = newton_Al(x_point, X, a)
AN_otn = AN/N[0][n]
N_x = round_(N_value, AN)

print(f"N (rounded) = {N_x},  AN = {AN:.3e},  AN_otn = {AN_otn:.3e}")
print(f"Время выполнения (Ньютон):  {(end-start)*1e6:.2f} мкс")
print(f"Занимаемая память (Ньютон):  {peak / 1024:.2f} КБ")

t = np.linspace(0, 6, 500)
N_vals = []
for x in t:
    N_graph= newton_value(x, X, Y)     
    N_vals.append(N_graph)

F_vals = [f_true(x) for x in t]

plt.figure(figsize=(10, 6))
plt.plot(t, F_vals, 'b-',  linewidth=1, label='Истинная f(x)')
plt.plot(t, N_vals, 'g--', linewidth=1, label='Многочлен Ньютона N(x)')
plt.plot(X, Y, 'ro', markersize=6, label='Узлы интерполирования')
plt.plot(x_point, N_value, 'm*', markersize=12, label=f'x = {x_point}')

plt.grid(True, alpha=0.3)
plt.xlabel('x'); plt.ylabel('y')
plt.title('Интерполяция Ньютона')
plt.legend()
plt.show()
