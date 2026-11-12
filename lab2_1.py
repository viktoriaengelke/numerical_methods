import matplotlib.pyplot as plt
import numpy as np
import math
import time
import tracemalloc

def lagrange(x, n, X, Y):
    L = 0
    A = 0
    for i in range(n):
        F = 1
        for j in range(n):
            if i != j:
                F *= (x - X[j]) / (X[i] - X[j])
        L += F * Y[i]
        A += abs(F)
    return L, A*a  

def f_true(x):
    return math.sin(x)-0.001*x**2

def round_(x, x_):
    k = 0
    m = -math.floor(math.log10(abs(x)))
    while x_ <= 0.5*10**k: k -= 1
    k += 1
    n = m - k + 1
    return round(x/10**m, n-1)*10**m


x_point = 3.053
X = [0.000, 1.500, 3.000, 4.500, 6.000]
Y = [0.000, 0.995, 0.132, -0.998, -0.315]
n = len(X)
a = 0.5 * 10**(-3)

tracemalloc.start()
start = time.perf_counter()
L1, A1 = lagrange(x_point, n, X, Y)
end = time.perf_counter()
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
A2 = A1/L1

t = np.linspace(0, 6, 500)
L_vals = []
for x in t:
    L, _ = lagrange(x, n, X, Y)     
    L_vals.append(L)

F_vals = [f_true(x) for x in t]

L_rounded = round_(L1, A1)   

print(f"L (rounded) = {L_rounded},  AL = {A1:.3e},  AL_otn = {A2:.3e}")
print(f"L = {L1}")

print(f"Время выполнения (Лагранж):  {(end-start)*1e6:.2f} мкс")
print(f"Занимаемая память (Лагранж):  {peak / 1024:.2f} КБ")

plt.figure(figsize=(10, 6))
plt.plot(t, F_vals, 'b-',  linewidth=1, label='Истинная f(x)')
plt.plot(t, L_vals, 'g--', linewidth=1, label='Многочлен Лагранжа L(x)')
plt.plot(X, Y, 'ro', markersize=6, label='Узлы интерполирования')
plt.plot(x_point, L1, 'm*', markersize=12, label=f'x = {x_point}')

plt.grid(True, alpha=0.3)
plt.xlabel('x'); plt.ylabel('y')
plt.title('Интерполяция Лагранжа')
plt.legend()
plt.show()