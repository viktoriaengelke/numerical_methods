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

def aitken(x, X, Y):
    n = len(X)
    L = [[0.0]*(n+1) for i in range(n)]
    for i in range(n):
        L[i][0] = X[i]
        L[i][1] = Y[i]
    for p in range(2, n + 1):
        for k in range(0, n - p + 1):
            L[k][p] = (L[k][p-1]*(L[k+p-1][0]-x) - (L[k][0]-x)*L[k+1][p-1])/(L[k+p-1][0] - L[k][0])
    return L

def aitken_value(x, X, Y):
    L = aitken(x, X, Y)
    return L[0][len(X)]

def print_aitken_table(L, n):
    for i in range(n):
        for j in range(n+1):
            print(f"{L[i][j]:<12.6f}", end="  ")
        print('\n')
    return 

def aitken_Al(x, X, a):
    n = len(X)
    A = 0.0
    for i in range(n):
        Y_test = [0.0]*n
        Y_test[i] = 1.0
        Li = aitken_value(x, X, Y_test)
        A += abs(Li)
    return A*a

x_point = 3.053
X = [0.000, 1.500, 3.000, 4.500, 6.000]
Y = [0.000, 0.995, 0.132, -0.998,-0.315]
n = len(X)
a = 0.5 * 10**(-3)

tracemalloc.start()
start = time.perf_counter()
L = aitken(x_point, X, Y)
end = time.perf_counter()
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print_aitken_table(L, n)
AL = aitken_Al(x_point, X, a)
AL_otn = AL/L[0][n]
L_x = round_(L[0][n], AL)

print(f"L (rounded) = {L_x},  AL = {AL:.3e},  AL_otn = {AL_otn:.3e}")
print(f"Время выполнения (Эйткен):  {(end-start)*1e6:.2f} мкс")
print(f"Занимаемая память (Эйткен):  {peak / 1024:.2f} КБ")