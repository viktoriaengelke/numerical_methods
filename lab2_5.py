import matplotlib.pyplot as plt
import numpy as np
import math

def f_true(x):
    return math.sin(x)-0.001*x**2

def f_deriv(x):
    return math.cos(x) - 0.002 * x

def lin_spline_a_b(X, Y):
    n = len(X)
    a = [0]*n
    b = [0]*n
    for i in range(1, n):
        a[i] = Y[i-1]
        b[i] = (Y[i]-Y[i-1])/(X[i]-X[i-1])
    return a, b

def lin_spline(x, X, a, b):
    P = 0
    for i in range(1, len(X)):
        if X[i-1] <= x <= X[i]:
            P = a[i] + b[i]*(x - X[i-1])
            break
    return P

def cub_spline_coeff(X, Y):
    n = len(X)
    a = [0]*n
    b = [0]*n
    c = [0]*(n+2)
    d = [0]*n
    g = [0]*(n+1)
    s = [0]*(n+1)
    t = [0]*(n+1)
    h = [0]*n
    for i in range(1, n):
        h[i] = X[i]-X[i-1]
    for i in range(1, n - 1):
        g[i+1] = 3*((Y[i+1] - Y[i])/h[i+1] - (Y[i] - Y[i-1])/h[i])
        s[i+1] = -h[i+1]/(h[i]*s[i]+2*(h[i]+h[i+1]))
        t[i+1] = (g[i+1]-h[i]*t[i])/(h[i]*s[i]+2*(h[i]+h[i+1]))
    for i in range(n-1, 0, -1):
        c[i+1] = s[i+1]*c[i+2] + t[i+1]
    for i in range(1, n):
        a[i] = Y[i-1]
        b[i] = (Y[i] - Y[i-1])/h[i] - (h[i]*(c[i+1] + 2*c[i]))/3
        d[i] = (c[i+1] - c[i])/(3*h[i])
    return a, b, c, d

def cub_spline(x, X, a, b, c, d):
    for i in range(1, len(X)):
        if X[i-1] <= x <= X[i]:
            P = a[i] + b[i]*(x - X[i-1]) + c[i]*(x - X[i-1])**2 + d[i]*(x - X[i-1])**3
            break
    return P

def hermite_spline_coeff(X, Y, Yd):
    n = len(X)
    h = [0.0] * n
    a = [0.0] * n
    b = [0.0] * n
    c = [0.0] * n
    d = [0.0] * n
    for i in range(1, n):
        h[i] = X[i] - X[i-1]
        a[i] = Y[i-1]
        b[i] = Yd[i-1]
        c[i] = (3*Y[i] - 3*Y[i-1] - 2*h[i]*Yd[i-1] - h[i]*Yd[i]) / h[i]**2
        d[i] = (2*Y[i-1] - 2*Y[i] + h[i]*Yd[i-1] + h[i]*Yd[i]) / h[i]**3
    return a, b, c, d

def hermite_spline(x, X, a, b, c, d):
    P = 0.0
    for i in range(1, len(X)):
        if X[i-1] <= x <= X[i]:
            dx = x - X[i-1]
            P = a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
            break
    return P

X = [0.000, 1.500, 3.000, 4.500, 6.000]
Y = [0.000, 0.995, 0.132, -0.998, -0.315]
n = len(X)

a1, b1 = lin_spline_a_b(X, Y)
a2, b2, c2, d2 = cub_spline_coeff(X, Y)

Yd = [f_deriv(x) for x in X]
a3, b3, c3, d3 = hermite_spline_coeff(X, Y, Yd)

t = np.linspace(min(X), max(X), 300)
lin_splain_vals = [lin_spline(x, X, a1, b1) for x in t]
cub_splain_vals = [cub_spline(x, X, a2, b2, c2, d2) for x in t]
hermite_vals = [hermite_spline(x, X, a3, b3, c3, d3) for x in t]

plt.figure(figsize=(10, 6))
plt.plot(t, lin_splain_vals, 'b-', linewidth=1, label='Линейный сплайн')
plt.plot(t, cub_splain_vals, 'r-', linewidth=1, label='Кубический сплайн')
plt.plot(t, hermite_vals, 'm-', linewidth=1, label='Эрмитов сплайн')
plt.plot(t, [f_true(x) for x in t], 'g--', linewidth=1, label='f_true(x)')
plt.plot(X, Y, 'ro', markersize=5, label='Данные')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

