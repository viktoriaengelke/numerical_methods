import matplotlib.pyplot as plt
import numpy as np
import math

def right_diff(Y, Ay):
    diff = [Y]
    n = len(Y)
    for j in range(1, n):
        prev = diff[-1] 
        if len(prev) < 2:        
            break
        curr = [prev[i+1] - prev[i] for i in range(len(prev)-1)]
        diff.append(curr)
    max_j = 0
    for j in range(1, n+1):
        max_val = max(abs(v) for v in diff[j])
        w = 2**j*Ay
        if max_val >= w:
            max_j = j
        else: break
    return diff, max_j

def gauss_2(z, X, Y, diff, max_j, Ay):

    h = X[1] - X[0]
    k = 0
    min_dist = float('inf')
    for i, xi in enumerate(X):
        dist = abs(z - xi)
        if dist < min_dist:
            min_dist = dist
            k = i
    xk = X[k]
    t = (z - xk) / h

    result = Y[k]
    term_t = 1.0
    fact = 1.0

    coefs_abs = [1.0]

    for j in range(1, max_j + 1):
        if j % 2 == 1:
            if j == 1:
                term_t *= t
            else:
                term_t *= (t**2 - ((j // 2)**2))
        else:       
            term_t *= (t + (j // 2))   
        fact *= j
        diff_idx = k - (j + 1) // 2
        diff_val = diff[j][diff_idx]

        term = (term_t / fact) * diff_val
        result += term   
        coefs_abs.append(abs(term_t / fact))  
    j_next = max_j + 1

    if j_next < len(diff):
        term_t_next = term_t
        if j_next % 2 == 1:
            if j_next == 1:
                term_t_next *= t
            else:
                term_t_next *= (t**2 - ((j_next // 2)**2))
        else:
            term_t_next *= (t + (j_next // 2))
        
        fact_next = fact * j_next
        diff_idx_next = k - (j_next + 1) // 2
        diff_val_next = diff[j_next][diff_idx_next]
        
        delta_M = abs((term_t_next / fact_next) * diff_val_next)
    else:
        delta_M = 0.0

    delta_eps = Ay 
    for j in range(1, max_j + 1):
        delta_eps += coefs_abs[j] * (2**j) * Ay
        
    delta_total = delta_M + delta_eps
    delta_rel = delta_total / abs(result)
    
    return result, delta_M, delta_eps, delta_total, delta_rel
  
Y = [10.824, 10.431, 9.918, 9.310, 8.631, 7.911, 7.173, 6.442, 5.735, 5.068, 4.451]
Ay = 0.5*10**(-3)
X = np.linspace(1.0, 2.0, 11) 
z = 1.68
diff, max_j = right_diff(Y, Ay)
answer, dM, dEps, dTotal, dRel = gauss_2(z, X, Y, diff, max_j, Ay)

if dTotal > 0:
    k_err = math.floor(math.log10(dTotal))
    d_rounded = math.ceil(dTotal / (10**k_err)) * (10**k_err)
else:
    k_err = 0
    d_rounded = 0.0

m_res = math.floor(math.log10(abs(answer)))
n_znakov = m_res - k_err + 1

if k_err < 0:
    answer_final = round(answer, -k_err)
else:
    answer_final = round(answer, 0)

print(f"Погрешность метода: {dM:.4e}")
print(f"Неустранимая погрешность: {dEps:.4e}")
print(f"Полная абсолютная: {dTotal:.4e}")
print(f"Относительная: {dRel:.4e}")
print(answer_final)
print(max_j)
