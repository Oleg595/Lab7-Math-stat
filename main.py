import numpy as np
import math
from scipy import integrate
import sympy
import matplotlib.pyplot as plt

def sum(x, n):
    result = 0.
    for i in range(n):
        result += x[i]
    return result

def sum_quad(x, n, m):
    result = 0.
    for i in range(n):
        result += (x[i] - m) * (x[i] - m)
    return result

def normal_params(x, n):
    result = np.zeros(2)
    s = sum(x, n)
    result[1] = s / n
    s_q = sum_quad(x, n, result[1])
    result[0] = math.sqrt(s_q / n)
    return result[0], result[1]

def get_cut(k):
    a = -1.01
    b = 1.56
    c = (b - a) / (k - 2)
    result = np.zeros(k + 1)
    result[0] = -1000
    for i in range(1, k):
        result[i] = a + (i - 1) * c
    result[k] = 1000
    return result

def rect_integral(f,xmin,xmax,n, m, sigma):
    dx=(xmax-xmin)/n
    area=0
    x=xmin
    for i in range(n):
        area+=dx*f(x, m, sigma)
        x+=dx
    return area

def func(x, m, sigma):
    result = 1 / (math.sqrt(2 * math.pi) * sigma) * math.exp(-((x - m) ** 2) / (2 * sigma * sigma))
    return result

def get_p(delta, k, sigma, m):
    result = np.zeros(k)
    result[0] = rect_integral(func, delta[0], delta[1], 1000000, m, sigma)
    sum = result[0]
    for i in range(1, k - 1):
        result[i] = rect_integral(func, delta[i], delta[i + 1], 100000, m, sigma)
        sum += result[i]
    result[k - 1] = 1 - sum
    return result

def get_n(delta, k, x, n):
    result = np.zeros(k)
    for i in range(n):
        j = 1
        while x[i] > delta[j] and j < k:
            j += 1
        result[j - 1] += 1
    return result

def div(n, p, k, size):
    result = np.zeros(k)
    for i in range(k):
        result[i] = n[i] - size * p[i]
    return result

def chi2(n, p, k, size):
    result = np.zeros(k)
    time = div(n, p, k, size)
    d = size * p
    for i in range(k):
        result[i] += time[i] * time[i] / d[i]
    return result

def sum(x, size):
    result = 0.
    for i in range(size):
        result += x[i]
    return result

print()
print("Normal, n = 100")
x = np.random.normal(0, 1, 100)
sig, m = normal_params(x, 100)
print(m, sig)
k = int(1 + 3.3 * math.log10(100))
delta = get_cut(k)
print("delta(i):", delta)
p = get_p(delta, k, sig, m)
print("p(i):", p, "sum =", sum(p, k))
n = get_n(delta, k, x, 100)
print("n(i):", n, "sum =", sum(n, k))
print("n * p(i):", 100 * p, "sum =", sum(100 * p, k))
d = div(n, p, k, 100)
print("n(i) - n * p(i):", d, "sum =", sum(d, k))
chi_2 = chi2(n, p, k, 100)
print("(n(i) - n * p(i)) ^ 2 / (n * p(i)):", chi_2, "sum =", sum(chi_2, k))

print()
print("Laplace, n = 20")
x = np.random.laplace(0, 1, 20)
sig, m = normal_params(x, 20)
print(m, sig)
k = int(1 + 3.3 * math.log10(20))
delta = get_cut(k)
print("delta(i):", delta)
p = get_p(delta, k, sig, m)
print("p(i):", p, "sum =", sum(p, k))
n = get_n(delta, k, x, 20)
print("n(i):", n, "sum =", sum(n, k))
print("n * p(i):", 20 * p, "sum =", sum(20 * p, k))
d = div(n, p, k, 20)
print("n(i) - n * p(i):", d, "sum =", sum(d, k))
chi_2 = chi2(n, p, k, 20)
print("(n(i) - n * p(i)) ^ 2 / (n * p(i)):", chi_2, "sum =", sum(chi_2, k))

print()
print("Uniform, n = 20")
x = np.random.uniform(0, 1, 20)
sig, m = normal_params(x, 20)
print(m, sig)
k = int(1 + 3.3 * math.log10(20))
delta = get_cut(k)
print("delta(i):", delta)
p = get_p(delta, k, sig, m)
print("p(i):", p, "sum =", sum(p, k))
n = get_n(delta, k, x, 20)
print("n(i):", n, "sum =", sum(n, k))
print("n * p(i):", 20 * p, "sum =", sum(20 * p, k))
d = div(n, p, k, 20)
print("n(i) - n * p(i):", d, "sum =", sum(d, k))
chi_2 = chi2(n, p, k, 20)
print("(n(i) - n * p(i)) ^ 2 / (n * p(i)):", chi_2, "sum =", sum(chi_2, k))
