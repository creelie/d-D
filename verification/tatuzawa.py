"""Theorem 1.4: Tatuzawa's bound with eps = 1/log(2^40) at |d| = 2^40."""
from math import log, pi, e
D = 2.0 ** 40
eps = 1 / log(D)
assert eps < 1 / 11.2 and 1 / eps >= 11.2
h = 0.655 * eps * D ** (0.5 - eps) / pi
assert abs(h - 0.655 * 2 ** 20 / (e * pi * log(D))) < 1e-6
print(f"eps = {eps:.6f}, lower bound for h(d), |d| > 2^40: {h:.3f}")
