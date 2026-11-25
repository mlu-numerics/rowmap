# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy"]
# ///
"""Symbolic check of the DIFKON / DIFFU2 source terms in docs/problems.md."""

import sympy as sp

t, x, y, eps = sp.symbols("t x y epsilon", real=True)
pi = sp.pi
s = sp.sin(pi * x) * sp.sin(pi * y)
cx = sp.cos(pi * x) * sp.sin(pi * y)
cy = sp.sin(pi * x) * sp.cos(pi * y)


def lap(u):
    return sp.diff(u, x, 2) + sp.diff(u, y, 2)


# DIFKON
E = sp.exp(-2 * eps * pi**2 * t)
u = s * (1 + 4 * x * y * E)
g_true = sp.diff(u, t) - eps * lap(u) + y * sp.diff(u, x) - x * sp.diff(u, y)
g0 = 2 * eps * pi**2 * s + pi * (y * cx - x * cy)
g1 = -8 * eps * pi * (y * cx + x * cy) + 4 * pi * x * y * (y * cx - x * cy) + 4 * s * (y**2 - x**2)
assert sp.simplify(g_true - (g0 + E * g1)) == 0, "DIFKON g"
assert sp.simplify(sp.diff(g0 + E * g1, t) - (-2 * eps * pi**2 * E * g1)) == 0, "DIFKON g_t"

# DIFFU2
u = s * (1 + 4 * x * y * sp.sin(t))
q_true = sp.diff(u, t) - lap(u)
q = 4 * x * y * s * sp.cos(t) + 2 * pi**2 * s * (1 + 4 * x * y * sp.sin(t)) - 8 * pi * sp.sin(t) * (y * cx + x * cy)
qt = -4 * x * y * s * sp.sin(t) + 8 * pi**2 * x * y * s * sp.cos(t) - 8 * pi * sp.cos(t) * (y * cx + x * cy)
assert sp.simplify(q_true - q) == 0, "DIFFU2 q"
assert sp.simplify(sp.diff(q, t) - qt) == 0, "DIFFU2 q_t"

print("source terms OK")
