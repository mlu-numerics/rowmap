"""Python (numpy) implementation of the test problems in docs/problems.md.

Must agree with fortran/experiments/problems.f90; checked by tools/crosscheck.py.
Grid arrays are indexed [i, j(, l)] with x index i fastest in the flat vector
(Fortran order).
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.sparse as sp

PI = np.pi
AUTONOMOUS = {"heat", "nilidi", "bruss", "food", "food3", "vdpl"}
PAPER_M = {"heat": 50, "nilidi": 69, "bruss": 50, "food": 50, "food3": 10,
           "vdpl": 20, "difkon": 69, "diffu2": 69}
# coarse grids for cargo test (AGENTS.md decision: small in tests, paper size in examples)
COARSE_M = {"heat": 10, "nilidi": 10, "bruss": 10, "food": 10, "food3": 5,
            "vdpl": 8, "difkon": 10, "diffu2": 10}


@dataclass
class Problem:
    name: str
    m: int
    n: int
    t0: float
    tend: float
    y0: np.ndarray
    f: Callable          # f(t, y) -> dy
    jv: Callable         # jv(t, y, v) -> J v
    ft: Callable         # ft(t, y) -> df/dt
    jac: Callable        # jac(t, y) -> scipy.sparse matrix (for Radau)


def _flat(a):
    return np.asarray(a).ravel(order="F")


def _grid(a, m, dim):
    return np.asarray(a).reshape((m,) * dim, order="F")


# ---------------------------------------------------------------- 1D operators
def _d2_dir(m, h):
    return sp.diags([np.ones(m - 1), -2 * np.ones(m), np.ones(m - 1)], [-1, 0, 1]) / h**2


def _d1_dir(m, h):
    return sp.diags([-np.ones(m - 1), np.ones(m - 1)], [-1, 1]) / (2 * h)


def _d2_neu(m, h):
    a = sp.diags([np.ones(m - 1), -2 * np.ones(m), np.ones(m - 1)], [-1, 0, 1]).tolil()
    a[0, 1] = 2.0      # ghost mirror u_0 = u_2
    a[m - 1, m - 2] = 2.0
    return a.tocsr() / h**2


def _kron_lap(d2, dim):
    m = d2.shape[0]
    eye = sp.identity(m, format="csr")
    if dim == 2:   # x fastest: index i + (j-1) m  ->  kron(I_y, D_x) + kron(D_y, I_x)
        return (sp.kron(eye, d2) + sp.kron(d2, eye)).tocsr()
    return (sp.kron(eye, sp.kron(eye, d2)) + sp.kron(eye, sp.kron(d2, eye))
            + sp.kron(d2, sp.kron(eye, eye))).tocsr()


def make(name: str, m: int | None = None) -> Problem:
    m = m or PAPER_M[name]
    zero = lambda t, y: np.zeros_like(y)

    if name in ("heat", "nilidi", "difkon", "diffu2"):
        lo, length = (0.0, PI / 3) if name == "nilidi" else (0.0, 1.0)
        h = length / (m + 1)
        xs = lo + h * np.arange(1, m + 1)
        X, Y = np.meshgrid(xs, xs, indexing="ij")
        L = _kron_lap(_d2_dir(m, h), 2)
        n = m * m
        eye = sp.identity(m, format="csr")
        Dx = sp.kron(eye, _d1_dir(m, h)).tocsr()
        Dy = sp.kron(_d1_dir(m, h), eye).tocsr()
        xf, yf = _flat(X), _flat(Y)

        if name == "heat":
            y0 = _flat(16 * X * (1 - X) * Y * (1 - Y))
            return Problem(name, m, n, 0.0, 10.24, y0, lambda t, y: L @ y,
                           lambda t, y, v: L @ v, zero, lambda t, y: L)
        if name == "nilidi":
            y0 = _flat(np.sin(3 * X) * np.sin(3 * Y))

            def f(t, y):
                eu = np.exp(y)
                return eu * (L @ y) + y * (18 * eu - 1)

            def jv(t, y, v):
                eu = np.exp(y)
                return eu * (L @ v + (L @ y) * v) + (18 * eu * (1 + y) - 1) * v

            def jac(t, y):
                eu = np.exp(y)
                return (sp.diags(eu) @ L + sp.diags(eu * (L @ y) + 18 * eu * (1 + y) - 1)).tocsr()

            return Problem(name, m, n, 0.0, 1.0, y0, f, jv, zero, jac)

        if name == "difkon":
            eps = 0.1
            A = (eps * L - sp.diags(yf) @ Dx + sp.diags(xf) @ Dy).tocsr()
            s = np.sin(PI * xf) * np.sin(PI * yf)
            cx = np.cos(PI * xf) * np.sin(PI * yf)
            cy = np.sin(PI * xf) * np.cos(PI * yf)
            g0 = 2 * eps * PI**2 * s + PI * (yf * cx - xf * cy)
            g1 = (-8 * eps * PI * (yf * cx + xf * cy) + 4 * PI * xf * yf * (yf * cx - xf * cy)
                  + 4 * s * (yf**2 - xf**2))
            E = lambda t: np.exp(-2 * eps * PI**2 * t)
            y0 = s * (1 + 4 * xf * yf)
            return Problem(name, m, n, 0.0, 1.0, y0, lambda t, y: A @ y + g0 + E(t) * g1,
                           lambda t, y, v: A @ v, lambda t, y: -2 * eps * PI**2 * E(t) * g1,
                           lambda t, y: A)

        # diffu2
        s = np.sin(PI * xf) * np.sin(PI * yf)
        cx = np.cos(PI * xf) * np.sin(PI * yf)
        cy = np.sin(PI * xf) * np.cos(PI * yf)
        w = yf * cx + xf * cy
        q = lambda t: (4 * xf * yf * s * np.cos(t) + 2 * PI**2 * s * (1 + 4 * xf * yf * np.sin(t))
                       - 8 * PI * np.sin(t) * w)
        qt = lambda t: (-4 * xf * yf * s * np.sin(t) + 8 * PI**2 * xf * yf * s * np.cos(t)
                        - 8 * PI * np.cos(t) * w)
        return Problem(name, m, n, 0.0, 1.0, s.copy(), lambda t, y: L @ y + q(t),
                       lambda t, y, v: L @ v, lambda t, y: qt(t), lambda t, y: L)

    # ------------------------------------------------ Neumann, two components
    dim = 3 if name == "food3" else 2
    h = 1.0 / (m - 1)
    xs = h * np.arange(m)
    L = _kron_lap(_d2_neu(m, h), dim)
    npts = m**dim
    n = 2 * npts
    G = np.meshgrid(*([xs] * dim), indexing="ij")
    G = [_flat(g) for g in G]

    def split(y):
        return y[:npts], y[npts:]

    def block(a11, a12, a21, a22):
        return sp.bmat([[a11, a12], [a21, a22]], format="csr")

    if name == "bruss":
        X, Y = G
        y0 = np.concatenate([0.5 + Y, 1 + 5 * X])
        al = 0.02

        def f(t, y):
            u, v = split(y)
            return np.concatenate([al * (L @ u) + 1 + u * u * v - 4 * u,
                                   al * (L @ v) + 3 * u - u * u * v])

        def jv(t, y, w):
            u, v = split(y)
            a, b = split(w)
            return np.concatenate([al * (L @ a) + (2 * u * v - 4) * a + u * u * b,
                                   al * (L @ b) + (3 - 2 * u * v) * a - u * u * b])

        def jac(t, y):
            u, v = split(y)
            return block(al * L + sp.diags(2 * u * v - 4), sp.diags(u * u),
                         sp.diags(3 - 2 * u * v), al * L - sp.diags(u * u))

        return Problem(name, m, n, 0.0, 1.0, y0, f, jv, zero, jac)

    if name == "food":
        X, Y = G
        y0 = np.concatenate([10 - 5 * np.cos(PI * X) * np.cos(10 * PI * Y),
                             17 + 5 * np.cos(10 * PI * X) * np.cos(PI * Y)])

        def f(t, y):
            u, v = split(y)
            return np.concatenate([0.05 * (L @ u) + u * (1 - 0.1 * v),
                                   L @ v + v * (-1000 + 100 * u)])

        def jv(t, y, w):
            u, v = split(y)
            a, b = split(w)
            return np.concatenate([0.05 * (L @ a) + (1 - 0.1 * v) * a - 0.1 * u * b,
                                   L @ b + 100 * v * a + (-1000 + 100 * u) * b])

        def jac(t, y):
            u, v = split(y)
            return block(0.05 * L + sp.diags(1 - 0.1 * v), sp.diags(-0.1 * u),
                         sp.diags(100 * v), L + sp.diags(-1000 + 100 * u))

        return Problem(name, m, n, 0.0, 3.0, y0, f, jv, zero, jac)

    if name == "food3":
        X, Y, Z = G
        b = (1 + 0.2 * X * Y * Z) * (1e6 - 1 + 1e-6)
        y0 = np.concatenate([500 + 250 * np.cos(PI * X) * np.cos(3 * PI * Y) * np.cos(10 * PI * Z),
                             200 + 150 * np.cos(10 * PI * X) * np.cos(PI * Y) * np.cos(3 * PI * Z)])

        def f(t, y):
            u, v = split(y)
            return np.concatenate([0.05 * (L @ u) + u * (b - 1e6 * u - v),
                                   L @ v + v * (b - (1e6 - 1) * u - 1e6 * v)])

        def jv(t, y, w):
            u, v = split(y)
            a, c = split(w)
            return np.concatenate([0.05 * (L @ a) + (b - 2e6 * u - v) * a - u * c,
                                   L @ c - (1e6 - 1) * v * a + (b - (1e6 - 1) * u - 2e6 * v) * c])

        def jac(t, y):
            u, v = split(y)
            return block(0.05 * L + sp.diags(b - 2e6 * u - v), sp.diags(-u),
                         sp.diags(-(1e6 - 1) * v), L + sp.diags(b - (1e6 - 1) * u - 2e6 * v))

        return Problem(name, m, n, 0.0, 10.0, y0, f, jv, zero, jac)

    if name == "vdpl":
        X, Y = G
        y0 = np.concatenate([1 - np.cos(PI * X) * np.cos(2 * PI * Y),
                             1 + 2 * np.cos(2 * PI * X) * np.cos(PI * Y)])

        def f(t, y):
            u, v = split(y)
            return np.concatenate([0.05 * (L @ u) + v, L @ v + 100 * (1 - u * u) * v - u])

        def jv(t, y, w):
            u, v = split(y)
            a, b = split(w)
            return np.concatenate([0.05 * (L @ a) + b,
                                   L @ b + (-200 * u * v - 1) * a + 100 * (1 - u * u) * b])

        def jac(t, y):
            u, v = split(y)
            return block(0.05 * L, sp.identity(npts), sp.diags(-200 * u * v - 1),
                         L + sp.diags(100 * (1 - u * u)))

        return Problem(name, m, n, 0.0, 1000.0, y0, f, jv, zero, jac)

    raise ValueError(f"unknown problem {name}")
