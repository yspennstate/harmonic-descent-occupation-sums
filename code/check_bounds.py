"""Numerical checks for "Occupation sums of the harmonic descent chain and slowly decaying weights".

The chain jumps from m >= 2 to m - d with probability 1/(d h_{m-1}), 1 <= d <= m-1.
a(n, i) is the probability that the chain started at n visits i, and a(i) = 6 h_{i-1} / (pi^2 (i-1)).

Part 1  all a(n, i), 2 <= i <= n, for one n: the identity of Proposition 1.3(a), the linear weight,
        a(n, i) >= a(i), Lemma 2.3 and Corollary 2.4(a).
Part 2  a(m, i) for fixed i and m = i, ..., i + K: monotonicity in m and Lemma 2.2.
Part 3  Table 1: expected visits to the window {n-L, ..., n-1}, with the bounds of Lemma 2.1 and Corollary 3.2.

Usage: python check_bounds.py [n]      (default n = 4000; needs numpy)
"""
import math
import sys

import numpy as np

EULER = 0.5772156649015329


def h(m):
    """Harmonic number h_m."""
    if m <= 0:
        return 0.0
    if m > 50:
        return math.log(m) + EULER + 1.0 / (2 * m) - 1.0 / (12.0 * m * m)
    return sum(1.0 / j for j in range(1, m + 1))


def all_visits(n):
    """a[i] = a(n, i) for 1 <= i <= n, and the harmonic numbers hh[0..n]."""
    hh = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, n + 1))))
    inv = np.concatenate(([0.0], 1.0 / np.arange(1, n + 1)))
    a = np.zeros(n + 1)
    a[n] = 1.0
    for i in range(n - 1, 0, -1):
        # first-step decomposition over the last state m > i visited before i
        a[i] = np.dot(a[i + 1:n + 1] / hh[i:n], inv[1:n - i + 1])
    return a, hh


def column(i, K):
    """a(m, i) for m = i, ..., i + K."""
    hh = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, i + K + 1))))
    col = np.zeros(K + 1)
    col[0] = 1.0
    for t in range(1, K + 1):
        m = i + t
        # a(m, i) = sum_{d=1}^{t} a(m-d, i) / (d h_{m-1})
        col[t] = np.dot(col[t - 1::-1], 1.0 / np.arange(1, t + 1)) / hh[m - 1]
    return col, hh


def band(n, L):
    """sum_{i=n-L}^{n-1} a(n, i), computed inside the window."""
    a = np.zeros(L + 1)
    a[0] = 1.0
    inv = 1.0 / np.arange(1, L + 1)
    for d in range(0, L):
        a[d + 1:] += (a[d] / h(n - d - 1)) * inv[:L - d]
    return float(a[1:].sum())


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    a, hh = all_visits(n)
    idx = np.arange(2, n + 1)
    lim = 6.0 * hh[idx - 1] / (math.pi ** 2 * (idx - 1))
    an = a[2:]

    print("Part 1, n = %d" % n)
    s = float(np.sum(an * (idx - 1) / hh[idx - 1]))
    sl = float(np.sum(lim * (idx - 1) / hh[idx - 1]))
    print("  sum a(n,i)(i-1)/h_{i-1} = %.10f   n-1 = %d" % (s, n - 1))
    print("  sum a(i)(i-1)/h_{i-1}   = %.10f   ratio = %.12f   pi^2/6 = %.12f" % (sl, s / sl, math.pi ** 2 / 6))
    s1 = float(np.sum(an * idx))
    s1l = float(np.sum(lim * idx))
    print("  sum a(n,i) i = %.2f   n log n = %.2f   sum a(i) i = %.2f   (6/pi^2) n log n = %.2f"
          % (s1, n * math.log(n), s1l, 6 / math.pi ** 2 * n * math.log(n)))
    print("  min a(n,i)/a(i) over 2 <= i <= n: %.6f (BCM monotonicity gives >= 1)" % float(np.min(an / lim)))
    far = idx + idx // 2 <= n
    print("  max a(n,i)/a(i) over n >= i + floor(i/2): %.4f   (Corollary 2.4(a): <= pi^2 = %.4f)"
          % (float(np.max(an[far] / lim[far])), math.pi ** 2))
    t = n - idx
    near = t < idx
    bound = hh[idx[near]] / ((t[near] + 1) * (hh[idx[near]] - hh[t[near]]))
    print("  max a(n,i) / [h_i/((t+1)(h_i-h_t))], t = n-i < i: %.6f   (Lemma 2.3: <= 1)"
          % float(np.max(an[near] / bound)))

    print("Part 2, a(m,i) for m = i..i+K")
    for i in (2, 3, 7, 50, 300, 1000):
        K = 3000
        col, hc = column(i, K)
        incr = int(np.sum(np.diff(col) > 1e-15))
        S = np.cumsum(col)
        tt = np.arange(0, min(i, K + 1))
        worst = float(np.max(S[tt] * (hc[i] - hc[tt]) / hc[i]))
        print("  i = %4d: increases in m: %d   max_t S_t (h_i - h_t)/h_i = %.6f (Lemma 2.2: <= 1)   a(i+K,i)/a(i) = %.6f"
              % (i, incr, worst, col[-1] / (6.0 * hc[i - 1] / (math.pi ** 2 * (i - 1)))))

    print("Part 3, expected visits to {n-L,...,n-1}")
    print("  %-8s %-10s %-10s %-10s %-10s" % ("n", "L", "visits", "lower", "upper"))
    rows = [(2 ** 10, 2 ** 6), (2 ** 14, 2 ** 10), (2 ** 18, 2 ** 14), (2 ** 22, 2 ** 14), (2 ** 30, 2 ** 14),
            (2 ** 40, 2 ** 14), (2 ** 9, 2 ** 9 // 9), (2 ** 16, 2 ** 16 // 16)]
    for nn, L in rows:
        v = band(nn, L)
        up = h(L) / (h(nn - L - 1) - h(L)) if nn >= 2 * L + 2 else float("nan")
        lam = math.log(8.0 * nn / L)
        lo = h(nn - 1) / (10 * lam) if (8 <= L <= nn // 2 and h(nn - 1) >= 2 * lam) else float("nan")
        print("  2^%-6d %-10d %-10.4f %-10.4f %-10.4f" % (int(math.log2(nn)), L, v, lo, up))


if __name__ == "__main__":
    main()
