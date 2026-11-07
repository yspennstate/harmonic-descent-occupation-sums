"""Second-route check of P0341: estimate the expected band visits by simulating the harmonic descent chain
directly (jump d from m with P(d) = 1/(d h_{m-1})), not by the forward recursion in band_visits.py."""
import math, random
G = 0.5772156649015329
H = [0.0]
for j in range(1, 201):
    H.append(H[-1] + 1.0 / j)


def h(x):
    return H[x] if x <= 200 else math.log(x) + G + 1 / (2 * x) - 1 / (12 * x * x) + 1 / (120 * x ** 4)


def jump(m, rng):
    t = rng.random() * h(m - 1)          # smallest d with h(d) >= t
    lo, hi = 1, m - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if h(mid) >= t:
            hi = mid
        else:
            lo = mid + 1
    return lo


def visits(n, L, chains, rng):
    tot = tot2 = 0
    for _ in range(chains):
        m, v = n, 0
        while True:
            m -= jump(m, rng)
            if m <= n - L - 1 or m < 1:
                break
            v += 1
        tot += v; tot2 += v * v
    mu = tot / chains
    return mu, math.sqrt((tot2 / chains - mu * mu) / chains)


rng = random.Random(20261108)
cases = [(2 ** 10, 64), (2 ** 14, 1024), (2 ** 18, 16384), (2 ** 22, 16384)]
cases += [(2 ** (j * j), int(2 ** (j * j) / j ** 4)) for j in (3, 4, 5)]
for n, L in cases:
    mu, se = visits(n, L, 200000, rng)
    print(n, L, "mc %.4f +- %.4f" % (mu, se), flush=True)
