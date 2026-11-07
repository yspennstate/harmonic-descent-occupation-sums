"""Expected visits of the harmonic descent chain started at n to the band (n-L, n), root n excluded.
a(n,n)=1; a(n,k) = sum_{m=k+1}^{n} a(n,m) / ((m-k) h_{m-1}); only band states feed band states."""
import numpy as np, math, sys
def h(m): return math.log(m) + 0.5772156649015329 + 1/(2*m) - 1/(12*m*m) if m > 50 else sum(1/j for j in range(1, m+1))
def band(n, L):
    a = np.zeros(L+1); a[0] = 1.0           # index d = n - k, d=0 is the root
    inv = 1.0/np.arange(1, L+1)
    for d in range(0, L):
        w = a[d] / h(n - d - 1)
        a[d+1:] += w * inv[:L-d]
    return a[1:].sum()
for n, L in [(2**10, 64), (2**14, 1024), (2**18, 16384), (2**22, 16384), (2**26, 16384), (2**30, 16384), (2**40, 16384)]:
    print(n, L, round(band(n, L), 4), "ln(n)/ln(n/L) bound-ish:", round(math.log(n)/max(1e-9, math.log(n/L)), 3), flush=True)
# the 0/1 example: n = 2^{j^2}, band n/j^4
for j in (3, 4, 5):
    n = 2**(j*j); L = int(n / j**4)
    print("j", j, "n", n, "L", L, "visits", round(band(n, L), 4), "j^2/ln j", round(j*j/math.log(j), 3), flush=True)
