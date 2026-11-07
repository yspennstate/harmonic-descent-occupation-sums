"""Check of the constants in Corollary 2.2: exact q = h_x / h_{N-1}, T = floor(j^2/log j), x = floor(L/T); report T q^T vs e^-9 T and T<=L."""
import math
g = 0.5772156649015329
def h(m): return math.log(m) + g + 1/(2*m) - 1/(12*m*m) if m > 50 else sum(1/k for k in range(1, m+1))
for j in list(range(3, 13)) + [20, 40, 80]:
    N = 2**(j*j); L = N // j**4; T = int(j*j/math.log(j)); x = L // T
    q = h(x)/h(N-1); u = (2+6*math.log(j))/(1+math.log(N))
    print(j, "T", T, "T<=L", T <= L, "q", round(q, 5), "1-u", round(1-u, 5), "q>=1-u", q >= 1-u, "q^T", f"{q**T:.3g}", ">=e^-9", q**T >= math.exp(-9), "Tu/(1-u)", round(T*u/(1-u), 3))
