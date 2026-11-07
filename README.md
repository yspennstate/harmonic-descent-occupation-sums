# Occupation sums of the harmonic descent chain and slowly decaying weights

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

The harmonic descent chain jumps from $m$ to $m-d$, $1\le d\le m-1$, with probability $1/(d\,h_{m-1})$. Let $a(n,i)$
be the probability that the chain started at $n$ visits $i$, and $a(i)=\lim_n a(n,i)$. Aldous and Janson
(arXiv:2303.02529, Ansatz 20) asked whether $\sum_{i=2}^n a(n,i)f(i)\to\sum_{i\ge2}a(i)f(i)$ for every non-negative $f$
of polynomial growth with a finite limit sum; Brandenberger, Chin and Mossel (arXiv:2505.24821) proved this for
$f(i)=O(i^{-\varepsilon})$. The note shows that for every $\delta>0$ there is an $f$ with $f(i)=O((\log i)^{-1+\delta})$
and $\sum_i a(i)f(i)<\infty$ whose occupation sums are unbounded along $n=2^{j^2}$, with $f(n)=0$ at those $n$. The
reason is that the chain started at $n$ makes at least of order $\log n/\log\log n$ expected visits to a window of width
$n/(\log n)^2$ just below $n$. Whether the statement holds for $f(i)=O((\log i)^{-1})$ is left open.

## Code

Python 3 with numpy.

- `band_visits.py`: the expected number of visits $\sum_{i=n-L}^{n-1}a(n,i)$ of the chain started at $n$ to the window
  just below $n$, computed exactly from $a(n,n)=1$ and $a(n,i)=\sum_{m=i+1}^{n}a(n,m)/((m-i)h_{m-1})$, in $O(L^2)$
  operations. Output `band_visits.out` (Table 1 of the paper).
- `cor_check.py`: the constants in the proof of Corollary 2.2 evaluated exactly for $3\le j\le 12$ and $j=20,40,80$:
  $T\le L_j$, $q\ge1-u$ and $q^T\ge e^{-9}$. Output `cor_check.out`.
- `audit_mc.py`: the same window visits estimated by simulating the chain (200,000 runs per window), with no use of
  the recursion. Output `audit_mc.out`; every estimate lies within 2.5 standard errors of `band_visits.out`.

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
