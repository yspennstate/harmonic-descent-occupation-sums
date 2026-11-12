# Occupation sums of the harmonic descent chain and slowly decaying weights

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

The harmonic descent chain jumps from a state $m\ge2$ to $m-d$ with probability proportional to $1/d$. Started at $n$,
it records the clade sizes on the path from the root to a uniform leaf of the critical beta-splitting tree with $n$
leaves. Let $a(n,i)$ be the probability that the chain started at $n$ visits $i$, and let $a(i)$ be its limit as
$n\to\infty$. Aldous and Janson asked whether the occupation sums $\sum_{i\le n}a(n,i)f(i)$ converge to
$\sum_i a(i)f(i)$ for every non-negative weight $f$ of polynomial growth for which the limit sum is finite.
Brandenberger, Chin and Mossel proved the convergence when $f(i)=O(i^{-\varepsilon})$. We prove it whenever
$f(i)=O(1/\log i)$, and we show that this rate is the threshold: for every $\varphi$ tending to infinity there is a
weight $f(i)\le\varphi(i)/\log i$, tending to zero, with a finite limit sum and unbounded occupation sums. The
companion statement for infinite limit sums fails for weights of linear growth, where the ratio of the two sides tends
to $\pi^2/6$, and for a bounded weight.

## Code

Python 3 with NumPy. Run `python code/check_bounds.py` from the repository root; it prints `code/check_bounds.out`.

`check_bounds.py` computes the occupation probabilities in double precision from $a(n,n)=1$ and
$a(n,i)=\sum_{m=i+1}^{n}a(n,m)/((m-i)\,h_{m-1})$, and checks the statements of the paper on them:

- Part 1, $n=4000$: the identities of Proposition 1.3(a) for $g(i)=(i-1)/h_{i-1}$, the sums for $f(i)=i$, and the
  bounds of Lemma 2.3 and Corollary 2.4(a).
- Part 2: the bound of Lemma 2.2 on $S_t(i)=\sum_{\ell=i}^{i+t}a(\ell,i)$ and the monotonicity of $a(m,i)$ in $m$.
- Part 3: the expected number of visits $V(n,L)$ to the window $\{n-L,\dots,n-1\}$ for $n$ up to $2^{40}$, with the
  lower bound of Corollary 3.2 (`nan` where the corollary does not apply) and the upper bound of Lemma 2.1
  (Table 1 of the paper).

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
