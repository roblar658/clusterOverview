#import "../template/ieee.typ": *

#alg_entry(
  id: 149,
  name: "STATPC: Statistical Projected Clustering",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2008",
  paradigm: "Statistical Non-Redundant Projected Clustering",
  complexity: $O(n_s d^2 n + |cal(H)|^2 n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

High-dimensional subspace and projected clustering algorithms frequently confront two coupled pathologies: arbitrary global density thresholds that degrade exponentially with subspace dimensionality @beyer1999nearest, and combinatorial redundancy whereby a single $r$-dimensional cluster spawns up to $2^r - 1$ overlapping lower-dimensional projections @kriegel2009subspacesurvey. STATPC (Statistical Projected Clustering) @moise2008statpc resolves both pathologies by casting projected clustering as a statistical optimization problem: discover a compact collection of *statistically significant, mutually non-redundant* axis-parallel hyper-rectangles whose point counts exhibit extreme upper-tail deviations from a uniform null hypothesis.

Let $cal(X) = \{bold(x)_1, dots.c, bold(x)_n\} subset.eq [0, 1]^d$ denote a finite dataset normalized coordinate-wise onto the unit hypercube $[0, 1]^d$. A *projected cluster region* is a pair $H = (S_H, B_H)$, where $S_H subset.eq \{1, dots.c, d\}$ is the relevant attribute subset of dimensionality $r_H = |S_H| >= 1$ and
$ B_H = product_(j in S_H) [l_(H, j), u_(H, j)] subset.eq [0, 1]^(r_H) $
is an axis-parallel hyper-rectangle spanning $[0, 1]$ on every irrelevant coordinate $j in.not S_H$. The covered point set (support) and normalized Lebesgue volume of $H$ are
$ C(H) = \{bold(x)_i in cal(X) : x_(i j) in [l_(H, j), u_(H, j)], forall j in S_H\}, $
$ m(H) = |C(H)|, quad p(H) = "vol"(B_H) = product_(j in S_H) (u_(H, j) - l_(H, j)). $

Under the null hypothesis $H_0$ that observations are drawn i.i.d. uniformly on $[0, 1]^d$, the random count $M(H) = |C(H)|$ follows the binomial law $M(H) tilde.op "Bin"(n, p(H))$ with expected support $mu(H) = n p(H)$. The exact upper-tail $p$-value and statistical surprisal (compactness score) of $H$ are
$ "pval"(H) = sum_(k = m(H))^n binom(n, k) p(H)^k (1 - p(H))^(n - k), $
$ F(H) = -log "pval"(H). $

#definition(title: "Statistical Significance and Non-Redundancy")[
  Fix a primary significance level $alpha_0 in (0, 1)$ and a residual significance level $alpha_H in [alpha_0, 1)$. A hyper-rectangle $H = (S_H, B_H)$ is *$alpha_0$-significant* if $m(H) > n p(H)$ and $"pval"(H) <= alpha_0$, i.e., $F(H) >= -log alpha_0$. Let $cal(H)_(alpha_0)$ denote the family of all $alpha_0$-significant regions. Given a set of already selected regions $cal(R)$, a candidate $H in cal(H)_(alpha_0)$ is *statistically redundant* with respect to $cal(R)$ if its unexplained residual support
  $ m_"res" (H | cal(R)) = |C(H) without union.big_(R in cal(R)) C(R)| $
  fails statistical significance at level $alpha_H$, i.e.,
  $ sum_(k = m_"res" (H | cal(R)))^n binom(n, k) p(H)^k (1 - p(H))^(n - k) > alpha_H. $
  Equivalently, in pairwise form, $(H_1, H_2) in E_"red"$ if either region renders the other non-significant at level $alpha_H$.
]

STATPC formalizes non-redundant projected clustering as maximizing aggregate statistical surprisal over a pairwise non-redundant subset of significant hyper-rectangles @moise2008statpc:
$ cal(R)^* = arg max_(cal(R) subset.eq cal(H)_(alpha_0)) sum_(H in cal(R)) F(H) \ "s.t." quad (H_a, H_b) in.not E_"red" quad forall H_a != H_b in cal(R). $
Because the continuous space of axis-parallel hyper-rectangles is infinite and discrete candidate selection is NP-hard, STATPC couples local 1D statistical neighborhood seeding with greedy subspace hill-climbing and descending-surprisal redundancy elimination.

#algbox(
  title: "STATPC (Statistical Projected Clustering)",
  inputs: [Normalized data $cal(X) subset.eq [0, 1]^d$; thresholds $alpha_0, alpha_H, alpha_K in (0, 1)$; seed budget $n_s <= n$.],
  outputs: [Non-redundant projected clusters $cal(R) = \{H_1, dots.c, H_K\}$ and point labels $ell_1, dots.c, ell_n in \{0, 1, dots.c, K\}$.],
)[
  1. *1D Statistical Interval Indexing*: For each coordinate $j in \{1, dots.c, d\}$, sort $\{x_(i j)\}_(i=1)^n$ and precompute minimal-width $alpha_K$-significant 1D intervals around candidate seeds.
  2. *Greedy Subspace Region Generation*: Sample $n_s$ seed points $cal(S) subset.eq cal(X)$. For each seed $bold(x)_s in cal(S)$, initialize $H$ from its highest-surprisal 1D interval and iteratively add (or refine) the coordinate $j in.not S_H$ maximizing $F(H union \{j\})$ until no single-axis addition or width adjustment strictly increases $F(H)$. Insert $H$ into $cal(H)_(alpha_0)$ if $F(H) >= -log alpha_0$.
  3. *Greedy Non-Redundant Selection*: Initialize $cal(R) <- emptyset$ and candidate pool $cal(Q) <- cal(H)_(alpha_0)$. While $cal(Q) != emptyset$:
     (a) Select $H^* <- arg max_(H in cal(Q)) F(H)$ and update $cal(R) <- cal(R) union \{H^*\}$.
     (b) Remove $H^*$ and every $H in cal(Q)$ rendered statistically redundant by $cal(R)$ at threshold $alpha_H$.
  4. *Optional Local Refinement & Labeling*: Re-optimize interval bounds of each $H_k in cal(R)$ on its assigned points. Assign each $bold(x)_i$ to $ell_i = arg max_(k : bold(x)_i in C(H_k)) F(H_k)$, or to outlier class $ell_i = 0$ if uncovered.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees characterizing STATPC: tight two-sided exponential bounds on statistical surprisal, an intrinsic noise-dimension self-regularization property, NP-hardness of optimal non-redundant selection, a neighborhood-independence approximation ratio for greedy pruning, and an explicit rotational failure counterexample.

#subsec("A", "Two-Sided Sanov-Chernoff Concentration of Surprisal")

Let $D_"KL" (q || p) = q log(q / p) + (1 - q) log((1 - q) / (1 - p))$ denote the binary Kullback-Leibler divergence for $p, q in (0, 1)$.

#lemma(title: "Exact Two-Sided Bounds on STATPC Surprisal", origin: "Original Proof")[
  Let $H = (S_H, B_H)$ have normalized volume $p = p(H) in (0, 1)$ and empirical support fraction $q = m(H) / n$ satisfying $p < q < 1$ (so $1 <= m(H) <= n - 1$). Then the statistical surprisal $F(H) = -log "pval"(H)$ satisfies
  $ n D_"KL" (q || p) <= F(H) <= n D_"KL" (q || p) + log(3 sqrt(n q (1 - q))). $
  When $m(H) = n$ ($q = 1$), $F(H) = n D_"KL" (1 || p) = n log(1 / p)$ holds with equality. Moreover, $F(H)$ is strictly increasing in $m(H)$ and strictly decreasing in $p(H)$.
]

#proof[
  For the upper bound on $"pval"(H)$, apply Markov's inequality to $e^(t M(H))$ for any $t > 0$:
  $ "pval"(H) = bb(P)(M(H) >= n q) <= e^(-t n q) (1 - p + p e^t)^n. $
  Because $q > p$, the unique minimizer satisfies $e^(t^*) = q (1 - p) / (p (1 - q)) > 1$. Substituting $t^*$ yields $1 - p + p e^(t^*) = (1 - p) / (1 - q)$ and hence
  $ "pval"(H) <= (q / p)^(-n q) ((1 - q) / (1 - p))^(-n (1 - q)) = exp(-n D_"KL" (q || p)). $
  Taking negative logarithms establishes the lower bound $F(H) >= n D_"KL" (q || p)$.

  For the upper bound on $F(H)$, retain the leading term $k = m = n q$ of the binomial sum:
  $ "pval"(H) >= binom(n, m) p^m (1 - p)^(n - m). $
  By Robbins' refinement of Stirling's formula, $sqrt(2 pi k) (k / e)^k e^(1 / (12 k + 1)) <= k! <= sqrt(2 pi k) (k / e)^k e^(1 / (12 k))$ for every integer $k >= 1$. Applying the lower bound to $n!$ and the upper bound to $m!$ and $(n - m)!$ yields
  $ binom(n, m) >= frac(e^(Delta_n), sqrt(2 pi n q (1 - q))) q^(-n q) (1 - q)^(-n (1 - q)), $
  where $Delta_n = 1 / (12 n + 1) - 1 / (12 m) - 1 / (12 (n - m))$. Since $m >= 1$ and $n - m >= 1$, we have $Delta_n > -1 / 6$, and numerical evaluation gives $e^(-1 / 6) / sqrt(2 pi) approx 0.3377 > 1 / 3$. Multiplying by $p^(n q) (1 - p)^(n (1 - q))$ gives
  $ "pval"(H) >= frac(1, 3 sqrt(n q (1 - q))) exp(-n D_"KL" (q || p)), $
  whose negative logarithm yields the upper bound. At $m = n$, $"pval"(H) = p^n = exp(-n D_"KL" (1 || p))$ directly. Finally, incrementing $m$ removes a strictly positive binomial probability mass term from $"pval"(H)$, while differentiating the incomplete beta representation gives $partial "pval"(H) / partial p = n binom(n - 1, m - 1) p^(m - 1) (1 - p)^(n - m) > 0$, proving strict monotonicity of $F(H)$ in both arguments.
]

#subsec("B", "Intrinsic Self-Regularization Against Noise Dimensions")

A crucial theoretical advantage of STATPC over grid-density methods such as CLIQUE @agrawal1998clique is that adding an irrelevant uniform coordinate automatically penalizes the statistical objective rather than inflating apparent density.

#proposition(title: "Strict Noise-Dimension Penalty Under Binomial Surprisal", origin: "Original Derivation")[
  Let $H = (S_H, B_H)$ have volume $p in (0, 1)$ and empirical fraction $q in (p, 1)$. Consider extending $H$ along a new coordinate $j in.not S_H$ with interval width $w_j in (0, 1)$, producing $H'$ with volume $p' = w_j p$ and point-retention ratio $gamma = m(H') / m(H) in (0, 1]$, so $q' = gamma q$.
  1. *Noise-axis degradation*: If coordinate $j$ is statistically independent and uniform on $[0, 1]$, the expected retention satisfies $bb(E)[gamma] = w_j$, and the asymptotic surprisal rate strictly decreases for every $w_j in (0, 1)$:
     $ D_"KL" (w_j q || w_j p) < D_"KL" (q || p). $
  2. *Relevant-axis gain*: Conversely, for any fixed retention $gamma in (0, 1]$ with $gamma q > w_j p$, there exists a unique critical width threshold $w_j^* (gamma) in (0, gamma)$ such that $D_"KL" (gamma q || w_j p) > D_"KL" (q || p)$ if and only if $w_j < w_j^* (gamma)$. In particular, if $gamma = 1$, every proper sub-interval $w_j < 1$ strictly increases $D_"KL"$.
]

#proof[
  Define $phi(t) = D_"KL" (t q || t p)$ on $t in (0, 1]$. Expanding and differentiating with respect to $t$ yields
  $ phi'(t) = q log(q / p) - q log(frac(1 - t q, 1 - t p)) - q + p frac(1 - t q, 1 - t p). $
  Introduce $u(t) = p (1 - t q) / (q (1 - t p))$. Since $0 < p < q < 1$ and $t in (0, 1]$, we have $p - t p q < q - t p q$, whence $0 < u(t) < 1$. Rewriting $phi'(t)$ in terms of $u(t)$ gives the exact identity
  $ phi'(t) = q (u(t) - 1 - log u(t)). $
  Because $f(u) = u - 1 - log u$ is strictly convex with unique global minimum $f(1) = 0$, we have $phi'(t) > 0$ strictly for all $t in (0, 1]$. Integrating from $t = w_j < 1$ to $t = 1$ proves $D_"KL" (w_j q || w_j p) = phi(w_j) < phi(1) = D_"KL" (q || p)$, establishing Part 1.

  For Part 2, fix $gamma in (0, 1]$ and consider $g(w) = D_"KL" (gamma q || w p)$ on $w in (0, gamma]$. Its derivative satisfies
  $ g'(w) = frac(w p - gamma q, w (1 - w p)) < 0 $
  because $w p <= gamma p < gamma q$. Moreover, $lim_(w -> 0^+) g(w) = +infinity$, whereas at the right endpoint $g(gamma) = phi(gamma) <= phi(1) = D_"KL" (q || p)$ (with strict inequality whenever $gamma < 1$). By strict monotonicity and the Intermediate Value Theorem, there exists a unique root $w_j^* (gamma) in (0, gamma)$ solving $g(w_j^*) = D_"KL" (q || p)$ when $gamma < 1$ (and $w_j^* (1) = 1$), such that $g(w_j) > D_"KL" (q || p)$ iff $w_j < w_j^* (gamma)$.
]

#subsec("C", "NP-Hardness of Non-Redundant Region Selection")

#theorem(title: "NP-Hardness of Optimal Non-Redundant Selection", origin: "Classical — Cited")[
  Given a dataset $cal(X) subset.eq [0, 1]^d$, fixed significance thresholds $alpha_0 = alpha_H in (0, 1)$, and a polynomial-size pool $cal(H)_(alpha_0)$ of $alpha_0$-significant hyper-rectangles, computing a pairwise non-redundant subset $cal(R)^* subset.eq cal(H)_(alpha_0)$ maximizing $sum_(H in cal(R)) F(H)$ is NP-hard @moise2008statpc.
]

#proof[
  We reduce from *Maximum Independent Set* (MIS) on an arbitrary simple undirected graph $G = (V, E)$ with $|V| = m >= 2$, maximum degree $Delta_G >= 1$, and edge set $E$. Set ambient dimension $d = m$, indexing coordinates by vertices $v in V$. For a width parameter $delta in (0, 1 / 2)$ chosen below, construct the dataset $cal(X)$ as follows:
  - For each edge $e = \{u, v\} in E$, insert one point $bold(x)_e$ with $x_(e, u) = x_(e, v) = delta / 2$ and $x_(e, w) = 1 - delta / 2$ for all $w in.not \{u, v\}$.
  - For each vertex $v in V$, insert $b_v = Delta_G - "deg"_G (v) >= 0$ private points placed at $x_v = delta / 2$ and $x_w = 1 - delta / 2$ for all $w != v$.
  The total dataset size $n = |E| + sum_(v in V) b_v <= m^2$ is polynomial in $m$. For each $v in V$, define the 1D projected candidate $H_v = (\{v\}, [0, delta])$ of volume $p(H_v) = delta$. By construction, $C(H_v)$ contains the $b_v$ private points of $v$ plus every edge point $bold(x)_e$ incident to $v$, yielding uniform support $m(H_v) = Delta_G$ for every $v in V$. Furthermore, for distinct $u, v in V$,
  $ |C(H_u) inter C(H_v)| = cases(1 &"if" \{u, v\} in E, 0 &"if" \{u, v\} in.not E.) $

  Let $T_k (p) = sum_(r=k)^n binom(n, r) p^r (1 - p)^(n - r)$. Because $T_(Delta_G) (p)$ is continuous and strictly increasing from $0$ to $1$ on $(0, 1)$, there is a unique $delta^* in (0, 1)$ with $T_(Delta_G) (delta^*) = alpha_0$. Since
  $ T_(Delta_G - 1) (delta^*) = alpha_0 + binom(n, Delta_G - 1) (delta^*)^(Delta_G - 1) (1 - delta^*)^(n - Delta_G + 1) > alpha_0, $
  continuity yields $delta in (0, min(delta^*, 1 / 2))$ satisfying $T_(Delta_G) (delta) < alpha_0 < T_(Delta_G - 1) (delta)$ (and $n delta < Delta_G$). Thus every $H_v$ is $alpha_0$-significant with identical score $F_0 = -log T_(Delta_G) (delta) > 0$. If $\{u, v\} in E$, then $|C(H_v) without C(H_u)| = Delta_G - 1$ yields residual $p$-value $T_(Delta_G - 1) (delta) > alpha_0$, so $(H_u, H_v) in E_"red"$; if $\{u, v\} in.not E$, the residual count remains $Delta_G$ with $p$-value $T_(Delta_G) (delta) < alpha_0$, so $(H_u, H_v) in.not E_"red"$. Hence $cal(R) subset.eq \{H_v : v in V\}$ is pairwise non-redundant iff $\{v : H_v in cal(R)\}$ is an independent set of $G$, and $sum_(H in cal(R)) F(H) = F_0 |cal(R)|$. Solving STATPC's selection problem therefore solves MIS in polynomial time.
]

#subsec("D", "Approximation Ratio of Greedy Redundancy Elimination")

Because exact optimization is NP-hard, Step 3 of STATPC greedily extracts the highest-surprisal remaining region and deletes its redundant neighbors. Let $G_"red" = (cal(H)_(alpha_0), E_"red")$ be the undirected redundancy conflict graph, let $Delta_"red"$ denote its maximum degree, and let $beta = max_(H in cal(H)_(alpha_0)) alpha_"ind" (N(H))$ denote its *neighborhood independence number* (the maximum number of mutually non-redundant candidates that can simultaneously conflict with a single region $H$).

#theorem(title: "Greedy Surprisal Approximation Guarantee", origin: "Original Proof")[
  Let $cal(R)_"greedy"$ be the non-redundant set produced by Step 3 of STATPC and let $cal(R)^*$ be any globally optimal non-redundant subset of $cal(H)_(alpha_0)$. Then
  $ sum_(H in cal(R)_"greedy") F(H) >= frac(1, max(1, beta)) sum_(H in cal(R)^*) F(H) >= frac(1, max(1, Delta_"red")) sum_(H in cal(R)^*) F(H). $
]

#proof[
  Enumerate $cal(R)_"greedy" = \{G_1, dots.c, G_k\}$ in the order selected by Step 3, so $F(G_1) >= F(G_2) >= dots.c >= F(G_k) > 0$. Because $cal(R)_"greedy"$ is a maximal independent set of $G_"red"$, every optimal region $H in cal(R)^*$ is either contained in $cal(R)_"greedy"$ or eliminated upon the selection of some earlier neighbor in $cal(R)_"greedy"$. Define a charging map $pi : cal(R)^* -> \{1, dots.c, k\}$ sending each $H in cal(R)^*$ to the smallest index $j$ such that $H in N[G_j] = \{G_j\} union N(G_j)$, and let $A_j = pi^(-1) (j) subset.eq cal(R)^*$.

  At iteration $j$, every $H in A_j$ remained in the active candidate pool $cal(Q)$ (otherwise $H$ would have conflicted with some $G_l$ for $l < j$, contradicting minimality of $j$). Greedy choice of $G_j$ therefore guarantees $F(H) <= F(G_j)$ for every $H in A_j$. Furthermore, $A_j$ is an independent set of $G_"red"$ contained in $N[G_j]$: if $G_j in A_j$, no neighbor of $G_j$ can appear in the independent set $cal(R)^*$, so $|A_j| = 1$; if $G_j in.not A_j$, then $A_j subset.eq N(G_j)$ implies $|A_j| <= alpha_"ind" (N(G_j)) <= beta <= Delta_"red"$. In both cases $|A_j| <= max(1, beta)$. Summing over all $j$ gives
  $ sum_(H in cal(R)^*) F(H) = sum_(j=1)^k sum_(H in A_j) F(H) <= sum_(j=1)^k |A_j| F(G_j) \ <= max(1, beta) sum_(j=1)^k F(G_j). $
  Dividing by $max(1, beta)$ and noting $beta <= Delta_"red"$ completes the proof.
]

#corollary(title: "Exact Recovery on Separated Mixtures", origin: "Original Proof")[
  If $G_"red"$ partitions into $K$ connected components $cal(H)_1, dots.c, cal(H)_K$ where each ground-truth region $H_k^* in cal(H)_k$ uniquely maximizes $F$ and conflicts with all $H in cal(H)_k without \{H_k^*\}$, then $cal(R)_"greedy" = \{H_1^*, dots.c, H_K^*\}$ achieves the exact global optimum.
]

#proof[
  Within each component $cal(H)_k$, Step 3 selects the unique maximizer $H_k^*$ first, pruning all remaining candidates in $cal(H)_k$ without affecting other components. Repeating across $k = 1, dots.c, K$ yields $cal(R)_"greedy" = \{H_1^*, dots.c, H_K^*\}$, which is globally optimal because any independent set selects at most one region per star component.
]

#subsec("E", "Coordinate Scaling Invariance and Rotational Failure")

#proposition(title: "Coordinate Affine Invariance and Diagonal Counterexample", origin: "Original Counterexample")[
  1. *Invariance*: After coordinate-wise min-max normalization onto $[0, 1]^d$, STATPC's volume $p(H)$, surprisal $F(H)$, and redundancy graph are strictly invariant under arbitrary positive coordinate-wise affine rescalings $x_(i j) |-> a_j x_(i j) + b_j$ ($a_j > 0$), point permutations, and feature permutations.
  2. *Rotational failure counterexample*: Consider $n >= 2$ points $bold(x)_i = (t_i, t_i)$ with $t_i = (i - 1) / (n - 1)$ on the diagonal $cal(L) = \{(t, t) : t in [0, 1]\} subset.eq [0, 1]^2$ (a 1D subspace rotated by $pi / 4$). Any axis-parallel box $B_H$ covering $m >= 2$ points satisfies $p(H) >= ((m - 1) / (n - 1))^2$. Consequently, the entire cluster ($m = n$) attains exact zero surprisal $F(H) = 0$ for every $n >= 2$, and for any fixed partial fraction $q = m / n in (0, 1)$ the asymptotic surprisal rate satisfies $lim sup_(n -> infinity) F(H) / n <= D_"KL" (q || q^2) <= e^(-1) < 0.368$, whereas rotating coordinates by $-pi / 4$ elevates the surprisal rate arbitrarily to $+infinity$.
]

#proof[
  Part 1 is immediate because for $a_j > 0$, $(a_j x_(i j) + b_j - (a_j min_l x_(l j) + b_j)) / (a_j (max_l x_(l j) - min_l x_(l j))) = (x_(i j) - min_l x_(l j)) / (max_l x_(l j) - min_l x_(l j))$ leaves all normalized coordinates identical; hyper-rectangle volumes and point counts are symmetric under coordinate and observation permutations.

  For Part 2, let $B_H = [l_1, u_1] times [l_2, u_2]$ cover a subset $C(H) subset.eq cal(L)$ of cardinality $m >= 2$. Every covered point $(t_i, t_i)$ obeys $t_i in [max(l_1, l_2), min(u_1, u_2)]$. Because consecutive points are spaced by $1 / (n - 1)$, covering $m$ points requires $min(u_1, u_2) - max(l_1, l_2) >= (m - 1) / (n - 1)$, forcing both side lengths $u_1 - l_1 >= (m - 1) / (n - 1)$ and $u_2 - l_2 >= (m - 1) / (n - 1)$, whence $p(H) = (u_1 - l_1)(u_2 - l_2) >= ((m - 1) / (n - 1))^2$. For $m = n$, this gives $p(H) = 1$, $"pval"(H) = 1$, and $F(H) = 0$ identically. For any sequence with $m / n -> q in (0, 1)$, we have $p_n = ((m - 1) / (n - 1))^2 -> q^2$, so Lemma 1 and monotonicity of $D_"KL" (q || p)$ in $p < q$ yield
  $ lim sup_(n -> infinity) frac(F(H), n) <= D_"KL" (q || q^2) = -q log q - (1 - q) log(1 + q) \ <= -q log q <= e^(-1), $
  since $(1 - q) log(1 + q) > 0$ and $-q log q$ attains its unique maximum $1 / e$ at $q = e^(-1)$. In contrast, rotating coordinates by $-pi / 4$ maps $cal(L)$ onto $u_2 = 0$, enclosed by a 1D interval of width $epsilon -> 0$ with $m = n$, giving $F(H) = n log(1 / epsilon) -> +infinity$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time Complexity.* Initial sorting along all $d$ coordinates costs $O(d n log n)$. Generating candidates from $n_s <= n$ seeds via greedy forward-backward subspace search evaluates at most $O(d)$ candidate axes per step over at most $r_"max" <= d$ steps; intersecting 1D membership sets costs $O(n)$ per evaluation (or $O(n / w)$ with $w$-bit hardware bitsets), yielding $O(n_s r_"max" d n) subset.eq O(n_s d^2 n)$. Sorting the $|cal(H)|$ significant candidates takes $O(|cal(H)| log |cal(H)|)$, and greedy redundancy elimination performs up to $O(|cal(H)|^2)$ residual support tests in $O(|cal(H)|^2 n)$ time (or $O(|cal(H)|^2 n / w)$ with bitset masks). Thus worst-case complexity is $O(n_s d^2 n + |cal(H)|^2 n)$, dropping sharply in practice when $n_s << n$ and clusters occupy low-dimensional subspaces $r_"max" << d$.

*Space and Cache Profile.* Storing the normalized data matrix requires $O(n d)$ memory; maintaining $|cal(H)|$ candidate bitsets and interval bounds requires $O(|cal(H)| (n / w + d))$ space. Bitset intersections execute sequential SIMD bitwise-AND operations with high L1/L2 cache hit rates, explaining STATPC's sub-millisecond runtime on moderate benchmarks.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n_s d^2 n + |cal(H)|^2 n)$; accelerated by $w$-bit bitsets ($1 / w$ factor)],
    [Space Complexity], [$O(n d + |cal(H)| (n / w + d))$ for data and candidate bitsets],
    [Cluster Geometry], [Axis-parallel hyper-rectangles in variable-dimension subspaces],
    [Cluster Count $K$], [Automatically inferred via significance ($alpha_0$) and non-redundancy ($alpha_H$)],
    [Noise / Outliers], [Explicitly supported (uncovered points assigned to background class $0$)],
    [Scale Invariance], [Strictly invariant to positive per-axis affine rescaling],
    [Rotation Invariance], [Not rotation-invariant (axis-aligned; fails on oblique subspaces)],
    [Determinism], [Deterministic if all $n$ seeds used; randomized if $n_s < n$ seeds sampled],
  ),
  caption: [Operational, geometric, and algorithmic profile of STATPC @moise2008statpc.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Dimensional Gene Expression and Functional Genomics.* In transcriptomic cohorts ($d tilde.op 10^3$–$10^4$ genes, $n tilde.op 10^2$–$10^3$ samples), disease subtypes manifest as tightly co-regulated expression intervals across small pathway subsets while irrelevant genes fluctuate randomly. By Proposition 1, STATPC's binomial surprisal rewards narrow co-expression boxes and penalizes noise genes, discovering interpretable gene-signature biclusters @cheng2000biclustering without presetting subtype counts.
2. *Multi-Channel Customer and Retail Behavioral Segmentation.* Enterprise transaction logs mix heterogeneous spending channels across disparate numerical scales. Because STATPC is invariant to coordinate-wise positive affine rescalings (Proposition 2) and outputs axis-parallel bounds $[l_(H, j), u_(H, j)]$, each discovered segment yields auditable business rules free of redundant sub-segment overlaps.
3. *Semiconductor Wafer-Test Yield and Root-Cause Binning.* Parametric electrical testing logs hundreds of voltage, leakage, and timing measurements per die, while fault mechanisms typically localize to 3–10 coupled parameters. STATPC's uniform null test isolates statistically over-represented failure boxes and suppresses nested sub-projections for rapid root-cause binning.
4. *Network Intrusion and Telemetry Anomaly Signatures.* Flow-level cybersecurity logs combine diffuse background traffic with coordinated scanning or botnet bursts tightly localized in a few protocol and packet-size attributes. Testing upper-tail count deviations relative to volume extracts low-false-positive signatures while relegating scattered benign traffic to the outlier class.
5. *Astronomical Photometric Catalogs and Multi-Band Surveys.* Multi-wavelength sky surveys record color indices across optical and infrared bands where stellar populations (white dwarfs, quasars, main-sequence branches) form sharp hyper-rectangles in distinct color subspaces, cleanly isolated by STATPC's dimension-adaptive surprisal without uniform density thresholds.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

STATPC's inductive bias breaks down under three concrete regimes: (i) *obliquely correlated or rotated manifolds*, where Proposition 2 shows axis-parallel hyper-rectangles shatter diagonal subspaces into fragmented stair-step boxes with capped surprisal rate $D_"KL" (q || q^2) <= e^(-1)$; (ii) *strongly non-uniform background marginals* (e.g., heavy-tailed power-law or skewed attributes), where the uniform null hypothesis $H_0 = "Unif"([0, 1]^d)$ flags benign high-density background regions as spurious clusters unless quantile/rank transformations are applied first; and (iii) *extreme sample sizes with dense candidate overlap*, where $O(|cal(H)|^2 n)$ redundancy checks grow costly without aggressive seed subsampling.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Gene Expression Subtypes], [Strong: sparse axis-aligned gene intervals; noise axes penalized], [Correlated factors: use ORCLUS @aggarwal2000orclus],
    [Retail & Tabular Segments], [Strong: scale-invariant, human-readable interval rules], [Heavy-tailed spend: apply rank normalization],
    [Semiconductor Fault Bins], [Strong: automatic $K$ and variable subspace dimensionality], [Curved drift: use DBSCAN @ester1996dbscan],
    [Correlated Sensor Arrays], [Weak: oblique linear correlations shatter into stair-steps], [Use SSC (Alg. 150) or ORCLUS @aggarwal2000orclus],
    [2D / 3D Curved Manifolds], [Weak: axis boxes cannot track crescents or rings], [Use SUBCLU @kailing2004subclu or spectral methods],
    [Massive Streaming Telemetry], [Moderate: $O(n_s d^2 n)$ offline hill-climbing overhead], [Subsample seeds $n_s$ or stream grids],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for STATPC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

STATPC occupies a central unifying position in Category VI, generalizing earlier axis-parallel subspace algorithms (CLIQUE @agrawal1998clique, ENCLUS @cheng1999enclus, PROCLUS @aggarwal1999proclus, DOC @procopiuc2002doc, FINDIT @woo2004findit, SUBCLU @kailing2004subclu, FIRES @kriegel2005fires, P3C @moise2006p3c, and INSCY @assent2008inscy) by replacing fixed density or support thresholds with dimension-adaptive binomial tail probabilities and formal redundancy elimination.

#proposition(title: "Asymptotic Dominance Over Fixed-Threshold Grid Subspace Methods", origin: "Original Derivation")[
  Consider CLIQUE @agrawal1998clique with fixed relative support threshold $tau in (0, 1)$ and fixed per-axis grid width $w in (0, 1)$, versus STATPC with significance level $alpha_0 in (0, 1)$. For any subspace dimensionality $r > log(tau) / log(w)$ (so cell volume $p_r = w^r < tau$) and any true cluster of volume $p_r$ and support fraction $q in (p_r, tau)$, CLIQUE provably misses the cluster for all $n$, whereas STATPC detects it with probability $1$ whenever
  $ n > frac(-log alpha_0, D_"KL" (q || p_r)). $
  Simultaneously, while CLIQUE's Apriori enumeration emits all $2^r - 1$ non-empty sub-projections of a dense $r$-dimensional cell, STATPC's residual redundancy filter collapses the entire chain to at most one non-redundant region.
]

#proof[
  By definition, CLIQUE retains a cell $H$ iff $m(H) / n >= tau$. Since $q < tau$, $H$ fails CLIQUE's criterion regardless of sample size $n$. Conversely, by Lemma 1, STATPC's surprisal obeys $F(H) >= n D_"KL" (q || p_r)$. Because $q > p_r$, $D_"KL" (q || p_r) > 0$, so choosing $n > (-log alpha_0) / D_"KL" (q || p_r)$ guarantees $F(H) > -log alpha_0$, i.e., $"pval"(H) < alpha_0$. For redundancy, every sub-projection $H' supset.eq H$ of an isolated cluster $C(H)$ satisfies $C(H') = C(H)$ and $p(H') = w^(r') > w^r = p(H)$ for $r' < r$, hence $F(H) > F(H')$ by Lemma 1 and $m_"res" (H' | \{H\}) = 0$, causing Step 3 to prune all $2^r - 2$ lower-dimensional projections immediately upon selecting $H$.
]

*Empirical Behavior on ClusterOverview Benchmarks.* On the 4-dataset ClusterOverview benchmark (Euclidean Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), STATPC delivers strong accuracy at sub-millisecond speed:
- *Iris* (`Sil = 0.406`, `ARI = 0.818`): Strong external species recovery because botanical classes separate cleanly along axis-aligned petal length and width while overlapping on sepals; STATPC's surprisal automatically selects the discriminating petal subspace.
- *Customers* (`Sil = 0.322`): Well-structured positive separation confirms that wholesale spending channels concentrate in axis-parallel boxes matched to STATPC's scale-invariant model.
- *Multi-Shape* (`Sil = 0.711`, `ARI = N/A`): High Silhouette cohesion reflects STATPC's ability to lock onto statistically dense cores and suppress overlapping sub-boxes, contrasting sharply with linear-subspace SSC (Alg. 150, `Sil = -0.011`).
- *High-D* (`Sil = 0.186`, `ARI = N/A`): Positive cohesion in high ambient dimension validates Proposition 1: STATPC discards uniform noise coordinates, avoiding the full-space distance concentration @beyer1999nearest that degrades global Euclidean partitions @lloyd1982least.
- *Mean Runtime* (`0.193 ms`): Fast 1D interval indexing, bitset intersection, and greedy pruning run $27 times$ faster than convex self-expressive optimization in SSC (`5.253 ms`).
