#import "../template/ieee.typ": *

#alg_entry(
  id: 16,
  name: "Kernel K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "1998 / 2004",
  paradigm: "RKHS Implicit Centroid Partitioning",
  complexity: $O(n^2 (d + T))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n}$ be a finite dataset drawn from an arbitrary input space $Omega$ (not necessarily Euclidean), equipped with positive observation weights $w_1, dots, w_n > 0$ (with standard unweighted case $w_i equiv 1$) and a symmetric *Mercer kernel* $kappa: Omega times Omega -> bb(R)$ inducing the $n times n$ Gram matrix $K_(i ell) = kappa(x_i, x_ell)$. By the Moore–Aronszajn theorem @scholkopf1998nonlinear @bishop2006pattern, whenever $K$ is positive semi-definite ($K succ.eq 0$), there exists a Reproducing Kernel Hilbert Space (RKHS) $(cal(H), chevron.l dot, dot chevron.r_cal(H))$ and nonlinear feature map $phi: Omega -> cal(H)$ such that $kappa(x, y) = chevron.l phi(x), phi(y) chevron.r_cal(H)$.

*Weighted Kernel K-Means* @dhillon2004kernel partitions ${1, dots, n}$ into $k$ disjoint non-empty clusters $cal(C) = {C_1, dots, C_k}$ by minimizing the weighted within-cluster sum of squared RKHS distances:
$ J_kappa (cal(C), bold(m)) = sum_(j=1)^k sum_(i in C_j) w_i norm(phi(x_i) - m_j)_cal(H)^2, $
where $m_j in cal(H)$ are implicit feature-space prototypes. Because $cal(H)$ may be infinite-dimensional (as for the Gaussian RBF kernel $kappa_gamma (x, y) = exp(-gamma norm(x - y)_2^2)$), $phi(x_i)$ and $m_j$ are never materialized explicitly; all distance calculations are performed purely via $K$.

#definition(title: "Implicit RKHS Centroid and Kernel Distance")[
  For a non-empty cluster $C_j subset.eq {1, dots, n}$ with total weight $W_j = sum_(ell in C_j) w_ell$, its *implicit RKHS centroid* is $m_j = W_j^(-1) sum_(ell in C_j) w_ell phi(x_ell) in cal(H)$. Defining the cluster cross-similarity and self-similarity terms
  $ beta_(i j) &:= 1 / W_j sum_(ell in C_j) w_ell K_(i ell), \
    gamma_j    &:= 1 / W_j^2 sum_(p, q in C_j) w_p w_q K_(p q), $
  the squared RKHS distance from $phi(x_i)$ to $m_j$ is
  $ D_kappa^2 (x_i, C_j) &:= norm(phi(x_i) - m_j)_cal(H)^2 \
                         &= K_(i i) - 2 beta_(i j) + gamma_j. $
]

Because $K_(i i)$ does not depend on cluster index $j$, the Voronoi assignment step in $cal(H)$ minimizes $gamma_j - 2 beta_(i j)$ over $j in {1, dots, k}$. Crucially, summing the cost of evaluating $beta_(i j)$ over all $j$ touches each pair $(i, ell)$ *exactly once* per iteration because ${C_1, dots, C_k}$ partitions ${1, dots, n}$, giving $O(n^2)$ work per sweep (Algorithm 1).

#algbox(
  title: "Weighted Kernel K-Means (Implicit Gram Formulation)",
  inputs: [Gram matrix $K in bb(R)^(n times n)$, positive weights $w in bb(R)_(>0)^n$, clusters $k >= 1$, initial partition $cal(C)^((0))$, shift $sigma >= 0$, max iterations $T$.],
  outputs: [Nonlinear partition $cal(C) = {C_1, dots, C_k}$.],
)[
  1. *Optional PSD shift:* If $lambda_"min" (W^(1/2) K W^(1/2)) < 0$, set $K <- K + sigma W^(-1)$ with $sigma >= -lambda_"min" (W^(1/2) K W^(1/2))$, where $W = "diag"(w_1, dots, w_n)$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) *Self-terms:* For each non-empty $C_j^((t))$, compute
     $ W_j     &<- sum_(ell in C_j^((t))) w_ell, \
       gamma_j &<- 1 / W_j^2 sum_(p in C_j^((t))) sum_(q in C_j^((t))) w_p w_q K_(p q). $
  4. #h(1em) *Cross-terms & Assignment:* For each $i = 1, dots, n$:
     - Accumulate $s_(i j) <- sum_(ell in C_j^((t))) w_ell K_(i ell)$ in one scan over $ell = 1, dots, n$, and set $beta_(i j) <- s_(i j) / W_j$.
     - Assign $c_i^((t+1)) <- arg min_(1 <= j <= k) (gamma_j - 2 beta_(i j))$, breaking ties by smallest cluster index $j$.
  5. #h(1em) Form $C_j^((t+1)) = {i : c_i^((t+1)) = j}$ for $j = 1, dots, k$.
  6. #h(1em) *If* $cal(C)^((t+1)) = cal(C)^((t))$, *break*.
  7. *Return* $cal(C)^((t+1))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results governing Kernel K-Means: (1) RKHS parallel-axis optimality and implicit distance identity; (2) monotonic convergence for PSD kernels contrasted with an explicit $3 times 3$ counterexample proving strict objective increase and infinite cycling for indefinite kernels, alongside exact diagonal-shift invariance; (3) exact algebraic equivalence to weighted Normalized Cut and Ratio Cut graph partitioning; (4) exact nonlinear separation of concentric shells that are provably inseparable by Euclidean K-Means; and (5) exact asymptotic degeneracies of the Gaussian RBF kernel at $gamma -> 0^+$ and $gamma -> +infinity$.

#subsec("A", "RKHS Descent, Indefinite-Kernel Divergence, and Shift Invariance")

#lemma(title: "RKHS Parallel-Axis Identity and Gram Expansion", origin: "Original Proof")[
  Let $K succ.eq 0$ with RKHS embedding $phi: Omega -> cal(H)$. For any non-empty cluster $C_j$, positive weights $w_i > 0$, and arbitrary prototype $m in cal(H)$,
  $ &sum_(i in C_j) w_i norm(phi(x_i) - m)_cal(H)^2 \
    &= sum_(i in C_j) w_i norm(phi(x_i) - m_j)_cal(H)^2 + W_j norm(m_j - m)_cal(H)^2, $
  where $m_j = W_j^(-1) sum_(ell in C_j) w_ell phi(x_ell)$ is the unique minimizer and $norm(phi(x_i) - m_j)_cal(H)^2 = K_(i i) - 2 beta_(i j) + gamma_j$.
]

#proof[
  Writing $phi(x_i) - m = (phi(x_i) - m_j) + (m_j - m)$ and expanding the squared Hilbert norm yields
  $ &sum_(i in C_j) w_i norm(phi(x_i) - m)_cal(H)^2 \
    &= sum_(i in C_j) w_i norm(phi(x_i) - m_j)_cal(H)^2 + W_j norm(m_j - m)_cal(H)^2 \
    &quad + 2 chevron.l sum_(i in C_j) w_i (phi(x_i) - m_j), space m_j - m chevron.r_cal(H). $
  By definition of $m_j$, $sum_(i in C_j) w_i (phi(x_i) - m_j) = sum_(i in C_j) w_i phi(x_i) - W_j m_j = 0_cal(H)$, so the inner product vanishes identically. Expanding $norm(phi(x_i) - m_j)_cal(H)^2$ via bilinearity of $chevron.l dot, dot chevron.r_cal(H)$ and $chevron.l phi(x_p), phi(x_q) chevron.r_cal(H) = K_(p q)$ directly yields $K_(i i) - 2 beta_(i j) + gamma_j$.
]

#theorem(title: "PSD Monotonicity, Indefinite Cycling, and Exact Diagonal Shift", origin: "Original Proof")[
  Let $J_K (cal(C)) := "Tr"(W K) - "Tr"(Z^T W^(1/2) K W^(1/2) Z)$ denote the kernel distortion for any symmetric $K in bb(R)^(n times n)$, where $Z_(i j) = w_i^(1/2) W_j^(-1/2) bb(1)_(i in C_j)$. Then:
  1. If $K succ.eq 0$, Algorithm 1 decreases $J_K (cal(C)^((t)))$ monotonically by at least $sum_(j=1)^k W_j^((t)) norm(m_j^((t+1)) - m_j^((t)))_cal(H)^2 >= 0$ per iteration and halts in at most $S(n, k)$ steps.
  2. If $K$ is symmetric *indefinite* ($lambda_"min" (K) < 0$), unshifted Algorithm 1 can *strictly increase* $J_K$ and enter an infinite cycle.
  3. Shifting $tilde(K) = K + sigma W^(-1)$ with $sigma >= -lambda_"min" (W^(1/2) K W^(1/2))$ guarantees $tilde(K) succ.eq 0$ (restoring monotonic convergence) while shifting $J_K (cal(C))$ by the partition-independent constant $sigma (n - k)$:
     $ J_(tilde(K)) (cal(C)) = J_K (cal(C)) + sigma (n - k). $
]

#proof[
  (1) When $K succ.eq 0$, $J_K (cal(C))$ equals $min_(bold(m)) J_kappa (cal(C), bold(m))$ in the induced RKHS $cal(H)$. At iteration $t$, Step 4 assigns each $i$ to minimize $D_kappa^2 (x_i, C_j^((t)))$, guaranteeing $J_kappa (cal(C)^((t+1)), bold(m)^((t))) <= J_kappa (cal(C)^((t)), bold(m)^((t)))$. Applying Lemma 1 to each new cluster $C_j^((t+1))$ with $m = m_j^((t))$ drops $J_kappa$ further by $sum_j W_j^((t+1)) norm(m_j^((t+1)) - m_j^((t)))_cal(H)^2$. Deterministic smallest-index tie-breaking eliminates zero-drop cycles exactly as in Lloyd's proof, ensuring finite termination in $<= S(n, k)$ steps.

  (2) Consider $n = 3$, $k = 2$, $w = (1, 1, 1)^T$, and the symmetric indefinite matrix
  $ K = mat(0, 1, 0; 1, 0, 2; 0, 2, 0), $
  whose eigenvalues are $+sqrt(5), 0, -sqrt(5)$. Initialize $C_1^((0)) = {2, 3}$ and $C_2^((0)) = {1}$. Because $"Tr"(K) = 0$, the initial distortion is $J_K (cal(C)^((0))) = 0 - (2 K_(2 3) / 2 + K_(1 1)) = -2$. Evaluating Step 3 gives $gamma_1 = 4/4 = 1$ and $gamma_2 = 0$. Evaluating scores $gamma_j - 2 beta_(i j)$ in Step 4:
  - $i = 1$: $1 - 2 (1/2) = 0$ (for $C_1$) vs. $0 - 0 = 0$ (for $C_2$) $=> c_1^((1)) = 1$ (smallest-index tie-break).
  - $i = 2$: $1 - 2 (2/2) = -1$ (for $C_1$) vs. $0 - 2 (1) = -2$ (for $C_2$) $=> c_2^((1)) = 2$.
  - $i = 3$: $1 - 2 (2/2) = -1$ (for $C_1$) vs. $0 - 2 (0) = 0$ (for $C_2$) $=> c_3^((1)) = 1$.
  Thus $C_1^((1)) = {1, 3}$ and $C_2^((1)) = {2}$, whose objective is
  $ J_K (cal(C)^((1))) = 0 - (2 K_(1 3) / 2 + K_(2 2)) = 0 > -2 = J_K (cal(C)^((0))), $
  a *strict increase* of $+2$! Running Step 4 from $cal(C)^((1))$ yields $gamma_1 = gamma_2 = 0$ and scores: $i = 1$: $0$ vs. $-2$ ($-> C_2$); $i = 2$: $-3$ vs. $0$ ($-> C_1$); $i = 3$: $0$ vs. $-4$ ($-> C_2$). Hence $C_1^((2)) = {2}, C_2^((2)) = {1, 3}$, and the next iteration swaps back to $C_1^((3)) = {1, 3}, C_2^((3)) = {2}$, cycling infinitely without terminating.

  (3) Observe that $Z^T Z = I_k$ for every $k$-partition with non-empty clusters, and setting $tilde(G) = W^(1/2) K W^(1/2)$ gives $W^(1/2) tilde(K) W^(1/2) = tilde(G) + sigma I_n succ.eq 0$ whenever $sigma >= -lambda_"min" (tilde(G))$. Computing $J_(tilde(K)) (cal(C))$ directly:
  $ J_(tilde(K)) (cal(C)) &= "Tr"(W K + sigma I_n) - "Tr"(Z^T (tilde(G) + sigma I_n) Z) \
                          &= "Tr"(W K) + sigma n - "Tr"(Z^T tilde(G) Z) - sigma "Tr"(I_k) \
                          &= J_K (cal(C)) + sigma (n - k), $
  which differs from $J_K (cal(C))$ solely by the constant $sigma (n - k)$.
]

#subsec("B", "Exact Algebraic Equivalence to Weighted Graph Cuts")

#theorem(title: "Dhillon–Guan–Kulis Graph Cut Equivalence", origin: "Classical — Cited")[
  Let $G = (V, E, A)$ be a weighted undirected graph with symmetric non-negative affinity matrix $A in bb(R)^(n times n)$, positive degree matrix $D = "diag"(d_1, dots, d_n)$ ($d_i = sum_ell A_(i ell) > 0$), and graph Laplacian $L = D - A$ @dhillon2004kernel.
  1. *(Normalized Cut)* Setting $W = D$ and $K = sigma D^(-1) + D^(-1) A D^(-1)$ renders minimizing $J_K (cal(C))$ algebraically identical to minimizing $k$-way Normalized Cut @shi2000ncut:
     $ "Ncut"(cal(C)) = sum_(j=1)^k ("cut"(C_j, overline(C)_j)) / ("vol"(C_j)). $
  2. *(Ratio Cut)* Setting $W = I_n$ and $K = sigma I_n - L$ renders minimizing $J_K (cal(C))$ algebraically identical to minimizing $k$-way Ratio Cut @hagen1992ratiocut:
     $ "Rcut"(cal(C)) = sum_(j=1)^k ("cut"(C_j, overline(C)_j)) / |C_j|. $
]

#proof[
  (1) With $W = D$, the indicator matrix has entries $Z_(i j) = d_i^(1/2) ("vol"(C_j))^(-1/2) bb(1)_(i in C_j)$ where $"vol"(C_j) = sum_(i in C_j) d_i$, satisfying $Z^T Z = I_k$. Since $W^(1/2) K W^(1/2) = sigma I_n + D^(-1/2) A D^(-1/2)$, each diagonal entry of $Z^T D^(-1/2) A D^(-1/2) Z$ equals
  $ (Z^T D^(-1/2) A D^(-1/2) Z)_(j j) &= 1 / ("vol"(C_j)) sum_(p, q in C_j) A_(p q) \
                                      &= ("assoc"(C_j, C_j)) / ("vol"(C_j)) = 1 - ("cut"(C_j, overline(C)_j)) / ("vol"(C_j)), $
  because $"assoc"(C_j, C_j) + "cut"(C_j, overline(C)_j) = "vol"(C_j)$. Summing over $j = 1, dots, k$ yields $"Tr"(Z^T D^(-1/2) A D^(-1/2) Z) = k - "Ncut"(cal(C))$, so
  $ J_K (cal(C)) = "Tr"(D K) - sigma k - k + "Ncut"(cal(C)). $
  Because $"Tr"(D K) - (sigma + 1) k$ is independent of $cal(C)$, minimizing $J_K$ is strictly equivalent to minimizing $"Ncut"$. Choosing $sigma >= -lambda_"min" (D^(-1/2) A D^(-1/2))$ guarantees $K succ.eq 0$ and monotonic descent of $"Ncut"$ without computing eigenvectors!

  (2) With $W = I_n$, $Z_(i j) = |C_j|^(-1/2) bb(1)_(i in C_j)$ and $W^(1/2) K W^(1/2) = sigma I_n - L$. Using the standard Laplacian quadratic identity $sum_(p, q in C_j) L_(p q) = "cut"(C_j, overline(C)_j)$,
  $ "Tr"(Z^T L Z) = sum_(j=1)^k 1 / |C_j| sum_(p, q in C_j) L_(p q) = "Rcut"(cal(C)). $
  Hence $J_K (cal(C)) = "Tr"(K) - sigma k + "Rcut"(cal(C))$, establishing exact algebraic equivalence to Ratio Cut.
]

#subsec("C", "Exact Nonlinear Separation of Concentric Shells")

#proposition(title: "Concentric-Shell Separation in Quadratic RKHS", origin: "Original Proof")[
  For any dimension $d >= 1$ and radii $0 < r_1 < r_2$, define the two concentric cross-polytope shells $S_1 = {plus.minus r_1 e_m : 1 <= m <= d}$ and $S_2 = {plus.minus r_2 e_m : 1 <= m <= d}$ in $bb(R)^d$, with $cal(X) = S_1 union S_2$ ($n = 4 d$).
  1. No Euclidean Voronoi 2-partition (Alg. 1) can separate $S_1$ from $S_2$ for any choice of centroids $mu_1, mu_2 in bb(R)^d$.
  2. Under the quadratic polynomial kernel $kappa(x, y) = (chevron.l x, y chevron.r + 1)^2$, $(S_1, S_2)$ is a *strict stationary RKHS Voronoi equilibrium* of unweighted Kernel K-Means with exact separation margin $(r_2^2 - r_1^2)^2 / d > 0$.
]

#proof[
  (1) Any Euclidean Voronoi 2-partition separates $S_1$ and $S_2$ by a hyperplane ${x : chevron.l v, x chevron.r = b}$ with $v != bold(0)$. Suppose $chevron.l v, x chevron.r <= b$ for all $x in S_1$ and $chevron.l v, y chevron.r >= b$ for all $y in S_2$. Averaging over the symmetric sets $S_1$ and $S_2$ (whose centroids both equal $bold(0)$) gives $0 <= b$ and $0 >= b$, forcing $b = 0$. Then $chevron.l v, plus.minus r_2 e_m chevron.r >= 0$ for every $m = 1, dots, d$, forcing $v_m = 0$ for all $m$, contradicting $v != bold(0)$.

  (2) Take $x = plus.minus r e_a$ ($r in {r_1, r_2}$) and cluster $S_j$ ($j in {1, 2}$, $|S_j| = 2 d$). As $y$ ranges over $S_j$, $chevron.l x, y chevron.r$ equals $+r r_j$ once, $-r r_j$ once, and $0$ for the remaining $2(d - 1)$ points. Therefore,
  $ beta(x, S_j) &= 1 / (2 d) ((r r_j + 1)^2 + (1 - r r_j)^2 + 2(d - 1)) \
                 &= 1 / (2 d) (2 r^2 r_j^2 + 2 d) = 1 + (r^2 r_j^2) / d. $
  Averaging over $x in S_j$ ($r = r_j$) gives $gamma(S_j) = 1 + r_j^4 / d$. Computing the Kernel K-Means assignment score $gamma(S_j) - 2 beta(x, S_j)$ for any point $x$ of radius $r$:
  $ gamma(S_j) - 2 beta(x, S_j) &= (1 + r_j^4 / d) - 2 (1 + (r^2 r_j^2) / d) \
                                &= -1 - r^4 / d + (r_j^2 - r^2)^2 / d. $
  Thus every $x in S_1$ ($r = r_1$) strictly prefers $S_1$ over $S_2$ by margin $(r_2^2 - r_1^2)^2 / d > 0$, and every $x in S_2$ ($r = r_2$) strictly prefers $S_2$ over $S_1$ by the exact same margin, proving $(S_1, S_2)$ is a strict fixed point of Algorithm 1.
]

#subsec("D", "Bandwidth Degeneracy Limits of the Gaussian RBF Kernel")

#proposition(title: "Small-Bandwidth Euclidean Recovery and Large-Bandwidth Freezing", origin: "Original Derivation")[
  Consider unweighted Kernel K-Means on distinct points $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ with Gaussian RBF kernel $kappa_gamma (x, y) = exp(-gamma norm(x - y)_2^2)$ ($gamma > 0$).
  1. As $gamma -> 0^+$, $J_(kappa_gamma) (cal(C)) = 2 gamma J_"Lloyd"^* (cal(C)) + O(gamma^2)$ uniformly over all partitions $cal(C)$, so small $gamma$ recovers Euclidean K-Means (Alg. 1).
  2. As $gamma -> +infinity$, $J_(kappa_gamma) (cal(C)) -> n - k$ for *every* $k$-partition $cal(C)$, and for sufficiently large $gamma$ every non-empty $k$-partition becomes a *strict fixed point* of Algorithm 1 (immediate freezing in one sweep).
]

#proof[
  (1) By Taylor's theorem, $1 - exp(-u) = u + O(u^2)$ uniformly on bounded intervals. Since $norm(phi(x_i) - phi(x_ell))_cal(H)^2 = 2 (1 - kappa_gamma (x_i, x_ell)) = 2 gamma norm(x_i - x_ell)_2^2 + O(gamma^2)$, substituting into the pairwise RKHS inertia representation gives
  $ J_(kappa_gamma) (cal(C)) &= sum_(j=1)^k 1 / (2 |C_j|) sum_(i, ell in C_j) norm(phi(x_i) - phi(x_ell))_cal(H)^2 \
                             &= 2 gamma sum_(j=1)^k 1 / (2 |C_j|) sum_(i, ell in C_j) norm(x_i - x_ell)_2^2 + O(gamma^2) \
                             &= 2 gamma J_"Lloyd"^* (cal(C)) + O(gamma^2). $

  (2) Let $delta_"min" = min_(i != ell) norm(x_i - x_ell)_2 > 0$. Then $K_(i i) = 1$ while $0 < K_(i ell) <= epsilon_gamma := exp(-gamma delta_"min"^2) -> 0$ for all $i != ell$. Hence $K = I_n + E_gamma$ with $norm(E_gamma)_F <= n epsilon_gamma$, giving $J_(kappa_gamma) (cal(C)) = "Tr"(I_n) - "Tr"(H^T I_n H) + O(n epsilon_gamma) = n - k + O(n epsilon_gamma) -> n - k$. Moreover, for any $i in C_j$, $gamma_j = |C_j|^(-1) + O(epsilon_gamma)$ and $beta_(i j) = |C_j|^(-1) + O(epsilon_gamma)$, whereas for any $ell != j$, $beta_(i ell) = O(epsilon_gamma)$. Therefore,
  $ gamma_j - 2 beta_(i j) &= -1 / |C_j| + O(epsilon_gamma) \
                           &< 1 / |C_ell| + O(epsilon_gamma) = gamma_ell - 2 beta_(i ell) $
  for all sufficiently large $gamma$ (specifically once $epsilon_gamma < 1 / (6 n)$), proving every datum strictly retains its initial cluster $C_j$ and freezing any initialization immediately.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Computing the dense $n times n$ Gram matrix $K$ requires $O(n^2 d)$ operations and $O(n^2)$ memory—the dominant bottleneck that limits exact Kernel K-Means to roughly $n <= 2 times 10^4$ on single-node RAM without Nyström or Random Fourier Feature approximations. Once $K$ is materialized, however, each iteration costs only $O(n^2)$ arithmetic operations *independent of $d$* and *independent of $k$*: evaluating $gamma_j$ visits each intra-cluster pair once ($sum_j |C_j|^2 <= n^2$), and accumulating $beta_(i j)$ visits every pair $(i, ell)$ exactly once per sweep. When $K$ is sparse (e.g., $k$-nearest-neighbor graph affinities with $|E|$ edges), per-iteration cost drops dramatically to $O(|E|)$, enabling multilevel graph-cut refinement on millions of nodes @dhillon2004kernel.

Geometrically, linear Voronoi hyperplanes in $cal(H)$ pull back to nonlinear curved hypersurfaces in $Omega$ (quadrics for degree-2 polynomial kernels; arbitrary smooth closed contours for Gaussian RBF kernels). Because centroids exist only implicitly as weighted sums in $cal(H)$, out-of-sample assignment of a new point $x_"new"$ requires $O(n)$ kernel evaluations against the training set via $gamma_j - 2 W_j^(-1) sum_(ell in C_j) w_ell kappa(x_"new", x_ell)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Kernel K-Means],
    [Time / Space Complexity], [Dense: $O(n^2 (d + T))$ time, $O(n^2)$ space; Sparse graph $K$: $O(|E| T)$ time, $O(|E|)$ space.],
    [Geometric Inductive Bias], [Linear Voronoi cells in RKHS $cal(H)$; arbitrary nonlinear boundaries in input space $Omega$.],
    [Convergence Condition], [Monotonic finite convergence iff $K succ.eq 0$ (or shifted $K + sigma W^(-1) succ.eq 0$, Thm. 1).],
    [Graph-Cut Equivalence], [Exact algebraic equivalence to Normalized Cut and Ratio Cut (Theorem 2) @dhillon2004kernel.],
    [Hyperparameter Sensitivity], [Severe RBF bandwidth sensitivity: $gamma -> 0^+$ recovers Lloyd; $gamma -> infinity$ freezes (Prop. 2).],
    [Recommended Companions], [Spectral Clustering (Algs. 103–105) for global eigen-initialization; Kernel FCM (Alg. 154) @zhang2004kfcm.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Kernel K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Kernel K-Means is the method of choice when cluster geometry is intrinsically nonlinear or when objects inhabit structured non-Euclidean domains endowed with domain-specific Mercer kernels.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Non-Euclidean Structured Data (Bioinformatics, Strings, & Graphs):* Protein sequences, DNA motifs, molecular graphs, and parse trees lack canonical vector coordinates, yet admit PSD string, spectrum, and Weisfeiler–Lehman graph kernels @bishop2006pattern. Because Algorithm 1 depends solely on the Gram matrix $K$, it clusters arbitrary structured objects without Euclidean vectorization.
2. *Eigenvector-Free Multilevel Graph Partitioning:* By Theorem 2 @dhillon2004kernel, weighted Kernel K-Means directly optimizes Normalized Cut and Ratio Cut on sparse affinity graphs in $O(|E| T)$ time without solving costly $O(n^3)$ or iterative Lanczos eigenproblems—forming the mathematical engine behind fast multilevel graph partitioners (Graclus) and community-refinement pipelines.
3. *Nonlinear Manifold and Concentric-Structure Segmentation:* In radar signal classification, astrophysical ring/halo separation, and biomedical image segmentation, classes wrap concentrically or along curved manifolds where linear Voronoi partitions fail (Proposition 1). RBF and polynomial kernels unfold these structures into linearly separable RKHS clusters.
4. *Outlier-Damped High-Dimensional Feature Clustering:* Under bounded translation-invariant kernels such as the Gaussian RBF ($0 < kappa_gamma (x, y) <= 1$) or Laplace kernel, every point maps onto the unit hypersphere $norm(phi(x_i))_cal(H) = 1$ in $cal(H)$. Extreme Euclidean outliers cannot pull RKHS centroids arbitrarily far outside the unit ball, naturally damping noise dimensions compared to unbounded $L_2^2$ distortion.
5. *Post-Refinement of Spectral Clustering Embeddings:* Spectral Clustering (Algs. 103–105) @ng2002spectral @vonluxburg2007tutorial solves a continuous eigenvector relaxation of Normalized Cut (Proposition 3 below) whose rounding step may not be locally optimal for discrete $"Ncut"$. Initializing weighted Kernel K-Means with the spectral partition guarantees strict monotonic refinement of exact discrete $"Ncut"$ @dhillon2004kernel.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Kernel K-Means fails in four main regimes: (a) *large dense datasets ($n > 5 times 10^4$)*, where $O(n^2)$ Gram storage exhausts memory (requiring Nyström approximation, Random Fourier Features, or sparse $k$NN graph kernels); (b) *miscalibrated RBF bandwidth $gamma$*, triggering Proposition 2's Euclidean collapse ($gamma -> 0^+$) or instant freezing ($gamma -> +infinity$); (c) *poor local minima under random initialization*, where non-convex RKHS basins trap Algorithm 1 far from global cuts (resolved by spectral initialization or Kernel K-Means++ seeding @arthur2007kmeanspp); and (d) *indefinite similarity matrices* (e.g., sigmoid kernels or cosine-shifted matrices) without diagonal shifting, triggering Theorem 1(2)'s divergence and cycling.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Kernel K-Means Fits / Mechanism], [Caution / Alternative],
    [Bio-Sequences & Graphs], [Needs only PSD Gram $K$; no vectors], [Ensure $K succ.eq 0$ (shift $sigma$)],
    [Sparse Graph Cuts], [Exact $"Ncut"$/$"Rcut"$ in $O(|E| T)$ (Thm. 2)], [Seed via Spectral (Alg. 103)],
    [Concentric / Curved Data], [RKHS hyperplanes unfold rings (Prop. 1)], [Tune RBF $gamma$ carefully],
    [Large Dense Tables ($n > 10^5$)], [$O(n^2)$ Gram memory bottleneck], [Use Hamerly (Alg. 14) / Nyström],
    [Soft Nonlinear Membership], [Hard RKHS Voronoi assignments only], [Use Kernel FCM (Alg. 154)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Kernel K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Kernel K-Means sits at the exact algebraic intersection of Prototype Partitioning (Cat. I), Kernel PCA @scholkopf1998nonlinear, and Spectral Graph Clustering (Cat. V) @ng2002spectral @vonluxburg2007tutorial @chung1997spectral.

#proposition(title: "Spectral Relaxation Bound via Kernel PCA and Stiefel Lifting", origin: "Original Unified Proof")[
  For any PSD Gram matrix $K succ.eq 0$ and positive weight matrix $W = "diag"(w_1, dots, w_n)$, define the weighted Gram operator $tilde(G) = W^(1/2) K W^(1/2)$ with eigenvalues $lambda_1 (tilde(G)) >= dots >= lambda_n (tilde(G)) >= 0$. Then the optimal Kernel K-Means objective satisfies the exact trace identity and spectral lower bound:
  $ min_(cal(C)) J_K (cal(C)) &= "Tr"(W K) - max_(Z in cal(Z)_"ind") "Tr"(Z^T tilde(G) Z) \
                              &>= "Tr"(W K) - sum_(r=1)^k lambda_r (tilde(G)), $
  where $cal(Z)_"ind"$ is the discrete set of weighted cluster indicators ($Z_(i j) = w_i^(1/2) W_j^(-1/2) bb(1)_(i in C_j)$). Relaxing $cal(Z)_"ind"$ to the continuous Stiefel manifold ${Z in bb(R)^(n times k) : Z^T Z = I_k}$ recovers top-$k$ Kernel PCA @scholkopf1998nonlinear (for $W = I_n$) and the normalized spectral clustering relaxation @shi2000ncut @ng2002spectral (for $W = D, K = D^(-1) A D^(-1)$).
]

#proof[
  Let $Phi in bb(R)^(n times p)$ (or Hilbert operator) represent the feature rows $phi(x_i)^T$ so $Phi Phi^T = K$. Observe that $(W^(-1/2) Z Z^T W^(1/2) Phi)_(i,:)$ equals $W_j^(-1) sum_(ell in C_j) w_ell phi(x_ell)^T = m_j^T$ for $i in C_j$. Using $Z^T Z = I_k$, the weighted RKHS distortion expands as
  $ J_K (cal(C)) &= norm(W^(1/2) Phi - Z Z^T W^(1/2) Phi)_F^2 \
                 &= "Tr"(Phi^T W^(1/2) (I_n - Z Z^T) W^(1/2) Phi) \
                 &= "Tr"(W K) - "Tr"(Z^T W^(1/2) K W^(1/2) Z). $
  Because every $Z in cal(Z)_"ind"$ satisfies $Z^T Z = I_k$, maximizing $"Tr"(Z^T tilde(G) Z)$ over the strictly larger Stiefel manifold ${Z : Z^T Z = I_k}$ yields $sum_(r=1)^k lambda_r (tilde(G))$ by the Ky Fan variational theorem, attained by the top-$k$ eigenvectors of $tilde(G)$—which are precisely the Kernel PCA coordinates when $W = I_n$ and the normalized graph eigenvectors when $tilde(G) = D^(-1/2) A D^(-1/2)$.
]

Setting $K = X X^T$ and $W = I_n$ reduces Algorithm 16 verbatim to Lloyd's K-Means (Alg. 1), while relaxing hard indicators to fuzzy memberships yields Kernel Fuzzy C-Means (Alg. 154) @zhang2004kfcm.

*Empirical Behavior on the ClusterOverview Benchmark Suite:* On the 4-dataset empirical benchmark @rousseeuw1987silhouettes @hubert1985comparing, Kernel K-Means delivers the strongest overall cluster-quality profile across the three partitioning algorithms: it achieves the highest *Iris ARI = 0.583* (with *Iris Sil = 0.447*, surpassing Hamerly's 0.563 and K-Medians's 0.440), the highest *Customers Sil = 0.361* (ahead of K-Medians 0.357 and Hamerly 0.341), strong *Multi-Shape Sil = 0.609* (ARI N/A), and strong *High-D Subspace+Noise Sil = 0.303* (ARI N/A, nearly doubling Euclidean Hamerly's 0.158 via bounded unit-sphere RKHS outlier damping). This nonlinear accuracy comes at a clear computational price: mean runtime is *1.052 ms*—over $21 times$ slower than Hamerly's K-Means (0.049 ms) and $7 times$ slower than K-Medians (0.147 ms)—directly manifesting the $O(n^2 (d + T))$ dense Gram-matrix complexity analyzed in Section III.
