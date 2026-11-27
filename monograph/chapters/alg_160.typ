#import "../template/ieee.typ": *

#alg_entry(
  id: 160,
  name: "ISA (Iterative Signature Algorithm)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2003",
  paradigm: "Iterative Thresholded Bipartite Power Iteration",
  complexity: $O(S T n m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The *Iterative Signature Algorithm* (ISA) @bergmann2003isa is a continuous, overlapping biclustering algorithm designed to extract coherent _transcription modules_ — pairs of row and column subsets $(I, J)$ exhibiting mutually reinforcing above-threshold profiles — from large real-valued matrices $E in bb(R)^(n times m)$ without discretizing the data into binary values as in Bimax @prelic2006bimax. Rather than optimizing a single global matrix-partitioning cost, ISA defines biclusters dynamically as self-consistent fixed points of an alternating linear projection and nonlinear z-score thresholding map.

Given a real data matrix $E in bb(R)^(n times m)$ with row index set $R = {1, dots, n}$ and column index set $C = {1, dots, m}$ (having non-constant rows and columns), ISA constructs *two distinct standardized representations* @bergmann2003isa:
$ (E_R)_(i j) = (E_(i j) - bar(E)_(i dot)) / (sigma_(R, i)), quad (E_C)_(i j) = (E_(i j) - bar(E)_(dot j)) / (sigma_(C, j)), $
where $bar(E)_(i dot) = (1 / m) sum_(j=1)^m E_(i j)$ and $sigma_(R, i)^2 = (1 / m) sum_(j=1)^m (E_(i j) - bar(E)_(i dot))^2 > 0$ standardize each row of $E_R$ to zero mean and unit variance, while $bar(E)_(dot j) = (1 / n) sum_(i=1)^n E_(i j)$ and $sigma_(C, j)^2 = (1 / n) sum_(i=1)^n (E_(i j) - bar(E)_(dot j))^2 > 0$ standardize each column of $E_C$. Matrix-form identities read $E_R = D_R^(-1) E P_m$ and $E_C = P_n E D_C^(-1)$, where $P_p = I_p - (1 / p) bold(1)_p bold(1)_p^top$ is the centering projector and $D_R = "diag"(sigma_(R, 1), dots, sigma_(R, n))$, $D_C = "diag"(sigma_(C, 1), dots, sigma_(C, m))$.

For any non-constant vector $v in bb(R)^p$ with empirical mean $mu(v) = (1 / p) bold(1)_p^top v$ and standard deviation $sigma(v) = (1 / sqrt(p)) norm(P_p v)_2 > 0$, define the standardized hard-thresholding operator $cal(T)_t : bb(R)^p -> bb(R)^p$ with threshold $t >= 0$ (one-sided for co-activation, or two-sided for signed co-regulation):
$ (cal(T)_t (v))_ell = cases(
  v_ell & "if" v_ell - mu(v) >= t sigma(v) " (or" |v_ell - mu(v)| >= t sigma(v) ")",
  0 & "otherwise."
) $
Let $hat(cal(T))_t (v) = cal(T)_t (v) / norm(cal(T)_t (v))_2$ whenever $cal(T)_t (v) != bold(0)$, and $bold(0)$ otherwise.

#definition(title: "ISA Transcription Module (Bicluster Fixed Point)")[
  Given dual-standardized matrices $(E_R, E_C)$ and threshold pair $(t_G, t_C) in bb(R)_(>= 0)^2$, define the coupled ISA iteration map $Phi_(t_G, t_C) : bb(R)^n -> bb(R)^n$ by
  $ c^((k)) &= hat(cal(T))_(t_C) (E_R^top g^((k))), \
    g^((k+1)) &= hat(cal(T))_(t_G) (E_C c^((k))) =: Phi_(t_G, t_C) (g^((k))). $
  A non-zero pair $(g^star, c^star) in bb(R)^n times bb(R)^m$ is an _ISA transcription module_ if $c^star = hat(cal(T))_(t_C) (E_R^top g^star)$ and $g^star = hat(cal(T))_(t_G) (E_C c^star)$. Its bicluster support is $(I^star, J^star) = ("supp"(g^star), "supp"(c^star))$, with continuous membership scores given by the non-zero entries of $(g^star, c^star)$.
]

#algbox(
  title: "Iterative Signature Algorithm (Multi-Seed ISA)",
  inputs: [Matrix $E in bb(R)^(n times m)$; thresholds $(t_G, t_C) > 0$; seed count $S$; tolerance $epsilon > 0$; max iterations $T$.],
  outputs: [Deduplicated collection $cal(M) = {(g^star_r, c^star_r)}$ of converged transcription modules.],
)[
  1. *Dual Normalization:* Compute row-standardized $E_R = D_R^(-1) E P_m$ and column-standardized $E_C = P_n E D_C^(-1)$.
  2. *Seed Loop:* For each seed $s = 1, dots, S$:
  3. *Initialization:* Draw a sparse random seed vector $g^((0)) in {0, 1}^n$ (normalized to $norm(g^((0)))_2 = 1$).
  4. *Alternating Signature Refinement:* For $k = 0, dots, T - 1$:
     a) Project onto columns: $tilde(c)^((k)) <- E_R^top g^((k))$ and threshold $c^((k)) <- hat(cal(T))_(t_C) (tilde(c)^((k)))$. If $c^((k)) = bold(0)$, abort seed.
     b) Project onto rows: $tilde(g)^((k+1)) <- E_C c^((k))$ and threshold $g^((k+1)) <- hat(cal(T))_(t_G) (tilde(g)^((k+1)))$. If $g^((k+1)) = bold(0)$, abort seed.
     c) Stop if $norm(g^((k+1)) - g^((k)))_2 < epsilon$ and $norm(c^((k)) - c^((k-1)))_2 < epsilon$.
  5. *Module Deduplication:* Retain non-empty converged pairs $(g^star, c^star)$ and merge pairs whose cosine similarity satisfies $|chevron.l g^star_a, g^star_b chevron.r| > 1 - delta$ and $|chevron.l c^star_a, c^star_b chevron.r| > 1 - delta$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Spectral Reduction in the Zero-Threshold Limit")

When thresholds vanish in the doubly-standardized symmetric regime ($E_R = E_C = M$), ISA collapses exactly to bipartite singular-vector power iteration, providing a clean spectral anchor for understanding how positive thresholds $(t_G, t_C)$ deform the singular spectrum into localized sparse attractors.

#theorem(title: "Exact Reduction of Unthresholded Symmetric ISA to Dominant SVD", origin: "Original Proof")[
  Suppose $E_R = E_C = M in bb(R)^(n times m)$ has singular value decomposition $M = sum_(r=1)^rho sigma_r u_r v_r^top$ with $sigma_1 > sigma_2 >= dots >= sigma_rho > 0$. Under two-sided thresholding with $t_G = t_C = 0$ (so $cal(T)_0 (v) = v$ identically), for any unit seed $g^((0))$ with $chevron.l g^((0)), u_1 chevron.r != 0$, the ISA iterates converge geometrically to the unique rank-$1$ singular bicluster $(plus.minus u_1, plus.minus v_1)$ with exact angular error bounds
  $ sin angle(g^((k)), u_1) <= (sigma_2 / sigma_1)^(2 k) tan angle(g^((0)), u_1), \
    sin angle(c^((k)), v_1) <= (sigma_2 / sigma_1)^(2 k + 1) tan angle(g^((0)), u_1). $
]

#proof[
  When $t_G = t_C = 0$, $c^((k)) = M^top g^((k)) / norm(M^top g^((k)))_2$ and $g^((k+1)) = M c^((k)) / norm(M c^((k)))_2$. By induction on $k$, positive scalar denominators fold into a single normalization factor, yielding
  $ g^((k)) = ((M M^top)^k g^((0))) / norm((M M^top)^k g^((0)))_2, quad c^((k)) = (M^top (M M^top)^k g^((0))) / norm(M^top (M M^top)^k g^((0)))_2. $
  Expand $g^((0)) = sum_(r=1)^rho alpha_r u_r + g_perp$ where $g_perp in ker(M^top)$ and $alpha_1 = chevron.l g^((0)), u_1 chevron.r != 0$. Since $(M M^top)^k u_r = sigma_r^(2 k) u_r$ and $M^top g_perp = bold(0)$, for every $k >= 1$ we obtain
  $ (M M^top)^k g^((0)) = alpha_1 sigma_1^(2 k) u_1 + sum_(r=2)^rho alpha_r sigma_r^(2 k) u_r. $
  Let $theta_k = angle(g^((k)), "span"(u_1)) in [0, pi / 2)$. By orthogonality of ${u_r}$,
  $ tan^2 theta_k = (sum_(r=2)^rho alpha_r^2 sigma_r^(4 k)) / (alpha_1^2 sigma_1^(4 k)) <= (sigma_2 / sigma_1)^(4 k) (sum_(r=2)^rho alpha_r^2 + norm(g_perp)_2^2) / (alpha_1^2) = (sigma_2 / sigma_1)^(4 k) tan^2 theta_0. $
  Since $sin theta_k <= tan theta_k$ for $theta_k in [0, pi / 2)$, taking square roots yields $sin theta_k <= (sigma_2 / sigma_1)^(2 k) tan theta_0$. Similarly, $M^top (M M^top)^k g^((0)) = sum_(r=1)^rho alpha_r sigma_r^(2 k + 1) v_r$, whose orthogonal decomposition along ${v_r}$ gives $tan^2 angle(c^((k)), v_1) <= (sigma_2 / sigma_1)^(4 k + 2) tan^2 theta_0$, completing the proof.
]

#subsec("B", "Variational Lyapunov Monotonicity of Cardinality-Thresholded ISA")

Adaptive z-score thresholds introduce non-variational coupling through $mu(v)$ and $sigma(v)$. However, replacing adaptive z-score truncation by fixed-cardinality hard thresholding $Pi_k (v) = arg max_(norm(u)_2 = 1, norm(u)_0 <= k) chevron.l u, v chevron.r$ (selecting the $k$ largest-magnitude entries and $ell_2$-normalizing) on symmetric $M$ yields an exact Lyapunov ascent system.

#theorem(title: "Lyapunov Ascent and Quadratic Progress of Symmetric Sparse ISA", origin: "Original Proof")[
  Let $M in bb(R)^(n times m)$ and consider the sparse bilinear objective $F(g, c) = g^top M c$ over $cal(S)_G times cal(S)_C$, where $cal(S)_G = {g in bb(R)^n : norm(g)_2 = 1, norm(g)_0 <= k_G}$ and $cal(S)_C = {c in bb(R)^m : norm(c)_2 = 1, norm(c)_0 <= k_C}$. Under alternating updates $c^((k)) = Pi_(k_C) (M^top g^((k)))$ and $g^((k+1)) = Pi_(k_G) (M c^((k)))$ (assuming non-zero projections), the objective is monotonically non-decreasing:
  $ F(g^((k)), c^((k-1))) <= F(g^((k)), c^((k))) <= F(g^((k+1)), c^((k))) <= sigma_1 (M). $
  Moreover, whenever supports stabilize ($"supp"(c^((k-1))) subset.eq "supp"(c^((k)))$), the step achieves exact quadratic progress
  $ F(g^((k)), c^((k))) - F(g^((k)), c^((k-1))) = (L_k / 2) norm(c^((k)) - c^((k-1)))_2^2, $
  where $L_k = F(g^((k)), c^((k))) > 0$, so $F(g^((k)), c^((k)))$ converges to a finite limit $F^star <= sigma_1 (M)$.
]

#proof[
  Fix $g = g^((k)) in cal(S)_G$ and set $v = M^top g$. For any $c in cal(S)_C$ with support $J = "supp"(c)$ ($|J| <= k_C$), Cauchy-Schwarz yields
  $ g^top M c = chevron.l v_J, c_J chevron.r <= norm(v_J)_2 norm(c_J)_2 = norm(v_J)_2, $
  with equality iff $c_J = v_J / norm(v_J)_2$. Maximizing $norm(v_J)_2^2 = sum_(j in J) v_j^2$ over $|J| <= k_C$ selects a subset $J^((k))$ of $k_C$ coordinates maximizing $|v_j|$. Therefore $c^((k)) = Pi_(k_C) (M^top g^((k)))$ is an exact global maximizer of $c |-> F(g^((k)), c)$ over $cal(S)_C$, proving $F(g^((k)), c^((k))) >= F(g^((k)), c^((k-1)))$. By complete symmetry between $(g, M)$ and $(c, M^top)$, $g^((k+1)) = Pi_(k_G) (M c^((k)))$ globally maximizes $g |-> F(g, c^((k)))$ over $cal(S)_G$, yielding $F(g^((k+1)), c^((k))) >= F(g^((k)), c^((k)))$. By definition of the spectral norm, $F(g, c) <= norm(M)_2 norm(g)_2 norm(c)_2 = sigma_1 (M)$ for all $(g, c) in cal(S)_G times cal(S)_C$.

  When $"supp"(c^((k-1))) subset.eq J^((k))$, we have $chevron.l v, c^((k-1)) chevron.r = chevron.l v_(J^((k))), c^((k-1)) chevron.r = L_k chevron.l c^((k)), c^((k-1)) chevron.r$ because $v_(J^((k))) = L_k c^((k))$ with $L_k = norm(v_(J^((k))))_2 = F(g^((k)), c^((k)))$. Expanding $norm(c^((k)) - c^((k-1)))_2^2 = 2 - 2 chevron.l c^((k)), c^((k-1)) chevron.r$ gives $1 - chevron.l c^((k)), c^((k-1)) chevron.r = (1 / 2) norm(c^((k)) - c^((k-1)))_2^2$, hence
  $ F(g^((k)), c^((k))) - F(g^((k)), c^((k-1))) = L_k (1 - chevron.l c^((k)), c^((k-1)) chevron.r) = (L_k / 2) norm(c^((k)) - c^((k-1)))_2^2. $
  Because $F(g^((k)), c^((k)))$ is bounded above and monotonically non-decreasing, it converges to a finite limit $F^star$.
]

#subsec("C", "Non-Normality, Negative Eigenvalues, and Rank-1 Nilpotency")

In the actual ISA algorithm @bergmann2003isa, row standardization and column standardization do *not* commute ($E_R != E_C$). We prove two sharp algebraic obstructions showing why dual standardization destroys positive semidefiniteness and makes nonlinear thresholding structurally indispensable.

#proposition(title: "Rank-1 Nilpotency and Negative Spectrum of Dual Standardization", origin: "Original Counterexample")[
  Let $A = E_C E_R^top in bb(R)^(n times n)$ denote the linear one-step ISA operator induced by dual standardization.
  1. *(Exact Rank-1 Nilpotency):* For any positive-scaled rank-$1$ matrix $E = u v^top$ with $u in bb(R)_(> 0)^n$ non-constant and $v in (bb(R) without {0})^m$ non-constant, $A$ is nilpotent of index $2$: $A^2 = bold(0)_(n times n)$. Consequently, unthresholded ISA ($t_G = t_C = 0$) annihilates a noiseless rank-$1$ multiplicative bicluster in two iterations.
  2. *(Indefinite Spectrum):* Consider the explicit $3 times 3$ family
     $ E(M) = mat(1, 1, -2; M, 0, -M; 0, M^2, -M^2), quad M > 0. $
     As $M -> infinity$, the two non-zero eigenvalues $(lambda_1, lambda_2)$ of $A(M) = E_C (M) E_R (M)^top$ satisfy $lambda_1 lambda_2 -> 9 / 2 - 3 sqrt(3) approx -0.696 < 0$. Hence for all sufficiently large $M$, $A(M)$ possesses a strictly negative real eigenvalue, ruling out any symmetric positive semidefinite quadratic Lyapunov energy $g^top A g$.
]

#proof[
  *Part 1:* Row $i$ of $E = u v^top$ is $u_i v^top$. Since $u_i > 0$, subtracting row mean $u_i bar(v)$ and dividing by row standard deviation $u_i sigma_v$ cancels $u_i$ identically, giving $(E_R)_(i, :) = hat(v)^top$ for all $i$, so $E_R = bold(1)_n hat(v)^top$ with $bold(1)_m^top hat(v) = 0$. Column $j$ of $E$ is $v_j u$; standardizing column $j$ cancels $|v_j|$ while retaining $"sgn"(v_j)$, yielding $E_C = hat(u) "sgn"(v)^top$ with $bold(1)_n^top hat(u) = 0$. Multiplying gives
  $ A = E_C E_R^top = hat(u) ("sgn"(v)^top hat(v)) bold(1)_n^top. $
  Squaring $A$ introduces the inner factor $bold(1)_n^top hat(u) = 0$, hence $A^2 = ("sgn"(v)^top hat(v))^2 hat(u) (bold(1)_n^top hat(u)) bold(1)_n^top = bold(0)_(n times n)$. Thus $g^((2)) prop A g^((1)) prop A^2 g^((0)) = bold(0)$ without thresholding; nonlinear thresholding $cal(T)_(t_G) (hat(u))$ rescues the signal precisely by breaking the zero-sum orthogonality $bold(1)_n^top g^((1)) = 0$.

  *Part 2:* Every row of $E(M)$ has zero mean. Dividing each row by its exact standard deviation ($sqrt(2)$, $M sqrt(2 / 3)$, and $M^2 sqrt(2 / 3)$) yields $M$-independent $E_R$:
  $ E_R = sqrt(3 / 2) mat(1 / sqrt(3), 1 / sqrt(3), -2 / sqrt(3); 1, 0, -1; 0, 1, -1). $
  For columns as $M -> infinity$, Column 1 $(1, M, 0)^top$ is dominated by $+M$ at row 2, Column 2 $(1, 0, M^2)^top$ by $+M^2$ at row 3, and Column 3 $(-2, -M, -M^2)^top$ by $-M^2$ at row 3. Standardizing a $3$-vector with one dominant entry $plus.minus M^p$ gives $plus.minus (2, -1, -1)^top / sqrt(2)$ in the limit, so
  $ lim_(M -> infinity) E_C (M) = (1 / sqrt(2)) mat(-1, -1, 1; 2, -1, 1; -1, 2, -2). $
  Computing $A_infinity = lim_(M -> infinity) E_C (M) E_R^top = (sqrt(3) / 2) W$ gives
  $ W = mat(-4 / sqrt(3), -2, -2; -1 / sqrt(3), 1, -2; 5 / sqrt(3), 1, 4). $
  Since $bold(1)_3^top W = bold(0)^top$, zero is an eigenvalue of $W$, and the product of the remaining two eigenvalues equals the sum of its three $2 times 2$ principal minors:
  $ (-4 / sqrt(3) - 2 / sqrt(3)) + (-16 / sqrt(3) + 10 / sqrt(3)) + (4 + 2) = 6 - 4 sqrt(3). $
  Scaling by $(sqrt(3) / 2)^2 = 3 / 4$ yields $lambda_1 lambda_2 -> (3 / 4)(6 - 4 sqrt(3)) = 9 / 2 - 3 sqrt(3) < 0$. By continuity of eigenvalues in $1 / M$, for all sufficiently large $M$ the two non-zero eigenvalues of $A(M)$ are real of opposite sign ($lambda_1 > 0 > lambda_2$).
]

#subsec("D", "Exact Resolution Capacity Bound and Planted Recovery")

A defining operational feature of ISA is its two-parameter multi-resolution hierarchy controlled by $(t_G, t_C)$ @bergmann2003isa. We derive an exact closed-form upper bound on the maximum fraction of rows and columns a module can occupy at given thresholds.

#theorem(title: "Exact Support-Fraction Capacity Bound of One-Sided ISA", origin: "Classical — Cited")[
  Let $(I^star, J^star)$ be any planted constant-level bicluster of size $n_0 times m_0$ ($0 < n_0 < n$, $0 < m_0 < m$) represented by binary-indicator fixed-point score profiles $g^star prop bold(1)_(I^star)$ and $c^star prop bold(1)_(J^star)$ under one-sided ISA thresholding @bergmann2003isa @prelic2006bimax. A necessary condition for $(I^star, J^star)$ to persist as an ISA fixed point is
  $ n_0 / n <= 1 / (1 + t_G^2), quad m_0 / m <= 1 / (1 + t_C^2). $
  Conversely, for an isolated block $E = mu bold(1)_(I^star) bold(1)_(J^star)^top$ ($mu > 0$) with $E_R, E_C$ standardized cleanly (no background noise), $(I^star, J^star)$ is an exact fixed point if and only if $t_G <= sqrt((n - n_0) / n_0)$ and $t_C <= sqrt((m - m_0) / m_0)$.
]

#proof[
  Let $w in {0, 1}^n$ be the indicator of $I^star$ with $n_0$ ones and $n - n_0$ zeros, and define fraction $p = n_0 / n in (0, 1)$. Its empirical mean and variance are $mu(w) = p$ and $sigma^2 (w) = p(1 - p)$. Therefore every active coordinate $i in I^star$ attains the standardized z-score
  $ z_"on" = (1 - p) / sqrt(p(1 - p)) = sqrt((1 - p) / p) = sqrt((n - n_0) / n_0), $
  whereas every inactive coordinate $i in.not I^star$ has $z_"off" = -sqrt(p / (1 - p)) < 0$. Under one-sided thresholding $z_i >= t_G$ with $t_G >= 0$, inactive rows are always zeroed, while active rows survive iff $t_G <= sqrt((n - n_0) / n_0)$. Squaring both sides yields $t_G^2 <= n / n_0 - 1$, equivalently $n_0 / n <= 1 / (1 + t_G^2)$. Applying the identical calculation to $c^star prop bold(1)_(J^star)$ gives $m_0 / m <= 1 / (1 + t_C^2)$.

  In the noiseless single-block case $E = mu bold(1)_(I^star) bold(1)_(J^star)^top$ (with infinitesimal tie-breaking noise on zero rows/columns, or projecting strictly on non-constant rows/columns), $E_R^top bold(1)_(I^star)$ is proportional to $bold(1)_(J^star)$ up to an additive constant, preserving z-scores $z_"on" = sqrt((m - m_0) / m_0) >= t_C$ on $J^star$ and $z_"off" < 0 < t_C$ off $J^star$. The subsequent row projection $E_C bold(1)_(J^star)$ recovers exact indicator z-scores $sqrt((n - n_0) / n_0) >= t_G$ on $I^star$, closing the fixed-point loop.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each ISA sweep executes two dense matrix-vector products ($E_R^top g$ and $E_C c$) plus two linear-time $O(n)$ and $O(m)$ z-score thresholding passes, costing $O(n m)$ floating-point operations per iteration. Over $S$ random seeds converging in at most $T$ iterations each (empirically $T approx 5"--"20$ when spectral gaps are sharp), total runtime is $O(S T n m)$, plus $O(n m)$ upfront for dual normalization and $O(S^2 (n + m))$ worst-case (or $O(K S (n + m))$ against $K$ unique modules) for cosine deduplication.

Critically, once $g^((k))$ and $c^((k))$ sparsify to supports $|I^((k))| = s_r << n$ and $|J^((k))| = s_c << m$, computing $E_R^top g^((k))$ only touches $s_r$ rows of $E_R$ ($O(s_r m)$), and $E_C c^((k))$ only touches $s_c$ columns of $E_C$ ($O(n s_c)$), reducing steady-state per-iteration cost to $O(s_r m + n s_c)$. Storing $E_R$ in row-major layout and $E_C$ in column-major layout (total memory $2 n m$ floats, i.e., $O(n m)$ space) guarantees contiguous unit-stride memory streaming in both half-steps, maximizing SIMD vectorization and L2/L3 cache line utilization, and trivially parallelizing across all $S$ seeds via embarrassingly parallel batched BLAS-3 GEMM ($E_R^top G$, $E_C C$).

Geometrically, ISA searches for sparse local singular-vector basins on the product of unit spheres $bb(S)^(n-1) times bb(S)^(m-1)$. By Theorem 4, $(t_G, t_C)$ sweeps a nested hierarchy of resolutions: low thresholds merge modules into broad global singular modes, while high thresholds carve out tight, strongly co-regulated micro-modules obeying $n_0 / n <= 1 / (1 + t_G^2)$. Dual standardization ensures exact invariance to positive row-wise affine transforms ($E_(i j) |-> a_i E_(i j) + b_i$ in $E_R$) and column-wise affine transforms ($E_(i j) |-> c_j E_(i j) + d_j$ in $E_C$), though *not* joint invariance to arbitrary heterogeneous rescaling (Proposition 1).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Input Domain], [Real continuous matrices $E in bb(R)^(n times m)$],
    [Time Complexity], [$O(S T n m)$ dense; $O(S T (s_r m + n s_c))$ sparse steady state],
    [Space Complexity], [$O(n m)$ (row-major $E_R$ + column-major $E_C$)],
    [Optimality], [Local fixed points; exact SVD at $t_G = t_C = 0$ (symmetric)],
    [Overlap Structure], [Continuous overlapping scores; multi-module membership],
    [Resolution Control], [Capacity $n_0 / n <= (1 + t_G^2)^(-1)$, $m_0 / m <= (1 + t_C^2)^(-1)$],
    [Convergence Behavior], [Monotone in symmetric fixed-$k$; cyclic/nilpotent risks if $E_R != E_C$],
    [Key Hyperparameters], [Thresholds $(t_G, t_C)$; random seed count $S$; one/two-sided mode],
  ),
  caption: [Operational and complexity profile of ISA @bergmann2003isa.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Transcriptional Regulatory Module Discovery:* Introduced specifically for yeast and human gene expression compendia @bergmann2003isa, ISA captures pleiotropic genes belonging to multiple overlapping pathways and ranks genes/conditions continuously by signature score. Sweeping $(t_G, t_C)$ naturally uncovers hierarchical biological organization from coarse cell-cycle programs down to specific regulons.
2. *Single-Cell Multi-Omics Co-Expression Modules:* In high-dimensional scRNA-seq matrices ($n approx 20\,000$ genes, $m approx 10^4"--"10^5$ cells), ISA's $O(s_r m + n s_c)$ sparse inner iterations scale linearly in $n$ and $m$, avoiding $O(n^2)$ pairwise gene-gene correlation matrices while tolerating heavy dropout noise far better than binary Bimax @prelic2006bimax.
3. *Financial Sector Co-Movement Regimes:* Applied to asset-return matrices (assets $times$ trading windows), two-sided ISA identifies transient regimes where subsets of equities move synchronously (positive and anti-correlated hedges) during macro shocks, weighting each asset continuously by its factor alignment strength.
4. *Recommender Systems and Latent Taste Niches:* On user-item rating matrices, ISA discovers overlapping micro-communities of users sharing strong preference signatures across niche item subsets without forcing every user into a single global cluster or requiring expensive MCMC inference.
5. *Neuroimaging Functional Subnetwork Extraction:* In voxel-by-time or ROI-by-task fMRI matrices, dual standardization equalizes baseline BOLD amplitude across brain regions while iterative z-thresholding isolates spatially sparse, task-locked functional assemblies.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four primary failure modes constrain ISA deployment: (i) *low-dimensional threshold starvation* (Theorem 4): when $m$ is tiny (e.g., $m = 4$), $m_0 / m <= 1 / (1 + t_C^2)$ leaves almost no room for multi-feature modules at standard thresholds $t_C >= 1$, collapsing modules to $1$ feature or empty sets; (ii) *asymmetric non-convergence and rank-1 collapse* (Proposition 1): severe cross-row scale disparity causes $E_C E_R^top$ to lose symmetry and positive semidefiniteness, creating limit cycles or zero-collapses; (iii) *seed coverage bias*: small modules with narrow basins of attraction can be missed if $S$ is too small, while dominant global modules absorb most random seeds; and (iv) *lack of explicit noise-variance modeling*, where heteroscedastic high-variance rows distort $E_C$ unless pre-filtered or modeled probabilistically via FABIA @hochreiter2010fabia.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Large Gene Expression Compendia], [Continuous overlapping modules; fast $O(S T n m)$ sweep], [Tune $(t_G, t_C)$ grid; deduplicate seeds],
    [Low-Feature Tabular ($m <= 6$)], [Poor; capacity bound $m_0 / m <= (1 + t_C^2)^(-1)$ starves columns], [Use $k$-Means, GMM, or spectral methods],
    [Heteroscedastic Multiplicative Data], [Moderate; dual z-scoring ignores row noise $psi_i$], [Prefer FABIA @hochreiter2010fabia],
    [Exact Binary Incidence Data], [Approximate; continuous averaging blurs exact $1$-blocks], [Use Bimax @prelic2006bimax for exactness],
    [Multi-Resolution Exploratory Mining], [Strong; $(t_G, t_C)$ traces nested attractor hierarchy], [Guard against asymmetric limit cycles],
  ),
  caption: [Domain suitability and failure-mode matrix for ISA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ISA bridges spectral biclustering @kluger2003spectral, bipartite co-clustering @dhillon2001concept, and sparse matrix factorization @lee1999nmf. The following proposition formalizes its precise algebraic relationship to Kluger's spectral biclustering and rank-$1$ sparse non-negative factorization.

#proposition(title: "Spectral Normalization and Sparse NMF Correspondence of ISA", origin: "Original Proof")[
  Let $E in bb(R)^(n times m)$ have strictly positive entries.
  1. *(Spectral Biclustering Limit):* Replace ISA's one-pass non-commuting standardization $(E_R, E_C)$ with Sinkhorn-Knopp bistochastic normalization $tilde(E) = D_1^(-1/2) E D_2^(-1/2)$ @kluger2003spectral @dhillon2001concept (so $E_R = E_C = tilde(E)$). Then unthresholded ISA ($t_G = t_C = 0$) recovers the exact leading singular pair $(u_1, v_1)$ of Kluger's normalized checkerboard operator, and projecting orthogonal to $(u_1, v_1)$ recovers the leading bicluster eigenvectors $(u_2, v_2)$.
  2. *(Rank-$1$ Cardinality-Constrained NMF Equivalence):* For $M >= 0$ entrywise and one-sided positive cardinality projection $Pi^+_k (v)$ onto ${u >= 0 : norm(u)_2 = 1, norm(u)_0 <= k}$, every fixed point $(g^star, c^star)$ of symmetric ISA with $sigma = (g^star)^top M c^star > 0$ is a block-coordinate stationary point of the sparse rank-$1$ Frobenius NMF objective @lee1999nmf
     $ min_(g >= 0, c >= 0, alpha >= 0) norm(M - alpha g c^top)_F^2 quad "s.t." norm(g)_2 = norm(c)_2 = 1, norm(g)_0 <= k_G, norm(c)_0 <= k_C. $
]

#proof[
  For (1), setting $E_R = E_C = tilde(E)$ satisfies the hypothesis of Theorem 1; thus unthresholded ISA performs exact power iteration on $tilde(E) tilde(E)^top$, converging to $(u_1, v_1)$ whenever $chevron.l g^((0)), u_1 chevron.r != 0$, and within $u_1^perp$ converging to $(u_2, v_2)$ at rate $(sigma_3 / sigma_2)^2$. For (2), expanding the Frobenius norm gives
  $ norm(M - alpha g c^top)_F^2 = norm(M)_F^2 - 2 alpha g^top M c + alpha^2 norm(g)_2^2 norm(c)_2^2. $
  Minimizing over $alpha >= 0$ at unit-norm $(g, c)$ yields $alpha^star = max(0, g^top M c) = g^top M c$ (since $M, g, c >= 0$), substituting back to give $norm(M)_F^2 - (g^top M c)^2$. Hence minimizing Frobenius error is strictly equivalent to maximizing $g^top M c$ over the non-negative sparse unit spheres. By Theorem 2 restricted to the non-negative orthant (where $v = M^top g >= 0$ automatically), the coordinate optimality conditions coincide verbatim with the fixed-point equations $c^star = Pi^+_(k_C) (M^top g^star)$ and $g^star = Pi^+_(k_G) (M c^star)$.
]

Unlike FABIA (Alg. 161 @hochreiter2010fabia), which infers $p$ multiplicative biclusters jointly via variational EM, ISA extracts modules sequentially from independent seeds; this grants ISA an order-of-magnitude speed advantage at the cost of no global probabilistic model for overlapping residual subtraction.

*Empirical Behavior on ClusterOverview Benchmarks:* Across the four-dataset ClusterOverview benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, ISA runs in a mean of *0.401 ms*, reflecting the lightweight $O(n m)$ matrix-vector sweeps and rapid convergence. On *Iris* (`Sil = N/A`, `ARI = 0.028`), the exact capacity bound of Theorem 4 directly explains the collapse: with only $m = 4$ features, any two-feature subset ($m_0 = 2$) achieves $z_"on" = 1.0$ and any three-feature subset ($m_0 = 3$) achieves $z_"on" = 1 / sqrt(3) approx 0.577$, so standard z-score thresholding starves column selection and collapses rows into degenerate or near-random labels (`ARI = 0.028`, `Sil = N/A`). Conversely, on *Multi-Shape* (`Sil = 0.773`, `ARI = N/A`), strong directional co-variation in the standardized feature space aligns sharply with ISA's bilinear fixed point, yielding the highest Silhouette score (`0.773`) of all four algorithms (Alg. 159–162). On *Customers* (`Sil = 0.135`) and *High-D* (`Sil = 0.075`, `ARI = N/A`), ISA successfully extracts non-degenerate partitions across coherent feature subsets, though global Euclidean Silhouette remains low because ISA optimizes localized bilinear coherence on $I^star times J^star$ rather than isotropic full-space compactness.
