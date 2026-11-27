#import "../template/ieee.typ": *

#alg_entry(
  id: 162,
  name: "k-Shape",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2015",
  paradigm: "Cross-Correlation Shift-Invariant Time-Series Partitioning",
  complexity: $O(T (n k m log m + n m^2 + k m^3))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard partitional clustering algorithms such as Lloyd's $k$-Means @lloyd1982least compare objects via coordinate-wise Euclidean distance, which fails catastrophically on temporal signals exhibiting amplitude scaling, baseline drift, or phase shifts. Conversely, elastic alignment distances such as Dynamic Time Warping (DTW) @sakoe1978dtw and DTW Barycenter Averaging (DBA) @petitjean2011dba incur $O(m^2)$ pairwise cost and lack closed-form global barycenters. *k-Shape* @paparrizos2015kshape bridges this gap via a shift- and scale-invariant partitional framework grounded in normalized cross-correlation and spectral Rayleigh-quotient centroid extraction.

Let $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^m$ be $n$ non-constant real time series of length $m >= 2$. Define the orthogonal centering projector $Q = I_m - (1 / m) bold(1)_m bold(1)_m^top$ ($Q^2 = Q = Q^top$, $"rank"(Q) = m - 1$). Each series is preprocessed via unit-norm z-normalization:
$ hat(bold(x))_i = (Q bold(x)_i) / norm(Q bold(x)_i)_2 in bb(S)^(m-1) inter bold(1)_m^perp, $
so $bold(1)_m^top hat(bold(x))_i = 0$ and $norm(hat(bold(x))_i)_2 = 1$ (equivalently standard z-scoring scaled by $1 / sqrt(m)$). For any shift lag $s in {-(m-1), dots, m-1}$, define the zero-padded shift operator $T_s : bb(R)^m -> bb(R)^m$ by $(T_s (bold(y)))_t = y_(t - s)$ if $1 <= t - s <= m$ and $0$ otherwise. Note the exact adjoint identity $chevron.l bold(x), T_s (bold(y)) chevron.r = chevron.l T_(-s) (bold(x)), bold(y) chevron.r$.

#definition(title: "Normalized Cross-Correlation and Shape-Based Distance (SBD)")[
  For non-zero series $bold(x), bold(y) in bb(R)^m$ and lag $s in {-(m-1), dots, m-1}$, define the linear cross-correlation $"CC"_s (bold(x), bold(y)) = chevron.l bold(x), T_s (bold(y)) chevron.r$ and its coefficient-normalized counterpart
  $ "NCC"_c (bold(x), bold(y))_s = ("CC"_s (bold(x), bold(y))) / (norm(bold(x))_2 norm(bold(y))_2). $
  The _Shape-Based Distance_ (SBD) @paparrizos2015kshape is defined by
  $ "SBD"(bold(x), bold(y)) = 1 - max_(-(m-1) <= s <= m-1) "NCC"_c (bold(x), bold(y))_s, $
  with optimal shift $s^star (bold(x), bold(y)) in arg max_s "NCC"_c (bold(x), bold(y))_s$ aligning $bold(y)$ onto $bold(x)$ via $T_(s^star) (bold(y))$.
]

By the Discrete Convolution Theorem, zero-padding $bold(x), bold(y)$ to length $L = 2^(ceil(log_2 (2 m - 1)))$ computes the entire cross-correlation vector ${"CC"_s (bold(x), bold(y))}_(s=-(m-1))^(m-1)$ simultaneously in $O(m log m)$ via $cal(F)_L^(-1) (cal(F)_L (bold(x)) dot.o overline(cal(F)_L (bold(y))))$, where $cal(F)_L$ denotes the $L$-point Fast Fourier Transform (FFT). Given a cluster $C_j$ whose members are shift-aligned to the previous centroid $bold(mu)_j$ as $tilde(bold(x))_i = T_(s_i^star) (hat(bold(x))_i)$, $k$-Shape extracts the updated z-normalized centroid by maximizing the sum of squared correlations:
$ bold(mu)_j^star = arg max_(bold(mu) in.not "span"(bold(1)_m)) sum_(i in C_j) (chevron.l Q bold(mu), tilde(bold(x))_i chevron.r / norm(Q bold(mu))_2)^2. $

#algbox(
  title: "k-Shape Time-Series Clustering Algorithm",
  inputs: [Time series $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^m$; cluster count $k$; max iterations $T$.],
  outputs: [Hard partition ${C_1, dots, C_k}$ and z-normalized shape centroids ${bold(mu)_1, dots, bold(mu)_k}$.],
)[
  1. *Z-Normalization:* Normalize each $hat(bold(x))_i <- (Q bold(x)_i) / norm(Q bold(x)_i)_2$ and precompute its zero-padded FFT $cal(F)_L (hat(bold(x))_i)$.
  2. *Initialization:* Assign each $hat(bold(x))_i$ randomly to an initial cluster $C_j^((0))$ (or seed $k$ centroids) and run Step 4 once.
  3. *Assignment & Alignment Step:* For each series $i = 1, dots, n$ and cluster $j = 1, dots, k$, evaluate $"SBD"(bold(mu)_j, hat(bold(x))_i)$ and optimal shift $s_(i j)^star$ in $O(m log m)$ via inverse FFT. Assign $i$ to $C_(j^star)$ where $j^star = arg min_j "SBD"(bold(mu)_j, hat(bold(x))_i)$ and record aligned series $tilde(bold(x))_i <- T_(s_(i j^star)^star) (hat(bold(x))_i)$.
  4. *Shape Extraction Step:* For each non-empty cluster $C_j$:
     a) Form scatter matrix $S_j <- sum_(i in C_j) tilde(bold(x))_i tilde(bold(x))_i^top$ and centered Rayleigh matrix $M_j <- Q S_j Q$.
     b) Compute the leading unit eigenvector $bold(v)_1$ of $M_j$.
     c) Resolve sign ambiguity: set $bold(mu)_j <- "sgn"(sum_(i in C_j) chevron.l bold(v)_1, tilde(bold(x))_i chevron.r) bold(v)_1$ and re-normalize $bold(mu)_j <- (Q bold(mu)_j) / norm(Q bold(mu)_j)_2$.
  5. *Termination:* Repeat Steps 3–4 until cluster assignments remain unchanged or $t = T$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Sharp Range, Affine Invariance, and Non-Metric Geometry of SBD")

Although the literature routinely states $"SBD"(bold(x), bold(y)) in [0, 2]$ by bounding each individual correlation in $[-1, 1]$ @paparrizos2015kshape, we prove a strictly sharper universal upper bound on z-normalized time series and construct an explicit minimal counterexample demonstrating that SBD violates the triangle inequality.

#proposition(title: "Sharp [0, 1) Range on Z-Normalized Series and Triangle Inequality Violation", origin: "Original Counterexample")[
  Let $hat(bold(x)), hat(bold(y)) in bb(S)^(m-1) inter bold(1)_m^perp$ be any two z-normalized non-constant time series ($m >= 2$). Then:
  1. *(Strict Unit Upper Bound):* $0 <= "SBD"(hat(bold(x)), hat(bold(y))) < 1$ strictly, and $"SBD"(hat(bold(x)), hat(bold(y))) = 0$ if and only if $hat(bold(x)) = T_s (hat(bold(y)))$ for some shift $s$ whose shifted-out boundary coordinates of $hat(bold(y))$ are identically zero.
  2. *(Positive-Affine Invariance):* For all $alpha_1, alpha_2 > 0$ and $beta_1, beta_2 in bb(R)$, z-normalizing $alpha_1 bold(x) + beta_1 bold(1)_m$ and $alpha_2 bold(y) + beta_2 bold(1)_m$ leaves $"SBD"$ strictly invariant.
  3. *(Non-Metric Counterexample):* In $bb(R)^4$, the three z-normalized unit vectors
     $ bold(a) = 1 / 2 (1, -1, 1, -1)^top, quad bold(b) = 1 / sqrt(2) (1, 0, 0, -1)^top, quad bold(c) = 1 / 2 (1, 1, -1, -1)^top $
     strictly violate the triangle inequality:
     $ "SBD"(bold(a), bold(b)) + "SBD"(bold(b), bold(c)) = 2 - sqrt(2) approx 0.5858 < 0.75 = "SBD"(bold(a), bold(c)). $
]

#proof[
  *Part 1:* By Cauchy-Schwarz, $"CC"_s (hat(bold(x)), hat(bold(y))) = chevron.l hat(bold(x)), T_s (hat(bold(y))) chevron.r <= norm(hat(bold(x)))_2 norm(T_s (hat(bold(y))))_2 <= 1$, giving $"SBD" >= 0$, with equality iff $norm(T_s (hat(bold(y))))_2 = 1$ (all shifted-out entries of $hat(bold(y))$ vanish) and $hat(bold(x)) = T_s (hat(bold(y)))$. For the upper bound, sum $"CC"_s$ across all $2 m - 1$ lags:
  $ sum_(s = -(m-1))^(m-1) "CC"_s (hat(bold(x)), hat(bold(y))) = sum_(t=1)^m sum_(u=1)^m hat(x)_t hat(y)_u = (bold(1)_m^top hat(bold(x))) (bold(1)_m^top hat(bold(y))) = 0. $
  Encode $P_(hat(bold(x))) (z) = sum_(t=1)^m hat(x)_t z^t$ and $Q_(hat(bold(y))) (z) = sum_(u=1)^m hat(y)_u z^(m - u)$ in the polynomial ring $bb(R)[z]$. The coefficient of $z^(m + s)$ in $P_(hat(bold(x))) (z) Q_(hat(bold(y))) (z)$ equals $"CC"_s (hat(bold(x)), hat(bold(y)))$. Since $hat(bold(x)), hat(bold(y)) != bold(0)$ and $bb(R)[z]$ is an integral domain, $P_(hat(bold(x))) Q_(hat(bold(y))) != 0$, so at least one lag satisfies $"CC"_s != 0$. Combined with $sum_s "CC"_s = 0$, at least one lag must be strictly positive, proving $max_s "NCC"_c (hat(bold(x)), hat(bold(y)))_s > 0$ and hence $"SBD"(hat(bold(x)), hat(bold(y))) < 1$.

  *Part 2:* Since $Q bold(1)_m = bold(0)$ and $alpha_1 > 0$, $Q(alpha_1 bold(x) + beta_1 bold(1)_m) / norm(Q(alpha_1 bold(x) + beta_1 bold(1)_m))_2 = (Q bold(x)) / norm(Q bold(x))_2 = hat(bold(x))$, leaving SBD unchanged.

  *Part 3:* Direct verification confirms $bold(1)_4^top bold(a) = bold(1)_4^top bold(b) = bold(1)_4^top bold(c) = 0$ and $norm(bold(a))_2 = norm(bold(b))_2 = norm(bold(c))_2 = 1$. Evaluating $4 "CC"_s (bold(a), bold(c))$ across all $7$ lags $s in {-3, dots, 3\}$ yields $(-1, 0, 1, 0, 1, 0, -1)$, whose maximum is $+1$ at $s = plus.minus 1$; thus $"SBD"(bold(a), bold(c)) = 1 - 1 / 4 = 3 / 4$. Evaluating $2 sqrt(2) "CC"_s (bold(a), bold(b))$ yields $(-1, 1, -1, 2, -1, 1, -1)$ (maximum $2$ at $s = 0$), and $2 sqrt(2) "CC"_s (bold(b), bold(c))$ yields $(-1, -1, 1, 2, 1, -1, -1)$ (maximum $2$ at $s = 0$). Hence $"SBD"(bold(a), bold(b)) = "SBD"(bold(b), bold(c)) = 1 - 1 / sqrt(2)$, and their sum $2 - sqrt(2)$ is strictly less than $3 / 4$ because $sqrt(2) > 5 / 4 <==> 32 > 25$.
]

#subsec("B", "Spectral Optimality and Automatic Zero-Mean Property of Shape Extraction")

A key mathematical insight of $k$-Shape is that the z-normalization constraint $bold(1)_m^top bold(mu) = 0$ does not require constrained Lagrange-multiplier projection during centroid extraction: it is enforced automatically by the nullspace geometry of $M_j = Q S_j Q$.

#theorem(title: "Rayleigh Optimality and Automatic Zero-Mean Eigenvector Property", origin: "Classical — Cited")[
  Let $S_j = sum_(i in C_j) tilde(bold(x))_i tilde(bold(x))_i^top succ.eq 0$ be the aligned scatter matrix of cluster $C_j$, and set $M_j = Q S_j Q$ @paparrizos2015kshape @vonluxburg2007tutorial. Assume $M_j != bold(0)$. Then:
  1. $max_(bold(mu) in.not "span"(bold(1)_m)) sum_(i in C_j) (chevron.l Q bold(mu), tilde(bold(x))_i chevron.r / norm(Q bold(mu))_2)^2 = lambda_"max" (M_j) > 0$.
  2. Every leading unit eigenvector $bold(v)_1$ of $M_j$ automatically satisfies $bold(1)_m^top bold(v)_1 = 0$ and $Q bold(v)_1 = bold(v)_1$, solving the constrained optimization exactly.
  3. Among the two optimal unit signs ${+bold(v)_1, -bold(v)_1}$, $bold(mu)_j = "sgn"(sum_(i in C_j) chevron.l bold(v)_1, tilde(bold(x))_i chevron.r) bold(v)_1$ uniquely minimizes $sum_(i in C_j) norm(tilde(bold(x))_i - bold(mu))_2^2$.
]

#proof[
  Set $bold(w) = Q bold(mu) / norm(Q bold(mu))_2$. Because $Q = Q^top = Q^2$, $bold(w)$ ranges over all unit vectors in $"range"(Q) = bold(1)_m^perp$, and
  $ sum_(i in C_j) (chevron.l Q bold(mu), tilde(bold(x))_i chevron.r / norm(Q bold(mu))_2)^2 = bold(w)^top S_j bold(w) = bold(w)^top (Q S_j Q) bold(w) = bold(w)^top M_j bold(w). $
  Because $M_j = Q S_j Q$ is symmetric positive semidefinite and non-zero, $lambda_"max" (M_j) > 0$. Furthermore, $M_j bold(1)_m = Q S_j (Q bold(1)_m) = bold(0)$, so $bold(1)_m$ lies in $ker(M_j)$. For any unit eigenvector $bold(v)_1$ satisfying $M_j bold(v)_1 = lambda_"max" (M_j) bold(v)_1$, left-multiplying by $bold(1)_m^top$ gives
  $ lambda_"max" (M_j) (bold(1)_m^top bold(v)_1) = (M_j bold(1)_m)^top bold(v)_1 = 0 ==> bold(1)_m^top bold(v)_1 = 0. $
  Thus $bold(v)_1 in bold(1)_m^perp$ ($Q bold(v)_1 = bold(v)_1$), so the unconstrained Rayleigh maximum $max_(norm(bold(v))_2 = 1) bold(v)^top M_j bold(v) = lambda_"max" (M_j)$ is attained inside $bold(1)_m^perp$, proving (1) and (2). For (3), expanding $sum_(i in C_j) norm(tilde(bold(x))_i - s bold(v)_1)_2^2 = sum_(i in C_j) (norm(tilde(bold(x))_i)_2^2 + 1) - 2 s sum_(i in C_j) chevron.l tilde(bold(x))_i, bold(v)_1 chevron.r$ for $s in {+1, -1}$ shows the minimum occurs at $s = "sgn"(sum_(i in C_j) chevron.l tilde(bold(x))_i, bold(v)_1 chevron.r)$.
]

#subsec("C", "Alternating Surrogate Monotonicity and Linear-SBD Sandwich Bound")

Because $k$-Shape assigns series by linear SBD ($1 - r_i$) while extracting centroids via squared correlation ($sum r_i^2$), we establish exact finite-step monotonicity for the squared-correlation objective and derive a tight sandwich inequality linking $lambda_"max" (M_j)$ directly to the linear SBD distortion.

#theorem(title: "Squared-Correlation Monotonicity and SBD Sandwich Inequality", origin: "Original Proof")[
  Define the squared-correlation objective over partitions ${C_j}$, shift lags ${s_i}$, and centroids $bold(mu)_j in bb(S)^(m-1) inter bold(1)_m^perp$:
  $ cal(J)_"sq" ({C_j}, {s_i}, {bold(mu)_j}) = sum_(j=1)^k sum_(i in C_j) chevron.l T_(s_i) (hat(bold(x))_i), bold(mu)_j chevron.r^2. $
  1. Under squared-correlation assignment $(j(i), s_i) = arg max_(j, s) chevron.l T_s (hat(bold(x))_i), bold(mu)_j chevron.r^2$ alternating with Shape Extraction (Step 4), $cal(J)_"sq"$ is monotonically non-decreasing at every half-step and terminates in finitely many iterations.
  2. For any cluster $C_j$ whose aligned correlations satisfy $r_i = chevron.l tilde(bold(x))_i, bold(mu)_j chevron.r in [0, 1]$ at $bold(mu)_j = bold(v)_1 (M_j)$, the linear SBD sum is tightly bracketed by $lambda_"max" (M_j)$:
     $ |C_j| - sqrt(|C_j| lambda_"max" (M_j)) <= sum_(i in C_j) (1 - r_i) <= |C_j| - lambda_"max" (M_j). $
]

#proof[
  *Part 1:* Fixing ${bold(mu)_j}$, independently selecting $(j(i), s_i) = arg max_(j, s) chevron.l T_s (hat(bold(x))_i), bold(mu)_j chevron.r^2$ for each $i$ maximizes every summand of $cal(J)_"sq"$ globally, so $cal(J)_"sq"$ cannot decrease. Fixing ${C_j, s_i}$, Theorem 1 guarantees that $bold(mu)_j = bold(v)_1 (M_j)$ globally maximizes $sum_(i in C_j) chevron.l tilde(bold(x))_i, bold(mu)_j chevron.r^2$ over $bb(S)^(m-1) inter bold(1)_m^perp$ for every $j$, so $cal(J)_"sq"$ again cannot decrease. Because the discrete configuration space of partition-shift pairs $({C_j}, {s_i})$ has finite cardinality at most $(k (2 m - 1))^n$ and deterministic tie-breaking prevents revisiting any state of strictly higher $cal(J)_"sq"$, the iteration terminates in finitely many steps.

  *Part 2:* Let $S_1 = sum_(i in C_j) r_i$ and $S_2 = sum_(i in C_j) r_i^2 = bold(v)_1^top M_j bold(v)_1 = lambda_"max" (M_j)$. Since $0 <= r_i <= 1$ implies $r_i^2 <= r_i$, summing over $i in C_j$ gives $S_2 <= S_1$. Conversely, Cauchy-Schwarz yields $S_1^2 <= |C_j| S_2$, hence $S_1 <= sqrt(|C_j| S_2)$. Substituting $sum_(i in C_j) (1 - r_i) = |C_j| - S_1$ into $S_2 <= S_1 <= sqrt(|C_j| S_2)$ directly yields $|C_j| - sqrt(|C_j| lambda_"max" (M_j)) <= sum_(i in C_j) (1 - r_i) <= |C_j| - lambda_"max" (M_j)$. Both inequalities are tight at perfect alignment ($r_i equiv 1$, where $lambda_"max" (M_j) = |C_j|$ and all three expressions vanish).
]

#subsec("D", "Spectral Perturbation Stability of Shape Extraction")

We quantify how robust $k$-Shape's spectral centroid is to additive waveform perturbations by proving an elementary self-contained Davis-Kahan $\sin Theta$ angular error bound @davis1970rotation.

#theorem(title: "Perturbation Stability of k-Shape Centroids", origin: "Original Derivation")[
  Suppose cluster $C_j$ of size $N = |C_j|$ consists of noisy aligned copies $tilde(bold(x))_i = bold(mu)^star + bold(eta)_i$, where $bold(mu)^star in bb(S)^(m-1) inter bold(1)_m^perp$ is the clean unit template. Let $M_0 = N bold(mu)^star (bold(mu)^star)^top$ and $M_j = Q (sum_(i=1)^N tilde(bold(x))_i tilde(bold(x))_i^top) Q = M_0 + E$. Define the normalized perturbation norm $rho = norm(E)_2 / N$. Whenever $rho < 1$, the extracted centroid $hat(bold(mu))_j = bold(v)_1 (M_j)$ satisfies
  $ sin angle(hat(bold(mu))_j, bold(mu)^star) <= rho / (1 - rho), $
  where $rho <= (2 / N) norm(sum_(i=1)^N Q bold(eta)_i)_2 + (1 / N) sum_(i=1)^N norm(Q bold(eta)_i)_2^2$.
]

#proof[
  Expanding $M_j = Q (sum_(i=1)^N (bold(mu)^star + bold(eta)_i)(bold(mu)^star + bold(eta)_i)^top) Q$ and using $Q bold(mu)^star = bold(mu)^star$ gives $M_j = M_0 + E$, where
  $ E = bold(mu)^star (sum_(i=1)^N Q bold(eta)_i)^top + (sum_(i=1)^N Q bold(eta)_i) (bold(mu)^star)^top + sum_(i=1)^N (Q bold(eta)_i)(Q bold(eta)_i)^top. $
  By the triangle inequality and $norm(bold(u) bold(v)^top)_2 = norm(bold(u))_2 norm(bold(v))_2$, we obtain $norm(E)_2 <= 2 norm(sum_i Q bold(eta)_i)_2 + sum_i norm(Q bold(eta)_i)_2^2$, verifying the upper bound on $rho = norm(E)_2 / N$.

  Now decompose the unit leading eigenvector $hat(bold(mu))_j$ of $M_j$ (with eigenvalue $hat(lambda)_1$) as $hat(bold(mu))_j = (cos theta) bold(mu)^star + (sin theta) bold(w)$, where $bold(w) perp bold(mu)^star$, $norm(bold(w))_2 = 1$, and $theta = angle(hat(bold(mu))_j, bold(mu)^star) in [0, pi / 2]$. Because $M_0 = N bold(mu)^star (bold(mu)^star)^top$ satisfies $bold(w)^top M_0 = bold(0)^top$, multiplying $M_j hat(bold(mu))_j = hat(lambda)_1 hat(bold(mu))_j$ on the left by $bold(w)^top$ yields
  $ hat(lambda)_1 sin theta = bold(w)^top (M_0 + E) hat(bold(mu))_j = bold(w)^top E hat(bold(mu))_j <= norm(E)_2. $
  Meanwhile, by Rayleigh-Ritz optimality of $hat(lambda)_1$,
  $ hat(lambda)_1 >= (bold(mu)^star)^top (M_0 + E) bold(mu)^star = N + (bold(mu)^star)^top E bold(mu)^star >= N - norm(E)_2 = N (1 - rho) > 0. $
  Dividing both sides by $hat(lambda)_1 >= N (1 - rho)$ yields $sin theta <= norm(E)_2 / (N - norm(E)_2) = rho / (1 - rho)$, matching the rank-$1$ Davis-Kahan perturbation rate @davis1970rotation.
]

Note that for zero-mean independent noise $bold(eta)_i$, the linear cross-term $(2 / N) norm(sum_i Q bold(eta)_i)_2$ concentrates at rate $O(sigma / sqrt(N))$, while for isotropic noise $bb(E)[(Q bold(eta)_i)(Q bold(eta)_i)^top] = sigma^2 Q$ adds only a multiple of the identity on $bold(1)_m^perp$ — which commutes with $M_0$ and leaves the leading eigenvector $bold(mu)^star$ _completely unperturbed_ in expectation!

#sec("III", "Complexity, Geometry and Operational Profile")

In each iteration of Algorithm 162, computing SBD between $n$ series and $k$ centroids requires $k$ forward FFTs of the centroids ($O(k m log m)$) plus $n k$ pointwise frequency-domain multiplications and inverse FFTs of length $L = Theta(m)$ ($O(n k m log m)$ total, since the $n$ data FFTs $cal(F)_L (hat(bold(x))_i)$ are precomputed once in $O(n m log m)$). During Shape Extraction, explicitly materializing $S_j = sum_(i in C_j) tilde(bold(x))_i tilde(bold(x))_i^top$ costs $O(n m^2)$ across all clusters, and dense symmetric eigendecomposition of the $k$ matrices $M_j in bb(R)^(m times m)$ costs $O(k m^3)$, giving total complexity $O(T(n k m log m + n m^2 + k m^3))$ @paparrizos2015kshape. Crucially, for long series ($m >= 512$), replacing dense eigendecomposition by $r$ matrix-free power iterations $bold(v) |-> Q (sum_(i in C_j) tilde(bold(x))_i chevron.l tilde(bold(x))_i, Q bold(v) chevron.r)$ avoids materializing $S_j$ entirely and reduces centroid extraction from $O(n m^2 + k m^3)$ to $O(r n m)$, making $O(T n k m log m)$ the sole bottleneck.

Space complexity is $O(n m + k m + m^2)$ (or $O(n m + k m)$ matrix-free), keeping FFT buffers cache-resident for typical lengths $m <= 4096$. Geometrically, $k$-Shape operates on the quotient of the centered unit sphere $bb(S)^(m-1) inter bold(1)_m^perp$ under rigid zero-padded translations $T_s$, conferring exact positive-affine invariance (Proposition 1) and global phase-shift tolerance without requiring any domain-specific distance hyperparameter (unlike DTW's Sakoe-Chiba warping window @sakoe1978dtw or Soft-DTW's smoothing temperature @cuturi2017softdtw).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Input Domain], [Equal-length ordered real time series $cal(X) subset bb(R)^m$],
    [Time Complexity], [$O(T(n k m log m + n m^2 + k m^3))$; $O(T n k m log m)$ power-iter],
    [Space Complexity], [$O(n m + k m + m^2)$ dense; $O(n m + k m)$ matrix-free],
    [Distance Geometry], [Non-metric SBD $in [0, 1)$ on z-normalized series (Prop. 1)],
    [Invariances], [Exact positive-affine ($alpha x + beta$) and rigid lag-shift],
    [Centroid Extraction], [Leading eigenvector of $M_j = Q S_j Q$ (automatically zero-mean)],
    [Noise Robustness], [Davis-Kahan $O(rho / (1 - rho))$; immune to isotropic noise mean],
    [Key Hyperparameters], [Cluster count $k$ only (parameter-free distance)],
  ),
  caption: [Operational and complexity profile of $k$-Shape @paparrizos2015kshape.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Biomedical Waveform Clustering (ECG, EEG, and PPG):* Cardiac cycles and neural evoked potentials vary markedly in sensor gain, baseline drift, and trigger-onset phase lag across patients. $k$-Shape's exact positive-affine invariance and rigid lag alignment cleanly group morphological arrhythmias without distorting diagnostic QRS durations via excessive elastic warping.
2. *Smart-Grid Load Profiling and Building Energy Analytics:* Daily household electricity consumption curves share characteristic morning/evening load shapes shifted by 1–3 hours depending on occupant schedules and scaled by household size. SBD's shift-invariant correlation recovers behavioral load archetypes while FFT acceleration scales easily to millions of smart meters.
3. *Industrial IoT and Rotating-Machinery Vibration Monitoring:* Periodic vibration transients and acoustic emission bursts sampled asynchronously from bearings or turbines exhibit arbitrary trigger phases; $k$-Shape extracts clean spectral fault templates immune to isotropic white sensor noise (Theorem 3).
4. *Financial Intraday Pattern and Regime Discovery:* Normalized intraday price and volume trajectories across equities exhibit delayed leader-follower responses; cross-correlation matching captures lead-lag shape synchronization across assets of disparate price scales.
5. *Astronomy Light-Curve Classification:* Variable stars and exoplanet transit light curves suffer from unknown orbital phase offsets and baseline flux differences; $k$-Shape groups phase-shifted photometric morphologies orders of magnitude faster than DTW/DBA @petitjean2011dba while maintaining comparable UCR benchmark accuracy @paparrizos2015kshape.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes govern $k$-Shape's domain boundaries: (i) *unordered tabular features*: shifting feature index $t$ onto index $t - s$ aligns semantically unrelated attributes (e.g., age onto income), injecting spurious similarity; (ii) *nonlinear time-axis stretching (variable speed)*: rigid global shift $T_s$ cannot align signals whose sub-events accelerate or decelerate locally, where DTW/DBA @sakoe1978dtw @petitjean2011dba or ROCKET representations @dempster2020rocket excel; (iii) *amplitude-driven discrimination*: z-normalization erases absolute magnitude and baseline offsets, failing when clusters differ purely by mean level or variance (e.g., constant low vs. high temperature sensors); and (iv) *boundary-clipping distortion under large shifts*, where zero-padding $T_s$ truncates large portions of the signal energy whenever $|s| approx m / 2$, shrinking $|"NCC"_c|$ artificially.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Phase-Shifted Waveforms (ECG, IoT)], [Exact affine + rigid shift invariance; fast $O(m log m)$ FFT], [Avoid huge lags $|s| > m / 2$],
    [Variable-Speed / Elastic Motion], [Weak; rigid shift $T_s$ cannot stretch local segments], [Use DTW/DBA @petitjean2011dba or Soft-DTW @cuturi2017softdtw],
    [Magnitude-Sensitive Signals], [Poor; z-normalization erases scale and DC offset], [Use Euclidean $k$-Means @lloyd1982least],
    [Unordered Static Tabular Data], [Misaligned; cross-feature shifts mix unrelated columns], [Use $k$-Means++ @arthur2007kmeanspp or GMM],
    [Large-Scale Time-Series Archives], [Strong; parameter-free SBD + matrix-free power iter], [Use multi-start seeds for local optima],
  ),
  caption: [Domain suitability and failure-mode matrix for $k$-Shape.],
)

#sec("V", "Structural Relations and Empirical Behavior")

$k$-Shape occupies a crisp structural position between Spherical $k$-Means @dhillon2001concept @banerjee2005vmf and elastic time-series barycenter clustering @petitjean2011dba.

#proposition(title: "Exact Zero-Lag Reduction to Spherical k-Means and Second-Order Barycenter Contrast", origin: "Original Proof")[
  Restrict the admissible lag set of $k$-Shape to $s = 0$ on unit-norm z-normalized series $hat(bold(x))_i in bb(S)^(m-1) inter bold(1)_m^perp$. Then:
  1. The assignment rule $arg min_j "SBD"_0 (bold(mu)_j, hat(bold(x))_i)$ is identically equivalent to both Spherical $k$-Means cosine assignment $arg max_j chevron.l bold(mu)_j, hat(bold(x))_i chevron.r$ @dhillon2001concept @banerjee2005vmf and Euclidean $k$-Means assignment $arg min_j norm(hat(bold(x))_i - bold(mu)_j)_2^2$ @lloyd1982least.
  2. Whereas Spherical $k$-Means updates $bold(mu)_j$ via the first-order Fréchet mean $bold(mu)_j^"sph" = (sum_(i in C_j) hat(bold(x))_i) / norm(sum_(i in C_j) hat(bold(x))_i)_2$, $k$-Shape replaces it with the second-order principal direction $bold(v)_1 (sum_(i in C_j) hat(bold(x))_i hat(bold(x))_i^top)$, which strictly dominates the first-order mean in quadratic coherence: $sum_(i in C_j) chevron.l bold(v)_1, hat(bold(x))_i chevron.r^2 >= sum_(i in C_j) chevron.l bold(mu)_j^"sph", hat(bold(x))_i chevron.r^2$.
]

#proof[
  For (1), when $s = 0$ and $norm(hat(bold(x))_i)_2 = norm(bold(mu)_j)_2 = 1$, $"SBD"_0 (bold(mu)_j, hat(bold(x))_i) = 1 - chevron.l bold(mu)_j, hat(bold(x))_i chevron.r = (1 / 2) norm(hat(bold(x))_i - bold(mu)_j)_2^2$. Minimizing SBD over $j$ therefore coincides identically with maximizing inner product (Spherical $k$-Means) and minimizing squared Euclidean distance ($k$-Means). For (2), because every $hat(bold(x))_i in bold(1)_m^perp$, $Q hat(bold(x))_i = hat(bold(x))_i$, so $M_j = sum_(i in C_j) hat(bold(x))_i hat(bold(x))_i^top$ and $bold(mu)_j^"sph" in bb(S)^(m-1) inter bold(1)_m^perp$. By Theorem 1, $bold(v)_1$ maximizes $bold(v)^top M_j bold(v)$ globally over all unit vectors, immediately yielding $bold(v)_1^top M_j bold(v)_1 >= (bold(mu)_j^"sph")^top M_j bold(mu)_j^"sph"$, with strict inequality whenever the cluster members do not share identical projections onto $bold(mu)_j^"sph"$.
]

*Empirical Behavior on ClusterOverview Benchmarks:* On the standardized four-dataset ClusterOverview benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, $k$-Shape records a mean runtime of *14.685 ms* — the highest among Algorithms 159–162 — directly reflecting the $O(n k m log m)$ FFT cross-correlations and $O(k m^3)$ Rayleigh eigensolves repeated over multiple iterations. Unlike the three biclustering algorithms (Bimax, ISA, FABIA), whose overlapping or thresholded formulations frequently produce `N/A` global Silhouette scores on low-dimensional static datasets, $k$-Shape's Lloyd-style alternating partitioning guarantees a complete, non-degenerate hard clustering across all four benchmarks. On *Iris* (`Sil = 0.237`, `ARI = 0.397`), treating the four ordered morphological measurements as a length-$4$ shape profile captures relative sepal-to-petal proportions cleanly enough to achieve `ARI = 0.397` and `Sil = 0.237` (the highest Iris Silhouette in this group). On *Customers* (`Sil = 0.181`) and *Multi-Shape* (`Sil = 0.200`, `ARI = N/A`), positive-affine z-normalization discards raw spending magnitudes and spatial coordinates, while cross-feature lag shifts $s != 0$ mix non-temporal dimensions, depressing Euclidean Silhouette relative to specialized geometric algorithms. Similarly, on *High-D* (`Sil = 0.055`, `ARI = N/A`), applying rigid 1D sequence shifts across unordered high-dimensional features blurs coordinate-specific cluster support, illustrating why $k$-Shape achieves peak accuracy specifically on genuinely sequential, phase-shifted time-series domains @paparrizos2015kshape.
