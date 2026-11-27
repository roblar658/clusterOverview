#import "../template/ieee.typ": *

#alg_entry(
  id: 29,
  name: "Minkowski Weighted K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2012",
  paradigm: "Cluster-Specific Weighted Lp Partitioning",
  complexity: $O(n d (k + T_"in") T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset with entries $x_(i v)$ ($v in {1, dots, d}$), $k in {2, dots, n}$ the target cluster count, and $p > 1$ a fixed *Minkowski exponent*. Standard Euclidean $k$-means (Alg. 1) weights all features identically, causing distance concentration when uninformative noise features proliferate @beyer1999nearest. Whereas W-K-Means (Alg. 28) @huang2005automated learns one global weight vector with an independent exponent $beta$, *Minkowski Weighted K-Means (MWK-Means)* @deamorim2012minkowski introduces two structural generalizations: (i) each cluster $C_j$ learns its own *local* weight vector $w_j = (w_(j 1), dots, w_(j d))^T$ on the probability simplex $Delta_(d-1) = { w in bb(R)_+^d : sum_(v=1)^d w_v = 1 }$, and (ii) the *same* exponent $p$ governs both distance and weight powers.

Thus the weighted Minkowski dissimilarity between $x_i$ and prototype $mu_j in bb(R)^d$ is an exact $L_p$ metric on cluster-locally rescaled coordinates:
$ d_p (x_i, mu_j; w_j) = ( sum_(v=1)^d w_(j v)^p |x_(i v) - mu_(j v)|^p )^(1/p) = norm( "diag"(w_j) (x_i - mu_j) )_p. $
Let $cal(C) = {C_1, dots, C_k}$ partition ${1, dots, n}$, $bold(mu) = (mu_1, dots, mu_k)$, $W = [w_(j v)] in (Delta_(d-1))^k$, and define the within-cluster $p$-dispersion $D_(j v) (C_j, mu_j) := sum_(i in C_j) |x_(i v) - mu_(j v)|^p$. Adding a non-negative dispersion floor $epsilon_0 >= 0$ (matching the dispersion shift in @deamorim2012minkowski) to prevent singular collapse on constant features yields the regularized *MWK-Means variational objective*:
$ tilde(J)_p (cal(C), bold(mu), W) &= sum_(j=1)^k sum_(i in C_j) sum_(v=1)^d w_(j v)^p |x_(i v) - mu_(j v)|^p + epsilon_0 sum_(j=1)^k norm(w_j)_p^p \
  &= sum_(j=1)^k sum_(v=1)^d w_(j v)^p ( D_(j v) (C_j, mu_j) + epsilon_0 ), $
where $epsilon_0 = 0$ recovers the unregularized objective $J_p (cal(C), bold(mu), W)$ of @deamorim2012minkowski.

#definition(title: "Weighted Minkowski Stationary Partition")[
  For $p > 1$ and $epsilon_0 >= 0$, a triple $(cal(C), bold(mu), W)$ is a *coordinate-wise MWK-Means equilibrium* if: (i) every $i in C_j$ satisfies $d_p (x_i, mu_j; w_j)^p <= d_p (x_i, mu_ell; w_ell)^p$ for all $ell != j$; (ii) each $mu_(j v)$ minimizes $c |-> sum_(i in C_j) |x_(i v) - c|^p$; and (iii) each $w_j in Delta_(d-1)$ minimizes $w |-> sum_(v=1)^d w_v^p (D_(j v) + epsilon_0)$.
]

MWK-Means minimizes $tilde(J)_p$ via exact three-block alternating minimization (Algorithm 1). In practice, de Amorim and Mirkin @deamorim2012minkowski recommend pairing Algorithm 1 with *Anomalous Pattern* (`iMWK-Means`) initialization, which extracts candidate clusters farthest from the data center before applying weighted refinement.

#algbox(
  title: "Minkowski Weighted K-Means (MWK-Means)",
  inputs: [Standardized dataset $cal(X) subset.eq bb(R)^d$, clusters $k >= 2$, exponent $p > 1$, floor $epsilon_0 > 0$, initial $(bold(mu)^((0)), W^((0)))$, max iterations $T$.],
  outputs: [Partition $cal(C)$, Minkowski centers $bold(mu)$, cluster-specific feature weights $W in (Delta_(d-1))^k$.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) *Assignment step:* For each $i in {1, dots, n}$, compute
     $ c_i^((t)) = arg min_(j in {1, dots, k}) sum_(v=1)^d (w_(j v)^((t)))^p |x_(i v) - mu_(j v)^((t))|^p, $
     breaking ties by smallest index $j$; set $C_j^((t)) = {i : c_i^((t)) = j}$.
  3. #h(1em) *Minkowski center update:* For each non-empty $C_j^((t))$ and $v in {1, dots, d}$, solve the 1D convex problem
     $ mu_(j v)^((t+1)) = arg min_(c in bb(R)) sum_(i in C_j^((t))) |x_(i v) - c|^p $
     (mean for $p = 2$; 1D Newton/steepest descent for $p != 2$).
  4. #h(1em) *Dispersion & simplex weight update:* Compute $tilde(D)_(j v)^((t+1)) = epsilon_0 + sum_(i in C_j^((t))) |x_(i v) - mu_(j v)^((t+1))|^p$ and set
     $ w_(j v)^((t+1)) = (tilde(D)_(j v)^((t+1)))^(-1/(p-1)) / ( sum_(u=1)^d (tilde(D)_(j u)^((t+1)))^(-1/(p-1)) ). $
  5. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$ (for $t >= 1$), *break*.
  6. *Return* $(cal(C)^((t)), bold(mu)^((t+1)), W^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five theoretical guarantees: (A) exact Hölder simplex optimality; (B) uniqueness of the 1D $L_p$ center and finite Lyapunov convergence; (C) noise-feature suppression rates; (D) scale collapse and non-convex Voronoi cells; and (E) outlier leverage decay and NP-hardness.

#subsec("A", "Hölder Simplex Optimality and Power-Mean Reduction")

#lemma(title: "Exact Hölder Minimization over the Weight Simplex", origin: "Original Proof")[
  Fix $p > 1$ and positive regularized dispersions $tilde(D)_(j v) > 0$ for $v = 1, dots, d$. Over the closed simplex $Delta_(d-1)$, the unique global minimizer of $w_j |-> sum_(v=1)^d w_(j v)^p tilde(D)_(j v)$ is
  $ w_(j v)^* = tilde(D)_(j v)^(-1/(p-1)) / ( sum_(u=1)^d tilde(D)_(j u)^(-1/(p-1)) ), $
  yielding the exact reduced cluster distortion
  $ min_(w_j in Delta_(d-1)) sum_(v=1)^d w_(j v)^p tilde(D)_(j v) = ( sum_(v=1)^d tilde(D)_(j v)^(-1/(p-1)) )^(1 - p). $
  Moreover, if $epsilon_0 = 0$ and $Z_j := {v : D_(j v) = 0} != emptyset$, the minimum is $0$ attained by any $w_j in Delta_(d-1)$ with $"supp"(w_j) subset.eq Z_j$; at the boundary $p = 1$, every simplex vertex $e_(v^*)$ with $v^* in arg min_v tilde(D)_(j v)$ is optimal.
]

#proof[
  Let $q = p / (p - 1)$ be the Hölder conjugate of $p$, so $1/p + 1/q = 1$, $q/p = 1/(p-1)$, and $p/q = p - 1$. For any $w_j in Delta_(d-1)$, set $a_v = w_(j v) tilde(D)_(j v)^(1/p)$ and $b_v = tilde(D)_(j v)^(-1/p)$. By Hölder's inequality,
  $ 1 = sum_(v=1)^d w_(j v) &= sum_(v=1)^d a_v b_v \
    &<= ( sum_(v=1)^d w_(j v)^p tilde(D)_(j v) )^(1/p) ( sum_(v=1)^d tilde(D)_(j v)^(-1/(p-1)) )^(1/q). $
  Raising both sides to the power $p > 1$ and dividing by $(sum_v tilde(D)_(j v)^(-1/(p-1)))^(p-1)$ establishes the lower bound
  $ sum_(v=1)^d w_(j v)^p tilde(D)_(j v) >= ( sum_(v=1)^d tilde(D)_(j v)^(-1/(p-1)) )^(1 - p). $
  By the strict equality condition of Hölder's inequality for $p > 1$, equality holds if and only if $a_v^p = kappa b_v^q$ for a positive constant $kappa > 0$ independent of $v$, i.e., $w_(j v)^p tilde(D)_(j v) = kappa tilde(D)_(j v)^(-1/(p-1))$, equivalently $w_(j v) prop tilde(D)_(j v)^(-1/(p-1))$. Normalizing to $sum_(v=1)^d w_(j v) = 1$ uniquely yields $w_(j v)^*$.

  When $epsilon_0 = 0$ and $Z_j != emptyset$, every $w_j in Delta_(d-1)$ supported inside $Z_j$ achieves $sum_v w_(j v)^p D_(j v) = 0$, showing why zero dispersion collapses weight onto constant axes unless $epsilon_0 > 0$. At $p = 1$, $w_j |-> sum_v w_(j v) tilde(D)_(j v)$ is linear on $Delta_(d-1)$ and attains $min_v tilde(D)_(j v)$ at extreme vertices $e_(v^*)$, so smooth multi-feature weighting requires $p > 1$. Notably, at $p = 2$, Lemma 1 yields the scaled *harmonic mean* $(sum_v tilde(D)_(j v)^(-1))^(-1) = d^(-1) H(tilde(D)_(j 1), dots, tilde(D)_(j d))$, which is dominated by the smallest coordinate variances and encodes axis-aligned subspace discovery.
]

#subsec("B", "Minkowski Center Uniqueness and Lyapunov Convergence")

#lemma(title: "Uniqueness & Monotone Root of the 1D Lp Center", origin: "Original Proof")[
  For any non-empty cluster $C_j$, feature $v$, and exponent $p > 1$, the univariate loss $g_(j v) (c) = sum_(i in C_j) |x_(i v) - c|^p$ is strictly convex and $C^1$-smooth with strictly increasing derivative
  $ g'_(j v) (c) = p sum_(i in C_j) |c - x_(i v)|^(p-1) "sign"(c - x_(i v)). $
  Hence $g_(j v)$ admits a *unique* global minimizer $mu_(j v)^* in [min_(i in C_j) x_(i v), max_(i in C_j) x_(i v)]$, independent of any positive weight $w_(j v) > 0$.
]

#proof[
  Since $p > 1$, the map $psi(t) = |t|^p$ is $C^1$ on $bb(R)$ with derivative $psi'(t) = p |t|^(p-1) "sign"(t)$, which is strictly increasing everywhere because $t |-> t^(p-1)$ is strictly increasing on $[0, infinity)$. For any $c_1 < c_2$ and every $i in C_j$, $c_1 - x_(i v) < c_2 - x_(i v)$ implies $psi'(c_1 - x_(i v)) < psi'(c_2 - x_(i v))$. Summing over $i in C_j$ ($|C_j| >= 1$) shows $g'_(j v) (c_1) < g'_(j v) (c_2)$, proving strict convexity of $g_(j v)$.

  Let $a = min_(i in C_j) x_(i v)$ and $b = max_(i in C_j) x_(i v)$. Every summand of $g'_(j v) (a)$ is non-positive, so $g'_(j v) (a) <= 0$, whereas every summand of $g'_(j v) (b)$ is non-negative, so $g'_(j v) (b) >= 0$. By the Intermediate Value Theorem and strict monotonicity of $g'_(j v)$, there exists a unique root $mu_(j v)^* in [a, b]$ solving $g'_(j v) (mu_(j v)^*) = 0$, which is the unique global minimizer. Multiplying $g_(j v)$ by any $w_(j v)^p > 0$ scales the gradient by a positive constant without altering its root.
]

#theorem(title: "Exact Three-Block Lyapunov Descent & Finite Termination", origin: "Original Proof")[
  For $p > 1$ and any dispersion floor $epsilon_0 > 0$, let ${(cal(C)^((t)), bold(mu)^((t)), W^((t)))}$ be generated by Algorithm 1 with deterministic tie-breaking and non-empty clusters. Then:
  1. At every iteration $t >= 0$, the regularized objective satisfies
     $ tilde(J)_p (cal(C)^((t+1)), bold(mu)^((t+1)), W^((t+1))) &<= tilde(J)_p (cal(C)^((t)), bold(mu)^((t+1)), W^((t+1))) \
       &<= tilde(J)_p (cal(C)^((t)), bold(mu)^((t+1)), W^((t))) \
       &<= tilde(J)_p (cal(C)^((t)), bold(mu)^((t)), W^((t))). $
  2. For any fixed partition $cal(C)$, the pair $(bold(mu)^*(cal(C)), W^*(cal(C)))$ computed by Steps 3–4 is the *unique* joint global minimizer of $(bold(mu), W) |-> tilde(J)_p (cal(C), bold(mu), W)$, achieving the reduced partition energy
     $ Phi_p (cal(C)) = sum_(j=1)^k ( sum_(v=1)^d ( D_(j v)^*(C_j) + epsilon_0 )^(-1/(p-1)) )^(1 - p). $
  3. Algorithm 1 terminates in at most $S(n, k) + 1 <= k^n / k! + 1$ iterations at a coordinate-wise MWK-Means equilibrium.
]

#proof[
  At iteration $t$, hold $cal(C)^((t))$ and $W^((t))$ fixed. Because $epsilon_0 sum_(j,v) (w_(j v)^((t)))^p$ is independent of $bold(mu)$ and $w_(j v)^((t)) > 0$, minimizing $tilde(J)_p$ over $bold(mu)$ decouples into $k d$ problems $min_c g_(j v) (c)$, uniquely solved by $bold(mu)^((t+1))$ (Lemma 2), giving the third inequality of Claim 1. Holding $(cal(C)^((t)), bold(mu)^((t+1)))$ fixed, $tilde(D)_(j v)^((t+1)) >= epsilon_0 > 0$, so Lemma 1 shows $W^((t+1))$ uniquely minimizes $tilde(J)_p$ over $W$, yielding the second inequality. Because $bold(mu)^((t+1))$ depends *only* on $cal(C)^((t))$ (Lemma 2), $(bold(mu)^((t+1)), W^((t+1)))$ is the *unique joint global minimizer* of $tilde(J)_p (cal(C)^((t)), dot, dot)$, establishing Claim 2.

  In Step 2 of iteration $t + 1$, nearest-cluster assignment minimizes the data term while leaving $epsilon_0 sum_(j,v) (w_(j v)^((t+1)))^p$ fixed, giving the first inequality of Claim 1. Chaining Claims 1–2 yields $Phi_p (cal(C)^((t+1))) <= tilde(J)_p (cal(C)^((t+1)), bold(mu)^*(cal(C)^((t))), W^*(cal(C)^((t)))) <= Phi_p (cal(C)^((t)))$. If equality holds, $(bold(mu)^*(cal(C)^((t))), W^*(cal(C)^((t))))$ uniquely minimizes $tilde(J)_p (cal(C)^((t+1)), dot, dot)$, forcing $(bold(mu)^((t+2)), W^((t+2))) = (bold(mu)^((t+1)), W^((t+1)))$; deterministic tie-breaking then yields $cal(C)^((t+2)) = cal(C)^((t+1))$ and immediate termination. Thus partitions strictly decrease $Phi_p$ prior to halting and cannot repeat; finiteness of $S(n, k)$ concludes the proof.
]

#subsec("C", "Quantitative Noise-Feature Suppression Rates")

#theorem(title: "Exact Feature-Share Identity & Noise Suppression Bound", origin: "Original Proof")[
  Let $p > 1$, $epsilon_0 = 0$, and suppose cluster $C_j$ has $D_(j v) > 0$ for all $v$. Partition ${1, dots, d} = cal(S)_j union.sq cal(N)_j$ into $d_s = |cal(S)_j|$ informative features with $D_(j v) <= sigma_s^p |C_j|$ and $d_n = |cal(N)_j|$ noise features with $D_(j u) >= sigma_n^p |C_j|$, where the scale-contrast ratio satisfies $rho := sigma_n / sigma_s > 1$. Then at the optimal weights $w_j^*$ of Lemma 1:
  1. *Distortion-Share Identity:* Each feature's fractional contribution to cluster $C_j$'s weighted distortion equals its weight identically:
     $ ( (w_(j v)^*)^p D_(j v) ) / ( sum_(r=1)^d (w_(j r)^*)^p D_(j r) ) = w_(j v)^* quad forall v in {1, dots, d}. $
  2. *Pairwise & Aggregate Noise Suppression:* For all $u in cal(N)_j$ and $v in cal(S)_j$, $w_(j u)^* / w_(j v)^* <= rho^(-p / (p - 1))$, and the total noise weight mass satisfies
     $ sum_(u in cal(N)_j) w_(j u)^* <= ( (d_n / d_s) rho^(-p/(p-1)) ) / ( 1 + (d_n / d_s) rho^(-p/(p-1)) ). $
]

#proof[
  Let $Z_j = sum_(r=1)^d D_(j r)^(-1/(p-1))$. Substituting $w_(j v)^* = D_(j v)^(-1/(p-1)) / Z_j$ directly into the $v$-th distortion summand gives
  $ (w_(j v)^*)^p D_(j v) &= Z_j^(-p) D_(j v)^(-p/(p-1) + 1) \
    &= Z_j^(-p) D_(j v)^(-1/(p-1)) = Z_j^(1-p) w_(j v)^*. $
  Summing over $v = 1, dots, d$ recovers $sum_(r=1)^d (w_(j r)^*)^p D_(j r) = Z_j^(1-p)$, and dividing the two identities yields Claim 1.

  For Claim 2, $D_(j u) / D_(j v) >= (sigma_n / sigma_s)^p = rho^p$ for every $u in cal(N)_j$, $v in cal(S)_j$. Raising to the negative power $-1/(p - 1)$ reverses the inequality:
  $ w_(j u)^* / w_(j v)^* = ( D_(j u) / D_(j v) )^(-1/(p-1)) <= rho^(-p / (p - 1)). $
  Summing $D_(j u)^(-1/(p-1)) <= |C_j|^(-1/(p-1)) sigma_n^(-p/(p-1))$ over $cal(N)_j$ and $D_(j v)^(-1/(p-1)) >= |C_j|^(-1/(p-1)) sigma_s^(-p/(p-1))$ over $cal(S)_j$, the ratio of total noise weight $W_"noise" := sum_(u in cal(N)_j) w_(j u)^*$ to signal weight $W_"sig" = 1 - W_"noise"$ satisfies
  $ W_"noise" / (1 - W_"noise") <= ( d_n sigma_n^(-p/(p-1)) ) / ( d_s sigma_s^(-p/(p-1)) ) = (d_n / d_s) rho^(-p / (p - 1)). $
  Solving for $W_"noise"$ produces the claimed rational upper bound. Because $-p/(p-1) -> -infinity$ as $p -> 1^+$, smaller exponents $p in (1, 2]$ yield dramatically sharper suppression of noisy features, whereas $p -> infinity$ weakens suppression to the inverse-linear limit $rho^(-1)$.
]

#subsec("D", "Scale Vulnerability and Non-Convex Voronoi Geometry")

#proposition(title: "Symmetry Breaking, Inverse Scale Collapse, and Non-Convex Cells", origin: "Original Proof")[
  Let $p > 1$. Then:
  1. *Symmetries & Inverse Scale Collapse:* MWK-Means is invariant under translations, coordinate permutations, sign flips, and uniform scalings $x |-> alpha x$ ($alpha > 0$), but is *not* invariant under orthogonal rotations $Q in O(d)$ nor anisotropic scalings $S(x) = "diag"(s_1, dots, s_d) x$. Specifically, shrinking a single coordinate $s_r -> 0^+$ while fixing $s_v = 1$ ($v != r$) forces $w_(j r)^* -> 1$ and $Phi_p -> 0$ for every cluster (when $epsilon_0 = 0$), even if coordinate $r$ is pure noise.
  2. *Non-Convex & Disconnected Voronoi Cells:* Even for $p = 2$, the cluster-specific Voronoi cells $V_j = {x in bb(R)^d : d_2 (x, mu_j; w_j) <= d_2 (x, mu_ell; w_ell) space forall ell}$ need neither be convex nor connected.
]

#proof[
  (1) Under uniform scaling $x |-> alpha x$ ($alpha > 0$), every dispersion scales homogeneously as $D'_(j v) = alpha^p D_(j v)$, so the common factor $alpha^(-p/(p-1))$ cancels in $w_(j v)^*$. Under anisotropic scaling with $s_r > 0$ and $s_v = 1$ ($v != r$), $D'_(j r) = s_r^p D_(j r)$, whence
  $ w'_(j r) = ( s_r^(-p/(p-1)) D_(j r)^(-1/(p-1)) ) / ( s_r^(-p/(p-1)) D_(j r)^(-1/(p-1)) + sum_(v != r) D_(j v)^(-1/(p-1)) ) -> 1 $
  as $s_r -> 0^+$, driving the reduced cluster objective $(s_r^(-p/(p-1)) D_(j r)^(-1/(p-1)) + C)^(1-p) = O(s_r^p) -> 0$. Because axis-aligned weights diagonalize the metric only in the canonical basis, rotating two axis-aligned subspaces by $45^degree$ destroys diagonal sparsity, breaking $O(d)$ invariance.

  (2) Construct a planar example ($d = 2, k = 2, p = 2$) with $mu_1 = (0, 0)$, $mu_2 = (2, 0)$ and weights $w_1 = (3/4, 1/4)$, $w_2 = (1/4, 3/4)$. For $x = (x_1, x_2) in bb(R)^2$,
  $ 16 ( d_2 (x, mu_1; w_1)^2 &- d_2 (x, mu_2; w_2)^2 ) \
    &= 9 x_1^2 + x_2^2 - (x_1 - 2)^2 - 9 x_2^2 \
    &= 4 (2 x_1^2 + x_1 - 1 - 2 x_2^2), $
  whose zero locus is a *hyperbola*. Verify two pairs of points:
  - $p_1 = (1, 2), p_2 = (1, -2)$ have $d_1^2 = 13/16 < 37/16 = d_2^2$ ($p_1, p_2 in V_1$), whereas their midpoint $(1, 0)$ has $d_1^2 = 9/16 > 1/16 = d_2^2$ (in $V_2 backslash V_1$), proving $V_1$ is *non-convex*.
  - $q_1 = (-2, 0), q_2 = (2, 0) in V_2$ ($d_2^2 <= 16/16 < 36/16 = d_1^2$), while their midpoint $(0, 0) = mu_1 in V_1$; $V_2$ splits into two disjoint hyperbolic branches.
]

#subsec("E", "Outlier Leverage Decay and Computational Complexity")

#proposition(title: "p-Dependent Outlier Leverage & Automatic Weight Damping", origin: "Original Proof")[
  Fix a cluster $C_j$ of size $m >= 3$ and send a single observation's $v$-th coordinate $x_(1 v) -> +infinity$ while holding all other data fixed. Then for every $p > 1$:
  1. The $L_p$ center $mu_(j v)^*$ diverges with asymptotic leverage rate
     $ lim_(x_(1 v) -> +infinity) mu_(j v)^* / x_(1 v) = 1 / ( 1 + (m - 1)^(1 / (p - 1)) ), $
     so the raw centroid breakdown point is $epsilon^* = 1 / n$, yet outlier leverage decays super-polynomially as $p -> 1^+$ (approaching the $L_1$ median's $50\%$ resistance) and worsens to $1/2$ as $p -> infinity$.
  2. Simultaneously, the contaminated coordinate's optimal weight vanishes as $w_(j v)^* = O(x_(1 v)^(-p / (p - 1))) -> 0$, insulating the remaining $d - 1$ coordinates from the corruption.
]

#proof[
  Let $m_0 = min_(i in C_j, i != 1) x_(i v)$ and $M_0 = max_(i in C_j, i != 1) x_(i v)$. For $x_(1 v)$ sufficiently large, the stationarity condition $g'_(j v) (mu) = 0$ of Lemma 2 forces $mu in (M_0, x_(1 v))$ and yields
  $ (x_(1 v) - mu)^(p-1) = sum_(i in C_j, i != 1) (mu - x_(i v))^(p-1). $
  Bounding each summand via $(mu - M_0)^(p-1) <= (mu - x_(i v))^(p-1) <= (mu - m_0)^(p-1)$ and taking the $(p-1)$-th root gives
  $ (m - 1)^(1/(p-1)) (mu - M_0) <= x_(1 v) - mu <= (m - 1)^(1/(p-1)) (mu - m_0). $
  Solving for $mu$ brackets $mu_(j v)^*$ within $O(1)$ of $x_(1 v) / (1 + (m - 1)^(1/(p-1)))$, proving Claim 1. Consequently $x_(1 v) - mu_(j v)^* = Theta(x_(1 v))$, so $D_(j v) >= (x_(1 v) - mu_(j v)^*)^p = Theta(x_(1 v)^p)$. Because all uncorrupted coordinates retain bounded dispersion $D_(j u) = O(1)$, Lemma 1 yields $w_(j v)^* prop (D_(j v) + epsilon_0)^(-1/(p-1)) = O(x_(1 v)^(-p/(p-1))) -> 0$.
]

#theorem(title: "NP-Hardness and Canonical Metric Specializations", origin: "Classical — Cited")[
  1. Globally minimizing $J_p (cal(C), bold(mu), W)$ (or $tilde(J)_p$) is NP-hard even for $k = 2$ @aloise2009np @deamorim2012minkowski.
  2. Fixing uniform weights $w_(j v) = 1/d$ recovers Lloyd's K-Means (Alg. 1) at $p = 2$ @lloyd1982least and K-Medians (Alg. 15) in the limit $p -> 1^+$ @bradley1997clustering; constraining weights to be globally shared ($w_1 = dots = w_k$) at $p = 2$ recovers W-K-Means (Alg. 28) at $beta = 2$ @huang2005automated.
]

#proof(label: "Proof Synthesis (Cited)")[
  By Proposition 1(1) and the reduction from standard $k$-means in @aloise2009np, embedding any 2-means instance along $d$ orthogonal copies with identical intra-cluster dispersion forces optimal uniform weights $w_(j v)^* = 1/d$ at $p = 2$, preserving NP-hardness of the joint partition problem @deamorim2012minkowski. Setting $w_(j v) = 1/d$ directly reduces $J_2$ to $d^(-2)$ times Euclidean WCSS @lloyd1982least, while letting $p -> 1^+$ drives the unique root of $g'_(j v) (c) = 0$ into the sample median interval @bradley1997clustering. Tying $w_(j v) equiv w_v$ at $p = 2$ reproduces the exact objective of @huang2005automated with weight exponent $beta = 2$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each outer iteration evaluates $n k$ weighted $L_p$ distances in $O(n k d)$ time, updates $k d$ dispersions and simplex weights in $O((n + k) d)$, and solves $k d$ 1D strictly convex center problems. Because every point belongs to one cluster, one sweep of 1D steepest descent or damped Newton across all $(j, v)$ costs $O(n d)$. With $T_"in"$ inner iterations ($T_"in" = 1$ closed-form at $p = 2$; $T_"in" = O(log(1/delta))$ via bracketed Newton/bisection for $p != 2$), total time is $O(n d (k + T_"in") T)$ and memory is strictly linear $O((n + k) d)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for MWK-Means (Alg. 29)],
    [Time / Space Complexity], [$O(n d (k + T_"in") T)$ time; $O((n + k) d)$ space; $T_"in" = 1$ at $p = 2$.],
    [Geometric Inductive Bias], [Cluster-specific axis-aligned $L_p$ subspaces with curved, potentially non-convex Voronoi cells (Prop. 1).],
    [Feature Selection Mechanism], [Automatic inverse-power weighting $w_(j v) prop tilde(D)_(j v)^(-1/(p-1))$; sharp near $p -> 1^+$ (Thm. 2).],
    [Invariances / Preprocessing], [Uniform-scale & permutation invariant; *must* range/Z-score standardize to avoid inverse scale collapse.],
    [Outlier Resistance], [Raw breakdown $1/n$, but leverage decays as $(m-1)^(-1/(p-1))$ near $p -> 1^+$ and outlier features self-downweight (Prop. 2).],
    [Recommended Companions], [Anomalous Pattern init (`iMWK-Means`) @deamorim2012minkowski; K-Means++ (Alg. 2); ORCLUS (Alg. 138) for rotated subspaces.],
  ),
  caption: [Theoretical and operational profile of Minkowski Weighted K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

MWK-Means excels when clusters inhabit distinct *axis-aligned* subspaces inside noisy tabular data and when $p in (1, 2)$ improves tail robustness.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Bioinformatics & Gene-Expression Subtyping:* Microarray and scRNA-seq cohorts embed small subtype biomarker panels inside noisy genes. By Theorem 2, MWK-Means suppresses background genes at rate $rho^(-p/(p-1))$ *per phenotype*, outputting interpretable weight vectors $w_j$ while preserving minority markers.
2. *Retail & Financial Segmentation:* Customer tables mix noisy demographics with coherent segment KPIs (spend frequency vs. tenure); local weights let cohorts form around distinct KPI subsets with auditable centroids.
3. *Industrial Sensor Fault Isolation:* Normal multi-sensor channels vary broadly while faults lock narrow sensor subsets; MWK-Means weights low-dispersion fault channels and dampens single-channel spikes (Prop. 2).
4. *Chemometrics & Spectral Fingerprinting:* Spectral profiles pair wide uninformative bands with sharp peaks; $p in [1.2, 1.8]$ blends $L_1$-like spike resistance with smooth multi-band weighting (Lemma 1).
5. *Anomaly-Initialized Discovery (`iMWK-Means`):* Paired with Anomalous Pattern seeding @deamorim2012minkowski, MWK-Means jointly recovers clusters and per-cluster feature relevance in noisy tabular analytics.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MWK-Means degrades under four regimes: (a) *unstandardized units*, where Proposition 1(1) collapses weights onto the smallest-scale raw axis; (b) *rotated/correlated subspaces*, where diagonal weights fail to align with oblique eigenvectors (use ORCLUS Alg. 138 @aggarwal2000orclus or GMM Alg. 81 @fraley2002mclust); (c) *$p -> 1^+$ over-sparsification*, discarding valid secondary features (Lemma 1); and (d) *non-convex manifolds* (use DBSCAN Alg. 31 or Spectral Alg. 103).

#figure(
  table(
    columns: (0.9fr, 1.15fr, 0.95fr),
    align: (left, left, left),
    [Domain], [Why MWK-Means Fits], [Caution / Alternative],
    [Genomics / Omics], [Per-cluster biomarker weights ($w_j$)], [Standardize; FABIA (Alg. 161)],
    [Tabular Marketing], [Noisy columns & cohort KPIs], [Mandatory scaling (Prop. 1)],
    [Sensor Fault Logs], [Outlier-channel damping ($p approx 1.5$)], [Trimmed (Alg. 26) for spikes],
    [Correlated Finance], [Weak fit on rotated subspaces], [ORCLUS (Alg. 138), GMM (Alg. 81)],
    [Spatial / LiDAR], [Distorts isotropic geometry], [Lloyd (Alg. 1), DBSCAN (Alg. 31)],
  ),
  caption: [Engineering suitability matrix and alternatives for MWK-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MWK-Means bridges centroid partitioning (Part I) and projected subspace clustering (Part VI: PROCLUS Alg. 136 @aggarwal1999proclus, COSA Alg. 145 @friedman2004cosa @kriegel2009subspacesurvey), with a strict variational advantage over global feature weighting:

#proposition(title: "Strict Subspace Separation over Global W-K-Means", origin: "Original Derivation")[
  Fix any partition $cal(C)$, centers $bold(mu)$, and $p > 1$ with $D_(j v) > 0$. Let $Phi_p^("MWK")$, $Phi_p^("WKM")$, and $Phi_p^("unif")$ be the optimal $J_p$ values under cluster-specific $W in (Delta_(d-1))^k$, globally tied $w_1 = dots = w_k in Delta_(d-1)$ (W-K-Means, Alg. 28), and uniform $w_(j v) = 1/d$. Then
  $ Phi_p^("MWK") <= Phi_p^("WKM") <= Phi_p^("unif") = d^(-p) sum_(j=1)^k sum_(v=1)^d D_(j v), $
  with $Phi_p^("MWK") = Phi_p^("WKM")$ iff $(D_(j 1), dots, D_(j d))$ are proportional across all $j$. When $k = d$ clusters occupy disjoint 1D subspaces (signal $1$, noise $R >> 1$), $Phi_p^("WKM") / Phi_p^("MWK") = Theta(R) -> infinity$.
]

#proof[
  Nested feasible sets ${(1/d) bold(1)} subset.eq Delta_"glob" subset.eq (Delta_(d-1))^k$ yield the chain. By Lemma 1, the global optimum $w^("glob")$ minimizes $sum_v w_v^p (sum_j D_(j v))$ and matches every local $w_j^*$ iff $D_(j v) / D_(j u) = D_(ell v) / D_(ell u)$ for all $j, ell, u, v$, i.e., proportional rows of $[D_(j v)]$. For $k = d$ disjoint subspaces ($D_(j j) = 1$, $D_(j v) = R$ for $v != j$), equal column sums $1 + (d - 1) R$ force $w_v^("glob") = 1/d$ and $Phi_p^("WKM") = d^(1-p)(1 + (d-1) R) = Theta(R)$, whereas Lemma 1 gives $Phi_p^("MWK") = d(1 + (d-1) R^(-1/(p-1)))^(1-p) <= d$, proving $Phi_p^("WKM") / Phi_p^("MWK") = Omega(R)$.
]

Empirically, on the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, MWK-Means achieves an outstanding *Iris ARI of 0.859* (unweighted Euclidean Silhouette *0.393*), far outperforming isotropic LBG (Alg. 30, Iris ARI 0.583). Confirming Theorem 2, MWK-Means concentrates weight on discriminative petal axes while down-weighting overlapping sepal variance, recovering species labels even though unweighted Euclidean Silhouette penalizes re-weighted geometry. On *Multi-Shape*, it reaches *Silhouette 0.718* (vs. 0.609 for LBG), exploiting the curved Voronoi boundaries of Proposition 1(2). On *Customers* and *High-D Subspace+Noise*, Euclidean Silhouettes are *0.229* and *0.229* (ARI: N/A), showing that ambient Euclidean Silhouette remains diluted by uninformative axes. Its mean runtime of *0.874 ms* confirms $O(n d (k + T_"in") T)$ scalability—roughly $5 times$ slower than LBG (0.171 ms) due to iterative Minkowski centers and per-cluster weight updates.
