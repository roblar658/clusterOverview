#import "../template/ieee.typ": *

#alg_entry(
  id: 8,
  name: "K-Prototypes",
  category: "Cat. I: Centroid & Partitioning",
  year: "1998",
  paradigm: "Mixed-Domain Alternating Prototype Minimization",
  complexity: $O(n k (d_r + d_c) T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Real-world enterprise, clinical, and telemetry datasets inhabit mixed product spaces $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^(d_r) times product_(m=1)^(d_c) cal(A)_m$, where each record $x_i = (x_i^((r)), x_i^((c)))$ pairs a numeric vector $x_i^((r)) in bb(R)^(d_r)$ with a categorical tuple $x_i^((c)) = (x_(i, 1)^((c)), dots, x_(i, d_c)^((c)))$ over finite unordered alphabets $cal(A)_m$ of size $v_m = |cal(A)_m| >= 2$. Standard K-Means (Alg. 1) cannot handle $cal(A)_m$ without artificial numeric coding, whereas K-Modes (Alg. 7) requires discretizing $bb(R)^(d_r)$ with quantization loss.

Huang @huang1998extensions unifies both paradigms by representing each cluster $C_j$ ($j = 1, dots, k$) via a *mixed prototype* $q_j = (q_j^((r)), q_j^((c))) in bb(R)^(d_r) times product_(m=1)^(d_c) cal(A)_m$ equipped with the mixed dissimilarity
$ d_gamma (x_i, q_j) = norm(x_i^((r)) - q_j^((r)))_2^2 + gamma H(x_i^((c)), q_j^((c))), $
where $gamma > 0$ balances numeric versus categorical distortion and $H(u, v) = sum_(m=1)^(d_c) delta(u_m, v_m)$ is the simple-matching Hamming distance ($delta(a, b) = bb(I)[a != b]$). The *K-Prototypes objective* seeks a disjoint partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ and prototypes $Q = (q_1, dots, q_k)$ minimizing
$ J_gamma (cal(C), Q) &= sum_(j=1)^k sum_(i in C_j) d_gamma (x_i, q_j) \
  &= J_"num" (cal(C), Q^((r))) + gamma J_"cat" (cal(C), Q^((c))), $
where $J_"num"$ is Euclidean WCSS and $J_"cat"$ is total categorical mismatch count. For non-empty $S subset.eq {1, dots, n}$ and $a in cal(A)_m$, let $f_m (a; S) = |{i in S : x_(i, m)^((c)) = a}|$ and $p_m (a; S) = f_m (a; S) / |S|$.

#definition(title: "Mixed-Domain Prototype Equilibrium")[
  A pair $(cal(C), Q)$ is a *mixed-domain prototype equilibrium* (coordinate-wise stationary point of $J_gamma$) if: (i) every $i in C_j$ satisfies $d_gamma (x_i, q_j) <= d_gamma (x_i, q_ell)$ for all $ell != j$; (ii) every non-empty cluster has numeric prototype equal to the arithmetic mean $q_j^((r)) = |C_j|^(-1) sum_(i in C_j) x_i^((r))$; and (iii) every categorical coordinate $q_(j, m)^((c))$ is a modal category maximizing $f_m (a; C_j)$ over $a in cal(A)_m$.
]

K-Prototypes minimizes $J_gamma$ by alternating hard Voronoi-like mixed assignments with decoupled mean-and-mode prototype projections (Algorithm 1), generalizing Lloyd's algorithm @lloyd1982least and K-Modes @huang1998extensions simultaneously.

#algbox(
  title: "Huang's K-Prototypes Algorithm",
  inputs: [Mixed dataset $cal(X) = {(x_i^((r)), x_i^((c)))}_(i=1)^n$, clusters $k >= 1$, balance weight $gamma > 0$, initial prototypes $Q^((0))$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and prototypes $Q = (q_1, dots, q_k)$.],
)[
  1. Fix a canonical total order on each alphabet $cal(A)_m$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) *Assignment step:* For each $i in {1, dots, n}$, compute
     $ c_i^((t)) = arg min_(1 <= j <= k) [ &norm(x_i^((r)) - q_j^((r, t)))_2^2 \
       &+ gamma H(x_i^((c)), q_j^((c, t))) ], $
     breaking ties by smallest cluster index $j$.
  4. #h(1em) Form $C_j^((t)) = {i : c_i^((t)) = j}$ for $j = 1, dots, k$ (reseed any empty cluster to the point maximizing $min_ell d_gamma (x_i, q_ell^((t)))$).
  5. #h(1em) *Prototype update:* For each non-empty cluster $C_j^((t))$:
     - Numeric mean: $q_j^((r, t+1)) = |C_j^((t))|^(-1) sum_(i in C_j^((t))) x_i^((r))$.
     - Categorical mode: for each $m = 1, dots, d_c$, set
       $ q_(j, m)^((c, t+1)) = arg max_(a in cal(A)_m) f_m (a; C_j^((t))), $
       retaining $q_(j, m)^((c, t))$ if tied for the maximum, otherwise choosing the canonically smallest maximizer.
  6. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$ and $Q^((t+1)) = Q^((t))$, *break*.
  7. *Return* $(cal(C)^((t)), Q^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of K-Prototypes: (1) exact separable prototype optimality and strict Lyapunov descent with finite-step termination; (2) an exact simplex-vertex isometric embedding revealing a universal factor-of-2 relation to scaled one-hot K-Means; (3) exact phase-transition thresholds for the balance parameter $gamma$; and (4) product-group symmetry invariances, asymmetric numeric-versus-categorical breakdown robustness, and NP-hardness.

#subsec("A", "Separable Prototype Optimality and Finite Convergence")

#lemma(title: "Exact Separable Prototype Excess Identity", origin: "Original Proof")[
  For any non-empty cluster $S subset.eq {1, dots, n}$, weight $gamma > 0$, and candidate prototype $q = (q^((r)), q^((c)))$, let $hat(q)_S = (overline(x)_S^((r)), hat(q)_S^((c)))$ be any prototype whose numeric block is the mean $overline(x)_S^((r))$ and whose categorical coordinates satisfy $hat(q)_(S, m)^((c)) in arg max_(a in cal(A)_m) f_m (a; S)$. Then $h_S (q) := sum_(i in S) d_gamma (x_i, q)$ satisfies
  $ h_S (q) - h_S (hat(q)_S) &= |S| norm(q^((r)) - overline(x)_S^((r)))_2^2 \
    &quad + gamma sum_(m=1)^(d_c) [ f_m (hat(q)_(S, m)^((c)); S) - f_m (q_m^((c)); S) ] >= 0. $
  Equality holds if and only if $q^((r)) = overline(x)_S^((r))$ and every $q_m^((c))$ attains the maximum frequency in $S$.
]

#proof[
  Decompose $h_S (q) = sum_(i in S) norm(x_i^((r)) - q^((r)))_2^2 + gamma sum_(m=1)^(d_c) sum_(i in S) delta(x_(i, m)^((c)), q_m^((c)))$. For the numeric sum, inserting $plus.minus overline(x)_S^((r))$ yields the parallel-axis identity
  $ sum_(i in S) norm(x_i^((r)) - q^((r)))_2^2 &= sum_(i in S) norm(x_i^((r)) - overline(x)_S^((r)))_2^2 \
    &quad + |S| norm(q^((r)) - overline(x)_S^((r)))_2^2, $
  because $chevron.l sum_(i in S) (x_i^((r)) - overline(x)_S^((r))), overline(x)_S^((r)) - q^((r)) chevron.r = 0$. For each categorical attribute $m in {1, dots, d_c}$ and symbol $a in cal(A)_m$, exact counting gives $sum_(i in S) delta(x_(i, m)^((c)), a) = sum_(i in S) (1 - bb(I)[x_(i, m)^((c)) = a]) = |S| - f_m (a; S)$. Subtracting the value at $hat(q)_(S, m)^((c))$ cancels $|S|$ identically:
  $ &sum_(i in S) [ delta(x_(i, m)^((c)), q_m^((c))) - delta(x_(i, m)^((c)), hat(q)_(S, m)^((c))) ] \
    &quad = f_m (hat(q)_(S, m)^((c)); S) - f_m (q_m^((c)); S). $
  Multiplying by $gamma > 0$ and summing over $m = 1, dots, d_c$ gives the exact excess formula. Because $|S| >= 1$, $gamma > 0$, and $f_m (hat(q)_(S, m)^((c)); S) >= f_m (q_m^((c)); S)$ for every $m$, every summand is non-negative and vanishes if and only if $q^((r)) = overline(x)_S^((r))$ and each $q_m^((c))$ is modal.
]

#theorem(title: "Monotone Lyapunov Descent and Finite-Step Convergence", origin: "Classical — Cited")[
  Let ${(cal(C)^((t)), Q^((t)))}$ be generated by Algorithm 1 with deterministic tie-breaking and non-empty clusters @huang1998extensions. Then:
  1. Every prototype update decreases $J_gamma$ by the exact non-negative margin
     $ Delta_"proto"^((t)) &= sum_(j=1)^k |C_j^((t))| norm(q_j^((r, t+1)) - q_j^((r, t)))_2^2 \
       &quad + gamma sum_(j=1)^k sum_(m=1)^(d_c) [ f_m (q_(j, m)^((c, t+1)); C_j^((t))) \
       &quad - f_m (q_(j, m)^((c, t)); C_j^((t))) ]. $
  2. Every assignment step satisfies $J_gamma (cal(C)^((t+1)), Q^((t+1))) <= J_gamma (cal(C)^((t)), Q^((t+1)))$.
  3. Algorithm 1 terminates in at most $T <= S(n, k) <= k^n / k!$ iterations at a mixed-domain prototype equilibrium.
]

#proof[
  Applying Lemma 1 independently to each cluster $S = C_j^((t))$ with $q = q_j^((t))$ and $hat(q)_S = q_j^((t+1))$ and summing over $j = 1, dots, k$ yields $J_gamma (cal(C)^((t)), Q^((t))) - J_gamma (cal(C)^((t)), Q^((t+1))) = Delta_"proto"^((t)) >= 0$, proving Claim 1. In the subsequent assignment step, choosing $c_i^((t+1)) = arg min_ell d_gamma (x_i, q_ell^((t+1)))$ guarantees $d_gamma (x_i, q_(c_i^((t+1)))^((t+1))) <= d_gamma (x_i, q_(c_i^((t)))^((t+1)))$ for every $i$, and summing over $i$ proves Claim 2.

  To exclude limit cycles when modes are non-unique, consider the tie-retention rule in Step 5 of Algorithm 1. Suppose $J_gamma (cal(C)^((t+1)), Q^((t+2))) = J_gamma (cal(C)^((t)), Q^((t+1)))$. Then both $Delta_"proto"^((t+1)) = 0$ and the assignment drop is zero. By Lemma 1, $Delta_"proto"^((t+1)) = 0$ forces $q_j^((r, t+2)) = q_j^((r, t+1))$ and $f_m (q_(j, m)^((c, t+2)); C_j^((t+1))) = f_m (q_(j, m)^((c, t+1)); C_j^((t+1)))$; hence the retention rule preserves $q_(j, m)^((c, t+2)) = q_(j, m)^((c, t+1))$, so $Q^((t+2)) = Q^((t+1))$. Because the assignment step achieved zero drop, every point had minimal dissimilarity already at its previous cluster, so deterministic smallest-index tie-breaking enforces $cal(C)^((t+1)) = cal(C)^((t))$ once prototypes stabilize, triggering immediate termination. Consequently, every non-terminal partition attains a strictly lower optimal objective $J_gamma^* (cal(C)) := min_Q J_gamma (cal(C), Q)$ and can never repeat. Since at most $sum_(m=1)^k S(n, m) <= k^n / k!$ partitions exist, Algorithm 1 halts finitely at a mixed-domain prototype equilibrium.
]

#subsec("B", "Simplex-Vertex Isometry and Factor-of-2 One-Hot Equivalence")

Practitioners frequently debate whether to cluster mixed data via K-Prototypes or via one-hot encoding followed by Euclidean K-Means. The next result establishes the exact geometric bridge between both methods.

#theorem(title: "One-Hot Euclidean Isometry and Universal Factor-of-2 Bound", origin: "Original Derivation")[
  Let $V = sum_(m=1)^(d_c) v_m$ and define the scaled one-hot embedding $phi_gamma : bb(R)^(d_r) times product_(m=1)^(d_c) cal(A)_m -> bb(R)^(d_r + V)$ by
  $ phi_gamma (x) = (x^((r)), sqrt(gamma / 2) e_(x_1^((c))), dots, sqrt(gamma / 2) e_(x_(d_c)^((c)))), $
  where $e_a in {0, 1}^(v_m)$ is the standard basis vector for $a in cal(A)_m$. Then:
  1. For every data point $x_i$ and discrete prototype $q_j$, $norm(phi_gamma (x_i) - phi_gamma (q_j))_2^2 = d_gamma (x_i, q_j)$ identically.
  2. For any non-empty cluster $S$, let $J_"cat"^"mode" (S)$ be its minimum K-Prototypes categorical cost and $J_"cat"^"mean" (S)$ be the Euclidean WCSS of its one-hot categorical subvectors around their continuous simplex mean $mu_S^((c))$. Then
     $ &J_"cat"^"mode" (S) - J_"cat"^"mean" (S) \
       &quad = (gamma |S|) / 2 sum_(m=1)^(d_c) [ 1 - 2 p_(m, max) (S) + sum_(a in cal(A)_m) p_m (a; S)^2 ] >= 0, $
     where $p_(m, max) (S) = max_(a in cal(A)_m) p_m (a; S)$.
  3. For every partition $cal(C)$, the relaxed one-hot WCSS $J_"OH"^* (cal(C))$ and discrete K-Prototypes cost $J_gamma^* (cal(C))$ obey the universal sandwich bound
     $ 1 / 2 J_gamma^* (cal(C)) &<= J_"num"^* (cal(C)) + 1 / 2 gamma J_"cat"^"mode" (cal(C)) \
       &<= J_"OH"^* (cal(C)) <= J_gamma^* (cal(C)). $
]

#proof[
  (1) Because $norm(e_a - e_b)_2^2 = 2 delta(a, b)$ for basis vectors $e_a, e_b in {0, 1}^(v_m)$, multiplying by $gamma / 2$ and summing over $m = 1, dots, d_c$ with $norm(x_i^((r)) - q_j^((r)))_2^2$ gives $norm(phi_gamma (x_i) - phi_gamma (q_j))_2^2 = norm(x_i^((r)) - q_j^((r)))_2^2 + sum_(m=1)^(d_c) (gamma / 2) dot 2 delta(x_(i, m)^((c)), q_(j, m)^((c))) = d_gamma (x_i, q_j)$.

  (2) Fix attribute $m$ and cluster $S$ of size $N = |S|$. Let $p_m = (p_m (a; S))_(a in cal(A)_m) in Delta_(v_m - 1)$ be the empirical probability vector, and let $a^* in arg max_a p_m (a; S)$ so $p_m (a^*; S) = p_(m, max) (S)$. The continuous mean of ${sqrt(gamma / 2) e_(x_(i, m)^((c))) : i in S}$ is $sqrt(gamma / 2) p_m$, whereas the discrete modal prototype embeds to $sqrt(gamma / 2) e_(a^*)$. By the parallel-axis identity on $bb(R)^(v_m)$,
  $ &sum_(i in S) gamma / 2 norm(e_(x_(i, m)^((c))) - e_(a^*))_2^2 - sum_(i in S) gamma / 2 norm(e_(x_(i, m)^((c))) - p_m)_2^2 \
    &quad = (gamma N) / 2 norm(e_(a^*) - p_m)_2^2. $
  Expanding $norm(e_(a^*) - p_m)_2^2 = (1 - p_(m, max) (S))^2 + sum_(a != a^*) p_m (a; S)^2 = 1 - 2 p_(m, max) (S) + sum_(a in cal(A)_m) p_m (a; S)^2$ and summing over $m = 1, dots, d_c$ yields the exact excess identity.

  (3) Directly evaluating the continuous one-hot inertia yields the scaled Gini impurity $J_"cat"^"mean" (S) = (gamma N / 2) sum_(m=1)^(d_c) (1 - sum_(a in cal(A)_m) p_m (a; S)^2)$. Because $p_m (a; S) <= p_(m, max) (S)$ for all $a in cal(A)_m$ and $sum_a p_m (a; S) = 1$, we have $sum_a p_m (a; S)^2 <= p_(m, max) (S)$, hence $1 - sum_a p_m (a; S)^2 >= 1 - p_(m, max) (S)$. Multiplying by $gamma N / 2$ yields
  $ 1 / 2 gamma J_"cat"^"mode" (S) &= (gamma N) / 2 sum_(m=1)^(d_c) (1 - p_(m, max) (S)) \
    &<= J_"cat"^"mean" (S) <= gamma J_"cat"^"mode" (S). $
  Adding $J_"num"^* (S)$ and summing over $S in cal(C)$ proves Claim 3, with equality on the lower categorical bound attained precisely by uniform distributions $p_m (a; S) = 1 / v_m$.
]

#remark(title: "Interpretability vs. Relaxation")[
  Theorem 2 shows K-Prototypes is vertex-constrained Euclidean K-Means on $bb(R)^(d_r + V)$: relaxing categorical vertices $e_a$ to simplex vectors $p_m in Delta_(v_m - 1)$ replaces mode matching with Gini minimization, lowering categorical distortion by at most $2 times$ while sacrificing single-symbol prototypes and inflating memory from $O(d_c)$ to $O(V)$.
]

#subsec("C", "Hyperparameter Phase Transitions and Scale Calibration")

The practical performance of K-Prototypes hinges critically on $gamma$. We derive exact dominance thresholds beyond which one data modality collapses into a mere tie-breaker.

#theorem(title: "Dominance Phase Transitions and Balanced Calibration of gamma", origin: "Original Proof")[
  Let $D_r = max_(i, ell) norm(x_i^((r)) - x_ell^((r)))_2 > 0$ be the numeric diameter of $cal(X)$, and assume numeric prototypes lie in $"conv"({x_1^((r)), dots, x_n^((r))})$ (satisfied at every iteration after data seeding).
  1. *Categorical-Dominance Phase:* If $gamma > D_r^2$, then whenever $H(x_i^((c)), q_j^((c))) < H(x_i^((c)), q_ell^((c)))$, one strictly has $d_gamma (x_i, q_j) < d_gamma (x_i, q_ell)$ irrespective of numeric values; numeric features act solely as tie-breakers on equal-Hamming prototypes.
  2. *Numeric-Dominance Phase:* For a given prototype set $Q^((r))$, let $Delta_min (Q^((r))) > 0$ be the smallest positive numeric squared-distance margin $norm(x_i^((r)) - q_ell^((r)))_2^2 - norm(x_i^((r)) - q_j^((r)))_2^2$. Whenever $0 < gamma < Delta_min (Q^((r))) / d_c$, categorical attributes cannot overturn any strict numeric Voronoi assignment.
  3. *Expected Pairwise Balance Calibration:* For independent draws $X, Y$ with numeric covariance $Sigma_r$ and independent categorical marginals $p_m$, equalizing expected pairwise dissimilarity contributions yields
     $ gamma_"bal" = (2 "Tr"(Sigma_r)) / (sum_(m=1)^(d_c) (1 - sum_(a in cal(A)_m) p_m (a)^2)). $
]

#proof[
  (1) Because $u |-> norm(x_i^((r)) - u)_2$ is convex on $bb(R)^(d_r)$, its maximum over $"conv"({x_1^((r)), dots, x_n^((r))})$ is attained at an extreme point, giving $0 <= norm(x_i^((r)) - q_j^((r)))_2^2 <= D_r^2$ for all $i, j$. Since Hamming distances are integer-valued, $H(x_i^((c)), q_j^((c))) < H(x_i^((c)), q_ell^((c)))$ implies $H(x_i^((c)), q_ell^((c))) - H(x_i^((c)), q_j^((c))) >= 1$, whence $d_gamma (x_i, q_ell) - d_gamma (x_i, q_j) >= 0 - D_r^2 + gamma dot 1 = gamma - D_r^2 > 0$.

  (2) Conversely, since $0 <= H(dot, dot) <= d_c$, the categorical difference satisfies $|H(x_i^((c)), q_j^((c))) - H(x_i^((c)), q_ell^((c)))| <= d_c$. For any pair $(j, ell)$ with strict numeric margin $norm(x_i^((r)) - q_ell^((r)))_2^2 - norm(x_i^((r)) - q_j^((r)))_2^2 >= Delta_min (Q^((r))) > gamma d_c$, we obtain $d_gamma (x_i, q_ell) - d_gamma (x_i, q_j) >= Delta_min (Q^((r))) - gamma d_c > 0$, preserving every strict numeric Voronoi preference identically.

  (3) For independent copies $X, Y$, linearity of expectation gives $bb(E)norm(X^((r)) - Y^((r)))_2^2 = 2 "Tr"(Sigma_r)$, while for each categorical attribute $m$, $bb(P)(X_m^((c)) != Y_m^((c))) = 1 - sum_(a in cal(A)_m) p_m (a)^2$. Setting $bb(E)norm(X^((r)) - Y^((r)))_2^2 = gamma bb(E)[H(X^((c)), Y^((c)))]$ yields $gamma_"bal"$. In particular, when numeric features are Z-standardized ($"Tr"(Sigma_r) = d_r$) and categorical alphabets are uniform of size $v$, $gamma_"bal" = 2 d_r / (d_c (1 - 1 / v))$, refining Huang's heuristic $gamma approx 0.5 overline(sigma)$ @huang1998extensions to account explicitly for modality dimensionality and alphabet entropy.
]

#subsec("D", "Product-Group Invariance, Breakdown, and Hardness")

#proposition(title: "Symmetry Group, Asymmetric Breakdown, and NP-Hardness", origin: "Original Proof")[
  For any mixed dataset $cal(X)$ and $gamma > 0$:
  1. $J_gamma$ is strictly invariant under the product group $G = E(d_r) times product_(m=1)^(d_c) S_(v_m)$ of rigid Euclidean isometries $x^((r)) |-> Q x^((r)) + b$ ($Q in O(d_r)$) paired with arbitrary symbol permutations $pi_m in S_(v_m)$. Uniform numeric scaling $x^((r)) |-> s x^((r))$ ($s > 0$) is exactly equivalent to rescaling $gamma |-> gamma / s^2$.
  2. The numeric prototype block has finite-sample breakdown point $epsilon_r^* = 1 / n$, whereas within any cluster $C_j$ whose top two attribute-$m$ counts satisfy $f_((1)) > f_((2))$, the modal coordinate $q_(j, m)^((c))$ is strictly preserved under up to $r_max = floor((f_((1)) - f_((2)) - 1) / 2)$ arbitrary categorical corruptions within $C_j$.
  3. Global minimization of $J_gamma$ is NP-hard even for $k = 2$.
]

#proof[
  (1) Orthogonal transforms preserve $norm(Q(x^((r)) - q^((r))))_2^2 = norm(x^((r)) - q^((r)))_2^2$ and commute with arithmetic means; bijections $pi_m in S_(v_m)$ satisfy $delta(pi_m (a), pi_m (b)) = delta(a, b)$ and preserve frequency rankings. Under $x^((r)) |-> s x^((r))$, $J_gamma$ maps to $s^2 (J_"num" + (gamma / s^2) J_"cat")$, yielding identical minimizers to $J_(gamma / s^2)$.

  (2) Sending a single point's numeric vector $norm(x_1^((r)))_2 -> infinity$ drives its cluster mean $norm(q_j^((r)))_2 >= |C_j|^(-1) norm(x_1^((r)))_2 - O(1) -> infinity$, proving $epsilon_r^* = 1 / n$. Conversely, arbitrarily altering $r$ categorical entries of $C_j$ on attribute $m$ decreases $f_((1))$ by at most $r$ and increases any rival count $f(b; C_j) <= f_((2))$ by at most $r$. Whenever $2 r < f_((1)) - f_((2))$, i.e., $r <= floor((f_((1)) - f_((2)) - 1) / 2)$, the perturbed counts obey $f'_((1)) >= f_((1)) - r > f_((2)) + r >= f'(b)$, preserving the unique mode strictly.

  (3) Embedding any Euclidean 2-means instance in $bb(R)^(d_r)$ with constant categorical coordinates $x_(i, m)^((c)) equiv a_0$ makes $J_"cat" equiv 0$, reducing global K-Prototypes directly to Euclidean 2-means, which is NP-hard @aloise2009np.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration evaluates $n k$ numeric squared distances in $O(n k d_r)$ flops and $n k$ Hamming distances in $O(n k d_c)$ integer comparisons. Updating prototypes requires $O(n d_r)$ additions for numeric means and $O(n d_c + k V)$ operations for categorical mode frequency tables (where $V = sum_m v_m$). Thus total time is $O(n k (d_r + d_c) T + k V T)$—strictly linear in $n$ and avoiding the $O(n^2)$ pairwise Gower-distance matrix required by PAM (Alg. 4 @kaufman1990finding). Memory consumption is $O(n(d_r + d_c) + k(d_r + V))$. Geometrically, each cluster region is a union of numeric convex Voronoi polyhedra stratified across discrete categorical Hamming shells, avoiding the curse of high-cardinality one-hot matrix blow-up while preserving exact single-symbol interpretability.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Prototypes],
    [Time / Space Complexity], [$O(n k (d_r + d_c) T)$ time; $O(n(d_r + d_c) + k(d_r + V))$ space; linear in $n$.],
    [Geometric Inductive Bias], [Spherical numeric Voronoi cells coupled with Hamming-ball modal shells.],
    [Hyperparameter Sensitivity], [Sharp $gamma$ phase transitions (Theorem 3); requires Z-scoring + $gamma_"bal"$.],
    [Invariances], [Invariant under $E(d_r) times product_m S_(v_m)$; sensitive to anisotropic numeric scale.],
    [Breakdown Robustness], [Asymmetric: $epsilon_r^* = 1 / n$ numeric vs. margin-tolerant categorical modes.],
    [Recommended Companions], [W-K-Means @huang2005automated for feature weights; PAM (Alg. 4) w/ Gower for numeric outliers; AutoClass @cheeseman1996autoclass for probabilistic mixed models.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of K-Prototypes.],
)

#sec("IV", "Engineering and Scientific Applications")

K-Prototypes is the standard linear-time partitioning engine for tabular records combining continuous measurements with high-cardinality discrete identifiers.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Electronic Health Records (EHR) & Clinical Phenotyping:* Patient cohorts blend continuous biomarkers (eGFR, HbA1c, blood pressure) with nominal ICD-10 codes, medication classes, and blood types; $G$-invariance (Proposition 1) preserves arbitrary coding while yielding auditable mean-and-modal clinical phenotypes.
2. *Retail & Banking Customer Segmentation:* CRM tables pair continuous spend, recency, and credit utilization with channel, region, and tier categories. Theorem 2 guarantees factor-of-2 fidelity to one-hot K-Means while keeping prototypes compact and directly actionable.
3. *Insurance Underwriting & Actuarial Risk Pooling:* Policy records combine continuous age, insured value, and mileage with vehicle class and coverage codes; calibrating $gamma_"bal"$ (Theorem 3) prevents high-variance monetary scales from drowning categorical risk tiers.
4. *Network Intrusion & Security Telemetry:* Flow logs pair packet rates and byte counts with protocol flags, ports, and TCP states; $O(n k (d_r + d_c))$ scaling handles millions of flows where $O(n^2)$ Gower-medoid matrices are infeasible.
5. *Manufacturing & Materials Informatics:* Process batches combine continuous temperatures and cycle times with alloy grades, shifts, and supplier lots, enabling modal prototypes to pinpoint discrete root-cause configurations behind continuous yield shifts.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-Prototypes fails systematically under four regimes: (a) *uncalibrated $gamma$ or unscaled numeric features*, triggering the phase collapse of Theorem 3(1)–(2); (b) *ordinal categorical attributes* (e.g., education levels or Likert scales), where the flat $0/1$ indicator $delta(a, b)$ ignores rank ordering; (c) *heavy-tailed numeric outliers*, which corrupt $q_j^((r))$ at breakdown point $1 / n$ (Proposition 1); and (d) *high-dimensional irrelevant noise features*, where unweighted Euclidean and Hamming sums dilute true cluster signal (mitigated by automated variable weighting @huang2005automated @deamorim2012minkowski).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Prototypes Fits / Prerequisite], [Caution / Alternative],
    [EHR Phenotyping], [Direct mean + modal clinical profiles], [Winsorize extreme lab outliers],
    [CRM / Banking], [Linear in $n$; avoids one-hot blow-up], [Calibrate $gamma$ via $gamma_"bal"$],
    [Insurance Risk], [Handles nominal tiers + continuous exposure], [Use Gower+PAM (Alg. 4) if small $n$],
    [NetFlow Telemetry], [Fast mixed-type anomaly grouping], [Use W-K-Means @huang2005automated for noise],
    [Purely Continuous], [Reduces identically to K-Means (Alg. 1)], [Use GMM (Alg. 81) for anisotropy],
  ),
  caption: [Application suitability and engineering guidance for K-Prototypes.],
)

#sec("V", "Structural Relations and Empirical Behavior")

K-Prototypes sits at the exact algebraic intersection of continuous centroid partitioning and discrete modal partitioning, admitting a clean limiting unification across Category I and probabilistic latent-variable families.

#proposition(title: "Limiting Reductions and Probabilistic Surrogate Relation", origin: "Original Unified Proof")[
  1. Setting $d_c = 0$ (or $gamma -> 0^+$ with distinct numeric margins) reduces Algorithm 1 identically to Lloyd's K-Means (Alg. 1 @lloyd1982least); setting $d_r = 0$ reduces Algorithm 1 identically to Huang's K-Modes (Alg. 7 @huang1998extensions).
  2. Consider a mixed generative model whose $j$-th component is a product of an isotropic Gaussian $cal(N)(q_j^((r)), sigma^2 I_(d_r))$ and $d_c$ independent categorical noise distributions $P_m (x_m | q_(j, m)^((c))) prop exp(-eta delta(x_m, q_(j, m)^((c))))$. Maximizing the joint hard-assignment log-likelihood over $(cal(C), Q)$ is strictly equivalent to minimizing $J_gamma (cal(C), Q)$ with $gamma = 2 sigma^2 eta$.
]

#proof[
  (1) When $d_c = 0$, $H equiv 0$ and the categorical update is vacuous, recovering Lloyd's K-Means verbatim; when $d_r = 0$, $norm(x_i^((r)) - q_j^((r)))_2^2 equiv 0$ and scaling by $gamma > 0$ preserves argmin, recovering K-Modes verbatim. By Theorem 3(2), $gamma -> 0^+$ yields identical strict numeric Voronoi assignments.

  (2) Taking negative log-likelihood of the conditional product density for $x_i$ assigned to component $j$ yields
  $ -log p(x_i | q_j) &= 1 / (2 sigma^2) norm(x_i^((r)) - q_j^((r)))_2^2 + eta H(x_i^((c)), q_j^((c))) + C_0 \
    &= 1 / (2 sigma^2) [ norm(x_i^((r)) - q_j^((r)))_2^2 + (2 sigma^2 eta) H(x_i^((c)), q_j^((c))) ] + C_0, $
  where $C_0$ is independent of $(cal(C), Q)$ because the normalization constant $Z_m (eta) = 1 + (v_m - 1) e^(-eta)$ does not depend on the modal symbol $q_(j, m)^((c))$. Summing over $i = 1, dots, n$ and multiplying by $2 sigma^2 > 0$ matches $J_gamma (cal(C), Q)$ identically for $gamma = 2 sigma^2 eta$.
]

This generative link connects K-Prototypes directly to Latent Class Analysis (Alg. 90 @lazarsfeld1968latent), AutoClass @cheeseman1996autoclass, and Bregman clustering @banerjee2005clustering, while automated feature-weighted extensions @huang2005automated replace scalar $gamma$ with learned coordinate weights $w_r, w_m$ to suppress irrelevant noise attributes.

On the ClusterOverview 4-dataset empirical benchmark, K-Prototypes executes with exceptional speed (mean runtime $0.154$ ms), achieving strong separation on Multi-Shape ($"Sil" = 0.712$, $"ARI" = "N/A"$) and moderate separation on Iris ($"Sil" = 0.474$, $"ARI" = 0.430$), while attaining $"Sil" = 0.342$ on Customers and $"Sil" = 0.303$ on High-D Subspace+Noise ($"ARI" = "N/A"$). Because these benchmark suites emphasize continuous and mixed-overlap structure with non-spherical or noisy subspaces, the empirical profile mirrors our theoretical bounds: fast linear-time Lyapunov convergence (Theorem 1), clean convex separation on compact blobs, and sensitivity to unweighted irrelevant dimensions and non-isotropic overlap (Theorem 3 and Proposition 1).
