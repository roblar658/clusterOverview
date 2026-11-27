#import "../template/ieee.typ": *

#alg_entry(
  id: 52,
  name: "SNN-DPC: Shared-Nearest-Neighbor Density Peaks Clustering",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2018",
  paradigm: "Shared-Nearest-Neighbor Metric & Two-Step Peak Assignment",
  complexity: $cal(O)(n^2 + n k^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Although Density Peaks Clustering (DPC, Alg. 51, @rodriguez2014dpc) recovers non-convex clusters in a single non-iterative sweep, Section II of Alg. 51 revealed three fundamental theoretical vulnerabilities: (i) sensitivity to the global cutoff parameter $d_c$, which suppresses sparse cluster peaks in multi-density mixtures; (ii) unscaled Euclidean higher-density distances $delta_i$, which allow dense-cluster shoulders to outrank well-separated sparse peaks under $gamma_i = rho_i delta_i$; and (iii) single-hop parent-pointer assignment, where one erroneous saddle edge triggers a domino misclassification cascade across an entire subtree.

Shared-Nearest-Neighbor Density Peaks Clustering (SNN-DPC) @liu2018snndpc resolves all three pathologies simultaneously by marrying DPC's decision-graph paradigm with the topological robustness of Shared Nearest Neighbor (SNN) similarity @ertoz2003snn. SNN-DPC replaces the scale-dependent cutoff $d_c$ with a single topological neighborhood size $k in bb(N)$, defines local density $rho_i$ from reciprocal $k$-NN overlap weighted by inverse two-hop path length, scales $delta_i$ by local neighborhood radii to equalize multi-density peaks, and replaces single-hop tree propagation with a two-step assignment architecture that separates high-confidence core links (*inevitably subordinate* points) from consensus-voted peripheral points (*possibly subordinate* points).

#definition(title: "SNN Similarity, Local Density, and Scale-Compensated Separation")[
  Let $cal(X) = {x_1, dots, x_n}$ be equipped with metric $d_(i j) = d(x_i, x_j)$ and integer $k >= 2$. Denote by $Gamma(i) = {i} union "kNN"_k (i)$ the closed $k$-nearest-neighbor set of $i$ (cardinality $k + 1$), and define the *shared-nearest-neighbor set* $"SNN"(i, j) = Gamma(i) inter Gamma(j)$. The *SNN similarity* between distinct points $i != j$ is @liu2018snndpc
  $
    "Sim"(i, j) = cases(
      (|"SNN"(i, j)|^2) / (sum_(p in "SNN"(i, j)) (d_(i p) + d_(j p))) & "if" i, j in "SNN"(i, j),
      0 & "otherwise".
    )
  $
  Note that $i, j in "SNN"(i, j)$ holds iff $i$ and $j$ are *mutual $k$-nearest neighbors* ($j in "kNN"_k (i)$ and $i in "kNN"_k (j)$). Let $L_k (i)$ denote the $k$ points $j != i$ maximizing $"Sim"(i, j)$ and let $S_i = sum_(p in Gamma(i)) d_(i p)$ be the local $k$-NN distance sum. SNN-DPC defines local density $rho_i$ and scale-compensated higher-density separation $delta_i$ by
  $
    rho_i &= sum_(j in L_k (i)) "Sim"(i, j), \
    delta_i &= min_(j : rho_j > rho_i) [ d_(i j) (S_i + S_j) ],
  $
  with $delta_(i^*) = max_(j != i^*) delta_j$ for the global density peak $i^* = arg max_i rho_i$, and decision score $gamma_i = rho_i delta_i$.
]

Because $"Sim"(i, j) > 0$ implies $j in "kNN"_k (i)$, at most $k$ points have strictly positive similarity to $i$; hence summing over $L_k (i)$ equals summing $"Sim"(i, j)$ over all mutual $k$-nearest neighbors of $i$. Once $K$ peaks $cal(S) = {c_1, dots, c_K}$ are selected from the decision graph ${(rho_i, delta_i)}$, SNN-DPC partitions non-center points into two tiers: *inevitably subordinate* points reached from $cal(S)$ by BFS across strong shared-neighbor edges satisfying $|"SNN"(p, q)| >= k / 2$, and *possibly subordinate* points assigned iteratively via a consensus $k$-NN recognition matrix $W_(u, c) = |{v in Gamma(u) : ell(v) = c}|$.

#algbox(
  title: "SNN-DPC: Shared-Nearest-Neighbor Density Peaks Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n}$ with metric $d$, neighbor parameter $k >= 2$, cluster count $K$.],
  outputs: [Cluster centers $cal(S) = {c_1, dots, c_K}$ and partition $ell: {1, dots, n} -> {1, dots, K}$.],
)[
  1. *Neighborhood & SNN Graph:* Compute pairwise distances $D = (d_(i j))$, closed $k$-NN sets $Gamma(i)$, and local radii sums $S_i = sum_(p in Gamma(i)) d_(i p)$. For every mutual $k$-NN pair $(i, j)$, compute $|"SNN"(i, j)|$ and $"Sim"(i, j)$.
  2. *Compute $rho_i$, $delta_i$, and Select Peaks:* Evaluate $rho_i = sum_(j in L_k (i)) "Sim"(i, j)$, scale-compensated $delta_i = min_(j : rho_j > rho_i) d_(i j) (S_i + S_j)$, and $gamma_i = rho_i delta_i$. Select top-$K$ peaks $cal(S) = {c_1, dots, c_K}$ and initialize $ell(c_m) = m$.
  3. *Step 1 (Inevitably Subordinate BFS):* Push $cal(S)$ into FIFO queue $Q$. While $Q != emptyset$, pop $p$; for every unassigned $q in Gamma(p)$ satisfying $|"SNN"(p, q)| >= k / 2$, assign $ell(q) <- ell(p)$ and push $q$ onto $Q$.
  4. *Step 2 (Possibly Subordinate Consensus Assignment):* Let $U$ be the remaining unassigned set and $tilde(k) <- k$. While $U != emptyset$:
     (a) For each $u in U$ and $c in {1, dots, K}$, compute $W_(u, c) = |{v in Gamma_(tilde(k)) (u) : ell(v) = c}|$.
     (b) If $max_(u, c) W_(u, c) = 0$, increment $tilde(k) <- tilde(k) + 1$ and repeat (a).
     (c) Pick $(u^*, c^*) = arg max_(u in U, c) W_(u, c)$, assign $ell(u^*) <- c^*$, remove $u^*$ from $U$, and update $W$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical results explaining SNN-DPC's empirical improvements over DPC: tight metric sandwich bounds on the SNN similarity kernel, an exact scale-cancellation theorem proving that SNN-DPC eliminates DPC's multi-density peak suppression across arbitrary density contrasts and dimensions, termination and domino-suppression guarantees for two-step assignment, and a sharp feasible-$k$ breakdown condition.

#subsec("A", "Metric Sandwich Bounds on SNN Similarity")

#lemma(title: "Reciprocity, Symmetry, and Geodesic Sandwich Bounds for SNN Similarity", origin: "Original Proof")[
  For any metric $d$ on $cal(X)$ and any distinct points $i != j$ with mutual $k$-NN relation ($i, j in "SNN"(i, j)$), let $r_k (i) = max_(p in Gamma(i)) d_(i p)$ be the $k$-NN radius of $i$. Then $"Sim"(i, j) = "Sim"(j, i) > 0$ satisfies the sharp sandwich inequality
  $
    (|"SNN"(i, j)|) / (r_k (i) + r_k (j)) <= "Sim"(i, j) <= (|"SNN"(i, j)|) / d_(i j) <= (k + 1) / d_(i j),
  $
  with the upper bound $"Sim"(i, j) = |"SNN"(i, j)| / d_(i j)$ attained iff every shared neighbor $p in "SNN"(i, j)$ lies on a metric geodesic segment between $x_i$ and $x_j$ ($d_(i p) + d_(j p) = d_(i j)$).
]
#proof[
  Symmetry of $"SNN"(i, j) = Gamma(i) inter Gamma(j)$ and $d_(i p) + d_(j p)$ immediately implies $"Sim"(i, j) = "Sim"(j, i)$. Since $i, j in "SNN"(i, j)$ and $d_(i j) > 0$, the denominator contains both $p = i$ (giving $d_(i i) + d_(j i) = d_(i j)$) and $p = j$ (giving $d_(i j) + d_(j j) = d_(i j)$), so the sum is strictly positive and $"Sim"(i, j) > 0$.

  For every $p in "SNN"(i, j) = Gamma(i) inter Gamma(j)$, the metric triangle inequality gives $d_(i p) + d_(j p) >= d_(i j)$, while definition of $r_k$ gives $d_(i p) <= r_k (i)$ and $d_(j p) <= r_k (j)$, hence
  $
    d_(i j) <= d_(i p) + d_(j p) <= r_k (i) + r_k (j).
  $
  Summing over all $|"SNN"(i, j)|$ shared neighbors yields
  $
    |"SNN"(i, j)| dot d_(i j) <= sum_(p in "SNN"(i, j)) (d_(i p) + d_(j p)) <= |"SNN"(i, j)| dot (r_k (i) + r_k (j)).
  $
  Dividing $|"SNN"(i, j)|^2$ by these bounds and noting $|"SNN"(i, j)| <= |Gamma(i)| = k + 1$ proves the sandwich inequality. Equality in the middle bound holds iff $d_(i p) + d_(j p) = d_(i j)$ for every $p in "SNN"(i, j)$.
]

Lemma 1 highlights the dual nature of SNN similarity: writing $bar(d)_"SNN" (i, j) = (1 / |"SNN"(i, j)|) sum_(p in "SNN"(i, j)) (d_(i p) + d_(j p))$, we have $"Sim"(i, j) = |"SNN"(i, j)| / bar(d)_"SNN" (i, j)$, simultaneously rewarding topological neighborhood agreement and penalizing metric stretch.

#subsec("B", "Scale Invariance and Exact Multi-Density Cancellation Theorem")

#theorem(title: "Global Scale Invariance and Multi-Density Contrast Cancellation", origin: "Original Derivation")[
  1. *(Global Similarity Invariance)* Under any Euclidean similarity $T(x) = alpha Q x + b$ ($alpha > 0$, $Q in O(d)$), SNN-DPC quantities scale homogeneously as
  $
    rho'_i = alpha^(-1) rho_i, quad delta'_i = alpha^2 delta_i, quad gamma'_i = alpha gamma_i.
  $
  Hence the exact rankings of $rho_i$, $delta_i$, and $gamma_i$, as well as both assignment steps, are strictly invariant to global scaling and rigid isometries with zero scale-parameter retuning.
  2. *(Exact Multi-Density Cancellation)* Let $P_0 subset bb(R)^d$ be any finite template of size $m > k$ with unique SNN-density peak $p_0$, peak density $rho_0 = rho(p_0)$, local distance sum $S_0 = S_(p_0)$, and maximum sum $S_"max" = max_(p in P_0) S_p$. Form a two-density mixture $cal(X) = C_1 union C_2$ where $C_1 = alpha_1 P_0$ and $C_2 = alpha_2 P_0 + v$ with arbitrary linear scale ratio $beta = alpha_2 / alpha_1 > 1$ (volume density contrast $beta^d$) and separation $D = "dist"(C_1, C_2) > alpha_2 "diam"(P_0)$. Then the sparse peak $c_2 in C_2$ satisfies the *scale-free* lower bound
  $
    gamma(c_2) >= D rho_0 S_0,
  $
  independently of $beta$ and dimension $d$. Moreover, whenever $D > 2 alpha_2 "diam"(P_0) (S_"max" / S_0)$, the two true peaks $c_1$ and $c_2$ strictly maximize $gamma_i$ over all $x in cal(X)$ for *every* density ratio $beta in (1, infinity)$.
]
#proof[
  *(1)* Because $d(T(x_i), T(x_j)) = alpha d_(i j)$, every $k$-NN set $Gamma(i)$ and shared set $"SNN"(i, j)$ is identical. Therefore $"Sim"'(i, j) = alpha^(-1) "Sim"(i, j)$, implying $rho'_i = alpha^(-1) rho_i$ and preserving the exact strict ordering $rho'_j > rho'_i <=> rho_j > rho_i$. Meanwhile $S'_i = alpha S_i$, so $d'_(i j) (S'_i + S'_j) = alpha^2 d_(i j) (S_i + S_j)$, giving $delta'_i = alpha^2 delta_i$ and $gamma'_i = (alpha^(-1) rho_i)(alpha^2 delta_i) = alpha gamma_i$. Because $alpha > 0$, all argmax/argmin rankings and neighborhood counts $|"SNN"(p, q)|$ are completely unchanged.

  *(2)* Because $D > alpha_2 "diam"(P_0) >= alpha_1 "diam"(P_0)$ and $|C_1| = |C_2| = m > k$, every point's closed $k$-NN set $Gamma(x)$ lies strictly inside its own cluster ($C_1$ or $C_2$). Applying Part (1) separately inside each isolated cluster yields, for every $x = alpha_1 tilde(x) in C_1$ and $y = alpha_2 tilde(y) + v in C_2$:
  $
    rho(x) = alpha_1^(-1) rho_(P_0) (tilde(x)), &quad S(x) = alpha_1 S_(P_0) (tilde(x)), \
    rho(y) = alpha_2^(-1) rho_(P_0) (tilde(y)), &quad S(y) = alpha_2 S_(P_0) (tilde(y)).
  $
  Since $beta = alpha_2 / alpha_1 > 1$, the peak $c_1 in C_1$ achieves $rho(c_1) = alpha_1^(-1) rho_0 > alpha_2^(-1) rho_0 = rho(c_2)$, and every $y in C_2$ satisfies $rho(y) <= rho(c_2)$. Consequently, the higher-density set ${z in cal(X) : rho(z) > rho(c_2)}$ is a nonempty subset of $C_1$ (containing $c_1$). Letting $j^* in C_1$ attain the minimum in $delta(c_2)$, we have $d(c_2, j^*) >= D$ and $S(j^*) > 0$, hence
  $
    delta(c_2) = d(c_2, j^*) (S(c_2) + S(j^*)) >= D S(c_2) = D alpha_2 S_0.
  $
  Multiplying by $rho(c_2) = alpha_2^(-1) rho_0$, the factor $alpha_2$ cancels identically:
  $
    gamma(c_2) = rho(c_2) delta(c_2) >= (alpha_2^(-1) rho_0) (D alpha_2 S_0) = D rho_0 S_0.
  $
  Now consider any non-peak point $u in C_1 without {c_1}$. Because $c_1 in C_1$ has strictly higher density than $u$ and $d(u, c_1) <= alpha_1 "diam"(P_0)$,
  $
    delta(u) <= d(u, c_1) (S(u) + S(c_1)) <= 2 alpha_1^2 "diam"(P_0) S_"max",
  $
  and $rho(u) <= alpha_1^(-1) rho_0$, so $gamma(u) <= 2 alpha_1 "diam"(P_0) rho_0 S_"max"$. By identical reasoning, every non-peak $w in C_2 without {c_2}$ satisfies $gamma(w) <= 2 alpha_2 "diam"(P_0) rho_0 S_"max"$. Since $alpha_1 < alpha_2$, the condition $D > 2 alpha_2 "diam"(P_0) (S_"max" / S_0)$ implies $gamma(c_2) > max_(z in.not {c_1, c_2}) gamma(z)$. Finally, $delta(c_1) = max_(j != c_1) delta(j) >= delta(c_2)$ and $rho(c_1) > rho(c_2)$ guarantee $gamma(c_1) > gamma(c_2)$, so ${c_1, c_2}$ are strictly the top-2 scores in $gamma$.
]

Theorem 2 reveals the precise algebraic mechanism of SNN-DPC's success on varied-density data: whereas standard DPC's density decays exponentially as $beta^(-d)$ while $delta$ remains $cal(O)(D)$, SNN-DPC's density decays only linearly as $beta^(-1)$ and is *exactly cancelled* by the linear factor $alpha_2$ inside $S(c_2) = alpha_2 S_0$!

#subsec("C", "Convergence and Domino-Resistance of Two-Step Assignment")

#theorem(title: "Completeness and Plurality Stability of Two-Step Assignment", origin: "Original Proof")[
  Let $cal(S) = {c_1, dots, c_K}$ be any nonempty center set ($1 <= K <= n$). Then:
  1. *(Finite Termination & Completeness)* Step 1 terminates in at most $n$ queue pops, and Step 2 assigns every remaining point in $U$ in exactly $|U|$ iterations (with neighborhood expansion $tilde(k) <= n - 1$), producing a complete partition of $cal(X)$.
  2. *(Core Bridge Resistance)* Step 1 never crosses between two disjoint subsets $A, B$ unless there exists an edge $(p, q) in A times B$ with $q in Gamma(p)$ and $|"SNN"(p, q)| >= k / 2$.
  3. *(Step 2 Plurality Robustness)* In Step 2, an unassigned point $u in U$ whose $k$-NN set contains $a$ points of its true cluster $c$ and at most $b$ points of any rival cluster $c' != c$ retains its correct assignment $ell(u) = c$ under any adversarial flip of up to $floor((a - b - 1) / 2)$ neighbor labels, preventing single-edge domino cascades.
]
#proof[
  *(1)* In Step 1, each point enters the FIFO queue $Q$ at most once (upon transitioning from unassigned to assigned), so at most $n$ pops occur. In Step 2, since $K >= 1$, at least one point is already assigned; if $max_(u, c) W_(u, c) = 0$ at radius $tilde(k)$, incrementing $tilde(k)$ strictly enlarges $Gamma_(tilde(k)) (u)$ until $tilde(k) <= n - 1$, where $Gamma_(n-1) (u) = cal(X)$ necessarily contains $cal(S)$ and yields $W_(u, c) >= 1$. Each assignment strictly decreases $|U|$ by $1$, terminating in $|U|$ steps with every vertex labeled.

  *(2)* Step 1 traverses exclusively directed edges $p -> q$ with $q in Gamma(p)$ and $|"SNN"(p, q)| >= k / 2$. By induction on BFS distance from $cal(S)$, a label originating in $A$ can reach $B$ only across a boundary edge satisfying this shared-neighbor majority criterion.

  *(3)* Flipping one neighbor of $u$ from cluster $c$ to a rival cluster $c'$ decreases $W_(u, c)$ by at most $1$ and increases $W_(u, c')$ by at most $1$, reducing the margin $W_(u, c) - W_(u, c')$ by at most $2$. Starting from margin $a - b >= 1$, flipping $t <= floor((a - b - 1) / 2)$ neighbors preserves $W'_(u, c) - W'_(u, c') >= (a - b) - 2 t >= 1 > 0$, guaranteeing that $c$ remains the strict plurality winner for $u$.
]

#subsec("D", "Neighborhood-Size Feasibility Window and Small-Cluster Breakdown")

#proposition(title: "Small-Cluster Spillover Breakdown When k Exceeds Cluster Cardinality", origin: "Original Counterexample")[
  Let $cal(X) = C_1 union C_2$ consist of two clusters of sizes $n_1 = |C_1|$ and $n_2 = |C_2|$ separated by arbitrary distance $D > 0$. Then:
  1. If $n_1 = k$ and $n_2 <= k / 2$, there exists a cross-cluster pair $(p, q^*) in C_1 times C_2$ with $q^* in Gamma(p)$ and $|"SNN"(p, q^*)| >= k + 2 - n_2 > k / 2$ for *every* $D > 0$, causing Step 1 BFS to flood unconditionally from $C_1$ into $C_2$.
  2. Even if $n_1 > k$ while $n_2 <= k / 2$, every point $q in C_2$ absorbs $k + 1 - n_2$ points of $C_1$ into $Gamma(q)$, collapsing its SNN local density as $rho(q) = cal(O)(n_2 k / D) -> 0$ as $D -> infinity$.
  Hence $k < min_m |C_m|$ is necessary to preserve intra-cluster $k$-NN purity.
]
#proof[
  Place all $n_1$ points of $C_1$ inside a ball of radius $eta << D$ around $bold(0)$ and all $n_2$ points of $C_2$ inside a ball of radius $eta$ around $(D, 0, dots, 0)$. Every cross-cluster pair $(p, q) in C_1 times C_2$ satisfies $D - 2 eta <= d_(p q) <= D + 2 eta$, whereas all intra-cluster distances are $<= 2 eta$.

  *(1)* When $n_1 = k$ and $n_2 <= k / 2$, every $p in C_1$ has only $k$ points inside $C_1$ (including itself), so its $(k + 1)$-element closed $k$-NN set $Gamma(p)$ equals $C_1 union {q^*}$ for the point $q^* in C_2$ nearest to $p$. Meanwhile, $q^*$ has only $n_2$ points in $C_2$, so $Gamma(q^*)$ must contain all of $C_2$ plus $k + 1 - n_2$ points of $C_1$. Intersecting the two sets gives
  $
    |"SNN"(p, q^*)| &= |Gamma(p) inter Gamma(q^*)| \
    &= (k + 1 - n_2) + |{q^*}| = k + 2 - n_2 > k / 2
  $
  since $n_2 <= k / 2$. Because $q^* in Gamma(p)$ and $|"SNN"(p, q^*)| > k / 2$, Step 1 BFS immediately propagates $C_1$'s label across the gap $D$ to $q^*$ regardless of how large $D$ is.

  *(2)* When $n_1 > k$ and $n_2 <= k / 2$, every $p in C_1$ satisfies $Gamma(p) subset C_1$, so no point of $C_1$ is a mutual $k$-NN of any $q in C_2$. Thus each $q in C_2$ has at most $n_2 - 1$ mutual $k$-NNs $q' in C_2 without {q}$. Because $"SNN"(q, q')$ contains $k + 1 - n_2$ shared neighbors in $C_1$ at distance $>= D - 2 eta$ from both $q$ and $q'$,
  $
    sum_(r in "SNN"(q, q')) (d_(q r) + d_(q' r)) >= 2 (k + 1 - n_2) (D - 2 eta).
  $
  Using $|"SNN"(q, q')| <= k + 1$ and $k + 1 - n_2 > k / 2$ yields $"Sim"(q, q') <= (k + 1)^2 / (k (D - 2 eta))$ and therefore $rho(q) <= (n_2 - 1) (k + 1)^2 / (k (D - 2 eta)) = cal(O)(n_2 k / D) -> 0$ as $D -> infinity$, severing $C_2$'s density from its true internal scale $eta$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Computing the pairwise distance matrix $D$ costs $cal(O)(n^2 d)$ time and $cal(O)(n^2)$ space; extracting $k$-NN lists $Gamma(i)$ takes $cal(O)(n^2 log k)$ (or $cal(O)(n^2)$ via linear selection). Because $"Sim"(i, j)$ is nonzero only on mutual $k$-NN pairs—of which there are at most $n k / 2$—computing set intersections $|Gamma(i) inter Gamma(j)|$ and sums costs only $cal(O)(n k^2)$ using sorted lists or hash sets (compared to $cal(O)(n^2 k)$ for naive all-pairs intersections @liu2018snndpc). Evaluating $rho_i$ and $S_i$ requires $cal(O)(n k)$; evaluating $delta_i = min_(j : rho_j > rho_i) d_(i j) (S_i + S_j)$ using precomputed $S_i$ takes exactly $cal(O)(n^2)$ time. Step 1 BFS checks at most $n k$ directed edges (with precomputed $|"SNN"(p, q)|$) in $cal(O)(n k)$ time, while Step 2 maintains the recognition matrix $W$ in $cal(O)(n k log n)$ time (worst-case $cal(O)(n^2)$ if fallback $tilde(k)$-expansion triggers). Thus overall complexity from a distance matrix is $cal(O)(n^2 + n k^2)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(n^2 d + n k^2)$ from coordinates; $cal(O)(n^2 + n k^2)$ from precomputed distance matrix.],
    [Space Complexity], [$cal(O)(n^2)$ for distance matrix $D$ plus $cal(O)(n k)$ for $k$-NN lists and sparse SNN graph.],
    [Inductive Bias], [Multi-density, non-convex unimodal basins connected by dense mutual-$k$-NN cores.],
    [Hyperparameters], [Neighborhood size $k$ (requiring $2 <= k < min_m |C_m|$) and cluster count $K$.],
    [Invariances], [Exact global Euclidean similarity invariance (Theorem 2); scale-free multi-density cancellation.],
    [Main Bottlenecks], [$cal(O)(n^2)$ memory ceiling for $D$, sensitivity when $k >= n_"min"$, and $k$-NN hubness in very high $d$.],
  ),
  caption: [Operational and complexity profile of SNN-DPC.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multi-Density Geospatial & Urban-Rural Segmentation:* Spatial point clouds combining dense metropolitan centers with sparse outlying towns defeat DPC's fixed $d_c$. Theorem 2 guarantees that SNN-DPC's local-$S_i$ scaling surfaces sparse town peaks cleanly on the decision graph alongside dense downtown peaks without manual multi-scale tuning @liu2018snndpc.
2. *Gene Expression & Single-Cell Subpopulation Discovery:* Transcriptomic profiles exhibit both heterogeneous cell-type abundances (multi-density) and moderate-to-high dimensionality where raw Euclidean radii concentrate. Mutual $k$-NN overlap captures local co-expression manifold topology far more reliably than global $d_c$ balls @ertoz2003snn.
3. *Hyperspectral Remote-Sensing Pixel Classification:* Land-cover spectra (water, urban concrete, varied vegetation canopies) form curved manifolds with vastly different spectral variance. SNN similarity suppresses spectral-intensity scale differences while Step 2 consensus voting reduces salt-and-pepper boundary misclassification.
4. *Complex Shape & Manifold Recognition in Computer Vision:* Benchmarks such as Aggregation, Jain, Spiral, and Flame combine touching multi-density blobs with intertwined spirals @liu2018snndpc. Step 1's $|"SNN"| >= k / 2$ BFS traces intertwined spirals reliably while preventing premature bridging across narrow inter-manifold gaps.
5. *Cybersecurity & Multi-Scale Behavioral Traffic Profiling:* Enterprise network flows mix ultra-dense routine protocol clusters with sparse low-and-slow exfiltration clusters; SNN-DPC elevates sparse behavioral peaks without splitting high-volume baseline traffic.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Key failure modes of SNN-DPC include: (i) *small-cluster $k$-spillover* (Proposition 4), where setting $k >= n_"min"$ contaminates small clusters' $k$-NN sets with neighbors from larger clusters; (ii) *quadratic memory exhaustion* above $n approx 5 times 10^4$ from storing $D$ and computing $delta_i$; (iii) *high-dimensional hubness*, where certain points appear in disproportionately many $k$-NN lists despite ambient concentration @beyer1999nearest, inflating spurious SNN intersections; and (iv) *multi-peak flat manifolds*, inherited from DPC whenever a single uniform ring has two weak local density maxima.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Multi-Density Non-Convex Manifolds], [Ideal: scale-free cancellation (Thm. 2) plus BFS core tracing.], [Keep $k < min_m |C_m|$ (Prop. 4).],
    [Moderate-Dimensional Omics / Spectra], [Strong: SNN reciprocity resists moderate distance concentration.], [Check for $k$-NN hubness when $d > 50$.],
    [Datasets With Tiny Rare Clusters ($n_"min" < 10$)], [Fragile: $k$-NN sets spill into adjacent dense clusters.], [Use HDBSCAN @campello2013hdbscan or reduce $k$ carefully.],
    [Massive Datasets ($n > 5 times 10^4$)], [Blocked by $cal(O)(n^2)$ distance and $delta_i$ sweep.], [Use MR-DBSCAN (Alg. 49) or ANN-accelerated SNN.],
    [Ultra-Low-Latency Homoscedastic Data], [Good quality, but $k$-NN intersection adds $~2.4 times$ latency.], [Use plain DPC (Alg. 51) if density is uniform.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SNN-DPC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Unification of SNN Connectivity Graphs and Density Peaks Mode Detection", origin: "Original Derivation")[
  SNN-DPC structurally bridges SNN-DBSCAN @ertoz2003snn, RNN-DBSCAN @bryant2018rnndbscan, and DPC @rodriguez2014dpc:
  1. Step 1 of SNN-DPC traverses connected components of the directed threshold graph $G_(k / 2)^"SNN"$ with edge criterion $q in Gamma(p)$ and $|"SNN"(p, q)| >= k / 2$, which is a seeded core-subgraph of SNN-DBSCAN's shared-neighbor graph restricted to components reachable from the decision-graph peaks $cal(S)$.
  2. The reciprocity gate $i, j in "SNN"(i, j)$ restricts positive SNN similarities exclusively to the symmetric mutual-$k$-NN graph ($j in "kNN"_k (i)$ and $i in "kNN"_k (j)$), sharing the exact reverse-nearest-neighbor density core of RNN-DBSCAN @bryant2018rnndbscan while replacing threshold-based component merging with DPC's peak-rooted separation $delta_i$.
]
#proof[
  For (1), SNN-DBSCAN @ertoz2003snn defines edge weight $w_(p q) = |"SNN"(p, q)|$ and extracts connected components of core points linked by edges with $w_(p q) >= "Eps"_"SNN"$. Setting $"Eps"_"SNN" = k / 2$ and initiating BFS solely from the $K$ density peaks $cal(S)$ yields Step 1's assigned set as the union of connected components of $G_(k/2)^"SNN"$ intersecting $cal(S)$ (with Step 2 absorbing all remaining points rather than marking them noise). For (2), $i, j in Gamma(i) inter Gamma(j)$ holds iff $i in Gamma(j)$ and $j in Gamma(i)$, which for $i != j$ is identical to mutual $k$-NN membership—precisely the bidirectional $k$-NN / RNN adjacency underlying RNN-DBSCAN @bryant2018rnndbscan.
]

Compared with VDBSCAN (Alg. 50, @liu2007vdbscan), OPTICS @ankerst1999optics, and HDBSCAN @campello2013hdbscan, SNN-DPC achieves multi-density adaptivity without iterative deflation passes or minimum spanning trees: local scale normalization is embedded algebraically into the pair $(rho_i, delta_i)$ itself (Theorem 2).

*Empirical Benchmark Analysis:* On the ClusterOverview four-dataset benchmark (Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), SNN-DPC records:
- *Iris:* Silhouette $= 0.476$, $"ARI" = 0.450$;
- *Customers:* Silhouette $= 0.297$;
- *Multi-Shape:* Silhouette $= 0.712$ ($"ARI" =$ N/A);
- *High-D:* Silhouette $= 0.297$ ($"ARI" =$ N/A);
- *Mean Runtime:* $2.136$ ms.

Direct head-to-head comparison with plain DPC (Alg. 51) vividly confirms Section II's theory. On *Iris*, SNN-DPC improves Silhouette by $+0.198$ ($0.476$ vs. $0.278$), demonstrating that two-step SNN assignment and consensus voting produce substantially cleaner, more compact cluster geometry than single-hop tree assignment, while the drop in $"ARI"$ ($0.450$ vs. $0.569$) illustrates the trade-off of Proposition 4 when a fixed $k$ bridges overlapping Versicolor/Virginica margins. On *Customers*, Silhouette increases from $0.247$ to $0.297$ ($+0.050$), directly reflecting Theorem 2's scale-compensated handling of varied spending densities, though still trailing pure density-reachability algorithms MR-DBSCAN ($0.663$) and VDBSCAN ($0.589$) that reject sparse customer outliers as noise rather than forcing complete assignment in Step 2. On *Multi-Shape* ($0.712$ vs. $0.718$) and *High-D* ($0.297$ vs. $0.296$), SNN-DPC matches DPC while completely eliminating $d_c$ tuning. Finally, its mean runtime ($2.136$ ms) incurs only a $2.4 times$ overhead over plain DPC ($0.888$ ms) from $k$-NN intersections and two-step assignment, remaining nearly $2 times$ faster than VDBSCAN ($4.156$ ms) and $4.4 times$ faster than MR-DBSCAN ($9.373$ ms).
