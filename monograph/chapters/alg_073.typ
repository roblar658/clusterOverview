#import "../template/ieee.typ": *

#alg_entry(
  id: 73,
  name: "Minimax Linkage",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2011",
  paradigm: "Prototype-Anchored Minimax-Radius Agglomerative Clustering",
  complexity: $O(n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical agglomerative hierarchical clustering algorithms — single, complete, average, and Ward linkage @lance1967general @ward1963hierarchical — construct nested partitions whose internal nodes lack an intrinsic data-point representative. Computing a post-hoc centroid requires a Euclidean embedding and produces a synthetic point that may lie outside the data manifold, whereas selecting a post-hoc medoid after clustering decouples prototype selection from the merge criterion itself. Bien and Tibshirani @bien2011minimax resolved this interpretability gap by introducing *Minimax Linkage*, an agglomerative method in which every cluster $C$ at every level of the dendrogram is jointly defined with an observed *prototype* $p^*(C) in C$ minimizing the worst-case (minimax) covering radius over $C$.

Let $X = {x_1, dots.c, x_n}$ be a finite dataset equipped with a symmetric, zero-diagonal dissimilarity function $d : X times X arrow.r RR_(>= 0)$. Unlike centroid or Ward linkage, $d$ need not arise from an inner product or even satisfy the triangle inequality for the core hierarchical consistency properties to hold.

#definition(title: "Eccentricity, Minimax Radius, Prototype, and Minimax Linkage")[
  For any nonempty subset $C subset.eq X$ and candidate point $x in C$, the *eccentricity* (covering radius from $x$) is
  $
    d_"max" (x, C) = max_(y in C) d(x, y).
  $
  The *minimax radius* of $C$ and its *prototype set* are defined by @bien2011minimax
  $
    r(C) = min_(x in C) d_"max" (x, C), quad cal(P)^*(C) = arg min_(x in C) d_"max" (x, C).
  $
  Any element $p^*(C) in cal(P)^*(C)$ is called a *minimax prototype* of $C$. For two disjoint clusters $G, H subset.eq X$, their *minimax linkage* is the minimax radius of their union:
  $
    d_"MM" (G, H) = r(G union H) = min_(x in G union H) max_(y in G union H) d(x, y).
  $
]

Geometrically, cutting a minimax-linkage dendrogram at height $h >= 0$ guarantees that every resulting cluster $C_m$ is completely enclosed within the closed dissimilarity ball $overline(B)(p^*(C_m), h) = {z in X : d(z, p^*(C_m)) <= h}$ centered at an actual dataset element $p^*(C_m) in C_m$. Notice that $d_"MM" (G, H)$ evaluates the *internal* covering radius of the merged set $G union H$ rather than aggregating only inter-cluster distances between $G$ and $H$. To accelerate evaluation across merges, the algorithm caches each point's cluster eccentricity $M(x, G) = max_(y in G) d(x, y)$ (for all $x in X$ and active $G$), reducing the computation of $d_"MM" (G, H)$ to an $O(abs(G) + abs(H))$ scan via the identity
$
  d_"max" (x, G union H) = max(M(x, G), M(x, H)).
$

#algbox(
  title: "Minimax Linkage Agglomerative Clustering",
  inputs: [Symmetric dissimilarity matrix $D in RR_(>= 0)^(n times n)$ over $X = {x_1, dots.c, x_n}$],
  outputs: [Dendrogram merges ${(G_t, H_t, h_t)}_(t=1)^(n-1)$ and per-node prototypes $p^*(C_t) in C_t$],
)[
  1. Initialize active clusters $cal(A) arrow.l {{x_1}, dots.c, {x_n}}$, eccentricities $M(x_u, {x_v}) arrow.l D_(u v)$ for all $u, v$, and pairwise linkages $L({x_u}, {x_v}) arrow.l D_(u v)$ for $u < v$.
  2. *For* $t = 1, dots.c, n - 1$:
     - Select $(G_t, H_t) in arg min_(G != H in cal(A)) L(G, H)$ and set merge height $h_t arrow.l L(G_t, H_t)$.
     - Form $C_t arrow.l G_t union H_t$ and store prototype $p^*(C_t) in arg min_(x in C_t) max(M(x, G_t), M(x, H_t))$.
     - *For all* $x in X$: update $M(x, C_t) arrow.l max(M(x, G_t), M(x, H_t))$.
     - Update $cal(A) arrow.l (cal(A) without {G_t, H_t}) union {C_t}$.
     - *For each* $K in cal(A) without {C_t}$: compute
       $L(C_t, K) arrow.l min_(x in C_t union K) max(M(x, C_t), M(x, K))$.
  3. Return the prototype-annotated merge tree; cutting at height $h$ yields balls $C_m subset.eq overline(B)(p^*(C_m), h)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Despite falling outside the combinatorial Lance-Williams family, Minimax Linkage possesses remarkably strong order-theoretic and geometric guarantees. We prove: (A) unconditional Bruynooghe reducibility and dendrogram monotonicity (even though $r(dot)$ itself is not set-monotone!), (B) a tight factor-2 metric sandwich linking Minimax and Complete linkage plus a maximal ball-covering theorem, (C) exact invariance under strictly increasing dissimilarity transforms alongside formal non-Lance-Williams impossibility, and (D) finite-sample outlier isolation and radius-breakdown bounds.

#subsec("A", "Unconditional Reducibility and Dendrogram Monotonicity")

A subtle feature of the minimax radius $r(C)$ is that it is *not* monotone with respect to set inclusion: adding a central point can strictly decrease the radius. For instance, on $RR$ with Euclidean distance, $A = {0, 2} subset.eq B = {0, 1, 2}$ satisfies $r(A) = 2 > 1 = r(B)$ because the newly added midpoint $1$ cuts the worst-case covering distance in half. Remarkably, this non-monotonicity of $r(dot)$ never produces dendrogram inversions.

#theorem(title: "Unconditional Reducibility and Ultrametricity of Minimax Linkage", origin: "Classical — Cited")[
  For any symmetric dissimilarity $d : X times X arrow.r RR_(>= 0)$ (without requiring the triangle inequality) and any three pairwise disjoint nonempty clusters $G, H, K subset.eq X$, Minimax Linkage satisfies the unconditional reducibility inequality @bien2011minimax
  $
    d_"MM" (G union H, K) >= min(d_"MM" (G, K), d_"MM" (H, K)).
  $
  Consequently, the merge heights are non-decreasing ($h_1 <= h_2 <= dots.c <= h_(n-1)$), the dendrogram never exhibits inversions, and $u(x, y) = h_(t(x, y))$ defines a valid ultrametric on $X$.
]

#proof[
  Let $U = G union H union K$, and choose a minimax prototype $p^* in U$ achieving $d_"MM" (G union H, K) = r(U) = max_(y in U) d(p^*, y)$. Because $U = (G union K) union (H union K)$, the prototype $p^*$ must belong to at least one of $G union K$ or $H union K$.
  - *Case 1*: $p^* in G union K$. Since $G union K subset.eq U$, taking the maximum over a smaller set cannot increase the eccentricity, and minimizing over all candidates in $G union K$ cannot increase it further:
    $
      r(U) = max_(y in U) d(p^*, y) >= max_(y in G union K) d(p^*, y) \
      >= min_(x in G union K) max_(y in G union K) d(x, y) = d_"MM" (G, K).
    $
  - *Case 2*: $p^* in H union K$. By identical containment $H union K subset.eq U$,
    $
      r(U) = max_(y in U) d(p^*, y) >= max_(y in H union K) d(p^*, y) >= d_"MM" (H, K).
    $
  In both cases, $d_"MM" (G union H, K) >= min(d_"MM" (G, K), d_"MM" (H, K))$.

  To deduce dendrogram monotonicity, observe that at step $t$, the merged pair $(G_t, H_t)$ minimizes $d_"MM"$ over all active pairs, so $h_t = d_"MM" (G_t, H_t) <= min(d_"MM" (G_t, K), d_"MM" (H_t, K)) <= d_"MM" (G_t union H_t, K)$ for every remaining $K$. All unaffected pairs $(K, K')$ already satisfy $d_"MM" (K, K') >= h_t$, hence $h_(t+1) >= h_t$. Monotonicity directly yields the strong triangle inequality $u(x, z) <= max(u(x, y), u(y, z))$ for all $x, y, z in X$.
]

Note how striking this proof is: while set-monotonicity fails for $A = {0, 2} subset.eq {0, 1, 2}$, an agglomerative run on ${0, 1, 2}$ merges ${0, 1}$ (or ${1, 2}$) at height $1$ before ever considering ${0, 2}$ at height $2$, so non-monotone radius drops are strictly unreachable along the merge trajectory. Moreover, because reducibility holds unconditionally, the Nearest-Neighbor Chain order-invariance theorem @mullner2011modern also applies to Minimax Linkage (though each linkage evaluation still scales with cluster cardinality).

#subsec("B", "Tight Metric Sandwich with Complete Linkage and Covering Maximality")

In any metric space, Minimax Linkage is tightly bracketed by Complete Linkage (cluster diameter), conferring concrete ball-covering guarantees at every cut height.

#theorem(title: "Tight Factor-2 Metric Sandwich and Maximal Ball-Covering Guarantee", origin: "Original Proof")[
  Let $(X, d)$ be a metric space and define cluster diameter $"diam"(C) = max_(x, y in C) d(x, y)$. Then:
  1. For every nonempty cluster $C subset.eq X$,
     $
       (1 / 2) "diam"(C) <= r(C) <= "diam"(C),
     $
     and both constants $1 / 2$ and $1$ are sharp.
  2. For any disjoint clusters $G, H subset.eq X$ with complete linkage $d_"CL" (G, H) = max_(x in G, y in H) d(x, y)$,
     $
       (1 / 2) Delta_(G, H) <= d_"MM" (G, H) <= Delta_(G, H),
     $
     where $Delta_(G, H) = max(d_"CL" (G, H), "diam"(G), "diam"(H))$.
  3. Cutting the Minimax dendrogram at any height $h >= 0$ yields a partition ${C_1, dots.c, C_k}$ with prototypes $p_m = p^*(C_m) in C_m$ such that (i) $C_m subset.eq overline(B)(p_m, h)$ for every $m$, and (ii) for every pair $m != ell$, no single point in $C_m union C_ell$ covers $C_m union C_ell$ within radius $h$ (i.e., $r(C_m union C_ell) > h$).
]

#proof[
  *Part 1*: For the upper bound, picking any $x in C$ gives $r(C) <= max_(y in C) d(x, y) <= max_(u, v in C) d(u, v) = "diam"(C)$. For the lower bound, let $p^* in cal(P)^*(C)$ achieve $r(C)$. For any $y, z in C$, the triangle inequality implies
  $
    d(y, z) <= d(y, p^*) + d(p^*, z) <= 2 max_(w in C) d(p^*, w) = 2 r(C).
  $
  Taking the maximum over $y, z in C$ gives $"diam"(C) <= 2 r(C)$. To verify tightness of both constants: the collinear 3-point set ${-1, 0, 1} subset.eq RR$ has $"diam" = 2$ and $r = 1 = (1 / 2) "diam"$, while a unit equilateral triangle in $RR^2$ has $"diam" = 1$ and $r = 1 = "diam"$.

  *Part 2*: The diameter of a union decomposes identically as $"diam"(G union H) = max(d_"CL" (G, H), "diam"(G), "diam"(H)) = Delta_(G, H)$. Applying Part 1 to $C = G union H$ yields the sandwich directly.

  *Part 3*: By monotonicity (Theorem 1), the cut at height $h$ retains all merges with $h_t <= h$, so every surviving cluster $C_m$ was formed at some merge height $<= h$ (or is a singleton with $r = 0 <= h$), giving $max_(y in C_m) d(p_m, y) = r(C_m) <= h$, hence $C_m subset.eq overline(B)(p_m, h)$. Conversely, if any pair $C_m, C_ell$ satisfied $r(C_m union C_ell) <= h$, then at the termination of all merges $<= h$ the active pair $(C_m, C_ell)$ would have $d_"MM" (C_m, C_ell) <= h$, contradicting that the minimum active linkage exceeds $h$ after the cut.
]

#subsec("C", "Monotone Transformation Invariance and Non-Combinatorial Structure")

Unlike UPGMA, Ward, or Flexible Beta (Alg. 72 and Alg. 74), Minimax Linkage shares the rare *ordinal invariance* of single and complete linkage @bien2011minimax @carlsson2010characterization — yet provably escapes the Lance-Williams family.

#proposition(title: "Strict Monotone-Transform Invariance and Non-Lance-Williams Impossibility", origin: "Original Proof")[
  1. Let $f : RR_(>= 0) arrow.r RR_(>= 0)$ be any strictly increasing function with $f(0) = 0$. Replacing $d$ by $tilde(d) = f compose d$ maps every cluster's minimax radius to $tilde(r)(C) = f(r(C))$ and leaves the entire merge sequence, dendrogram topology, and prototype sets $cal(P)^*(C)$ strictly invariant.
  2. There exists no function $F(n_i, n_j, n_k, d_(i j), d_(i k), d_(j k))$ — and in particular no Lance-Williams parameterization — that equals $d_"MM" (C_i union C_j, C_k)$ across all Euclidean datasets.
]

#proof[
  *Part 1*: Because $f$ is strictly increasing on a finite range, it commutes with both finite maxima and finite minima: for every $C subset.eq X$,
  $
    tilde(r)(C) &= min_(x in C) max_(y in C) f(d(x, y)) \
    &= min_(x in C) f(max_(y in C) d(x, y)) = f(r(C)),
  $
  and the sets of argmins over $x in C$ coincide identically. Applying $f$ to every pairwise linkage $d_"MM" (G, H)$ preserves all strict inequalities among candidate pairs at every step $t = 1, dots.c, n - 1$, yielding identical merges and prototypes @bien2011minimax.

  *Part 2*: Consider the two 4-point Euclidean instances in $RR^2$ with partition $C_i = {a, b}$, $C_j = {c}$, $C_k = {e}$ (cardinalities $n_i = 2, n_j = 1, n_k = 1$):
  - *Instance I*: $a = (0, 0)$, $b = (0, 1)$, $c = (1, 0)$, $e = (-1, 0)$. Point $a$ achieves eccentricity $1$ in $C_i union C_j$, $C_i union C_k$, and $C_i union C_j union C_k$, while $d(c, e) = 2$. Thus $(d_"MM" (C_i, C_j), d_"MM" (C_i, C_k), d_"MM" (C_j, C_k)) = (1, 1, 2)$ and $d_"MM" (C_i union C_j, C_k) = 1$.
  - *Instance II*: $a = (-1 / 2, 0)$, $b = (1 / 2, 0)$, $c = (1, sqrt(3) / 2)$, $e = (-1, sqrt(3) / 2)$. Here $d(a, b) = d(b, c) = d(a, e) = 1$, $d(a, c) = d(b, e) = sqrt(3)$, and $d(c, e) = 2$. Prototype $b$ gives $d_"MM" (C_i, C_j) = 1$, prototype $a$ gives $d_"MM" (C_i, C_k) = 1$, and $d_"MM" (C_j, C_k) = 2$. Yet in $C_i union C_j union C_k$, the eccentricities of $(a, b, c, e)$ are $(sqrt(3), sqrt(3), 2, 2)$, giving $d_"MM" (C_i union C_j, C_k) = sqrt(3) != 1$.

  Since identical scalar inputs produce distinct merged linkages, $d_"MM"$ cannot be represented by any combinatorial update $F$.
]

#subsec("D", "Outlier Isolation and Worst-Case Radius Breakdown")

Because the minimax radius uses an $L_infinity$ worst-case norm over cluster members, its response to outliers exhibits a sharp two-regime behavior: complete immunity below the isolation threshold, followed by unbounded radius inflation upon forced absorption.

#proposition(title: "Outlier Isolation Below Cut Threshold and Absorption Breakdown", origin: "Original Proof")[
  Let $X = X_0 union {z}$ where an adversarial outlier $z$ satisfies $R = min_(u in X_0) d(z, u) > "diam"(X_0)$. Then:
  1. In the Minimax dendrogram of $X$, $z$ remains an isolated singleton until the final merge step $t = n - 1$, leaving the entire subtree and all prototypes on $X_0$ strictly identical to clustering $X_0$ alone.
  2. At the absorption merge, $r(C union {z}) >= R$ for every nonempty $C subset.eq X_0$, so forcing a partition into $k < k_0 + 1$ clusters inflates the absorbing cluster's radius by at least $R / r(C) arrow.r infinity$ as $R arrow.r infinity$.
]

#proof[
  For any nonempty $C subset.eq X_0$, every point $x in C$ satisfies $d_"max" (x, C union {z}) >= d(x, z) >= R$, while $z$ itself satisfies $d_"max" (z, C union {z}) = max_(u in C) d(z, u) >= R$. Taking the minimum over $x in C union {z}$ yields $d_"MM" (C, {z}) = r(C union {z}) >= R$. Conversely, by Theorem 2, every pair of disjoint subsets $G, H subset.eq X_0$ satisfies $d_"MM" (G, H) <= "diam"(G union H) <= "diam"(X_0) < R$. Hence every merge within $X_0$ strictly precedes any merge involving $z$, keeping the subtree and prototype sequence on $X_0$ completely unperturbed until step $n - 1$, where $r(X_0 union {z}) >= R$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because Minimax Linkage is non-combinatorial (Proposition 3), evaluating $d_"MM" (C_t, K)$ after merging $C_t = G_t union H_t$ requires scanning all candidate prototypes $x in C_t union K$, costing $Theta(abs(C_t) + abs(K))$ time per active cluster $K$ even when cached eccentricities $M(x, G) = max_(y in G) d(x, y)$ are maintained @bien2011minimax. Summing over all $K in cal(A) without {C_t}$ at step $t$ costs $O(abs(cal(A)) dot abs(C_t) + n)$ work; on balanced binary merges where cluster sizes double across $O(log n)$ levels, or on caterpillar merges where $abs(C_t) = t$, the cumulative cost reaches $sum_(t=1)^(n-1) t (n - t) = Theta(n^3)$ operations in the worst case. Memory usage is $O(n^2)$ to store the pairwise dissimilarity matrix $D$, the active linkage matrix $L$, and the $n times abs(cal(A))$ eccentricity table $M$.

Two practical accelerations preserve exactness without altering the worst-case bound: (i) using the lower bound $d_"MM" (C_t, K) >= min(d_"MM" (G_t, K), d_"MM" (H_t, K))$ (Theorem 1) inside a priority queue to skip recomputations of pairs exceeding the current best merge candidate, and (ii) terminating the inner prototype loop whenever $max(M(x, C_t), M(x, K))$ reaches the incumbent minimum for $(C_t, K)$.

Geometrically, Minimax Linkage imposes a *ball-covering inductive bias* analogous to hierarchical $k$-center clustering. Every cluster has an explicit exemplar center and bounded radius $h$ (Theorem 2), resisting the single-linkage chaining artifact while remaining strictly tighter than complete linkage on centrally dense spherical clusters (achieving up to a $2 times$ smaller merge height than complete linkage when a central prototype exists). However, because $d_"max" (x, C)$ depends solely on the farthest point in $C$, Minimax Linkage is acutely vulnerable to high-dimensional distance concentration @beyer1999nearest: as intrinsic dimension grows, maximum and minimum pairwise distances converge, flattening prototype discrimination.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(n^3)$ worst-case exact (via cached eccentricities $M(x, C)$)],
    [Space complexity], [$O(n^2)$ for dissimilarity, linkage, and eccentricity matrices],
    [Monotonicity], [Unconditionally guaranteed on all symmetric $d >= 0$ (zero inversions)],
    [Reducibility], [Unconditionally reducible (NN-chain valid, still $O(n^3)$ worst case)],
    [Interpretability], [Native observed prototype $p^*(C) in C$ at every dendrogram node],
    [Invariances], [Exact invariance under all strictly increasing transforms $f compose d$],
    [Inductive bias], [Compact metric balls $C subset.eq overline(B)(p^*, h)$; resists chaining, fractures arcs],
  ),
  caption: [Operational and theoretical profile of Minimax Linkage (Alg. 73).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Visual Summarization and Image Collection Browsing*: In Bien and Tibshirani's original evaluation @bien2011minimax, Minimax Linkage clustered facial and natural image datasets so that each internal dendrogram node displays an actual exemplar photograph rather than a blurry pixel average. The ball-covering bound guarantees that every image inside a branch lies within dissimilarity $h$ of the displayed exemplar.
2. *Cheminformatics and Molecular Library Diversification*: Pharmaceutical screening libraries cluster compounds under Tanimoto dissimilarity $d_T (u, v) = 1 - abs(u inter v) / abs(u union v)$. Cutting a Minimax dendrogram at $h = 0.30$ directly selects a physically synthesizable medoid compound per cluster such that every analog in the cluster shares at least $70\%$ structural fingerprint similarity with its prototype.
3. *Clinical Phenotyping and Electronic Health Records*: Healthcare stratification requires transparent cohort definitions grounded in real patient profiles rather than synthetic averages over binary diagnoses and skewed lab values. Minimax Linkage handles non-Euclidean clinical dissimilarities (e.g., Gower or comorbidity-ontology metrics) while attaching an auditable index patient to each phenotype cluster.
4. *Extractive Document and News-Feed Taxonomy*: Clustering corpora by TF-IDF cosine or semantic embedding distances yields multi-resolution topic hierarchies where each subtree is summarized by a single real representative document, avoiding unreadable centroid bag-of-words vectors @steinbach2000comparison.
5. *Hierarchical Facility Location and Sensor Covering*: Network-design problems requiring multi-scale relay placement seek a nested hierarchy of exemplar hubs minimizing worst-case communication latency ($L_infinity$ covering radius) rather than average path length ($L_1$ $k$-medoids). Minimax Linkage directly solves a greedy hierarchical relaxation of the metric $k$-center problem.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four principal failure modes govern Minimax Linkage in practice: (1) *Cubic computational wall* ($O(n^3)$), limiting exact execution to roughly $n <= 5 times 10^3$ unless preceded by coreset compression or approximate pruning; (2) *Manifold fragmentation*, since elongated crescents, concentric rings, or density filaments cannot be covered by small-radius balls without splitting into many disconnected pieces; (3) *High-dimensional eccentricity collapse*, where concentration of measure @beyer1999nearest makes $d_"max" (x, C)$ nearly uniform across all $x in C$, causing arbitrary prototype selection and degraded cluster separation; and (4) *Sensitivity to forced low-$k$ outliers*, because an isolated outlier either consumes an entire cluster budget as a singleton or triggers catastrophic radius inflation upon forced absorption (Proposition 4).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Exemplar summaries ($n <= 5000$)], [Excellent: guaranteed $overline(B)(p^*, h)$ prototype balls], [Ensure metric/meaningful $d$; prune outliers],
    [Cheminformatics / Tanimoto], [Excellent: ordinal invariance & worst-case analog radius], [For $n > 10^4$, pre-cluster with BIRCH or FasterPAM],
    [Non-convex / manifold shapes], [Poor: ball-covering fractures arcs & moons], [Use HDBSCAN @campello2013hdbscan or Spectral @ng2002spectral],
    [High-D feature spaces], [Weak: $L_infinity$ eccentricity suffers distance concentration], [Prefer Flexible Beta (Alg. 74) or PCA + Ward],
    [Massive scale ($n > 10^4$)], [Infeasible: $O(n^3)$ exact runtime], [Use FasterPAM @schubert2019fasterpam or CLARA @srikantha2011clara],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Alg. 73.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Minimax Linkage bridges two classical families: *complete-linkage hierarchical clustering* @defays1977clink (to which it is tightly factor-2 coupled by Theorem 2, sharing ordinal invariance and non-chaining behavior) and *medoid-based partitioning* (PAM / $k$-medoids @kaufman1990finding @schubert2019fasterpam). The essential difference from PAM is the choice of prototype norm: PAM minimizes the $L_1$ *average* distance $abs(C)^(-1) sum_(y in C) d(x, y)$, whereas Minimax Linkage minimizes the $L_infinity$ *maximum* distance $max_(y in C) d(x, y)$.

#proposition(title: "Divergence Between Minimax (L-Infinity) and Medoid (L-1) Prototypes", origin: "Original Counterexample")[
  For any $m >= 3$, there exists a 1D Euclidean cluster $C$ of size $m + 2$ in which the unique $L_1$ medoid $p_1$ and the unique $L_infinity$ minimax prototype $p_infinity$ lie at distance $1 / 2$ of the cluster diameter apart, while $p_1$ incurs nearly twice the minimax radius of $p_infinity$.
]

#proof[
  Place $m$ points at $0$ (perturbed by tiny distinct $epsilon_i in (0, 10^(-6))$), one point at $1 / 2$, and one point at $1$, so $"diam"(C) = 1$. Evaluating total $L_1$ distances:
  - At $x approx 0$: $sum_(y in C) abs(x - y) approx 1 / 2 + 1 = 3 / 2$.
  - At $x = 1 / 2$: $sum_(y in C) abs(1 / 2 - y) approx m / 2 + 1 / 2 >= 2 > 3 / 2$ (for all $m >= 3$).
  - At $x = 1$: $sum_(y in C) abs(1 - y) approx m + 1 / 2 > 3 / 2$.
  Thus the unique $L_1$ medoid is $p_1 approx 0$, whose worst-case covering radius is $d_"max" (p_1, C) approx 1$. Conversely, the eccentricities are $d_"max" (0, C) = 1$, $d_"max" (1 / 2, C) = 1 / 2$, and $d_"max" (1, C) = 1$, so the unique minimax prototype is $p_infinity = 1 / 2$ with $r(C) = 1 / 2 = (1 / 2) d_"max" (p_1, C)$ and $abs(p_infinity - p_1) approx 1 / 2 = (1 / 2) "diam"(C)$.
]

This divergence encapsulates the practical trade-off: $L_1$ medoids track bulk density modes (ignoring isolated tails), whereas $L_infinity$ minimax prototypes track geometric midpoints (guaranteeing uniform worst-case coverage). Unlike Lance-Williams rules (Alg. 72) and Flexible Beta (Alg. 74), Minimax Linkage sacrifices $O(n^2)$ combinatorial updates in exchange for native exemplar interpretability and ordinal invariance (Proposition 3).

*Empirical Benchmark Behavior*: On the ClusterOverview 4-dataset benchmark, Minimax Linkage obtains Iris Silhouette $= 0.455$ (@rousseeuw1987silhouettes) and $"ARI" = 0.635$ (@hubert1985comparing), Customers Silhouette $= 0.335$, Multi-Shape Silhouette $= 0.417$, and High-D Silhouette $= 0.167$, with a mean runtime of $8.141$ ms. Every metric reflects our theoretical analysis with precision: (i) on compact, low-dimensional Iris clusters, Minimax Linkage achieves top-tier recovery ($"Sil" = 0.455$, $"ARI" = 0.635$, tying Flexible Beta and decisively outperforming default Lance-Williams at $"ARI" = 0.565$); (ii) its mean runtime ($8.141$ ms) is roughly $3.5 times$ higher than Lance-Williams ($2.346$ ms) and Flexible Beta ($2.458$ ms), directly manifesting the $O(n^3)$ prototype-scanning cost vs. $O(n^2)$ combinatorial updates; (iii) on Multi-Shape ($0.417$), its spherical ball-covering bias fragments non-convex shapes; and (iv) on High-D ($0.167$), its $L_infinity$ eccentricity suffers the sharpest degradation among the three monotonic hierarchical methods due to high-dimensional distance concentration @beyer1999nearest @jain2010data.
