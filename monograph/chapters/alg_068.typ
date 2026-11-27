#import "../template/ieee.typ": *

#alg_entry(
  id: 68,
  name: "MONA: Monothetic Analysis for Binary Data",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1990",
  paradigm: "Monothetic Divisive Binary Hierarchical Clustering",
  complexity: $O(n d^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Most hierarchical clustering methods are _polythetic_: splits or merges
depend simultaneously on all $d$ attributes through an aggregate
dissimilarity matrix @kaufman1990finding @xu2005survey. In contrast,
*Monothetic Analysis* (*MONA*), formalized by Kaufman and Rousseeuw
@kaufman1990finding as a robust refinement of Williams and Lambert's
(1959) association analysis, constructs a top-down divisive hierarchy
over binary observations $X in {0, 1}^(n times d)$ by splitting each
active cluster along a *single* binary coordinate at every step. Every
resulting cluster is an exact logical conjunction of attribute literals,
yielding direct rule-based interpretability without post-hoc surrogate
models.

Let $I = {1, dots, n}$ index the objects and $F = {1, dots, d}$ index
the binary attributes, where $x_(i f) in {0, 1}$ denotes the value of
attribute $f$ on object $i$. For any active cluster $C subset.eq I$ of
cardinality $m = abs(C) >= 2$ and any pair of attributes $f, g in F$,
define the $2 times 2$ contingency table counts on $C$:
$
  a_(f g) (C) &= |{i in C : x_(i f) = 1, x_(i g) = 1}|, \
  b_(f g) (C) &= |{i in C : x_(i f) = 1, x_(i g) = 0}|, \
  c_(f g) (C) &= |{i in C : x_(i f) = 0, x_(i g) = 1}|, \
  d_(f g) (C) &= |{i in C : x_(i f) = 0, x_(i g) = 0}|,
$
satisfying $a_(f g) (C) + b_(f g) (C) + c_(f g) (C) + d_(f g) (C) = m$.
Let $s_f (C) = a_(f g) (C) + b_(f g) (C) = sum_(i in C) x_(i f)$ denote
the column sum of attribute $f$ inside $C$, and let $p_f (C) = s_f (C)
/ m$ denote its empirical prevalence.

#definition(title: "MONA Total Association Criterion")[
  For a cluster $C subset.eq I$ and distinct binary attributes $f, g in
  F$, the *pairwise determinant association* is
  $
    alpha_C (f, g) = |a_(f g) (C) d_(f g) (C) - b_(f g) (C) c_(f g) (C)|.
  $
  The *total association* of attribute $f$ within cluster $C$ is
  $
    A_C (f) = sum_(g in F without {f}) alpha_C (f, g).
  $
  The monothetic splitting variable for $C$ is selected as
  $f^* (C) in arg max_(f in F) A_C (f)$, partitioning $C$ into
  $C_1 = {i in C : x_(i, f^*) = 1}$ and
  $C_0 = {i in C : x_(i, f^*) = 0}$ whenever $A_C (f^*) > 0$.
]

Classical association analysis maximizes the sum of Pearson $chi^2$
statistics $sum_(g != f) m (a_(f g) d_(f g) - b_(f g) c_(f g))^2 /
(s_f (m - s_f) s_g (m - s_g))$, which suffers division-by-zero
singularities whenever an attribute becomes constant within a subcluster
and strongly favors extreme minority splits. Replacing $chi^2$ by the
unnormalized $2 times 2$ determinant magnitude $alpha_C (f, g)$
simultaneously eliminates division by zero and introduces a built-in
variance-weighted damping against outlier splits, as proven in
Section II.A.

Before recursive splitting begins, MONA handles missing entries in
$X in {0, 1, "NA"}^(n times d)$ via a single-pass deterministic proxy
imputation @kaufman1990finding. For each attribute $f$ possessing
missing entries, MONA computes $alpha_I (f, g)$ on the jointly observed
subset $I_(f g) = {i : x_(i f) != "NA", x_(i g) != "NA"}$, identifies
the maximally associated proxy $g^* = arg max_(g != f) alpha_I (f, g)$,
and sets $x_(i f) <- x_(i g^*)$ if $a_(f g^*) d_(f g^*) > b_(f g^*)
c_(f g^*)$ or $x_(i f) <- 1 - x_(i g^*)$ otherwise.

#algbox(
  title: "MONA (Monothetic Analysis)",
  inputs: [Binary matrix $X in {0, 1, "NA"}^(n times d)$; optional depth
    cap $h_"max"$],
  outputs: [Divisive binary tree $cal(T)$ with splitting variables and
    separation-step matrix $S in bb(N)^(n times n)$],
)[
  + *Imputation:* For each $f in F$ with missing entries, find
    $g^* = arg max_(g != f) alpha_(I_(f g)) (f, g)$ and impute missing
    $x_(i f)$ from $x_(i g^*)$ using the sign of
    $a_(f g^*) d_(f g^*) - b_(f g^*) c_(f g^*)$.
  + *Initialization:* Set active frontier $cal(Q) <- {(I, 1)}$ and
    record level $ell = 1$.
  + *While* $cal(Q)$ contains a cluster $(C, ell)$ with $abs(C) >= 2$
    and $ell <= h_"max"$:
    - Compute the $d times d$ co-occurrence Gram matrix
      $G_C = X_C^top X_C$ and column marginals $s_C = "diag"(G_C)$.
    - For all $f != g$, form $a_(f g) = (G_C)_(f g)$,
      $b_(f g) = s_f - a_(f g)$, $c_(f g) = s_g - a_(f g)$, and
      $d_(f g) = abs(C) - s_f - s_g + a_(f g)$.
    - Evaluate $A_C (f) = sum_(g != f) |a_(f g) d_(f g) - b_(f g)
      c_(f g)|$ and select $f^* = arg max_(f in F) A_C (f)$.
    - *If* $A_C (f^*) = 0$ (all active pairs uncorrelated or at most one
      non-constant attribute remains): mark $C$ as a terminal leaf.
    - *Else:* split $C$ into $C_1 = {i in C : x_(i, f^*) = 1}$ and
      $C_0 = C without C_1$; record separation step $S_(i j) = ell$ for
      all $(i, j) in C_1 times C_0$; push $(C_1, ell + 1)$ and
      $(C_0, ell + 1)$ onto $cal(Q)$.
  + *Return* the binary decision hierarchy $cal(T)$ and matrix $S$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four formal guarantees for MONA: an exact empirical
covariance representation showing why MONA suppresses outlier splits
relative to $chi^2$ association analysis, a geometric characterization
of subcube termination and bit-flip invariance, a constructive parity
counterexample exposing second-order blindness, and a classical analysis
of proxy imputation together with an accelerated Gram-subtraction
complexity theorem.

#subsec("A", "Covariance Equivalence and Outlier-Split Damping")

#lemma(
  title: "Exact Covariance Identity and Outlier Damping",
  origin: "Original Derivation",
)[
  Let $C subset.eq I$ have size $m = abs(C) >= 2$, empirical attribute
  means $p_f = s_f (C) / m$, and empirical covariance matrix
  $hat(Sigma)_C in bb(R)^(d times d)$ with entries
  $(hat(Sigma)_C)_(f g) = m^(-1) sum_(i in C) (x_(i f) - p_f)
  (x_(i g) - p_g)$. Then for every $f, g in F$:
  $
    alpha_C (f, g) = m^2 |"Cov"_C (X_f, X_g)|,
  $
  and consequently
  $
    A_C (f) = m^2 (norm(hat(Sigma)_C e_f)_1 - p_f (1 - p_f)).
  $
  Furthermore, $alpha_C (f, g) <= m^2 sqrt(p_f (1 - p_f) p_g (1 - p_g))
  <= (m^2 / 2) sqrt(p_f (1 - p_f))$, whereas any singleton-support
  attribute ($s_f = 1$) satisfies $alpha_C (f, g) <= m - 1$, yielding a
  factor-$Theta(m)$ suppression relative to balanced correlated splits
  where Pearson's $chi^2$ assigns identical maximal score $m$.
]

#proof[
  Using $d_(f g) = m - a_(f g) - b_(f g) - c_(f g)$ together with
  $s_f = a_(f g) + b_(f g)$ and $s_g = a_(f g) + c_(f g)$, direct
  algebraic expansion gives
  $
    a_(f g) d_(f g) - b_(f g) c_(f g)
      &= a_(f g) (m - a_(f g) - b_(f g) - c_(f g)) - b_(f g) c_(f g) \
      &= m a_(f g) - (a_(f g) + b_(f g)) (a_(f g) + c_(f g)) \
      &= m^2 (a_(f g) / m - p_f p_g).
  $
  Because $x_(i f), x_(i g) in {0, 1}$, we have
  $m^(-1) sum_(i in C) x_(i f) x_(i g) = a_(f g) / m$, so
  $a_(f g) / m - p_f p_g = "Cov"_C (X_f, X_g)$. Taking absolute values
  proves $alpha_C (f, g) = m^2 |"Cov"_C (X_f, X_g)|$. Summing over all
  $g != f$ and noting $("Cov"_C (X_f, X_f)) = p_f (1 - p_f) >= 0$
  yields $A_C (f) = m^2 (norm(hat(Sigma)_C e_f)_1 - p_f (1 - p_f))$.

  For the variance bound, apply Cauchy--Schwarz to the centered vectors
  $u_f = (x_(i f) - p_f)_(i in C)$ and $u_g = (x_(i g) - p_g)_(i in C)$,
  whose squared norms equal $m p_f (1 - p_f)$ and $m p_g (1 - p_g)$:
  $
    alpha_C (f, g) = m |chevron.l u_f, u_g chevron.r|
      <= m^2 sqrt(p_f (1 - p_f) p_g (1 - p_g)).
  $
  Since $p_g (1 - p_g) <= 1 / 4$, this is at most
  $(m^2 / 2) sqrt(p_f (1 - p_f))$, with global maximum $m^2 / 4$
  attained when $p_f = p_g = 1 / 2$ and $X_f = X_g$ (where
  $a_(f g) = d_(f g) = m / 2$, $b_(f g) = c_(f g) = 0$, giving
  $chi^2 = m$). Conversely, if $s_f = 1$, then $a_(f g) in {0, 1}$:
  if $a_(f g) = 1$, then $a_(f g) d_(f g) - b_(f g) c_(f g) = m - s_g
  <= m - 1$; if $a_(f g) = 0$, then $|a_(f g) d_(f g) - b_(f g)
  c_(f g)| = s_g <= m - 1$. In the tied rare-outlier case $s_f = s_g =
  1$ with $a_(f g) = 1$, Pearson's $chi^2 = m (m - 1)^2 / ((m - 1)^2)
  = m$ ties the perfectly balanced split, whereas MONA scores $m - 1$
  versus $m^2 / 4$, an asymptotic damping ratio of $m / 4$.
]

#subsec("B", "Subcube Termination, Ultrametricity, and Symmetry")

#theorem(
  title: "Finite Subcube Termination, Ultrametricity, and Bit-Flip Invariance",
  origin: "Original Proof",
)[
  Let $X in {0, 1}^(n times d)$ be complete. Then:
  + Every attribute $f in F$ is selected as a splitting variable at most
    once along any root-to-leaf path of $cal(T)$; consequently, MONA
    terminates in depth $h <= min(d - 1, n - 1)$ with at most
    $min(2^(d - 1), n)$ leaves, and each node at depth $ell$ is an
    exact Boolean subcube of dimension $d - ell$.
  + Defining $u (i, j) = 2^(-S_(i j))$ when $i, j$ are separated at step
    $S_(i j) >= 1$ and $u (i, j) = 0$ when they remain in the same
    terminal leaf yields an exact ultrametric on the leaf-quotient
    space.
  + MONA's partition hierarchy is strictly invariant to object
    permutations, attribute permutations, and arbitrary attribute-wise
    bit-complementations $x_(i f) mapsto 1 - x_(i f)$, but is *not*
    invariant to attribute duplication.
]

#proof[
  (1) Suppose node $C$ splits on $f^*$. Then $x_(i, f^*) = 1$ for all
  $i in C_1$ and $x_(i, f^*) = 0$ for all $i in C_0$. For any descendant
  $D subset.eq C_1$, $p_(f^*) (D) = 1$, so by Lemma 1,
  $alpha_D (f^*, g) <= abs(D)^2 sqrt(1 dot 0) = 0$ for all $g in F$,
  whence $A_D (f^*) = 0$; symmetrically $A_D (f^*) = 0$ for all
  $D subset.eq C_0$. Thus $f^*$ can never attain positive total
  association below $C$. Moreover, if only one non-constant variable $f$
  remains in $D$, every $g != f$ is constant so $A_D (f) = 0$ and $D$
  terminates. Hence along every root-to-leaf path at most $d - 1$ splits
  occur; since each split yields two nonempty children, depth is also
  bounded by $n - 1$, giving $h <= min(d - 1, n - 1)$ and at most
  $min(2^(d - 1), n)$ leaves. Intersecting the $ell$ axis-aligned
  half-spaces $x_(i, f_r) = v_r$ ($r = 1, dots, ell$) identifies each
  depth-$ell$ cluster with a $(d - ell)$-dimensional subcube of
  ${0, 1}^d$.

  (2) Symmetry $u (i, j) = u (j, i) >= 0$ is immediate. For any three
  objects $i, j, k$, let $C$ be the deepest cluster containing all
  three, at level $ell$. If $C$ is a leaf, $u (i, j) = u (j, k) =
  u (i, k) = 0$. Otherwise $C$ partitions into $(C_1, C_0)$, so by the
  pigeonhole principle at least two objects—say $i, j$—fall in the same
  child and the third $k$ in the opposite child (or all three enter the
  same child, contradicting maximality of $ell$). Thus
  $S_(i k) = S_(j k) = ell$ while $S_(i j) > ell$ (or $i, j$ never
  separate), giving $u (i, j) <= 2^(-(ell + 1)) < 2^(-ell) =
  max(u (i, k), u (j, k))$. Hence $u$ satisfies the strong triangle
  inequality.

  (3) Object and attribute permutations merely relabel rows and columns
  of $G_C$, leaving $A_C (f)$ invariant. Complementing attribute $f$ via
  $tilde(x)_(i f) = 1 - x_(i f)$ swaps $(a_(f g), b_(f g), c_(f g),
  d_(f g)) mapsto (c_(f g), d_(f g), a_(f g), b_(f g))$ for every
  $g != f$, whence
  $
    |tilde(a)_(f g) tilde(d)_(f g) - tilde(b)_(f g) tilde(c)_(f g)|
      = |c_(f g) b_(f g) - d_(f g) a_(f g)| = alpha_C (f, g).
  $
  Thus every pairwise score and total association $A_C (f)$ is strictly
  preserved. Finally, duplicating attribute $f$ into a copy $f'$ adds
  $alpha_C (f, f') = m^2 p_f (1 - p_f) > 0$ to $A_C (f)$ and $A_C (f')$
  plus an extra $alpha_C (f, g)$ to every other $g$, so duplicating a
  weakly correlated balanced variable enough times forces MONA to split
  on $f$, disproving duplication invariance.
]

#subsec("C", "Higher-Order Parity Blindness Counterexample")

#theorem(
  title: "Exact Pairwise-Independence Stall on Parity Structures",
  origin: "Original Counterexample",
)[
  For every integer $d >= 3$ and every cluster multiplicity $r >= 1$,
  there exists a dataset $X in {0, 1}^(n times d)$ with $n = 2^(d - 1)
  r$ objects partitioned into $2^(d - 1)$ perfectly pure, well-separated
  clusters of equal size $r$ and minimum inter-cluster Hamming distance
  $d_H >= 2$ such that $A_I (f) = 0$ for all $f in F$, causing MONA to
  terminate immediately at the root with $K = 1$ cluster.
]

#proof[
  Let $E_d = {v in {0, 1}^d : sum_(f = 1)^d v_f equiv 0 med (mod 2)}$ be
  the even-parity Boolean code of length $d >= 3$, which has cardinality
  $|E_d| = 2^(d - 1)$. Construct $X$ by placing exactly $r$ identical
  objects at each codeword $v in E_d$, so $n = 2^(d - 1) r$. Distinct
  codewords in $E_d$ differ in an even positive number of coordinates,
  so their pairwise Hamming distance satisfies $d_H (v, v') >= 2$,
  whereas every intra-cluster distance is $0$.

  Fix any two distinct coordinates $f, g in {1, dots, d}$ and any bit
  pair $(u_1, u_2) in {0, 1}^2$. Since $d >= 3$, choose a third
  coordinate $h in F without {f, g}$. For every assignment to the
  remaining $d - 3$ coordinates, there exists a unique value of
  $v_h in {0, 1}$ making the total sum even. Hence exactly $2^(d - 3)$
  codewords in $E_d$ satisfy $(v_f, v_g) = (u_1, u_2)$, giving
  $
    a_(f g) (I) = b_(f g) (I) = c_(f g) (I) = d_(f g) (I) = 2^(d - 3) r.
  $
  Consequently,
  $
    alpha_I (f, g) = |(2^(d - 3) r)^2 - (2^(d - 3) r)^2| = 0
    quad forall f != g.
  $
  Therefore $A_I (f) = 0$ for every $f in F$, and MONA halts at level
  $ell = 1$ without executing a single split, despite the joint mutual
  information $I (X_d ; X_1, dots, X_(d - 1)) = 1$ bit being maximal.
]

#subsec("D", "Proxy Imputation Mechanics and Gram-Subtraction Acceleration")

#theorem(
  title: "Conditional-Lift Imputation and Accelerated Gram Subtraction",
  origin: "Classical — Cited",
)[
  Under Kaufman and Rousseeuw's missing-value imputation rule
  @kaufman1990finding:
  + Positive association $a_(f g) d_(f g) > b_(f g) c_(f g)$ on
    $I_(f g)$ is equivalent to strictly positive conditional lift
    $hat("Pr")(X_f = 1 | X_g = 1) > hat("Pr")(X_f = 1 | X_g = 0)$, and
    imputing $k_1$ rows with $x_(i g^*) = 1$ and $k_0$ rows with
    $x_(i g^*) = 0$ strictly preserves the association sign while
    increasing $alpha (f, g^*)$ by $k_1 d_(f g^*) + k_0 a_(f g^*) +
    k_1 k_0 >= 0$. However, positive association does *not* guarantee
    that $x_(i f) = x_(i g^*)$ maximizes agreement count over
    $1 - x_(i g^*)$.
  + By computing $G_(C_"small") = X_(C_"small")^top X_(C_"small")$ only
    on the smaller child ($abs(C_"small") <= abs(C) / 2$) and forming
    $G_(C_"large") = G_C - G_(C_"small")$, the full MONA hierarchy is
    computed identically in $O(n d^2 min(d, log n))$ worst-case time,
    strictly improving the naive $O(n d^3)$ bound @kaufman1990finding.
]

#proof[
  (1) Dividing $a_(f g) d_(f g) - b_(f g) c_(f g) > 0$ by
  $(a_(f g) + c_(f g)) (b_(f g) + d_(f g)) > 0$ yields
  $
    a_(f g) / (a_(f g) + c_(f g)) - b_(f g) / (b_(f g) + d_(f g)) > 0,
  $
  which is precisely $hat("Pr")(X_f = 1 | X_g = 1) > hat("Pr")(X_f = 1
  | X_g = 0)$. Imputing $x_(i f) = x_(i g^*)$ adds $k_1$ to $a_(f g^*)$
  and $k_0$ to $d_(f g^*)$, changing the determinant to
  $(a_(f g^*) + k_1)(d_(f g^*) + k_0) - b_(f g^*) c_(f g^*) =
  (a_(f g^*) d_(f g^*) - b_(f g^*) c_(f g^*)) + k_1 d_(f g^*) +
  k_0 a_(f g^*) + k_1 k_0$, preserving sign and inflating magnitude.
  Conversely, the table $(a, b, c, d) = (2, 1, 15, 10)$ has
  $a d - b c = 20 - 15 = 5 > 0$, yet agreement for $x_f = x_g$ is
  $a + d = 12 < 16 = b + c$, proving that conditional lift does not
  imply majority agreement under asymmetric marginals.

  (2) Let $M$ be the number of internal splits in $cal(T)$; every
  internal split separates a cluster of size $>= 2$ into two nonempty
  children, so $M <= n - 1$. At each internal node $C$, evaluating
  $A_C (f)$ from $G_C$ takes $O(d^2)$ time, totaling $O(M d^2) subset.eq
  O(n d^2)$, and partitioning $C$ takes $O(abs(C))$, which sums to
  $O(n h) subset.eq O(n d)$ across the $h <= d - 1$ levels. Updating
  $G_(C_"small")$ costs $O(abs(C_"small") d^2)$ and subtracting
  $G_(C_"large") = G_C - G_(C_"small")$ costs $O(d^2)$ per split. Each
  object $i in I$ enters a smaller child at most $floor(log_2 n)$ times
  (because cluster cardinality halves upon every entry into a smaller
  child) and at most $h <= d - 1$ times (by Theorem 1). Summing over
  all internal nodes gives
  $
    sum_(C "internal") abs(C_"small") d^2
      <= n d^2 min(d - 1, floor(log_2 n)),
  $
  yielding total exact worst-case runtime $O(n d^2 min(d, log n))$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Classical MONA recomputes all $binom(d, 2)$ contingency tables from
scratch at every node, incurring $O(abs(C) d^2)$ per cluster and
$O(n d^2)$ per tree level, or $O(n d^3)$ across a worst-case path of
depth $d - 1$ @kaufman1990finding. With the smaller-half Gram
subtraction of Theorem 3, worst-case runtime drops to
$O(n d^2 min(d, log n))$, and bit-packed $w$-bit SIMD popcount kernels
(`x_f & x_g`) further divide the Gram-formation constant by $w = 64$.
Space complexity is $O(n d + d^2 h)$ if storing Gram matrices along the
active DFS stack, or $O(n d + d^2)$ by recomputing at most one stack
frame.

Geometrically, MONA partitions the Boolean hypercube ${0, 1}^d$ into
axis-aligned orthogonal subcubes (equivalently, an unsupervised decision
tree). Unlike polythetic divisive algorithms such as DIANA
@kaufman1990finding—which compute an $O(n^2)$ dissimilarity matrix and
split along oblique boundaries—MONA never materializes an $n times n$
matrix, scaling linearly in $n$ (up to $log n$) and excelling when
$n >> d$. Its operational bottleneck is quadratic scaling in $d$: when
$d$ exceeds several thousand binary indicators, evaluating $A_C (f)$
dominates runtime and sparse inverted-list intersections are required.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Binary matrix ${0, 1}^(n times d)$ (with NA support)],
    [Time (Classical / Accel.)], [$O(n d^3)$ worst-case / $O(n d^2 min(d, log n))$],
    [Space Complexity], [$O(n d + d^2)$ (DFS with single active Gram frame)],
    [Cluster Geometry], [Axis-aligned Boolean subcubes (monothetic rules)],
    [Invariances], [Row/column permutations, bit flips $x mapsto 1 - x$],
    [Non-Invariances], [Attribute duplication, monotone real transforms],
    [Tree Height / Leaves], [$h <= min(d - 1, n - 1)$; leaves $<= min(2^(d-1), n)$],
    [Determinism], [Fully deterministic (up to exact $A_C (f)$ ties)],
  ),
  caption: [Operational and theoretical profile of MONA.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

MONA is uniquely suited to domains where native observations are binary
presence/absence vectors and where stakeholders require exact,
auditable, single-attribute decision rules:

+ *Ecological Floristics and Species Assemblage Classification:*
  Quadrat survey matrices record binary presence/absence of $d$ plant
  or marine species across $n$ sites. Indicator species strongly
  co-occur with environmental guilds, producing high pairwise covariance
  $|hat("Cov")(X_f, X_g)|$ and allowing MONA to extract hierarchical
  taxonomic keys (each split naming a single diagnostic indicator
  species) without arbitrary distance metric choices @kaufman1990finding.
+ *Clinical Phenotyping and Comorbidity Stratification:* Electronic
  health records encode ICD diagnostic flags and binary symptom panels.
  Because medical protocols demand transparent stratification trees,
  MONA's subcube leaves map directly to conjunctive clinical criteria
  while Lemma 1's outlier-damping property prevents ultra-rare syndromes
  ($s_f = 1$) from hijacking top-level splits.
+ *Bioinformatics Presence--Absence Pan-Genome Analysis:* Bacterial
  pan-genomes represent strains by binary accessory-gene presence
  profiles. Co-inherited operons and mobile genetic elements induce
  block-structured pairwise covariances; MONA isolates functional
  pathotypes along single marker genes and tolerates bit-polarity
  conventions via Theorem 1(3).
+ *Market-Basket and Software Feature-Flag Telemetry:* Binary
  subscription bundles or feature-flag activation logs exhibit dense
  co-adoption blocks. Since MONA scales as $O(n d^2 log n)$ without
  forming an $n times n$ Jaccard matrix, it segments millions of user
  profiles into human-readable product-tier trees directly usable by
  downstream rules engines.
+ *Psychometric Questionnaire and Survey Typology:* Yes/no item banks
  measuring latent attitudes @lazarsfeld1968latent often contain
  item-level missingness and reverse-coded questions. MONA's bit-flip
  invariance eliminates manual polarity alignment, while its native
  proxy imputation handles scattered non-responses without listwise
  deletion @kaufman1990finding.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes govern MONA's practical envelope. First,
*binarization information loss*: when applied to continuous data via
median or quantile thresholding, MONA discards metric magnitude and
intra-quadrant geometry, forcing axis-aligned cuts through oblique or
curved manifolds. Second, *higher-order parity blindness* (Theorem 2):
whenever cluster membership depends on XOR-like or three-way
interactions with balanced two-way marginals, $A_C (f) = 0$ stalls the
split. Third, *irrevocable divisive errors*: a suboptimal early split on
a noisy indicator cannot be undone at deeper levels. Fourth,
*redundancy sensitivity* (Theorem 1(3)): blocks of duplicated or
near-collinear binary flags inflate $A_C (f)$ linearly with copy count,
biasing top splits toward over-represented feature groups.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Binary ecology & taxonomy], [Exact indicator-key subcubes; strong co-occurrence], [Avoid collinear clones; check $d >> 10^3$],
    [Clinical symptom panels], [Auditable rule leaves; robust to rare flags], [Validate missingness; Latent Class @lazarsfeld1968latent],
    [Continuous Euclidean data], [Poor: thresholding destroys metric geometry], [Use Ward, Genie, or BUBBLE],
    [Sparse high-$d$ baskets], [Good if $d$ moderate; $O(d^2)$ bottleneck if huge], [Use ROCK @guha2000rock or $k$-Modes @huang1998extensions],
    [Non-linear / XOR patterns], [Fails: zero pairwise covariance stalls root], [Use spectral or density clustering],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for MONA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MONA bridges hierarchical clustering, unsupervised decision trees, and
categorical partitioning. Unlike polythetic binary/categorical methods
such as ROCK @guha2000rock (which links objects via shared-neighbor
Jaccard graphs) or $k$-Modes @huang1998extensions (which minimizes
Hamming distortion iteratively), MONA optimizes an explicit second-order
covariance norm over single coordinates, paralleling COBWEB's
@fisher1987cobweb category utility but restricting each split to one
attribute.

#proposition(
  title: "Reduction of One-Step MONA to Weighted Cut-Covariance Maximization",
  origin: "Original Proof",
)[
  Splitting cluster $C$ on binary variable $f$ partitions $C$ into
  $(C_1, C_0)$ with sizes $(m_1, m_0)$. Let $mu_1, mu_0 in [0, 1]^d$ be
  the binary centroid vectors of $C_1$ and $C_0$. Then MONA's objective
  satisfies the exact geometric identity
  $
    A_C (f) = m_1 m_0 (norm(mu_1 - mu_0)_1 - 1).
  $
]

#proof[
  For any attribute $g != f$, the within-child means are
  $(mu_1)_g = a_(f g) / m_1$ and $(mu_0)_g = c_(f g) / m_0$. Using
  $m_1 = a_(f g) + b_(f g)$ and $m_0 = c_(f g) + d_(f g)$:
  $
    m_1 m_0 |(mu_1)_g - (mu_0)_g|
      &= |a_(f g) (c_(f g) + d_(f g)) - c_(f g) (a_(f g) + b_(f g))| \
      &= |a_(f g) d_(f g) - b_(f g) c_(f g)| = alpha_C (f, g).
  $
  Summing over all $g != f$ and noting $|(mu_1)_f - (mu_0)_f| = |1 - 0|
  = 1$ yields
  $A_C (f) = m_1 m_0 (sum_(g = 1)^d |(mu_1)_g - (mu_0)_g| - 1) =
  m_1 m_0 (norm(mu_1 - mu_0)_1 - 1)$.
]

Proposition 1 gives a striking geometric interpretation: whereas Ward's
method maximizes the $L_2^2$ between-cluster variance
$(m_1 m_0 / m) norm(mu_1 - mu_0)_2^2$ over arbitrary partitions, MONA
maximizes the $L_1$ between-cluster separation $m_1 m_0 (norm(mu_1 -
mu_0)_1 - 1)$ restricted to axis-aligned binary hyperplanes—explaining
both its computational speed and its loss of fidelity on continuous
non-axis-aligned geometries @dasgupta2016cost @kleinberg2002impossibility.

This theoretical trade-off is confirmed directly on the ClusterOverview
4-dataset continuous benchmark (using silhouette @rousseeuw1987silhouettes
and adjusted Rand index @hubert1985comparing after binarization). MONA
executes at an exceptionally low mean runtime of *0.176 ms*—over an
order of magnitude faster than metric-tree and MST-based hierarchical
rivals—reflecting its $O(n d^2)$ bitwise co-occurrence kernel with no
pairwise $n times n$ distance matrix. Conversely, quantization onto
axis-aligned binary subcubes caps geometric cluster quality on
continuous benchmarks: Iris achieves $"Sil" = 0.233$ and $"ARI" =
0.354$ (oblivious to oblique sepal/petal covariance structure),
Customers reaches $"Sil" = 0.328$, Multi-Shape attains $"Sil" = 0.304$
($"ARI" = "N/A"$, as orthogonal hypercube cuts fracture curved rings and
moons), and High-D scores $"Sil" = 0.292$ ($"ARI" = "N/A"$), where
binary thresholding actually exhibits slightly gentler degradation than
Euclidean metric balls in high dimensions. Consequently, MONA should be
deployed strictly where binary semantics and monothetic rule clarity are
primary objectives, rather than as a substitute for metric hierarchical
clustering on continuous manifolds @jain2010data.
