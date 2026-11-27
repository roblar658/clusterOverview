#import "../template/ieee.typ": *

#alg_entry(
  id: 139,
  name: "DOC: Density-Based Optimal Projective Clustering",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2002",
  paradigm: "Monte Carlo Axis-Aligned Projective Clustering",
  complexity: $cal(O)(n d dot.op alpha^(-(r + 1)))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Projected clustering algorithms face an inherent combinatorial
trade-off: exhaustively enumerating axis-parallel subspaces
(CLIQUE @agrawal1998clique, SUBCLU @kailing2004subclu) scales
exponentially in cluster dimensionality, whereas medoid-based
heuristics (PROCLUS @aggarwal1999proclus) require the user to
specify a fixed average subspace dimensionality $l$ and offer no
approximation guarantee against an optimal subspace cluster.
DOC (Density-Based Optimal Projective Clustering), introduced by
Procopiuc, Jones, Agarwal, and Murali @procopiuc2002doc, bridges
this gap by formulating an explicit mathematical objective over
hypercube-bounded projective clusters and solving it with a
provably accurate Monte Carlo algorithm whose runtime is strictly
linear in $n$ and polynomial in $d$ for fixed density parameters.

Let $cal(X) = \{bold(x)_1, ..., bold(x)_n\} subset bb(R)^d$ be
a dataset over attribute set $A = \{1, ..., d\}$, and let
$w > 0$ be a fixed coordinate side width, $alpha in (0, 1)$ a
minimum cluster density fraction, and $beta in (0, 1 / 2)$ a
dimension-cardinality balance parameter (canonically $beta = 0.25$).
For any subset $C subset.eq cal(X)$ and coordinate $i in A$,
define the projected coordinate span
$"diam"_i (C) = max_(bold(q) in C) q_i - min_(bold(q) in C) q_i$.

#definition(title: "Alpha-Beta Projective Cluster and Quality")[
  A pair $(C, D)$ with $C subset.eq cal(X)$ and
  $emptyset != D subset.eq A$ is an _$(alpha, beta)$-projective
  cluster_ of width $w$ if it satisfies:
  1. _Support_: $|C| >= alpha n$;
  2. _Bounded projection_: $"diam"_i (C) <= w$ for every
     $i in D$;
  3. _$beta$-dispersion outside $D$_: for every $i in A without D$
     and every closed interval $I subset bb(R)$ of length $w$,
     $|\{bold(q) in C : q_i in I\}| <= beta |C|$.
  To balance cluster cardinality $a = |C|$ against subspace
  dimensionality $b = |D|$, DOC defines the multiplicative
  quality function @procopiuc2002doc
  $
    mu(a, b) = a dot.op (1 / beta)^b.
  $
  An _optimal projective cluster_ maximizes $mu(|C|, |D|)$ over
  all pairs $(C, D)$ satisfying Conditions 1 and 2.
]

Because $beta < 1 / 2$, adding one bounded coordinate to $D$
multiplies $mu$ by $1 / beta > 2$; hence the objective prefers a
higher-dimensional subspace whenever doing so retains more than
a fraction $beta$ of the cluster's points. Searching directly
over all $(2^d - 1)$ subspaces and all box placements is NP-hard
in general, so DOC leverages a clean randomized witness idea: if
a random seed (medoid) $bold(p)$ lands inside an optimal cluster
$C^*$, then *every* point $bold(q) in C^*$ lies inside the
$L_infinity$ half-width-$w$ box
$
  B_w (bold(p), D)
  = \{bold(q) in cal(X) : |q_i - p_i| <= w, forall i in D\}
$
for $D = D^*$. To identify $D^*$ without scanning all subspaces,
DOC draws a small _discriminating set_ $X subset cal(X)$ of size
$r$ and keeps precisely the coordinates on which every point of
$X$ stays within distance $w$ of $bold(p)$.

#algbox(
  title: "DOC Monte Carlo Optimal Projective Clustering",
  inputs: [Dataset $cal(X) subset bb(R)^d$, width $w > 0$,
    parameters $alpha in (0, 1)$, $beta in (0, 1 / 2)$, failure
    tolerance $delta in (0, 1)$, cluster count $k$.],
  outputs: [Sequence of up to $k$ projective clusters
    $\{(hat(C)_t, hat(D)_t)\}$.],
)[
  1. *Sample Sizes*: Set
     $r <- ceil(ln(2 d) / ln(1 / (2 beta)))$,
     outer trials $T_"out" <- ceil(alpha^(-1) ln(2 / delta))$,
     and inner trials $m <- ceil(2 alpha^(-r) ln(2 / delta))$
     (capped by `"MAXITER"` in FastDOC @procopiuc2002doc).
  2. *Extract One Cluster*: Initialize $mu_"best" <- 0$,
     $(hat(C), hat(D)) <- (emptyset, emptyset)$.
  3. *Outer Medoid Loop*: Repeat $T_"out"$ times:
     Sample $bold(p) in cal(X)$ uniformly at random.
  4. *Inner Discriminating Loop*: Repeat $m$ times:
     Draw $X = \{bold(q)_1, ..., bold(q)_r\}$ uniformly at random
     (with replacement) from $cal(X)$.
  5. *Subspace Inference*: Set
     $D <- \{i in A : |q_i - p_i| <= w, forall bold(q) in X\}$.
  6. *Box Evaluation*: Form $C <- B_w (bold(p), D)$. If
     $|C| >= alpha n$ and $mu(|C|, |D|) > mu_"best"$, update
     $(hat(C), hat(D)) <- (C, D)$ and
     $mu_"best" <- mu(|C|, |D|)$.
  7. *Iterative Extraction*: Output $(hat(C), hat(D))$, remove
     $hat(C)$ from $cal(X)$, and repeat Steps 2--6 up to $k$
     times (or until no cluster meets $|C| >= alpha n$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees: the sharp Monte Carlo
recovery and width-2 approximation theorem of DOC, an original
closed-form derivation of its quasi-polynomial dimension exponent
and phase boundary, an original log-utility proof showing when
optimality automatically implies $beta$-dispersion, a tight
truncation failure bound explaining FastDOC's degradation in high
dimensions, and a complete symmetry/invariance characterization.

#subsec("A", "Monte Carlo Recovery and Width-2 Approximation")

#theorem(
  title: "High-Probability Recovery and Width-2 Approximation",
  origin: "Classical — Cited",
)[
  Let $(C^*, D^*)$ be an $(alpha, beta)$-projective cluster of
  width $w$ with $alpha in (0, 1)$ and $beta in (0, 1 / 2)$.
  Running Steps 1--6 of DOC with
  $r = ceil(ln(2 d) / ln(1 / (2 beta)))$,
  $T_"out" = ceil(alpha^(-1) ln(2 / delta))$, and
  $m = ceil(2 alpha^(-r) ln(2 / delta))$ guarantees that with
  probability at least $1 - delta$, the returned pair
  $(hat(C), hat(D))$ satisfies:
  1. $mu(|hat(C)|, |hat(D)|) >= mu(|C^*|, |D^*|)$;
  2. Every coordinate $i in hat(D)$ has diameter
     $"diam"_i (hat(C)) <= 2 w$;
  3. If $(C^*, D^*)$ is uniquely optimal, some trial recovers
     $D = D^*$ with $C^* subset.eq B_w (bold(p), D^*)$
     @procopiuc2002doc.
]

#proof[
  Consider two failure events. First, let $E_"out"$ be the event
  that at least one of the $T_"out"$ outer medoid samples falls
  inside $C^*$. Since $|C^*| / n >= alpha$ and outer draws are
  independent,
  $
    "Pr"[E_"out"^c]
    = (1 - |C^*| / n)^(T_"out")
    <= (1 - alpha)^(T_"out")
    <= exp(-alpha T_"out")
    <= delta / 2.
  $
  Second, condition on an outer draw with $bold(p) in C^*$, and
  examine a single inner trial drawing $X = \{bold(q)_1, ...,
  bold(q)_r\}$ i.i.d. uniformly from $cal(X)$. Let $A_1$ be the
  event $X subset.eq C^*$. Clearly
  $"Pr"[A_1] = (|C^*| / n)^r >= alpha^r$. Whenever $A_1$ occurs,
  Condition 2 of Definition 1 ensures $|q_i - p_i| <= "diam"_i (C^*)
  <= w$ for every $bold(q) in X$ and every $i in D^*$; therefore
  $D^* subset.eq D$ deterministically.

  Conditioned on $A_1$, the $r$ points of $X$ are i.i.d. uniform
  on $C^*$. Fix any irrelevant coordinate $i in A without D^*$.
  The symmetric window $[p_i - w, p_i + w]$ is the union of two
  closed intervals of length $w$, namely $[p_i - w, p_i]$ and
  $[p_i, p_i + w]$. Applying Condition 3 ($beta$-dispersion) to
  each half-interval gives
  $
    |\{bold(q) in C^* : |q_i - p_i| <= w\}|
    <= 2 beta |C^*|.
  $
  Hence a single random point of $C^*$ satisfies $|q_i - p_i| <= w$
  with probability at most $2 beta$, and all $r$ independent
  points of $X$ do so with probability at most $(2 beta)^r$:
  $
    "Pr"[i in D | A_1] <= (2 beta)^r.
  $
  By a union bound over all $|A without D^*| <= d$ irrelevant
  dimensions and our choice $r >= ln(2 d) / ln(1 / (2 beta))$,
  $
    "Pr"[D != D^* | A_1]
    <= d dot.op (2 beta)^r
    <= d dot.op (1 / (2 d))
    = 1 / 2.
  $
  Therefore, each single inner trial recovers $D = D^*$ exactly
  with probability
  $
    "Pr"[D = D^*]
    >= "Pr"[A_1] dot.op "Pr"[D = D^* | A_1]
    >= (1 / 2) alpha^r.
  $
  Repeating $m >= 2 alpha^(-r) ln(2 / delta)$ independent inner
  trials bounds the conditional miss probability by
  $
    (1 - (1 / 2) alpha^r)^m
    <= exp(-(m / 2) alpha^r)
    <= delta / 2.
  $
  By the union bound, with probability at least $1 - delta$, some
  trial realizes both $bold(p) in C^*$ and $D = D^*$. On that
  trial, Condition 2 implies $C^* subset.eq B_w (bold(p), D^*)$,
  so $|B_w (bold(p), D^*)| >= |C^*| >= alpha n$ and
  $mu(|B_w (bold(p), D^*)|, |D^*|) >= mu(|C^*|, |D^*|)$. Because
  DOC keeps the maximum-$mu$ feasible box across all trials and
  every box $B_w (bold(p)', D')$ satisfies
  $|q_i - q'_i| <= |q_i - p'_i| + |p'_i - q'_i| <= 2 w$ for all
  $i in D'$, all three conclusions follow.
]

Note that our exact proof achieves $m = ceil(2 alpha^(-r) ln(2 / delta))$,
sharpening the original $(2 / alpha)^r ln 4$ trial bound in
@procopiuc2002doc by a factor of $2^(r - 1)$ through conditioning
directly on $X subset.eq C^*$ before invoking the union bound.

#subsec("B", "Exact Polynomial Dimension Exponent and Phase Boundary")

#proposition(
  title: "Dimension-Scaling Exponent and Complexity Phase Boundary",
  origin: "Original Derivation",
)[
  Define the scaling exponent
  $
    gamma(alpha, beta)
    = (ln(1 / alpha)) / (ln(1 / (2 beta))) > 0.
  $
  For any fixed $(alpha, beta, delta) in (0, 1) times (0, 1 / 2)
  times (0, 1)$, the worst-case runtime of one DOC extraction
  satisfies
  $
    T_"DOC" (n, d)
    = cal(O)(n d dot.op alpha^(-2) (2 d)^(gamma(alpha, beta))
      ln^2 (2 / delta)),
  $
  which is strictly linear in $n$ and polynomial in $d$ of exact
  degree $1 + gamma(alpha, beta)$. Furthermore, the critical
  phase boundary $gamma(alpha, beta) = 1$ occurs along the line
  $alpha = 2 beta$: for $alpha < 2 beta$ the dimension degree
  exceeds $2$, and as $beta -> (1 / 2)^-$ the exponent diverges
  as $gamma(alpha, beta) = Theta((1 - 2 beta)^(-1) ln(1 / alpha))$.
]

#proof[
  Since $r = ceil(ln(2 d) / ln(1 / (2 beta))) <
  1 + ln(2 d) / ln(1 / (2 beta))$, exponentiating base $1 / alpha$
  gives
  $
    alpha^(-r)
    < alpha^(-1) exp((ln(2 d) dot.op ln(1 / alpha)) /
      ln(1 / (2 beta)))
    = alpha^(-1) (2 d)^(gamma(alpha, beta)).
  $
  Multiplying $T_"out" = cal(O)(alpha^(-1) ln(2 / delta))$ by
  $m = cal(O)(alpha^(-r) ln(2 / delta))$ and the $cal(O)(n d)$
  scan required to compute $D$ and evaluate $B_w (bold(p), D)$
  (noting $r << n$ so inferring $D$ costs $cal(O)(r d) <= cal(O)(n d)$)
  yields the stated complexity bound. Setting $gamma(alpha, beta) = 1$
  equates $ln(1 / alpha) = ln(1 / (2 beta))$, i.e., $alpha = 2 beta$.
  Finally, writing $1 - 2 beta = eta -> 0^+$ gives
  $ln(1 / (2 beta)) = -ln(1 - eta) = eta + cal(O)(eta^2)$, proving
  $gamma(alpha, beta) = Theta(eta^(-1) ln(1 / alpha))$.
]

Proposition 1 quantifies the practical computational cliff of
DOC: at canonical defaults $alpha = 0.1$ and $beta = 0.25$,
$gamma = ln 10 / ln 2 approx 3.32$, so full DOC scales as
$cal(O)(n d^(4.32))$. Raising $beta$ to $0.45$ pushes the degree
to $1 + 21.85 approx 22.9$, rendering exact sampling infeasible
even for $d = 20$.

#subsec("C", "Log-Utility Equivalence and Automatic Dispersion")

#lemma(
  title: "Log-Linear Objective and Self-Implied Beta-Dispersion",
  origin: "Original Proof",
)[
  Maximizing $mu(|C|, |D|) = |C| (1 / beta)^(|D|)$ is equivalent
  to maximizing the additive log-utility
  $Phi(C, D) = ln |C| + |D| ln(1 / beta)$. Moreover, if
  $(C^*, D^*)$ maximizes $mu$ over all width-$w$ pairs with
  $|C| >= alpha n$ and satisfies the support margin
  $|C^*| >= (alpha / beta) n$, then Condition 3 ($beta$-dispersion
  outside $D^*$) holds *automatically* without being assumed a
  priori.
]

#proof[
  Monotonicity of $t |-> ln t$ establishes the equivalence
  $ln mu(|C|, |D|) = ln |C| + |D| ln(1 / beta) = Phi(C, D)$.
  Now suppose $(C^*, D^*)$ is $mu$-optimal with
  $|C^*| >= (alpha / beta) n$, and assume toward contradiction
  that some coordinate $i in A without D^*$ and closed interval
  $I$ of length $w$ satisfy
  $|C'| > beta |C^*|$, where $C' = \{bold(q) in C^* : q_i in I\}$.
  Then:
  1. $|C'| > beta |C^*| >= beta dot.op (alpha / beta) n = alpha n$,
     so $(C', D^* union \{i\})$ satisfies the support constraint;
  2. $"diam"_j (C') <= "diam"_j (C^*) <= w$ for all $j in D^*$
     and $"diam"_i (C') <= |I| = w$, so $(C', D^* union \{i\})$
     satisfies the width constraint;
  3. Its objective strictly exceeds that of $(C^*, D^*)$:
     $
       mu(|C'|, |D^*| + 1)
       = |C'| (1 / beta)^(|D^*| + 1)
       > beta |C^*| dot.op (1 / beta)^(|D^*| + 1)
       = mu(|C^*|, |D^*|),
     $
     contradicting optimality of $(C^*, D^*)$. Thus
     $|C'| <= beta |C^*|$ for all $i in.not D^*$ and length-$w$
     intervals $I$.
]

#subsec("D", "FastDOC Truncation Failure and Oblique-Box Mismatch")

To escape the $d^(gamma(alpha, beta))$ explosion of Proposition 1,
practical implementations (FastDOC @procopiuc2002doc) cap the
inner loop at a fixed constant `MAXITER` (e.g., $m <= 200$). We
prove a tight lower bound on the resulting miss probability.

#theorem(
  title: "FastDOC Truncation Lower Bound and Oblique Failure",
  origin: "Original Counterexample",
)[
  1. (_Sampling Truncation Miss Bound_.) Let $(C^*, D^*)$ be a
     unique $(alpha, beta)$-cluster of exact size $|C^*| = alpha n$
     embedded in uniform background noise such that any
     discriminating point outside $C^*$ eliminates at least one
     dimension of $D^*$. For any fixed inner budget $m$ and outer
     budget $T_"out"$, the probability that DOC ever draws a pure
     discriminating set $X subset.eq C^*$ is bounded above by
     $
       "Pr"[exists "pure trial"]
       <= T_"out" m dot.op alpha^(r + 1).
     $
     Consequently, if $m$ and $T_"out"$ are held constant while
     $d -> infinity$ (so $r -> infinity$), recovery fails with
     probability $1 - cal(O)(d^(-gamma(alpha, beta))) -> 1$.
  2. (_Diagonal Manifold Collapse_.) In $bb(R)^2$, a 1D line
     segment cluster of length $L > 2 w sqrt(2)$ oriented at
     $45^degree$ has zero intrinsic width orthogonal to its axis,
     yet any axis-aligned box $B_w (bold(p), \{1, 2\})$ captures
     at most a fraction $2 w sqrt(2) / L$ of the cluster.
]

#proof[
  (1) In any single outer-inner pair $(t, s)$, the seed $bold(p)$
  and $r$ discriminating points $bold(q)_1, ..., bold(q)_r$ are
  drawn independently with probability $|C^*| / n = alpha$ each
  of hitting $C^*$. Thus the joint probability that
  $\{bold(p)\} union X subset.eq C^*$ is exactly $alpha^(r + 1)$.
  Applying the union bound across all $T_"out" m$ evaluations
  gives $"Pr"[exists "pure trial"] <= T_"out" m alpha^(r + 1)$.
  Since $alpha^r <= (2 d)^(-gamma(alpha, beta))$, holding
  $T_"out"$ and $m$ fixed drives the upper bound to zero at rate
  $cal(O)(d^(-gamma(alpha, beta)))$.

  (2) Parameterize the $45^degree$ segment as
  $bold(x)(s) = (s / sqrt(2), s / sqrt(2))$ for $s in [0, L]$
  with uniform density $1 / L$. For any seed
  $bold(p) = bold(x)(s_0)$, the $L_infinity$ condition
  $|x_1 (s) - x_1 (s_0)| <= w$ and $|x_2 (s) - x_2 (s_0)| <= w$
  holds iff $|s - s_0| <= w sqrt(2)$, which covers an interval of
  length at most $2 w sqrt(2)$, capturing fraction at most
  $2 w sqrt(2) / L < 1$ of the points.
]

#subsec("E", "Exact Symmetry and Metric Invariance Profile")

#proposition(
  title: "Hyperoctahedral Invariance and Non-Invariance Under Scaling",
  origin: "Original Proof",
)[
  The output distribution of DOC is invariant under: (i) all
  translations $bold(x) |-> bold(x) + bold(t)$, (ii) all signed
  coordinate permutations in the hyperoctahedral group $B_d$ of
  order $2^d d!$, and (iii) joint uniform scaling
  $(bold(x), w) |-> (c bold(x), c w)$ for $c > 0$. It is *not*
  invariant under $O(d)$ rotations or uncompensated anisotropic
  coordinate rescaling.
]

#proof[
  DOC accesses the data exclusively through coordinatewise
  differences $|q_i - p_i| <= w$. Translations leave
  $q_i - p_i$ invariant; sign flips $x_i |-> -x_i$ leave
  $|q_i - p_i|$ invariant; coordinate permutations relabel $A$
  bijectively without changing $|D|$, $|B_w (bold(p), D)|$, or
  $mu(|C|, |D|)$; and scaling both data and $w$ by $c > 0$
  scales $|c q_i - c p_i| <= c w$ equivalently. Non-invariance
  under $45^degree$ rotation follows directly from Theorem 2(2):
  an axis-aligned segment of length $L$ along $x_1$ is captured
  completely in $D = \{2\}$ with $|C| = n$, whereas rotating by
  $45^degree$ destroys boundedness along both coordinates. Under
  anisotropic scaling $x_1 |-> (2 w / epsilon) x_1$ with $x_2$
  fixed, any non-degenerate cluster along axis $1$ loses
  coordinate $1$ from $D$ while $x_2$ is unaffected.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact DOC requires $T_"out" m = cal(O)(alpha^(-2) (2 d)^(gamma(alpha, beta))
ln^2 (2 / delta))$ box evaluations per cluster, each costing
$cal(O)(n d)$ to scan coordinates against the $r$ discriminating
points and count points in $B_w (bold(p), D)$, yielding total
serial cost $cal(O)(k n d dot.op alpha^(-(r + 1)))$ (or
$cal(O)(k n d dot.op T_"out" "MAXITER")$ in FastDOC). Notice the
complete absence of pairwise distance matrices ($cal(O)(n^2)$)
and matrix eigendecompositions ($cal(O)(d^3)$): DOC needs only
$cal(O)(n d)$ RAM for the dataset plus $cal(O)(n + d)$ scratch
arrays, making every inner trial a tight, vectorizable streaming
scan with exemplary cache locality.

Geometrically, DOC extracts axis-parallel hyper-rectangular
slabs ($L_infinity$ balls of side $2 w$ in $D$, unbounded in
$A without D$). Because clusters are extracted greedily one at a
time and removed from $cal(X)$, DOC naturally accommodates
clusters with *distinct* subspace dimensionalities $|D_1| != |D_2|$,
unlike PROCLUS and ORCLUS which enforce a common target
dimension $l$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(n d dot.op alpha^(-(r + 1)))$ exact ($r = cal(O)(log d)$);
    $cal(O)(k n d dot.op T_"out" "MAXITER")$ for FastDOC.],
    [Space complexity],
    [$cal(O)(n d)$ total memory; no $cal(O)(n^2)$ affinity matrix
    or $cal(O)(d^3)$ eigensolver workspace.],
    [Cluster geometry],
    [Axis-aligned hyper-rectangular boxes of uniform width $2 w$
    in cluster-specific subspaces $D_t subset.eq A$.],
    [Noise & overlap],
    [Explicit background noise (unassigned points); greedy
    removal yields disjoint clusters (overlap optional).],
    [Invariance],
    [Invariant to translations, signed axis permutations $B_d$,
    and joint scaling; sensitive to rotations and feature units.],
    [Key hyperparameters],
    [Box width $w > 0$, density fraction $alpha in (0, 1)$,
    balance $beta in (0, 1 / 2)$, iteration cap `"MAXITER"`.],
  ),
  caption: [Operational, geometric, and complexity profile of DOC.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Throughput Tabular Customer and Cohort Segmentation*:
   Retail, telecom, and SaaS behavioral tables contain tens of
   normalized KPIs where actionable segments are defined by crisp
   interval bounds on a subset of metrics (e.g., tenure, spend,
   recency). DOC's $L_infinity$ box representation translates
   directly into human-auditable conjunctive business rules
   ($|x_i - p_i| <= w$ for $i in D$) at sub-millisecond speed.
2. *Large-Scale Streaming Log and Telemetry Summarization*: In
   cloud infrastructure monitoring with millions of event vectors
   ($n >> 10^6$), quadratic pairwise methods (SUBCLU) are
   infeasible. DOC's strictly $cal(O)(n)$-linear streaming scan
   and constant-memory inner loop enable fast extraction of
   dominant fault signatures defined over subsets of system
   counters.
3. *Network Intrusion and Fraud Pattern Discovery*: Attack
   campaigns (port scans, botnet floods, credential stuffing)
   concentrate sharply on a few packet or session features while
   randomizing irrelevant headers to evade full-space detectors.
   DOC's $(alpha, beta)$-cluster model detects dense axis-aligned
   signatures with provable statistical recall (Theorem 1).
4. *E-Commerce Attribute-Subspace Recommendation*: Product
   catalogs span heterogeneous specifications where buyers group
   tightly on a small subset of attributes (price band, weight,
   power rating) while ignoring dozens of irrelevant fields.
   DOC's variable subspace cardinality $|D_t|$ adapts naturally
   to categories with differing numbers of relevant specs.
5. *Bioinformatics Microarray Co-Expression Box Mining*:
   Finding subsets of samples exhibiting tightly bounded gene
   expression within a width-$w$ tolerance across a maximal set
   of genes corresponds directly to maximizing $mu(|C|, |D|)$,
   linking DOC to constant-value biclustering @aggarwal2013data.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DOC suffers four principal failure modes: (i) _High-$d$ FastDOC
undersampling_ (Theorem 2(1)), where fixed `MAXITER` drops
exponentially below $alpha^(-r)$ as $d$ grows; (ii) _Scale
sensitivity_, because a single global width $w$ fails if
features are not standardized to comparable variance; (iii)
_Oblique and non-convex geometry_ (Theorem 2(2)), where rotated
or curved clusters are fragmented by $L_infinity$ hypercubes;
and (iv) _Small minority clusters_ ($|C^*| < alpha n$), which are
invisible to uniform medoid sampling and require lowering $alpha$
at a steep $alpha^(-(r + 1))$ computational penalty.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Standardized tabular rule mining],
    [Strong: fast $cal(O)(n)$ scan; interpretable box rules.],
    [Standardize all features before setting global width $w$.],
    [Variable-dim axis subspaces],
    [Strong: $mu(|C|, |D|)$ adapts $|D_t|$ per cluster.],
    [Tune $beta in [0.15, 0.35]$ to prevent Prop. 1 blowup.],
    [High-$d$ subtle minority clusters],
    [Poor under FastDOC cap (Thm. 2(1) miss bound).],
    [Use P3C @moise2006p3c, STATPC @moise2008statpc, or ORCLUS.],
    [Rotated correlated manifolds],
    [Poor: $L_infinity$ boxes clip oblique segments (Thm. 2(2)).],
    [Use ORCLUS @aggarwal2000orclus (Alg. 138).],
    [Curved / multi-shape densities],
    [Weak: hypercubes cut crescent and ring geometries.],
    [Use SUBCLU @kailing2004subclu (Alg. 137).],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for DOC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DOC occupies a unique theoretical niche among subspace clustering
algorithms @kriegel2009subspacesurvey: whereas CLIQUE
@agrawal1998clique and MAFIA @goil1999mafia quantize space into
fixed or adaptive grids and PROCLUS @aggarwal1999proclus
optimizes $k$-medoid average dispersion, DOC directly optimizes a
probabilistic hypercube density-dimension objective via random
witnesses, directly inspiring deterministic statistical
successors P3C @moise2006p3c and STATPC @moise2008statpc.

#proposition(
  title: "Exact Reduction Between DOC Boxes and L-Infinity DBSCAN Cores",
  origin: "Original Proof",
)[
  Let $d_D^infinity (bold(x), bold(y)) = max_(i in D) |x_i - y_i|$
  be the $L_infinity$ subspace metric on $D subset.eq A$. For any
  seed $bold(p) in cal(X)$ and subspace $D$, DOC's candidate
  cluster $B_w (bold(p), D)$ coincides identically with the
  $L_infinity$ $w$-neighborhood $N_(w, infinity)^D (bold(p))$ in
  SUBCLU. Consequently, DOC's support feasibility test
  $|B_w (bold(p), D)| >= alpha n$ holds if and only if $bold(p)$
  is a SUBCLU core point in subspace $D$ under $L_infinity$ radius
  $epsilon = w$ and density threshold $mu = ceil(alpha n)$.
]

#proof[
  By definition, $bold(q) in B_w (bold(p), D)$ holds iff
  $|q_i - p_i| <= w$ for all $i in D$, which is equivalent to
  $max_(i in D) |q_i - p_i| = d_D^infinity (bold(p), bold(q)) <= w$,
  i.e., $bold(q) in N_(w, infinity)^D (bold(p))$. Hence
  $|B_w (bold(p), D)| >= alpha n <=> |N_(w, infinity)^D (bold(p))|
  >= ceil(alpha n)$, which is precisely the $L_infinity$ core-point
  predicate of SUBCLU (Definition 1 of Alg. 137). Unlike SUBCLU,
  however, DOC does *not* chain overlapping core balls into
  non-convex connected components; it stops at a single maximal
  core ball maximizing $mu(|C|, |D|)$.
]

The ClusterOverview 4-dataset empirical benchmark
(@rousseeuw1987silhouettes, @hubert1985comparing) aligns closely
with our theoretical bounds. On *Iris*, DOC achieves
$"Sil" = 0.350$ (higher Silhouette than both SUBCLU $0.315$ and
ORCLUS $0.292$) yet modest ground-truth agreement
($"ARI" = 0.398$ vs. ORCLUS $0.681$ and SUBCLU $0.516$): its
axis-aligned hypercube cleanly isolates the compact _Setosa_ box
but partitions the obliquely correlated _Versicolor_ / _Virginica_
pair imperfectly (Theorem 2(2)). On *Customers* ($"Sil" = 0.214$),
DOC outperforms SUBCLU ($0.170$) and MAFIA ($"N/A"$) while
trailing ORCLUS ($0.318$). On *Multi-Shape* ($"Sil" = 0.312$,
$"ARI" = "N/A"$), stopping at single $L_infinity$ core balls
without density chaining (Proposition 3) leaves DOC well behind
SUBCLU ($0.584$). Most tellingly, on *High-D* ($"Sil" = 0.068$,
$"ARI" = "N/A"$), DOC scores lowest among all four algorithms—a
direct empirical manifestation of Theorem 2(1), where practical
inner-loop truncation under-samples the exponentially growing
discriminating-set space as $d$ increases. In exchange, DOC
delivers ultra-fast execution at $0.112 "ms"$ mean runtime (25x
faster than ORCLUS and 111x faster than SUBCLU), trailing only
MAFIA's histogram pass ($0.058 "ms"$).
