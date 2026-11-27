#import "../template/ieee.typ": *

#alg_entry(
  id: 137,
  name: "SUBCLU: Density-Connected Subspace Clustering",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2004",
  paradigm: "Bottom-Up Density-Connected Subspace Search",
  complexity: $cal(O)(|cal(S)_"cand"| dot.op n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

High-dimensional datasets routinely harbor clusters that exist
only within low-dimensional axis-parallel projections, while
the remaining attributes behave as irrelevant background noise
that destroys full-space distance contrast @beyer1999nearest.
Grid-based subspace algorithms such as CLIQUE @agrawal1998clique
and MAFIA @goil1999mafia mitigate this curse of dimensionality
by discretizing each axis into bins, yet rigid hyper-rectangular
partitions inevitably sever non-convex clusters and lose density
resolution near bin boundaries. SUBCLU (Density-Connected
Subspace Clustering), introduced by Kailing, Kriegel, and
Kröger @kailing2004subclu, resolves this tension by lifting the
density-connected cluster model of DBSCAN @ester1996dbscan
directly to the power set of coordinate subspaces while
preserving Apriori-style bottom-up pruning.

Let $cal(X) = \{bold(x)_1, ..., bold(x)_n\} subset bb(R)^d$
denote a finite dataset over attribute index set
$A = \{1, ..., d\}$. For any non-empty subspace
$S subset.eq A$ with cardinality $m = |S|$, define the
canonical orthogonal projection
$pi_S : bb(R)^d -> bb(R)^m$ retaining only coordinates indexed
by $S$, and induce the subspace $L_p$ metric ($p in [1, infinity]$)
$
  d_S (bold(x), bold(y))
  &= ||pi_S (bold(x)) - pi_S (bold(y))||_p \
  &= (sum_(j in S) |x_j - y_j|^p)^(1 / p),
$
taking $p = 2$ (Euclidean) as the canonical setting. Given
global neighborhood radius $epsilon > 0$ and minimum cardinality
threshold $mu in bb(N)_(>= 2)$ (denoted $"MinPts"$ in
@ester1996dbscan), the subspace $epsilon$-neighborhood of
$bold(x) in cal(X)$ in projection $S$ is
$
  N_epsilon^S (bold(x))
  = \{bold(y) in cal(X) : d_S (bold(x), bold(y)) <= epsilon\}.
$

#definition(title: "Density-Connected Subspace Cluster")[
  Fix $S subset.eq A$ ($S != emptyset$), $epsilon > 0$, and
  $mu >= 2$. An object $bold(x) in cal(X)$ is a _core object_
  in $S$, written $"Core"_(epsilon, mu)^S (bold(x))$, if
  $|N_epsilon^S (bold(x))| >= mu$. An object $bold(y)$ is
  _directly density-reachable_ from $bold(x)$ in $S$
  ($bold(y) <-_S bold(x)$) if $"Core"_(epsilon, mu)^S (bold(x))$
  and $bold(y) in N_epsilon^S (bold(x))$. It is
  _density-reachable_ ($bold(y) <-_S^* bold(x)$) if there exists
  a finite chain $bold(x) = bold(z)_0, bold(z)_1, ..., bold(z)_r
  = bold(y)$ in $cal(X)$ with $bold(z)_t <-_S bold(z)_(t-1)$ for
  all $t = 1, ..., r$. Two points $bold(x), bold(y)$ are
  _density-connected_ in $S$ ($bold(x) <->_S bold(y)$) if some
  $bold(o) in cal(X)$ satisfies both $bold(x) <-_S^* bold(o)$
  and $bold(y) <-_S^* bold(o)$. A non-empty subset
  $C subset.eq cal(X)$ is a _density-connected subspace cluster_
  in $S$ if it satisfies:
  1. _Maximality_: $bold(x) in C$ and $bold(y) <-_S^* bold(x)$
     imply $bold(y) in C$;
  2. _Connectivity_: $bold(x) <->_S bold(y)$ for every
     $bold(x), bold(y) in C$.
  The collection of all such clusters in $S$ is denoted
  $cal(C)(S)$, and its support union is
  $U(S) = union.big_(C in cal(C)(S)) C$.
]

SUBCLU seeks the complete family
$\{(S, cal(C)(S)) : emptyset != S subset.eq A, cal(C)(S) != emptyset\}$
of all non-empty subspaces containing at least one
density-connected cluster. Because enumerating all $2^d - 1$
subspaces is intractable for moderate $d$, SUBCLU exploits a
monotonicity theorem analogous to CLIQUE's downward closure:
every cluster in a $(k + 1)$-dimensional subspace $S$ must be
contained inside a cluster of every $k$-dimensional projection
$T subset S$. Crucially, SUBCLU does not merely prune the
subspace lattice; it also restricts the subsequent DBSCAN call
on candidate $S$ to the points $U(T^*)$ surviving in the
minimum-support parent $T^* subset S$, eliminating noise points
that cannot possibly join any cluster in $S$.

#algbox(
  title: "SUBCLU Bottom-Up Subspace Search",
  inputs: [Dataset $cal(X) subset bb(R)^d$, radius $epsilon > 0$,
    density threshold $mu >= 2$.],
  outputs: [All subspaces $S subset.eq A$ with $cal(C)(S) != emptyset$
    and their cluster sets $cal(C)(S)$.],
)[
  1. *1D Initialization*: For each $j in \{1, ..., d\}$, run
     $"DBSCAN"(cal(X), \{j\}, epsilon, mu)$ to obtain
     $cal(C)(\{j\})$. Set
     $cal(S)_1 <- \{\{j\} : cal(C)(\{j\}) != emptyset\}$ and
     $k <- 1$.
  2. *Level-Wise Loop*: While $cal(S)_k != emptyset$:
  3. *Apriori Candidate Generation*: Form
     $"Cand"_(k+1)$ by joining pairs $S_1, S_2 in cal(S)_k$ that
     agree on their first $k - 1$ sorted indices, then prune any
     candidate $S$ possessing a $k$-subset $T subset S$ with
     $T in.not cal(S)_k$.
  4. *Minimum-Support Parent Selection*: For each
     $S in "Cand"_(k+1)$, select
     $T^* = arg min_(T subset S, |T| = k) |U(T)|$.
  5. *Restricted Density Clustering*: For each cluster
     $C in cal(C)(T^*)$, invoke $"DBSCAN"(C, S, epsilon, mu)$
     and collect all discovered clusters into $cal(C)(S)$.
  6. *Update Active Lattice*: If $cal(C)(S) != emptyset$, insert
     $S$ into $cal(S)_(k+1)$. Increment $k <- k + 1$.
  7. *Return*: Union of $(S, cal(C)(S))$ across all levels
     $1 <= k <= d$.
]

Notice step 5 executes DBSCAN separately on each individual
cluster $C in cal(C)(T^*)$ rather than on the full dataset
$cal(X)$. As proved in Section II, this decomposition is exact:
no cluster in $S$ can straddle two distinct clusters of $T^*$ or
include points outside $U(T^*)$.

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees governing SUBCLU: the
foundational monotonicity and exactness theorems of
@kailing2004subclu, an original dimensional volume-collapse
bound proving the finite subspace horizon under fixed $epsilon$,
a constructive if-and-only-if counterexample showing why naive
dimension-adaptive radii break Apriori pruning, and an original
analysis of multi-parent intersection pruning.

#subsec("A", "Monotonicity and Downward Closure of Density Connectivity")

#lemma(
  title: "Subspace Metric Monotonicity and Neighborhood Inclusion",
  origin: "Classical — Cited",
)[
  Let $d_S$ be induced by any $L_p$ norm ($p in [1, infinity]$).
  For all non-empty subspaces $T subset.eq S subset.eq A$, all
  $bold(x), bold(y) in cal(X)$, and every $epsilon > 0$
  @kailing2004subclu:
  $
    d_T (bold(x), bold(y)) <= d_S (bold(x), bold(y)),
    quad N_epsilon^S (bold(x)) subset.eq N_epsilon^T (bold(x)).
  $
]

#proof[
  For finite $p >= 1$, since every summand $|x_j - y_j|^p >= 0$
  and $T subset.eq S$,
  $
    d_T (bold(x), bold(y))^p
    = sum_(j in T) |x_j - y_j|^p
    <= sum_(j in S) |x_j - y_j|^p
    = d_S (bold(x), bold(y))^p.
  $
  Taking the strictly increasing $(1/p)$-th root yields
  $d_T (bold(x), bold(y)) <= d_S (bold(x), bold(y))$. For
  $p = infinity$, taking the maximum over a smaller index set
  gives $max_(j in T) |x_j - y_j| <= max_(j in S) |x_j - y_j|$.
  Consequently, any $bold(y) in N_epsilon^S (bold(x))$ satisfies
  $d_T (bold(x), bold(y)) <= d_S (bold(x), bold(y)) <= epsilon$,
  hence $bold(y) in N_epsilon^T (bold(x))$.
]

#theorem(
  title: "Downward Closure and Exactness of Restricted DBSCAN",
  origin: "Classical — Cited",
)[
  For fixed $(epsilon, mu)$ and any non-empty
  $T subset.eq S subset.eq A$, if $C in cal(C)(S)$ is a
  density-connected subspace cluster in $S$, then there exists a
  unique cluster $C' in cal(C)(T)$ such that $C subset.eq C'$
  @kailing2004subclu. Moreover, running $"DBSCAN"(C', S, epsilon,
  mu)$ independently over each $C' in cal(C)(T)$ recovers the
  exact cluster collection $cal(C)(S)$ produced by
  $"DBSCAN"(cal(X), S, epsilon, mu)$ (modulo border-point
  tie-breaking between overlapping core hulls).
]

#proof[
  By Lemma 1, $|N_epsilon^T (bold(x))| >= |N_epsilon^S (bold(x))|$
  for every $bold(x) in cal(X)$, so
  $"Core"_(epsilon, mu)^S (bold(x)) => "Core"_(epsilon, mu)^T (bold(x))$.
  Furthermore, $bold(y) <-_S bold(x)$ requires both
  $"Core"_(epsilon, mu)^S (bold(x))$ and
  $bold(y) in N_epsilon^S (bold(x)) subset.eq N_epsilon^T (bold(x))$,
  which implies $bold(y) <-_T bold(x)$. By induction on chain
  length $r$, $bold(y) <-_S^* bold(x) => bold(y) <-_T^* bold(x)$,
  and therefore $bold(x) <->_S bold(y) => bold(x) <->_T bold(y)$.

  Now let $C in cal(C)(S)$. Since $C$ is non-empty and
  density-connected in $S$, it contains at least one core point
  $bold(o)$ in $S$ (hence also a core point in $T$) such that
  every $bold(x) in C$ satisfies $bold(x) <-_S^* bold(o)$, thus
  $bold(x) <-_T^* bold(o)$. Let $C' in cal(C)(T)$ be the maximal
  density-connected cluster in $T$ containing the core point
  $bold(o)$. By maximality of $C'$ in $T$, every
  $bold(x) in C$ belongs to $C'$, proving $C subset.eq C'$.
  Uniqueness of $C'$ follows because two distinct clusters in
  $cal(C)(T)$ cannot share the core point $bold(o)$.

  For exactness of restricted DBSCAN on $C'$, observe that for
  every $S$-core point $bold(x) in C subset.eq C'$, every
  neighbor $bold(y) in N_epsilon^S (bold(x))$ satisfies
  $bold(y) <-_S bold(x) => bold(y) <-_T bold(x)$, so
  $bold(y) in C'$ by maximality of $C'$. Hence
  $N_epsilon^S (bold(x)) inter C' = N_epsilon^S (bold(x))$:
  restricting queries to $C'$ loses zero neighbors of any
  $S$-core point in $C'$. Conversely, restricting to $C'$ can
  never increase any neighborhood cardinality, so no non-core
  point in $cal(X)$ becomes a core point inside $C'$. Thus the
  $S$-core points and their direct-reachability graph inside $C'$
  are identical to those in $cal(X)$, proving exact recovery of
  all core-induced clusters in $cal(C)(S)$.
]

#subsec("B", "Volume-Collapse Dimensional Horizon Under Fixed Radius")

While Theorem 1 makes SUBCLU's lattice pruning exact, keeping
$epsilon$ constant across all dimensions $m = |S|$ exposes the
algorithm to rapid spherical volume collapse. We quantify the
exact rate at which fixed-$epsilon$ core points vanish as
subspace dimensionality grows.

#proposition(
  title: "Super-Exponential Fixed-Radius Dimensional Decay",
  origin: "Original Proof",
)[
  Suppose a cluster $C subset.eq cal(X)$ of size $n_C = |C|$ has
  i.i.d. isotropic Gaussian projections
  $pi_S (bold(X)_i) tilde.op cal(N)(bold(mu)_S, sigma^2 bold(I)_m)$
  in an $m$-dimensional subspace $S$ with $L_2$ metric. Then for
  any fixed radius $epsilon > 0$ and threshold $mu >= 2$:
  $
    bb(E)[|N_epsilon^S (bold(X)_1) inter C| - 1]
    <= (n_C - 1) dot.op
    ((epsilon / (2 sigma))^m) / Gamma(m / 2 + 1).
  $
  Moreover, applying Stirling's bound
  $Gamma(m / 2 + 1) >= sqrt(pi m) (m / (2 e))^(m / 2)$, the
  probability that $C$ contains even a single core point in $S$
  satisfies
  $
    "Pr"[exists bold(x) in C : |N_epsilon^S (bold(x)) inter C| >= mu]
    <= (n_C (n_C - 1)) / ((mu - 1) sqrt(pi m))
    dot.op ((e epsilon^2) / (2 m sigma^2))^(m / 2),
  $
  which decays super-exponentially to zero once
  $m > e epsilon^2 / (2 sigma^2)$.
]

#proof[
  For any distinct pair $i != j$ in $C$, let
  $bold(Z) = pi_S (bold(X)_i) - pi_S (bold(X)_j)$. By independence,
  $bold(Z) tilde.op cal(N)(bold(0), 2 sigma^2 bold(I)_m)$, whose
  probability density is uniformly upper-bounded by its peak at
  the origin:
  $
    p(bold(z))
    = (4 pi sigma^2)^(-m / 2) exp(-||bold(z)||_2^2 / (4 sigma^2))
    <= (4 pi sigma^2)^(-m / 2).
  $
  Integrating over the $m$-dimensional Euclidean ball
  $B_m (bold(0), epsilon)$ of Lebesgue volume
  $V_m epsilon^m = pi^(m / 2) epsilon^m / Gamma(m / 2 + 1)$ yields
  $
    "Pr"[d_S (bold(X)_i, bold(X)_j) <= epsilon]
    &= integral_(B_m (bold(0), epsilon)) p(bold(z)) dif bold(z) \
    &<= (4 pi sigma^2)^(-m / 2) dot.op
    (pi^(m / 2) epsilon^m) / Gamma(m / 2 + 1) \
    &= ((epsilon / (2 sigma))^m) / Gamma(m / 2 + 1),
  $
  where the factors $pi^(m/2)$ and $pi^(-m/2)$ cancel identically.
  Summing over the $n_C - 1$ indices $j != 1$ by linearity of
  expectation proves the first inequality.

  For the second bound, let
  $Y_i = |N_epsilon^S (bold(X)_i) inter C| - 1 >= 0$ be the
  number of distinct cluster neighbors of $bold(X)_i$. A point
  $bold(X)_i$ is a core point inside $C$ iff $Y_i >= mu - 1$.
  Applying a union bound over $i in \{1, ..., n_C\}$ followed by
  Markov's inequality gives
  $
    "Pr"[exists i : Y_i >= mu - 1]
    <= sum_(i=1)^(n_C) (bb(E)[Y_i]) / (mu - 1)
    <= (n_C (n_C - 1)) / (mu - 1) dot.op
    ((epsilon / (2 sigma))^m) / Gamma(m / 2 + 1).
  $
  Substituting Stirling's lower bound
  $Gamma(m / 2 + 1) >= sqrt(pi m) (m / (2 e))^(m / 2)$ (valid for
  all $m >= 1$) simplifies the factor
  $
    (epsilon / (2 sigma))^m / (m / (2 e))^(m / 2)
    = ((epsilon^2 / (4 sigma^2)) dot.op (2 e / m))^(m / 2)
    = ((e epsilon^2) / (2 m sigma^2))^(m / 2),
  $
  completing the proof. Whenever $m > e epsilon^2 / (2 sigma^2)$,
  the base is strictly less than $1$ and shrinks as $cal(O)(m^(-1/2))$,
  forcing the core-existence probability to collapse faster than
  any exponential in $m$.
]

Proposition 1 formalizes a fundamental dilemma of SUBCLU: a
small $epsilon$ tailored to clean 1D and 2D projections destroys
all high-dimensional clusters above the critical horizon
$m^* = Theta(epsilon^2 / sigma^2)$, whereas choosing $epsilon$
large enough to retain $m$-dimensional clusters collapses
low-dimensional projections into giant merged noise blobs.

#subsec("C", "Impossibility of Monotone Dimension-Adaptive Radii")

A natural engineering impulse is to counteract Proposition 1 by
equipping SUBCLU with a dimension-dependent radius schedule
$epsilon_1, epsilon_2, ..., epsilon_d$ (for instance
$epsilon_m = epsilon_1 sqrt(m)$ to match expected Euclidean norm
growth). The next theorem proves that *any* strict increase in
$epsilon_m$ destroys downward closure and invalidates Apriori
pruning.

#theorem(
  title: "Exact Characterization of Downward-Closed Radius Schedules",
  origin: "Original Counterexample",
)[
  Let $epsilon_1, ..., epsilon_d > 0$ be a radius schedule indexed
  by subspace cardinality $|S|$, with fixed $mu >= 2$ and any
  $L_p$ metric ($p in [1, infinity]$). Downward closure of
  density-connected clusters holds universally across all finite
  datasets $cal(X) subset bb(R)^d$ if and only if the schedule
  is non-increasing:
  $
    epsilon_1 >= epsilon_2 >= ... >= epsilon_d.
  $
]

#proof[
  (_Sufficiency_.) Suppose $epsilon_1 >= ... >= epsilon_d$, and
  let $emptyset != T subset.eq S subset.eq A$. Then
  $|T| <= |S|$, so $epsilon_(|S|) <= epsilon_(|T|)$. Combined with
  Lemma 1, any $bold(y) in N_(epsilon_(|S|))^S (bold(x))$ obeys
  $d_T (bold(x), bold(y)) <= d_S (bold(x), bold(y)) <= epsilon_(|S|)
  <= epsilon_(|T|)$, whence
  $N_(epsilon_(|S|))^S (bold(x)) subset.eq N_(epsilon_(|T|))^T (bold(x))$.
  Repeating the proof of Theorem 1 verbatim establishes downward
  closure for all datasets.

  (_Necessity via Constructive Counterexample_.) Suppose to the
  contrary that $epsilon_(k+1) > epsilon_k$ for some
  $1 <= k <= d - 1$. Define the positive gap parameter
  $eta = (epsilon_(k+1) - epsilon_k) / 4 > 0$, and let
  $m_1 = floor(mu / 2)$ and $m_2 = ceil(mu / 2)$ so that
  $m_1 + m_2 = mu$. Since $mu >= 2$, both $m_1 <= mu - 1$ and
  $m_2 <= mu - 1$. Construct a dataset $cal(X)$ of $n = mu$
  distinct points embedded along coordinate $1$ (setting
  coordinates $2, ..., d$ identically to $0$):
  - Place $m_1$ distinct points in $[0, eta]$ along axis $1$;
  - Place $m_2$ distinct points in
    $[epsilon_k + 2 eta, epsilon_k + 3 eta]$ along axis $1$.
  Set $T = \{1, ..., k\}$ and $S = \{1, ..., k + 1\}$. Because
  all coordinates beyond index $1$ vanish,
  $d_T (bold(x), bold(y)) = d_S (bold(x), bold(y)) = |x_1 - y_1|$
  for every $p in [1, infinity]$.

  In subspace $T$ (radius $epsilon_k$), the separation between the
  two groups is at least
  $(epsilon_k + 2 eta) - eta = epsilon_k + eta > epsilon_k$.
  Hence every point in the first group has neighborhood cardinality
  $m_1 < mu$, and every point in the second group has cardinality
  $m_2 < mu$. No core point exists in $T$, so
  $cal(C)(T) = emptyset$, and Apriori pruning permanently discards
  the $(k + 1)$-subspace $S$.

  In subspace $S$ (radius $epsilon_(k+1)$), however, the maximum
  pairwise distance across the entire dataset is at most
  $
    epsilon_k + 3 eta
    = (epsilon_k + 3 epsilon_(k+1)) / 4
    < epsilon_(k+1).
  $
  Consequently, every point in $cal(X)$ has all $mu$ points inside
  its $epsilon_(k+1)$-neighborhood in $S$, forming a valid
  density-connected cluster $cal(C)(S) = \{cal(X)\} != emptyset$.
  Downward closure fails, proving necessity.
]

Theorem 2 shows that SUBCLU's fixed-radius bottleneck is not an
accidental oversight, but a structural barrier inherent to
Apriori pruning over density-connected graphs. Overcoming it
requires abandoning strict downward closure via top-down or
statistical approaches such as INSCY @assent2008inscy or
FIRES @kriegel2005fires.

#subsec("D", "Multi-Parent Support Intersection and Tight Pruning")

Standard SUBCLU (step 4) restricts DBSCAN on candidate $S$ to
the single smallest parent $T^* subset S$. We prove that
intersecting *all* $k + 1$ parent supports remains exact and can
yield an $Omega(n)$ reduction in the DBSCAN input size.

#proposition(
  title: "Exactness and Asymptotic Gain of All-Parent Intersection",
  origin: "Original Derivation",
)[
  For any candidate $(k + 1)$-subspace $S in "Cand"_(k+1)$,
  define the multi-parent intersection
  $I(S) = inter.big_(T subset S, |T| = k) U(T)$. Then:
  1. Running $"DBSCAN"(I(S), S, epsilon, mu)$ returns the exact
     cluster collection $cal(C)(S)$ (modulo border tie-breaking);
  2. There exist families of datasets of arbitrary size
     $n = 2 m + mu$ where every single parent satisfies
     $|U(T)| = (n + mu) / 2 = Theta(n)$, yet $|I(S)| = mu = cal(O)(1)$,
     reducing the exact DBSCAN cost on $S$ from $Theta(n^2)$ to
     $cal(O)(n + mu^2)$.
]

#proof[
  (1) By Theorem 1, every cluster $C in cal(C)(S)$ satisfies
  $C subset.eq U(T)$ for every $k$-subset $T subset S$, hence
  $C subset.eq I(S)$. In particular, if $bold(x) in C$ is an
  $S$-core point in $cal(X)$, then every neighbor
  $bold(y) in N_epsilon^S (bold(x))$ is directly density-reachable
  from $bold(x)$ in $S$ and therefore belongs to $C subset.eq I(S)$.
  Thus $N_epsilon^S (bold(x)) inter I(S) = N_epsilon^S (bold(x))$
  for every $S$-core point, while shrinking the point set from
  $cal(X)$ to $I(S)$ cannot turn any non-core point into a core
  point. Therefore $"DBSCAN"(I(S), S, epsilon, mu)$ recovers
  $cal(C)(S)$ exactly.

  (2) Let $k = 1$, $S = \{1, 2\}$, and fix $mu >= 2$ and
  $m in bb(N)$. Place $mu$ core points at the origin $(0, 0)$,
  $m$ points at $(0, 3 j epsilon)$ for $j = 1, ..., m$, and
  $m$ points at $(3 j epsilon, 0)$ for $j = 1, ..., m$, giving
  $n = 2 m + mu$. In $T_1 = \{1\}$, the origin group plus the
  $m$ points on the vertical axis all project to $x_1 = 0$, so
  $|U(\{1\})| = m + mu = (n + mu) / 2$; the remaining $m$ points
  are separated by $3 epsilon > epsilon$ and form noise.
  Symmetrically, $|U(\{2\})| = m + mu = (n + mu) / 2$. Single-parent
  SUBCLU therefore runs DBSCAN on $(n + mu)/2 = Theta(n)$ points
  at $Theta(n^2)$ pairwise cost. However, the intersection
  $I(\{1, 2\}) = U(\{1\}) inter U(\{2\})$ contains solely the
  $mu$ origin points, computable in $cal(O)(n)$ bitset operations,
  slashing the subsequent DBSCAN pass to $cal(O)(mu^2)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

SUBCLU's computational profile is governed by the product of the
number of candidate subspaces $|cal(S)_"cand"|$ generated across
all levels $1 <= k <= k_"max"$ and the cost of invoking DBSCAN on
the pruned point set $U(T^*)$ (or $I(S)$). Without spatial
indexing, each DBSCAN pass requires $cal(O)(|U(T^*)|^2 dot.op k)$
time, giving total runtime
$cal(O)(sum_(S in cal(S)_"cand") |U(T^*(S))|^2 |S|)$, which is
bounded above by $cal(O)(|cal(S)_"cand"| n^2 d)$ and in the
pathological worst case (when every projection is dense) reaches
$cal(O)(2^d n^2 d)$. When low-dimensional projections ($|S| <= 6$)
are indexed by an R\*-tree or $k$-d tree, neighborhood queries
accelerate to $cal(O)(n log n)$ per subspace; in higher
subspaces, tree indices degenerate to linear scans
@schubert2017dbscanrevisited.

Space complexity is dominated by caching cluster assignments
across active level-$k$ subspaces to support step 4 and step 5,
requiring $cal(O)(|cal(S)_k| n)$ memory per level (at most
$cal(O)(binom(d, k) n)$ in the worst case), plus $cal(O)(n d)$
for the raw data matrix. Unlike grid methods that inspect $m^k$
cells independent of $n$, SUBCLU's memory scales with $n$ and the
number of *non-empty* subspaces rather than grid resolution, but
redundant sub-clusters appear across every projection of a true
$m$-dimensional cluster ($2^m - 1$ nested projections), creating
heavy output redundancy unless post-filtered for maximality or
depth-first pruned as in INSCY @assent2008inscy.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(|cal(S)_"cand"| n^2)$ naive; $cal(O)(|cal(S)_"cand"| n log n)$
    with spatial index in low $|S|$; worst-case $cal(O)(2^d n^2)$.],
    [Space complexity],
    [$cal(O)(n d + |cal(S)_k| n)$ active level cache; avoids
    $cal(O)(m^k)$ grid allocation.],
    [Cluster geometry],
    [Arbitrary non-convex, curved, and elongated manifolds in
    any axis-parallel subspace $S subset.eq A$.],
    [Noise & overlap],
    [Explicit background noise rejection; allows overlapping
    clusters across distinct subspaces.],
    [Invariance],
    [Invariant to translations and coordinate permutations;
    sensitive to anisotropic scaling and rotations.],
    [Key hyperparameters],
    [Global radius $epsilon > 0$ and density count $mu >= 2$
    (fixed across all subspace dimensions $|S|$).],
  ),
  caption: [Operational, geometric, and complexity profile of SUBCLU.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multi-Omics and Gene Expression Subspace Discovery*: In
   microarray and single-cell transcriptomics, co-regulated
   pathways manifest as non-linear density manifolds across small
   subsets of genes ($|S| in [3, 15]$) out of thousands of
   uninformative background probes. SUBCLU's exact density
   connectivity captures curved regulatory trajectories without
   the rectangular quantization artifacts of CLIQUE
   @agrawal1998clique.
2. *Sensor-Network and Telemetry Anomaly Isolation*: Industrial
   turbines and avionics buses stream dozens of correlated
   channels where incipient fault modes induce dense operating
   regimes in 2D--5D subsets of sensors (e.g., pressure vs.
   vibration harmonics). SUBCLU isolates localized non-convex
   operating regimes while treating transient spikes as explicit
   DBSCAN noise.
3. *Financial Microstructure and Multi-Factor Regime Mining*:
   Cross-asset order-book descriptors (spread, depth imbalance,
   realized volatility, skew) exhibit heavy-tailed, non-ellipsoidal
   regimes active on subsets of liquidity factors. Because SUBCLU
   imposes no Gaussian covariance assumption, it recovers
   asymmetric liquidity pockets cleanly without requiring prior
   knowledge of $k$.
4. *Astronomical Photometric and Spectroscopic Surveys*: Catalogs
   combining multi-band magnitudes, proper motions, and line
   indices contain stellar streams and kinematic substructures
   curving through specific photometric subspaces. Density
   connectivity traces filamentary stellar streams that
   axis-aligned hypercubes (DOC, CLIQUE) fragment.
5. *Cheminformatics and Structure-Activity Subspace Profiling*:
   Molecular descriptor tables contain hundreds of topological and
   physicochemical features where a given pharmacological
   mechanism depends on only 3--6 descriptors. SUBCLU uncovers
   overlapping activity-cliff neighborhoods across multiple
   descriptor subspaces simultaneously, matching polypharmacology
   where one molecule participates in distinct subspace clusters.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SUBCLU exhibits four primary failure modes rooted directly in
Section II: (i) _Fixed-$epsilon$ dimensional collapse_
(Proposition 1), where no single $epsilon$ simultaneously resolves
dense 1D projections and sparse 8D+ clusters; (ii) _Exponential
subspace explosion_, where a single true 20D cluster spawns
$2^(20) - 1 > 10^6$ redundant projections; (iii) _Varying
intra-subspace densities_, inherited from DBSCAN
@ester1996dbscan, causing dense clusters to merge or sparse ones
to vanish; and (iv) _Oblique correlated manifolds_, where
clusters live on rotated linear planes not aligned with coordinate
axes, favoring ORCLUS @aggarwal2000orclus.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Non-convex low-$|S|$ manifolds],
    [Exact density connectivity avoids grid boundary slicing.],
    [Tune $epsilon$ via 1D/2D $k$-dist plots; HDBSCAN @campello2013hdbscan
    in full space.],
    [Multi-factor overlapping regimes],
    [Naturally outputs overlapping memberships across subspaces.],
    [Post-prune redundant lower-D projections; INSCY @assent2008inscy.],
    [High-$d$ dense clusters ($|S| > 12$)],
    [Poor: $2^(|S|)$ lattice blowup and Prop. 1 volume decay.],
    [Use FIRES @kriegel2005fires, PROCLUS @aggarwal1999proclus, or DOC
    @procopiuc2002doc.],
    [Rotated / correlated hyperplanes],
    [Poor: axis-parallel projections smear oblique correlations.],
    [Use ORCLUS @aggarwal2000orclus (Alg. 138).],
    [Heterogeneous cluster densities],
    [Moderate: single $(epsilon, mu)$ misses sparse clusters.],
    [Use OPTICS @ankerst1999optics or multi-scale density methods.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SUBCLU.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SUBCLU sits at the exact theoretical intersection of full-space
density-based clustering (DBSCAN @ester1996dbscan, GDBSCAN
@sander1998gdbscan) and bottom-up subspace search (CLIQUE
@agrawal1998clique, MAFIA @goil1999mafia). Removing grid
discretization strictly generalizes CLIQUE's density model:

#proposition(
  title: "Grid-Refinement Embedding Between CLIQUE and SUBCLU",
  origin: "Original Proof",
)[
  Let $cal(G)_h$ be a CLIQUE grid of uniform cell width $h$ in an
  $m$-dimensional subspace $S$ with cell density threshold $mu$.
  Under Euclidean radius $epsilon = h sqrt(m + 3)$ and density
  threshold $mu$, every point inside a dense CLIQUE cell
  $u in cal(G)_h$ satisfies $"Core"_(epsilon, mu)^S (bold(x))$,
  and every connected component of face-adjacent dense CLIQUE
  cells in $S$ is wholly contained inside a single SUBCLU cluster
  in $cal(C)(S)$.
]

#proof[
  For any cell $u in cal(G)_h$ of edge length $h$ in $bb(R)^m$,
  its Euclidean diameter is $h sqrt(m) < h sqrt(m + 3) = epsilon$.
  If $u$ is dense ($|u inter cal(X)| >= mu$), every point
  $bold(x) in u inter cal(X)$ satisfies
  $|N_epsilon^S (bold(x))| >= |u inter cal(X)| >= mu$, so every
  point of $u$ is a SUBCLU core point and all points of $u$ are
  mutually directly density-reachable. Furthermore, if two dense
  cells $u, v in cal(G)_h$ share an $(m - 1)$-dimensional face,
  any pair $bold(x) in u$ and $bold(y) in v$ differs by at most
  $2 h$ along the normal coordinate separating $u$ and $v$ and
  by at most $h$ along each of the remaining $m - 1$ coordinates;
  hence
  $
    d_S (bold(x), bold(y))^2
    <= (2 h)^2 + (m - 1) h^2
    = (m + 3) h^2 = epsilon^2.
  $
  Thus every point in $u$ is directly density-reachable from
  every point in $v$ (and vice versa). By induction along paths
  in CLIQUE's cell-adjacency graph, every connected component of
  dense CLIQUE cells lies inside a single SUBCLU cluster at
  radius $h sqrt(m + 3)$, whereas SUBCLU strictly avoids CLIQUE's
  boundary-alignment misses on shifted non-convex shapes.
]

On the ClusterOverview 4-dataset empirical benchmark, SUBCLU's
performance mirrors these theoretical properties quantitatively
(@rousseeuw1987silhouettes, @hubert1985comparing). Its standout
result occurs on *Multi-Shape* ($"Sil" = 0.584$, $"ARI" = "N/A"$),
by far the highest Silhouette score among the four subspace
algorithms (Algs. 137--140), directly reflecting its ability to
trace curved, non-convex density contours where axis-aligned box
methods (DOC: $0.312$) and linear-subspace projections (ORCLUS:
$0.411$) impose rigid geometric bias, and where adaptive grids
(MAFIA: $"N/A"$) fail to separate shapes. On *Iris*
($"Sil" = 0.315$, $"ARI" = 0.516$), SUBCLU outperforms both DOC
($"ARI" = 0.398$) and MAFIA ($"ARI" = 0.388$) in true label
recovery, though it trails ORCLUS ($"ARI" = 0.681$) because the
partially overlapping _Versicolor_ / _Virginica_ classes sit on
an oblique linear discriminant plane rather than axis-parallel
projections. On *Customers* ($"Sil" = 0.170$) and *High-D*
($"Sil" = 0.212$, $"ARI" = "N/A"$), fixed-$epsilon$ volume decay
(Proposition 1) limits cluster separation in higher subspaces
compared to ORCLUS ($0.318$ and $0.254$). Finally, SUBCLU's mean
runtime of $12.412 "ms"$ is the heaviest of the quartet—4.4x
slower than ORCLUS ($2.792 "ms"$), 111x slower than DOC
($0.112 "ms"$), and 214x slower than MAFIA ($0.058 "ms"$)—a
direct empirical signature of repeated quadratic DBSCAN passes
across the bottom-up subspace lattice.
