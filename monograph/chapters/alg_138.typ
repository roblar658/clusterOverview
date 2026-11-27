#import "../template/ieee.typ": *

#alg_entry(
  id: 138,
  name: "ORCLUS: Arbitrarily Oriented Projected Clustering",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2000",
  paradigm: "Hierarchical PCA-Projected Subspace Clustering",
  complexity: $cal(O)(k_0^3 + k_0 n d + k_0^2 d^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Axis-parallel subspace clustering methods—including CLIQUE
@agrawal1998clique, PROCLUS @aggarwal1999proclus, DOC
@procopiuc2002doc, and SUBCLU @kailing2004subclu—implicitly assume
that relevant cluster attributes align with the canonical
coordinate basis of $bb(R)^d$. In physical, econometric, and
biomedical datasets, however, strong inter-feature correlations
confine clusters to *oblique* low-dimensional affine subspaces
whose orientation is rotated relative to the measurement axes.
ORCLUS (Arbitrarily Oriented Projected Cluster Generation),
proposed by Aggarwal and Yu @aggarwal2000orclus, generalizes
PROCLUS by replacing canonical coordinate subsets with locally
adapted orthonormal bases learned via Principal Component
Analysis (PCA) inside an agglomerative hierarchical schedule.

Let $cal(X) = \{bold(x)_1, ..., bold(x)_n\} subset bb(R)^d$ be
a dataset, $k in bb(N)$ the target number of clusters, and
$l in \{1, ..., d\}$ the target subspace dimensionality. Let
$cal(V)_(d, l) = \{bold(E) in bb(R)^(d times l) : bold(E)^top
bold(E) = bold(I)_l\}$ denote the Stiefel manifold of
orthonormal $l$-frames in $bb(R)^d$. Each cluster $C_i$ carries
a centroid $bold(mu)_i = (1 / |C_i|) sum_(bold(x) in C_i) bold(x)$,
sample covariance matrix
$
  bold(Sigma)_i
  = (1 / |C_i|) sum_(bold(x) in C_i)
  (bold(x) - bold(mu)_i)(bold(x) - bold(mu)_i)^top,
$
and local orthonormal basis $bold(E)_i in cal(V)_(d, l)$ spanning
the $l$-dimensional projected subspace $cal(U)_i = "span"(bold(E)_i)$.

#definition(title: "Generalized Projected Cluster and Energy")[
  A _generalized projected cluster_ of target dimension $l$ is a
  pair $(C_i, bold(E)_i)$ where $C_i subset.eq cal(X)$ is a
  non-empty point subset and $bold(E)_i in cal(V)_(d, l)$ is an
  orthonormal $l$-frame. The _projected Euclidean distance_ of
  $bold(x) in bb(R)^d$ to cluster $(C_i, bold(E)_i)$ is
  $
    d(bold(x), bold(mu)_i ; bold(E)_i)
    = ||bold(E)_i^top (bold(x) - bold(mu)_i)||_2
    = (sum_(r=1)^l chevron.l bold(x) - bold(mu)_i,
      bold(e)_(i, r) chevron.r^2)^(1 / 2).
  $
  The _projected energy_ (mean squared projected spread) of
  $C_i$ along $bold(E)_i$ is
  $
    E(C_i, bold(E)_i)
    = (1 / |C_i|) sum_(bold(x) in C_i)
      ||bold(E)_i^top (bold(x) - bold(mu)_i)||_2^2
    = "Tr"(bold(E)_i^top bold(Sigma)_i bold(E)_i).
  $
]

A critical distinction separates ORCLUS from standard subspace
PCA (such as sparse subspace clustering @elhamifar2013ssc):
whereas representation-oriented PCA retains the *largest*-variance
eigenvectors, projected clustering seeks directions along which
cluster members are *tightly concentrated* while background
noise spreads along the complementary $d - l$ directions.
Accordingly, ORCLUS selects the $l$ orthonormal eigenvectors of
$bold(Sigma)_i$ corresponding to its *$l$ smallest eigenvalues*
$0 <= lambda_1 (bold(Sigma)_i) <= ... <= lambda_l (bold(Sigma)_i)$.
Because estimating $bold(E)_i$ directly in $l << d$ dimensions
from an initial poor partition traps local search in spurious
orientations, ORCLUS couples $k$-means-style alternating
assignment @lloyd1982least with top-down hierarchical merging
starting from $k_0 > k$ overspecified seeds in full dimension
$l_0 = d$, geometrically annealing both $(k_c, l_c)$ toward
$(k, l)$.

#algbox(
  title: "ORCLUS Iterative Hierarchical Subspace Reduction",
  inputs: [Dataset $cal(X) subset bb(R)^d$, target clusters $k$,
    target subspace dimension $l <= d$, initial seeds $k_0 > k$,
    decay factor $alpha in (0, 1)$ (default $alpha = 0.5$).],
  outputs: [Partition $\{C_1, ..., C_k\}$ (plus outliers $cal(O)$)
    and orthonormal bases $\{bold(E)_1, ..., bold(E)_k\}$.],
)[
  1. *Initialization*: Sample $k_c <- k_0$ seeds
     $\{bold(mu)_1, ..., bold(mu)_(k_0)\} subset cal(X)$; set
     $l_c <- d$ and $bold(E)_i <- bold(I)_d$ for all $i$. Set
     $beta <- exp((ln(l / d) dot.op ln alpha) / ln(k / k_0))$ so
     $k_c$ and $l_c$ reach $(k, l)$ in the same number of rounds.
  2. *Hierarchical Loop*: While $k_c > k$:
  3. *Assign*: Assign each $bold(x) in cal(X)$ to
     $i^* = arg min_(1 <= i <= k_c) ||bold(E)_i^top
     (bold(x) - bold(mu)_i)||_2$, forming clusters
     $C_1, ..., C_(k_c)$, and update centroids $bold(mu)_i$.
  4. *FindVectors*: For each cluster $C_i$, compute
     $bold(Sigma)_i$ and set $bold(E)_i in cal(V)_(d, l_c)$ to
     the $l_c$ orthonormal eigenvectors of $bold(Sigma)_i$ with
     smallest eigenvalues.
  5. *Schedule Update*: Set
     $k_"new" <- max(k, floor(alpha k_c))$ and
     $l_"new" <- max(l, floor(beta l_c))$.
  6. *Merge*: For every pair $1 <= i < j <= k_c$, form
     $C_(i j) = C_i union C_j$, compute covariance
     $bold(Sigma)_(i j)$ and minimum projected energy
     $r_(i j) = sum_(r=1)^(l_"new") lambda_r (bold(Sigma)_(i j))$.
     Greedily merge the pair $(i^*, j^*)$ minimizing $r_(i j)$
     and update affected pairs until $k_"new"$ clusters remain;
     set $k_c <- k_"new"$, $l_c <- l_"new"$.
  7. *Final Pass and Outlier Pruning*: Perform one final
     `FindVectors` and `Assign` pass at $(k, l)$; optionally flag
     $bold(x) in C_i$ as outlier if
     $d(bold(x), bold(mu)_i ; bold(E)_i) > min_(j != i)
     d(bold(mu)_i, bold(mu)_j ; bold(E)_i)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five formal mathematical results: Ky Fan spectral
optimality of the `FindVectors` step, an original proof of
orthogonal invariance coupled with an unbounded approximation gap
over axis-aligned PROCLUS, an original proof of fixed-$(k_c, l_c)$
coordinate-descent monotonicity, an exact algebraic merge-energy
decomposition, and a self-contained Davis-Kahan finite-sample
subspace stability bound.

#subsec("A", "Spectral Optimality of Minimum-Variance Subspaces")

#theorem(
  title: "Ky Fan Minimum-Energy Subspace Optimality",
  origin: "Classical — Cited",
)[
  Let $bold(Sigma)_i in bb(R)^(d times d)$ be the symmetric
  positive semidefinite covariance matrix of cluster $C_i$, with
  ordered eigenvalues
  $0 <= lambda_1 <= lambda_2 <= ... <= lambda_d$ and associated
  orthonormal eigenvectors $bold(v)_1, ..., bold(v)_d$
  @aggarwal2000orclus @bishop2006pattern. Then for every
  $l in \{1, ..., d\}$,
  $
    min_(bold(E) in cal(V)_(d, l)) E(C_i, bold(E))
    = min_(bold(E) in cal(V)_(d, l))
      "Tr"(bold(E)^top bold(Sigma)_i bold(E))
    = sum_(r=1)^l lambda_r,
  $
  attained by $bold(E)_i = [bold(v)_1, ..., bold(v)_l]$, and
  uniquely (up to right-multiplication by $O(l)$) whenever
  $lambda_l < lambda_(l+1)$.
]

#proof[
  Write the spectral decomposition
  $bold(Sigma)_i = bold(V) bold(Lambda) bold(V)^top$ with
  $bold(V) in O(d)$ and
  $bold(Lambda) = "diag"(lambda_1, ..., lambda_d)$. For any
  $bold(E) in cal(V)_(d, l)$, define
  $bold(U) = bold(V)^top bold(E) in cal(V)_(d, l)$, and set
  $w_t = sum_(r=1)^l U_(t, r)^2$ for $t = 1, ..., d$. Because the
  columns of $bold(U)$ are orthonormal, $bold(U)$ can be extended
  to an orthogonal matrix in $O(d)$, implying both
  $0 <= w_t <= 1$ for every $t$ and
  $sum_(t=1)^d w_t = "Tr"(bold(U)^top bold(U)) = l$.
  Expanding the projected energy yields
  $
    "Tr"(bold(E)^top bold(Sigma)_i bold(E))
    = "Tr"(bold(U)^top bold(Lambda) bold(U))
    = sum_(t=1)^d lambda_t w_t.
  $
  Now subtract $sum_(t=1)^l lambda_t$ using
  $sum_(t=1)^l (1 - w_t) = l - sum_(t=1)^l w_t = sum_(t=l+1)^d w_t$:
  $
    sum_(t=1)^d lambda_t w_t - sum_(t=1)^l lambda_t
    &= sum_(t=l+1)^d lambda_t w_t
      - sum_(t=1)^l lambda_t (1 - w_t) \
    &>= lambda_(l+1) sum_(t=l+1)^d w_t
      - lambda_l sum_(t=1)^l (1 - w_t) \
    &= (lambda_(l+1) - lambda_l) sum_(t=l+1)^d w_t
    >= 0,
  $
  since $w_t >= 0$, $1 - w_t >= 0$, and $lambda_(l+1) >= lambda_l$.
  Choosing $bold(E)_i = [bold(v)_1, ..., bold(v)_l]$ gives
  $w_1 = ... = w_l = 1$ and $w_(l+1) = ... = w_d = 0$, attaining
  the lower bound. When $lambda_l < lambda_(l+1)$, equality
  forces $w_t = 0$ for all $t > l$, i.e.,
  $"span"(bold(E)) = "span"(bold(v)_1, ..., bold(v)_l)$.
]

#subsec("B", "Rigid-Motion Invariance and Gap Over Axis-Aligned PROCLUS")

#theorem(
  title: "O(d) Equivariance and Unbounded Gap Over PROCLUS",
  origin: "Original Proof",
)[
  1. (_Exact Rigid-Motion Invariance_.) Under any Euclidean
     isometry $T(bold(x)) = bold(Q) bold(x) + bold(t)$ with
     $bold(Q) in O(d)$ and $bold(t) in bb(R)^d$, ORCLUS
     transforms centroids as $bold(mu)_i |-> bold(Q) bold(mu)_i
     + bold(t)$ and bases as $bold(E)_i |-> bold(Q) bold(E)_i$,
     leaving all projected distances, energies $E(C_i, bold(E)_i)$,
     merge costs $r_(i j)$, and cluster assignments identically
     invariant (given transformed initial seeds).
  2. (_Unbounded Energy Separation from PROCLUS_.) For every
     dimension $d >= 2$ and variance ratio
     $rho = sigma_"out"^2 / sigma_"in"^2 > 1$, there exists a
     cluster covariance $bold(Sigma) succ.eq bold(0)$ such that
     for $l = 1$ (and for every $1 <= l < d$ whenever $d = 2^m$),
     ORCLUS achieves optimal projected energy
     $E_"ORCLUS" = l sigma_"in"^2$, whereas *every* $l$-dimensional
     axis-aligned PROCLUS subspace incurs energy
     $
       E_"PROCLUS"
       = l ((l / d) sigma_"in"^2 + (1 - l / d) sigma_"out"^2),
     $
     yielding an unbounded energy-inflation factor
     $E_"PROCLUS" / E_"ORCLUS" = (l / d) + (1 - l / d) rho$.
]

#proof[
  (1) Under $bold(x)' = bold(Q) bold(x) + bold(t)$, the sample
  mean shifts to $bold(mu)'_i = bold(Q) bold(mu)_i + bold(t)$
  and the covariance transforms by orthogonal similarity:
  $bold(Sigma)'_i = bold(Q) bold(Sigma)_i bold(Q)^top$. Thus
  $bold(Sigma)'_i$ has identical eigenvalues to $bold(Sigma)_i$,
  and its bottom-$l$ eigenvectors are $bold(E)'_i = bold(Q) bold(E)_i$.
  Consequently,
  $
    ||bold(E)'_i^top (bold(x)' - bold(mu)'_i)||_2
    = ||bold(E)_i^top bold(Q)^top bold(Q)(bold(x) - bold(mu)_i)||_2
    = ||bold(E)_i^top (bold(x) - bold(mu)_i)||_2,
  $
  preserving every assignment, merge cost $r_(i j)$, and final
  partition identically.

  (2) First take $l = 1$ and arbitrary $d >= 2$. Let
  $bold(u) = d^(-1/2) (1, ..., 1)^top$ ($||bold(u)||_2 = 1$), and
  define the rank-one-deflated covariance
  $
    bold(Sigma)
    = sigma_"out"^2 bold(I)_d
    - (sigma_"out"^2 - sigma_"in"^2) bold(u) bold(u)^top.
  $
  Clearly $bold(Sigma) bold(u) = sigma_"in"^2 bold(u)$ while
  $bold(Sigma) bold(w) = sigma_"out"^2 bold(w)$ for all
  $bold(w) perp bold(u)$, so Theorem 1 gives
  $E_"ORCLUS" = lambda_1 (bold(Sigma)) = sigma_"in"^2$. For any
  canonical basis vector $bold(e)_j$, however,
  $
    bold(e)_j^top bold(Sigma) bold(e)_j
    = sigma_"out"^2 - (sigma_"out"^2 - sigma_"in"^2) / d
    = (1 / d) sigma_"in"^2 + (1 - 1 / d) sigma_"out"^2,
  $
  giving $E_"PROCLUS" / E_"ORCLUS" = (1/d) + (1 - 1/d) rho$.
  For general $1 <= l < d$ with $d = 2^m$, let
  $bold(H) in O(d)$ be the normalized Walsh-Hadamard matrix
  satisfying $H_(r c)^2 = 1 / d$ for all $r, c$, and set
  $bold(Sigma) = bold(H) bold(Lambda) bold(H)^top$ where
  $bold(Lambda)$ has $l$ entries equal to $sigma_"in"^2$ and
  $d - l$ entries equal to $sigma_"out"^2$. Then
  $E_"ORCLUS" = sum_(r=1)^l lambda_r = l sigma_"in"^2$, while
  every diagonal entry of $bold(Sigma)$ equals
  $
    Sigma_(j j)
    = sum_(c=1)^d H_(j c)^2 Lambda_(c c)
    = (1 / d) "Tr"(bold(Lambda))
    = (l / d) sigma_"in"^2 + (1 - l / d) sigma_"out"^2.
  $
  Summing over any $l$-subset of coordinate axes yields
  $E_"PROCLUS" = l Sigma_(1 1)$, establishing the claimed exact
  ratio, which diverges linearly with $rho = sigma_"out"^2 / sigma_"in"^2$.
]

#subsec("C", "Fixed-Dimension Monotonicity and Merge Decomposition")

#proposition(
  title: "Lyapunov Monotonicity of Fixed-(k, l) Alternating Passes",
  origin: "Original Derivation",
)[
  At fixed $(k_c, l_c)$, consider the global projected distortion
  functional
  $
    J(\{C_i, bold(mu)_i, bold(E)_i\})
    = sum_(i=1)^(k_c) sum_(bold(x) in C_i)
      ||bold(E)_i^top (bold(x) - bold(mu)_i)||_2^2.
  $
  Iterating (i) nearest-projected-centroid assignment, (ii) mean
  update $bold(mu)_i = (1 / |C_i|) sum_(bold(x) in C_i) bold(x)$,
  and (iii) Ky Fan basis update `FindVectors` decreases $J$
  monotonically at every sub-step and terminates in finitely many
  iterations.
]

#proof[
  Step (i) fixes $\{bold(mu)_i, bold(E)_i\}$ and assigns each
  $bold(x)$ independently to
  $arg min_i ||bold(E)_i^top (bold(x) - bold(mu)_i)||_2^2$, which
  cannot increase $J$. In step (ii), for fixed $C_i$ and
  $bold(E)_i$, expand the cluster objective around the sample
  mean $overline(bold(x))_i = (1 / |C_i|) sum_(bold(x) in C_i) bold(x)$:
  $
    sum_(bold(x) in C_i) ||bold(E)_i^top (bold(x) - bold(mu))||_2^2
    &= sum_(bold(x) in C_i)
      ||bold(E)_i^top (bold(x) - overline(bold(x))_i)||_2^2 \
    &+ |C_i| dot.op ||bold(E)_i^top (overline(bold(x))_i - bold(mu))||_2^2,
  $
  where the cross-term vanishes identically because
  $sum_(bold(x) in C_i) (bold(x) - overline(bold(x))_i) = bold(0)$.
  Setting $bold(mu)_i = overline(bold(x))_i$ sets the second term
  to zero, globally minimizing $J$ over $bold(mu)_i$. In step
  (iii), with $bold(mu)_i = overline(bold(x))_i$, the $i$-th
  summand equals $|C_i| E(C_i, bold(E)_i)$, which is globally
  minimized over $cal(V)_(d, l_c)$ by Theorem 1. Thus $J$ is
  non-increasing across all three sub-steps. Because there are at
  most $k_c^n$ distinct partitions and each partition maps to a
  unique optimal value of $J$, strict descent terminates in
  finitely many steps.
]

Remarkably, across ORCLUS's *hierarchical merge* steps, $J$ is
not monotone because $k_c$ shrinks (raising distortion) while
$l_c$ simultaneously shrinks (lowering distortion). The next
lemma dissects exactly what the greedy merge criterion rewards.

#lemma(
  title: "Parallel-Axis Covariance Identity and Merge Energy Lower Bound",
  origin: "Original Proof",
)[
  Let $C_i, C_j subset.eq cal(X)$ be disjoint clusters with sizes
  $n_i, n_j$, means $bold(mu)_i, bold(mu)_j$, covariances
  $bold(Sigma)_i, bold(Sigma)_j$, and union $C_(i j) = C_i union C_j$
  of size $n_(i j) = n_i + n_j$. Write
  $bold(Delta)_(i j) = bold(mu)_i - bold(mu)_j$. Then the merged
  covariance obeys the exact rank-one parallel-axis identity
  $
    bold(Sigma)_(i j)
    = (n_i / n_(i j)) bold(Sigma)_i
    + (n_j / n_(i j)) bold(Sigma)_j
    + (n_i n_j / n_(i j)^2) bold(Delta)_(i j) bold(Delta)_(i j)^top.
  $
  Consequently, for any dimension $l_"new"$, the optimal merged
  energy $r_(i j) = E(C_(i j), bold(E)_(i j))$ satisfies
  $
    r_(i j)
    >= (n_i / n_(i j)) E_(l_"new")(C_i)
    + (n_j / n_(i j)) E_(l_"new")(C_j)
    + (n_i n_j / n_(i j)^2)
      ||bold(E)_(i j)^top bold(Delta)_(i j)||_2^2,
  $
  where $E_(l_"new")(C_i) = sum_(r=1)^(l_"new") lambda_r (bold(Sigma)_i)$,
  with equality if and only if $bold(E)_(i j)$ simultaneously
  minimizes both constituent energies.
]

#proof[
  The merged mean is
  $bold(mu)_(i j) = (n_i bold(mu)_i + n_j bold(mu)_j) / n_(i j)$,
  so $bold(mu)_i - bold(mu)_(i j) = (n_j / n_(i j)) bold(Delta)_(i j)$
  and $bold(mu)_j - bold(mu)_(i j) = -(n_i / n_(i j)) bold(Delta)_(i j)$.
  Writing $bold(x) - bold(mu)_(i j) = (bold(x) - bold(mu)_i)
  + (bold(mu)_i - bold(mu)_(i j))$ for $bold(x) in C_i$, the
  cross-terms vanish on summation over $C_i$, leaving
  $
    sum_(bold(x) in C_i)
      (bold(x) - bold(mu)_(i j))(bold(x) - bold(mu)_(i j))^top
    = n_i bold(Sigma)_i
      + (n_i n_j^2 / n_(i j)^2)
        bold(Delta)_(i j) bold(Delta)_(i j)^top.
  $
  Adding the symmetric expression for $C_j$ (whose rank-one
  weight is $n_j n_i^2 / n_(i j)^2$) and dividing by $n_(i j)$
  gives $(n_i n_j (n_j + n_i)) / n_(i j)^3 = n_i n_j / n_(i j)^2$,
  proving the covariance identity. Taking
  $"Tr"(bold(E)_(i j)^top bold(Sigma)_(i j) bold(E)_(i j))$ and
  applying Theorem 1 to bound
  $"Tr"(bold(E)_(i j)^top bold(Sigma)_i bold(E)_(i j)) >= E_(l_"new")(C_i)$
  and similarly for $C_j$ yields the exact lower bound.
]

Lemma 1 reveals a clean geometric decomposition of ORCLUS's merge
objective: $r_(i j)$ penalizes (a) intrinsic spread of $C_i$ and
$C_j$, (b) subspace misalignment via the nonnegative gap
$"Tr"(bold(E)_(i j)^top bold(Sigma)_i bold(E)_(i j)) - E_(l_"new")(C_i)$,
and (c) centroid displacement projected onto the merged subspace
via the Ward-like term $(n_i n_j / n_(i j)^2) ||bold(E)_(i j)^top
bold(Delta)_(i j)||_2^2$ @ward1963hierarchical.

#subsec("D", "Finite-Sample Subspace Estimation Stability")

#theorem(
  title: "Davis-Kahan Subspace Recovery Bound for ORCLUS",
  origin: "Classical — Cited",
)[
  Suppose cluster $C_i$ is sampled from a population covariance
  $bold(Sigma)^*$ whose eigenvalues satisfy the separation gap
  $gamma = lambda_(l+1)(bold(Sigma)^*) - lambda_l (bold(Sigma)^*) > 0$,
  with true bottom-$l$ frame $bold(E)^* in cal(V)_(d, l)$ and
  empirical frame $hat(bold(E))_i in cal(V)_(d, l)$ computed from
  $hat(bold(Sigma))_i$. Let $bold(Delta) = hat(bold(Sigma))_i - bold(Sigma)^*$.
  Then the principal-angle sine matrix between $"span"(hat(bold(E))_i)$
  and $"span"(bold(E)^*)$ obeys the unconditional Davis-Kahan bound
  @davis1970rotation:
  $
    ||sin bold(Theta)(hat(bold(E))_i, bold(E)^*)||_F
    <= (2 sqrt(l) dot.op ||hat(bold(Sigma))_i - bold(Sigma)^*||_2)
       / gamma.
  $
]

#proof[
  Let $bold(E)_perp^* in cal(V)_(d, d - l)$ span the top $d - l$
  eigenvectors of $bold(Sigma)^*$ with diagonal eigenvalue matrix
  $bold(Lambda)_perp^* succ.eq lambda_(l+1)^* bold(I)_(d-l)$, and
  let $hat(bold(Lambda))_l prec.eq hat(lambda)_l bold(I)_l$ contain
  the bottom $l$ eigenvalues of $hat(bold(Sigma))_i$. Recall that
  $||sin bold(Theta)(hat(bold(E))_i, bold(E)^*)||_F
  = ||bold(E)_perp^{* top} hat(bold(E))_i||_F <= sqrt(l)$.
  If $||bold(Delta)||_2 >= gamma / 2$, the right-hand side is at
  least $sqrt(l)$ and the bound holds trivially.

  Assume therefore $||bold(Delta)||_2 < gamma / 2$. By Weyl's
  eigenvalue perturbation inequality,
  $hat(lambda)_l <= lambda_l^* + ||bold(Delta)||_2$, so the
  spectral separation satisfies
  $lambda_(l+1)^* - hat(lambda)_l >= gamma - ||bold(Delta)||_2
  > gamma / 2 > 0$. Using
  $bold(E)_perp^{* top} bold(Sigma)^* = bold(Lambda)_perp^*
  bold(E)_perp^{* top}$ and
  $hat(bold(Sigma))_i hat(bold(E))_i = hat(bold(E))_i
  hat(bold(Lambda))_l$, form the Sylvester identity:
  $
    bold(Lambda)_perp^* (bold(E)_perp^{* top} hat(bold(E))_i)
    - (bold(E)_perp^{* top} hat(bold(E))_i) hat(bold(Lambda))_l
    = -bold(E)_perp^{* top} bold(Delta) hat(bold(E))_i.
  $
  Denoting $bold(M) = bold(E)_perp^{* top} hat(bold(E))_i$, each
  $(r, s)$ entry of the left-hand side equals
  $(lambda_(l+r)^* - hat(lambda)_s) M_(r s)$ with
  $lambda_(l+r)^* - hat(lambda)_s >= lambda_(l+1)^* - hat(lambda)_l
  > gamma / 2$. Taking Frobenius norms on both sides yields
  $
    (gamma / 2) ||bold(M)||_F
    <= ||bold(E)_perp^{* top} bold(Delta) hat(bold(E))_i||_F
    <= ||bold(Delta)||_2 dot.op ||hat(bold(E))_i||_F
    = sqrt(l) ||bold(Delta)||_2,
  $
  since $||bold(E)_perp^{* top}||_2 = 1$ and
  $||hat(bold(E))_i||_F = sqrt(l)$. Dividing by $gamma / 2$
  completes the proof.
]

Theorem 3 immediately exposes ORCLUS's sample-size floor: when
$|C_i| <= d$, $hat(bold(Sigma))_i$ is rank-deficient and exhibits
at least $d - |C_i| + 1$ spurious zero eigenvalues regardless of
$bold(Sigma)^*$, destroying the eigengap and causing
`FindVectors` to overfit arbitrary noise hyperplanes.

#sec("III", "Complexity, Geometry and Operational Profile")

Across the hierarchical annealing schedule
$k_c <- max(k, floor(alpha k_c))$, the number of outer rounds is
$R = ceil(log_(1 / alpha)(k_0 / k))$. In round $t$ with $k_t$
clusters, the `Assign` step projects $n$ points onto $k_t$ bases
of dimension $l_t <= d$ in $cal(O)(k_t n d)$ time; `FindVectors`
forms and eigendecomposes $k_t$ covariance matrices of size
$d times d$ in $cal(O)(k_t n d + k_t d^3)$ time; and `Merge`
evaluates up to $binom(k_t, 2)$ pairwise merged covariances
(using the Lemma 1 rank-one formula in $cal(O)(d^2)$ per pair,
avoiding an $cal(O)(n)$ rescan) and eigendecomposes them in
$cal(O)(k_t^2 d^3)$ time, with $cal(O)(k_t^3)$ bookkeeping for
repeated merge selection @aggarwal2000orclus. Because $k_t$
decays geometrically ($sum_t k_t^p = cal(O)(k_0^p)$ for
$p in \{1, 2, 3\}$), total time complexity simplifies to
$cal(O)(k_0^3 + k_0 n d + k_0^2 d^3)$.

Memory complexity is $cal(O)(n d + k_0 d^2 + k_0^2)$ to hold the
dataset, active covariance matrices and eigenvectors, and the
pairwise merge cost table. Note the steep cubic dependence on
$d$: for moderate dimensionality ($d <= 100$) and small $k_0$,
ORCLUS is fast and cache-friendly (dense BLAS matrix products),
whereas for $d > 1000$, the $k_0^2 d^3$ merge eigendecompositions
become prohibitive unless preceded by randomized sketching or
truncated SVD.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity],
    [$cal(O)(k_0^3 + k_0 n d + k_0^2 d^3)$ via geometric decay
    of $k_c$; linear in $n$, cubic in $d$ and $k_0$.],
    [Space complexity],
    [$cal(O)(n d + k_0 d^2 + k_0^2)$ for data, cluster
    covariances, bases, and merge heap.],
    [Cluster geometry],
    [Arbitrarily oriented affine linear subspaces (rotated
    hyperplanes and correlated ellipsoids) of equal dim $l$.],
    [Noise & overlap],
    [Disjoint hard partition plus optional distance-threshold
    outlier set $cal(O)$; no cluster overlap.],
    [Invariance],
    [Strictly invariant under rigid motions $O(d) plus.o bb(R)^d$
    and uniform scaling; sensitive to anisotropic rescaling.],
    [Key hyperparameters],
    [Cluster count $k$, uniform target subspace dim $l$, initial
    seed count $k_0 > k$, decay rate $alpha$.],
  ),
  caption: [Operational, geometric, and complexity profile of ORCLUS.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multi-Asset Quantitative Factor and Yield-Curve Regimes*:
   Sovereign bond yields and equity factor returns move along
   oblique combinations (level, slope, curvature) rather than
   individual maturities. ORCLUS discovers regime clusters defined
   by distinct co-integration hyperplanes (where linear combinations
   of rates exhibit minimal variance), directly matching the
   minimum-eigenvalue inductive bias of Theorem 1.
2. *Hyperspectral Remote Sensing and Spectral Unmixing*: Each
   material class (minerals, vegetation canopies, urban surfaces)
   occupies a rotated low-dimensional affine manifold across
   narrowband reflectance channels. Because atmospheric scattering
   and sensor calibration rotate spectral signatures relative to
   wavelength axes, ORCLUS's $O(d)$ invariance (Theorem 2) cleanly
   outperforms axis-aligned subspace grids.
3. *Robotics Kinematics and Multi-Modal Sensor Calibration*:
   Joint-torque and inertial measurement streams obey nonlinear
   mechanical constraints that locally linearize into oblique
   low-variance constraint subspaces corresponding to distinct
   contact or grasp modes. Extracting minimum-variance normal
   frames isolates each physical contact manifold automatically.
4. *Chemometrics and Spectroscopic Process Control*: Near-infrared
   (NIR) and Raman spectra exhibit broad overlapping absorption
   bands where chemical concentrations govern linear subspace
   constraints across correlated wavenumbers. ORCLUS separates
   batch operating regimes defined by different stoichiometric
   invariants without requiring peak picking.
5. *Biomedical Morphometry and Clinical Phenotyping*: In
   multivariate biomarker panels (e.g., metabolic ratios, cardiac
   hemodynamics), physiological homeostasis constrains healthy and
   pathological subpopulations onto oblique linear trade-off
   surfaces. ORCLUS captures these correlated phenotypes even when
   every single marginal biomarker overlaps across groups.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ORCLUS fails predictably in four regimes: (i) _Small-cluster rank
deficiency_ ($|C_i| <= d$), where Theorem 3 breaks down and
spurious zero eigenvalues trap random hyperplanes; (ii)
_Heterogeneous subspace dimensionalities_, because ORCLUS forces
a uniform final dimension $l$ across all $k$ clusters; (iii)
_Non-linear or curved manifolds_ (concentric rings, moons), which
violate the linear covariance model; and (iv) _Ultra-high
dimensionality_ ($d >> 10^3$), where the $k_0^2 d^3$ merge cost
dominates runtime.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Correlated oblique affine clusters],
    [Optimal: $O(d)$-invariant Ky Fan energy minimization.],
    [Ensure $|C_i| >> d$ and tune $l$ via scree plot.],
    [Axis-aligned interpretable features],
    [Moderate: rotated eigenvectors mix all $d$ features.],
    [Use PROCLUS @aggarwal1999proclus or DOC @procopiuc2002doc
    for axis interpretability.],
    [Non-convex / curved shapes],
    [Poor: linear subspace distance cuts curved manifolds.],
    [Use SUBCLU @kailing2004subclu (Alg. 137) or spectral methods.],
    [Ultra-high-$d$ sparse data ($d > 10^3$)],
    [Poor: $cal(O)(k_0^2 d^3)$ eigendecomposition bottleneck.],
    [Use SSC @elhamifar2013ssc or random projections first.],
    [Variable cluster subspace dims],
    [Weak: enforces common $l$ across all $k$ clusters.],
    [Use correlation clustering or adaptive-rank mixtures.],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for ORCLUS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ORCLUS forms the rotated generalization of PROCLUS
@aggarwal1999proclus and bridges projected clustering with
subspace-segmentation frameworks such as Sparse Subspace
Clustering @elhamifar2013ssc @soltanolkotabi2012geometric and
hierarchical Ward clustering @ward1963hierarchical.

#proposition(
  title: "Exact Reduction of ORCLUS to PROCLUS and Full-Space k-Means",
  origin: "Original Proof",
)[
  1. Restricting the feasible set of `FindVectors` from the
     Stiefel manifold $cal(V)_(d, l)$ to the finite discrete
     subset $cal(V)_(d, l)^"axis" = \{[bold(e)_(j_1), ...,
     bold(e)_(j_l)] : 1 <= j_1 < ... < j_l <= d\}$ reduces the
     projected energy to the sum of the $l$ smallest coordinate
     variances $sum_(t=1)^l Sigma_(i, j_t j_t)$, recovering
     PROCLUS's axis-selection criterion @aggarwal1999proclus
     under Euclidean spread.
  2. Setting $l = d$ fixes $bold(E)_i bold(E)_i^top = bold(I)_d$,
     collapsing ORCLUS's `Assign` pass to standard Lloyd $k$-means
     @lloyd1982least and reducing the merge energy $r_(i j)$ of
     Lemma 1 to the trace-normalized Ward variance criterion
     $(1 / n_(i j)) sum_(bold(x) in C_(i j)) ||bold(x) - bold(mu)_(i j)||_2^2$.
]

#proof[
  (1) For $bold(E) = [bold(e)_(j_1), ..., bold(e)_(j_l)]$,
  $"Tr"(bold(E)^top bold(Sigma)_i bold(E)) = sum_(t=1)^l
  bold(e)_(j_t)^top bold(Sigma)_i bold(e)_(j_t) = sum_(t=1)^l
  Sigma_(i, j_t j_t)$, minimized by choosing the $l$ coordinates
  with smallest within-cluster variance—the $L_2$ counterpart of
  PROCLUS's $L_1$ Manhattan segmental selection.
  (2) When $l = d$, every $bold(E)_i in O(d)$ obeys
  $||bold(E)_i^top (bold(x) - bold(mu)_i)||_2 = ||bold(x) - bold(mu)_i||_2$,
  recovering Euclidean $k$-means assignment identically, while
  $r_(i j) = "Tr"(bold(Sigma)_(i j))$ equals the mean squared
  within-cluster scatter of $C_(i j)$.
]

Empirical results on the ClusterOverview 4-dataset benchmark
(@rousseeuw1987silhouettes, @hubert1985comparing) vividly confirm
ORCLUS's geometric advantages over its axis-aligned peers. On
*Iris*, ORCLUS achieves $"ARI" = 0.681$ ($"Sil" = 0.292$), the
highest ground-truth recovery among Algs. 137--140 by a wide
margin (beating SUBCLU $0.516$, DOC $0.398$, and MAFIA $0.388$):
because _Iris_ petal length and petal width are strongly
correlated along an oblique axis, ORCLUS's rotated PCA frames
(Theorem 2) separate _Versicolor_ and _Virginica_ cleanly even
while full-space Euclidean Silhouette remains moderate. ORCLUS
likewise leads the quartet on both *Customers* ($"Sil" = 0.318$
vs. DOC $0.214$, SUBCLU $0.170$, MAFIA $"N/A"$) and *High-D*
($"Sil" = 0.254$ vs. SUBCLU $0.212$, MAFIA $0.177$, DOC $0.068$),
reflecting the stability of covariance projections in moderate
dimensionality when $n_i > d$. Conversely, on *Multi-Shape*
($"Sil" = 0.411$, $"ARI" = "N/A"$), ORCLUS trails density-based
SUBCLU ($0.584$) because linear subspaces cannot bend around
non-convex contours. Its mean runtime of $2.792 "ms"$ is 4.4x
faster than SUBCLU ($12.412 "ms"$), while trailing lightweight
sampling/histogram algorithms DOC ($0.112 "ms"$) and MAFIA
($0.058 "ms"$) due to repeated $d times d$ eigendecompositions
inside the hierarchical merge loop.
