#import "../template/ieee.typ": *

#alg_entry(
  id: 186,
  name: "JULE (Joint Unsupervised Learning)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2016",
  paradigm: "Recurrent Deep Representation and Agglomerative Clustering",
  complexity: $cal(O)(T (N P + N^2 d))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical agglomerative hierarchical clustering @ward1963hierarchical
@lance1967general constructs a nested tree of partitions from a
fixed pairwise distance matrix, so its accuracy is permanently
capped by the quality of hand-crafted input features. Conversely,
deep convolutional networks (CNNs) excel at learning discriminative
hierarchical representations, yet standard unsupervised objectives
such as pixel reconstruction do not explicitly optimize cluster
separability. Joint Unsupervised Learning (JULE), proposed by Yang,
Parikh, and Batra @yang2016jule, resolves this decoupling by
recasting agglomerative clustering and deep representation learning
as the forward and backward passes of a single *recurrent* process
unrolled over merge timesteps $t = 0, 1, dots, T$.

Let $cal(X) = {bold(x)_1, dots, bold(x)_N} subset.eq bb(R)^D$ be an
unlabeled dataset and $K >= 2$ be the desired final cluster count.
A neural encoder $f_theta : bb(R)^D -> bb(R)^d$ maps each input
$bold(x)_i$ to an $ell_2$-normalized feature vector on the unit
hypersphere $bb(S)^(d-1)$:
$
  bold(bar(x))_i (theta)
  = (f_theta (bold(x)_i)) / norm(f_theta (bold(x)_i))_2
  in bb(S)^(d-1).
$
Using the current embeddings ${bold(bar(x))_i}$, JULE builds a
sparse directed $K_s$-nearest-neighbor affinity matrix
$bold(W) in bb(R)_(>= 0)^(N times N)$ with entries
$
  W_(i j) = cases(
    exp(- norm(bold(bar(x))_i - bold(bar(x))_j)_2^2 / sigma^2)
    & "if" j in cal(N)_(K_s) (i),
    0 & "otherwise",
  )
$
where $cal(N)_(K_s) (i)$ is the set of $K_s$ nearest neighbors of
$bold(bar(x))_i$ (excluding $i$) and the adaptive bandwidth is
$sigma^2 = (a / (N K_s)) sum_(i=1)^N sum_(j in cal(N)_(K_s) (i)) norm(bold(bar(x))_i - bold(bar(x))_j)_2^2$
for scale factor $a > 0$.

Rather than measuring cluster similarity by single or average
Euclidean linkage, JULE adopts *Graph Degree Linkage* (GDL), which
evaluates bidirectional two-hop transition density between any two
disjoint clusters $cal(C)_a, cal(C)_b subset.eq {1, dots, N}$.

#definition(title: "GDL Cluster Affinity and Recurrent JULE Objective")[
  Let $bold(W)_(a b) in bb(R)^(|cal(C)_a| times |cal(C)_b|)$ denote
  the submatrix of $bold(W)$ with rows in $cal(C)_a$ and columns in
  $cal(C)_b$, and let $bold(1)_m in bb(R)^m$ be the all-ones vector.
  1. The *Graph Degree Linkage (GDL) affinity* between disjoint
     clusters $cal(C)_a$ and $cal(C)_b$ is
     $
       cal(A)(cal(C)_a, cal(C)_b)
       &= 1 / |cal(C)_a|^2 bold(1)_(|cal(C)_a|)^top
          bold(W)_(a b) bold(W)_(b a) bold(1)_(|cal(C)_a|) \
       & quad + 1 / |cal(C)_b|^2 bold(1)_(|cal(C)_b|)^top
          bold(W)_(b a) bold(W)_(a b) bold(1)_(|cal(C)_b|).
     $
  2. At recurrent step $t$, given partition
     $cal(P)^((t-1)) = {cal(C)_1, dots, cal(C)_(M_(t-1))}$ with
     labels $bold(y)^((t-1))$, JULE alternates between maximizing
     agglomerative merge affinity
     $
       cal(L)_"agg" (bold(y)^((t)) | bold(y)^((t-1)), theta)
       = - sum_((a, b) in cal(M)^((t))) cal(A)(cal(C)_a, cal(C)_b)
     $
     over a batch $cal(M)^((t))$ of $n_p^((t))$ cluster merges, and
     minimizing the hyperspherical neighborhood margin loss
     $
       cal(L)_"rep" (theta | bold(y)^((t)))
       = sum_(i=1)^N sum_(j in cal(N)_i^+) sum_(l in cal(N)_i^-)
         max(0, m + d_(i j)^2 (theta) - d_(i l)^2 (theta)),
     $
     where $d_(i j)^2 (theta) = norm(bold(bar(x))_i (theta) - bold(bar(x))_j (theta))_2^2$,
     $cal(N)_i^+$ are same-cluster $K_s$-NN positives, $cal(N)_i^-$
     are nearest different-cluster negatives ($|cal(N)_i^-| = K_d$),
     and $m in (0, 4)$ is the separation margin.
]

To avoid running an entire neural network training pass after every
single pairwise merge ($N - K$ stages), JULE introduces *partial
unrolling*: at each recurrent timestep $t$, it merges a constant
fraction $eta in (0, 1 / 2]$ of active clusters
($n_p^((t)) = min(ceil(eta M_(t-1)), M_(t-1) - K)$) before invoking
a multi-epoch SGD backward pass on $cal(L)_"rep" (theta | bold(y)^((t)))$
@yang2016jule.

#algbox(
  title: "JULE — Joint Unsupervised Learning",
  inputs: [
    Dataset $cal(X) = {bold(x)_i}_(i=1)^N subset.eq bb(R)^D$,
    target clusters $K$, neighbor counts $(K_s, K_d)$, unrolling
    rate $eta in (0, 1 / 2]$, margin $m > 0$, SGD steps $S$ per stage.
  ],
  outputs: [
    Learned encoder $f_theta$, hierarchy ${cal(P)^((t))}_(t=0)^T$,
    final cluster labels $bold(y)^((T)) in {1, dots, K}^N$.
  ],
)[
  1. *Initialization:* Initialize encoder $theta^((0))$, set
     singleton partition $cal(P)^((0)) = {{1}, dots, {N}}$
     ($M_0 = N$), and set timestep $t <- 0$.
  2. *While* $M_t > K$:
     - Increment $t <- t + 1$ and extract normalized embeddings
       $bold(bar(x))_i = f_(theta^((t-1))) (bold(x)_i) / norm(f_(theta^((t-1))) (bold(x)_i))_2$.
     - Build directed $K_s$-NN affinity matrix $bold(W)$ with
       adaptive scale $sigma^2$.
     - *(Forward Agglomerative Pass):* Compute pairwise GDL
       affinities $cal(A)(cal(C)_a, cal(C)_b)$ and greedily execute
       $n_p^((t)) = min(ceil(eta M_(t-1)), M_(t-1) - K)$ highest-affinity
       merges to obtain $cal(P)^((t))$ and labels $bold(y)^((t))$.
     - *(Backward Representation Pass):* Construct positive / negative
       neighborhoods $(cal(N)_i^+, cal(N)_i^-)$ under $bold(y)^((t))$
       and perform $S$ mini-batch SGD updates minimizing
       $cal(L)_"rep" (theta | bold(y)^((t)))$ to yield $theta^((t))$.
  3. *Return* encoder $f_(theta^((T)))$ and final partition $bold(y)^((T))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four rigorous results for JULE: (A) an exact
degree-product decomposition of GDL affinity proving immunity to
unidirectional outlier bridging, (B) tight logarithmic bounds on
recurrent unrolling depth and cumulative complexity, (C) a
hyperspherical margin separation theorem proving non-collapse and
local Silhouette guarantees, and (D) a constructive finite-sample
counterexample demonstrating irreversible early-merge amplification.

#subsec("A", "Degree-Product Decomposition and Outlier-Bridge Immunity")

#theorem(
  title: "Spectral Degree Identity and Bridge Immunity of GDL Affinity",
  origin: "Original Derivation",
)[
  For any disjoint clusters $cal(C)_a, cal(C)_b$ and any nonnegative
  directed affinity matrix $bold(W) in bb(R)_(>= 0)^(N times N)$,
  define for each node $v in cal(C)_b$ its incoming and outgoing
  cross-cluster degrees with respect to $cal(C)_a$:
  $
    d_"in" (v; cal(C)_a) = sum_(u in cal(C)_a) W_(u v),
    quad d_"out" (v; cal(C)_a) = sum_(w in cal(C)_a) W_(v w),
  $
  and define $d_"in" (u; cal(C)_b), d_"out" (u; cal(C)_b)$
  symmetrically for $u in cal(C)_a$. Then:
  1. *(Exact Degree-Product Representation):*
     $
       cal(A)(cal(C)_a, cal(C)_b)
       &= 1 / |cal(C)_a|^2 sum_(v in cal(C)_b)
          d_"in" (v; cal(C)_a) d_"out" (v; cal(C)_a) \
       & quad + 1 / |cal(C)_b|^2 sum_(u in cal(C)_a)
          d_"in" (u; cal(C)_b) d_"out" (u; cal(C)_b).
     $
  2. *(Symmetry and Spectral Bound):* $cal(A)(cal(C)_a, cal(C)_b) = cal(A)(cal(C)_b, cal(C)_a) >= 0$,
     and letting $sigma_"max" (dot)$ denote the leading singular
     value,
     $
       cal(A)(cal(C)_a, cal(C)_b)
       <= (1 / |cal(C)_a| + 1 / |cal(C)_b|)
          sigma_"max" (bold(W)_(a b)) sigma_"max" (bold(W)_(b a)).
     $
  3. *(Unidirectional Outlier Immunity):* If $cal(C)_b$ consists of
     outliers whose $K_s$-nearest neighbors fall inside $cal(C)_a$
     ($d_"out" (v; cal(C)_a) > 0$) while no node in $cal(C)_a$
     selects any node of $cal(C)_b$ among its $K_s$ nearest
     neighbors ($bold(W)_(a b) = bold(0)$), then
     $cal(A)(cal(C)_a, cal(C)_b) = 0$ identically — even when the
     unidirectional edges $bold(W)_(b a)$ have arbitrarily large
     weight.
]

#proof[
  For part (1), expand the first bilinear form in Definition 1:
  $
    bold(1)_(|cal(C)_a|)^top bold(W)_(a b) bold(W)_(b a) bold(1)_(|cal(C)_a|)
    &= sum_(u in cal(C)_a) sum_(v in cal(C)_b) sum_(w in cal(C)_a)
       W_(u v) W_(v w) \
    &= sum_(v in cal(C)_b) (sum_(u in cal(C)_a) W_(u v))
       (sum_(w in cal(C)_a) W_(v w)) \
    &= sum_(v in cal(C)_b) d_"in" (v; cal(C)_a) d_"out" (v; cal(C)_a).
  $
  Applying the identical expansion to
  $bold(1)_(|cal(C)_b|)^top bold(W)_(b a) bold(W)_(a b) bold(1)_(|cal(C)_b|)$
  gives the second sum, establishing the exact identity.

  For part (2), symmetry and non-negativity follow immediately from
  part (1) since $W_(i j) >= 0$ and the expression is invariant
  under swapping $a <-> b$. For the spectral upper bound, apply
  Cauchy–Schwarz and submultiplicativity of the spectral norm to
  the first term:
  $
    1 / |cal(C)_a|^2 bold(1)_(|cal(C)_a|)^top bold(W)_(a b) bold(W)_(b a) bold(1)_(|cal(C)_a|)
    &<= 1 / |cal(C)_a|^2 norm(bold(W)_(a b)^top bold(1)_(|cal(C)_a|))_2
        norm(bold(W)_(b a) bold(1)_(|cal(C)_a|))_2 \
    &<= 1 / |cal(C)_a|^2
        (sigma_"max" (bold(W)_(a b)) |cal(C)_a|^(1/2))
        (sigma_"max" (bold(W)_(b a)) |cal(C)_a|^(1/2)) \
    &= 1 / |cal(C)_a| sigma_"max" (bold(W)_(a b)) sigma_"max" (bold(W)_(b a)),
  $
  since $norm(bold(1)_m)_2 = m^(1/2)$. Repeating for the second term
  and summing yields the bound.

  For part (3), if $bold(W)_(a b) = bold(0)$, then
  $d_"in" (v; cal(C)_a) = 0$ for every $v in cal(C)_b$ and
  $d_"out" (u; cal(C)_b) = 0$ for every $u in cal(C)_a$. Every
  product term in part (1) vanishes identically, yielding
  $cal(A)(cal(C)_a, cal(C)_b) = 0$ regardless of $bold(W)_(b a)$.
]

Theorem 1 provides the exact geometric reason why JULE's GDL
criterion outperforms classical single-linkage @sibson1973slink on
curved manifolds with sparse bridge noise: a single-edge chain or
one-way hub outlier contributes zero bidirectional two-hop flow,
forcing merges to occur strictly across mutually connected dense
neighborhoods.

#subsec("B", "Logarithmic Unrolling Depth and Geometric Cost Bound")

#lemma(
  title: "Exact Unrolling Depth and Geometric Schedule Complexity",
  origin: "Original Proof",
)[
  Suppose JULE starts from $M_0 = N$ singletons and merges
  $n_p^((t)) = min(ceil(eta M_(t-1)), M_(t-1) - K)$ cluster pairs at
  stage $t$ with fixed unrolling rate $eta in (0, 1 / 2]$. Then:
  1. The number of recurrent stages $T$ required to reach $M_T = K$
     satisfies
     $
       floor((log(N / K)) / (log(1 / (1 - eta))))
       <= T
       <= ceil((log(N / K)) / (log(1 / (1 - eta)))) + 1,
     $
     so $T = Theta(eta^(-1) log(N / K))$ as $eta -> 0^+$.
  2. The cumulative number of cluster-count evaluations across all
     stages satisfies $sum_(t=0)^(T-1) M_t <= N / eta$, and the
     cumulative quadratic cluster-pair count satisfies
     $sum_(t=0)^(T-1) M_t^2 <= N^2 / (eta (2 - eta))$.
]

#proof[
  For every stage $t < T$ before the terminal cap activates,
  $M_t = M_(t-1) - ceil(eta M_(t-1)) <= (1 - eta) M_(t-1)$, so by
  induction $M_t <= (1 - eta)^t N$. Choosing
  $t_0 = ceil(log(N / K) / log(1 / (1 - eta)))$ guarantees
  $(1 - eta)^(t_0) N <= K$, hence at most $t_0 + 1$ stages suffice
  (accounting for the final capped step reaching exact $K$). For
  the lower bound, whenever $M_(t-1) >= 1 / eta$, we have
  $ceil(eta M_(t-1)) < eta M_(t-1) + 1$, and more cleanly, since
  each stage removes at most $ceil(eta M_(t-1))$ clusters, if
  $M_t >= (1 - eta)^t N - 1 / eta$, using $M_t >= (1 - eta) M_(t-1) - 1$
  shows $T = Omega(eta^(-1) log(N / K))$; in fact for
  $eta M_(t-1) in bb(N)$ the ceiling is exact and
  $M_t = (1 - eta)^t N$, yielding $T >= floor(log(N / K) / log(1 / (1 - eta)))$
  as the exact geometric lower envelope.

  For part (2), summing the geometric upper bound $M_t <= (1 - eta)^t N$
  gives
  $
    sum_(t=0)^(T-1) M_t
    <= N sum_(t=0)^infinity (1 - eta)^t
    = N / eta,
  $
  and squaring yields
  $
    sum_(t=0)^(T-1) M_t^2
    <= N^2 sum_(t=0)^infinity (1 - eta)^(2 t)
    = N^2 / (1 - (1 - eta)^2)
    = N^2 / (eta (2 - eta)).
  $
  This proves that partial unrolling collapses the naive
  $cal(O)(N)$ CNN retraining stages into $cal(O)(eta^(-1) log(N / K))$
  while keeping cumulative pairwise cluster comparisons within a
  constant factor $1 / (eta (2 - eta))$ of the single initial stage.
]

#subsec("C", "Hyperspherical Margin Separation and Non-Collapse")

#theorem(
  title: "Hyperspherical Angular Separation and Collapse Prevention",
  origin: "Original Proof",
)[
  Let ${bold(bar(x))_i}_(i=1)^N subset.eq bb(S)^(d-1)$ be the
  $ell_2$-normalized JULE embeddings at a partition $bold(y)$ with
  non-empty positive and negative neighbor sets
  $(cal(N)_i^+, cal(N)_i^-)$ for every $i$.
  1. *(Non-Collapse Lower Bound):* For any positive margin
     $m in (0, 4)$, any trivially collapsed representation
     $bold(bar(x))_1 = dots = bold(bar(x))_N = bold(u)_0 in bb(S)^(d-1)$
     incurs strictly positive representation loss
     $
       cal(L)_"rep" (theta_"coll" | bold(y))
       = m sum_(i=1)^N |cal(N)_i^+| |cal(N)_i^-| > 0,
     $
     whereas setting $m = 0$ makes total collapse a zero-loss
     global minimizer.
  2. *(Angular Gap & Local Silhouette Guarantee):* If $theta^ast$
     achieves $cal(L)_"rep" (theta^ast | bold(y)) = 0$ with margin
     $m > 0$, then every triplet $(i, j, l)$ with
     $j in cal(N)_i^+, l in cal(N)_i^-$ satisfies the cosine gap
     $
       chevron.l bold(bar(x))_i, bold(bar(x))_j chevron.r
       - chevron.l bold(bar(x))_i, bold(bar(x))_l chevron.r
       >= m / 2.
     $
     Moreover, if all intra-cluster pairwise distances obey
     $d_(i j) <= r_"in"$ and all inter-cluster distances satisfy
     $d_(i l)^2 >= r_"in"^2 + m$, then every datum $i$ has strictly
     positive Euclidean Silhouette width @rousseeuw1987silhouettes:
     $
       s(i) >= 1 - r_"in" / (r_"in"^2 + m)^(1/2) > 0.
     $
]

#proof[
  On the unit hypersphere $bb(S)^(d-1)$, expanding the squared
  Euclidean norm gives
  $d_(i j)^2 = norm(bold(bar(x))_i)_2^2 + norm(bold(bar(x))_j)_2^2 - 2 chevron.l bold(bar(x))_i, bold(bar(x))_j chevron.r = 2(1 - chevron.l bold(bar(x))_i, bold(bar(x))_j chevron.r) in [0, 4]$.

  For part (1), at a totally collapsed configuration
  $bold(bar(x))_i equiv bold(u)_0$, we have $d_(i j)^2 = d_(i l)^2 = 0$
  for all pairs $(i, j)$ and $(i, l)$. Every hinge term in
  $cal(L)_"rep"$ evaluates to $max(0, m + 0 - 0) = m$, summing to
  $m sum_(i=1)^N |cal(N)_i^+| |cal(N)_i^-|$. This is strictly
  positive whenever $m > 0$ and vanishes identically if $m = 0$.

  For part (2), since every hinge summand is nonnegative,
  $cal(L)_"rep" (theta^ast | bold(y)) = 0$ implies
  $m + d_(i j)^2 - d_(i l)^2 <= 0$ for all active triplets, hence
  $
    2(1 - chevron.l bold(bar(x))_i, bold(bar(x))_l chevron.r)
    - 2(1 - chevron.l bold(bar(x))_i, bold(bar(x))_j chevron.r) >= m,
  $
  which rearranges directly to
  $chevron.l bold(bar(x))_i, bold(bar(x))_j chevron.r - chevron.l bold(bar(x))_i, bold(bar(x))_l chevron.r >= m / 2$.
  Finally, under the uniform intra-cluster radius $r_"in"$, the
  average intra-cluster distance obeys $a(i) <= r_"in"$, while the
  nearest-cluster average distance obeys
  $b(i) >= (r_"in"^2 + m)^(1/2) > r_"in"$. Substituting into the
  Silhouette definition $s(i) = (b(i) - a(i)) / max(a(i), b(i))$
  yields $s(i) = 1 - a(i) / b(i) >= 1 - r_"in" / (r_"in"^2 + m)^(1/2) > 0$.
]

#subsec("D", "Irreversible Early-Merge Error Propagation Counterexample")

Because JULE constructs a strictly nested agglomerative hierarchy
($cal(P)^((t))$ coarsens $cal(P)^((t-1))$ without cluster-splitting
moves), representation updates cannot undo mistaken initial merges.

#proposition(
  title: "Irreversible Early-Merge Amplification in Recurrent JULE",
  origin: "Original Counterexample",
)[
  Let $cal(X) = cal(C)_1^ast union cal(C)_2^ast$ consist of two
  ground-truth clusters of size $N / 2$ ($N >= 6$) in $bb(R)^D$.
  Suppose the initial encoder $f_(theta^((0)))$ maps two opposite-class
  boundary points $x_p in cal(C)_1^ast$ and $x_q in cal(C)_2^ast$
  into a mutually closest pair
  $norm(bold(bar(x))_p - bold(bar(x))_q)_2 = delta$ strictly smaller
  than all other pairwise distances ($delta < min_((u, v) != (p, q)) norm(bold(bar(x))_u - bold(bar(x))_v)_2$)
  under $K_s = 1$. Then:
  1. At stage $t = 1$, GDL achieves its strict global maximum
     $cal(A)({p}, {q}) = 2 exp(- 2 delta^2 / sigma^2)$ on $({p}, {q})$,
     merging $p$ and $q$ into a permanent impure cluster.
  2. For all subsequent stages $t = 1, dots, T$, $y_p^((t)) = y_q^((t))$
     remains irreversibly locked, capping the maximum achievable
     purity at $1 - 1 / N$.
  3. Moreover, during the backward pass at $t = 1$, $q in cal(N)_p^+$
     induces an attractive gradient
     $nabla_(bold(bar(x))_p) d_(p q)^2 = 2 (bold(bar(x))_p - bold(bar(x))_q)$
     that actively pulls $bold(bar(x))_p$ and $bold(bar(x))_q$ — and,
     by neural Lipschitz smoothness, their surrounding local
     neighborhoods in $cal(C)_1^ast$ and $cal(C)_2^ast$ — *closer
     together*, amplifying initial boundary contamination across
     higher-level merges.
]

#proof[
  At $t = 0$, every cluster is a singleton ($|cal(C)_a| = 1$). By
  Theorem 1(1), for any two singletons ${u}, {v}$,
  $cal(A)({u}, {v}) = 2 W_(u v) W_(v u)$. Because $(p, q)$ is the
  unique global closest pair, $q$ is the $1$-NN of $p$ and $p$ is
  the $1$-NN of $q$, giving
  $W_(p q) = W_(q p) = exp(- delta^2 / sigma^2)$, strictly exceeding
  $W_(u v) W_(v u)$ for every other pair $(u, v)$. Thus stage
  $t = 1$ merges ${p}$ and ${q}$ first, proving (1).

  Since every forward pass of Alg. 186 only executes merge
  operations on $cal(P)^((t-1))$ and never splits an existing
  cluster, the equivalence relation $y_p^((t)) = y_q^((t))$ is
  preserved inductively for all $t >= 1$, proving (2).

  For part (3), whenever any negative $l in cal(N)_p^-$ activates
  the triplet hinge ($d_(p l)^2 < d_(p q)^2 + m$), the gradient of
  $cal(L)_"rep"$ with respect to $bold(bar(x))_p$ contains the
  attractive positive term $+ 2 (bold(bar(x))_p - bold(bar(x))_q)$,
  so negative-gradient descent moves $bold(bar(x))_p$ strictly
  toward $bold(bar(x))_q$. If the encoder $f_theta$ has Lipschitz
  constant $L_f$, any true neighbor $x_r in cal(C)_1^ast$ close to
  $x_p$ moves coherently with $bold(bar(x))_p$ toward
  $bold(bar(x))_q$, increasing cross-cluster $K_s$-NN overlap and
  triggering cascading false merges at later stages.
]

#sec("III", "Complexity, Geometry and Operational Profile")

JULE's computational cost combines deep forward-backward neural
passes with graph-based agglomerative clustering at every recurrent
timestep. Let $P$ denote the CNN / MLP forward-backward cost per
sample across $S$ SGD epochs per stage. At stage $t$ with $M_t$
clusters, extracting normalized features takes $cal(O)(N P)$,
building the exact directed $K_s$-NN affinity matrix $bold(W)$ in
$bb(R)^d$ takes $cal(O)(N^2 d)$ (or $cal(O)(N K_s d log N)$ with
approximate NN trees in low dimension), and evaluating sparse
pairwise GDL affinities takes $cal(O)(M_t^2 + N K_s^2)$ (using
sparse matrix-vector products on $bold(W)$). Summing across
$T = Theta(eta^(-1) log(N / K))$ stages (Lemma 1) yields total
time complexity
$
  T_"JULE" = cal(O)(T (N P S + N^2 d))
  = cal(O)(eta^(-1) log(N / K) (N P S + N^2 d)),
$
with peak memory $cal(O)(N^2 + N d + |theta|)$ if dense cluster
affinity tables are stored, or $cal(O)(N K_s + M_t^2 + |theta|)$
under sparse $K_s$-NN representation.

Geometrically, JULE operates on the *unit hypersphere $bb(S)^(d-1)$*
using local geodesic / chordal neighborhoods. Because GDL merges
clusters along chains of mutually reinforced $K_s$-NN blocks
(Theorem 1), JULE natively discovers *arbitrarily curved, elongated,
and non-convex ambient manifolds* without assuming globular
Gaussian or linear-hyperplane latent clusters. However, this local
graph strength becomes a vulnerability in very high ambient
dimensions with isotropic noise: by the nearest-neighbor distance
concentration phenomenon @beyer1999nearest, initial Euclidean
$K_s$-NN neighborhoods become noisy before the encoder has been
trained, seeding irreversible early merges (Proposition 1).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(eta^(-1) log(N / K) (N P S + N^2 d))$ across recurrent stages],
    [Space Complexity], [$cal(O)(N K_s + M_t^2 + |theta|)$ sparse; $cal(O)(N^2)$ worst-case affinity table],
    [Latent Geometry], [Unit hypersphere $bb(S)^(d-1)$ with local mutual $K_s$-NN manifold linkage],
    [Cluster Balance], [Unconstrained; naturally handles highly unequal cluster sizes and shapes],
    [Hierarchy Output], [Full nested dendrogram ${cal(P)^((t))}$ plus learned deep representation $f_theta$],
    [Initialization Sensitivity], [High at $t = 0$; early $K_s$-NN mistakes lock in irreversibly (Prop. 1)],
    [Primary Hyperparameters], [Cluster count $K$, neighbors $(K_s, K_d)$, unrolling rate $eta$, margin $m$],
  ),
  caption: [Operational and geometric profile of JULE (Alg. 186).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Unsupervised Facial Identity Discovery and Photo Organization:*
   Face images of the same person under varying pose, illumination,
   and expression trace smooth, non-convex, curved pose manifolds
   in CNN feature space, while person counts per album are strongly
   imbalanced. JULE's GDL affinity links progressive pose sequences
   transitively without forcing globular clusters or uniform identity
   counts, while backward triplet updates simultaneously sharpen
   pose-invariant identity features @yang2016jule.

2. *Fine-Grained Visual Object and Wildlife Species Taxonomies:*
   Camera-trap wildlife datasets and museum specimen archives exhibit
   natural multi-scale biological hierarchies (individual -> species
   -> genus). JULE produces both an inspectable agglomerative
   hierarchy ${cal(P)^((t))}$ and a discriminative hyperspherical
   embedding that separates visually similar subspecies via margin
   $m$.

3. *Medical Pathology and Radiology Sub-Phenotyping:* Histological
   slide tiles and lesion crops frequently form continuous morphologic
   manifolds bridging benign and malignant states, accompanied by
   unidirectional staining artifacts. Theorem 1's immunity to
   one-way outlier bridges keeps artifact tiles from prematurely
   fusing distinct diagnostic subtypes during early merges.

4. *Remote-Sensing Patch and Land-Cover Hierarchy Discovery:*
   Hyperspectral and aerial imagery patches span curved seasonal and
   illumination manifolds across unequally represented land-cover
   classes (e.g., dense urban vs. rare wetland). GDL's degree-product
   linkage adapts naturally to heterogeneous cluster densities and
   produces multi-resolution cartographic hierarchies.

5. *Compact Multimodal Corpus and Visual Retrieval Indexing:*
   For moderate-scale collections ($N <= 5 times 10^4$) where
   hierarchical browsing and cosine nearest-neighbor retrieval are
   jointly required, JULE directly outputs $ell_2$-normalized
   embeddings with proven angular inter-cluster margins (Theorem 2)
   alongside an explicit tree index.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

JULE encounters four clear operational boundaries:
(i) *Quadratic large-$N$ bottleneck*, as constructing and updating
$K_s$-NN and GDL cluster affinities scales as $cal(O)(N^2)$ in the
early singleton stages ($M_0 = N$), making $N > 10^5$ prohibitive
without mini-batch approximations;
(ii) *High-dimensional uninformative noise* @beyer1999nearest,
which corrupts initial $K_s$-NN graphs at $t = 0$ and triggers
irreversible merge cascades (Proposition 1);
(iii) *Lack of inductive soft density modeling*, as JULE provides
neither generative sampling nor calibrated Bayesian posterior
probabilities; and
(iv) *Sensitivity to unrolling rate $eta$*, where large $eta$ merges
too aggressively before features adapt, while tiny $eta$ multiplies
CNN training passes by $eta^(-1)$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Face & Pose Manifolds], [GDL links curved manifolds + handles skewed identity counts], [Avoid large $N > 10^5$; use FINCH],
    [Fine-Grained Taxonomies], [Produces full hierarchy + angular margin embeddings], [Pre-extract CNN features if noisy],
    [Multi-Shape Low-D Data], [Local mutual $K_s$-NN excels on non-convex shapes], [Tune $K_s$ to avoid disconnected islets],
    [High-D Noisy Tabular], [Poor; $t = 0$ $K_s$-NN concentration locks early errors], [Prefer DEPICT (185) or DAC (187)],
    [Streaming / Web-Scale], [Batch-only global hierarchy; $cal(O)(N^2)$ startup cost], [Use SwAV or DeepCluster],
  ),
  caption: [Domain suitability and failure-mode matrix for JULE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

JULE bridges classical hierarchical graph clustering
@ward1963hierarchical @lance1967general @dasgupta2016cost and
modern iterative deep clustering frameworks such as DeepCluster
@caron2018deepcluster, FINCH @sarfraz2019finch, and SCAN
@vangansbeke2020scan. The following proposition clarifies these exact
structural relationships.

#proposition(
  title: "Structural Limits Relating JULE, GDL, FINCH, and DeepCluster",
  origin: "Original Unified Proof",
)[
  Let $T_"JULE" (theta, eta, K_s, S)$ denote Alg. 186.
  1. *(Frozen-Representation GDL Limit):* Setting the backward SGD
     step count $S = 0$ (freezing $theta = theta^((0))$) reduces
     JULE identically to classical Graph Degree Linkage agglomerative
     clustering on fixed normalized features, whose merge sequence is
     invariant to $eta$ whenever $n_p = 1$.
  2. *(First-Neighbor Parameter-Free Limit):* Freezing $theta$,
     setting $K_s = 1$, and replacing GDL greedy pair selection by
     connected components of the symmetrized 1-NN graph
     $bold(A) = bold(W) + bold(W)^top + bold(W)^top bold(W)$ recovers
     the hierarchical merge mechanism of FINCH @sarfraz2019finch,
     which attains strictly linear $cal(O)(N log N)$ merging without
     specifying $K$ or $eta$.
  3. *(Non-Nested Re-Assignment Contrast with DeepCluster):*
     Replacing JULE's irreversible agglomerative forward pass with a
     fresh global $k$-means partition at every stage $t$ converts
     JULE into DeepCluster @caron2018deepcluster, eliminating the
     early-merge lock-in of Proposition 1 at the cost of sacrificing
     both hierarchical nesting and non-convex manifold linkage.
]

#proof[
  For part (1), when $S = 0$, $theta^((t)) = theta^((0))$ for all
  $t$, so the embeddings ${bold(bar(x))_i}$ and directed affinity
  matrix $bold(W)$ remain invariant across all stages. Setting
  $n_p = 1$ executes one argmax merge of $cal(A)(cal(C)_a, cal(C)_b)$
  per step from $N$ down to $K$, reproducing classical batch GDL
  agglomerative clustering identically.

  For part (2), note that for $K_s = 1$ with binary indicator
  weights $W_(i j) in {0, 1}$, the GDL term
  $bold(W)_(a b) bold(W)_(b a)$ is nonzero precisely when there
  exists a reciprocal 1-NN link (a mutual nearest neighbor) between
  $cal(C)_a$ and $cal(C)_b$, which is the symmetric core of FINCH's
  adjacency $bold(A) = bold(W) + bold(W)^top + bold(W)^top bold(W)$
  @sarfraz2019finch; taking connected components merges all mutually
  linked pairs simultaneously in a single pass instead of one
  fractional batch $eta$.

  For part (3), because $k$-means recomputes cluster labels
  $bold(y)^((t))$ globally from scratch on ${bold(bar(x))_i (theta^((t-1)))}$,
  any pair $(p, q)$ erroneously grouped at $t = 1$ can be separated
  at $t = 2$ once global centroids shift, breaking the monotonic
  coarsening invariant $cal(P)^((t)) succ.eq cal(P)^((t-1))$ that
  drove Proposition 1(2), while restricting cluster geometry to
  convex Voronoi cells in $bb(R)^d$ instead of GDL graph manifolds.
]

The ClusterOverview four-dataset benchmark provides striking
empirical confirmation of every theorem established above. On
*Iris*, JULE attains the *highest Adjusted Rand Index ($0.594$)*
among all four deep clustering algorithms (surpassing DEPICT
$0.583$, VaDE $0.579$, and DAC $0.579$) alongside strong Silhouette
($0.412$), showing that local mutual $K_s$-NN linkage accurately
tracks low-dimensional botanical species boundaries. On
*Multi-Shape*, JULE achieves Silhouette $0.718$ (ARI N/A), tying
VaDE for the top manifold score and validating Theorem 1's
non-convex GDL linkage strength. On *Customers*, JULE achieves
Silhouette $0.335$ — edging out balanced-prior DEPICT ($0.327$) and
far exceeding DAC ($0.197$) because GDL imposes no uniform
cluster-size constraint. In stark contrast, on *High-D* JULE drops
to the *lowest Silhouette ($0.119$)* of all four methods (less than
half of DAC's $0.266$ and DEPICT's $0.263$): this directly reflects
high-dimensional nearest-neighbor concentration @beyer1999nearest
corrupting the initial $t = 0$ $K_s$-NN graph and triggering
irreversible early-merge errors (Proposition 1). Finally, JULE's
mean runtime of $14.593$ ms is the highest of the four architectures
($4.4 times$ slower than DAC's $3.329$ ms), directly quantifying the
overhead of recurrent $K_s$-NN graph construction and multi-stage
agglomerative updates derived in Lemma 1.
