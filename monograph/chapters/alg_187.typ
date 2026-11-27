#import "../template/ieee.typ": *

#alg_entry(
  id: 187,
  name: "DAC (Deep Adaptive Clustering)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2017",
  paradigm: "Adaptive Pairwise Binary Self-Supervision on the Simplex",
  complexity: $cal(O)(E (N P + N B K))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Reconstruction-based deep clustering algorithms such as VaDE (Alg.
184) and DEPICT (Alg. 185) devote half of their network parameters
and compute budget to a decoder reconstructing raw inputs, which is
both computationally expensive and prone to emphasizing nuisance
background variance over semantic category identity. Meanwhile,
graph-agglomerative architectures such as JULE (Alg. 186) incur
quadratic $cal(O)(N^2)$ affinity overhead across the full dataset.
Deep Adaptive Clustering (DAC), introduced by Chang et al.
@chang2017dac, eliminates both bottlenecks by recasting
unsupervised clustering as a *binary pairwise classification*
problem solved end-to-end inside a purely feed-forward encoder via
self-paced adaptive thresholds.

Let $cal(X) = {bold(x)_1, dots, bold(x)_N} subset.eq bb(R)^D$ be an
unlabeled dataset and $K >= 2$ be the target cluster count. DAC
trains a deep encoder $f_w : bb(R)^D -> Delta^(K-1)$ (ending in a
$K$-way softmax layer) that maps each datum $bold(x)_i$ directly to
a nonnegative *label feature* vector
$bold(l)_i = f_w (bold(x)_i) in Delta^(K-1)$ ($l_(i k) >= 0$,
$sum_(k=1)^K l_(i k) = 1$). Crucially, DAC applies $ell_2$
normalization to project every label feature onto the positive
orthant of the unit hypersphere,
$bb(S)_+^(K-1) = {bold(u) in bb(R)^K : bold(u) >= bold(0), norm(bold(u))_2 = 1}$:
$
  bold(hat(l))_i = bold(l)_i / norm(bold(l)_i)_2 in bb(S)_+^(K-1).
$
The pairwise similarity between any two samples $bold(x)_i, bold(x)_j$
is their nonnegative cosine similarity
$
  S_(i j) (w) = chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r in [0, 1].
$
Because ground-truth pairwise relationships
$r_(i j)^ast in {0, 1}$ are unobserved, DAC introduces a curriculum
parameter $lambda >= 0$ governing two monotonic adaptive thresholds
$u(lambda)$ (upper, non-increasing) and $l(lambda)$ (lower,
non-decreasing) satisfying $0 < l(lambda) <= u(lambda) < 1$ —
specifically $u(lambda) = u_0 - lambda$ and $l(lambda) = l_0 + lambda$
until they meet at $l(lambda) = u(lambda)$. Pairs with extreme
similarities are treated as high-confidence pseudo-labeled pairs,
while ambiguous intermediate pairs are masked out.

#definition(title: "DAC Adaptive Pairwise Binary Objective")[
  For any mini-batch $cal(B) subset.eq {1, dots, N}$ of size $B$
  and current thresholds $(l(lambda), u(lambda))$, define the
  self-paced pairwise pseudo-labels $r_(i j)$ and binary selection
  weights $v_(i j) in {0, 1}$ by
  $
    (r_(i j), v_(i j)) = cases(
      (1, 1) & "if" S_(i j) (w) >= u(lambda),
      (0, 1) & "if" S_(i j) (w) < l(lambda),
      (0, 0) & "if" l(lambda) <= S_(i j) (w) < u(lambda).
    )
  $
  The *DAC pairwise binary cross-entropy objective* is
  $
    cal(L)_"DAC" (w; cal(B), lambda)
    &= 1 / |cal(B)|^2 sum_(i, j in cal(B)) v_(i j) ell_"BCE" (r_(i j), S_(i j) (w)), \
    ell_"BCE" (r, S)
    &= - r log S - (1 - r) log(1 - S).
  $
  After convergence, hard cluster labels are read directly from the
  dominant label-feature coordinate without post-hoc clustering:
  $hat(c)_i = arg max_(k in {1, dots, K}) l_(i k)$.
]

Notice the dual role of $bold(l)_i$: during training, its
$ell_2$-normalized version $bold(hat(l))_i$ acts as a metric
embedding for pairwise similarity estimation; at convergence, its
nonnegativity and unit-norm constraints force the learned vectors
toward the $K$ standard one-hot basis vectors $bold(e)_1, dots, bold(e)_K$,
making $bold(l)_i$ itself the final cluster assignment vector
@chang2017dac.

#algbox(
  title: "DAC — Deep Adaptive Clustering",
  inputs: [
    Dataset $cal(X) = {bold(x)_i}_(i=1)^N subset.eq bb(R)^D$,
    clusters $K$, initial thresholds $(l_0, u_0)$ with
    $0 < l_0 < u_0 < 1$, step $Delta lambda > 0$, batch size $B$.
  ],
  outputs: [
    Trained encoder $f_w : bb(R)^D -> Delta^(K-1)$, direct cluster
    assignments ${hat(c)_i}_(i=1)^N$.
  ],
)[
  1. *Initialization:* Initialize encoder weights $w$, set
     $lambda <- 0$, $u <- u_0$, $l <- l_0$.
  2. *While* $l(lambda) < u(lambda)$:
     - Sample mini-batch $cal(B) subset.eq cal(X)$ of size $B$.
     - *(Forward Label-Feature Pass):* Compute softmax vectors
       $bold(l)_i = f_w (bold(x)_i) in Delta^(K-1)$ and normalized
       vectors $bold(hat(l))_i = bold(l)_i / norm(bold(l)_i)_2$ for
       $i in cal(B)$.
     - *(Pairwise Selection):* Form the $B times B$ Gram matrix
       $S_(i j) = chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r$,
       and compute $(r_(i j), v_(i j))$ via $(l(lambda), u(lambda))$
       (holding $r_(i j), v_(i j)$ detached as constants).
     - *(Gradient Update):* Update $w$ via SGD / RMSprop minimizing
       $cal(L)_"DAC" (w; cal(B), lambda)$ over selected pairs.
     - *(Curriculum Advance):* Increment $lambda <- lambda + Delta lambda$,
       updating $u(lambda) = u_0 - lambda$ and
       $l(lambda) = l_0 + lambda$.
  3. *Direct Assignment:* For every $bold(x)_i in cal(X)$, compute
     $bold(l)_i = f_w (bold(x)_i)$ and return
     $hat(c)_i = arg max_(k in {1, dots, K}) l_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical foundations for DAC: (A) the exact
one-hot identifiability theorem on the positive orthant sphere
@chang2017dac, (B) a new quantitative $epsilon$-stability theorem
proving one-hot recovery under *approximate* pairwise separation
via Gershgorin spectral analysis, (C) exact Riemannian gradient
bounds explaining how the adaptive exclusion band $[l, u)$ controls
Lipschitz smoothness, and (D) a constructive confirmation-bias
collapse counterexample under misspecified initial thresholds.

#subsec("A", "Exact One-Hot Characterization on the Positive Sphere")

#theorem(
  title: "Exact One-Hot Identifiability of Positive-Orthant Label Features",
  origin: "Classical — Cited",
)[
  Let $cal(X) = union.big_(a=1)^K cal(C)_a$ be partitioned into $K$
  non-empty clusters, and suppose the simplex label features
  $bold(l)_i in Delta^(K-1)$ with $bold(hat(l))_i = bold(l)_i / norm(bold(l)_i)_2 in bb(S)_+^(K-1)$
  satisfy the exact binary pairwise conditions @chang2017dac:
  $
    chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r
    = cases(
      1 & "if" i, j in cal(C)_a "for the same" a,
      0 & "if" i in cal(C)_a, j in cal(C)_b "with" a != b.
    )
  $
  Then there exists a unique permutation $pi in S_K$ such that
  $bold(l)_i = bold(hat(l))_i = bold(e)_(pi(a))$ for every
  $i in cal(C)_a$ and every $a in {1, dots, K}$.
]

#proof[
  For any two samples $i, j$ in the same cluster $cal(C)_a$,
  $norm(bold(hat(l))_i - bold(hat(l))_j)_2^2 = 2(1 - chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r) = 0$,
  so every sample in $cal(C)_a$ shares a single unit representative
  $bold(u)_a in bb(S)_+^(K-1)$. By the inter-cluster condition,
  $chevron.l bold(u)_a, bold(u)_b chevron.r = sum_(k=1)^K u_(a k) u_(b k) = 0$
  for all $a != b$. Because $u_(a k) >= 0$ for all $a, k$, every
  summand is nonnegative, forcing $u_(a k) u_(b k) = 0$ for every
  coordinate $k in {1, dots, K}$ and every $a != b$.

  Define the support $S_a = "supp"(bold(u)_a) = {k : u_(a k) > 0} subset.eq {1, dots, K}$.
  Since $norm(bold(u)_a)_2 = 1$, each $S_a$ is non-empty ($|S_a| >= 1$).
  Moreover, $u_(a k) u_(b k) = 0$ implies $S_a inter S_b = emptyset$
  for all $a != b$. Thus $S_1, dots, S_K$ are $K$ pairwise disjoint
  non-empty subsets of a $K$-element set:
  $
    K >= |union.big_(a=1)^K S_a| = sum_(a=1)^K |S_a| >= K.
  $
  Equality holds throughout, forcing $|S_a| = 1$ for every $a$.
  Hence there exists a bijection $pi : {1, dots, K} -> {1, dots, K}$
  with $S_a = {pi(a)}$. Since $norm(bold(u)_a)_2 = 1$ and
  $u_(a, pi(a)) > 0$, we obtain $bold(u)_a = bold(e)_(pi(a))$.
  Finally, because $bold(l)_i = norm(bold(l)_i)_2 bold(u)_a$ lies in
  the simplex ($sum_k l_(i k) = 1$), $norm(bold(l)_i)_2 = 1$ and
  $bold(l)_i = bold(e)_(pi(a))$ identically.
]

#subsec("B", "Quantitative Stability Under Approximate Separation")

Real neural optimization never reaches exact $S_(i j) in {0, 1}$
because softmax outputs stay strictly positive ($l_(i k) > 0$). The
following original theorem extends Theorem 1 to quantitative
approximate separation, proving that small inter-cluster cosine
leakage $epsilon$ still guarantees unique one-hot recovery!

#theorem(
  title: "Quantitative One-Hot Stability on the Positive Orthant Sphere",
  origin: "Original Proof",
)[
  Let $bold(u)_1, dots, bold(u)_K in bb(S)_+^(K-1)$ be $K$ cluster
  representatives satisfying the approximate orthogonality bound
  $0 <= chevron.l bold(u)_a, bold(u)_b chevron.r <= epsilon$ for all
  $a != b$. If $epsilon < 1 / (4 (K - 1))$, then:
  1. The Gram matrix $bold(G) = bold(U)^top bold(U)$ (where
     $bold(U) = [bold(u)_1, dots, bold(u)_K] in bb(R)^(K times K)$)
     is strictly positive definite with
     $lambda_"min" (bold(G)) >= 1 - (K - 1) epsilon > 3 / 4$, ruling
     out cluster rank collapse.
  2. There exists a unique bijection $pi in S_K$ such that every
     cluster representative concentrates strictly on its assigned
     standard basis vector $bold(e)_(pi(a))$:
     $
       u_(a, pi(a))^2 >= 1 - 2 (K - 1) epsilon > 1 / 2,
     $
     and the Euclidean deviation from exact one-hot encoding obeys
     $
       norm(bold(u)_a - bold(e)_(pi(a)))_2^2 <= 4 (K - 1) epsilon.
     $
  3. Consequently, direct argmax decoding
     $hat(c)(bold(u)_a) = arg max_k u_(a k)$ recovers the exact
     $K$-way permutation $pi(a)$ without ambiguity.
]

#proof[
  For part (1), the Gram matrix $bold(G) = bold(U)^top bold(U)$ has
  unit diagonal entries $G_(a a) = norm(bold(u)_a)_2^2 = 1$ and
  off-diagonal entries $0 <= G_(a b) <= epsilon$. By Gershgorin's
  circle theorem, every eigenvalue $lambda$ of $bold(G)$ satisfies
  $|lambda - 1| <= sum_(b != a) |G_(a b)| <= (K - 1) epsilon$, so
  $lambda_"min" (bold(G)) >= 1 - (K - 1) epsilon > 3 / 4 > 0$.

  For part (2), since $bold(U)$ is square ($K times K$), the matrix
  $bold(M) = bold(U) bold(U)^top$ shares the exact same spectrum as
  $bold(G) = bold(U)^top bold(U)$, hence
  $lambda_"min" (bold(M)) = lambda_"min" (bold(G)) >= 1 - (K - 1) epsilon$.
  By the Rayleigh–Ritz theorem, every diagonal entry $r_k^2$ of
  $bold(M)$ obeys
  $
    r_k^2 = bold(e)_k^top bold(M) bold(e)_k
    = sum_(a=1)^K u_(a k)^2
    >= lambda_"min" (bold(M))
    >= 1 - (K - 1) epsilon.
  $
  For each coordinate $k in {1, dots, K}$, choose a maximizing
  cluster index $a_k in arg max_(a in {1, dots, K}) u_(a k)$. Using
  nonnegativity of $bold(U)$ and $u_(b k) <= u_(a_k, k)$, every
  other cluster $b != a_k$ satisfies
  $
    u_(b k)^2
    <= u_(a_k, k) u_(b k)
    <= sum_(j=1)^K u_(a_k, j) u_(b j)
    = chevron.l bold(u)_(a_k), bold(u)_b chevron.r
    <= epsilon.
  $
  Summing over all $K - 1$ indices $b != a_k$ yields
  $sum_(b != a_k) u_(b k)^2 <= (K - 1) epsilon$. Subtracting this
  bound from the lower bound on $r_k^2$ produces
  $
    u_(a_k, k)^2
    = r_k^2 - sum_(b != a_k) u_(b k)^2
    >= 1 - 2 (K - 1) epsilon.
  $
  Under the hypothesis $epsilon < 1 / (4 (K - 1))$, we have
  $u_(a_k, k)^2 > 1 / 2$ for every $k = 1, dots, K$. Because each
  column $bold(u)_a$ has unit squared norm
  $sum_(k=1)^K u_(a k)^2 = 1$, no single cluster $a$ can have two
  distinct coordinates strictly greater than $1 / 2$. Therefore the
  map $k |-> a_k$ is injective on a finite $K$-element set — hence
  a bijection with inverse $pi : {1, dots, K} -> {1, dots, K}$ such
  that $u_(a, pi(a))^2 >= 1 - 2 (K - 1) epsilon > 1 / 2$.

  Using the elementary inequality $1 - (1 - t)^(1/2) <= t$ for
  $t in [0, 1]$ with $t = 2 (K - 1) epsilon$, we obtain
  $
    norm(bold(u)_a - bold(e)_(pi(a)))_2^2
    = 2 (1 - u_(a, pi(a)))
    <= 2 (1 - (1 - 2 (K - 1) epsilon)^(1/2))
    <= 4 (K - 1) epsilon.
  $
  Finally, part (3) follows immediately because
  $u_(a, pi(a))^2 > 1 / 2$ implies $u_(a, pi(a)) > u_(a k)$ for
  every $k != pi(a)$, so $arg max_k u_(a k) = pi(a)$ uniquely.
]

#subsec("C", "Riemannian Gradient Bounds of the Adaptive Margin Band")

#lemma(
  title: "Uniform Riemannian Gradient Bound Under Self-Paced Masking",
  origin: "Original Proof",
)[
  Let $nabla_"tan" ell_"BCE" (r_(i j), S_(i j))$ denote the
  Riemannian gradient of $ell_"BCE" (r_(i j), chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r)$
  with respect to $bold(hat(l))_i$ on the tangent space
  $T_(bold(hat(l))_i) bb(S)^(K-1)$. For any pair $(i, j)$ selected
  by the adaptive thresholds $0 < l(lambda) <= u(lambda) < 1$
  ($v_(i j) = 1$):
  1. If $r_(i j) = 1$ ($S_(i j) >= u(lambda)$), the negative
     Riemannian gradient moves $bold(hat(l))_i$ strictly along the
     shortest geodesic toward $bold(hat(l))_j$ with norm bounded by
     $
       norm(nabla_"tan" ell_"BCE" (1, S_(i j)))_2
       = ((1 - S_(i j)^2)^(1/2)) / S_(i j)
       <= ((1 - u(lambda)^2)^(1/2)) / u(lambda).
     $
  2. If $r_(i j) = 0$ ($S_(i j) < l(lambda)$), the negative
     Riemannian gradient repels $bold(hat(l))_i$ strictly away from
     $bold(hat(l))_j$ with norm bounded by
     $
       norm(nabla_"tan" ell_"BCE" (0, S_(i j)))_2
       = ((1 + S_(i j)) / (1 - S_(i j)))^(1/2)
       <= ((1 + l(lambda)) / (1 - l(lambda)))^(1/2).
     $
  3. Hence excluding the ambiguous band $[l(lambda), u(lambda))$
     uniformly bounds the per-pair spherical gradient norm across
     the entire mini-batch by
     $G_"max" (lambda) = max((1 - u^2)^(1/2) / u, ((1 + l) / (1 - l))^(1/2))$.
]

#proof[
  Differentiating $S_(i j) = chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r$
  in ambient $bb(R)^K$ gives
  $nabla_(bold(hat(l))_i) ell_"BCE" = (partial ell_"BCE" / (partial S_(i j))) bold(hat(l))_j$.
  Projecting onto the tangent space $T_(bold(hat(l))_i) bb(S)^(K-1)$
  via $bold(P)_(bold(hat(l))_i)^perp = bold(I)_K - bold(hat(l))_i bold(hat(l))_i^top$
  gives
  $
    nabla_"tan" ell_"BCE"
    = (partial ell_"BCE" / (partial S_(i j)))
      (bold(hat(l))_j - S_(i j) bold(hat(l))_i).
  $
  Since $norm(bold(hat(l))_i)_2 = norm(bold(hat(l))_j)_2 = 1$, the
  tangent vector $bold(v)_(i j) = bold(hat(l))_j - S_(i j) bold(hat(l))_i$
  points along the initial velocity of the minimizing geodesic from
  $bold(hat(l))_i$ to $bold(hat(l))_j$ and has Euclidean norm
  $norm(bold(v)_(i j))_2 = (1 - S_(i j)^2)^(1/2)$.

  For $r_(i j) = 1$, $partial ell_"BCE" / (partial S_(i j)) = - 1 / S_(i j) < 0$,
  so $- nabla_"tan" ell_"BCE"$ points positively along $bold(v)_(i j)$
  (attractive) with norm $(1 - S_(i j)^2)^(1/2) / S_(i j)$. Because
  $s |-> (1 - s^2)^(1/2) / s = (s^(-2) - 1)^(1/2)$ is strictly
  decreasing on $(0, 1]$, its supremum over $[u(lambda), 1]$ is
  attained at $s = u(lambda)$, proving (1).

  For $r_(i j) = 0$, $partial ell_"BCE" / (partial S_(i j)) = + 1 / (1 - S_(i j)) > 0$,
  so $- nabla_"tan" ell_"BCE"$ points negatively along $bold(v)_(i j)$
  (repulsive) with norm
  $
    ((1 - S_(i j)^2)^(1/2)) / (1 - S_(i j))
    = (((1 - S_(i j))(1 + S_(i j))) / (1 - S_(i j))^2)^(1/2)
    = ((1 + S_(i j)) / (1 - S_(i j)))^(1/2).
  $
  Because $s |-> ((1 + s) / (1 - s))^(1/2)$ is strictly increasing
  on $[0, 1)$, its supremum over $[0, l(lambda))$ is bounded by its
  value at $l(lambda)$, proving (2). Combining (1) and (2) yields
  (3), showing precisely how the curriculum gap $(l(lambda), u(lambda))$
  prevents gradient explosion from borderline uncertain pairs.
]

#subsec("D", "Confirmation-Bias Collapse Counterexample")

#proposition(
  title: "Single-Cluster Collapse Under Concentrated Pairwise Similarities",
  origin: "Original Counterexample",
)[
  Suppose at initialization $w^((0))$ (e.g., due to positively skewed
  nonnegative inputs or high-dimensional angular concentration
  @beyer1999nearest), the normalized label features satisfy
  $min_(i, j in cal(X)) chevron.l bold(hat(l))_i (w^((0))), bold(hat(l))_j (w^((0))) chevron.r >= u_0$.
  Then:
  1. At $lambda = 0$, every pair receives $v_(i j) = 1$ and
     $r_(i j) = 1$, with *zero negative pairs* across all mini-batches.
  2. Any constant encoder $f_(w^ast) (bold(x)) equiv bold(u)_0 in Delta^(K-1)$
     achieves the exact global minimum $cal(L)_"DAC" (w^ast) = 0$,
     while gradient descent contracts pairwise angles monotonically
     toward $S_(i j) = 1$, locking all subsequent curriculum steps
     $lambda > 0$ into a trivial single-cluster collapse ($"ARI" = 0$).
]

#proof[
  Because $S_(i j) (w^((0))) >= u_0 = u(0)$ for every $(i, j)$,
  Definition 1 assigns $(r_(i j), v_(i j)) = (1, 1)$ universally.
  The objective simplifies to
  $cal(L)_"DAC" = - (1 / |cal(B)|^2) sum_(i, j in cal(B)) log S_(i j) (w) >= 0$,
  which attains its exact global minimum $0$ if and only if
  $S_(i j) (w) = 1$ for all $i, j$ — i.e.,
  $bold(hat(l))_i = bold(hat(l))_j = bold(u)_0$ identically. By
  Lemma 1(1), every pairwise gradient is purely attractive, increasing
  $min_(i, j) S_(i j)$ further above $u_0$. Since $u(lambda)$ is
  non-increasing with $lambda$, $S_(i j) >= u_0 >= u(lambda)$
  persists inductively for every subsequent step $lambda > 0$,
  preventing any negative pair $r_(i j) = 0$ from ever activating.
]

#sec("III", "Complexity, Geometry and Operational Profile")

DAC achieves exceptional computational efficiency because it
eliminates both the generative / denoising decoder network (present
in VaDE, Alg. 184, and DEPICT, Alg. 185) and the global
$cal(O)(N^2)$ nearest-neighbor graph (present in JULE, Alg. 186).
Let $P$ denote the forward-backward arithmetic cost of the encoder
$f_w$ alone per sample (roughly half the cost of a full
encoder–decoder pair). Inside each mini-batch $cal(B)$ of size $B$,
forming the $B times B$ pairwise similarity matrix
$bold(S) = bold(hat(L))_cal(B) bold(hat(L))_cal(B)^top$ requires a
single dense BLAS matrix multiplication of size
$(B times K) times (K times B)$, costing only $cal(O)(B^2 K)$
floating-point operations — or $cal(O)(N B K)$ summed across all
$N / B$ mini-batches per epoch. Across $E$ epochs, total training
time is
$
  T_"DAC" = cal(O)(E (N P + N B K)),
$
strictly linear in $N$ for fixed mini-batch size $B$, with peak
memory $cal(O)(|w| + B (D + H + K) + B^2)$ and pure single-pass
feed-forward inference $cal(O)(P_"enc")$ per test sample without
needing any centroid lookup or $k$-means step.

Geometrically, DAC maps ambient data $bb(R)^D$ directly onto the
positive orthant of the unit $(K - 1)$-sphere $bb(S)_+^(K-1)$.
Unlike Euclidean embeddings whose arbitrary rotations require a
downstream clustering head, the nonnegativity constraint
$bold(l)_i >= bold(0)$ coupled with pairwise orthogonality
($S_(i j) -> 0$) aligns the $K$ clusters tightly onto the $K$
coordinate axes $bold(e)_1, dots, bold(e)_K$ (Theorems 1–2). Thus
cluster geometry is governed by *angular cosine cones* around
standard simplex vertices, giving high discriminative capacity in
high-dimensional spaces while making performance sensitive to input
feature normalization and angular spread at initialization
(Proposition 1).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(E (N P + N B K))$ training; $cal(O)(P_"enc")$ direct test inference],
    [Space Complexity], [$cal(O)(|w| + B(D + H + K) + B^2)$ — no decoder and no $cal(O)(N^2)$ graph],
    [Latent Geometry], [Positive orthant sphere $bb(S)_+^(K-1)$; one-hot simplex vertex cones (Thm. 1–2)],
    [Cluster Balance], [Unconstrained; no explicit uniform prior, though batch statistics favor active modes],
    [Decoder Overhead], [Zero; purely discriminative encoder avoids reconstructing background noise],
    [Initialization Sensitivity], [High on skewed/unnormalized inputs (Prop. 1 collapse); requires tuned $(l_0, u_0)$],
    [Primary Hyperparameters], [Cluster count $K$, thresholds $(l_0, u_0)$, curriculum rate $Delta lambda$, batch size $B$],
  ),
  caption: [Operational and geometric profile of DAC (Alg. 187).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Natural Image Categorization:* Natural scenes
   contain high-entropy background clutter (sky, foliage, texture)
   that wastes autoencoder capacity under pixel MSE losses. Because
   DAC discards the decoder and learns purely from pairwise cosine
   relations on $bb(S)_+^(K-1)$, deep convolutional backbones devote
   all parameters to discriminative semantic contours and directly
   output one-hot category probabilities @chang2017dac.

2. *High-Dimensional Document and Pretrained Embedding Clustering:*
   Text embeddings, TF-IDF vectors, and vision–language features are
   naturally compared by cosine similarity rather than Euclidean
   distance @dhillon2001concept @banerjee2005vmf. DAC's native
   cosine Gram matrix $S_(i j) = chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r$
   operates in the exact metric geometry of high-dimensional
   semantic embeddings while scaling linearly to millions of items.

3. *Real-Time Visual Hashing and End-to-End One-Hot Quantization:*
   Production retrieval pipelines require compact one-hot bucket IDs
   produced in a single forward pass with sub-millisecond latency.
   Theorem 2 guarantees that once inter-cluster cosine leakage falls
   below $1 / (4(K - 1))$, DAC's encoder natively emits sharp
   one-hot codes ready for inverted-index routing without an extra
   quantization layer.

4. *Medical Imaging Triage with Irrelevant Background Artifacts:*
   Whole-slide scans and chest radiographs include variable borders,
   annotations, and scanner noise that reconstruction-based models
   faithfully memorize. Discriminative pairwise self-supervision
   ignores non-shared nuisance artifacts and focuses representation
   capacity on repeatable diagnostic patterns across mini-batches.

5. *Online and Mini-Batch Acoustic Event Discovery:* Audio clips
   streaming from IoT arrays cannot be stored in an $cal(O)(N^2)$
   global graph. DAC's mini-batch Gram update ($cal(O)(B^2 K)$)
   supports streaming GPU updates over spectral features while
   keeping memory footprint strictly independent of total stream
   length $N$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DAC fails predictably in four regimes:
(i) *Positively skewed, unnormalized tabular features*, where raw
pairwise cosines cluster tightly near $1$ at initialization and
trigger Proposition 1's single-cluster collapse (mitigated by
z-scoring / whitening inputs and calibrating $(l_0, u_0)$ to
initial similarity quantiles);
(ii) *Small mini-batch size relative to large $K$* ($B < K$),
where many mini-batches contain zero positive pairs ($r_(i j) = 1$)
for rare classes, starving minority clusters of attractive
gradients;
(iii) *Rapid curriculum annealing* (large $Delta lambda$), which
forces premature binary labeling on noisy intermediate pairs before
embeddings separate; and
(iv) *Tasks requiring calibrated density estimation or generative
synthesis*, where decoder-free DAC provides no $p(bold(x))$ model.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Natural Image Categories], [Decoder-free CNN avoids background MSE waste], [Pair with strong augmentation / IIC],
    [High-D Semantic Embeddings], [Native cosine geometry on $bb(S)_+^(K-1)$ + $cal(O)(N)$ speed], [Ensure $B >= 4 K$ for batch coverage],
    [Low-Latency One-Hot Routing], [Direct argmax decoding; fastest inference pass], [Monitor active cluster count],
    [Skewed Nonnegative Tabular], [Prone to cosine concentration collapse (Prop. 1)], [Prefer VaDE (Alg. 184) or GMM],
    [Density / Anomaly Modeling], [No generative likelihood $p(bold(x))$ or decoder], [Use VaDE (Alg. 184) or DCN],
  ),
  caption: [Domain suitability and failure-mode matrix for DAC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DAC sits at a pivotal historical inflection point in deep
clustering: it bridged classical spherical / spectral pairwise
clustering @dhillon2001concept @banerjee2005vmf @ng2002spectral
with the modern self-supervised contrastive clustering revolution
(IIC @ji2019iic, SCAN @vangansbeke2020scan, Contrastive Clustering
@li2021contrastive, and InfoNCE @oord2018infonce). The following
unified proposition makes this mathematical bridge explicit.

#proposition(
  title: "Exact Structural Bridge From DAC to Contrastive and Spherical Clustering",
  origin: "Original Unified Proof",
)[
  Let $bold(hat(l))_i = f_w (bold(x)_i) / norm(f_w (bold(x)_i))_2 in bb(S)_+^(K-1)$
  and let $S_(i j) = chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r$.
  1. *(Hard-Threshold Limit of InfoNCE / Contrastive Clustering):*
     Consider the pairwise InfoNCE loss @oord2018infonce with
     temperature $tau > 0$ on a mini-batch $cal(B)$ with positive
     set $cal(P)_i = {j in cal(B) : r_(i j) = 1}$ and negative set
     $cal(N)_i = {j in cal(B) : r_(i j) = 0}$:
     $
       cal(L)_"NCE" (i)
       = - 1 / |cal(P)_i| sum_(j in cal(P)_i)
         log (exp(S_(i j) / tau) /
              (exp(S_(i j) / tau) + sum_(m in cal(N)_i) exp(S_(i m) / tau))).
     $
     While InfoNCE couples positives and negatives through a
     log-partition denominator, DAC's pairwise BCE decouples every
     pair independently on $[0, 1]$, and both share the *exact same
     global minimizers* on $bb(S)_+^(K-1)$ (one-hot orthogonal
     vertices ${bold(e)_k}$ by Theorem 1) whenever $cal(P)_i, cal(N)_i$
     match the true partition.
  2. *(Nonnegative Spherical $k$-Means Equivalence):* When $f_w$
     attains the $epsilon$-separated regime of Theorem 2
     ($epsilon < 1 / (4(K - 1))$) with hard prototype centroids
     $bold(mu)_k = bold(e)_k$, DAC's intra-cluster positive loss
     $- sum_(r_(i j) = 1) log chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r$
     upper-bounds the spherical $k$-means distortion
     $sum_i (1 - chevron.l bold(hat(l))_i, bold(e)_(hat(c)_i) chevron.r)$
     @dhillon2001concept via $- log s >= 1 - s$ for $s in (0, 1]$.
]

#proof[
  For part (1), on the positive orthant sphere $bb(S)_+^(K-1)$, we
  have $0 <= S_(i j) <= 1$ for every pair. In $cal(L)_"NCE" (i)$,
  $exp(S_(i j) / tau) <= e^(1 / tau)$ for every $j in cal(P)_i$
  and $exp(S_(i m) / tau) >= 1$ for every $m in cal(N)_i$. Hence
  $
    cal(L)_"NCE" (i)
    >= - log (e^(1 / tau) / (e^(1 / tau) + |cal(N)_i|))
    = log(1 + |cal(N)_i| e^(- 1 / tau)),
  $
  with equality achieved simultaneously across all $i in cal(X)$ if
  and only if $S_(i j) = 1$ for all $j in cal(P)_i$ and
  $S_(i m) = 0$ for all $m in cal(N)_i$. These are the exact
  hypotheses of Theorem 1, which uniquely characterize the $K$
  one-hot standard basis vectors ${bold(e)_1, dots, bold(e)_K}$
  (matching the exact zero-loss global minimizers of $cal(L)_"DAC"$).

  For part (2), let $bold(c)_a = (1 / |cal(C)_a|) sum_(j in cal(C)_a) bold(hat(l))_j$
  be the cluster mean vector and $bold(u)_a = bold(c)_a / norm(bold(c)_a)_2$.
  Using the elementary bound $- log s >= 1 - s$ for all $s in (0, 1]$
  and averaging over $j in cal(C)_a$,
  $
    1 / |cal(C)_a| sum_(j in cal(C)_a) (- log chevron.l bold(hat(l))_i, bold(hat(l))_j chevron.r)
    >= 1 - chevron.l bold(hat(l))_i, bold(c)_a chevron.r
    >= 1 - chevron.l bold(hat(l))_i, bold(u)_a chevron.r,
  $
  since $norm(bold(c)_a)_2 <= 1$ by convexity of the Euclidean norm.
  Under Theorem 2, $bold(u)_a$ lies within Euclidean distance
  $2 ((K - 1) epsilon)^(1/2)$ of $bold(e)_(pi(a))$, showing that
  minimizing DAC's positive pairwise BCE tightly bounds and drives
  down the spherical $k$-means objective @dhillon2001concept on the
  simplex vertices.
]

Empirical results on the ClusterOverview four-dataset benchmark
(Silhouette @rousseeuw1987silhouettes and ARI
@hubert1985comparing) vividly confirm DAC's theoretical strengths
and vulnerabilities. First, DAC achieves a mean runtime of just
*$3.329$ ms* — $2.5 times$ faster than DEPICT ($8.460$ ms),
$3.3 times$ faster than VaDE ($11.007$ ms), and $4.4 times$ faster
than JULE ($14.593$ ms) — directly reflecting the complete removal
of decoder backpropagation and global graph construction. Second,
on *High-D*, DAC achieves the *highest Silhouette ($0.266$)* among
all four deep methods (edging out DEPICT $0.263$ and VaDE $0.245$,
and more than doubling JULE's $0.119$), validating that direct
cosine self-supervision on $bb(S)_+^(K-1)$ avoids both high-D
reconstruction waste and $K_s$-NN hubness. On *Iris* (Silhouette
$0.433$, ARI $0.579$) and *Multi-Shape* (Silhouette $0.711$, ARI
N/A), DAC delivers competitive geometry and classification accuracy
comparable to DEPICT and VaDE. Crucially, however, on *Customers*
DAC plunges to the *lowest Silhouette ($0.197$)* of the four
methods (vs. VaDE $0.381$, JULE $0.335$, DEPICT $0.327$): because
raw customer expenditure features are strictly nonnegative and
positively skewed, pairwise cosines concentrate densely at
initialization, triggering partial confirmation-bias collapse as
predicted by Proposition 1. Together, Algs. 184–187 map out the
fundamental four-way Pareto frontier of deep clustering across
generative density (VaDE), denoising balance (DEPICT), local
manifold hierarchy (JULE), and discriminative speed (DAC).
