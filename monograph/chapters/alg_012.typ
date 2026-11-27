#import "../template/ieee.typ": *

#alg_entry(
  id: 12,
  name: "Spherical K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2001",
  paradigm: "Unit-Hypersphere Cosine / Directional Partitioning",
  complexity: $O("nnz"(X) k T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

In text mining, information retrieval, and deep neural representation learning, feature vectors frequently exhibit wide variation in Euclidean length—for example, document length in TF-IDF bags of words—while semantic similarity is encoded primarily by *direction*. The *Spherical K-Means* algorithm, formalized by Dhillon and Modha @dhillon2001concept, restricts both data representations and cluster prototypes to the unit $(d - 1)$-hypersphere $bb(S)^(d-1) = {x in bb(R)^d : norm(x)_2 = 1}$ and replaces Euclidean distortion with *cosine similarity*.

Given raw non-zero vectors ${tilde(x)_1, dots, tilde(x)_n} subset.eq bb(R)^d backslash {bold(0)}$, define the $L_2$-normalized dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(S)^(d-1)$ via $x_i = tilde(x)_i / norm(tilde(x)_i)_2$. For a target cluster count $k in {1, dots, n}$, Spherical K-Means seeks a disjoint partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ and unit-norm *concept vectors* $bold(mu) = (mu_1, dots, mu_k) in (bb(S)^(d-1))^k$ maximizing total within-cluster cosine coherence @dhillon2001concept:
$ Q(cal(C), bold(mu)) &= sum_(j=1)^k sum_(i in C_j) chevron.l x_i, mu_j chevron.r, \
                      &quad "subject to" quad norm(mu_j)_2 = 1. $
Equivalently, Spherical K-Means minimizes total *cosine dissimilarity* $J_"cos" (cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) (1 - chevron.l x_i, mu_j chevron.r) = n - Q(cal(C), bold(mu))$. Letting $s_j = sum_(i in C_j) x_i$ denote the unnormalized resultant vector of cluster $C_j$, maximizing $Q$ over $bold(mu) in (bb(S)^(d-1))^k$ for a fixed partition eliminates the prototypes in closed form (Lemma 1 below) and yields the purely combinatorial *resultant-length* objective:
$ Q^*(cal(C)) = sum_(j=1)^k norm(s_j)_2 = sum_(j=1)^k norm(sum_(i in C_j) x_i)_2. $
Conversely, holding $bold(mu)$ fixed partitions $bb(S)^(d-1)$ into *spherical Voronoi cells*
$ V_j^bb(S) (bold(mu)) = {&x in bb(S)^(d-1) : chevron.l x, mu_j chevron.r >= chevron.l x, mu_ell chevron.r \
                          &quad forall ell != j}, $
which are intersections of closed central polyhedral cones with $bb(S)^(d-1)$ bounded by great hyperspheres.

#definition(title: "Spherical Centroidal Voronoi Equilibrium")[
  A pair $(cal(C), bold(mu))$ with $bold(mu) in (bb(S)^(d-1))^k$ is a *spherical centroidal Voronoi tessellation (SCVT)*—a coordinate-wise stationary point of $Q$—if: (i) every $i in C_j$ satisfies $x_i in V_j^bb(S) (bold(mu))$; and (ii) for every cluster $C_j$ with non-zero resultant $s_j = sum_(i in C_j) x_i != bold(0)$, the prototype equals the normalized resultant $mu_j = s_j / norm(s_j)_2$.
]

Algorithm 1 implements alternating maximization of $Q(cal(C), bold(mu))$, which is directly compatible with sparse matrix representations of high-dimensional corpora.

#algbox(
  title: "Spherical K-Means Algorithm",
  inputs: [Non-zero vectors $tilde(cal(X)) subset.eq bb(R)^d$, integer $k >= 1$, initial unit prototypes $bold(mu)^((0)) in (bb(S)^(d-1))^k$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and unit concept vectors $bold(mu)$.],
)[
  1. *Normalize inputs:* Set $x_i <- tilde(x)_i / norm(tilde(x)_i)_2$ for $i = 1, dots, n$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) *Directional assignment:* For each $i in {1, dots, n}$, assign
     $ c_i^((t)) = arg max_(j in {1, dots, k}) chevron.l x_i, mu_j^((t)) chevron.r, $
     breaking ties deterministically by smallest cluster index $j$.
  4. #h(1em) Form $C_j^((t)) = {i : c_i^((t)) = j}$ and resultants $s_j^((t)) = sum_(i in C_j^((t))) x_i$.
  5. #h(1em) *Spherical mean update:* For each $j = 1, dots, k$: if $norm(s_j^((t)))_2 > 0$, set
     $ mu_j^((t+1)) = s_j^((t)) / norm(s_j^((t)))_2; $
     otherwise reseed $mu_j^((t+1))$ to the point $x_i$ minimizing $max_ell chevron.l x_i, mu_ell^((t)) chevron.r$.
  6. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$ (or $bold(mu)^((t+1)) = bold(mu)^((t))$), *break*.
  7. *Return* $(cal(C)^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of Spherical K-Means: (1) an exact spherical analogue of the parallel-axis identity yielding strict quantitative ascent and finite-step termination; (2) chordal assignment equivalence combined with a 4-point constructive proof that global Spherical and Euclidean K-Means minimizers *strictly disagree* even on unit-normalized data; (3) exact equivalence to the zero-temperature hard-EM limit of von Mises–Fisher (vMF) mixture models; (4) orthogonal and per-point radial invariance together with an exact tight angular corruption bound against adversarial outliers; and (5) structural connections to concept decompositions and spectral clustering.

#subsec("A", "Spherical Pythagorean Identity and Finite Convergence")

#lemma(title: "Exact Resultant-Direction Optimality Gap", origin: "Original Proof")[
  Let $S subset.eq bb(S)^(d-1)$ be a finite subset with non-zero resultant $s = sum_(x in S) x != bold(0)$, and let $hat(mu) = s / norm(s)_2 in bb(S)^(d-1)$. Then for every unit vector $u in bb(S)^(d-1)$,
  $ norm(s)_2 - sum_(x in S) chevron.l x, u chevron.r = 1 / 2 norm(s)_2 norm(u - hat(mu))_2^2 >= 0, $
  with equality if and only if $u = hat(mu)$. Hence $hat(mu)$ is the unique maximizer of directional coherence over $bb(S)^(d-1)$.
]

#proof[
  By linearity of the inner product, $sum_(x in S) chevron.l x, u chevron.r = chevron.l s, u chevron.r = norm(s)_2 chevron.l hat(mu), u chevron.r$. Since $norm(u)_2 = norm(hat(mu))_2 = 1$, expanding the squared chordal distance yields
  $ norm(u - hat(mu))_2^2 &= norm(u)_2^2 + norm(hat(mu))_2^2 - 2 chevron.l u, hat(mu) chevron.r \
                          &= 2 (1 - chevron.l hat(mu), u chevron.r). $
  Multiplying both sides by $(1 / 2) norm(s)_2 > 0$ gives
  $ 1 / 2 norm(s)_2 norm(u - hat(mu))_2^2 &= norm(s)_2 (1 - chevron.l hat(mu), u chevron.r) \
                                          &= norm(s)_2 - sum_(x in S) chevron.l x, u chevron.r. $
  Because $norm(s)_2 > 0$ and $norm(u - hat(mu))_2^2 >= 0$ vanishes if and only if $u = hat(mu)$, the normalized resultant uniquely maximizes coherence on $bb(S)^(d-1)$.
]

#theorem(title: "Quantitative Monotone Ascent and Finite Termination", origin: "Original Proof")[
  Let ${(cal(C)^((t)), bold(mu)^((t)))}$ be generated by Algorithm 1 with deterministic smallest-index tie-breaking and non-vanishing cluster resultants $norm(s_j^((t)))_2 > 0$. Then at every iteration $t >= 0$:
  1. $Q(cal(C)^((t)), bold(mu)^((t+1))) = Q(cal(C)^((t)), bold(mu)^((t))) + 1 / 2 sum_(j=1)^k norm(s_j^((t)))_2 norm(mu_j^((t+1)) - mu_j^((t)))_2^2$.
  2. $Q(cal(C)^((t+1)), bold(mu)^((t+1))) >= Q(cal(C)^((t)), bold(mu)^((t+1)))$, with strict inequality whenever any point switches to a prototype of strictly higher cosine similarity.
  3. Algorithm 1 terminates in at most $T <= S(n, k) <= k^n / k!$ iterations at a spherical CVT equilibrium (strengthening @dhillon2001concept with an exact step-gain identity).
]

#proof[
  Holding $cal(C)^((t))$ fixed and summing Lemma 1 over $j = 1, dots, k$ with $u = mu_j^((t))$ and $hat(mu) = mu_j^((t+1))$ immediately yields Claim 1. In the assignment step at iteration $t + 1$, every point $x_i$ satisfies $chevron.l x_i, mu_(c_i^((t+1)))^((t+1)) chevron.r >= chevron.l x_i, mu_(c_i^((t)))^((t+1)) chevron.r$; summing over $i = 1, dots, n$ proves Claim 2. Thus $Q_t := Q(cal(C)^((t)), bold(mu)^((t+1))) = Q^*(cal(C)^((t)))$ is monotonically non-decreasing.

  To exclude cycling, suppose $Q_(t+1) = Q_t$ while $cal(C)^((t+1)) != cal(C)^((t))$. Then Claim 1 enforces $mu_j^((t+2)) = mu_j^((t+1))$ for all $j$, and Claim 2 forces every reassigned point $i$ to achieve equal inner product with $mu_(c_i^((t)))^((t+1))$ and $mu_(c_i^((t+1)))^((t+1))$. Under smallest-index deterministic tie-breaking, stationary prototypes produce identical assignments, contradicting $cal(C)^((t+1)) != cal(C)^((t))$. Consequently, every partition visited prior to halting achieves a strictly higher value of $Q^*(cal(C))$ and cannot repeat. Finiteness of the Stirling partition count $S(n, k)$ guarantees finite-step convergence to an SCVT.
]

#subsec("B", "Non-Equivalence to Euclidean K-Means on Unit Vectors")

A pervasive misconception in applied machine learning asserts that Spherical K-Means is equivalent to running Lloyd's Euclidean K-Means (Alg. 1) on $L_2$-normalized data. We prove that while *assignment* coincides for unit prototypes, the *global partition objectives* are strictly inequivalent.

#theorem(title: "Chordal Equivalence vs. Global Partition Disagreement", origin: "Original Counterexample")[
  Let $cal(X) subset.eq bb(S)^(d-1)$ be a unit-normalized dataset.
  1. For any unit prototypes $bold(mu) in (bb(S)^(d-1))^k$, $norm(x_i - mu_j)_2^2 = 2 (1 - chevron.l x_i, mu_j chevron.r)$, so minimizing squared chordal distance is identical to maximizing cosine similarity.
  2. However, unconstrained Euclidean K-Means places centroids in the *interior* ball ($norm(overline(x)_j)_2 <= 1$) and maximizes $B_"Euc" (cal(C)) := n - J^*(cal(C)) = sum_(j=1)^k norm(s_j)_2^2 / |C_j|$, whereas Spherical K-Means maximizes $Q^*(cal(C)) = sum_(j=1)^k norm(s_j)_2$. There exists a 4-point dataset in $bb(S)^2 subset.eq bb(R)^3$ with $k = 2$ whose unique global Spherical K-Means optimum and unique global Euclidean K-Means optimum are *strictly different partitions*.
]

#proof[
  (1) Direct expansion gives $norm(x_i - mu_j)_2^2 = norm(x_i)_2^2 + norm(mu_j)_2^2 - 2 chevron.l x_i, mu_j chevron.r = 2(1 - chevron.l x_i, mu_j chevron.r)$.

  (2) For Euclidean K-Means on unit vectors, $sum_(i in C_j) norm(x_i - overline(x)_j)_2^2 = |C_j| - |C_j| norm(overline(x)_j)_2^2 = |C_j| - norm(s_j)_2^2 / |C_j|$, so $J^*(cal(C)) = n - B_"Euc" (cal(C))$. Now construct $cal(X) = {x_1, x_2, x_3, x_4} subset.eq bb(S)^2$ by
  $ &x_1 = (3/5, 4/5, 0), quad x_2 = (3/5, -4/5, 0), \
    &x_3 = (0, 0, 1), quad x_4 = (-1, 0, 0), $
  with $k = 2$. Exhaustively evaluating all $2^3 - 1 = 7$ non-trivial 2-partitions yields:
  - $P_1 = {{1, 2}, {3, 4}}$: $s_1 = (6/5, 0, 0)$, $s_2 = (-1, 0, 1)$, giving $Q^*(P_1) = 6/5 + sqrt(2) approx 2.6142$ and $B_"Euc" (P_1) = 18/25 + 1 = 1.7200$.
  - $P_2 = {{1, 2, 3}, {4}}$: $s_1 = (6/5, 0, 1)$ ($norm(s_1)_2^2 = 61/25$), $s_2 = (-1, 0, 0)$, giving $Q^*(P_2) = sqrt(61)/5 + 1 approx 2.5620$ and $B_"Euc" (P_2) = 61/75 + 1 approx 1.8133$.
  - $P_3 = {{1, 2, 4}, {3}}$: $Q^*(P_3) = 1.2000$, $B_"Euc" (P_3) = 76/75 approx 1.0133$.
  - $P_4 = {{1, 3}, {2, 4}}$ and symmetric $P_5 = {{2, 3}, {1, 4}}$: $norm(s_1)_2^2 = 2$, $norm(s_2)_2^2 = 4/5$, giving $Q^* = sqrt(2) + 2/sqrt(5) approx 2.3086$ and $B_"Euc" = 1 + 2/5 = 1.4000$.
  - $P_6 = {{1}, {2, 3, 4}}$ and symmetric $P_7 = {{2}, {1, 3, 4}}$: $norm(s_1)_2 = 1$, $norm(s_2)_2^2 = 9/5$, giving $Q^* = 1 + 3/sqrt(5) approx 2.3416$ and $B_"Euc" = 1 + 3/5 = 1.6000$.

  Therefore, $P_1$ is the *unique* global maximizer of $Q^*$ ($2.6142 > 2.5620$), whereas $P_2$ is the *unique* global minimizer of Euclidean WCSS ($B_"Euc" = 1.8133 > 1.7200$). Intuitively, because $overline(x)_({1,2}) = (3/5, 0, 0)$ shrinks inward toward the origin, Euclidean K-Means shifts the bisector toward the smaller-norm interior centroid, pulling the orthogonal point $x_3$ into ${1, 2}$, whereas Spherical K-Means projects $mu_({1,2}) = (1, 0, 0)$ back onto $bb(S)^2$, preserving equal angular radii.
]

#subsec("C", "Hard-EM Limit of von Mises–Fisher Mixture Models")

#theorem(title: "Generative Equivalence to Homoscedastic vMF Hard-EM", origin: "Classical — Cited")[
  Consider a $k$-component mixture of von Mises–Fisher (vMF) distributions on $bb(S)^(d-1)$ with densities @banerjee2005vmf
  $ p(x | mu_j, kappa_j) = C_d (kappa_j) exp(kappa_j chevron.l x, mu_j chevron.r), $
  where $C_d (kappa) = kappa^(d/2 - 1) / ((2 pi)^(d/2) I_(d/2 - 1) (kappa))$ and $I_nu$ is the modified Bessel function of the first kind. Under uniform mixing weights $pi_j = 1 / k$ and shared concentration $kappa_1 = dots = kappa_k = kappa$, the Expectation-Maximization algorithm @dempster1977em reduces identically to Algorithm 1 in the zero-temperature limit $kappa -> +infinity$ @banerjee2005vmf.
]

#proof[
  In the E-step with $pi_j = 1 / k$ and shared $kappa_j = kappa$, the normalizing constants $C_d (kappa)$ cancel identically, giving posterior responsibilities
  $ gamma_(i j)^((kappa)) = exp(kappa chevron.l x_i, mu_j chevron.r) / (sum_(ell=1)^k exp(kappa chevron.l x_i, mu_ell chevron.r)). $
  Let $m_i = max_ell chevron.l x_i, mu_ell chevron.r$ with unique maximizer $j^*(i)$ and gap $Delta_(i ell) := m_i - chevron.l x_i, mu_ell chevron.r >= 0$. Dividing numerator and denominator by $exp(kappa m_i)$ yields
  $ gamma_(i j)^((kappa)) = exp(-kappa Delta_(i j)) / (1 + sum_(ell != j^*(i)) exp(-kappa Delta_(i ell))). $
  As $kappa -> +infinity$, every strictly suboptimal exponent tends to $-infinity$, yielding $lim_(kappa -> infinity) gamma_(i j)^((kappa)) = bb(1)_(j = j^*(i))$, which is precisely Step 3 of Algorithm 1.

  In the M-step, maximizing the expected complete-data log-likelihood over $mu_j in bb(S)^(d-1)$ reduces to maximizing the Lagrangian $cal(L)(mu_j, lambda_j) = kappa sum_(i=1)^n gamma_(i j) chevron.l x_i, mu_j chevron.r - lambda_j (norm(mu_j)_2^2 - 1)$. Setting $nabla_(mu_j) cal(L) = bold(0)$ yields $2 lambda_j mu_j = kappa sum_(i=1)^n gamma_(i j) x_i$. Enforcing $norm(mu_j)_2 = 1$ with $lambda_j > 0$ (required for maximization) yields $mu_j = (sum_i gamma_(i j) x_i) / norm(sum_i gamma_(i j) x_i)_2$. Substituting the hard limit $gamma_(i j) = bb(1)_(i in C_j)$ recovers Step 5 identically.
]

#subsec("D", "Invariances and Exact Angular Outlier-Robustness Bound")

#proposition(title: "Directional Invariances and Tight Single-Outlier Deflection Bound", origin: "Original Proof")[
  1. *Orthogonal and Positive-Ray Invariance:* The normalized problem and all SCVT equilibria are strictly invariant under any global orthogonal transformation $tilde(x)_i |-> Q tilde(x)_i$ ($Q in O(d)$) and under arbitrary *pointwise* positive radial rescalings $tilde(x)_i |-> alpha_i tilde(x)_i$ ($alpha_i > 0$).
  2. *Translation Non-Invariance:* Shifting $tilde(x)_i |-> tilde(x)_i + b$ ($b != bold(0)$) is not an invariance of Spherical K-Means.
  3. *Exact Bounded Angular Corruption:* Let cluster $C_j$ have uncorrupted resultant $s_j$ of magnitude $R_j = norm(s_j)_2 > 1$ and prototype $mu_j = s_j / R_j$. Inserting any single adversarial point $tilde(x)_"adv" != bold(0)$ into $C_j$ deflects the updated prototype $tilde(mu)_j = (s_j + x_"adv") / norm(s_j + x_"adv")_2$ by an angle $theta = arccos chevron.l mu_j, tilde(mu)_j chevron.r$ satisfying the *exact tight bound*
     $ theta <= arcsin (1 / R_j), $
     with equality achieved whenever $chevron.l s_j, x_"adv" chevron.r = -1$.
]

#proof[
  (1) For $Q in O(d)$ and $alpha_i > 0$, $(Q (alpha_i tilde(x)_i)) / norm(Q (alpha_i tilde(x)_i))_2 = Q x_i$, and $chevron.l Q x_i, Q mu_j chevron.r = chevron.l x_i, mu_j chevron.r$, so both assignment and normalized-resultant updates commute identically with $Q$.

  (2) Take $tilde(x)_1 = (1, 1)$, $tilde(x)_2 = (1, -1)$ (orthogonal, cosine $0$); translating by $b = (100, 0)$ yields nearly collinear vectors with cosine $> 0.9998$, altering angular distances and optimal partitions.

  (3) Regardless of $norm(tilde(x)_"adv")_2$, pre-normalization maps it to a unit vector $x_"adv" in bb(S)^(d-1)$. Consider the triangle formed by the origin $O$, $S = s_j$ (length $R_j > 1$), and $V = s_j + x_"adv"$ (with side $S V = norm(x_"adv")_2 = 1$). Let $theta = arccos chevron.l mu_j, tilde(mu)_j chevron.r in [0, pi/2)$ be the angle at $O$ and $psi$ the interior angle at $V$. By the Law of Sines,
  $ (sin theta) / 1 = (sin psi) / R_j <= 1 / R_j, $
  hence $theta <= arcsin(1 / R_j)$. Equality holds when $psi = pi / 2$, i.e., $x_"adv" perp (s_j + x_"adv")$, which simplifies to $chevron.l s_j, x_"adv" chevron.r + norm(x_"adv")_2^2 = 0$, or $chevron.l s_j, x_"adv" chevron.r = -1$ (realizable in $d >= 2$ whenever $R_j > 1$). In contrast to Euclidean K-Means (where a single outlier drives $norm(mu_j)_2 -> infinity$), coherent clusters with $R_j = Omega(|C_j|)$ sustain only $O(1 / |C_j|)$ angular perturbation, and overturning a prototype to its antipode requires at least $ceil(R_j)$ adversarial points.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Spherical K-Means attains exceptional computational efficiency on sparse high-dimensional data. Storing $cal(X)$ in Compressed Sparse Row (CSR) format with $"nnz"(X)$ non-zero entries, evaluating $chevron.l x_i, mu_j chevron.r$ for all $i, j$ costs $O("nnz"(X) k)$ fused multiply-adds per iteration, while accumulating resultants $s_j$ and $L_2$-normalizing costs $O("nnz"(X) + k d)$. Across $T$ iterations, total work is $O(("nnz"(X) + d) k T)$ (simplifying to $O("nnz"(X) k T)$ when $"nnz"(X) >= d$, and $O(n k d T)$ on dense inputs) with memory $O("nnz"(X) + k d)$. Note that while sparse document vectors remain sparse, dense concept vectors $mu_j$ generally fill in unless explicitly pruned or top-$m$ truncated. In high dimensions, unit normalization also mitigates radial norm explosion in distance concentration @beyer1999nearest, concentrating discriminative signal into pairwise angular structure.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Spherical K-Means],
    [Time / Space Complexity], [$O("nnz"(X) k T)$ sparse ($O(n k d T)$ dense); $O("nnz"(X) + k d)$ space.],
    [Geometric Inductive Bias], [Equal-concentration spherical caps (central Voronoi cones) on $bb(S)^(d-1)$.],
    [Generative Counterpart], [Hard-EM limit ($kappa -> infinity$) of equal-concentration vMF mixtures @banerjee2005vmf.],
    [Invariances], [Invariant under $Q in O(d)$ and pointwise scaling $alpha_i > 0$; *not* translation-invariant.],
    [Outlier Robustness], [Norm-immune; single outlier deflects prototype by at most $arcsin(1 / R_j)$ (Prop. 1).],
    [Recommended Upgrades], [vMF-EM @banerjee2005vmf for unequal spreads; Cosine K-Means++ seeding; Kernel K-Means (Alg. 16) for non-linear arcs.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Spherical K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Spherical K-Means is the canonical partitioning workhorse across modern directional, textual, and embedding-centric domains where vector length reflects nuisance volume rather than latent class identity.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Text, News, and Patent Clustering:* In TF-IDF vector spaces, doubling a document's length doubles its Euclidean norm without changing its topic. Pointwise radial invariance (Proposition 1(1)) guarantees identical clustering for short abstracts and long articles, while CSR matrix-vector inner products scale effortlessly to millions of documents @dhillon2001concept @steinbach2000comparison.
2. *Dense Neural Embeddings and Retrieval Indexing:* Self-supervised contrastive encoders (InfoNCE @oord2018infonce, SwAV @caron2020swav, CLIP, Sentence-BERT) explicitly optimize cosine similarity on the unit hypersphere $bb(S)^(d-1)$. Spherical K-Means aligns exactly with the training geometry for semantic deduplication, topic discovery, and cosine IVF vector-database partitioning.
3. *Biometric Speaker and Face Template Clustering:* Deep speaker embeddings (x-vectors/ECAPA) and angular-margin face representations (ArcFace) reside on $bb(S)^(d-1)$ by construction; Spherical K-Means recovers identity prototypes while capping single-frame corruption via Proposition 1(3).
4. *Genomic Co-Expression Profile Grouping:* Clustering genes across experimental conditions often prioritizes *relative up/down regulation patterns* rather than absolute transcript abundance. After per-gene centering and $L_2$ normalization, cosine similarity equals the Pearson correlation coefficient, making Spherical K-Means an exact correlation-clustering primitive.
5. *Recommender-System Item and User Cohorting:* Sparse user-item interaction vectors suffer severe activity-bias (power users vs. casual users); radial normalization isolates taste direction from purchase volume, enabling fast cosine neighborhood partitioning.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Spherical K-Means degrades under four conditions: (a) *magnitude-informative data*, where radial distance from the origin distinguishes classes (e.g., low- vs. high-spending customers along the same spending mix) and normalization erases class separation; (b) *zero or near-zero vectors*, where $x_i / norm(x_i)_2$ is undefined or dominated by measurement noise; (c) *antipodal cancellation*, where opposing points within a cluster cause $s_j approx bold(0)$ and erratic prototype flips; and (d) *unequal angular concentrations*, where equal-$kappa$ bisectors cause compact clusters to lose boundary points to wide diffuse clusters (resolved by soft vMF-EM @banerjee2005vmf).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Spherical Fits / Prerequisite], [Caution / Alternative],
    [TF-IDF / NLP], [Document-length invariant; sparse speed], [Filter empty/stopword docs],
    [Neural Embeddings], [Matches contrastive cosine geometry], [Use vMF-EM for varied $kappa_j$],
    [Gene Expression], [Equals Pearson correlation after centering], [Exclude flat zero-variance genes],
    [Customer Spend], [Captures basket *mix* only], [Use Euclidean (Alg. 1/13) for volume],
    [Geospatial / Physical], [Origin-dependent angles distort layout], [Use Lloyd (Alg. 1) or DBSCAN],
  ),
  caption: [Application suitability and domain engineering guide for Spherical K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Spherical K-Means connects directly to low-rank matrix factorization, normalized kernel methods, and spectral graph partitioning via the following structural identity.

#proposition(title: "Concept Decomposition and Normalized-Kernel Equivalence", origin: "Original Derivation")[
  Let $X in bb(R)^(n times d)$ have unit-norm rows $x_i^T$, let $M = [mu_1, dots, mu_k] in bb(R)^(d times k)$ hold the unit concept vectors, and let $G in {0, 1}^(n times k)$ be the hard assignment indicator matrix.
  1. Maximizing $Q(cal(C), bold(mu)) = "Tr"(G^T X M)$ is strictly equivalent to minimizing the constrained rank-$k$ matrix approximation error $norm(X - G M^T)_F^2 = 2 n - 2 "Tr"(G^T X M)$ over $G$ and column-normalized $M$.
  2. Unlike Dhillon's Kernel K-Means @dhillon2004kernel with normalized cosine kernel $K = X X^T$ (which minimizes $n - sum_j (bold(1)_(C_j)^T K bold(1)_(C_j)) / |C_j|$), Spherical K-Means maximizes $sum_j sqrt(bold(1)_(C_j)^T K bold(1)_(C_j))$—replacing each cluster's mean Gram energy by its *square-root total Gram energy*.
]

#proof[
  (1) Because every row of $X$ and every row of $G M^T$ (namely $mu_(c_i)^T$) has unit Euclidean norm, $norm(X)_F^2 = norm(G M^T)_F^2 = n$. Expanding the Frobenius norm yields
  $ norm(X - G M^T)_F^2 &= norm(X)_F^2 + norm(G M^T)_F^2 - 2 "Tr"(G^T X M) \
                        &= 2 n - 2 Q(cal(C), bold(mu)), $
  linking Spherical K-Means directly to Dhillon and Modha's *concept decomposition* $X approx Z M^T$ @dhillon2001concept (obtained by subsequently solving unconstrained least squares for $Z$ given $M$).

  (2) Writing $norm(s_j)_2^2 = norm(sum_(i in C_j) x_i)_2^2 = sum_(i, ell in C_j) chevron.l x_i, x_ell chevron.r = bold(1)_(C_j)^T K bold(1)_(C_j)$, Spherical K-Means maximizes $Q^*(cal(C)) = sum_(j=1)^k (bold(1)_(C_j)^T K bold(1)_(C_j))^(1/2)$, whereas unweighted Kernel K-Means on the same Gram matrix $K$ maximizes $sum_(j=1)^k |C_j|^(-1) (bold(1)_(C_j)^T K bold(1)_(C_j))$, explaining algebraically the global partition divergence demonstrated in Theorem 2. Similarly, Ng–Jordan–Weiss spectral clustering @ng2002spectral row-normalizes Laplacian eigenvectors onto $bb(S)^(k-1)$, where Spherical K-Means often replaces Euclidean K-Means to respect angular community separation.
]

Empirical results on the ClusterOverview 4-dataset benchmark suite (measured by Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing) closely corroborate this theory. On *Iris*, Spherical K-Means achieves *Sil = 0.444* and *ARI = 0.589*, edging past both Elkan's Euclidean K-Means (ARI = 0.583) and G-Means (ARI = 0.567) because botanical petal-to-sepal aspect ratios act as directional signatures. On *Customers*, however, its Silhouette drops to *Sil = 0.319* (below Elkan's 0.334 and G-Means's 0.342), directly illustrating the loss of radial spend-volume information predicted by Proposition 1(1). On *Multi-Shape*, it records *Sil = 0.711* (ARI = N/A), essentially matching Euclidean Voronoi partitioning (0.712). Most notably, on *High-D Subspace+Noise*, Spherical K-Means attains *Sil = 0.303* (ARI = N/A)—substantially outperforming Euclidean Elkan (0.239) by $+0.064$—demonstrating the directional concentration stability and bounded outlier deflection established in Proposition 1(3). Finally, with a mean runtime of just *0.048 ms*, it is the fastest method among Algorithms 11–13, benefiting from branch-free unit inner products and minimal per-iteration bookkeeping.
