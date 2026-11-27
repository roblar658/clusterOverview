#import "../template/ieee.typ": *

#alg_entry(
  id: 192,
  name: "Contrastive Clustering (CC)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2021",
  paradigm: "Dual-Level Instance & Cluster Contrastive Deep Learning",
  complexity: $O(E N (C_"enc" + N m + K^2))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Contrastive Clustering (CC) @li2021contrastive casts deep clustering as simultaneous representation learning and cluster assignment via a single end-to-end contrastive objective. Let $cal(X) subset.eq RR^d$ denote a data space and let $cal(T)$ be a stochastic family of label-preserving augmentations (e.g., random cropping, color jitter, or Gaussian perturbation). Given a mini-batch $X = {x_1, dots, x_N} subset cal(X)$ of $N$ instances and a target cluster count $K >= 2$, CC draws two independent augmentations $t_i^a, t_i^b tilde.op cal(T)$ per instance to form paired views $x_i^a = t_i^a (x_i)$ and $x_i^b = t_i^b (x_i)$.

A shared neural backbone $f_theta : RR^d -> RR^p$ encodes both views into latent features $h_i^v = f_theta (x_i^v)$ for $v in {a, b}$. Two nonlinear MLP projection heads operate in parallel on $h_i^v$:
1. An *instance-level head* $g_I : RR^p -> SS^(m-1)$ producing $ell_2$-normalized embeddings $z_i^v = g_I (h_i^v) / norm(g_I (h_i^v))_2$, stacked row-wise into $bold(Z)^a, bold(Z)^b in RR^(N times m)$.
2. A *cluster-level head* $g_C : RR^p -> Delta^(K-1)$ ending in a `"softmax"` layer, producing soft cluster assignment rows $y_i^v = g_C (h_i^v)$ stacked into row-stochastic matrices $bold(Y)^a, bold(Y)^b in [0, 1]^(N times K)$.

CC's central structural insight is *row-column duality*: row $i$ of $bold(Y)^v$ is the soft cluster probability of instance $i$, whereas column $k$, denoted $hat(y)_k^v in [0, 1]^N$, is a representation of cluster $k$ across the mini-batch. Writing cosine similarity as $s(u, v) = chevron.l u, v chevron.r / (norm(u)_2 norm(v)_2)$, the instance-level InfoNCE loss @oord2018infonce at temperature $tau_I > 0$ is
$
  ell_i^a = -log (exp(s(z_i^a, z_i^b) / tau_I)) / (sum_(j=1)^N [bb(1)_(j != i) exp(s(z_i^a, z_j^a) / tau_I) + exp(s(z_i^a, z_j^b) / tau_I)]) ,
$
with $ell_i^b$ defined symmetrically, yielding $L_"ins" = (1 / (2N)) sum_(i=1)^N (ell_i^a + ell_i^b)$. Analogously, at cluster temperature $tau_C > 0$, treating $(hat(y)_k^a, hat(y)_k^b)$ as a positive cluster pair and the remaining $2K - 2$ columns as negatives gives
$
  hat(ell)_k^a = -log (exp(s(hat(y)_k^a, hat(y)_k^b) / tau_C)) / (sum_(l=1)^K [bb(1)_(l != k) exp(s(hat(y)_k^a, hat(y)_l^a) / tau_C) + exp(s(hat(y)_k^a, hat(y)_l^b) / tau_C)]) .
$
To avoid trivial assignment collapse, CC subtracts the marginal cluster assignment entropy $H(bold(Y)) = -sum_(k=1)^K [P_k^a log P_k^a + P_k^b log P_k^b]$, where $P_k^v = (1 / N) sum_(i=1)^N Y_(i k)^v$. The joint objective is
$
  L_"CC" (theta, g_I, g_C) = L_"ins" + 1 / (2K) sum_(k=1)^K (hat(ell)_k^a + hat(ell)_k^b) - lambda_H H(bold(Y)) ,
$
with default entropy weight $lambda_H = 1$.

#definition(title: "Dual-Level Contrastive Partition")[
  A parameterization $(theta, g_I, g_C)$ induces a *dual-level contrastive partition* on $cal(X)$ if instance embeddings maximize cross-view cosine agreement on $SS^(m-1)$ under $L_"ins"$ while the cluster head $g_C$ maps each $x in cal(X)$ to a hard label $hat(c)(x) = arg max_(k in {1, dots, K}) [g_C (f_theta (x))]_k$ without requiring a post-hoc clustering stage.
]

#algbox(
  title: "Contrastive Clustering (CC)",
  inputs: [Dataset $cal(D)$, cluster count $K$, batch size $N$, temperatures $(tau_I, tau_C)$, entropy weight $lambda_H$, augmentation family $cal(T)$, epochs $E$.],
  outputs: [Trained encoder $f_theta$, cluster head $g_C$, and cluster assignments ${hat(c)(x_i)}_(i=1)^(|cal(D)|)$.],
)[
  1. Initialize backbone $f_theta$, instance MLP $g_I$, and cluster MLP $g_C$.
  2. *For* epoch $e = 1, dots, E$ *do*:
  3. #h(1em) Sample mini-batch $X = {x_1, dots, x_N} subset cal(D)$ and independent augmentations $t_i^a, t_i^b tilde.op cal(T)$.
  4. #h(1em) Compute shared features $h_i^v = f_theta (t_i^v (x_i))$, normalized embeddings $z_i^v = g_I (h_i^v) / norm(g_I (h_i^v))_2$, and soft assignments $y_i^v = g_C (h_i^v)$ for $v in {a, b}$.
  5. #h(1em) Compute pairwise instance similarities $bold(S)^I in [-1, 1]^(2N times 2N)$ and instance loss $L_"ins"$.
  6. #h(1em) Extract columns $hat(y)_k^v$ of $bold(Y)^v$, compute cluster similarities $bold(S)^C in [0, 1]^(2K times 2K)$, marginal probabilities $P_k^v = (1/N) norm(hat(y)_k^v)_1$, and $L_"clu" = (1/(2K)) sum_k (hat(ell)_k^a + hat(ell)_k^b) - lambda_H H(bold(Y))$.
  7. #h(1em) Update $(theta, g_I, g_C)$ via stochastic gradient descent (Adam) on $nabla (L_"ins" + L_"clu")$.
  8. *Return* $hat(c)(x) = arg max_k [g_C (f_theta (x))]_k$ evaluated on unaugmented inputs $x in cal(D)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Global Lower Bound and Exact Optimality Characterization")

Unlike unconstrained unit-sphere embeddings whose cosine similarities span $[-1, 1]$, the cluster-head columns lie in the nonnegative orthant $[0, 1]^N$, imposing a sharp geometric floor on $L_"clu"$.

#lemma(title: "Exact Minimum and Hard-Assignment Optimality of Cluster Contrastive Loss", origin: "Original Derivation")[
  Let $bold(Y)^a, bold(Y)^b in [0, 1]^(N times K)$ be row-stochastic matrices with strictly nonzero columns ($hat(y)_k^v != bold(0)$ for all $k, v$), $tau_C > 0$, and $lambda_H > 0$. Then:
  1. The cluster-level objective satisfies the exact lower bound
  $
    L_"clu" (bold(Y)^a, bold(Y)^b) >= log(1 + (2K - 2) e^(-1 / tau_C)) - 2 lambda_H log K .
  $
  2. When $K$ divides $N$, equality is attained *if and only if* $bold(Y)^a = bold(Y)^b in {0, 1}^(N times K)$ is a view-consistent, one-hot, exactly balanced hard partition with $sum_(i=1)^N Y_(i k)^v = N / K$ for every $k in {1, dots, K}$.
]

#proof[
  Dividing numerator and denominator of $hat(ell)_k^a$ by $exp(s(hat(y)_k^a, hat(y)_k^b) / tau_C)$ yields
  $
    hat(ell)_k^a = log (1 &+ sum_(l != k) exp((s(hat(y)_k^a, hat(y)_l^a) - s(hat(y)_k^a, hat(y)_k^b)) / tau_C) \
                          &+ sum_(l != k) exp((s(hat(y)_k^a, hat(y)_l^b) - s(hat(y)_k^a, hat(y)_k^b)) / tau_C)) .
  $
  Because $bold(Y)^a, bold(Y)^b$ have nonnegative entries, Cauchy--Schwarz implies $0 <= s(hat(y)_k^u, hat(y)_l^v) <= 1$ for all $k, l$ and $u, v in {a, b}$. Hence every exponent in the $2(K - 1)$ negative terms satisfies
  $
    (s(hat(y)_k^a, hat(y)_l^v) - s(hat(y)_k^a, hat(y)_k^b)) / tau_C >= (0 - 1) / tau_C = -1 / tau_C ,
  $
  so $hat(ell)_k^a >= log(1 + (2K - 2) e^(-1 / tau_C))$, and identically for $hat(ell)_k^b$. Averaging over $k$ establishes the contrastive floor. For the entropy term, strict concavity of $t |-> -t log t$ on $[0, 1]$ and $sum_(k=1)^K P_k^v = (1/N) sum_(i,k) Y_(i k)^v = 1$ give by Jensen's inequality
  $
    -sum_(k=1)^K P_k^v log P_k^v <= log K , quad v in {a, b} ,
  $
  with equality iff $P_k^v = 1 / K$ for all $k$. Combining the two inequalities yields part 1.

  For part 2, equality holds iff every contrastive bound and both entropy bounds are tight simultaneously. Tightness of the contrastive bound forces $s(hat(y)_k^a, hat(y)_k^b) = 1$ for all $k$ and $chevron.l hat(y)_k^v, hat(y)_l^u chevron.r = 0$ for all $l != k$. Since $Y_(i k)^v >= 0$, orthogonality $sum_(i=1)^N Y_(i k)^v Y_(i l)^v = 0$ forces $Y_(i k)^v Y_(i l)^v = 0$ for every row $i$ and $l != k$, so each row of $bold(Y)^v$ has at most one positive entry. Row-stochasticity $sum_k Y_(i k)^v = 1$ then compels every row to be a standard basis vector, i.e., $bold(Y)^v in {0, 1}^(N times K)$. Because $s(hat(y)_k^a, hat(y)_k^b) = 1$ implies $hat(y)_k^a = alpha_k hat(y)_k^b$ for $alpha_k > 0$, binary entries and nonzero columns force $alpha_k = 1$, hence $bold(Y)^a = bold(Y)^b$. Finally, $P_k^v = 1/K$ yields $sum_i Y_(i k)^v = N / K$. Conversely, any such balanced one-hot matrix achieves the bound identically.
]

#subsec("B", "Finite-Sample Mutual Information Lower Bound")

The instance-level head anchors the shared encoder $f_theta$ by maximizing cross-view mutual information @oord2018infonce @li2021contrastive.

#theorem(title: "Exact Finite-Sample Mutual Information Bound for Instance Contrast", origin: "Classical — Cited")[
  Let $(x_1^a, x_1^b), dots, (x_N^a, x_N^b)$ be i.i.d. pairs drawn from the joint augmentation law $p(x^a, x^b)$, inducing embeddings $(Z_i^a, Z_i^b)$ with mutual information $I(Z^a; Z^b)$. For every finite batch size $N >= 2$ and temperature $tau_I > 0$,
  $
    I(Z^a; Z^b) >= log N - EE[L_"ins"] .
  $
]

#proof[
  Fix index $i$ and define the cross-view-only loss
  $
    ell_(i,"cross")^a = -log (exp(s(Z_i^a, Z_i^b) / tau_I)) / (sum_(j=1)^N exp(s(Z_i^a, Z_j^b) / tau_I)) .
  $
  Because the within-view negative terms $sum_(j != i) exp(s(Z_i^a, Z_j^a) / tau_I) > 0$ are strictly positive, we have $ell_i^a > ell_(i,"cross")^a$ pointwise (and symmetrically $ell_i^b > ell_(i,"cross")^b$). Consider the $N$-way identification experiment: draw $J tilde.op "Unif"{1, dots, N}$, sample $V_1, dots, V_N tilde.op^( "i.i.d." ) p(z^b)$, and sample $U tilde.op p(z^a | Z^b = V_J)$. By exchangeability of the mini-batch, $EE[ell_(i,"cross")^a] = EE[-log q(J | U, V_(1:N))]$ for the softmax classifier $q(j | U, V_(1:N)) prop exp(s(U, V_j) / tau_I)$.

  Let $p_"post" (j | U, V_(1:N))$ be the Bayes posterior of $J$. By non-negativity of Kullback--Leibler divergence,
  $
    EE[-log q(J | U, V_(1:N))] &= H(J | U, V_(1:N)) + EE["KL"(p_"post" || q)] \
                               &>= H(J | U, V_(1:N)) = log N - I(J; U, V_(1:N)) .
  $
  Since $V_(1:N)$ has identical product law $product_(r=1)^N p(v_r)$ under every $J = j$, we have $I(J; V_(1:N)) = 0$. Applying the chain rule and noting that conditioning reduces entropy yields
  $
    I(J; U, V_(1:N)) &= I(J; U | V_(1:N)) = H(U | V_(1:N)) - H(U | J, V_(1:N)) \
                     &<= H(U) - H(U | J, V_J) = H(Z^a) - H(Z^a | Z^b) = I(Z^a; Z^b) ,
  $
  where $H(U | J, V_(1:N)) = H(U | J, V_J)$ holds because $U$ is conditionally independent of $V_(-J)$ given $(J, V_J)$. Combining both bounds gives $EE[ell_i^a] > EE[ell_(i,"cross")^a] >= log N - I(Z^a; Z^b)$, and averaging over $i$ and both views completes the proof.
]

#subsec("C", "Riemannian Gradient Dynamics and Hardness-Aware Repulsion")

The temperature parameters $(tau_I, tau_C)$ and $ell_2$-normalization govern both gradient stability and hard-negative focus.

#theorem(title: "Exact Spherical Gradient Norm and Hard-Negative Concentration", origin: "Original Derivation")[
  Let $tilde(z)_i^a in RR^m without {bold(0)}$ be the pre-normalized instance output with $z_i^a = tilde(z)_i^a / norm(tilde(z)_i^a)_2 in SS^(m-1)$, and let $p_(i j)^(a v) = exp(s(z_i^a, z_j^v) / tau_I) / Z_i^a$ denote the softmax probabilities over the $2N - 1$ denominator terms of $ell_i^a$. Then:
  1. The gradient with respect to $tilde(z)_i^a$ lies in the tangent space $T_(z_i^a) SS^(m-1)$ and satisfies
  $
    nabla_(tilde(z)_i^a) ell_i^a &= (bold(I) - z_i^a (z_i^a)^top) / (tau_I norm(tilde(z)_i^a)_2) [ -(1 - p_(i i)^(a b)) z_i^b \
                                 &+ sum_(j != i) p_(i j)^(a b) z_j^b + sum_(j != i) p_(i j)^(a a) z_j^a ] .
  $
  2. The Riemannian gradient norm on $SS^(m-1)$ obeys the tight self-modulated upper bound
  $
    norm(nabla_(z_i^a)^SS ell_i^a)_2 <= 2 / tau_I (1 - p_(i i)^(a b)) .
  $
  3. Analogously, $norm(nabla_(hat(y)_k^a) hat(ell)_k^a)_2 prop norm(hat(y)_k^a)_2^(-1)$, so vanishing clusters induce unbounded gradient amplification.
]

#proof[
  On $SS^(m-1)$, cosine similarity reduces to the Euclidean inner product $s(z_i^a, u) = chevron.l z_i^a, u chevron.r$. Differentiating $ell_i^a (z) = -chevron.l z, z_i^b chevron.r / tau_I + log Z_i^a (z)$ with respect to $z in RR^m$ yields
  $
    nabla_z ell_i^a (z_i^a) = 1 / tau_I [ -z_i^b + sum_(j=1)^N p_(i j)^(a b) z_j^b + sum_(j != i) p_(i j)^(a a) z_j^a ] .
  $
  Grouping the $j = i$ positive term gives the bracketed vector $bold(v)_i = -(1 - p_(i i)^(a b)) z_i^b + sum_(j != i) p_(i j)^(a b) z_j^b + sum_(j != i) p_(i j)^(a a) z_j^a$. The Jacobian of $tilde(z) |-> tilde(z) / norm(tilde(z))_2$ is $bold(J) = (bold(I) - z_i^a (z_i^a)^top) / norm(tilde(z)_i^a)_2$. By the chain rule, $nabla_(tilde(z)_i^a) ell_i^a = bold(J) nabla_z ell_i^a (z_i^a)$, proving part 1, and the Riemannian gradient is $nabla_(z_i^a)^SS ell_i^a = (bold(I) - z_i^a (z_i^a)^top) bold(v)_i / tau_I$.

  Since $bold(P) = bold(I) - z_i^a (z_i^a)^top$ is an orthogonal projection, $norm(bold(P))_"op" = 1$. Moreover, all embeddings have unit norm and the negative weights sum to $sum_(j != i) (p_(i j)^(a b) + p_(i j)^(a a)) = 1 - p_(i i)^(a b)$. By the triangle inequality,
  $
    norm(nabla_(z_i^a)^SS ell_i^a)_2 <= 1 / tau_I [(1 - p_(i i)^(a b)) + sum_(j != i) (p_(i j)^(a b) + p_(i j)^(a a))] = 2 / tau_I (1 - p_(i i)^(a b)) .
  $
  Conditional on negative selection, $p_(i j)^(a v) / (1 - p_(i i)^(a b)) prop exp(s(z_i^a, z_j^v) / tau_I)$ concentrates exponentially on the hardest negatives as $tau_I -> 0^+$. Applying the identical spherical chain rule to the column map $hat(y)_k^a |-> hat(y)_k^a / norm(hat(y)_k^a)_2$ introduces the factor $norm(hat(y)_k^a)_2^(-1)$, proving part 3.
]

#subsec("D", "Scale-Invariance Degeneracy and Necessity of Entropy Regularization")

Without the marginal entropy penalty, column-wise cosine normalization makes the cluster contrastive loss blind to catastrophic cluster imbalance.

#proposition(title: "Degenerate Collapse Counterexample without Entropy Penalty", origin: "Original Counterexample")[
  Fix any $K >= 2$, $N >= K$, and $tau_C > 0$. When $lambda_H = 0$, the severely collapsed hard partition assigning $N - K + 1$ instances to cluster $1$ and exactly $1$ instance to each cluster $k in {2, dots, K}$ achieves the *identical* global minimum of $L_"clu"$ as a perfectly balanced ground-truth partition. For any $lambda_H > 0$, entropy regularization strictly separates the two partitions by the exact objective gap
  $
    Delta L_"clu" = 2 lambda_H [ log K + (1 - (K - 1) / N) log(1 - (K - 1) / N) + (K - 1) / N log(1 / N) ] ,
  $
  which converges to $2 lambda_H log K > 0$ as $N -> oo$.
]

#proof[
  Construct $bold(Y)_"col"^a = bold(Y)_"col"^b in {0, 1}^(N times K)$ by setting rows $1, dots, N - K + 1$ to $e_1^top$ and row $N - K + k$ to $e_k^top$ for $k = 2, dots, K$. Every column $hat(y)_k$ is nonzero ($norm(hat(y)_1)_2 = sqrt(N - K + 1)$ and $norm(hat(y)_k)_2 = 1$ for $k >= 2$), and columns have mutually disjoint supports. Because cosine similarity $s(hat(y)_k, hat(y)_l)$ is strictly invariant to positive rescaling of column norms, $s(hat(y)_k^a, hat(y)_k^b) = 1$ for all $k$ and $s(hat(y)_k^u, hat(y)_l^v) = 0$ for all $l != k$. Thus $bold(Y)_"col"$ achieves $hat(ell)_k^v = log(1 + (2K - 2) e^(-1 / tau_C))$, matching the exact global lower bound of Lemma 1 when $lambda_H = 0$, while its marginal distribution is $P = (1 - (K-1)/N, 1/N, dots, 1/N)$.

  When $lambda_H > 0$, subtracting the balanced partition's objective $L_"clu" (bold(Y)_"bal")$ from $L_"clu" (bold(Y)_"col")$ cancels the identical contrastive terms and yields
  $
    Delta L_"clu" = lambda_H (H(bold(Y)_"bal") - H(bold(Y)_"col")) = 2 lambda_H (log K - H_1 (P)) ,
  $
  where $H_1 (P) = -(1 - (K-1)/N) log(1 - (K-1)/N) - ((K-1)/N) log(1/N)$. As $N -> oo$, $(1 - (K-1)/N) log(1 - (K-1)/N) -> 0$ and $((K-1)/N) log N -> 0$, so $H_1 (P) -> 0$ and $Delta L_"clu" -> 2 lambda_H log K > 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each mini-batch iteration of CC executes two encoder forward passes of cost $O(N C_"enc")$, instance similarity matrix construction in $O(N^2 m)$, and cluster-column similarity evaluation in $O(K^2 N)$ plus $O(N K)$ entropy calculation. Across $E$ epochs over a dataset of size $n$ partitioned into $n / N$ mini-batches, total training time is $O(E n (C_"enc" + N m + K^2))$. Crucially, unlike offline graph-based deep methods such as JULE @yang2016jule or spectral clustering @ng2002spectral, CC never materializes an $n times n$ affinity matrix; peak GPU memory scales only with the mini-batch as $O(N C_"act" + N^2 + K^2)$ where $C_"act"$ is backbone activation memory. At inference time, cluster prediction requires a single unaugmented forward pass through $g_C compose f_theta$ in $O(C_"enc" + p K)$ per sample, supporting streaming online inference once trained.

Geometrically, CC enforces *dual spherical separation*: instance embeddings distribute uniformly over $SS^(m-1)$ while preserving augmentation neighborhoods, and cluster probability columns become mutually orthogonal rays in the nonnegative orthant $RR_+^N$. Because the cluster head couples all $N$ batch samples simultaneously through column normalization, large batch sizes ($N >= 256$ in vision) stabilize marginal entropy estimates $P_k^v$ and supply diverse negatives, whereas very small $N$ increases variance in $H(bold(Y))$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Training Time], [$O(E n (C_"enc" + N m + K^2))$ across $E$ epochs, batch size $N$.],
    [Inference Time], [$O(C_"enc" + p K)$ per point; single forward pass, no centroid lookup.],
    [Space Complexity], [$O(abs(theta) + N C_"act" + N^2 + K^2)$ GPU memory; $O(1)$ in $n$.],
    [Inductive Bias], [Augmentation-invariant manifolds with approximately balanced cluster priors.],
    [Key Hyperparameters], [Cluster count $K$, temperatures $(tau_I, tau_C)$, batch size $N$, augmentation $cal(T)$.],
    [Scalability], [High via mini-batch SGD; supports multi-GPU distributed batch gathering.],
    [Differentiability], [Fully end-to-end differentiable; no alternating EM or offline $k$-Means step.],
  ),
  caption: [Operational and computational profile of Contrastive Clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Unsupervised Visual Category Discovery:* Natural images undergo nuisance transformations (crop, scale, illumination, color shift) that form a near-invariance group $cal(T)$. By Theorem 2, $L_"ins"$ extracts semantically invariant features on $SS^(m-1)$ while $L_"clu"$ directly separates semantic categories without collapsing into texture-dominated prototypes @li2021contrastive @vangansbeke2020scan.
2. *Single-Cell Transcriptomics (scRNA-seq) Phenotyping:* High-dimensional gene expression profiles suffer from dropout noise and technical batch effects. Modeling dropout masking and mild Gaussian jitter as $cal(T)$ enables CC to learn denoised cell embeddings while jointly discovering cell-type assignments in a single pass.
3. *Short-Text and Document Topic Clustering:* Pairing pretrained language-model backbones with token-deletion or back-translation augmentations aligns semantically equivalent sentences in $SS^(m-1)$ while column orthogonality in $bold(Y)$ sharpens topical separation, overcoming the sparsity of bag-of-words metrics.
4. *Industrial Wafer-Map and Medical Histopathology Stratification:* Defect patterns and tissue morphology exhibit rotational and stain-color symmetries. Encoding these symmetries directly into $cal(T)$ aligns CC's invariance bias with domain physics, enabling end-to-end subtyping of unlabeled scans.
5. *Tabular and Multi-Modal Telemetry Representation:* When domain-aware tabular corruptions (feature masking, subspace mixup) define $cal(T)$, CC maps heterogeneous continuous features onto a normalized hypersphere where cosine distance mitigates Euclidean scale disparity across sensors.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CC exhibits three mathematically predictable failure regimes. First, *augmentation mismatch*: if $cal(T)$ violates true class boundaries (e.g., aggressive cropping or heavy tabular noise that moves a point across manifolds), $L_"ins"$ forces distinct classes together. Second, *severe class imbalance*: by Lemma 1 and Proposition 4, maximizing $H(bold(Y))$ biases CC toward uniform cluster proportions $P_k approx 1/K$, degrading performance on long-tailed datasets unless $H(bold(Y))$ is replaced by a KL penalty toward a known prior. Third, *small-sample / low-dimensional tabular overhead*: without rich invariances or large $n$, neural contrastive training can overfit and incurs higher latency than lightweight geometric methods.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Natural Image Collections], [Strong: geometric / color $cal(T)$ maximizes semantic MI.], [Tune $tau_I, tau_C$; compare SwAV @caron2020swav / SCAN @vangansbeke2020scan.],
    [High-D Omics & Text], [Strong: spherical projection mitigates distance concentration.], [Requires domain-valid $cal(T)$; else IDEC @guo2017idec / VaDE @jiang2017vade.],
    [Long-Tailed / Rare-Class Data], [Weak: $-H(bold(Y))$ penalizes skewed cluster priors.], [Use density methods (HDBSCAN @campello2013hdbscan) or weighted prior KL.],
    [Low-D Streaming Telemetry], [Weak: iterative SGD too slow for sub-ms streams.], [Use StreamKM++ @ackermann2012streamkmpp or DenStream @cao2006denstream.],
    [Arbitrary Non-Convex 2D Shapes], [Moderate: MLP can warp manifolds if $cal(T)$ is local.], [Local noise scale critical; compare Spectral @ng2002spectral / DBSCAN @ester1996dbscan.],
  ),
  caption: [Domain suitability and failure-mode matrix for Contrastive Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CC occupies a pivotal position in deep clustering, unifying self-supervised instance discrimination @oord2018infonce with one-stage clustering methods such as IIC @ji2019iic, SwAV @caron2020swav, and DAC @chang2017dac, while eliminating the alternating offline $k$-Means re-clustering of DeepCluster @caron2018deepcluster and the reconstruction bottleneck of DEC / IDEC @xie2016dec @guo2017idec and DCN @yang2017dcn. Its high-temperature limit reveals an exact algebraic bridge to spherical clustering @dhillon2001concept @banerjee2005vmf and normalized spectral relaxations @shi2000ncut.

#proposition(title: "High-Temperature Reduction to Spherical Dispersion and Normalized Cut Gram Alignment", origin: "Original Derivation")[
  In the zero-augmentation-noise limit ($z_i^a = z_i^b = z_i in SS^(m-1)$) as $tau_I -> oo$, the instance contrastive loss admits the exact asymptotic expansion
  $
    L_"ins" = log(2N - 1) + (2N) / (tau_I (2N - 1)) (norm(overline(z))_2^2 - 1) + O(tau_I^(-2)) ,
  $
  where $overline(z) = (1 / N) sum_(i=1)^N z_i$ is the batch mean embedding. Simultaneously, for $bold(Y)^a = bold(Y)^b = bold(Y)$, the cluster contrastive loss minimizes the off-diagonal entries of the normalized Gram matrix $bold(D)^(-1/2) bold(Y)^top bold(Y) bold(D)^(-1/2)$ with $bold(D) = "diag"(norm(hat(y)_1)_2^2, dots, norm(hat(y)_K)_2^2)$.
]

#proof[
  Set $epsilon = 1 / tau_I$ and $s_(i j) = chevron.l z_i, z_j chevron.r$, noting $s_(i i) = 1$. Since $z_i^a = z_i^b = z_i$,
  $
    ell_i^a &= log(1 + 2 sum_(j != i) exp(epsilon (s_(i j) - 1))) \
            &= log((2N - 1) [1 + (2 epsilon) / (2N - 1) sum_(j != i) (s_(i j) - 1) + O(epsilon^2)]) .
  $
  Expanding $log(1 + u) = u + O(u^2)$ and averaging over $i = 1, dots, N$ yields
  $
    L_"ins" = log(2N - 1) + (2 epsilon) / (2N - 1) [ 1 / N sum_(i=1)^N sum_(j != i) s_(i j) - (N - 1) ] + O(epsilon^2) .
  $
  Using $sum_(i=1)^N sum_(j != i) chevron.l z_i, z_j chevron.r = norm(sum_(i=1)^N z_i)_2^2 - sum_(i=1)^N norm(z_i)_2^2 = N^2 norm(overline(z))_2^2 - N$ gives $(1/N)(N^2 norm(overline(z))_2^2 - N) - (N - 1) = N(norm(overline(z))_2^2 - 1)$, proving the expansion. Minimizing $norm(overline(z))_2^2$ drives the resultant mean vector to zero—matching maximum spherical dispersion on $SS^(m-1)$ @banerjee2005vmf. For the cluster head, $s(hat(y)_k, hat(y)_l) = [bold(D)^(-1/2) bold(Y)^top bold(Y) bold(D)^(-1/2)]_(k l)$, so minimizing $hat(ell)_k$ directly penalizes off-diagonal mass of the column-normalized Gram matrix, enforcing the soft orthogonality constraint of normalized spectral clustering @ng2002spectral @shi2000ncut.
]

*Empirical Benchmark Analysis:* On the standardized ClusterOverview four-dataset benchmark (evaluating Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), Contrastive Clustering achieves Iris Sil = 0.440, ARI = 0.613; Customers Sil = 0.299; Multi-Shape Sil = 0.611 (ARI = N/A); and High-D Sil = 0.286 (ARI = N/A), at a mean runtime of 9.230 ms. These numbers reflect the theory developed above: CC delivers solid high-dimensional separation (High-D Sil = 0.286, outperforming CluStream's 0.253 and avoiding DenStream's high-dimensional collapse) because spherical projection and nonlinear MLP encoding counter Euclidean distance concentration @beyer1999nearest, while achieving strong Iris ARI (0.613) and Multi-Shape Sil (0.611). Conversely, its Customers Silhouette (0.299) is more modest due to the uniform-prior bias of $-H(bold(Y))$ on skewed commercial segments, and its 9.230 ms mean runtime is 4--17$times$ slower than single-pass streaming algorithms (Algs. 193--195) owing to iterative dual-head gradient updates.
