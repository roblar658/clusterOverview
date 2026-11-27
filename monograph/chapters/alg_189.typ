#import "../template/ieee.typ": *

#alg_entry(
  id: 189,
  name: "SCAN (Semantic Clustering by Adopting Nearest neighbors)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2020",
  paradigm: "Decoupled Self-Supervised Neighbor-Consistency Clustering",
  complexity: $O(N k log N + E N (C_"net" + k K))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

End-to-end deep clustering methods such as DeepCluster @caron2018deepcluster, DEC @xie2016dec, and JULE @yang2016jule suffer from circular coupling: noisy early features produce corrupted pseudo-labels that steer representations toward low-level shortcuts. *Semantic Clustering by Adopting Nearest neighbors (SCAN)* @vangansbeke2020scan breaks this coupling via a decoupled multi-stage pipeline.

Given $cal(X) = {bold(x)_1, ..., bold(x)_N} subset.eq bb(R)^d$, Stage 1 trains a self-supervised pretext encoder $Phi_phi : bb(R)^d -> bb(S)^(r-1)$ (e.g., InfoNCE @oord2018infonce) invariant to augmentation group $cal(T)$, yielding $bold(z)_i = Phi_phi (bold(x)_i)$. Rather than running Euclidean K-Means globally, SCAN mines the directed $k$-nearest neighbor set $cal(N)_k (i) subset.eq {1, ..., N} without {i}$ of each point in embedding space, leveraging high local semantic purity even on non-convex manifolds. Stage 2 then trains a softmax clustering head $Phi_theta : bb(R)^d -> Delta^(K-1)$ ($bold(p)_i = Phi_theta (bold(x)_i)$) enforcing consistent, confident cluster predictions across mined neighbors while avoiding collapse through marginal entropy regularization.

#definition(title: "Semantic Neighbor-Consistency Objective and Self-Labeling Stage")[
  Let $cal(N)_k (i)$ be the Stage-1 $k$-NN set of $bold(x)_i$, $bold(p)_i = Phi_theta (bold(x)_i) in Delta^(K-1)$, and $bar(bold(p)) = (1 / N) sum_(i=1)^N bold(p)_i$. For entropy weight $lambda > 0$, the normalized *SCAN Stage-2 objective* is
  $
    cal(L)_"SCAN" (theta) &= - 1 / (N k) sum_(i=1)^N sum_(j in cal(N)_k (i)) ln chevron.l bold(p)_i, bold(p)_j chevron.r \
      &quad - lambda H(bar(bold(p))),
  $
  where $H(bar(bold(p))) = - sum_(c=1)^K bar(p)_c ln bar(p)_c$ (absorbing $1 / k$ into $lambda$ recovers the unnormalized form of @vangansbeke2020scan). In *Stage-3 self-labeling*, confident points $cal(S)_tau = {i : max_c p_(i, c) > tau}$ receive hard pseudo-labels $hat(y)_i = arg max_c p_(i, c)$ to fine-tune $Phi_theta$ via inverse-frequency weighted cross-entropy under strong augmentations $cal(T)_"strong"$.
]

Because $chevron.l bold(p)_i, bold(p)_j chevron.r$ equals the probability that independent draws from $bold(p)_i$ and $bold(p)_j$ agree, minimizing $-ln chevron.l bold(p)_i, bold(p)_j chevron.r$ jointly drives *neighbor agreement* ($bold(p)_i approx bold(p)_j$) and *one-hot sharpness* ($||bold(p)_i||_2 -> 1$).

#algbox(
  title: "SCAN: Decoupled Semantic Neighbor-Consistency Clustering",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$; clusters $K$; neighbors $k$; weight $lambda > 0$; threshold $tau in (1 / K, 1)$; augmentations $cal(T), cal(T)_"strong"$; epochs $E_1, E_2, E_3$.],
  outputs: [Clustering head $Phi_theta$ and hard labels $hat(y)_i = arg max_c [Phi_theta (bold(x)_i)]_c$.],
)[
  + *Stage 1 (Pretext & $k$-NN Mining):* Train $Phi_phi : bb(R)^d -> bb(S)^(r-1)$ for $E_1$ epochs under $cal(T)$ @oord2018infonce; embed $bold(z)_i = Phi_phi (bold(x)_i)$ and index $cal(N)_k (i)$.
  + *Stage 2 (Neighbor-Consistency Optimization):* Initialize $Phi_theta$ from $phi$ with a $K$-way softmax head. *For* $e = 1, ..., E_2$:
    + For each minibatch $cal(B)$, sample $j in cal(N)_k (i)$ (or traverse all $k$ neighbors), compute $bold(p)_i = Phi_theta (g(bold(x)_i))$, $bold(p)_j = Phi_theta (g'(bold(x)_j))$, and batch prior $bar(bold(p))_cal(B)$.
    + Descend $-(1 / |cal(B)|) sum_(i in cal(B)) ln chevron.l bold(p)_i, bold(p)_j chevron.r - lambda H(bar(bold(p))_cal(B))$ via SGD.
  + *Stage 3 (Self-Labeling):* *For* $e = 1, ..., E_3$, form $cal(S)_tau = {i : max_c p_(i, c) > tau}$ with $hat(y)_i = arg max_c p_(i, c)$ and inverse-frequency weights $w_c$, and minimize weighted CE on $(g(bold(x)_i), hat(y)_i)$ ($g tilde.op cal(T)_"strong"$).
  + *Return* $Phi_theta$ and $hat(y)_i = arg max_c [Phi_theta (bold(x)_i)]_c$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Fano Collapse Threshold and Intrinsic Confidence Calibration")

We determine the exact entropy weight $lambda$ preventing single-cluster collapse ($bold(p)_i equiv bold(e)_1$, yielding $cal(L)_"SCAN" = 0$) when fraction $eta > 0$ of $k$-NN edges cross true semantic class boundaries. Remarkably, symmetric optimization admits an exact closed-form solution equal to the Fano conditional entropy bound.

#theorem(title: "Exact Fano Collapse Threshold and Stage-2 Uncertainty Calibration", origin: "Original Proof")[
  Let $cal(X)$ partition into $K >= 2$ balanced classes with $k$-NN impurity rate $eta = (1 / (N k)) |{(i, j) : j in cal(N)_k (i), y_i != y_j}| in (0, 1 - 1 / K)$. Over symmetric class-aligned assignments $p_(i, y_i) = alpha in [1 / K, 1]$ and $p_(i, c) = (1 - alpha) / (K - 1)$ ($c != y_i$):
  + The unique optimal confidence $alpha^*$ minimizing $cal(L)_"SCAN"$ is strictly interior ($alpha^* < 1$ when $eta > 0$):
    $
      alpha^* = 1 / K + 1 / K sqrt((K - 1) (K (1 - eta) - 1)).
    $
  + The minimal Stage-2 loss equals the Fano-regularized value
    $
      cal(L)_"SCAN"^* (eta) = h_2 (eta) + eta ln(K - 1) - lambda ln K,
    $
    where $h_2 (eta) = -eta ln eta - (1 - eta) ln(1 - eta)$.
  + The true partition strictly beats single-cluster collapse ($cal(L)_"SCAN" = 0$) iff
    $
      lambda > lambda_"crit" (eta, K) := (h_2 (eta) + eta ln(K - 1)) / (ln K) in (0, 1).
    $
]

#proof[
  *Step 1 (Inner-product parameterization).* Set $alpha = 1 / K + delta$ and $beta = 1 / K - delta / (K - 1)$ for $delta in [0, 1 - 1 / K]$. Balance and symmetry give $bar(p)_c = 1 / K$ and $-lambda H(bar(bold(p))) = -lambda ln K$. Same-class pairs satisfy
  $
    s_+ (delta) &= alpha^2 + (K - 1) beta^2 = 1 / K + K / (K - 1) delta^2,
  $
  while $s_+ (delta) + (K - 1) s_- (delta) = (alpha + (K - 1) beta)^2 = 1$ gives cross-class inner products
  $
    s_- (delta) = (1 - s_+ (delta)) / (K - 1) = 1 / K - K / (K - 1)^2 delta^2.
  $
  Setting $u = (K^2 / (K - 1)) delta^2 in [0, K - 1]$ yields $s_+ = (1 + u) / K$ and $s_- = (1 - u / (K - 1)) / K$.

  *Step 2 (Exact minimization).* Substituting into $cal(L)_"SCAN"$ gives
  $
    cal(L)_"SCAN" (u) &= (1 - lambda) ln K - (1 - eta) ln(1 + u) \
      &quad - eta ln(1 - u / (K - 1)),
  $
  which is strictly convex on $[0, K - 1)$ and diverges as $u -> K - 1$ for $eta > 0$. Setting the derivative to zero gives
  $
    (1 - eta) / (1 + u) = eta / (K - 1 - u) <=> u^* = K(1 - eta) - 1,
  $
  which lies in $(0, K - 1)$ for $eta in (0, 1 - 1 / K)$. Solving $(K^2 / (K - 1))(delta^*)^2 = u^*$ yields $alpha^* = 1 / K + (1 / K) sqrt((K - 1)(K(1 - eta) - 1)) < 1$.

  *Step 3 (Collapse threshold).* At $u^*$, $s_+^* = 1 - eta$ and $s_-^* = eta / (K - 1)$, yielding
  $
    cal(L)_"SCAN"^* (eta) &= -(1 - eta) ln(1 - eta) - eta ln(eta / (K - 1)) - lambda ln K \
      &= h_2 (eta) + eta ln(K - 1) - lambda ln K.
  $
  Since collapse ($bold(p)_i equiv bold(e)_1$) yields $cal(L)_"SCAN" = 0$, the true partition strictly beats collapse iff $lambda > lambda_"crit" (eta, K)$. Concavity of entropy ensures $h_2 (eta) + eta ln(K - 1) < ln K$ for $eta < 1 - 1 / K$, hence $lambda_"crit" in (0, 1)$.
]

Theorem 1 proves both that any $lambda >= 1$ globally prevents collapse whenever $eta < 1 - 1 / K$, and that Stage 2 intrinsically softens confidence to $alpha^* < 1$ under neighbor noise—explaining why Stage-3 self-labeling is needed to sharpen predictions.

#subsec("B", "Second-Order Spectral Relaxation to Graph Laplacian Clustering")

Expanding $cal(L)_"SCAN"$ around the uniform barycenter exposes an exact link to spectral graph cuts @shi2000ncut @vonluxburg2007tutorial and a local symmetry-breaking threshold.

#lemma(title: "Spectral Laplacian Relaxation and Local Symmetry-Breaking Threshold", origin: "Original Proof")[
  Let $bold(p)_i = (1 / K) bold(1)_K + bold(v)_i$ ($bold(1)_K^top bold(v)_i = 0$, $max_i ||bold(v)_i||_2 <= epsilon$), and stack $bold(v)_i^top$ into $bold(V) in bb(R)^(N times K)$. Let $bold(A) = (bold(W) + bold(W)^top) / 2$ with $W_(i j) = (1 / k) bb(1)_(j in cal(N)_k (i))$, $bold(D) = "diag"(bold(A) bold(1)_N)$, and $bold(L) = bold(D) - bold(A)$. Then
  $
    cal(L)_"SCAN" (bold(V)) &= (1 - lambda) ln K + K / N "Tr"(bold(V)^top bold(L) bold(V)) \
      &quad - K / N "Tr"(bold(V)^top bold(D) bold(V)) \
      &quad + (lambda K) / (2 N^2) ||bold(V)^top bold(1)_N||_2^2 + O(epsilon^3).
  $
  On a $k$-regular mutual graph ($bold(D) = bold(I)_N$) with Fiedler eigenvalue $mu_1 = lambda_2 (bold(L)) > 0$, the Fiedler clustering mode descends strictly faster than constant collapse $bold(u)_0 = N^(-1 / 2) bold(1)_N$ iff $lambda > 2 mu_1$, and collapse is locally strictly repulsive iff $lambda > 2$.
]

#proof[
  *Step 1 (Consistency expansion).* Since $bold(1)_K^top bold(v)_i = 0$, $chevron.l bold(p)_i, bold(p)_j chevron.r = 1 / K + chevron.l bold(v)_i, bold(v)_j chevron.r$, so $-ln chevron.l bold(p)_i, bold(p)_j chevron.r = ln K - K chevron.l bold(v)_i, bold(v)_j chevron.r + O(epsilon^4)$. Averaging over edges and using symmetry of $chevron.l bold(v)_i, bold(v)_j chevron.r$ replaces $bold(W)$ with $bold(A) = bold(D) - bold(L)$:
  $
    - 1 / (N k) sum_(i=1)^N sum_(j in cal(N)_k (i)) ln chevron.l bold(p)_i, bold(p)_j chevron.r
    &= ln K + K / N "Tr"(bold(V)^top bold(L) bold(V)) \
    &quad - K / N "Tr"(bold(V)^top bold(D) bold(V)) + O(epsilon^4).
  $

  *Step 2 (Entropy expansion).* With $bar(bold(v)) = (1 / N) bold(V)^top bold(1)_N$ and $sum_c bar(v)_c = 0$, expanding $(1 / K + bar(v)_c) ln(1 / K + bar(v)_c) = -(1 / K) ln K + bar(v)_c (1 - ln K) + (K / 2) bar(v)_c^2 + O(epsilon^3)$ and summing over $c$ yields
  $
    -lambda H(bar(bold(p))) &= -lambda ln K + (lambda K) / 2 ||bar(bold(v))||_2^2 + O(epsilon^3),
  $
  which equals $-lambda ln K + (lambda K / (2 N^2)) ||bold(V)^top bold(1)_N||_2^2 + O(epsilon^3)$. Adding Steps 1 and 2 gives the expansion.

  *Step 3 (Curvature threshold).* For $bold(D) = bold(I)_N$ and rank-1 perturbation $bold(V) = epsilon bold(u) bold(d)^top$ ($||bold(u)||_2 = ||bold(d)||_2 = 1, bold(1)_K^top bold(d) = 0$), the collapse mode $bold(u)_0 = N^(-1 / 2) bold(1)_N$ has $bold(L) bold(u)_0 = bold(0)$ and $||bold(V)^top bold(1)_N||_2^2 = epsilon^2 N$, giving quadratic coefficient $(K / N)(lambda / 2 - 1)$. The Fiedler mode $bold(u)_1$ ($bold(u)_1^top bold(1)_N = 0, bold(L) bold(u)_1 = mu_1 bold(u)_1$) has coefficient $(K / N)(mu_1 - 1)$. Thus $mu_1 - 1 < lambda / 2 - 1 <=> lambda > 2 mu_1$, while collapse has positive curvature iff $lambda > 2$.
]

#subsec("C", "Pretext Representation Bayes Error and k-NN Purity Guarantees")

Because Stage 2 optimizes over fixed neighborhoods $cal(N)_k (i)$, its recoverability depends directly on the asymptotic $k$-NN impurity of the Stage-1 pretext space @vangansbeke2020scan @beyer1999nearest.

#proposition(title: "Asymptotic Neighbor Impurity Bound via Pretext Bayes Risk", origin: "Classical — Cited")[
  Let $p_c (bold(z)) = "Pr"(Y = c | bold(z))$ be the class posterior in Stage-1 embedding space and $R_"Bayes" = bb(E)[1 - max_c p_c (bold(z))]$ the Bayes error. Under density regularity as $N -> infinity$ with $k = o(N)$, the limiting $k$-NN impurity $eta_infinity = lim_(N -> infinity) bb(E)[eta]$ satisfies @vangansbeke2020scan @beyer1999nearest
  $
    R_"Bayes" <= eta_infinity <= 2 R_"Bayes" - K / (K - 1) R_"Bayes"^2.
  $
  Under balanced isotropic Gaussian class embeddings $cal(N)(bold(mu)_c, sigma^2 bold(I)_r)$ with minimum separation $Delta = min_(c != c') ||bold(mu)_c - bold(mu)_(c')||_2$, we obtain $eta_infinity <= 2 (K - 1) exp(-Delta^2 / (8 sigma^2))$.
]

#proof[
  As $N -> infinity$ with $k = o(N)$, neighbors of $bold(z)$ converge to $bold(z)$ in probability, so neighbor labels draw asymptotically independently from $P(Y | bold(z))$, giving $eta_infinity = bb(E)[1 - sum_(c=1)^K p_c (bold(z))^2]$. Setting $e(bold(z)) = 1 - p_"max" (bold(z))$, we have $sum_c p_c^2 <= p_"max" sum_c p_c = p_"max"$, yielding $eta_infinity >= R_"Bayes"$. Conversely, Cauchy--Schwarz on the $K - 1$ non-maximal posteriors gives $sum_(c != c^*) p_c^2 >= e(bold(z))^2 / (K - 1)$, whence
  $
    1 - sum_(c=1)^K p_c (bold(z))^2 &<= 1 - (1 - e(bold(z)))^2 - e(bold(z))^2 / (K - 1) \
      &= 2 e(bold(z)) - K / (K - 1) e(bold(z))^2.
  $
  Taking expectations and invoking Jensen's inequality $bb(E)[e^2] >= R_"Bayes"^2$ yields the upper bound. Under balanced Gaussian embeddings, pairwise Bayes error is bounded by $exp(-Delta^2 / (8 sigma^2))$ via the Chernoff bound; a union bound over $K - 1$ competing classes gives $R_"Bayes" <= (K - 1) exp(-Delta^2 / (8 sigma^2))$.
]

#subsec("D", "Stage-3 Error Contraction and Pretext-Misalignment Failure")

While Stage-3 self-labeling sharpens calibrated Stage-2 probabilities into accurate hard clusters, it irreversibly amplifies errors whenever the Stage-1 pretext task misaligns with target semantics.

#theorem(title: "Stage-3 Orbit Contraction and Pretext-Misalignment Counterexample", origin: "Original Counterexample")[
  + *(Stage-3 Contraction).* If $cal(S)_tau$ has class-symmetric pseudo-label error $epsilon_tau = "Pr"(hat(y) != Y | bold(x) in cal(S)_tau) < 1 - 1 / K$ and class-preserving augmentations $cal(T)_"strong"$ mix orbits across $cal(S)_tau$, then the Stage-3 population minimizer $p^*$ satisfies $arg max_c p^* (c | bold(x)) = Y(bold(x))$ identically.
  + *(Irreversible Pretext Misalignment).* Let $K = 2$ with balanced $(y, s) in {-1, +1}^2$ ($y$ semantic, $s$ nuisance). If Stage-1 embeds $Phi_phi (bold(x)) = (delta y, M s) + bold(xi)$ with $bold(xi) tilde.op cal(N)(bold(0), sigma^2 bold(I)_2)$ and $M >> sigma >> delta$, then $eta_"sem" -> 1 / 2$ while $eta_"nuis" <= 2 exp(-M^2 / (2 sigma^2))$. For every $lambda > 0$, Stage 2 uniquely prefers the nuisance partition by margin $ln 2$ with confidence $alpha^* -> 1$, causing Stage 3 to lock in $50\%$ semantic error.
]

#proof[
  *Part 1.* Minimizing $bb(E)_(bold(x) in cal(S)_tau, g tilde.op cal(T)_"strong") [-ln p(hat(y)(bold(x)) | g(bold(x)))]$ over $cal(T)_"strong"$-invariant posteriors yields by Gibbs' inequality $p^*(c | bold(x)') = "Pr"(hat(y) = c | bold(x)' in cal(O)(bold(x)), bold(x) in cal(S)_tau)$. Because $cal(T)_"strong"$ preserves class $y$, $p^*(y | bold(x)') = 1 - epsilon_tau > epsilon_tau / (K - 1) = p^*(c | bold(x)')$ ($c != y$) whenever $epsilon_tau < 1 - 1 / K$, giving zero Bayes error.

  *Part 2.* Under $Phi_phi (bold(x)) = (delta y, M s) + bold(xi)$, nuisance classes $s = plus.minus 1$ have separation $2 M$, so Proposition 1 gives $eta_"nuis" <= 2 exp(-M^2 / (2 sigma^2)) -> 0$. Semantic classes $y = plus.minus 1$ have separation $2 delta << sigma$, so $p_(y = +1) (bold(z)) -> 1 / 2$ a.s. and $eta_"sem" -> 1 / 2$. By Theorem 1, $cal(L)^*_"SCAN" (eta_"nuis") -> -lambda ln 2$, whereas $cal(L)^*_"SCAN" (eta_"sem") -> (1 - lambda) ln 2$, strictly favoring the nuisance clustering by margin $ln 2$ with confidence $alpha^* = (1 / 2)(1 + sqrt(1 - 2 eta_"nuis")) -> 1 > tau$. Every sample enters $cal(S)_tau$ carrying nuisance label $s$, and Stage 3 locks in $50\%$ semantic error.
]

#sec("III", "Complexity, Geometry and Operational Profile")

SCAN executes in three phases: Stage 1 trains the pretext encoder in $O(E_1 N C_"net")$ time and builds the $k$-NN graph in $O(N k log N)$ time via approximate indexing (HNSW/FAISS) or $O(N^2 r)$ exact search; Stage 2 optimizes neighbor consistency in $O(E_2 N (C_"net" + k K))$ time; Stage 3 fine-tunes on $|cal(S)_tau| <= N$ confident samples in $O(E_3 N C_"net")$ time, totaling $O(N k log N + E N (C_"net" + k K))$ ($E = E_1 + E_2 + E_3$). Space complexity is $O(|theta| + N r + N k)$ during Stage 1 and $O(|theta| + N k + B K)$ during Stage 2. Geometrically, Lemma 1 shows SCAN is a parametric neural generalization of normalized spectral clustering @vonluxburg2007tutorial on the pretext $k$-NN manifold, scaling linearly in $K$ (avoiding IIC's $O(K^2)$ joint matrix) yet remaining sensitive to high-dimensional distance concentration @beyer1999nearest.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N k log N + E N (C_"net" + k K))$ ($O(N^2 r)$ exact $k$-NN); $O(N C_"net")$ inference.],
    [Space Complexity], [$O(|theta| + N r + N k)$ in Stage 1; $O(|theta| + N k + B K)$ in Stage 2.],
    [Inductive Bias], [Manifold smoothness on the Stage-1 $k$-NN graph with balanced cluster priors.],
    [Collapse Resistance], [Global for $lambda > lambda_"crit" (eta, K)$ (Theorem II.A); local for $lambda > 2 mu_1$ (Lemma II.B).],
    [Confidence Calibration], [Stage 2 self-calibrates to $alpha^* < 1$ under $eta > 0$; Stage 3 sharpens via $tau$.],
    [Scalability in $K$], [Linear $O(K)$ per edge; no $O(K^2)$ matrix or $B >> K^2$ batch constraint.],
    [Hyperparameters], [Pretext $Phi_phi$, neighbors $k$, weight $lambda$, threshold $tau$, clusters $K$.],
    [Out-of-Sample Inference], [Direct parametric forward pass $Phi_theta (bold(x))$; no test-time $k$-NN lookup.],
  ),
  caption: [Operational and geometric profile of SCAN (Alg. 189).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Large-Scale Unsupervised Visual Classification:* On CIFAR-10/100, STL-10, and ImageNet subsets, SimCLR/MoCo pretext embeddings @oord2018infonce achieve high local neighbor purity ($eta < 0.15$ at $k = 20$), allowing SCAN to outperform end-to-end DeepCluster and DEC substantially @vangansbeke2020scan.
+ *Fine-Grained Radiology and Pathology Phenotyping:* Self-supervised medical backbones organize tissue tiles on curved scanner-invariant manifolds; Stage-2 neighbor adoption respects non-convex boundaries while Stage-3 thresholding filters ambiguous lesions.
+ *Single-Cell Transcriptomics Cell Typing:* Contrastive gene-expression embeddings place developmental trajectories on low-dimensional curved graphs where local $k$-NN agreement preserves cell identity better than Euclidean centroids while offering calibrated Stage-2 transition uncertainty.
+ *Unsupervised Document and Multimodal Intent Discovery:* Pretrained sentence and vision-language encoders already supply semantically pure Stage-1 neighborhoods, which SCAN converts into discrete intent clusters without destabilizing pretrained weights.
+ *Industrial Optical Defect Taxonomy Discovery:* Decoupling illumination-invariant self-supervised pretraining from neighbor adoption separates subtle surface defects from lighting variation on millions of unlabeled inspection crops.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Three failure modes dominate: (1) *pretext-semantic misalignment* (Theorem II.D) locks Stages 2 and 3 onto nuisance clusters irreversibly; (2) *high-dimensional distance concentration and hubness* @beyer1999nearest raise $eta$ when $r$ is large relative to $N$; and (3) *extreme class imbalance* conflicts with $-lambda H(bar(bold(p)))$ unless replaced by a non-uniform prior KL penalty @vangansbeke2020scan.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Natural Images with SSL Backbone], [Excellent: high Stage-1 $k$-NN purity yields strong semantic clustering.], [Verify pretext invariances align with labels.],
    [Pretrained Text / Multimodal Spaces], [Strong: exploits frozen semantic neighborhood geometry directly.], [Choose $k$ smaller than rarest class size.],
    [Tabular / Low-D Non-Convex Data], [Strong: local $k$-NN graph captures non-convex shapes cleanly.], [Higher runtime than IIC from $k$-NN indexing.],
    [High-D Noisy Spaces ($r >> log N$)], [Moderate-Weak: $k$-NN hubness and concentration inflate $eta$.], [Reduce dimension; compare IIC or SwAV.],
    [Imbalanced / Streaming Data], [Weak: static $k$-NN graph and uniform entropy penalize rare modes.], [Use non-uniform prior KL or online SwAV.],
  ),
  caption: [Domain suitability and failure-mode matrix for SCAN (Alg. 189).],
)

#sec("V", "Structural Relations and Empirical Behavior")

SCAN links deep clustering to Normalized Cuts @shi2000ncut, Spectral Clustering @ng2002spectral @vonluxburg2007tutorial, and FINCH @sarfraz2019finch, while generalizing IIC @ji2019iic from augmentation pairs to mined $k$-NN pairs.

#proposition(title: "Variational Reduction Between SCAN, IIC, and Soft Spectral Cuts", origin: "Original Derivation")[
  Define the $k$-NN paired cluster matrix $bold(P)_"knn" = (1 / (N k)) sum_(i=1)^N sum_(j in cal(N)_k (i)) bold(p)_i bold(p)_j^top$.
  + By Jensen's inequality, Stage-2 consistency upper-bounds the negative log-trace of $bold(P)_"knn"$:
    $
      - 1 / (N k) sum_(i=1)^N sum_(j in cal(N)_k (i)) ln chevron.l bold(p)_i, bold(p)_j chevron.r
      >= - ln "Tr"(bold(P)_"knn").
    $
  + Maximizing $"Tr"(bold(P)_"knn")$ with $-lambda H(bar(bold(p)))$ penalizes off-diagonal cut mass while driving rows to simplex vertices; replacing $k$-NN pairs by augmentation pairs and substituting $-I_lambda (bold(P)_"aug")$ recovers IIC @ji2019iic.
]

#proof[
  Convexity of $t |-> -ln t$ on $(0, infinity)$ applied to the uniform average over all $N k$ edges yields
  $
    - 1 / (N k) sum_(i, j) ln chevron.l bold(p)_i, bold(p)_j chevron.r
    >= - ln ( 1 / (N k) sum_(i, j) chevron.l bold(p)_i, bold(p)_j chevron.r )
    = - ln "Tr"(bold(P)_"knn"),
  $
  using $"Tr"(bold(P)_"knn") = (1 / (N k)) sum_(i, j) chevron.l bold(p)_i, bold(p)_j chevron.r$. Since $chevron.l bold(p)_i, bold(p)_j chevron.r <= ||bold(p)_i||_2 ||bold(p)_j||_2 <= 1$, driving $"Tr"(bold(P)_"knn") -> 1$ forces $bold(p)_i = bold(p)_j$ on edges and one-hot sharpness $||bold(p)_i||_2 -> 1$, while $-lambda H(bar(bold(p)))$ prevents collapse just as IIC's marginal entropy does on augmentation pairs @ji2019iic.
]

*Empirical Behavior on the ClusterOverview Benchmark.* On the 4-dataset ClusterOverview benchmark @rousseeuw1987silhouettes @hubert1985comparing, SCAN achieves *Iris* $"Sil" = 0.430$, $"ARI" = 0.632$; *Customers* $"Sil" = 0.319$; *Multi-Shape* $"Sil" = 0.709$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.211$ ($"ARI" =$ N/A), at mean runtime *7.005 ms*. SCAN achieves the *highest Iris ARI ($0.632$)*, *highest Iris Silhouette ($0.430$)*, and *highest Customers Silhouette ($0.319$)* across Algs. 188--191, validating decoupled neighbor mining on low- and moderate-dimensional manifolds. Conversely, its *High-D* Silhouette drops to $0.211$ (behind IIC's $0.298$ and SwAV's $0.292$), corroborating Proposition II.C and @beyer1999nearest as $k$-NN distance concentration elevates $eta$. Its $7.005$ ms runtime is the slowest of the quartet due to explicit $k$-NN indexing and multi-stage training.
