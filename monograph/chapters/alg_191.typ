#import "../template/ieee.typ": *

#alg_entry(
  id: 191,
  name: "DeepCluster",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2018",
  paradigm: "Alternating Deep Feature K-Means and Pseudo-Label Training",
  complexity: $O(E (N C_"net" + I_"km" N K r'))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Self-supervised representation learning historically relied on domain-specific pretext tasks. *DeepCluster* @caron2018deepcluster demonstrated that coupling a standard neural encoder with classical K-Means @lloyd1982least learns transferable visual representations at ImageNet scale without handcrafted pretext objectives, by alternating *offline* global feature clustering with *discriminative* training on the resulting cluster indices as pseudo-labels.

Let $cal(X) = {bold(x)_1, ..., bold(x)_N} subset.eq bb(R)^d$, encoder $f_theta : bb(R)^d -> bb(R)^r$, and linear classifier head $(bold(W), bold(b)) in bb(R)^(K times r) times bb(R)^K$ giving $p_(theta, bold(W)) (k | bold(x)_i) = ["softmax"(bold(W) f_theta (bold(x)_i) + bold(b))]_k$. At epoch $t$, DeepCluster extracts $bold(z)_i = f_(theta^((t-1))) (bold(x)_i)$, applies PCA whitening $bold(U)_(r') = bold(V)_(r') bold(Lambda)_(r')^(-1 / 2)$ ($r' <= r$), and $ell_2$-normalizes onto $bb(S)^(r'-1)$:
$
  tilde(bold(z))_i = (bold(U)_(r')^top (bold(z)_i - bar(bold(z)))) / (||bold(U)_(r')^top (bold(z)_i - bar(bold(z)))||_2).
$

#definition(title: "Alternating K-Means Pseudo-Labeling and Balanced Discriminative Objective")[
  Across epochs $t = 1, ..., E$, DeepCluster alternates two subproblems:
  + *Offline Clustering Step:* Fixing $theta^((t-1))$, solve Euclidean K-Means on ${tilde(bold(z))_i}$ for centroids $bold(C)^((t)) in bb(R)^(r' times K)$ and hard labels $bold(Y)^((t)) = {bold(y)_i^((t))} subset.eq {0, 1}^K$ ($bold(1)_K^top bold(y)_i^((t)) = 1$):
    $
      min_(bold(C), bold(Y)) cal(J)_"KM" (bold(C), bold(Y); theta^((t-1)))
      := 1 / N sum_(i=1)^N ||tilde(bold(z))_i - bold(C) bold(y)_i||_2^2,
    $
    splitting a non-empty cluster $bold(c)_(k^*) plus.minus bold(delta)$ whenever any cluster empties.
  + *Balanced Discriminative Step:* Reset $(bold(W), bold(b))$ and update $(theta, bold(W), bold(b))$ via SGD minimizing *cluster-uniformly reweighted* cross-entropy:
    $
      cal(J)_"CE"^"bal" (theta, bold(W), bold(b); bold(Y)^((t)))
      &= 1 / K sum_(k=1)^K 1 / (|cal(C)_k|) \
      &quad times sum_(i in cal(C)_k) ell_"CE" (bold(W) f_theta (bold(x)_i) + bold(b), k),
    $
    implemented by sampling batches uniformly across clusters @caron2018deepcluster.
]

Empty-cluster splitting, uniform cluster resampling, and epoch-wise classifier resets are mathematically indispensable safeguards, as proved below.

#algbox(
  title: "DeepCluster: Alternating K-Means and Discriminative Pseudo-Labeling",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$; clusters $K$ (typically $K > K_"true"$); PCA dim $r'$; Lloyd iterations $I_"km"$; batch $B$; step $eta$; epochs $E$.],
  outputs: [Encoder $f_theta$, centroids $bold(C)^((E))$, hard assignments $hat(y)_i in {1, ..., K}$.],
)[
  + Initialize encoder $f_theta$ (with Sobel preprocessing for raw images @caron2018deepcluster).
  + *For* epoch $t = 1, ..., E$ *do*:
    + Extract $bold(z)_i <- f_(theta^((t-1))) (bold(x)_i)$ ($i = 1, ..., N$), fit PCA whitening $bold(U)_(r')$, and form $tilde(bold(z))_i in bb(S)^(r'-1)$.
    + Run $I_"km"$ Lloyd steps @lloyd1982least (with K-Means++ seeding @arthur2007kmeanspp and empty-cluster splitting) to obtain $bold(C)^((t))$ and ${cal(C)_k^((t))}_(k=1)^K$.
    + Re-initialize linear classifier $(bold(W), bold(b))$ to break inter-epoch index permutation mismatch.
    + Draw minibatches uniformly over clusters $k in {1, ..., K}$ and update $(theta, bold(W), bold(b))$ via SGD on $cal(J)_"CE"^"bal"$ under augmentations $cal(T)$.
  + *Return* $f_(theta^((E)))$, $bold(C)^((E))$, and final pseudo-labels $hat(y)_i$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Absence of a Joint Lyapunov Functional and Whitening Degeneracy")

Unlike DCN @yang2017dcn or Lloyd's algorithm @lloyd1982least @bottou1995convergence, DeepCluster alternates incompatible objectives ($cal(J)_"KM"$ on PCA-whitened normalized features vs. $cal(J)_"CE"^"bal"$ on raw features). We prove constructively that no joint monotonic Lyapunov descent exists and that PCA whitening erases linear discriminative scaling.

#theorem(title: "Non-Monotonicity of K-Means Distortion and Whitening-Induced Partition Flip", origin: "Original Counterexample")[
  Consider $cal(X) = {(plus.minus 1, plus.minus 1)} subset.eq bb(R)^2$ with $K = 2$ and diagonal encoder $f_theta (bold(x)) = "diag"(a_1, a_2) bold(x)$ ($a_1, a_2 > 0$).
  + *(Strict CE Descent with Strict K-Means Increase).* Under $ell_2$ normalization without whitening, start at $(a_1, a_2) = (1, 1)$ with horizontal partition $cal(C)_1 = {x_1 = +1}, cal(C)_2 = {x_1 = -1}$ and $bold(W) = ((w, 0), (-w, 0))$ ($w > 0$). Updating to $(a'_1, a'_2) = (2, 4)$ strictly decreases cross-entropy from $ln(1 + e^(-2 w))$ to $ln(1 + e^(-4 w))$, yet strictly increases horizontal K-Means distortion from $2.0$ to $3.2$ and flips the optimal K-Means partition to vertical (distortion $0.8 < 3.2$).
  + *(Whitening Erasure).* Under PCA whitening + $ell_2$ normalization, for *every* diagonal encoder $(a_1, a_2) > 0$—even as $a_1 / a_2 -> infinity$—the whitened embeddings stay fixed at $tilde(bold(z))_i = (plus.minus 1 / sqrt(2), plus.minus 1 / sqrt(2))$, erasing all discriminative scaling and leaving horizontal and vertical partitions identically tied at distortion $2.0$.
]

#proof[
  *Part 1.* At $(1, 1)$, normalized features are $(plus.minus 1 / sqrt(2), plus.minus 1 / sqrt(2))$. Horizontal centroids $(plus.minus 1 / sqrt(2), 0)$ incur squared distance $1 / 2$ per point and total distortion $4 times (1 / 2) = 2.0$. Classifier $bold(W) = ((w, 0), (-w, 0))$ gives margin $2 w a_1$ and per-sample loss $ell_"CE" = ln(1 + e^(-2 w a_1))$. At $(a'_1, a'_2) = (2, 4)$, $ell_"CE"$ drops strictly to $ln(1 + e^(-4 w))$ because $ell_"CE"$ is blind to the orthogonal axis $a_2$. Normalized features become $(plus.minus 1 / sqrt(5), plus.minus 2 / sqrt(5))$: horizontal centroids $(plus.minus 1 / sqrt(5), 0)$ incur distortion $4 times (4 / 5) = 3.2 > 2.0$, while vertical centroids $(0, plus.minus 2 / sqrt(5))$ achieve distortion $4 times (1 / 5) = 0.8 < 3.2$, forcing K-Means to flip labels orthogonally (zero epoch-to-epoch ARI).

  *Part 2.* For $bold(z)_i = (plus.minus a_1, plus.minus a_2)$, empirical covariance is $bold(Sigma) = "diag"(a_1^2, a_2^2)$. PCA whitening multiplies by $bold(Sigma)^(-1 / 2) = "diag"(a_1^(-1), a_2^(-1))$, mapping all four points to $(plus.minus 1, plus.minus 1)$ and after normalization to $(plus.minus 1 / sqrt(2), plus.minus 1 / sqrt(2))$ identically, independent of $(a_1, a_2)$.
]

Theorem 1 reveals a core trade-off: PCA whitening blocks rank-1 collapse onto a single dominant axis, yet erases linear boundary scaling so improvement must come from nonlinear manifold folding.

#subsec("B", "Gradient-Domination Collapse and Uniform-Resampling Stabilization")

Without uniform cluster resampling, unweighted ERM induces a self-reinforcing majority-collapse loop through classifier bias and gradient imbalance.

#lemma(title: "Bias-Domination Collapse under ERM and Exact Uniform-Resampling Equilibrium", origin: "Original Proof")[
  Let clusters ${cal(C)_k}$ be non-empty with proportions $pi_k = |cal(C)_k| / N$ and logits $s_(i, k) = bold(w)_k^top f_theta (bold(x)_i) + b_k$.
  + Under *unweighted ERM* $cal(J)_"CE"^"unw"$, every stationary point satisfies $bar(p)_k := (1 / N) sum_i p_(i, k) = pi_k$; at $bold(w)_k = bold(0)$, $b_k - b_(k') = ln(pi_k / pi_(k'))$ and minority clusters receive only $O(pi_k)$ encoder gradient weight.
  + Under *uniform cluster sampling* $cal(J)_"CE"^"bal"$, every stationary point obeys $(1 / K) sum_(c=1)^K bar(p)_(k | c) = 1 / K$ (where $bar(p)_(k | c) = (1 / |cal(C)_c|) sum_(i in cal(C)_c) p_(i, k)$). At $bold(w)_k = bold(0)$, the unique optimum is strictly uniform ($p_k = 1 / K$, $b_1 = ... = b_K$) for *every* imbalanced partition $bold(pi) in "int"(Delta^(K-1))$.
]

#proof[
  *Part 1.* Using $partial ell_"CE" (bold(s)_i, y_i) / (partial b_k) = p_(i, k) - bb(1)_(y_i = k)$ and averaging over $i$ gives
  $
    (partial cal(J)_"CE"^"unw") / (partial b_k) &= 1 / N sum_(i=1)^N (p_(i, k) - bb(1)_(y_i = k)) \
      &= bar(p)_k - pi_k = 0 <=> bar(p)_k = pi_k.
  $
  At $bold(w)_k = bold(0)$, $p_(i, k) = p_k$ is independent of $i$, forcing $exp(b_k) / exp(b_(k')) = pi_k / pi_(k')$; dominant clusters ($pi_1 -> 1$) monopolize $nabla_theta cal(J)_"CE"^"unw"$ and absorb boundary points.

  *Part 2.* Differentiating $cal(J)_"CE"^"bal"$ with respect to $b_k$ yields
  $
    (partial cal(J)_"CE"^"bal") / (partial b_k) &= 1 / K sum_(c=1)^K 1 / (|cal(C)_c|) sum_(i in cal(C)_c) (p_(i, k) - delta_(k c)) \
      &= 1 / K sum_(c=1)^K bar(p)_(k | c) - 1 / K = 0.
  $
  At $bold(w)_k = bold(0)$, $bar(p)_(k | c) = p_k$, uniquely solving to $p_k = 1 / K$ ($b_1 = ... = b_K$) regardless of $bold(pi)$, while each cluster exerts equal $1 / K$ gradient weight on $theta$.
]

#subsec("C", "Spherical Resultant Identity, K-Means Convergence, and Overclustering")

Within each offline step, $ell_2$ normalization endows DeepCluster's K-Means subproblem with an exact hyperspherical resultant geometry @dhillon2001concept @banerjee2005vmf and finite-step Lloyd guarantees @lloyd1982least @arthur2007kmeanspp.

#proposition(title: "Hyperspherical Resultant Equivalence and K-Means Step Guarantees", origin: "Classical — Cited")[
  Let ${tilde(bold(z))_i}_(i=1)^N subset.eq bb(S)^(r'-1)$ be $ell_2$-normalized whitened features @caron2018deepcluster @lloyd1982least @arthur2007kmeanspp @dhillon2001concept.
  + For optimal centroids $bold(c)_k^* = (1 / |cal(C)_k|) sum_(i in cal(C)_k) tilde(bold(z))_i$, Euclidean K-Means distortion satisfies
    $
      cal(J)_"KM" (bold(C)^*, bold(Y)) = 1 - 1 / N sum_(k=1)^K |cal(C)_k| ||bold(c)_k^*||_2^2.
    $
  + Although global minimization is NP-hard @aloise2009np, Lloyd's algorithm strictly decreases $cal(J)_"KM"$ and terminates finitely @bottou1995convergence, with K-Means++ seeding guaranteeing $bb(E)[cal(J)_"KM"] <= 8 (ln K + 2) cal(J)_"KM"^*$ @arthur2007kmeanspp.
]

#proof[
  *Part 1.* Since $||tilde(bold(z))_i||_2^2 = 1$, expanding within-cluster distortion gives
  $
    sum_(i in cal(C)_k) ||tilde(bold(z))_i - bold(c)_k^*||_2^2
    &= sum_(i in cal(C)_k) ( 1 - 2 chevron.l tilde(bold(z))_i, bold(c)_k^* chevron.r + ||bold(c)_k^*||_2^2 ) \
    &= |cal(C)_k| - 2 |cal(C)_k| ||bold(c)_k^*||_2^2 + |cal(C)_k| ||bold(c)_k^*||_2^2 \
    &= |cal(C)_k| (1 - ||bold(c)_k^*||_2^2).
  $
  Summing over $k$ and dividing by $N$ gives the identity, closely paralleling spherical K-Means @dhillon2001concept (which maximizes $sum_k |cal(C)_k| ||bold(c)_k^*||_2$).

  *Part 2.* Alternating exact Voronoi assignment and barycenter updates strictly decreases $cal(J)_"KM"$ on every partition change over the finite set of at most $K^N$ partitions @lloyd1982least @bottou1995convergence; the $O(ln K)$ seeding bound follows from @arthur2007kmeanspp.
]

Overclustering ($K approx 10 K_"true"$ in @caron2018deepcluster) shrinks Voronoi cell diameters on $bb(S)^(r'-1)$, preventing linear Voronoi cells from cutting across curved multi-modal semantic classes.

#subsec("D", "Permutation-Shock Instability and Necessity of Classifier-Head Reset")

Because K-Means cluster labels are defined only up to arbitrary permutations $sigma in S_K$, retaining $bold(W)^((t-1))$ across epochs causes severe loss spikes and spurious backbone gradients even when the geometric partition is unchanged.

#theorem(title: "Permutation-Induced Loss Spike and Gradient Shock Without Head Reset", origin: "Original Proof")[
  Suppose at epoch $t - 1$, $bold(W)^((t-1)) = [bold(w)_1, ..., bold(w)_K]^top$ achieves margin $gamma > 0$ on balanced clusters ($bold(w)_(y_i)^top bold(z)_i - bold(w)_k^top bold(z)_i >= gamma$ for $k != y_i$) with $||bold(w)_k||_2 = R$, $max_(k != k') chevron.l bold(w)_k, bold(w)_(k') chevron.r <= 0$, and $p_(i, y_i) >= 1 - delta$. If epoch-$t$ K-Means returns the *same partition* permuted by $sigma in S_K$ moving fraction $alpha_"perm" = (1 / K) |{k : sigma(k) != k}|$ of clusters:
  + Retaining $bold(W)^((t-1))$ causes initial epoch-$t$ loss to spike to $cal(J)_"CE"^"bal" >= alpha_"perm" ln(1 + e^gamma) > alpha_"perm" gamma$.
  + Every sample $i$ with $sigma(y_i) != y_i$ suffers spurious feature gradient norm
    $
      ||nabla_(bold(z)_i) ell_"CE" (bold(W)^((t-1)) bold(z)_i, sigma(y_i))||_2 >= (sqrt(2) - 2 delta) R.
    $
  + Resetting $bold(W)^((t)) = bold(0)$ restores exact $S_K$-permutation invariance with initial loss $ln K$.
]

#proof[
  *Part 1.* For any $i$ with $sigma(y_i) != y_i$, term $k = y_i$ lies in the denominator of $ell_"CE"$:
  $
    ell_"CE" (bold(W)^((t-1)) bold(z)_i, sigma(y_i))
    &= ln ( 1 + sum_(k != sigma(y_i)) e^((bold(w)_k - bold(w)_(sigma(y_i)))^top bold(z)_i) ) \
    &> ln ( 1 + e^((bold(w)_(y_i) - bold(w)_(sigma(y_i)))^top bold(z)_i) ) \
    &>= ln (1 + e^gamma) > gamma.
  $
  Averaging across the fraction $alpha_"perm"$ of permuted clusters gives $cal(J)_"CE"^"bal" >= alpha_"perm" ln(1 + e^gamma) > alpha_"perm" gamma$.

  *Part 2.* Rewriting the exact feature gradient $nabla_(bold(z)_i) ell_"CE" = sum_k p_(i, k) bold(w)_k - bold(w)_(sigma(y_i))$ gives
  $
    nabla_(bold(z)_i) ell_"CE" &= (bold(w)_(y_i) - bold(w)_(sigma(y_i))) \
      &quad + sum_(k != y_i) p_(i, k) (bold(w)_k - bold(w)_(y_i)).
  $
  Since $chevron.l bold(w)_(y_i), bold(w)_(sigma(y_i)) chevron.r <= 0$ and $||bold(w)_k||_2 = R$, we have $||bold(w)_(y_i) - bold(w)_(sigma(y_i))||_2 = sqrt(2 R^2 - 2 bold(w)_(y_i)^top bold(w)_(sigma(y_i))) >= sqrt(2) R$. Meanwhile $||bold(w)_k - bold(w)_(y_i)||_2 <= 2 R$ and $sum_(k != y_i) p_(i, k) <= delta$, so the second sum has norm at most $2 delta R$. The reverse triangle inequality yields $||nabla_(bold(z)_i) ell_"CE"||_2 >= (sqrt(2) - 2 delta) R$.

  *Part 3.* Setting $bold(W)^((t)) = bold(0)$ and $bold(b)^((t)) = bold(0)$ yields $p_(i, k) = 1 / K$, initial loss $ln K$, and identical SGD trajectories up to column permutation for every $sigma in S_K$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each outer epoch performs full-dataset feature extraction ($O(N C_"net")$), PCA whitening ($O(N r^2 + r^3)$), $I_"km"$ Lloyd iterations ($O(I_"km" N K r')$), and one epoch of balanced SGD ($O(N C_"net")$), totaling $O(E (N C_"net" + I_"km" N K r'))$ (with subsampled PCA) and $O(N (C_"net" + K r'))$ at inference. Space is $O(|theta| + N r' + K r')$ due to caching the full $N times r'$ feature matrix offline @caron2018deepcluster—the chief scalability bottleneck resolved by online SwAV (Alg. 190 @caron2020swav). Geometrically, DeepCluster induces spherical Voronoi cells in PCA-whitened space, making it fast on low-dimensional manifolds yet vulnerable to distance concentration @beyer1999nearest and whitening flips (Theorem II.A) in high dimension.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(E (N C_"net" + I_"km" N K r'))$ ($+ O(E N r^2)$ full PCA); $O(N (C_"net" + K r'))$ inference.],
    [Space Complexity], [$O(|theta| + N r' + K r')$; requires offline full-dataset feature cache each epoch.],
    [Inductive Bias], [Spherical Voronoi cells in PCA-whitened space reinforced by discriminative margins.],
    [Collapse Resistance], [Relies on uniform resampling (Lemma II.B), empty-cluster splitting, and head reset.],
    [Lyapunov Stability], [No joint monotonic Lyapunov function; whitening permits flips (Theorem II.A).],
    [Overclustering Gain], [Peaks at $K approx 10 K_"true"$; exact $K = K_"true"$ degrades curved-manifold fit.],
    [Hyperparameters], [Clusters $K$, PCA dim $r'$, Lloyd steps $I_"km"$, augmentation policy $cal(T)$.],
    [Streaming / Online], [Offline epoch-synchronized design; not natively streaming (prefer SwAV, Alg. 190).],
  ),
  caption: [Operational and geometric profile of DeepCluster (Alg. 191).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Pretext-Free Self-Supervised Vision Pretraining:* Random CNNs possess strong spatial-filter inductive bias at initialization @caron2018deepcluster; DeepCluster bootstraps this signal into transferable representations on ImageNet and YFCC100M without handcrafted pretext losses.
+ *Audio-Visual and Acoustic Event Discovery:* Pairing spectrogram CNNs with overclustered K-Means ($K >> K_"true"$) separates acoustic event modes across pitch and timbral shifts without complex audio augmentations.
+ *Histopathology Representation Bootstrapping:* Where aggressive color augmentations risk erasing diagnostic staining cues, DeepCluster leverages texture sensitivity with minimal geometric transformations.
+ *Non-Convex Tabular Manifold Warping:* On moderate-dimensional tabular benchmarks, alternating MLP discriminative steps with spherical K-Means rapidly unfolds curved shapes into compact Voronoi clusters without $k$-NN indexing.
+ *Unsupervised Domain Adaptation:* Periodic target-domain K-Means re-clustering adapts pretrained backbones to covariate shift while uniform cluster sampling prevents adaptation collapse onto majority target classes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes bound deployment: (1) *high-dimensional concentration* @beyer1999nearest and *PCA whitening equalization* (Theorem II.A) amplify noise axes and destabilize K-Means; (2) *epoch-to-epoch label oscillation* under aggressive learning rates; (3) *low-level RGB shortcuts* without Sobel filtering @caron2018deepcluster; and (4) *streaming / billion-scale memory limits* preventing full-dataset feature caching.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Visual Pretraining (Medium $N$)], [Strong: pretext-free bootstrapping with $K approx 10 K_"true"$ and Sobel.], [Overcluster $K$; for web scale use SwAV (Alg. 190).],
    [Non-Convex Low-D Manifolds], [Strong: discriminative MLP warps curves into compact Voronoi cells.], [Track epoch-to-epoch NMI for oscillation.],
    [High-D Noisy Spaces], [Poor: PCA whitening boosts noise axes and K-Means concentrates.], [Prefer IIC (Alg. 188) or SwAV (Alg. 190).],
    [Online / Billion-Scale Streams], [Poor: offline feature cache and global K-Means stall scaling.], [Use online SwAV (Alg. 190) or Mini-Batch KM.],
    [Skewed Class Clustering], [Moderate: uniform resampling blocks collapse but biases priors.], [Tune sampling weights or use SCAN (Alg. 189).],
  ),
  caption: [Domain suitability and failure-mode matrix for DeepCluster (Alg. 191).],
)

#sec("V", "Structural Relations and Empirical Behavior")

DeepCluster sits at the center of Category VIII between reconstruction-coupled deep K-Means (DCN @yang2017dcn, DEC @xie2016dec, IDEC @guo2017idec) and modern online multi-view clustering (SwAV @caron2020swav, SCAN @vangansbeke2020scan).

#proposition(title: "Unifying Hierarchy Across DeepCluster, DCN, and SwAV", origin: "Original Derivation")[
  Let $f_theta : bb(R)^d -> bb(S)^(r-1)$ and unit centroids $bold(C) = [bold(c)_1, ..., bold(c)_K]$.
  + Replacing DeepCluster's decoupled cross-entropy $cal(J)_"CE"^"bal"$ with joint reconstruction plus distortion $cal(L)_"rec" (theta) + (beta / 2) cal(J)_"KM" (bold(C), bold(Y); theta)$ recovers *DCN* @yang2017dcn, restoring monotonic Lyapunov descent at the cost of reconstruction capacity.
  + Replacing offline K-Means and reset head $bold(W)$ with persistent prototypes $bold(W) = bold(C)^top / tau$, online minibatch Sinkhorn equipartition $bold(Q) in cal(U)(B, K)$, and swapped view targets recovers *SwAV* @caron2020swav (Alg. 190), simultaneously removing the $O(N r')$ cache, whitening degeneracy (Theorem II.A), and permutation head reset (Theorem II.D).
]

#proof[
  In DCN @yang2017dcn, both alternating steps minimize the identical distortion $cal(J)_"KM"$ coupled to an autoencoder loss $cal(L)_"rec"$, guaranteeing monotonic non-increase across iterations, whereas DeepCluster drops the decoder and switches to $cal(J)_"CE"^"bal"$. Tying classifier rows persistently to $bold(w)_k = bold(c)_k / tau$ eliminates inter-epoch $sigma$-permutations (Theorem II.D), while online Sinkhorn transport on $cal(U)(B, K)$ replaces both global K-Means and uniform resampling (Lemma II.B), recovering SwAV @caron2020swav.
]

*Empirical Behavior on the ClusterOverview Benchmark.* On the 4-dataset ClusterOverview benchmark @rousseeuw1987silhouettes @hubert1985comparing, DeepCluster attains *Iris* $"Sil" = 0.360$, $"ARI" = 0.387$; *Customers* $"Sil" = 0.259$; *Multi-Shape* $"Sil" = 0.716$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.109$ ($"ARI" =$ N/A), at mean runtime *5.374 ms*. DeepCluster achieves the *highest Multi-Shape Silhouette ($0.716$)* across Algs. 188--191 while running faster than SwAV ($6.479$ ms) and SCAN ($7.005$ ms), confirming that iterative cross-entropy warping excels on low-dimensional non-convex geometries. Conversely, its sharp drop on *High-D* ($"Sil" = 0.109$, under half of IIC's $0.298$ and SwAV's $0.292$) directly reflects Theorem II.A and @beyer1999nearest: PCA whitening equalizes noise dimensions and high-dimensional K-Means distance concentration destabilizes pseudo-labels, while its moderate *Iris* ARI ($0.387$) trails SCAN ($0.632$) and IIC ($0.531$).
