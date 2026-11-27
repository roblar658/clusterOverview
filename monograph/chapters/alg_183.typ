#import "../template/ieee.typ": *

#alg_entry(
  id: 183,
  name: "DCN (Deep Clustering Network)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2017",
  paradigm: "Joint Autoencoder and Hard K-Means Optimization",
  complexity: $O(E N (W_"AE" + K p))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The Deep Clustering Network (DCN) @yang2017dcn formulates deep clustering as a direct joint optimization of nonlinear dimensionality reduction and hard Euclidean K-Means clustering @lloyd1982least. Rather than relying on heuristic soft self-training distributions (as in DEC @xie2016dec and IDEC @guo2017idec), DCN couples an autoencoder's reconstruction error directly with the discrete K-Means objective in latent space, yielding a "K-Means-friendly" representation whose latent clusters are explicitly shaped into compact, linearly separable Voronoi cells.

#definition(title: "K-Means-Friendly Latent Representation Objective")[
  Let $X = [x_1, ..., x_N] in bb(R)^(d times N)$ be the data matrix, let $f_theta : bb(R)^d -> bb(R)^p$ and $g_psi : bb(R)^p -> bb(R)^d$ be encoder and decoder networks, let $M = [mu_1, ..., mu_K] in bb(R)^(p times K)$ be latent centroids, and let $S = [s_1, ..., s_N] in {0, 1}^(K times N)$ be hard one-hot assignment vectors satisfying $bold(1)^top s_i = 1$. For regularization weight $lambda > 0$, the *DCN optimization problem* is
  $ min_(theta, psi, M, S) L_"DCN" (theta, psi, M, S) = 1 / N sum_(i=1)^N ( norm(x_i - g_psi (f_theta (x_i)))_2^2 \
    + lambda / 2 norm(f_theta (x_i) - M s_i)_2^2 ). $
]

Because $S$ is discrete and $(theta, psi)$ are non-convex neural parameters, Yang et al. @yang2017dcn solve $L_"DCN"$ via alternating block coordinate optimization after pretraining $(f_theta, g_psi)$ on reconstruction loss and initializing $(M, S)$ with K-Means++ @arthur2007kmeanspp:
(1) *Network step:* holding $(M, S)$ fixed, update $(theta, psi)$ by mini-batch SGD on the differentiable sum of reconstruction error and squared attraction toward each sample's assigned centroid $M s_i$;
(2) *Assignment step:* holding $(theta, psi, M)$ fixed, assign each embedding $z_i = f_theta (x_i)$ to its exact nearest centroid $s_(i j) = bb(I)[j = arg min_k norm(z_i - mu_k)_2^2]$;
(3) *Centroid step:* holding $(theta, psi, S)$ fixed, update each centroid either via exact batch mean $mu_k = abs(C_k)^(-1) sum_(i in C_k) z_i$ or via the streaming Robbins--Monro count update @bottou1995convergence @yang2017dcn: increment $c_k <- c_k + 1$ and set
$ mu_k <- mu_k - (1 / c_k) (mu_k - f_theta (x_i)) quad (i in C_k). $

#algbox(
  title: "Deep Clustering Network (DCN)",
  inputs: [Dataset $X = {x_1, ..., x_N} subset bb(R)^d$, clusters $K$, latent dim $p$, weight $lambda > 0$, epochs $E$, learning rate $eta$.],
  outputs: [Autoencoder $(f_theta, g_psi)$, centroids $M$, hard cluster assignments $S$.],
)[
  + Pretrain $(f_theta, g_psi)$ on $X$ minimizing $N^(-1) sum_i norm(x_i - g_psi (f_theta (x_i)))_2^2$.
  + Compute latent codes $z_i = f_theta (x_i)$ and initialize $(M, S)$ via K-Means++ on ${z_i}$; set counts $c_k <- 0$.
  + For epoch $t = 1, ..., E$:
    - *Network update:* for mini-batches $cal(B)$, step $(theta, psi)$ by SGD on $abs(cal(B))^(-1) sum_(i in cal(B)) (norm(x_i - g_psi (f_theta (x_i)))_2^2 + (lambda / 2) norm(f_theta (x_i) - M s_i)_2^2)$.
    - *Assignment update:* set $j^star (i) = arg min_k norm(f_theta (x_i) - mu_k)_2^2$ and $s_i <- e_(j^star(i))$ for all $i$.
    - *Centroid update:* compute exact cluster means $mu_k <- abs(C_k)^(-1) sum_(i in C_k) f_theta (x_i)$ (or online $c_k <- c_k + 1$, $mu_k <- mu_k - c_k^(-1) (mu_k - f_theta (x_i))$).
  + Return $(f_theta, g_psi, M, S)$ with labels $hat(y)_i = arg max_k s_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Unlike DEC and IDEC, DCN optimizes a single unified objective across all variable blocks. We prove exact Lyapunov monotonicity of batch alternating minimization, establish NP-hardness and an exact spectral PCA--LDA trace reduction for linear DCN, prove an ambient quantization upper bound alongside a scale-shrinkage counterexample, and derive DCN's finite-sample breakdown point.

#subsec("A", "Monotone Block-Coordinate Descent and Vanishing Residuals")

#theorem(title: "Exact Lyapunov Monotonicity of Batch Alternating DCN", origin: "Original Proof")[
  Assume $(theta, psi) |-> L_"DCN" (theta, psi, M, S)$ has $L_"lip"$-Lipschitz gradient for fixed $(M, S)$, and let batch DCN alternate gradient steps on $(theta, psi)$ with step size $eta in (0, 1 / L_"lip"]$, exact nearest-centroid updates for $S$, and exact cluster-mean updates for $M$. Then:
  (i) $L_"DCN"^((t+1)) <= L_"DCN"^((t))$ for every iteration $t >= 0$, and $L_"DCN"^((t))$ converges to a finite limit $L^star >= 0$;
  (ii) the network gradients and centroid displacements are square-summable:
  $ sum_(t=0)^infinity ( eta / 2 norm(nabla_(theta, psi) L_"DCN"^((t)))_2^2 + lambda / (2 N) sum_(k=1)^K abs(C_k^((t+1))) norm(mu_k^((t+1)) - mu_k^((t)))_2^2 ) \
    <= L_"DCN"^((0)) - L^star < infinity, $
  so $nabla_(theta, psi) L_"DCN"^((t)) -> 0$ and $norm(mu_k^((t+1)) - mu_k^((t)))_2 -> 0$ for every non-empty cluster $k$.
]

#proof[
  At iteration $t$, decompose the transition $(theta^((t)), psi^((t)), M^((t)), S^((t))) -> (theta^((t+1)), psi^((t+1)), M^((t+1)), S^((t+1)))$ into its three block updates.

  *Step 1 (Network gradient step).* By the standard descent lemma for $L_"lip"$-smooth objectives with $eta <= 1 / L_"lip"$,
  $ L_"DCN" (theta^((t+1)), psi^((t+1)), M^((t)), S^((t))) <= L_"DCN"^((t)) - (eta / 2) norm(nabla_(theta, psi) L_"DCN"^((t)))_2^2. $

  *Step 2 (Exact assignment step).* Let $z_i^((t+1)) = f_(theta^((t+1))) (x_i)$. Since $s_i^((t+1))$ minimizes $norm(z_i^((t+1)) - M^((t)) s_i)_2^2$ independently for every $i$, while $L_r$ does not depend on $S$,
  $ L_"DCN" (theta^((t+1)), psi^((t+1)), M^((t)), S^((t+1))) <= L_"DCN" (theta^((t+1)), psi^((t+1)), M^((t)), S^((t))). $

  *Step 3 (Exact centroid step).* Setting $mu_k^((t+1)) = abs(C_k^((t+1)))^(-1) sum_(i in C_k^((t+1))) z_i^((t+1))$ satisfies the parallel-axis identity for every non-empty cluster $C_k^((t+1))$:
  $ sum_(i in C_k^((t+1))) norm(z_i^((t+1)) - mu_k^((t)))_2^2 - sum_(i in C_k^((t+1))) norm(z_i^((t+1)) - mu_k^((t+1)))_2^2 \
    = abs(C_k^((t+1))) norm(mu_k^((t+1)) - mu_k^((t)))_2^2. $
  Multiplying by $lambda / (2 N)$ and combining all three steps yields
  $ L_"DCN"^((t)) - L_"DCN"^((t+1)) >= eta / 2 norm(nabla_(theta, psi) L_"DCN"^((t)))_2^2 + lambda / (2 N) sum_(k=1)^K abs(C_k^((t+1))) norm(mu_k^((t+1)) - mu_k^((t)))_2^2. $
  Since $L_"DCN" >= 0$, ${L_"DCN"^((t))}$ is monotone non-increasing and bounded below, hence converges to $L^star >= 0$. Telescoping the inequality from $t = 0$ to $T$ and letting $T -> infinity$ yields the finite sum bound and vanishing residuals.
]

#subsec("B", "NP-Hardness and Exact Linear-DCN Spectral Reduction")

#proposition(title: "NP-Hardness via Exact Reduction from Euclidean K-Means", origin: "Classical — Cited")[
  Globally minimizing $L_"DCN" (theta, psi, M, S)$ is NP-hard even when $f_theta$ and $g_psi$ are restricted to linear orthonormal maps ($p = d$) @aloise2009np @yang2017dcn.
]

#proof[
  Restrict the encoder and decoder to $f_W (x) = W^top x$ and $g_W (z) = W z$ with $W in O(d)$ ($p = d$). Then $g_W (f_W (x_i)) = W W^top x_i = x_i$, so the reconstruction term vanishes identically ($L_r = 0$). By orthogonal invariance of the Euclidean norm, setting $tilde(mu)_k = W mu_k$ transforms the remaining objective into $(lambda / (2 N)) sum_(i=1)^N norm(x_i - tilde(M) s_i)_2^2$, which is the exact Euclidean K-Means objective—proven NP-hard for general $d$ even at $K = 2$ @aloise2009np and for general $K$ even at $d = 2$ @mahajan2012planar.
]

#theorem(title: "Exact Trace Reduction of Orthonormal Linear DCN", origin: "Original Derivation")[
  Consider tied orthonormal linear DCN on centered data ($X bold(1) = 0$) with $f_W (x) = W^top x$, $g_W (z) = W z$, and $W in bb(R)^(d times p)$ satisfying $W^top W = I_p$. Let $S_T = N^(-1) X X^top$ be total scatter, let $H_S = S^top (S S^top)^(-1) S$ be the orthogonal cluster-mean projector for partition $S$ (with all $K$ clusters non-empty), and define between-cluster scatter $S_B (S) = N^(-1) X H_S X^top$ and within-cluster scatter $S_W (S) = S_T - S_B (S)$. Minimizing out $M$ analytically gives
  $ min_M L_"DCN" (W, M, S) = "Tr"(S_T) - "Tr"(W^top (S_T - lambda / 2 S_W (S)) W) \
    = "Tr"(S_T) - "Tr"(W^top (S_B (S) + (1 - lambda / 2) S_W (S)) W). $
  Consequently, for any fixed $S$, the globally optimal subspace $W^star (S)$ is spanned by the top-$p$ eigenvectors of $S_T - (lambda / 2) S_W (S)$, interpolating smoothly between exact PCA ($lambda = 0$), pure between-cluster discriminant projection ($lambda = 2$), and additive Fisher-criterion maximization ($lambda > 2$).
]

#proof[
  With $Z = W^top X in bb(R)^(p times N)$ and $W^top W = I_p$, the reconstruction term simplifies via $W W^top$ being an orthogonal projector:
  $ 1 / N norm(X - W W^top X)_F^2 = "Tr"(S_T) - "Tr"(W^top S_T W). $
  For fixed $W$ and $S$, the unique optimal centroid matrix is $M^star = Z S^top (S S^top)^(-1)$, so $M^star S = Z H_S$. Because $I_N - H_S$ is symmetric and idempotent,
  $ 1 / N norm(Z - M^star S)_F^2 &= 1 / N "Tr"(Z (I_N - H_S)^2 Z^top) \
    &= 1 / N "Tr"(W^top X (I_N - H_S) X^top W) = "Tr"(W^top S_W (S) W). $
  Adding the reconstruction loss and $(lambda / 2)$ times the clustering loss yields
  $ min_M L_"DCN" (W, M, S) = "Tr"(S_T) - "Tr"(W^top (S_T - lambda / 2 S_W (S)) W). $
  Substituting $S_T = S_B (S) + S_W (S)$ gives $S_T - (lambda / 2) S_W (S) = S_B (S) + (1 - lambda / 2) S_W (S)$. By Ky Fan's variational principle, maximizing $"Tr"(W^top A W)$ over $W^top W = I_p$ is achieved by the top-$p$ eigenvectors of symmetric $A = S_T - (lambda / 2) S_W (S)$. At $lambda = 0$, $A = S_T$ recovers exact PCA; at $lambda = 2$, $A = S_B (S)$ retains exclusively between-cluster variation; for $lambda > 2$, $A = S_B (S) - (lambda / 2 - 1) S_W (S)$ directly penalizes within-cluster spread in additive Fisher form.
]

#subsec("C", "Ambient Quantization Guarantee vs. Scale-Shrinkage Degeneracy")

#theorem(title: "Ambient Quantization Upper Bound and Unconstrained Scale Shrinkage", origin: "Original Counterexample")[
  (i) *Ambient quantization guarantee:* If the decoder $g_psi$ is $L_g$-Lipschitz, then the decoded cluster prototypes ${g_psi (mu_k)}_(k=1)^K$ bound the *ambient-space* K-Means distortion directly in terms of $L_"DCN"$:
  $ 1 / N sum_(i=1)^N norm(x_i - g_psi (M s_i))_2^2 <= max(2, (4 L_g^2) / lambda) L_"DCN" (theta, psi, M, S). $
  (ii) *Zero-scaling degeneracy:* Conversely, in an unconstrained autoencoder with affine bottleneck layer $z = W_"enc" h + b_"enc"$ and decoder input $W_"dec" z + b_"dec"$, scaling $(W_"enc", b_"enc", M) |-> (c W_"enc", c b_"enc", c M)$ and $W_"dec" |-> c^(-1) W_"dec"$ preserves $L_r$ identically while shrinking the K-Means term by $c^2 -> 0$ as $c -> 0^+$, so $inf_(theta, psi, M, S) L_"DCN" = L_r^star$ for *every* partition $S$ unless weight decay or decoder Lipschitz regularization is enforced.
]

#proof[
  (i) For each sample $i$, apply the inequality $norm(a + b)_2^2 <= 2 norm(a)_2^2 + 2 norm(b)_2^2$ together with $L_g$-Lipschitz continuity of $g_psi$:
  $ norm(x_i - g_psi (M s_i))_2^2 &<= 2 norm(x_i - g_psi (z_i))_2^2 + 2 norm(g_psi (z_i) - g_psi (M s_i))_2^2 \
    &<= 2 norm(x_i - g_psi (z_i))_2^2 + 2 L_g^2 norm(z_i - M s_i)_2^2. $
  Averaging over $i = 1, ..., N$ and writing $2 L_g^2 = (4 L_g^2 / lambda) (lambda / 2)$ yields
  $ 1 / N sum_(i=1)^N norm(x_i - g_psi (M s_i))_2^2 &<= 2 L_r + (4 L_g^2) / lambda (lambda / (2 N) sum_(i=1)^N norm(z_i - M s_i)_2^2) \
    &<= max(2, (4 L_g^2) / lambda) L_"DCN". $
  (ii) Under the bottleneck rescaling by $c > 0$, $W'_"dec" z'_i = (c^(-1) W_"dec")(c z_i) = W_"dec" z_i$, leaving $g_(psi') (f_(theta') (x_i)) = g_psi (f_theta (x_i))$ and $L_r$ invariant. Meanwhile, the K-Means penalty is homogeneous of degree $2$ in $(z, M)$:
  $ (lambda / (2 N)) sum_(i=1)^N norm(c z_i - c M s_i)_2^2 = c^2 dot.op (lambda / (2 N)) sum_(i=1)^N norm(z_i - M s_i)_2^2 -> 0 $
  as $c -> 0^+$, regardless of $S$. Enforcing bounded $L_g$ (e.g., via $L_2$ weight decay or spectral normalization) prohibits $norm(W'_"dec")_2 = c^(-1) norm(W_"dec")_2 -> infinity$ and activates part (i).
]

#subsec("D", "Outlier Sensitivity and Breakdown Point")

#proposition(title: "Linear Centroid Influence and Zero Breakdown Point", origin: "Original Proof")[
  Assume the encoder is coercive (or has bounded decoder Lipschitz constant $L_g < infinity$ with $L_r <= C_0$). Then under exact centroid updates, DCN has finite-sample adversarial breakdown point $1 / N$ (asymptotic $0\%$), in contrast to the uniformly bounded redescending influence of DEC and IDEC.
]

#proof[
  Fix $N - 1$ inliers inside $B(0, 1)$ and let a single corrupted sample $x_"out"$ satisfy $norm(x_"out")_2 = R -> infinity$. Since $L_r <= C_0$ implies $norm(x_"out" - g_psi (z_"out"))_2 <= sqrt(N C_0)$, the $L_g$-Lipschitz condition forces $norm(z_"out")_2 >= (R - sqrt(N C_0) - norm(g_psi (0))_2) / L_g -> infinity$. In the centroid update for the cluster $C_k$ containing $x_"out"$,
  $ mu_k = 1 / abs(C_k) (z_"out" + sum_(i in C_k without {"out"}) z_i), $
  so $norm(mu_k)_2 -> infinity$ linearly at rate $norm(z_"out")_2 / abs(C_k)$ (whether $x_"out"$ forms a singleton cluster or drags an entire inlier cluster's centroid to infinity). Thus a single sample arbitrarily corrupts $M$, proving breakdown point $1 / N$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

DCN shares the $O(E N (W_"AE" + K p))$ asymptotic training complexity of IDEC, but avoids computing and storing an $O(N K)$ global target matrix $P = cal(T)(Q)$. Both hard nearest-centroid assignments ($O(K p)$ per sample) and online running-count centroid updates ($O(p)$ per sample) execute strictly inside mini-batches, matching mini-batch K-Means @sculley2010webscale in scalability and reducing peak memory to $O(B sum_ell d_ell^"AE" + W_"AE" + K p)$.

Geometrically, DCN shapes the latent space so that clusters occupy tight Euclidean Voronoi balls while retaining a decodable embedding manifold (Theorem 3). Theorem 2 shows why this works even linearly: as $lambda$ increases, the encoder shifts spectral energy away from within-cluster scatter $S_W$ and concentrates it on between-cluster separation $S_B$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O((E_"pre" + E_"ft") N W_"AE" + E_"ft" N K p)$],
    [Space Complexity], [$O(B sum_ell d_ell^"AE" + W_"AE" + K p)$ (no $N K$ table)],
    [Cluster Geometry], [K-Means-friendly Voronoi polytopes in latent space],
    [Cluster Count $K$], [Fixed a priori; NP-hard globally @aloise2009np],
    [Hyperparameters], [Architecture, $K, p, lambda$, learning rate $eta$],
    [Convergence], [Exact monotone Lyapunov descent in batch mode],
    [Outlier Breakdown], [$1 / N$ (unbounded linear centroid pull)],
    [Streaming Fit], [High; native mini-batch running-count centroid step],
  ),
  caption: [Operational and geometric profile of Deep Clustering Network.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Large-Scale Streaming Document and Web-Log Clustering.* Because DCN's online centroid update requires no global $N times K$ target matrix, it scales directly to millions of TF-IDF documents (e.g., Reuters, RCV1 @yang2017dcn) with constant memory in $N$ while learning K-Means-friendly latent topics.

*2. Hyperspectral Remote Sensing and Spectral Unmixing.* High-dimensional reflectance spectra exhibit nonlinear atmospheric mixing alongside discrete land-cover classes. Theorem 2's additive Fisher-like weighting suppresses within-material illumination scatter while preserving material separability.

*3. Neural Vector Quantization and Learned Codebook Compression.* Theorem 3(i) proves that minimizing $L_"DCN"$ with a Lipschitz decoder directly upper-bounds the ambient quantization distortion of decoded codewords ${g_psi (mu_k)}$, providing a principled foundation for neural codebook learning.

*4. Genomics and Proteomics Expression Stratification.* When biological cohorts differ along nonlinear gene-regulatory combinations, DCN's joint reconstruction--discriminant projection boosts between-class variance over PCA+K-Means while retaining reconstructive gene profiles for interpretability.

*5. Healthcare and Customer Cohort Segmentation.* Clinical and commercial registries often require crisp, mutually exclusive cohort assignments with interpretable representative profiles; hard one-hot assignments $S$ and decoded centroids $g_psi (mu_k)$ satisfy both requirements natively.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes require mitigation: (i) *unconstrained scale shrinkage* (Theorem 3(ii)), where $c -> 0^+$ nullifies the K-Means penalty unless weight decay or batch/spectral normalization is applied; (ii) *outlier sensitivity* (Proposition 2), where heavy-tailed noise pulls centroids away from inliers; (iii) *empty clusters* during hard nearest-centroid reassignment when $lambda$ is large; and (iv) *discrete assignment chattering* under noisy mini-batch gradients if $eta$ is not decayed.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Large-Scale Text], [Mini-batch $O(K p)$ updates, no $N K$ matrix], [Apply $L_2$ weight decay (Thm. 3)],
    [Neural Codebooks], [Ambient quantization upper bound (Thm. 3)], [Handle empty codewords by reseeding],
    [Overlapping Gaussians], [Fisher-like $S_B$ amplification (Thm. 2)], [Tune $lambda$ via elbow/validation],
    [Heavy-Tailed Noise], [Quadratic penalty has $1/N$ breakdown], [Use IDEC, DEC, or Huber-DCN],
    [Continuous Manifolds], [Hard partitioning splits smooth arcs], [Use GCS or soft VaDE @jiang2017vade],
  ),
  caption: [Domain suitability and failure-mode matrix for DCN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Unified Comparison: DCN, K-Means, IDEC, DEC, and DeepCluster", origin: "Original Unified Proof")[
  Within deep clustering, DCN bridges classical Bregman quantization and modern self-supervised partitioning:
  (i) Fixing $(theta, psi)$ to the identity recovers Lloyd's K-Means @lloyd1982least (generalizable to any Bregman divergence @banerjee2005clustering);
  (ii) replacing DCN's hard quadratic penalty $(lambda / 2) norm(z_i - M s_i)_2^2$ with soft Student-$t$ KL divergence $gamma D_"KL" (P || Q)$ yields IDEC @guo2017idec, and further dropping $L_r$ yields DEC @xie2016dec—trading DCN's exact Lyapunov monotonicity (Theorem 1) for bounded redescending outlier influence;
  (iii) replacing DCN's autoencoder reconstruction anchor $L_r$ by a cross-entropy classifier head over periodic K-Means pseudo-labels yields DeepCluster @caron2018deepcluster.
]

#proof[
  (i) Immediate from setting $f_theta = g_psi = "Id"$, under which $L_r = 0$ and the alternating steps reduce verbatim to Lloyd's assignment and centroid-mean updates @lloyd1982least.
  (ii) Comparing objectives directly, $L_"DCN"$ and $L_"IDEC"$ share $L_r (theta, psi)$ identically; DCN replaces IDEC's non-Lyapunov target-sharpening functional $D_"KL" (cal(T)(Q) || Q)$ with a joint quadratic functional over $(theta, M, S)$ that decreases monotonically at every block update (Theorem 1) at the expense of linear outlier sensitivity (Proposition 2).
  (iii) DeepCluster retains DCN's alternating K-Means assignment step in latent space but removes the decoder $g_psi$, substituting a softmax classification head and uniform pseudo-label sampling to avoid trivial collapse @caron2018deepcluster.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark, DCN achieves $"Sil" = 0.441$ @rousseeuw1987silhouettes and $"ARI" = 0.648$ @hubert1985comparing on *Iris*, $"Sil" = 0.326$ on *Customers*, $"Sil" = 0.716$ ($"ARI"$ N/A) on *Multi-Shape*, and $"Sil" = 0.298$ ($"ARI"$ N/A) on *High-D*, with a mean runtime of *10.429 ms*. This empirical profile aligns remarkably well with our theorems. Most strikingly, DCN attains by far the highest ground-truth recovery on *Iris* ($"ARI" = 0.648$ vs. IDEC $0.584$, DEC $0.574$, and GCS $0.563$): Theorem 2's additive Fisher-like suppression of within-cluster scatter $S_W$ separates the overlapping Versicolor and Virginica manifolds much more cleanly than either raw quantization or frequency-balanced KL self-training. Simultaneously, on *Multi-Shape* ($"Sil" = 0.716$), DCN virtually matches DEC ($0.718$) and greatly outperforms IDEC ($0.563$), showing that direct quadratic attraction toward hard centroids reshapes non-convex clusters into compact Voronoi balls even while retaining the decoder. On *Customers* ($"Sil" = 0.326$), DCN trails DEC ($0.344$) slightly—consistent with Proposition 2's higher sensitivity to heavy-tailed spending outliers—while on *High-D* ($"Sil" = 0.298$) it matches IDEC ($0.300$) and outperforms DEC ($0.295$) against distance concentration @beyer1999nearest. Finally, DCN is the fastest of the three deep methods (*10.429 ms* vs. IDEC $10.549$ ms and DEC $10.574$ ms), confirming the lower overhead of $O(K p)$ hard updates over full-batch Student-$t$ target passes.
