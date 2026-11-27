#import "../template/ieee.typ": *

#alg_entry(
  id: 182,
  name: "IDEC (Improved DEC)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2017",
  paradigm: "Structure-Preserving Deep Embedded Clustering",
  complexity: $O(E N (W_"AE" + K p))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Improved Deep Embedded Clustering (IDEC) @guo2017idec directly addresses the fundamental theoretical vulnerability of DEC @xie2016dec: once DEC discards its decoder after pretraining, minimizing the clustering KL divergence alone can distort or collapse the latent manifold (destroying local geometric structure that is not aligned with initial pseudo-labels). IDEC resolves this pathology by keeping the full undercomplete autoencoder $(f_theta, g_psi)$ active throughout clustering and jointly optimizing reconstruction fidelity alongside KL self-training.

#definition(title: "Structure-Preserving Deep Embedded Clustering Objective")[
  Let $f_theta : bb(R)^d -> bb(R)^p$ be an encoder with latent codes $z_i = f_theta (x_i)$, let $g_psi : bb(R)^p -> bb(R)^d$ be a decoder, and let $M = {mu_j}_(j=1)^K subset bb(R)^p$ be $K$ cluster centroids. Define the empirical reconstruction loss and Student-$t$ KL clustering loss ($alpha = 1$) by
  $ L_r (theta, psi) &= 1 / N sum_(i=1)^N norm(x_i - g_psi (f_theta (x_i)))_2^2, \
    L_c (theta, M) &= 1 / N sum_(i=1)^N sum_(j=1)^K p_(i j) log (p_(i j) / q_(i j)), $
  where soft assignments $q_(i j)$ and sharpened balanced targets $p_(i j)$ are
  $ q_(i j) = (1 + norm(z_i - mu_j)_2^2)^(-1) / (sum_(j'=1)^K (1 + norm(z_i - mu_(j'))_2^2)^(-1)), \
    p_(i j) = (q_(i j)^2 / f_j) / (sum_(j'=1)^K q_(i j')^2 / f_(j')), quad f_j = sum_(i=1)^N q_(i j). $
  For a trade-off hyperparameter $gamma > 0$, the *IDEC objective* is
  $ L_"IDEC" (theta, psi, M) = L_r (theta, psi) + gamma L_c (theta, M). $
]

Note that Guo et al. @guo2017idec originally omit the $1 / N$ normalization on $L_c$, which is algebraically equivalent to rescaling $gamma <- N gamma$; using per-sample averages keeps both terms $O(1)$ across dataset sizes.

After pretraining $(f_theta, g_psi)$ on $L_r$ and initializing $M$ via K-Means @lloyd1982least on ${f_theta (x_i)}$, IDEC alternates between computing the target distribution $P = cal(T)(Q)$ every $T_"upd"$ steps and updating all three parameter blocks $(theta, psi, M)$ simultaneously via mini-batch stochastic gradient descent with $P$ held fixed. Because the encoder receives the coupled gradient $nabla_theta L_r + gamma nabla_theta L_c$, the decoder acts as an active geometric anchor: any latent movement that improves cluster sharpness must remain invertible back to input space $X$ without incurring excessive reconstruction distortion.

#algbox(
  title: "Improved Deep Embedded Clustering (IDEC)",
  inputs: [Dataset $X = {x_1, ..., x_N} subset bb(R)^d$, clusters $K$, latent dim $p$, weight $gamma > 0$, update interval $T_"upd"$, tolerance $delta$, learning rate $eta$.],
  outputs: [Encoder--decoder $(f_theta, g_psi)$, centroids $M$, cluster assignments $hat(y)_1, ..., hat(y)_N$.],
)[
  + Pretrain full autoencoder $(f_theta, g_psi)$ to minimize $L_r (theta, psi) = N^(-1) sum_i norm(x_i - g_psi (f_theta (x_i)))_2^2$; *retain* $g_psi$.
  + Compute latent codes $z_i = f_theta (x_i)$ and initialize centroids $M = {mu_j}_(j=1)^K$ via K-Means on ${z_i}$.
  + Initialize soft assignments $Q$ and hard labels $hat(y)_i^((0)) = arg max_j q_(i j)$.
  + For $t = 1, 2, ...$:
    - If $t mod T_"upd" = 0$: recompute $Q$, update $P <- cal(T)(Q)$, set $hat(y)_i^((t)) = arg max_j q_(i j)$, and terminate if $N^(-1) sum_i bb(I)[hat(y)_i^((t)) != hat(y)_i^((t-1))] < delta$.
    - Draw mini-batch $cal(B) subset.eq {1, ..., N}$ and update $(theta, psi, M)$ by SGD on $L_r (cal(B)) + gamma L_c (cal(B); P)$ with $P$ frozen.
  + Output $(f_theta, g_psi, M)$ and labels $hat(y)_i = arg max_j q_(i j)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove a quantitative bi-Lipschitz non-collapse theorem showing how reconstruction prevents latent collapse, characterize the exact Pareto trade-off in $gamma$, uncover a fundamental bottleneck scale-bypass degeneracy in unconstrained autoencoders together with its weight-decay cure, and establish coupled-gradient dynamics.

#subsec("A", "Quantitative Non-Collapse and Local Structure Preservation")

#theorem(title: "Lipschitz Non-Collapse and Pairwise Separation Bound", origin: "Original Proof")[
  Suppose the decoder $g_psi : bb(R)^p -> bb(R)^d$ is $L_g$-Lipschitz: $norm(g_psi (z) - g_psi (z'))_2 <= L_g norm(z - z')_2$ for all $z, z' in bb(R)^p$. Let $e_i = norm(x_i - g_psi (z_i))_2$ with empirical reconstruction MSE $epsilon_r = L_r (theta, psi) = N^(-1) sum_(i=1)^N e_i^2$. Denote the ambient and latent standard deviations by
  $ sigma_X = (1 / N sum_(i=1)^N norm(x_i - bar(x))_2^2)^(1/2), quad
    sigma_Z = (1 / N sum_(i=1)^N norm(z_i - bar(z))_2^2)^(1/2). $
  Then every IDEC iterate satisfies:
  (i) *Pairwise distance lower bound:* for all $i, j in {1, ..., N}$,
  $ norm(z_i - z_j)_2 >= 1 / L_g (norm(x_i - x_j)_2 - e_i - e_j); $
  (ii) *Global non-collapse bound:*
  $ sigma_Z >= max(0, sigma_X - sqrt(epsilon_r)) / L_g; $
  (iii) *Anti-point-collapse bound:* if $z_i = mu_(c(i))$ collapses to $K$ exact centroids, then $epsilon_r$ is lower-bounded by the optimal ambient $K$-Means error of $X$:
  $ epsilon_r >= min_(abs(cal(C)) <= K) 1 / N sum_(i=1)^N min_(c in cal(C)) norm(x_i - c)_2^2. $
]

#proof[
  (i) By the triangle inequality and $L_g$-Lipschitz continuity of $g_psi$,
  $ norm(x_i - x_j)_2 &<= norm(x_i - g_psi (z_i))_2 + norm(g_psi (z_i) - g_psi (z_j))_2 \
    &quad + norm(g_psi (z_j) - x_j)_2 \
    &<= e_i + L_g norm(z_i - z_j)_2 + e_j. $
  Rearranging yields the pairwise separation lower bound.

  (ii) Define centered sample matrices $X_c, G_c, E_c in bb(R)^(d times N)$ with columns $x_i - bar(x)$, $g_psi (z_i) - bar(g)$, and $(x_i - g_psi (z_i)) - (bar(x) - bar(g))$, respectively, so $X_c = G_c + E_c$. Because mean-centering is an orthogonal projection across columns, $norm(E_c)_F <= norm(E)_F = sqrt(N epsilon_r)$. By the Frobenius triangle inequality,
  $ sqrt(N) sigma_X = norm(X_c)_F <= norm(G_c)_F + norm(E_c)_F <= norm(G_c)_F + sqrt(N epsilon_r). $
  Moreover, the pairwise variance identity and $L_g$-Lipschitzness give
  $ norm(G_c)_F^2 &= 1 / (2 N) sum_(i, j=1)^N norm(g_psi (z_i) - g_psi (z_j))_2^2 \
    &<= L_g^2 / (2 N) sum_(i, j=1)^N norm(z_i - z_j)_2^2 = N L_g^2 sigma_Z^2. $
  Substituting $norm(G_c)_F <= sqrt(N) L_g sigma_Z$ and dividing by $sqrt(N)$ yields $sigma_X <= L_g sigma_Z + sqrt(epsilon_r)$, proving (ii).

  (iii) If $z_i = mu_(c(i))$ takes at most $K$ distinct values, then ${g_psi (z_i)}_(i=1)^N subset.eq cal(C)_psi := {g_psi (mu_1), ..., g_psi (mu_K)}$ has cardinality at most $K$. Hence
  $ epsilon_r = 1 / N sum_(i=1)^N norm(x_i - g_psi (mu_(c(i))))_2^2 >= 1 / N sum_(i=1)^N min_(c in cal(C)_psi) norm(x_i - c)_2^2, $
  which is bounded below by the optimal ambient $K$-Means distortion. Thus whenever $epsilon_r$ remains strictly below the ambient $K$-Means error (typically achieved by any undercomplete autoencoder with $p >= 2$ on continuous manifolds), total collapse to $K$ points (as occurs in DEC) is mathematically impossible!
]

#subsec("B", "Pareto Monotonicity and Explicit Distortion Budget")

#lemma(title: "Concavity of Value Function and Distortion Budget in gamma", origin: "Original Derivation")[
  Fix a target matrix $P$ (or consider exact joint minimizers) and define the optimal value function $Phi(gamma) = inf_(theta, psi, M) {L_r (theta, psi) + gamma L_c (theta, M)}$ for $gamma >= 0$, attained at $(theta_gamma, psi_gamma, M_gamma)$. Let $L_r^star = inf_(theta, psi) L_r (theta, psi)$ and $L_c^0 = L_c (theta_0, M_0)$ at the pretrained initialization. Then:
  (i) $Phi(gamma)$ is concave and non-decreasing on $[0, infinity)$;
  (ii) $L_c (theta_gamma, M_gamma)$ is non-increasing in $gamma$, while $L_r (theta_gamma, psi_gamma)$ is non-decreasing in $gamma$;
  (iii) any solution improving upon initialization obeys the explicit reconstruction distortion budget
  $ L_r (theta_gamma, psi_gamma) <= L_r^star + gamma L_c^0 <= L_r^star + gamma log (N K / f_"min"^((0))), $
  where $f_"min"^((0)) = min_j f_j^((0)) > 0$ (reducing to $L_r^star + gamma log K$ when initial soft frequencies are balanced).
]

#proof[
  (i) For each fixed $(theta, psi, M)$, the map $gamma |-> L_r (theta, psi) + gamma L_c (theta, M)$ is affine with non-negative slope $L_c >= 0$. As the pointwise infimum of non-decreasing affine functions, $Phi(gamma)$ is concave and non-decreasing.

  (ii) Let $0 < gamma_1 < gamma_2$ with optimal solutions $u_1, u_2$. Optimality yields
  $ L_r (u_1) + gamma_1 L_c (u_1) &<= L_r (u_2) + gamma_1 L_c (u_2), \
    L_r (u_2) + gamma_2 L_c (u_2) &<= L_r (u_1) + gamma_2 L_c (u_1). $
  Adding the two inequalities cancels $L_r$ and gives $(gamma_2 - gamma_1)(L_c (u_2) - L_c (u_1)) <= 0$, so $L_c (u_2) <= L_c (u_1)$. Substituting back into the first inequality gives $L_r (u_1) - L_r (u_2) <= gamma_1 (L_c (u_2) - L_c (u_1)) <= 0$, so $L_r (u_1) <= L_r (u_2)$.

  (iii) Because $L_c >= 0$, $L_r (theta_gamma, psi_gamma) <= L_"IDEC" (theta_gamma, psi_gamma, M_gamma) <= L_r^star + gamma L_c^0$. To bound $L_c^0$, note that by Cauchy--Schwarz, $sum_(j'=1)^K q_(i j')^2 >= (sum_(j') q_(i j'))^2 / K = 1 / K$. Since every $f_(j') <= N$,
  $ sum_(j'=1)^K q_(i j')^2 / f_(j') >= 1 / N sum_(j'=1)^K q_(i j')^2 >= 1 / (N K). $
  Therefore, using $q_(i j) <= 1$ and $f_j >= f_"min"^((0))$,
  $ p_(i j) / q_(i j) = (q_(i j) / f_j) / (sum_(j'=1)^K q_(i j')^2 / f_(j')) <= (N K) / f_"min"^((0)). $
  Taking logarithms and averaging with respect to $p_(i j)$ gives $L_c^0 <= log (N K / f_"min"^((0)))$. When $f_j^((0)) = N / K$ for all $j$, $f_j$ cancels directly to give $p_(i j) / q_(i j) = q_(i j) / sum_(j') q_(i j')^2 <= K$, yielding $L_c^0 <= log K$.
]

#subsec("C", "Bottleneck Scale-Bypass Degeneracy and Coercive Regularization")

#theorem(title: "Bottleneck Scale-Bypass Counterexample and Weight-Decay Cure", origin: "Original Counterexample")[
  Consider any deep autoencoder $(f_theta, g_psi)$ whose final encoder layer is affine ($z = W_"enc" h + b_"enc"$) and whose first decoder layer is affine ($tilde(h) = W_"dec" z + b_"dec"$), with arbitrary nonlinear hidden layers and biases elsewhere.
  (i) *Unconstrained scale bypass:* For every $c > 0$, rescaling $(W_"enc", b_"enc", M) |-> (c W_"enc", c b_"enc", c M)$ and $W_"dec" |-> c^(-1) W_"dec"$ leaves $g_psi (f_theta (x))$ and hence $L_r (theta, psi)$ *strictly invariant*, while driving $L_c (theta, M) -> 0$ as $c -> 0^+$. Thus without norm constraints on $g_psi$, $inf_(theta, psi, M) L_"IDEC" = L_r^star$ is achieved trivially with uniform assignments $q_(i j) -> 1 / K$!
  (ii) *Coercive cure:* Adding $L_2$ weight decay $beta (norm(theta)_2^2 + norm(psi)_2^2)$ ($beta > 0$) bounds $L_g <= C(beta)$ uniformly along any sublevel set, restoring Theorem 1's non-collapse bound and ruling out the $c -> 0^+$ bypass.
]

#proof[
  (i) Under the transformation $W'_"enc" = c W_"enc"$, $b'_"enc" = c b_"enc"$, $M' = c M$, and $W'_"dec" = c^(-1) W_"dec"$, the latent code becomes $z'_i = c z_i$, and the first decoder pre-activation is
  $ W'_"dec" z'_i + b_"dec" = (c^(-1) W_"dec")(c z_i) + b_"dec" = W_"dec" z_i + b_"dec". $
  Every subsequent decoder layer receives identical inputs, so $g_(psi') (f_(theta') (x_i)) = g_psi (f_theta (x_i))$ for all $x_i$ and $L_r (theta', psi') = L_r (theta, psi)$. Meanwhile, $norm(z'_i - mu'_j)_2^2 = c^2 norm(z_i - mu_j)_2^2 -> 0$ as $c -> 0^+$, whence $q_(i j) (c) -> 1 / K$, $p_(i j) (c) -> 1 / K$, and $L_c (c) -> 0$. Consequently, an unregularized autoencoder can eliminate the clustering loss entirely without modifying its reconstruction map!

  (ii) With weight decay $beta (norm(theta)_2^2 + norm(psi)_2^2)$, the rescaled penalty satisfies
  $ beta (c^2 norm(W_"enc")_F^2 + c^(-2) norm(W_"dec")_F^2) >= 2 beta norm(W_"enc")_F norm(W_"dec")_F $
  and diverges to $+infinity$ as either $c -> 0^+$ or $c -> infinity$ whenever $W_"dec" != 0$ and $W_"enc" != 0$. On any sublevel set $L_"IDEC"^"reg" <= C_0$, every decoder weight matrix obeys $norm(W_ell)_2 <= norm(W_ell)_F <= sqrt(C_0 / beta)$, bounding $L_g$ above by $(C_0 / beta)^(D_"dec" / 2)$ for 1-Lipschitz activations with decoder depth $D_"dec"$. Applying Theorem 1(ii) then guarantees $sigma_Z >= max(0, sigma_X - sqrt(C_0)) / L_g > 0$ whenever $C_0 < sigma_X^2$.
]

#subsec("D", "Coupled Gradient Structure and Alternating Convergence")

#proposition(title: "Joint Encoder--Decoder Gradient Alignment and Step-Size Stability", origin: "Classical — Cited")[
  For fixed target $P$, the IDEC parameter gradients satisfy @guo2017idec @xie2016dec
  $ nabla_psi L_"IDEC" &= nabla_psi L_r (theta, psi), \
    nabla_(mu_j) L_"IDEC" &= - (2 gamma) / N sum_(i=1)^N (1 + norm(z_i - mu_j)_2^2)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j), \
    nabla_theta L_"IDEC" &= 1 / N sum_(i=1)^N J_(f_theta) (x_i)^top ( -2 J_(g_psi) (z_i)^top (x_i - hat(x)_i) + gamma g_i^"clus" ), $
  where $hat(x)_i = g_psi (z_i)$ and $g_i^"clus" = 2 sum_j (1 + norm(z_i - mu_j)_2^2)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j)$. Under bounded activations and diminishing SGD step sizes $sum_t eta_t = infinity, sum_t eta_t^2 < infinity$, fixed-$P$ updates converge almost surely to a stationary point of $L_r + gamma L_c$ @bottou1995convergence.
]

#proof[
  Because $L_c$ is independent of $psi$, $nabla_psi L_"IDEC" = nabla_psi L_r$. Applying the chain rule through $z_i = f_theta (x_i)$ to $L_r = N^(-1) sum_i norm(x_i - g_psi (z_i))_2^2$ yields $partial L_r / partial z_i = -(2 / N) J_(g_psi) (z_i)^top (x_i - hat(x)_i)$. Combining with the $alpha = 1$ Student-$t$ gradient of DEC @xie2016dec gives the stated formula for $nabla_theta L_"IDEC"$. Geometrically, $nabla_theta L_r$ projects reconstruction residuals into the encoder's tangent space, opposing any component of $gamma g_i^"clus"$ that pushes $z_i$ off the reconstructible data manifold. Since $L_r + gamma L_c >= 0$ has Lipschitz gradients on compact sets, classical Robbins--Monro stochastic approximation @bottou1995convergence yields almost-sure convergence of $norm(nabla L_"IDEC")_2 -> 0$ for fixed $P$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Unlike DEC, IDEC evaluates both encoder and decoder forward and backward passes in every fine-tuning epoch, raising Phase 2 cost from $O(E_"ft" N (W_"enc" + K p))$ to $O(E_"ft" N (W_"AE" + K p))$ where $W_"AE" = W_"enc" + W_"dec" approx 2 W_"enc"$ for symmetric architectures. Because $L_r$ stabilizes the latent coordinate chart, IDEC often requires smaller $gamma$ (Guo et al. recommend $gamma = 0.1$) and converges cleanly without aggressive early stopping. Memory overhead is $O(B sum_ell d_ell^"AE" + W_"AE" + N K)$.

Geometrically, IDEC induces a *constrained diffeomorphic deformation*: rather than tearing non-convex manifolds into disconnected point masses (as DEC can), IDEC preserves ambient topology whenever $L_r$ stays small (Theorem 1), trading maximal latent compactness for faithful metric and neighborhood preservation.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O((E_"pre" + E_"ft") N W_"AE" + E_"ft" N K p)$],
    [Space Complexity], [$O(B sum_ell d_ell^"AE" + W_"AE" + N K)$],
    [Cluster Geometry], [Manifold-preserving soft Student-$t$ latent clusters],
    [Cluster Count $K$], [Fixed a priori; initialized by K-Means],
    [Hyperparameters], [Architecture, $K, p, gamma$ (default $0.1$), $T_"upd", delta$],
    [Collapse Protection], [Provable $sigma_Z$ floor when $L_g$ bounded (Thm. 1)],
    [Scale Regularization], [Requires weight decay or tied weights (Thm. 2)],
    [Dual Utility], [Simultaneous cluster labels + reconstruction residual],
  ),
  caption: [Operational and geometric profile of Improved DEC (IDEC).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Medical Imaging and Histopathology Patch Clustering.* Diagnostic tissue slides contain continuous morphological gradients alongside discrete phenotypes. Retaining $L_r$ ensures via Theorem 1(i) that subtle cellular variations are never collapsed away during clustering, allowing pathologists to decode and inspect cluster centroids $g_psi (mu_j)$.

*2. Single-Cell Genomics with Continuous Differentiation.* Biological tissues mix discrete cell types with continuous developmental trajectories. IDEC's reconstruction anchor preserves local cell-to-cell distances along trajectories while gently tightening cluster boundaries around mature populations.

*3. Joint Industrial Clustering and Anomaly Scoring.* Because $g_psi$ stays synchronized with the evolving latent space, the reconstruction error $norm(x - g_psi (f_theta (x)))_2^2$ serves directly as an out-of-distribution anomaly score alongside the Student-$t$ cluster membership $q_(i j)$, supporting simultaneous regime clustering and fault detection.

*4. Visual Retrieval and Reusable Representation Learning.* Embeddings collapsed by DEC lose intra-cluster ranking fidelity, hurting nearest-neighbor retrieval. By Lemma 1, tuning $gamma$ controls the exact Pareto point between retrieval fidelity ($L_r$) and semantic clustering ($L_c$).

*5. Speech and Acoustic Representation Discovery.* Phonetic frames vary smoothly in pitch and speaker timbre; preserving invertibility prevents speaker nuisance factors from distorting phoneme clusters while retaining signal synthesis capability.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Key failure modes include: (i) *over-regularization* ($gamma$ too small), where $L_r$ pins the latent space near PCA/autoencoder geometry and blocks unwrapping of intertwined non-convex manifolds; (ii) *bottleneck scale bypass* (Theorem 2), if weight decay or norm clipping is omitted; (iii) *pixel-space MSE bias*, where squared Euclidean reconstruction in raw image space prioritizes low-frequency background illumination over semantic foreground objects; and (iv) *class-imbalance distortion* inherited from DEC's $1 / f_j$ frequency normalization.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Medical Imaging], [Preserves subtle morphology + decodable centroids], [Use perceptual loss if raw MSE blurs edges],
    [Single-Cell Omics], [Balances discrete types and trajectories], [Tune $gamma$; or use VaDE @jiang2017vade],
    [Sensor Diagnostics], [Cluster labels + OOD reconstruction score], [Apply weight decay (Thm. 2)],
    [Strongly Non-Convex], [Decoder resists tearing linked shapes], [Increase $gamma$ or use DCN / DEC],
    [Natural RGB Images], [Pixel MSE focuses on color/background], [Use contrastive methods @caron2018deepcluster],
  ),
  caption: [Domain suitability and failure-mode matrix for IDEC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Homotopy Embedding of Autoencoders, IDEC, DEC, and DCN", origin: "Original Unified Proof")[
  Let $cal(F) = {(theta, psi, M)}$ be a weight-decayed autoencoder--centroid family. Then IDEC defines a continuous homotopy indexed by $gamma in [0, infinity)$:
  (i) at $gamma = 0$, IDEC reduces to pure undercomplete representation learning (and to exact top-$p$ PCA for linear tied weights);
  (ii) rescaling $gamma^(-1) L_"IDEC" = gamma^(-1) L_r + L_c$ and letting $gamma -> infinity$ recovers unregularized DEC @xie2016dec;
  (iii) replacing $L_c = D_"KL" (P || Q)$ by hard latent quantization $(lambda / 2) N^(-1) sum_i norm(z_i - mu_(c(i)))_2^2$ transforms IDEC directly into DCN @yang2017dcn, while replacing $L_r$ by a Gaussian-mixture ELBO yields VaDE @jiang2017vade @kingma2014vae.
]

#proof[
  (i) Setting $gamma = 0$ removes $M$ from $L_"IDEC"$, leaving $min_(theta, psi) L_r$; for $f_theta (x) = W^top x$ and $g_psi (z) = W z$ ($W^top W = I_p$), the Eckart--Young--Mirsky theorem identifies the global minimizer with the top-$p$ principal subspace of $X X^top$.
  (ii) Dividing $L_"IDEC"$ by $gamma > 0$ preserves all minimizers and yields $gamma^(-1) L_r + L_c$; as $gamma -> infinity$ on a bounded set, $gamma^(-1) L_r -> 0$ uniformly, recovering $L_"DEC"$.
  (iii) Both IDEC and DCN retain identical $L_r (theta, psi)$ and differ strictly in their clustering functional: soft Student-$t$ KL self-training with $1 / f_j$ balancing vs. hard Euclidean K-Means assignment.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark, IDEC achieves $"Sil" = 0.444$ @rousseeuw1987silhouettes and $"ARI" = 0.584$ @hubert1985comparing on *Iris*, $"Sil" = 0.297$ on *Customers*, $"Sil" = 0.563$ ($"ARI"$ N/A) on *Multi-Shape*, and $"Sil" = 0.300$ ($"ARI"$ N/A) on *High-D*, with a mean runtime of *10.549 ms*. These results expose the exact mathematical trade-off proved in Theorem 1 and Lemma 1. On *Iris* and *High-D*, IDEC surpasses DEC in both label recovery ($"ARI"$ $0.584$ vs. $0.574$) and high-dimensional separation ($"Sil"$ $0.300$ vs. $0.295$), confirming that retaining the decoder prevents spurious feature distortion under high-dimensional noise @beyer1999nearest and overlapping Gaussian species. Conversely, on *Multi-Shape* ($"Sil" = 0.563$ vs. DEC $0.718$ and DCN $0.716$) and *Customers* ($"Sil" = 0.297$ vs. DEC $0.344$), the reconstruction anchor constrains the encoder from aggressively collapsing non-convex arcs or skewed spending tails into tight spherical balls—precisely the pairwise separation floor of Theorem 1(i) in action. Its $10.549$ ms runtime matches DEC ($10.574$ ms) closely, demonstrating that keeping the decoder imposes negligible practical latency overhead when convergence stabilizes earlier @ghasedi2017depict.
