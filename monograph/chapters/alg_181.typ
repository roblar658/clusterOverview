#import "../template/ieee.typ": *

#alg_entry(
  id: 181,
  name: "DEC (Deep Embedded Clustering)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2016",
  paradigm: "Autoencoder Self-Training via KL Divergence",
  complexity: $O(E N (W_"enc" + K p))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Deep Embedded Clustering (DEC) @xie2016dec is a seminal deep clustering paradigm that jointly optimizes a nonlinear neural embedding $f_theta : bb(R)^d -> bb(R)^p$ ($p << d$) and $K$ latent cluster centroids $M = [mu_1, ..., mu_K] in bb(R)^(p times K)$ without requiring labeled supervision. Rather than clustering raw high-dimensional inputs directly with K-Means @lloyd1982least, DEC first learns a low-dimensional manifold via a stacked denoising autoencoder, discards the decoder, and then iteratively sharpens soft cluster assignments by minimizing a Kullback--Leibler (KL) divergence against a self-generated auxiliary target distribution.

#definition(title: "Student-t Soft Assignment and Auxiliary Target Operator")[
  Given latent embeddings $z_i = f_theta (x_i) in bb(R)^p$ ($i = 1, ..., N$) and centroids $M = {mu_j}_(j=1)^K subset bb(R)^p$, the *Student-$t$ soft assignment matrix* $Q in (0, 1)^(N times K)$ with degrees of freedom $alpha > 0$ (default $alpha = 1$) has entries
  $ q_(i j) = (1 + norm(z_i - mu_j)_2^2 / alpha)^(-(alpha + 1) / 2) / (sum_(j'=1)^K (1 + norm(z_i - mu_(j'))_2^2 / alpha)^(-(alpha + 1) / 2)). $
  Let $f_j = sum_(i=1)^N q_(i j)$ denote the soft cluster frequency of cluster $j$. The *auxiliary target operator* $cal(T) : Q |-> P$ constructs $P in (0, 1)^(N times K)$ by squaring and normalizing:
  $ p_(i j) = (q_(i j)^2 / f_j) / (sum_(j'=1)^K q_(i j')^2 / f_(j')). $
]

DEC minimizes the KL divergence between target distribution $P$ and model distribution $Q$:
$ L_"DEC" (theta, M) = D_"KL" (P || Q) = sum_(i=1)^N sum_(j=1)^K p_(i j) log (p_(i j) / q_(i j)). $
Because $P$ is treated as a fixed pseudo-label target during each inner optimization stage (stop-gradient on $P$), the active clustering loss over mini-batches equals the cross-entropy $-sum_(i, j) p_(i j) log q_(i j)$ up to the constant $sum_(i, j) p_(i j) log p_(i j)$.

Training proceeds in two distinct phases. In *Phase 1 (Parameter Initialization)*, an autoencoder $(f_theta, g_psi)$ is pretrained (optionally layer-wise as a denoising autoencoder) to minimize reconstruction error $N^(-1) sum_(i=1)^N norm(x_i - g_psi (f_theta (x_i)))_2^2$. The decoder $g_psi$ is then discarded, and $K$-Means++ @arthur2007kmeanspp initializes $M^((0))$ in the pretrained latent space ${f_(theta^((0))) (x_i)}$. In *Phase 2 (Self-Training)*, DEC alternates between recomputing $P = cal(T)(Q)$ every $T_"upd"$ iterations and performing stochastic gradient descent (SGD with momentum) on $(theta, M)$ with $P$ held fixed. Optimization halts when the fraction of hard-assignment flips $N^(-1) sum_(i=1)^N bb(I)[hat(y)_i^((t)) != hat(y)_i^((t-1))]$ falls below a tolerance $delta$ (e.g., $delta = 10^(-3)$), where $hat(y)_i = arg max_j q_(i j)$.

#algbox(
  title: "Deep Embedded Clustering (DEC)",
  inputs: [Dataset $X = {x_1, ..., x_N} subset bb(R)^d$, clusters $K$, latent dim $p$, Student-$t$ parameter $alpha$, update interval $T_"upd"$, stopping threshold $delta$, step size $eta$.],
  outputs: [Trained encoder $f_theta$, latent centroids $M$, hard cluster labels $hat(y)_1, ..., hat(y)_N$.],
)[
  + Pretrain autoencoder $(f_theta, g_psi)$ on $X$ by minimizing $N^(-1) sum_(i=1)^N norm(x_i - g_psi (f_theta (x_i)))_2^2$; discard decoder $g_psi$.
  + Compute initial embeddings $z_i = f_theta (x_i)$ and initialize $M = {mu_j}_(j=1)^K$ via K-Means on ${z_i}_(i=1)^N$.
  + Compute soft assignments $Q$ via the Student-$t$ kernel and store initial hard labels $hat(y)_i^((0)) = arg max_j q_(i j)$.
  + Repeat for $t = 1, 2, ...$:
    - If $t mod T_"upd" = 0$: recompute $Q$, update $P <- cal(T)(Q)$, and set $hat(y)_i^((t)) = arg max_j q_(i j)$; stop if $N^(-1) sum_i bb(I)[hat(y)_i^((t)) != hat(y)_i^((t-1))] < delta$.
    - Sample mini-batch $cal(B)$ and update $theta <- theta - eta nabla_theta L_"DEC"$ and $M <- M - eta nabla_M L_"DEC"$ with $P$ fixed.
  + Return $(f_theta, M)$ and final labels $hat(y)_i = arg max_j q_(i j)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We derive the exact heavy-tailed gradients and redescending influence bound of the Student-$t$ objective, prove strict entropy sharpening and fixed-point stability of the target operator $cal(T)$, construct an exact zero-loss latent collapse counterexample, and characterize DEC's scale non-invariance.

#subsec("A", "Exact Gradients and Redescending Influence Bound")

#lemma(title: "Closed-Form DEC Gradients and Bounded Outlier Influence", origin: "Classical — Cited")[
  For fixed target $P$, the gradients of $L_"DEC" (theta, M) = sum_(i=1)^N L_i$ with respect to $z_i$ and $mu_j$ are @xie2016dec
  $ (partial L_"DEC") / (partial z_i) &= (alpha + 1) / alpha sum_(j=1)^K (1 + norm(z_i - mu_j)_2^2 / alpha)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j), \
    (partial L_"DEC") / (partial mu_j) &= - (alpha + 1) / alpha sum_(i=1)^N (1 + norm(z_i - mu_j)_2^2 / alpha)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j). $
  Furthermore, the gradient influence of any single embedding $z_i$ on centroid $mu_j$ is uniformly bounded and redescending as $r = norm(z_i - mu_j)_2 -> infinity$:
  $ norm((partial L_i) / (partial mu_j))_2 <= (alpha + 1) / alpha dot.op r / (1 + r^2 / alpha) <= (alpha + 1) / (2 sqrt(alpha)). $
]

#proof[
  Set $s_(i j) = 1 + norm(z_i - mu_j)_2^2 / alpha$, $w_(i j) = s_(i j)^(-(alpha+1)/2)$, and $Z_i = sum_(k=1)^K w_(i k)$ so $q_(i j) = w_(i j) / Z_i$. Because $sum_(j=1)^K p_(i j) = 1$, we write $L_i = "const" - sum_(j=1)^K p_(i j) log w_(i j) + log Z_i$. Differentiating with respect to $z_i$ gives
  $ nabla_(z_i) log w_(i j) = - (alpha + 1) / alpha s_(i j)^(-1) (z_i - mu_j), $
  and by the chain rule,
  $ nabla_(z_i) log Z_i = sum_(j=1)^K q_(i j) nabla_(z_i) log w_(i j). $
  Combining both terms yields
  $ (partial L_i) / (partial z_i) = sum_(j=1)^K (q_(i j) - p_(i j)) nabla_(z_i) log w_(i j) = (alpha + 1) / alpha sum_(j=1)^K s_(i j)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j). $
  For $mu_j$, only $w_(i j)$ depends on $mu_j$, with $nabla_(mu_j) log w_(i j) = -nabla_(z_i) log w_(i j)$, giving the stated expression for $partial L_"DEC" / partial mu_j$.

  For the influence bound, observe $abs(p_(i j) - q_(i j)) <= 1$ since $p_(i j), q_(i j) in [0, 1]$. Setting $r = norm(z_i - mu_j)_2$,
  $ norm((partial L_i) / (partial mu_j))_2 <= (alpha + 1) / alpha dot.op r / (1 + r^2 / alpha). $
  By the AM--GM inequality, $1 + r^2 / alpha >= 2 r / sqrt(alpha)$ with equality at $r = sqrt(alpha)$, so $r / (1 + r^2 / alpha) <= sqrt(alpha) / 2$, yielding the uniform cap $(alpha + 1) / (2 sqrt(alpha))$ (equal to $1$ for $alpha = 1$) and $O(1 / r)$ decay as $r -> infinity$, matching robust Student-$t$ M-estimation @peel2000robust.
]

#subsec("B", "Entropy Sharpening and Stability of the Target Map")

#theorem(title: "Strict Entropy Contraction and Vertex Stability of Target Operator", origin: "Original Proof")[
  Suppose soft cluster frequencies are balanced, $f_1 = ... = f_K = N / K$, so that $cal(T)$ acts row-wise as $p_(i j) = q_(i j)^2 / sum_(k=1)^K q_(i k)^2$ on the simplex $Delta_(K-1)$. Let $H(q) = -sum_(j=1)^K q_j log q_j$ be the Shannon entropy. Then for every $q in Delta_(K-1)$:
  (i) $H(cal(T)(q)) <= H(q)$, with strict inequality unless $q$ is uniform on its support $"supp"(q)$;
  (ii) the fixed points of $cal(T)$ are precisely the uniform distributions on non-empty subsets $S subset.eq {1, ..., K\}$;
  (iii) the $K$ vertices $e_1, ..., e_K$ are quadratically attracting fixed points, whereas every non-vertex fixed point ($abs(S) >= 2$) is linearly unstable with an eigenvalue of $2$.
]

#proof[
  (i) Fix $q in Delta_(K-1)$ with support $S = {j : q_j > 0}$, and define the exponential-family interpolation for $beta > 0$:
  $ q_j^((beta)) = q_j^beta / Z(beta), quad Z(beta) = sum_(k in S) q_k^beta, $
  so $q^((1)) = q$ and $q^((2)) = cal(T)(q)$. Put $u_j = log q_j$ and $psi(beta) = log Z(beta)$. Standard exponential-family identities give $psi'(beta) = bb(E)_(q^((beta))) [u_J]$ and $psi''(beta) = "Var"_(q^((beta))) (u_J) >= 0$. The entropy $h(beta) = H(q^((beta)))$ equals
  $ h(beta) = -sum_(j in S) q_j^((beta)) (beta u_j - psi(beta)) = psi(beta) - beta psi'(beta). $
  Differentiating with respect to $beta$ yields
  $ h'(beta) = psi'(beta) - psi'(beta) - beta psi''(beta) = -beta "Var"_(q^((beta))) (log q_J) <= 0. $
  Integrating over $[1, 2]$ gives
  $ H(cal(T)(q)) - H(q) = - integral_1^2 beta "Var"_(q^((beta))) (log q_J) dif beta <= 0, $
  which is strictly negative unless $log q_J$ is constant on $S$, i.e., $q_j = 1 / abs(S)$ for all $j in S$.

  (ii) Clearly $cal(T)(q) = q$ requires $q_j (q_j - sum_k q_k^2) = 0$ for all $j$, so every $q_j > 0$ equals the common value $sum_k q_k^2 = 1 / abs(S)$.

  (iii) Near vertex $e_a$, let $q_a = 1 - epsilon$ and $q_b = epsilon_b >= 0$ ($b != a$) with $sum_(b != a) epsilon_b = epsilon$. Then for $b != a$,
  $ (cal(T)(q))_b = epsilon_b^2 / ((1 - epsilon)^2 + sum_(c != a) epsilon_c^2) <= epsilon^2 / (1 - epsilon)^2 = O(epsilon^2), $
  proving quadratic local attraction to $e_a$. Conversely, at a face fixed point $q^star$ with $abs(S) = m >= 2$ and $q_j^star = 1 / m$ ($j in S$), consider the tangent perturbation $q_1 = 1/m + delta$, $q_2 = 1/m - delta$ ($1, 2 in S$), and $q_j = 1/m$ otherwise. Then
  $ (cal(T)(q))_1 - (cal(T)(q))_2 = (4 delta / m) / (1 / m + 2 delta^2) = 2 (q_1 - q_2) / (1 + 2 m delta^2). $
  Taking $delta -> 0$ shows the directional derivative along $e_1 - e_2$ equals $2 > 1$, so all interior and face fixed points with $m >= 2$ are exponentially repelling.
]

#subsec("C", "Unregularized Latent Collapse Counterexample")

#theorem(title: "Zero-Loss Representational Collapse Without Decoder Regularization", origin: "Original Counterexample")[
  Let $X = {x_1, ..., x_N} subset bb(R)^d$ consist of distinct points, and let $cal(C) = {C_1, ..., C_K}$ be *any* non-empty $K$-partition of $X$ (even one completely unrelated to the true cluster structure). If the encoder class ${f_theta}$ can interpolate $X$ onto $K$ arbitrary vectors (e.g., a 2-layer `"ReLU"` network of width $Omega(N)$), then
  $ inf_(theta, M) L_"DEC" (theta, M) = 0 $
  is attained asymptotically by embeddings realizing the arbitrary partition $cal(C)$ with zero within-cluster latent variance.
]

#proof[
  Pick $K$ distinct unit vectors $u_1, ..., u_K in bb(R)^p$ with pairwise separation $Delta_"min" = min_(j != j') norm(u_j - u_(j'))_2 > 0$. For any scale $R > 0$, set centroids $mu_j (R) = R u_j$ and choose encoder parameters $theta(R)$ (by scaling the final linear layer by $R$) such that $f_(theta(R)) (x_i) = R u_(c(i))$ where $c(i)$ is the index satisfying $x_i in C_(c(i))$. Then for every $i$:
  $ norm(z_i - mu_(c(i)))_2 = 0, quad norm(z_i - mu_j)_2 >= R Delta_"min" quad (j != c(i)). $
  Hence $w_(i, c(i)) = 1$, while for $j != c(i)$, $w_(i j) <= (1 + R^2 Delta_"min"^2 / alpha)^(-(alpha+1)/2) =: epsilon(R) -> 0$ as $R -> infinity$. Consequently,
  $ 1 / (1 + (K - 1) epsilon(R)) <= q_(i, c(i)) <= 1, quad 0 < q_(i j) <= epsilon(R). $
  Because each $C_j$ is non-empty, $f_j = sum_(i=1)^N q_(i j) -> abs(C_j) >= 1$, whence $p_(i, c(i)) -> 1$ and, for $j != c(i)$, $p_(i j) / q_(i j) = (q_(i j) / f_j) / (sum_k q_(i k)^2 / f_k) <= epsilon(R) (abs(C_"max") / abs(C_"min") + o(1)) -> 0$. Writing $p_(i j) log (p_(i j) / q_(i j)) = q_(i j) dot.op (p_(i j) / q_(i j)) log (p_(i j) / q_(i j))$ with $q_(i j) <= 1$ and applying $lim_(t -> 0^+) t log t = 0$,
  $ L_"DEC" (theta(R), M(R)) = sum_(i=1)^N sum_(j=1)^K p_(i j) log (p_(i j) / q_(i j)) -> 0 quad (R -> infinity). $
  Thus DEC's objective admits zero-loss global minima for arbitrary partitions without preserving ambient geometry, formally explaining why discarding the decoder risks feature-space distortion @guo2017idec.
]

#subsec("D", "Scale Non-Invariance and Non-Lyapunov Alternating Dynamics")

#proposition(title: "Scale Non-Invariance, Zero-Loss Collapse, and Frozen-Target Stationarity", origin: "Original Proof")[
  (i) $L_"DEC"$ is *not* invariant to isotropic latent scaling $(z_i, mu_j) |-> (c z_i, c mu_j)$ ($c > 0$). For every configuration, $lim_(c -> 0^+) L_"DEC" (c Z, c M) = 0$ (uniform collapse), whereas as $c -> infinity$, $q_(i j) (c)$ approaches the scale-free power kernel $d_(i j)^(-(alpha+1)) / sum_k d_(i k)^(-(alpha+1))$ on generic points ($d_(i j) > 0$) and hard 1-NN assignment (with $L_"DEC" -> 0$) whenever every point coincides with its centroid.
  (ii) For fixed $P^((t))$, $L_"DEC"$ is $L$-smooth on compact parameter domains with smooth activations, so diminishing-step SGD reaches an $epsilon$-stationary point in $O(1 / epsilon^2)$ iterations; however, the full alternating map $P <- cal(T)(Q)$ lacks a universal monotonic Lyapunov functional.
]

#proof[
  (i) Under scaling by $c > 0$, $w_(i j) (c) = (1 + c^2 d_(i j)^2 / alpha)^(-(alpha+1)/2)$ where $d_(i j) = norm(z_i - mu_j)_2$. As $c -> 0^+$, $w_(i j) (c) -> 1$ uniformly for all $(i, j)$, so $q_(i j) (c) -> 1/K$, $f_j (c) -> N/K$, $p_(i j) (c) -> 1/K$, and $D_"KL" (P || Q) -> 0$ unconditionally. As $c -> infinity$, if $d_(i j) > 0$ for all $j$, factoring out $(c^2 / alpha)^(-(alpha+1)/2)$ from numerator and denominator yields
  $ lim_(c -> infinity) q_(i j) (c) = d_(i j)^(-(alpha + 1)) / (sum_(k=1)^K d_(i k)^(-(alpha + 1))), $
  which differs strictly from $q_(i j) (1)$ whenever distances vary, proving scale non-invariance at all finite $c > 0$. If instead $d_(i, j^star(i)) = 0$ while $d_(i j) > 0$ ($j != j^star(i)$), then $w_(i, j^star(i)) = 1$ and $w_(i j) (c) -> 0$, driving $q_(i j) (c) -> bb(I)[j = j^star(i)]$ and $L_"DEC" (c Z, c M) -> 0$.

  (ii) For fixed $P^((t))$, $L_"DEC" >= 0$ is twice continuously differentiable in $(z_i, mu_j)$ with bounded Hessian (by Lemma 1), so SGD with step sizes $eta_k = O(1 / sqrt(k))$ guarantees $min_(k <= T) bb(E)[norm(nabla L_"DEC")_2^2] = O(1 / sqrt(T))$. However, replacing $P^((t))$ by $P^((t+1)) = cal(T)(Q^((t+1)))$ changes the first argument of $D_"KL" (P || Q)$: because $Q^((t+1))$ typically has $D_"KL" (cal(T)(Q) || Q) > 0$ (Theorem 1), updating $P$ strictly *increases* the KL loss whenever $Q^((t+1))$ was close to $P^((t))$ yet remains non-degenerate, precluding monotone descent of $L_"DEC"$ across outer updates (unlike exact EM @dempster1977em).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $W_"enc" = sum_(ell=1)^L d_(ell-1) d_ell$ be the encoder forward-pass cost per sample ($d_0 = d, d_L = p$) and $W_"AE" = W_"enc" + W_"dec"$ the full autoencoder cost. Pretraining takes $O(E_"pre" N W_"AE")$, K-Means initialization takes $O(I_"km" N K p)$, and Phase 2 self-training takes $O(E_"ft" N (W_"enc" + K p))$ since the decoder is dropped. Computing the global target $P = cal(T)(Q)$ every $T_"upd"$ steps requires $O(N K p)$ time and $O(N K)$ memory (or streaming frequency estimation for large $N$). Peak training memory for mini-batch size $B$ is $O(B (d + sum_ell d_ell) + W_"enc" + N K)$.

Geometrically, DEC acts as a learned metric warper: pretraining extracts a nonlinear coordinate chart on the data manifold, while KL self-training pulls high-confidence points toward their assigned centroids ($p_(i j) > q_(i j)$) and repels ambiguous points with heavy-tailed, redescending forces (Lemma 1). This lets DEC untangle strongly non-convex shapes into spherical latent clusters, provided pretraining preserves cluster separability.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(E_"pre" N W_"AE" + E_"ft" N (W_"enc" + K p))$],
    [Space Complexity], [$O(B sum_ell d_ell + W_"enc" + N K)$],
    [Cluster Geometry], [Nonlinear latent warp to spherical Student-$t$ balls],
    [Cluster Count $K$], [Fixed a priori (NP-hard init via K-Means @aloise2009np)],
    [Hyperparameters], [Architecture, $K, p, alpha, T_"upd", delta, eta$],
    [Outlier Influence], [Uniformly bounded $<= (alpha+1)/(2 sqrt(alpha))$, $O(1/r)$ decay],
    [Collapse Risk], [High without decoder / early stopping (Thm. 2)],
    [Imbalance Bias], [Moderate; $1 / f_j$ favors balanced cluster sizes],
  ),
  caption: [Operational and geometric profile of Deep Embedded Clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. High-Dimensional Image Clustering.* On handwritten digits and object benchmarks (MNIST, USPS, STL-10) @xie2016dec, pixel-space Euclidean distances fail due to translation and stroke variation. Pretrained autoencoders map images onto compact semantic manifolds where DEC's KL sharpening yields substantial accuracy gains over raw K-Means.

*2. Short-Text and Document Topic Clustering.* Sparse TF-IDF document vectors suffer severe distance concentration @beyer1999nearest. DEC compresses thousands of lexical dimensions into $p approx 10$ dense latent axes, while the heavy-tailed Student-$t$ kernel prevents atypical outlier documents from destabilizing topic centroids (Lemma 1).

*3. Single-Cell Transcriptomics (scRNA-seq).* Gene expression profiles exhibit dropout noise and nonlinear cell-state manifolds. DEC-family models (such as scDeepCluster) leverage heavy-tailed latent kernels to absorb technical expression spikes while sharpening discrete cell-type assignments.

*4. Customer and Behavioral Segmentation.* Heterogeneous purchasing and telemetry features rarely follow Gaussian clusters in raw coordinates. Nonlinear encoding extracts latent behavioral factors, and the $1 / f_j$ target normalization prevents a single majority customer segment from absorbing smaller niches.

*5. Industrial Acoustic and Vibration Diagnostics.* High-dimensional spectrogram frames from rotating machinery lie near low-dimensional fault manifolds; DEC separates nonlinear operating modes into well-isolated latent clusters without fault labels.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four mathematically grounded failure modes govern DEC deployment: (i) *latent distortion / collapse* (Theorem 2), where excessive fine-tuning destroys fine-grained structure; (ii) *severe class imbalance*, because dividing by $f_j$ upweights tiny clusters and splits large skewed modes; (iii) *initialization lock-in*, since self-training amplifies initial K-Means mistakes quadratically ($p_(i j) / p_(i k) prop (q_(i j) / q_(i k))^2$); and (iv) *scale drift* (Proposition 1), where unconstrained weight norms either collapse $Z -> 0$ or saturate $Q$ prematurely.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Image Features], [Non-convex manifold unwrapping], [Keep decoder (IDEC) or use SwAV @caron2020swav],
    [TF-IDF Text], [Mitigates high-D concentration], [Imbalanced topics; use VaDE @jiang2017vade],
    [scRNA-seq Data], [Redescending Student-$t$ resists spikes], [Continuous trajectories; use IDEC],
    [Tabular Skewed], [Latent compression + soft sharpening], [Monitor $delta$ early stopping],
    [Small-Sample Low-D], [Over-parameterized; risk of collapse], [Use GCS, GMM, or K-Means],
  ),
  caption: [Domain suitability and failure-mode matrix for DEC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Formal Link Between DEC, Student-t EM, and Deep Clustering Variants", origin: "Original Unified Proof")[
  DEC sits at the intersection of robust mixture modeling and deep representation learning:
  (i) Freezing $theta$ and replacing $p_(i j) = cal(T)(Q)_(i j)$ by the exact Bayesian posterior $r_(i j)$ recovers an M-step of EM for an isotropic Student-$t$ mixture @peel2000robust @dempster1977em with uniform mixing weights; DEC modifies EM by inserting quadratic sharpening ($q_(i j)^2$), frequency balancing ($1 / f_j$), and joint backpropagation into $f_theta$.
  (ii) Adding reconstruction penalty $L_r (theta, psi)$ yields IDEC @guo2017idec, eliminating the zero-loss collapse path of Theorem 2 under bounded decoder Lipschitz constants.
  (iii) Replacing $D_"KL" (P || Q)$ by hard latent K-Means $norm(z_i - mu_(c(i)))_2^2$ alongside $L_r$ yields DCN @yang2017dcn, while placing a Gaussian mixture prior inside a variational autoencoder @kingma2014vae yields VaDE @jiang2017vade.
]

#proof[
  For an isotropic Student-$t$ mixture with density $p(z | mu_j) prop (1 + norm(z - mu_j)_2^2 / alpha)^(-(alpha+p)/2)$, the EM complete-data log-likelihood differentiates with respect to $mu_j$ to give $sum_i r_(i j) u_(i j) (z_i - mu_j) = 0$ with Mahalanobis weights $u_(i j) = (alpha + p) / (alpha + norm(z_i - mu_j)_2^2)$ @peel2000robust @bishop2006pattern. Comparing with Lemma 1 shows that DEC shares the identical redescending weight factor $(1 + norm(z_i - mu_j)_2^2 / alpha)^(-1)$, differing in three explicit operations: fixing exponent $(alpha + 1)/2$ independent of $p$ (avoiding high-D underflow), replacing $r_(i j)$ by the sharpened residual $p_(i j) - q_(i j)$, and propagating gradients into $theta$. Setting $L_"IDEC" = L_r + gamma L_"DEC"$ recovers DEC in the formal limit $gamma -> infinity$ (rescaling by $1/gamma$), while replacing soft KL by hard assignment $s_i in {0, 1}^K$ and squared latent error yields $L_"DCN"$ directly.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark, DEC obtains $"Sil" = 0.445$ @rousseeuw1987silhouettes and $"ARI" = 0.574$ @hubert1985comparing on *Iris*, $"Sil" = 0.344$ on *Customers*, $"Sil" = 0.718$ ($"ARI"$ N/A) on *Multi-Shape*, and $"Sil" = 0.295$ ($"ARI"$ N/A) on *High-D*, with a mean runtime of *10.574 ms*. These numbers vividly confirm DEC's theoretical geometry. On *Multi-Shape*, DEC achieves the highest Silhouette ($0.718$) across Algorithms 180--183, substantially exceeding GCS ($0.609$) and IDEC ($0.563$): discarding the decoder frees $f_theta$ to aggressively unwrap non-convex rings and crescents into tightly separated Student-$t$ spheres (Theorem 1). Likewise, DEC leads on *Customers* ($"Sil" = 0.344$), benefiting from both nonlinear feature compression and redescending robustness to spending outliers (Lemma 1). Conversely, on *Iris* ($"ARI" = 0.574$) and *High-D* ($"Sil" = 0.295$), DEC slightly trails reconstruction-regularized IDEC ($"ARI" = 0.584$, *High-D* $"Sil" = 0.300$) and DCN ($"ARI" = 0.648$, *High-D* $"Sil" = 0.298$), reflecting the mild feature-distortion risk proved in Theorem 2 when no decoder constrains latent updates. Finally, its $10.574$ ms runtime reflects the two-phase pretraining and backpropagation overhead shared across autoencoder clustering architectures @caron2018deepcluster @ghasedi2017depict.
