#import "../template/ieee.typ": *

#alg_entry(
  id: 86,
  name: "CLASSIT",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1989",
  paradigm: "Incremental Continuous Conceptual Clustering",
  complexity: $cal(O)(n B^2 d h)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

CLASSIT @gennari1989classit extends COBWEB @fisher1987cobweb from discrete nominal streams to real-valued Euclidean observations while preserving incremental, unsupervised concept-hierarchy formation. Let $cal(D) = (bold(x)_1, dots, bold(x)_n)$ be an ordered stream of $n$ instances in $bb(R)^d$, with $bold(x)_t = (x_(t, 1), dots, x_(t, d))^T$. Rather than storing discrete frequency histograms, each node $C$ in the concept tree $cal(T)$ models its $d$ continuous attributes as conditionally independent univariate Gaussians $A_i | C ~ cal(N)(mu_(C, i), sigma_(C, i)^2)$, summarized solely by the online sufficient statistics $(n_C, bold(mu)_C, bold(M)_(2, C))$, where $bold(mu)_C = n_C^(-1) sum_(bold(x) in C) bold(x)$, $M_(2, C, i) = sum_(bold(x) in C) (x_i - mu_(C, i))^2$, and empirical standard deviation $sigma_(C, i) = sqrt(M_(2, C, i) / n_C)$.

To generalize Gluck and Corter's discrete Category Utility to continuous densities, Gennari, Langley, and Fisher @gennari1989classit replace each discrete sum of squared probabilities $sum_j P(A_i = v_(i, j) | C)^2$ by the $L^2$-norm integral of the univariate Gaussian density $f_(C, i) (u) = (2 pi sigma_(C, i)^2)^(-1/2) exp(-(u - mu_(C, i))^2 / (2 sigma_(C, i)^2))$:
$ integral_(-infinity)^(+infinity) f_(C, i) (u)^2 dif u = 1 / (2 sqrt(pi) sigma_(C, i)). $
Because a singleton node satisfies $sigma_(C, i) = 0$ and causes $1 / sigma_(C, i)$ to diverge, CLASSIT introduces a positive perceptual _acuity_ floor $sigma_"ac" > 0$ and a _cutoff_ pruning threshold $theta_"cut" >= 0$.

#definition(title: "Regularized Continuous Category Utility")[
  Given acuity parameter $sigma_"ac" > 0$, define the acuity-regularized standard deviation $tilde(sigma)_(C, i) = max(sigma_(C, i), sigma_"ac")$. For an internal node $C$ with parent deviations $tilde(sigma)_(0, i)$ and child partition $cal(P) = {C_1, dots, C_K}$ with weights $P(C_k) = n_(C_k) / n_C$, the _Continuous Category Utility_ of $cal(P)$ is
  $ "CU"_"c"(cal(P)) = 1 / (2 sqrt(pi) K) sum_(k=1)^K P(C_k) sum_(i=1)^d ( 1 / tilde(sigma)_(k, i) - 1 / tilde(sigma)_(0, i) ). $
]

Because the constant $(2 sqrt(pi))^(-1) > 0$ is shared across all candidate partitions of $C$, implementations often omit it when ranking operators, while retaining it preserves the exact continuous-limit interpretation proved below. CLASSIT uses the same four local hill-climbing operators as COBWEB—*Incorporate*, *Create*, *Merge*, and *Split*—augmented by a *Cutoff* rule that halts recursion whenever incorporating $bold(x)$ into the best host yields negligible utility change relative to $theta_"cut"$, thereby pruning redundant noisy subtrees.

#algbox(
  title: "CLASSIT Incremental Gaussian Concept Tree",
  inputs: [Stream $cal(D) = (bold(x)_1, dots, bold(x)_n) subset bb(R)^d$, acuity $sigma_"ac" > 0$, cutoff $theta_"cut" >= 0$.],
  outputs: [Concept tree $cal(T)$ with Gaussian sufficient statistics $(n_C, bold(mu)_C, bold(M)_(2, C))$ per node.],
)[
  + Initialize root $R$ with $(1, bold(x)_1, bold(0))$; for $t = 2, dots, n$, invoke `Insert`$(R, bold(x)_t)$.
  + *Procedure* `Insert`$(C, bold(x))$: if $C$ is a leaf, either absorb $bold(x)$ into $C$ if $max_i |x_i - mu_(C, i)| <= theta_"cut"$ or expand $C$ into two leaf children holding the prior node statistics and ${bold(x)}$, update $C$, and return.
  + Otherwise ($C$ is internal with children $cal(P) = {C_1, dots, C_K}$), update $(n_C, bold(mu)_C, bold(M)_(2, C))$ with $bold(x)$.
  + *Incorporate*: for each $k in {1, dots, K}$, compute $"CU"_k = "CU"_"c"(cal(P)_(k <- bold(x)))$ using updated sufficient statistics; let $k_1, k_2$ index the top two scores.
  + *Create*: evaluate $"CU"_"new" = "CU"_"c"(cal(P) union {{bold(x)}})$ where the new singleton has $tilde(sigma)_(K+1, i) = sigma_"ac"$.
  + *Merge*: if $K >= 2$, form $C_m = C_(k_1) union C_(k_2)$ via exact sufficient-statistic fusion and evaluate $"CU"_"merge"$ with $bold(x)$ added to $C_m$.
  + *Split*: if $C_(k_1)$ is internal, promote its children into $C$ and evaluate $"CU"_"split"$ after re-incorporating $bold(x)$.
  + Execute the highest-scoring operator. If *Incorporate* wins and the normalized score gain falls below $theta_"cut"$, halt recursion (pruning); otherwise recurse on the selected host (or re-enter `Insert`$(C, bold(x))$ after *Split*).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results for CLASSIT: the exact Riemann-limit derivation from discrete COBWEB and Jensen non-negativity in the unclipped regime, a constructive proof of singleton divergence and acuity-induced sign reversal, an anisotropic scale-sensitivity theorem proving inverse-scale dominance, exact $cal(O)(d)$ sufficient-statistic fusion identities, and a rigorous comparison against diagonal Gaussian differential entropy.

#subsec("A", "Continuous Limit and Unclipped Non-Negativity")

#theorem(title: "Exact Continuum Limit and Jensen Non-Negativity", origin: "Original Proof")[
  1. For any Gaussian density $f(u) = (2 pi sigma^2)^(-1/2) exp(-(u - mu)^2 / (2 sigma^2))$ with $sigma > 0$, uniform histogram discretization into bins $I_j = [j Delta, (j + 1) Delta)$ satisfies
     $ lim_(Delta -> 0^+) 1 / Delta sum_(j=-infinity)^(+infinity) (integral_(I_j) f(u) dif u)^2 &= integral_(-infinity)^(+infinity) f(u)^2 dif u \
       &= 1 / (2 sqrt(pi) sigma). $
     Consequently, for any continuous Gaussian mixture partition $cal(P)$ with positive component deviations $sigma_(k, i) > 0$, the rescaled discrete Category Utility converges exactly: $lim_(Delta -> 0^+) Delta^(-1) "CU"_"disc"(cal(P)_Delta) = "CU"_"c"(cal(P))|_(sigma_"ac" = 0)$.
  2. Whenever $sigma_(k, i) > 0$ for all $k, i$ and $sigma_"ac" = 0$, the unclipped Continuous Category Utility satisfies $"CU"_"c"(cal(P)) >= 0$, with equality if and only if $mu_(k, i) = mu_(0, i)$ and $sigma_(k, i) = sigma_(0, i)$ for every $k, i$.
]

#proof[
  For Claim 1, by the mean-value theorem on each bin $I_j$, there exists $xi_j in I_j$ such that $integral_(I_j) f(u) dif u = Delta f(xi_j)$. Therefore
  $ 1 / Delta sum_(j=-infinity)^(+infinity) (integral_(I_j) f(u) dif u)^2 = sum_(j=-infinity)^(+infinity) Delta f(xi_j)^2. $
  Because $f^2$ is smooth, bounded, and rapidly decaying, this Riemann sum converges to $integral_(-infinity)^(+infinity) f(u)^2 dif u$ as $Delta -> 0^+$. Evaluating the Gaussian integral via the substitution $z = (u - mu) / sigma$ gives
  $ integral_(-infinity)^(+infinity) f(u)^2 dif u &= 1 / (2 pi sigma^2) integral_(-infinity)^(+infinity) exp(-(u - mu)^2 / sigma^2) dif u \
    &= 1 / (2 pi sigma) integral_(-infinity)^(+infinity) exp(-z^2) dif z = 1 / (2 sqrt(pi) sigma). $
  Applying this identity to each child and parent attribute density yields the exact rescaled limit of COBWEB's Category Utility.

  For Claim 2, fix attribute $i$ and denote $w_k = P(C_k)$. By the law of total variance,
  $ sigma_(0, i)^2 &= sum_(k=1)^K w_k sigma_(k, i)^2 + sum_(k=1)^K w_k (mu_(k, i) - mu_(0, i))^2 \
    &>= sum_(k=1)^K w_k sigma_(k, i)^2. $
  Consider $g(v) = v^(-1/2)$ on $(0, infinity)$. Since $g'(v) = -(1/2) v^(-3/2) < 0$ and $g''(v) = (3/4) v^(-5/2) > 0$, $g$ is strictly decreasing and strictly convex. Applying Jensen's inequality followed by monotonicity of $g$ yields
  $ sum_(k=1)^K w_k 1 / sigma_(k, i) &= sum_(k=1)^K w_k g(sigma_(k, i)^2) \
    &>= g(sum_(k=1)^K w_k sigma_(k, i)^2) >= g(sigma_(0, i)^2) = 1 / sigma_(0, i). $
  The first inequality is tight if and only if all $sigma_(k, i)$ are equal, and the second is tight if and only if $sum_k w_k (mu_(k, i) - mu_(0, i))^2 = 0$, i.e., $mu_(k, i) = mu_(0, i)$ for all $k$. Summing over $i = 1, dots, d$ and multiplying by $(2 sqrt(pi) K)^(-1)$ proves Claim 2.
]

#subsec("B", "Singleton Singularity and Acuity-Induced Negativity")

Although Theorem 1 guarantees non-negativity without acuity clipping, setting $sigma_"ac" = 0$ is catastrophic algorithmically, whereas imposing $sigma_"ac" > 0$ subtly breaks the convexity needed for non-negativity.

#proposition(title: "Singularity Without Acuity and Acuity-Induced Sign Reversal", origin: "Original Counterexample")[
  1. If $sigma_"ac" = 0$, any partition $cal(P)$ containing at least one singleton (or zero-variance) cluster on any attribute yields $"CU"_"c"(cal(P)) = +infinity$, causing CLASSIT to degenerate into a flat star of singletons under *Create*.
  2. Conversely, for any $sigma_"ac" > 0$, the clipped map $tilde(g)(v) = max(v, sigma_"ac"^2)^(-1/2)$ loses global convexity, and $"CU"_"c"(cal(P))$ can become strictly negative. Specifically, on the 1D dataset $cal(D) = {0, 0, -sqrt(2), sqrt(2)}$ with $sigma_"ac" = 1$, the partition $cal(P) = {{0, 0}, {-sqrt(2), sqrt(2)}}$ achieves $"CU"_"c"(cal(P)) = (sqrt(2) - 2) / (16 sqrt(pi)) < 0$.
]

#proof[
  Claim 1 is immediate: for a singleton child $C_k = {bold(x)}$, $M_(2, k, i) = 0$, so $sigma_(k, i) = 0$ and $P(C_k) / sigma_(k, i) = +infinity$ whenever $P(C_k) > 0$. Since *Create* always forms a singleton with weight $1 / n_C > 0$, its unclipped Category Utility is $+infinity$, dominating any non-degenerate *Incorporate* or *Merge* partition.

  For Claim 2, consider $cal(D) = {0, 0, -sqrt(2), sqrt(2)} subset bb(R)$ ($d = 1$) partitioned into $C_1 = {0, 0}$ and $C_2 = {-sqrt(2), sqrt(2)}$ with equal weights $w_1 = w_2 = 1/2$. The parent has mean $mu_0 = 0$ and empirical variance $sigma_0^2 = (0 + 0 + 2 + 2) / 4 = 1$, so $sigma_0 = 1$. Meanwhile, $sigma_1 = 0$ and $sigma_2 = sqrt((2 + 2) / 2) = sqrt(2)$. Setting $sigma_"ac" = 1$ gives clipped deviations $tilde(sigma)_0 = max(1, 1) = 1$, $tilde(sigma)_1 = max(0, 1) = 1$, and $tilde(sigma)_2 = max(sqrt(2), 1) = sqrt(2)$. Substituting into Definition 1 with $K = 2$ yields
  $ "CU"_"c"(cal(P)) &= 1 / (4 sqrt(pi)) [ 1 / 2 (1 / 1 - 1 / 1) + 1 / 2 (1 / sqrt(2) - 1 / 1) ] \
    &= 1 / (8 sqrt(pi)) (sqrt(2) / 2 - 1) = (sqrt(2) - 2) / (16 sqrt(pi)) approx -0.0206 < 0. $
  Geometrically, clipping replaces $v^(-1/2)$ by a flat cap $sigma_"ac"^(-1)$ on $[0, sigma_"ac"^2]$, creating a non-convex elbow at $v = sigma_"ac"^2$ where high-variance clusters are penalized while zero-variance clusters receive no compensatory reward beyond $sigma_"ac"^(-1)$.
]

#subsec("C", "Anisotropic Scale Sensitivity and Inverse-Scale Dominance")

A crucial geometric property separates CLASSIT from Euclidean variance-minimization methods such as $k$-means and Ward clustering @ward1963hierarchical.

#theorem(title: "Inverse-Scale Dominance and Standardization Invariance", origin: "Original Proof")[
  Consider an invertible diagonal affine transformation $x_i' = a_i x_i + b_i$ with scales $a_i > 0$ and shifts $b_i in bb(R)$ for $i = 1, dots, d$.
  1. CLASSIT is strictly invariant to all translations $(b_1, dots, b_d)$.
  2. In the unclipped regime ($min(a_i, 1) sigma_(k, i) >= sigma_"ac"$), each attribute's contribution to $"CU"_"c"(cal(P))$ scales as $a_i^(-1)$:
     $ "CU"_"c"(cal(P)') = 1 / (2 sqrt(pi) K) sum_(k=1)^K P(C_k) sum_(i=1)^d 1 / a_i ( 1 / sigma_(k, i) - 1 / sigma_(0, i) ). $
     Consequently, shrinking the numerical units of attribute $i$ by $a_i << 1$ amplifies its weight in $"CU"_"c"$ by $a_i^(-1) >> 1$—the exact inverse of Euclidean SSE, where large-scale features dominate quadratically as $a_i^2$.
  3. Standardizing each feature via $hat(x)_i = (x_i - mu_(0, i)) / sigma_(0, i)$ with dimensionless acuity $alpha_"ac" > 0$ (equivalently setting per-feature acuity $sigma_("ac", i) = alpha_"ac" sigma_(0, i)$) renders CLASSIT strictly invariant to all positive diagonal rescalings.
]

#proof[
  Under $x_i' = a_i x_i + b_i$, every cluster mean shifts as $mu_(k, i)' = a_i mu_(k, i) + b_i$, so $(x_i' - mu_(k, i)')^2 = a_i^2 (x_i - mu_(k, i))^2$ and $sigma_(k, i)' = a_i sigma_(k, i)$, independent of $b_i$; this proves Claim 1. In the unclipped regime where $tilde(sigma)_(k, i)' = a_i sigma_(k, i)$ and $tilde(sigma)_(0, i)' = a_i sigma_(0, i)$, substituting directly into Definition 1 factors $1 / a_i$ out of the $i$-th summand, proving Claim 2. Notice that even under isotropic scaling $a_i = a$ with fixed $sigma_"ac"$, varying $a$ shifts the effective dimensionless acuity $sigma_"ac" / (a sigma_(0, i))$, altering which nodes hit the acuity floor! Finally, defining $sigma_("ac", i) = alpha_"ac" sigma_(0, i)$ yields
  $ 1 / max(sigma_(k, i)', alpha_"ac" sigma_(0, i)') - 1 / sigma_(0, i)' \
    = 1 / (a_i sigma_(0, i)) ( 1 / max(sigma_(k, i) / sigma_(0, i), alpha_"ac") - 1 ), $
  and multiplying each term by the normalized feature scale $a_i sigma_(0, i)$ cancels $a_i$ identically across all $i$, proving Claim 3.
]

#subsec("D", "Exact Online Sufficient-Statistic Fusion and Complexity")

#lemma(title: "Numerically Stable O(d) Merge and Insertion Identities", origin: "Original Derivation")[
  Represent each node $C$ by the triple $(n_C, bold(mu)_C, bold(M)_(2, C))$, where $M_(2, C, i) = sum_(bold(x) in C) (x_i - mu_(C, i))^2$. Then merging any two disjoint nodes $C_a, C_b$ into $C_(a b) = C_a union C_b$ satisfies the exact identities
  $ n_(a b) = n_a + n_b, quad bold(mu)_(a b) = (n_a bold(mu)_a + n_b bold(mu)_b) / n_(a b), $
  $ M_(2, a b, i) = M_(2, a, i) + M_(2, b, i) + (n_a n_b) / n_(a b) (mu_(a, i) - mu_(b, i))^2. $
  Consequently, every single-point incorporation ($n_b = 1, M_(2, b, i) = 0$) and every sibling merge executes in $cal(O)(d)$ time without raw-instance access, yielding total stream complexity $cal(O)(n B^2 d h)$ naively and $cal(O)(n B d h)$ when caching child sums $S_k = sum_(i=1)^d tilde(sigma)_(k, i)^(-1)$ @gennari1989classit.
]

#proof[
  The identities for $n_(a b)$ and $bold(mu)_(a b)$ follow immediately from linearity of summation. For each coordinate $i$, decompose $M_(2, a b, i)$ by adding and subtracting $mu_(a, i)$ on $C_a$ and $mu_(b, i)$ on $C_b$:
  $ sum_(bold(x) in C_a) (x_i - mu_(a b, i))^2 = M_(2, a, i) + n_a (mu_(a, i) - mu_(a b, i))^2, $
  because $sum_(bold(x) in C_a) (x_i - mu_(a, i)) = 0$. Similarly on $C_b$. Substituting $mu_(a, i) - mu_(a b, i) = (n_b / n_(a b))(mu_(a, i) - mu_(b, i))$ and $mu_(b, i) - mu_(a b, i) = -(n_a / n_(a b))(mu_(a, i) - mu_(b, i))$ gives
  $ n_a (mu_(a, i) - mu_(a b, i))^2 + n_b (mu_(b, i) - mu_(a b, i))^2 &= (n_a n_b^2 + n_b n_a^2) / n_(a b)^2 (mu_(a, i) - mu_(b, i))^2 \
    &= (n_a n_b) / n_(a b) (mu_(a, i) - mu_(b, i))^2. $
  Summing the $C_a$ and $C_b$ parts establishes the exact merge law, avoiding catastrophic cancellation inherent in naive $sum x_i^2 - n mu_i^2$ formulas. At each of at most $h$ levels, evaluating *Incorporate* across $K <= B$ hosts naively recomputes $K$ sums of $d$ inverse deviations in $cal(O)(B^2 d)$ time; caching $S_k$ and $sum_(k=1)^K w_k S_k$ updates only the modified host in $cal(O)(d)$ time, lowering per-level cost to $cal(O)(B d)$ and total stream complexity over $n$ points to $cal(O)(n B d h)$.
]

#subsec("E", "Comparison with Diagonal Gaussian Differential Entropy")

#proposition(title: "Polynomial vs. Logarithmic Precision Bias Relative to Diagonal GMM", origin: "Original Proof")[
  For a fixed $K$-way hard partition $cal(P)$ with positive deviations $sigma_(k, i) > 0$, the empirical diagonal-Gaussian log-likelihood gain over a single Gaussian equals the Shannon differential-entropy reduction
  $ Delta cal(L)_"GMM"(cal(P)) = sum_(k=1)^K w_k sum_(i=1)^d log (sigma_(0, i) / sigma_(k, i)), $
  whereas CLASSIT maximizes $K "CU"_"c"(cal(P)) = (2 sqrt(pi))^(-1) sum_(k=1)^K w_k sum_(i=1)^d (sigma_(k, i)^(-1) - sigma_(0, i)^(-1))$. Although both functions $psi_"G"(v) = -(1/2) log v$ and $psi_"CL"(v) = (2 sqrt(pi))^(-1) v^(-1/2)$ are strictly decreasing and strictly convex in variance $v = sigma^2$, their ratio $psi_"CL"(v) / psi_"G"(v) -> +infinity$ as $v -> 0^+$. Hence CLASSIT exerts an asymptotically unbounded polynomial preference for ultra-tight micro-clusters relative to diagonal Gaussian MLE.
]

#proof[
  Maximizing the diagonal Gaussian log-likelihood of cluster $C_k$ over $(mu_(k, i), sigma_(k, i)^2)$ gives $hat(ell)_k = -(n_(C_k) / 2) sum_(i=1)^d (1 + log(2 pi sigma_(k, i)^2))$. Dividing $sum_k hat(ell)_k - hat(ell)_0$ by $n$ cancels the constants $1 + log(2 pi)$ and yields $Delta cal(L)_"GMM"(cal(P)) = sum_k w_k sum_i log(sigma_(0, i) / sigma_(k, i))$. Both $psi_"G"(v)$ and $psi_"CL"(v)$ have negative first and positive second derivatives on $(0, infinity)$, confirming identical directionality on homogeneous-variance partitions. However, by L'Hôpital's rule as $v -> 0^+$,
  $ lim_(v -> 0^+) (v^(-1/2)) / (-(1/2) log v) &= lim_(v -> 0^+) (-(1/2) v^(-3/2)) / (-(1/2) v^(-1)) \
    &= lim_(v -> 0^+) v^(-1/2) = +infinity. $
  Thus, whenever one subset of points lies densely around a mode while another is diffuse, CLASSIT assigns disproportionately greater reward to isolating the dense micro-cluster near $sigma_"ac"$ than Gaussian log-likelihood does, explaining CLASSIT's heightened sensitivity to $sigma_"ac"$ and tendency to over-fragment heteroscedastic data.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Replacing COBWEB's $V$-bin frequency tables by two scalar sufficient statistics $(mu_(C, i), M_(2, C, i))$ per dimension yields major computational and memory gains: space shrinks from $cal(O)(|cal(T)| d V)$ to $cal(O)(|cal(T)| d)$ and per-instance time drops by a factor of $V$ to $cal(O)(B d h)$ with caching (Lemma 1). Furthermore, every update relies strictly on contiguous float arrays, achieving high CPU cache locality and eliminating hash-map overhead.

Geometrically, CLASSIT induces axis-aligned hyper-ellipsoidal concepts parameterized by diagonal Gaussian contours. Because the acuity floor $sigma_"ac"$ acts globally across coordinates, unstandardized features with differing physical units distort the geometry via the $a_i^(-1)$ law of Theorem 2; z-score normalization (or online running standardization) is therefore mandatory in practice. As with COBWEB, local greedy tree restructuring inherits presentation-order dependence and lacks global optimality guarantees.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Continuous real vectors $bold(x) in bb(R)^d$],
    [Objective Function], [Acuity-regularized Continuous Category Utility $"CU"_"c"(cal(P))$],
    [Time Complexity], [$cal(O)(n B^2 d h)$ naive; $cal(O)(n B d h)$ with cached inverse-deviation sums],
    [Space Complexity], [$cal(O)(|cal(T)| d)$ storing $(n_C, bold(mu)_C, bold(M)_(2, C))$ per tree node],
    [Hyperparameters], [Perceptual acuity $sigma_"ac" > 0$ and cutoff threshold $theta_"cut" >= 0$],
    [Invariances], [Translation-invariant; *not* scale-invariant unless standardized; order-sensitive],
    [Cluster Geometry], [Axis-aligned Gaussian concepts (diagonal covariance inductive bias)],
    [Streaming Support], [Native single-pass online insertion with $cal(O)(d)$ sufficient-statistic fusion],
  ),
  caption: [Operational and geometric profile of CLASSIT.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Continuous Sensor-Stream Concept Discovery:* Industrial IoT and environmental monitoring arrays emit continuous multivariate streams (temperature, vibration RMS, pressure, acoustic energy) that cannot be buffered in full. CLASSIT builds an anytime hierarchy of operating regimes using $cal(O)(d)$ sufficient statistics per node, while the cutoff parameter $theta_"cut"$ filters sub-sensor-noise fluctuations directly at the acuity level of the instrumentation.

*2. Cognitive Modeling of Perceptual Categorization:* Originally validated on continuous biological and geometric stimuli @gennari1989classit, CLASSIT captures perceptual acuity limits directly through $sigma_"ac"$—modeling the psychophysical just-noticeable difference (JND) below which human observers treat distinct physical measurements as perceptually equivalent.

*3. Robotic Skill and Tactile Regime Taxonomies:* Manipulation policies stream low-dimensional continuous joint torques, contact forces, and end-effector velocities. CLASSIT incrementally organizes contact states into hierarchical motion primitives with explicit per-joint Gaussian envelopes $(mu_(C, i), sigma_(C, i))$, allowing fast online regime recognition and anomalous-contact detection by checking Mahalanobis z-scores along the tree path.

*4. Real-Time Telemetry and Cloud Workload Summarization:* Datacenter monitors track continuous CPU utilization, memory pressure, I/O latency, and queue depth across thousands of microservices. CLASSIT's lightweight $cal(O)(B d h)$ updates and hierarchical abstraction let operators inspect coarse load archetypes near the root or drill into granular microservice profiles near the leaves without re-clustering historical logs.

*5. Adaptive Database Query-Plan Indexing:* Continuous workload descriptors (selectivity estimates, join cardinalities, buffer hit rates) arrive sequentially as queries execute. CLASSIT groups similar query profiles into a dynamic concept tree whose leaf Gaussian summaries predict runtime and memory footprints for incoming plans via top-down tree traversal.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CLASSIT exhibits four principal failure modes: (i) *unscaled heterogeneous features*, where Theorem 2 shows tiny-variance features hijack $"CU"_"c"$ inversely as $a_i^(-1)$; (ii) *acuity miscalibration*, where too small $sigma_"ac"$ triggers singleton fragmentation (Proposition 1) and too large $sigma_"ac"$ causes acuity clipping and negative utility scores; (iii) *correlated or oblique Gaussian clusters*, where diagonal independence cannot rotate axes along covariance eigenvectors; and (iv) *high-dimensional isotropic noise* @beyer1999nearest, where additive $1 / tilde(sigma)_i$ terms across uninformative coordinates drown out low-dimensional cluster separation.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Standardized low-D sensor streams], [High: fast $cal(O)(d)$ updates, natural acuity matched to sensor noise], [Standardize online; tune $sigma_"ac"$ to noise floor],
    [Psychophysical / cognitive studies], [High: direct formalization of JND acuity and hierarchical concepts], [Average over shuffled orders to verify stability],
    [Correlated ellipsoidal tabular data], [Low: diagonal bias splits tilted clusters; order-sensitive], [Use AutoClass (Alg. 87), MCLUST, or TMM (Alg. 88)],
    [Heavy-tailed / outlier-rich telemetry], [Poor: squared sufficient statistics $M_2$ explode under outliers], [Use Student's $t$-mixtures (Alg. 88) or trimmed methods],
    [High-D embeddings / gene expression], [Poor: additive inverse-std dilution and $h$ tree skew], [Use subspace, spectral, or spherical clustering],
  ),
  caption: [Domain suitability matrix and failure mitigations for CLASSIT.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CLASSIT forms the continuous bridge between symbolic conceptual trees and modern streaming micro-cluster hierarchies. Its sufficient-statistic node summary $(n_C, bold(mu)_C, bold(M)_(2, C))$ directly foreshadowed the Clustering Feature triple $"CF" = (n, bold("LS"), "SS")$ of BIRCH @zhang1996birch and the anytime hierarchical index of ClusTree @kranen2011clustree—differing primarily in that BIRCH tracks isotropic Euclidean radius ($"SS" = sum norm(bold(x))_2^2$) whereas CLASSIT tracks coordinate-wise variances $M_(2, C, i)$ and evaluates probabilistic Category Utility rather than radius thresholds.

Relative to COBWEB @fisher1987cobweb (Alg. 85), Theorem 1 establishes CLASSIT as the exact $Delta -> 0^+$ continuum limit under Gaussian densities. Relative to AutoClass @cheeseman1996autoclass (Alg. 87) and diagonal Gaussian mixtures @fraley2002mclust @mclachlan2000finite, Proposition 2 shows that CLASSIT trades global EM optimization and Shannon log-likelihood for greedy single-pass speed and a sharper $sigma^(-1)$ variance penalty. Online tree methods such as PERCH @kobren2017perch and GRINCH @monath2019grinch address CLASSIT's order sensitivity using masking and tree-rotation invariants.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* On the standardized 4-dataset suite, CLASSIT achieves Iris Silhouette $= 0.162$ and $"ARI" = 0.066$ @hubert1985comparing @rousseeuw1987silhouettes, Customers Silhouette $= 0.318$, Multi-Shape Silhouette $= 0.402$ ($"ARI" =$ N/A), and High-D Silhouette $= 0.144$ ($"ARI" =$ N/A), with a mean runtime of just $2.460$ ms. These empirical numbers align closely with our proofs:
1. *Speed Advantage:* At $2.460$ ms, CLASSIT is $4.6 times$ faster than COBWEB ($11.350$ ms), $7.8 times$ faster than AutoClass ($19.192$ ms), and $9.3 times$ faster than TMM ($22.958$ ms), directly validating the lightweight $cal(O)(B d h)$ sufficient-statistic machinery of Lemma 1.
2. *Continuous vs. Discretized Geometry:* Removing nominal binning raises Iris Silhouette nearly eightfold over COBWEB ($0.162$ vs. $0.021$) and edges out COBWEB on Customers ($0.318$ vs. $0.312$). Yet Iris $"ARI"$ remains low ($0.066$) because diagonal Gaussian independence fails to model the strong sepal-petal correlation, while Proposition 2's $sigma^(-1)$ micro-cluster bias and greedy order sensitivity fragment the overlapping _versicolor_ / _virginica_ manifold.
3. *High-Dimensional and Non-Convex Degradation:* On High-D ($"Sil" = 0.144$) and Multi-Shape ($"Sil" = 0.402$), additive coordinate-wise $tilde(sigma)_(k, i)^(-1)$ scoring suffers from inverse-scale distortion (Theorem 2) and inability to represent curved or elongated manifolds, motivating full Bayesian and heavy-tailed mixture estimators in Algs. 87 and 88.
