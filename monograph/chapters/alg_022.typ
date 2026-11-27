#import "../template/ieee.typ": *

#alg_entry(
  id: 22,
  name: "Forgy's K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "1965",
  paradigm: "Batch Alternating Partitioning with Random Data Seeding",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset with at least $k$ distinct observations, and let $k in {1, dots, n}$ be the target number of clusters. Independently of Lloyd's unpublished 1957 Bell Labs report @lloyd1982least, Edward W. Forgy @forgy1965cluster introduced the canonical *batch alternating $k$-means* algorithm coupled with *uniform without-replacement data-point seeding*—now universally termed *Forgy initialization*.

The optimization target is the Euclidean within-cluster sum of squares (WCSS) functional over disjoint partitions $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ and prototype centroids $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$:
$ J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2. $
Unlike MacQueen's sequential variant (Alg. 21), Forgy holds all $k$ centroids *frozen* across a complete synchronous pass over $cal(X)$ and updates every centroid simultaneously once all $n$ assignments are finalized.

#definition(title: "Forgy Seeding Distribution & Synchronous Batch Operator")[
  The *Forgy seeding distribution* $bb(P)_"Forgy"$ draws a $k$-element subset $I = {I_1, dots, I_k} subset.eq {1, dots, n}$ uniformly without replacement from all $binom(n, k)$ choices and initializes $mu_j^((0)) = x_(I_j)$ for $j = 1, dots, k$. Given $bold(mu)^((t))$, the *synchronous batch operator* $cal(T)_"batch"$ maps $bold(mu)^((t)) |-> bold(mu)^((t+1))$ by forming Voronoi cells $C_j^((t)) = {i : j = arg min_ell norm(x_i - mu_ell^((t)))_2^2}$ (smallest-index tie-break) and setting $mu_j^((t+1)) = |C_j^((t))|^(-1) sum_(i in C_j^((t))) x_i$.
]

Algorithm 1 formalizes Forgy's complete batch procedure with $O(k)$ Fisher–Yates seeding and synchronous updates @forgy1965cluster.

#algbox(
  title: "Forgy's Batch K-Means Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, integer $k >= 1$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. *Forgy Seeding:* Sample $I = {I_1, dots, I_k} subset.eq {1, dots, n}$ uniformly without replacement; set $mu_j^((0)) <- x_(I_j)$ for $j = 1, dots, k$.
  2. *For* $t = 0, 1, dots, T - 1$ *do*:
  3. #h(1em) *Batch Assignment:* For all $i = 1, dots, n$ in parallel,
     $ c_i^((t)) <- arg min_(1 <= j <= k) norm(x_i - mu_j^((t)))_2^2. $
  4. #h(1em) Form $C_j^((t)) = {i : c_i^((t)) = j}$ for $j = 1, dots, k$.
  5. #h(1em) *Batch Centroid Update:* For each non-empty $C_j^((t))$,
     $ mu_j^((t+1)) <- 1 / (|C_j^((t))|) sum_(i in C_j^((t))) x_i. $
     (If any $C_j^((t)) = emptyset$ for $t >= 1$, reseed $mu_j^((t+1))$ to a farthest unassigned point.)
  6. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$, *break*.
  7. *Return* $(cal(C)^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Because Forgy's batch loop shares Lloyd's Voronoi–centroid fixed-point mechanics (Alg. 1), the defining mathematical theory of Forgy's K-Means centers on the *probabilistic geometry of uniform data-point seeding* and its contrast with Random-Partition seeding, MacQueen online updating, and K-Means++ $D^2$-seeding.

#subsec("A", "Exact Well-Separated Recovery Probability & Exponential Decay")

#theorem(title: "Exact Forgy Hit Probability and Exponential Decay", origin: "Original Proof")[
  Suppose $cal(X)$ partitions into $k$ well-separated ground-truth clusters $C_1^*, dots, C_k^*$ with cardinalities $n_j = |C_j^*|$ and proportions $p_j = n_j / n$ ($sum_(j=1)^k p_j = 1$). The exact probability $P_"hit"$ that Forgy seeding places exactly one initial centroid into every cluster $C_j^*$ is
  $ P_"hit" &= (product_(j=1)^k n_j) / binom(n, k) \
    &= (k! product_(j=1)^k n_j) / (n (n - 1) dots.c (n - k + 1)). $
  Consequently, for $n >= 2 k^2$:
  1. $P_"hit" <= (n / (n - k + 1))^k (k! / k^k) <= sqrt(2 pi e k) space e^(-k + 1 / (12 k)) = O(sqrt(k) e^(-k))$.
  2. If a minority cluster satisfies $p_"min" = epsilon << 1 / k$, then $P_"hit" <= O(k epsilon space e^(-k))$, collapsing linearly in $epsilon$.
]

#proof[
  Under uniform sampling without replacement, all $binom(n, k)$ unordered $k$-subsets of ${1, dots, n}$ are equally likely. The number of $k$-subsets containing exactly one element from each disjoint cluster $C_1^*, dots, C_k^*$ is $product_(j=1)^k n_j$ by the multiplication principle, yielding the exact formula for $P_"hit"$.

  To upper-bound $P_"hit"$, note $n (n - 1) dots.c (n - k + 1) >= (n - k + 1)^k$, hence
  $ P_"hit" <= (n / (n - k + 1))^k k! product_(j=1)^k p_j. $
  By the AM–GM inequality, $(product_(j=1)^k p_j)^(1/k) <= k^(-1) sum_(j=1)^k p_j = 1 / k$, so $product_(j=1)^k p_j <= k^(-k)$ with equality if and only if $p_1 = dots = p_k = 1 / k$. Applying Stirling's upper bound $k! <= sqrt(2 pi k) (k / e)^k e^(1 / (12 k))$ yields $k! / k^k <= sqrt(2 pi k) e^(-k + 1/(12 k))$. Finally, for $n >= 2 k^2$, $(n / (n - k + 1))^k = (1 + (k - 1)/(n - k + 1))^k <= exp(k (k - 1) / (n - k + 1)) <= e^(1/2)$, establishing Claim 1. For Claim 2, fix $p_1 = epsilon$ and apply AM–GM to the remaining $k - 1$ proportions summing to $1 - epsilon$:
  $ k! product_(j=1)^k p_j &<= k! space epsilon ((1 - epsilon) / (k - 1))^(k - 1) \
    &<= k epsilon space (k - 1)! / (k - 1)^(k - 1) = O(k epsilon space e^(-k)), $
  proving simultaneous exponential decay in $k$ and linear collapse in minority cluster mass $epsilon$.
]

#subsec("B", "Initial Dispersion: Forgy vs. Random-Partition Seeding")

A classical question in prototype clustering is why Forgy's data-point seeding outperforms *Random Partition* (RP) seeding—assigning each point uniformly to ${1, dots, k}$ and averaging—in iteration count and early symmetry breaking.

#theorem(title: "Exact Initial Centroid Spread: Forgy vs. Random Partition", origin: "Original Proof")[
  Let $overline(x) = n^(-1) sum_(i=1)^n x_i$ and $sigma_X^2 = n^(-1) sum_(i=1)^n norm(x_i - overline(x))_2^2$ denote the global empirical mean and total variance of $cal(X)$.
  1. Under Forgy seeding, each seed satisfies $bb(E)[norm(mu_j^((0)) - overline(x))_2^2] = sigma_X^2$, and every distinct seed pair ($j != ell$) satisfies
     $ bb(E)[norm(mu_j^((0)) - mu_ell^((0)))_2^2] = (2 n) / (n - 1) sigma_X^2. $
  2. Under balanced Random-Partition seeding (partitioning $cal(X)$ uniformly into $k$ groups of size $m = n / k$), each initial centroid $tilde(mu)_j^((0))$ satisfies
     $ bb(E)[norm(tilde(mu)_j^((0)) - overline(x))_2^2] = ((k - 1) / (n - 1)) sigma_X^2. $
  Thus Forgy's initial centroids exhibit an exact $(n - 1) / (k - 1) = Theta(n / k)$ factor larger expected squared radius from $overline(x)$ than Random-Partition centroids.
]

#proof[
  (1) Marginally, each $mu_j^((0)) = x_(I_j)$ is uniform over $cal(X)$, so $bb(E)[norm(mu_j^((0)) - overline(x))_2^2] = n^(-1) sum_(i=1)^n norm(x_i - overline(x))_2^2 = sigma_X^2$. For distinct $j != ell$, $(I_j, I_ell)$ is uniform over the $n (n - 1)$ ordered pairs $(i, r)$ with $i != r$. Using the identity $sum_(i=1)^n sum_(r=1)^n norm(x_i - x_r)_2^2 = 2 n sum_(i=1)^n norm(x_i - overline(x))_2^2 = 2 n^2 sigma_X^2$ and noting diagonal terms ($i = r$) vanish,
  $ bb(E)[norm(mu_j^((0)) - mu_ell^((0)))_2^2] &= (2 n^2 sigma_X^2) / (n (n - 1)) \
    &= (2 n) / (n - 1) sigma_X^2. $

  (2) Let $S_j subset.eq {1, dots, n}$ be a uniform random subset of size $m = n / k$ without replacement, with centroid $tilde(mu)_j^((0)) = m^(-1) sum_(i in S_j) x_i$. Define centered vectors $z_i = x_i - overline(x)$ (so $sum_(i=1)^n z_i = bold(0)$). Expanding the squared norm:
  $ m^2 bb(E)[norm(tilde(mu)_j^((0)) - overline(x))_2^2] &= bb(E)[sum_(i in S_j) norm(z_i)_2^2] \
    &quad + bb(E)[sum_(i != r in S_j) chevron.l z_i, z_r chevron.r]. $
  Every index $i$ belongs to $S_j$ with probability $m / n$, and every ordered distinct pair $(i, r)$ belongs to $S_j$ with probability $m (m - 1) / (n (n - 1))$. Because $sum_(i != r) chevron.l z_i, z_r chevron.r = norm(sum_i z_i)_2^2 - sum_i norm(z_i)_2^2 = - n sigma_X^2$, we obtain
  $ bb(E)[norm(tilde(mu)_j^((0)) - overline(x))_2^2] &= 1 / m^2 (m sigma_X^2 - (m (m - 1)) / (n - 1) sigma_X^2) \
    &= ((n - m) / (m (n - 1))) sigma_X^2. $
  Substituting $m = n / k$ yields $((k - 1) / (n - 1)) sigma_X^2$, proving that RP centroids collapse tightly around the unstable global mean saddle $overline(x)$ at rate $O(k / n)$, whereas Forgy seeds span the full data distribution.
]

#subsec("C", "Non-Empty Guarantee, Expected Initial Cost & Convergence")

#theorem(title: "Expected Initial Cost Bound & Batch Convergence", origin: "Original Proof")[
  Let $cal(X)$ contain at least $k$ distinct points and let $J_1 = sum_(i=1)^n norm(x_i - overline(x))_2^2 = n sigma_X^2$. Under Forgy's K-Means (Algorithm 1):
  1. *(Non-Empty Guarantee at $t = 0$):* Every initial Voronoi cluster satisfies $|C_j^((0))| >= 1$.
  2. *(Expected Initial Distortion):* The expected initial WCSS satisfies
     $ bb(E)_(I tilde bb(P)_"Forgy")[J(cal(C)^((0)), bold(mu)^((0)))] <= (2 (n - k)) / (n - 1) J_1. $
  3. *(Order-Invariant Finite Descent):* Given $bold(mu)^((0))$, every batch iteration is strictly invariant under dataset permutation $pi in S_n$, decreases $J$ by $sum_(j=1)^k |C_j^((t))| norm(mu_j^((t+1)) - mu_j^((t)))_2^2$ on each update, and terminates in at most $S(n, k)$ iterations at a discrete CVT.
]

#proof[
  (1) Because $I_1, dots, I_k$ index distinct points $x_(I_1), dots, x_(I_k)$, for each $j$ we have $norm(x_(I_j) - mu_j^((0)))_2^2 = 0 < norm(x_(I_j) - mu_ell^((0)))_2^2$ for all $ell != j$, ensuring $I_j in C_j^((0))$.

  (2) Every seed index $i in I$ incurs zero initial distortion, while every $i in.not I$ satisfies $min_j norm(x_i - x_(I_j))_2^2 <= norm(x_i - x_(I_1))_2^2$. Conditioning on $I_1$, the remaining $k - 1$ seeds are uniform among the $n - 1$ indices $i != I_1$, so $bb(P)(i in.not I | I_1) = (n - k) / (n - 1)$. Taking expectations,
  $ &bb(E)[J(cal(C)^((0)), bold(mu)^((0)))] \
    &quad <= bb(E)_(I_1)[sum_(i != I_1) bb(P)(i in.not I | I_1) norm(x_i - x_(I_1))_2^2] \
    &quad = ((n - k) / (n - 1)) 1 / n sum_(r=1)^n sum_(i=1)^n norm(x_i - x_r)_2^2 \
    &quad = (2 (n - k)) / (n - 1) J_1. $

  (3) Synchronous assignment and arithmetic averaging depend only on the multiset $cal(X)$ (independent of index order $pi$). Exact Lyapunov descent and finite termination in at most $S(n, k)$ steps follow from the parallel-axis identity identically to Lloyd's algorithm (Alg. 1 @lloyd1982least).
]

#subsec("D", "Unbounded Expected Approximation Ratio & Classical Bounds")

#theorem(title: "Unbounded Expected Approximation Ratio of Forgy's K-Means", origin: "Original Counterexample")[
  For any target factor $M > 0$, there exists a 4-point dataset $cal(X) subset.eq bb(R)$ with $k = 3$ on which the *expected* WCSS of Forgy's K-Means over its random initialization $bb(P)_"Forgy"$ satisfies
  $ bb(E)[J(cal(C)_"Forgy")] / J(cal(C)^*) > M. $
]

#proof[
  Set $k = 3$, $D = sqrt(2 M) + 1$, and $L = 2 D + 3$, and construct $cal(X) = {x_1, x_2, x_3, x_4} = {0, 1, L, L + D} subset.eq bb(R)$. The unique optimal partition $cal(C)^* = {{0, 1}, {L}, {L + D}}$ has cost $J(cal(C)^*) = (1/2)^2 + (1/2)^2 = 1 / 2$. Forgy sampling selects one of the $binom(4, 3) = 4$ seed triples uniformly with probability $1/4$ each:
  - *Triples ${0, 1, L}$ and ${0, 1, L + D}$:* Since $L - 1 > 2 D + 2 > D$, both $L$ and $L + D$ assign to the third seed, producing updated centroids $(0, 1, L + D / 2)$. Because $D / 2 < L - 1$, partition ${{0}, {1}, {L, L + D}}$ is a strict CVT fixed point with WCSS $0 + 0 + (D / 2)^2 + (D / 2)^2 = D^2 / 2$.
  - *Triples ${0, L, L + D}$ and ${1, L, L + D}$:* Points $0$ and $1$ both assign to the first seed (since $1 < L - 1$), converging immediately to the optimal partition ${{0, 1}, {L}, {L + D}}$ with WCSS $1 / 2$.

  Averaging over all four equally likely Forgy initializations yields
  $ bb(E)[J(cal(C)_"Forgy")] = 1 / 2 (D^2 / 2) + 1 / 2 (1 / 2) = (D^2 + 1) / 4. $
  Dividing by $J(cal(C)^*) = 1 / 2$ gives
  $ bb(E)[J(cal(C)_"Forgy")] / J(cal(C)^*) = (D^2 + 1) / 2 > (2 M + 1) / 2 > M. $
  (In contrast, K-Means++ @arthur2007kmeanspp selects the distant cluster with probability approaching $1$ and maintains an $O(log k)$ expected ratio uniformly over all $D, L > 0$.)
]

#proposition(title: "Classical Complexity, Smoothed Runtime & Breakdown", origin: "Classical — Cited")[
  1. Global WCSS minimization is NP-hard for $k = 2$ in general $d$ @aloise2009np and for $d = 2$ with arbitrary $k$ @mahajan2012planar; worst-case batch iterations can reach $2^(Omega(sqrt(n)))$ @vattani2011kmeans, whereas Gaussian-perturbed smoothed iteration count is $"poly"(n, k, d, 1 / sigma)$ @arthur2011smoothed.
  2. Empirical global minimizers are strongly consistent under finite second moment @pollard1981strong.
  3. Forgy's centroid breakdown point is $epsilon^* = 1 / n$, and any extreme outlier has probability $k / n$ of being chosen *directly* as an initial seed.
]

#proof[
  Claims (1)–(2) follow directly from @aloise2009np @mahajan2012planar @vattani2011kmeans @arthur2011smoothed @pollard1981strong because Forgy executes exact synchronous Lloyd iterations after seeding. For Claim (3), any corrupted point $norm(x_1)_2 -> infinity$ drives its cluster centroid to infinity (breakdown $1/n$), and under uniform without-replacement sampling of $k$ indices out of $n$, the marginal inclusion probability of index $1$ is exactly $binom(n - 1, k - 1) / binom(n, k) = k / n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Forgy's initialization requires only *$O(k)$* time via partial Fisher–Yates sampling—eliminating K-Means++'s $O(n k d)$ seeding scan—and each synchronous batch iteration costs $O(n k d)$ time and $O((n + k) d)$ memory. Because all $n$ assignments share a fixed centroid matrix $M in bb(R)^(k times d)$, squared distances expand as $norm(x_i - mu_j)_2^2 = norm(x_i)_2^2 + norm(mu_j)_2^2 - 2 x_i^T mu_j$, turning the assignment step into a pure BLAS-3 GEMM ($X M^T$) with maximal SIMD/GPU arithmetic intensity and zero serial data dependencies across points. When paired with $R$ independent multi-starts to mitigate Theorem 1's seeding decay, total complexity is $O(R n k d T)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Forgy's K-Means],
    [Time / Space Complexity], [$O(k)$ seeding; $O(n k d T)$ batch time; $O((n+k)d)$ space; full BLAS-3 parallelism.],
    [Seeding Recovery Rate], [$P_"hit" = O(sqrt(k) e^(-k))$ on balanced clusters; $O(k epsilon e^(-k))$ on minority mass $epsilon$ (Thm. 1).],
    [Approximation Guarantee], [Unbounded expected ratio $bb(E)[J_"Forgy"] / J^* -> infinity$ (Theorem 4).],
    [Initial Centroid Spread], [$Theta(n / k)$ larger initial variance than Random-Partition seeding (Theorem 2).],
    [Invariances / Breakdown], [Permutation-invariant (given seeds); isometry-equivariant; $epsilon^* = 1/n$ ($k/n$ seed hit).],
    [Recommended Upgrades], [K-Means++ (Alg. 2) for $O(log k)$ seeding; Elkan/Hamerly (Algs. 13–14) for exact distance pruning.],
  ),
  caption: [Theoretical guarantees and operational profile of Forgy's K-Means (Alg. 22).],
)

#sec("IV", "Engineering and Scientific Applications")

Forgy's K-Means remains a workhorse primitive whenever ultra-low initialization overhead, full batch vectorization, or deliberate stochastic ensemble diversity outweighs single-run worst-case guarantees.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Throughput GPU / SIMD Multi-Start Clustering:* When clustering small-to-moderate $k$ ($k <= 8$) across thousands of independent batched sub-problems on GPUs, Forgy's $O(k)$ seed draw avoids K-Means++'s sequential $k$-step kernel synchronization barrier, while batched GEMM exploits full hardware tensor throughput.
2. *Consensus Clustering and Co-Association Ensembles:* Consensus Clustering @monti2003consensus and cluster ensembles @strehl2002cluster intentionally seek diverse local equilibria across $R$ stochastic runs to construct stability matrices; Forgy's uniform data seeding naturally explores multiple CVT basins at minimal computational cost.
3. *Real-Time Color Palette Quantization & Preview Rendering:* Interactive image thumbnails and video-frame palette quantizers require sub-millisecond latency on dense, relatively balanced RGB/CIELAB clouds where Forgy's $O(k)$ seeding and fast convergence from dispersed seeds (Theorem 2) outperform heavier initializers.
4. *Balanced Customer and Operations Tiering:* When marketing or logistics cohorts are roughly balanced ($p_j approx 1 / k$) and $k$ is small ($k in {2, 3, 4}$), Forgy's hit probability ($k! / k^k = 0.50, 0.22, 0.09$) allows a handful of fast parallel restarts ($R approx 10$–$30$) to recover optimal convex partitions reliably.
5. *Baseline Vector Quantizer Codebook Initialization:* In classical LBG vector quantization @linde1980algorithm and large-scale data-clustering surveys @jain2010data, Forgy provides the canonical zero-overhead baseline against which adaptive seeding ($D^2$ sampling) and hierarchical splitting are benchmarked.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Forgy's K-Means fails predictably on (a) *large $k$ or imbalanced minority classes* (Theorem 1), where $P_"hit"$ vanishes exponentially; (b) *widely separated multi-scale clusters* (Theorem 4), where two seeds merge on one dense group while stranding distant clusters together; (c) *heavy outlier contamination*, where outliers are selected directly as initial centroids with probability $k / n$; and (d) *non-convex manifolds*, precluded by linear Voronoi boundaries.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit / Prerequisite], [Caution / Alternative],
    [Batched GPU Small-$k$], [$O(k)$ seeding + full BLAS-3 GEMM], [Run $R$ parallel multi-starts],
    [Consensus Ensembles], [Diverse CVT basin exploration], [Pair with stability scoring],
    [Imbalanced / Large $k$], [Avoid: $P_"hit" = O(k epsilon e^(-k))$ (Thm. 1)], [Use K-Means++ (Alg. 2)],
    [Outlier-Rich Finance], [Avoid: $k/n$ outlier seed probability], [Use PAM (Alg. 4) / Trimmed],
    [Non-Convex Manifolds], [Avoid: linear Voronoi hyperplanes], [Use Spectral / DBSCAN (Alg. 31)],
  ),
  caption: [Application suitability and engineering matrix for Forgy's K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Forgy's K-Means is the stochastic-seeding counterpart of Lloyd's deterministic batch iteration (Alg. 1), and its multi-start amplification law quantifies the exact computational price of uniform seeding relative to $D^2$ adaptive seeding.

#proposition(title: "Multi-Start Amplification Law & Exact Accelerator Equivalence", origin: "Original Derivation")[
  1. Executing $R$ independent Forgy restarts and selecting the lowest-WCSS output reduces the probability of failing to place one seed in every well-separated ground-truth cluster $C_1^*, dots, C_k^*$ to $(1 - P_"hit")^R <= exp(- R P_"hit")$. Thus $R >= ceil(P_"hit"^(-1) ln(1 / delta))$ restarts guarantee recovery with probability at least $1 - delta$—requiring $R = Omega(k^(-1/2) e^k ln(1 / delta))$ restarts on balanced clusters.
  2. Given any Forgy seed tuple $bold(mu)^((0))$, Elkan's (Alg. 13 @elkan2003using) and Hamerly's (Alg. 14 @hamerly2010making) triangle-inequality accelerators produce *bit-for-bit identical* iterate trajectories $(cal(C)^((t)), bold(mu)^((t)))$ to Algorithm 1.
]

#proof[
  (1) By independence across restarts, the joint miss probability is $(1 - P_"hit")^R$. Applying the elementary inequality $1 - u <= e^(-u)$ for $u in [0, 1]$ yields $(1 - P_"hit")^R <= exp(- R P_"hit") <= delta$ whenever $R >= P_"hit"^(-1) ln(1 / delta)$. Substituting Theorem 1's upper bound $P_"hit" = O(sqrt(k) e^(-k))$ proves the exponential restart lower bound $R = Omega(k^(-1/2) e^k ln(1 / delta))$.
  (2) Both Elkan @elkan2003using and Hamerly @hamerly2010making prune only cluster pairs provably excluded by the Euclidean triangle inequality, computing the exact same argmin assignments and synchronous means at every iteration $t >= 0$.
]

Furthermore, Forgy's batch update arises identically as the zero-temperature limit ($sigma -> 0^+$) of isotropic Gaussian Mixture EM (Alg. 81 @dempster1977em) initialized at $k$ random data means, and forms the $b = n$ synchronous endpoint contrasting with MacQueen's $b = 1$ online update (Alg. 21).

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark suite, Forgy's K-Means achieves *Iris* $"Sil" = 0.473$, $"ARI" = 0.431$; *Customers* $"Sil" = 0.330$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.238$ ($"ARI" = "N/A"$), with a mean runtime of *0.063 ms*. These results showcase every facet of our theory: Forgy is the fastest of the three variants ($2.2 times$ faster than MacQueen's $0.140$ ms and $16.4 times$ faster than COP-KMeans's $1.031$ ms), reflecting its $O(k)$ seed draw and fully vectorizable batch geometry. On *Multi-Shape* ($0.609$) and *High-D Subspace+Noise* ($0.238$), Forgy's dispersed data-point seeds (Theorem 2) and order-invariant batch means markedly outperform MacQueen's sequential drift ($0.357$ and $0.112$). Most tellingly, on *Iris*, Forgy yields the highest Euclidean Silhouette ($0.473$) yet the lowest ARI ($0.431$ vs. MacQueen's $0.586$ and COP-KMeans's $0.564$)—a direct empirical signature of Theorems 1 and 4, where single-draw uniform seeding ($k = 3, P_"hit" approx 22.2\%$) frequently lands in a compact, high-Silhouette local CVT that splits one species while merging two others.
