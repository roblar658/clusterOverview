#import "../template/ieee.typ": *

#alg_entry(
  id: 2,
  name: "K-Means++",
  category: "Cat. I: Centroid & Partitioning",
  year: "2007",
  paradigm: "Adaptive D²-Weighted Seeding",
  complexity: $O(n k d (1 + T))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset in Euclidean space and $k in {1, dots, n}$ the desired cluster count. For any non-empty center set $C subset.eq bb(R)^d$, denote the distance from $x in cal(X)$ to its nearest center by
$ D_C (x) := min_(c in C) norm(x - c)_2, $
and define the *within-cluster potential* of any subset $S subset.eq cal(X)$ relative to $C$ by
$ phi_C (S) := sum_(x in S) D_C (x)^2, quad phi(C) := phi_C (cal(X)). $
Let $cal(C)^* = {A_1, dots, A_k}$ denote a globally optimal $k$-means partition of $cal(X)$ with centroids $mu_r^* = |A_r|^(-1) sum_(x in A_r) x$ and optimal cost $phi_"OPT" := sum_(r=1)^k phi_"OPT" (A_r)$, where $phi_"OPT" (A_r) = sum_(x in A_r) norm(x - mu_r^*)_2^2$.

Whereas uniform Forgy initialization (Alg. 22) @forgy1965cluster suffers an unbounded worst-case approximation ratio even after Lloyd refinement (Alg. 1, Theorem 3), *K-Means++* @arthur2007kmeanspp replaces uniform seeding with an adaptive sequential sampling process that weights each candidate point proportionally to its current squared residual distance $D_C (x)^2$.

#definition(title: "D²-Weighted Adaptive Seeding Process")[
  Given $cal(X) subset.eq bb(R)^d$ and $k >= 1$, the *$D^2$-seeding Markov chain* constructs nested center sets $C_1 subset.eq C_2 subset.eq dots.c subset.eq C_k subset.eq cal(X)$ as follows:
  1. Draw $c_1 ~ "Unif"(cal(X))$ and set $C_1 = {c_1}$.
  2. For $t = 2, dots, k$, if $phi(C_(t-1)) > 0$, sample $c_t in cal(X)$ from the adaptive probability mass function
     $ p_(C_(t-1)) (x) := (D_(C_(t-1)) (x)^2) / (phi(C_(t-1))), quad x in cal(X), $
     and set $C_t = C_(t-1) union {c_t}$ (if $phi(C_(t-1)) = 0$, draw uniformly from $cal(X) backslash C_(t-1)$).
  The full *K-Means++* procedure outputs $C_k$ as the initial centroids $bold(mu)^((0))$ for Lloyd's algorithm (Alg. 1) @lloyd1982least.
]

Intuitively, any optimal cluster $A_r in cal(C)^*$ that has not yet been hit by $C_(t-1)$ contributes disproportionately large squared distances to $phi(C_(t-1))$, concentrating the $D^2$ probability mass onto uncovered regions of $cal(X)$. Algorithm 1 implements $D^2$-seeding in $O(n k d)$ time by maintaining a running distance cache $d_i = D_(C_t) (x_i)^2$.

#algbox(
  title: "K-Means++ Seeding and Lloyd Refinement",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, cluster count $k >= 1$, max Lloyd iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. Choose $c_1 in cal(X)$ uniformly at random; set $C <- {c_1}$.
  2. For each $i in {1, dots, n}$, cache $d_i <- norm(x_i - c_1)_2^2$.
  3. *For* $t = 2, dots, k$ *do*:
  4. #h(1em) Compute total potential $Phi <- sum_(i=1)^n d_i$.
  5. #h(1em) If $Phi = 0$, pick $c_t in cal(X) backslash C$ uniformly; else sample $c_t = x_i$ with probability $d_i / Phi$.
  6. #h(1em) Update $C <- C union {c_t}$.
  7. #h(1em) For each $i = 1, dots, n$, update $d_i <- min(d_i, norm(x_i - c_t)_2^2)$.
  8. Initialize $bold(mu)^((0)) <- C$ and execute up to $T$ iterations of Lloyd's K-Means (Alg. 1) until convergence.
  9. *Return* $(cal(C), bold(mu))$.
]

In production libraries, step 5 is often augmented with *greedy multi-trial $D^2$-seeding*, drawing $L = floor(2 + ln k)$ independent candidates from $p_(C_(t-1))$ and retaining the candidate minimizing the resulting potential $phi(C_(t-1) union {c})$.

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results characterizing K-Means++: (1) an exact factor-2 expectation identity for the uniform first seed; (2) the factor-8 potential bound for any $D^2$-covered optimal cluster; (3) the full $O(log k)$ expected approximation theorem via backward potential induction; (4) exponential high-probability boosting under restarts together with similarity-group invariance; and (5) an explicit adversarial counterexample exposing outlier-induced center starvation.

#subsec("A", "Single-Cluster Uniform and D²-Weighted Potential Bounds")

#lemma(title: "Exact Factor-2 Identity for Uniform First-Center Seeding", origin: "Original Proof")[
  Let $A subset.eq cal(X)$ be any finite non-empty cluster with centroid $overline(x)_A = |A|^(-1) sum_(x in A) x$ and optimal inertia $phi_"OPT" (A) = sum_(x in A) norm(x - overline(x)_A)_2^2$. If $c$ is drawn uniformly from $A$, then
  $ bb(E)_(c ~ "Unif"(A)) [phi_({c}) (A)] = 2 phi_"OPT" (A). $
]

#proof[
  For any fixed $c in A$, expanding $norm(x - c)_2^2 = norm((x - overline(x)_A) + (overline(x)_A - c))_2^2$ and summing over $x in A$ cancels the linear cross-term because $sum_(x in A) (x - overline(x)_A) = bold(0)$:
  $ phi_({c}) (A) &= sum_(x in A) norm(x - c)_2^2 \
                  &= phi_"OPT" (A) + |A| norm(c - overline(x)_A)_2^2. $
  Taking the expectation over $c ~ "Unif"(A)$ yields
  $ bb(E)[phi_({c}) (A)] &= phi_"OPT" (A) + sum_(c in A) norm(c - overline(x)_A)_2^2 \
                         &= phi_"OPT" (A) + phi_"OPT" (A) = 2 phi_"OPT" (A). $
  Notice that the pairwise sum satisfies the useful corollary identity $sum_(a_0 in A) sum_(a in A) norm(a_0 - a)_2^2 = 2 |A| phi_"OPT" (A)$.
]

#lemma(title: "Factor-8 Bound for a D²-Sampled Optimal Cluster", origin: "Classical — Cited")[
  Let $A in cal(C)^*$ be an arbitrary optimal cluster and $C subset.eq cal(X)$ any existing center set with $phi_C (A) > 0$. Conditioned on the next $D^2$-sampled center $c$ landing in $A$, the expected potential of $A$ after adding $c$ satisfies @arthur2007kmeanspp
  $ bb(E)[phi_(C union {c}) (A) | c in A] <= 8 phi_"OPT" (A). $
]

#proof[
  Conditioned on $c in A$, each point $a_0 in A$ is chosen with conditional probability $D_C (a_0)^2 / phi_C (A)$, and the resulting cluster potential is $sum_(a in A) min(D_C (a)^2, norm(a - a_0)_2^2)$. Thus
  $ &bb(E)[phi_(C union {c}) (A) | c in A] \
    &= sum_(a_0 in A) (D_C (a_0)^2) / (phi_C (A)) sum_(a in A) min(D_C (a)^2, norm(a - a_0)_2^2). $
  For any pair $a_0, a in A$, the triangle inequality $D_C (a_0) <= D_C (a) + norm(a_0 - a)_2$ combined with $(u + v)^2 <= 2 u^2 + 2 v^2$ yields
  $ D_C (a_0)^2 <= 2 D_C (a)^2 + 2 norm(a_0 - a)_2^2. $
  Averaging this inequality over all $a in A$ gives
  $ D_C (a_0)^2 &<= 2 / (|A|) sum_(a in A) D_C (a)^2 \
                &quad + 2 / (|A|) sum_(a in A) norm(a_0 - a)_2^2. $
  Substituting this upper bound for $D_C (a_0)^2$ into the conditional expectation splits it into two sums $T_1 + T_2$. In $T_1$, we bound $min(D_C (a)^2, norm(a - a_0)_2^2) <= norm(a - a_0)_2^2$; in $T_2$, we bound $min(D_C (a)^2, norm(a - a_0)_2^2) <= D_C (a)^2$. Using $sum_(a in A) D_C (a)^2 = phi_C (A)$ and the pairwise identity from Lemma 1:
  $ T_1 &<= 2 / (|A| phi_C (A)) phi_C (A) sum_(a_0, a' in A) norm(a' - a_0)_2^2 \
        &= 2 / (|A|) dot 2 |A| phi_"OPT" (A) = 4 phi_"OPT" (A), $
  and symmetrically for the second term,
  $ T_2 &<= 2 / (|A| phi_C (A)) (sum_(a_0, a in A) norm(a_0 - a)_2^2) phi_C (A) \
        &= 2 / (|A|) dot 2 |A| phi_"OPT" (A) = 4 phi_"OPT" (A). $
  Summing $T_1 + T_2 <= 8 phi_"OPT" (A)$ completes the proof.
]

#subsec("B", "The O(log k) Expected Approximation Guarantee")

#theorem(title: "Expected 8(ln k + 2)-Approximation of K-Means++", origin: "Classical — Cited")[
  Let $C_k$ be the center set returned by $D^2$-seeding and $C_"final"$ the centroids after any subsequent iterations of Lloyd's algorithm. Then @arthur2007kmeanspp
  $ bb(E)[phi(C_"final")] <= bb(E)[phi(C_k)] <= 8 (ln k + 2) phi_"OPT". $
]

#proof[
  Let $H_t := sum_(j=1)^t 1/j$ denote the $t$-th harmonic number ($H_0 = 0$). Consider any fixed partial center set $C subset.eq cal(X)$ that intersects $k - u$ optimal clusters (whose union is $X_c$, the *covered* points) and misses $u >= 1$ optimal clusters (whose union is $X_u$, the *uncovered* points). Let $C'$ be obtained by appending $t <= u$ additional centers via $D^2$-seeding. We prove by induction on $(u, t)$ that
  $ bb(E)[phi(C')] &<= (phi_C (X_c) + 8 phi_"OPT" (X_u))(1 + H_t) \
                   &quad + (u - t) / u phi_C (X_u). $
  *Base cases:* For $t = 0$, $H_0 = 0$, so the right-hand side equals $phi_C (X_c) + 8 phi_"OPT" (X_u) + phi_C (X_u) >= phi(C)$. For $t = 1$ and $u = 1$, either the drawn center $c$ lands in $X_c$ (probability $p_c = phi_C (X_c) / phi(C)$, leaving potential $<= phi(C)$) or in $X_u$ (probability $p_u = phi_C (X_u) / phi(C)$, yielding expected potential $<= phi_C (X_c) + 8 phi_"OPT" (X_u)$ by Lemma 2). Thus
  $ bb(E)[phi(C')] &<= p_c (phi_C (X_c) + phi_C (X_u)) \
                   &quad + p_u (phi_C (X_c) + 8 phi_"OPT" (X_u)) \
                   &= phi_C (X_c) + p_u phi_C (X_c) + 8 p_u phi_"OPT" (X_u) \
                   &<= 2 (phi_C (X_c) + 8 phi_"OPT" (X_u)), $
  matching $1 + H_1 = 2$ and $(u - t)/u = 0$.

  *Inductive step:* Assume the claim holds for $(u, t - 1)$ and $(u - 1, t - 1)$ with $1 <= t <= u$. Condition on the first of the $t$ centers, $c$:
  - *Case 1 ($c in X_c$, prob. $phi_C (X_c) / phi(C)$):* Adding $c$ cannot increase any point's distance, and $(u, t - 1)$ remain uncovered/unselected. By the induction hypothesis,
    $ bb(E)[phi(C') | c in X_c] &<= (phi_C (X_c) + 8 phi_"OPT" (X_u))(1 + H_(t-1)) \
                                &quad + (u - t + 1) / u phi_C (X_u). $
  - *Case 2 ($c in A subset.eq X_u$, prob. $phi_C (A) / phi(C)$):* Cluster $A$ becomes covered, leaving $u - 1$ uncovered clusters and $t - 1$ draws. Combining the induction hypothesis on $(u - 1, t - 1)$ with Lemma 2 ($bb(E)[phi_(C union {c}) (A) | c in A] <= 8 phi_"OPT" (A)$) yields
    $ bb(E)[phi(C') | c in A] &<= (phi_C (X_c) + 8 phi_"OPT" (X_u))(1 + H_(t-1)) \
                              &quad + (u - t) / (u - 1) (phi_C (X_u) - phi_C (A)). $

  Weighting Case 2 across all $u$ uncovered clusters $A subset.eq X_u$ and invoking the Cauchy–Schwarz inequality $sum_(A subset.eq X_u) phi_C (A)^2 >= u^(-1) phi_C (X_u)^2$ gives
  $ &sum_(A subset.eq X_u) (phi_C (A)) / (phi(C)) (phi_C (X_u) - phi_C (A)) \
    &quad <= (u - 1) / u (phi_C (X_u)^2) / (phi(C)). $
  Multiplying by $(u - t)/(u - 1)$ and combining with Case 1 yields
  $ bb(E)[phi(C')] &<= (phi_C (X_c) + 8 phi_"OPT" (X_u))(1 + H_(t-1)) \
                   &quad + (u - t) / u phi_C (X_u) (phi_C (X_c) + phi_C (X_u)) / (phi(C)) \
                   &quad + 1 / u (phi_C (X_c) phi_C (X_u)) / (phi(C)). $
  Since $phi_C (X_c) + phi_C (X_u) = phi(C)$, $phi_C (X_u) <= phi(C)$, and $1 / u <= 1 / t$, the final term is bounded by $t^(-1) phi_C (X_c) <= t^(-1) (phi_C (X_c) + 8 phi_"OPT" (X_u))$, completing the induction because $1 + H_(t-1) + 1/t = 1 + H_t$.

  Finally, apply the bound immediately after the uniform first center $c_1$ lands in some optimal cluster $A_1$, where $u = t = k - 1$ and $X_c = A_1$. By Lemma 1, $bb(E)[phi_({c_1}) (A_1) | c_1 in A_1] = 2 phi_"OPT" (A_1) <= 8 phi_"OPT" (A_1)$. Taking total expectation over $A_1$ yields
  $ bb(E)[phi(C_k)] <= 8 (1 + H_(k-1)) phi_"OPT" <= 8 (ln k + 2) phi_"OPT", $
  using $H_(k-1) <= 1 + ln k$ for $k >= 1$. Because Lloyd iterations monotonically non-increase $phi$ (Alg. 1, Theorem 1), $bb(E)[phi(C_"final")] <= bb(E)[phi(C_k)]$.
]

#remark(title: "Tightness of the Logarithmic Factor")[
  Arthur and Vassilvitskii @arthur2007kmeanspp construct a $(k - 1)$-dimensional regular simplex of $k$ equal-size clusters separated by $Delta >> 1$ on which $bb(E)[phi(C_k)] >= 2 ln k dot phi_"OPT"$ asymptotically, proving that the $O(log k)$ factor in Theorem 1 is tight up to constants.
]

#subsec("C", "High-Probability Amplification and Similarity Invariance")

#proposition(title: "Exponential Tail Boosting and Similarity Invariance", origin: "Original Proof")[
  1. Let $C^((1)), dots, C^((R))$ be $R = ceil(log_2 (1 / delta))$ independent $D^2$-seedings and $C_"min" = arg min_(1 <= r <= R) phi(C^((r)))$. Then with probability at least $1 - delta$,
     $ phi(C_"min") <= 16 (ln k + 2) phi_"OPT". $
  2. The $D^2$-seeding distribution is strictly invariant under all Euclidean similarity maps $T(x) = alpha Q x + b$ ($alpha > 0$, $Q in O(d)$, $b in bb(R)^d$), yet is not invariant under anisotropic diagonal scaling $S(x) = "diag"(s_1, dots, s_d) x$ ($s_p != s_q > 0$).
]

#proof[
  (1) Because $phi(C^((r))) >= 0$, Markov's inequality applied to Theorem 1 yields
  $ &bb(P)[phi(C^((r))) > 16 (ln k + 2) phi_"OPT"] \
    &quad <= (bb(E)[phi(C^((r)))]) / (16 (ln k + 2) phi_"OPT") <= 1 / 2. $
  By independence of the $R$ trials, the probability that every trial exceeds the threshold is at most $2^(-R) <= delta$, giving the high-probability guarantee for $C_"min"$.

  (2) Under any similarity $T(x) = alpha Q x + b$ with $Q^T Q = I_d$, squared distances scale uniformly by $alpha^2$: $D_(T(C)) (T(x))^2 = alpha^2 D_C (x)^2$, and $phi(T(C)) = alpha^2 phi(C)$. Hence the $D^2$ transition kernel satisfies
  $ p_(T(C)) (T(x)) = (alpha^2 D_C (x)^2) / (alpha^2 phi(C)) = p_C (x) $
  identically at every step $t = 1, dots, k$. Conversely, on $cal(X) = {(0,0), (1,0), (0,10)}$ with $C_1 = {(0,0)}$, $S = "diag"(100, 1)$ reverses the sampling ratio from $1:100$ to $100:1$, altering the seeding law.
]

#subsec("D", "Adversarial Outlier Capture and Center Starvation")

#theorem(title: "Outlier-Driven Center Starvation in D²-Seeding", origin: "Original Counterexample")[
  For every $k >= 3$, every cluster size $m >= 2$, and every target distortion $M > 0$, there exists a dataset $cal(X) = cal(X)_"in" union {x_"out"} subset.eq bb(R)$ with $n = k m$ well-separated inliers and one single outlier $x_"out"$ such that with probability at least $1 - 2 / n$, $D^2$-seeding selects $x_"out" in C_k$ and suffers inlier clustering cost exceeding $M$ times the optimal $k$-cluster inlier cost.
]

#proof[
  Place $k$ inlier clusters $A_1, dots, A_k$ of $m$ points each along $bb(R)$: for $r = 1, dots, k$, put $floor(m / 2)$ points at $(r - 1) Delta - 1$ and $floor(m / 2)$ points at $(r - 1) Delta + 1$ (if $m$ is odd, place the remaining point at $(r - 1) Delta$), with inter-cluster separation $Delta > 4$. Then $phi_"OPT" (cal(X)_"in") <= k m$. Append one outlier $x_"out" = R$ with $R > (n + 1) k Delta$.

  At step $t = 1$, $c_1 in cal(X)_"in"$ occurs with probability $n / (n + 1) > 1 - 1 / n$. Conditioned on $c_1 in cal(X)_"in"$, every inlier satisfies $D_(C_1) (x) <= (k - 1) Delta + 2 < k Delta$, so $phi_(C_1) (cal(X)_"in") < n k^2 Delta^2$, whereas $D_(C_1) (x_"out") > R - k Delta > n k Delta$. Therefore at step $t = 2$,
  $ &bb(P)(c_2 = x_"out" | c_1 in cal(X)_"in") \
    &quad = (D_(C_1) (x_"out")^2) / (D_(C_1) (x_"out")^2 + phi_(C_1) (cal(X)_"in")) \
    &quad >= (n^2 k^2 Delta^2) / (n^2 k^2 Delta^2 + n k^2 Delta^2) = n / (n + 1). $
  By the union bound, $bb(P)(x_"out" in C_2) >= 1 - 2 / (n + 1) > 1 - 2 / n$. Once $x_"out" in C_k$, at most $k - 1$ centers remain in $C_k$ to cover the $k$ distinct clusters $A_1, dots, A_k$. By the pigeonhole principle, at least one pair of distinct clusters separated by at least $Delta - 2$ shares a single Voronoi centroid, forcing inlier inertia at least $m ((Delta - 2)/2)^2$. Choosing $Delta > 2 + 2 sqrt(k M)$ yields inlier cost $> M k m >= M phi_"OPT" (cal(X)_"in")$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Caching $d_i = D_(C_t) (x_i)^2$ reduces $D^2$-seeding from naive $O(n k^2 d)$ to $O(n k d)$—identical in asymptotic cost to a single iteration of Lloyd's algorithm—while requiring only $O((n + k) d)$ memory. Because each draw $c_t$ depends strictly on the updated distances from $c_(t-1)$, classic K-Means++ imposes $k$ sequential passes over $cal(X)$, creating an I/O synchronization barrier when $k$ is large ($k >= 10^3$) or data are distributed across cluster nodes. This bottleneck directly motivated parallel *K-Means||* (Alg. 23) @bahmani2012scalable, which oversamples $ell = Theta(k)$ centers per pass to achieve an $O(1)$-approximation in $O(log phi_1)$ rounds.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Means++],
    [Time / Space Complexity], [$O(n k d)$ seeding + $O(n k d T)$ Lloyd; $O((n+k)d)$ space; $k$ sequential passes.],
    [Approximation Guarantee], [$bb(E)[phi] <= 8(ln k + 2) phi_"OPT"$ (tight $Omega(log k)$); $1 - delta$ via $O(log(1/delta))$ restarts.],
    [Determinism / Variance], [Stochastic seeding; dramatically lower variance & faster convergence than Forgy (Alg. 22).],
    [Invariances], [Invariant under similarities $alpha Q x + b$; sensitive to anisotropic feature scaling.],
    [Outlier Breakdown Point], [$epsilon^* = 1/n$; squared weights actively attract seeds to extreme outliers (Theorem 2).],
    [Recommended Upgrades], [K-Means|| (Alg. 23) for large $k$ / MapReduce; Elkan/Hamerly (Algs. 13–14) for Lloyd; PAM (Alg. 4) for outliers.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of K-Means++.],
)

#sec("IV", "Engineering and Scientific Applications")

K-Means++ is the default initialization standard across virtually all scientific computing ecosystems (`scikit-learn`, `FAISS`, `Spark MLlib`), converting Lloyd's heuristic into a provably competitive approximation algorithm at negligible computational overhead.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Vector Database Coarse Quantizers (IVF / ScaNN):* Partitioning $10^7$–$10^9$ neural embeddings into thousands of Voronoi buckets demands uniform spatial coverage across semantic modes. $D^2$-seeding prevents redundant centroids from collapsing onto high-density hubs, minimizing worst-case inverted-list imbalance.
2. *Initialization Backbone for Mixture and Spectral Pipelines:* Gaussian Mixture Models (Alg. 81), Spherical K-Means (Alg. 12), and Spectral Clustering (Algs. 103–105) all exhibit non-convex likelihood or trace landscapes; K-Means++ seeding provably places initial prototypes inside well-separated modes before EM or k-means refinement.
3. *Color Palette & Image Superpixel Seeding:* In lossy image compression and visual codebook learning, rare but visually salient chromatic regions contribute large squared RGB/CIELAB residuals; $D^2$-sampling allocates palette entries to salient minority colors that uniform seeding misses.
4. *Coreset Construction for Streaming Analytics:* Sensitivity-weighted sampling derived from $D^2$-seeding forms the theoretical core of StreamKM++ (Alg. 195) @ackermann2012streamkmpp and coreset compression @harpeled2004coresets, preserving $k$-means costs within $(1 plus.minus epsilon)$ on compact weighted sketches.
5. *Industrial Telemetry and Smart-Grid Load Profiling:* Daily electricity consumption curves exhibit distinct seasonal and peak-tariff regimes; $D^2$-seeding reliably discovers rare peak-demand profiles prior to Euclidean prototype refinement.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-Means++ fails under (a) *extreme uncleaned outliers*, where $D^2$-weighting actively starves genuine clusters by placing seeds directly on remote anomalies (Theorem 2; mitigate via Trimmed K-Means Alg. 26 or PAM Alg. 4); (b) *massive $k$ in distributed storage*, where $k$ sequential passes dominate wall-clock time (mitigate via K-Means|| Alg. 23); (c) *non-convex manifold geometries* inherited from Euclidean Voronoi cells (mitigate via DBSCAN Alg. 31 or Spectral Clustering Algs. 103–105); and (d) *unscaled heterogeneous features* (Proposition 1(2)).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Means++ Fits / Prerequisite], [Caution / Alternative],
    [ANN Codebooks], [Prevents mode collapse across $k$ buckets], [Use K-Means|| (Alg. 23) for huge $k$],
    [GMM / Spectral Init], [Provable $O(log k)$ warm-start basin], [Standardize axes first (Prop. 1)],
    [Color / Image Quant.], [Captures rare high-contrast chromatic modes], [Use LBG (Alg. 30) for deterministic tree],
    [Anomaly-Heavy Logs], [Fast coarse partition after clipping], [Outliers steal seeds; use PAM (Alg. 4)],
    [Streaming Sketches], [Exact sensitivity basis for coresets], [Use StreamKM++ (Alg. 195) / BICO (Alg. 78)],
  ),
  caption: [Application suitability and engineering matrix for K-Means++.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Beyond initialization, the $D^2$ distribution is mathematically equivalent to an optimal *importance-sampling proposal* for estimating $k$-means objectives over compressed subsets:

#proposition(title: "Unbiased Importance-Weighted Cost Estimation via D²-Sampling", origin: "Original Derivation")[
  Fix a pilot center set $C$ with $phi(C) > 0$, and let $S = {y_1, dots, y_m}$ consist of $m$ i.i.d. draws from the $D^2$ distribution $p_C (x) = D_C (x)^2 / phi(C)$ equipped with inverse-probability weights $w(y_r) = phi(C) / (m D_C (y_r)^2)$. Then for every candidate center set $Q subset.eq bb(R)^d$, the weighted coreset estimator
  $ hat(phi)_S (Q) := sum_(r=1)^m w(y_r) D_Q (y_r)^2 $
  is strictly unbiased on $cal(X)_+ = {x in cal(X) : D_C (x) > 0}$, i.e., $bb(E)[hat(phi)_S (Q)] = phi_Q (cal(X)_+)$, with variance bounded by $m^(-1) phi(C) sum_(x in cal(X)_+) D_Q (x)^4 / D_C (x)^2$.
]

#proof[
  By linearity of expectation and i.i.d. draws $y_r ~ p_C$:
  $ bb(E)[hat(phi)_S (Q)] &= m sum_(x in cal(X)_+) (D_C (x)^2) / (phi(C)) dot (phi(C)) / (m D_C (x)^2) D_Q (x)^2 \
                          &= sum_(x in cal(X)_+) D_Q (x)^2 = phi_Q (cal(X)_+). $
  Similarly, for a single draw $Z = w(y_1) D_Q (y_1)^2$, dropping the non-negative squared mean yields
  $ "Var"(hat(phi)_S (Q)) &<= 1 / m bb(E)[(phi(C) D_Q (y_1)^2 / D_C (y_1)^2)^2] \
                          &= (phi(C)) / m sum_(x in cal(X)_+) (D_Q (x)^4) / (D_C (x)^2), $
  demonstrating why $D^2$-sampling suppresses estimator variance whenever $D_Q (x) = O(D_C (x))$ across distant points.
]

This importance-sampling bridge directly unifies K-Means++ with coresets @harpeled2004coresets, StreamKM++ (Alg. 195) @ackermann2012streamkmpp, BICO (Alg. 78) @fichtenberger2013bico, and parallel K-Means|| (Alg. 23) @bahmani2012scalable, while exact geometric accelerators such as Elkan (Alg. 13) @elkan2003using and Hamerly (Alg. 14) @hamerly2010making accelerate the subsequent Lloyd phase without altering a single iterate.

On the ClusterOverview 4-dataset empirical benchmark, K-Means++ achieves *Iris Silhouette = 0.447* and *Iris ARI = 0.583*, *Customers Silhouette = 0.334*, *Multi-Shape Silhouette = 0.712* (ARI N/A), and *High-D Subspace+Noise Silhouette = 0.303* (ARI N/A), at a remarkably fast *mean runtime of 0.229 ms*. Its sub-millisecond runtime reflects rapid Lloyd convergence from near-optimal initial seeds, whereas its sharp drop to Silhouette = 0.303 on the noisy high-dimensional benchmark empirically mirrors the outlier-capture and squared-distance concentration vulnerability established in Theorem 2.
