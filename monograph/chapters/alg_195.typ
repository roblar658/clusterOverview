#import "../template/ieee.typ": *

#alg_entry(
  id: 195,
  name: "StreamKM++",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2012",
  paradigm: "Coreset-Tree Merge-and-Reduce Streaming k-Means++",
  complexity: $O(n m d + I_"off" m k d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

StreamKM++ @ackermann2012streamkmpp solves the Euclidean $k$-Means problem over a massive single-pass data stream $P = (x_1, x_2, dots, x_n) subset RR^d$ with strict sublinear memory and provable multiplicative approximation guarantees. For any weighted point set $(Q, w)$ with $w : Q -> RR_(>0)$ and any $k$-center set $C subset RR^d$ ($abs(C) = k$), define the weighted $k$-Means objective
$
  Phi(Q, C) = sum_(x in Q) w(x) min_(c in C) norm(x - c)_2^2 ,
$
and let $"OPT"_k (Q) = min_(abs(C) = k) Phi(Q, C)$. Because exact $k$-Means is NP-hard @aloise2009np and storing $P$ is infeasible when $n$ is huge, StreamKM++ reduces $P$ online to a small weighted *$(k, epsilon)$-coreset* $(S, w_S)$ of size $abs(S) = m << n$ whose objective approximates $Phi(P, C)$ simultaneously for *every* candidate solution $C$.

#definition(title: "Strong (k, epsilon)-Coreset for Euclidean k-Means")[
  Let $epsilon in (0, 1)$. A weighted point set $(S, w_S)$ is a *strong $(k, epsilon)$-coreset* of $(P, w_P)$ for the $k$-Means problem @harpeled2004coresets if for every $k$-element subset $C subset RR^d$,
  $
    (1 - epsilon) Phi(P, C) <= Phi(S, C) <= (1 + epsilon) Phi(P, C) .
  $
]

StreamKM++ achieves single-pass coreset construction by combining two mechanisms:
1. *The Coreset Tree (`Reduce`):* Given a weighted set $Q$ of size $2m$, `Reduce` compresses $Q$ to size $m$ by constructing a binary partition tree $T$ with $m$ leaves via $m - 1$ randomized splits. Starting from a single root leaf $v_0$ holding $Q$ with a randomly drawn representative $q(v_0) in Q$, each split traverses from root to a leaf $v$ by moving at every internal node $u$ to left child $u_L$ with probability $"SSE"(u_L) / ("SSE"(u_L) + "SSE"(u_R))$ and to right child $u_R$ otherwise, where $"SSE"(v) = sum_(x in Q(v)) w(x) norm(x - q(v))_2^2$. At leaf $v$, a new center $q_2 in Q(v)$ is drawn via $D^2$-sampling @arthur2007kmeanspp with probability proportional to $w(x) norm(x - q(v))_2^2$. Then $Q(v)$ is partitioned into $Q(v_L)$ (points closer to $q(v)$) and $Q(v_R)$ (points closer to $q_2$), and ancestor costs are updated bottom-up. After $m - 1$ splits, each leaf $v$ emits its weighted centroid $c(v) = (1 / W(v)) sum_(x in Q(v)) w(x) x$ with total weight $W(v) = sum_(x in Q(v)) w(x)$.
2. *Merge-and-Reduce Streaming Architecture:* Inspired by Bentley--Saxe dynamization and Har-Peled--Mazumdar @harpeled2004coresets, StreamKM++ maintains a logarithmic array of buckets $B_0, B_1, dots, B_L$ ($L = ceil(log_2 (n / m)) + 1$), each holding either $emptyset$ or a weighted coreset of size $m$. Incoming points fill $B_0$; whenever $abs(B_0) = m$, the bucket cascades upward like binary carry propagation—merging with occupied buckets $B_i$ and applying `Reduce` to restore size $m$ until reaching the first empty slot. At query time, all nonempty buckets are merged and reduced to a single coreset $S$ of size $m$, on which weighted $k$-Means++ @arthur2007kmeanspp and Lloyd refinement @lloyd1982least produce the final $k$ centers.

#algbox(
  title: "StreamKM++ (Merge-and-Reduce with Coreset Tree)",
  inputs: [Stream $P = (x_1, dots, x_n)$, cluster count $k$, coreset size $m >= k$.],
  outputs: [Final $k$-Means center set $hat(C) subset RR^d$ of size $k$.],
)[
  1. *Init:* Initialize empty bucket array $B_0, B_1, dots <- emptyset$ (allocated dynamically as needed).
  2. *Online Stream Ingestion:* For each arriving point $x in P$:
  3. #h(1em) Insert $(x, 1)$ into $B_0$. If $abs(B_0) < m$, continue.
  4. #h(1em) Set $Q <- B_0$, $B_0 <- emptyset$, and $i <- 1$.
  5. #h(1em) *While* $B_i != emptyset$: set $Q <- "Reduce"(Q union B_i, m)$ via the Coreset Tree ($m - 1$ randomized $D^2$-splits), clear $B_i <- emptyset$, and set $i <- i + 1$.
  6. #h(1em) Store $B_i <- Q$.
  7. *Offline Extraction:* Collect all nonempty buckets $B_i$ ($i >= 0$), iteratively apply `"Reduce"` to obtain a single weighted coreset $S$ of size $m$.
  8. Run weighted $k$-Means++ seeding followed by weighted Lloyd iterations on $S$; *return* $hat(C)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Distributional Equivalence of Coreset-Tree Sampling")

A central algorithmic contribution of StreamKM++ @ackermann2012streamkmpp is avoiding full $O(m^2 d)$ global distance updates per split while preserving exact $D^2$-sampling probabilities relative to the tree partition.

#lemma(title: "Telescoping Equivalence of Tree Walk and Global Partition D2-Sampling", origin: "Original Proof")[
  At step $t in {1, dots, m - 1}$ of `"Reduce"(Q, m)` with current leaf set $cal(L)_t$ and maintained internal sums $"SSE"(u) = "SSE"(u_L) + "SSE"(u_R)$:
  1. The root-to-leaf random walk followed by intra-leaf $D^2$-sampling selects point $x in Q$ with *exact* probability
  $
    "Pr"("select " x) = (w(x) norm(x - q(v(x)))_2^2) / (sum_(y in Q) w(y) norm(y - q(v(y)))_2^2) ,
  $
  where $v(x) in cal(L)_t$ is the unique leaf containing $x$.
  2. Emitting each leaf's weighted centroid $c(v)$ instead of its sampled representative $q(v)$ strictly reduces quantization error by the exact nonnegative deficit
  $
    sum_(v in cal(L)_m) sum_(x in Q(v)) w(x) norm(x - c(v))_2^2 = "SSE"(u_0) - sum_(v in cal(L)_m) W(v) norm(q(v) - c(v))_2^2 .
  $
]

#proof[
  Let $u_0 -> u_1 -> dots -> u_h = v$ be the unique path from root $u_0$ to leaf $v in cal(L)_t$. Assuming $"SSE"(u_0) > 0$ (otherwise every point coincides with its leaf representative and the coreset is exact), conditional independence of branch choices yields the telescoping product
  $
    "Pr"("reach " v) = product_(r=0)^(h-1) ("SSE"(u_(r+1))) / ("SSE"(u_r)) = ("SSE"(v)) / ("SSE"(u_0)) ,
  $
  since $"SSE"(u_r) = "SSE"(u_(r, L)) + "SSE"(u_(r, R))$ at every internal node. Conditioned on reaching $v$, point $x in Q(v)$ is drawn with probability $w(x) norm(x - q(v))_2^2 / "SSE"(v)$. Multiplying cancels $"SSE"(v)$ identically, and noting $"SSE"(u_0) = sum_(u in cal(L)_t) "SSE"(u) = sum_(y in Q) w(y) norm(y - q(v(y)))_2^2$ establishes part 1.

  For part 2, for each final leaf $v in cal(L)_m$, write $x - q(v) = (x - c(v)) + (c(v) - q(v))$. Since $sum_(x in Q(v)) w(x)(x - c(v)) = bold(0)$, expanding the squared norm yields
  $
    "SSE"(v) = sum_(x in Q(v)) w(x) norm(x - c(v))_2^2 + W(v) norm(q(v) - c(v))_2^2 .
  $
  Summing over all $v in cal(L)_m$ and rearranging proves part 2.
]

#subsec("B", "Exact Amortized Carry Count of Merge-and-Reduce")

The binary bucket cascade guarantees both logarithmic space and $O(1)$ amortized reductions per chunk, independently of stream length $n$.

#lemma(title: "Exact Reduction Count and Amortized Streaming Cost", origin: "Original Derivation")[
  Let a stream of $n$ points be ingested in chunks of size $m$, with $M = floor(n / m) >= 1$ completed chunks and binary popcount (number of $1$-bits in the binary expansion of $M$) denoted $s_2 (M)$. Then:
  1. The number of active buckets never exceeds $floor(log_2 M) + 2$, requiring $O(m d log(n / m))$ working space.
  2. Across the entire stream, `"Reduce"` is invoked *exactly*
  $
    N_"red" (n) = M - s_2 (M) < n / m
  $
  times, and no stream point undergoes more than $h_"max" = floor(log_2 M)$ sequential reductions.
]

#proof[
  At any time after $M$ full chunks have been inserted, bucket $B_i$ ($i >= 1$) is occupied iff the $(i - 1)$-th bit of the binary representation of $M$ is $1$, while $B_0$ buffers at most $m - 1$ residual points. Since $M$ uses $floor(log_2 M) + 1$ bits, at most $floor(log_2 M) + 2$ buckets are allocated, each storing at most $m$ points in $RR^d$, proving part 1.

  For part 2, define the potential $Psi_r$ as the number of occupied buckets among $B_1, B_2, dots$ immediately after chunk $r$ finishes cascading ($Psi_0 = 0$). Each full chunk introduces $1$ new bucket ($+1$) and each `"Reduce"` call merges $2$ full buckets into $1$ ($-1$). Summing over $r = 1, dots, M$ gives $Psi_M = M - N_"red" (n)$. Because $Psi_M$ equals the exact number of $1$-bits in $M$, namely $s_2 (M) >= 1$, we obtain $N_"red" (n) = M - s_2 (M) < M <= n / m$. Finally, a bucket at level $B_i$ is formed by combining two buckets from level $B_(i-1)$, so points in $B_i$ have undergone exactly $i - 1 <= floor(log_2 M)$ reductions.
]

#subsec("C", "Merge-and-Reduce Coreset Composability Theorem")

The soundness of streaming compression rests on the algebraic closure of strong $(k, epsilon)$-coresets under disjoint union and sequential composition @harpeled2004coresets @ackermann2012streamkmpp.

#theorem(title: "Exact Coreset Closure and Height-h Distortion Amplification", origin: "Classical — Cited")[
  Let $epsilon in (0, 1]$ and $h >= 1$ be the maximum Merge-and-Reduce tree height.
  1. *Merge closure:* If $S_1, S_2$ are $(k, epsilon)$-coresets of disjoint sets $P_1, P_2$, then $S_1 union S_2$ is a $(k, epsilon)$-coreset of $P_1 union P_2$.
  2. *Transitivity:* If $S'$ is a $(k, delta)$-coreset of $S$ and $S$ is a $(k, epsilon)$-coreset of $P$, then $S'$ is a $(k, epsilon + delta + epsilon delta)$-coreset of $P$.
  3. *Tree error budget:* If every `"Reduce"` node in a Merge-and-Reduce tree of height $h <= ceil(log_2 (n / m))$ produces a $(k, epsilon')$-coreset with $epsilon' = epsilon / (2 h)$ and failure probability $delta' = delta m / n$, then with probability at least $1 - delta$ the final coreset $S$ is a strong $(k, epsilon)$-coreset of $P$.
]

#proof[
  For part 1, for every $k$-center set $C subset RR^d$, disjointness implies $Phi(P_1 union P_2, C) = Phi(P_1, C) + Phi(P_2, C)$ and $Phi(S_1 union S_2, C) = Phi(S_1, C) + Phi(S_2, C)$. Summing $(1 - epsilon) Phi(P_r, C) <= Phi(S_r, C) <= (1 + epsilon) Phi(P_r, C)$ over $r in {1, 2}$ proves Merge closure.

  For part 2, chaining the two-sided inequalities gives
  $
    (1 - delta)(1 - epsilon) Phi(P, C) <= Phi(S', C) <= (1 + delta)(1 + epsilon) Phi(P, C) .
  $
  Noting $(1 + delta)(1 + epsilon) = 1 + (epsilon + delta + epsilon delta)$ and $(1 - delta)(1 - epsilon) = 1 - (epsilon + delta - epsilon delta) >= 1 - (epsilon + delta + epsilon delta)$ proves transitivity.

  For part 3, by Lemma 2 there are strictly fewer than $n / m$ calls to `"Reduce"`, so a union bound ensures every call succeeds simultaneously with probability at least $1 - (n / m) delta' = 1 - delta$. Inducting on tree height $h$ via parts 1 and 2, the root coreset $S$ satisfies
  $
    (1 - epsilon')^h Phi(P, C) <= Phi(S, C) <= (1 + epsilon')^h Phi(P, C) .
  $
  For the lower bound, Bernoulli's inequality yields $(1 - epsilon / (2 h))^h >= 1 - epsilon / 2 > 1 - epsilon$. For the upper bound, $1 + x <= e^x$ gives $(1 + epsilon / (2 h))^h <= e^(epsilon / 2)$. Because $g(epsilon) = e^(epsilon / 2) - (1 + epsilon)$ is convex on $[0, 1]$ with $g(0) = 0$ and $g(1) = sqrt(e) - 2 < 0$, we have $e^(epsilon / 2) <= 1 + epsilon$ for all $epsilon in (0, 1]$, completing the proof.
]

#subsec("D", "End-to-End Streaming O(log k) Approximation Guarantee")

Combining Theorem 3 with weighted $k$-Means++ seeding @arthur2007kmeanspp yields a rigorous end-to-end streaming approximation bound @ackermann2012streamkmpp.

#theorem(title: "End-to-End Expected Approximation of Coreset-Backed Weighted k-Means++", origin: "Original Unified Proof")[
  Conditioned on the event that the Merge-and-Reduce root $S$ is a strong $(k, epsilon)$-coreset of $P$ with $epsilon in (0, 1 / 2]$, let $hat(C)$ be the output of weighted $k$-Means++ seeding followed by any number of weighted Lloyd iterations on $S$. Then
  $
    EE[Phi(P, hat(C))] <= (1 + epsilon) / (1 - epsilon) dot.op 8 (log k + 2) "OPT"_k (P) .
  $
  In particular, for $epsilon <= 1 / 3$, $EE[Phi(P, hat(C))] <= 16 (log k + 2) "OPT"_k (P)$.
]

#proof[
  Let $C^((0))$ be the initial center set produced by weighted $k$-Means++ on $S$. For rational weights, replicating each weighted point recovers unweighted $k$-Means++ @arthur2007kmeanspp; by continuity of sampling probabilities and $Phi(S, C)$ in weights $w$, the Arthur--Vassilvitskii expectation bound holds for arbitrary $w : S -> RR_(>0)$:
  $
    EE[Phi(S, C^((0)))] <= 8 (log k + 2) "OPT"_k (S) .
  $
  Each subsequent weighted Lloyd iteration monotonically decreases the objective on $S$ (nearest-center reassignment and weighted-centroid updates never increase $Phi$), hence $Phi(S, hat(C)) <= Phi(S, C^((0)))$ pointwise.

  Now let $C_P^* = arg min_(abs(C) = k) Phi(P, C)$ be an exact optimal solution for the full stream $P$. By the upper coreset inequality applied to $C_P^*$,
  $
    "OPT"_k (S) <= Phi(S, C_P^*) <= (1 + epsilon) Phi(P, C_P^*) = (1 + epsilon) "OPT"_k (P) .
  $
  Conversely, applying the lower coreset inequality to the random output $hat(C)$ gives $Phi(P, hat(C)) <= (1 - epsilon)^(-1) Phi(S, hat(C))$. Taking expectations and chaining the three inequalities yields
  $
    EE[Phi(P, hat(C))] &<= 1 / (1 - epsilon) EE[Phi(S, hat(C))] <= 1 / (1 - epsilon) EE[Phi(S, C^((0)))] \
                       &<= (8 (log k + 2)) / (1 - epsilon) "OPT"_k (S) <= (1 + epsilon) / (1 - epsilon) dot.op 8 (log k + 2) "OPT"_k (P) .
  $
  Substituting $epsilon <= 1 / 3$ gives $(1 + epsilon) / (1 - epsilon) <= 2$, establishing the corollary.
]

#sec("III", "Complexity, Geometry and Operational Profile")

A single call to `"Reduce"(Q, m)` on $2m$ weighted points performs $m - 1$ leaf splits. If splits are balanced, tree depth is $O(log m)$ and total point-distance re-evaluations across all splits take $O(m d log m)$; in the worst case of a maximally skewed tree, `"Reduce"` takes $O(m^2 d)$ time. By Lemma 2, `"Reduce"` runs strictly fewer than $n / m$ times across the entire stream, yielding total streaming time $O(n m d)$ worst-case and $O(n d log m)$ expected under balanced splits—an $O(m d)$ (or $O(d log m)$) amortized cost per stream point, with practical recommended size $m = 200 k$ @ackermann2012streamkmpp. Working RAM is $O(m d log(n / m))$ across the logarithmic bucket cascade. Offline extraction merges at most $ceil(log_2 (n / m))$ buckets and executes $I_"off"$ weighted Lloyd iterations on $m$ points in $O(I_"off" m k d)$ time.

Geometrically, StreamKM++ preserves the *global Euclidean $k$-Means Voronoi geometry* uniformly across all dimensions via $D^2$-importance sampling: unlike fixed-radius micro-clusters (CluStream, DenStream), coreset tree splits automatically adapt resolution to local variance, placing denser prototypes in high-cost regions and remaining immune to fixed-radius high-dimensional collapse.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Online Time (Total)], [$O(n m d)$ worst-case; $O(n d log m)$ balanced Coreset Tree.],
    [Offline Query Time], [$O(m^2 d log(n / m) + I_"off" m k d)$ via bucket merge + weighted $k$-Means++.],
    [Working Memory], [$O(m d log(n / m))$; sublinear in stream length $n$.],
    [Approximation Guarantee], [$O(log k)$ expected $k$-Means cost with $(1 plus.minus epsilon)$ coreset (Thm. 4).],
    [Inductive Bias], [Convex Euclidean Voronoi cells minimizing global SSE.],
    [Key Hyperparameters], [Cluster count $k$ and coreset bucket size $m$ (default $m approx 200 k$).],
    [Concept Drift Handling], [None natively: infinite-window uniform weighting over all past data.],
  ),
  caption: [Operational and computational profile of StreamKM++.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Web-Scale Embedding Quantization and Vector Search Indexing:* Large retrieval systems construct IVF (Inverted File) $k$-Means codebooks over billions of high-dimensional neural embeddings ($d in [128, 1536]$). StreamKM++ streams shards from disk in one pass with $O(m d log(n / m))$ RAM while preserving the exact Euclidean quantization SSE within $(1 plus.minus epsilon)$ (Theorem 4).
2. *High-Throughput Astronomy and Particle Physics Catalogs:* Sky surveys and collider triggers emit stationary multi-terabyte feature streams that exceed node memory. Merge-and-Reduce produces mergeable coresets across distributed compute partitions (by Merge closure, Theorem 3), enabling exact parallel reduction before final clustering.
3. *Large-Scale Customer Lifetime Value and Retail Segmentation:* Enterprise transaction tables with tens of millions of customer records require rapid prototyping across multiple candidate $k$ values. Because a single $(k_"max", epsilon)$-coreset serves queries for any $k <= k_"max"$, analysts sweep $k$ offline in milliseconds.
4. *Genomic $k$-mer and Metagenomic Sequence Sketching:* High-dimensional $k$-mer frequency vectors suffer severe distance concentration where radius-threshold stream algorithms fail. StreamKM++'s variance-adaptive $D^2$-sampling concentrates coreset mass on informative cluster modes without requiring a brittle Euclidean radius $epsilon$.
5. *Edge-to-Cloud Distributed Telemetry Aggregation:* Because strong coresets satisfy exact Merge closure (Theorem 3), edge devices can independently maintain local StreamKM++ coresets and transmit only $m$ weighted centroids upstream for lossless cloud-level global clustering.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

StreamKM++ has three clear limitations. First, *concept drift blindness*: as proved in Proposition 5 below, infinite-window unweighted aggregation lets stale historical mass overwhelm recent distribution shifts. Second, *non-convex manifolds*: because the downstream objective is Euclidean $k$-Means, intertwined crescents or nested rings are partitioned into convex cells. Third, *adversarial extreme outliers*: because $D^2$-sampling weights points quadratically by distance, a tiny fraction of corrupted mega-outliers can absorb disproportionate coreset slots unless pre-filtered.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Stationary Massive Euclidean Data], [Excellent: single pass, sublinear RAM, $O(log k)$ bound.], [Set $m >= 200 k$; compare BICO @fichtenberger2013bico / Mini-Batch @sculley2010webscale.],
    [High-Dimensional Feature Streams], [Strong: $D^2$-sampling avoids fixed-radius collapse.], [Sanitize extreme outliers before $D^2$ seeding.],
    [Drifting / Non-Stationary Streams], [Weak: stale mass dominates coreset weights (Prop. 5).], [Use CluStream @aggarwal2003clustream or DenStream @cao2006denstream.],
    [Non-Convex Noisy Manifolds], [Weak: convex $k$-Means Voronoi bias cuts curves.], [Use DenStream @cao2006denstream or HDBSCAN @campello2013hdbscan.],
    [Distributed MapReduce Streams], [Excellent: Merge closure enables tree aggregation.], [Ensure per-node $epsilon' = epsilon / (2 h)$ budget.],
  ),
  caption: [Domain suitability and failure-mode matrix for StreamKM++.],
)

#sec("V", "Structural Relations and Empirical Behavior")

StreamKM++ bridges classical $k$-Means++ @arthur2007kmeanspp and geometric coreset theory @harpeled2004coresets, directly motivating BICO @fichtenberger2013bico (which accelerates the online phase by fusing BIRCH's CF-tree @zhang1996birch with coreset projections) and scalable parallel seeding $k$-Means$||$ @bahmani2012scalable. Its sharpest contrast with CluStream @aggarwal2003clustream and DenStream @cao2006denstream lies in its treatment of temporal non-stationarity.

#proposition(title: "Quantitative Stale-Mass Dominance under Abrupt Concept Drift", origin: "Original Derivation")[
  Suppose a stream consists of $N_1$ points $P_1$ drawn from an obsolete regime $cal(D)_1$ followed by $N_2$ points $P_2$ drawn from an active regime $cal(D)_2$, and let $(S, w_S)$ be StreamKM++'s final root coreset. Then:
  1. Total coreset weight and first moment are conserved identically: $sum_(x in S) w_S (x) = N_1 + N_2$, so the optimal $1$-Means coreset centroid deviates from the active mean $mu_2$ by
  $
    norm(c_S^* - mu_2)_2 = (N_1) / (N_1 + N_2) norm(mu_1 - mu_2)_2 -> norm(mu_1 - mu_2)_2 .
  $
  2. If $(S, w_S)$ is a $(k, epsilon)$-coreset of $P = P_1 union P_2$, then for any $k$-center set $C$ incurring positive per-point MSE $overline(phi)_1 (C) > 0$ on $P_1$ and bounded MSE $overline(phi)_2 (C)$ on $P_2$, the active regime's relative contribution to $Phi(P, C)$ satisfies
  $
    (Phi(P_2, C)) / (Phi(P, C)) = (N_2 overline(phi)_2 (C)) / (N_1 overline(phi)_1 (C) + N_2 overline(phi)_2 (C)) = O(N_2 / N_1) -> 0 ,
  $
  forcing $Phi(S, C)$ within $(1 plus.minus epsilon)$ of the obsolete objective $Phi(P_1, C)$ as $N_1 / N_2 -> oo$.
]

#proof[
  At every leaf $v$ of every `"Reduce"` call, the emitted centroid carries weight $W(v) = sum_(x in Q(v)) w(x)$ and weighted sum $W(v) c(v) = sum_(x in Q(v)) w(x) x$, conserving both zeroth and first moments identically across merges and reductions. By induction over the Merge-and-Reduce tree, $sum_(x in S) w_S (x) = N_1 + N_2$ and the optimal $1$-Means centroid of $S$ satisfies
  $
    c_S^* = (1 / (N_1 + N_2)) sum_(x in S) w_S (x) x = (N_1 mu_1 + N_2 mu_2) / (N_1 + N_2) ,
  $
  where $mu_1, mu_2$ are the sample means of $P_1, P_2$. Subtracting $mu_2$ yields $norm(c_S^* - mu_2)_2 = (N_1 / (N_1 + N_2)) norm(mu_1 - mu_2)_2 -> norm(mu_1 - mu_2)_2$ as $N_1 / N_2 -> oo$, proving part 1. For part 2, additivity $Phi(P, C) = N_1 overline(phi)_1 (C) + N_2 overline(phi)_2 (C)$ gives the exact ratio above; applying the coreset sandwich $(1 - epsilon) Phi(P, C) <= Phi(S, C) <= (1 + epsilon) Phi(P, C)$ shows that $(Phi(S, C)) / (Phi(P_1, C)) in [(1 - epsilon), (1 + epsilon)(1 + O(N_2 / N_1))]$, proving that stale historical mass dominates $Phi(S, C)$ unless explicit decay or windowing is applied.
]

*Empirical Benchmark Analysis:* On the stationary ClusterOverview four-dataset benchmark, StreamKM++ delivers the strongest overall balance of quality and speed across Algs. 192--195: Iris Sil = 0.455, ARI = 0.618 (highest ARI among all four); Customers Sil = 0.361 (highest Customers Silhouette); Multi-Shape Sil = 0.611 (ARI = N/A); High-D Sil = 0.302 (highest High-D Silhouette, surpassing Contrastive Clustering's 0.286, CluStream's 0.253, and DenStream's N/A collapse); and a remarkably low mean runtime of 0.550 ms (16.8$times$ faster than Contrastive Clustering, 4.1$times$ faster than CluStream, and 3.5$times$ faster than DenStream). These outcomes align precisely with our theory: on stationary datasets without concept drift, variance-adaptive $D^2$-coreset trees preserve Euclidean structure across arbitrary dimension (Theorem 4) in $O(m d)$ amortized time (Lemma 2), yielding only to DenStream (0.793 vs. 0.611) on strongly non-convex Multi-Shape geometry.
