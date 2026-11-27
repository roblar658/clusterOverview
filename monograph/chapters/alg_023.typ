#import "../template/ieee.typ": *

#alg_entry(
  id: 23,
  name: "K-Means|| (Scalable K-Means++)",
  category: "Cat. I: Centroid & Partitioning",
  year: "2012",
  paradigm: "Parallel Oversampled D²-Seeding & Coreset Reclustering",
  complexity: $O(n d r ell + r ell k d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset equipped with the Euclidean metric, and let $k in {1, dots, n}$ be the target number of clusters. For any non-empty center set $C subset.eq bb(R)^d$, define the point-to-set distance $D_C (x) := min_(c in C) norm(x - c)_2$ and the *clustering potential* (within-cluster inertia)
$ phi_X (C) := sum_(x in X) D_C (x)^2, quad phi^* := min_(|C| = k) phi_X (C). $
Classical K-Means++ @arthur2007kmeanspp (Alg. 2) achieves an $O(log k)$-approximation in expectation by drawing $k$ centers sequentially, sampling each center $c_t$ with probability $D_(C^((t-1))) (x)^2 / phi_X (C^((t-1)))$. However, this sequential dependence forces $k$ full passes over $X$ with $k$ global synchronization barriers—a severe bottleneck on distributed clusters (MapReduce, Spark) when $k in [10^3, 10^5]$.

*K-Means||* (Scalable K-Means++), introduced by Bahmani et al. @bahmani2012scalable, breaks this sequential bottleneck by oversampling $ell = Theta(k)$ centers *independently in parallel* during each of $r = O(log psi)$ rounds, where $psi = phi_X ({c_1})$ is the initial one-center potential. The resulting candidate set $C$ of size $O(r ell) << n$ is then compressed into a weighted coreset $(C, w)$ and reclustered into $k$ final centers on a single coordinator node.

#definition(title: "Parallel Oversampled D²-Seeding & Weighted Coreset")[
  Given oversampling parameter $ell > 0$ and round budget $r >= 1$, the *K-Means|| process* starts from a uniformly drawn $c_1 in X$ ($C^((0)) = {c_1}$) and, at round $t in {1, dots, r}$, draws a random subset $C'_t subset.eq X$ by including each $x in X$ independently as a Bernoulli random variable with probability
  $ p_t (x) := min(1, space (ell dot D_(C^((t-1))) (x)^2) / (phi_X (C^((t-1))))), $
  setting $C^((t)) = C^((t-1)) union C'_t$. After $r$ rounds, each $c in C := C^((r))$ receives Voronoi multiplicity
  $ w_c := |{x in X : arg min_(c' in C) norm(x - c')_2 = c}|, $
  inducing the weighted $k$-means objective $phi_((C, w)) (M) = sum_(c in C) w_c D_M (c)^2$ over $|M| = k$.
]

#algbox(
  title: "K-Means|| (Scalable K-Means++)",
  inputs: [Dataset $X = {x_1, dots, x_n} subset.eq bb(R)^d$, cluster count $k >= 1$, oversampling factor $ell > 0$ (typically $ell in [0.5 k, 2 k]$), rounds $r >= 1$ (typically $r = ceil(log psi)$ or $r in [5, 8]$).],
  outputs: [Initial $k$-center set $M subset.eq bb(R)^d$ (optionally refined by Lloyd iterations).],
)[
  1. Sample $c_1 in X$ uniformly at random; set $C <- {c_1}$, $C'_0 <- {c_1}$.
  2. Initialize $d^2 [x] <- norm(x - c_1)_2^2$ for all $x in X$ in parallel; aggregate $psi <- sum_(x in X) d^2 [x]$.
  3. *For* $t = 1, dots, r$ *do*:
  4. #h(1em) If $t > 1$, update $d^2 [x] <- min(d^2 [x], min_(c in C'_(t-1)) norm(x - c)_2^2)$ for all $x in X$ in parallel.
  5. #h(1em) Compute $phi_(t-1) <- sum_(x in X) d^2 [x]$ via a single all-reduce; if $phi_(t-1) = 0$, *break*.
  6. #h(1em) In parallel, form $C'_t subset.eq X$ by sampling each $x in X$ independently with probability $p_t (x) = min(1, ell dot d^2 [x] / phi_(t-1))$.
  7. #h(1em) Broadcast $C'_t$ and set $C <- C union C'_t$.
  8. *Weighting pass:* For each $x in X$ in parallel, find $pi(x) = arg min_(c in C) norm(x - c)_2$; aggregate weights $w_c = |pi^(-1) (c)|$ for all $c in C$.
  9. *Weighted reclustering:* Run weighted K-Means++ (Alg. 2) or local search on $(C, w)$ to obtain $k$ centers $M$.
  10. *Return* $M$ (and optionally run Lloyd's K-Means, Alg. 1, initialized at $M$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational results governing K-Means||: (1) exact adaptive-martingale Chernoff concentration of candidate set size $|C|$; (2) a sharp $ell_2$-Minkowski distortion transfer bound from the weighted coreset $(C, w)$ back to $X$; (3) single-round geometric potential contraction yielding an $O(1)$ bicriteria and $O(alpha)$ end-to-end guarantee; and (4) an explicit multiscale counterexample proving that fixed-round seeding fails when the aspect ratio grows faster than $r$ rounds can resolve.

#subsec("A", "Adaptive Supermartingale Concentration of Candidate Set Size")

Because round-$t$ inclusion probabilities $p_t (x)$ depend adaptively on the random centers $C^((t-1))$ sampled in earlier rounds, the total candidate count $|C^((r))|$ is not a sum of independent variables. Nevertheless, a deterministic normalization invariant yields an exact Chernoff tail bound.

#lemma(title: "Adaptive Chernoff Bound on Oversampled Cardinality", origin: "Original Proof")[
  Let $N_t = |C'_t|$ be the number of points sampled in round $t$ of Algorithm 1, and let $cal(F)_t = sigma(c_1, C'_1, dots, C'_t)$ be the natural filtration. Then $bb(E)[|C^((r))|] <= 1 + r ell$, and for every $delta > 0$,
  $ bb(P)[sum_(t=1)^r N_t >= (1 + delta) r ell] &<= (e^delta / ((1 + delta)^(1 + delta)))^(r ell) \
    &<= exp(- (delta^2) / (2 + delta) r ell). $
]

#proof[
  Conditioned on $cal(F)_(t-1)$, the inclusion indicators $Z_(t, x) = bb(1){x in C'_t}$ are mutually independent Bernoulli random variables with parameters $p_t (x) = min(1, ell D_(C^((t-1))) (x)^2 / phi_(t-1))$. Summing over $X$ yields the *deterministic* conditional mean bound
  $ mu_t := sum_(x in X) p_t (x) <= sum_(x in X) (ell D_(C^((t-1))) (x)^2) / (phi_(t-1)) = ell, $
  on every sample path (with $mu_t = 0 <= ell$ if $phi_(t-1) = 0$). Taking expectations gives $bb(E)[|C^((r))|] <= 1 + sum_(t=1)^r bb(E)[N_t] <= 1 + r ell$.

  For the upper tail, fix any $lambda > 0$. By conditional independence across $x in X$ and the elementary inequality $1 + u <= e^u$,
  $ bb(E)[e^(lambda N_t) | cal(F)_(t-1)] &= product_(x in X) (1 + p_t (x)(e^lambda - 1)) \
    &<= exp((e^lambda - 1) sum_(x in X) p_t (x)) <= exp(ell (e^lambda - 1)). $
  Define $M_t := exp(lambda sum_(s=1)^t N_s - t ell (e^lambda - 1))$ with $M_0 = 1$. The inequality above states $bb(E)[M_t | cal(F)_(t-1)] <= M_(t-1)$, so ${M_t}$ is a non-negative supermartingale. By induction, $bb(E)[M_r] <= 1$, i.e., $bb(E)[exp(lambda sum_(t=1)^r N_t)] <= exp(r ell (e^lambda - 1))$. Applying Markov's inequality with $lambda = ln(1 + delta) > 0$ yields
  $ bb(P)[sum_(t=1)^r N_t >= (1 + delta) r ell] &<= exp(r ell (e^lambda - 1 - lambda (1 + delta))) \
    &= (e^delta / ((1 + delta)^(1 + delta)))^(r ell). $
  Finally, define $g(delta) = (1 + delta) ln(1 + delta) - delta - delta^2 / (2 + delta)$. We have $g(0) = 0$, $g'(0) = 0$, and for all $delta > 0$,
  $ g''(delta) &= 1 / (1 + delta) - 8 / ((2 + delta)^3) \
    &= (4 delta + 6 delta^2 + delta^3) / ((1 + delta)(2 + delta)^3) > 0. $
  Hence $g(delta) > 0$ for all $delta > 0$, establishing the simplified exponential tail bound.
]

#subsec("B", "Sharp L2-Minkowski Coreset Distortion Transfer Bound")

Step 9 of Algorithm 1 replaces the full dataset $X$ by the weighted Voronoi projection $(C, w)$. Standard analyses apply $(a + b)^2 <= 2 a^2 + 2 b^2$ twice, incurring a leading factor of $4 alpha$. Formulating the transfer directly in $ell_2 (bb(R)^n)$ via Minkowski's inequality tightens the leading constant by a factor of two.

#theorem(title: "Tight Minkowski Coreset Projection Bound", origin: "Original Proof")[
  Let $C subset.eq X$ be any candidate center set with weighted projection $(C, w)$ defined by $w_c = |pi^(-1) (c)|$, and let $M$ be any $k$-center set achieving an $alpha$-approximation on $(C, w)$, i.e., $phi_((C, w)) (M) <= alpha phi_((C, w))^*$. Then deterministically,
  $ sqrt(phi_X (M)) <= (1 + sqrt(alpha)) sqrt(phi_X (C)) + sqrt(alpha) sqrt(phi^*). $
  In particular, squaring and applying AM–GM yields the explicit linear bound
  $ phi_X (M) <= (2 alpha + 3 sqrt(alpha) + 1) phi_X (C) + (2 alpha + sqrt(alpha)) phi^*. $
]

#proof[
  Define three non-negative vectors $u, v, z in bb(R)^n$ indexed by $x in X$:
  $ v_x := norm(x - pi(x))_2, quad u_x := D_M (pi(x)), quad z_x := D_(M^*) (x), $
  where $M^*$ attains $phi^* = phi_X (M^*)$. By definition of $w_c$,
  $ &norm(v)_2 = sqrt(phi_X (C)), quad norm(z)_2 = sqrt(phi^*), \
    &norm(u)_2 = sqrt(phi_((C, w)) (M)). $
  For every $x in X$, the Euclidean triangle inequality implies $D_M (x) <= norm(x - pi(x))_2 + D_M (pi(x)) = v_x + u_x$. Taking the $ell_2$-norm across $x in X$ and invoking Minkowski's triangle inequality gives
  $ sqrt(phi_X (M)) &<= norm(v + u)_2 <= norm(v)_2 + norm(u)_2 \
    &= sqrt(phi_X (C)) + sqrt(phi_((C, w)) (M)). $
  Applying the exact same triangle inequality to the optimal center set $M^*$, for every $x in X$ we have $D_(M^*) (pi(x)) <= norm(pi(x) - x)_2 + D_(M^*) (x) = v_x + z_x$, hence by Minkowski's inequality and optimality of $phi_((C, w))^*$,
  $ sqrt(phi_((C, w))^*) &<= sqrt(phi_((C, w)) (M^*)) <= norm(v)_2 + norm(z)_2 \
    &= sqrt(phi_X (C)) + sqrt(phi^*). $
  Substituting $sqrt(phi_((C, w)) (M)) <= sqrt(alpha) sqrt(phi_((C, w))^*)$ into the first inequality yields
  $ sqrt(phi_X (M)) <= (1 + sqrt(alpha)) sqrt(phi_X (C)) + sqrt(alpha) sqrt(phi^*). $
  Setting $A = sqrt(phi_X (C))$, $B = sqrt(phi^*)$, expanding the square, and applying $2 A B <= A^2 + B^2$ to the cross term $2 sqrt(alpha)(1 + sqrt(alpha)) A B$ produces
  $ phi_X (M) &<= (1 + sqrt(alpha))^2 A^2 + sqrt(alpha)(1 + sqrt(alpha))(A^2 + B^2) + alpha B^2 \
    &= (2 alpha + 3 sqrt(alpha) + 1) phi_X (C) + (2 alpha + sqrt(alpha)) phi^*, $
  completing the proof.
]

#subsec("C", "Geometric Potential Contraction and End-to-End Approximation")

We now prove the hallmark bicriteria contraction property of parallel $D^2$-sampling established by Bahmani et al. @bahmani2012scalable, via a concise self-contained core-hitting argument.

#theorem(title: "Single-Round Contraction & O(alpha) End-to-End Bound", origin: "Classical — Cited")[
  Let $C = C^((t-1))$ be any non-empty center set with potential $phi = phi_X (C) > 0$, and let $C' = C'_t$ be sampled in one round of Algorithm 1 with oversampling factor $ell > 0$ @bahmani2012scalable. Then
  $ bb(E)[phi_X (C union C') | C] <= 12 phi^* + (8 k) / (e ell) phi_X (C). $
  Consequently, whenever $ell > 8 k / e$ (so $rho := 8 k / (e ell) < 1$), after $r >= ceil(log_(1/rho) (psi / phi^*))$ rounds the candidate set $C^((r))$ satisfies $bb(E)[phi_X (C^((r)))] = O(phi^*)$, and Algorithm 1 with an $alpha$-approximate weighted solver yields $bb(E)[phi_X (M)] = O(alpha) phi^*$ (specifically $O(log k) phi^*$ with weighted K-Means++ @arthur2007kmeanspp and $O(1) phi^*$ with local search @kanungo2004local).
]

#proof[
  Let ${X_1, dots, X_k}$ be an optimal $k$-means partition of $X$ with centroids $mu_j$ and optimal costs $phi_j^* = sum_(x in X_j) norm(x - mu_j)_2^2$, so $phi^* = sum_(j=1)^k phi_j^*$. Denote current cluster inertia by $phi_j = sum_(x in X_j) D_C (x)^2$. For each $j$, define the *core subset*
  $ G_j := {y in X_j : norm(y - mu_j)_2^2 <= (2 phi_j^*) / (|X_j|)}. $
  By Markov's inequality, $|X_j backslash G_j| < |X_j| / 2$ (trivially $G_j = X_j$ if $phi_j^* = 0$), so $|G_j| >= |X_j| / 2$. Whenever any point $y in G_j$ is sampled into $C'$, the parallel-axis identity gives
  $ phi_(X_j) (C union C') &<= sum_(x in X_j) norm(x - y)_2^2 \
    &= phi_j^* + |X_j| norm(y - mu_j)_2^2 <= 3 phi_j^*. $
  Now split clusters into two regimes:
  - *Cheap clusters ($phi_j <= 12 phi_j^*$):* Monotonicity of potential immediately gives $phi_(X_j) (C union C') <= phi_j <= 12 phi_j^*$.
  - *Expensive clusters ($phi_j > 12 phi_j^*$):* For any $x in X_j$ and $y in G_j$, $D_C (x)^2 <= 2 D_C (y)^2 + 2 norm(x - y)_2^2$. Summing over $x in X_j$ and using $sum_(x in X_j) norm(x - y)_2^2 <= 3 phi_j^*$ gives $phi_j <= 2 |X_j| D_C (y)^2 + 6 phi_j^*$. Averaging over $y in G_j$ and using $|X_j| / |G_j| <= 2$ yields
    $ phi_j <= 4 phi_(G_j) (C) + 6 phi_j^* \
      ==> phi_(G_j) (C) >= (phi_j - 6 phi_j^*) / 4 > phi_j / 8. $
    Let $E_j = {C' inter G_j = emptyset}$ be the event that $G_j$ is missed. If some $y in G_j$ has $p_t (y) = 1$, then $bb(P)[E_j | C] = 0$; otherwise $p_t (y) = ell D_C (y)^2 / phi$ for all $y in G_j$, so
    $ bb(P)[E_j | C] &= product_(y in G_j) (1 - p_t (y)) \
      &<= exp(- (ell phi_(G_j) (C)) / phi) <= exp(- (ell phi_j) / (8 phi)). $
    Conditioned on $E_j^c$, $phi_(X_j) (C union C') <= 3 phi_j^*$; conditioned on $E_j$, it is at most $phi_j$. Using $u e^(-u) <= e^(-1)$ for $u = ell phi_j / (8 phi) >= 0$,
    $ bb(E)[phi_(X_j) (C union C') | C] &<= 3 phi_j^* + phi_j exp(- (ell phi_j) / (8 phi)) \
      &<= 3 phi_j^* + (8 / (e ell)) phi. $
  Summing over all $k$ clusters (bounding $3 phi_j^* <= 12 phi_j^*$) proves $bb(E)[phi_X (C union C') | C] <= 12 phi^* + (8 k / (e ell)) phi$. Unrolling this affine recursion over $r$ rounds with $rho = 8 k / (e ell) < 1$ gives
  $ bb(E)[phi_X (C^((r)))] <= rho^r psi + 12 / (1 - rho) phi^*. $
  For $r >= ceil(log_(1/rho) (psi / phi^*))$, $rho^r psi <= phi^*$, so $bb(E)[phi_X (C^((r)))] <= (1 + 12/(1 - rho)) phi^*$. Taking expectations in Theorem 1 (which is linear in $phi_X (C^((r)))$) yields $bb(E)[phi_X (M)] = O(alpha) phi^*$.
]

#subsec("D", "Multiscale Lower Bound: Necessity of Logarithmic Rounds")

Practitioners frequently run K-Means|| for a small fixed number of rounds $r in [2, 5]$. The following counterexample proves rigorously that when the data exhibits a geometric hierarchy of cluster scales, *each* round of parallel $D^2$-sampling can resolve at most one scale, making $r = Omega(log psi)$ unavoidable in the worst case.

#theorem(title: "Unbounded Ratio Under Fixed-Round Multiscale Cascades", origin: "Original Counterexample")[
  For any fixed round count $r >= 1$, any finite oversampling factor $ell > 0$, and any target factor $B > 0$, there exists a one-dimensional dataset $X subset.eq bb(R)$ with $k = r + 2$ clusters ($n = 2(r + 2)$ points) such that after $r$ rounds of Algorithm 1,
  $ (bb(E)[phi_X (C^((r)))]) / (phi^*) > B. $
]

#proof[
  Set $R := ceil(3 sqrt(r (r + 1) ell)) + 4$ and let $epsilon in (0, 1/4)$ be a tiny intra-cluster radius chosen below. Construct $k = r + 2$ two-point clusters $X_0, X_1, dots, X_(r+1) subset.eq bb(R)$ centered at $a_0 = 0$ and $a_j = R^(j - 1)$ for $j = 1, dots, r + 1$, with $X_j = {a_j - epsilon, a_j + epsilon}$. Placing one center at each $a_j$ gives optimal cost $phi^* = 2(r + 2) epsilon^2$.

  Condition on the initial uniform draw landing in $X_0$, which occurs with probability $1 / (r + 2)$. For $t in {0, 1, dots, r}$, define the event
  $ A_t := {C^((t)) inter union.big_(j=1)^(r + 1 - t) X_j = emptyset}. $
  Clearly $bb(P)[A_0 | c_1 in X_0] = 1$. Now fix $t in {1, dots, r}$ and condition on $A_(t-1)$. Because no point in $X_1, dots, X_(r + 2 - t)$ has been sampled yet, the largest uncovered cluster index $m_t := max{j >= 1 : C^((t-1)) inter X_j = emptyset}$ satisfies $m_t >= r + 2 - t$. For any point $x in X_(m_t)$, every center in $C^((t-1))$ lies either in $X_0$ (distance $>= R^(m_t - 1) - 2 epsilon$) or in $X_j$ with $j >= m_t + 1$ (distance $>= (R - 1) R^(m_t - 1) - 2 epsilon > R^(m_t - 1)$ since $R >= 4$). Because $2(1 - 1/8)^2 > 1$, the current potential satisfies
  $ phi_(t-1) >= 2 (R^(m_t - 1) - 2 epsilon)^2 >= R^(2(m_t - 1)). $
  Meanwhile, every point $x in union.big_(j=1)^(r + 1 - t) X_j$ satisfies $j <= m_t - 1$, so its distance to $c_1 in X_0$ is at most $R^(m_t - 2) + 2 epsilon <= (3/2) R^(m_t - 2)$. Hence every such point has sampling probability
  $ p_t (x) <= (ell dot (9/4) R^(2(m_t - 2))) / (R^(2(m_t - 1))) = (9 ell) / (4 R^2). $
  By a union bound over the $2(r + 1 - t)$ points in $union.big_(j=1)^(r + 1 - t) X_j$,
  $ bb(P)[A_t^c | A_(t-1)] <= 2(r + 1 - t) (9 ell) / (4 R^2). $
  Chaining across $t = 1, dots, r$ yields
  $ bb(P)[A_r^c | c_1 in X_0] &<= sum_(t=1)^r (9 ell (r + 1 - t)) / (2 R^2) \
    &= (9 r (r + 1) ell) / (4 R^2) <= 1 / 4, $
  by our choice of $R$. On the event $A_r$, cluster $X_1$ (at $a_1 = 1$) remains completely uncovered, and every center in $C^((r))$ lies at distance $>= 1 - 2 epsilon >= 1/2$ from $X_1$, forcing $phi_X (C^((r))) >= 2 (1/2)^2 = 1/2$. Therefore,
  $ (bb(E)[phi_X (C^((r)))]) / (phi^*) &>= (1 / (r + 2) dot 3 / 4 dot 1 / 2) / (2 (r + 2) epsilon^2) \
    &= 3 / (16 (r + 2)^2 epsilon^2). $
  Selecting $epsilon < min(1/4, (4(r + 2) sqrt(B))^(-1))$ forces the expected approximation ratio strictly above $B$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Using incremental distance caching (Step 4 of Algorithm 1), round $t$ computes squared Euclidean distances only between $X$ and the *newly added* batch $C'_(t-1)$, requiring $O(n d |C'_(t-1)|)$ floating-point operations. By Lemma 1, $sum_(t=1)^r |C'_t| = O(r ell)$ with high probability, so the total distributed work across all $r$ rounds is $O(n d r ell)$—identical in asymptotic FLOPs to sequential K-Means++ ($O(n d k)$) when $r ell = O(k log psi)$, yet requiring only $r = O(log psi)$ global MapReduce/All-Reduce synchronization rounds instead of $k$.

Across $P$ distributed workers, per-worker memory is $O((n / P + r ell) d)$ and total broadcast/gather network traffic is $O(r ell d P)$. Once the weighted coreset $(C, w)$ of cardinality $|C| = O(r ell)$ is gathered at the coordinator, weighted K-Means++ reclustering costs $O(r ell k d)$—entirely independent of $n$. In practice, setting $ell = 2 k$ and $r = 5$ (the Apache Spark MLlib default) achieves $O(1)$ empirical approximation while keeping $|C| <= 10 k$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Means||],
    [Time / Sync Complexity], [$O(n d r ell + r ell k d)$ FLOPs; only $r = O(log psi)$ parallel barriers (vs. $k$ for Alg. 2).],
    [Space / Communication], [$O((n/P + r ell) d)$ per worker; $O(r ell d P)$ network broadcast volume.],
    [Approximation Guarantee], [$O(1)$ bicriteria after $r = O(log psi)$ rounds ($ell = Omega(k)$); $O(log k)$ or $O(1)$ after Step 9.],
    [Geometric Invariances], [Invariant under orthogonal $Q in O(d)$, translation, and uniform scaling; sensitive to anisotropic feature scales.],
    [Outlier Sensitivity], [High $D^2$-probability draws outliers into $C$, but Step 8 assigns them unit weight $w_c = 1$, mitigating centroid pull during reclustering.],
    [Recommended Companions], [Lloyd's K-Means (Alg. 1) or Elkan/Hamerly (Algs. 13–14) for refinement; StreamKM++ (Alg. 195) or BICO (Alg. 78) for one-pass streams.],
  ),
  caption: [Theoretical guarantees and distributed operational profile of K-Means||.],
)

#sec("IV", "Engineering and Scientific Applications")

K-Means|| is the de facto initialization standard in distributed data-processing frameworks (Apache Spark MLlib, Dask, Ray, distributed FAISS) whenever both $n$ and $k$ exceed single-machine latency budgets.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Distributed Enterprise Telemetry & Cloud Log Partitioning:* Multi-terabyte behavioral and telemetry tables stored across HDFS or object storage cannot afford $k = 10^3$ sequential passes. K-Means||'s $r approx 5$ bulk-synchronous MapReduce rounds overlap network all-reduce latency with parallel BLAS-3 distance scans, delivering near-linear horizontal speedup across hundreds of executors.
2. *Billion-Scale Approximate Nearest Neighbor (IVF) Codebook Initialization:* Inverted-file (IVF) vector indexes require $k in [10^4, 2^(16)]$ coarse Voronoi centroids over $10^8$–$10^9$ neural embeddings. Sequential K-Means++ would demand tens of thousands of full corpus scans; K-Means|| initializes well-spread codebooks in $< 10$ distributed passes, preventing empty or imbalanced inverted lists.
3. *Single-Cell Atlas & Population-Scale Genomics:* Biobank-scale PCA embeddings ($n > 10^6$ cells, $d approx 50$, $k approx 200$ cell subtypes) exhibit multi-order-of-magnitude abundance disparities. Because $D^2$-oversampling gives rare distant subpopulations disproportionate inclusion probability in $C$ without overwhelming memory, K-Means|| captures rare cell populations missed by uniform Forgy seeding.
4. *High-Resolution Satellite Imagery & Remote Sensing Quantization:* Continental multispectral rasters contain billions of spectral pixels partitioned into $k$ land-cover prototypes. K-Means|| maps directly onto tiled raster partitions, broadcasting only compact $O(r ell d)$ spectral center tables between tiles.
5. *Large-Scale Cosmological & Molecular Dynamics State Discretization:* Markov State Models (MSMs) in molecular dynamics discretize $10^7+$ conformation frames into $k in [500, 5000]$ microstates; K-Means||'s parallel oversampling discovers metastable transition basins rapidly before local Lloyd refinement.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-Means|| exhibits four principal failure regimes: (a) *extreme multiscale hierarchies under fixed small $r$* (Theorem 3), where a single distant cluster absorbs nearly all $D^2$ probability mass in early rounds; (b) *severe outlier contamination*, where isolated extreme points consume candidate budget and—if $k$ is large relative to $r ell$—may still claim final cluster slots during weighted K-Means++; (c) *non-convex manifolds and anisotropic covariances* inherited from the squared Euclidean $k$-means objective; and (d) *very small $n$ or $k$*, where oversampling adds stochastic variance without parallelization payoff.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Means|| Fits], [Caution / Alternative],
    [Spark / Distributed ETL], [$O(log psi)$ BSP rounds; linear $P$-worker scaling], [Use Mini-Batch for 1-pass streams],
    [IVF Vector Search ($k >= 10^4$)], [Eliminates $O(k)$ serial scan barrier], [Standardize / L2-normalize first],
    [Single-Cell Atlases], [Oversampling captures rare distant states], [Increase $r$ if multiscale (Thm. 3)],
    [Heavy-Tailed Finance], [Fast coarse seeding after clipping], [Use CLARA-MED (Alg. 25) or Trimmed (Alg. 26)],
    [Non-Convex Shapes], [Poor fit (Voronoi hyperplanes)], [Use DBSCAN (Alg. 31) or Spectral (Algs. 103–105)],
  ),
  caption: [Domain suitability matrix and mitigation guide for K-Means||.],
)

#sec("V", "Structural Relations and Empirical Behavior")

K-Means|| bridges sequential $D^2$-seeding (K-Means++, Alg. 2) and distributed coreset methods (BICO Alg. 78 @fichtenberger2013bico, StreamKM++ Alg. 195 @ackermann2012streamkmpp, and Har-Peled–Mazumdar coresets @harpeled2004coresets). The following exact non-asymptotic limit shows that K-Means|| continuously interpolates back to exact K-Means++ as the oversampling rate vanishes.

#proposition(title: "Exact Vanishing-Oversampling Convergence to K-Means++", origin: "Original Derivation")[
  Fix any non-empty center set $C subset.eq X$ with $phi_X (C) > 0$, let $q(x) := D_C (x)^2 / phi_X (C)$ denote the exact K-Means++ distribution, and let $C'$ be a single K-Means|| round with $ell in (0, 1]$ satisfying $ell max_(x in X) q(x) <= 1$. Then conditioned on drawing at least one center ($|C'| >= 1$):
  1. $bb(P)[|C'| >= 2 | |C'| >= 1] <= ell / (1 - ell / 2) <= 2 ell$, and
  2. for every $x in X$, $|bb(P)[C' = {x} | |C'| >= 1] - q(x)| <= ell q(x)$, so the total variation error relative to exact K-Means++ decays at rate $O(ell)$ as $ell -> 0^+$.
]

#proof[
  Because $ell q(x) <= 1$ for all $x$, $p(x) = ell q(x)$ with $sum_(x in X) p(x) = ell$. By Bernoulli's inequality and $1 - u <= e^(-u) <= 1 - u + u^2 / 2$ (for $u >= 0$),
  $ 1 - ell <= product_(y in X) (1 - ell q(y)) <= e^(-ell) <= 1 - ell + ell^2 / 2. $
  Subtracting from $1$ gives tight two-sided bounds on the non-empty draw probability:
  $ ell (1 - ell / 2) <= bb(P)[|C'| >= 1] <= ell. $
  For any single point $x in X$, $bb(P)[C' = {x}] = ell q(x) product_(y != x) (1 - ell q(y))$, which similarly satisfies
  $ ell q(x) (1 - ell) <= bb(P)[C' = {x}] <= ell q(x). $
  Summing the lower bound over $x in X$ gives $bb(P)[|C'| = 1] >= ell (1 - ell)$, hence $bb(P)[|C'| >= 2] = bb(P)[|C'| >= 1] - bb(P)[|C'| = 1] <= ell^2$, proving Claim 1:
  $ bb(P)[|C'| >= 2 | |C'| >= 1] &<= (ell^2) / (ell (1 - ell / 2)) \
    &= ell / (1 - ell / 2) <= 2 ell. $
  For Claim 2, dividing the bounds on $bb(P)[C' = {x}]$ by those on $bb(P)[|C'| >= 1]$ yields
  $ q(x)(1 - ell) &<= bb(P)[C' = {x} | |C'| >= 1] \
    &<= (q(x)) / (1 - ell / 2) <= q(x)(1 + ell), $
  where $(1 - ell / 2)^(-1) <= 1 + ell$ holds for all $ell in (0, 1]$ since $(1 - ell / 2)(1 + ell) = 1 + ell(1 - ell)/2 >= 1$. Thus $|bb(P)[C' = {x} | |C'| >= 1] - q(x)| <= ell q(x)$, and summing over $x in X$ bounds the singleton $L_1$ discrepancy by $ell$.
]

*Empirical Behavior on ClusterOverview Benchmarks:* Across the 4-dataset ClusterOverview evaluation suite, K-Means|| achieves the fastest mean runtime in its immediate cohort at *0.105 ms* (reflecting lightweight vectorized $D^2$ passes without pairwise matrix construction). On *Iris*, it records Silhouette = *0.473* but ARI = *0.431*—noticeably lower label agreement than deterministic medoid search (FastPAM ARI = 0.663, Alg. 24), illustrating that oversampling on small $n$ with overlapping Versicolor/Virginica boundaries increases seeding variance and favors spread-out boundary prototypes. On *Customers* (Sil = *0.301*), *Multi-Shape* (Sil = *0.609*, ARI = N/A), and *High-D Subspace+Noise* (Sil = *0.249*, ARI = N/A), its scores trail $L_1$-medoid methods (Algs. 24–25), directly corroborating the squared-$L_2$ sensitivity of $D^2$-seeding to non-convex geometry and high-dimensional noise directions.
