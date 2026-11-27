#import "../template/ieee.typ": *

#alg_entry(
  id: 178,
  name: "Adaptive Resonance Theory (ART / Fuzzy ART)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1987",
  paradigm: "Vigilance-Gated Resonant Prototype Learning",
  complexity: $cal(O)(T K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Adaptive Resonance Theory (ART), originated by Carpenter and Grossberg for binary patterns (ART1) @carpenter1987art and generalized to real-valued spaces via fuzzy set operations (Fuzzy ART) @carpenter1991fuzzyart, was explicitly engineered to resolve the *stability-plasticity dilemma*: how an online neural clustering system can remain plastic enough to form new categories for novel patterns while remaining stable enough to avoid catastrophic forgetting of previously learned clusters.

Let raw data points $bold(a)_1, dots, bold(a)_n$ be normalized into the unit hypercube $[0, 1]^d$. To prevent unbounded weight erosion (proven inevitable without normalization in Section II.C), Fuzzy ART applies *complement coding*, mapping each input $bold(a) in [0, 1]^d$ to the $2d$-dimensional vector
$
  bold(I) := (bold(a), bold(1) - bold(a)) in [0, 1]^(2 d),
$
which enforces constant $L_1$ norm $norm(bold(I))_1 = sum_(r=1)^d (a_r + 1 - a_r) = d$. The network maintains a dynamically growing bank of $K$ committed category templates $bold(w)_1, dots, bold(w)_K in [0, 1]^(2 d)$. Denote the fuzzy AND (componentwise minimum) of two vectors $bold(x), bold(y) in [0, 1]^(2 d)$ by $(bold(x) and bold(y))_r := min(x_r, y_r)$.

#definition(title: "Fuzzy ART Choice, Vigilance Match, and Resonance Rules")[
  Fix choice parameter $alpha > 0$, vigilance threshold $rho in [0, 1]$, and learning rate $beta in (0, 1]$. For a complement-coded input $bold(I) in [0, 1]^(2 d)$ and committed prototype $bold(w)_j in [0, 1]^(2 d)$, define the bottom-up *choice activation* and top-down *match ratio* by
  $
    T_j (bold(I)) := norm(bold(I) and bold(w)_j)_1 / (alpha + norm(bold(w)_j)_1),
    quad M_j (bold(I)) := norm(bold(I) and bold(w)_j)_1 / norm(bold(I))_1.
  $
  Category $j$ achieves *resonance* with $bold(I)$ if and only if $M_j (bold(I)) >= rho$. Upon resonance at the highest-choice qualifying category $J$, its template updates via
  $
    bold(w)_J^"new" = beta (bold(I) and bold(w)_J^"old")
      + (1 - beta) bold(w)_J^"old".
  $
  If every committed category fails vigilance ($M_j (bold(I)) < rho$ for all $j in [K]$), a new category $K + 1$ is committed with $bold(w)_(K+1) = bold(I)$.
]

Intuitively, $T_j (bold(I))$ measures how strongly input $bold(I)$ is contained within template $bold(w)_j$ (favoring tighter, more specific templates as $alpha -> 0^+$), while $M_j (bold(I)) >= rho$ enforces a hard top-down fidelity constraint preventing excessive prototype generalization.

#algbox(
  title: "Complement-Coded Fuzzy ART Clustering",
  inputs: [Scaled data ${bold(a)_m}_(m=1)^n subset [0, 1]^d$; vigilance $rho in [0, 1]$; choice $alpha > 0$; learning rate $beta in (0, 1]$ (fast learning: $beta = 1$).],
  outputs: [Committed count $K$, templates ${bold(w)_j}_(j=1)^K subset [0, 1]^(2 d)$ (boxes $R_j = [bold(u)_j, bold(v)_j]$), and assignments $c(m) in [K]$.],
)[
  + *Complement Code:* Map each $bold(a)_m |-> bold(I)_m = (bold(a)_m, bold(1) - bold(a)_m) in [0, 1]^(2 d)$. Initialize $K <- 0$.
  + *For each input $bold(I)_m$ (across epochs or online stream) do:*
    - Compute $T_j (bold(I)_m) = norm(bold(I)_m and bold(w)_j)_1 / (alpha + norm(bold(w)_j)_1)$ for $j = 1, dots, K$.
    - Sort candidate indices $(j_1, dots, j_K)$ in descending order of $T_j (bold(I)_m)$.
    - *Vigilance Search:* Set $"found" <- "false"$. For $r = 1, dots, K$:
      - Let $j = j_r$; evaluate $M_j (bold(I)_m) = norm(bold(I)_m and bold(w)_j)_1 / d$.
      - If $M_j (bold(I)_m) >= rho$: update $bold(w)_j <- beta (bold(I)_m and bold(w)_j) + (1 - beta) bold(w)_j$; assign $c(m) <- j$; set $"found" <- "true"$; break.
    - *Novelty Commitment:* If not $"found"$: increment $K <- K + 1$, initialize $bold(w)_K <- bold(I)_m$, and assign $c(m) <- K$.
  + *Return:* Learned hyperrectangle prototypes ${bold(w)_j}_(j=1)^K$ and cluster assignments $c$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present four rigorous results characterizing Fuzzy ART: an exact hyperrectangle $L_1$-perimeter invariant and vigilance size bound, finite-sample stability and epoch convergence under fast learning, a constructive counterexample proving unbounded forgetting and proliferation without complement coding, and a high-dimensional Hoeffding concentration theorem explaining Fuzzy ART's empirical degradation in high dimensions.

#subsec("A", "Hyperrectangle Geometry and Exact Vigilance Size Bound")

#theorem(title: "Exact L1-Perimeter Identity and Universal Box-Size Bound", origin: "Original Proof")[
  Under complement coding $bold(I) = (bold(a), bold(1) - bold(a)) in [0, 1]^(2 d)$ and any learning rate $beta in (0, 1]$, every committed template $bold(w)_j (t)$ admits the unique representation
  $
    bold(w)_j (t) = (bold(u)_j (t), bold(1) - bold(v)_j (t)),
    quad bold(0) <= bold(u)_j (t) <= bold(v)_j (t) <= bold(1),
  $
  defining an axis-aligned hyperrectangle $R_j (t) := product_(r=1)^d [u_(j r) (t), v_(j r) (t)] subset.eq [0, 1]^d$ whose $L_1$ side-length sum $|R_j (t)|_1 := sum_(r=1)^d (v_(j r) (t) - u_(j r) (t))$ satisfies the exact norm identity and uniform vigilance bound
  $
    norm(bold(w)_j (t))_1 = d - |R_j (t)|_1,
    quad |R_j (t)|_1 <= d (1 - rho).
  $
]

#proof[
  We argue by induction on the number of resonant updates to category $j$. Upon commitment at $bold(I) = (bold(a), bold(1) - bold(a))$, $bold(u)_j = bold(v)_j = bold(a)$, so $bold(0) <= bold(u)_j <= bold(v)_j <= bold(1)$, $|R_j|_1 = 0 <= d(1 - rho)$, and $norm(bold(w)_j)_1 = d = d - |R_j|_1$.

  Suppose the claim holds prior to a resonant update with $bold(I) = (bold(a), bold(1) - bold(a))$. Compute the two $d$-dimensional halves of $bold(I) and bold(w)_j$: for every coordinate $r in [d]$,
  $
    (bold(I) and bold(w)_j)_r &= min(a_r, u_(j r)), \
    (bold(I) and bold(w)_j)_(d + r) &= min(1 - a_r, 1 - v_(j r))
      = 1 - max(a_r, v_(j r)).
  $
  Substituting into $bold(w)_j^"new" = beta (bold(I) and bold(w)_j) + (1 - beta) bold(w)_j$ yields $bold(w)_j^"new" = (bold(u)_j^"new", bold(1) - bold(v)_j^"new")$ with
  $
    u_(j r)^"new" &= beta min(a_r, u_(j r)) + (1 - beta) u_(j r) <= u_(j r), \
    v_(j r)^"new" &= beta max(a_r, v_(j r)) + (1 - beta) v_(j r) >= v_(j r).
  $
  Hence $bold(0) <= bold(u)_j^"new" <= bold(u)_j <= bold(v)_j <= bold(v)_j^"new" <= bold(1)$, proving $R_j (t) subset.eq R_j (t + 1)$ and preserving the representation. Summing all $2d$ coordinates gives
  $
    norm(bold(w)_j^"new")_1
      = sum_(r=1)^d (u_(j r)^"new" + 1 - v_(j r)^"new")
      = d - sum_(r=1)^d (v_(j r)^"new" - u_(j r)^"new")
      = d - |R_j^"new"|_1.
  $
  Finally, because $bold(I) and bold(w)_j <= bold(w)_j$ coordinatewise and $beta in (0, 1]$, we have $bold(w)_j^"new" >= bold(I) and bold(w)_j$ coordinatewise. Combined with the resonance condition $norm(bold(I) and bold(w)_j)_1 = d M_j (bold(I)) >= rho d$, we obtain
  $
    d - |R_j^"new"|_1 = norm(bold(w)_j^"new")_1
      >= norm(bold(I) and bold(w)_j)_1 >= rho d,
  $
  which rearranges immediately to $|R_j^"new"|_1 <= d (1 - rho)$.
]

Under fast learning ($beta = 1$), $R_j (t)$ is precisely the minimal axis-aligned bounding box enclosing every input that has resonated with category $j$, and its expansion halts strictly once its $L_1$ perimeter reaches $d(1 - rho)$.

#subsec("B", "Stability-Plasticity and Exact Finite-Update Termination")

#theorem(title: "Monotone Weight Stability and Finite-Epoch Convergence", origin: "Classical — Cited")[
  Let $X = {bold(a)_1, dots, bold(a)_n} subset [0, 1]^d$ be any finite dataset presented repeatedly in arbitrary order to complement-coded Fuzzy ART @carpenter1991fuzzyart @carpenter1987art.
  + For any $beta in (0, 1]$ and any infinite input stream, every committed weight vector satisfies $bold(0) <= bold(w)_j (t + 1) <= bold(w)_j (t)$ coordinatewise and converges monotonically to a limit $bold(w)_j^* in [0, 1]^(2 d)$ with $norm(bold(w)_j^*)_1 >= rho d$.
  + Under fast learning ($beta = 1$), the total number of committed categories across all epochs satisfies $K <= n$, each category undergoes at most $n - 1$ strict weight updates, and the network reaches a permanently fixed weight state in at most $1 + K(n - 1)$ epochs.
]

#proof[
  Part (1) follows immediately from $(bold(I) and bold(w)_j)_r = min(I_r, w_(j r)) <= w_(j r)$, so $bold(w)_j (t + 1) = bold(w)_j (t) - beta (bold(w)_j (t) - bold(I) and bold(w)_j (t)) <= bold(w)_j (t)$. Being bounded below by $bold(0)$ and coordinatewise non-increasing, every sequence $bold(w)_j (t)$ converges monotonically, while Theorem II.A maintains $norm(bold(w)_j (t))_1 >= rho d$.

  For Part (2) with $beta = 1$, observe that if $bold(a)_m in R_j (t)$ then $u_(j r) <= a_(m r) <= v_(j r)$ for all $r in [d]$, whence $bold(I)_m and bold(w)_j = bold(w)_j$. Consequently, $M_j (bold(I)_m) = norm(bold(w)_j)_1 / d >= rho$, and resonance with $j$ leaves $bold(w)_j$ completely unchanged. On its very first presentation in epoch $1$, every point $bold(a)_m$ either resonates with an existing category (entering $R_j^"new"$) or commits a new category $R_k = {bold(a)_m}$. Since boxes only expand ($R_j (t) subset.eq R_j (t + 1)$), every point $bold(a)_m$ lies inside at least one committed box for all subsequent presentations and therefore always has at least one category satisfying $M_j (bold(I)_m) >= rho$. Thus no new categories can ever be created after epoch $1$, proving $K <= n$.

  Moreover, a strict update $bold(w)_j^"new" != bold(w)_j$ at input $bold(a)_m$ occurs only if $bold(a)_m in.not R_j$ prior to the update, after which $bold(a)_m in R_j^"new"$ permanently. Because $|X| = n$ and $R_j$ initially contains one point upon commitment, category $j$ can absorb at most $n - 1$ new points across all epochs, giving at most $K(n - 1)$ strict updates after epoch $1$. Any epoch with zero strict updates fixes all weights permanently, bounding convergence by $1 + K(n - 1)$ epochs.
]

#subsec("C", "Necessity of Complement Coding: Constructive Proliferation Proof")

#theorem(title: "Unbounded Erosion and Category Proliferation Without Complement Coding", origin: "Original Counterexample")[
  Consider uncoded Fuzzy ART operating directly on $bold(x) in (0, 1]^d$ with $beta = 1$, arbitrary $alpha > 0$, and any vigilance $rho in (0, 1)$. Fix any $gamma in (0, rho)$ and integer $L >= 2$, and define the $L$-point dataset $bold(x)_k := gamma^k bold(1) in (0, 1]^d$ for $k = 0, 1, dots, L - 1$, presented in decreasing order $(bold(x)_0, dots, bold(x)_(L-1))$ each epoch. Then:
  + Without complement coding, epoch $1$ collapses all $L$ points into a single eroded prototype $bold(w)_1 = gamma^(L-1) bold(1)$ that forgets every prior point ($M_1 (bold(x)_k) < rho$ for all $k < L - 1$), and $L$ epochs proliferate to $K = L$ singleton categories.
  + With complement coding on the exact same dataset and order, zero forgetting occurs and the total category count is uniformly bounded by $K <= ceil(1 / (1 - rho))$, independent of $L$.
]

#proof[
  For (1), we prove by induction on epoch $e in {1, dots, L}$ that at the end of epoch $e$, uncoded Fuzzy ART holds exactly $e$ committed categories with weights $bold(w)_m = gamma^(L - m) bold(1)$ for $m = 1, dots, e$.
  - *Base case ($e = 1$):* $bold(x)_0 = bold(1)$ commits $bold(w)_1 = bold(1)$. For each $k = 1, dots, L - 1$, prior weight is $gamma^(k-1) bold(1) > gamma^k bold(1)$, so $bold(x)_k and bold(w)_1 = bold(x)_k$, giving $M_1 (bold(x)_k) = norm(bold(x)_k)_1 / norm(bold(x)_k)_1 = 1 >= rho$. Thus $bold(x)_k$ resonates and erodes $bold(w)_1 <- gamma^k bold(1)$, finishing epoch $1$ at $bold(w)_1 = gamma^(L-1) bold(1)$. For every $k <= L - 2$, $M_1 (bold(x)_k) = gamma^(L - 1 - k) <= gamma < rho$, so $bold(w)_1$ rejects all earlier inputs.
  - *Inductive step ($e -> e + 1$, $e < L$):* At input $bold(x)_0$, every existing category $m <= e$ has $M_m (bold(x)_0) = gamma^(L - m) <= gamma < rho$, forcing a new category $bold(w)_(e+1) = bold(1)$. For each $1 <= k <= L - e - 1$, every $m <= e$ satisfies $L - m > k$, so $M_m (bold(x)_k) = gamma^(L - m - k) <= gamma < rho$ and resets; category $e + 1$ satisfies $M_(e+1) (bold(x)_k) = 1 >= rho$ and erodes sequentially down to $bold(w)_(e+1) = gamma^(L - e - 1) bold(1)$. For each remaining $k >= L - e$, category $m^* := L - k <= e$ matches $bold(x)_k$ exactly ($bold(w)_(m^*) = bold(x)_k$, $M_(m^*) = 1$), strictly dominates all $m > m^*$ in choice because $norm(bold(w)_(m^*))_1 < norm(bold(w)_m)_1$, while all $m < m^*$ fail vigilance ($M_m <= gamma < rho$). Thus $bold(x)_k$ resonates with $m^*$ without modifying any weights, completing the induction and yielding $K = L$ after $L$ epochs.

  For (2), under complement coding, $bold(a)_0$ starts $R_1 = [gamma^0, gamma^0]^d$, which absorbs every consecutive $bold(a)_k$ satisfying $|R_1|_1 = d(1 - gamma^k) <= d(1 - rho)$, i.e., $gamma^k >= 1 - (1 - rho)$. Each newly spawned box $R_j$ starts at a scalar $r_j < r_(j-1) - (1 - rho)$. Since all $r_j in [0, 1]$, at most $ceil(1 / (1 - rho))$ boxes can be committed regardless of $L$, and by Theorem II.B no point is ever forgotten.
]

#subsec("D", "High-Dimensional Concentration of Vigilance Match")

#proposition(title: "Order Sensitivity and High-Dimensional Vigilance Concentration", origin: "Original Proof")[
  + *Symmetry & Order Sensitivity:* Complement-coded Fuzzy ART is invariant to coordinate permutations and reflections $a_r |-> 1 - a_r$, not invariant to orthogonal rotations, and strictly dependent on sample presentation order.
  + *High-Dimensional Match Concentration:* Let $bold(a), bold(b)$ be independent draws from $"Unif"([0, 1]^d)$ with complement codes $bold(I)_a, bold(I)_b$. Then $bb(E)[M(bold(I)_a, bold(I)_b)] = 2 / 3$, and for every $epsilon > 0$,
    $
      bb(P)(|M(bold(I)_a, bold(I)_b) - 2 / 3| >= epsilon)
        <= 2 exp(-2 d epsilon^2).
    $
]

#proof[
  For (1), coordinate permutations and reflections merely permute the $2d$ entries of $bold(I)$, preserving componentwise minima and $L_1$ norms; axis-aligned bounding boxes clearly break rotational equivariance. To prove order sensitivity constructively in $d = 1$ with $rho = 0.5$ ($1 - rho = 0.5$), present ${0.0, 0.4, 0.8}$: order $(0.0, 0.4, 0.8)$ creates boxes $[0.0, 0.4]$ and $[0.8, 0.8]$, whereas order $(0.8, 0.4, 0.0)$ creates $[0.4, 0.8]$ and $[0.0, 0.0]$, assigning $0.4$ to opposite clusters.

  For (2), note that for each coordinate $r in [d]$,
  $
    min(a_r, b_r) + min(1 - a_r, 1 - b_r)
      = min(a_r, b_r) + 1 - max(a_r, b_r)
      = 1 - |a_r - b_r|.
  $
  Summing over $r = 1, dots, d$ and dividing by $d$ yields $M(bold(I)_a, bold(I)_b) = 1 - (1 / d) sum_(r=1)^d Z_r$, where $Z_r := |a_r - b_r| in [0, 1]$ are i.i.d. random variables with mean
  $
    bb(E)[Z_r] = 2 integral_0^1 integral_0^x (x - y) dif y dif x
      = integral_0^1 x^2 dif x = 1 / 3.
  $
  Thus $bb(E)[M(bold(I)_a, bold(I)_b)] = 2 / 3$. Applying Hoeffding's inequality to the average of $d$ independent $[0, 1]$-valued variables gives $bb(P)(|(1 / d) sum_(r=1)^d Z_r - 1 / 3| >= epsilon) <= 2 exp(-2 d epsilon^2)$.
]

Proposition II.D exposes an acute high-dimensional phase transition @beyer1999nearest: as $d$ grows, pairwise match scores collapse into a narrow band $2 / 3 +- cal(O)(d^(-1 / 2))$. Setting $rho$ slightly above $2 / 3$ triggers massive category proliferation, while setting $rho$ slightly below $2 / 3$ merges unrelated clusters into overlapping axis-aligned boxes that cut across Euclidean cluster boundaries.

#sec("III", "Complexity, Geometry and Operational Profile")

For each input presentation with $K$ active templates, Fuzzy ART computes $K$ choice activations via $2d$-element componentwise minima and $L_1$ sums in $cal(O)(K d)$ operations, sorts $K$ scalars in $cal(O)(K log K)$ time, and tests vigilance until the first resonant match in $cal(O)(r d)$ operations ($r <= K$ resets). Across $T$ presentations, total time is $cal(O)(T K (d + log K))$, simplifying to $cal(O)(T K d)$ whenever $d >= log K$. Because Fuzzy ART relies purely on additions, comparisons, and scalar divisions—requiring no square roots, exponentials, or matrix inversions—its per-operation constant factor is exceptionally small and easily vectorized or synthesized in hardware.

Memory complexity is $cal(O)(K d)$, storing only the $K$ $2d$-dimensional templates without retaining past observations. When $rho$ bounds $K$ independently of stream length $n$, Fuzzy ART is a true single-pass $cal(O)(1)$-memory streaming clusterer.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T K (d + log K))$; arithmetic-light (min/add only)],
    [Space Complexity], [$cal(O)(K d)$ templates; stream-length independent when $K$ bounded],
    [Cluster Geometry], [Axis-aligned $L_1$-bounded hyperrectangles $R_j subset [0, 1]^d$],
    [Hyperparameters], [Vigilance $rho in [0, 1]$, choice $alpha > 0$, learning rate $beta in (0, 1]$],
    [Invariance], [Axis-permutation & reflection invariant; not rotation-invariant],
    [Stability], [Zero catastrophic forgetting; monotone template convergence],
    [Outlier Robustness], [High against prototype drift (outliers spawn isolated categories)],
  ),
  caption: [Operational and geometric profile of ART / Fuzzy ART (Alg. 178).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Online Intrusion Detection and Zero-Day Novelty Discovery:* Network-security monitors must recognize known attack families without forgetting historical signatures while immediately flagging unseen zero-day traffic. ART's vigilance gate mathematically separates known resonant traffic ($M_j >= rho$) from novel anomalies that commit isolated new categories without corrupting existing prototypes @xu2005survey.
+ *Aerospace and Industrial Fault Diagnostics:* Flight-control and turbine diagnostic systems operate under strict safety certification forbidding catastrophic forgetting of rare failure modes during lifelong adaptation. Fuzzy ART's coordinatewise monotonic box bounds (Theorem II.A) guarantee that previously certified fault templates remain intact as new normal operating regimes are learned.
+ *Continual Customer-Behavior Segmentation:* Retail and financial transaction streams introduce emerging consumer cohorts over time. Fuzzy ART adapts in single-pass streaming mode without retraining on historical logs, while axis-aligned box boundaries $[u_(j r), v_(j r)]$ translate directly into human-auditable IF-THEN business rules per feature.
+ *Low-Power Embedded and FPGA Pattern Classification:* Edge sensors and radar emitter recognizers face strict latency and power envelopes. Because Fuzzy ART requires only min-comparisons and $L_1$ summations on fixed-point $[0, 1]$ features, it maps directly onto parallel hardware comparators at microsecond latencies.
+ *Continuous Biomedical Telemetry and Drift-Resilient ECG Monitoring:* Wearable cardiac monitors must adapt to patient-specific baseline shifts while never overwriting rare arrhythmia templates. Vigilance-controlled commitment isolates transient arrhythmias into dedicated templates while allowing normal sinus templates to track gradual physiological drift via $beta < 1$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Fuzzy ART exhibits five documented failure modes: (i) *high-dimensional vigilance collapse* (Proposition II.D), where $M(bold(I)_a, bold(I)_b) -> 2 / 3$ makes tuning $rho$ knife-edge unstable; (ii) *presentation-order sensitivity*, where early outlier or boundary points dictate box placement; (iii) *category proliferation* when $rho -> 1$ or when complement coding is omitted (Theorem II.C); (iv) *min-max normalization fragility*, since a single extreme outlier compresses all normal coordinates into a tiny sub-interval of $[0, 1]^d$, artificially inflating $M_j$; and (v) *axis-alignment mismatch*, where diagonal ellipsoids or curved manifolds require staircased collections of boxes.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Lifelong / Continual Streams], [Excellent: strict stability-plasticity], [Order bias; ensemble or slow $beta$],
    [Rule-Based Tabular Data], [Strong: axis-aligned boxes $[bold(u)_j, bold(v)_j]$], [Sensitive to min-max outliers],
    [High-Dim. Continuous ($d >> 10$)], [Poor: $L_1$ match concentrates at $2 / 3$], [Use Neural Gas (Alg. 179) or spherical $k$-means],
    [Curved / Non-Convex Shapes], [Weak: staircased box tiling], [Use GNG (Alg. 177) or DBSCAN @ester1996dbscan],
    [Binary Feature Patterns], [Excellent: exact ART1 Boolean intersection], [Tune $rho$ to avoid proliferation],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for ART / Fuzzy ART.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Exact Reduction of Fuzzy ART to Binary ART1 and Comparison to Metric Thresholding", origin: "Original Derivation")[
  When raw inputs are restricted to binary vectors $bold(a) in {0, 1}^d$ (without complement coding) and $beta = 1$, componentwise minimum $bold(a) and bold(w)_j$ coincides identically with bitwise logical intersection $bold(a) inter bold(w)_j$, $norm(bold(a))_1 = |bold(a)|$ counts active bits, and Definition I's choice, match, and fast-learning rules reduce identically to Carpenter-Grossberg ART1 @carpenter1987art. Furthermore, under complement coding, the vigilance test $M_j (bold(I)) >= rho$ is algebraically equivalent to the $L_1$ bounding-box growth constraint $|"box"(R_j union {bold(a)})|_1 <= d(1 - rho)$.
]

#proof[
  For $u, v in {0, 1}$, $min(u, v) = u dot v = bold(1)_(u = 1 "and" v = 1)$, so $bold(a) and bold(w)_j$ equals the bitwise AND of binary templates, and its $L_1$ norm equals the Hamming weight $|bold(a) inter bold(w)_j|$. Substituting into $T_j$ and $M_j$ yields the exact ART1 equations of @carpenter1987art. Under complement coding, Theorem II.A proved $norm(bold(I) and bold(w)_j)_1 = d - |"box"(R_j union {bold(a)})|_1$, whence $M_j (bold(I)) >= rho <=> d - |"box"(R_j union {bold(a)})|_1 >= rho d <=> |"box"(R_j union {bold(a)})|_1 <= d(1 - rho)$.
]

Structurally, Fuzzy ART replaces the Euclidean spherical Cluster Feature balls of BIRCH @zhang1996birch with $L_1$-perimeter-constrained hyperrectangles, and contrasts sharply with SOM (Alg. 176), GNG (Alg. 177), and Neural Gas (Alg. 179), which rely on unconstrained $L_2$ prototype pull and therefore lack hard vigilance guarantees against catastrophic forgetting.

*ClusterOverview Benchmark Analysis:* On the 4-dataset ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), ART achieves Iris Sil = $0.357$, ARI = $0.793$; Customers Sil = $0.281$; Multi-Shape Sil = $0.505$ (ARI = N/A); and High-D Sil = $-0.222$ (ARI = N/A), at a mean runtime of $1.700$ ms. Every single metric is illuminated by our formal theory. First, ART attains the *highest Iris ARI* ($0.793$) among all four neural algorithms (surpassing GNG $0.712$, Neural Gas $0.583$, and SOM $0.562$), because the three Iris species align remarkably well with axis-aligned petal length and petal width intervals—precisely Fuzzy ART's hyperrectangle inductive bias (Theorem II.A). Second, Multi-Shape Silhouette is lower ($0.505$ vs. $0.718$ for Neural Gas) because curved shapes force staircased axis-aligned partitions. Third, High-D Silhouette collapses to $-0.222$, providing striking empirical verification of Proposition II.D: in high dimensions, complement-coded match scores concentrate tightly around $2 / 3$, causing hyperrectangles to overlap across Euclidean cluster boundaries and yield negative Silhouette margins.
