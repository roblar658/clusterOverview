#import "../template/ieee.typ": *

#alg_entry(
  id: 169,
  name: "TRACLUS (Trajectory Clustering: Partition-and-Group)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2007",
  paradigm: "MDL Sub-Trajectory Partitioning & Line-Segment Density Clustering",
  complexity: $O(M^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Spatio-temporal tracking systems—cyclone telemetry, GPS vehicle fleets, wildlife collars, and maritime AIS transponders—produce ordered polylines in $bb(R)^d$ ($d in \{2, 3\}$). Whole-trajectory metrics such as Dynamic Time Warping @sakoe1978dtw compare paths end-to-end and miss *common sub-trajectories* when objects share an intermediate corridor yet diverge at their origins or destinations. Lee, Han, and Whang @lee2007traclus introduced TRACLUS, a three-phase *partition-and-group* framework that partitions trajectories into directed segments via Minimum Description Length (MDL), clusters segments via a density-based DBSCAN adaptation, and synthesizes sweep-line representative trajectories.

Let $cal(D) = \{"TR"_1, dots, "TR"_m\}$ be a database of trajectories $"TR"_i = chevron.l bold(p)_1, dots, bold(p)_(n_i) chevron.r subset bb(R)^d$. Selecting characteristic indices $1 = c_1 < dots < c_q = n_i$ splits $"TR"_i$ into directed segments $L_r = overline(bold(p)_(c_r) bold(p)_(c_(r+1)))$. Let $cal(L) = \{L_1, dots, L_M\}$ be the pooled set of $M$ segments across $cal(D)$ with source-trajectory mapping $"tid"(L_u) in \{1, dots, m\}$.

To compare directed segments $L_i = overline(bold(s)_i bold(e)_i)$ and $L_j = overline(bold(s)_j bold(e)_j)$, adopt the convention $norm(L_i)_2 >= norm(L_j)_2$ (swapping if needed) and orthogonally project $bold(s)_j, bold(e)_j$ onto the affine line of $L_i$ to get $bold(p)_s, bold(p)_e$. Set perpendicular offsets $l_(perp 1) = norm(bold(s)_j - bold(p)_s)_2$, $l_(perp 2) = norm(bold(e)_j - bold(p)_e)_2$, parallel offsets $l_(parallel 1) = min(norm(bold(p)_s - bold(s)_i)_2, norm(bold(p)_s - bold(e)_i)_2)$, $l_(parallel 2) = min(norm(bold(p)_e - bold(s)_i)_2, norm(bold(p)_e - bold(e)_i)_2)$, and directed heading angle $theta in [0, pi]$ between $bold(e)_i - bold(s)_i$ and $bold(e)_j - bold(s)_j$.

#definition(title: "Composite Line-Segment Dissimilarity and MDL Partitioning Cost")[
  For $norm(L_i)_2 >= norm(L_j)_2$, the three TRACLUS component dissimilarities are @lee2007traclus:
  $ d_perp (L_i, L_j) = (l_(perp 1)^2 + l_(perp 2)^2) / (l_(perp 1) + l_(perp 2)), quad d_parallel (L_i, L_j) = min(l_(parallel 1), l_(parallel 2)), $
  $ d_theta (L_i, L_j) = cases(norm(L_j)_2 sin theta & "if" 0 <= theta < pi \/ 2, norm(L_j)_2 & "if" pi \/ 2 <= theta <= pi,) $
  with $d_perp = 0$ if $l_(perp 1) = l_(perp 2) = 0$. Given weights $w_perp, w_parallel, w_theta > 0$ (default $1$), total dissimilarity is $"dist"(L_i, L_j) = w_perp d_perp + w_parallel d_parallel + w_theta d_theta$. In Phase 1, approximating $chevron.l bold(p)_a, dots, bold(p)_b chevron.r$ by chord $L_(a b) = overline(bold(p)_a bold(p)_b)$ versus preserving elementary edges $e_k = overline(bold(p)_k bold(p)_(k+1))$ incurs MDL bit costs
  $ "MDL"_"par" (a, b) &= log_2 norm(L_(a b))_2 + sum_(k=a)^(b-1) (log_2 d_perp (L_(a b), e_k) + log_2 d_theta (L_(a b), e_k)), \
    "MDL"_"nopar" (a, b) &= sum_(k=a)^(b-1) log_2 norm(e_k)_2, $
  where arguments inside $log_2$ are floored at $1$ to keep bit lengths non-negative.
]

In Phase 2, TRACLUS clusters segments using $N_epsilon (L_i) = \{L_j in cal(L) : "dist"(L_i, L_j) <= epsilon\}$ and threshold `"MinLns"`, retaining density-connected segment clusters $C$ iff distinct trajectory count $|"PTR"(C)| = |\{"tid"(L) : L in C\}| >= "MinLns"$.

#algbox(
  title: "TRACLUS: Partition-and-Group Trajectory Clustering",
  inputs: [Trajectories $cal(D) = \{"TR"_1, dots, "TR"_m\}$, radius $epsilon > 0$, threshold $"MinLns"$, weights $(w_perp, w_parallel, w_theta)$, smoothing step $gamma > 0$.],
  outputs: [Sub-trajectory clusters $\{C_r\}$ and representative trajectories $\{R(C_r)\}$.],
)[
  + *(Phase 1: Greedy MDL Partitioning)* For each $"TR"_i$: start $a <- 1$; advance $b <- a + 1, dots, n_i$ while $"MDL"_"par" (a, b) <= "MDL"_"nopar" (a, b)$; upon first violation, emit segment $overline(bold(p)_a bold(p)_(b-1))$ into $cal(L)$ and reset $a <- b - 1$.
  + *(Phase 2: Segment Density Grouping)* Execute DBSCAN over $(cal(L), "dist")$ with $(epsilon, "MinLns")$ to form maximal density-connected segment clusters $\{C_r\}$.
  + *(Cardinality Filter)* Discard any $C_r$ with $|"PTR"(C_r)| < "MinLns"$ to prevent single looping trajectories from forming clusters.
  + *(Phase 3: Representative Synthesis)* For each surviving $C_r$: align the $x'$-axis with average direction $bar(bold(v)) = |C_r|^(-1) sum_(L in C_r) (bold(e)_L - bold(s)_L)$, sort the $2 |C_r|$ projected endpoints, and sweep a perpendicular hyperplane; whenever active segment count $>= "MinLns"$ and spacing from the prior vertex $>= gamma$, append the back-rotated mean transverse point to $R(C_r)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

TRACLUS's composite segment dissimilarity exhibits sharp non-metric pathologies that invalidate metric trees, yet satisfies an exact Euclidean MBR lower bound certifying R-tree pruning alongside tight Lehmer-mean, MDL-optimality, and sweep-line guarantees.

#subsec("A", "Violation of Metric Space Axioms and Boundary Discontinuity")

#theorem(title: "Violation of Identity, Triangle Inequality, Symmetry, and Continuity", origin: "Original Counterexample")[
  With unit weights $w_perp = w_parallel = w_theta = 1$, $"dist"(dot.c, dot.c)$ on directed segments in $bb(R)^2$ is neither a metric nor continuous:
  + *(Non-Identity)* Distinct segments $L_1 = overline((0,0)(4,0))$ and $L_2 = overline((0,0)(2,0))$ satisfy $"dist"(L_1, L_2) = 0$.
  + *(Infinite-Ratio Triangle Violation)* Collinear segments $L_1 = overline((0,0)(10,0))$, $L_2 = overline((0,0)(1,0))$, and $L_3 = overline((9,0)(10,0))$ satisfy $"dist"(L_2, L_3) = 8 > 0 = "dist"(L_2, L_1) + "dist"(L_1, L_3)$.
  + *(Equal-Length Asymmetry & Jump Discontinuity)* For equal-length orthogonal segments $L_a = overline((0,0)(2,0))$ and $L_b = overline((0,1)(0,3))$, projecting $L_b$ onto $L_a$ gives $"dist" = 4.5$, whereas projecting $L_a$ onto $L_b$ gives $"dist" = 5.0$, creating a $0.5$ jump discontinuity across $norm(L_a)_2 = norm(L_b)_2$.
]
#proof[
  *(Part 1)* Since $norm(L_1)_2 = 4 > 2 = norm(L_2)_2$, projecting $L_2$ onto the $x$-axis gives $l_(perp 1) = l_(perp 2) = 0 => d_perp = 0$, $l_(parallel 1) = 0 => d_parallel = 0$, and $theta = 0 => d_theta = 0$, so $"dist"(L_1, L_2) = 0$ despite $L_1 != L_2$.

  *(Part 2)* All pairs have $d_perp = d_theta = 0$. Comparing $L_2$ and $L_3$ individually to $L_1$ ($norm(L_1)_2 = 10 > 1$) yields $l_(parallel 1) = 0$ at $(0,0)$ and $l_(parallel 2) = 0$ at $(10,0)$, giving $"dist"(L_2, L_1) = "dist"(L_1, L_3) = 0$. Yet comparing $L_2$ and $L_3$ directly yields $l_(parallel 1) = 8$, $l_(parallel 2) = 9$, and $"dist"(L_2, L_3) = 8 > 0 + 0$.

  *(Part 3)* Let $L_a^delta = overline((0,0)(2 + delta, 0))$ and $L_b = overline((0,1)(0,3))$ ($theta = pi \/ 2$). For $delta > 0$, $L_b$ projects onto the $x$-axis at $(0,0)$: $l_(perp 1) = 1, l_(perp 2) = 3 => d_perp = 10 \/ 4 = 2.5$; $d_parallel = 0$; $d_theta = 2$, so $lim_(delta -> 0^+) "dist" = 4.5$. For $delta < 0$, $L_a^delta$ projects onto the $y$-axis at $(0,0)$: $l_(perp 1) = 0, l_(perp 2) = 2 + delta => d_perp = 2 + delta$; projection $(0,0)$ lies distance $1$ from $(0,1)$, so $d_parallel = 1$; $d_theta = 2 + delta$, giving $lim_(delta -> 0^-) "dist" = 5.0 != 4.5$.
]

Theorem 1 proves that triangle-inequality metric indices (Ball-Trees, VP-Trees) are *unsound* for TRACLUS and can silently miss true neighbors.

#subsec("B", "Lehmer-Mean Bounds and Exact Soundness of R-Tree MBR Pruning")

#lemma(title: "Tight Lehmer-Mean Constants and Exact Soundness of R-Tree MBR Pruning", origin: "Original Proof")[
  For directed segments $L_i, L_j subset bb(R)^d$ and $w_"min" = min(w_perp, w_parallel, w_theta) > 0$:
  + *(Lehmer-Mean Dominance)* With $l_"max" = max(l_(perp 1), l_(perp 2))$,
    $ (l_(perp 1) + l_(perp 2)) / 2 <= d_perp (L_i, L_j) <= l_"max", quad d_perp (L_i, L_j) >= 2 (sqrt(2) - 1) l_"max", $
    where $2 (sqrt(2) - 1) approx 0.8284$ is tight.
  + *(Rigid Equivariance)* $"dist"$ is $E(d)$-invariant under $R bold(p) + bold(t)$ ($R in O(d)$) and positively $1$-homogeneous: $"dist"(alpha L_i, alpha L_j) = alpha "dist"(L_i, L_j)$ ($alpha > 0$).
  + *(Sound MBR Lower Bound)* The minimum Euclidean distance $d_"MBR" (L_i, L_j)$ between axis-aligned bounding boxes $"MBR"(L_i), "MBR"(L_j)$ satisfies
    $ d_"MBR" (L_i, L_j) <= d_"seg" (L_i, L_j) <= "dist"(L_i, L_j) / w_"min", $
    guaranteeing zero false dismissals for R-tree queries expanded by $epsilon \/ w_"min"$.
]
#proof[
  *(Part 1)* For $l_"max" > 0$, put $r = min(l_(perp 1), l_(perp 2)) \/ l_"max" in [0, 1]$, so $d_perp \/ l_"max" = g(r) = (1 + r^2) \/ (1 + r) <= 1$. Also $2(1 + r^2) - (1 + r)^2 = (1 - r)^2 >= 0 => d_perp >= (l_(perp 1) + l_(perp 2)) \/ 2$. Setting $g'(r) = (r^2 + 2 r - 1) \/ (1 + r)^2 = 0$ yields unique interior minimum $r^* = sqrt(2) - 1$ with $g(r^*) = 2 (sqrt(2) - 1)$.

  *(Part 2)* Orthogonal projections, norms, and inner-product angles are preserved by $E(d)$ and scale linearly by $alpha > 0$.

  *(Part 3)* Assume $norm(L_i)_2 >= norm(L_j)_2$. Projecting $bold(e)_j - bold(s)_j$ orthogonally to $L_i$ gives $norm((bold(e)_j - bold(p)_e) - (bold(s)_j - bold(p)_s))_2 = norm(L_j)_2 sin theta <= d_theta (L_i, L_j)$. By the reverse triangle inequality, $|l_(perp 1) - l_(perp 2)| <= d_theta$, and since $min(l_(perp 1), l_(perp 2)) <= d_perp$, we get $l_"max" <= d_perp + d_theta$. Let $bold(q) in \{bold(s)_j, bold(e)_j\}$ attain $d_parallel$ with projection $bold(p)_q$ at distance $d_parallel$ from endpoint $bold(u) in \{bold(s)_i, bold(e)_i\}$. By orthogonality,
  $ d_"MBR" (L_i, L_j) &<= d_"seg" (L_i, L_j) <= norm(bold(q) - bold(u))_2 \
    &= sqrt(norm(bold(q) - bold(p)_q)_2^2 + d_parallel^2) <= l_"max" + d_parallel \
    &<= d_perp + d_theta + d_parallel <= "dist"(L_i, L_j) / w_"min". $
  Hence $"dist"(L_i, L_j) <= epsilon => d_"MBR" (L_i, L_j) <= epsilon \/ w_"min"$, certifying exact R-tree pruning.
]

#subsec("C", "Global Dynamic Programming vs. Greedy MDL Optimality Gap")

#theorem(title: "Exact DAG Dynamic Program and Greedy MDL Suboptimality", origin: "Original Derivation")[
  For trajectory $"TR" = chevron.l bold(p)_1, dots, bold(p)_n chevron.r$:
  + Global MDL minimization over all $2^(n - 2)$ characteristic-point subsets reduces to a shortest path on the DAG $(\{1, dots, n\}, \{(a, b) : a < b\})$ with arc cost $w_"MDL" (a, b) = "MDL"_"par" (a, b)$ ($w_"MDL" (a, a+1) = log_2 norm(e_a)_2$), solvable exactly in $O(n^3)$ time and $O(n)$ space.
  + Phase 1's greedy scan runs in $O(n bar(ell))$ time but is *strictly suboptimal*: on the 4-point integer polyline $(0,0), (16,0), (32,3), (160,27)$, greedy outputs a single chord costing $13.99$ bits versus the exact DP optimum of $11.20$ bits ($Delta_"MDL" = 2.79$ bits).
]
#proof[
  *(Part 1)* Total MDL cost decomposes additively over consecutive intervals $[c_r, c_(r+1)]$. The recurrence $"DP"[1] = 0$, $"DP"[b] = min_(1 <= a < b) ("DP"[a] + w_"MDL" (a, b))$ computes the exact global optimum across $O(n^2)$ arcs, each evaluated in $O(n)$ time ($O(n^3)$ total).

  *(Part 2)* Evaluate $bold(p)_1 = (0,0), bold(p)_2 = (16,0), bold(p)_3 = (32,3), bold(p)_4 = (160,27)$, where $e_1 = (16,0)$ ($norm(e_1) = 16$), $e_2 = (16,3)$ ($norm(e_2) = sqrt(265) approx 16.28$), and $e_3 = (128,24) = 8 e_2$ is strictly collinear with $e_2$. Starting at $a = 1$:
  - At $b = 3$: chord $L_(1 3) = (32,3)$ has length $sqrt(1033) approx 32.14$ ($log_2 = 5.006$). Both $d_perp$ and $d_theta$ for $e_1, e_2$ equal $48 \/ sqrt(1033) approx 1.493$ ($log_2 = 0.579$), so $"MDL"_"par" (1, 3) = 5.006 + 4(0.579) = 7.320 <= 8.025 = "MDL"_"nopar" (1, 3)$.
  - At $b = 4$: chord $L_(1 4) = (160,27)$ has length $sqrt(26329) approx 162.26$ ($log_2 = 7.342$). Perpendicular distances of $bold(p)_2, bold(p)_3$ to $L_(1 4)$ are $432 \/ 162.26 = 2.662$ and $384 \/ 162.26 = 2.367$. Computing exact components gives $d_perp (e_1) = d_theta (e_1) = 2.662$, $d_perp (e_2) = 2.523$, $d_theta (e_2) = 0.296 < 1$ ($log_2 = 0$), and $d_perp (e_3) = d_theta (e_3) = 2.367$, summing to $"MDL"_"par" (1, 4) = 7.342 + 1.413 + 1.413 + 1.335 + 0 + 1.243 + 1.243 = 13.989 <= 15.050 = "MDL"_"nopar" (1, 4)$.
  Thus greedy accepts $(1, 4)$ at $13.989$ bits. By contrast, splitting at $bold(p)_2$ merges collinear $(2, 4)$ at zero distortion ($norm(L_(2 4)) = 9 sqrt(265) approx 146.51$), achieving DP cost $log_2 16 + log_2 (9 sqrt(265)) = 4.000 + 7.195 = 11.195$ bits—saving $2.794$ bits.
]

#subsec("D", "Sweep-Line Representative Trajectory Optimality")

#theorem(title: "Monotonicity, Non-Self-Intersection, and L2 Optimality of R(C)", origin: "Original Proof")[
  Let cluster $C subset.eq cal(L)$ ($s = |C|$) have average direction $bar(bold(v)) != bold(0)$, and assume every $L in C$ satisfies $chevron.l bold(e)_L - bold(s)_L, bar(bold(v)) chevron.r > 0$. Then in the rotated frame aligned with $hat(bold(e))_(x') = bar(bold(v)) \/ norm(bar(bold(v)))_2$: (i) every segment is strictly $x'$-monotone ($x'_s (L) < x'_e (L)$) and meets each sweep line $x' = t$ at most once at $bold(y)'_L (t) in bb(R)^(d-1)$; (ii) the emitted vertex $bar(bold(y))'(t) = |A(t)|^(-1) sum_(L in A(t)) bold(y)'_L (t)$ uniquely minimizes $sum_(L in A(t)) norm(bold(y)' - bold(y)'_L (t))_2^2$; and (iii) $R(C)$ is strictly non-self-intersecting and computable in $O(s log s)$ time.
]
#proof[
  Positive inner product gives $x'_e (L) - x'_s (L) > 0$, making each segment an affine graph $bold(y)'_L (t) = bold(a)_L + t bold(b)_L$ over $[x'_s (L), x'_e (L)]$. Strict convexity of $F(bold(y)') = sum_(L in A(t)) norm(bold(y)' - bold(y)'_L (t))_2^2$ (Hessian $2 |A(t)| I_(d-1) succ 0$) makes $bar(bold(y))'(t)$ the unique global minimizer. Because sweep coordinates strictly increase ($t_1 < dots < t_r$), $R(C)$ is a single-valued graph over $x'$ and hence free of self-intersections under rigid back-rotation. Sorting $2 s$ events takes $O(s log s)$ time, and maintaining running affine sums $sum_(L in A(t)) bold(a)_L$ and $sum_(L in A(t)) bold(b)_L$ in $O(1)$ per insertion/deletion evaluates each sweep mean in $O(1)$, yielding $O(s log s)$ total time.
]

#subsec("E", "Entropy Heuristic and Total Framework Complexity")

#theorem(title: "Entropy-Guided Radius Selection and Partition-and-Group Complexity", origin: "Classical — Cited")[
  Let $N_"pts" = sum_(i=1)^m n_i$ be total vertices and $M = |cal(L)| <= N_"pts"$ the partitioned segments @lee2007traclus. Minimizing neighborhood entropy $H(epsilon) = - sum_(i=1)^M p_epsilon (L_i) log_2 p_epsilon (L_i)$ with $p_epsilon (L_i) = |N_epsilon (L_i)| \/ sum_j |N_epsilon (L_j)|$ (and $"MinLns" = ceil("avg"_i |N_epsilon (L_i)|) + 1$) selects the scale maximizing corridor-vs-noise contrast. Total runtime is $O(N_"pts" + M^2)$ naively and $O(N_"pts" + M log M)$ with R-tree MBR pruning (Lemma 1) under bounded overlap.
]
#proof[
  As $epsilon -> 0$ or $epsilon -> infinity$, $p_epsilon$ approaches uniform $1 \/ M$ maximizing $H(epsilon) = log_2 M$; intermediate scales concentrating mass on dense corridors sharply reduce $H(epsilon)$ @lee2007traclus. Phase 1 runs in $O(N_"pts")$ time, Phase 2 executes $M$ range queries in $O(M^2)$ or $O(M log M)$ via Lemma 1's exact MBR filter, and Phase 3 runs in $O(sum_r |C_r| log |C_r|) = O(M log M)$ by Theorem 3.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Phase 1 MDL compression decouples clustering cost from GPS sampling rate, routinely compressing $n_i approx 10^3$ raw vertices into $O(10)$ segments ($M << N_"pts"$). Phase 2 requires $O(M^2)$ time naively and $O(M log M)$ with R-tree MBR pruning certified by Lemma 1(3), using $O(N_"pts" + M)$ space while avoiding the unsound metric-tree indexing exposed in Theorem 1. Geometrically, operating in $(2 d)$-dimensional directed segment space isolates co-directional movement corridors regardless of divergent trajectory endpoints.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N_"pts" + M^2)$ naive; $O(N_"pts" + M log M)$ with R-tree MBR pruning (Lemma 1)],
    [Space Complexity], [$O(N_"pts" + M)$ for trajectories, $M$ segments, and DBSCAN structures],
    [Dissimilarity Space], [Non-metric directed segment score $w_perp d_perp + w_parallel d_parallel + w_theta d_theta$],
    [Metric Pathologies], [Violates identity, triangle inequality, equal-length symmetry, continuity (Thm. 1)],
    [Inductive Bias], [Co-directional sub-trajectory corridors supported by $>= "MinLns"$ distinct tracks],
    [Key Hyperparameters], [Radius $epsilon$, threshold $"MinLns"$, weights $(w_perp, w_parallel, w_theta)$, sweep step $gamma$],
    [Invariance], [Exact $E(d)$ rigid invariance; positively $1$-homogeneous under scaling],
    [Output Summary], [Segment clusters plus $L_2$-optimal sweep-line representative polylines $R(C)$],
  ),
  caption: [Operational and theoretical profile of TRACLUS (Alg. 169).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Meteorological Hurricane & Cyclone Track Analysis:* Atlantic storms originate across wide longitudes yet recurve along shared Gulf Stream steering corridors; MDL partitioning isolates turning points and clusters shared recurvature sub-tracks missed by whole-track metrics @lee2007traclus.
+ *Wildlife Migration & Ecological Telemetry:* GPS-collared animals exhibit local foraging loops interrupted by shared migration bottlenecks; the distinct-ID filter $|"PTR"(C)| >= "MinLns"$ blocks single looping individuals from forming false corridors while extracting multi-animal highways.
+ *Maritime AIS & Aviation Corridor Monitoring:* Because $d_theta$ penalizes opposing headings ($theta >= pi \/ 2 => d_theta = norm(L_j)_2$), TRACLUS cleanly separates inbound and outbound shipping lanes sharing the same strait and flags off-corridor deviations.
+ *Urban Mobility & Ride-Hailing Route Mining:* Phase 1 MDL strips collinear GPS points along straight avenues while retaining turns, enabling $O(M log M)$ arterial corridor discovery across large taxi fleets.
+ *Computer Vision Surveillance & Crowd Flow Parsing:* Sweep-line representative polylines (Theorem 3) extract clean, non-self-intersecting dominant pedestrian flow axes across crowded plazas for anomaly and counter-flow detection.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure modes bound deployment: (1) *Static Non-Trajectory Data*, where synthetic point orderings produce erratic segment headings $d_theta$; (2) *Invalid Metric Trees* causing silent false dismissals (Theorem 1); (3) *Multi-Scale Heterogeneity* where one global $epsilon$ cannot span dense junctions and wide highways; and (4) *High-Curvature Reversals* within a single cluster violating Theorem 3's $< pi \/ 2$ alignment assumption.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Cyclone / Wildlife / AIS Tracks], [Excellent: sub-trajectory MDL + directional corridors], [Use R-tree (Lemma 1), never metric trees],
    [Whole-Sequence Shape Matching], [Low: partitioning discards global temporal alignment], [Use DTW @sakoe1978dtw, DBA @petitjean2011dba, or k-Shape @paparrizos2015kshape],
    [Spatio-Temporal Event Streams], [Moderate: discards timestamps once ordered into polylines], [Use ST-DBSCAN @birant2007stdbscan for time sync],
    [Multi-Density Road Networks], [Moderate: global $epsilon$ misses sparse rural corridors], [Swap Phase 2 to HDBSCAN @campello2013hdbscan on segments],
    [Static / High-D Tabular Data], [Poor: artificial segments distort isotropic proximity], [Use Watershed (Alg. 167), DBSCAN, or GMM],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for TRACLUS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

TRACLUS merges MDL curve encoding @wallace1968snob with Generalized DBSCAN (GDBSCAN) @sander1998gdbscan, contrasting sharply with whole-sequence elastic alignment @sakoe1978dtw @petitjean2011dba @cuturi2017softdtw.

#proposition(title: "GDBSCAN Reduction and Incomparability with Whole-Trajectory DTW", origin: "Original Proof")[
  Define predicate $"NPred"(L_i, L_j) <==> "dist"(L_i, L_j) <= epsilon$ with unit weight $w(L) = 1$ (under deterministic equal-length tie-breaking). Then: (i) Phase 2 is an exact instance of GDBSCAN @sander1998gdbscan augmented by post-cluster filter $|"PTR"(C)| >= "MinLns"$; and (ii) for general multi-segment tracks, TRACLUS sub-trajectory similarity is strictly incomparable to whole-trajectory DTW @sakoe1978dtw (tracks of arbitrarily large DTW distance can share zero-distance TRACLUS segments).
]
#proof[
  With deterministic tie-breaking, `"NPred"` is reflexive and symmetric, so $N_epsilon (L_i)$ and $|N_epsilon (L_i)| >= "MinLns"$ directly instantiate GDBSCAN @sander1998gdbscan, followed by $|"PTR"(C)| >= "MinLns"$. For (ii), two trajectories sharing $overline((0,0)(1,0))$ while diverging to $(1, R)$ and $(1, -R)$ achieve $"dist" = 0$ on their shared sub-segment yet incur whole-track DTW distance $Omega(R) -> infinity$.
]

On the *ClusterOverview* four-dataset benchmark (where static point clouds are mapped to synthetic local segments), TRACLUS records *Iris* Silhouette @rousseeuw1987silhouettes $0.391$ and ARI @hubert1985comparing $0.467$, *Customers* Silhouette $0.321$, *Multi-Shape* Silhouette $-0.454$ (ARI N/A), and *High-D* Silhouette $-0.107$ (ARI N/A), with mean runtime $1.644$ ms—slowest among Algs. 166--169 due to segment extraction and $O(M^2)$ composite scoring. The *negative* Silhouettes on *Multi-Shape* ($-0.454$) and *High-D* ($-0.107$) cleanly expose the mismatch between directed segment dissimilarity and isotropic Euclidean point geometry: on unordered curved clouds, artificial tangents assign large $d_theta$ penalties between adjacent points while Theorem 1's $d_parallel$ triangle violation chains distant collinear segments, inverting $a(i) > b(i)$ ($s(i) < 0$). Conversely, monotonic low-D trends in *Iris* ($0.467$ ARI) and *Customers* ($0.321$ Sil) retain partial segment alignment, confirming that TRACLUS excels strictly on genuine sequential trajectory manifolds.
