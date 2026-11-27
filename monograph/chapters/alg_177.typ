#import "../template/ieee.typ": *

#alg_entry(
  id: 177,
  name: "Growing Neural Gas (GNG)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1995",
  paradigm: "Incremental Topology-Learning Neural Graph",
  complexity: $cal(O)(T M_"max" d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Growing Neural Gas (GNG), introduced by Fritzke @fritzke1995gng, eliminates two structural restrictions of fixed-lattice Self-Organizing Maps (Alg. 176): a pre-specified neuron count $M$ and a rigid low-dimensional grid topology. Instead, GNG maintains a dynamically evolving undirected graph $G_t = (V_t, E_t)$ embedded in $bb(R)^d$, incrementally inserting prototypes where local quantization distortion is largest and learning neighborhood edges directly from data via Competitive Hebbian Learning (CHL) @martinetz1993neuralgas.

Each vertex $i in V_t$ holds a prototype vector $bold(w)_i in bb(R)^d$ and a nonnegative cumulative local error variable $e_i >= 0$. Each undirected edge $e = {i, j} in E_t$ carries an integer age $a_(i j) in bb(N)_0$ measuring how many winner events have elapsed at its endpoints without refreshing that connection. Let $N_t (i) := {j in V_t : {i, j} in E_t}$ denote the 1-hop topological neighborhood of vertex $i$ at step $t$.

#definition(title: "Voronoi Regions and Induced Delaunay Triangulation")[
  Let $P$ be a Borel probability measure on $bb(R)^d$ with support $S := "supp"(P)$ and let $W = {bold(w)_i}_(i in V)$ be finite prototypes in general position. The first-order Voronoi cell of $i in V$ and second-order Voronoi cell of an unordered pair ${i, j} subset V$ are
  $
    V_i (W) &:= {bold(x) in bb(R)^d : norm(bold(x) - bold(w)_i)_2
      <= norm(bold(x) - bold(w)_k)_2, forall k in V}, \
    V_(i j) (W) &:= {bold(x) in bb(R)^d : max(norm(bold(x) - bold(w)_i)_2,
      norm(bold(x) - bold(w)_j)_2) \
      &quad <= norm(bold(x) - bold(w)_k)_2, forall k in V without {i, j}}.
  $
  The *induced Delaunay triangulation* of $W$ restricted to $S$ is the graph $D_P (W) := (V, E_P (W))$ with edge set
  $
    E_P (W) := {{i, j} subset V : V_(i j) (W) inter S != emptyset}.
  $
]

Notice that $E_P (W)$ automatically adapts to the local intrinsic dimension of $S$: on a 1D curve it forms chains, on a 2D surface it forms triangulated meshes, and across disconnected cluster components with empty separating regions $V_(i j) (W) inter S = emptyset$ it forms disconnected subgraphs whose connected components directly recover non-convex clusters.

#algbox(
  title: "Growing Neural Gas (GNG)",
  inputs: [Input stream / dataset $X subset bb(R)^d$; learning rates $0 < epsilon_n < epsilon_b < 1$; max edge age $a_"max" in bb(N)$; insertion cadence $lambda in bb(N)$; error split $alpha in (0, 1)$ and decay $beta in (0, 1)$; node cap $M_"max"$.],
  outputs: [Neural graph $G_T = (V_T, E_T)$, prototypes ${bold(w)_i}_(i in V_T)$, and cluster labels (connected components or prototype assignments).],
)[
  + *Initialize:* Set $V <- {1, 2}$ with $bold(w)_1, bold(w)_2$ sampled from $X$, $E <- emptyset$, $e_1 <- 0$, $e_2 <- 0$.
  + *For $t = 1, dots, T$ do:*
    - Sample $bold(x)_t in X$; find winner $s_1$ and runner-up $s_2$:
      $
        s_1 &= arg min_(i in V) norm(bold(x)_t - bold(w)_i)_2^2, \
        s_2 &= arg min_(i in V without {s_1}) norm(bold(x)_t - bold(w)_i)_2^2.
      $
    - Age incident edges: $a_(s_1 j) <- a_(s_1 j) + 1$ for all $j in N(s_1)$.
    - Accumulate error: $e_(s_1) <- e_(s_1) + norm(bold(x)_t - bold(w)_(s_1))_2^2$.
    - Adapt winner and 1-hop neighbors:
      $
        bold(w)_(s_1) &<- bold(w)_(s_1) + epsilon_b (bold(x)_t - bold(w)_(s_1)), \
        bold(w)_j &<- bold(w)_j + epsilon_n (bold(x)_t - bold(w)_j),
          quad forall j in N(s_1).
      $
    - *Competitive Hebbian Edge:* Add ${s_1, s_2}$ to $E$ if absent and reset $a_(s_1 s_2) <- 0$.
    - *Pruning:* Delete every edge with $a_(i j) > a_"max"$; delete any vertex $i$ that becomes isolated ($N(i) = emptyset$).
    - *Growth:* If $t equiv 0 ("mod" lambda)$ and $|V| < M_"max"$: select $q = arg max_(i in V) e_i$ and $f = arg max_(j in N(q)) e_j$; insert vertex $r$ at $bold(w)_r = (bold(w)_q + bold(w)_f) / 2$; replace ${q, f}$ by ${q, r}, {r, f}$ of age $0$; rescale $e_q <- alpha e_q$, $e_f <- alpha e_f$, $e_r <- e_q$.
    - *Global Decay:* Set $e_i <- (1 - beta) e_i$ for all $i in V$.
  + *Return:* Graph $G_T$, prototypes $W$, and nearest-prototype or connected-component cluster partition.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four analytical results for GNG: exponential-rate recovery of the induced Delaunay triangulation, an exact renewal-theoretic phase separation between stale and valid edge lifetimes under age pruning, steady-state distortion tracking of the exponentially decayed local error variables, and monotone global distortion reduction upon midpoint node insertion.

#subsec("A", "Induced Delaunay Recovery Under Competitive Hebbian Learning")

#theorem(title: "Exact Recovery of the Induced Delaunay Subgraph", origin: "Classical — Cited")[
  Let $W = {bold(w)_i}_(i in V)$ be fixed prototypes in general position and let $bold(x)_1, dots, bold(x)_T$ be i.i.d. draws from $P$ with support $S$. Assume $p_"min" := min_(e in E_P (W)) P(V_e (W)) > 0$. Under Competitive Hebbian edge creation without premature deletion ($a_"max" = infinity$), the learned edge set $E_T$ satisfies $E_T subset.eq E_P (W)$ almost surely for all $T >= 1$, and
  $
    bb(P)(E_T = E_P (W)) >= 1 - |E_P (W)| (1 - p_"min")^T.
  $
  In particular, $E_T = E_P (W)$ with probability at least $1 - delta$ whenever $T >= (1 / p_"min") ln(|E_P (W)| / delta)$ @martinetz1993neuralgas @fritzke1995gng.
]

#proof[
  At any step $t$, the pair ${s_1, s_2} = {i, j}$ is selected if and only if $bold(x)_t in V_(i j) (W)$ (up to tie boundaries of measure zero under general position). Because $bold(x)_t in S$ almost surely, every inserted edge satisfies $V_(i j) (W) inter S != emptyset$, hence ${i, j} in E_P (W)$. Thus $E_T subset.eq E_P (W)$ almost surely for every $T >= 1$.

  For completeness, an edge $e in E_P (W)$ is absent from $E_T$ if and only if none of the $T$ independent draws falls inside $V_e (W)$, an event of probability $(1 - P(V_e (W)))^T <= (1 - p_"min")^T$. Applying the union bound across all $e in E_P (W)$ yields
  $
    bb(P)(E_T != E_P (W))
      <= sum_(e in E_P (W)) (1 - P(V_e (W)))^T
      <= |E_P (W)| (1 - p_"min")^T.
  $
  Using $1 - p_"min" <= exp(-p_"min")$, setting $|E_P (W)| exp(-p_"min" T) <= delta$ directly yields the sample complexity bound.
]

#subsec("B", "Exact Renewal Lifetime of Valid vs. Stale Edges")

Because GNG continuously moves prototypes, previously valid edges can become geometrically invalid as Voronoi cells shift. The age threshold $a_"max"$ must prune stale edges rapidly while preserving valid manifold edges. We derive the exact expected survival time of any edge under quasi-stationary prototypes.

#theorem(title: "Exact Expected Edge Lifetime and Exponential Phase Separation", origin: "Original Proof")[
  Under quasi-stationary prototypes $W$ and i.i.d. draws $bold(x)_t ~ P$, consider an existing edge $e = {i, j} in E_0$ initialized at age $a_(i j) = 0$. Let $p_i := P(V_i (W))$, $p_j := P(V_j (W))$ with $p_i + p_j > 0$, and $p_(i j) := P(V_(i j) (W))$. Define the conditional refresh probability $rho_e := p_(i j) / (p_i + p_j) in [0, 1)$. Then the expected number of input presentations $T_"del" (e)$ until edge $e$ is pruned ($a_(i j) = a_"max" + 1$) is given exactly by
  $
    bb(E)[T_"del" (e)] = cases(
      (a_"max" + 1) / (p_i + p_j) & "if" p_(i j) = 0 " (stale edge),",
      ((1 - rho_e)^(-(a_"max" + 1)) - 1) / p_(i j)
        & "if" p_(i j) > 0 " (valid edge)."
    )
  $
]

#proof[
  Observe that $a_(i j)$ changes only on steps where $s_1 in {i, j}$, i.e., when $bold(x)_t in V_i (W) union V_j (W)$. Because $V_i (W)$ and $V_j (W)$ intersect only on a hyperplane of $P$-measure zero, each draw triggers an endpoint-hit event independently with probability $p_e := p_i + p_j$. Conditional on an endpoint hit, $a_(i j)$ resets to $0$ if $bold(x)_t in V_(i j) (W)$ (probability $rho_e = p_(i j) / p_e$) and increments by $1$ otherwise (probability $1 - rho_e$).

  Let $m := a_"max" + 1$, and let $K$ be the number of endpoint hits until the absorbing state $a_(i j) = m$ is first reached from state $0$. Denote $u_k := bb(E)[K | a_(i j) = k]$ for $0 <= k <= m$, with $u_m = 0$. One-step Markov conditioning gives the linear system
  $
    u_k = 1 + rho_e u_0 + (1 - rho_e) u_(k+1),
      quad k = 0, 1, dots, m - 1.
  $
  Unrolling this recurrence from $k = 0$ down to $u_m = 0$ yields
  $
    u_0 = (1 + rho_e u_0) sum_(r=0)^(m-1) (1 - rho_e)^r.
  $
  If $p_(i j) = 0$ (so $rho_e = 0$), the sum equals $m$ and $u_0 = m = a_"max" + 1$. If $p_(i j) > 0$ (so $rho_e in (0, 1)$), the geometric series equals $(1 - (1 - rho_e)^m) / rho_e$, giving
  $
    u_0 (1 - (1 - (1 - rho_e)^m))
      = (1 - (1 - rho_e)^m) / rho_e,
  $
  and solving for $u_0$ yields
  $
    u_0 = ((1 - rho_e)^(-m) - 1) / rho_e.
  $
  Finally, the total number of samples $T_"del" (e)$ is the sum of $K$ independent $"Geom"(p_e)$ inter-arrival times between endpoint hits. Because $K$ is a stopping time for the i.i.d. sequence with $bb(E)[K] = u_0 < infinity$, Wald's identity yields $bb(E)[T_"del" (e)] = u_0 / p_e$. Substituting $rho_e p_e = p_(i j)$ gives both branches of the theorem.
]

Theorem II.B explains why GNG's simple age counter succeeds robustly: stale edges ($p_(i j) = 0$) expire in *linear* time $cal(O)(a_"max" / (p_i + p_j))$, whereas valid Delaunay edges with modest refresh ratio $rho_e$ survive for an *exponentially* long lifetime $Omega((1 - rho_e)^(-a_"max"))$.

#subsec("C", "Exponentially Decayed Error Dynamics and Minimax Targeting")

#lemma(title: "Uniform Error Boundedness and Stationary Distortion Tracking", origin: "Original Derivation")[
  Suppose $S = "supp"(P)$ has finite diameter $D := sup_(bold(x), bold(y) in S) norm(bold(x) - bold(y))_2 < infinity$, prototypes are initialized in $"conv"(S)$, and $epsilon_b, epsilon_n, alpha, beta in (0, 1)$. Then:
  + For every node $i in V_t$ and every step $t >= 0$, the accumulated error satisfies the uniform deterministic bound
    $
      0 <= e_i (t) <= ((1 - beta) / beta) D^2.
    $
  + Between node insertions with fixed $W$ and i.i.d. draws $bold(x)_t ~ P$, the expected error $mu_i (t) := bb(E)[e_i (t)]$ converges geometrically at rate $(1 - beta)^t$ to the exact local distortion
    $
      mu_i^* = ((1 - beta) / beta) Q_i (W),
      quad Q_i (W) := integral_(V_i (W)) norm(bold(x) - bold(w)_i)_2^2 dif P(bold(x)).
    $
]

#proof[
  First, observe by induction that every prototype stays inside $"conv"(S)$: both winner/neighbor adaptation $bold(w)' = (1 - epsilon) bold(w) + epsilon bold(x)_t$ ($epsilon in (0, 1)$) and midpoint insertion $bold(w)_r = (bold(w)_q + bold(w)_f) / 2$ form convex combinations of elements in $"conv"(S)$. Hence $norm(bold(x)_t - bold(w)_i)_2^2 <= D^2$ for all $i$ and $t$.

  Let $B := ((1 - beta) / beta) D^2$. At $t = 0$, $e_i (0) = 0 <= B$. Assume $e_i (t) <= B$ for all $i in V_t$. After error increment and optional node insertion, every node error before global decay is at most $e_i (t) + D^2 <= B + D^2$ (note that insertion scales $e_q$ and $e_f$ by $alpha < 1$ and sets $e_r = alpha e_q < B + D^2$). Applying global decay $(1 - beta)$ gives
  $
    e_i (t + 1) <= (1 - beta)(B + D^2)
      = (1 - beta) (D^2 / beta) = B,
  $
  establishing (1) by induction.

  For (2), at each non-insertion step, $e_i$ receives random increment $bold(1)_(bold(x)_t in V_i (W)) norm(bold(x)_t - bold(w)_i)_2^2$, whose expectation is precisely $Q_i (W)$, followed by multiplication by $(1 - beta)$. Taking expectations yields the affine contraction
  $
    mu_i (t + 1) = (1 - beta)(mu_i (t) + Q_i (W)).
  $
  Its unique fixed point satisfies $mu_i^* = (1 - beta)(mu_i^* + Q_i (W))$, i.e., $mu_i^* = ((1 - beta) / beta) Q_i (W)$, and subtracting yields $mu_i (t) - mu_i^* = (1 - beta)^t (mu_i (0) - mu_i^*)$.
]

Lemma II.C proves that GNG's error-selection heuristic $q = arg max_i e_i$ asymptotically tracks the Voronoi cell with *maximum local quantization distortion* $Q_i (W)$, driving greedy minimax error equalization across the manifold. Furthermore, setting $alpha = 1 / 2$ strictly conserves the error mass of the split cell ($e'_q + e'_r = 2 alpha e_q = e_q$), preventing runaway consecutive insertions at a single node.

#subsec("D", "Similarity Equivariance and Monotone Insertion Gain")

#proposition(title: "Similarity Equivariance and Strict Quantization Gain Under Midpoint Insertion", origin: "Original Proof")[
  Let $Q(W) := integral_(bb(R)^d) min_(i in V) norm(bold(x) - bold(w)_i)_2^2 dif P(bold(x))$ denote the global quantization distortion.
  + Under any similarity map $bold(x) |-> s Q bold(x) + bold(t)$ ($s > 0$, $Q in O(d)$), GNG prototypes transform as $bold(w)_i |-> s Q bold(w)_i + bold(t)$, node errors scale homogeneously as $e_i |-> s^2 e_i$, and the graph topology $G_t$ is strictly invariant.
  + Inserting vertex $r$ at $bold(w)_r = (bold(w)_q + bold(w)_f) / 2$ (with $bold(w)_q != bold(w)_f$) yields $Q(W union {bold(w)_r}) <= Q(W)$, with strict inequality whenever $bold(w)_r in.not W$ and $P(V_r (W union {bold(w)_r})) > 0$.
]

#proof[
  For (1), similarity maps scale all squared pairwise distances uniformly by $s^2$, preserving winner/runner-up indices $(s_1, s_2)$, scaling every error increment by $s^2$, and therefore preserving $arg max_i e_i$, $arg max_(j in N(q)) e_j$, and all edge ages $a_(i j)$.

  For (2), define $d_W (bold(x)) := min_(i in V) norm(bold(x) - bold(w)_i)_2^2$ and $d_(W') (bold(x)) := min(d_W (bold(x)), norm(bold(x) - bold(w)_r)_2^2)$ for $W' = W union {bold(w)_r}$. Clearly $d_(W') (bold(x)) <= d_W (bold(x))$ pointwise, so $Q(W') <= Q(W)$. On the interior of $V_r (W')$, $norm(bold(x) - bold(w)_r)_2^2 < d_W (bold(x))$ strictly whenever $bold(w)_r in.not W$. Since boundary hyperplanes carry zero measure relative to an interior of positive $P$-mass, integrating over $V_r (W')$ gives $Q(W) - Q(W') = integral_(V_r (W')) (d_W (bold(x)) - norm(bold(x) - bold(w)_r)_2^2) dif P(bold(x)) > 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

At step $t$ with current node count $M_t <= M_"max"$, finding $(s_1, s_2)$ scans $M_t$ prototypes in $cal(O)(M_t d)$ time without sorting. Updating the winner, its 1-hop neighbors $N_t (s_1)$, and incident edge ages costs only $cal(O)((1 + "deg"_(G_t) (s_1)) d)$, which is $cal(O)(d)$ on bounded-degree manifold graphs. Node insertion occurs once every $lambda$ steps in $cal(O)(M_t)$ time to scan $e_i$. Over $T$ total presentations, the worst-case runtime is $cal(O)(T M_"max" d)$, and when $T$ is comparable to $lambda M_"max"$, linear growth from $M_0 = 2$ halves the effective constant factor to approximately $T M_"max" d / 2$.

Space complexity is $cal(O)(M_"max" d + |E|)$. On manifolds of low intrinsic dimension $q <= 3$, the induced Delaunay subgraph is sparse ($|E| = cal(O)(M_"max")$), keeping adjacency-list memory negligible; in high-dimensional ambient spaces, however, Delaunay degree can approach $cal(O)(M_"max"^2)$. Unlike batch SOM or Neural Gas, GNG processes one point at a time in $cal(O)(1)$ working memory relative to stream length, making it naturally suited to embedded and online settings.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T M_"max" d)$ total; $cal(O)(M_t d)$ per step without sorting],
    [Space Complexity], [$cal(O)(M_"max" d + |E|)$; true $cal(O)(1)$ w.r.t. stream length $n$],
    [Cluster Geometry], [Adaptive induced Delaunay graph; arbitrary non-convex manifolds],
    [Hyperparameters], [$epsilon_b, epsilon_n$, $a_"max"$, $lambda$, $alpha, beta$, node cap $M_"max"$],
    [Invariance], [Similarity-equivariant; anisotropic-scale & stream-order sensitive],
    [Topology Learning], [Online CHL edges; variable intrinsic dimension per component],
    [Outlier Robustness], [Moderate-low (squared-error inflation can trigger outlier node spawns)],
  ),
  caption: [Operational and geometric profile of Growing Neural Gas (Alg. 177).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *3D Point-Cloud Surface Reconstruction and LiDAR Mapping:* Unorganized 3D laser scans sample complex non-convex surfaces with holes, bifurcations, and varying curvature. GNG's Competitive Hebbian edges recover the induced Delaunay mesh directly on the 2D manifold embedded in $bb(R)^3$ without assuming planar topology or fixed mesh resolution @fritzke1995gng.
+ *Non-Stationary Streaming Telemetry and Concept-Drift Tracking:* Industrial IoT and network-traffic streams evolve over time as new operating regimes appear and old ones vanish. Because GNG continuously inserts nodes in high-error regions while pruning stale edges in linear time (Theorem II.B) and deleting isolated vertices, its graph tracks drifting distributions without requiring full retraining @aggarwal2003clustream.
+ *Robotic Motion-Planning Roadmaps and Free-Space Topology:* Mobile robots exploring unknown environments sample collision-free configurations sequentially. GNG incrementally builds a sparse topological navigation graph whose edges remain confined to sampled free space ($V_(i j) inter S != emptyset$), directly providing connectivity roadmaps for path planning.
+ *Adaptive Color and Feature-Space Vector Quantization:* Multimodal image and video palettes often consist of sharply separated chromatic islands of unequal volume. By equalizing local distortion via Lemma II.C, GNG allocates prototypes proportionally to high-variance color manifolds without wasting dead nodes in empty chromatic space.
+ *Multi-Regime Biomedical and Sensor Manifold Discovery:* Physiological state spaces (e.g., sleep stages, gait dynamics) combine 1D transition paths with higher-dimensional attractor basins. GNG captures mixed-dimensional topology natively in a single graph, enabling downstream spectral or connected-component segmentation @vonluxburg2007tutorial.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GNG suffers from five practical failure modes: (i) *hyperparameter coupling*, where setting $a_"max"$ too small relative to $M_t$ prunes valid edges (Theorem II.B) while setting $a_"max"$ too large retains stale bridges across non-convex gaps; (ii) *outlier-driven node proliferation*, because squared-error accumulation $norm(bold(x)_t - bold(w)_(s_1))_2^2$ disproportionately inflates $e_i$ near extreme outliers; (iii) *curse of dimensionality* in high-$d$ spaces @beyer1999nearest, where concentration of distances makes second-nearest neighbors noisy and inflates graph degree; (iv) *lack of a direct cluster-count parameter $k$*, requiring post-hoc connected-component or graph-cut extraction; and (v) *premature stopping* if $T$ is small relative to $lambda M_"max"$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [3D Surfaces & Manifolds], [Excellent: recovers induced Delaunay mesh], [Tune $a_"max"$ to avoid bridges],
    [Drifting Data Streams], [Strong: online growth + stale edge pruning], [For density micro-clusters: DenStream @cao2006denstream],
    [Non-Convex Shapes], [Strong: disconnected CHL components], [Outlier bridges; compare HDBSCAN @campello2013hdbscan],
    [High-Dim. Sparse Data], [Moderate-weak: noisy $s_2$, dense graph], [Use Neural Gas (Alg. 179) or spectral cuts],
    [Exact Fixed-$k$ Partition], [Indirect: requires post-hoc graph merging], [Use $k$-means++ @arthur2007kmeanspp or Neural Gas],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for GNG.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Structural Unification of GNG with GCS, Neural Gas, and SOM", origin: "Original Derivation")[
  Growing Neural Gas occupies the exact structural intersection of three neural quantization families:
  + Relative to *Growing Cell Structures (GCS)* @fritzke1994gcs, GNG replaces fixed-$q$-dimensional simplicial complexes with unconstrained CHL edges, generalizing topological representation from uniform dimension $q$ to heterogeneous local intrinsic dimensions across components.
  + Relative to *Neural Gas (Alg. 179)* @martinetz1993neuralgas, GNG replaces $cal(O)(M log M)$ global distance-rank sorting and all-prototype updates with $cal(O)("deg"_G (s_1))$ 1-hop graph updates, trading NG's exact global energy potential for sparse dynamic growth and faster per-step adaptation.
  + Relative to *SOM (Alg. 176)* @kohonen1982som, GNG replaces static lattice distances $d_G (i, j)$ with data-induced Delaunay connectivity $E_P (W)$, eliminating topological folding and fixed-grid dead units.
]

#proof[
  In GCS @fritzke1994gcs, every inserted node must preserve a pure $q$-simplex manifold structure, constraining every vertex neighborhood to dimension $q$. Replacing simplex-preserving edge rewiring by Competitive Hebbian insertion ${s_1, s_2}$ removes the fixed-$q$ constraint while retaining error-driven midpoint insertion between $(q, f)$. Conversely, fixing $V_t equiv [M]$ and replacing 1-hop unweighted adaptation ($epsilon_b$ on $s_1$, $epsilon_n$ on $N_t (s_1)$) with continuous rank weights $h_lambda (k_i (bold(x)_t, W))$ recovers online Neural Gas @martinetz1993neuralgas; replacing $N_t (s_1)$ with a static grid kernel $h_sigma (s_1, j)$ recovers online SOM @kohonen1982som.
]

*ClusterOverview Benchmark Analysis:* Across the 4-dataset ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), GNG achieves Iris Sil = $0.441$, ARI = $0.712$; Customers Sil = $0.367$; Multi-Shape Sil = $0.609$ (ARI = N/A); and High-D Sil = $0.249$ (ARI = N/A), with a remarkable mean runtime of $0.304$ ms. GNG is by far the fastest of the four neural algorithms—over $8 times$ faster than SOM ($2.477$ ms), $5.6 times$ faster than ART ($1.700$ ms), and $4.4 times$ faster than Neural Gas ($1.327$ ms)—directly reflecting our Section III analysis: GNG avoids both sorting and full-codebook updates, touching only $(s_1, N(s_1))$ while starting from $M_0 = 2$ nodes. On Iris and Customers, adaptive error-driven node placement yields strong Silhouette ($0.441$ and $0.367$, the highest Customers score among Algs. 176--179) and high Iris ARI ($0.712$), substantially outperforming SOM ($0.290$ / $0.562$). On Multi-Shape ($0.609$) and High-D ($0.249$), GNG trails Neural Gas ($0.718$ and $0.303$) slightly, consistent with our theoretical analysis: in higher dimensions, noisy second-nearest neighbors degrade CHL edge precision, while finite-horizon insertion leaves slight residual quantization variance compared to annealed all-node Neural Gas updates.
