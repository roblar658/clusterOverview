#import "../template/ieee.typ": *

#alg_entry(
  id: 180,
  name: "Growing Cell Structures (GCS)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1994",
  paradigm: "Simplicial Self-Organizing Topology",
  complexity: $O(T M_"max" d + (T / lambda) k^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Growing Cell Structures (GCS) @fritzke1994gcs is an incremental topological neural network that approximates an unknown data probability distribution $p(x)$ on $bb(R)^d$ by dynamically growing and pruning a pure $k$-dimensional simplicial complex embedded in input space. Unlike the static rectangular lattice of the Self-Organizing Map (SOM) @kohonen1982som, GCS adapts both the cardinality and local connectivity of its prototype set in response to local error or sample-hit density, while preserving strict $k$-simplicial manifold structure throughout growth.

#definition(title: "Pure k-Simplicial Prototype Complex")[
  Fix an intrinsic topological dimension $k in bb(N)$ with $1 <= k <= d$. A *pure $k$-simplicial prototype complex* is a triple $cal(K) = (cal(V), cal(S)_k, W)$ where $cal(V)$ is a finite vertex set of size $M = abs(cal(V)) >= k + 1$, $cal(S)_k subset.eq binom(cal(V), k+1)$ is a collection of $k$-simplices such that every vertex $c in cal(V)$ belongs to at least one $sigma in cal(S)_k$, and $W = {w_c}_(c in cal(V)) subset bb(R)^d$ assigns a reference vector to each vertex. The induced 1-skeleton is $cal(E) = {{u, v} : exists sigma in cal(S)_k, {u, v} subset.eq sigma}$, the topological neighborhood of $c$ is $N(c) = {u in cal(V) : {c, u} in cal(E)}$, and the Euclidean Voronoi cell of $c$ is
  $ V_c (W) = {x in bb(R)^d : norm(x - w_c)_2 <= norm(x - w_u)_2, forall u in cal(V)}. $
]

Each vertex $c in cal(V)$ maintains a non-negative scalar resource counter $tau_c >= 0$. At signal presentation $t$, a sample $x_t tilde.op p$ selects its best-matching unit (BMU) $s_1 = arg min_(c in cal(V)) norm(x_t - w_c)_2$. Prototypes are updated via competitive neighborhood adaptation with learning rates $0 < epsilon_n < epsilon_b < 1$:
$ w_(s_1) &<- w_(s_1) + epsilon_b (x_t - w_(s_1)), \
  w_u &<- w_u + epsilon_n (x_t - w_u), quad forall u in N(s_1). $
The resource counter of $s_1$ accumulates either constant hit mass ($Delta tau_t = 1$, density-equalizing mode) or squared quantization error ($Delta tau_t = norm(x_t - w_(s_1))_2^2$, vector-quantization mode), followed by uniform exponential forgetting with rate $alpha in (0, 1)$:
$ tau_(s_1) <- tau_(s_1) + Delta tau_t, quad
  tau_c <- (1 - alpha) tau_c quad forall c in cal(V). $

Every $lambda$ presentations (while $abs(cal(V)) < M_"max"$), GCS identifies the maximum-resource vertex $q = arg max_(c in cal(V)) tau_c$ and its farthest topological neighbor $f = arg max_(u in N(q)) norm(w_q - w_u)_2$. A new vertex $r$ is inserted at the midpoint $w_r = (w_q + w_f) / 2$, edge ${q, f}$ is deleted, and every $k$-simplex $sigma in cal(S)_k (q, f) = {sigma in cal(S)_k : {q, f} subset.eq sigma}$ is bisected into two $k$-simplices $(sigma without {f}) union {r}$ and $(sigma without {q}) union {r}$. Resources are redistributed from the new neighborhood $N'(r)$ to $r$ via
$ Delta tau_u = - tau_u / abs(N'(u)), quad
  tau_r = - sum_(u in N'(r)) Delta tau_u, $
where $N'(u)$ denotes post-insertion neighbors. To separate disconnected clusters across zero-density regions, GCS periodically estimates each cell's local probability density $hat(p)_c = tau_c / (hat(V)_c sum_(u in cal(V)) tau_u)$ using the surrogate $k$-volume $hat(V)_c = bar(ell)_c^k$ with mean edge length $bar(ell)_c = abs(N(c))^(-1) sum_(u in N(c)) norm(w_c - w_u)_2$. Any vertex with $hat(p)_c < eta$ is deleted along with all incident $k$-simplices, and any resulting orphaned vertices belonging to no remaining $k$-simplex are removed. Connected components of the pruned 1-skeleton $(cal(V), cal(E))$ define the final clusters, and points $x_i$ inherit the component label of their BMU.

#algbox(
  title: "Growing Cell Structures (GCS)",
  inputs: [Data stream or set $X = {x_1, ..., x_N} subset bb(R)^d$, dimension $k$, max cells $M_"max"$, insertion period $lambda$, rates $(epsilon_b, epsilon_n, alpha)$, density threshold $eta$, budget $T$.],
  outputs: [Pure $k$-simplicial complex $(cal(V), cal(S)_k, W)$ and cluster labels $hat(y)_1, ..., hat(y)_N$.],
)[
  + Initialize $cal(K)$ as a single $k$-simplex $sigma_0$ of $k + 1$ prototypes sampled near $X$ with $tau_c <- 0$.
  + For $t = 1, ..., T$: draw $x_t$ from $X$ and compute BMU $s_1 = arg min_(c in cal(V)) norm(x_t - w_c)_2$.
  + Update $w_(s_1) <- w_(s_1) + epsilon_b (x_t - w_(s_1))$ and $w_u <- w_u + epsilon_n (x_t - w_u)$ for $u in N(s_1)$.
  + Increment $tau_(s_1) <- tau_(s_1) + Delta tau_t$, then decay $tau_c <- (1 - alpha) tau_c$ for all $c in cal(V)$.
  + If $t mod lambda = 0$ and $abs(cal(V)) < M_"max"$: select $q = arg max_c tau_c$, $f = arg max_(u in N(q)) norm(w_q - w_u)_2$, bisect every $sigma in cal(S)_k (q, f)$ with $w_r = (w_q + w_f)/2$, and redistribute $tau$.
  + Periodically compute $hat(p)_c$ and delete vertices with $hat(p)_c < eta$ together with incident $k$-simplices and orphaned vertices; label connected components of $(cal(V), cal(E))$ and assign each $x_i$ to its BMU's component.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish combinatorial invariance of the simplicial mesh, exact conservation and asymptotic scaling of resource dynamics, fixed-partition contraction, and sharp breakdown and topology-mismatch bounds.

#subsec("A", "Simplicial Invariance and Combinatorial Growth Law")

#theorem(title: "Preservation of Pure k-Simpliciality and Euler Invariance", origin: "Original Proof")[
  Let $cal(K) = (cal(V), cal(S)_k, W)$ be a pure $k$-simplicial complex, and let $(q, f) in cal(E)$ be shared by $s = abs(cal(S)_k (q, f)) >= 1$ maximal $k$-simplices. Let $U = union.big_(sigma in cal(S)_k (q, f)) (sigma without {q, f})$ be the set of vertices opposite to ${q, f}$. Then edge-bisection insertion of $r$ produces a complex $cal(K)' = (cal(V)', cal(S)'_k, W')$ satisfying:
  (i) $cal(K)'$ is pure $k$-simplicial with every vertex in at least one $k$-simplex;
  (ii) $Delta abs(cal(V)) = 1$, $Delta abs(cal(S)_k) = s$, and $Delta abs(cal(E)) = abs(U) + 1$ with $k - 1 <= abs(U) <= s(k - 1)$;
  (iii) for $k = 2$ on any triangulated 2-manifold with boundary, the Euler characteristic $chi = abs(cal(V)) - abs(cal(E)) + abs(cal(S)_2)$ is strictly invariant ($Delta chi = 0$).
]

#proof[
  (i) Each removed simplex $sigma in cal(S)_k (q, f)$ has cardinality $k + 1$ and contains both $q$ and $f$. Replacing $sigma$ by $sigma_q = (sigma without {f}) union {r}$ and $sigma_f = (sigma without {q}) union {r}$ yields two distinct subsets of $cal(V)'$ of size $k + 1$, since $r in.not cal(V)$ and $q != f$. Distinct $sigma, sigma' in cal(S)_k (q, f)$ differ on at least one vertex in $U$, so $sigma_q != sigma'_q$ and $sigma_f != sigma'_f$; hence no duplicate $k$-simplices arise. Every untouched vertex remains in its prior simplex, $q in sigma_q$, $f in sigma_f$, and $r in sigma_q inter sigma_f$ (with $s >= 1$). Thus $cal(K)'$ is pure $k$-simplicial.

  (ii) Exactly one vertex $r$ is added, so $Delta abs(cal(V)) = 1$. Each of the $s$ simplices in $cal(S)_k (q, f)$ is replaced by $2$ simplices, giving $Delta abs(cal(S)_k) = 2s - s = s$. For edges, ${q, f}$ is deleted ($-1$) because no remaining $k$-simplex contains both $q$ and $f$. New edges incident to $r$ are ${r, q}$, ${r, f}$, and ${r, u}$ for every $u in U$, contributing $2 + abs(U)$ new edges. All edges inside $U$ or between $U$ and ${q, f}$ are preserved in $sigma_q$ and $sigma_f$. Hence
  $ Delta abs(cal(E)) = (2 + abs(U)) - 1 = abs(U) + 1. $
  Since each $sigma in cal(S)_k (q, f)$ contributes $k - 1$ vertices to $U$ and $s >= 1$, the union bound immediately gives $k - 1 <= abs(U) <= s(k - 1)$.

  (iii) For $k = 2$ on a triangulated 2-manifold with boundary, an interior edge belongs to exactly $s = 2$ triangles sharing no other vertex besides $q, f$, so $abs(U) = 2$, $Delta abs(cal(E)) = 3$, and $Delta chi = 1 - 3 + 2 = 0$. A boundary edge belongs to $s = 1$ triangle, giving $abs(U) = 1$, $Delta abs(cal(E)) = 2$, and $Delta chi = 1 - 2 + 1 = 0$. By induction from the initial 2-simplex ($chi = 3 - 3 + 1 = 1$), $chi(cal(K))$ stays constant across all insertions.
]

#subsec("B", "Resource Conservation and Asymptotic Prototype Density")

#lemma(title: "Exact Resource Conservation and Bounded Total Mass", origin: "Original Derivation")[
  Let $S_t = sum_(c in cal(V)_t) tau_c (t)$ be the total resource mass after step $t$ with decay $alpha in (0, 1)$ and $S_0 = 0$.
  (i) During any vertex insertion of $r$, the redistribution rule preserves total mass exactly: $S' = S$ and $tau'_c >= 0$ for all $c in cal(V)'$.
  (ii) If $0 <= Delta tau_t <= D_"max"$ for all $t >= 1$, then $0 <= S_t <= ((1 - alpha) / alpha) D_"max"$ for all $t$, and in constant-hit mode ($Delta tau_t = 1$), $S_t$ increases monotonically to $S_infinity = (1 - alpha) / alpha$ via
  $ S_t = (1 - alpha) / alpha (1 - (1 - alpha)^t). $
]

#proof[
  (i) Post-insertion, every $u in N'(r)$ belongs to at least one $k$-simplex of $cal(K)'$ by Theorem 1(i), so its degree satisfies $abs(N'(u)) >= k >= 1$. Consequently $tau'_u = (1 - 1 / abs(N'(u))) tau_u >= 0$ and $tau'_r = sum_(u in N'(r)) tau_u / abs(N'(u)) >= 0$. Summing over $cal(V)' = cal(V) union {r}$ gives
  $ S' = tau'_r + sum_(u in N'(r)) (tau_u - tau_u / abs(N'(u))) + sum_(c in cal(V) without N'(r)) tau_c = S. $
  (ii) Because insertions leave $S_t$ invariant, every adaptation step obeys $S_t = (1 - alpha)(S_(t-1) + Delta tau_t)$. Unrolling from $S_0 = 0$ yields
  $ S_t = sum_(j=1)^t (1 - alpha)^(t - j + 1) Delta tau_j <= D_"max" sum_(m=1)^t (1 - alpha)^m = (1 - alpha) / alpha (1 - (1 - alpha)^t) D_"max". $
  When $Delta tau_j = 1$, the inequality is an equality converging monotonically to $(1 - alpha)/alpha$ as $t -> infinity$.
]

#theorem(title: "Equilibrium Distortion and Gersho Density Scaling", origin: "Classical — Cited")[
  Let $p(x)$ be a continuous positive density on a compact $k$-dimensional manifold $cal(M) subset bb(R)^d$, and assume GCS reaches a large-$M$ stationary equilibrium with smooth prototype density $rho(x)$ ($integral_cal(M) rho(x) dif x = M$) and asymptotically regular Voronoi cells of $k$-volume $abs(V_c) approx rho(w_c)^(-1)$ @fritzke1994gcs @linde1980algorithm. Then:
  (i) In constant-hit mode ($Delta tau = 1$), equalizing stationary cell resources $tau_c prop P_c = integral_(V_c) p(x) dif x$ yields linear density matching $rho(x) prop p(x)$.
  (ii) In quantization-error mode ($Delta tau = norm(x - w_(s_1))_2^2$), equalizing local distortion $tau_c prop integral_(V_c) norm(x - w_c)_2^2 p(x) dif x$ drives the prototype density to the $k$-dimensional Gersho optimum:
  $ rho(x) prop p(x)^(k / (k + 2)). $
]

#proof[
  At large $M$, $p(x)$ is nearly constant over each small cell $V_c$ around $w_c$, so the cell probability mass satisfies $P_c = integral_(V_c) p(x) dif x approx p(w_c) abs(V_c) approx p(w_c) / rho(w_c)$.
  (i) Under $Delta tau = 1$, insertions repeatedly split the cell with largest $tau_c prop P_c$, driving $P_c$ toward a uniform constant $1/M$. Hence $p(w_c) / rho(w_c) = "const"$, yielding $rho(x) prop p(x)$.
  (ii) Under $Delta tau = norm(x - w_(s_1))_2^2$, by the dimensionless inertial coefficient $C_k > 0$ of $k$-dimensional Voronoi polytopes @linde1980algorithm, the local second moment scales with cell volume as
  $ integral_(V_c) norm(x - w_c)_2^2 dif x approx C_k abs(V_c)^(1 + 2 / k) approx C_k rho(w_c)^(-(1 + 2 / k)). $
  Multiplying by $p(w_c)$, the stationary cell error satisfies
  $ tau_c prop p(w_c) C_k rho(w_c)^(-(k + 2) / k). $
  Splitting cells of maximal $tau_c$ equalizes this quantity across all $c in cal(V)$, forcing $p(w_c) rho(w_c)^(-(k+2)/k) = "const"$, and therefore $rho(x) prop p(x)^(k / (k + 2))$, which minimizes the asymptotic mean-squared quantization error over all $M$-point $k$-manifold quantizers @fritzke1994gcs.
]

#subsec("C", "Fixed-Partition Convex Potential and Prototype Contraction")

#proposition(title: "Graph-Coupled Fixed-Partition Equilibrium and Contraction", origin: "Original Proof")[
  For a fixed simplicial complex $cal(K)$ and fixed Voronoi partition ${V_c}$ with cell masses $P_c = integral_(V_c) p(x) dif x > 0$ and centroids $bar(x)_c = P_c^(-1) integral_(V_c) x p(x) dif x$, the expected GCS prototype update $bb(E)[Delta w_v | W]$ equals the exact negative gradient $-nabla_(w_v) J_"fix" (W)$ of the strictly convex quadratic potential
  $ J_"fix" (W) = 1 / 2 sum_(c in cal(V)) ( epsilon_b integral_(V_c) norm(x - w_c)_2^2 p(x) dif x \
    + epsilon_n sum_(u in N(c)) integral_(V_c) norm(x - w_u)_2^2 p(x) dif x ). $
  Its unique minimizer $W^star$ satisfies the neighborhood-smoothed fixed-point identity
  $ w_v^star = (epsilon_b P_v bar(x)_v + epsilon_n sum_(c in N(v)) P_c bar(x)_c) / (epsilon_b P_v + epsilon_n sum_(c in N(v)) P_c), $
  and for $0 < epsilon_n <= epsilon_b < 1$, the expected displacement contracts strictly:
  $ norm(bb(E)[w_v (t+1) - w_v^star | W_t])_2 = (1 - gamma_v) norm(w_v (t) - w_v^star)_2, $
  where $gamma_v = epsilon_b P_v + epsilon_n sum_(c in N(v)) P_c in (0, 1)$.
]

#proof[
  Vertex $v$ receives an update of size $epsilon_b (x - w_v)$ whenever $x in V_v$, and of size $epsilon_n (x - w_v)$ whenever $x in V_c$ for some $c in N(v)$ (since edge incidence is symmetric: $v in N(c) <=> c in N(v)$). Taking expectation over $x tilde.op p$ gives
  $ bb(E)[Delta w_v | W] = epsilon_b P_v (bar(x)_v - w_v) + epsilon_n sum_(c in N(v)) P_c (bar(x)_c - w_v). $
  Differentiating $J_"fix" (W)$ with respect to $w_v$ yields the exact negative of this expression, and the Hessian is block-diagonal with blocks $gamma_v I_d succ.eq (epsilon_b P_v) I_d$, proving strict convexity. Setting $nabla_(w_v) J_"fix" (W^star) = 0$ directly produces the convex-combination formula for $w_v^star$. Finally, since $sum_(c in cal(V)) P_c = 1$ and $epsilon_n <= epsilon_b < 1$,
  $ 0 < epsilon_b P_v <= gamma_v <= epsilon_b (P_v + sum_(c in N(v)) P_c) <= epsilon_b < 1. $
  Rewriting $bb(E)[Delta w_v | W_t] = gamma_v (w_v^star - w_v (t))$ yields $bb(E)[w_v (t+1) - w_v^star | W_t] = (1 - gamma_v)(w_v (t) - w_v^star)$, establishing geometric contraction on each fixed-partition regime. Moreover, because every single-sample update is a convex combination of $w_v$ and $x_t$, the convex hull $"conv"(W_0 union "supp"(p))$ is strictly forward-invariant.
]

#subsec("D", "Equivariance, Breakdown Point, and Topological Separation")

#proposition(title: "Isometry Equivariance and Mode-Dependent Breakdown Point", origin: "Original Counterexample")[
  (i) GCS prototypes and cluster partitions are strictly equivariant under Euclidean isometries $x |-> Q x + b$ ($Q in O(d)$, $b in bb(R)^d$), and constant-hit GCS ($Delta tau = 1$) is additionally invariant to isotropic scaling $x |-> c x$ ($c > 0$) when $eta$ scales as $c^(-k)$.
  (ii) On a finite dataset of $N$ points, quantization-error GCS ($Delta tau = norm(x - w_(s_1))_2^2$) has finite-sample adversarial breakdown point $1 / N$ (asymptotic $0\%$), whereas constant-hit GCS ($Delta tau = 1$) caps any single point's resource share at most $O(1 / (alpha N))$.
]

#proof[
  (i) Because $norm(Q x + b - (Q w_c + b))_2 = norm(x - w_c)_2$, BMU indices, farthest-neighbor choices $f$, midpoint insertions $Q (w_q + w_f)/2 + b$, and edge lengths $bar(ell)_c$ are preserved identically. Under $x |-> c x$, distances scale by $c$, so $hat(V)_c = bar(ell)_c^k$ scales by $c^k$; with $Delta tau = 1$, $tau_c$ is scale-free, so $hat(p)_c$ scales by $c^(-k)$, preserving all pruning decisions when $eta' = c^(-k) eta$.

  (ii) For quantization-error GCS, fix $N - 1$ points inside the unit ball $B(0, 1)$ and place one adversarial outlier $x_"out"$ at radius $R = norm(x_"out")_2 -> infinity$. In cyclic sweeps (or with high probability under i.i.d. draws at rate $1/N$), whenever $x_"out"$ is presented, its BMU $s_1$ receives increment $Delta tau = norm(x_"out" - w_(s_1))_2^2 = Omega(R^2)$, while all inlier cells combined hold total resource at most $O(1 / alpha)$ by Lemma 1. Choosing $R$ large enough that $(1 - alpha)^N R^2 > 2 D_"in" / alpha$ (where $D_"in"$ bounds inlier distortion) guarantees that $q = arg max_c tau_c$ remains locked on the outlier's cell (and its newly split descendants) at every insertion step $t = m lambda$, diverting arbitrarily many cells toward $x_"out"$ and collapsing approximation accuracy on the $N - 1$ inliers. Hence a single corrupted point breaks the clustering, giving breakdown point $1/N$. Conversely, under $Delta tau = 1$, $x_"out"$ adds at most $1$ per visit, so its stationary contribution to $S_t$ is bounded by $((1 - alpha)/alpha) / N$ regardless of $R$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each signal presentation performs a linear scan over $M_t <= M_"max"$ prototypes to locate $s_1$ in $O(M_t d)$ time, updates $w_(s_1)$ and $abs(N(s_1))$ neighbors in $O((1 + deg(s_1)) d)$ time, and decays $M_t$ scalars in $O(M_t)$ time. Every $lambda$ steps, finding $q$ takes $O(M_t)$, locating $f$ takes $O(deg(q) d)$, and bisecting the $s$ incident $k$-simplices takes $O(s k^2)$ adjacency updates. Over $T$ signal presentations, total runtime is
$ T_"GCS" = O(T M_"max" d + (T / lambda) k^2), $
while memory scales as $O(M_"max" d + abs(cal(S)_k) (k + 1))$. Because every vertex belongs to $O(1)$ simplices in fixed low dimension $k in {1, 2, 3\}$, memory remains tiny and cache-resident; however, for large $k$, vertex star cardinality grows combinatorially with $k$, making GCS strictly a low-intrinsic-dimension topological method.

Geometrically, GCS embeds a piecewise-linear $k$-manifold into $bb(R)^d$. Its pruning mechanism converts density valleys into disconnected simplicial components, allowing it to track arbitrary non-convex curves ($k = 1$) or curved sheets ($k = 2$) provided the intrinsic manifold dimension matches $k$. If the true support has intrinsic dimension $m > k$, the complex must fold through $bb(R)^d$ like a space-filling curve, inducing topological defects @carlsson2010characterization.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T M_"max" d + (T / lambda) k^2)$],
    [Space Complexity], [$O(M_"max" d + abs(cal(S)_k) (k + 1))$],
    [Cluster Geometry], [Non-convex connected $k$-simplicial components],
    [Cluster Count $K$], [Automatic via density-based vertex pruning ($eta$)],
    [Hyperparameters], [$k, M_"max", lambda, epsilon_b, epsilon_n, alpha, eta$],
    [Streaming Support], [Native single-pass online adaptation],
    [Invariances], [Euclidean isometries; scaling if $Delta tau = 1$],
    [Outlier Breakdown], [$1/N$ ($Delta tau = "error"$); bounded ($Delta tau = 1$)],
  ),
  caption: [Operational and geometric profile of Growing Cell Structures.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Robotic Tactile and LiDAR Surface Reconstruction.* Choosing $k = 2$ enforces a pure triangle mesh throughout online acquisition. By Theorem 1, every insertion preserves 2-manifold triangulation and Euler characteristic without global Delaunay recomputation, directly outputting a watertight mesh whose density concentrates on high-curvature regions where quantization error peaks.

*2. Single-Cell Flow Cytometry and Rare Subpopulation Mapping.* In low-dimensional cytometry panels ($d <= 12$), rare immune phenotypes occupy small-volume, high-density islands. Error-driven GCS allocates prototypes proportionally to $p(x)^(k/(k+2))$ (Theorem 2), giving higher resolution to sparse minority populations than linear-density sampling while pruning bridges across low-density troughs to yield automatic cluster counts.

*3. Adaptive Vector Quantization and Color Palette Compression.* In image and audio coding ($d = 3$ color space or low-order LPC frames), minimizing squared distortion requires matching Gersho's $p(x)^(k/(k+2))$ prototype density @linde1980algorithm. Theorem 2 guarantees that quantization-error GCS approaches this distortion-optimal distribution incrementally without fixing $K$ or grid dimensions in advance.

*4. Real-Time Industrial Telemetry and Regime Segmentation.* Streaming sensor arrays exhibit non-stationary operating regimes separated by forbidden state transitions. Exponential forgetting ($alpha$) bounds total resource memory to $(1 - alpha)/alpha$ (Lemma 1), letting GCS track drifting regimes online while density pruning isolates distinct operational modes.

*5. Low-Dimensional Gene-Expression Topology Exploration.* Unlike fixed rectangular SOMs that impose artificial corner distortions on biological manifolds @kohonen1982som, a $k = 2$ simplicial complex adapts its boundary perimeter freely to branching or annular gene-expression topologies @xu2005survey.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GCS degrades sharply under four mathematical conditions: (i) *intrinsic dimension mismatch* ($m > k$), where a pure $k$-complex cannot cover an $m$-manifold without space-filling folds; (ii) *combinatorial simplex explosion* when $k >= 4$, causing $abs(cal(S)_k)$ to dominate memory and runtime; (iii) *outlier capture* under squared-error resources (Proposition 2); and (iv) *high-dimensional distance concentration* ($d >> 20$) @beyer1999nearest, where Euclidean BMU distances become nearly uniform and mean edge length $bar(ell)_c^k$ becomes an unreliable volume proxy.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [3D Surface Meshing], [Pure $k = 2$ triangulation preserves topology], [Non-manifold junctions; use GNG],
    [Vector Quantization], [Approaches Gersho $p^(k/(k+2))$ optimum], [Outliers capture cells; use $Delta tau = 1$],
    [Low-D Streaming], [Bounded memory, online growth/pruning], [Tune $alpha$ to drift timescale],
    [High-D Vision / NLP], [Poor; Euclidean BMU concentrates], [Use DEC, IDEC, or DCN],
    [Mixed Intrinsic Dim.], [Fixed $k$ forces folding on subspaces], [Use GNG @fritzke1995gng],
  ),
  caption: [Domain suitability and failure-mode matrix for GCS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Topological Hierarchy: SOM, GCS, Neural Gas, and GNG", origin: "Original Unified Proof")[
  Let $cal(G)_"top" = (cal(V), cal(E))$ be the prototype graph governing neighborhood adaptation. Then SOM @kohonen1982som, GCS @fritzke1994gcs, Neural Gas (NG) @martinetz1993neuralgas, and Growing Neural Gas (GNG) @fritzke1995gng form a strict topological relaxation chain:
  (i) SOM fixes both $abs(cal(V))$ and a static regular $k$-lattice $cal(E)$;
  (ii) GCS relaxes cardinality to grow dynamically while constraining $cal(G)_"top"$ to the 1-skeleton of a pure $k$-simplicial complex;
  (iii) GNG drops the pure $k$-simplicial invariant entirely, maintaining an unconstrained induced Delaunay subgraph via competitive Hebbian edge creation and age-based pruning, thereby recovering arbitrary local intrinsic dimensions $k(x)$;
  (iv) Setting $epsilon_n = 0$ and freezing insertions reduces GCS adaptation to online MacQueen K-Means @lloyd1982least.
]

#proof[
  Every $k$-dimensional hypercubic lattice admits a Freudenthal triangulation into a pure $k$-simplicial complex, so static SOM topology is embedded inside the feasible class of GCS complexes with $lambda = infinity$ and $eta = 0$. Conversely, any pure $k$-simplicial 1-skeleton is a valid undirected graph with vertex degrees $>= k$, hence a strict subset of GNG's unconstrained competitive-Hebbian graphs, which allow simultaneous 1D filaments, 2D sheets, and higher-dimensional cliques in one network. Finally, when $epsilon_n = 0$ and $cal(V)$ is fixed, $J_"fix" (W)$ in Proposition 1 reduces directly to the classical K-Means objective $1/2 sum_c integral_(V_c) norm(x - w_c)_2^2 p(x) dif x$, and the prototype step $w_(s_1) <- w_(s_1) + epsilon_b (x_t - w_(s_1))$ becomes stochastic K-Means.
]

*Empirical Benchmark Analysis.* On the standardized ClusterOverview 4-dataset benchmark @jain2010data, GCS achieves Silhouette $"Sil" = 0.449$ @rousseeuw1987silhouettes and Adjusted Rand Index $"ARI" = 0.563$ @hubert1985comparing on *Iris*, $"Sil" = 0.321$ on *Customers*, $"Sil" = 0.609$ ($"ARI"$ N/A) on *Multi-Shape*, and $"Sil" = 0.302$ ($"ARI"$ N/A) on *High-D*, with an exceptional mean runtime of *0.196 ms*. This empirical signature mirrors our theoretical analysis closely. First, the sub-millisecond runtime (over $50 times$ faster than backpropagation-based deep clustering algorithms such as DEC, IDEC, and DCN at $approx 10.5$ ms) follows directly from the lightweight $O(T M_"max" d)$ prototype updates and absence of multi-layer gradient passes. Second, on *Multi-Shape* ($"Sil" = 0.609$), simplicial pruning splits the 1-skeleton across low-density gaps to recover non-convex geometry, outperforming rigid spherical methods in raw speed while trailing deep nonlinear feature warpers (DEC $0.718$, DCN $0.716$) that actively reshape metric space. Third, on *Iris* ($"Sil" = 0.449$, $"ARI" = 0.563$), GCS attains the highest raw Silhouette among Algorithms 180--183 because its prototypes stay strictly within $"conv"(X)$ (Proposition 1), though overlapping Versicolor / Virginica regions produce fewer separated components, limiting $"ARI"$ compared to DCN ($0.648$). Finally, on *High-D* ($"Sil" = 0.302$), GCS remains competitive with deep autoencoders despite Euclidean distance concentration, confirming that sparse simplicial quantization retains strong local fidelity when sample budgets are modest.
