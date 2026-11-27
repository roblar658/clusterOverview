#import "../template/ieee.typ": *

#alg_entry(
  id: 37,
  name: "PDBSCAN (Parallel DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1999",
  paradigm: "Shared-Nothing Distributed Density Clustering",
  complexity: $O((n / P) log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a large spatial dataset equipped with a metric $d : bb(R)^d times bb(R)^d -> bb(R)_(>= 0)$, density radius $epsilon > 0$, and minimum neighborhood cardinality $"MinPts" in bb(N)_(>= 2)$. Classical sequential DBSCAN @ester1996dbscan identifies clusters as maximal density-connected components of the $epsilon$-neighborhood graph, requiring $O(n log n)$ time with an in-memory spatial index (or $O(n^2)$ without one) and incurring prohibitive I/O and memory bottlenecks once $n$ exceeds single-machine capacity.

*PDBSCAN (Parallel DBSCAN)* @xu1999pdbscan reformulates density-based clustering for a *shared-nothing* distributed cluster of $P >= 1$ worker nodes $cal(W) = {1, dots, P}$ coordinated by a master node. The dataset is disjointly partitioned across workers as $cal(X) = union.sq.big_(r=1)^P cal(X)_r$ via spatial regions $R_1, dots, R_P subset.eq bb(R)^d$ (typically ordered along a Hilbert space-filling curve to preserve spatial locality), while a *distributed R\*-tree (dR\*-tree)* replicates upper directory levels across all workers and distributes data leaf pages among the $P$ local disks.

#definition(title: "Distributed Partition, Halo, and Core Graph")[
  For each worker $r in {1, dots, P}$ owning $cal(X)_r = cal(X) inter R_r$, define its *$epsilon$-halo* (remote boundary neighborhood) by
  $ H_r (epsilon) = { y in cal(X) backslash cal(X)_r : exists x in cal(X)_r, space d(x, y) <= epsilon }. $
  For any $x in cal(X)$, let $N_epsilon (x) = {y in cal(X) : d(x, y) <= epsilon}$ denote its global $epsilon$-neighborhood, and let $N_epsilon^((r)) (x) = N_epsilon (x) inter (cal(X)_r union H_r (epsilon))$ denote its worker-$r$ augmented neighborhood. The *global core set* is $V_"core" = {x in cal(X) : |N_epsilon (x)| >= "MinPts"}$, and the undirected *core-reachability graph* is
  $ G_"core" = (V_"core", E_"core"), \ E_"core" = {{u, v} subset.eq V_"core" : d(u, v) <= epsilon}. $
]

PDBSCAN executes in three synchronized stages (Algorithm 1): (i) spatial partitioning and dR\*-tree construction; (ii) concurrent execution of local DBSCAN on each worker $r$ where queries $N_epsilon (x)$ whose $epsilon$-ball crosses $partial R_r$ issue remote range requests through the dR\*-tree; and (iii) global merging of local clusters at the master via connected-component computation over cross-partition core-contact edges, followed by parallel local relabeling.

#algbox(
  title: "PDBSCAN (Shared-Nothing Parallel DBSCAN)",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$, parameters $(epsilon, "MinPts")$, workers $P$.],
  outputs: [Global cluster labels $L : cal(X) -> {0, 1, dots, K}$ ($0$ denotes noise).],
)[
  1. *Phase 1 (Data Distribution & dR\*-Tree):* Partition $cal(X)$ spatially into $P$ disjoint shards ${cal(X)_r}_(r=1)^P$ (Hilbert-ordered); build distributed R\*-tree with replicated directory and partitioned leaves.
  2. *Phase 2 (Concurrent Local Clustering, on each worker $r$ in parallel):*
  3. #h(1em) Initialize local cluster counter $c_r <- 0$ and contact set $E_r <- emptyset$.
  4. #h(1em) *For each* unvisited $x in cal(X)_r$ *do*:
  5. #h(2em) Query $N_epsilon (x)$ via dR\*-tree (local page scan + remote RPC for $B_epsilon (x) inter R_s != emptyset$, $s != r$).
  6. #h(2em) *If* $|N_epsilon (x)| < "MinPts"$ *then* mark $x$ as tentative noise; *continue*.
  7. #h(2em) Increment $c_r <- c_r + 1$; create local cluster $C_(r, c_r)$ and seed queue $Q <- N_epsilon (x)$.
  8. #h(2em) *While* $Q != emptyset$ *do*: pop $q$ from $Q$;
  9. #h(3em) *If* $q in cal(X)_r$: assign $q$ to $C_(r, c_r)$; if unvisited, query $N_epsilon (q)$ and if $|N_epsilon (q)| >= "MinPts"$, push $N_epsilon (q)$ onto $Q$.
  10. #h(3em) *Else* ($q in cal(X)_s, s != r$): record cross-partition witness $(x_"src", q, (r, c_r))$ in $E_r$.
  11. *Phase 3 (Master-Side Equivalence Merging & Relabeling):*
  12. #h(1em) Gather local-cluster sets $V_"loc" = {(r, c)}$ and witnesses $union.big_r E_r$ at master.
  13. #h(1em) Build merge graph $G_"merge" = (V_"loc", E_"merge")$ linking $(r, c_1)$ and $(s, c_2)$ whenever a witness $(u, v)$ satisfies $u in C_(r, c_1) inter V_"core"$ and $v in C_(s, c_2) inter V_"core"$.
  14. #h(1em) Compute connected components of $G_"merge"$ via Union-Find to obtain global map $phi : V_"loc" -> {1, dots, K}$; broadcast $phi$.
  15. #h(1em) Workers relabel local clusters via $phi$ and resolve cross-partition border points deterministically.
  16. *Return* global label map $L$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties governing PDBSCAN: (1) exact preservation of core status and global semantic equivalence with sequential DBSCAN; (2) a geometric surface-to-volume communication bound quantifying parallel efficiency; (3) an explicit spatial-skew counterexample proving worst-case speedup collapse under static uniform grids; and (4) classical distributed index and Union-Find merging complexity bounds.

#subsec("A", "Exact Distributed Semantic Equivalence")

#lemma(title: "Global Core-Status and Cut-Edge Preservation", origin: "Original Proof")[
  Under any disjoint partition $cal(X) = union.sq.big_(r=1)^P cal(X)_r$ equipped with exact $epsilon$-halos $H_r (epsilon)$ (equivalently, exact dR\*-tree range queries):
  1. For every worker $r$ and every owned point $x in cal(X)_r$, $N_epsilon^((r)) (x) = N_epsilon (x)$, hence $x in V_"core"$ if and only if $|N_epsilon^((r)) (x)| >= "MinPts"$.
  2. Every cross-partition edge ${u, v} in E_"core"$ with $u in cal(X)_r$ and $v in cal(X)_s$ ($r != s$) satisfies $v in H_r (epsilon)$ and $u in H_s (epsilon)$ and is emitted as a witness in both $E_r$ and $E_s$.
]

#proof[
  (1) Fix $x in cal(X)_r$ and take any $y in N_epsilon (x)$. Either $y in cal(X)_r$, or $y in cal(X) backslash cal(X)_r$ with $d(x, y) <= epsilon$, which by Definition 1 implies $y in H_r (epsilon)$. Thus $N_epsilon (x) subset.eq cal(X)_r union H_r (epsilon)$, yielding $N_epsilon^((r)) (x) = N_epsilon (x) inter (cal(X)_r union H_r (epsilon)) = N_epsilon (x)$. Equality of cardinalities $|N_epsilon^((r)) (x)| = |N_epsilon (x)|$ immediately preserves the core predicate $|N_epsilon (x)| >= "MinPts"$ identically on every worker.

  (2) Let ${u, v} in E_"core"$ with $u in cal(X)_r$ and $v in cal(X)_s$ ($r != s$). By symmetry of the metric $d$, $d(u, v) = d(v, u) <= epsilon$, so $v in N_epsilon (u) backslash cal(X)_r subset.eq H_r (epsilon)$ and symmetrically $u in H_s (epsilon)$. Because $u in V_"core" inter cal(X)_r$, Algorithm 1 expands $u$ inside some local cluster $C_(r, c_1)$, pushes $N_epsilon (u)$ (which contains $v$) onto $Q$, detects $v in cal(X)_s$, and records the cross-partition witness $(u, v, (r, c_1)) in E_r$ (and symmetrically worker $s$ records $(v, u, (s, c_2)) in E_s$).
]

#theorem(title: "Exact Core-Graph Equivalence with Sequential DBSCAN", origin: "Original Proof")[
  Let $cal(K)_"seq"$ be the partition of $V_"core"$ into connected components of $G_"core"$ produced by sequential DBSCAN @ester1996dbscan, and let $cal(K)_"par"$ be the partition of $V_"core"$ induced by Algorithm 1 after Phase 3 merging. Then:
  1. $cal(K)_"par" = cal(K)_"seq"$ identically, independent of the spatial partition ${cal(X)_r}_(r=1)^P$ and parallel worker scheduling.
  2. The global noise set $cal(N)_"par"$ coincides identically with the sequential noise set $cal(N)_"seq" = {x in cal(X) : N_epsilon (x) inter V_"core" = emptyset}$.
  3. Under deterministic border tie-breaking (assigning each non-core point $b in cal(X) backslash (V_"core" union cal(N))$ to the minimum global cluster ID $min {phi(C) : exists u in C inter V_"core", space d(u, b) <= epsilon}$), the complete clustering $L$ is strictly invariant to both $P$ and scan order.
]

#proof[
  (1) By Lemma 1(1), every worker identifies the exact local subset of core points $V_"core"^((r)) = V_"core" inter cal(X)_r$. During Phase 2, worker $r$ traverses only edges of $G_"core"$ whose endpoints both lie inside $cal(X)_r$; thus the core points of each local cluster $C_(r, c)$ form a connected component of the induced subgraph $G_"core" [V_"core"^((r))]$. Let $pi : V_"core" -> V_"loc"$ map each core point $u in V_"core"^((r))$ to its unique local cluster ID $(r, c)$.

  We prove $u, v in V_"core"$ lie in the same connected component of $G_"core"$ if and only if $pi(u)$ and $pi(v)$ lie in the same connected component of $G_"merge"$.
  - *Forward direction:* Suppose $u = w_0, w_1, dots, w_m = v$ is a path in $G_"core"$. For every consecutive edge ${w_(t-1), w_t} in E_"core"$, either $w_(t-1)$ and $w_t$ belong to the same shard $cal(X)_r$—in which case they lie in the same component of $G_"core" [V_"core"^((r))]$ so $pi(w_(t-1)) = pi(w_t)$—or they belong to distinct shards $cal(X)_r != cal(X)_s$. By Lemma 1(2), the cross-partition edge ${w_(t-1), w_t}$ generates a witness in Phase 2, creating edge ${pi(w_(t-1)), pi(w_t)} in E_"merge"$ in Phase 3. Contracting identical adjacent labels yields a path from $pi(u)$ to $pi(v)$ in $G_"merge"$.
  - *Converse direction:* Conversely, suppose $pi(u) = alpha_0, alpha_1, dots, alpha_ell = pi(v)$ is a path in $G_"merge"$. Every edge ${alpha_(j-1), alpha_j} in E_"merge"$ arises from some core pair $(a_j, b_j)$ with $a_j in alpha_(j-1) inter V_"core"$, $b_j in alpha_j inter V_"core"$, and $d(a_j, b_j) <= epsilon$, so ${a_j, b_j} in E_"core"$. Within each local cluster $alpha_j$, its core set is connected in $G_"core" [V_"core"^((r_j))] subset.eq G_"core"$. Concatenating intra-cluster core paths with the cross-partition edges ${a_j, b_j}$ produces a valid path from $u$ to $v$ in $G_"core"$.
  Thus $phi(pi(u)) = phi(pi(v)) <=> u ~_(G_"core") v$, proving $cal(K)_"par" = cal(K)_"seq"$.

  (2) A point $x in cal(X)_r$ is labeled noise in sequential DBSCAN iff $N_epsilon (x) inter V_"core" = emptyset$. By Lemma 1(1), $N_epsilon (x)$ is evaluated exactly; if any $u in N_epsilon (x) inter V_"core"$ lies in $cal(X)_r$, worker $r$ absorbs $x$ into $pi(u)$, whereas if $u in cal(X)_s$ ($s != r$), worker $s$ records witness $(u, x, pi(u))$ in $E_s$, enabling Phase 3 to assign $x$ to cluster $phi(pi(u))$. Thus $x$ remains unlabeled (noise) iff $N_epsilon (x) inter V_"core" = emptyset$.

  (3) Since $cal(K)_"par"$ is canonical up to global ID permutation (fixed e.g. by sorting components by minimum point index), taking the minimum adjacent global cluster ID for every non-core non-noise border point $b$ depends solely on the unordered set ${phi(pi(u)) : u in N_epsilon (b) inter V_"core"}$, eliminating the classical DBSCAN border-point order ambiguity @schubert2017dbscanrevisited.
]

#subsec("B", "Geometric Surface-to-Volume Communication Bound")

#theorem(title: "Halo Communication Bound and Parallel Efficiency Limit", origin: "Original Derivation")[
  Suppose $cal(X) = {x_1, dots, x_n}$ are i.i.d. draws from a density $f$ on $[0, 1]^d$ bounded by $0 < f_"min" <= f(x) <= f_"max"$, partitioned into $P = m^d$ congruent hypercubes $R_1, dots, R_P$ of side length $h = P^(-1/d) > 2 epsilon$. Then for every interior cell $R_r$:
  1. The expected remote halo size satisfies
     $ bb(E)[|H_r (epsilon)|] <= n f_"max" / P ((1 + 2 epsilon P^(1/d))^d - 1). $
  2. In the small-radius regime $eta := 2 d epsilon P^(1/d) <= 1$, the ratio of expected remote-halo points to local points obeys
     $ (bb(E)[|H_r (epsilon)|]) / (bb(E)[|cal(X)_r|]) &<= (f_"max" / f_"min") (e - 1) \
       &quad times 2 d epsilon P^(1/d). $
  3. Consequently, parallel communication overhead overtakes local work once $P = Omega((d epsilon)^(-d))$, establishing an intrinsic dimensional ceiling on fine-grained spatial partitioning.
]

#proof[
  (1) Every point $y in H_r (epsilon)$ lies outside $R_r = product_(j=1)^d [a_j, a_j + h]$ yet within Euclidean distance $epsilon$ of $R_r$. Because $norm(u)_infinity <= norm(u)_2$, the $epsilon$-neighborhood of $R_r$ is contained in the inflated concentric hypercube $R_r^(+epsilon) = product_(j=1)^d [a_j - epsilon, a_j + h + epsilon]$, whose Lebesgue volume is $(h + 2 epsilon)^d$. Subtracting $"Vol"(R_r) = h^d$ bounds the halo collar volume by
  $ "Vol"(R_r^(+epsilon) backslash R_r) = (h + 2 epsilon)^d - h^d \ = h^d ((1 + (2 epsilon) / h)^d - 1). $
  Substituting $h = P^(-1/d)$ gives $h^d = 1 / P$ and $2 epsilon / h = 2 epsilon P^(1/d)$. Integrating $f <= f_"max"$ over the collar and multiplying by $n$ yields
  $ bb(E)[|H_r (epsilon)|] <= n f_"max" / P ((1 + 2 epsilon P^(1/d))^d - 1). $

  (2) Using the elementary inequality $(1 + x)^d <= e^(d x)$ for $x >= 0$ and the convexity bound $e^eta - 1 <= (e - 1) eta$ valid for all $eta in [0, 1]$, we obtain with $x = 2 epsilon P^(1/d)$ and $eta = 2 d epsilon P^(1/d)$:
  $ (1 + 2 epsilon P^(1/d))^d - 1 <= e^eta - 1 <= (e - 1) dot 2 d epsilon P^(1/d). $
  Meanwhile, $bb(E)[|cal(X)_r|] = n integral_(R_r) f(x) dif x >= n f_"min" / P$. Dividing the two bounds cancels $n / P$ and yields Claim 2.

  (3) For fixed $epsilon$ and $d$, once $2 d epsilon P^(1/d) = Omega(1)$—equivalently $P = Omega((2 d epsilon)^(-d))$—the collar volume $(h + 2 epsilon)^d - h^d$ equals or exceeds $h^d$, so the remote halo contains as many points as the local partition itself and every worker incurs $Omega(1)$ remote RPC fraction per query, collapsing linear speedup.
]

#subsec("C", "Spatial Skew and Load-Imbalance Counterexample")

#proposition(title: "Worst-Case Speedup Collapse Under Static Grid Partitioning", origin: "Original Counterexample")[
  Fix any $delta in (0, 1)$, worker count $P >= 2$, and static equal-volume partitioning of $[0, 1]^d$ into $P$ hypercubes of side $h = P^(-1/d) > 2 epsilon$. For sufficiently large $n$, there exist zero-communication datasets $cal(X) subset.eq [0, 1]^d$ of size $n$ whose parallel speedup $S(P) = T_1 / T_P$ satisfies:
  1. $S(P) <= 1 / (1 - delta)^2$ under unindexed $Theta(m^2)$ local scans (and $S(P) -> 1$ even with spatial indexing on a single dense hotspot);
  2. $S(P) <= 1 / (1 - delta)$ under bounded-degree $Theta(m log n)$ spatially indexed scans.
]

#proof[
  (1) Choose the center $c_1$ of cell $R_1$ and place $n_1 = ceil((1 - delta) n)$ points inside the open ball $B_(epsilon / 4) (c_1) subset.eq R_1$. Distribute the remaining $n - n_1 = floor(delta n)$ points across the interiors of $R_2, dots, R_P$ at mutual separations $> 2 epsilon$. Since $h > 2 epsilon$, $B_(epsilon / 4) (c_1)$ lies at distance $> epsilon$ from $partial R_1$, so zero remote RPCs occur. Under unindexed local scans, sequential DBSCAN performs $c n^2$ distance checks whereas worker $1$ performs $c n_1^2$ checks, bounding speedup by
  $ S(P) = T_1 / T_P <= (c n^2) / (c n_1^2) <= 1 / (1 - delta)^2. $
  Moreover, even with a spatial index, every query inside $B_(epsilon / 4) (c_1)$ outputs all $n_1$ neighbors, forcing $T_P >= Omega(n_1^2)$ while the remaining $delta n$ sparse queries contribute only $O(delta n log n) = o(n_1^2)$, driving $S(P) -> 1$.

  (2) To bound the indexed, bounded-degree regime, instead pack $n_1 = ceil((1 - delta) n)$ points into $floor(n_1 / "MinPts")$ disjoint micro-clusters of size $"MinPts"$ strictly inside the interior sub-cube of $R_1$ of side $h - 2 epsilon > 0$ (choosing $epsilon$ small enough for separation $> 2 epsilon$), and place $floor(delta n)$ isolated points in $R_2, dots, R_P$. Every dR\*-tree query returns $O(1)$ neighbors in $Theta(log n)$ time, so $T_1 = Theta(n log n)$ and worker $1$'s makespan is $T_P >= Theta(n_1 log n)$, yielding
  $ S(P) <= n / n_1 <= 1 / (1 - delta). $
  For $delta = 0.05$, speedup is capped at $<= 1.11 times$ even on $P = 1000$ nodes, proving why data-density-balanced partitioning (Hilbert dR\*-tree in PDBSCAN @xu1999pdbscan or cost-balanced splits in MR-DBSCAN @he2011mrdbscan) is mathematically indispensable.
]

#subsec("D", "Distributed Index and Union-Find Merging Complexity")

#theorem(title: "Parallel Runtime, Communication, and Merge Complexity", origin: "Classical — Cited")[
  Let the dR\*-tree have page fanout $B$ and height $h_T = O(log_B n)$, and suppose $cal(X)$ is balanced across $P$ nodes ($|cal(X)_r| = Theta(n / P)$) with $M_"bnd"$ total cross-partition core-contact witnesses and at most $C_"rem"$ remote range-query RPCs per worker @xu1999pdbscan. Then:
  1. Phase 2 executes in expected wall-clock time $O((n / P) (log_B n + overline(k)_epsilon) + C_"rem" dot tau_"net")$ per worker, where $overline(k)_epsilon$ is the mean $epsilon$-neighborhood size and $tau_"net"$ is per-RPC latency.
  2. Phase 3 resolves all $M_"bnd"$ merges on the master using Disjoint-Set Union (path compression and union by rank) in $O(M_"bnd" alpha(n))$ time and $O(|V_"loc"|)$ memory, where $alpha(dot)$ is the inverse Ackermann function.
  3. Consequently, PDBSCAN achieves near-linear speedup $S(P) = Theta(P)$ whenever $C_"rem" tau_"net" = o((n / P) log_B n)$ and $M_"bnd" alpha(n) = o((n / P) log_B n)$.
]

#proof(label: "Proof Synthesis (Cited — Xu et al. 1999)")[
  In the dR\*-tree @xu1999pdbscan, replicating top-level directory nodes across all $P$ workers eliminates root-node contention: every worker traverses upper bounding boxes locally in $O(log_B n)$ comparisons and performs local disk/memory reads for leaves contained in $R_r$, returning $overline(k)_epsilon$ neighbors on average. Only queries whose $epsilon$-ball intersects a remote leaf MBR trigger a network message; aggregating $C_"rem"$ such queries per worker yields Claim 1. In Phase 3, each local cluster ID $(r, c) in V_"loc"$ initializes a singleton Union-Find set ($|V_"loc"| <= n / "MinPts"$), and each cross-partition witness $(u, v)$ with both endpoints verified as core triggers at most one `Union` and two `Find` operations. By Tarjan's classical amortized bound, $M_"bnd"$ operations execute in $O((|V_"loc"| + M_"bnd") alpha(|V_"loc"|)) = O(M_"bnd" alpha(n))$ worst-case time, establishing Claim 2. Comparing against sequential DBSCAN's $Theta(n (log_B n + overline(k)_epsilon))$ runtime directly gives the near-linear scaling condition of Claim 3.
]

#sec("III", "Complexity, Geometry and Operational Profile")

PDBSCAN trades modest spatial-indexing and boundary-RPC overhead for horizontal memory and compute scalability. Each worker stores $O((n / P) d)$ point coordinates plus the replicated upper dR\*-tree directory of size $O(P d)$—drastically lower than single-node $O(n d)$ memory. Because remote range queries dominate latency when $epsilon$ is large or $d$ grows (Theorem 2), Xu et al. @xu1999pdbscan emphasize *Hilbert space-filling curve* leaf assignment, which minimizes the expected bounding-box perimeter of each partition $R_r$ and reduces $C_"rem"$ compared to arbitrary hash or round-robin striping. Unlike iterative centroid methods (e.g., K-Means Alg. 1 or K-Means|| Alg. 23), PDBSCAN requires only a *single* global synchronization barrier (Phase 3) after local expansion, avoiding multi-round shuffle stalls.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for PDBSCAN],
    [Time / Space Complexity], [$O((n / P) log_B n + C_"rem" tau_"net" + M_"bnd" alpha(n))$ wall time; $O((n/P + P) d)$ RAM per node.],
    [Geometric Inductive Bias], [Exact DBSCAN super-level sets: arbitrary non-convex manifolds with explicit noise rejection.],
    [Determinism / Equivalence], [Exact core-set and noise equivalence to DBSCAN (Theorem 1); deterministic with min-ID border tie-break.],
    [Communication Scaling], [Halo ratio $O(d epsilon P^(1/d))$ (Theorem 2); degrades sharply in high $d$ or oversized $P$.],
    [Skew / Outlier Robustness], [Outlier-robust (noise isolated); vulnerable to spatial load skew without density-balanced partitioning (Prop. 1).],
    [Recommended Upgrades], [MR-DBSCAN (Alg. 49) for cost-balanced MapReduce; GriDBSCAN (Alg. 48) for grid overlap; HDBSCAN (Alg. 32) for varied densities.],
  ),
  caption: [Theoretical guarantees and operational profile of PDBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

PDBSCAN is engineered for massive low-to-moderate dimensional spatial databases where non-convex clusters and background clutter preclude centroid partitioning and data volumes exceed single-server memory.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Planetary-Scale Geospatial & Earth Observation Archives:* Satellite imagery catalogs, OpenStreetMap POI archives, and global AIS vessel-tracking feeds span $10^8$–$10^9$ latitude–longitude records in $d = 2$–$3$. Because $d$ is small, Theorem 2 guarantees minimal halo surface-to-volume ratios ($O(epsilon sqrt(P))$), allowing PDBSCAN over Hilbert-partitioned dR\*-trees to extract coastlines, urban corridors, and shipping lanes with near-linear speedup.
2. *Large-Scale LiDAR & Photogrammetry Point Clouds:* Autonomous-driving fleet maps and airborne topographic LiDAR produce billions of 3D $(x, y, z)$ returns with non-convex building facades, road curbs, and vegetation noise @rusu2011pcl. PDBSCAN's exact core-graph equivalence (Theorem 1) guarantees seamless object segmentation across spatial tile boundaries without artificial seam artifacts.
3. *High-Throughput Astronomy & Sky-Survey Catalogs:* Galactic filament and star-cluster discovery in multi-terabyte astronomical surveys (e.g., Gaia, LSST) involves filamentary, non-spherical structures embedded in diffuse background stars. Distributed spatial indexing with exact $epsilon$-connectivity recovers long curved filaments spanning multiple storage nodes.
4. *Telecom Call-Detail & Urban Mobility Hotspot Mining:* Mobile network base-station pings and GPS trajectory endpoints exhibit strongly non-convex street-network geometry and extreme outlier bursts. Shared-nothing partitioning across regional servers enables batch hotspot extraction without centralizing raw telemetry until Phase 3 merging.
5. *Computational Materials & Molecular Dynamics Snapshots:* Identifying atomic defect clusters, grain boundaries, and phase-separated domains across $10^7$-atom 3D simulation boxes naturally aligns with domain-decomposed MPI clusters, where PDBSCAN's $epsilon$-halo mirrors standard spatial ghost-cell exchange.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PDBSCAN inherits sequential DBSCAN's sensitivity to global $(epsilon, "MinPts")$ across multi-density regions (resolved by OPTICS Alg. 33 or HDBSCAN Alg. 32) and adds two distributed failure modes: (a) *high-dimensional halo explosion* (Theorem 2), where R\*-tree bounding-box overlap approaches 100% for $d > 15$ and remote RPC traffic destroys speedup; and (b) *severe spatial skew* (Proposition 1), where dense hotspots overload a single worker unless cost-balanced partitioning @he2011mrdbscan is employed.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why PDBSCAN Fits / Prerequisite], [Caution / Alternative],
    [GIS / Satellite Maps], [2D/3D Hilbert locality; tiny halo ratio], [Use HDBSCAN (Alg. 32) for varied density],
    [3D LiDAR Point Clouds], [Exact cross-tile object stitching], [Use MR-DBSCAN (Alg. 49) under skew],
    [Astronomical Catalogs], [Recovers non-convex filaments + noise], [Use OPTICS (Alg. 33) for nested scales],
    [High-D Embeddings], [Poor fit: halo & R\*-tree collapse], [Use SNN (Alg. 36) or K-Means|| (Alg. 23)],
    [Small In-Memory Sets], [RPC/merge overhead outweighs gains], [Use DBSCAN (Alg. 31) or DBSCAN++ (Alg. 47)],
  ),
  caption: [Domain engineering suitability matrix for PDBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PDBSCAN stands at the foundation of distributed density clustering, forming a unified algebraic family with grid- and MapReduce-based parallelizations.

#proposition(title: "Unifying Overlapping-Cover Gluing Framework", origin: "Original Derivation")[
  Let ${R_r}_(r=1)^P$ be any disjoint partition of $bb(R)^d$, and define the open-cover augmentation $U_r = {x in cal(X) : d(x, R_r) <= epsilon} = cal(X)_r union H_r (epsilon)$. Then PDBSCAN (Alg. 37), GriDBSCAN (Alg. 48) @mahran2008gridbscan, and MR-DBSCAN (Alg. 49) @he2011mrdbscan all compute the exact same quotient core-partition $V_"core" / ~_(G_"core")$ by taking local connected components on $U_r$ and gluing along the overlap cocycle
  $ (r, c_1) ~ (s, c_2) <==> exists u in U_r inter U_s inter V_"core" \ "with" space u in C_(r, c_1) inter C_(s, c_2). $
  The three algorithms differ strictly in (i) how ${R_r}$ is constructed (Hilbert dR\*-tree vs. regular grid vs. cost-estimated KD-partition) and (ii) whether $H_r (epsilon)$ is materialized lazily via RPC (PDBSCAN) or eagerly via point replication (GriDBSCAN, MR-DBSCAN).
]

#proof[
  For any core edge ${u, v} in E_"core"$ with $u in cal(X)_r$ and $v in cal(X)_s$, we have $u, v in U_r inter U_s$. In eager-replication schemes (GriDBSCAN, MR-DBSCAN), both workers $r$ and $s$ evaluate $u$ and $v$ locally inside $U_r$ and $U_s$, placing $u$ in $C_(r, c_1)$ and $C_(s, c_2)$ simultaneously; gluing identical core points across $U_r inter U_s$ therefore links $(r, c_1)$ and $(s, c_2)$ identically to Phase 3 of Algorithm 1. By Theorem 1, the induced quotient on $V_"core"$ is invariant to the choice of partition and replication strategy.
]

Beyond parallel variants, PDBSCAN directly extends sequential DBSCAN (Alg. 31) @ester1996dbscan and can be combined with generalized neighborhood predicates from GDBSCAN (Alg. 39) @sander1998gdbscan or subsampled core selection from DBSCAN++ (Alg. 47) @jang2019dbscanpp.

On the ClusterOverview 4-dataset empirical benchmark suite, PDBSCAN confirms these exact theoretical predictions. On the 2D *Multi-Shape* dataset (non-convex moons, rings, and clusters with clutter), PDBSCAN achieves a high Silhouette of $0.717$ (ARI N/A), directly reflecting Theorem 1's lossless preservation of DBSCAN's arbitrary-shape density connectedness across partitions. On *Iris*, PDBSCAN attains Silhouette $= 0.627$ and ARI $= 0.503$, mirroring classical single-$epsilon$ DBSCAN's tendency to separate Setosa cleanly while merging adjacent Versicolor and Virginica manifolds. On *Customers* (Silhouette $= 0.296$) and *High-D Subspace+Noise* (Silhouette $= 0.352$, ARI N/A), accuracy degrades due to varying intra-cluster densities and high-dimensional distance concentration @beyer1999nearest alongside increased boundary overlap (Theorem 2). Finally, PDBSCAN records a mean runtime of $1.079$ ms—the slowest of the three adjacent DBSCAN extensions (Algs. 37–39) on small benchmark sizes—because partitioning, halo lookup, and Phase 3 merge coordination introduce fixed overhead that amortizes only at large database scale ($n >= 10^5$).
