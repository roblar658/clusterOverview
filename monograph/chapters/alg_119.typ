#import "../template/ieee.typ": *

#alg_entry(
  id: 119,
  name: "Fluid Communities (FluidC)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2018",
  paradigm: "Density-Inhibited Asynchronous Propagation",
  complexity: $cal(O)(T (|V| + |E|))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Community detection partitions a graph into densely internally connected, sparsely externally connected groups. While modularity @newman2004finding @blondel2008louvain and flow compression @rosvall2008infomap choose the community count $k$ implicitly, many applications require an explicit $k$-partition at linear cost. Unconstrained asynchronous Label Propagation (LPA) @raghavan2007lpa runs in $cal(O)(|V| + |E|)$ per pass but often collapses into a single "monster community" on graphs with weak bottlenecks. *Fluid Communities* (FluidC) @pares2018fluidc resolves this degeneracy by modeling $k$ competing incompressible fluids that expand and contract across the graph until local neighborhood pressures equilibrate.

Let $G = (V, E)$ be a simple, connected, undirected graph with $n = |V|$ vertices, $m = |E|$ edges, and prescribed $k in {1, dots.c, n}$. For $v in V$, let $N(v) = {u in V : {u, v} in E}$, $N[v] = N(v) union {v}$, and $"deg"(v) = |N(v)|$. At asynchronous micro-step $t >= 0$, the labeling $ell^((t)) : V -> {0, 1, dots.c, k}$ ($0$ = uncolonized) induces fluid supports $C_c^((t)) = {u in V : ell^((t)) (u) = c}$.

#definition(title: "Fluid Density and Neighborhood Update Potential")[
  Each active fluid $c in {1, dots.c, k}$ with non-empty support $C_c^((t)) != emptyset$ carries a unit total mass distributed uniformly over its vertices, defining its *inverse-cardinality density*
  $ rho_c^((t)) = 1 / (|C_c^((t))|) in (0, 1], quad "with" quad rho_0^((t)) = 0. $
  Consequently, total fluid mass satisfies the conservation law $sum_(v in C_c^((t))) rho_c^((t)) = 1$ for every active $c$, and $sum_(c=1)^k sum_(v in C_c^((t))) rho_c^((t)) = k$. For any vertex $v in V$ and fluid $c in {1, dots.c, k}$, the *closed-neighborhood density support* is
  $ S^((t)) (v, c) &= sum_(u in N[v]) rho_c^((t)) bb(I)[u in C_c^((t))] \
  &= (|N[v] inter C_c^((t))|) / (|C_c^((t))|) in [0, 1]. $
]

Unlike LPA's raw count $|N(v) inter C_c|$, FluidC normalizes neighborhood support by global size $|C_c|$: expanding fluids dilute ($rho_c = 1 / |C_c|$ decreases), while contracting fluids concentrate and resist encroachment. Each asynchronous superstep traverses $V$ in a uniformly random permutation $pi$, updating $v$ to a maximizer in $M^((t)) (v) = arg max_(c > 0) S^((t)) (v, c)$ while retaining $ell^((t)) (v)$ whenever $ell^((t)) (v) in M^((t)) (v)$.

#algbox(
  title: "Fluid Communities (FluidC)",
  inputs: [Connected graph $G = (V, E)$, community count $k in {1, dots.c, n}$, superstep cap $T_"max"$.],
  outputs: [$k$-partition $cal(C) = {C_1, dots.c, C_k}$ of $V$.],
)[
  1. Sample $k$ distinct seed vertices $s_1, dots.c, s_k in V$ uniformly without replacement; set $ell(s_c) <- c$, $|C_c| <- 1$, $rho_c <- 1$ for $c = 1, dots.c, k$, and $ell(v) <- 0$ for all remaining $v in V$.
  2. For superstep $r = 1, dots.c, T_"max"$:
     a. Draw a uniform random permutation $pi$ of $V$; set `"changed" <- "false"`.
     b. For each $v$ in order $pi$:
        - Compute $S(v, c) = |N[v] inter C_c| / |C_c|$ for every active fluid $c > 0$ appearing in $N[v]$.
        - If no active fluid appears in $N[v]$, leave $ell(v) = 0$ unchanged and continue.
        - Let $M(v) = arg max_(c > 0) S(v, c)$.
        - If $ell(v) in.not M(v)$: sample $c' tilde.op "Unif"(M(v))$; if $ell(v) > 0$, decrement $|C_(ell(v))|$ and update $rho_(ell(v))$; increment $|C_(c')|$, update $rho_(c')$, set $ell(v) <- c'$ and `"changed" <- "true"`.
     c. If `"changed" = "false"`, terminate and return $cal(C) = {C_1, dots.c, C_k}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal properties governing FluidC: exact non-extinction and coverage time, its exact variational relation (and non-monotonicity) relative to Ratio Association, a complete equilibrium and disconnected-community pathology on star graphs, a sharp Chernoff–Hellinger stability threshold under the Stochastic Block Model, and linear per-superstep complexity.

#subsec("A", "Non-Extinction Invariant and Diameter-Bounded Coverage")

#proposition(title: "Strict Non-Extinction and Exact k-Partition Coverage", origin: "Original Proof")[
  Let $G = (V, E)$ be connected and $1 <= k <= n$. Under the tie-retaining closed-neighborhood FluidC dynamics of Algorithm 1:
  1. Every community remains strictly non-empty at every micro-step: $|C_c^((t))| >= 1$ for all $c in {1, dots.c, k}$ and all $t >= 0$.
  2. Every vertex acquires a non-zero fluid label within at most $"diam"(G)$ supersteps, producing a valid $k$-partition of $V$.
]

#proof[
  *(1) Non-extinction.* We induct on micro-step $t$. At $t = 0$, the $k$ distinct seeds satisfy $|C_c^((0))| = 1$ for every $c$. Assume $|C_c^((t))| >= 1$ for all $c$. Only an update of a vertex $v in C_c^((t))$ can reduce $|C_c^((t))|$, and extinction requires $|C_c^((t))| = 1$, i.e., $C_c^((t)) = {v}$. Because the support is evaluated over the *closed* neighborhood $N[v] = N(v) union {v}$, we have $v in N[v] inter C_c^((t))$, yielding
  $ S^((t)) (v, c) = (|N[v] inter {v}|) / (|{v}|) = 1 / 1 = 1. $
  For any competing fluid $c' != c$, we have $v in.not C_(c')^((t))$, hence $N[v] inter C_(c')^((t)) = N(v) inter C_(c')^((t)) subset.eq C_(c')^((t))$, which implies
  $ S^((t)) (v, c') &= (|N(v) inter C_(c')^((t))|) / (|C_(c')^((t))|) \
  &<= (|C_(c')^((t))|) / (|C_(c')^((t))|) = 1 = S^((t)) (v, c). $
  Therefore $c in M^((t)) (v) = arg max_(c') S^((t)) (v, c')$. By the current-label retention rule, $v$ retains label $c$, so $|C_c^((t+1))| = 1$. By induction, no fluid ever vanishes.

  *(2) Coverage time.* Let $U_r = {v in V : ell(v) > 0}$ after superstep $r$, with $U_0 = {s_1, dots.c, s_k}$. Once $ell(v) = c > 0$, $S(v, c) >= 1 / n > 0 = S(v, 0)$ precludes reverting to $0$, so $U_0 subset.eq U_1 subset.eq dots.c$. We induct on $r >= 0$ that $B_r = {v in V : "dist"_G (v, U_0) <= r} subset.eq U_r$. The base $r = 0$ is immediate. If $B_(r-1) subset.eq U_(r-1)$ and $"dist"_G (v, U_0) = r$, then $v$ has a neighbor $u in B_(r-1) subset.eq U_(r-1)$ holding an active label throughout superstep $r$. When $v$ is visited in superstep $r$, $N[v]$ contains $u$, ensuring $max_(c > 0) S(v, c) >= 1 / n > 0$ and forcing $v in U_r$. Since $G$ is connected, $B_("diam"(G)) = V$, hence $U_("diam"(G)) = V$.
]

#subsec("B", "Exact Relation to Ratio Association and Non-Monotonicity")

A natural question is whether FluidC performs coordinate ascent on a global graph-partitioning potential. We connect its local transition criterion exactly to the *Ratio Association* objective @dhillon2004kernel (the algebraic dual of Ratio Cut @hagen1992ratiocut) and prove that FluidC omits a global recoil term, causing strict non-monotonicity.

#theorem(title: "Ratio-Association Recoil Decomposition and Non-Monotonicity", origin: "Original Derivation")[
  For a $k$-partition $cal(C) = {C_1, dots.c, C_k}$, define the *Ratio Association* functional and internal average degree
  $ Psi(cal(C)) = sum_(c=1)^k (2 |E(C_c)|) / (|C_c|), quad bar(d)_c = (2 |E(C_c)|) / (|C_c|), $
  where $E(C_c) = {{u, w} in E : u, w in C_c}$. When vertex $v in C_a$ (with $n_a = |C_a| >= 2$) moves to $C_b$ ($n_b = |C_b| >= 1$) with $d_c (v) = |N(v) inter C_c|$, the exact change $Delta Psi = Psi(cal(C)') - Psi(cal(C))$ satisfies
  $ Delta Psi &= 2 ( (d_b (v)) / (n_b + 1) - (d_a (v)) / (n_a - 1) ) \
  &quad + (bar(d)_a) / (n_a - 1) - (bar(d)_b) / (n_b + 1). $
  Moreover:
  1. If $n_a = n_b + 1$ and $bar(d)_a >= bar(d)_b$, every FluidC move strictly satisfies $Delta Psi > 0$.
  2. In general, FluidC does *not* monotonically increase $Psi$: there exist graphs and states where a valid FluidC move yields $Delta Psi < 0$.
]

#proof[
  Let $e_c = |E(C_c)|$. Moving $v$ from $C_a$ to $C_b$ transforms $(n_a, e_a) |-> (n_a - 1, e_a - d_a (v))$ and $(n_b, e_b) |-> (n_b + 1, e_b + d_b (v))$, leaving all other communities unchanged. Expanding $Delta Psi / 2$ directly gives
  $ (Delta Psi) / 2 &= (e_a - d_a (v)) / (n_a - 1) - (e_a) / (n_a) + (e_b + d_b (v)) / (n_b + 1) - (e_b) / (n_b) \
  &= (e_a - n_a d_a (v)) / (n_a (n_a - 1)) + (n_b d_b (v) - e_b) / (n_b (n_b + 1)) \
  &= (d_b (v)) / (n_b + 1) - (d_a (v)) / (n_a - 1) + (e_a) / (n_a (n_a - 1)) - (e_b) / (n_b (n_b + 1)). $
  Multiplying by $2$ and substituting $2 e_c / n_c = bar(d)_c$ proves the exact identity.

  *(1)* Suppose $n_a = n_b + 1$ and $bar(d)_a >= bar(d)_b$. The FluidC move criterion $S(v, b) > S(v, a)$ reads $d_b (v) / n_b > (d_a (v) + 1) / (n_b + 1)$, i.e.,
  $ n_b (d_b (v) - d_a (v)) > n_b - d_b (v) >= 0, $
  where $d_b (v) <= n_b$. By integrality, $d_b (v) >= d_a (v) + 1$. Substituting $n_a - 1 = n_b$ into the first bracket gives
  $ (d_b (v)) / (n_b + 1) - (d_a (v)) / (n_b) &>= (d_a (v) + 1) / (n_b + 1) - (d_a (v)) / (n_b) \
  &= (n_b - d_a (v)) / (n_b (n_b + 1)) >= 0, $
  since $d_a (v) <= |C_a without {v}| = n_b$. Meanwhile, $bar(d)_a / (n_a - 1) - bar(d)_b / (n_b + 1) >= bar(d)_b (1 / n_b - 1 / (n_b + 1)) > 0$ (or $d_a (v) < n_b$ if $bar(d)_b = 0$ since $d_b (v) >= 1$), yielding $Delta Psi > 0$.

  *(2) Constructive counterexample.* Take $13$ vertices split into a path $C_a = P_9$ ($n_a = 9$, $e_a = 8$, $bar(d)_a = 16 / 9$, endpoint $v$ with $d_a (v) = 1$) and a clique $C_b = K_4$ ($n_b = 4$, $e_b = 6$, $bar(d)_b = 3$), joined by one bridge edge from $v$ ($d_b (v) = 1$). Then $S(v, b) = 1 / 4 > 2 / 9 = S(v, a)$, so FluidC moves $v$ to $C_b$, yet
  $ Delta Psi &= 2 (1 / 5 - 1 / 8) + (16 / 9) / 8 - 3 / 5 \
  &= 3 / 20 + 2 / 9 - 12 / 20 = - 41 / 180 < 0. $
  FluidC dilutes $C_b$ by ignoring the global recoil term $-bar(d)_b / (n_b + 1)$.
]

#subsec("C", "Star-Graph Equilibrium and Disconnected Communities")

Although FluidC prevents single-community collapse, density equalization over closed neighborhoods can produce completely disconnected communities around high-degree hubs.

#theorem(title: "Exact Star-Graph Equilibria and Disconnected-Leaf Pathology", origin: "Original Counterexample")[
  Consider the star graph $K_(1, n-1)$ ($n >= 4$) with hub $v_0$ and $n - 1$ leaves, partitioned by FluidC into $k >= 2$ communities. Then:
  1. A $k$-partition $cal(C)$ with $v_0 in C_(c^*)$ is a stationary fixed point of FluidC if and only if $|C_(c^*)| >= max_(j != c^*) |C_j|$.
  2. In every fixed point, all $k - 1$ non-hub communities $C_j$ ($j != c^*$) induce edgeless subgraphs ($|E(C_j)| = 0$, disconnected whenever $|C_j| >= 2$) while $|C_(c^*)| >= ceil(n / k)$.
  3. From any fully labeled $k$-partition, FluidC converges in at most $1$ superstep (and from $k$-seed initialization in at most $2$ supersteps).
]

#proof[
  Let $v_0 in C_(c^*)$ with $r = |C_(c^*)|$.
  - *Hub invariance:* Because $N[v_0] = V$, every active fluid $c in {1, dots.c, k}$ attains $S(v_0, c) = |C_c| / |C_c| = 1$. By tie-retention, $v_0$ never leaves $C_(c^*)$.
  - *Hub-community leaves:* For any leaf $u in C_(c^*) without {v_0}$, $N[u] = {v_0, u} subset.eq C_(c^*)$, so $S(u, c^*) = 2 / r > 0 = S(u, j)$ for all $j != c^*$. Thus no vertex ever leaves $C_(c^*)$, meaning $r$ is monotonically non-decreasing.
  - *External leaves:* For any leaf $w in C_j$ ($j != c^*$), $N[w] = {v_0, w\}$, giving $S(w, c^*) = 1 / r$, $S(w, j) = 1 / |C_j|$, and $S(w, c') = 0$ for $c' in.not {c^*, j}$. Hence $w$ switches to $C_(c^*)$ if and only if $1 / r > 1 / |C_j|$, i.e., $|C_j| > r$.

  *(1)–(2)* Consequently, no vertex can move if and only if $|C_j| <= r$ for all $j != c^*$. Summing over all $k$ parts gives $n = r + sum_(j != c^*) |C_j| <= k r$, forcing $r >= ceil(n / k)$. Because every edge of $K_(1, n-1)$ is incident to $v_0 in C_(c^*)$, every non-hub community $C_j$ contains zero internal edges.

  *(3)* Starting from any fully labeled $k$-partition, suppose for contradiction that after one full superstep some $j != c^*$ still satisfies $|C_j^("end")| > r^("end")$. Because $r$ is non-decreasing and $|C_j|$ is non-increasing during the superstep, at the moment the last vertex $w in C_j$ was visited we had $|C_j^("visit")| >= |C_j^("end")| > r^("end") >= r^("visit")$, forcing $w$ to transfer to $C_(c^*)$, contradiction. From $k$-seed initialization, $v_0$ is labeled in superstep $1$ (if not already a seed, its visit finds all $k$ singleton seeds tied at $1$); in superstep $2$, every remaining unlabeled leaf has only $v_0$ as a labeled neighbor and joins $C_(c^*)$, leaving every other seed as a singleton of size $1 <= r$, which is already stationary.
]

#subsec("D", "Sharp SBM Stability Threshold and Linear Complexity")

Despite pathological hubs, FluidC exhibits provably sharp statistical stability on assortative modular networks matching the information-theoretic Chernoff–Hellinger recovery scaling of the Stochastic Block Model (SBM) @holland1983sbm @abbe2018sbmreview.

#theorem(title: "Planted-Partition Stability under the Balanced SBM", origin: "Original Proof")[
  Let $G tilde.op "SBM"(n, k, p, q)$ be a balanced $k$-block Stochastic Block Model with equal community sizes $s = n / k$, intra-block edge probability $p$, and inter-block probability $q < p$. Whenever $(s - 1) p > s q$, the planted partition $cal(C)^*$ is an exact stationary fixed point of FluidC with probability at least
  $ 1 - n (k - 1) exp( - (sqrt((s - 1) p) - sqrt(s q))^2 ). $
  In the logarithmic connectivity regime $p = a (ln n) / n$, $q = b (ln n) / n$, $cal(C)^*$ is asymptotically stable with probability $1 - o(1)$ whenever $sqrt(a) - sqrt(b) > sqrt(k)$.
]

#proof[
  Under $cal(C)^*$, all blocks have equal cardinality $|C_c^*| = s$. Fix a vertex $v in C_a^*$ and any rival block $b != a$. Vertex $v$ prefers $C_b^*$ over $C_a^*$ if and only if
  $ S(v, b) > S(v, a) <==> (d_b (v)) / s > (d_a (v) + 1) / s <==> Z >= 2, $
  where $Z = d_b (v) - d_a (v) = sum_(u in C_b^*) X_u - sum_(w in C_a^* without {v}) Y_w$, with independent $X_u tilde.op "Ber"(q)$ and $Y_w tilde.op "Ber"(p)$. Let $A = (s - 1) p$ and $B = s q$, with $A > B > 0$. For any $lambda > 0$, independence and $1 + x <= e^x$ yield the moment-generating bound
  $ bb(E)[e^(lambda Z)] &= (1 - q + q e^lambda)^s (1 - p + p e^(-lambda))^(s - 1) \
  &<= exp( B (e^lambda - 1) + A (e^(-lambda) - 1) ). $
  Choosing the exact minimizer $e^lambda = sqrt(A / B) > 1$ evaluates the exponent to
  $ B sqrt(A / B) - B + A sqrt(B / A) - A = - (sqrt(A) - sqrt(B))^2. $
  By Markov's inequality and $e^(-2 lambda) = B / A < 1$,
  $ "Pr"[Z >= 2] <= e^(-2 lambda) bb(E)[e^(lambda Z)] <= exp( - (sqrt(A) - sqrt(B))^2 ). $
  Applying a union bound over all $n$ vertices and $k - 1$ rival blocks gives the finite-sample guarantee. Setting $p = a (ln n) / n$, $q = b (ln n) / n$, and $s = n / k$ yields $(sqrt(A) - sqrt(B))^2 = ((sqrt(a) - sqrt(b))^2 / k) ln n - cal(O)((ln n) / n)$, so the failure probability scales as $(k - 1) n^(1 - (sqrt(a) - sqrt(b))^2 / k + o(1)) -> 0$ strictly when $sqrt(a) - sqrt(b) > sqrt(k)$.
]

#proposition(title: "Exact Linear Per-Superstep Complexity", origin: "Classical — Cited")[
  Using an adjacency-list representation and an epoch-stamped density accumulator array of length $k$, each superstep of FluidC executes in $cal(O)(|V| + |E|)$ time and $cal(O)(|V| + k)$ auxiliary space, matching asynchronous LPA @raghavan2007lpa @pares2018fluidc.
]

#proof[
  Maintain integer arrays `count[1..k]`, `stamp[1..k]`, and `size[1..k]` initialized once in $cal(O)(k)$ time. When processing vertex $v$ at micro-step $t$, iterate over $u in N[v]$: for each $c = ell(u) > 0$, if `stamp[c] != t`, set `stamp[c] <- t`, `count[c] <- 1`, and append $c$ to a local active list $L_v$; otherwise increment `count[c]`. Computing $S(v, c) = "count"[c] / "size"[c]$ and finding $M(v)$ over $c in L_v$ takes $cal(O)(|L_v|) = cal(O)("deg"(v) + 1)$ time, independent of $k$. Updating `size` takes $cal(O)(1)$. Summing across all $v in V$ yields $sum_(v in V) cal(O)("deg"(v) + 1) = cal(O)(|V| + |E|)$ per superstep, with $cal(O)(|V| + k)$ auxiliary memory for labels, permutation buffer, and epoch arrays.
]

#sec("III", "Complexity, Geometry and Operational Profile")

FluidC achieves one of the lowest computational footprints in graph clustering: $cal(O)(T (|V| + |E|))$ time over $T$ supersteps (typically $T in [5, 25]$ on sparse modular graphs @pares2018fluidc) and $cal(O)(|V| + k)$ auxiliary space beyond $cal(O)(|V| + |E|)$ CSR storage. Random vertex sweeps cause irregular neighbor reads, making cache bandwidth—not arithmetic—the primary bottleneck above L3 cache. Geometrically, FluidC drives *equalizing geodesic wavefronts* from $k$ seeds whose expansion speed scales inversely with volume $1 / |C_c|$, migrating colliding boundaries toward sparse topological cuts without eigensolves or matrix factorizations. Its output is stochastic (depending on seeds, permutations, and tie breaks) and requires connected input graphs.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(T (|V| + |E|))$ total; $cal(O)(|V| + |E|)$ per superstep],
    [Space Complexity], [$cal(O)(|V| + |E| + k)$ including CSR graph and epoch buffers],
    [Community Count $k$], [Explicitly user-specified ($1 <= k <= |V|$); guaranteed exact],
    [Deterministic], [No; depends on seeds, sweep order $pi$, and tie breaks],
    [Connectivity Req.], [Connected $G$ required; disconnected components need sub-budgets],
    [Objective Geometry], [Local approximation of Ratio Association; density-inhibited fronts],
    [Collapse Resistance], [Provably non-extinct ($|C_c| >= 1$); immune to LPA monster cluster],
  ),
  caption: [Operational and complexity profile of Fluid Communities (FluidC).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Social and Collaboration Graphs:* Small-world networks have diameter $cal(O)(log n)$, so Proposition 1 guarantees full colonization in $cal(O)(log n)$ linear passes at fixed $k$—avoiding LPA collapse and modularity's resolution limit @fortunato2007resolution.
2. *Manifold Clustering on Sparse $k$-NN Graphs:* On non-convex or intertwined manifolds, symmetric $k$-NN graphs convert geometry into local connectivity; fluid fronts propagate smoothly along manifold geodesics without Euclidean centroid assumptions.
3. *Decentralized Sensor and Mesh Partitioning:* Each update requires only $1$-hop neighbor labels plus $k$ global cardinality scalars $|C_c|$, fitting distributed wireless networks lacking coordinates while maintaining self-regulating territory balance.
4. *Fixed-Granularity Biological Module Extraction:* On PPI and co-expression networks queried at prescribed pathway counts $k$, non-extinction ensures exactly $k$ non-empty modules at linear cost, enabling rapid multi-seed consensus clustering @monti2003consensus.
5. *Interactive $k$-Sweeps and Warm Initialization:* Sub-millisecond linear passes support sweeping $k$ via silhouette/ARI curves @rousseeuw1987silhouettes @hubert1985comparing and fast warm starts for Louvain @blondel2008louvain, Leiden @traag2019leiden, or Kernighan–Lin (`Alg. 120`).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

FluidC exhibits four primary failure modes: (i) *Hub-induced disconnection* on star-like or extreme power-law networks (Theorem 2); (ii) *Size-imbalance distortion*, since $rho_c = 1 / |C_c|$ over-equalizes disparate cluster sizes and ignores global recoil (Theorem 1); (iii) *Seed collision variance*, where two seeds in one dense cluster split it permanently unless mitigated by multi-start consensus; and (iv) *Disconnected graphs and unknown $k$*, requiring component-wise $k$ allocation or external model selection.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Sparse $k$-NN Manifolds], [Geodesic fronts track non-convex shapes], [Use multi-seed consensus; spectral @ng2002spectral if tiny $n$],
    [Assortative SBM / Modular Nets], [Sharp stability at $sqrt(a) - sqrt(b) > sqrt(k)$ (Thm. 3)], [Unknown $k$: Leiden @traag2019leiden or Infomap @rosvall2008infomap],
    [Hub / Star-Like Scale-Free Nets], [Poor; edgeless leaf communities (Thm. 2)], [Leiden @traag2019leiden or DC-SBM @karrer2011dcsbm],
    [Skewed Cluster Cardinalities], [Density $1 / |C_c|$ over-equalizes cluster sizes], [Louvain @blondel2008louvain or HDBSCAN @campello2013hdbscan],
    [Distributed / Linear Sharding], [Local $1$-hop rule + $k$ scalars; $cal(O)(|E|)$ pass], [For minimum cut use METIS (`Alg. 121`)],
  ),
  caption: [Domain suitability and failure-mode matrix for FluidC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

FluidC sits midway between unconstrained label propagation @raghavan2007lpa and balanced graph partitioning (`Alg. 120`–`Alg. 122`).

#proposition(title: "Exact Reduction of FluidC to Self-Looped Label Propagation", origin: "Original Proof")[
  On any balanced $k$-partition $cal(C)$ with $|C_1| = dots.c = |C_k| = n / k$, FluidC's decision set coincides identically with self-looped asynchronous Label Propagation @raghavan2007lpa:
  $ arg max_(c in {1, dots.c, k}) S(v, c) = arg max_(c in {1, dots.c, k}) |N[v] inter C_c|. $
  For unequal sizes, $S(v, c)$ acts as an inverse-size multiplier transferring boundary vertices from oversized to undersized communities whenever $|N[v] inter C_b| / |N[v] inter C_a| > |C_b| / |C_a|$.
]

#proof[
  When $|C_c| = n / k$ for all $c$, the positive constant $k / n$ factors out of $arg max_c (k / n) |N[v] inter C_c|$. For general sizes, $S(v, b) > S(v, a)$ holds iff $|N[v] inter C_b| / |C_b| > |N[v] inter C_a| / |C_a|$, yielding the stated ratio criterion when $|N[v] inter C_a| > 0$.
]

Relative to Chinese Whispers @biemann2006chinese and LPA @raghavan2007lpa, FluidC trades automatic $k$-selection for guaranteed non-extinction (Proposition 1); relative to Kernighan–Lin (`Alg. 120` @kernighan1970efficient), METIS (`Alg. 121` @karypis1998metis), and Chaco (`Alg. 122` @hendrickson1995chaco), it replaces hard balance constraints and eigensolves with a soft $1 / |C_c|$ restoring pressure.

*Empirical Behavior on the ClusterOverview Benchmark.* Across the four standardized ClusterOverview datasets (evaluated on mutual $k$-NN graphs via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), FluidC achieves: *Iris* $"Sil" = 0.476$, $"ARI" = 0.450$; *Customers* $"Sil" = 0.370$; *Multi-Shape* $"Sil" = 0.708$ ($"ARI" =$ N/A); *High-D* $"Sil" = 0.300$ ($"ARI" =$ N/A); mean runtime *0.131 ms*. These outcomes closely match Sections II–III: at *0.131 ms*, FluidC is $4.0 times$ faster than METIS (*0.524 ms*), $22.2 times$ faster than Chaco (*2.906 ms*), and $42.3 times$ faster than KL (*5.536 ms*), confirming its low-constant linear pass (Proposition 2). It leads all four partitioners in silhouette on *Iris* ($0.476$) and *Customers* ($0.370$) and achieves strong non-convex separation on *Multi-Shape* ($0.708$), while trailing METIS in *Iris* $"ARI"$ ($0.450$ vs. $0.678$) due to stochastic seeding and non-monotonic boundary drift (Theorem 1). On *High-D* ($"Sil" = 0.300$), $k$-NN hubness partially triggers the leaf-fragmentation mechanism of Theorem 2 while still outperforming METIS ($0.158$) and Chaco ($0.221$).
