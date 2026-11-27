#import "../template/ieee.typ": *

#alg_entry(
  id: 117,
  name: "Spinglass Community Detection",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2006",
  paradigm: "Statistical Mechanics Potts Hamiltonian Optimization",
  complexity: $cal(O)(N_("sw") (m + n q))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Reichardt and Bornholdt @reichardt2006spinglass recast community detection as finding the zero-temperature ground state of an infinite-range $q$-state Potts spin glass defined on an undirected graph $G = (V, E, w)$ with $n = |V|$ vertices, $m = |E|$ edges (total weight $m = 1 / 2 sum_(i=1)^n k_i$), adjacency matrix $A in bb(R)^(n times n)$ ($A_(i i) = 0$), and vertex degrees $k_i = sum_j A_(i j)$. Each vertex $i in V$ carries a discrete spin state $sigma_i in {1, dots, q}$ indexing its community assignment, forming a global configuration $bold(sigma) = (sigma_1, dots, sigma_n) in {1, dots, q}^n$.

Guided by four first-principles axioms—reward existing intra-community edges, penalize missing intra-community edges, penalize inter-community edges, and reward missing inter-community edges—Reichardt and Bornholdt @reichardt2006spinglass showed that up to an additive constant, any pairwise spin Hamiltonian reduces to
$ cal(H)(bold(sigma)) = -sum_(i < j) (a_(i j) A_(i j) - b_(i j) (1 - A_(i j))) delta(sigma_i, sigma_j), $
where $a_(i j), b_(i j) >= 0$ specify ferromagnetic attraction along present edges and antiferromagnetic repulsion across absent edges. Choosing the configuration-model coupling $a_(i j) = 1 - gamma p_(i j)$ and $b_(i j) = gamma p_(i j)$ with null probability $p_(i j) = k_i k_j / (2 m)$ and tunable resolution parameter $gamma > 0$ yields the canonical *Reichardt-Bornholdt (RB) Potts Hamiltonian*
$ cal(H)_(gamma) (bold(sigma)) = -1 / 2 sum_(i != j) (A_(i j) - gamma (k_i k_j) / (2 m)) delta(sigma_i, sigma_j). $

#definition(title: "Potts Spin-Glass Gibbs Measure and Local Spin Field")[
  Given inverse temperature $beta = 1 / T > 0$ and resolution $gamma > 0$, the equilibrium *Gibbs-Boltzmann distribution* over $bold(sigma) in {1, dots, q}^n$ is
  $ pi_(beta) (bold(sigma)) = 1 / (Z(beta)) exp(-beta cal(H)_(gamma) (bold(sigma))), \
    Z(beta) = sum_(bold(sigma) in {1, dots, q}^n) exp(-beta cal(H)_(gamma) (bold(sigma))). $
  For vertex $i$, let $k_(i, s) = sum_(j != i, sigma_j = s) A_(i j)$ denote the edge weight linking $i$ to spin state $s$, and let $K_s = sum_(j : sigma_j = s) k_j$ be the total degree of spin $s$. The *local effective energy* of placing spin $i$ in state $s$ (holding $bold(sigma)_(-i)$ fixed) is
  $ h_i (s | bold(sigma)_(-i)) = -k_(i, s) + gamma (k_i) / (2 m) (K_s - k_i delta(sigma_i, s)). $
]

Because minimizing $cal(H)_(gamma) (bold(sigma))$ is NP-hard @brandes2008modularity, Spinglass community detection locates low-energy states via *simulated annealing*: single-spin Gibbs heat-bath updates sample $sigma_i$ from its conditional Boltzmann law while slowly lowering temperature $T_t -> 0$ ($beta_t -> infinity$).

#algbox(
  title: "Reichardt-Bornholdt Potts Spinglass Simulated Annealing",
  inputs: [Connected graph $G = (V, E, w)$; max spin states $q$; resolution $gamma > 0$; cooling schedule $T_0 > T_("stop") > 0$, cooling factor $alpha in (0, 1)$, sweeps per temperature $L$.],
  outputs: [Low-energy spin partition $hat(bold(sigma)) in {1, dots, q}^n$ and occupied community clusters.],
)[
  1. *Random spin initialization:* Sample $sigma_i ~ "Unif"{1, dots, q}$ independently for $i = 1, dots, n$; precompute $K_s <- sum_(j : sigma_j = s) k_j$ for $s = 1, dots, q$; set $T <- T_0$.
  2. *Annealing loop:* While $T > T_("stop")$ and spin flips remain active:
     - Set $beta <- 1 / T$ and repeat for $L$ full sweeps:
       - Visit vertices $i in V$ in random order. Let $r <- sigma_i$.
       - Accumulate neighbor weights $k_(i, s)$ for all spin states $s$ adjacent to $i$ in $cal(O)(k_i)$ time.
       - Compute local energies $h_i (s) = -k_(i, s) + gamma k_i (K_s - k_i delta(r, s)) / (2 m)$ for $s = 1, dots, q$.
       - Sample new spin $sigma'_i$ from the heat-bath conditional distribution
         $ "Pr"(sigma'_i = s | bold(sigma)_(-i)) = (exp(-beta h_i (s))) / (sum_(u=1)^q exp(-beta h_i (u))). $
       - If $sigma'_i != r$, update $K_r <- K_r - k_i$, $K_(sigma'_i) <- K_(sigma'_i) + k_i$, and $sigma_i <- sigma'_i$.
     - Cool temperature: $T <- alpha T$.
  3. *Zero-temperature quench & return:* Perform greedy deterministic moves $sigma_i <- arg min_s h_i (s)$ until no single-spin move lowers $cal(H)_(gamma)$; return non-empty spin classes.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Hamiltonian-Modularity Isomorphism and Ground-State Cohesion")

We establish the exact algebraic equivalence between the RB Potts Hamiltonian and multi-resolution modularity, and prove necessary structural inequalities obeyed by every ground state.

#theorem(title: "Exact Hamiltonian-Modularity Isomorphism and Ground-State Cohesion", origin: "Original Proof")[
  Define the resolution-$gamma$ modularity functional
  $ Q_(gamma) (bold(sigma)) = 1 / (2 m) sum_(i=1)^n sum_(j=1)^n (A_(i j) - gamma (k_i k_j) / (2 m)) delta(sigma_i, sigma_j). $
  Then for every simple graph ($A_(i i) = 0$) and every spin configuration $bold(sigma) in {1, dots, q}^n$:
  (i) $cal(H)_(gamma) (bold(sigma)) = -m Q_(gamma) (bold(sigma)) - (gamma / (4 m)) sum_(i=1)^n k_i^2$, so minimizing $cal(H)_(gamma)$ is strictly equivalent to maximizing $Q_(gamma)$ (reducing to Newman-Girvan modularity @newman2004finding at $gamma = 1$);
  (ii) If $bold(sigma)^*$ is a global ground state with $K < q$ occupied spin states (so at least one unused spin exists), then every occupied community $C_r$ satisfies:
  - *Internal cohesion:* For every proper bipartition $C_r = S union T$,
    $ e(S, T) >= gamma ("vol"(S) "vol"(T)) / (2 m), $
    where $e(S, T) = sum_(u in S, v in T) A_(u v)$; in particular, every $G[C_r]$ is connected when $gamma > 0$ and all $k_i > 0$;
  - *External separation:* For every pair of distinct communities $C_r != C_s$,
    $ e(C_r, C_s) <= gamma ("vol"(C_r) "vol"(C_s)) / (2 m). $
]

#proof[
  For (i), split the full double sum in $-m Q_(gamma) (bold(sigma))$ into diagonal ($i = j$) and off-diagonal ($i != j$) terms. Since $A_(i i) = 0$ and $delta(sigma_i, sigma_i) = 1$,
  $ -m Q_(gamma) (bold(sigma)) &= -1 / 2 sum_(i != j) (A_(i j) - gamma (k_i k_j) / (2 m)) delta(sigma_i, sigma_j) \
                              &quad + 1 / 2 sum_(i=1)^n gamma (k_i^2) / (2 m) \
                              &= cal(H)_(gamma) (bold(sigma)) + gamma / (4 m) sum_(i=1)^n k_i^2. $
  Because the last term is independent of $bold(sigma)$, $arg min cal(H)_(gamma) = arg max Q_(gamma)$.

  For (ii), suppose $bold(sigma)^*$ minimizes $cal(H)_(gamma)$ and has an unused spin $u in {1, dots, q}$. Take any community $C_r$ and any proper bipartition $C_r = S union T$ ($S, T != emptyset$). Reassigning all vertices in $T$ from spin $r$ to the unused spin $u$ produces a valid configuration $bold(sigma)'$ whose energy change equals the negation of the removed cross terms:
  $ cal(H)_(gamma) (bold(sigma)') - cal(H)_(gamma) (bold(sigma)^*) = sum_(i in S, j in T) (A_(i j) - gamma (k_i k_j) / (2 m)) \
    = e(S, T) - gamma ("vol"(S) "vol"(T)) / (2 m). $
  Optimality of $bold(sigma)^*$ forces $cal(H)_(gamma) (bold(sigma)') - cal(H)_(gamma) (bold(sigma)^*) >= 0$, yielding the internal cohesion inequality. If $G[C_r]$ were disconnected into $S$ and $T$ with positive degrees, then $e(S, T) = 0 < gamma "vol"(S) "vol"(T) / (2 m)$, a contradiction; hence $G[C_r]$ is connected. Similarly, merging distinct communities $C_r$ and $C_s$ changes energy by
  $ cal(H)_(gamma) (bold(sigma)'') - cal(H)_(gamma) (bold(sigma)^*) &= -(e(C_r, C_s) - gamma ("vol"(C_r) "vol"(C_s)) / (2 m)) \
                                                                    &>= 0, $
  yielding the external separation bound (cf. @reichardt2006spinglass).
]

#subsec("B", "Exact Multi-Scale Phase Window Overcoming the Resolution Limit")

Standard modularity ($gamma = 1$) fails on rings of $R$ cliques once $R > 2(l + 1)$ @fortunato2007resolution. We derive the exact lower and upper critical resolution thresholds $[gamma_("merge"), gamma_("split")]$ in which the RB Potts ground state resolves all $R$ cliques faithfully without either merging pairs or shattering cliques.

#theorem(title: "Exact Resolution Phase Window on a Ring of Cliques", origin: "Original Derivation")[
  Consider a ring of $R >= 3$ identical cliques $C_1, dots, C_R tilde.op K_p$ ($p >= 3$, internal edges $l = binom(p, 2)$) linked consecutively by single bridge edges, so $m = R(l + 1)$. Assume $q > R$. Then the exact $R$-clique spin state $bold(sigma)^*$ satisfies both ground-state cohesion and separation (Theorem II.A) if and only if $gamma$ lies inside the closed phase window
  $ gamma_("merge") <= gamma <= gamma_("split"), $
  where the critical phase boundaries are exactly
  $ gamma_("merge") &= R / (2 (l + 1)) = R / (p^2 - p + 2), \
    gamma_("split") &= (R (p^2 - p + 2)(p - 1)) / (p ((p - 1)^2 + 1)). $
  Furthermore, the phase-window ratio scales as $gamma_("split") / gamma_("merge") = Theta(p^2) > 1$ for every $p >= 3$, proving that tuning $gamma$ completely eliminates the Fortunato-Barthélemy resolution limit for arbitrarily large $R$.
]

#proof[
  In the ring of cliques, each clique $C_r$ has $p - 2$ internal vertices of degree $p - 1$ and $2$ bridge vertices of degree $p$, so $"vol"(C_r) = p(p - 1) + 2 = 2(l + 1)$ and $2 m = 2 R(l + 1)$.

  *Lower boundary (clique merging):* Adjacent cliques $C_r, C_(r+1)$ share $e(C_r, C_(r+1)) = 1$ edge (non-adjacent cliques share $0$). By external separation, no pair merges iff
  $ 1 <= gamma ("vol"(C_r) "vol"(C_(r+1))) / (2 m) = gamma (4 (l + 1)^2) / (2 R (l + 1)) = (2 gamma (l + 1)) / R, $
  which rearranges identically to $gamma >= R / (2(l + 1)) = R / (p^2 - p + 2) = gamma_("merge")$. At $gamma = 1$, this recovers the classical limit $R <= 2(l + 1)$ @fortunato2007resolution.

  *Upper boundary (clique shattering):* Partition a single clique $C_r$ into non-empty subsets $S union T$ with $|S| = s$ and $|T| = p - s$ ($1 <= s <= p - 1$). Because $C_r$ is complete, $e(S, T) = s(p - s)$. Internal cohesion requires $gamma <= 2 m e(S, T) / ("vol"(S) "vol"(T))$ for all splits, so $gamma_("split")$ is governed by maximizing $"vol"(S) "vol"(T) / e(S, T)$. There are two placements of the two bridge vertices (each adding $+1$ to its side's degree sum):
  - *Case 1 (one bridge in $S$, one in $T$):* $"vol"(S) = s(p - 1) + 1$ and $"vol"(T) = (p - s)(p - 1) + 1$, giving
    $ ("vol"(S) "vol"(T)) / (e(S, T)) = (p - 1)^2 + (p(p - 1) + 1) / (s(p - s)). $
    Since $s(p - s)$ is minimized on ${1, dots, p - 1}$ at $s = 1$ (value $p - 1$), this expression attains its exact maximum at $s = 1$ (detaching a single bridge vertex):
    $ M_1 = (p - 1)^2 + p + 1 / (p - 1) = (p((p - 1)^2 + 1)) / (p - 1). $
  - *Case 2 (both bridges in $S$, $s >= 2$):* $"vol"(S) = s(p - 1) + 2$, $"vol"(T) = (p - s)(p - 1)$, giving
    $ ("vol"(S) "vol"(T)) / (e(S, T)) = (p - 1)^2 + (2(p - 1)) / s <= (p - 1)^2 + p - 1 < M_1. $
  Thus Case 1 at $s = 1$ is the unique global maximizer over all bipartitions of $C_r$. Substituting $M_1$ into $gamma <= 2 m / M_1$ with $2 m = R(p^2 - p + 2)$ gives the exact formula for $gamma_("split")$. Finally, dividing gives
  $ (gamma_("split")) / (gamma_("merge")) = ((p^2 - p + 2)^2 (p - 1)) / (p ((p - 1)^2 + 1)) >= 128 / 15 approx 8.53 $
  already at $p = 3$ (and growing as $p^2$), confirming a wide, non-empty stable phase for every $R >= 3$.
]

#subsec("C", "Detailed Balance and Simulated Annealing Convergence")

#theorem(title: "Detailed Balance of Heat-Bath Dynamics and Annealing Convergence", origin: "Classical — Cited")[
  At any fixed inverse temperature $beta > 0$, the single-spin Gibbs heat-bath kernel selecting vertex $i$ uniformly at random and resampling $sigma_i$ from $"Pr"(sigma'_i = s | bold(sigma)_(-i)) prop exp(-beta h_i (s))$ is irreducible, aperiodic, and satisfies strict detailed balance with respect to the Gibbs-Boltzmann measure $pi_(beta) (bold(sigma))$ @reichardt2006spinglass. Moreover, under logarithmic cooling $T_t >= Delta^* / ln(t + 2)$, where $Delta^*$ is the maximum energy-barrier depth of any non-global local minimum of $cal(H)_(gamma)$, the inhomogeneous Markov chain converges in probability to the set of global ground states: $lim_(t -> infinity) "Pr"(bold(sigma)^((t)) in arg min cal(H)_(gamma)) = 1$.
]

#proof[
  Let $bold(sigma)$ and $bold(sigma)'$ differ only at coordinate $i$ ($sigma_i = r != s = sigma'_i$). The transition probability $P_(beta) (bold(sigma) -> bold(sigma)')$ equals
  $ P_(beta) (bold(sigma) -> bold(sigma)') = 1 / n dot (exp(-beta h_i (s | bold(sigma)_(-i)))) / (sum_(u=1)^q exp(-beta h_i (u | bold(sigma)_(-i)))). $
  Crucially, $cal(H)_(gamma) (bold(sigma)) = cal(H)_(-i) (bold(sigma)_(-i)) + h_i (r | bold(sigma)_(-i))$ and $cal(H)_(gamma) (bold(sigma)') = cal(H)_(-i) (bold(sigma)_(-i)) + h_i (s | bold(sigma)_(-i))$, where $cal(H)_(-i)$ omits all terms incident to $i$. Substituting into $pi_(beta) (bold(sigma)) P_(beta) (bold(sigma) -> bold(sigma)')$ gives
  $ pi_(beta) (bold(sigma)) P_(beta) (bold(sigma) -> bold(sigma)') \
    quad = (exp(-beta [cal(H)_(-i) + h_i (r) + h_i (s)])) / (n Z(beta) sum_(u=1)^q exp(-beta h_i (u))), $
  which is manifestly symmetric in $r$ and $s$, proving $pi_(beta) (bold(sigma)) P_(beta) (bold(sigma) -> bold(sigma)') = pi_(beta) (bold(sigma)') P_(beta) (bold(sigma)' -> bold(sigma))$. Because $exp(-beta h_i (s)) > 0$ for all finite $beta$ and every state $s in {1, dots, q}$, every single-spin transition (including self-transitions $s = r$) has strictly positive probability, guaranteeing irreducibility and aperiodicity; hence $pi_(beta)$ is the unique stationary distribution. As $beta -> infinity$, $pi_(beta)$ concentrates uniformly on $arg min cal(H)_(gamma)$, and the Hajek-Geman theorem establishes global convergence under $T_t >= Delta^* / ln(t + 2)$ (whereas practical geometric schedules $T_t = T_0 alpha^t$ run in polynomial time at the risk of metastable spin-glass freezing).
]

#subsec("D", "High-Temperature Susceptibility and Spectral Link")

#proposition(title: "Exact High-Temperature Susceptibility Identity and Zero-Temperature Termination", origin: "Original Proof")[
  (i) In the high-temperature paramagnetic limit $beta -> 0$, the exact Gibbs co-membership correlation between any distinct vertices $i != j$ satisfies
  $ bb(E)_(beta) [delta(sigma_i, sigma_j)] &= 1 / q + beta ((q - 1) / q^2) (A_(i j) - gamma (k_i k_j) / (2 m)) \
                                           &quad + cal(O)(beta^2), $
  so the first-order spin susceptibility matrix is directly proportional to the generalized modularity matrix $B_(gamma) = A - gamma (2 m)^(-1) k k^T$;
  (ii) In the zero-temperature quench ($beta -> infinity$) on an unweighted graph with rational $gamma = a / b$ ($a, b in bb(N)_+$), every strictly improving single-spin move lowers $cal(H)_(gamma)$ by at least $1 / (2 b m)$, guaranteeing finite termination in at most $cal(O)(b m^2)$ flips.
]

#proof[
  For (i), write $cal(H)_(gamma) (bold(sigma)) = -sum_(u < v) J_(u v) delta(sigma_u, sigma_v)$ with $J_(u v) = A_(u v) - gamma k_u k_v / (2 m)$. At $beta = 0$, spins are i.i.d. uniform on ${1, dots, q}$, so $bb(E)_0 [delta(sigma_i, sigma_j)] = 1 / q$. Differentiating $bb(E)_(beta) [delta(sigma_i, sigma_j)]$ at $beta = 0$ yields
  $ (dif / (dif beta) bb(E)_(beta) [delta(sigma_i, sigma_j)])_(beta = 0) &= -"Cov"_0 (delta(sigma_i, sigma_j), cal(H)_(gamma)) \
    &= sum_(u < v) J_(u v) "Cov"_0 (delta(sigma_i, sigma_j), delta(sigma_u, sigma_v)). $
  Under i.i.d. uniform spins, if ${u, v} = {i, j}$ then $"Var"_0 (delta(sigma_i, sigma_j)) = (1 / q)(1 - 1 / q) = (q - 1) / q^2$. If ${u, v} != {i, j}$—even when the two pairs share a single vertex, say $u = i, v != j$—the indicators $delta(sigma_i, sigma_j)$ and $delta(sigma_i, sigma_v)$ are pairwise independent because $"Pr"(sigma_i = sigma_j = sigma_v) = q / q^3 = 1 / q^2 = (1 / q)^2$, so $"Cov"_0 = 0$. Every cross-term vanishes identically, leaving $((q - 1) / q^2) J_(i j)$ as the exact linear coefficient.

  For (ii), since $A_(i j) in {0, 1}$ and $k_i, k_j in bb(Z)$, the value $2 b m cal(H)_(gamma) (bold(sigma))$ is an integer bounded in $[-2 b m^2, 2 a m^2]$. Every strictly energy-lowering spin flip reduces $2 b m cal(H)_(gamma)$ by at least $1$, bounding the total number of improving flips by $2 (a + b) m^2 = cal(O)(b m^2)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Maintaining the community degree totals $K_s$ ($s = 1, dots, q$) reduces each vertex's local field computation to scanning its $k_i$ incident edges (accumulating $k_(i, s)$) plus $q$ scalar evaluations of $h_i (s)$. One complete heat-bath sweep over all $n$ vertices therefore costs $cal(O)(m + n q)$ time (and $cal(O)(m + n)$ for Metropolis single-proposal updates). Across $N_("sw") = L ceil(log(T_("stop") / T_0) / log alpha)$ total sweeps, total time is $cal(O)(N_("sw") (m + n q))$ and memory is strictly $cal(O)(m + n + q)$. In standard implementations on moderate graphs with $q = 25$--$50$ and small $n$ ($n <= 10^3$), Spinglass executes in sub-millisecond to millisecond time; on large graphs ($n > 10^4$), long relaxation times near the spin-glass freezing temperature $T_g$ inflate $N_("sw")$ sharply @reichardt2006spinglass.

*Energy Landscape and Geometry.* Unlike CNM (deterministic agglomerative) and Leading Eigenvector (deterministic spectral bisection), Spinglass explores a rugged non-convex energy landscape stochastically. Finite temperature allows spins to jump out of shallow local traps—yielding high accuracy on compact, mildly overlapping clusters—while the resolution parameter $gamma$ smoothly sweeps cluster granularity from global ($gamma < 1$) to micro-scale ($gamma > 1$). However, the configuration null $k_i k_j / (2 m)$ assumes a single connected component: on disconnected or sparsely bridged elongated $k$NN manifolds, antiferromagnetic repulsion across distant parts of a chain creates glassy domain walls that fragment or merge shapes erratically.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(N_("sw") (m + n q))$ for $N_("sw")$ heat-bath sweeps and $q$ spins],
    [Space Complexity], [$cal(O)(m + n + q)$ via adjacency lists and degree cache $K_s$],
    [Cluster Count $K$], [Emergent ($K <= q$), tuned continuously via resolution $gamma$],
    [Determinism], [Stochastic (simulated annealing); consensus across runs advisable],
    [Connectedness], [Guaranteed in true ground state ($K < q, gamma > 0$, Thm. II.A)],
    [Invariances], [Positive weight-scale invariant; permutation equivariant in law],
    [Key Limitation], [Requires connected graph; slow annealing on large $n$; glassy traps],
  ),
  caption: [Operational and complexity profile of Spinglass Community Detection.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multi-Scale Brain Connectomics (fMRI / Structural Networks):* Cortical parcellation spans hierarchical scales from hemispheres down to functional sub-nuclei. Sweeping $gamma$ across the stability plateaus of $cal(H)_(gamma)$ (Theorem II.B) reveals multi-resolution functional modules while naturally accommodating signed positive/negative functional correlations via separate ferromagnetic and antiferromagnetic couplings.
2. *Compact Overlapping Biological & Phenotypic Similarity Graphs:* Gene co-expression modules and phenotypic $k$NN graphs (e.g., botanical or clinical cohorts of $10^2$--$10^3$ samples) often exhibit noisy boundary overlap where greedy agglomeration makes irreversible mistakes. Thermal fluctuations during annealing escape shallow local minima and achieve superior partition recovery.
3. *Signed Social Balance and Polarization Networks:* Because the general RB Hamiltonian accepts both positive ($A_(i j) > 0$, friend) and negative ($A_(i j) < 0$, foe) edges directly via $a_(i j)$ and $b_(i j)$, Spinglass naturally partitions signed trust/distrust networks into mutually cohesive, mutually antagonistic factions without ad-hoc edge thresholding.
4. *Financial Portfolio Correlation Clustering:* Asset correlation matrices contain strong global market modes and overlapping sector hierarchies. Tuning $gamma$ isolates fine-grained industry sub-sectors while stochastic restarts quantify cluster assignment uncertainty via co-occurrence matrices.
5. *Benchmarking Ground-State Modularity Bounds:* On small-to-medium scientific networks ($n <= 5000$), slow-cooled Spinglass frequently attains strictly higher modularity $Q$ than deterministic heuristics (CNM, Leading Eigenvector, Walktrap @pons2006walktrap), serving as a high-precision reference baseline.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Spinglass exhibits four primary failure modes: (i) *disconnected graphs* where the configuration null $2 m$ couples independent components artificially (standard implementations refuse disconnected inputs outright); (ii) *elongated, sparse non-convex manifolds* where low-conductance geometric chains trap thermal domain walls and suffer severe fragmentation or random cross-merging; (iii) *large-scale networks ($n > 10^5$)* where annealing convergence slows super-linearly; and (iv) *underspecified spin cap $q < K^*$* where the pigeonhole principle forces distinct true communities to share spin labels.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Compact / Noisy $k$NN Data], [Very high: escapes greedy traps], [Run multiple seeds; verify $q > K$],
    [Multi-Scale / Signed Graphs], [Very high: exact $gamma$ phase control], [Tune $gamma$ via stability plateaus],
    [Multi-Shape Elongated Clouds], [Very low: glassy domain pinning], [Use CNM (Alg. 115) or HDBSCAN @campello2013hdbscan],
    [Disconnected Networks], [Inapplicable globally without split], [Apply per connected component],
    [Million-Node Web Graphs], [Low: prohibitive $N_("sw")$ sweeps], [Use Leiden @traag2019leiden or LPA @raghavan2007lpa],
  ),
  caption: [Domain suitability matrix for Spinglass Community Detection.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Spinglass generalizes virtually all modularity-family algorithms within a single Hamiltonian umbrella, connecting CNM (Alg. 115), Leading Eigenvector (Alg. 116), Louvain @blondel2008louvain, and Leiden @traag2019leiden through its choice of null matrix $p_(i j)$ and optimization dynamics.

#proposition(title: "Unification of Modularity, Constant Potts Model, and Leiden Objectives", origin: "Original Proof")[
  Under the generalized RB Potts Hamiltonian $cal(H)(bold(sigma)) = -1 / 2 sum_(i != j) (A_(i j) - gamma p_(i j)) delta(sigma_i, sigma_j)$:
  (i) Setting $p_(i j) = k_i k_j / (2 m)$ and $gamma = 1$ yields an exact affine equivalent of Newman-Girvan modularity maximized by CNM (Alg. 115), Leading Eigenvector (Alg. 116), and Louvain @blondel2008louvain;
  (ii) Setting $p_(i j) = 1$ (uniform Erdős-Rényi null) yields the Constant Potts Model (CPM) objective $cal(H)_("CPM") (bold(sigma)) = -sum_(r=1)^q (e_r - gamma binom(|C_r|, 2))$ optimized by Leiden @traag2019leiden, which is strictly invariant to adding disjoint components and therefore completely free of the global resolution limit.
]

#proof[
  Part (i) follows directly from Theorem II.A(i). For (ii), substituting $p_(i j) = 1$ into $cal(H)(bold(sigma))$ gives
  $ cal(H)_("CPM") (bold(sigma)) &= -sum_(r=1)^q sum_(i < j in C_r) (A_(i j) - gamma) \
                                 &= -sum_(r=1)^q (e_r - gamma binom(|C_r|, 2)), $
  where $e_r$ is the internal edge count of $C_r$. Because the penalty $gamma binom(|C_r|, 2)$ depends exclusively on $|C_r|$ and $gamma$—with zero dependence on the global edge count $m$—embedding any subgraph $H$ inside a larger network $G$ leaves the relative energy of all partitions of $H$ strictly unchanged, eliminating the global $sqrt(m / 2)$ resolution limit @fortunato2007resolution @traag2019leiden.
]

*Empirical Behavior on ClusterOverview Benchmarks.* On the ClusterOverview 4-dataset benchmark, Spinglass achieves *Iris* $"Sil" = 0.411$, $"ARI" = 0.636$; *Customers* $"Sil" = 0.356$; *Multi-Shape* $"Sil" = 0.088$; and *High-D* $"Sil" = 0.303$, with the fastest mean runtime of all four chapters at $0.546$ ms. This striking contrast—best-in-group *Iris* ARI ($0.636$ vs. $0.373$ for CNM, $0.495$ for Leading Eigenvector, and $0.564$ for HCS) paired with dramatic collapse on *Multi-Shape* ($0.088$ vs. $0.782$ for CNM)—directly illustrates our theoretical findings. On the compact, connected 150-point *Iris* $k$NN graph, finite-temperature heat-bath sampling escapes the greedy bridge traps of CNM and the binary-bisection constraints of Leading Eigenvector while exploiting cached $cal(O)(m + n q)$ sweeps for sub-millisecond speed. Conversely, on *Multi-Shape*, elongated low-conductance manifolds with heterogeneous shape densities trigger glassy spin-domain freezing under fixed $gamma = 1$ and short geometric cooling, causing arbitrary splits along thin geometric necks and severely degrading Silhouette. Spinglass is thus the premier modularity optimizer for compact, connected, noisy graphs where thermal escape matters, but should be replaced by CNM, Leiden, or density-based methods on sparse elongated manifolds.
