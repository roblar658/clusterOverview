#import "../template/ieee.typ": *

#alg_entry(
  id: 6,
  name: "CLARANS (Randomized Graph Search)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1994 / 2002",
  paradigm: "Stochastic Hill Climbing on the Johnson Graph",
  complexity: $O(r dot.c m_"max" dot.c n dot.c T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(cal(X), d)$ be a finite metric space of cardinality $|cal(X)| = n$ with symmetric dissimilarity $d : cal(X) times cal(X) -> [0, infinity)$, and let $k in {1, dots, n - 1}$ be the desired number of medoids. As in PAM (Alg. 4) and CLARA (Alg. 5), the objective is to select $M subset.eq cal(X)$ with $|M| = k$ minimizing the total assignment dissimilarity:
$ L(M) = sum_(i=1)^n min_(m in M) d(x_i, m). $

Ng and Han @ng2002clarans reformulate $k$-medoids optimization explicitly as a discrete search problem on a structured combinatorial graph. Define the *Johnson graph* $cal(G)_(n, k) = (cal(V), cal(E))$ whose vertex set $cal(V) = binom(cal(X), k)$ consists of all $binom(n, k)$ distinct $k$-element subsets of $cal(X)$, with an undirected edge ${M, M'} in cal(E)$ if and only if $|M inter M'| = k - 1$. Equivalently, every neighbor $M' in N(M)$ is obtained from $M$ via a single medoid swap $M' = (M backslash {m_j}) union {x_c}$ for some $m_j in M$ and $x_c in cal(X) backslash M$. The graph $cal(G)_(n, k)$ is $k(n - k)$-regular, vertex-transitive, and has diameter $min(k, n - k)$.

Under this lens, PAM @kaufman1990finding performs exhaustive steepest descent inspecting all $k(n - k)$ neighbors at every visited vertex, whereas CLARA @srikantha2011clara restricts the search permanently to a fixed induced subgraph $cal(G)_(s, k) subset.eq cal(G)_(n, k)$. *CLARANS* (Clustering Large Applications Based on Randomized Search) @ng2002clarans retains access to the *entire* vertex set $cal(V)$ while sampling edges dynamically: from each current vertex $M$, it draws neighbors uniformly at random and transitions immediately to the first strictly improving neighbor, halting a search trajectory after $m_"max"$ (`maxneighbor`) consecutive rejected trials and repeating for $r$ (`numlocal`) independent restarts.

#definition(title: "Improving Neighborhood & Stopping State")[
  At any vertex $M in cal(V)$ of $cal(G)_(n, k)$, define the *strictly improving neighborhood* and its *improving density* by
  $ N^- (M) &= {M' in N(M) : L(M') < L(M)}, \
    rho(M) &= (|N^- (M)|) / (k(n - k)) in [0, 1]. $
  A vertex $M$ is a *true 1-swap local minimum* if $rho(M) = 0$, and an *$epsilon$-approximate local minimum* if $rho(M) <= epsilon$. CLARANS terminates a restart whenever $m_"max"$ i.i.d. uniform draws from $N(M)$ fall in $N(M) backslash N^- (M)$.
]

Crucially, CLARANS avoids storing an $n times n$ matrix by caching only three scalars per point—nearest medoid index $c_1 (i)$, nearest distance $d_1 (i)$, and second-nearest distance $d_2 (i)$—evaluating each random swap in $O(n)$ time (Algorithm 1).

#algbox(
  title: "CLARANS (Randomized Johnson-Graph Search)",
  inputs: [Dataset $(cal(X), d)$ of size $n$, integer $k >= 1$, restart count $r >= 1$ (`numlocal`), trial cap $m_"max" >= 1$ (`maxneighbor`).],
  outputs: [Best medoid set $M_"best" subset.eq cal(X)$ and partition $cal(C)(M_"best")$.],
)[
  1. Set $L_"best" <- +infinity$ and $M_"best" <- emptyset$.
  2. *For* $ell = 1, dots, r$ *do*:
  3. #h(1em) Sample initial vertex $M tilde "Unif"(binom(cal(X), k))$ and cache $(c_1 (i), d_1 (i), d_2 (i))$ for all $i = 1, dots, n$; set $t <- 0$.
  4. #h(1em) *While* $t < m_"max"$ *do*:
  5. #h(2em) Sample $m_j tilde "Unif"(M)$ and $x_c tilde "Unif"(cal(X) backslash M)$.
  6. #h(2em) Compute exact incremental swap delta $Delta L(m_j -> x_c)$ in $O(n)$ time (Lemma 1).
  7. #h(2em) *If* $Delta L(m_j -> x_c) < 0$:
  8. #h(3em) Update $M <- (M backslash {m_j}) union {x_c}$, refresh $(c_1, d_1, d_2)$, and reset $t <- 0$.
  9. #h(2em) *Else*: increment $t <- t + 1$.
  10. #h(1em) *If* $L(M) < L_"best"$, set $L_"best" <- L(M)$ and $M_"best" <- M$.
  11. Assign each $x_i in cal(X)$ to $c_i = arg min_(m in M_"best") d(x_i, m)$.
  12. *Return* $(M_"best", cal(C)(M_"best"))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical properties of CLARANS: (1) an exact $O(n)$ incremental swap-cost identity and an original stopping-time certificate bounding the residual improving density $rho(M_tau)$; (2) a self-contained proof of the classical factor-5 approximation of true 1-swap local optima on metric spaces; (3) the exact geometric-distribution phase transition governing `maxneighbor` scaling; (4) an explicit Euclidean simplex counterexample proving unbounded stalling error when $m_"max"$ is held fixed; and (5) the $O(k)$ amortized acceleration of FastCLARANS.

#subsec("A", "Exact O(n) Delta Cost and Stopping-Time Optimality Bound")

#lemma(title: "Exact O(n) Incremental Swap-Cost Identity", origin: "Original Derivation")[
  For $M = {m_1, dots, m_k}$ with cached $(c_1 (i), d_1 (i), d_2 (i))$, swapping $m_j in M$ for $x_c in cal(X) backslash M$ changes $L(M)$ by exactly
  $ Delta L(m_j -> x_c) &= sum_(i : c_1 (i) = j) (min(d_2 (i), d(x_i, x_c)) - d_1 (i)) \
    &quad + sum_(i : c_1 (i) != j) min(0, space d(x_i, x_c) - d_1 (i)). $
  Hence evaluating a random neighbor in $cal(G)_(n, k)$ requires exactly $n$ distance queries and $O(n)$ arithmetic operations, independent of $k$.
]

#proof[
  Partition ${1, dots, n}$ into $C_j = {i : c_1 (i) = j}$ and its complement. For any $i in C_j$, deleting $m_j$ leaves $d_2 (i) = min_(ell != j) d(x_i, m_ell)$ as the distance to the nearest surviving medoid; inserting $x_c$ therefore yields new distance $min(d_2 (i), d(x_i, x_c))$, contributing $min(d_2 (i), d(x_i, x_c)) - d_1 (i)$. For any $i in.not C_j$, its current nearest medoid $m_(c_1 (i))$ remains in $M'$, so the new distance is $min(d_1 (i), d(x_i, x_c)) = d_1 (i) + min(0, d(x_i, x_c) - d_1 (i))$. Summing over all $i = 1, dots, n$ proves the identity identically.
]

#theorem(title: "Acyclic Termination and Stopping-Time Certificate", origin: "Original Proof")[
  Let $M_0, M_1, dots, M_tau$ be the sequence of vertices visited during a single restart of Algorithm 1, where $tau$ is the random number of accepted improving swaps before halting. Then:
  1. $L(M_0) > L(M_1) > dots > L(M_tau)$, so the trajectory is strictly acyclic and terminates almost surely with $tau <= binom(n, k) - 1$.
  2. For any threshold $epsilon in (0, 1)$, the stopping vertex $M_tau$ satisfies the exact unconditional stopping-time tail bound
     $ bb(P)(rho(M_tau) > epsilon) <= (bb(E)[tau] + 1) e^(-epsilon m_"max"). $
     In particular, setting $m_"max" >= epsilon^(-1) ln((bb(E)[tau] + 1) / delta)$ certifies $rho(M_tau) <= epsilon$ with probability at least $1 - delta$.
]

#proof[
  (1) Step 7 transitions only when $Delta L(m_j -> x_c) < 0$, enforcing strict monotone decrease of $L(M_t)$. Hence no vertex in $cal(V)$ can be revisited; since $|cal(V)| = binom(n, k) < infinity$, every trajectory halts after at most $binom(n, k) - 1$ transitions.

  (2) Let $cal(F)_t$ be the $sigma$-algebra generated by the initial vertex $M_0$ and all random draws up to the instant the walk enters vertex $M_t$. On the $cal(F)_t$-measurable event ${tau >= t}$, the vertex $M_t$ and its improving density $rho(M_t)$ are $cal(F)_t$-measurable, whereas the next $m_"max"$ neighbor trials at $M_t$ are drawn freshly and independently of $cal(F)_t$. Therefore, halting precisely at step $t$ has exact conditional probability
  $ bb(P)(tau = t | cal(F)_t) = (1 - rho(M_t))^(m_"max") quad "on" {tau >= t}. $
  Multiplying by the $cal(F)_t$-measurable indicator $bold(1)_{rho(M_t) > epsilon}$ and taking expectations gives, for every $t >= 0$,
  $ bb(P)(tau = t, space rho(M_t) > epsilon) &= bb(E)[bold(1)_{tau >= t, space rho(M_t) > epsilon} \
    &quad dot.c (1 - rho(M_t))^(m_"max")] \
    &<= (1 - epsilon)^(m_"max") bb(P)(tau >= t) \
    &<= e^(-epsilon m_"max") bb(P)(tau >= t). $
  Because $tau < infinity$ almost surely, summing over all $t = 0, 1, 2, dots$ and applying the integer tail-sum identity $sum_(t=0)^infinity bb(P)(tau >= t) = 1 + sum_(t=1)^infinity bb(P)(tau >= t) = bb(E)[tau] + 1$ yields
  $ bb(P)(rho(M_tau) > epsilon) &= sum_(t=0)^infinity bb(P)(tau = t, space rho(M_t) > epsilon) \
    &<= (bb(E)[tau] + 1) e^(-epsilon m_"max"), $
  completing the proof.
]

#subsec("B", "Factor-5 Approximation of Metric Local Optima")

Whenever $m_"max"$ is large enough that CLARANS reaches a true local minimum ($rho(M_tau) = 0$), metric local-search theory guarantees constant-factor global quality (cf. @kanungo2004local @ng2002clarans).

#theorem(title: "Factor-5 Bound for 1-Swap Metric Local Optima", origin: "Classical — Cited")[
  Let $(cal(X), d)$ be a metric space and let $M subset.eq cal(X)$ ($|M| = k$) be any true 1-swap local minimum of $cal(G)_(n, k)$ ($rho(M) = 0$). Then for any globally optimal $k$-medoid set $M^*$,
  $ L(M) <= 5 L(M^*). $
]

#proof[
  For each $x in cal(X)$, let $m(x) in M$ and $m^*(x) in M^*$ be its closest medoids, with $s_x = d(x, m(x))$ and $o_x = d(x, m^*(x))$, so $L(M) = sum_x s_x$ and $L(M^*) = sum_x o_x$. Define the nearest-medoid map $phi : M^* -> M$ by $phi(m^*) in arg min_(m in M) d(m^*, m)$, and partition $M$ by preimage degree: $M_0 = {m : |phi^(-1)(m)| = 0}$, $M_1 = {m : |phi^(-1)(m)| = 1}$, and $M_(>= 2) = {m : |phi^(-1)(m)| >= 2}$.

  Construct $k$ test swaps $(m -> m^*)$ pairing every $m^* in M^*$ with a medoid $m in M$ as follows: if $m in M_1$ with $phi^(-1)(m) = {m^*}$, pair $(m, m^*)$; let $R^* = union.big_(m in M_(>= 2)) phi^(-1)(m)$ denote the remaining optimal medoids. Since $|M_0| + |M_1| + |M_(>= 2)| = k = |M_1| + |R^*|$ and $|R^*| >= 2 |M_(>= 2)|$, we have $2 |M_0| = 2(|R^*| - |M_(>= 2)|) >= |R^*|$. Thus every $m^* in R^*$ can be paired with some $m in M_0$ so that each $m in M_0$ appears at most twice, each $m in M_1$ once, and each $m in M_(>= 2)$ zero times.

  Because $rho(M) = 0$, every pair $(m -> m^*)$ satisfies $Delta L(m -> m^*) >= 0$. Upper-bound $Delta L(m -> m^*)$ via a feasible reassignment: send every $x in C^*(m^*)$ to the newly opened $m^*$ (net change $o_x - s_x$), and send every $x in C(m) backslash C^*(m^*)$ to $phi(m^*(x))$. Observe that $phi(m^*(x)) != m$ always holds: if $m in M_0$ then $phi^(-1)(m) = emptyset$, and if $m in M_1$ then $m^*(x) != m^*$ implies $phi(m^*(x)) != m$. By the triangle inequality and minimality of $phi$,
  $ d(x, phi(m^*(x))) - s_x &<= o_x + d(m^*(x), phi(m^*(x))) - s_x \
    &<= o_x + d(m^*(x), m) - s_x \
    &<= o_x + (o_x + s_x) - s_x = 2 o_x. $
  Therefore, for every test swap $(m -> m^*)$,
  $ 0 <= Delta L(m -> m^*) <= sum_(x in C^*(m^*)) (o_x - s_x) + sum_(x in C(m)) 2 o_x. $
  Summing over all $k$ pairs, each $m^* in M^*$ occurs exactly once while each $m in M$ occurs at most twice, yielding
  $ 0 <= (L(M^*) - L(M)) + 2 sum_(x in cal(X)) 2 o_x = 5 L(M^*) - L(M), $
  which rearranges to $L(M) <= 5 L(M^*)$.
]

#subsec("C", "Expected Step Cost and Phase Transition of maxneighbor")

#theorem(title: "Exact Expected Trials and Scaling Phase Transition", origin: "Original Derivation")[
  At any vertex $M in cal(V)$ with improving density $rho = rho(M) > 0$, let $N_"eval" (M)$ be the random number of neighbor evaluations performed by CLARANS before either accepting a swap or exhausting $m_"max"$ trials. Then:
  1. $bb(E)[N_"eval" (M)] = (1 - (1 - rho)^(m_"max")) / rho <= min(m_"max", space 1 / rho)$.
  2. Escaping a bottleneck vertex having only $c = O(1)$ improving swaps ($rho = c / [k(n - k)]$) with probability at least $p_0 in (0, 1)$ requires $m_"max" >= p_0 k(n - k) / c = Omega(k n / c)$.
]

#proof[
  Let $G tilde "Geom"(rho)$ be the trial index of the first improving neighbor draw ($bb(P)(G > t) = (1 - rho)^t$). Algorithm 1 executes $N_"eval" (M) = min(G, m_"max")$ evaluations. Applying the tail-sum formula for nonnegative integers gives
  $ bb(E)[N_"eval" (M)] &= sum_(t=0)^(m_"max" - 1) bb(P)(G > t) = sum_(t=0)^(m_"max" - 1) (1 - rho)^t \
    &= (1 - (1 - rho)^(m_"max")) / rho. $
  Since $1 - (1 - rho)^(m_"max") <= 1$ and each summand $(1 - rho)^t <= 1$, Claim 1 follows immediately. For Claim 2, the escape probability satisfies $p_"esc" (M) = 1 - (1 - rho)^(m_"max") <= rho m_"max"$ by Bernoulli's inequality. Substituting $rho = c / [k(n - k)]$ yields $p_0 <= p_"esc" (M) <= c m_"max" / [k(n - k)]$, hence $m_"max" >= p_0 k(n - k) / c$.
]

Theorem 3 explains the sharp computational phase transition in Ng and Han's @ng2002clarans parameterization: during early descent when $rho(M) = Omega(1 / k)$, CLARANS advances in $O(k)$ evaluations—saving a factor of $Theta(n)$ over PAM—whereas guaranteeing escape from narrow bottlenecks forces $m_"max" = Theta(k n)$ (Ng and Han's default $m_"max" = max(250, 0.0125 k(n - k))$), transitioning per-step cost from $O(n)$ to $Theta(k n^2)$.

#subsec("D", "Adversarial Stalling Under Constant Trial Budget")

Theorem 2 requires a *true* local optimum ($rho = 0$). Holding $m_"max"$ constant breaks this guarantee arbitrarily badly.

#theorem(title: "Unbounded Stalling Error of Fixed-Budget CLARANS", origin: "Original Counterexample")[
  For any fixed trial cap $m_"max" >= 1$, restart count $r >= 1$, target ratio $M > 0$, and tolerance $eta in (0, 1)$, there exists a finite metric space $(cal(X), d)$ of $n$ distinct points (embeddable in $bb(R)^(n-1)$) with $k = 2$ on which CLARANS stalls with probability at least $1 - eta$ at a vertex $hat(M)$ satisfying $L(hat(M)) / L(M^*) > M$.
]

#proof[
  Choose $n > 2 + r(4 + 2 m_"max") / eta$ and let $cal(X) = C_1 union {y_1, y_2}$ with $|C_1| = n - 2$ distinct points. Define the metric $d$ by $d(u, u') = 1$ for distinct $u, u' in C_1$, $d(y_1, y_2) = 1$, and $d(u, y_a) = D > 1$ for all $u in C_1, a in {1, 2}$. The triangle inequality holds everywhere because $1 <= 1 + 1$, $D <= 1 + D$, and $1 <= 2 D$ (indeed, $C_1$ embeds isometrically as a regular $(n - 3)$-simplex of edge $1$ separated orthogonally from a unit segment ${y_1, y_2}$).

  The global optimum $M^* = {u_1, y_1}$ ($u_1 in C_1$) attains $L(M^*) = (n - 3) dot.c 1 + 1 = n - 2$. By contrast, *every* pair $M = {u_1, u_2} subset.eq C_1$ incurs identical cost $L(M) = (n - 4) dot.c 1 + 2 D$. Consequently, swapping either medoid of $M$ to any $u' in C_1 backslash M$ yields $Delta L = 0$ (non-improving), whereas only the $4$ swaps replacing $u_1$ or $u_2$ with $y_1$ or $y_2$ strictly improve $L(M)$. Thus $rho(M) = 4 / [2(n - 2)] = 2 / (n - 2)$.

  In a single restart, the initial vertex $M_0 tilde "Unif"(binom(cal(X), 2))$ lands inside $C_1$ with probability
  $ bb(P)(M_0 subset.eq C_1) &= binom(n - 2, 2) / binom(n, 2) \
    &= 1 - (4 n - 6) / (n(n - 1)) >= 1 - 4 / n. $
  Conditioned on $M_0 subset.eq C_1$, all $m_"max"$ trials miss the $4$ improving edges with probability $(1 - 2 / (n - 2))^(m_"max") >= 1 - 2 m_"max" / (n - 2)$. By union bound across all $r$ restarts, the probability that every restart stalls immediately at $M_0 subset.eq C_1$ is at least
  $ 1 - r (4 / n + (2 m_"max") / (n - 2)) > 1 - eta. $
  On this event, $L(hat(M)) / L(M^*) = (n - 4 + 2 D) / (n - 2) > 2 D / n$. Choosing $D = n M / 2 + 1$ guarantees $L(hat(M)) / L(M^*) > M$ with probability at least $1 - eta$.
]

#subsec("E", "FastCLARANS: Amortized O(k) Multi-Swap Scanning")

#proposition(title: "O(k) Amortized Speedup of FastCLARANS", origin: "Classical — Cited")[
  By sampling a single non-medoid $x_c tilde "Unif"(cal(X) backslash M)$ and evaluating $Delta L(m_j -> x_c)$ simultaneously for all $j = 1, dots, k$ during one $O(n)$ scan over $(c_1, d_1, d_2)$, FastCLARANS @schubert2019fasterpam @schubert2021fast inspects $k$ Johnson-graph neighbors per pass at the same $n$-distance cost as evaluating a single pair $(m_j, x_c)$ in classic CLARANS, yielding an exact $O(k)$ reduction in distance evaluations per neighborhood coverage level.
]

#proof(label: "Proof Synthesis (Cited)")[
  Precompute $Delta L^(-m_j) = sum_(i : c_1 (i) = j) (d_2 (i) - d_1 (i))$ for all $j = 1, dots, k$ in one $O(n)$ pass upon entering $M$. For a sampled candidate $x_c$, initialize an array $A[1..k] <- Delta L^(-m_j)$ and scalar $Delta_"add" <- 0$. Scanning each $x_i in cal(X)$ once computes $d_(i c) = d(x_i, x_c)$ once: if $d_(i c) < d_1 (i)$, add $(d_(i c) - d_1 (i))$ to $Delta_"add"$ and add $(d_1 (i) - d_2 (i))$ to $A[c_1 (i)]$; else if $d_(i c) < d_2 (i)$, add $(d_(i c) - d_2 (i))$ to $A[c_1 (i)]$. After one $O(n)$ sweep, $Delta L(m_j -> x_c) = A[j] + Delta_"add"$ is exact for all $j = 1, dots, k$, matching Lemma 1 while amortizing the $n$ distance computations across $k$ neighbors @schubert2021fast.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Using Lemma 1's cached $(c_1, d_1, d_2)$ structure, each neighbor trial takes $O(n d)$ time in $bb(R)^d$ (or $n$ metric oracle calls), and each accepted swap refreshes $(c_1, d_1, d_2)$ in $O(n k)$ time only for affected points (worst-case $O(n k)$ when a medoid is removed). Thus $r$ restarts with $T$ accepted swaps per restart cost $O(r dot.c T dot.c m_"max" dot.c n)$ distance evaluations and strictly $O(n d + n)$ memory—completely eliminating PAM's $O(n^2)$ matrix while retaining full access to all $n$ candidate medoids. Because the trial cap $m_"max"$ smoothly interpolates between $O(n)$ linear runtime (fixed $m_"max"$) and $O(k n^2)$ quadratic convergence (proportional $m_"max" = 0.0125 k(n - k)$ @ng2002clarans), CLARANS exposes an explicit compute-versus-quality dial. In spatial databases, Ng and Han @ng2002clarans further integrate R\*-tree bounding queries (`SD(CLARANS)` and `IR-CLARANS`) to prune I/O pages for non-adjacent clusters.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for CLARANS],
    [Time / Space Complexity], [$O(r dot.c m_"max" dot.c n dot.c T)$ evaluations ($O(n)$ per trial via Lemma 1); $O(n d + n)$ space.],
    [Geometric Inductive Bias], [Full-dataset metric Voronoi cells; dynamic medoid eligibility over all $cal(X)$.],
    [Approximation Guarantee], [Factor-$5$ global bound at true 1-swap optima (Thm. 2); stopping bound on $rho(M_tau)$ (Thm. 1).],
    [Invariances], [Strictly invariant under metric isometries and positive distance scaling $c dot.c d$.],
    [Outlier Robustness], [Within-cluster breakdown $1/2$; global $k >= 2$ hijacking possible if outlier is sampled.],
    [Recommended Companions], [FastPAM/FasterPAM (Alg. 24) for eager full-graph $O(n^2)$; CLARA (Alg. 5) for massive $n$; DBSCAN (Alg. 31) for non-convex spatial shapes.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of CLARANS.],
)

#sec("IV", "Engineering and Scientific Applications")

CLARANS was designed for large spatial and metric databases where CLARA's static subsamples sacrifice boundary accuracy yet PAM's quadratic distance table cannot fit in memory.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Spatial Data Mining & GIS Hotspot Discovery:* Ng and Han's @ng2002clarans motivating domain—identifying urban service centers, crime hotspots, and ecological zones in large spatial databases. Coupling CLARANS with R\*-trees retrieves candidate points on disk while ensuring every cluster center is a genuine geographic parcel.
2. *Graph & Road-Network Facility Siting:* On large transportation networks where pairwise shortest-path distances $d_G (u, v)$ cannot be precomputed in $O(n^2)$ space, CLARANS evaluates each candidate medoid $x_c$ via a single-source Dijkstra pass (especially potent in FastCLARANS, which amortizes one shortest-path tree across all $k$ swap tests).
3. *Chemical-Graph & Molecular Library Diversity Selection:* Pharmaceutical libraries of $10^5$–$10^6$ compounds scored by Tanimoto fingerprint dissimilarity lack meaningful vector averages. Randomized Johnson-graph search selects structurally diverse exemplar compounds across the entire chemical space without pre-filtering to a tiny subsample.
4. *Medium-to-Large Customer & Behavioral Segmentation:* On $n in [5 times 10^3, 2 times 10^5]$ customer tables with mixed Manhattan/Gower distances and overlapping cohorts, dynamic medoid search avoids CLARA's subsample projection error while staying within linear memory.
5. *Non-Euclidean Metric Pivot & Medoid-Shift Initialization:* CLARANS provides high-quality starting exemplars for metric trees and non-parametric mode-seeking methods such as Medoidshift (Alg. 42 @sheikh2007medoidshift), dramatically shortening downstream convergence.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CLARANS suffers under (a) *under-budgeted `maxneighbor` on imbalanced clusters*, where tiny improving density $rho(M) = O(1 / (k n))$ causes premature stalling (Theorem 4); (b) *massive $n > 10^6$ with proportional $m_"max" = 0.0125 k(n - k)$*, where $O(k n^2)$ distance calls become slow compared to CLARA or coreset compression @harpeled2004coresets; (c) *non-convex filaments or varying-density spatial clusters* (resolved by DBSCAN Alg. 31 @ester1996dbscan or HDBSCAN Alg. 32 @campello2013hdbscan); and (d) *high-dimensional noise* inducing distance concentration @beyer1999nearest.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why CLARANS Fits / Prerequisite], [Caution / Alternative],
    [GIS / Spatial DBs], [Full-dataset medoids; R\*-tree friendly], [Use DBSCAN (Alg. 31) for arcs],
    [Graph Logistics], [Single-source shortest paths per trial], [Use FastCLARANS (Prop. 1)],
    [Cheminformatics], [Tanimoto metric; real exemplar drugs], [Tune $m_"max"$ for rare scaffolds],
    [Massive $n > 10^6$], [Proportional $m_"max"$ turns $O(k n^2)$], [Use CLARA (Alg. 5) or coresets],
    [Skewed Minority], [Stalls if $m_"max" << k n$ (Thm. 4)], [Use FasterPAM (Alg. 24) / $D^2$],
  ),
  caption: [Application suitability and engineering matrix for CLARANS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PAM (Alg. 4), CLARA (Alg. 5), and CLARANS (Alg. 6) represent three distinct search policies over the exact same Johnson graph $cal(G)_(n, k)$, unified by the following structural characterization.

#proposition(title: "Unified Johnson-Graph Hierarchy of PAM, CLARA, and CLARANS", origin: "Original Unified Proof")[
  Let $cal(T)_"steep"$, $cal(T)_"sub"$, and $cal(T)_"rand"$ denote the transition operators of PAM, CLARA, and CLARANS on $cal(G)_(n, k) = (cal(V), cal(E))$. Then:
  1. PAM restricts neither vertices nor edges, traversing $arg min_(M' in N(M)) L(M')$ over all $k(n - k)$ incident edges;
  2. CLARA restricts the search space to static induced subgraphs $cal(G)_(S_r, k)$ of vertex cardinality $binom(s, k) << binom(n, k)$ while searching exhaustively inside each subgraph;
  3. CLARANS retains the full vertex set $cal(V)$ while sampling random edge subsets of size $<= m_"max"$, and in the limit $m_"max" -> infinity$ (with $r = 1$) converges almost surely to a true 1-swap local optimum satisfying $L(M_tau) <= 5 L(M^*)$ (Theorem 2), whereas CLARA with fixed $s < n$ and finite $R$ remains permanently excluded from all vertices $M in.not union.big_(r=1)^R binom(S_r, k)$.
]

#proof[
  Claims 1 and 2 follow directly from the definitions of PAM and CLARA. For Claim 3, set $epsilon = 1 / [2 k(n - k)] > 0$ in Theorem 1(2). Because $rho(M)$ takes values in the discrete lattice ${0, 1/[k(n-k)], dots, 1}$, the condition $rho(M_tau) <= epsilon$ implies $rho(M_tau) = 0$ identically. Since $tau <= binom(n, k) - 1$ uniformly regardless of $m_"max"$, Theorem 1(2) gives
  $ bb(P)(rho(M_tau) > 0) <= binom(n, k) exp(- m_"max" / (2 k(n - k))) -> 0 $
  as $m_"max" -> infinity$. Hence $M_tau$ is a true 1-swap local minimum with probability $1$ in the limit, immediately inheriting the factor-$5$ guarantee of Theorem 2. By contrast, CLARA can only output $M in union.big_(r=1)^R binom(S_r, k)$, which has total cardinality at most $R binom(s, k) << binom(n, k)$.
]

In empirical benchmarks on the ClusterOverview 4-dataset suite (scored by Silhouette @rousseeuw1987silhouettes and ARI), CLARANS achieves `Sil = 0.439` and `ARI = 0.676` on Iris, `Sil = 0.327` on Customers, `Sil = 0.712` (`ARI = N/A`) on Multi-Shape, and `Sil = 0.295` (`ARI = N/A`) on High-D Subspace+Noise, with a mean runtime of `0.377 ms`. Compared directly against CLARA (Alg. 5), CLARANS matches Silhouette across Iris (`0.439`), Customers (`0.327`), and Multi-Shape (`0.712`) while delivering a $+0.071$ jump in Iris ARI (`0.676` vs. `0.605`) at $2.24 times$ lower runtime (`0.377 ms` vs. `0.845 ms`). This empirical gain confirms Proposition 2: dynamic access to the full vertex set $cal(V)$ locates superior species-boundary medoids on Iris without CLARA's subsample restriction, while randomized first-improvement hill climbing reaches quality basins far faster than running multiple exhaustive PAM subproblems.
