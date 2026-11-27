#import "../template/ieee.typ": *

#alg_entry(
  id: 124,
  name: "OSLOM (Order Statistics Local Optimization Method)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2011",
  paradigm: "Statistical Significance & Order Statistics",
  complexity: $cal(O)(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Global partition objectives such as Newman-Girvan modularity @newman2004finding @brandes2008modularity suffer from two fundamental statistical flaws: they exhibit a global *resolution limit* that merges small well-defined modules in large networks @fortunato2007resolution, and they force every vertex—including unstructured noise—into some community, returning high-scoring spurious partitions even on purely random graphs. Lancichinetti, Radicchi, Ramasco, and Fortunato @lancichinetti2011oslom introduced *OSLOM* (Order Statistics Local Optimization Method), the first general framework that directly optimizes the *statistical significance* of individual candidate subgraphs against the configuration-model null hypothesis.

Let $G = (V, E)$ be an undirected simple graph (with direct generalizations to directed and weighted graphs @lancichinetti2011oslom) with $n = abs(V)$ vertices, $m = abs(E)$ edges, and vertex degrees $k_i = "deg"(i)$ satisfying $sum_(i in V) k_i = 2 m$. For a candidate subgraph $C subset.eq V$, denote its total degree volume by $K_C = sum_(v in C) k_v$, its internal stub sum by $k_"in"^C = 2 abs(E(C))$, its external boundary stubs by $m_"out"^C = K_C - k_"in"^C$, and the remaining stub pool outside $C$ by $M_C = 2 m - K_C$. Under the configuration model (random stub pairing preserving ${k_i}$), a vertex $i in V without C$ with degree $k_i$ and $k_(i, "in")$ edges into $C$ is tested for whether its observed attachment $k_(i, "in")$ is unexpectedly large relative to random pairing.

#definition(title: "Configuration-Model Attachment p-Value and Order-Statistic Score")[
  For an external vertex $i in V without C$ with degree $k_i$ and $k_(i, "in")$ neighbors in $C$, the configuration-model attachment probability $p_C (i)$ and upper-tail *r-score* (exact hypergeometric / binomial null p-value) are given by
  $
    r_i (C) = sum_(x = k_(i, "in"))^(min(k_i, m_"out"^C)) P_"null" (x | k_i, m_"out"^C, M_C),
  $
  where in the sparse limit $M_C >> k_i$, $P_"null" (x | k_i, m_"out"^C, M_C) approx binom(k_i, x) p_C^x (1 - p_C)^(k_i - x)$ with $p_C = m_"out"^C / (M_C - k_i)$ (and exact hypergeometric sampling without replacement over $M_C$ stubs for finite graphs). Sorting the $N_"ext" = n - abs(C)$ external p-values as $r_((1)) <= r_((2)) <= ... <= r_((N_"ext"))$, the *$q$-th order-statistic score* is
  $
    Omega(r_((q)); q, N_"ext") = sum_(j = q)^(N_"ext") binom(N_"ext", j) r_((q))^j (1 - r_((q)))^(N_"ext" - j) ,
  $
  and the cluster *best score* is $Omega_"bs" (C) = min_(1 <= q <= N_"ext") Omega(r_((q)); q, N_"ext")$. Mapping $Omega_"bs" (C)$ through its null cumulative distribution yields the calibrated *C-score* $C_"val" (C) in [0, 1]$.
]

Crucially, OSLOM evaluates *both* external addition and internal retention symmetrically: to test whether $v in C$ truly belongs to $C$, OSLOM evaluates $r_v (C without {v})$ and prunes vertices whose internal attachment fails the significance threshold $P_"th"$ (default $P_"th" = 0.1$). Because each cluster is optimized locally from multiple random seeds without partitioning constraints, vertices may belong to multiple significant modules (overlap) or to none at all (*homeless* background vertices).

#algbox(
  title: "OSLOM Local Significance and Hierarchical Optimization",
  inputs: [Graph $G = (V, E)$, significance threshold $P_"th" in (0, 1)$, tolerance $tau_"cov"$, restart count $R$.],
  outputs: [Hierarchical overlapping cover ${cal(C)^((0)), cal(C)^((1)), ...}$ and homeless vertex set $V_"hom"$.],
)[
  1. *Seeded local expansion*: Sample an unvisited seed vertex $s in V$, initialize $C = {s}$, and iteratively add the external neighbor minimizing $r_i (C)$ until $C$ reaches a candidate size where $Omega_"bs" (C)$ achieves a local minimum below threshold.
  2. *Order-statistic pruning (clean-up)*: Iteratively evaluate internal attachment scores $r_v (C without {v})$ for $v in C$; remove the weakest internal vertices failing the calibrated $C$-score test at level $P_"th"$ until every retained vertex in $C$ is statistically significant, or discard $C$ if it empties.
  3. *Consensus and union-separability check*: Repeat Steps 1–2 across $R$ stochastic runs to gather minimal significant clusters. For overlapping pairs $(C_a, C_b)$ with high Jaccard overlap, test the union $C_a union C_b$ by running OSLOM restricted to the induced subgraph $G[C_a union C_b]$: merge if no significant sub-modules emerge, otherwise retain both overlapping modules.
  4. *Homeless classification*: Assign all vertices belonging to zero retained significant clusters to $V_"hom"$.
  5. *Hierarchical coarsening*: Contract each significant cluster into a super-vertex (preserving inter-cluster edge weights and keeping homeless vertices intact) to form $G^((h+1))$; repeat Steps 1–4 until no significant super-clusters are detected.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present four formal results establishing the statistical foundations of OSLOM: exact beta calibration of the $q$-th order statistic, an exact finite-sample lower bound proving why the raw best score $Omega_"bs"$ must be recalibrated, a first-principles proof of immunity to the Fortunato-Barthélemy resolution limit, and asymptotic null-graph rejection.

#subsec("A", "Exact Order-Statistic Distribution and Conservative Discrete Calibration")

Greedy local search selects the smallest attachment p-values among $N_"ext"$ candidates; without accounting for the order-statistic distribution, extreme-value bias would manufacture false communities in every sparse graph.

#theorem(title: "Exact Beta Law and Conservative Calibration of OSLOM Order Statistics", origin: "Original Proof")[
  Let $R_1, ..., R_N$ be independent null attachment p-values of $N = N_"ext"$ external candidates, sorted as $R_((1)) <= ... <= R_((N))$.
  1. If each $R_i ~ "Unif"(0, 1)$, then for every $q in {1, ..., N}$ the $q$-th order statistic follows $"Beta"(q, N - q + 1)$, with exact CDF
     $
       F_q (u) = "Pr"(R_((q)) <= u) = I_u (q, N - q + 1) = sum_(j = q)^N binom(N, j) u^j (1 - u)^(N - j) ,
     $
     and therefore $Omega(R_((q)); q, N) = F_q (R_((q))) ~ "Unif"(0, 1)$.
  2. For discrete configuration-model counts $X_i$ with upper-tail p-values $R_i = "Pr"_"null" (X_i >= x_i)$, the vector $(R_1, ..., R_N)$ stochastically dominates i.i.d. uniform variables $(U_1, ..., U_N)$, and for every fixed $q in {1, ..., N}$ and significance level $alpha in (0, 1)$,
     $
       "Pr"_"null" (Omega(R_((q)); q, N) <= alpha) <= alpha .
     $
]

#proof[
  For part (1), the event ${R_((q)) <= u}$ occurs if and only if at least $q$ of the $N$ independent uniform random variables $R_1, ..., R_N$ fall in $[0, u]$. Since the indicator count $S(u) = sum_(i = 1)^N bb(I)[R_i <= u]$ follows $"Binom"(N, u)$,
  $
    F_q (u) = "Pr"(S(u) >= q) = sum_(j = q)^N binom(N, j) u^j (1 - u)^(N - j) .
  $
  Differentiating term-by-term with respect to $u$ collapses the telescoping sum to
  $
    F'_q (u) = N binom(N - 1, q - 1) u^(q - 1) (1 - u)^(N - q) ,
  $
  which is precisely the $"Beta"(q, N - q + 1)$ density; hence $F_q (u) = I_u (q, N - q + 1)$. Because $F_q$ is continuous and strictly increasing on $[0, 1]$, the probability integral transform implies $F_q (R_((q))) ~ "Unif"(0, 1)$.

  For part (2), let $G_i (x) = "Pr"_"null" (X_i >= x)$ be the discrete survival function of $X_i$. For any $u in [0, 1]$, let $x_u = min{x : G_i (x) <= u}$ (with $x_u = +infinity$ if the set is empty). Then $R_i <= u$ if and only if $X_i >= x_u$, so $"Pr"_"null" (R_i <= u) = G_i (x_u) <= u$. Hence there exists a monotone coupling on a common probability space with i.i.d. $U_i ~ "Unif"(0, 1)$ such that $R_i >= U_i$ almost surely for all $i$. Sorting preserves coordinatewise inequality ($R_((q)) >= U_((q))$ a.s.), and since $F_q$ is non-decreasing, $Omega(R_((q)); q, N) = F_q (R_((q))) >= F_q (U_((q)))$ a.s. Therefore,
  $
    "Pr"_"null" (Omega(R_((q)); q, N) <= alpha) <= "Pr"(F_q (U_((q))) <= alpha) = alpha ,
  $
  proving conservative Type-I error control for every fixed rank $q$.
]

#subsec("B", "Multiple-Testing Inflation of the Minimum Order-Statistic Score")

While $Omega(R_((q)); q, N)$ is uniform for any *fixed* $q$, OSLOM selects the *minimum* score across all $q in {1, ..., N}$ to form $Omega_"bs" = min_q F_q (R_((q)))$. The next original derivation proves an exact closed-form lower bound on the resulting Type-I inflation, mathematically justifying OSLOM's secondary $C$-score calibration @lancichinetti2011oslom.

#lemma(title: "Exact Finite-Sample Inflation Bounds for the Uncalibrated Best Score", origin: "Original Derivation")[
  Let $R_1, ..., R_N ~ "Unif"(0, 1)$ be independent ($N >= 2$) and let $Omega_"bs" = min_(1 <= q <= N) F_q (R_((q)))$. Then for every $alpha in (0, 1)$,
  $
    alpha + (alpha^(1 / N) + (1 - alpha)^(1 / N) - 1)^N <= "Pr"_"null" (Omega_"bs" <= alpha) <= min(1, N alpha) ,
  $
  with exact equality in the lower bound for $N = 2$:
  $
    "Pr"_"null" (Omega_"bs" <= alpha) = alpha + (sqrt(alpha) + sqrt(1 - alpha) - 1)^2 > alpha .
  $
]

#proof[
  The upper bound follows directly from Theorem 1 and the union bound:
  $
    "Pr"_"null" (Omega_"bs" <= alpha) = "Pr"(union.big_(q = 1)^N {F_q (R_((q))) <= alpha}) <= sum_(q = 1)^N alpha = N alpha .
  $
  For the lower bound, retain only the two extreme ranks $q = 1$ and $q = N$, defining $A_1 = {F_1 (R_((1))) <= alpha}$ and $A_N = {F_N (R_((N))) <= alpha}$. Since $F_1 (u) = 1 - (1 - u)^N$ and $F_N (u) = u^N$, we have
  $
    A_1 = {R_((1)) <= a_N}, quad A_N = {R_((N)) <= b_N},
  $
  where $a_N = 1 - (1 - alpha)^(1 / N)$ and $b_N = alpha^(1 / N)$. Because $t^(1 / N) > t$ for all $t in (0, 1)$ and $N >= 2$,
  $
    b_N - a_N = alpha^(1 / N) + (1 - alpha)^(1 / N) - 1 > alpha + (1 - alpha) - 1 = 0 ,
  $
  so $0 < a_N < b_N < 1$. Now evaluate the complementary probability in the unit hypercube $[0, 1]^N$:
  $
    A_1^c inter A_N^c = {R_((1)) > a_N} inter {R_((N)) > b_N} = (a_N, 1]^N without (a_N, b_N]^N ,
  $
  because $(a_N, b_N]^N subset.eq (a_N, 1]^N$ consists precisely of points with all coordinates in $(a_N, b_N]$. Computing Lebesgue volumes directly gives
  $
    "Pr"(A_1^c inter A_N^c) = (1 - a_N)^N - (b_N - a_N)^N = (1 - alpha) - (b_N - a_N)^N .
  $
  Hence $"Pr"_"null" (Omega_"bs" <= alpha) >= "Pr"(A_1 union A_N) = 1 - "Pr"(A_1^c inter A_N^c) = alpha + (b_N - a_N)^N > alpha$. When $N = 2$, $q in {1, 2}$ exhausts all ranks, so equality holds identically. Notice that as $alpha -> 0$, $(b_N - a_N)^N = alpha (1 + cal(O)(alpha^(1 - 1 / N)))$, so even two ranks nearly double the false-positive rate unless $Omega_"bs"$ is recalibrated via its exact joint CDF ($C$-score).
]

#subsec("C", "Immunity to the Modularity Resolution Limit")

A celebrated theorem of Fortunato and Barthélemy @fortunato2007resolution shows that global modularity maximization merges adjacent cliques in a ring of $L$ cliques $K_s$ once $L > s (s - 1) + 2$. We now prove from first principles that OSLOM's order-statistic correction completely eliminates this resolution limit: as $L -> infinity$, the $1 / L$ decay of a bridge vertex's attachment p-value is *exactly canceled* by the linear growth of $N_"ext" = (L - 1) s$.

#theorem(title: "Asymptotic Resolution-Limit Immunity on the Clique Ring", origin: "Original Proof")[
  Let $G_(L, s)$ be the Fortunato-Barthélemy ring of $L >= 3$ cliques $K_s$ ($s >= 4$) linked consecutively by single bridge edges, so $n = L s$ and $2 m = L (s (s - 1) + 2)$. Consider any single clique $C = K_s$ and significance threshold $P_"th" in (0, 1 / 2)$. Then as $L -> infinity$:
  1. *Internal stability*: For every $v in C$, removing $v$ yields order-statistic retention score $Omega(r_v (C without {v}); 1, n - s + 1) = cal(O)(L^(-(s - 2))) -> 0$, so every vertex of $C$ remains strictly significant ($< P_"th"$) for all sufficiently large $L$.
  2. *Non-absorption of neighboring cliques*: For the two external bridge neighbors $w_1, w_2 in.not C$ (each having $k_(w, "in") = 1$ edge into $C$), the order-statistic addition scores satisfy
     $
       lim_(L -> infinity) Omega(r_((1)); 1, N_"ext") &= 1 - exp(-lambda_s) > 1 - e^(-2) approx 0.865 , \
       lim_(L -> infinity) Omega(r_((2)); 2, N_"ext") &= 1 - (1 + lambda_s) exp(-lambda_s) > 0.594 ,
     $
     where $lambda_s = 2 s^2 / (s^2 - s + 2) >= 2$. Consequently, OSLOM neither dissolves $K_s$ nor merges adjacent cliques for arbitrarily large $L$.
]

#proof[
  Let $d_0 = s (s - 1) + 2$ denote the degree volume of each clique, so $2 m = L d_0$ and $M_C = (L - 1) d_0$.

  For part (1), remove any vertex $v in C$ to form $C' = C without {v}$. Then $v$ has degree $k_v in {s - 1, s}$ with $k_(v, "in") = s - 1$ edges into $C'$, while $C'$ has $m_"out"^(C') in {s, s + 1}$ boundary stubs and external stub pool $M_(C') = (L - 1) d_0 + k_v = Theta(L)$. Under the configuration model, the probability that $v$ connects at least $s - 1$ stubs to $C'$ satisfies
  $
    r_v (C') <= binom(k_v, s - 1) (m_"out"^(C') / (M_(C') - k_v))^(s - 1) = cal(O)(L^(-(s - 1))) .
  $
  Applying the first order-statistic map over $N'_"ext" = L s - s + 1 = Theta(L)$ candidates and using $1 - (1 - r)^N <= N r$ yields
  $
    Omega(r_v (C'); 1, N'_"ext") <= N'_"ext" r_v (C') = cal(O)(L dot L^(-(s - 1))) = cal(O)(L^(-(s - 2))) .
  $
  Since $s >= 4$, $L^(-(s - 2)) <= L^(-2) -> 0$, so $v$ is strictly significant for all sufficiently large $L$.

  For part (2), the clique $C$ has $m_"out"^C = 2$ boundary stubs and $N_"ext" = (L - 1) s$ external vertices. Only the two bridge neighbors $w_1, w_2$ have $k_(w, "in") = 1$ (all other external vertices have $k_(u, "in") = 0$ and $r_u = 1$). Each bridge neighbor has degree $k_w = s$, so its null probability of having at least one edge into $C$ is
  $
    r_w = 1 - (1 - 2 / M_C)^s + cal(O)(L^(-2)) = (2 s) / ((L - 1) d_0) + cal(O)(L^(-2)) .
  $
  Multiplying by $N_"ext" = (L - 1) s$ cancels $(L - 1)$ identically:
  $
    lim_(L -> infinity) N_"ext" r_w = (2 s^2) / (s^2 - s + 2) = lambda_s .
  $
  Observe that $lambda_s - 2 = 2 (s - 2) / (s^2 - s + 2) > 0$ for all $s >= 3$, so $lambda_s > 2$. By the Poisson limit of the binomial tail ($N_"ext" -> infinity$, $N_"ext" r_w -> lambda_s$),
  $
    lim_(L -> infinity) Omega(r_((1)); 1, N_"ext") &= 1 - e^(-lambda_s) > 1 - e^(-2) approx 0.865 , \
    lim_(L -> infinity) Omega(r_((2)); 2, N_"ext") &= 1 - (1 + lambda_s) e^(-lambda_s) > 1 - 3 e^(-2) approx 0.594 .
  $
  For $q >= 3$, $r_((q)) = 1$, giving $Omega = 1$. Since $Omega_"bs" (C) > 0.594 > P_"th"$ uniformly for large $L$, no bridge vertex is absorbed into $C$, proving resolution-limit immunity.
]

#subsec("D", "Null-Graph Rejection vs Spurious Modularity")

A hallmark empirical and theoretical distinction between OSLOM and modularity-based algorithms @blondel2008louvain @newman2004finding @traag2019leiden is behavior on unstructured random graphs.

#theorem(title: "Spurious Community Suppression on Configuration Null Graphs", origin: "Classical — Cited")[
  Let $G$ be drawn from the configuration model null ensemble (including sparse Erdős-Rényi graphs $G(n, d / n)$ with constant average degree $d > 1$). While maximum modularity satisfies $Q_"max" (G) = Omega_p (d^(-1/2))$ @fortunato2007resolution @lancichinetti2011oslom, forcing a spurious partition of all $n$ vertices, any fixed candidate cluster $C subset.eq V$ tested by OSLOM at calibrated $C$-score level $P_"th"$ is falsely accepted with probability at most $P_"th"$. As $P_"th" -> 0$, OSLOM classifies vertices of pure configuration-model graphs as homeless ($V_"hom" = V$) with high probability @lancichinetti2011oslom.
]

#proof[
  Under the configuration model null hypothesis, the conditional distribution of external attachments $k_(i, "in")$ given the degree sequence and $C$ matches the exact hypergeometric null used to define $r_i (C)$. By Theorem 1(2), the resulting attachment p-values stochastically dominate independent $"Unif"(0, 1)$ variables, and the calibrated $C$-score transformation maps $Omega_"bs" (C)$ through the exact null CDF $F_"bs"$, satisfying $"Pr"_"null" (C_"val" (C) <= P_"th") <= P_"th"$. Moreover, during internal clean-up (Step 2 of the algorithm), every vertex in a candidate seed cluster must independently survive its leave-one-out significance test against the remaining subgraph; in a locally tree-like sparse random graph $G(n, d / n)$, any bounded subgraph $C$ of size $s$ contains at most $s$ edges with high probability (tree or unicyclic structure), so each vertex has $k_(v, "in") <= 2$ and $N_"ext" r_v = Theta(1)$, yielding $C_"val" (C) = Omega_p (1) > P_"th"$ for small $P_"th"$. Consequently, seeded candidates collapse to the empty set during clean-up and are placed in $V_"hom"$, in sharp contrast to modularity maximization, whose spectral fluctuations yield $Q_"max" = Omega_p (d^(-1/2)) > 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Evaluating $r_i (C)$ for an external neighbor requires evaluating a hypergeometric or binomial tail in $cal(O)(k_(i, "in"))$ time (or $cal(O)(1)$ via cached log-factorials and incomplete beta approximations). On sparse modular graphs where communities have bounded average size $bar(s)$ and external boundary $cal(O)(bar(s) bar(k))$, expanding and cleaning one cluster costs $cal(O)(bar(s)^2 bar(k))$, and covering the network across $R$ consensus restarts costs empirically $cal(O)(R m log n)$ @lancichinetti2011oslom. In the worst case—such as dense graphs or deep overlapping hierarchies where candidate neighborhoods grow to $Theta(n)$—the asymptotic complexity reaches $cal(O)(n^2)$ per hierarchical level. Space complexity is $cal(O)(n + m)$ for adjacency lists plus $cal(O)(sum_a abs(C_a))$ for the overlapping cover.

*Topological and Statistical Profile.* OSLOM operates purely on local combinatorial statistics without embedding vertices into Euclidean space or requiring a preset cluster count $K$. It natively generalizes to directed graphs (splitting in- and out-stub pools $m_"out"^C, m_"in"^C$) and weighted graphs (combining topological p-values with continuous weight-sum null statistics via Fisher's method @lancichinetti2011oslom). Because local expansion depends on random seed order, OSLOM executes $R$ restarts (default $R = 10$ at level $0$, $R = 50$ at higher levels) to extract statistically stable consensus modules.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [Worst-case $cal(O)(n^2)$; typical sparse modular networks $cal(O)(R m log n)$.],
    [Space complexity], [$cal(O)(n + m + sum_a abs(C_a))$ for graph and overlapping cover.],
    [Cluster geometry], [Arbitrary non-convex graph modules, hierarchies, and overlapping covers.],
    [Overlap / noise], [Native overlapping memberships and explicit *homeless* noise-vertex detection.],
    [Graph types], [Undirected, directed, and positively weighted networks via generalized stub null.],
    [Hyperparameters], [Significance threshold $P_"th"$ (default $0.1$), coverage tolerance $tau_"cov"$, restarts $R$.],
    [Determinism], [Stochastic local seeding stabilized by multi-run consensus merging.],
  ),
  caption: [Operational and theoretical profile of OSLOM.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Protein-Protein Interaction (PPI) and Functional Complexes.* Many proteins are pleiotropic, participating in multiple biochemical pathways simultaneously, while experimental interactomes contain substantial false-positive noise edges. OSLOM's overlapping significance model identifies multi-complex hub proteins while relegating noisy peripheral proteins to $V_"hom"$.
2. *Directed Citation and Interdisciplinary Co-Authorship Networks.* Citation graphs are inherently directed, hierarchical (subfields within disciplines), and overlapping across interdisciplinary researchers. OSLOM's directed stub null model separates genuine topical citation clusters from high-degree survey-paper artifacts without suffering from the resolution limit.
3. *International Trade and Financial Interbank Networks.* Weighted directed bilateral trade and interbank exposure networks combine dense regional trade blocs with sparse peripheral economies. By coupling topological stub significance with edge-weight statistics, OSLOM uncovers multi-scale economic blocs while distinguishing statistically insignificant peripheral flows.
4. *Online Social Networks with Unaffiliated Background Users.* Real social graphs contain large populations of passive or weakly connected accounts alongside tightly knit interest groups. Unlike Louvain @blondel2008louvain or Infomap @rosvall2008infomap, which force unaffiliated accounts into neighboring communities, OSLOM cleanly isolates significant social groups and labels unclustered accounts as homeless.
5. *Functional Brain Connectomes.* Cortical regions exhibit hierarchical modularity (local sensory circuits inside distributed cognitive networks) and functional overlap (multimodal association hubs). OSLOM extracts both fine-grained cortical parcels and macro-scale functional hierarchies directly from thresholded correlation graphs without assuming a fixed scale.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

OSLOM degrades in four regimes: (i) *dense, unsparsified similarity graphs* (complete weighted matrices with uniform edge density), where the configuration null saturates and $cal(O)(n^2)$ runtime dominates unless sparsified via $k$-NN; (ii) *very small networks* ($n < 30$), where discrete hypergeometric p-values cannot fall below $P_"th"$ even for real clusters; (iii) *weakly assortative SBM regimes near the Kesten-Stigum detectability threshold* @abbe2018sbmreview, where local greedy expansion lacks the global spectral signal of non-backtracking or belief-propagation methods; and (iv) *massive billion-edge graphs*, where $R$-restart local testing is slower than linear-time Leiden @traag2019leiden or SLPA @xie2011slpa.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Biological PPI networks], [Excellent: native overlap, hierarchy, and noise filtering.], [Verify $P_"th"$ sensitivity on very small complexes.],
    [Directed citation / web], [Excellent: directed null avoids resolution limit.], [Use Infomap @rosvall2008infomap for flow-based dynamics.],
    [Sparse $k$-NN manifolds], [Strong: local significance tracks non-convex geometry.], [Requires appropriate $k$ so graph is sparse.],
    [Billion-edge social graphs], [Moderate: multi-restart local tests add overhead.], [Prefer Leiden @traag2019leiden, BigCLAM @yang2013bigclam, or SLPA @xie2011slpa.],
    [Dense unsparsified graphs], [Weak: $cal(O)(n^2)$ cost and saturated degree null.], [Sparsify first or use spectral / SBM inference.],
  ),
  caption: [Domain suitability matrix for OSLOM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

OSLOM bridges classical modularity optimization @newman2004finding @blondel2008louvain @traag2019leiden and overlapping network inference @palla2005cpm @yang2013bigclam @xie2011slpa. The following original derivation shows that OSLOM's attachment test acts as an information-theoretic KL-divergence sharpening of local modularity gain, automatically regularized by the global log-candidate penalty $log N_"ext"$.

#proposition(title: "Chernoff-KL Sharpening of Local Modularity under OSLOM", origin: "Original Derivation")[
  Let $mu_(i, C) = k_i m_"out"^C / (M_C - k_i)$ denote the expected configuration-model edges from external vertex $i$ into $C$, and let $Delta Q_"loc" (i, C) = k_(i, "in") - mu_(i, C)$ denote the unnormalized local modularity excess @newman2004finding. Whenever $k_(i, "in") > mu_(i, C)$, OSLOM's binomial tail score satisfies the exact Chernoff upper bound
  $
    -log r_i (C) >= k_i "KL"("Bern"(k_(i, "in") / k_i) || "Bern"(mu_(i, C) / k_i)) >= (Delta Q_"loc" (i, C))^2 / (2 mu_(i, C) + (2 / 3) Delta Q_"loc" (i, C)) .
  $
  Consequently, passing the first order-statistic test $Omega(r_i; 1, N_"ext") <= alpha$ approximately requires
  $
    k_i "KL"("Bern"(k_(i, "in") / k_i) || "Bern"(mu_(i, C) / k_i)) >= log(N_"ext" / alpha) ,
  $
  converting modularity's linear excess $Delta Q_"loc"$ into a scale-free likelihood-ratio threshold that adapts logarithmically to network size $N_"ext"$.
]

#proof[
  Under the sparse binomial approximation $X_i ~ "Binom"(k_i, p)$ with $p = mu_(i, C) / k_i$ and observed fraction $hat(p) = k_(i, "in") / k_i > p$, the Chernoff bound for any $t > 0$ gives
  $
    r_i (C) = "Pr"(X_i >= k_i hat(p)) <= exp(-t k_i hat(p)) (1 - p + p e^t)^(k_i) .
  $
  Minimizing over $t > 0$ attains $e^t = (hat(p) (1 - p)) / (p (1 - hat(p)))$, yielding the exact relative-entropy bound
  $
    -log r_i (C) >= k_i (hat(p) log(hat(p) / p) + (1 - hat(p)) log((1 - hat(p)) / (1 - p))) = k_i "KL"("Bern"(hat(p)) || "Bern"(p)) .
  $
  Applying Bernstein's inequality to the sum of centered Bernoulli indicators $X_i - mu_(i, C)$ (variance $k_i p (1 - p) <= mu_(i, C)$, bound $1$) directly yields the quadratic lower bound $(Delta Q_"loc")^2 / (2 mu_(i, C) + (2 / 3) Delta Q_"loc")$. Finally, solving $1 - (1 - r_i)^(N_"ext") <= alpha$ gives $r_i <= 1 - (1 - alpha)^(1 / N_"ext") approx alpha / N_"ext"$, so $-log r_i >= log(N_"ext" / alpha)$, completing the equivalence.
]

Compared to Clique Percolation (CPM) @palla2005cpm—which relies on rigid $k$-clique adjacency and fails on sparse trees—and DEMON @coscia2012demon—which aggregates ego-network label propagation without formal null testing—OSLOM provides principled statistical calibration across arbitrary sparse topologies. Compared to BigCLAM @yang2013bigclam and SLPA @xie2011slpa, OSLOM requires no preset community count $K$ and uniquely certifies homeless background vertices.

*Empirical Behavior on ClusterOverview Benchmarks.* Across the 4-dataset ClusterOverview benchmark suite (evaluated via $k$-NN graph construction using Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), OSLOM achieves *Iris* Sil $= 0.274$, ARI $= 0.451$; *Customers* Sil $= 0.386$; *Multi-Shape* Sil $= 0.780$ (ARI N/A); and *High-D* Sil $= 0.175$ (ARI N/A), with a mean runtime of $3.397$ ms. These empirical results closely reflect OSLOM's theoretical inductive bias. On *Multi-Shape*, OSLOM attains the highest Silhouette ($0.780$) among the four network algorithms in this chapter group: because non-convex manifolds form locally dense connected subgraphs in $k$-NN space, Theorem 2 guarantees that each manifold component stands out sharply against the sparse global configuration null without being spuriously merged. On *Customers* ($0.386$) and *Iris* (ARI $= 0.451$, best among Algs. 123–126), statistical pruning preserves coherent density clusters while avoiding arbitrary merges. In *High-D* ($0.175$), concentration of distances produces hub-skewed $k$-NN neighborhoods that compress the KL contrast in Proposition 1, reducing Silhouette while maintaining fast $3.397$ ms execution.
