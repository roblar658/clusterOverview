#import "../template/ieee.typ": *

#alg_entry(
  id: 39,
  name: "GDBSCAN (Generalized DBSCAN)",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1998",
  paradigm: "Predicate- & Weight-Generalized Density Connectedness",
  complexity: $O(n T_"pred")$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard DBSCAN @ester1996dbscan is restricted to unweighted point objects in a metric space using Euclidean (or metric) $epsilon$-balls and raw point counts. Real-world spatial and scientific databases, however, contain *spatially extended objects* (polygons, polylines, trajectories, raster tiles) and *non-spatial attributes* (weights, spectral intensities, population mass, timestamps) for which metric point-counting fails to capture domain semantics.

*GDBSCAN (Generalized Density-Based Spatial Clustering of Applications with Noise)* @sander1998gdbscan abstracts DBSCAN to an arbitrary finite object database $cal(D) = {o_1, dots, o_n}$ by introducing two orthogonal axiomatic generalizations: (1) an arbitrary binary *neighborhood predicate* $"NPred" : cal(D) times cal(D) -> {"true", "false"}$ that is reflexive and symmetric, defining each object's neighborhood $N_"NPred" (o) = {q in cal(D) : "NPred"(o, q)}$; and (2) an arbitrary non-negative *weighted cardinality functional* $"wCard" : 2^(cal(D)) -> bb(R)_(>= 0)$ (most commonly additive weights $"wCard"(S) = sum_(q in S) w(q)$ with $w : cal(D) -> bb(R)_(>= 0)$) paired with a density threshold $"MinCard" > 0$.

#definition(title: "Generalized Density-Reachability and Clusters")[
  Given $(cal(D), "NPred", "wCard", "MinCard")$ with reflexive, symmetric $"NPred"$ @sander1998gdbscan:
  1. An object $o in cal(D)$ is a *core object* ($o in cal(D)_"core"$) if $"wCard"(N_"NPred" (o)) >= "MinCard"$.
  2. Object $q$ is *directly density-reachable* from $o$ ($o ->_D q$) iff $q in N_"NPred" (o)$ and $o in cal(D)_"core"$.
  3. Object $q$ is *density-reachable* from $o$ ($o ->_D^* q$) iff there exists a chain $o = p_0, p_1, dots, p_m = q$ with $p_(t-1) ->_D p_t$ for all $t = 1, dots, m$.
  4. Objects $o, q$ are *density-connected* ($o <->_D q$) iff there exists $p in cal(D)$ such that $p ->_D^* o$ and $p ->_D^* q$.
  5. A *generalized density-based cluster* $C subset.eq cal(D)$ is a non-empty subset satisfying *maximality* ($o in C and o ->_D^* q => q in C$) and *connectivity* ($forall o, q in C, space o <->_D q$). Noise is $cal(N) = cal(D) backslash union.big_j C_j$.
]

GDBSCAN discovers all generalized clusters in a single linear pass over $cal(D)$ by expanding BFS queues from each unvisited core object (Algorithm 1), evaluating $N_"NPred" (o)$ via a spatial or relational index (e.g., R\*-tree for polygon intersection or compound spatial–attribute queries).

#algbox(
  title: "GDBSCAN (Generalized Density-Based Clustering)",
  inputs: [Object database $cal(D)$, reflexive–symmetric predicate $"NPred"$, functional $"wCard"$, threshold $"MinCard" > 0$.],
  outputs: [Partition of $cal(D)$ into clusters ${C_1, dots, C_K}$ and noise $cal(N)$.],
)[
  1. Initialize label $L(o) <- "UNCLASSIFIED"$ for all $o in cal(D)$; $c <- 0$.
  2. *For each* $o in cal(D)$ *do*:
  3. #h(1em) *If* $L(o) != "UNCLASSIFIED"$ *then* *continue*.
  4. #h(1em) Retrieve $S <- N_"NPred" (o) = {q in cal(D) : "NPred"(o, q)}$.
  5. #h(1em) *If* $"wCard"(S) < "MinCard"$ *then*:
  6. #h(2em) Set $L(o) <- "NOISE"$; *continue*.
  7. #h(1em) Increment $c <- c + 1$; set $L(o) <- c$; initialize queue $Q <- emptyset$; for each $q in S backslash {o}$ with $L(q) in {"UNCLASSIFIED", "NOISE"}$, push $q$ into $Q$ if $L(q) = "UNCLASSIFIED"$ and set $L(q) <- c$.
  8. #h(1em) *While* $Q != emptyset$ *do*:
  9. #h(2em) Pop $q$ from $Q$; retrieve $S' <- N_"NPred" (q)$.
  10. #h(2em) *If* $"wCard"(S') >= "MinCard"$ *then*:
  11. #h(3em) *For each* $p in S'$ *do*:
  12. #h(4em) *If* $L(p) in {"UNCLASSIFIED", "NOISE"}$ *then*:
  13. #h(5em) *If* $L(p) = "UNCLASSIFIED"$ *then* push $p$ into $Q$.
  14. #h(5em) Set $L(p) <- c$.
  15. *Return* clusters $C_j = {o : L(o) = j}$ and noise $cal(N) = {o : L(o) = "NOISE"}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical properties of GDBSCAN: (1) an equivalence-relation theorem showing reflexivity and symmetry of $"NPred"$ suffice to make core clusters canonical for *any* `"wCard"` functional; (2) an explicit counterexample proving breakdown under asymmetric predicates; (3) an order-theoretic nesting (functoriality) theorem across refined predicates and weight functionals; (4) classical correctness and border-object order-dependence bounds; and (5) a weighted kernel density level-set equivalence with dual spatial/weight breakdown analysis.

#subsec("A", "Core Equivalence Relation and Necessity of Symmetry")

#theorem(title: "Canonical Core Equivalence Under Reflexive–Symmetric Predicates", origin: "Original Proof")[
  Let $"NPred"$ be reflexive and symmetric on $cal(D)$, and let $"wCard" : 2^(cal(D)) -> bb(R)_(>= 0)$ be an *arbitrary* set functional. Then:
  1. Restricted to $cal(D)_"core"$, density-reachability $->_D^*$ is reflexive, symmetric, and transitive—hence an *equivalence relation* coinciding identically with $<->_D$ on $cal(D)_"core"$.
  2. For every $o in cal(D)_"core"$, the set $C(o) := {q in cal(D) : o ->_D^* q}$ is a valid generalized density-based cluster, and for any $o_1, o_2 in cal(D)_"core"$, $C(o_1) = C(o_2)$ iff $o_1 <->_D o_2$.
]

#proof[
  (1) *Reflexivity:* For $o in cal(D)_"core"$, reflexivity of $"NPred"$ gives $o in N_"NPred" (o)$, hence $o ->_D o$ and $o ->_D^* o$. *Transitivity:* Concatenating two chains $o ->_D^* p$ and $p ->_D^* q$ yields $o ->_D^* q$. *Symmetry on $cal(D)_"core"$:* Let $u, v in cal(D)_"core"$ satisfy $u ->_D^* v$ via chain $u = p_0 ->_D p_1 ->_D dots ->_D p_m = v$. By Definition 1(2), every intermediate object $p_0, dots, p_(m-1)$ must belong to $cal(D)_"core"$, and by hypothesis $p_m = v in cal(D)_"core"$ as well. Because $"NPred"$ is symmetric, $"NPred"(p_(t-1), p_t)$ implies $"NPred"(p_t, p_(t-1))$ for each $t = 1, dots, m$; combined with $p_t in cal(D)_"core"$, this gives $p_t ->_D p_(t-1)$. Reversing the chain yields $v = p_m ->_D p_(m-1) ->_D dots ->_D p_0 = u$, so $v ->_D^* u$. Finally, if $u, v in cal(D)_"core"$ satisfy $u <->_D v$ via witness $p in cal(D)$, then $p ->_D^* u$ and $p ->_D^* v$ force $p in cal(D)_"core"$ (otherwise $p$ could not reach anything unless $p = u = v$). By symmetry on $cal(D)_"core"$, $u ->_D^* p$, and transitivity yields $u ->_D^* v$. Thus $->_D^*$ and $<->_D$ coincide on $cal(D)_"core"$.

  (2) Fix $o in cal(D)_"core"$. Clearly $o in C(o)$ ($C(o) != emptyset$). For *connectivity*, any $q_1, q_2 in C(o)$ satisfy $o ->_D^* q_1$ and $o ->_D^* q_2$, so $o$ witnesses $q_1 <->_D q_2$. For *maximality*, if $q_1 in C(o)$ and $q_1 ->_D^* q_2$, transitivity yields $o ->_D^* q_2$, so $q_2 in C(o)$. Finally, if $o_1, o_2 in cal(D)_"core"$ satisfy $o_1 <->_D o_2$, then by (1) both $o_1 ->_D^* o_2$ and $o_2 ->_D^* o_1$ hold; transitivity therefore forces $C(o_1) = C(o_2)$ (including border objects). Conversely, $C(o_1) = C(o_2)$ implies $o_2 in C(o_1)$, so $o_1 ->_D^* o_2$ and $o_1 <->_D o_2$. Note that no additivity or monotonicity of $"wCard"$ is required for this theorem.
]

#proposition(title: "Breakdown of Cluster Well-Definedness Under Asymmetric Predicates", origin: "Original Counterexample")[
  Drop symmetry of $"NPred"$ while preserving reflexivity. Specifically, take $cal(D) = {o_1, o_2, o_3} subset.eq bb(R)$ at coordinates $o_1 = 0, o_2 = 2, o_3 = 3$, with the asymmetric directed 1-NN predicate $"NPred"(o, q) := [q = o or q = arg min_(p != o) |o - p|]$, unit weights $"wCard"(S) = |S|$, and $"MinCard" = 2$. Then every object is core, yet:
  1. $o_1 ->_D^* o_3$ while $not (o_3 ->_D^* o_1)$, so $->_D^*$ fails to be an equivalence relation on $cal(D)_"core"$.
  2. Algorithm 1 returns *one* cluster ${o_1, o_2, o_3}$ if $o_1$ is visited first, but *two* disjoint clusters ${o_2, o_3}$ and ${o_1}$ if $o_2$ is visited first.
]

#proof[
  Compute neighborhoods directly: $N_"NPred" (o_1) = {o_1, o_2}$ (nearest neighbor of $0$ is $2$); $N_"NPred" (o_2) = {o_2, o_3}$ ($|2 - 3| = 1 < |2 - 0| = 2$); and $N_"NPred" (o_3) = {o_3, o_2}$. Every object has $"wCard"(N_"NPred" (o_i)) = 2 >= "MinCard"$, so $cal(D)_"core" = {o_1, o_2, o_3}$.

  (1) We have $o_1 ->_D o_2$ and $o_2 ->_D o_3$, hence $o_1 ->_D^* o_3$. However, $N_"NPred" (o_2) = N_"NPred" (o_3) = {o_2, o_3}$ contains no outgoing edge to $o_1$, so $not (o_2 ->_D^* o_1)$ and $not (o_3 ->_D^* o_1)$, violating symmetry.

  (2) If Algorithm 1 starts at $o_1$, line 7 sets $L(o_1) = L(o_2) = 1$ and pushes $o_2$; expanding $o_2$ absorbs $o_3$ into $C_1$, returning a *single* cluster $C_1 = {o_1, o_2, o_3}$. Conversely, if Algorithm 1 starts at $o_2$, line 7 sets $L(o_2) = L(o_3) = 1$ and halts expansion because $N_"NPred" (o_2) union N_"NPred" (o_3) = {o_2, o_3}$. Subsequently visiting $o_1$ evaluates $"wCard"(N_"NPred" (o_1)) = 2 >= "MinCard"$, creates cluster $C_2$ with $L(o_1) = 2$, and leaves $o_2$ at $L(o_2) = 1$ since $L(o_2) in.not {"UNCLASSIFIED", "NOISE"}$, returning *two* clusters ${o_2, o_3}$ and ${o_1}$. Symmetrizing $"NPred"$ via mutual-$k$NN ($"NPred"(o, q) and "NPred"(q, o)$) or symmetric-$k$NN ($"or"$) immediately restores Theorem 1.
]

#subsec("B", "Functorial Nesting Under Predicate & Weight Monotonicity")

#theorem(title: "Hierarchical Refinement in Predicate–Weight Space", origin: "Original Proof")[
  Consider two GDBSCAN parameterizations $theta_1 = ("NPred"_1, "wCard"_1, tau_1)$ and $theta_2 = ("NPred"_2, "wCard"_2, tau_2)$ on $cal(D)$ with reflexive–symmetric predicates satisfying:
  (i) $"NPred"_1 (o, q) => "NPred"_2 (o, q)$ for all $o, q in cal(D)$;
  (ii) $"wCard"_1 (S) <= "wCard"_2 (S)$ for all $S subset.eq cal(D)$, with $"wCard"_2$ monotone ($S_1 subset.eq S_2 => "wCard"_2 (S_1) <= "wCard"_2 (S_2)$); and
  (iii) $tau_1 >= tau_2$.
  Then:
  1. $cal(D)_"core"^((1)) subset.eq cal(D)_"core"^((2))$ and $cal(N)^((2)) subset.eq cal(N)^((1))$.
  2. Every core equivalence class $[o]_1$ of $theta_1$ is wholly contained in a unique core equivalence class $[o]_2$ of $theta_2$, forming a strict *hierarchical nesting* across parameter chains $theta_1 prec.eq theta_2 prec.eq dots prec.eq theta_L$.
  3. If $"wCard"$ is positively homogeneous of degree $1$ (e.g., additive weights $sum_(q in S) w(q)$), scaling $("wCard", "MinCard") |-> (c dot "wCard", c dot "MinCard")$ for any $c > 0$ leaves all clusters and noise strictly invariant.
]

#proof[
  (1) By (i), $N_("NPred"_1) (o) subset.eq N_("NPred"_2) (o)$ for every $o in cal(D)$. Applying (ii) and (iii) for any $o in cal(D)_"core"^((1))$ gives
  $ "wCard"_2 (N_("NPred"_2) (o)) &>= "wCard"_2 (N_("NPred"_1) (o)) \
    &>= "wCard"_1 (N_("NPred"_1) (o)) >= tau_1 >= tau_2, $
  hence $o in cal(D)_"core"^((2))$. Furthermore, any non-noise point $q in cal(D) backslash cal(N)^((1))$ has some $u in cal(D)_"core"^((1))$ with $"NPred"_1 (u, q)$; then $u in cal(D)_"core"^((2))$ and $"NPred"_2 (u, q)$ hold, so $q in.not cal(N)^((2))$, proving $cal(N)^((2)) subset.eq cal(N)^((1))$.

  (2) Let $u, v in cal(D)_"core"^((1))$ satisfy $u ->_(D, 1)^* v$ via a core chain $u = p_0, dots, p_m = v$. By (1), every $p_t in cal(D)_"core"^((2))$, and by (i), $"NPred"_2 (p_(t-1), p_t)$ holds for all $t$, giving $u ->_(D, 2)^* v$. Thus $[u]_1 subset.eq [u]_2$. Because distinct equivalence classes of $theta_2$ are disjoint, $[u]_2$ is unique, establishing strict hierarchical refinement.

  (3) Positive homogeneity gives $c dot "wCard"(S) >= c dot "MinCard" <=> "wCard"(S) >= "MinCard"$ for all $c > 0$, preserving $cal(D)_"core"$ and all direct-reachability edges identically.
]

#subsec("C", "Correctness and Border-Object Order Dependence")

#theorem(title: "Algorithmic Correctness Up to Border Ambiguity", origin: "Classical — Cited")[
  For any reflexive–symmetric predicate $"NPred"$ and threshold $"MinCard" > 0$, Algorithm 1 terminates after evaluating $N_"NPred" (o)$ at most once per object $o in cal(D)$ @sander1998gdbscan. Its output clusters ${C_1, dots, C_K}$ satisfy:
  1. The collection of core subsets ${C_j inter cal(D)_"core"}_(j=1)^K$ equals the exact partition of $cal(D)_"core"$ into equivalence classes of $<->_D$, independent of object scan order.
  2. The noise set $cal(N)$ is unique and independent of scan order.
  3. Only *border objects* $b in cal(D) backslash (cal(D)_"core" union cal(N))$ adjacent to two or more distinct core classes can change label across database permutations @schubert2017dbscanrevisited.
]

#proof(label: "Proof Synthesis (Cited — Sander et al. 1998; Schubert et al. 2017)")[
  Each object $o in cal(D)$ transitions from `"UNCLASSIFIED"` to either `"NOISE"` or a cluster ID $c >= 1$ at most once; if labeled `"NOISE"` in line 6, it has already been verified non-core and is never pushed onto $Q$ again even if later absorbed as a border object in line 14. Thus $N_"NPred" (o)$ is queried at most once per object, ensuring finite termination. When a new cluster $c$ starts at an unvisited core object $o$, the BFS loop explores every object reachable from $o$ along core edges, recovering $[o]_(<->_D)$ completely by Theorem 1. Any non-core object $b$ with $N_"NPred" (b) inter cal(D)_"core" = emptyset$ can never be reached from a core object and remains `"NOISE"`. Conversely, a border object $b in N_"NPred" (u_1) inter N_"NPred" (u_2)$ with $u_1 in [o_1], u_2 in [o_2]$ ($[o_1] != [o_2]$) is permanently claimed by whichever core class expands first in line 14, isolating order dependence strictly to border objects.
]

#subsec("D", "Weighted Kernel Level-Set Equivalence and Breakdown Point")

#proposition(title: "Weighted KDE Equivalence & Dual Breakdown Analysis", origin: "Original Derivation")[
  Let $cal(D) subset.eq bb(R)^d$ carry positive weights $w_i > 0$ normalized to $sum_(i=1)^n w_i = 1$, with $"NPred"(x, y) := [norm(x - y)_2 <= epsilon]$ and $"wCard"(S) := sum_(y in S) w(y)$. Then:
  1. The core set $cal(D)_"core"$ equals the empirical super-level set ${x_i in cal(D) : hat(f)_(w, epsilon) (x_i) >= tau}$ of the weighted spherical Parzen–Rosenblatt KDE
     $ hat(f)_(w, epsilon) (x) = 1 / (V_d epsilon^d) sum_(j=1)^n w_j bb(I)(norm(x - x_j)_2 <= epsilon) $
     at level $tau = "MinCard" / (V_d epsilon^d)$.
  2. *Spatial vs. Weight Breakdown:* Moving $m$ points to infinity ($norm(x_j)_2 -> infinity$) cannot alter core status of any remaining point outside their $epsilon$-balls unless the removed weight mass reaches $"MinCard"$. Conversely, corrupting the *weight* of a single object ($w_1 -> infinity$) turns $o_1$ into a spurious core point absorbing all of $N_"NPred" (o_1)$, yielding a weight breakdown point of $epsilon^*_w = 1 / n$ unless weights are upper-bounded by $w_"max" < "MinCard"$.
]

#proof[
  (1) Dividing the GDBSCAN core inequality $sum_(j : norm(x_i - x_j)_2 <= epsilon) w_j >= "MinCard"$ by $V_d epsilon^d > 0$ gives $hat(f)_(w, epsilon) (x_i) >= "MinCard" / (V_d epsilon^d) = tau$, proving exact algebraic identity.

  (2) Let $O subset.eq cal(D)$ be spatially corrupted to infinity. Any uncorrupted core point $x_i$ with $N_epsilon (x_i) inter O = emptyset$ retains identical neighborhood mass, and even if $N_epsilon (x_i) inter O != emptyset$, its core status survives provided $sum_(j in N_epsilon (x_i) backslash O) w_j >= "MinCard"$. By contrast, fixing positions and sending a single weight $w_1 -> infinity$ guarantees $"wCard"(N_"NPred" (o_1)) >= w_1 >= "MinCard"$ (making $o_1$ core regardless of isolation) and forces $"wCard"(N_"NPred" (q)) >= w_1 >= "MinCard"$ for every $q in N_"NPred" (o_1)$, bridging otherwise separated clusters touching $o_1$. Enforcing $w_j <= w_"max" < "MinCard"$ restores a positive weight breakdown threshold of at least $ceil("MinCard" / w_"max") / n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

GDBSCAN executes at most $n$ evaluations of $N_"NPred" (o)$ and `"wCard"` (Theorem 3), yielding total time $O(n T_"pred" (n) + n T_"card")$ and space $O(n)$. When $"NPred"$ is supported by a spatial access method (e.g., an R\*-tree for polygon intersection, MBR overlap, or $epsilon$-distance queries in low $d$) and `"wCard"` is additive ($T_"card" = O(|N_"NPred" (o)|)$), average query time is $T_"pred" (n) = O(log n)$, giving $O(n log n)$ overall complexity @sander1998gdbscan. Without indexing—or for opaque relational/semantic predicates—$T_"pred" (n) = O(n)$, yielding $O(n^2)$ runtime. Notice that unlike centroid methods, GDBSCAN never requires a vector space or metric triangle inequality: any reflexive–symmetric graph or geometric incidence relation suffices.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for GDBSCAN],
    [Time / Space Complexity], [$O(n T_"pred" (n))$ time ($O(n log n)$ with R\*-tree, $O(n^2)$ unindexed); $O(n)$ space.],
    [Domain Generality], [Handles polygons, polylines, trajectories, graphs, and weighted non-spatial attributes.],
    [Core Invariance], [Exact core equivalence & unique noise for any reflexive–symmetric `"NPred"` (Theorem 1).],
    [Nesting / Scale], [Hierarchical refinement under predicate/weight monotonicity; homogeneous scaling invariant (Theorem 2).],
    [Breakdown Robustness], [Spatially outlier-robust; weight breakdown $1/n$ unless clipped ($w_"max" < "MinCard"$, Prop. 2).],
    [Recommended Upgrades], [ST-DBSCAN (Alg. 40) for explicit space–time; HDBSCAN (Alg. 32) for multi-scale extraction; PDBSCAN (Alg. 37) for scale.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of GDBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

By decoupling density connectedness from Euclidean point sets, GDBSCAN unlocks density clustering for complex geometric objects and heterogeneous scientific records.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Extended Polygon & Cadastral GIS Clustering:* Sander et al. @sander1998gdbscan originally demonstrated GDBSCAN on regional economic and geographic polygon databases using $"NPred"(P_1, P_2) := [P_1 inter P_2 != emptyset and |a(P_1) - a(P_2)| <= delta]$ (shared boundary plus similar socioeconomic index) and $"wCard"(S) = sum_(P in S) "Area"(P)$. Centroid-based Euclidean metrics fail on irregular administrative polygons, whereas GDBSCAN captures topological contiguity natively.
2. *Satellite Remote Sensing & Multispectral Segmentation:* In 5D+ Earth observation imagery (2D pixel coordinates plus spectral bands), setting $"NPred"$ to spatial adjacency conjoined with spectral similarity ($norm(s_i - s_j) <= delta_s$) and weighting pixels by confidence or area recovers non-convex agricultural fields, cloud formations, and coastlines while isolating sensor speckles as noise @sander1998gdbscan.
3. *Molecular Biology & 3D Protein Surface Binding Pockets:* In structural bioinformatics, atoms or surface vertices carry physicochemical attributes (electrostatic potential, hydrophobicity) and surface-area weights. GDBSCAN clusters spatially contiguous patches sharing compatible chemical predicates into candidate ligand-binding sites regardless of pocket concavity.
4. *Spatiotemporal Event & Trajectory Analysis:* Combining spatial proximity, temporal window overlap, and heading/velocity similarity into a compound symmetric predicate enables clustering of seismic swarms, epidemic outbreaks, and GPS trajectory segments @lee2007traclus, directly spawning ST-DBSCAN (Alg. 40) @birant2007stdbscan.
5. *Weighted Customer, Retail & Network Anomaly Profiling:* Assigning transaction revenue, population count, or edge bandwidth as `"wCard"` weights lets analysts discover high-mass behavioral cohorts—detecting dense pockets of high-value activity that unweighted point counts would dismiss as sparse noise.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GDBSCAN's flexibility carries three engineering hazards: (a) *accidentally specifying an asymmetric predicate* (such as directed $k$-NN or directional temporal precedence) breaks core equivalence and makes clusters traversal-order dependent (Proposition 1); (b) *unbounded heavy-tailed weights* enable a single extreme record to forge spurious core hubs (Proposition 2); and (c) *non-indexable compound predicates* force $O(n^2)$ pairwise evaluation on large datasets.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why GDBSCAN Fits / Prerequisite], [Caution / Alternative],
    [Polygon / Vector GIS], [Exact topological intersection + mass weights], [Ensure R\*-tree MBR filter exists],
    [Remote Sensing], [Joint spatial–spectral symmetric predicates], [Normalize spectral scales carefully],
    [Weighted Customer Data], [Captures revenue/volume density mass], [Clip extreme weights ($w < "MinCard"$)],
    [Directed Graphs / kNN], [Asymmetric `"NPred"` breaks equivalence (Prop. 1)], [Symmetrize or use SNN (Alg. 36)],
    [Multi-Density Streams], [Fixed `"MinCard"` misses sparse modes], [Use HDBSCAN (Alg. 32) or OPTICS (Alg. 33)],
  ),
  caption: [Domain suitability and engineering risk matrix for GDBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

GDBSCAN acts as the *universal algebraic abstraction* unifying an entire family of density-based clustering algorithms across Category II.

#proposition(title: "Universal Subsumption of Density Clustering Subfamilies", origin: "Original Unified Proof")[
  Every algorithm below is obtained as an exact specialization of GDBSCAN $(cal(D), "NPred", "wCard", "MinCard")$ with reflexive–symmetric `"NPred"`:
  1. *DBSCAN (Alg. 31)* @ester1996dbscan: $"NPred"(x, y) = [d(x, y) <= epsilon]$, $"wCard"(S) = |S|$, $"MinCard" = "MinPts"$.
  2. *ST-DBSCAN (Alg. 40)* @birant2007stdbscan: $"NPred"(x, y) = [d_"spat" (x, y) <= epsilon_1 and d_"temp" (x, y) <= epsilon_2]$, $"wCard"(S) = |S|$ (with cluster-mean drift check as an extension).
  3. *SNN Clustering (Alg. 36)* @ertoz2003snn: $"NPred"(x, y) = [x in "kNN"_k (y) and y in "kNN"_k (x) and |"kNN"_k (x) inter "kNN"_k (y)| >= tau_"SNN"]$, $"wCard"(S) = |S|$.
  4. *Weighted Spatial/Polygon Clustering* @sander1998gdbscan: $"NPred"(P_1, P_2) = [P_1 inter P_2 != emptyset and |a(P_1) - a(P_2)| <= delta]$, $"wCard"(S) = sum_(P in S) w(P)$.
  Moreover, sweeping a totally ordered chain of predicates or thresholds via Theorem 2 generates the exact hierarchical cluster tree materialized by OPTICS (Alg. 33) @ankerst1999optics and HDBSCAN (Alg. 32) @campello2013hdbscan.
]

#proof[
  In (1), $d(x, x) = 0 <= epsilon$ and $d(x, y) = d(y, x)$ guarantee reflexivity and symmetry, and $"wCard"(N_"NPred" (x)) = |N_epsilon (x)| >= "MinPts"$ reproduces Definition 1 of @ester1996dbscan verbatim. In (2), conjunction of two reflexive–symmetric inequalities remains reflexive and symmetric, reproducing ST-DBSCAN's core spatiotemporal adjacency predicate @birant2007stdbscan. In (3), mutual $k$-NN membership together with intersection cardinality $|"kNN"_k (x) inter "kNN"_k (y)| >= tau_"SNN"$ is explicitly symmetric in $(x, y)$ and reflexive (since $|"kNN"_k (x)| = k >= tau_"SNN"$), reproducing the SNN density core-connectivity graph @ertoz2003snn. In (4), non-empty set intersection and absolute attribute difference are reflexive and symmetric, with additive mass $"wCard"$. Finally, fixing $"NPred"_epsilon (x, y) = [d(x, y) <= epsilon]$ and varying $epsilon_1 < epsilon_2$ satisfies Theorem 2's refinement hypotheses, yielding the nested cluster hierarchy extracted by OPTICS and HDBSCAN.
]

On the ClusterOverview 4-dataset empirical benchmark suite, GDBSCAN demonstrates the strongest overall balance across heterogeneous data geometries among Algorithms 37–39. On *Iris*, GDBSCAN achieves the highest ARI ($0.547$, vs. PDBSCAN $0.503$ and DBCLASD $0.314$) alongside a strong Silhouette of $0.615$. On *Customers*, GDBSCAN leads all three methods with Silhouette $= 0.435$ (vs. DBCLASD $0.414$ and PDBSCAN $0.296$), reflecting the benefit of weighted/attribute-aware neighborhood predicates when density mass varies across commercial segments. On *High-D Subspace+Noise*, GDBSCAN likewise attains the top Silhouette of $0.356$ (ARI N/A), edging out PDBSCAN ($0.352$) and DBCLASD ($0.317$). On purely geometric 2D *Multi-Shape*, GDBSCAN achieves Silhouette $= 0.577$ (ARI N/A)—substantially outperforming DBCLASD ($0.334$) while trailing crisp unweighted PDBSCAN ($0.717$) where uniform point weighting offers no extra discriminative signal. Its mean runtime of $0.911$ ms sits cleanly between DBCLASD ($0.381$ ms) and PDBSCAN ($1.079$ ms), quantifying the modest overhead of generalized predicate and weight accumulation.
