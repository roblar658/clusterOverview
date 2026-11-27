#import "../template/ieee.typ": *

#cmp_entry(
  roman: "XI",
  title: "Axiomatic Foundations, Impossibility Theorems, Concentration of Measure, and Robustness",
  scope_label: "Cross-Paradigm Meta-Theorems Governing All 200 Clustering Algorithms",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Comparison Chapters I--II unified all 200 clustering algorithms through six constructive variational, topological, and spectral lenses, while Chapters III--X benchmarked competing methods inside each of Categories I--VIII. Yet an overarching theoretical question remains: _why does no single algorithm dominate every geometry, dimension, contamination regime, and validation metric?_ The obstruction is not optimization heuristics or model capacity, but six universal mathematical laws—impossibility, functorial rigidity, measure concentration, stability-symmetry bifurcation, finite-sample breakdown ceilings, and validation convexity bias—that constrain *every* clustering operator on finite data @kleinberg2002impossibility, @ackerman2008measures, @carlsson2010characterization, @beyer1999nearest, @ben-david2006sober, @jain2010data.

This chapter proves all six cross-paradigm meta-theorems in complete, self-contained detail and maps their exact consequences across Categories I--VIII (Algorithms 1--200): Kleinberg's impossibility and an 8-family axiom classification (Section I); Ackerman--Ben-David quality measures and Carlsson--Mémoli functorial uniqueness of Single Linkage (Section II); Beyer et al.'s high-dimensional distance concentration alongside collapse and survival theorems for subspace, cosine/vMF, and SNN geometry (Section III); Ben-David et al.'s stability-symmetry dichotomy (Section IV); an original five-tier Donoho--Huber finite-sample breakdown hierarchy (Section V); closed-form convexity-bias proofs for Silhouette, CHI, DBI, and Dunn contrasted with chance-corrected ARI (Section VI); and a unified meta-law governance matrix (Section VII).

// =============================================================================
// SECTION I: KLEINBERG'S IMPOSSIBILITY THEOREM
// =============================================================================
#sec("I", "Kleinberg's Impossibility Theorem and Axiomatic Fault Lines")

At its most abstract, any flat clustering algorithm is a deterministic map from pairwise dissimilarities to set partitions, stripped of initialization heuristics and optimization dynamics @kleinberg2002impossibility.

#definition(title: "Dissimilarity Space, Partition Lattice, and Clustering Function")[
  Let $S$ be a finite set with $|S| = n >= 2$.
  1. A *distance function* (semi-metric) $d: S times S -> bb(R)_(>= 0)$ satisfies $d(x, y) = 0 <=> x = y$ and $d(x, y) = d(y, x)$. Denote by $cal(D)(S)$ all distance functions on $S$, by $cal(M)(S)$ metrics, by $cal(U)(S)$ ultrametrics ($d(x, z) <= max(d(x, y), d(y, z))$), and by $cal(E)(S)$ metrics isometrically embeddable in some Euclidean space $(bb(R)^m, norm(dot)_2)$.
  2. Denote by $cal(P)(S)$ the lattice of partitions $Gamma = {C_1, dots, C_k}$ of $S$. We say $Gamma_0$ *refines* $Gamma_1$ ($Gamma_0 prec.eq Gamma_1$) if every block of $Gamma_0$ lies in a block of $Gamma_1$, and write $Gamma_0 prec Gamma_1$ when $Gamma_0 != Gamma_1$. Write $x tilde.op_Gamma y$ if $x, y$ share a block of $Gamma$, and $x tilde.not_Gamma y$ otherwise. A *clustering function* is any map $f: cal(D)_0 (S) -> cal(P)(S)$ on a domain $cal(D)_0 (S) subset.eq cal(D)(S)$.
]

#definition(title: "Kleinberg's Three Clustering Axioms")[
  A clustering function $f: cal(D)_0 (S) -> cal(P)(S)$ @kleinberg2002impossibility satisfies:
  1. *(Scale-Invariance):* $f(alpha dot d) = f(d)$ for every $d in cal(D)_0 (S)$ and every $alpha > 0$, where $(alpha dot d)(x, y) := alpha d(x, y)$.
  2. *(Richness):* $"Range"(f) := {f(d) : d in cal(D)_0 (S)} = cal(P)(S)$. (*Near-Richness* requires $"Range"(f)$ to contain every non-trivial partition $Gamma in cal(P)^* (S) := cal(P)(S) without {{S}, {{x} : x in S}}$.)
  3. *(Consistency):* Letting $Gamma = f(d)$, if $d' in cal(D)_0 (S)$ is a *$Gamma$-transformation* of $d$—meaning $d'(x, y) <= d(x, y)$ whenever $x tilde.op_Gamma y$ and $d'(x, y) >= d(x, y)$ whenever $x tilde.not_Gamma y$—then $f(d') = Gamma$.
]

#subsec("A", "Complete Proof of Kleinberg's Impossibility Theorem")

The core obstruction is order-theoretic: Scale-Invariance and Consistency jointly force $"Range"(f)$ to be an *antichain* in $(cal(P)(S), prec.eq)$. We prove a strengthened Antichain Lemma whose witnessing metric $d^*$ is simultaneously an *ultrametric* and *isometrically Euclidean*, showing that restricting to Euclidean point clouds provides no escape.

#lemma(title: "Ultrametric–Euclidean Antichain Lemma", origin: "Classical — Cited")[
  Let $|S| = n >= 2$ and $cal(D)_0 (S) supset.eq cal(U)(S) inter cal(E)(S)$. If $f: cal(D)_0 (S) -> cal(P)(S)$ satisfies *Scale-Invariance* and *Consistency*, then $"Range"(f)$ is an antichain in $(cal(P)(S), prec.eq)$: no $Gamma_0, Gamma_1 in "Range"(f)$ satisfy $Gamma_0 prec Gamma_1$ @kleinberg2002impossibility.
]

#proof[
  Suppose for contradiction that $Gamma_0, Gamma_1 in "Range"(f)$ satisfy $Gamma_0 prec Gamma_1$, with $f(d_0) = Gamma_0$ and $f(d_1) = Gamma_1$ for $d_0, d_1 in cal(D)_0 (S)$. Let $m_r := min_(x != y) d_r (x, y) > 0$ and $M_r := max_(x, y) d_r (x, y) >= m_r$ ($r in {0, 1}$). Partition all unordered distinct pairs ${x, y} subset.eq S$ into three disjoint sets:
  - $P_1 := {{x, y} : x tilde.op_(Gamma_0) y}$ (intra-$Gamma_0$, hence intra-$Gamma_1$);
  - $P_2 := {{x, y} : x tilde.not_(Gamma_0) y, x tilde.op_(Gamma_1) y}$ (non-empty since $Gamma_0 prec Gamma_1$);
  - $P_3 := {{x, y} : x tilde.not_(Gamma_1) y}$ (inter-$Gamma_1$, hence inter-$Gamma_0$).

  Set $alpha_1 := (2 M_0) / m_1 > 0$ and three strictly ordered positive constants
  $ a := m_0 / 2, quad b := (3 / 2) M_0, \
    c := 2 alpha_1 M_1 = (4 M_0 M_1) / m_1. $
  Since $0 < m_r <= M_r$, we have $0 < a <= M_0 / 2 < b < 2 M_0 <= c$. Define $d^*: S times S -> bb(R)_(>= 0)$ by $d^*(x, x) = 0$ and, for $x != y$,
  $ d^*(x, y) := cases(
      a & "if " {x, y} in P_1,
      b & "if " {x, y} in P_2,
      c & "if " {x, y} in P_3.
    ) $

  *Step 1 ($d^* in cal(U)(S) inter cal(E)(S)$):* For any distinct $x, y, z in S$: if ${x, z} in P_1$, then $d^*(x, z) = a <= max(d^*(x, y), d^*(y, z))$; if ${x, z} in P_2$, transitivity of $tilde.op_(Gamma_0)$ forbids both ${x, y}, {y, z} in P_1$, so $max(d^*(x, y), d^*(y, z)) >= b = d^*(x, z)$; if ${x, z} in P_3$, transitivity of $tilde.op_(Gamma_1)$ forbids both ${x, y}, {y, z} in P_1 union P_2$, so $max(d^*(x, y), d^*(y, z)) = c = d^*(x, z)$. Thus $d^* in cal(U)(S)$. Indexing blocks $B_1, dots, B_(k_1)$ of $Gamma_1$, blocks $A_1, dots, A_(k_0)$ of $Gamma_0$, and elements $x_1, dots, x_n in S$ with orthonormal basis ${u_r} union {v_q} union {w_i}$ of $bb(R)^(k_1 + k_0 + n)$, map $x_i in A_q subset.eq B_r$ to
  $ phi(x_i) &:= sqrt((c^2 - b^2) / 2) u_r \
    &quad + sqrt((b^2 - a^2) / 2) v_q + sqrt(a^2 / 2) w_i. $
  All radicands are positive ($0 < a < b < c$), and orthogonality yields $norm(phi(x_i) - phi(x_j))_2^2 = a^2$ on $P_1$, $(b^2 - a^2) + a^2 = b^2$ on $P_2$, and $(c^2 - b^2) + (b^2 - a^2) + a^2 = c^2$ on $P_3$. Thus $d^* in cal(U)(S) inter cal(E)(S) subset.eq cal(D)_0 (S)$.

  *Step 2 ($Gamma_0$-transformation of $d_0$):* On $P_1$, $d^*(x, y) = m_0 / 2 < m_0 <= d_0 (x, y)$; on $P_2 union P_3$, $d^*(x, y) >= (3 / 2) M_0 > M_0 >= d_0 (x, y)$. Consistency gives $f(d^*) = f(d_0) = Gamma_0$.

  *Step 3 ($Gamma_1$-transformation of $alpha_1 d_1$):* By Scale-Invariance, $f(alpha_1 d_1) = Gamma_1$. On $P_1 union P_2$, $d^*(x, y) <= (3 / 2) M_0 < 2 M_0 = alpha_1 m_1 <= alpha_1 d_1 (x, y)$; on $P_3$, $d^*(x, y) = 2 alpha_1 M_1 > alpha_1 M_1 >= alpha_1 d_1 (x, y)$. Consistency gives $f(d^*) = f(alpha_1 d_1) = Gamma_1$, contradicting $Gamma_0 != Gamma_1$.
]

#theorem(title: "Kleinberg's Impossibility Theorem — Euclidean Strengthening", origin: "Classical — Cited")[
  For every $|S| = n >= 2$ and every domain $cal(D)_0 (S) supset.eq cal(U)(S) inter cal(E)(S)$ (including $cal(D)(S)$, $cal(M)(S)$, $cal(U)(S)$, and $n$-point subsets of $bb(R)^d$ for $d >= n - 1$), *no* clustering function $f: cal(D)_0 (S) -> cal(P)(S)$ satisfies *Scale-Invariance*, *Richness*, and *Consistency* simultaneously @kleinberg2002impossibility. For all $n >= 4$, impossibility persists even under *Near-Richness*.
]

#proof[
  If $f$ satisfies Scale-Invariance and Consistency, Lemma 1 forces $"Range"(f)$ to be an antichain. For $n >= 2$, $Gamma_"bot" := {{x} : x in S} prec {S} =: Gamma_"top"$, refuting Richness. For $n >= 4$, pick distinct $x_1, x_2 in S$ and set $Gamma_0 := {{x_1}, {x_2}, S without {x_1, x_2}}$ and $Gamma_1 := {{x_1, x_2}, S without {x_1, x_2}}$. Since $|S without {x_1, x_2}| = n - 2 >= 2$, both are non-trivial with $Gamma_0 prec Gamma_1$, refuting Near-Richness. (At $n = 3$, non-trivial partitions form a 3-element antichain, realized by 2-cluster Single Linkage.) Finally, projecting Lemma 1's $n$ embedded vertices onto their $(n - 1)$-dimensional affine hull yields an exact isometric embedding in $bb(R)^(n - 1)$.
]

#subsec("B", "Pairwise Realizability and the Failure of k-Means Consistency")

Kleinberg's axioms are *minimally incompatible*: dropping any one axiom renders the other two jointly satisfiable on metric spaces via Single Linkage (Alg. 56 @sibson1973slink).

#proposition(title: "Pairwise Realizability via Single-Linkage Stopping Rules", origin: "Classical — Cited")[
  On $cal(M)(S)$ with $|S| = n >= 3$ @kleinberg2002impossibility: (1) *Fixed-$k$ Single Linkage* (Alg. 56 halted at $k$ blocks) satisfies Scale-Invariance + Consistency (sacrificing Richness); (2) *Fixed-threshold Single Linkage* $f_r (d)$ (edges $d(x, y) <= r$, i.e., DBSCAN Alg. 31 with $"minPts" = 1, epsilon = r$) satisfies Richness + Consistency (sacrificing Scale-Invariance); and (3) *Scale-$alpha$ Single Linkage* $f_"rel" (d)$ at $alpha = 1 / 2$ (edges $d(x, y) <= M / 2$, $M = max d$) satisfies Scale-Invariance + full Richness (sacrificing Consistency).
]

#proof[
  *(1)* Scaling preserves MST edge order (Scale-Invariance). Under a $Gamma$-transformation $d'$ of $d$, all $n - k$ intra-$Gamma$ MST edges remain $<= t^*$ while all inter-$Gamma$ edges remain $>= t^*$, so Kruskal outputs $Gamma$ (Consistency).
  *(2)* Every $Gamma in cal(P)(S)$ is realized by $d_(Gamma, r / 2, 2 r)$ (Richness); $Gamma$-transformations preserve intra-$Gamma$ edges $<= r$ and keep inter-$Gamma$ pairs $> r$ (Consistency); scaling breaks $f_r$.
  *(3)* Scaling multiplies $d$ and $M$ equally (Scale-Invariance). For Richness ($n >= 3$): $Gamma_"bot"$ is attained by uniform $d(x, y) = 1$ ($M / 2 = 1 / 2 < 1$); $Gamma_"top"$ by $x_1 = 0, x_2 = 1, x_3 = 2, x_i = 1 + (i - 3) / (2 n) in bb(R)$ ($M = 2$, star edges to $x_2$ all $<= 1$); every intermediate $Gamma$ by $d_(Gamma, 1, 3)$ ($M / 2 = 3 / 2 in (1, 3)$). Consistency fails at $n >= 5$: take blocks $C_1, C_2, C_3$ with ultrametric distances $1$ internally, $3$ between $C_1, C_2$, and $3.5$ from $C_3$ ($M / 2 = 1.75 => {C_1, C_2, C_3}$); raising $d(C_3, C_1 union C_2)$ to $10$ is a $Gamma$-transformation yet raises $M' / 2 = 5 > 3$, merging $C_1$ and $C_2$.
]

Crucially, fixing $k$—as in Category I—does *not* rescue Consistency: global $k$-Means (Alg. 1) violates *both* Richness and Consistency @kleinberg2002impossibility. Below we give an elementary 13-point Euclidean proof in $bb(R)^3$ where a $Gamma$-transformation shrinks a single intra-cluster edge while leaving all 77 remaining pairwise distances *identically invariant*.

#proposition(title: "Global 2-Means Violates Consistency Under Axial Rotation", origin: "Original Counterexample")[
  Let $f_"KM,2"$ return the global minimizer of 2-Means $"SSE"(C_1, C_2) = sum_(j=1)^2 sum_(x in C_j) norm(x - mu_j)_2^2$ (Alg. 1 @lloyd1982least). On $|S| = 13$ distinct points in $bb(R)^3$, there exist a configuration $X$ with unique optimum $Gamma = {C_1, C_2}$ and a Euclidean $Gamma$-transformation $X -> X'$ shrinking strictly one intra-$C_1$ edge while preserving all other 77 pairwise distances identically, yet flipping $f_"KM,2" (X') != Gamma$.
]

#proof[
  Fix $0 < epsilon, eta < 10^(-4)$ and place 13 distinct points in $bb(R)^3$: $B := {(0, 0, r epsilon) : r in {-4, dots, 5}}$ ($10$ points on the $z$-axis), $A := {a_1 = (1, 0, 0), a_2 = (-1, 0, 0)}$, and $c_2 := (0, 0, sqrt(3 / 2))$ ($D^2 = 3 / 2$). Any optimal Voronoi 2-partition keeps $B$ intact; checking all five non-trivial topologies (up to $a_1 <-> a_2$ symmetry):
  1. $Gamma = {A union B, {c_2}}$: $"SSE"_X (Gamma) = 1^2 + 1^2 + O(epsilon^2) = 2.000 + O(epsilon^2)$.
  2. $Gamma' = {A, B union {c_2}}$: $"SSE"_X (Gamma') = 2 + (10 / 11)(3 / 2) + O(epsilon) approx 3.364$.
  3. ${{a_1}, {a_2} union B union {c_2}}$: $(11 / 12)(1 + 3 / 2) + O(epsilon) = 55 / 24 + O(epsilon) approx 2.292$.
  4. ${{a_1, c_2}, {a_2} union B}$: $5 / 4 + 10 / 11 + O(epsilon) = 95 / 44 + O(epsilon) approx 2.159$.
  5. ${A union {c_2}, B}$: $2 + (2 / 3)(3 / 2) + O(epsilon) = 3.000$.

  Hence $f_"KM,2" (X) = Gamma$ uniquely by margin $> 0.15$. Now rotate *only* $a_2$ in the $x y$-plane to $a'_2 := (cos eta, sin eta, 0)$. For every $p = (0, 0, z_p) in B union {c_2}$, axial symmetry gives
  $ norm(a'_2 - p)_2 = sqrt(cos^2 eta + sin^2 eta + z_p^2) = norm(a_2 - p)_2 $
  identically, while $norm(a_1 - a'_2)_2 = 2 sin(eta / 2) = O(eta) < 2 = norm(a_1 - a_2)_2$. Thus $X -> X'$ is an exact Euclidean $Gamma$-transformation. Yet in $X'$, using $"SSE"(U) = (2 |U|)^(-1) sum_(u, v in U) norm(u - v)_2^2$:
  $ "SSE"_(X') (Gamma) &= (2 dot 10) / 12 + O(epsilon + eta^2) approx 1.667, \
    "SSE"_(X') (Gamma') &= (10 / 11)(3 / 2) + O(epsilon + eta^2) approx 1.364, $
  so $"SSE"_(X') (Gamma') < "SSE"_(X') (Gamma)$ and $f_"KM,2" (X') != Gamma$.
]

The algebraic cause is the inverse-cardinality weight $1 / (2 |U|)$ in $"SSE"(U)$: shrinking distances inside a sub-clump $A subset.eq C_1$ lowers the SSE of splinter block $A$ with leverage $1 / (2 |A|) > 1 / (2 |C_1|)$, breaking Consistency across *all* average-distortion and likelihood partitioners—including $k$-Medians (\#15), PAM (\#4), Ward (\#59), GMM (\#81), and FCM (\#151).

#corollary(title: "Exact Extension to L1 K-Medians and Ward's Linkage", origin: "Original Proof")[
  Under the exact 13-point axial-rotation $Gamma$-transformation $X -> X'$ of Proposition 2 (rescaling $z$-height to $D = 1.15$), *both* global 2-Medians (Alg. 15, minimizing $L_2$-spatial-median sum $sum_j sum_(x in C_j) norm(x - m_j)_2$) and Ward's agglomerative linkage (Alg. 59 @ward1963hierarchical) strictly violate Kleinberg Consistency.
]

#proof[
  In $X$ with $D = 1.15$ and $epsilon < 10^(-4)$: for $Gamma = {A union B, {c_2}}$, the spatial median of $A union B$ sits at $(0, 0, O(epsilon))$ (cost $1 + 1 = 2.000$); for $Gamma' = {A, B union {c_2}}$, the spatial median of $A$ is $(0, 0, 0)$ (cost $2$) while the spatial median of $B union {c_2}$ stays inside $B$ (since $|B| = 10 > 1$), costing $D = 1.15$, total $3.15$; splitting $a_1$ away costs $1 + sqrt(1 + D^2) approx 2.523$ or $1 + D = 2.15 > 2.000$. Thus $Gamma$ is uniquely optimal in $X$. After rotating $a_2$ to $a'_2 = (cos eta, sin eta, 0)$, the median of $A union B$ remains at $B$ (cost $2.000$), whereas the median of $A$ moves to $(1, O(eta), 0)$ (cost $O(eta)$), dropping the total cost of $Gamma'$ to $0 + 1.15 + O(eta) = 1.15 < 2.000$! For Ward's method (\#59), merging two clusters $U, V$ incurs exact SSE increment $Delta_"Ward" (U, V) = ((abs(U) dot abs(V)) / (abs(U) + abs(V))) norm(mu_U - mu_V)_2^2$. In $X'$, $a_1$ and $a'_2$ merge first at height $O(eta^2)$, $B$ merges at height $O(epsilon^2)$, and then $B$ merges with ${c_2}$ at height $(10 / 11)(3 / 2) = 1.364$ strictly before merging with $A$ at height $(20 / 12)(1) = 1.667$, yielding $Gamma' != Gamma$.
]

Table I classifies all eight families (Algorithms 1--200) against Kleinberg's three axioms.

#wide_figure(
  caption: [*Table I:* Kleinberg Axiomatic Classification Across All Eight Clustering Families (Algorithms 1--200). *S-Inv.* = Scale-Invariance, *Rich.* = Richness, *Cons.* = Consistency; "(Cond.)" denotes satisfaction conditional on scale-adaptive kernel/bandwidth selection.],
  table(
    columns: (1.28fr, 0.44fr, 0.40fr, 0.44fr, 1.05fr, 1.39fr),
    align: (left, center, center, center, left, left),
    [Algorithm Family (Algs.)], [S-Inv.], [Rich.], [Cons.], [Sacrificed Axiom(s)], [Exact Mechanism & Representative Algorithms],
    [*Cat. I:* Centroid / Medoid (1--30)], [Yes], [No], [No], [Richness & Consistency], [Fixed $k$ restricts range to $k$-block antichain; Prop. 2 leverage breaks Consistency ($k$-Means \#1, PAM \#4, K-Medians \#15, Trimmed \#26); X/G-Means (\#10--11) vary $k$ via BIC/normality but violate Cons.],
    [*Cat. II:* Density / Mode (31--55)], [No (fixed $epsilon, h$)], [Yes (w/ noise)], [Cond. / No], [Scale-Inv. (or Cons. if adaptive)], [Fixed-$epsilon$ DBSCAN (\#31, $"minPts"=1$) is Rich+Cons. but breaks S-Inv.; $"minPts">1$ breaks Rich. on singletons; HDBSCAN (\#32) & SNN (\#36) restore S-Inv. via mutual reachability / ranks but break Cons.],
    [*Cat. III:* Hierarchical (56--80)], [Yes (fixed $k$)], [No (fixed $k$)], [SL: Yes; others: No], [Richness (SL); Rich. & Cons. (CL, UPGMA, Ward)], [Fixed-$k$ Single Linkage (\#56) is S-Inv.+Cons. (Prop. 1); Complete (\#57), UPGMA (\#58), Ward (\#59) violate Cons.; fixed-height cuts sacrifice S-Inv.; BIRCH (\#64) fixed threshold $T$ breaks S-Inv.],
    [*Cat. IV:* Mixture / Topic (81--102)], [Yes (full $Sigma$)], [No], [No], [Richness & Consistency], [Parametric $k$-mixtures (GMM \#81, TMM \#88, vMF \#89) fix $k$ and fail Prop. 2; DPGMM/IGMM (\#82, \#99) & BIC-MCLUST (\#93) vary $k$ but fixed prior hyperparameters break S-Inv. & Cons.],
    [*Cat. V:* Graph / Spectral (103--132)], [Cond. ($k$-NN)], [No], [No], [Richness & Consistency], [Fixed-$sigma$ RBF affinity breaks S-Inv. (self-tuning/$k$-NN graphs restore S-Inv.); Spectral NJW/Ncut (\#103--104) fix $k$ and fail Cons.; Louvain/Leiden (\#107--108) suffer modularity resolution limit @fortunato2007resolution.],
    [*Cat. VI:* Grid / Subspace (133--150)], [Cond.], [No], [No], [Scale-Inv. (fixed grid) & Consistency], [Fixed-cell grids (STING \#133, CLIQUE \#134, WaveCluster \#135) break S-Inv.; PROCLUS/ORCLUS (\#136, \#138) fix $(k, ell)$ and break Rich.+Cons.; SSC (\#150) is S-Inv. on rays but violates Cons. across subspaces.],
    [*Cat. VII:* Fuzzy / Domain (151--175)], [Yes (FCM / GK)], [No], [No], [Richness & Consistency], [FCM/GK/GG (\#151, \#153--154) inherit $k$-Means non-Consistency; PCM (\#152) fixes $eta_j$ (breaking S-Inv.); k-Shape (\#162) & DTW (\#163--164) achieve shift/warp invariance but fix $k$ and violate Cons.],
    [*Cat. VIII:* Neural / Stream (176--200)], [Cond. (BN/norm)], [No], [No], [Richness & Consistency], [DEC/IDEC/SwAV (\#181--182, \#190) fix $k$ and depend on encoder Lipschitz constants; ART (\#178) vigilance $rho$ & DenStream (\#194) radius $epsilon$ break S-Inv.; GKA/PSO/ACO (\#198--200) optimize non-consistent SSE.],
  ),
)

// =============================================================================
// SECTION II: ACKERMAN-BEN-DAVID RELAXATION & CARLSSON-MEMOLI FUNCTORIALITY
// =============================================================================
#sec("II", "Axiomatic Relaxations: Quality Measures and Dendrogram Functors")

Kleinberg's theorem identifies two overly rigid requirements in a flat clustering map $f: cal(D)(S) -> cal(P)(S)$: (i) committing to a single winning partition rather than scoring candidate partitions, and (ii) collapsing multi-scale structure into one flat cut. Relaxing either yields a rich positive axiomatic theory @ackerman2008measures, @carlsson2010characterization.

#subsec("A", "Ackerman–Ben-David Clustering Quality Measures (CQMs)")

Ackerman and Ben-David @ackerman2008measures noted that shrinking distances inside one block $C_1$ improves $Gamma$ *while potentially making a finer sub-partition of $C_1$ even better* (as in Proposition 2). Reformulating Kleinberg's axioms for real-valued *Clustering Quality Measures (CQMs)* resolves the contradiction.

#definition(title: "Clustering Quality Measure and CQM Axioms")[
  Let $cal(P)^* (S) := {Gamma in cal(P)(S) : 2 <= |Gamma| <= n - 1}$ ($|S| = n >= 3$). A *Clustering Quality Measure* $m: cal(P)^* (S) times cal(D)(S) -> bb(R)_(>= 0)$ @ackerman2008measures satisfies the four *Ackerman--Ben-David axioms* if: (1) *Scale-Invariance*: $m(Gamma, alpha d) = m(Gamma, d)$ for all $alpha > 0$; (2) *Consistency*: $m(Gamma, d') >= m(Gamma, d)$ for every $Gamma$-transformation $d'$ of $d$; (3) *Richness*: for every $Gamma in cal(P)^* (S)$ there exists $d in cal(D)(S)$ with $Gamma in arg max_(Gamma' in cal(P)^* (S)) m(Gamma', d)$; and (4) *Isomorphism-Invariance*: $m(pi(Gamma), pi(d)) = m(Gamma, d)$ for every bijection $pi: S -> S$.
]

#theorem(title: "Ackerman–Ben-David Possibility Theorem for Quality Measures", origin: "Classical — Cited")[
  For every $|S| = n >= 3$, Scale-Invariance, Consistency, Richness, and Isomorphism-Invariance are *jointly satisfiable* on $cal(U)(S) inter cal(E)(S)$ @ackerman2008measures. Letting $E_"in" (Gamma)$ and $E_"out" (Gamma)$ be the intra- and inter-cluster unordered pairs of $Gamma in cal(P)^* (S)$, both the *Dunn ratio* @dunn1974well and the *mean-ratio CQM* satisfy all four axioms simultaneously:
  $ m_"Dunn" (Gamma, d) &:= (min_(x tilde.not_Gamma y) d(x, y)) / (max_(x tilde.op_Gamma y) d(x, y)), \
    m_"avg" (Gamma, d) &:= (|E_"out" (Gamma)|^(-1) sum_(e in E_"out" (Gamma)) d(e)) \
    &quad dot (|E_"in" (Gamma)|^(-1) sum_(e in E_"in" (Gamma)) d(e))^(-1). $
  Moreover, taking $arg max_(Gamma') m_"Dunn" (Gamma', d)$ does *not* contradict Theorem 1 because CQM Consistency allows a competing partition $Gamma' != Gamma$ to increase even faster under a $Gamma$-transformation.
]

#proof[
  Since $2 <= |Gamma| <= n - 1$, both $E_"in" (Gamma)$ and $E_"out" (Gamma)$ are non-empty, giving positive denominators. Scaling $d |-> alpha d$ cancels $alpha$ in numerator and denominator, and both depend only on the multisets of intra- and inter-cluster distances (Scale- and Isomorphism-Invariance). Under any $Gamma$-transformation $d'$, pairs in $E_"out" (Gamma)$ weakly increase and pairs in $E_"in" (Gamma)$ weakly decrease, so $m_"Dunn" (Gamma, d') >= m_"Dunn" (Gamma, d)$ and $m_"avg" (Gamma, d') >= m_"avg" (Gamma, d)$ (Consistency).

  For *Richness*, given any $Gamma in cal(P)^* (S)$ and $R > 1$, let $d = d_(Gamma, 1, R) in cal(U)(S) inter cal(E)(S)$ ($1$ on $E_"in" (Gamma)$, $R$ on $E_"out" (Gamma)$), giving $m_"Dunn" (Gamma, d) = m_"avg" (Gamma, d) = R$. For any $Gamma' in cal(P)^* (S)$ with $Gamma' != Gamma$:
  - *Case 1 ($E_"in" (Gamma') inter E_"out" (Gamma) != emptyset$):* Then $max_(E_"in" (Gamma')) d = R$, so $m_"Dunn" (Gamma', d) <= R / R = 1 < R$; letting $p_"in" := abs(E_"in" (Gamma') inter E_"out" (Gamma)) / abs(E_"in" (Gamma')) > 0$, $m_"avg" (Gamma', d) <= R / (1 + p_"in" (R - 1)) < R$.
  - *Case 2 ($E_"in" (Gamma') subset.neq E_"in" (Gamma)$):* Then $min_(E_"out" (Gamma')) d = 1$ and $max_(E_"in" (Gamma')) d = 1$, giving $m_"Dunn" (Gamma', d) = 1 < R$; letting $q_"out" := abs(E_"out" (Gamma') inter E_"in" (Gamma)) / abs(E_"out" (Gamma')) > 0$, $m_"avg" (Gamma', d) = R - q_"out" (R - 1) < R$.
  Thus $Gamma$ uniquely maximizes both CQMs.

  *Resolution with Theorem 1:* On $|S| = 4$ with $Gamma = {{x_1, x_2, x_3}, {x_4}}$, intra-distance $1$, and inter-distance $3$, $m_"Dunn" (Gamma, d) = 3 > 1$ is uniquely optimal. Shrinking *only* $d'(x_1, x_2) = 0.1$ preserves $m_"Dunn" (Gamma, d') = 3 >= 3$, yet lifts $Gamma' = {{x_1, x_2}, {x_3}, {x_4}}$ to $m_"Dunn" (Gamma', d') = 1 / 0.1 = 10 > 3$.
]

#subsec("B", "Carlsson–Mémoli Functoriality and Uniqueness of Single Linkage")

The second escape replaces flat cuts with *multi-scale hierarchies* (Category III, Algs. 56--80), recasting hierarchical clustering as a functor between metric and ultrametric categories @carlsson2010characterization.

#definition(title: "Categories, Dendrograms, and Carlsson–Mémoli Axioms")[
  Let $bold("Met")^"gen"$ be the category of finite metric spaces $(X, d_X)$ with *distance-non-increasing* (1-Lipschitz) morphisms $phi: (X, d_X) -> (Y, d_Y)$ ($d_Y (phi(x), phi(x')) <= d_X (x, x')$), and let $bold("Ult")^"gen" subset.eq bold("Met")^"gen"$ be the subcategory of finite ultrametrics $(X, u_X)$ (bijectively identified with dendrograms $theta_X: [0, infinity) -> cal(P)(X)$ via $x tilde.op_(theta_X (t)) x' <=> u_X (x, x') <= t$). Write $"sep"(X, d_X) := min_(a != b) d_X (a, b)$ and $Delta_2 (delta) := ({p, q}, d(p, q) = delta)$. A hierarchical method $frak(H): (X, d_X) |-> (X, u_X)$ satisfies the *Carlsson--Mémoli axioms* @carlsson2010characterization if:
  - *(Axiom I — Two-Point Normalization):* $u_(Delta_2 (delta)) (p, q) = delta$ for all $delta > 0$.
  - *(Axiom II — Functoriality):* Every distance-non-increasing $phi: (X, d_X) -> (Y, d_Y)$ is ultrametric-non-increasing: $u_Y (phi(x), phi(x')) <= u_X (x, x')$.
  - *(Axiom III — Separation):* For $|X| >= 2$ and $x != x'$, $u_X (x, x') >= "sep"(X, d_X)$.
]

Functoriality unifies contraction monotonicity ($phi = "id"_X$), isomorphism invariance (bijective isometries), subspace coherence (injections), and coarse-graining compatibility (quotient surjections). Single Linkage (Alg. 56 @sibson1973slink, @gower1969minimum) outputs the maximal sub-dominant ultrametric
$ u_"SL"^X (x, x') := min_(x = z_0, dots, z_m = x') max_(0 <= r < m) d_X (z_r, z_(r + 1)). $

#theorem(title: "Carlsson–Mémoli Functorial Uniqueness of Single Linkage", origin: "Classical — Cited")[
  Single Linkage $frak(H)_"SL": (X, d_X) |-> (X, u_"SL"^X)$ (Alg. 56) is the *unique* hierarchical clustering method satisfying Axioms I--III @carlsson2010characterization. Furthermore, $frak(H)_"SL"$ is *excisive* ($u_"SL"^C = u_"SL"^X |_(C times C)$ for every $C in theta_X (t)$) and 1-Lipschitz in Gromov--Hausdorff distance:
  $ d_"GH" (frak(H)_"SL" (X, d_X), frak(H)_"SL" (Y, d_Y)) <= d_"GH" (X, Y). $
]

#proof[
  *Part 1 (Verification):* Concatenating minimax chains proves $u_"SL"^X$ is an ultrametric. Axiom I is immediate on $Delta_2 (delta)$. For Axiom II, pushing any chain $x = z_0, dots, z_m = x'$ forward along 1-Lipschitz $phi$ gives $d_Y (phi(z_r), phi(z_(r + 1))) <= d_X (z_r, z_(r + 1))$, whence $u_"SL"^Y (phi(x), phi(x')) <= u_"SL"^X (x, x')$. For Axiom III, any chain from $x$ to $x' != x$ contains a distinct edge of weight $>= "sep"(X, d_X)$.

  *Part 2 (Upper Bound $u_X <= u_"SL"^X$):* Let $frak(H): (X, d_X) |-> (X, u_X)$ satisfy Axioms I--III. For any distinct $a, b in X$ with $delta := d_X (a, b)$, the isometric inclusion $iota_(a b): Delta_2 (delta) -> (X, d_X)$ sending $(p, q) |-> (a, b)$ gives, by Axioms II and I,
  $ u_X (a, b) &= u_X (iota_(a b) (p), iota_(a b) (q)) \
    &<= u_(Delta_2 (delta)) (p, q) = d_X (a, b). $
  Iterating the ultrametric inequality along any chain $x = z_0, dots, z_m = x'$ gives $u_X (x, x') <= max_r u_X (z_r, z_(r + 1)) <= max_r d_X (z_r, z_(r + 1))$; minimizing over chains yields $u_X (x, x') <= u_"SL"^X (x, x')$.

  *Part 3 (Lower Bound $u_X >= u_"SL"^X$):* Fix distinct $x, x' in X$ and set $epsilon := u_"SL"^X (x, x') > 0$. Define $a tilde.op b <=> u_"SL"^X (a, b) < epsilon$ (an equivalence relation because $u_"SL"^X$ is an ultrametric). Let $Y := X slash tilde.op$, let $pi: X -> Y$ be the canonical quotient map, and equip $Y$ with the constant metric $d_Y ([a], [b]) := epsilon dot bb(1)[[a] != [b]]$. For $[a] != [b]$, $d_X (a, b) >= u_"SL"^X (a, b) >= epsilon = d_Y ([a], [b])$, so $pi$ is distance-non-increasing. Since $u_"SL"^X (x, x') = epsilon$, we have $pi(x) != pi(x')$ and $"sep"(Y, d_Y) = epsilon$. Axioms II and III yield
  $ u_X (x, x') >= u_Y (pi(x), pi(x')) >= "sep"(Y, d_Y) = epsilon = u_"SL"^X (x, x'), $
  completing $u_X = u_"SL"^X$.

  *Part 4 (Excisiveness & Stability):* For $x, x' in C in theta_X (t)$, any chain achieving $u_"SL"^X (x, x') <= t$ stays entirely inside $C$ (every intermediate vertex $z_r$ satisfies $u_"SL"^X (x, z_r) <= t$), proving $u_"SL"^C = u_"SL"^X |_(C times C)$. For any correspondence $R subset.eq X times Y$ of distortion $"dis"(R)$, lifting minimax chains across $R$ changes edge weights by at most $"dis"(R)$, giving $|u_"SL"^X (x, x') - u_"SL"^Y (y, y')| <= "dis"(R)$ and hence $d_"GH" (frak(H)_"SL" (X), frak(H)_"SL" (Y)) <= d_"GH" (X, Y)$.
]

#proposition(title: "Functorial Failure of CLINK, UPGMA, and Ward & Tightness of SLINK Chaining", origin: "Original Proof")[
  1. *(Violation of Axiom II by Category III Linkages):* Complete Linkage (CLINK \#57 @defays1977clink), Average Linkage (UPGMA \#58 @sokal1958statistical), and Ward's method (\#59 @ward1963hierarchical) all violate Carlsson--Mémoli Functoriality even for isometric inclusions $Delta_2 (delta) -> Y$ on 4-point subsets of $bb(R)$.
  2. *(Chaining Saturates Gromov--Hausdorff Stability):* Adding a midpoint $L / 2$ to $X = {0, L}$ ($Y = {0, L / 2, L}$) incurs exact perturbation $d_"GH" (X, Y) = L / 4$ while dropping Single Linkage's merge height from $L$ to $L / 2$, achieving exact equality $d_"GH" (frak(H)_"SL" (X), frak(H)_"SL" (Y)) = d_"GH" (X, Y)$—showing that chaining is the unavoidable topological price of 1-Lipschitz functoriality.
]

#proof[
  *(1)* On $Y = {0, 1, 2.1, 3.1} subset.eq bb(R)$, all three methods first merge ${0, 1}$ and ${2.1, 3.1}$ at height $1$, then merge the two pairs at $u_"CL"^Y (1, 2.1) = 3.1$, $u_"UPGMA"^Y (1, 2.1) = (1.1 + 2.1 + 2.1 + 3.1) / 4 = 2.1$, and $u_"Ward"^Y (1, 2.1) = sqrt(2 dot 2 / 4) |2.6 - 0.5| = 2.1$. Each strictly exceeds $d_Y (1, 2.1) = 1.1$, violating the sub-dominance inequality $u_Y <= d_Y$ forced by Axioms I + II (Theorem 3, Part 2) along the isometric inclusion $Delta_2 (1.1) -> Y$ sending $(p, q) |-> (1, 2.1)$.
  *(2)* Any correspondence $R subset.eq X times Y$ must pair $L / 2$ with $0$ or $L$, yielding minimum distortion $"dis"(R) = L / 2$ and $d_"GH" (X, Y) = (1 / 2) "dis"(R) = L / 4$. Meanwhile $u_"SL"^X (0, L) = L$ drops to $u_"SL"^Y (0, L) = L / 2$, giving $d_"GH" (frak(H)_"SL" (X), frak(H)_"SL" (Y)) = L / 4 = d_"GH" (X, Y)$.
]

#remark(title: "Density-Robustified Functoriality in HDBSCAN (Alg. 32)")[
  Restricting morphisms to bijective contractions restores UPGMA and CLINK only as tie-dependent set-valued functors without Gromov--Hausdorff stability @carlsson2010characterization. HDBSCAN (Alg. 32 @campello2013hdbscan) achieves the decisive practical synthesis by applying the exact Single-Linkage functor $frak(H)_"SL"$ to the density-robustified mutual reachability metric $d_"mreach-k"$, suppressing single-point chaining bridges while preserving MST functoriality.
]

// =============================================================================
// SECTION III: HIGH-DIMENSIONAL DISTANCE CONCENTRATION
// =============================================================================
#sec("III", "High-Dimensional Distance Concentration and Survival Geometry")

As ambient dimension $d -> infinity$, pairwise distances among independent random vectors lose relative contrast, undermining algorithms reliant on raw Euclidean neighborhoods or Voronoi balls. We prove Beyer et al.'s concentration theorem @beyer1999nearest, derive collapse thresholds for ambient $k$-Means (\#1) and fixed-$epsilon$ DBSCAN (\#31), and prove why Subspace (Cat. VI), Cosine/vMF (Algs. 12, 89), and SNN (Algs. 36, 52) methods survive.

#subsec("A", "Chebyshev Proof of Relative Distance Concentration")

In ambient spaces $bb(R)^d$ ($d >= 1$), let $Q_d in bb(R)^d$ be a query and $X_(1, d), dots, X_(n, d) in bb(R)^d$ ($n >= 2$) be random points with $L_p$ distances ($p > 0$) $D_i^((d)) := norm(X_(i, d) - Q_d)_p$, nearest distance $D_"min"^((d)) := min_i D_i^((d))$, and farthest distance $D_"max"^((d)) := max_i D_i^((d))$.

#theorem(title: "High-Dimensional Distance Concentration", origin: "Classical — Cited")[
  Fix $n >= 2, p > 0$. Suppose $W_i^((d)) := (D_i^((d)))^p$ ($i = 1, dots, n$) are identically distributed with $mu_d := bb(E)[W_1^((d))] > 0$ and $sigma_d^2 := "Var"(W_1^((d)))$ satisfying $lim_(d -> infinity) sigma_d^2 / mu_d^2 = 0$ @beyer1999nearest. Then
  $ (D_"max"^((d)) - D_"min"^((d))) / D_"min"^((d)) attach(->, b: p) 0 quad "as " d -> infinity. $
  If coordinate differences $Z_(i, r) := X_(i, d, r) - Q_(d, r)$ ($r = 1, dots, d$) are i.i.d. with $m_p := bb(E)[|Z_(1, 1)|^p] > 0$ and $m_(2 p) := bb(E)[|Z_(1, 1)|^(2 p)] < infinity$, then for all $epsilon > 0$ and $d >= 1$,
  $ bb(P)( (D_"max"^((d)) - D_"min"^((d))) / D_"min"^((d)) >= epsilon ) &<= (n (m_(2 p) - m_p^2)) / (d m_p^2 delta(epsilon, p)^2) \
    &= O( n / (d epsilon^2) ), $
  where $delta(epsilon, p) := ((1 + epsilon)^p - 1) / ((1 + epsilon)^p + 1) in (0, 1)$.
]

#proof[
  Set $Y_i^((d)) := W_i^((d)) / mu_d$, so $bb(E)[Y_i^((d))] = 1$ and $"Var"(Y_i^((d))) = sigma_d^2 / mu_d^2 =: v_d -> 0$. For any $epsilon > 0$, define $delta := ((1 + epsilon)^p - 1) / ((1 + epsilon)^p + 1) in (0, 1)$; direct algebra confirms
  $ (1 + delta) / (1 - delta) = (2 (1 + epsilon)^p) / 2 = (1 + epsilon)^p. $
  On the event $E_d := {forall i <= n : |Y_i^((d)) - 1| < delta}$, every $Y_i^((d))$ lies in $(1 - delta, 1 + delta)$, so $Y_"max"^((d)) / Y_"min"^((d)) < (1 + delta) / (1 - delta) = (1 + epsilon)^p$. Taking $p$-th roots (canceling $mu_d > 0$) gives $D_"max"^((d)) / D_"min"^((d)) < 1 + epsilon$, i.e., $(D_"max"^((d)) - D_"min"^((d))) / D_"min"^((d)) < epsilon$. By union bound and Chebyshev's inequality,
  $ bb(P)( (D_"max"^((d)) - D_"min"^((d))) / D_"min"^((d)) >= epsilon ) &<= sum_(i=1)^n bb(P)(|Y_i^((d)) - 1| >= delta) \
    &<= (n v_d) / delta^2 -> 0. $
  Under i.i.d. coordinates, $mu_d = d m_p$ and $sigma_d^2 = d (m_(2 p) - m_p^2)$, giving $v_d = (m_(2 p) - m_p^2) / (d m_p^2)$ and the finite-$d$ rate $O(n / (d epsilon^2))$ since $delta(epsilon, p) = (p / 2) epsilon + O(epsilon^2)$.
]

For standard Gaussian coordinates, $sigma_D / mu_D approx (1 / p) sqrt(v_d)$ evaluates to $0.756 / sqrt(d)$ for $L_1$ versus $0.707 / sqrt(d)$ for $L_2$: Manhattan K-Medians (\#15) and Minkowski-weighted $k$-Means (\#29 @deamorim2012minkowski) retain marginally better contrast than $L_2$, yet every global $L_p$ norm collapses as $O(d^(-1 / 2))$. Moreover, because radial norms concentrate in a thin shell $sqrt(d) sigma plus.minus O(sigma)$, points slightly closer to the population mean $bold(0)$ become nearest neighbors to disproportionately many queries (*hubness*), skewing directed $k$-NN in-degrees and motivating symmetric rank-overlap metrics such as SNN (\#36, \#52).

#subsec("B", "Ambient Collapse vs. Subspace, Directional, and SNN Survival")

#theorem(title: "High-Dimensional Collapse vs. Subspace, Directional, and SNN Survival", origin: "Original Unified Proof")[
  Let ambient dimension $d -> infinity$.
  1. *(Ambient $k$-Means \#1 & Fixed-$epsilon$ DBSCAN \#31 Collapse):* Consider a two-cluster sparse-signal Gaussian mixture supported on an unknown $s$-dimensional subspace $S subset.eq {1, dots, d}$ ($s << d$):
     $ X | (Y = plus.minus 1) tilde.op cal(N)(plus.minus (Delta / 2) u_S, sigma^2 I_d), \
       norm(u_S)_2 = 1, quad "supp"(u_S) subset.eq S. $
     - With $n$ training points per cluster, the empirical centroid rule $hat(w) := hat(mu)_+ - hat(mu)_-$ suffers asymptotic test error
       $ "err"(hat(w)) = Phi( - (Delta^2 / 2) / sqrt(sigma^2 Delta^2 + (2 d sigma^4) / n) ) -> 1 / 2 $
       whenever $Delta^4 / sigma^4 = o(d / n)$.
     - For any fixed $epsilon > 0$ and $d > epsilon^2 / (2 sigma^2)$, same-cluster neighbor probability vanishes exponentially:
       $ bb(P)(norm(X - X')_2 <= epsilon | "same") \
         quad <= exp( - (2 d sigma^2 - epsilon^2)^2 / (16 d sigma^4) ) -> 0, $
       labeling $100\%$ of points as noise. Even with oracle-retuned $epsilon_d$, whenever $Delta^2 = o(sqrt(d) sigma^2)$ the feasible window collapses: $sup_t |F_"same" (t) - F_"diff" (t)| = O(Delta^2 / (sqrt(d) sigma^2)) -> 0$.
  2. *(Subspace Survival — Cat. VI: CLIQUE \#134, PROCLUS \#136, SUBCLU \#137, ORCLUS \#138, SSC \#150):* Restricting distortion or density to $S$ eliminates $d - s$ noise axes, lifting pairwise SNR by $sqrt(d / s)$ and shifting the $k$-Means error exponent from $O(n Delta^4 / (d sigma^4))$ to $O(n Delta^4 / (s sigma^4))$, independent of $d$.
  3. *(Cosine / vMF Concentration Blessing — Algs. 12, 89 @dhillon2001concept, @banerjee2005vmf):* For radial-nuisance directional clusters $X_i = R_i tilde(X)_i$ ($R_i > 0$) with $tilde(X)_i = gamma mu_(y_i) + sqrt(1 - gamma^2) U_i in bb(S)^(d - 1)$ ($chevron.l mu_a, mu_b chevron.r = rho < 1$, $gamma in (0, 1)$, $U_i tilde.op "Unif"(bb(S)^(d - 2) inter mu_(y_i)^perp)$), normalizing to $bb(S)^(d - 1)$ eliminates $R_i$ identically while cosine distance $d_"cos" := 1 - chevron.l tilde(X), tilde(X)' chevron.r$ obeys
     $ bb(E)[d_"cos" | "diff"] - bb(E)[d_"cos" | "same"] = gamma^2 (1 - rho), \
       max("Var"(d_"cos" | "same"), "Var"(d_"cos" | "diff")) <= 1 / (d - 1), $
     turning high-$d$ concentration into a *blessing*: $"SNR"_"cos" >= sqrt(d - 1) gamma^2 (1 - rho) -> +infinity$ whenever $gamma^2 (1 - rho) >> d^(-1 / 2)$.
  4. *(SNN Rank Amplification — Algs. 36, 52 @ertoz2003snn, @liu2018snndpc):* Suppose $k >= 2$ balanced clusters of size $m = n / k$ suffer raw distance concentration $(D_"max" - D_"min") / D_"min" -> 0$, yet each $K$-NN set $N_K (x)$ retains excess same-cluster purity $p_"in" = 1 / k + eta$ ($eta > 0$) with uniform leakage $p_"out" = (1 - p_"in") / (k - 1)$ drawn independently across pairs. Then $s_"SNN" (x, y) := |N_K (x) inter N_K (y)|$ satisfies
     $ bb(E)[s_"SNN" | "same"] - bb(E)[s_"SNN" | "diff"] \
       quad = (K^2 / m) (p_"in" - p_"out")^2 = (K^2 / m) ((k eta) / (k - 1))^2, $
     achieving exact dimension-free relative contrast $(8 eta^2) / (1 - 4 eta^2)$ at $k = 2$ even as $(D_"max" - D_"min") / D_"min" = O(d^(-1 / 2)) -> 0$.
]

#proof[
  *Part 1:* Write $hat(w) = Delta u_S + xi$ with $xi tilde.op cal(N)(bold(0), (2 sigma^2 / n) I_d)$ and $hat(m) = (hat(mu)_+ + hat(mu)_-) / 2 tilde.op cal(N)(bold(0), (sigma^2 / (2 n)) I_d)$. For test draw $X = (Delta / 2) u_S + Z$ ($Z tilde.op cal(N)(bold(0), sigma^2 I_d)$), $chevron.l X - hat(m), hat(w) chevron.r$ has mean $Delta^2 / 2$ and variance
  $ sigma^2 (1 + 1 / (2 n))(Delta^2 + (2 d sigma^2) / n) + (Delta^2 sigma^2) / (2 n) \
    quad >= sigma^2 Delta^2 + (2 d sigma^4) / n, $
  so by the CLT as $d -> infinity$, $"err"(hat(w)) -> 1 / 2$ whenever $Delta^4 / sigma^4 = o(d / n)$. For DBSCAN, $norm(X - X')_2^2 / (2 sigma^2) | "same" tilde.op chi_d^2$; applying the standard $chi_d^2$ lower-tail bound $bb(P)(W <= d - u) <= exp(-u^2 / (4 d))$ at $u = d - epsilon^2 / (2 sigma^2) > 0$ yields the exponential noise-collapse bound. Under retuned $epsilon_d$, intra-cluster squared norms follow $2 sigma^2 chi_d^2$ while inter-cluster ones follow non-central $2 sigma^2 chi_d^2 (Delta^2 / (2 sigma^2))$ whose mean shift $Delta^2$ is dwarfed by std $2 sqrt(2 d) sigma^2$, giving $sup_t |F_"same" (t) - F_"diff" (t)| = O(Delta^2 / (sqrt(d) sigma^2)) -> 0$.

  *Part 2:* Projecting onto $S$ replaces $I_d$ by $I_s$, shrinking variance $(2 d sigma^4) / n$ to $(2 s sigma^4) / n$ (an exact $sqrt(d / s)$ SNR gain). Category VI recovers $S$ via downward-closed density pruning (CLIQUE \#134 @agrawal1998clique, SUBCLU \#137 @kailing2004subclu), local spread minimization (PROCLUS \#136 @aggarwal1999proclus, ORCLUS \#138 @aggarwal2000orclus), or $ell_1$ self-expression (SSC \#150 @elhamifar2013ssc, @soltanolkotabi2012geometric).

  *Part 3:* Normalizing $X_i / norm(X_i)_2 = tilde(X)_i$ cancels $R_i$. With $bb(E)[U] = bb(E)[V] = bold(0)$ and $"Cov"(U) = (I_d - mu_a mu_a^top) / (d - 1)$, expanding $chevron.l tilde(X), tilde(Y) chevron.r$ gives
  $ chevron.l tilde(X), tilde(Y) chevron.r &= gamma^2 chevron.l mu_a, mu_b chevron.r + (1 - gamma^2) chevron.l U, V chevron.r \
    &quad + gamma sqrt(1 - gamma^2) chevron.l mu_a, V chevron.r \
    &quad + gamma sqrt(1 - gamma^2) chevron.l U, mu_b chevron.r. $
  Taking expectation gives $gamma^2$ ($a = b$) vs. $gamma^2 rho$ ($a != b$), proving mean gap $gamma^2 (1 - rho)$. The three zero-mean fluctuation terms $A_1, A_2, A_3$ are pairwise uncorrelated. For $a = b$, $"Var"(d_"cos") = (1 - gamma^2)^2 "tr"("Cov"(U)^2) = (1 - gamma^2)^2 / (d - 1)$. For $a != b$, $"Var"(A_1) = "Var"(A_2) = gamma^2 (1 - gamma^2)(1 - rho^2) / (d - 1)$ and
  $ "Var"(A_3) &= (1 - gamma^2)^2 (d - 2 + rho^2) / (d - 1)^2 \
    &<= (1 - gamma^2)^2 / (d - 1). $
  Summing and using $1 - rho^2 <= 1$ gives $"Var"(d_"cos" | "diff") <= (1 - gamma^4) / (d - 1) <= 1 / (d - 1)$, proving $"SNR"_"cos" >= sqrt(d - 1) gamma^2 (1 - rho)$.

  *Part 4:* SNN depends solely on neighbor ranks, invariant under any monotone distance warp $g(d)$. Setting $q_"in" = (K p_"in") / m$ and $q_"out" = (K p_"out") / m$, linearity of expectation over all $k m$ points gives
  $ bb(E)[s_"SNN" | "same"] &= (K^2 / m)[p_"in"^2 + (k - 1) p_"out"^2], \
    bb(E)[s_"SNN" | "diff"] &= (K^2 / m)[2 p_"in" p_"out" + (k - 2) p_"out"^2]. $
  Subtracting collapses into $(K^2 / m)(p_"in" - p_"out")^2 = (K^2 / m)(k eta / (k - 1))^2$. At $k = 2$, $bb(E)[s_"SNN" | "diff"] = (K^2 / m)(1 / 2 - 2 eta^2)$ and the gap is $(K^2 / m)(4 eta^2)$, giving exact relative contrast $(8 eta^2) / (1 - 4 eta^2)$ independent of $d$.
]

#wide_figure(
  caption: [*Table II:* High-Dimensional Concentration & Survival Mechanisms Across Ambient, Subspace, Directional, and Secondary-Rank Regimes.],
  table(
    columns: (1.08fr, 1.05fr, 0.90fr, 1.05fr, 1.25fr),
    align: (left, left, center, left, left),
    [Geometry / Mechanism], [Effective Contrast / SNR], [Ambient $d$ Scaling], [Key Structural Prior], [Representative Algorithms],
    [Ambient $L_2$ Euclidean], [$Delta^2 / (2 sqrt(2 d) sigma^2)$], [Collapses $O(d^(-1 / 2))$], [Full-space isotropy], [$k$-Means (\#1), DBSCAN (\#31), Ward (\#59)],
    [Ambient $L_1$ / Minkowski], [$0.756 / sqrt(d)$ vs. $0.707 / sqrt(d)$], [Collapses $O(d^(-1 / 2))$], [Sparse coordinate noise], [K-Medians (\#15), MWK-Means (\#29)],
    [Subspace / Projected], [$Delta^2 / (2 sqrt(2 s) sigma^2)$], [*Dimension-free* ($s << d$)], [Cluster-local subspace $S$], [CLIQUE (\#134), PROCLUS (\#136), SUBCLU (\#137), ORCLUS (\#138)],
    [$ell_1$ Self-Expression], [Inradius-coherence ratio], [*Dimension-free*], [Union of linear subspaces], [SSC (\#150) @elhamifar2013ssc, @soltanolkotabi2012geometric],
    [Cosine / vMF ($bb(S)^(d - 1)$)], [$sqrt(d - 1) gamma^2 (1 - rho)$], [*Blessing* $O(sqrt(d))$], [Angular signal $gamma >> d^(-1 / 4)$], [Spherical $k$-Means (\#12), vMF (\#89), SwAV (\#190)],
    [Shared $K$-NN (SNN)], [$(8 eta^2) / (1 - 4 eta^2)$ ($k = 2$)], [*Dimension-free*], [Local rank purity $eta > 0$], [SNN (\#36), SNN-DPC (\#52), Spectral $k$-NN (\#103)],
  ),
)

#remark(title: "Blessing vs. Curse in Deep & Directional Clustering")[
  Theorem 5(3) explains why $ell_2$-normalized hyperspherical representations dominate modern deep clustering (SwAV \#190 @caron2020swav, Contrastive Clustering \#192 @li2021contrastive, SCAN \#189 @vangansbeke2020scan) and text corpora (Algs. 12, 89): whenever semantic alignment $gamma > 0$ stays bounded away from zero, high dimension *sharpens* angular separation as $O(sqrt(d))$. Conversely, feeding raw tabular data with $d - s$ pure-noise columns gives $gamma = O(sqrt(s / d))$, requiring Category VI subspace recovery.
]

// =============================================================================
// SECTION IV: STABILITY THEORY & MODEL SELECTION PITFALLS
// =============================================================================
#sec("IV", "Stability Theory and Model Selection Pitfalls")

Practitioners widely select the cluster count $k$—including via Consensus Clustering (Alg. 170 @monti2003consensus)—by drawing independent subsamples $S_n, S'_n$ and choosing $k$ to minimize disagreement between $cal(A)_k (S_n)$ and $cal(A)_k (S'_n)$. Ben-David, von Luxburg, and Pál @ben-david2006sober proved that *asymptotic stability diagnoses minimizer uniqueness and symmetry rather than correctness of $k$*.

#definition(title: "Clustering Distance, Risk Consistency, and Instability")[
  Let $P$ be a Borel probability measure on bounded $cal(X) subset.eq bb(R)^d$, and let $cal(F)_k$ be a compact class of measurable $k$-clusterings $C: cal(X) -> {1, dots, k}$ equipped with the *$P$-disagreement metric*
  $ d_P (C, C') := min_(pi in S_k) P({x in cal(X) : C(x) != pi(C'(x))}). $
  For algorithm $cal(A)_k$ acting on i.i.d. draws $S_n, S'_n tilde.op P^n$, the *expected instability* is
  $ "instab"(cal(A)_k, P, n) := bb(E)[d_P (cal(A)_k (S_n), cal(A)_k (S'_n))]. $
  Given population risk $R_P (C) := bb(E)_P [ell(X, C)]$ with infimum $R_P^*$ and minimizer set $cal(C)_k^* (P) := arg min_(C in cal(F)_k) R_P (C)$, algorithm $cal(A)_k$ is *risk-consistent* if $R_P (cal(A)_k (S_n)) attach(->, b: p) R_P^*$ (holding for empirical $k$-Means by Pollard @pollard1981strong).
]

#theorem(title: "Stability Convergence & Symmetry Dichotomy", origin: "Classical — Cited")[
  Let $cal(A)_k$ be risk-consistent over compact $(cal(F)_k, d_P)$ with continuous $R_P$ @ben-david2006sober.
  1. *(Unique Minimizer $==>$ Perfect Stability):* If $cal(C)_k^* (P) = {C^*}$ is unique (up to label permutation), then $lim_(n -> infinity) "instab"(cal(A)_k, P, n) = 0$.
  2. *(Symmetric Minimizer Orbit $==>$ Permanent Instability):* If a finite isometry group $G$ preserves $P$, $cal(A)_k$ is $G$-equivariant, and $cal(C)_k^* (P) = {C_1^*, dots, C_m^*}$ ($m >= 2$) is a single transitive $G$-orbit with $delta := min_(a != b) d_P (C_a^*, C_b^*) > 0$, then
     $ lim inf_(n -> infinity) "instab"(cal(A)_k, P, n) >= (m - 1) / m dot delta > 0. $
]

#proof[
  *Part 1:* For any $epsilon > 0$, $K_epsilon := {C in cal(F)_k : d_P (C, C^*) >= epsilon}$ is compact. If $K_epsilon != emptyset$, continuity of $R_P$ and uniqueness of $C^*$ give a strict risk gap $eta_epsilon := inf_(C in K_epsilon) R_P (C) - R_P^* > 0$. Risk consistency implies $bb(P)(d_P (C_n, C^*) >= epsilon) <= bb(P)(R_P (C_n) - R_P^* >= eta_epsilon) -> 0$. By triangle inequality and $d_P <= 1$,
  $ "instab"(cal(A)_k, P, n) <= 2 epsilon + 2 bb(P)(d_P (C_n, C^*) >= epsilon), $
  so $lim sup_(n -> infinity) "instab" <= 2 epsilon$; letting $epsilon -> 0^+$ yields $0$.

  *Part 2:* Fix $r in (0, delta / 2)$ and open balls $B_a := {C : d_P (C, C_a^*) < r}$. For $a != b$, $C in B_a$, $C' in B_b$, $d_P (C, C') > delta - 2 r > 0$, so $B_1, dots, B_m$ are disjoint. Compactness of $cal(F)_k without union_a B_a$ yields $eta_r > 0$ and $p_n := bb(P)(C_n in union_a B_a) -> 1$. By $G$-invariance of $P$, $G$-equivariance of $cal(A)_k$, and transitivity of $G$ on ${B_a}$, each ball has equal mass $bb(P)(C_n in B_a) = p_n / m$. By independence of $S_n, S'_n$,
  $ bb(P)(exists a != b : C_n in B_a, C'_n in B_b) = p_n^2 (m - 1) / m, $
  whence $"instab"(cal(A)_k, P, n) >= p_n^2 ((m - 1) / m)(delta - 2 r)$. Letting $n -> infinity$ ($p_n -> 1$) and then $r -> 0^+$ establishes the bound.
]

#proposition(title: "Dual Inversion of Stability Model Selection", origin: "Original Proof")[
  For global Euclidean $k$-Means $cal(A)_k$ (Alg. 1 @lloyd1982least):
  1. *(False Stability at Wrong $k$):* On the clearly separated $k^* = 2$ mixture $P = (1 / 2) "Unif"([-D - L, -D + L] times [-1, 1]) + (1 / 2) "Unif"(B((D, 0), r_0))$ ($D > 10 L$, $L > 1$, $r_0 < 1 / 2$), global 3-Means ($k = 3 != k^*$) has a *unique* population minimizer (cutting the rectangle vertically at $x_1 = -D$ and leaving the ball intact), so $"instab"(cal(A)_3, P, n) -> 0$.
  2. *(Permanent Instability at True $k^*$):* On two concentric circles $P = (1 / 2) "Unif"(S^1 (r_1)) + (1 / 2) "Unif"(S^1 (r_2))$ ($0 < r_1 < r_2$, true $k^* = 2$), $"SO"(2)$ invariance yields an entire circle of optimal diameter cuts $C_theta^*$ ($theta in [0, pi)$) and exact permanent instability $lim_(n -> infinity) "instab"(cal(A)_2, P, n) = 1 / 4 > 0$.
]

#proof[
  *Part 1:* Because $D > 10 L$, merging the two macro-clusters costs $>= D^2 / 4 >> L^2$, so the third centroid splits either the rectangle $R_L$ or the ball. Vertically bisecting $R_L$ at $x_1 = -D$ reduces horizontal variance from $L^2 / 3$ to $L^2 / 12$—saving $(1 / 2)(L^2 / 4) = L^2 / 8$, strictly beating horizontal bisection ($1 / 8 < L^2 / 8$) and any split of the ball ($<= r_0^2 / 4 < L^2 / 8$). Uniqueness of 1D uniform two-quantization yields $|cal(C)_3^* (P)| = 1$, so Theorem 6(1) forces $"instab" -> 0$.

  *Part 2:* By central radial symmetry, every optimal 2-Means separating line passes through the origin, placing centroids at $plus.minus ((r_1 + r_2) / pi) u_theta$ for any angle $theta in [0, pi)$. Rotational equivariance makes the sample cut angle $hat(theta)_n tilde.op "Unif"[0, pi)$. Two independent uniform diameters meet at acute angle $Delta := min(|theta - theta'|, pi - |theta - theta'|) tilde.op "Unif"[0, pi / 2]$, disagreeing on two opposite wedges of total arc $2 Delta$ out of $2 pi$, so $d_P (C_theta^*, C_(theta')^*) = Delta / pi$ and $lim_(n -> infinity) "instab" = bb(E)[Delta / pi] = 1 / 4$.
]

#remark(title: "Consensus Ensembles (Alg. 170) & Principled k-Selection")[
  Proposition 4 shows why Monti's Consensus Clustering (\#170 @monti2003consensus) endorses $k = 3$ on asymmetric two-cluster data while rejecting $k = 2$ on concentric rings: resampling stability measures *basin uniqueness*, not model order. Principled $k$-selection requires complexity-penalized likelihood (BIC in X-Means \#10 @pelleg2000xmeans and MCLUST \#93 @fraley2002mclust, @schwarz1978bic; MML in SNOB \#102 @wallace1968snob), Dirichlet process priors (DPGMM \#82, IGMM \#99 @rasmussen2000igmm), or density persistence (HDBSCAN \#32 @campello2013hdbscan).
]

// =============================================================================
// SECTION V: COMPARATIVE FINITE-SAMPLE BREAKDOWN POINT HIERARCHY
// =============================================================================
#sec("V", "Comparative Finite-Sample Breakdown Point Hierarchy")

To rank outlier resistance rigorously across all eight families, we adapt the Donoho--Huber finite-sample replacement breakdown point @cuesta1997trimmed to both prototype estimators (Categories I, IV, VII, VIII) and density-connectivity clusterers (Categories II, III, V).

#definition(title: "Prototype, Local-Basin, and Structural Breakdown Points")[
  Let $X = union_(j=1)^k C_j subset.eq bb(R)^d$ ($k >= 2$, $|C_j| = n_j$, $n = sum_j n_j$, $n_"min" := min_j n_j$) be a clean dataset, and let $cal(X)_m (X)$ be all $n$-point datasets $tilde(X)$ formed by replacing at most $m$ points of $X$ arbitrarily.
  1. *(Global Prototype Breakdown):* For $k$-prototype estimator $Theta(tilde(X)) = {mu_1 (tilde(X)), dots, mu_k (tilde(X))}$,
     $ epsilon_n^* (Theta; X) := min { m / n : &sup_(tilde(X) in cal(X)_m (X)) \
       &quad max_j norm(mu_j (tilde(X)))_2 = +infinity }. $
  2. *(Within-Cluster Local Breakdown $epsilon_n^"loc"$):* The smallest fraction $m / n$ capable of sending a prototype to infinity when each cluster retains its clean basin (or $k = 1$, excluding singleton-cluster hijacking).
  3. *(Density Structural Breakdown):* For noise-isolating methods, $epsilon_n^"disp"$ is the smallest replacement fraction of mutually isolated outliers ($d > epsilon$ from each other and clean clusters) that destroys a clean core cluster, and $epsilon_n^"bridge"$ is the smallest fraction of collusive outliers capable of merging two clean clusters separated by gap $D > 3 epsilon$.
]

#theorem(title: "Five-Tier Comparative Finite-Sample Breakdown Point Hierarchy", origin: "Original Proof")[
  Suppose clean clusters obey $norm(x - c_j^*)_2 <= M$ ($x in C_j$) with separation $Delta_"min" := min_(a != b) d(C_a, C_b) >> M$. Then:
  1. *(Tier 1 — Centroid, Ward, GMM, Simplex FCM):* For $k$-Means (\#1 @lloyd1982least), Ward (\#59 @ward1963hierarchical), GMM (\#81 @dempster1977em), and FCM (\#151 @bezdek1981fcm), $epsilon_n^* = epsilon_n^"loc" = 1 / n$.
  2. *(Tier 2 — $L_1$ K-Medians \#15 & PAM / FasterPAM \#4, \#24 @kaufman1990finding, @schubert2019fasterpam):* Within any fixed cluster $C_j$, replacing $m < n_j / 2$ points keeps the prototype inside radius $(2 (n_j - m) / (n_j - 2 m)) M$, giving exact local breakdown
     $ epsilon_n^"loc" = ceil(n_"min" / 2) / n approx 1 / (2 k) quad ("balanced " n_j = n / k). $
     In unrestricted global $k >= 2$ optimization, one outlier at $R > n_"min" Delta_"min" + 2 n M$ hijacks a prototype ($epsilon_n^* = 1 / n$).
  3. *(Tier 3 — Trimmed $k$-Means \#26 @cuesta1997trimmed):* For $alpha in (0, n_"min" / (2 n))$ and $(n_"min" - 2 floor(alpha n)) Delta_"min"^2 > 2 n M^2$, global breakdown is exact:
     $ epsilon_n^* ("TrimKM"_alpha; X) = (floor(alpha n) + 1) / n -> alpha < n_"min" / (2 n). $
     (For $alpha in [n_"min" / (2 n), n_"min" / n)$, tight-clump hijacking caps $epsilon_n^*$ at $ceil(n_"min" / 2) / n$.)
  4. *(Tier 4 — Redescending M-Estimators: TMM \#88 @peel2000robust & PCM \#152 @krishnapuram1993pcm):*
     - Isotropic Student-$t$ mixtures ($Sigma_j = sigma^2 I_d, nu > 0$) have bounded redescending score
       $ psi_"TMM" (delta) = ((nu + d) delta) / (nu + delta^2) <= (nu + d) / (2 sqrt(nu)) -> 0, $
       giving $O(R^(-1))$ local outlier pull and an *exponential* global hijack barrier
       $ R_"hijack"^"TMM" &>= sigma sqrt(nu) \
         &quad dot [ (1 + Delta_"min"^2 / (4 nu sigma^2))^(n_"min" / 2) - 1 ]^(1 / 2). $
     - PCM (\#152) has uniformly bounded loss $rho_(m, eta_j) (r) < eta_j$ and redescending pull $O(r^(-(m + 1) / (m - 1)))$, preventing single-point global hijacking at *any* $R -> +infinity$ (at the trade-off of coincident collapse).
  5. *(Tier 5 — Density Noise Isolation: DBSCAN \#31 @ester1996dbscan & HDBSCAN \#32 @campello2013hdbscan):* On clean core cliques of size $>= "minPts"$ separated by $D > 3 epsilon$, DBSCAN achieves dispersed breakdown $epsilon_n^"disp" = 1 - "minPts" / n$ (or $1 - k dot "minPts" / n$ for $k$ cliques), while adversarial collusive bridging obeys the sharp Vitali bracket
     $ ("minPts" / n) (ceil(D / (3 epsilon)) - 1) &<= epsilon_n^"bridge" ("DBSCAN") \
       &<= ("minPts" / n) (ceil(D / epsilon) - 1). $
]

#proof[
  *Tier 1:* Replace $x_1$ by $x_1^((R)) = R e_1$ ($R -> +infinity$), with clean bound $B := max_i norm(x_i)_2$. In $k$-Means, whichever cluster $C_j$ ($1 <= |C_j| <= n$) absorbs $x_1^((R))$ obeys $norm(mu_j)_2 >= (R - (abs(C_j) - 1) B) / abs(C_j) >= R / n - B -> +infinity$ (covering both absorption and singleton capture). In FCM ($m > 1$), $sum_j u_(1 j) = 1$ gives some $u_(1 j^*) >= 1 / k$, forcing $norm(v_(j^*))_2 >= (k^(-m) R - (n - 1) B) / (n - 1 + k^(-m)) -> +infinity$. In GMM, $x_1^((R))$ either captures a component or explodes its mean/covariance.

  *Tier 2:* Translate $c_j^* = bold(0)$ so clean $x_i in C_j$ satisfy $norm(x_i)_2 <= M$, and replace $m < n_j / 2$ points by arbitrary $z_1, dots, z_m$. For $F(mu) := sum_(i=1)^(n_j - m) norm(x_i - mu)_2 + sum_(r=1)^m norm(z_r - mu)_2$ and any $mu$ with $norm(mu)_2 = R$: clean terms satisfy $norm(x_i - mu)_2 - norm(x_i)_2 >= R - 2 M$, while triangle inequality gives $norm(z_r - mu)_2 - norm(z_r)_2 >= -R$. Summing yields
  $ F(mu) - F(bold(0)) >= (n_j - 2 m) R - 2 (n_j - m) M > 0 $
  whenever $R > (2 (n_j - m) / (n_j - 2 m)) M$ (for PAM, comparing to a clean point $x_1$ replaces $2 M$ by $3 M$; coordinate-wise $L_1$ applies identically per axis). At $m >= ceil(n_j / 2)$, moving $m$ points to $R e_1$ pulls the median to infinity, proving $epsilon_n^"loc" = ceil(n_"min" / 2) / n$. For *global* $k >= 2$ hijacking, dedicating one prototype to $x_1^((R))$ and merging two closest clean clusters costs at most $n_"min" (Delta_"min" + 2 M) + (n - 1) M$, strictly below $R - B$ for large $R$.

  *Tier 3:* Let $h := n - floor(alpha n)$. At $m = floor(alpha n) + 1$, every $h$-subset $H$ contains an outlier (pigeonhole), forcing $max_j norm(mu_j)_2 >= R / h - B -> +infinity$. Conversely, for $m <= floor(alpha n)$: *Option A* discards all $m$ escaping outliers, keeping cost $<= J_"clean" := "SSE"_k (X) <= n M^2$; *Option B* retains an escaping outlier in a bounded cluster, costing $Omega(R^2) -> +infinity$; *Option C* retains $m_1 <= m$ escaping outliers in $H^*$ using $r >= 1$ escaping prototypes, leaving $<= k - 1$ bounded prototypes for $H_0 := H^* inter X_"clean"$. Since $m$ clean points were replaced and $floor(alpha n) - (m - m_1)$ more were trimmed, at most $floor(alpha n) + m_1 <= 2 floor(alpha n) < n_"min"$ clean points are absent from $H_0$, so all $k$ clean clusters retain $>= n_"min" - 2 floor(alpha n) >= 1$ points. By pigeonhole, two distinct cluster remnants share one centroid, costing at least $(n_"min" - 2 floor(alpha n)) Delta_"min"^2 / 2 > n M^2 >= J_"clean"$, ruling out Option C.

  *Tier 4:* For isotropic TMM, $rho_"TMM" (delta) = ((nu + d) / 2) log(1 + delta^2 / nu)$ ($delta = norm(x - mu_j)_2 / sigma$) has score $psi_"TMM" (delta) = ((nu + d) delta) / (nu + delta^2)$, uniquely maximized at $delta^* = sqrt(nu)$ with peak $(nu + d) / (2 sqrt(nu))$ and $O(delta^(-1))$ decay. In EM, one outlier at radius $R$ exerts pull $<= sigma^2 (nu + d) / R -> 0$. Global hijacking requires $((nu + d) / 2) log(1 + R^2 / (nu sigma^2)) >= n_"min" ((nu + d) / 2) log(1 + Delta_"min"^2 / (4 nu sigma^2))$, yielding $R_"hijack"^"TMM" = exp(Omega(n_"min"))$. In PCM (\#152), profiled loss $rho_(m, eta_j) (r) = eta_j [1 + (eta_j / r^2)^(1 / (m - 1))]^(1 - m) < eta_j$ caps any single outlier's global penalty at $eta_j$, preventing single-point hijacking whenever $n_"min" rho_(m, eta) (Delta_"min" / 2) > eta_j$.

  *Tier 5:* Under dispersed noise, every outlier $z$ has $|B(z, epsilon) inter tilde(X)| = 1 < "minPts"$ and lands in $C_(-1)$, preserving a clean $"minPts"$-clique intact until $m = n - "minPts" + 1$ ($epsilon_n^"disp" = 1 - "minPts" / n$). For adversarial bridging ($D > 3 epsilon$): placing $"minPts"$ points at each of $ceil(D / epsilon) - 1$ evenly spaced hubs along a shortest segment proves the upper bound. Conversely, let $z_0 in C_a, z_1, dots, z_L, z_(L + 1) in C_b$ be any core chain ($norm(z_t - z_(t + 1))_2 <= epsilon$). Define a greedy *farthest-exit* subsequence $w_1, dots, w_R$:
  1. Set $t_1 := max {t <= L : d(z_t, C_a) <= epsilon} + 1$ and $w_1 := z_(t_1)$, so $d(w_1, C_a) <= 2 epsilon$ while $d(z_u, C_a) > epsilon$ for all $u >= t_1$.
  2. While $d(w_r, C_b) > 4 epsilon$, set $w_(r + 1) := z_(t_(r + 1))$ with
     $ t_(r + 1) := max { t <= L : norm(w_r - z_t)_2 <= 2 epsilon } + 1; $
     otherwise stop with $R := r$. Each step obeys $norm(w_r - w_(r + 1))_2 <= 2 epsilon + epsilon = 3 epsilon$, guaranteeing $d(w_(r + 1), C_b) >= d(w_r, C_b) - 3 epsilon > epsilon$.
  3. By maximality of $t_(r + 1) - 1$, *every* $r' > r$ satisfies $norm(w_r - w_(r'))_2 > 2 epsilon$ strictly (ruling out path self-intersection!). Hence the $R$ closed balls $B(w_r, epsilon)$ are *pairwise disjoint* and disjoint from $C_a union C_b$.
  4. By the triangle inequality, $D <= d(C_a, w_1) + sum_(r=1)^(R - 1) norm(w_r - w_(r + 1))_2 + d(w_R, C_b) <= 2 epsilon + 3 (R - 1) epsilon + 4 epsilon = (3 R + 3) epsilon$, forcing $R >= ceil(D / (3 epsilon)) - 1$. Since each $w_r$ is a core point, each disjoint ball contains $>= "minPts"$ contaminated points, giving $m >= "minPts" (ceil(D / (3 epsilon)) - 1)$.
]

#corollary(title: "Bregman–Robustness Conservation Law", origin: "Original Proof")[
  Every clustering objective minimizing an unweighted strictly convex Bregman divergence $sum_(i in C_j) D_phi (x_i || mu_j)$ over an unbounded convex domain $Omega = bb(R)^d$ (including all exponential-family hard/soft MLEs of Lens 1 in Chapter I @banerjee2005clustering) has exact breakdown point $epsilon_n^* = epsilon_n^"loc" = 1 / n$. Lifting breakdown above $1 / n$ *necessarily* sacrifices exact arithmetic-mean optimality (and hence lossless additive sufficient-statistic streaming compression as in BIRCH \#64 and CluStream \#193).
]

#proof[
  By Banerjee et al. @banerjee2005clustering (Chapter I, Theorem 1), the unique right-hand minimizer of $sum_(i in C_j) D_phi (x_i || mu_j)$ is identically the arithmetic mean $mu_j^* = |C_j|^(-1) sum_(i in C_j) x_i$, irrespective of the strictly convex generator $phi$. Replacing one point $x_1$ by $R e_1$ ($R -> +infinity$) therefore forces $norm(mu_j^*)_2 >= R / n - B -> +infinity$ verbatim as in Tier 1 of Theorem 7.
]

#remark(title: "HDBSCAN Bridge Immunity")[
  HDBSCAN (\#32 @campello2013hdbscan) elevates bridge resistance still further: narrow adversarial chains incur inflated mutual reachability core radii and negligible Excess-of-Mass persistence, causing HDBSCAN to prune bridges as noise rather than merging parent clusters. Table III summarizes the complete five-tier hierarchy.
]

#wide_figure(
  caption: [*Table III:* Comparative Finite-Sample Breakdown Point Hierarchy Across All Eight Clustering Families (Algorithms 1--200).],
  table(
    columns: (0.76fr, 1.10fr, 0.78fr, 0.86fr, 0.86fr, 1.14fr),
    align: (center, left, center, center, left, left),
    [Tier], [Representative Methods (Algs.)], [Global $epsilon_n^*$], [Local $epsilon_n^"loc"$], [Influence $psi(r)$], [Dominant Failure / Hijack Mechanism],
    [*Tier 1*], [$k$-Means (\#1), K-Means++ (\#2), Ward (\#59), GMM (\#81), FCM (\#151), BIRCH (\#64)], [$1 / n$], [$1 / n$], [Unbounded $O(r)$], [Single outlier pulls centroid/mean or explodes $Sigma_j$ at any $R -> infinity$],
    [*Tier 2*], [PAM (\#4), K-Medians (\#15), FasterPAM (\#24), CLARA (\#5), Medoidshift (\#42)], [$1 / n$ (hijack)], [$ceil(n_"min" / 2) / n approx 1 / (2 k)$], [Bounded $O(1)$], [Immune inside fixed cluster up to $50\%$; linear global hijack $R > O(n_"min" Delta)$],
    [*Tier 3*], [Trimmed $k$-Means (\#26), Robust CURE (\#65), Genie (\#69)], [$approx alpha$ ($< n_"min" / (2 n)$)], [$approx alpha$], [Hard-trimmed $0$ ($r > r_alpha$)], [Immune to up to $floor(alpha n)$ outliers; breaks at $floor(alpha n) + 1$ or bad $alpha$ choice],
    [*Tier 4*], [Student-$t$ Mix. (\#88), PCM (\#152), Gaussian KFCM-K (\#155), Power $k$-Means (\#27)], [Exp. / Bounded], [High (local basin)], [Redescending $-> 0$], [TMM: exp. hijack $exp(Omega(n_"min"))$ & $det Sigma -> 0$; PCM: coincident collapse],
    [*Tier 5*], [DBSCAN (\#31), HDBSCAN (\#32), OPTICS (\#33), DENCLUE (\#35), SNN (\#36), DenStream (\#194)], [$1 - "minPts" / n$ (dispersed)], [Up to $1 - k "minPts" / n$], [Zero (noise $C_(-1)$)], [Quarantines arbitrary dispersed noise; vulnerable to collusive bridge $Omega("minPts" D / epsilon)$],
  ),
)

// =============================================================================
// SECTION VI: INTERNAL & EXTERNAL VALIDATION GEOMETRY
// =============================================================================
#sec("VI", "Internal and External Validation Geometry: Convexity Bias")

Practitioners routinely compare algorithms across Categories I--VIII using classical internal validity indices—Mean Silhouette @rousseeuw1987silhouettes, Caliński--Harabasz (CHI) @calinski1974dendrite, Davies--Bouldin (DBI) @davies1979cluster, and Dunn @dunn1974well—assuming they act as geometry-neutral referees. We prove that *all four encode severe Euclidean-convexity biases*: on concentric non-convex manifolds, CHI, DBI, and Silhouette strictly prefer a geometrically *wrong* linear $k$-Means cut over *exact* manifold recovery.

#definition(title: "Classical Internal Validity Functionals")[
  For $Gamma = {C_1, dots, C_k}$ ($2 <= k <= n - 1$) on $X subset.eq bb(R)^d$ with centroids $mu_j$ and global mean $mu$:
  1. *(Caliński--Harabasz @calinski1974dendrite, higher better):* $"CH"(Gamma) := ("tr"(B_Gamma) / (k - 1)) / ("tr"(W_Gamma) / (n - k))$, where $"tr"(B_Gamma) = sum_j |C_j| norm(mu_j - mu)_2^2$ and $"tr"(W_Gamma) = sum_j sum_(x in C_j) norm(x - mu_j)_2^2$.
  2. *(Davies--Bouldin @davies1979cluster, lower better):* With $S_j := |C_j|^(-1) sum_(x in C_j) norm(x - mu_j)_2$,
     $ "DB"(Gamma) := 1 / k sum_(j=1)^k max_(ell != j) (S_j + S_ell) / norm(mu_j - mu_ell)_2. $
  3. *(Mean Silhouette @rousseeuw1987silhouettes, higher better):* For $x in C_j$, let $a(x) := (|C_j| - 1)^(-1) sum_(y in C_j without {x}) norm(x - y)_2$ and $b(x) := min_(ell != j) |C_ell|^(-1) sum_(y in C_ell) norm(x - y)_2$; then $s(x) := (b(x) - a(x)) / max(a(x), b(x))$ and $overline(s)(Gamma) := n^(-1) sum_i s(x_i)$.
  4. *(Dunn @dunn1974well, higher better):*
     $ "DI"(Gamma) := (min_(a != b) d(C_a, C_b)) / (max_(1 <= j <= k) "diam"(C_j)). $
]

#subsec("A", "Exact Convexity-Bias Theorem on Non-Convex Manifolds")

Consider two concentric circles $C_1 = S^1 (r_1), C_2 = S^1 (r_2) subset.eq bb(R)^2$ ($0 < r_1 < r_2$, ratio $rho := r_1 / r_2 in (0, 1)$, continuum limit of $n / 2$ uniform samples per ring). Compare exact recovery $Gamma_"true" = {C_1, C_2}$ (DBSCAN \#31, HDBSCAN \#32, Spectral \#103) against the false linear cut $Gamma_"lin" = {{x_1 >= 0}, {x_1 < 0}}$ ($k$-Means \#1).

#theorem(title: "Spherical-Convexity Bias of CHI, DBI, Silhouette, and Dunn", origin: "Original Proof")[
  On $C_1 = S^1 (r_1), C_2 = S^1 (r_2)$ with ratio $rho = r_1 / r_2 in (0, 1)$:
  1. *(CHI Collapse):* For every $rho in (0, 1)$, $"CH"(Gamma_"true") = 0$, whereas
     $ "CH"(Gamma_"lin") &= (n - 2) dot ((1 + rho)^2 / pi^2) / ((1 + rho^2) / 2 - (1 + rho)^2 / pi^2) \
       &= Theta(n) -> +infinity. $
  2. *(DBI Singularity):* For every $rho in (0, 1)$, $"DB"(Gamma_"true") = +infinity$, whereas
     $ "DB"(Gamma_"lin") &<= pi / (1 + rho) sqrt((1 + rho^2) / 2 - (1 + rho)^2 / pi^2) \
       &< pi / sqrt(2) approx 2.221. $
  3. *(Universal Negative Outer-Ring Silhouette & Inversion):* Letting $f(rho) := (2 pi)^(-1) integral_0^(2 pi) sqrt(1 + rho^2 - 2 rho cos theta) dif theta in [1, 4 / pi)$, for *every* $rho in (0, 1)$, $100\%$ of outer-ring points $y in C_2$ satisfy
     $ s_"true" (y) = pi / 4 f(rho) - 1 < 0, $
     approaching $pi / 4 - 1 approx -0.2146$ as $rho -> 0^+$, and the global true Silhouette satisfies $overline(s)(Gamma_"true") <= (pi / 8)(1 - rho)$ ($pi / 8 approx 0.3927$). As $rho -> 1^-$, the false cut attains
     $ lim_(rho -> 1^-) overline(s)(Gamma_"lin") &= 2 - (4 sqrt(2)) / pi ln(1 + sqrt(2)) approx 0.4130 \
       &>> 0 = lim_(rho -> 1^-) overline(s)(Gamma_"true") $
     (and at $rho = 1 / 2$, $overline(s)(Gamma_"lin") approx 0.3537 > 3 dot overline(s)(Gamma_"true") approx 0.1166$).
  4. *(Dunn Diameter Penalty):* $"DI"(Gamma_"true") = (1 - rho) / 2 -> 0$ on concentric circles and $h / (2 L) -> 0$ on parallel clean segments $[-L, L] times {0, h}$.
]

#proof[
  Scale by $r_2 = 1$ ($r_1 = rho in (0, 1)$).
  *Part 1 (CHI):* Under $Gamma_"true"$, $mu_1 = mu_2 = mu = bold(0)$, so $"tr"(B_(Gamma_"true")) = 0$ and $"CH"(Gamma_"true") = 0$. Under $Gamma_"lin"$, right/left semicircle centroids sit at $mu_(plus.minus) = (plus.minus (1 + rho) / pi, 0)$, giving $"tr"(B_(Gamma_"lin")) = n (1 + rho)^2 / pi^2$ and
  $ "tr"(W_(Gamma_"lin")) = n [ (1 + rho^2) / 2 - (1 + rho)^2 / pi^2 ] > 0, $
  yielding $"CH"(Gamma_"lin") = Theta(n)$.

  *Part 2 (DBI):* Under $Gamma_"true"$, $norm(mu_1 - mu_2)_2 = 0$ and $S_1 + S_2 = rho + 1 > 0$, giving $"DB"(Gamma_"true") = +infinity$. Under $Gamma_"lin"$, $norm(mu_+ - mu_-)_2 = 2 (1 + rho) / pi$ and Cauchy--Schwarz gives $S_+ = S_- <= sqrt(2 "tr"(W_(Gamma_"lin")) / n)$, bounding $"DB"(Gamma_"lin") < pi / sqrt(2) approx 2.221$.

  *Part 3 (Silhouette):* On $S^1 (r)$, mean chord distance is
  $ I(r) := 1 / (2 pi) integral_0^(2 pi) 2 r sin(theta / 2) dif theta = (4 r) / pi. $
  Differentiating $f(rho)$ twice on $[0, 1)$ gives
  $ f'(rho) &= 1 / (2 pi) integral_0^(2 pi) (rho - cos theta) (1 + rho^2 - 2 rho cos theta)^(-1 / 2) dif theta, \
    f''(rho) &= 1 / (2 pi) integral_0^(2 pi) (sin^2 theta) (1 + rho^2 - 2 rho cos theta)^(-3 / 2) dif theta > 0. $
  Since $f'(0) = 0$ and $f'' > 0$, $f$ is strictly increasing and strictly convex from $f(0) = 1$ to $f(1) = 4 / pi$, so $f(rho) <= 1 + (4 / pi - 1) rho$.
  - On $C_2$: $a_"true" (y) = 4 / pi > f(rho) = b_"true" (y)$ for all $rho in (0, 1)$, proving $s_"true" (y) = (pi / 4) f(rho) - 1 < 0$.
  - On $C_1$: $a_"true" (x) = (4 / pi) rho$ and $b_"true" (x) = f(rho)$. Because
    $ (1 + rho^2 - 2 rho cos theta) - rho^2 (2 - 2 cos theta) \
      quad = (1 - rho)[(1 - rho) + 2 rho (1 - cos theta)] >= (1 - rho)^2 > 0, $
    taking square roots and integrating proves $f(rho) > (4 / pi) rho$ everywhere on $(0, 1)$, so $s_"true" (x) = 1 - (4 rho) / (pi f(rho))$.
  - Averaging both circles cancels $1$:
    $ overline(s)(Gamma_"true") &= 1 / 2 [ pi / 4 f(rho) - (4 rho) / (pi f(rho)) ] \
      &<= 1 / 2 [ pi / 4 + (1 - pi / 4) rho - rho ] = (pi / 8)(1 - rho). $

  For $Gamma_"lin"$ as $rho -> 1^-$, at angle $phi in [-pi / 2, pi / 2]$ on $S^1 (1)$, splitting $|sin((theta - phi) / 2)|$ at $theta = phi$ yields
  $ a_"lin" (phi) &= 4 / pi [ 2 - cos((pi / 2 + phi) / 2) - cos((pi / 2 - phi) / 2) ] \
    &= 4 / pi [ 2 - sqrt(2) cos(phi / 2) ]. $
  Since $(a_"lin" + b_"lin") / 2 = I(1) = 4 / pi$, we get $b_"lin" (phi) = (4 sqrt(2) / pi) cos(phi / 2) > a_"lin" (phi)$ on $(-pi / 2, pi / 2)$ and $s_"lin" (phi) = 2 - sqrt(2) "sec"(phi / 2)$. Substituting $u = phi / 2$ yields
  $ lim_(rho -> 1^-) overline(s)(Gamma_"lin") &= 2 - (4 sqrt(2)) / pi integral_0^(pi / 4) "sec"(u) dif u \
    &= 2 - (4 sqrt(2)) / pi ln(1 + sqrt(2)) approx 0.4130, $
  strictly exceeding $lim_(rho -> 1^-) overline(s)(Gamma_"true") = 0$.

  *Part 4:* Concentric gap $r_2 - r_1$ over diameter $2 r_2$ gives $(1 - rho) / 2 -> 0$; parallel segments give $h / (2 L) -> 0$.
]

#corollary(title: "Failure of CHI and DBI Even on Linearly Separable Convex Bars", origin: "Original Proof")[
  Even when two clusters are *strictly convex and linearly separable*—two parallel uniform bars $C_1 = [-L, L] times {0}$ and $C_2 = [-L, L] times {h}$ ($L, h > 0$)—CHI strictly prefers the *wrong* transverse cut $Gamma_"vert" = {{x_1 >= 0}, {x_1 < 0}}$ over $Gamma_"true" = {C_1, C_2}$ for *every* aspect ratio $L / h > 1$, and DBI prefers $Gamma_"vert"$ for every $L / h > sqrt((1 + sqrt(37)) / 6) approx 1.087$:
  $ "CH"(Gamma_"vert") / "CH"(Gamma_"true") &= (4 L^4) / (h^2 (L^2 + 3 h^2)) > 1 quad (L > h), \
    "DB"(Gamma_"vert") &<= sqrt(1 / 3 + h^2 / L^2) < L / h = "DB"(Gamma_"true"). $
  Thus CHI and DBI mathematically enforce *isotropic sphericity*, not merely linear separability or convexity.
]

#proof[
  Under $Gamma_"true"$, centroids $(0, 0)$ and $(0, h)$ give $"tr"(B) = n h^2 / 4$, $"tr"(W) = n L^2 / 3$, and $S_1 = S_2 = L / 2$, yielding $"CH"(Gamma_"true") = (n - 2)(3 h^2) / (4 L^2)$ and $"DB"(Gamma_"true") = L / h$. Under $Gamma_"vert"$, centroids $(plus.minus L / 2, h / 2)$ give $"tr"(B) = n L^2 / 4$ and $"tr"(W) = n (L^2 / 12 + h^2 / 4)$, yielding $"CH"(Gamma_"vert") = (n - 2)(3 L^2) / (L^2 + 3 h^2)$ and, via Cauchy--Schwarz, $"DB"(Gamma_"vert") <= 2 sqrt(L^2 / 12 + h^2 / 4) / L = sqrt(1 / 3 + h^2 / L^2)$. For CHI, $4 L^4 - h^2 (L^2 + 3 h^2) = (L^2 - h^2)(4 L^2 + 3 h^2) > 0 <=> L > h$. For DBI, writing $r := L / h$, $sqrt(1 / 3 + r^(-2)) < r <=> 3 r^4 - r^2 - 3 > 0 <=> r > sqrt((1 + sqrt(37)) / 6) approx 1.0865$.
]

#wide_figure(
  caption: [*Table IV:* Geometric Bias & Invariance Audit of Classical Internal and External Validation Metrics.],
  table(
    columns: (1.05fr, 0.80fr, 1.00fr, 0.82fr, 1.25fr),
    align: (left, center, left, center, left),
    [Validation Index], [Type], [Concentric Circles Verdict], [Convexity Bias], [Recommended Usage / Failure Mode],
    [Silhouette @rousseeuw1987silhouettes], [Internal], [Outer ring $< 0$; favors cut $3 times$], [*Strong* (chord)], [Spherical/compact clusters; fails rings & rings-in-rings],
    [Caliński--Harabasz @calinski1974dendrite], [Internal], [$"CH"(Gamma_"true") = 0$ vs. $Theta(n)$], [*Severe* ($L > h$)], [Isotropic $k$-Means only; fails elongated convex bars too],
    [Davies--Bouldin @davies1979cluster], [Internal], [$"DB"(Gamma_"true") = +infinity$], [*Severe* ($L > 1.09 h$)], [Compact separated balls; explodes on shared centroids],
    [Dunn Index @dunn1974well], [Internal], [$(1 - rho) / 2 -> 0$], [Diameter-biased], [Satisfies CQM axioms (Thm. 2) but penalizes filaments],
    [DBCV / MST-Validity], [Internal], [Rewards dense rings], [*None* (geodesic)], [Density/manifold algorithms (Cat. II & Cat. V)],
    [ARI @hubert1985comparing], [External], [Exact $1.000$], [*None* (exact)], [Chance-corrected benchmark; quadratic minority weight],
    [NMI @strehl2002cluster], [External], [Exact $1.000$], [*None* (exact)], [Entropy-balanced; pair with ARI on imbalanced classes],
  ),
)

Theorem 8 explains why maximizing Silhouette, CHI, or DBI systematically selects spherical $k$-Means over density (Cat. II) and spectral (Cat. V) methods on non-convex manifolds (confirmed empirically in Chapter XIII). Manifold validation instead requires density-geodesic indices such as DBCV built on mutual reachability MSTs.

#subsec("B", "External Validation Geometry: Hypergeometric Chance Correction of ARI")

Given ground truth $U = {U_1, dots, U_r}$ and prediction $V = {V_1, dots, V_c}$ with contingency counts $n_(i j) := |U_i inter V_j|$ and marginals $a_i, b_j$, raw Rand Index $"RI"(U, V)$ drifts toward $1$ as $k$ grows. Hubert and Arabie @hubert1985comparing cured this bias via the *Adjusted Rand Index (ARI)*.

#proposition(title: "Exact Hypergeometric Expectation & Homeomorphism Invariance of ARI", origin: "Classical — Cited")[
  Let $P := sum_i binom(a_i, 2)$, $Q := sum_j binom(b_j, 2)$, $T := sum_(i, j) binom(n_(i j), 2)$, and $N := binom(n, 2)$ @hubert1985comparing.
  1. Under the fixed-marginal random permutation model ($V$ uniform with block sizes $b_j$, independent of $U$), $bb(E)[T] = (P Q) / N$ and $bb(E)["RI"] = 1 - (P + Q) / N + (2 P Q) / N^2$.
  2. Defining
     $ "ARI"(U, V) &:= ("RI" - bb(E)["RI"]) / (1 - bb(E)["RI"]) \
       &= (T - (P Q) / N) / (1 / 2 (P + Q) - (P Q) / N), $
     yields $bb(E)["ARI"] = 0$, $"ARI" <= 1$ with equality iff $U = V$ up to relabeling, and complete invariance under arbitrary ambient homeomorphisms (zero convexity bias).
  3. *(Majority-Cluster Dominance):* With $k - 1$ clusters of size $m$ and one minority cluster of size $m_"min" << m$, absorbing the minority cluster changes $"ARI"$ by only $O(m_"min" / (k m))$, and minority pairs carry quadratic weight $O(m_"min"^2 / (k m^2))$.
]

#proof[
  *Part 1:* Letting $I_e^U := bb(1)[x tilde.op_U y]$ and $I_e^V := bb(1)[x tilde.op_V y]$ over all $N$ unordered pairs $e = {x, y}$, transitivity of random permutations yields $bb(P)(I_e^V = 1) = Q / N$, so $bb(E)[T] = sum_e I_e^U (Q / N) = (P Q) / N$. Since $"RI" = 1 - (P + Q - 2 T) / N$, taking expectation gives $bb(E)["RI"]$.
  *Part 2:* Both $"RI" - bb(E)["RI"] = (2 / N)(T - (P Q) / N)$ and $1 - bb(E)["RI"] = (2 / N)((P + Q) / 2 - (P Q) / N)$ share factor $2 / N$, giving the stated ratio with $bb(E)["ARI"] = 0$. Since $T <= min(P, Q) <= (P + Q) / 2$ with equality iff $T = P = Q$ ($U = V$), $"ARI" <= 1$ characterizes exact recovery and depends solely on $[n_(i j)]$, independent of ambient geometry.
  *Part 3:* Merging a minority block into one size-$m$ cluster increases $Q$ by $m m_"min"$ while $P = Theta(k m^2)$, depressing $"ARI"$ by $O(m_"min" / (k m))$ while $binom(m_"min", 2) / P = O(m_"min"^2 / (k m^2))$.
]

Consequently, imbalanced benchmarks should always pair $"ARI"$ with Normalized Mutual Information ($"NMI"$) or macro-averaged cluster F1.

// =============================================================================
// SECTION VII: UNIFIED META-SYNTHESIS ACROSS ALL 200 ALGORITHMS
// =============================================================================
#sec("VII", "Unified Meta-Synthesis: Six Laws Across All 200 Algorithms")

The six meta-theorems of Sections I--VI interlock as a rigid conservation law governing the entire 200-algorithm design space.

#proposition(title: "No-Free-Lunch Conservation Principle of Unsupervised Clustering", origin: "Original Unified Proof")[
  No single flat metric clustering operator without an explicit inductive prior can simultaneously achieve: (1) Scale-Invariance, Richness, and Consistency (Theorem 1); (2) non-vanishing contrast under isotropic i.i.d. coordinates as $d -> infinity$ (Theorems 4--5); (3) assumption-free $k$-selection via empirical stability (Theorem 6, Proposition 4); (4) exact Gaussian/Bregman mean optimality alongside Donoho--Huber breakdown $> 1 / n$ (@banerjee2005clustering, Theorem 7); and (5) exact non-convex manifold recovery maximizing Euclidean Silhouette, CHI, or DBI (Theorem 8).
]

#proof[
  Each obstruction is immediate: (1) Lemma 1 forbids chains in $"Range"(f)$; (2) Theorem 4 drives relative $L_p$ contrast to $0$, forcing subspace, spherical, or rank priors (Theorem 5); (3) Proposition 4 constructs $"instab" -> 0$ at $k != k^*$ and $"instab" -> 1 / 4$ at $k = k^*$; (4) Bregman mean uniqueness @banerjee2005clustering fixes arithmetic means as unique quadratic minimizers while Theorem 7(1) pins their breakdown to $1 / n$; and (5) Theorem 8 yields $"CH" = 0$, $"DB" = +infinity$, and negative outer-ring Silhouette on exact concentric manifolds.
]

Table V synthesizes all six meta-laws across Categories I--VIII, bridging directly to the 200-algorithm master matrices (Chapter XII), empirical multi-geometry validation (Chapter XIII), and practitioner selection architecture (Chapter XIV).

#wide_figure(
  caption: [*Table V:* Grand Unified Governance Matrix Mapping All Six Meta-Theorems Across the Eight Clustering Families (Algorithms 1--200).],
  table(
    columns: (0.98fr, 0.88fr, 0.88fr, 0.90fr, 0.82fr, 0.84fr, 1.00fr),
    align: (left, left, left, left, left, left, left),
    [Category (Algs.)], [Kleinberg Axiom Trade-off], [Functor / CQM Regime], [High-$d$ Behavior ($d -> infinity$)], [Stability & $k$-Selection], [Breakdown Tier ($epsilon_n^*$)], [Validation Geometry Alignment],
    [*I.* Centroid / Medoid (1--30)], [S-Inv. only; sacrifices Rich. & Cons. (Prop. 2)], [Fixed-$k$ rel. CQM; non-functorial], [Euclidean collapses; Sph. \#12 thrives on $bb(S)^(d - 1)$], [Prone to Prop. 4 false stability; use BIC (\#10)], [Tier 1 ($1 / n$) to Tier 3 ($approx alpha$, \#26)], [Maximizes CHI / Sil. / DBI (convex bias)],
    [*II.* Density / Mode (31--55)], [DBSCAN: Rich.+Cons., breaks S-Inv.; HDBSCAN: S-Inv.], [HDBSCAN lifts $frak(H)_"SL"$ to $d_"mreach"$], [Fixed-$epsilon$ collapses; SNN (\#36, \#52) survives], [EOM persistence (\#32) avoids resampling trap], [Tier 5: $1 - "minPts" / n$ disp.; Vitali bridge], [Penalized by CHI/DBI/Sil.; use DBCV / ARI],
    [*III.* Hierarchical (56--80)], [Fixed-$k$ SL: S-Inv.+Cons.; CL/UPGMA/Ward: S-Inv. only], [*SL \#56 is unique C--M functor* (Thm. 3)], [Global metric collapses; MST chains degrade], [Multi-scale tree defers flat $k$ cut], [Tier 1 (Ward) to Tier 3 (CURE \#65, Genie \#69)], [Dendrogram cophenetic & Dasgupta cost @dasgupta2016cost],
    [*IV.* Mixture / Topic (81--102)], [Full-$Sigma$ GMM: S-Inv.; DPGMM/BIC break S-Inv. & Cons.], [Likelihood ratio CQM; non-functorial], [Full $Sigma$ $O(d^2)$ overfits; vMF \#89 & LDA survive], [BIC/ICL (\#93), DP prior (\#82, \#99), MML (\#102)], [Tier 1 (GMM) to Tier 4 (TMM \#88, exp. hijack)], [Ellipsoidal BIC / log-lik.; penalized on rings],
    [*V.* Graph / Spectral (103--132)], [$k$-NN graph S-Inv.; modularity resolution limit], [Graph cut / Ncut CQM; non-functorial], [$k$-NN / SNN graph survives; fixed-$sigma$ RBF fails], [Eigengap (\#103), SBM likelihood (\#130--131)], [Robust to isolated nodes; bridge-sensitive], [Recovers manifolds; penalized by Euclidean Sil./CHI],
    [*VI.* Grid / Subspace (133--150)], [Fixed grid breaks S-Inv.; SSC ray S-Inv., non-Cons.], [Subspace density / $ell_1$ CQM], [*Built to beat Thm. 4*: $sqrt(d / s)$ SNR gain], [MDL / entropy (\#143), eigengap (\#150)], [Tier 3--5 (noise cells / $ell_1$ sparse outliers)], [Requires subspace-projected validation / ARI],
    [*VII.* Fuzzy / Domain (151--175)], [FCM/GK S-Inv., non-Cons.; PCM breaks S-Inv.], [Fuzzy / hypergraph CQM; non-functorial], [Domain-specific (DTW, bicluster, multi-view)], [Monti \#170 hits Thm. 6 symmetry trap!], [Tier 1 (FCM) to Tier 4 (PCM \#152, Gauss KFCM)], [Xie--Beni / fuzzy Sil. (still convex-biased)],
    [*VIII.* Neural / Stream (176--200)], [Encoder/vigilance breaks S-Inv. & Cons.], [Variational / OT / contrastive CQM], [Sph. contrastive (\#190, \#192) uses Thm. 5(3) blessing], [Self-supervised OT / Sinkhorn (\#190) prevents collapse], [Tier 1 (DEC) to Tier 5 (DenStream \#194)], [Latent-space Sil. valid only if encoder unfolds manifold],
  ),
)

The six meta-theorems proved in this chapter complete the theoretical foundation of Part IX. Rather than pursuing an impossible universal clustering algorithm—ruled out jointly by Theorems 1, 4, 6, 7, and 8—rigorous unsupervised modeling matches an algorithm's *inductive invariants* to the observable invariants of the data.

#subsec("A", "Operational Design Laws Governing Chapters XII–XIV")

Each negative boundary in Sections I--VI directly prescribes a concrete architectural rule for navigating the 200 algorithms cataloged in this monograph:
1. *Axiomatic Choice Over Impossible Ideals (Theorems 1--3):* Commit explicitly to which Kleinberg axiom to release. When multi-scale structural consistency under perturbations is paramount, prefer the unique Gromov--Hausdorff 1-Lipschitz functor Single Linkage (\#56) or its density-regularized lifting HDBSCAN (\#32) over fixed-$k$ centroid partitioners; when partition ranking across variable $k$ is required, evaluate scale-invariant Clustering Quality Measures rather than unnormalized distortion.
2. *Dimensionality Regime Matching (Theorems 4--5):* Never deploy ambient Euclidean $k$-Means (\#1) or fixed-$epsilon$ DBSCAN (\#31) blindly when $d >> s$ without testing relative contrast $(D_"max" - D_"min") / D_"min"$. When relevant structure lies on low-dimensional subspaces, deploy Category VI (CLIQUE \#134, PROCLUS \#136, SUBCLU \#137, SSC \#150); when directional alignment persists on $bb(S)^(d - 1)$ ($gamma >> d^(-1 / 4)$), exploit the spherical concentration blessing via Spherical $k$-Means (\#12), vMF Mixtures (\#89), or hyperspherical deep models (SwAV \#190, CC \#192); when raw metrics concentrate yet local neighborhood purity remains positive, lift affinities to Shared Nearest Neighbor ranks (SNN \#36, SNN-DPC \#52).
3. *Model-Order Selection Beyond Resampling Stability (Theorem 6, Proposition 4):* Treat resampling consensus (Monti \#170) as a diagnostic of *objective uniqueness and basin curvature*, never as a standalone oracle for $k^*$. Pair stability checks with explicit complexity penalties (BIC in X-Means \#10 and MCLUST \#93; MML in SNOB \#102), nonparametric Bayesian priors (DPGMM \#82, IGMM \#99), or excess-of-mass topological persistence (HDBSCAN \#32).
4. *Robustness Tier Alignment (Theorem 7, Corollary 2):* Whenever data contains heavy tails, measurement spikes, or adversarial contamination, depart explicitly from Tier 1 Bregman arithmetic means ($epsilon_n^* = 1 / n$): select Tier 2 $L_1$ Medians/PAM (\#4, \#15, \#24) for bounded local contamination below $1 / (2 k)$, Tier 3 Trimmed $k$-Means (\#26) when an explicit contamination cap $alpha < n_"min" / (2 n)$ is known, Tier 4 Student-$t$ Mixtures (\#88) or calibrated PCM (\#152) for redescending heavy-tail immunity, and Tier 5 DBSCAN/HDBSCAN (\#31, \#32) when dispersed background noise must be isolated into $C_(-1)$ while guarding inter-cluster gaps against collusive bridging.
5. *Geometry-Faithful Validation (Theorem 8, Corollary 3, Proposition 5):* Never adjudicate non-convex or elongated manifold recoveries using Euclidean Silhouette, CHI, or DBI, which strictly enforce isotropic sphericity and invert rankings by over $3 times$ on concentric rings and even linearly separable bars ($L > 1.09 h$). Evaluate non-convex models internally with geodesic/density-MST criteria (DBCV) and externally with chance-corrected ARI paired with entropy-balanced NMI on imbalanced classes.

Comparison Chapter XII codifies these theoretical invariants exhaustively across all 200 individual algorithms in master reference matrices; Chapter XIII validates our exact analytical predictions—including high-dimensional contrast decay and Theorem 8's Silhouette inversion on concentric manifolds—empirically on the ClusterOverview benchmark engine; and Chapter XIV distills all six meta-laws into an actionable engineering selection architecture.
