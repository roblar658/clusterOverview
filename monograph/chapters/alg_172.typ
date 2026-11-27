#import "../template/ieee.typ": *

#alg_entry(
  id: 172,
  name: "HGPA (HyperGraph Partitioning Algorithm)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2002",
  paradigm: "Direct k-Way Hypergraph Minimum-Cut Partitioning",
  complexity: $O(n k r)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Given $n$ objects indexed by $V = {1, dots, n}$ and a cluster ensemble $Lambda = {lambda^((1)), dots, lambda^((r))}$ of $r$ base partitions with cluster counts $k^((1)), dots, k^((r))$, Strehl and Ghosh @strehl2002cluster observe that constructing the dense $n times n$ co-association matrix of CSPA (Alg. 171) incurs quadratic $O(n^2)$ memory and expands every cluster into a clique of pairwise edges. The _HyperGraph Partitioning Algorithm (HGPA)_ avoids pairwise clique expansion altogether by representing the entire ensemble directly as a single hypergraph and formulating consensus clustering as a direct $k$-way minimum hyperedge-cut problem.

Let $H = [H^((1)), dots, H^((r))] in {0, 1}^(n times m)$ be the concatenated cluster-membership incidence matrix with $m = sum_(q=1)^r k^((q))$ columns. Each column $h_a in {0, 1}^n$ defines a _hyperedge_ $e_a = {i in V : H_(i a) = 1}$ corresponding to one cluster from one base clustering. This induces the hypergraph $cal(H) = (V, cal(E), w)$ with vertex set $V$, hyperedge multiset $cal(E) = {e_1, dots, e_m}$, and uniform hyperedge weights $w(e_a) = 1$ (generalizable to $w(e_a) > 0$ for weighted ensembles).

#definition(title: "Hyperedge Span, Unit Cut, and Balanced HGPA Objective")[
  For any $k$-way partition $cal(C) = {C_1, dots, C_k}$ of $V$ and any hyperedge $e in cal(E)$, define its _span_ (or connectivity) across $cal(C)$ by
  $ lambda(e; cal(C)) = abs({c in {1, dots, k} : e inter C_c != emptyset}). $
  A hyperedge $e$ is _cut_ by $cal(C)$ iff $lambda(e; cal(C)) >= 2$. The _unit hyperedge cut_ and the _$(lambda - 1)$-connectivity cut_ are defined respectively by
  $ "hcut"(cal(C); cal(H)) = sum_(e in cal(E)) w(e) dot bb(1)(lambda(e; cal(C)) >= 2), $
  $ "hcut"_(lambda - 1) (cal(C); cal(H)) = sum_(e in cal(E)) w(e) (lambda(e; cal(C)) - 1). $
  HGPA seeks a consensus $k$-partition minimizing $"hcut"(cal(C); cal(H))$ subject to strict vertex balance:
  $ (1 - epsilon_"bal") n / k <= abs(C_c) <= (1 + epsilon_"bal") n / k, quad forall c in {1, dots, k}. $
]

Because cutting a hyperedge $e_a$ corresponds to splitting a base cluster across two or more consensus clusters, minimizing $"hcut"$ seeks a balanced partition that keeps as many base clusters wholly intact as possible. HGPA solves this combinatorial program via multilevel hypergraph partitioning (HMETIS, based on the multilevel coarsening and Fiduccia–Mattheyses / Kernighan–Lin refinement paradigm @karypis1998metis @kernighan1970efficient).

#algbox(
  title: "HGPA (HyperGraph Partitioning Algorithm)",
  inputs: [Cluster ensemble $Lambda = {lambda^((1)), dots, lambda^((r))}$ on $n$ objects, target cluster count $k$, balance tolerance $epsilon_"bal" >= 0$.],
  outputs: [Consensus $k$-partition $lambda^* : V -> {1, dots, k}$ and hyperedge-cut score $"hcut"$.],
)[
  + *Hypergraph construction*: initialize vertex set $V = {1, dots, n}$ with unit vertex weights $mu(i) = 1$; for every base clustering $q = 1, dots, r$ and non-empty cluster $C^((q))_c$, add hyperedge $e_a = C^((q))_c$ with weight $w(e_a) = 1$ (total pins $sum_a abs(e_a) <= n r$).
  + *Multilevel coarsening*: iteratively contract high-overlap vertex subsets via heavy-edge / hyperedge-coarsening matchings $phi_t : V_t -> V_(t+1)$, summing vertex weights and merging duplicate hyperedges while pruning singleton hyperedges ($abs(e') = 1$).
  + *Initial $k$-way partitioning*: on the coarsest hypergraph $cal(H)_L$ ($abs(V_L) = O(k)$), compute a balanced $k$-way partition minimizing $"hcut"$ via repeated randomized recursive bisection or greedy growing.
  + *Uncoarsening and FM refinement*: project the partition back through $cal(H)_(L-1), dots, cal(H)_0$; at each level, execute boundary Fiduccia–Mattheyses vertex-move passes @kernighan1970efficient that decrease $"hcut"$ while preserving balance $abs(C_c) <= (1 + epsilon_"bal") n / k$.
  + *Consensus output*: return the refined labeling $lambda^*$ on $V$ and evaluate $Phi^("ANMI") (Lambda, lambda^*)$ for supra-consensus selection @strehl2002cluster.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results revealing both HGPA's information-theoretic rationale and its acute mathematical vulnerability to boundary noise: an exact mutual-information lower bound, a constructive proof of single-vertex cut saturation, submodularity and clique-expansion distortion bounds, NP-hardness plus multilevel coarsening cut-preservation, and heterogeneous-resolution bias.

#subsec("A", "Information-Theoretic Lower Bound on Average Mutual Information")

Strehl and Ghosh @strehl2002cluster motivate HGPA heuristically as an approximate maximizer of Average Normalized Mutual Information. The following theorem makes this link mathematically precise by bounding the ensemble average mutual information directly from below in terms of $"hcut"$.

#theorem(title: "Hyperedge-Cut Lower Bound on Ensemble Mutual Information", origin: "Original Proof")[
  Let $cal(H)$ be the unweighted ensemble hypergraph ($w(e) equiv 1$) of a complete ensemble $Lambda$, and let $s_"max" = max_(e in cal(E)) abs(e)$ denote the maximum base-cluster cardinality. For any consensus $k$-partition $cal(C)$ with Shannon entropy $H(cal(C)) = -sum_(c=1)^k (abs(C_c) / n) log_2 (abs(C_c) / n)$, the average mutual information obeys
  $ 1 / r sum_(q=1)^r I(cal(C); lambda^((q))) >= H(cal(C)) - (s_"max" log_2 k) / (r n) "hcut"(cal(C); cal(H)), $
  and under the connectivity cut,
  $ 1 / r sum_(q=1)^r I(cal(C); lambda^((q))) >= H(cal(C)) - s_"max" / (r n) "hcut"_(lambda - 1) (cal(C); cal(H)). $
]
#proof[
  Fix a base clustering $lambda^((q))$ with clusters $e_a = C^((q))_a$ ($a = 1, dots, k^((q))$). Conditioning $cal(C)$ on the event ${lambda^((q)) = a}$ restricts the consensus labels to the $s_a = lambda(e_a; cal(C)) in {1, dots, k}$ blocks of $cal(C)$ that intersect $e_a$. Because any discrete distribution supported on at most $s_a$ atoms has binary entropy at most $log_2 s_a$,
  $ H(cal(C) | lambda^((q)) = a) <= log_2 lambda(e_a; cal(C)). $
  Observe two elementary inequalities for all integers $s in {1, dots, k}$: first, $log_2 s <= (log_2 k) dot bb(1)(s >= 2)$ (both sides vanish at $s = 1$); second, $log_2 s <= s - 1$ (since $2^(s-1) >= s$ for all $s >= 1$, with equality at $s in {1, 2}$). Weighting by the cluster probability $abs(e_a) / n <= s_"max" / n$ and summing over $a = 1, dots, k^((q))$ gives
  $ H(cal(C) | lambda^((q))) & = sum_(a=1)^(k^((q))) abs(e_a) / n H(cal(C) | lambda^((q)) = a) \
    & <= (s_"max" log_2 k) / n sum_(a=1)^(k^((q))) bb(1)(lambda(e_a; cal(C)) >= 2), $
  and likewise $H(cal(C) | lambda^((q))) <= (s_"max" / n) sum_a (lambda(e_a; cal(C)) - 1)$. Averaging over $q = 1, dots, r$ and substituting into $I(cal(C); lambda^((q))) = H(cal(C)) - H(cal(C) | lambda^((q)))$ yields both claimed inequalities directly.
]

Theorem 1 reveals the exact two-part logic of HGPA: enforcing strict cluster balance pushes $H(cal(C))$ toward its maximum $log_2 k$, while minimizing $"hcut"(cal(C); cal(H))$ minimizes the conditional-entropy upper bound; when $"hcut" = 0$, the bound is tight and $I(cal(C); lambda^((q))) = H(cal(C))$ for all $q$.

#subsec("B", "Single-Vertex Noise Saturation and Objective Collapse")

Despite Theorem 1, the unit hyperedge cut suffers from a severe structural fragility: cutting a single vertex out of a cluster of size $n / k$ incurs the exact same unit penalty ($1$) as bisecting the cluster cleanly in half. We prove constructively that $O(1 / n)$ label noise completely blinds HGPA to the true partition.

#theorem(title: "Zero-Gradient Cut Saturation Under Single-Point Base Noise", origin: "Original Counterexample")[
  Let $V$ have $n = 2 t$ vertices ($t >= 4$ even) partitioned into two balanced ground-truth clusters $C_1^*, C_2^*$ ($abs(C_1^*) = abs(C_2^*) = t$). Construct an ensemble $Lambda$ of $r >= 1$ base bipartitions where each $lambda^((q))$ agrees with $cal(C)^*$ on $n - 2$ vertices and swaps exactly one pair $(u_q, v_q) in C_1^* times C_2^*$ (error rate $2 / n = O(1 / n)$). Let $cal(C)_"rand"$ be an orthogonal balanced bipartition splitting each $C_c^*$ exactly in half ($"ARI"(cal(C)_"rand", cal(C)^*) = 0$). Then:
  (a) Under HGPA's unit hyperedge cut, $cal(C)^*$ and $cal(C)_"rand"$ are *completely indistinguishable*:
  $ "hcut"(cal(C)^*; cal(H)) = "hcut"(cal(C)_"rand"; cal(H)) = 2 r = m. $
  (b) Under CSPA's co-association cut (Alg. 171) on the identical ensemble, the ratio of cut costs diverges linearly with $n$:
  $ "Cut"(cal(C)_"rand"; S) / "Cut"(cal(C)^*; S) >= (t / 2 - 1)^2 / (t - 1) = Omega(n). $
]
#proof[
  *Part (a):* In replicate $q$, the two hyperedges are $e^((q))_1 = (C_1^* without {u_q}) union {v_q}$ and $e^((q))_2 = (C_2^* without {v_q}) union {u_q}$. Because $v_q in C_2^*$ and $C_1^* without {u_q} != emptyset$, $e^((q))_1$ intersects both $C_1^*$ and $C_2^*$, so $lambda(e^((q))_1; cal(C)^*) = 2$; identically $lambda(e^((q))_2; cal(C)^*) = 2$. Hence every single one of the $2 r$ hyperedges is cut by $cal(C)^*$, giving $"hcut"(cal(C)^*; cal(H)) = 2 r$. Now let $cal(C)_"rand" = {A_1 union B_1, A_2 union B_2}$ where $C_1^* = A_1 union A_2$ and $C_2^* = B_1 union B_2$ with $abs(A_1) = abs(A_2) = abs(B_1) = abs(B_2) = t / 2 >= 2$. Since removing one element $u_q$ from $C_1^*$ leaves at least $t / 2 - 1 >= 1$ vertices in each of $A_1$ and $A_2$, both $e^((q))_1$ and $e^((q))_2$ intersect both blocks of $cal(C)_"rand"$, giving $"hcut"(cal(C)_"rand"; cal(H)) = 2 r$.

  *Part (b):* For CSPA, in run $q$ the only cross-pairs between $C_1^*$ and $C_2^*$ that co-cluster are ${u_q} times (C_2^* without {v_q})$ and ${v_q} times (C_1^* without {u_q})$, totaling $2 (t - 1)$ pairs per run. Averaging over $r$ yields $"Cut"(cal(C)^*; S) = 2 (t - 1)$. Under $cal(C)_"rand"$, in run $q$ at least $t / 2 - 1$ elements of $A_1$ co-cluster with at least $t / 2 - 1$ elements of $A_2$, and likewise for $B_1$ and $B_2$, contributing at least $2 (t / 2 - 1)^2$ cut weight per run. Thus $"Cut"(cal(C)_"rand"; S) >= 2 (t / 2 - 1)^2$, and the cut ratio satisfies
  $ "Cut"(cal(C)_"rand"; S) / "Cut"(cal(C)^*; S) >= (t / 2 - 1)^2 / (t - 1) = Omega(t) = Omega(n). $
]

Theorem 2 provides a rigorous explanation for HGPA's well-known empirical fragility: whenever every base partition contains even minimal boundary scatter, the hypergraph cut landscape collapses into a flat plateau at $"hcut" = m$, depriving HMETIS of gradient information to distinguish coherent clusters from noise.

#subsec("C", "Submodularity and Clique-Expansion Distortion Bounds")

#lemma(title: "Submodularity of Hypergraph Cut and Clique Distortion Bounds", origin: "Original Derivation")[
  Let $cal(H) = (V, cal(E), w)$ with $w(e) >= 0$ and hyperedge cardinalities $2 <= s_"min" <= abs(e) <= s_"max"$.
  (a) The cut function $f(A) = "hcut"({A, V without A}; cal(H))$ is symmetric ($f(A) = f(V without A)$) and submodular: for all $A, B subset.eq V$,
  $ f(A) + f(B) >= f(A union B) + f(A inter B). $
  (b) Let $G_"clq"$ be the weighted clique-expansion graph assigning weight $w(e) / abs(e)$ to every unordered pair ${u, v} subset.eq e$ (summed over $e in cal(E)$). Then for every cut $A subset.eq V$,
  $ (1 - 1 / s_"min") f(A) <= "cut"(A, V without A; G_"clq") <= s_"max" / 4 f(A). $
]
#proof[
  *Part (a):* Symmetry is immediate from the definition of cut hyperedges. Since $f(A) = sum_(e in cal(E)) w(e) f_e (A)$ with $w(e) >= 0$, it suffices to verify submodularity for a single hyperedge indicator $f_e (A) = bb(1)(e inter A != emptyset and e inter A^c != emptyset) in {0, 1}$.
  - If $f_e (A union B) + f_e (A inter B) = 0$, the inequality is trivial.
  - If $f_e (A union B) = 1$, there exist $x in e inter (A union B)$ and $y in e without (A union B)$. Then $y in.not A$ and $y in.not B$, while $x$ belongs to at least one of $A$ or $B$, forcing $f_e (A) + f_e (B) >= 1$.
  - If $f_e (A inter B) = 1$, there exist $x in e inter (A inter B)$ (so $x in A$ and $x in B$) and $y in e without (A inter B)$ (so $y$ misses at least one of $A$ or $B$), forcing $f_e (A) + f_e (B) >= 1$.
  - If both $f_e (A union B) = 1$ and $f_e (A inter B) = 1$, then $e$ contains both some $x in A inter B$ and some $y in.not A union B$; hence $x in A, y in.not A => f_e (A) = 1$ and $x in B, y in.not B => f_e (B) = 1$, giving $2 >= 2$.

  *Part (b):* Fix a hyperedge $e$ of size $s = abs(e) >= 2$, and let $a = abs(e inter A)$. If $a in {0, s}$, both $f_e (A)$ and the clique contribution vanish. If $1 <= a <= s - 1$, then $f_e (A) = 1$, while $e$ contributes $(w(e) / s) a (s - a)$ to $"cut"(A, V without A; G_"clq")$. Over integers $1 <= a <= s - 1$, the quadratic $a(s - a)$ is minimized at $a in {1, s - 1}$ with value $s - 1$ and maximized at $a = s / 2$ with value $s^2 / 4$. Dividing by $s$ yields $(1 - 1 / s) w(e) <= (w(e) / s) a (s - a) <= (s / 4) w(e)$. Summing over all $e in cal(E)$ establishes the two-sided bound.
]

Because $s_"max" = Theta(n / k)$ in balanced ensembles, Lemma 3 proves that the multiplicative distortion between hypergraph cut and clique-expansion cut (as used by CSPA and spectral hypergraph relaxations @chung1997spectral @shi2000ncut) grows as $Theta(n / k)$, formalizing why HGPA and CSPA behave as fundamentally distinct objectives.

#subsec("D", "NP-Hardness and Exact Multilevel Coarsening Invariance")

#theorem(title: "NP-Hardness and Exact Coarsening Cut-Preservation Invariant", origin: "Classical — Cited")[
  (a) Minimizing $"hcut"(cal(C); cal(H))$ subject to $abs(C_c) = n / k$ is NP-hard even for $k = 2$ and uniform hyperedge size $abs(e) = 2$ @karypis1998metis @strehl2002cluster.
  (b) Let $phi : V -> V'$ be any surjective coarsening contraction mapping $cal(H)$ to the coarsened hypergraph $cal(H)' = (V', cal(E)', w)$ with supernode weights $mu(v') = abs(phi^(-1) (v'))$ and projected hyperedges $e' = phi(e)$ (omitting collapsed singletons $abs(phi(e)) = 1$). For *every* $k$-partition $cal(C)' = {C'_1, dots, C'_k}$ of $V'$, its pullback $cal(C) = {phi^(-1) (C'_1), dots, phi^(-1) (C'_k)}$ on $V$ preserves both exact balance and exact hyperedge cuts:
  $ abs(phi^(-1) (C'_c)) = sum_(v' in C'_c) mu(v'), quad "hcut"(cal(C); cal(H)) = "hcut"(cal(C)'; cal(H)'), $
  and identically $"hcut"_(lambda - 1) (cal(C); cal(H)) = "hcut"_(lambda - 1) (cal(C)'; cal(H)')$.
]
#proof[
  *Part (a):* Restricting to $2$-uniform hypergraphs ($abs(e) = 2$ for all $e$) and $k = 2$ with $epsilon_"bal" = 0$ reduces HGPA identically to the classical NP-hard Minimum Graph Bisection problem @kernighan1970efficient @karypis1998metis.

  *Part (b):* Because ${phi^(-1) (v') : v' in V'}$ partitions $V$, we immediately obtain $abs(phi^(-1) (C'_c)) = sum_(v' in C'_c) abs(phi^(-1) (v')) = sum_(v' in C'_c) mu(v')$. For any hyperedge $e in cal(E)$ and any block $C'_c subset.eq V'$,
  $ e inter phi^(-1) (C'_c) != emptyset <=> exists x in e : phi(x) in C'_c <=> phi(e) inter C'_c != emptyset. $
  Consequently, the span is invariant for every hyperedge: $lambda(e; cal(C)) = lambda(phi(e); cal(C)')$. In particular, any collapsed hyperedge with $abs(phi(e)) = 1$ lies inside a single supernode and thus satisfies $lambda(e; cal(C)) = 1$, contributing $0$ to both $"hcut"$ and $"hcut"_(lambda - 1)$. Summing over all $e in cal(E)$ establishes both cut-preservation identities. Because every FM uncoarsening move accepts only non-increasing changes in $"hcut"$, multilevel projection guarantees monotonic cut non-increase from the coarsest level down to $cal(H)$.
]

#subsec("E", "Resolution Bias Under Heterogeneous Base Cluster Counts")

#proposition(title: "Hyperedge-Count Domination by Fine-Grained Base Partitions", origin: "Original Proof")[
  Under uniform hyperedge weights $w(e) equiv 1$, if an ensemble combines $r_1$ coarse partitions with $k_1$ clusters and $r_2$ fine partitions with $k_2 > k_1$ clusters, the maximum total cut weight assignable to the fine group exceeds that of the coarse group by the factor $(r_2 k_2) / (r_1 k_1)$. Even for $r_1 = r_2 = 1$, splitting a single micro-cluster of size $2$ incurs the same unit penalty as splitting a macro-cluster of size $n / 2$, biasing HGPA toward preserving fine-grained partitions.
]
#proof[
  Each partition $q$ contributes exactly $k^((q))$ unit-weight hyperedges to $cal(E)$, so the fine sub-ensemble contributes $r_2 k_2$ total weight versus $r_1 k_1$ for the coarse sub-ensemble, irrespective of cluster sizes. Moreover, a candidate partition that aligns with a coarse $k_1$-partition may cut up to $k_2$ fine hyperedges (penalty $k_2$), whereas aligning with a fine $k_2$-partition (or a balanced merge of fine clusters) can cut at most $k_1$ coarse hyperedges (penalty $k_1 << k_2$). Weighting hyperedges inversely by $k^((q))$ or proportionally to $abs(e_a) / n$ is required to equalize voting influence across resolutions.
]

#sec("III", "Complexity, Geometry and Operational Profile")

A defining advantage of HGPA over CSPA is its *linear scaling in $n$*. The hypergraph $cal(H)$ contains $n$ vertices, $m = sum_(q=1)^r k^((q)) = O(k r)$ hyperedges, and total pin count $P = sum_(a=1)^m abs(e_a) = n r$ (for complete ensembles). In multilevel HMETIS @karypis1998metis, each coarsening level shrinks $abs(V)$ geometrically while touching each pin $O(1)$ times, and Fiduccia–Mattheyses gain-bucket refinement updates each cut hyperedge in amortized $O(1)$ time per pin, yielding $O(n r log k)$ runtime for recursive bisection and $O(n k r)$ worst-case complexity for direct $k$-way partitioning @strehl2002cluster. Space complexity is strictly $O(n r)$ to store vertex-to-hyperedge and hyperedge-to-vertex adjacency lists, avoiding CSPA's $Theta(n^2)$ memory barrier entirely.

Despite superior asymptotic scaling, pointer-based hypergraph coarsening, gain-bucket updates, and repeated randomized multi-start runs incur larger constant factors than cache-friendly dense matrix operations on small-to-medium datasets ($n <= 2000$), explaining why HGPA can run slower than CSPA and MCLA on compact benchmarks while dominating CSPA at $n >= 10^5$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n k r)$ worst-case @strehl2002cluster; $O(n r log k)$ multilevel bisection],
    [Space Complexity], [$O(n r)$ linear in pin count (no $n times n$ matrix)],
    [Geometric Bias], [Strongly favors equal-sized clusters; zero-gradient under boundary noise],
    [Hyperparameters], [Base ensemble $Lambda$, target $k$, balance $epsilon_"bal"$, HMETIS trials],
    [Determinism], [Stochastic (randomized coarsening & initial bisection); multi-start needed],
    [Invariances], [Label-permutation invariant; feature-free; non-invariant to $k^((q))$ mix],
    [Cut Geometry], [Submodular unit cut; up to $Theta(n / k)$ distortion vs. clique cut (Lem. 3)],
  ),
  caption: [Operational and geometric profile of HGPA (Alg. 172).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. VLSI Netlist and Parallel-Computing Partitioning Ensembles:* HGPA inherits its core solver directly from VLSI circuit partitioning @karypis1998metis @hendrickson1995chaco. When combining multiple coarse-grained logical groupings of gates or mesh elements into $k$ balanced processor blocks, strict balance $abs(C_c) approx n / k$ is a hard physical requirement rather than a liability, and minimizing cut nets directly minimizes inter-chip communication.

*2. Large-Scale Document and Web-Graph Taxonomy Fusion:* For $n > 10^5$ documents clustered independently across multiple topic models or community detectors @blondel2008louvain, forming CSPA's $n times n$ matrix is impossible. HGPA's $O(n r)$ pin representation processes millions of memberships natively, especially when base clusters are clean, fine-grained micro-clusters ($k^((q)) >> k$) that rarely suffer boundary-swap saturation.

*3. Distributed Federated Clustering With Minimal Bandwidth:* When remote data silos transmit sparse cluster-membership lists to a central coordinator, each cluster is already a hyperedge; HGPA ingests the streamed incidence lists directly into sparse CSR / CSC structures without expanding cross-silo pairwise affinities.

*4. Balanced Workload Sharding Across Heterogeneous Partitions:* In cloud database sharding and distributed graph databases, combining historical query co-access clusters across $r$ time windows into $k$ equally loaded shards matches HGPA's strict load-balance objective and minimizes cross-shard transaction cuts.

*5. Supra-Consensus Portfolio Diversification:* Within Strehl and Ghosh's supra-consensus architecture @strehl2002cluster, HGPA serves as a fast, structurally orthogonal candidate alongside CSPA and MCLA. On ensembles composed of high-purity, fine-resolution partitions where hyperedges remain largely uncut, HGPA occasionally discovers cleaner balanced cuts than pairwise methods, and $"ANMI"$ validation safely discards it whenever saturation occurs.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

HGPA exhibits the narrowest safe operating envelope of the Strehl–Ghosh triad due to three mathematical pathologies:
(i) _Noise-induced cut saturation_ (Thm. 2): even $1\%$ random point swaps across base clusters cause almost every large hyperedge to span multiple blocks, flattening the objective landscape and yielding near-random partitions.
(ii) _Hard balance distortion_: like METIS-based CSPA, HGPA fails severely on skewed cluster sizes (e.g., power-law communities or outlier-rich cohorts).
(iii) _Coarse-ensemble degeneracy and resolution bias_ (Prop. 5): when $k^((q)) = k$ and $r$ is small, the hypergraph has very few hyperedges ($m = k r$) of huge cardinality ($approx n / k$), producing sparse, uninformative gain buckets; HGPA performs noticeably better when base clusterings over-segment the data ($k^((q)) >> k$).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Balanced Sharding / VLSI], [High: strict balance and net-cut match physical cost], [Verify hyperedges are not saturated],
    [Large-Scale Fine Ensembles ($k^((q)) >> k$)], [Medium-High: $O(n r)$ memory; small hyperedges resist saturation], [Combine in supra-consensus with MCLA],
    [Noisy General Data Clustering], [Low: unit cut saturates under minor boundary errors (Thm. 2)], [Prefer CSPA (Alg. 171) or MCLA (Alg. 173)],
    [Skewed / Imbalanced Clusters], [Poor: hard balance constraint slices dominant clusters], [Use Monti (Alg. 170) or modularity methods],
    [Soft / Fuzzy Membership Needs], [Unsuitable: yields hard vertex partition only], [Switch to MCLA (Alg. 173)],
  ),
  caption: [Domain suitability and failure-mode matrix for HGPA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

HGPA occupies the _hyperedge-cut_ vertex of the ensemble triangle alongside CSPA's _object-clique_ view (Alg. 171) and MCLA's _cluster-meta-graph_ view (Alg. 173). Whereas CSPA expands each hyperedge $e_a$ into an all-pairs clique (introducing quadratic storage but smooth $O(n)$ cut gradients), HGPA retains the exact hypergraph topology at the cost of a non-smooth binary indicator $bb(1)(lambda(e; cal(C)) >= 2)$.

#proposition(title: "Exact Zero-Cut Coincidence Across Ensemble Algorithms", origin: "Original Proof")[
  Suppose an ensemble $Lambda$ admits a balanced $k$-partition $cal(C)^*$ that is a coarsening of every base partition $lambda^((q))$ (i.e., every base cluster $C^((q))_a$ is wholly contained in a single block of $cal(C)^*$). Then $cal(C)^*$ simultaneously achieves the exact global minimum $"hcut"(cal(C)^*; cal(H)) = 0$ of HGPA and $"Cut"(cal(C)^*; S) = 0$ of CSPA, and maximizes $1 / r sum_q I(cal(C); lambda^((q))) = H(cal(C)^*) = log_2 k$.
]
#proof[
  If every $e_a = C^((q))_a$ satisfies $e_a subset.eq C^*_(c(a))$ for some $c(a) in {1, dots, k}$, then $lambda(e_a; cal(C)^*) = 1$ for all $a = 1, dots, m$, giving $"hcut"(cal(C)^*; cal(H)) = 0$. Furthermore, for any pair $(i, j)$ in distinct blocks of $cal(C)^*$, no hyperedge can contain both $i$ and $j$, so $S_(i j) = 0$ and $"Cut"(cal(C)^*; S) = 0$. By Theorem 1, $H(cal(C)^* | lambda^((q))) = 0$ for every $q$, so $(1 / r) sum_q I(cal(C)^*; lambda^((q))) = H(cal(C)^*) = log_2 k$ under exact balance.
]

Proposition 6 clarifies the sharp contrast between noiseless and noisy regimes: HGPA is exact when base clusters cleanly refine the true partition, yet degrades precipitously the moment boundary contamination violates pure containment (Theorem 2).

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* On the ClusterOverview 4-dataset benchmark, HGPA records: *Iris* $"Sil" = 0.239$, $"ARI" = 0.492$; *Customers* $"Sil" = 0.294$; *Multi-Shape* $"Sil" = 0.563$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.259$ ($"ARI" = "N/A"$); with mean runtime $4.296$ ms. Consistent with Strehl and Ghosh's own empirical comparisons @strehl2002cluster and directly predicted by Theorem 2, HGPA scores lowest among Algs. 170–173 on *Iris* (trailing CSPA and MCLA by $-0.211$ in Silhouette and Monti by $-0.142$ in ARI), *Customers* ($-0.073$ behind CSPA), and *Multi-Shape* ($-0.149$ behind CSPA), while exhibiting the highest mean runtime ($4.296$ ms) due to multilevel hypergraph coarsening overhead at moderate $n$. Only on *High-D* ($"Sil" = 0.259$) does HGPA tie MCLA ($0.259$) and marginally edge Monti ($0.257$), reflecting that high-dimensional dispersion already weakens base-partition structure across all methods. In practice, HGPA should rarely be deployed standalone; rather, it functions best as an $O(n k r)$ candidate inside a supra-consensus selector audited by $"ANMI"$ or Silhouette @rousseeuw1987silhouettes.
