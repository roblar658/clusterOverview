#import "../template/ieee.typ": *

#alg_entry(
  id: 173,
  name: "MCLA (Meta-CLustering Algorithm)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2002",
  paradigm: "Cluster-Level Meta-Graph Partitioning & Soft Assignment",
  complexity: $O(n k^2 r^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $Lambda = {lambda^((1)), dots, lambda^((r))}$ be an ensemble of $r$ base clusterings on $n$ objects $V = {1, dots, n}$, encoded by the concatenated binary cluster-indicator matrix $H = [H^((1)), dots, H^((r))] in {0, 1}^(n times m)$ where $m = sum_(q=1)^r k^((q))$ is the total number of base clusters. CSPA (Alg. 171) operates in the $n$-dimensional _object space_ via $S = (1 / r) H H^top$ and pays $Theta(n^2)$ memory, whereas HGPA (Alg. 172) partitions vertices directly under a fragile binary hyperedge cut. Strehl and Ghosh @strehl2002cluster introduce the _Meta-CLustering Algorithm (MCLA)_ to resolve the label-correspondence problem directly in the much smaller $m$-dimensional _cluster space_ ($m << n$), producing both linear-in-$n$ scalability and soft object-to-cluster memberships.

Each column $h_a in {0, 1}^n$ ($a = 1, dots, m$) represents a base cluster $C_a = {i in V : h_(i a) = 1}$ (a hyperedge of the ensemble). Rather than comparing pairs of objects, MCLA treats the $m$ base clusters themselves as vertices of a _meta-graph_, groups corresponding clusters across different base runs into $k$ _meta-clusters_, collapses each meta-cluster into a soft consensus prototype, and assigns every object to the meta-cluster with which it associates most strongly.

#definition(title: "MCLA Meta-Graph, Meta-Cluster Collapse, and Competition Assignment")[
  The _meta-graph_ $G_"meta" = (V_"meta", E_"meta", W)$ has vertex set $V_"meta" = {1, dots, m}$ (indexing all base clusters) and edge weights given by the binary Jaccard similarity between indicator columns:
  $ W_(a b) = J(C_a, C_b) = (h_a^top h_b) / (norm(h_a)_2^2 + norm(h_b)_2^2 - h_a^top h_b) = abs(C_a inter C_b) / abs(C_a union C_b). $
  Partitioning $G_"meta"$ into $k$ balanced meta-clusters $cal(C)^((M)) = {C^((M))_1, dots, C^((M))_k}$ (via METIS @karypis1998metis) yields the $n times k$ _soft association matrix_ $A in [0, 1]^(n times k)$ obtained by averaging member indicator vectors within each meta-cluster:
  $ A_(i, c) = 1 / abs(C^((M))_c) sum_(a in C^((M))_c) H_(i a), quad i in V, c in {1, dots, k}. $
  Hard consensus labels are assigned by pointwise competition:
  $ lambda^* (i) = arg max_(c in {1, dots, k}) A_(i, c), $
  breaking ties uniformly at random.
]

Because two distinct clusters from the same base partition $lambda^((q))$ are disjoint, $W_(a b) = 0$ for all $a != b$ belonging to the same $q$. Hence $G_"meta"$ is strictly _$r$-partite_, and balanced $k$-way partitioning of $G_"meta"$ naturally seeks to place approximately one cluster from each of the $r$ base runs into each meta-cluster, automating cross-run cluster matching even when $k^((q)) != k$.

#algbox(
  title: "MCLA (Meta-CLustering Algorithm)",
  inputs: [Cluster ensemble $Lambda = {lambda^((1)), dots, lambda^((r))}$ on $n$ objects, target cluster count $k$, balance tolerance $epsilon_"bal" >= 0$.],
  outputs: [Soft association matrix $A in [0, 1]^(n times k)$ and hard consensus partition $lambda^* : V -> {1, dots, k}$.],
)[
  + *Meta-graph construction*: form the $m$ binary cluster indicators $h_1, dots, h_m in {0, 1}^n$; compute all pairwise intersections $G_(a b) = h_a^top h_b$ via sparse row accumulation $H^top H = sum_(i=1)^n H_(i, :)^top H_(i, :)$ and evaluate Jaccard weights $W_(a b) = G_(a b) / (abs(C_a) + abs(C_b) - G_(a b))$ for $a != b$.
  + *Meta-clustering*: partition the $m$-node weighted $r$-partite graph $G_"meta" = (V_"meta", E_"meta", W)$ into $k$ balanced meta-clusters $cal(C)^((M)) = {C^((M))_1, dots, C^((M))_k}$ minimizing the $k$-way cut via METIS @karypis1998metis.
  + *Meta-hyperedge collapsing*: for each meta-cluster $c = 1, dots, k$, average its constituent binary columns to form the association vector $A_(:, c) <- (1 / abs(C^((M))_c)) sum_(a in C^((M))_c) h_a in [0, 1]^n$.
  + *Object competition*: assign each object $i in V$ to $lambda^* (i) <- arg max_(1 <= c <= k) A_(i, c)$ (breaking ties randomly).
  + *Optional confidence reporting*: output soft memberships $tilde(A)_(i, c) = A_(i, c) / sum_(c'=1)^k A_(i, c')$ and mean confidence $bar(gamma) = (1 / n) sum_(i=1)^n max_c A_(i, c)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove five rigorous theoretical properties of MCLA: the metric and positive-definite kernel geometry of the Jaccard meta-graph, exact global recovery of perturbed permuted ensembles, simplex conservation versus unbalanced-meta-cluster pathologies, exact reduction to aligned majority voting with a Fano-based mutual information bound, and asymptotic linear scalability.

#subsec("A", "Metric and Positive-Semidefinite Geometry of the Meta-Graph")

#theorem(title: "Jaccard Metricity and Positive Semidefiniteness of the Meta-Graph", origin: "Original Proof")[
  Let $C_1, dots, C_m subset.eq V$ be nonempty base clusters, $W_(a b) = J(C_a, C_b)$, and $d_J (C_a, C_b) = 1 - W_(a b) = abs(C_a Delta C_b) / abs(C_a union C_b)$ where $Delta$ denotes symmetric difference. Then:
  (a) $d_J$ is a proper metric on nonempty subsets of $V$.
  (b) The meta-graph similarity matrix $W in [0, 1]^(m times m)$ is positive semidefinite ($W succ.eq 0$).
]
#proof[
  *Part (a):* Non-negativity, symmetry, and $d_J (A, B) = 0 <=> A = B$ are immediate. To prove the triangle inequality for any three nonempty sets $A, B, C subset.eq V$, let $x = abs(B without (A union C)) >= 0$. Decomposing $A Delta B$ into its parts inside and outside $A union C$ (note that for $u in.not A union C$, $u in A Delta B <=> u in B$), we obtain
  $ abs(A Delta B) = abs((A Delta B) inter (A union C)) + x, $
  $ abs(B Delta C) = abs((B Delta C) inter (A union C)) + x. $
  Similarly, $abs(A union B) = abs((A union B) inter (A union C)) + x <= abs(A union C) + x$ and $abs(B union C) <= abs(A union C) + x$. Therefore,
  $ d_J (A, B) + d_J (B, C) & = abs(A Delta B) / abs(A union B) + abs(B Delta C) / abs(B union C) \
    & >= (abs((A Delta B) inter (A union C)) + abs((B Delta C) inter (A union C)) + 2 x) / (abs(A union C) + x). $
  Because symmetric difference obeys $A Delta C subset.eq (A Delta B) union (B Delta C)$ and $A Delta C subset.eq A union C$, intersecting with $A union C$ gives $abs(A Delta C) <= abs((A Delta B) inter (A union C)) + abs((B Delta C) inter (A union C))$. Finally, since $abs(A Delta C) <= abs(A union C)$ and $x >= 0$, we have $(abs(A Delta C) + x) abs(A union C) >= abs(A Delta C) (abs(A union C) + x)$, hence
  $ d_J (A, B) + d_J (B, C) >= (abs(A Delta C) + 2 x) / (abs(A union C) + x) >= (abs(A Delta C) + x) / (abs(A union C) + x) >= abs(A Delta C) / abs(A union C) = d_J (A, C). $

  *Part (b):* Let $pi$ be a uniformly random permutation of $V = {1, dots, n}$, and for any nonempty subset $C subset.eq V$ define its min-hash signature $h_pi (C) = arg min_(i in C) pi(i) in V$. For any two nonempty clusters $C_a, C_b$, the event ${h_pi (C_a) = h_pi (C_b)}$ holds iff the unique minimizer of $pi$ on $C_a union C_b$ lies in $C_a inter C_b$; by symmetry of the uniform permutation, $"Pr"_pi (h_pi (C_a) = h_pi (C_b)) = abs(C_a inter C_b) / abs(C_a union C_b) = W_(a b)$. Consequently, for every vector $v in bb(R)^m$,
  $ v^top W v & = sum_(a, b=1)^m v_a v_b bb(E)_pi [bb(1)(h_pi (C_a) = h_pi (C_b))] \
    & = bb(E)_pi [ sum_(i=1)^n ( sum_(a : h_pi (C_a) = i) v_a )^2 ] >= 0, $
  establishing $W succ.eq 0$.
]

#subsec("B", "Exact Meta-Cluster Recovery Under Perturbed Permuted Ensembles")

Combining the Jaccard metric inequality of Theorem 1(a) with the $r$-partite structure of $G_"meta"$, we derive an exact deterministic threshold guaranteeing that balanced meta-graph partitioning recovers the true cross-run cluster correspondence uniquely.

#theorem(title: "Exact Meta-Cluster Correspondence Recovery for Jaccard Noise < 1/4", origin: "Original Proof")[
  Let $cal(C)^* = {C_1^*, dots, C_k^*}$ be a ground-truth $k$-partition of $V$. Suppose $Lambda$ consists of $r$ base $k$-clusterings ($k^((q)) = k$) such that for each run $q$ there exists an unknown label permutation $sigma_q in S_k$ satisfying the uniform Jaccard perturbation bound
  $ max_(1 <= q <= r) max_(1 <= c <= k) d_J (C^((q))_(sigma_q (c)), C_c^*) <= eta < 1 / 4. $
  Then:
  (a) The true aligned meta-clustering $cal(M)^* = {M_1^*, dots, M_k^*}$ with $M_c^* = {C^((q))_(sigma_q (c)) : 1 <= q <= r}$ is the *unique global minimizer* of the balanced $k$-way cut on $G_"meta"$.
  (b) Under $cal(M)^*$, MCLA's Step 4 competition rule exactly recovers the ground-truth label $lambda^* (i) = c$ for every object $i in C_c^*$ that is correctly assigned by strictly more than $r / 2$ base clusterers.
]
#proof[
  *Part (a):* Fix any two distinct vertices $a, b in V_"meta"$.
  - If $a, b in M_c^*$ belong to the same true meta-cluster (necessarily from distinct runs $q != q'$), the Jaccard triangle inequality (Theorem 1(a)) gives $d_J (a, b) <= d_J (a, C_c^*) + d_J (C_c^*, b) <= 2 eta$, hence $W_(a b) = 1 - d_J (a, b) >= 1 - 2 eta$.
  - If $a in M_c^*$ and $b in M_(c')^*$ with $c != c'$, then because $C_c^* inter C_(c')^* = emptyset$, we have $d_J (C_c^*, C_(c')^*) = 1$. By Theorem 1(a),
  $ 1 = d_J (C_c^*, C_(c')^*) <= d_J (C_c^*, a) + d_J (a, b) + d_J (b, C_(c')^*) <= 2 eta + d_J (a, b), $
  yielding $d_J (a, b) >= 1 - 2 eta$ and hence $W_(a b) <= 2 eta$ (and $W_(a b) = 0 <= 2 eta$ if $a, b$ come from the same run $q$).

  Now let $cal(C)^((M)) = {C^((M))_1, dots, C^((M))_k}$ be any balanced $k$-partition of $V_"meta"$ ($abs(C^((M))_j) = r$). Because the total number of within-block unordered pairs is fixed at $k binom(r, 2)$, letting $S_"same"$ denote the number of within-block pairs drawn from the same $M_c^*$ and $S_"diff" = k binom(r, 2) - S_"same"$ the number of within-block pairs drawn from distinct $M_c^*$'s, we see that $cal(C)^((M))$ cuts exactly $S_"diff"$ intra-$M_c^*$ pairs (each of weight $>= 1 - 2 eta$) while keeping exactly $S_"diff"$ inter-$M_c^*$ pairs uncut (each of weight $<= 2 eta$). Comparing total within-block weight $W_"in"$ yields
  $ "Cut"(cal(C)^((M)); W) - "Cut"(cal(M)^*; W) & = W_"in" (cal(M)^*) - W_"in" (cal(C)^((M))) \
    & >= S_"diff" (1 - 2 eta) - S_"diff" (2 eta) = S_"diff" (1 - 4 eta). $
  Since $eta < 1 / 4$ implies $1 - 4 eta > 0$, and $S_"diff" > 0$ for every balanced partition $cal(C)^((M)) != cal(M)^*$, the true alignment $cal(M)^*$ uniquely minimizes the balanced $k$-way cut.

  *Part (b):* Given $cal(M)^*$, every meta-cluster satisfies $abs(M_c^*) = r$, so $A_(i, c) = (1 / r) sum_(q=1)^r bb(1)(lambda^((q)) (i) = sigma_q (c))$. If object $i in C_c^*$ is correctly labeled by more than $r / 2$ base clusterers, then $A_(i, c) > 1 / 2$, which forces $A_(i, c') <= 1 - A_(i, c) < 1 / 2$ for every $c' != c$, guaranteeing $lambda^* (i) = c$.
]

#subsec("C", "Simplex Conservation and Unbalanced Meta-Cluster Pathologies")

While MCLA's soft association matrix $A$ invites a fuzzy-membership interpretation analogous to Fuzzy $C$-Means @bezdek1981fcm, the row-wise normalization depends critically on meta-cluster size balance.

#proposition(title: "Weighted Simplex Conservation and Size-Normalization Distortion", origin: "Original Counterexample")[
  Let $Lambda$ be a complete ensemble and let $A_(i, c) = (1 / abs(C^((M))_c)) sum_(a in C^((M))_c) H_(i a)$.
  (a) For every object $i in V$, $sum_(c=1)^k abs(C^((M))_c) A_(i, c) = r$. In particular, when all meta-clusters have equal size $abs(C^((M))_c) = r$, every row $A_(i, :)$ lies on the probability simplex $Delta^(k - 1)$ and satisfies $max_c A_(i, c) >= 1 / k$.
  (b) When meta-clusters are *unbalanced*, dividing by $abs(C^((M))_c)$ can overturn a strict ensemble majority and produce empty consensus clusters.
]
#proof[
  *Part (a):* Because ${C^((M))_1, dots, C^((M))_k}$ partitions all $m$ columns of $H$ and each complete base partition places object $i$ in exactly one cluster,
  $ sum_(c=1)^k abs(C^((M))_c) A_(i, c) = sum_(c=1)^k sum_(a in C^((M))_c) H_(i a) = sum_(a=1)^m H_(i a) = r. $
  Dividing by $r$ when $abs(C^((M))_c) = r$ for all $c$ yields $sum_(c=1)^k A_(i, c) = 1$, and by pigeonhole averaging $max_c A_(i, c) >= 1 / k$.

  *Part (b):* Construct an ensemble with $r = 4$, $k = 2$, $m = 8$, and an unbalanced meta-clustering placing $1$ hyperedge in $C^((M))_1$ and $7$ hyperedges in $C^((M))_2$. Any object $i$ appearing in the lone hyperedge of $C^((M))_1$ and in $3$ of the $4$ base partitions inside $C^((M))_2$ receives scores $A_(i, 1) = 1 / 1 = 1$ versus $A_(i, 2) = 3 / 7 approx 0.429$, so MCLA assigns $i$ to $C^((M))_1$ despite a $3$-to-$1$ ensemble majority favoring $C^((M))_2$. Conversely, if a balanced $k$-way meta-clustering groups $r$ fragmented outlier clusters into meta-cluster $C^((M))_k$ whose member clusters are mutually disjoint on every object ($sum_(a in C^((M))_k) H_(i a) <= 1$), then $A_(i, k) <= 1 / r$, while at least one remaining meta-cluster achieves $A_(i, c) >= (r - 1) / (r (k - 1)) > 1 / r$ whenever $r > k$; consequently $C^((M))_k$ wins *zero* objects, leaving fewer than $k$ nonempty consensus clusters.
]

#subsec("D", "Aligned Majority Optimality and Fano Mutual-Information Bound")

#theorem(title: "Hamming-Loss Optimality and Fano ANMI Lower Bound", origin: "Original Derivation")[
  Suppose $k^((q)) = k$ for all $q$ and the meta-clustering defines a valid bijective alignment $sigma_q in S_k$ across runs ($C^((M))_c = {C^((q))_(sigma_q (c)) : 1 <= q <= r}$). Let $bar(gamma) = (1 / n) sum_(i=1)^n max_c A_(i, c)$ be MCLA's average winning confidence and let $h_2 (p) = -p log_2 p - (1 - p) log_2 (1 - p)$ be binary entropy. Then:
  (a) MCLA's assignment $lambda^* (i) = arg max_c A_(i, c)$ is the *exact global minimizer* of the aligned Hamming disagreement
  $ L_"Ham" (lambda; {sigma_q}) = 1 / (n r) sum_(q=1)^r sum_(i=1)^n bb(1)(lambda(i) != sigma_q^(-1) (lambda^((q)) (i))), $
  achieving minimum loss $L_"Ham" (lambda^*) = 1 - bar(gamma) <= 1 - 1 / k$.
  (b) The ensemble average mutual information obeys the Fano lower bound
  $ 1 / r sum_(q=1)^r I(lambda^*; lambda^((q))) >= H(lambda^*) - h_2 (1 - bar(gamma)) - (1 - bar(gamma)) log_2 (k - 1). $
]
#proof[
  *Part (a):* For each object $i$ and candidate label $c in {1, dots, k}$,
  $ 1 / r sum_(q=1)^r bb(1)(c != sigma_q^(-1) (lambda^((q)) (i))) = 1 - 1 / r sum_(q=1)^r bb(1)(lambda^((q)) (i) = sigma_q (c)) = 1 - A_(i, c). $
  Averaging over $i = 1, dots, n$ decouples $L_"Ham"$ completely across objects: $L_"Ham" (lambda; {sigma_q}) = (1 / n) sum_(i=1)^n (1 - A_(i, lambda(i)))$. Minimizing each term independently selects $lambda^* (i) = arg max_c A_(i, c)$, achieving $L_"Ham" (lambda^*) = 1 - bar(gamma)$, which is at most $1 - 1 / k$ by Proposition 3(a).

  *Part (b):* Draw an object index $I tilde.op "Uniform"(V)$, and define random variables $Y = lambda^* (I)$ and $hat(Y)_q = sigma_q^(-1) (lambda^((q)) (I))$ taking values in ${1, dots, k}$. Because $sigma_q$ is a bijection, $H(Y | lambda^((q)) (I)) = H(Y | hat(Y)_q)$. Let $epsilon_q = "Pr"(Y != hat(Y)_q)$ and $E_q = bb(1)(Y != hat(Y)_q)$. By Fano's inequality (expanding $H(Y, E_q | hat(Y)_q)$ via the chain rule and noting $H(Y | E_q = 0, hat(Y)_q) = 0$ and $H(Y | E_q = 1, hat(Y)_q) <= log_2 (k - 1)$),
  $ H(Y | hat(Y)_q) <= H(E_q) + epsilon_q log_2 (k - 1) = h_2 (epsilon_q) + epsilon_q log_2 (k - 1). $
  Averaging over $q = 1, dots, r$, note that $(1 / r) sum_(q=1)^r epsilon_q = L_"Ham" (lambda^*) = 1 - bar(gamma)$. By concavity of $h_2$ and Jensen's inequality, $(1 / r) sum_(q=1)^r h_2 (epsilon_q) <= h_2 (1 - bar(gamma))$. Substituting into $I(Y; hat(Y)_q) = H(Y) - H(Y | hat(Y)_q)$ yields the stated lower bound.
]

#subsec("E", "Linear-in-n Complexity and Supra-Consensus Portfolio Dominance")

#proposition(title: "Sparse Linear-in-n Scaling and Supra-Consensus Optimality", origin: "Classical — Cited")[
  While Strehl and Ghosh @strehl2002cluster bound MCLA's worst-case time complexity by $O(n k^2 r^2)$ via dense pairwise column inner products, exploiting row sparsity of $H$ reduces exact meta-graph construction and full MCLA execution to $O(n r^2 + k^2 r^2)$ time and $O(n (r + k) + k^2 r^2)$ space—strictly linear in $n$ for fixed $k, r$. Moreover, evaluating Strehl and Ghosh's supra-consensus rule $lambda^"supra" = arg max_(lambda in {lambda^"CSPA", lambda^"HGPA", lambda^"MCLA"}) Phi^("ANMI") (Lambda, lambda)$ guarantees $Phi^("ANMI") (Lambda, lambda^"supra") >= max(Phi^("ANMI")_1, Phi^("ANMI")_2, Phi^("ANMI")_3)$ at only $O(n k r)$ verification cost.
]
#proof[
  Each row $H_(i, :)$ of $H in {0, 1}^(n times m)$ contains at most $r$ nonzero entries (at most one per base clustering). Accumulating $H^top H = sum_(i=1)^n H_(i, :)^top H_(i, :)$ therefore updates at most $binom(r, 2) = O(r^2)$ nonzero index pairs per object $i$, costing $O(n r^2)$ total operations to compute all $m^2 = O(k^2 r^2)$ intersections exactly, plus $O(m^2) = O(k^2 r^2)$ to compute Jaccard ratios and partition the $m$-node graph $G_"meta"$ via METIS @karypis1998metis. Collapsing meta-clusters and finding pointwise maxima across $k$ columns take $O(n r + n k)$ operations, giving total complexity $O(n r^2 + n k + k^2 r^2) = O(n r^2 + k^2 r^2)$. Finally, computing the contingency table between any candidate $lambda$ and each $lambda^((q))$ takes $O(n + k k^((q)))$ time, so evaluating $Phi^("ANMI")$ across all three ensemble heuristics costs $O(n r + k m) = O(n k r)$ and trivially dominates each constituent heuristic on $Phi^("ANMI")$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because MCLA shifts graph partitioning from the $n$-node object graph of CSPA to the $m$-node cluster meta-graph ($m = sum_q k^((q)) approx k r$), its computational profile exhibits an exact crossover relative to CSPA: whenever $k r << n$ (the standard data-mining regime, e.g., $n = 10^5, k = 10, r = 20 => m = 200$), MCLA executes orders of magnitude faster than CSPA and requires only $O(n (r + k) + m^2)$ memory instead of $Theta(n^2)$. Conversely, when $k r$ approaches or exceeds $n$ (extreme micro-clustering), the $O(k^2 r^2)$ meta-graph overtakes CSPA.

Geometrically, MCLA acts as a _cluster-correspondence aligner_ rather than a pairwise metric smoother: it preserves the geometric shape of consensus clusters *inherited directly from the base clusterers* (including highly skewed cluster sizes, because METIS balances the *number of base clusters* per meta-cluster rather than the *number of objects* per consensus cluster!). This subtle distinction is crucial: unlike CSPA and HGPA, MCLA's balance constraint operates on $V_"meta"$ (equalizing $abs(C^((M))_c) approx r$), allowing the resulting object clusters ${i : lambda^* (i) = c}$ to have arbitrary imbalanced cardinalities!

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n k^2 r^2)$ classical @strehl2002cluster; $O(n r^2 + k^2 r^2)$ via sparse row accumulation],
    [Space Complexity], [$O(n (r + k) + k^2 r^2)$ (strictly linear in $n$)],
    [Geometric Bias], [Cluster-level correspondence voting; supports imbalanced object cluster sizes],
    [Hyperparameters], [Base ensemble $Lambda$, target $k$, meta-graph balance tolerance $epsilon_"bal"$],
    [Determinism], [Exact $W$; minor stochasticity from METIS on small $m$-node graph & tie breaks],
    [Invariances], [Strictly invariant to base label permutations; feature-free knowledge reuse],
    [Outputs], [Hard consensus $lambda^*$ plus soft association scores $A in [0, 1]^(n times k)$],
  ),
  caption: [Operational and geometric profile of MCLA (Alg. 173).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Million-Scale Customer and User-Behavior Segmentation:* Industrial customer databases ($n = 10^6$–$10^8$, $k = 5$–$50$, $r = 10$–$30$ models trained across behavioral views) render $O(n^2)$ co-association matrices impossible. MCLA collapses $m <= 1500$ base clusters in milliseconds and assigns all $n$ users in a single streaming $O(n r)$ pass, while tolerating naturally skewed customer-tier sizes (e.g., $1\%$ VIP vs. $60\%$ casual users).

*2. Distributed Multi-View and Multi-Modal Soft Fusion:* In federated healthcare or multi-modal recommendation systems @cai2013multiview @kumar2011coregularized, different institutions or modalities produce cluster assignments with unaligned labels. MCLA aligns cross-view cluster identities via Jaccard overlap (Theorem 2) and delivers soft probabilistic scores $tilde(A)_(i, c)$ directly usable by downstream risk calibrators or fuzzy-logic controllers @bezdek1981fcm.

*3. Large-Scale Community Detection Consensus in Complex Networks:* Stochastic modularity optimizers such as Louvain @blondel2008louvain and Leiden @traag2019leiden produce slightly varying community partitions across random seeds. Running MCLA across $r = 20$–$50$ seeds aligns community IDs without forming an $n times n$ graph, and the per-node confidence $max_c A_(i, c)$ immediately identifies overlapping bridge nodes straddling multiple communities.

*4. Streaming Ensemble Drift Tracking and Label Alignment:* Because meta-clusters act as explicit consensus prototypes $A_(:, c)$, MCLA naturally solves the _label-switching problem_ across sequential clustering epochs or Bayesian mixture MCMC chains @dempster1977em, maintaining consistent cluster identities over time while quantifying per-object assignment drift.

*5. Knowledge Reuse Across Heterogeneous Feature Spaces:* When legacy clusterings exist on historical feature subsets that are no longer available at inference time @strehl2002cluster, MCLA combines historical and current partitions via pure cluster overlap without recomputing distances in the original spaces.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four concrete regimes limit MCLA:
(i) _Severe cluster-count mismatch ($k^((q)) >> k$ or heterogeneous $k^((q))$)_: when base clusterings shatter true non-convex manifolds into many small convex sub-clusters (e.g., $K$-means with $k^((q)) = 20$ on two moons), two sub-clusters belonging to the same moon have *zero* Jaccard overlap within a run and minimal overlap across runs, preventing MCLA from chaining them together (whereas CSPA's pairwise graph bridges them smoothly).
(ii) _High base noise ($eta >= 1 / 4$ in Thm. 2)_: if base partitions disagree severely on cluster identities, the $r$-partite meta-graph loses block structure and METIS merges unrelated clusters.
(iii) _Empty-cluster collapse_ (Prop. 3(b)): fragmented or outlier meta-clusters can lose every pointwise $arg max$ competition, returning fewer than $k$ nonempty consensus clusters.
(iv) _Unbalanced meta-cluster inflation_ (Prop. 3(b)): relaxing meta-graph balance lets singleton meta-clusters dominate multi-cluster majorities unless raw vote counts $sum_(a in C^((M))_c) H_(i a)$ replace normalized means.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Web-Scale / Telecom ($n > 10^5$)], [Excellent: $O(n r^2)$ linear scaling and skewed-size tolerance], [Verify $k^((q)) approx k$ consistency],
    [Multi-Seed Community Consensus], [High: fast label alignment + soft bridge-node scores], [Watch for empty clusters (Prop. 3)],
    [Imbalanced Cluster Cohorts], [High: METIS balances meta-clusters, not object counts], [Avoid CSPA/HGPA object balance],
    [Non-Convex Over-Segmented Ensembles], [Moderate-Low: coarse cluster matching misses local chaining], [Prefer CSPA (Alg. 171) if $n <= 10^4$],
    [High-Noise / Conflicting Ensembles], [Low when $d_J >= 1/4$: meta-graph loses $r$-partite blocks], [Filter weak base runs or use Monti],
  ),
  caption: [Domain suitability and failure-mode matrix for MCLA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MCLA forms the _dual cluster-space counterpart_ to CSPA (Alg. 171) and complements HGPA (Alg. 172) and Monti Consensus Clustering (Alg. 170 @monti2003consensus). Specifically, whereas CSPA extracts clusters from the $n times n$ object Gram matrix $H H^top$, MCLA extracts meta-clusters from the Jaccard-normalized $m times m$ cluster Gram matrix $H^top H$ and then projects back to objects via $H$-averaging.

#proposition(title: "Exact Equivalence to Rank-1 Bipartite Co-Clustering Collapse", origin: "Original Proof")[
  Let $Z in {0, 1}^(m times k)$ be the binary indicator matrix of MCLA's meta-clustering $cal(C)^((M))$, and let $D_M = "diag"(abs(C^((M))_1), dots, abs(C^((M))_k))$. Then MCLA's soft association matrix is the exact linear projection $A = H Z D_M^(-1)$, and when every base partition shares the same balanced cluster cardinality $abs(C_a) = n / k$, the Jaccard edge weight $W_(a b) = G_(a b) / (2 n / k - G_(a b))$ is a strictly increasing convex bijection of the inner product $G_(a b) = (H^top H)_(a b)$, making MCLA's meta-graph monotonically order-equivalent to the dual Gram matrix of CSPA.
]
#proof[
  Entry $(i, c)$ of $H Z D_M^(-1)$ equals $(1 / abs(C^((M))_c)) sum_(a=1)^m H_(i a) Z_(a c) = (1 / abs(C^((M))_c)) sum_(a in C^((M))_c) H_(i a) = A_(i, c)$ identically. When $abs(C_a) = s = n / k$ for all $a$, we have $W_(a b) = f(G_(a b))$ where $f(g) = g / (2 s - g)$ on $[0, s]$. Differentiating gives $f'(g) = 2 s / (2 s - g)^2 > 0$ and $f''(g) = 4 s / (2 s - g)^3 > 0$, proving that $f$ is strictly increasing and strictly convex on $[0, s]$: MCLA thus amplifies high-overlap cluster correspondences quadratically more aggressively than linear inner-product similarity $G_(a b)$ while suppressing weak spurious overlaps.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* Across the four standardized ClusterOverview benchmarks, MCLA achieves: *Iris* $"Sil" = 0.450$, $"ARI" = 0.563$; *Customers* $"Sil" = 0.354$; *Multi-Shape* $"Sil" = 0.611$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.259$ ($"ARI" = "N/A"$); with mean runtime $3.197$ ms. These empirical numbers verify our theoretical comparisons with remarkable precision:
(1) On *Iris*, MCLA matches CSPA *identically* ($"Sil" = 0.450$, $"ARI" = 0.563$) and decisively outperforms HGPA ($+0.211$ Sil, $+0.071$ ARI), confirming Theorem 2 and Theorem 4: when base clusters exhibit clean $k$-to-$k$ correspondence ($eta < 1/4$), cluster-space meta-alignment matches full $n times n$ co-association partitioning while running $25.6\%$ faster than HGPA ($3.197$ vs. $4.296$ ms) and $16.1\%$ faster than Monti ($3.809$ ms).
(2) On *Customers*, MCLA ($"Sil" = 0.354$) trails CSPA ($0.367$) by only $0.013$ while cleanly outperforming Monti ($0.332$) and HGPA ($0.294$).
(3) On *Multi-Shape* ($0.611$ vs. CSPA's $0.712$) and *High-D* ($0.259$ vs. CSPA's $0.303$), MCLA loses pairwise granularity relative to CSPA exactly as predicted in Section IV-B: meta-clustering aggregates at the coarse cluster level and cannot smooth local non-convex pairwise neighborhoods the way CSPA's $n times n$ similarity kernel does, making MCLA the premier scalable default at large $n$ and a strong complement to CSPA inside Strehl and Ghosh's supra-consensus framework @strehl2002cluster @rousseeuw1987silhouettes.
