#import "../template/ieee.typ": *

#alg_entry(
  id: 85,
  name: "COBWEB",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1987",
  paradigm: "Incremental Conceptual Clustering via Category Utility",
  complexity: $cal(O)(n B^2 d V h)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

COBWEB @fisher1987cobweb formulates unsupervised learning as _incremental conceptual clustering_: rather than partitioning a fixed batch of vectors in metric space, it constructs a probabilistic taxonomy tree $cal(T)$ online, one nominal instance at a time, such that every internal level maximizes predictive inference over all attributes simultaneously. Let $cal(D) = (x_1, dots, x_n)$ be an ordered stream of $n$ instances over $d$ nominal attributes $A_1, dots, A_d$, where attribute $A_i$ takes values in a finite domain $cal(V)_i = {v_(i, 1), dots, v_(i, m_i)}$ with $m_i = |cal(V)_i| <= V$.

Each node $C in cal(T)$ represents a probabilistic concept summarizing a multiset of $n_C >= 1$ incorporated instances through empirical frequency tables $n_(C, i, j) = |{x in C : x_i = v_(i, j)}|$. These sufficient statistics induce the empirical conditional probabilities
$ P(A_i = v_(i, j) | C) = n_(C, i, j) / n_C, quad sum_(j=1)^(m_i) P(A_i = v_(i, j) | C) = 1. $
For an internal node $C$ whose direct children induce a partition $cal(P) = {C_1, dots, C_K}$ of $C$, the child base rates are $P(C_k | C) = n_(C_k) / n_C$ (written $P(C_k)$ when $C$ is clear from context). Adapting the psychological categorization metric of Gluck and Corter, Fisher @fisher1987cobweb evaluates every candidate partition by its _Category Utility_.

#definition(title: "Category Utility")[
  Let $cal(P) = {C_1, dots, C_K}$ be a $K$-way partition of a parent node $C$ with attribute marginals $P(A_i = v_(i, j)) = P(A_i = v_(i, j) | C)$. The _Category Utility_ of $cal(P)$ is
  $ "CU"(cal(P)) = 1 / K sum_(k=1)^K P(C_k) [& sum_(i=1)^d sum_(j=1)^(m_i) P(A_i = v_(i, j) | C_k)^2 \
    &- sum_(i=1)^d sum_(j=1)^(m_i) P(A_i = v_(i, j))^2 ]. $
]

Probabilistically, under a probability-matching strategy that guesses $A_i = v_(i, j)$ with probability $P(A_i = v_(i, j) | C_k)$ inside cluster $C_k$, the inner sum equals the expected gain in correctly predicted attribute values over guessing from the parent marginals alone. Dividing by $K$ amortizes this gain across the branching factor, balancing within-cluster predictive specificity against taxonomic parsimony.

To maintain $cal(T)$ incrementally without storing or re-clustering past instances, COBWEB performs greedy hill-climbing at each visited internal node using four local tree-restructuring operators: _incorporating_ $x$ into an existing child, _creating_ a new singleton child, _merging_ the two best host children before incorporating $x$, and _splitting_ the best host child by promoting its grandchildren.

#algbox(
  title: "COBWEB Incremental Concept Formation",
  inputs: [Ordered stream $cal(D) = (x_1, dots, x_n)$ over $d$ nominal attributes $A_1, dots, A_d$.],
  outputs: [Probabilistic concept hierarchy $cal(T)$ with sufficient statistics $(n_C, {n_(C, i, j)})$ at every node $C in cal(T)$.],
)[
  + Initialize root $R$ with $x_1$; for each $t = 2, dots, n$, invoke `Insert`$(R, x_t)$.
  + *Procedure* `Insert`$(C, x)$: if $C$ is a leaf and $x$ matches $C$'s unique profile, increment $n_C$ and return; if $C$ is a leaf differing from $x$, expand $C$ into two leaf children holding the prior multiset and ${x}$, update counts at $C$, and return.
  + Otherwise ($C$ is internal with children $cal(P) = {C_1, dots, C_K}$), update $n_C$ and ${n_(C, i, j)}$ with $x$.
  + *Incorporate*: for each $k in {1, dots, K}$, compute $"CU"_k = "CU"(cal(P)_(k <- x))$ when adding $x$ to $C_k$; let $k_1, k_2$ index the highest and second-highest values.
  + *Create*: compute $"CU"_"new" = "CU"(cal(P) union {{x}})$ by adding a new singleton child.
  + *Merge*: if $K >= 2$, form $C_m = C_(k_1) union C_(k_2)$ (with $C_(k_1), C_(k_2)$ as children of $C_m$) and compute $"CU"_"merge" = "CU"((cal(P) without {C_(k_1), C_(k_2)}) union {C_m union {x}})$.
  + *Split*: if $C_(k_1)$ is internal with children $cal(S)_(k_1)$, form $cal(P)_"sp" = (cal(P) without {C_(k_1)}) union cal(S)_(k_1)$ and evaluate $"CU"_"split"$ after re-incorporating $x$ into $cal(P)_"sp"$.
  + Execute the operator attaining $max("CU"_(k_1), "CU"_"new", "CU"_"merge", "CU"_"split")$. If *Incorporate* or *Merge* wins, recurse on the chosen host; if *Split* wins, recurse on `Insert`$(C, x)$; if *Create* wins, terminate.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal characterizations of COBWEB: an exact equivalence to amortized Gini-impurity reduction, closed-form Ward-type merge and split decision laws, a minimal finite-sample counterexample proving order sensitivity despite bidirectional operators, tight incremental complexity bounds, and an exact algebraic characterization of COBWEB's conditional-independence inductive bias.

#subsec("A", "Gini-Impurity Equivalence and Exact Bounds")

#lemma(title: "Category Utility as Amortized Gini-Impurity Reduction", origin: "Original Proof")[
  For any node $C$ and attribute $A_i$, define the empirical Gini impurity $G_i (C) = 1 - sum_(j=1)^(m_i) P(A_i = v_(i, j) | C)^2$. Then for every $K$-way partition $cal(P) = {C_1, dots, C_K}$ of $C$,
  $ "CU"(cal(P)) = 1 / K sum_(i=1)^d ( G_i (C) - sum_(k=1)^K P(C_k) G_i (C_k) ). $
  Moreover, $"CU"(cal(P))$ satisfies the sharp two-sided bounds
  $ 0 <= "CU"(cal(P)) <= 1 / K sum_(i=1)^d G_i (C) <= d / K (1 - 1 / V), $
  with $"CU"(cal(P)) = 0$ if and only if $P(A_i = v_(i, j) | C_k) = P(A_i = v_(i, j))$ for every $i, j, k$, and with the middle upper bound attained if and only if every leaf distribution $P(dot.c | C_k)$ is pure (a point mass on every attribute).
]

#proof[
  Since $sum_(k=1)^K P(C_k) = 1$, rewriting the definition of $"CU"(cal(P))$ yields
  $ "CU"(cal(P)) &= 1 / K sum_(i=1)^d [ sum_(k=1)^K P(C_k) sum_(j=1)^(m_i) P(A_i = v_(i, j) | C_k)^2 \
    &quad - sum_(j=1)^(m_i) P(A_i = v_(i, j))^2 ] \
    &= 1 / K sum_(i=1)^d [ (1 - sum_(j=1)^(m_i) P(A_i = v_(i, j))^2) \
    &quad - sum_(k=1)^K P(C_k) (1 - sum_(j=1)^(m_i) P(A_i = v_(i, j) | C_k)^2) ], $
  which directly gives the Gini-impurity reduction identity.

  For the lower bound, observe that $P(A_i = v_(i, j)) = sum_(k=1)^K P(C_k) P(A_i = v_(i, j) | C_k)$ by the law of total probability. Because $phi(u) = u^2$ is strictly convex on $[0, 1]$, Jensen's inequality gives, for every $(i, j)$,
  $ sum_(k=1)^K P(C_k) P(A_i = v_(i, j) | C_k)^2 >= P(A_i = v_(i, j))^2, $
  with equality if and only if $P(A_i = v_(i, j) | C_k)$ is constant over all $k$ with $P(C_k) > 0$, i.e., equal to $P(A_i = v_(i, j))$. Summing over $i, j$ and dividing by $K$ yields $"CU"(cal(P)) >= 0$ with the stated equality condition.

  For the upper bounds, $G_i (C_k) >= 0$ for all $i, k$, with equality if and only if $P(A_i = v_(i, j) | C_k) in {0, 1}$ for all $j$, yielding $"CU"(cal(P)) <= K^(-1) sum_(i=1)^d G_i (C)$. Finally, by Cauchy-Schwarz on the probability simplex of dimension $m_i <= V$,
  $ sum_(j=1)^(m_i) P(A_i = v_(i, j))^2 >= 1 / m_i >= 1 / V, $
  so $G_i (C) <= 1 - 1 / V$, and summing over $i = 1, dots, d$ proves $"CU"(cal(P)) <= (d / K)(1 - 1 / V)$.
]

#subsec("B", "Exact Ward-Space Merge and Split Laws")

Let $M = sum_(i=1)^d m_i$ and associate each node $C_k$ with its concatenated empirical probability vector $bold(p)_k in [0, 1]^M$ whose $(i, j)$ entry is $P(A_i = v_(i, j) | C_k)$, and weight $w_k = P(C_k)$. Then $"CU"(cal(P)) = K^(-1) (sum_(k=1)^K w_k norm(bold(p)_k)_2^2 - norm(bold(p)_0)_2^2)$, where $bold(p)_0 = sum_(k=1)^K w_k bold(p)_k$.

#theorem(title: "Exact Algebraic Merge and Split Decision Criteria", origin: "Original Derivation")[
  Let $cal(P) = {C_1, dots, C_K}$ be a $K$-way partition ($K >= 2$).
  1. Merging two siblings $C_a, C_b in cal(P)$ into $C_(a b) = C_a union C_b$ produces partition $cal(P)_"mg"$ with exact utility change
     $ "CU"(cal(P)_"mg") - "CU"(cal(P)) = 1 / (K - 1) [ "CU"(cal(P)) - Delta_(a b) ], $
     where $Delta_(a b) = (w_a w_b) / (w_a + w_b) norm(bold(p)_a - bold(p)_b)_2^2$ is the exact Ward linkage divergence in probability-vector space. Consequently, merging strictly improves Category Utility if and only if $Delta_(a b) < "CU"(cal(P))$.
  2. Splitting a cluster $C_a in cal(P)$ into its $r >= 2$ subclusters ${C_(a, 1), dots, C_(a, r)}$ (having local Category Utility $"CU"_"loc"(C_a)$ inside $C_a$) yields partition $cal(P)_"sp"$ of size $K + r - 1$ with exact utility change
     $ "CU"(cal(P)_"sp") - "CU"(cal(P)) = (r w_a "CU"_"loc"(C_a) - (r - 1) "CU"(cal(P))) / (K + r - 1). $
     Consequently, splitting strictly improves Category Utility if and only if $w_a "CU"_"loc"(C_a) > ((r - 1) / r) "CU"(cal(P))$.
]

#proof[
  Let $U = sum_(k=1)^K w_k norm(bold(p)_k)_2^2 - norm(bold(p)_0)_2^2 = K "CU"(cal(P))$. Under a merge of $C_a$ and $C_b$, $w_(a b) = w_a + w_b$ and $bold(p)_(a b) = (w_a bold(p)_a + w_b bold(p)_b) / w_(a b)$. Multiplying the numerator loss $L_(a b) = w_a norm(bold(p)_a)_2^2 + w_b norm(bold(p)_b)_2^2 - w_(a b) norm(bold(p)_(a b))_2^2$ by $w_(a b)$ gives
  $ w_(a b) L_(a b) &= w_a (w_a + w_b) norm(bold(p)_a)_2^2 + w_b (w_a + w_b) norm(bold(p)_b)_2^2 \
    &quad - w_a^2 norm(bold(p)_a)_2^2 - w_b^2 norm(bold(p)_b)_2^2 - 2 w_a w_b chevron.l bold(p)_a, bold(p)_b chevron.r \
    &= w_a w_b norm(bold(p)_a - bold(p)_b)_2^2, $
  so $L_(a b) = Delta_(a b)$. Since the merged partition has $K - 1$ clusters and numerator $U' = U - Delta_(a b)$,
  $ "CU"(cal(P)_"mg") - "CU"(cal(P)) &= (K "CU"(cal(P)) - Delta_(a b)) / (K - 1) - "CU"(cal(P)) \
    &= ("CU"(cal(P)) - Delta_(a b)) / (K - 1), $
  proving Claim 1.

  For Claim 2, let $w_(a, l) = P(C_(a, l))$ so $sum_(l=1)^r w_(a, l) = w_a$ and $sum_(l=1)^r w_(a, l) bold(p)_(a, l) = w_a bold(p)_a$. Expanding around $bold(p)_a$ yields
  $ sum_(l=1)^r w_(a, l) norm(bold(p)_(a, l))_2^2 - w_a norm(bold(p)_a)_2^2 \
    = sum_(l=1)^r w_(a, l) norm(bold(p)_(a, l) - bold(p)_a)_2^2 = r w_a "CU"_"loc"(C_a), $
  using $P(C_(a, l) | C_a) = w_(a, l) / w_a$. Hence the split numerator is $U'' = K "CU"(cal(P)) + r w_a "CU"_"loc"(C_a)$ over $K + r - 1$ clusters, giving
  $ "CU"(cal(P)_"sp") - "CU"(cal(P)) &= (K "CU"(cal(P)) + r w_a "CU"_"loc"(C_a)) / (K + r - 1) - "CU"(cal(P)) \
    &= (r w_a "CU"_"loc"(C_a) - (r - 1) "CU"(cal(P))) / (K + r - 1), $
  which completes the proof.
]

Theorem 1 reveals a sharp structural duality: COBWEB's merge criterion is mathematically identical to Ward's minimum-variance criterion @ward1963hierarchical in $M$-dimensional probability space, thresholded adaptively against the current global $"CU"(cal(P))$.

#subsec("C", "Order Sensitivity Despite Bidirectional Restructuring")

Fisher @fisher1987cobweb introduced *Merge* and *Split* specifically to mitigate order sensitivity in incremental tree learning. Nevertheless, because *Merge* only evaluates the two best hosts for the *current* instance rather than all $binom(K, 2)$ sibling pairs, COBWEB remains strictly order-dependent even on tiny symmetric datasets.

#proposition(title: "Minimal Finite-Sample Counterexample to Order Independence", origin: "Original Counterexample")[
  Consider the $n = 4$, $d = 4$ binary dataset $cal(D) = {x_1, x_2, x_3, x_4} subset {0, 1}^4$ given by
  $ x_1 = (0, 0, 0, 0), quad x_2 = (0, 0, 0, 1), \
    x_3 = (1, 1, 1, 0), quad x_4 = (1, 1, 1, 1). $
  Under presentation order $pi_A = (x_1, x_3, x_2, x_4)$, COBWEB recovers the global Category Utility optimum $cal(P)^* = {{x_1, x_2}, {x_3, x_4}}$ with $"CU"(cal(P)^*) = 3 / 4$. Under presentation order $pi_B = (x_1, x_2, x_3, x_4)$, COBWEB converges to the strictly suboptimal three-cluster partition $cal(P)_B = {{x_1}, {x_2}, {x_3, x_4}}$ with $"CU"(cal(P)_B) = 7 / 12 < 3 / 4$.
]

#proof[
  On the full dataset $cal(D)$, every attribute has balanced marginal $(1/2, 1/2)$, so $norm(bold(p)_0)_2^2 = 4(1/4 + 1/4) = 2$. In $cal(P)^*$, both clusters agree on $A_1, A_2, A_3$ and split $A_4$ evenly, yielding $norm(bold(p)_k)_2^2 = 3 + 1/2 = 7/2$ and $"CU"(cal(P)^*) = (7/2 - 2)/2 = 3/4$. Any 4-singleton partition gives $(4 - 2)/4 = 1/2$, any 3-cluster partition with one matching pair gives $((1/2)(7/2) + 1 + 1 - 2)/3 = 7/12$, and any 3+1 partition gives $((3/4)(20/9) + 1 - 2)/2 = 1/3$; hence $cal(P)^*$ is the unique global optimum.

  *Trace of $pi_A = (x_1, x_3, x_2, x_4)$*: After $(x_1, x_3)$, the root holds leaves $C_1 = {x_1}, C_2 = {x_3}$. On $x_2$, each root attribute has counts $(2, 1)$, so $norm(bold(p)_0)_2^2 = 20/9$. Incorporating $x_2$ into $C_1$ yields $"CU" = ((2/3)(7/2) + 4/3 - 20/9)/2 = 13/18$, beating Incorporate into $C_2$ ($2/9$), Create ($16/27$), and Merge ($0$), forming ${{x_1, x_2}, {x_3}}$. On $x_4$, incorporating into $C_2 = {x_3}$ attains $"CU" = 3/4$, strictly beating Create ($7/12$), Incorporate into $C_1$ ($1/3$), Split of $C_1$ ($7/12$), and Merge ($0$), recovering $cal(P)^*$.

  *Trace of $pi_B = (x_1, x_2, x_3, x_4)$*: After $(x_1, x_2)$, the root holds leaves $C_1 = {x_1}, C_2 = {x_2}$. On $x_3$ ($norm(bold(p)_0)_2^2 = 20/9$), incorporating into $C_1$ agrees only on $A_4 = 0$ ($"CU" = 7/18$) and into $C_2$ gives $2/9$, so Create wins with $"CU" = 16/27$, forming ${{x_1}, {x_2}, {x_3}}$. On $x_4$ ($norm(bold(p)_0)_2^2 = 2$), the three host incorporations score $7/12$ ($C_3$, best host $k_1 = 3$), $5/12$ ($C_2$, agreeing on $A_4 = 1$, second-best host $k_2 = 2$), and $1/3$ ($C_1$). Create yields $1/2$, leaf $C_3$ cannot split, and merging the *two best hosts* $(C_3, C_2)$ with $x_4$ gives ${{x_1}, {x_2, x_3, x_4}}$ with $"CU" = 1/3$. Because COBWEB never evaluates merging non-hosts $(C_1, C_2)$ while assigning $x_4$ to $C_3$, Incorporate into $C_3$ wins, yielding $cal(P)_B$ with $"CU" = 7/12 < 3/4$.
]

#subsec("D", "Incremental Complexity and Cached Sufficient-Statistic Bounds")

#theorem(title: "Per-Instance and Stream Complexity of COBWEB", origin: "Classical — Cited")[
  Let $cal(T)$ have height $h$, maximum branching factor $B$, and maximum attribute cardinality $V = max_i |cal(V)_i|$. Following Fisher @fisher1987cobweb and Gennari et al. @gennari1989classit, naively re-evaluating $"CU"$ from scratch for every candidate operator requires $cal(O)(B^2 d V h)$ time per instance, or $cal(O)(n B^2 d V log_B n)$ over a balanced stream of $n$ instances. Moreover, caching the per-child squared probability norms $q_k = norm(bold(p)_k)_2^2$ and their weighted sum $W = sum_(k=1)^K w_k q_k$ reduces per-instance cost to $cal(O)(B d V h)$ while producing identical decisions.
]

#proof[
  At a single internal node with $K <= B$ children, naively computing $"CU"$ for one partition sums $sum_(i=1)^d m_i <= d V$ squared probabilities across all $K$ children in $cal(O)(B d V)$ time. Evaluating *Incorporate* across all $K$ hosts requires $K$ evaluations, costing $cal(O)(B^2 d V)$. *Create*, *Merge*, and *Split* add at most $cal(O)(B d V)$ work (with *Split* bounded by $B$ promoted children). Descending at most $h$ levels yields $cal(O)(B^2 d V h)$ per instance; when $cal(T)$ remains balanced with $h = cal(O)(log_B n)$, summing over $n$ insertions gives $cal(O)(n B^2 d V log_B n)$ (and $cal(O)(n^2 B^2 d V)$ in the adversarial degenerate-chain worst case).

  With caching, updating a single child $C_k$ by adding $x$ modifies only $d$ of the $M <= d V$ frequency entries of $C_k$. Recomputing $q_k'$ from $q_k$ and $n_(C_k)$ requires only $cal(O)(d)$ arithmetic operations after a single $cal(O)(d V)$ pass to rescale $q_k$ when $n_(C_k)$ increments by $1$ — in fact, storing integer sum-of-squares $Q_k = sum_(i=1)^d sum_(j=1)^(m_i) n_(C_k, i, j)^2$ updates exactly via
  $ Q_k' = Q_k + sum_(i=1)^d (2 n_(C_k, i, x_i) + 1) $
  in $cal(O)(d)$ time! Then $w_k' q_k' = Q_k' / (n_C (n_(C_k) + 1))$, while unchanged siblings scale trivially by $(n_C - 1) / n_C$. Evaluating all $K <= B$ incorporations therefore costs $cal(O)(B d)$ (and $cal(O)(d V)$ for a merge table addition), bounding each level by $cal(O)(B d + d V) subset.eq cal(O)(B d V)$ and total per-instance time by $cal(O)(B d V h)$.
]

#subsec("E", "Implicit Conditional-Independence Inductive Bias")

#proposition(title: "Blindness to Higher-Order Attribute Interactions", origin: "Original Proof")[
  Category Utility is completely invariant to within-cluster statistical dependence among attributes. Specifically, there exist two-cluster joint distributions over binary attributes $(A_1, A_2)$ with maximal cluster-attribute mutual information $I(C; A_1, A_2) = log 2$ bits for which $"CU"({C_1, C_2}) = 0$.
]

#proof[
  Let $C_1 = {(0, 0), (1, 1)}$ (even parity) and $C_2 = {(0, 1), (1, 0)}$ (odd parity) with equal weights $P(C_1) = P(C_2) = 1/2$. Given $(A_1, A_2)$, the cluster label $C$ is deterministically $A_1 xor A_2$, so $H(C | A_1, A_2) = 0$ and $I(C; A_1, A_2) = H(C) = log 2$. Yet for both $k in {1, 2}$ and $i in {1, 2}$, every univariate conditional marginal satisfies $P(A_i = 0 | C_k) = P(A_i = 1 | C_k) = 1/2 = P(A_i = 0)$. By Lemma 1, $"CU"({C_1, C_2}) = 0$, rendering the clusters mathematically indistinguishable from a random partition under COBWEB.
]

#sec("III", "Complexity, Geometry and Operational Profile")

COBWEB operates strictly in the discrete probability-simplex geometry of attribute marginals rather than Euclidean space. Its memory footprint scales as $cal(O)(|cal(T)| d V)$ for storing frequency tables across $|cal(T)| <= 2 n - 1$ tree nodes, and because every node summarizes its descendant stream purely via additive integer counts, COBWEB is genuinely single-pass and streaming-capable without retaining raw instances. Using the integer sum-of-squares trick in Theorem 2 reduces per-level evaluation from $cal(O)(B^2 d V)$ to $cal(O)(B d + d V)$, substantially accelerating wide-branching trees and eliminating floating-point drift.

However, COBWEB exhibits two operational bottlenecks: (i) in high-cardinality nominal domains ($V >> 10$), sparse contingency tables inflate memory and slow merge operations unless hash-backed sparse maps are used; (ii) because trees are never globally rebalanced, sorted or drifting input streams can induce skewed trees of height $h = Omega(n)$ and trigger excessive *Split* / *Merge* thrashing.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Discrete nominal attributes; continuous features require prior discretization],
    [Objective Function], [Local parent-child Category Utility $"CU"(cal(P))$ (amortized Gini reduction)],
    [Time Complexity], [$cal(O)(n B^2 d V h)$ naive; $cal(O)(n B d V h)$ with integer sum-of-squares caching],
    [Space Complexity], [$cal(O)(|cal(T)| d V)$ with $|cal(T)| <= 2 n - 1$; no raw-instance retention required],
    [Number of Clusters $K$], [Determined dynamically at every internal tree level via $1 / K$ amortization],
    [Invariances], [Invariant to attribute-value relabeling and attribute ordering; *not* order-invariant in $n$],
    [Missing Values], [Handled natively by excluding missing attributes from node count denominators],
    [Inductive Bias], [Axis-aligned conditional independence (Naive Bayes marginal factorization)],
  ),
  caption: [Operational and geometric profile of COBWEB.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Cognitive Modeling and Human Concept Formation:* COBWEB was directly motivated by Rosch's basic-level category theory and Gluck-Corter human learning experiments @fisher1987cobweb. The $1 / K$ normalization naturally produces intermediate "basic-level" concepts that maximize predictive cue validity without collapsing into overly general roots or rote leaf memorization, explaining human typicality gradients and incremental category acquisition.

*2. Electronic Health Record (EHR) Phenotyping:* Clinical registries record predominantly discrete diagnostic codes (ICD-10), medication classes, and categorical comorbidity flags with pervasive missingness. COBWEB's native handling of nominal tables and skipped missing attributes constructs interpretable clinical taxonomies where each internal node specifies explicit conditional symptom probabilities $P(A_i = v_(i, j) | C_k)$, supporting differential diagnosis and missing-feature prediction.

*3. Software Defect and IT Incident Taxonomy Induction:* Continuous integration pipelines emit streams of structured categorical fault tickets (subsystem, compiler flag, error code, severity, OS target). Because incidents arrive online and new failure families emerge over time, COBWEB's *Create* and *Split* operators dynamically spawn novel bug taxonomies while *Merge* consolidates redundant historical categories without batch re-training.

*4. E-Commerce Catalog and Clickstream Behavioral Typing:* Retail product schemas and user sessions consist of multi-valued categorical facets (brand tier, department, device class, payment mode). COBWEB builds navigable hierarchical facet trees directly from streaming logs, enabling cold-start attribute imputation by descending the tree along observed facets and reading off $arg max_j P(A_i = v_(i, j) | C)$ at the matched leaf.

*5. Network Intrusion Signature Grouping:* Packet-level firewall and DNS alerts exhibit categorical protocol, flag, port-class, and rule-category signatures. COBWEB incrementally partitions alert streams into interpretable behavioral concepts, flagging zero-day signatures whenever the *Create* operator isolates a low-support anomalous branch at high tree depth.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

COBWEB degrades severely under four mathematical regimes: (i) *continuous metric data*, where pre-discretization introduces arbitrary bin-edge artifacts and destroys inter-bin distance ordering; (ii) *adversarially ordered or concept-drifting streams*, where Proposition 1 shows local *Merge* / *Split* cannot escape suboptimal early splits; (iii) *high-dimensional noisy attributes*, since $"CU"$ sums linearly over all $d$ attributes and uniform noise dilutes signal gain without feature selection; and (iv) *XOR or strong pairwise interaction structures*, where Proposition 2 proves $"CU"$ vanishes identically.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Categorical EHR / clinical codes], [High: native nominal tables, missingness tolerance, interpretable posteriors], [Prune deep leaves; compare AutoClass (Alg. 87)],
    [Streaming symbolic logs / tickets], [High: single-pass online taxonomy with dynamic $K$], [Shuffle initial warm-up batch to curb order bias],
    [Continuous sensor / metric data], [Poor: discretization loses ordinality and covariance], [Use CLASSIT (Alg. 86), AutoClass, or TMM],
    [High-D sparse text / genomics], [Low: $cal(O)(d V)$ table blowup and additive noise dilution], [Use spherical $k$-means, LDA, or sparse NMF],
    [Relational / parity interactions], [Unsuitable: marginal $"CU"$ blind to joint interactions], [Use full-covariance GMM or kernel methods],
  ),
  caption: [Domain suitability matrix and failure mitigations for COBWEB.],
)

#sec("V", "Structural Relations and Empirical Behavior")

COBWEB sits at the intersection of probabilistic mixture modeling and online hierarchical clustering.

#proposition(title: "Structural Embedding of COBWEB in Categorical Mixture and Tree Families", origin: "Original Derivation")[
  Under a fixed partition size $K$, maximizing COBWEB's Category Utility over hard partitions $cal(P)$ is equivalent to maximizing the Tsallis-$2$ (quadratic) mutual information $I_2 (C; A_1, dots, A_d) = sum_(i=1)^d (G_i (C) - sum_k P(C_k) G_i (C_k))$ under an independent-attribute Naive Bayes model, whereas maximizing the hard-assignment log-likelihood of Latent Class Analysis @lazarsfeld1968latent maximizes the Shannon mutual information $I_1 (C; A_1, dots, A_d) = sum_(i=1)^d (H_i (C) - sum_k P(C_k) H_i (C_k))$.
]

#proof[
  In Latent Class Analysis @lazarsfeld1968latent with hard assignments $cal(P)$, the maximum-likelihood empirical log-likelihood of $cal(D)$ is $sum_(k=1)^K n_(C_k) sum_(i=1)^d sum_(j=1)^(m_i) P(A_i = v_(i, j) | C_k) log P(A_i = v_(i, j) | C_k) = -n sum_(i=1)^d sum_(k=1)^K P(C_k) H_i (C_k)$. Subtracting the constant one-cluster baseline $-n sum_(i=1)^d H_i (C)$ and dividing by $n$ yields precisely $I_1 (C; A_1, dots, A_d)$. Replacing Shannon entropy $H_i (C) = -sum_j p_(i j) log p_(i j)$ by the quadratic Tsallis-$2$ entropy $G_i (C) = 1 - sum_j p_(i j)^2$ and dividing by $K$ yields Lemma 1's exact formula for $"CU"(cal(P))$.
]

This linkage clarifies COBWEB's place in the clustering catalog: CLASSIT @gennari1989classit (Alg. 86) replaces COBWEB's discrete quadratic probabilities with continuous $L^2$-integrated Gaussian densities; AutoClass @cheeseman1996autoclass (Alg. 87) replaces greedy online quadratic utility with global Bayesian MAP-EM under Shannon log-likelihood and conjugate priors; $k$-Modes @huang1998extensions and ROCK @guha2000rock offer batch categorical alternatives; and modern online hierarchical algorithms such as BIRCH @zhang1996birch, PERCH @kobren2017perch, and GRINCH @monath2019grinch generalize COBWEB's local tree operators with rotation and grafting moves that improve theoretical cost guarantees @dasgupta2016cost.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized continuous evaluation suite (requiring nominal discretization for COBWEB), COBWEB records Iris Silhouette $= 0.021$ and $"ARI" = 0.066$ @hubert1985comparing @rousseeuw1987silhouettes, Customers Silhouette $= 0.312$, Multi-Shape Silhouette $= 0.501$ ($"ARI" =$ N/A), and High-D Silhouette $= 0.205$ ($"ARI" =$ N/A), with a mean runtime of $11.350$ ms. These measurements corroborate our theoretical analysis directly:
1. On *Iris*, severe degradation ($"Sil" = 0.021$, $"ARI" = 0.066$) stems from the compounding effects of metric-destroying discretization, conditional-independence blindness to strong petal-length / petal-width correlation (Proposition 2), and incremental order sensitivity (Proposition 1) that fragments the overlapping _versicolor_ and _virginica_ species across bin boundaries.
2. On *Multi-Shape*, well-separated clusters yield moderate axis-aligned bin purity ($"Sil" = 0.501$), whereas on *Customers* ($"Sil" = 0.312$) and *High-D* ($"Sil" = 0.205$), additive unweighted Gini aggregation dilutes discriminative axes across noisy dimensions.
3. The $11.350$ ms runtime—nearly $4.6 times$ slower than CLASSIT's $2.460$ ms on the same datasets—reflects the $cal(O)(d V)$ discrete contingency-table overhead established in Theorem 2.
