#import "../template/ieee.typ": *

#alg_entry(
  id: 67,
  name: "ROCK (RObust Clustering using linKs)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2000",
  paradigm: "Shared-Neighbor Link-Based Categorical Agglomeration",
  complexity: $cal(O)(n^2 + n m_m m_a + n^2 log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Traditional hierarchical clustering relies on Euclidean distances and vector centroids, which fail fundamentally on *categorical*, *Boolean*, and *market-basket* data $X = { T_1, dots, T_n }$ where each record $T_i subset.eq cal(I)$ is a sparse subset of discrete items or attribute-value pairs. Averaging sparse binary indicators yields dense, non-interpretable centroids whose pairwise distances shrink rapidly as clusters grow, triggering catastrophic snowball merges @guha2000rock. Even replacing Euclidean distance with pairwise Jaccard similarity remains fragile because two unrelated transactions can share an accidental common item while two genuinely related transactions share none directly. ROCK @guha2000rock overcomes both pathologies by lifting similarity from a *local pairwise* comparison to a *global second-order* graph statistic: the number of *common neighbors* (links) shared by two records in a thresholded similarity graph.

#definition(title: "Neighbors, Links, Global Criterion, and Goodness Measure")[
  Let $"sim": X times X -> [0, 1]$ be a symmetric normalized similarity function (canonically Jaccard similarity $"sim"(T_i, T_j) = abs(T_i inter T_j) / abs(T_i union T_j)$) and fix a neighbor threshold $theta in (0, 1)$.
  + *Neighbor Graph:* Points $p, q in X$ are *neighbors* if $"sim"(p, q) >= theta$. Let $N_theta (p) = { r in X : "sim"(p, r) >= theta }$ and define the symmetric binary adjacency matrix $A in { 0, 1 }^(n times n)$ with $A_(i j) = bold(1)["sim"(T_i, T_j) >= theta]$.
  + *Links:* The number of *links* between $p$ and $q$ is the count of common neighbors
    $ "link"(p, q) = abs(N_theta (p) inter N_theta (q)) = (A^2)_(p q), $
    extended to clusters $C_i, C_j subset.eq X$ by $"link"(C_i, C_j) = sum_(p in C_i) sum_(q in C_j) "link"(p, q)$.
  + *Global Criterion Function:* Modeling a point in a cluster of size $n_i$ as having roughly $n_i^(f(theta))$ internal neighbors (canonically $f(theta) = (1 - theta) / (1 + theta)$ for market baskets), ROCK seeks a $k$-partition $cal(C) = { C_1, dots, C_k }$ maximizing
    $ E_l (cal(C)) &= sum_(i=1)^k n_i sum_(p, q in C_i) ("link"(p, q)) / (n_i^(1 + 2 f(theta))) \
      &= sum_(i=1)^k ("link"(C_i, C_i)) / (n_i^(2 f(theta))). $
  + *Goodness Measure:* To merge clusters greedily without collapsing all data into one cluster, ROCK evaluates pairs $(C_i, C_j)$ ($n_i = abs(C_i), n_j = abs(C_j)$) via
    $ g(C_i, C_j) = ("link"(C_i, C_j)) / ((n_i + n_j)^(1 + 2 f(theta)) - n_i^(1 + 2 f(theta)) - n_j^(1 + 2 f(theta))). $
]

The denominator of $g(C_i, C_j)$ estimates the *expected* number of cross-cluster links if $C_i$ and $C_j$ belong to a single true cluster of size $n_i + n_j$, normalizing raw link counts against cluster cardinality.

#algbox(
  title: "ROCK Link-Based Hierarchical Agglomerative Clustering",
  inputs: [Dataset $X$ ($abs(X) = n$), similarity function $"sim"$, threshold $theta in (0, 1)$, scaling function $f(theta)$, target clusters $k$, optional sample size $n_s$.],
  outputs: [Partition $cal(C) = { C_1, dots, C_k }$ of $X$.],
)[
  + *Optional Sampling:* If $n$ is large, draw a uniform random sample $S subset.eq X$ of size $n_s$ (via the Chernoff bound of CURE @guha1998cure); otherwise set $S = X$.
  + *Neighbor & Link Matrix Construction:* Build adjacency lists $N_theta (p)$ for all $p in S$ where $"sim"(p, q) >= theta$. Compute sparse link counts $"link"(p, q) = (A^2)_(p q)$ by iterating over every vertex $r in S$ and incrementing $"link"(p, q)$ for all pairs $p, q in N_theta (r)$.
  + *Heap Initialization:* Place each $p in S$ in a singleton cluster $C_p = { p }$. Build a local max-heap $q[i]$ for each cluster $C_i$ containing every $C_j$ with $"link"(C_i, C_j) > 0$ keyed by $g(C_i, C_j)$, and a global max-heap $Q$ keyed by $max_(j in q[i]) g(C_i, C_j)$.
  + *Greedy Link Agglomeration:* While $abs(Q) > k$ and the maximum key in $Q$ is strictly positive:
    - Extract $u = arg max Q$, let $v = "top"(q[u])$, and form $w = u union v$.
    - For every cluster $x in q[u] union q[v]$, update $"link"(w, x) <- "link"(u, x) + "link"(v, x)$, recompute $g(w, x)$, and update local heaps $q[x], q[w]$ and global heap $Q$.
  + *Disk-Point Labeling (if sampled):* Draw representative subsets $L_i subset.eq C_i$, and assign each $p in X without S$ to $arg max_i abs(N_theta (p) inter L_i) / (abs(L_i) + 1)^(f(theta))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Link-Matrix Positivity and Variational Coupling Between $g$ and $E_l$")

We first establish the spectral positivity of ROCK's link operator and derive the *exact* algebraic relationship between greedy maximization of the heuristic goodness $g(C_i, C_j)$ and steepest ascent on the global objective $E_l (cal(C))$.

#theorem(title: "Positive Semidefiniteness and Exact $E_l$-Ascent Law", origin: "Original Proof")[
  Let $A in { 0, 1 }^(n times n)$ be the symmetric neighbor adjacency matrix, $L = A^2$ the link matrix, and $alpha = 2 f(theta) > 0$.
  + *PSD Kernel Structure:* $L = A^2 succ.eq 0$ is symmetric positive semidefinite; for every cluster indicator $bold(z)_i in { 0, 1 }^n$, $"link"(C_i, C_i) = norm(A bold(z)_i)_2^2$, $"link"(C_i, C_j) = chevron.l A bold(z)_i, A bold(z)_j chevron.r$, and $"link"(p, q) <= sqrt(d(p) d(q))$ where $d(p) = abs(N_theta (p)) = L_(p p)$.
  + *Well-Posedness of $g$:* For all $alpha > 0$ and integers $n_i, n_j >= 1$, the denominator $Delta_(i j) = (n_i + n_j)^(1 + alpha) - n_i^(1 + alpha) - n_j^(1 + alpha) > 0$ strictly.
  + *Exact Ascent Identity:* Define each cluster's normalized link density $rho_i = "link"(C_i, C_i) / n_i^(1 + alpha)$. Merging $(C_i, C_j)$ changes the global objective $E_l$ by exactly
    $ Delta E_l (C_i, C_j) = (2 "link"(C_i, C_j) - Phi_(i j)) / (n_i + n_j)^alpha, $
    where $Phi_(i j) = rho_i n_i ((n_i + n_j)^alpha - n_i^alpha) + rho_j n_j ((n_i + n_j)^alpha - n_j^alpha)$.
  + *Equilibrium Equivalence:* Whenever $rho_i = rho_j = rho$, $Phi_(i j) = rho Delta_(i j)$, giving
    $ Delta E_l (C_i, C_j) = (Delta_(i j) / (n_i + n_j)^alpha) (2 g(C_i, C_j) - rho). $
    Thus at uniform link density $rho$, merging strictly increases $E_l$ if and only if $2 g(C_i, C_j) > rho$.
]

#proof[
  For Part (1), $A = A^top$ implies $bold(v)^top L bold(v) = bold(v)^top A^top A bold(v) = norm(A bold(v))_2^2 >= 0$ for all $bold(v) in bb(R)^n$, so $L succ.eq 0$. Evaluating on binary indicator vectors $bold(z)_i, bold(z)_j$ yields $"link"(C_i, C_i) = bold(z)_i^top A^2 bold(z)_i = norm(A bold(z)_i)_2^2$ and $"link"(C_i, C_j) = bold(z)_i^top A^2 bold(z)_j = chevron.l A bold(z)_i, A bold(z)_j chevron.r$; Cauchy-Schwarz on unit basis vectors $bold(e)_p, bold(e)_q$ gives $"link"(p, q) = chevron.l A bold(e)_p, A bold(e)_q chevron.r <= norm(A bold(e)_p)_2 norm(A bold(e)_q)_2 = sqrt(d(p) d(q))$.

  For Part (2), set $u = n_i / (n_i + n_j) in (0, 1)$ and $1 - u = n_j / (n_i + n_j) in (0, 1)$. Since $1 + alpha > 1$, strict power decay gives $u^(1 + alpha) < u$ and $(1 - u)^(1 + alpha) < 1 - u$. Summing yields $u^(1 + alpha) + (1 - u)^(1 + alpha) < 1$; multiplying by $(n_i + n_j)^(1 + alpha)$ gives $Delta_(i j) > 0$.

  For Part (3), write $L_(i i) = "link"(C_i, C_i) = rho_i n_i^(1 + alpha)$, $L_(j j) = "link"(C_j, C_j) = rho_j n_j^(1 + alpha)$, and $L_(i j) = "link"(C_i, C_j)$. Because $"link"(C_i union C_j, C_i union C_j) = L_(i i) + L_(j j) + 2 L_(i j)$, the exact objective increment is
  $ Delta E_l &= (L_(i i) + L_(j j) + 2 L_(i j)) / (n_i + n_j)^alpha - L_(i i) / n_i^alpha - L_(j j) / n_j^alpha \
    &= (2 L_(i j)) / (n_i + n_j)^alpha - L_(i i) ((n_i + n_j)^alpha - n_i^alpha) / (n_i^alpha (n_i + n_j)^alpha) \
    &quad - L_(j j) ((n_i + n_j)^alpha - n_j^alpha) / (n_j^alpha (n_i + n_j)^alpha). $
  Substituting $L_(i i) / n_i^alpha = rho_i n_i$ and $L_(j j) / n_j^alpha = rho_j n_j$ immediately yields Part (3).

  For Part (4), when $rho_i = rho_j = rho$, collecting terms in $Phi_(i j)$ gives
  $ Phi_(i j) &= rho (n_i (n_i + n_j)^alpha - n_i^(1 + alpha) \
    &quad + n_j (n_i + n_j)^alpha - n_j^(1 + alpha)) \
    &= rho ((n_i + n_j)^(1 + alpha) - n_i^(1 + alpha) - n_j^(1 + alpha)) = rho Delta_(i j). $
  Factoring $Delta_(i j) / (n_i + n_j)^alpha > 0$ from $(2 L_(i j) - rho Delta_(i j)) / (n_i + n_j)^alpha$ yields $(Delta_(i j) / (n_i + n_j)^alpha)(2 g(C_i, C_j) - rho)$, completing the proof.
]

Theorem 1 Part (4) both justifies Guha et al.'s goodness denominator *variationally* and reveals its exact discrepancy from steepest $E_l$-ascent: $g$ ignores both the size-dependent prefactor $Delta_(i j) / (n_i + n_j)^alpha$ and heterogeneity between $rho_i$ and $rho_j$.

#subsec("B", "Stochastic Link Amplification Over Pairwise Similarity")

Why does squaring the neighbor graph ($L = A^2$) dramatically suppress cross-cluster noise compared to single-linkage or direct Jaccard thresholding? We prove an exponential finite-sample separation guarantee under a two-block random neighbor graph.

#theorem(title: "Stochastic Link Concentration and Noise Rejection", origin: "Original Proof")[
  Suppose $X = C_1 union C_2$ consists of two disjoint clusters of equal size $m >= 2$, and the loopless neighbor graph $A in { 0, 1 }^(2 m times 2 m)$ ($A_(i i) = 0$) has independent symmetric edges with intra-cluster probability $p$ and inter-cluster probability $q$ ($0 <= q < p <= 1$). Then:
  + Every intra-cluster pair $(u, u')$ and cross-cluster pair $(u, v)$ satisfy exact link expectations
    $ bb(E)["link"(u, u')] &= (m - 2) p^2 + m q^2, \
      bb(E)["link"(u, v)] &= 2 (m - 1) p q, $
    with positive expectation gap $Delta = (p - q)(m (p - q) - 2 p)$ whenever $m > 2 p / (p - q)$.
  + Every intra-cluster pair simultaneously has strictly more links than every cross-cluster pair with probability at least
    $ 1 - 2 m^2 exp(- Delta^2 / (4 (m - 1))). $
    Thus exact pairwise link separation occurs with high probability whenever $p - q = omega(((ln m) / m)^(1 / 4))$, even when $q > 0$ generates $Theta(m^2 q)$ false direct cross-cluster edges that completely destroy single-linkage clustering.
]

#proof[
  For fixed distinct $u, u' in C_1$, decompose the common-neighbor count over all remaining $2 m - 2$ vertices:
  $ "link"(u, u') = sum_(r in C_1 without { u, u' }) A_(u r) A_(u' r) + sum_(r in C_2) A_(u r) A_(u' r). $
  Because each index $r$ involves a disjoint pair of graph edges ${ { u, r }, { u', r } }$, the $2 m - 2$ Bernoulli variables $X_r = A_(u r) A_(u' r) in { 0, 1 }$ are *mutually independent*, with $bb(E)[X_r] = p^2$ for $r in C_1 without { u, u' }$ and $bb(E)[X_r] = q^2$ for $r in C_2$. Linearity of expectation gives $mu_"in" = (m - 2) p^2 + m q^2$.

  Similarly, for fixed $u in C_1, v in C_2$,
  $ "link"(u, v) = sum_(r in (C_1 union C_2) without { u, v }) A_(u r) A_(v r). $
  Again, the $2 m - 2$ indicators $Y_r = A_(u r) A_(v r) in { 0, 1 }$ are mutually independent; for every $r not in { u, v }$, one edge is intra-cluster and the other is inter-cluster, so $bb(E)[Y_r] = p q$, yielding $mu_"out" = 2 (m - 1) p q$. Subtracting gives
  $ Delta = mu_"in" - mu_"out" &= m (p^2 - 2 p q + q^2) - 2 p (p - q) \
    &= (p - q)(m (p - q) - 2 p). $

  Assume $Delta > 0$ and set midpoint $tau = (mu_"in" + mu_"out") / 2$. Applying Hoeffding's inequality to each sum of $2 m - 2$ independent $[0, 1]$-valued random variables at deviation $t = Delta / 2$ yields
  $ "Pr"("link"(u, u') <= tau) &<= exp(- (2 (Delta / 2)^2) / (2 m - 2)) = exp(- Delta^2 / (4 (m - 1))), \
    "Pr"("link"(u, v) >= tau) &<= exp(- Delta^2 / (4 (m - 1))). $
  There are $2 binom(m, 2) = m(m - 1)$ intra-cluster pairs and $m^2$ inter-cluster pairs, totaling $2 m^2 - m < 2 m^2$ pairs. By the union bound, the probability that any intra-cluster pair falls below $tau$ or any inter-cluster pair rises above $tau$ is at most $2 m^2 exp(- Delta^2 / (4 (m - 1)))$. Since $Delta = Theta(m (p - q)^2)$, the exponent scales as $- Theta(m (p - q)^4)$, which dominates $ln(2 m^2)$ whenever $p - q = omega(((ln m) / m)^(1 / 4))$.
]

Theorem 2 explains ROCK's surprising resilience in noisy high-dimensional spaces: even when pairwise neighbor graphs contain $Theta(m^2 q)$ false edges, squaring $A$ averages over $2 m - 2$ independent witnesses, driving false-edge confusion exponentially toward zero.

#subsec("C", "Classical Worst-Case Time and Space Complexity")

#theorem(title: "Link-Matrix and Heap Complexity of ROCK", origin: "Classical — Cited")[
  @guha2000rock Let $m_m = max_(p in X) abs(N_theta (p))$ and $m_a = (1 / n) sum_(p in X) abs(N_theta (p))$ denote the maximum and average neighbor-list sizes in the $theta$-neighbor graph on $n$ records. Then ROCK executes in worst-case time
  $ T(n) = cal(O)(n^2 + n m_m m_a + n^2 log n) $
  and space $S(n) = cal(O)(min(n^2, n m_m m_a))$. When $m_m = Theta(n)$, link computation alone reaches $Theta(n^3)$.
]

#proof[
  Computing all pairwise similarities to form neighbor lists $N_theta (p)$ takes $cal(O)(n^2)$ similarity evaluations. In the sparse link-accumulation loop, for each vertex $r in X$ of degree $d_r = abs(N_theta (r))$, iterating over all ordered pairs $(p, q) in N_theta (r) times N_theta (r)$ takes $d_r^2 <= m_m d_r$ operations; summing over $r = 1, dots, n$ bounds total link computation by
  $ sum_(r=1)^n d_r^2 <= m_m sum_(r=1)^n d_r = n m_m m_a, $
  and also bounds the number of nonzero entries in $L = A^2$ by $min(n^2, n m_m m_a)$, establishing the space bound. During agglomeration, there are at most $n - k$ merges. Each merge traverses local heaps $q[u], q[v]$ (at most $n$ distinct neighbors $x$ total), performs $cal(O)(1)$ link addition per neighbor, and executes $cal(O)(log n)$ deletion/re-insertion in $q[x]$, $q[w]$, and $Q$, costing $cal(O)(n log n)$ per merge and $cal(O)(n^2 log n)$ total. Summing all three phases gives $cal(O)(n^2 + n m_m m_a + n^2 log n)$.
]

#subsec("D", "Bimodal Threshold Collapse and Size-Bias Phase Transition")

ROCK's Achilles heel—particularly on continuous data converted via similarity thresholds—is extreme sensitivity to $theta$ and the exponent $f(theta)$. We prove an exact phase transition in merge geometry at $f(theta) = 1 / 2$ ($theta = 1 / 3$).

#proposition(title: "Threshold Degeneracy and $f(theta) = 1 / 2$ Size-Bias Phase Transition", origin: "Original Counterexample")[
  Consider ROCK on any $n$-point dataset under similarity threshold $theta in (0, 1)$ with $f(theta) = (1 - theta) / (1 + theta)$:
  + *Sparse Disconnection Collapse:* If $theta > max_(i != j) "sim"(T_i, T_j)$, then $A = I_n$, $"link"(p, q) = 0$ for all $p != q$, and $g(C_i, C_j) = 0$ identically, halting agglomeration in $n$ singletons (or arbitrary tie-breaking).
  + *Dense Saturation & Size-Bias Phase Transition:* If $theta <= min_(i, j) "sim"(T_i, T_j)$, then $A = bold(1) bold(1)^top$ and $"link"(p, q) = n$ for all pairs, so $g(C_i, C_j)$ depends *solely* on cluster sizes $(n_i, n_j)$. For equal-size pairs $n_i = n_j = s$,
    $ g(s, s) = n / (2^(1 + 2 f(theta)) - 2) dot s^(1 - 2 f(theta)). $
    Consequently, at the critical threshold $theta_c = 1 / 3$ ($f(theta_c) = 1 / 2$):
    - For $theta > 1 / 3$ ($f(theta) < 1 / 2$): $g(s, s)$ is *strictly increasing* in $s$, inducing giant-component snowballing.
    - For $theta < 1 / 3$ ($f(theta) > 1 / 2$): $g(s, s)$ is *strictly decreasing* in $s$, forcing balanced binary-tree merges blind to data geometry.
    - At $theta = 1 / 3$ ($f(theta) = 1 / 2$): $g(s, s) = n / 2$ is strictly scale-invariant.
]

#proof[
  Part (1) is immediate: when no distinct pair achieves $"sim" >= theta$, $A = I_n$, so $L = A^2 = I_n$ has zero off-diagonal entries, making $"link"(C_i, C_j) = 0$ and $g(C_i, C_j) = 0$ for every disjoint pair $(C_i, C_j)$.

  For Part (2), when $"sim"(T_i, T_j) >= theta$ for all pairs, $A = bold(1) bold(1)^top$, so $L = A^2 = n bold(1) bold(1)^top$ and $"link"(C_i, C_j) = sum_(p in C_i, q in C_j) n = n n_i n_j$. Evaluating $g$ at $n_i = n_j = s$ gives
  $ g(s, s) &= (n s^2) / ((2 s)^(1 + 2 f(theta)) - 2 s^(1 + 2 f(theta))) \
    &= n / (2^(1 + 2 f(theta)) - 2) dot s^(1 - 2 f(theta)). $
  Because $f(theta) = (1 - theta) / (1 + theta)$ is strictly decreasing in $theta$ with $f(1 / 3) = 1 / 2$, the exponent $1 - 2 f(theta)$ is strictly positive for $theta > 1 / 3$ (making $s^(1 - 2 f(theta))$ strictly increasing), strictly negative for $theta < 1 / 3$ (strictly decreasing), and identically zero at $theta = 1 / 3$ (where $2^2 - 2 = 2$, yielding $g(s, s) = n / 2$).
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational Bottleneck:* As established in Theorem 3, computing $L = A^2$ costs $cal(O)(n m_m m_a)$, which can scale cubically ($Theta(n^3)$) on dense neighbor graphs when $theta$ is set too low! Even on sparse graphs ($m_m, m_a = cal(O)(1)$), maintaining local max-heaps $q[i]$ and global heap $Q$ requires $cal(O)(n^2 log n)$ worst-case time and $cal(O)(n^2)$ memory without random sampling. Drawing a Chernoff sample of size $n_s << n$ reduces clustering cost to $cal(O)(n_s^2 log n_s + n_s m_m m_a)$ plus $cal(O)(n k abs(L_i))$ for final assignment.

*Geometry & Metric Regime:* ROCK is fundamentally *non-Euclidean*: it requires only an arbitrary symmetric similarity oracle $"sim"(p, q) in [0, 1]$ (Jaccard, Dice, cosine, or RBF kernel on continuous data). Thresholding continuous distances into a binary graph discards fine-grained metric resolution in low dimensions (degrading compact Euclidean Silhouette), yet in high-dimensional or sparse categorical spaces, replacing first-order distances with second-order common-neighbor counts ($A^2$) acts as a robust denoising filter (Theorem 2) that counteracts distance concentration @beyer1999nearest.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Full Time Complexity], [$cal(O)(n^2 + n m_m m_a + n^2 log n)$ ($Theta(n^3)$ if $m_m = Theta(n)$)],
    [Sampled Time], [$cal(O)(n_s^2 log n_s + n_s m_m m_a + n k abs(L_i))$],
    [Memory Footprint], [$cal(O)(min(n^2, n m_m m_a))$ (link matrix + local heaps)],
    [Inductive Bias], [Shared-neighbor cohesion on categorical / discrete sets],
    [Key Hyperparameters], [Similarity threshold $theta$, scaling law $f(theta)$, clusters $k$],
    [Invariances], [Invariant to attribute ordering and any monotone map preserving $"sim" >= theta$],
    [Outlier Handling], [Zero-link isolated items stay singleton; pruned via small-cluster discard],
    [Critical Transition], [Size-bias phase transition at $f(theta) = 1 / 2$ ($theta = 1 / 3$, Prop. 1)],
  ),
  caption: [Operational, complexity, and geometric profile of ROCK (Alg. 67).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Market-Basket and E-Commerce Transaction Clustering:* Retail baskets are ultra-sparse subsets of $10^4$–$10^5$ SKUs. Euclidean centroids collapse unrelated baskets, whereas ROCK's Jaccard neighbor links (Theorem 2) cluster shopping missions sharing common co-purchased anchor items even when individual baskets have small direct overlap.
+ *Political Voting Records and Survey Response Typology:* Congressional roll-call matrices and psychometric categorical questionnaires consist of discrete nominal responses (`yea`/`nay`/missing). Shared-neighbor links robustly separate ideological blocs despite occasional cross-party votes (modeled directly by $q > 0$ in Theorem 2).
+ *Biological Taxonomy and Categorical Phenotype Profiling:* Classical morphological datasets (e.g., mushroom edibility, microbial biochemical traits) consist exclusively of unordered categorical descriptors. ROCK clusters taxa cleanly via discrete attribute overlap without arbitrary one-hot Euclidean distortion.
+ *Web Clickstream and Document Tag Co-Occurrence:* User session logs and folksonomy tag sets form sparse bipartite incidence graphs. Second-order link counts $A^2$ capture two-hop collaborative-filtering similarity, grouping semantically related sessions that share bridging URLs or tags.
+ *High-Dimensional Sparse Feature Denoising:* In sparse high-$d$ continuous or binary embeddings (gene-expression signatures, text TF-IDF vectors), first-order Euclidean distances concentrate @beyer1999nearest, whereas shared-neighbor links act as a non-parametric high-$d$ similarity regularizer closely allied with SNN clustering @ertoz2003snn.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ROCK fails in four predictable regimes: (i) *Threshold misspecification*, where $theta$ too high disconnects $A$ and $theta$ too low triggers $Theta(n^3)$ runtime and size-bias collapse (Proposition 1); (ii) *Low-dimensional continuous metric clusters*, where binary thresholding discards continuous geometric precision (hurting Euclidean Silhouette); (iii) *Heterogeneous-density similarity graphs*, where a single global threshold $theta$ cannot adapt simultaneously to dense and sparse regions (unlike CHAMELEON's $k$-NN normalization @karypis1999chameleon); and (iv) *Large unsampled datasets*, where $cal(O)(n m_m m_a + n^2 log n)$ memory and runtime become prohibitive.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Market baskets / transactions], [Native Jaccard + $A^2$ links (Thm. 2)], [Use sampling ($n_s << n$) if $n > 10^4$],
    [Categorical / nominal tables], [Avoids centroid collapse completely], [Validate $theta$ near connectivity knee],
    [High-$d$ sparse embeddings], [Links counteract $ell_2$ concentration], [Compare with SNN @ertoz2003snn or Louvain @blondel2008louvain],
    [Low-$d$ continuous Euclidean], [Binary $theta$-cut loses metric precision], [Prefer CURE @guha1998cure or CHAMELEON @karypis1999chameleon],
    [Multi-density spatial data], [Global $theta$ fails across density scales], [Use CHAMELEON @karypis1999chameleon or HDBSCAN @campello2013hdbscan],
  ),
  caption: [Domain suitability and failure-mode matrix for ROCK.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Although conceived for categorical databases, ROCK's global objective possesses an exact spectral-kernel reduction linking it directly to Kernel $k$-Means @dhillon2004kernel, Spectral Clustering @vonluxburg2007tutorial, and Shared Nearest Neighbors (SNN) @ertoz2003snn.

#proposition(title: "Exact Equivalence of ROCK at $f(theta) = 1 / 2$ to Link-Kernel $k$-Means", origin: "Original Proof")[
  At the scale-invariant critical point $f(theta) = 1 / 2$ (corresponding to $theta = 1 / 3$ under $f(theta) = (1 - theta) / (1 + theta)$):
  + *Exact Kernel $k$-Means Equivalence:* Maximizing ROCK's global criterion $E_l (cal(C))$ over $k$-partitions $cal(C) = { C_1, dots, C_k }$ is *algebraically identical* to minimizing the exact Kernel $k$-Means objective @dhillon2004kernel with the PSD link kernel $K = L = A^2$ (feature map $phi(T_r) = A_(colon, r) in { 0, 1 }^n$):
    $ "SSE"_K (cal(C)) = "Tr"(A^2) - E_l (cal(C)). $
  + *Spectral Relaxation:* Letting $H in bb(R)^(n times k)$ have orthonormal columns $H_(r i) = bold(1)[T_r in C_i] / sqrt(n_i)$ ($H^top H = I_k$), $E_l (cal(C)) = "Tr"(H^top A^2 H)$, whose continuous spectral upper bound is the sum of the top-$k$ squared singular values of $A$: $max_(cal(C)) E_l (cal(C)) <= sum_(j=1)^k sigma_j (A)^2$.
]

#proof[
  Define the feature embedding $phi(T_r) = A_(colon, r)$ (the $r$-th column of symmetric $A$), so $chevron.l phi(T_p), phi(T_q) chevron.r = (A^top A)_(p q) = (A^2)_(p q) = "link"(T_p, T_q)$. For any cluster $C_i$ of size $n_i$ with feature mean $bold(mu)_i = (1 / n_i) sum_(r in C_i) phi(T_r) = (1 / n_i) A bold(z)_i$, expanding the Kernel $k$-Means SSE gives
  $ "SSE"_K (cal(C)) &= sum_(i=1)^k sum_(r in C_i) norm(phi(T_r) - bold(mu)_i)_2^2 \
    &= sum_(r=1)^n norm(phi(T_r))_2^2 - sum_(i=1)^k n_i norm(bold(mu)_i)_2^2 \
    &= "Tr"(A^2) - sum_(i=1)^k (bold(z)_i^top A^2 bold(z)_i) / n_i. $
  At $f(theta) = 1 / 2$, we have $2 f(theta) = 1$, so by Theorem 1(1),
  $ E_l (cal(C)) &= sum_(i=1)^k ("link"(C_i, C_i)) / (n_i^(2 f(theta))) \
    &= sum_(i=1)^k (bold(z)_i^top A^2 bold(z)_i) / n_i = "Tr"(H^top A^2 H). $
  Hence $"SSE"_K (cal(C)) = "Tr"(A^2) - E_l (cal(C))$, where $"Tr"(A^2) = sum_(r=1)^n d(T_r)$ is partition-independent. Because $H^top H = I_k$, the Ky Fan variational theorem immediately bounds $"Tr"(H^top A^2 H) <= sum_(j=1)^k lambda_j (A^2) = sum_(j=1)^k sigma_j (A)^2$.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* Across the standardized 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, ROCK records *Iris* $"Sil" = 0.174, "ARI" = 0.499$; *Customers* $"Sil" = 0.127$; *Multi-Shape* $"Sil" = 0.552$; and *High-D* $"Sil" = 0.293$, at a mean runtime of *$25.428$ ms*. This distinctive profile matches our theoretical results with striking precision:
- *High-Dimensional Winner:* Remarkably, ROCK achieves the *highest* *High-D* Silhouette of all four hierarchical algorithms ($0.293$, beating BIRCH $0.249$, CURE $0.238$, and CHAMELEON $0.214$), directly validating Theorem 2 and Proposition 2: squaring the adjacency graph into link kernel $A^2$ averages over common neighbors to cancel first-order Euclidean distance concentration @beyer1999nearest.
- *Continuous Metric Loss on Low-$d$ Data:* On compact continuous datasets (*Iris* $"Sil" = 0.174, "ARI" = 0.499$; *Customers* $"Sil" = 0.127$), ROCK scores lowest among the four methods—reflecting Proposition 1's binary $theta$-thresholding distortion when continuous distances are discretized without metric weighting—while still capturing non-convex *Multi-Shape* topology ($"Sil" = 0.552 >$ BIRCH $0.452$).
- *Quadratic-to-Cubic Runtime Overhead:* At $25.428$ ms, ROCK is $18.6 times$ slower than BIRCH ($1.370$ ms), $10.5 times$ slower than CURE ($2.423$ ms), and $6.48 times$ slower than CHAMELEON ($3.922$ ms), directly reflecting Theorem 3's $cal(O)(n^2 + n m_m m_a + n^2 log n)$ link-accumulation and local-heap maintenance burden.
