#import "../template/ieee.typ": *

#alg_entry(
  id: 7,
  name: "K-Modes (Categorical Mode Partitioning)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1998",
  paradigm: "Alternating Hamming-Mode Minimization",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(A) = A_1 times A_2 times dots times A_d$ be a finite categorical product space, where each attribute $r in {1, dots, d}$ takes values in a finite unordered nominal alphabet $A_r = {a_(r, 1), dots, a_(r, m_r)}$ of cardinality $m_r = |A_r| >= 2$ (total category count $M_"cat" = sum_(r=1)^d m_r$). Given a dataset $cal(X) = {x_1, dots, x_n} subset.eq cal(A)$ with $x_i = (x_(i, 1), dots, x_(i, d))$, Euclidean arithmetic centroids and squared $L_2$ norms are mathematically undefined because nominal symbols possess neither ordering nor algebraic addition.

Huang @huang1998extensions equips $cal(A)$ with the *simple matching (Hamming) metric*:
$ d_H (x, q) = sum_(r=1)^d delta(x_r, q_r), quad delta(a, b) = cases(0 &"if" a = b, 1 &"if" a != b.) $
The *K-Modes clustering problem* seeks a disjoint $k$-partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ together with categorical prototype vectors (*modes*) $Q = (q_1, dots, q_k) in cal(A)^k$ minimizing the total within-cluster Hamming mismatch functional @huang1998extensions:
$ J_H (cal(C), Q) &= sum_(j=1)^k sum_(i in C_j) d_H (x_i, q_j) \
  &= sum_(j=1)^k sum_(i in C_j) sum_(r=1)^d delta(x_(i, r), q_(j, r)). $

#definition(title: "Coordinate Mode & Hamming–Voronoi Equilibrium")[
  For any non-empty subset $S subset.eq {1, dots, n}$ and attribute $r in {1, dots, d}$, let $f_r (a | S) = |{i in S : x_(i, r) = a}|$ denote the frequency of symbol $a in A_r$ in $S$. A vector $q in cal(A)$ is a *mode* of $S$ if
  $ f_r (q_r | S) = max_(a in A_r) f_r (a | S) quad forall r in {1, dots, d}. $
  A pair $(cal(C), Q)$ is a *Hamming–Voronoi mode equilibrium* if (i) every $i in C_j$ satisfies $d_H (x_i, q_j) <= d_H (x_i, q_ell)$ for all $ell != j$, and (ii) each $q_j$ is a mode of $C_j$.
]

Unlike Euclidean centroids, a mode $q_j$ is generally *non-unique* (due to frequency ties) and may be a *synthetic* attribute tuple not present in $cal(X)$. Huang's K-Modes algorithm solves the discrete non-convex optimization via alternating assignment and coordinate-wise frequency maximization (Algorithm 1).

#algbox(
  title: "Huang's K-Modes Algorithm",
  inputs: [Categorical dataset $cal(X) subset.eq cal(A)$ of size $n$, integer $k >= 1$, initial modes $Q^((0)) in cal(A)^k$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and modes $Q = (q_1, dots, q_k)$.],
)[
  1. Initialize $c_i^((-1)) <- 0$ for all $i$; set $t <- 0$.
  2. *While* $t < T$ *do*:
  3. #h(1em) *Assignment step (Hamming projection):* For each $i in {1, dots, n}$, compute $d_(i j) = d_H (x_i, q_j^((t)))$ for $j = 1, dots, k$. Retain $c_i^((t)) <- c_i^((t-1))$ if it already achieves $min_ell d_(i ell)$; otherwise set $c_i^((t)) <- arg min_j d_(i j)$ (smallest index).
  4. #h(1em) Form $C_j^((t)) = {i : c_i^((t)) = j}$ for $j = 1, dots, k$.
  5. #h(1em) *Mode update step:* For each non-empty $C_j^((t))$ and each attribute $r in {1, dots, d}$, compute $f_r (a | C_j^((t)))$ for all $a in A_r$. Retain $q_(j, r)^((t+1)) <- q_(j, r)^((t))$ if it achieves $max_a f_r (a | C_j^((t)))$; else pick the lexicographically first maximizer.
  6. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$ and $Q^((t+1)) = Q^((t))$, *break*.
  7. #h(1em) Increment $t <- t + 1$.
  8. *Return* $(cal(C)^((t)), Q^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of K-Modes: (1) separable mode optimality and a strictly polynomial $O(n d)$ worst-case iteration bound—sharply beating Lloyd's exponential lower bound; (2) an exact one-hot isometric embedding proving a tight factor-2 sandwich inequality with Euclidean K-Means and an induced $O(log k)$ approximation; (3) exact Hard-EM equivalence to symmetric categorical mixtures and a first-principles derivation of cardinality bias; (4) nominal relabeling invariance, coordinate breakdown bounds, and a tie-independent unbounded-suboptimality counterexample; and (5) global NP-hardness.

#subsec("A", "Separable Mode Optimality and Polynomial Iteration Bound")

#lemma(title: "Exact Mode Optimality and Frequency-Gap Identity", origin: "Original Proof")[
  For any non-empty cluster $S subset.eq {1, dots, n}$ and any candidate prototype $q in cal(A)$, the within-cluster Hamming cost satisfies
  $ sum_(i in S) d_H (x_i, q) = d |S| - sum_(r=1)^d f_r (q_r | S). $
  Consequently, $q^*$ minimizes $q |-> sum_(i in S) d_H (x_i, q)$ if and only if $q^*$ is a coordinate-wise mode of $S$, and replacing any $q$ by $q^*$ reduces cost by the exact nonnegative integer gap
  $ sum_(i in S) d_H (x_i, q) - sum_(i in S) d_H (x_i, q^*) \
    = sum_(r=1)^d (f_r (q_r^* | S) - f_r (q_r | S)). $
]

#proof[
  Since $delta(x_(i, r), q_r) = 1 - bold(1)_{x_(i, r) = q_r}$, exchanging the order of summation gives
  $ sum_(i in S) d_H (x_i, q) &= sum_(r=1)^d sum_(i in S) (1 - bold(1)_{x_(i, r) = q_r}) \
    &= d |S| - sum_(r=1)^d f_r (q_r | S). $
  Because the coordinates $q_1, dots, q_d$ decouple completely across $A_1, dots, A_d$, minimizing the left-hand side is equivalent to maximizing each $f_r (q_r | S)$ independently over $A_r$ (recovering and quantifying Huang's mode theorem @huang1998extensions). Subtracting the identities for $q$ and $q^*$ gives the exact integer frequency-gap drop.
]

#theorem(title: "Strict Integer Descent and Polynomial Iteration Bound", origin: "Original Proof")[
  Under Algorithm 1's conservative tie-breaking (retaining $c_i^((t-1))$ and $q_(j, r)^((t))$ upon ties), every non-terminating iteration $t >= 1$ strictly decreases $J_H$ by at least $1$:
  $ J_H (cal(C)^((t)), Q^((t+1))) <= J_H (cal(C)^((t-1)), Q^((t))) - 1. $
  Consequently, Algorithm 1 terminates at a Hamming–Voronoi mode equilibrium in at most
  $ T <= 1 + floor(n sum_(r=1)^d (1 - 1 / m_r)) <= n d $
  iterations—strictly polynomial in $n$ and $d$, in contrast to the $2^(Omega(sqrt(n)))$ worst-case iterations of Euclidean K-Means @vattani2011kmeans.
]

#proof[
  Consider any iteration $t >= 1$ that does not break at Step 6. Then either (i) at least one point $i$ changes cluster ($c_i^((t)) != c_i^((t-1))$), or (ii) assignments remain unchanged and at least one mode coordinate changes ($q_(j, r)^((t+1)) != q_(j, r)^((t))$).

  In case (i), conservative tie-breaking in Step 3 switches $c_i$ only when $d_H (x_i, q_(c_i^((t)))^((t))) < d_H (x_i, q_(c_i^((t-1)))^((t)))$. Because $d_H$ takes values in $bb(Z)_(>= 0)$, every reassigned point lowers $J_H$ by at least $1$, and the subsequent mode update cannot increase $J_H$ by Lemma 1. In case (ii), conservative tie-breaking in Step 5 alters $q_(j, r)$ only when $f_r (q_(j, r)^((t+1)) | C_j^((t))) > f_r (q_(j, r)^((t)) | C_j^((t)))$, which by integrality of frequencies and Lemma 1 lowers $J_H$ by at least $1$. Thus every non-terminating iteration reduces the integer objective $J_H >= 0$ by at least $1$.

  Finally, after the initial mode update $Q^((1))$, the pigeonhole principle ensures $max_(a in A_r) f_r (a | C_j^((0))) >= |C_j^((0))| / m_r$ for every $(j, r)$. Summing over $j$ and $r$ via Lemma 1 bounds the post-first-pass cost by
  $ J_H (cal(C)^((0)), Q^((1))) &<= n d - sum_(r=1)^d n / m_r \
    &= n sum_(r=1)^d (1 - 1 / m_r) < n d. $
  Since $J_H$ drops by at least $1$ per subsequent iteration and remains nonnegative, termination occurs in at most $1 + floor(n sum_(r=1)^d (1 - 1 / m_r)) <= n d$ iterations.
]

#subsec("B", "One-Hot Isometric Embedding and Factor-2 K-Means Bridge")

#theorem(title: "Exact One-Hot Isometry & Factor-2 K-Means Sandwich", origin: "Original Derivation")[
  Define the one-hot embedding $Phi : cal(A) -> {0, 1}^(M_"cat")$ mapping $x = (x_1, dots, x_d)$ to the concatenation of standard basis vectors $Phi(x) = (e^((1))_(x_1), dots, e^((d))_(x_d))$. Then:
  1. For all $x, y in cal(A)$, $d_H (x, y) = 1 / 2 norm(Phi(x) - Phi(y))_1 = 1 / 2 norm(Phi(x) - Phi(y))_2^2$.
  2. For every partition $cal(C)$, the optimal K-Modes cost $J_H^* (cal(C)) = min_Q J_H (cal(C), Q)$ and the optimal one-hot Euclidean K-Means inertia $J_"KM"^* (cal(C)) = sum_(j=1)^k sum_(i in C_j) norm(Phi(x_i) - overline(Phi)_(C_j))_2^2$ satisfy the tight universal factor-2 sandwich inequality:
     $ 1 / 2 J_"KM"^* (cal(C)) <= J_H^* (cal(C)) <= J_"KM"^* (cal(C)). $
]

#proof[
  (1) On each block $r$, $e^((r))_(x_r) - e^((r))_(y_r)$ equals $bold(0)$ if $x_r = y_r$, and has exactly one $+1$ and one $-1$ if $x_r != y_r$. Thus both its $ell_1$ norm and squared $ell_2$ norm equal $2 delta(x_r, y_r)$. Summing over $r = 1, dots, d$ proves Claim 1.

  (2) Fix any cluster $C_j$ and let $hat(p)_(j, r) (a) = f_r (a | C_j) / |C_j|$ be its empirical symbol distribution on attribute $r$, with maximum probability $p^*_(j, r) = max_(a in A_r) hat(p)_(j, r) (a)$. The arithmetic mean $overline(Phi)_(C_j)$ has $r$-th block equal to the probability vector $hat(p)_(j, r)$. Expanding its squared $ell_2$ inertia on block $r$ yields the Gini impurity:
  $ sum_(i in C_j) norm(e^((r))_(x_(i, r)) - hat(p)_(j, r))_2^2 \
    = |C_j| (1 - sum_(a in A_r) hat(p)_(j, r) (a)^2), $
  whereas Lemma 1 gives the K-Modes block cost $|C_j| (1 - p^*_(j, r))$. For any discrete probability distribution $(p_a)$ with maximum entry $p^*$, two elementary inequalities hold:
  - *Upper bound:* $sum_a p_a^2 <= p^* sum_a p_a = p^*$, hence $1 - p^* <= 1 - sum_a p_a^2$.
  - *Lower bound:* $2(1 - p^*) - (1 - sum_a p_a^2) = 1 - 2 p^* + sum_a p_a^2 >= 1 - 2 p^* + (p^*)^2 = (1 - p^*)^2 >= 0$, hence $1 / 2 (1 - sum_a p_a^2) <= 1 - p^*$.

  Multiplying by $|C_j|$ and summing over all $j = 1, dots, k$ and $r = 1, dots, d$ yields $1 / 2 J_"KM"^* (cal(C)) <= J_H^* (cal(C)) <= J_"KM"^* (cal(C))$. Tightness of the upper bound is attained when $hat(p)_(j, r)$ is uniform on $A_r$ with $m_r = 2$, and tightness of the lower bound is approached as $p^*_(j, r) -> 1^-$.
]

#corollary(title: "Induced O(log k) Approximation via Hamming D^2-Seeding", origin: "Original Derivation")[
  Because $norm(Phi(x) - Phi(y))_2^2 = 2 d_H (x, y)$, applying K-Means++ seeding @arthur2007kmeanspp with sampling probability proportional to $min_(q in Q) d_H (x, q)$ (equivalent to squared Euclidean $D^2$-seeding on $Phi(cal(X))$) followed by K-Modes refinement produces a partition satisfying $bb(E)[J_H (cal(C), Q)] <= 16 (ln k + 2) J_H^* (cal(C)^*)$.
]

#proof[
  By @arthur2007kmeanspp, $D^2$-seeding on $Phi(cal(X))$ achieves $bb(E)[J_"KM"^* (cal(C)^((0)))] <= 8(ln k + 2) J_"KM"^* (cal(C)^*_"KM")$. Applying Theorem 2(2) twice and Theorem 1's monotonicity of K-Modes refinement gives
  $ bb(E)[J_H (cal(C), Q)] &<= bb(E)[J_H^* (cal(C)^((0)))] <= bb(E)[J_"KM"^* (cal(C)^((0)))] \
    &<= 8(ln k + 2) J_"KM"^* (cal(C)^*_H) <= 16(ln k + 2) J_H^* (cal(C)^*_H). $
]

#subsec("C", "Probabilistic Hard-EM Equivalence and Cardinality Bias")

#theorem(title: "Symmetric Mixture Equivalence and Cardinality Bias", origin: "Original Derivation")[
  Consider a uniform-prior factorized categorical mixture model (restricted Latent Class Analysis, Alg. 90 @lazarsfeld1968latent) where, given cluster $Z = j$, each attribute $X_r$ equals mode $q_(j, r)$ with probability $1 - eta$ and takes each remaining symbol $a in A_r backslash {q_(j, r)}$ uniformly with probability $eta / (m_r - 1)$, for noise level $eta in (0, 1/2)$. Then the complete-data negative log-likelihood equals
  $ -ln p(cal(X), cal(C) | Q) \
    = C_0 + sum_(j=1)^k sum_(i in C_j) sum_(r=1)^d w_r delta(x_(i, r), q_(j, r)), $
  where $C_0 = n ln k - n d ln(1 - eta)$ and $w_r = ln(((m_r - 1)(1 - eta)) / eta) > 0$. Consequently, unweighted K-Modes is *exact Hard-EM* if and only if all attributes share equal cardinality $m_r equiv m$; otherwise, unweighted Hamming distance systematically under-weights high-cardinality attributes relative to exact maximum likelihood by a factor proportional to $ln(m_r - 1)$.
]

#proof[
  For each observation $i in C_j$ and attribute $r$, the conditional probability mass is
  $ p(x_(i, r) | q_(j, r)) = (1 - eta) (eta / ((m_r - 1)(1 - eta)))^(delta(x_(i, r), q_(j, r))). $
  Because $m_r >= 2$ and $eta < 1 / 2$, we have $(m_r - 1)(1 - eta) / eta > 1$, so $w_r > 0$. Taking $-ln$ and summing over $i = 1, dots, n$, $r = 1, dots, d$ (with $-ln(1/k) = ln k$ per object) yields the stated identity. When $m_r equiv m$, $w_r equiv w > 0$ factors out identically, leaving $C_0 + w J_H (cal(C), Q)$; when $m_r$ varies across attributes, exact log-likelihood requires weighting mismatch on coordinate $r$ by $w_r$, growing logarithmically in $m_r - 1$.
]

#subsec("D", "Invariances, Breakdown Robustness, and Adversarial Traps")

#proposition(title: "Invariance, Breakdown Margin, and Unbounded Suboptimality", origin: "Original Counterexample")[
  1. *Relabeling Invariance:* $J_H$ and the set of K-Modes equilibria are strictly invariant under arbitrary per-attribute symbol bijections $pi_r : A_r -> A_r$ and feature permutations $sigma in S_d$.
  2. *Coordinate Breakdown Point:* In a fixed cluster $C_j$ of size $m$ whose mode on attribute $r$ has frequency $f_1$ and runner-up frequency $f_2$, corrupting $b < ceil((f_1 - f_2 + 1) / 2)$ records cannot alter $q_(j, r)$ (reaching $ceil(m / 2) / m -> 1 / 2$ for a pure cluster).
  3. *Unbounded Suboptimality:* For any $M > 0$, there exists a 4-point binary dataset $cal(X) subset.eq {0, 1}^d$ ($k = 3$) and a distinct-object initialization such that K-Modes halts immediately at a strict fixed point with $J_H (cal(C), Q) / J_H^* > M$ under *every* mode tie-breaking rule.
]

#proof[
  (1) Every bijection $pi_r$ satisfies $delta(pi_r (a), pi_r (b)) = delta(a, b)$, and summation over $r$ is permutation-invariant.

  (2) Changing at most $b$ entries in $C_j$ decreases $f_1$ by at most $b$ and increases $f_2$ by at most $b$. The mode strictly persists whenever $f_1 - b > f_2 + b$, i.e., $2 b < f_1 - f_2$, equivalently $b < ceil((f_1 - f_2 + 1) / 2)$.

  (3) Set $L = ceil(M / 2) + 1$ and $d = 4 L + 1$, partitioning coordinates into four blocks of sizes $(1, 2 L, L, L)$. Define four binary strings:
  $ x_1 &= (0, 0^(2L), 0^L, 0^L), & x_2 &= (1, 0^(2L), 0^L, 0^L), \
    x_3 &= (0, 1^(2L), 1^L, 0^L), & x_4 &= (0, 1^(2L), 0^L, 1^L). $
  Verify exact pairwise Hamming distances: $d_H (x_1, x_2) = 1$, $d_H (x_3, x_4) = 2 L$, $d_H (x_1, x_3) = d_H (x_1, x_4) = 3 L$, and $d_H (x_2, x_3) = d_H (x_2, x_4) = 3 L + 1$. Initialize $k = 3$ at $Q^((0)) = (x_1, x_2, x_3)$. Then $x_1 -> q_1$, $x_2 -> q_2$, $x_3 -> q_3$, and $x_4 -> q_3$ (since $2 L < 3 L$), giving $cal(C) = {{x_1}, {x_2}, {x_3, x_4}}$. Any valid mode $q_3$ of ${x_3, x_4}$ has $0$ in block 1, $1^(2L)$ in block 2, and arbitrary bits in blocks 3–4, so $d_H (x_3, q_3) + d_H (x_4, q_3) = 2 L$ and $max(d_H (x_3, q_3), d_H (x_4, q_3)) <= 2 L < 3 L$. Thus $cal(C)$ is a strict fixed point of cost $2 L$ regardless of tie-breaking. By contrast, $cal(C)^* = {{x_1, x_2}, {x_3}, {x_4}}$ attains cost $J_H^* = 1$, yielding $J_H / J_H^* = 2 L > M$.
]

#subsec("E", "NP-Hardness of Global K-Modes Optimization")

#theorem(title: "NP-Hardness of Optimal Hamming K-Modes", origin: "Classical — Cited")[
  Minimizing $J_H (cal(C), Q)$ over $k$-partitions of categorical vectors is NP-hard even for $k = 2$ on binary alphabets $cal(A) = {0, 1}^d$ @huang1998extensions @aloise2009np.
]

#proof(label: "Proof Synthesis (Cited)")[
  By Theorem 2(1), for binary attributes $m_r = 2$, minimizing $J_H$ over $k = 2$ clusters with discrete prototypes in ${0, 1}^d$ is equivalent to minimizing $ell_1$-median / hypercube consensus clustering, which reduces directly from the NP-hard *Densest Cut* / *2-Means* problem on binary incidence matrices @aloise2009np (and from Exact Cover / Consensus String for general $k$ @huang1998extensions).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Algorithm 1 evaluates $n k$ Hamming distances in $O(n k d)$ integer comparisons and updates $k$ per-attribute frequency tables in $O(n d + k M_"cat")$ time, giving total runtime $O((n k d + k M_"cat") T)$—strictly linear in $n$ when $M_"cat" = O(n d)$—with guaranteed polynomial iteration ceiling $T <= n d$ (Theorem 1) and typically $T in [5, 20]$ in practice. Memory consumption is $O(n d + k M_"cat")$ using compact integer indices (never materializing sparse $M_"cat"$-dimensional one-hot matrices). Geometrically, K-Modes partitions vertices of the Hamming graph $K_(m_1) square dots square K_(m_d)$; because Hamming distances take only $d + 1$ discrete values, Voronoi ties occur with non-negligible probability—making Huang's frequency-dispersed seeding @huang1998extensions or Hamming $D^2$-seeding (Corollary 1) essential to prevent empty clusters and poor traps (Proposition 1(3)).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for K-Modes],
    [Time / Space Complexity], [$O((n k d + k M_"cat") T)$ time with proven $T <= n d$ (Thm. 1); $O(n d + k M_"cat")$ space.],
    [Geometric Inductive Bias], [Hamming balls on nominal product spaces; coordinate-wise independent marginal modes.],
    [Approximation Guarantee], [Arbitrary traps under random seeds (Prop. 1); $O(log k)$ expected ratio via $D^2$-seeding (Cor. 1).],
    [Invariances], [Strictly invariant to arbitrary nominal relabeling $pi_r$ and feature permutations $sigma$.],
    [Outlier Robustness], [Coordinate breakdown up to $1/2$ for pure modes; bounded mismatch $d_H <= d$.],
    [Recommended Companions], [K-Prototypes (Alg. 8) for mixed data; ROCK (Alg. 67) for sparse links; LCA (Alg. 90) for soft probabilistic classes; W-K-Means (Alg. 28) @huang2005automated for learned weights.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of K-Modes.],
)

#sec("IV", "Engineering and Scientific Applications")

K-Modes is the canonical linear-time workhorse for large purely nominal relational tables where arithmetic averaging is invalid and interpretable symbolic prototypes are required.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Electronic Health Records (EHR) & Diagnostic Phenotyping:* Patient encounters encode nominal ICD-10 diagnosis codes, drug classes, admission types, and demographic categories. K-Modes runs in linear time over millions of visits while yielding clinically readable modal patient profiles.
2. *Cybersecurity & Network Intrusion Log Triage:* Firewall and SIEM events consist of categorical protocols, HTTP verbs, TLS cipher suites, and status codes. Integer-only Hamming comparisons and Theorem 1's rapid convergence support high-throughput log clustering without one-hot memory explosion.
3. *E-Commerce & Categorical Customer Profiling:* Session attributes (device OS, channel, payment method, product category) are purely nominal; modal prototypes directly translate into human-auditable marketing personas without distorting labels via arbitrary numeric encoding.
4. *Genomics & Multi-Locus SNP Haplotyping:* Diploid genotypes ${A A, A a, a a}$ across $d$ polymorphic loci share uniform alphabet size $m_r = 3$, matching the exact equal-cardinality regime of Theorem 3 where K-Modes coincides with Hard-EM under a symmetric mutation model.
5. *Psychometric & Survey Typology Discovery:* Multiple-choice questionnaires with unordered categorical responses resist Euclidean K-Means; K-Modes preserves exact relabeling invariance (Proposition 1(1)) while summarizing respondent archetypes via modal answers.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

K-Modes fails on (a) *continuous or ordinal numerical features*, where discretization discards metric ordering and treats adjacent bins identically to distant extremes; (b) *heterogeneous-cardinality tables* where unweighted Hamming distance under-weights rare matches in large alphabets (Theorem 3, mitigated by frequency-weighted matching or automated variable weighting @huang2005automated); (c) *ultra-sparse binary transaction baskets*, where the coordinate-wise mode collapses to the all-zero vector (resolved by ROCK Alg. 67 @guha2000rock or Bernoulli Mixtures Alg. 96); and (d) *strongly correlated categorical attributes* violating conditional independence (resolved by LCA Alg. 90 or AutoClass Alg. 87 @cheeseman1996autoclass).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Modes Fits / Prerequisite], [Caution / Alternative],
    [Clinical EHR], [Nominal ICD/drug codes; modal profiles], [Use K-Prototypes (Alg. 8) if labs],
    [SIEM / NetLogs], [Integer Hamming speed; linear RAM], [Weight high-cardinality ports],
    [SNP Haplotypes], [Uniform $m_r = 3$ matches Thm. 3 MLE], [Use LCA (Alg. 90) for soft calls],
    [Market Baskets], [Sparse 0/1 collapses to all-zero mode], [Use ROCK (Alg. 67) or Jaccard],
    [Continuous Data], [Binning destroys metric geometry], [Use K-Means (Alg. 1) or PAM (Alg. 4)],
  ),
  caption: [Application suitability and engineering matrix for K-Modes.],
)

#sec("V", "Structural Relations and Empirical Behavior")

K-Modes serves as the discrete categorical counterpart to Lloyd's K-Means (Alg. 1) and anchors a continuous-to-discrete spectrum spanning partitioning, hierarchical, and probabilistic paradigms.

#proposition(title: "Embedding Inside K-Prototypes, LCA, and COBWEB", origin: "Original Derivation")[
  1. *K-Prototypes Unification:* In Huang's mixed-type K-Prototypes objective (Alg. 8 @huang1998extensions), $J_gamma = sum_(j=1)^k sum_(i in C_j) (norm(x_i^"num" - mu_j)_2^2 + gamma d_H (x_i^"cat", q_j))$, setting $d_"num" = 0$ recovers K-Modes identically while setting $d_"cat" = 0$ recovers Lloyd's K-Means (Alg. 1).
  2. *Probabilistic & Conceptual Generalizations:* Relaxing hard assignments and equal-dispersion constraints in Theorem 3 to unconstrained conditional tables $theta_(j, r, a) = bb(P)(X_r = a | Z = j)$ yields Latent Class Analysis (Alg. 90 @lazarsfeld1968latent), while replacing Hamming mismatch with expected category-utility gain in a concept hierarchy yields COBWEB (Alg. 85 @fisher1987cobweb).
]

#proof[
  Claim 1 holds verbatim by substituting $d_"num" = 0$ (with $gamma = 1$) or $d_"cat" = 0$ into $J_gamma$. For Claim 2, allowing unrestricted $theta_(j, r, a) in Delta^(m_r - 1)$ in the log-likelihood of Theorem 3 replaces $-w_r delta(x_(i, r), q_(j, r))$ by $ln theta_(j, r, x_(i, r))$, whose soft EM maximization over mixing weights $pi_j$ and tables $theta_(j, r, dot.c)$ is precisely standard LCA @lazarsfeld1968latent.
]

Empirical results on the ClusterOverview 4-dataset benchmark suite vividly corroborate the mathematical boundary between continuous and nominal spaces: K-Modes records `Sil = 0.332` and `ARI = 0.486` on Iris, `Sil = 0.277` on Customers, `Sil = 0.441` (`ARI = N/A`) on Multi-Shape, and `Sil = 0.072` (`ARI = N/A`) on High-D Subspace+Noise, while achieving an exceptionally fast mean runtime of `0.136 ms` ($6.2 times$ faster than CLARA and $2.8 times$ faster than CLARANS). The ultra-low latency reflects integer Hamming arithmetic and Theorem 1's strict polynomial convergence, whereas the reduced Silhouette and ARI stem directly from discretizing intrinsically continuous benchmark features into unordered bins—erasing inter-bin metric geometry and amplifying uniform Hamming noise in High-D Subspace+Noise (`Sil = 0.072`). These metrics underscore a cardinal engineering rule: deploy K-Modes natively on genuinely nominal categorical data (or K-Prototypes Alg. 8 on mixed tables), never via lossy discretization of continuous Euclidean manifolds.
