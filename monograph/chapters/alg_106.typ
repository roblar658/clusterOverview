#import "../template/ieee.typ": *

#alg_entry(
  id: 106,
  name: "Affinity Propagation",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2007",
  paradigm: "Factor-Graph Max-Sum Message Passing",
  complexity: $O(T n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Unlike centroid methods ($K$-means @lloyd1982least) or classical $K$-medoids (PAM @kaufman1990finding), which require fixing the cluster count $k$ a priori and suffer sensitivity to random initial seed placement, *Affinity Propagation* (AP), introduced by Frey and Dueck (@frey2007affinity), treats *every* data point simultaneously as a candidate exemplar and determines both the optimal exemplar subset and cluster assignments deterministically via max-sum message passing on a factor graph. Crucially, AP operates on an arbitrary real-valued pairwise similarity matrix $S = [s(i, k)] in bb(R)^(n times n)$ that need not be symmetric, nonnegative, or metric.

The diagonal entries $p_k = s(k, k) in bb(R)$ are called *preferences*: larger $p_k$ increases the propensity of point $k$ to become an exemplar (commonly set uniformly to the median or minimum of off-diagonal similarities). Clustering is encoded by an assignment vector $c = (c_1, dots, c_n) in {1, dots, n}^n$, where $c_i = k$ indicates that point $i$ selects point $k$ as its exemplar (and $c_k = k$ indicates that $k$ is an exemplar).

#definition(title: "Exemplar-Consistent Net Similarity Objective")[
  Given similarity matrix $S in bb(R)^(n times n)$, Affinity Propagation seeks an exemplar assignment $c in {1, dots, n}^n$ maximizing the constrained net similarity
  $ cal(E)(c) = sum_(i=1)^n s(i, c_i) + sum_(k=1)^n delta_k (c), $
  where the exemplar-consistency barrier $delta_k (c)$ enforces that any point selected as an exemplar by another point must also select itself:
  $ delta_k (c) = cases(-infinity & "if" c_k != k "and" exists i "with" c_i = k, 0 & "otherwise.") $
  Equivalently, with binary indicators $C_(i k) = bb(I)(c_i = k) in {0, 1}$, AP maximizes $sum_(i,k=1)^n C_(i k) s(i, k)$ subject to the row simplex constraint $sum_(k=1)^n C_(i k) = 1$ ($forall i$) and consistency constraint $C_(i k) <= C_(k k)$ ($forall i, k$).
]

Instead of solving this NP-hard Uncapacitated Facility Location integer program directly, AP runs damped max-sum belief propagation between two scalar message matrices: *responsibility* $r(i, k)$, quantifying how well-suited $k$ is to serve as exemplar for $i$ relative to competing candidates $k' != k$, and *availability* $a(i, k)$, quantifying accumulated evidence from other points $i' != i$ supporting $k$ as an exemplar (@frey2007affinity).

#algbox(
  title: "Affinity Propagation (Damped Max-Sum Message Passing)",
  inputs: [Similarity matrix $S in bb(R)^(n times n)$ with preferences $p_k = s(k, k)$, damping $lambda in [0.5, 1)$, max iterations $T$.],
  outputs: [Exemplar set $E = {k : hat(c)_k = k}$ and assignments $hat(c)_i in E$.],
)[
  1. *Initialization:* Set $a^((0))(i, k) = 0$ and $r^((0))(i, k) = 0$ for all $i, k in {1, dots, n}$.
  2. *Responsibility Update:* For each $i, k$, compute the raw responsibility and damp:
     $ tilde(r)(i, k) &= s(i, k) - max_(k' != k) {a^((t-1))(i, k') + s(i, k')}, \
       r^((t))(i, k) &= lambda r^((t-1))(i, k) + (1 - lambda) tilde(r)(i, k). $
  3. *Availability Update:* Compute raw self- and off-diagonal availabilities and damp:
     $ tilde(a)(k, k) &= sum_(i' != k) max{0, r^((t))(i', k)}, \
       tilde(a)(i, k) &= min{0, r^((t))(k, k) + sum_(i' in.not {i, k}) max{0, r^((t))(i', k)}} quad (i != k), \
       a^((t))(i, k) &= lambda a^((t-1))(i, k) + (1 - lambda) tilde(a)(i, k). $
  4. *Convergence Check & Decoding:* Stop if $t = T$ or exemplar assignments $hat(c)_i = arg max_k {a^((t))(i, k) + r^((t))(i, k)}$ remain unchanged for $T_"conv"$ consecutive steps (e.g., $15$); output $hat(c)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Scalar Reduction of Binary Factor-Graph Max-Sum BP")

A central mathematical achievement of Affinity Propagation is collapsing $2 n^2$ vector-valued binary belief-propagation messages over $n$-ary constraints into two simple scalar recursions without approximation.

#theorem(title: "Exact Scalar Reduction of Binary Max-Sum BP", origin: "Original Derivation")[
  Represent Definition 1 as a factor graph over binary variables $C_(i k) in {0, 1}$ with unary potentials $psi_(i k)(C_(i k)) = C_(i k) s(i, k)$, row factors $I_i (C_(i,:)) = 0$ if $sum_k C_(i k) = 1$ (else $-infinity$), and column factors $E_k (C_(:,k)) = 0$ if $C_(k k) = 1$ or $sum_i C_(i k) = 0$ (else $-infinity$). Let $rho_(i k)(b)$ and $alpha_(i k)(b)$ ($b in {0, 1}$) denote the binary max-sum messages from $I_i$ to $C_(i k)$ and from $E_k$ to $C_(i k)$, respectively. Defining scalar log-odds differences $r(i, k) = rho_(i k)(1) - rho_(i k)(0)$ and $a(i, k) = alpha_(i k)(1) - alpha_(i k)(0)$ yields *identically* the AP recursions in Steps 2--3 of the algorithm box.
]

#proof[
  In standard max-sum belief propagation (@bishop2006pattern), variable $C_(i k)$ sends incoming message $beta_(i k)(b) = alpha_(i k)(b) + b dot s(i, k)$ to row factor $I_i$ and $eta_(i k)(b) = rho_(i k)(b) + b dot s(i, k)$ to column factor $E_k$. Wait: notice the unary term $b dot s(i, k)$ enters once in the factor graph; following the exact symmetric convention of @frey2007affinity, attach $s(i, k)$ to the edge entering $I_i$, so the message from $C_(i k)$ to $I_i$ is $beta_(i k)(b) = alpha_(i k)(b) + b dot s(i, k)$ while the message from $C_(i k)$ to $E_k$ is $rho_(i k)(b)$ (shifting $r(i, k)$ by $s(i, k)$).

  *Row-factor message $I_i -> C_(i k)$:* Enforcing $sum_(k') C_(i k') = 1$:
  - If $C_(i k) = 1$, every $k' != k$ must have $C_(i k') = 0$, so
    $ rho_(i k)(1) = s(i, k) + sum_(k' != k) alpha_(i k')(0). $
  - If $C_(i k) = 0$, exactly one $k' != k$ must equal $1$ and all others $0$, so
    $ rho_(i k)(0) = max_(k' != k) (s(i, k') + alpha_(i k')(1) + sum_(j in.not {k, k'}) alpha_(i j)(0)). $
  Subtracting $rho_(i k)(0)$ from $rho_(i k)(1)$ cancels $sum_(j != k) alpha_(i j)(0)$ identically:
  $ r(i, k) = rho_(i k)(1) - rho_(i k)(0) = s(i, k) - max_(k' != k) {a(i, k') + s(i, k')}. $

  *Column-factor message $E_k -> C_(i k)$:*
  - Case $i = k$ (self-availability): If $C_(k k) = 1$, every $i' != k$ may freely take $0$ or $1$, giving $alpha_(k k)(1) = sum_(i' != k) max{rho_(i' k)(0), rho_(i' k)(1)}$. If $C_(k k) = 0$, exemplar consistency forces $C_(i' k) = 0$ for all $i' != k$, giving $alpha_(k k)(0) = sum_(i' != k) rho_(i' k)(0)$. Subtracting yields
    $ a(k, k) = sum_(i' != k) max{0, rho_(i' k)(1) - rho_(i' k)(0)} = sum_(i' != k) max{0, r(i', k)}. $
  - Case $i != k$: If $C_(i k) = 1$, consistency forces $C_(k k) = 1$ while all $i' in.not {i, k}$ are free, so
    $ alpha_(i k)(1) = rho_(k k)(1) + sum_(i' in.not {i, k}) max{rho_(i' k)(0), rho_(i' k)(1)}. $
    If $C_(i k) = 0$, we choose the better of $C_(k k) = 1$ (all $i' in.not {i, k}$ free) or $C_(k k) = 0$ (all $i' in.not {i, k}$ forced to $0$), giving
    $ alpha_(i k)(0) = max(rho_(k k)(1) + sum_(i' in.not {i, k}) max{rho_(i' k)(0), rho_(i' k)(1)}, \ rho_(k k)(0) + sum_(i' in.not {i, k}) rho_(i' k)(0)). $
    Subtracting $alpha_(i k)(0)$ from $alpha_(i k)(1)$ cancels the first branch to $0$ and subtracts the second branch to $r(k, k) + sum_(i' in.not {i, k}) max{0, r(i', k)}$, yielding
    $ a(i, k) = min{0, r(k, k) + sum_(i' in.not {i, k}) max{0, r(i', k)}}. $
  Finally, the max-marginal difference for $C_(i k) = 1$ vs. $0$ is $a(i, k) + r(i, k)$, confirming the decoding rule $hat(c)_i = arg max_k {a(i, k) + r(i, k)}$.
]

#subsec("B", "Preference Monotonicity and Exact Phase Transitions")

#theorem(title: "Global Preference Monotonicity and Exemplar Thresholds", origin: "Original Proof")[
  Let $S$ have fixed off-diagonal entries $s(i, j)$ ($i != j$) with $s_"min" = min_(i != j) s(i, j)$ and $s_"max" = max_(i != j) s(i, j)$, and uniform preference $s(k, k) = p in bb(R)$ for all $k$. Let $K^*(p)$ denote the number of exemplars in a global maximizer $c^*(p)$ of $cal(E)_p (c)$ for $n >= 2$. Then:
  (i) *Exact Monotonicity:* for any $p_1 < p_2$, every pair of global maximizers satisfies $K^*(p_1) <= K^*(p_2)$;
  (ii) *Singleton Collapse:* whenever $p > s_"max"$, the unique global maximizer has $K^*(p) = n$ ($c_i^* = i$ for all $i$);
  (iii) *Single-Cluster Regime:* whenever $p < s_"min" - (n - 2)(s_"max" - s_"min")$, every global maximizer has $K^*(p) = 1$, centered at the $1$-medoid $k^* = arg max_k sum_(i != k) s(i, k)$.
]

#proof[
  Write $cal(E)_p (c) = A(c) + K(c) p$, where $K(c) = |{k : c_k = k}|$ is the exemplar count and $A(c) = sum_(i : c_i != i) s(i, c_i)$ is the off-diagonal assignment sum.

  For (i), let $c^((1))$ maximize $cal(E)_(p_1)$ with $K_1 = K(c^((1)))$ and $c^((2))$ maximize $cal(E)_(p_2)$ with $K_2 = K(c^((2)))$. By global optimality of each assignment at its respective preference,
  $ A(c^((1))) + K_1 p_1 &>= A(c^((2))) + K_2 p_1, \
    A(c^((2))) + K_2 p_2 &>= A(c^((1))) + K_1 p_2. $
  Adding the two inequalities cancels $A(c^((1))) + A(c^((2)))$ and yields $(K_2 - K_1)(p_2 - p_1) >= 0$. Since $p_2 - p_1 > 0$, we obtain $K_2 >= K_1$.

  For (ii), if $p > s_"max"$, then for any assignment $c$ with $K < n$, we have $cal(E)_p (c) = K p + sum_(i : c_i != i) s(i, c_i) <= K p + (n - K) s_"max" < n p = cal(E)_p ("id")$, proving $K^*(p) = n$ uniquely.

  For (iii), let $c^((1))$ be the optimal $1$-exemplar solution with medoid $k^*$, achieving $cal(E)_p (c^((1))) = p + sum_(i != k^*) s(i, k^*) >= p + (n - 1) s_"min"$. Any valid assignment $c$ with $K >= 2$ exemplars has $n - K$ non-exemplars, so $cal(E)_p (c) <= K p + (n - K) s_"max"$. Thus
  $ cal(E)_p (c^((1))) - cal(E)_p (c)
    &>= (n - 1) s_"min" - (n - K) s_"max" - (K - 1) p \
    &= (K - 1)(s_"min" - p) - (n - K)(s_"max" - s_"min"). $
  For every $K in {2, dots, n}$ and $n >= 2$, note that $(K - 1)(n - 2) - (n - K) = (K - 2)(n - 1) >= 0$, so $(n - K) <= (K - 1)(n - 2)$. Substituting $p < s_"min" - (n - 2)(s_"max" - s_"min")$ yields $cal(E)_p (c^((1))) - cal(E)_p (c) > (K - 1)(n - 2)(s_"max" - s_"min") - (n - K)(s_"max" - s_"min") >= 0$, forcing $K^*(p) = 1$.
]

#subsec("C", "Undamped Limit-Cycle Oscillation and Damping Contraction")

Because AP's factor graph contains dense 4-cycles, undamped belief propagation ($lambda = 0$) need not converge even on trivial symmetric graphs. We prove both an exact period-2 oscillation counterexample and the exact Schur-Cohn contraction law explaining why damping $lambda in (0, 1)$ stabilizes the dynamics.

#theorem(title: "Exact 2-Node Oscillation and Damping Spectral Radius", origin: "Original Counterexample")[
  Consider $n = 2$ vertices with symmetric similarity $s(1, 2) = s(2, 1) = 0$ and preference $s(1, 1) = s(2, 2) = p < 0$, initialized at $a^((0)) = 0$.
  (i) With zero damping ($lambda = 0$), AP enters an exact perpetual period-2 limit cycle: $a^((2 t)) = 0$ and $a^((2 t + 1)) = mat(-p, p; p, -p)$ for all $t >= 0$.
  (ii) For any $lambda in (0, 1)$, the unique symmetric fixed point $r_"diag"^* = a_"off"^* = p / 2$, $r_"off"^* = a_"diag"^* = -p / 2$ is locally exponentially stable with exact linearized Jacobian $M(lambda)$ having characteristic polynomial $q(mu) = mu^2 - (2 lambda - (1 - lambda)^2) mu + lambda^2$; its spectral radius satisfies $rho(M(lambda)) < 1$ for all $lambda in (0, 1)$, and equals $rho(M(lambda)) = lambda$ for all $lambda in (3 - 2 sqrt(2), 1)$.
]

#proof[
  For (i) with $lambda = 0$, starting from $a^((0)) = 0$:
  - Iteration $t = 1$: $r^((1))(1, 1) = p - (0 + 0) = p < 0$, $r^((1))(1, 2) = 0 - (0 + p) = -p > 0$ (and symmetrically $r^((1))(2, 2) = p$, $r^((1))(2, 1) = -p$). Updating availabilities gives $a^((1))(1, 1) = max(0, -p) = -p > 0$ and $a^((1))(1, 2) = min(0, p) = p < 0$ (since no $i' in.not {1, 2}$ exists), and symmetrically $a^((1))(2, 2) = -p$, $a^((1))(2, 1) = p$.
  - Iteration $t = 2$: $r^((2))(1, 1) = p - (p + 0) = 0$ and $r^((2))(1, 2) = 0 - (-p + p) = 0$, so $r^((2)) = 0$. Consequently $a^((2))(1, 1) = max(0, 0) = 0$ and $a^((2))(1, 2) = min(0, 0) = 0$, returning identically to $a^((2)) = a^((0)) = 0$ and repeating with period $2$ forever.

  For (ii), inside the cone $r_"diag" < 0, r_"off" > 0$, the symmetric updates decouple into two identical 2D linear systems. For $(r_"diag", a_"off")$, setting $e_r^((t)) = r_"diag"^((t)) - p / 2$ and $e_a^((t)) = a_"off"^((t)) - p / 2$ yields
  $ e_r^((t)) &= lambda e_r^((t-1)) - (1 - lambda) e_a^((t-1)), \
    e_a^((t)) &= lambda e_a^((t-1)) + (1 - lambda) e_r^((t)), $
  giving the $2 times 2$ iteration matrix
  $ M(lambda) = mat(lambda, -(1 - lambda); lambda (1 - lambda), lambda - (1 - lambda)^2). $
  Its trace and determinant are $"Tr"(M) = 2 lambda - (1 - lambda)^2$ and $det(M) = lambda^2$, yielding $q(mu) = mu^2 - (2 lambda - (1 - lambda)^2) mu + lambda^2$. By the Jury / Schur-Cohn stability criterion for quadratic polynomials, both roots lie strictly inside the unit disk iff: (a) $det(M) = lambda^2 < 1$, (b) $q(1) = 2(1 - lambda)^2 > 0$, and (c) $q(-1) = 4 lambda > 0$, all of which hold strictly for every $lambda in (0, 1)$ (while $q(-1) = 0$ at $lambda = 0$ confirms the period-2 bifurcation $mu = -1$). Moreover, the discriminant $Delta = -(1 - lambda)^2 (4 lambda - (1 - lambda)^2)$ is strictly negative for $lambda in (3 - 2 sqrt(2), 1)$, so the roots form a complex conjugate pair with exact modulus $|mu_(1,2)| = sqrt(det(M)) = lambda < 1$.
]

#subsec("D", "Affine Invariance and LP-Relaxation Optimality")

#theorem(title: "Affine Invariance and Integrality Optimality", origin: "Classical — Cited")[
  (i) For any positive scale $alpha > 0$, constant shift $beta in bb(R)$, and permutation matrix $P$, transforming $S |-> alpha P S P^top + beta bb(1) bb(1)^top$ leaves all optimal exemplar sets of Definition 1 and all damped AP iterations (up to scaling messages by $alpha$ and permuting indices by $P$) strictly invariant.
  (ii) Whenever the linear programming relaxation of Definition 1 (relaxing $C_(i k) in [0, 1]$) has a unique integral optimum $C^*$, any fixed point of AP message passing with unique max-marginals $hat(c)_i = arg max_k (a(i, k) + r(i, k))$ recovers the exact global optimum $C^*$ (@frey2007affinity, @bishop2006pattern).
]

#proof[
  Under $S' = alpha S + beta bb(1) bb(1)^top$, $cal(E)'(c) = sum_(i=1)^n (alpha s(i, c_i) + beta) = alpha cal(E)(c) + n beta$, preserving all global optima. By induction on $t$, assume $a'^((t-1)) = alpha a^((t-1))$ and $r'^((t-1)) = alpha r^((t-1))$. In Step 2, the shift $beta$ cancels identically inside $s'(i, k) - max_(k' != k) (a'(i, k') + s'(i, k')) = alpha (s(i, k) - max_(k' != k) (a(i, k') + s(i, k')))$, so $r'^((t)) = alpha r^((t))$. Because $max{0, alpha x} = alpha max{0, x}$ and $min{0, alpha x} = alpha min{0, x}$ for $alpha > 0$, Step 3 preserves $a'^((t)) = alpha a^((t))$. Simultaneous relabeling $P S P^top$ commutes with all row maxima and column sums. Statement (ii) follows from the exact correspondence between AP max-sum fixed points and dual-feasible reparameterizations of the Uncapacitated Facility Location LP relaxation (@frey2007affinity).
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each AP iteration updates all $n^2$ responsibilities and $n^2$ availabilities in $O(n^2)$ time by precomputing for each row $i$ the largest and second-largest values of $a(i, k') + s(i, k')$ (avoiding an $O(n^3)$ naive max loop) and for each column $k$ the positive sum $sum_(i' != k) max{0, r(i', k)}$. Total runtime over $T$ iterations is $O(T n^2)$ with $O(n^2)$ memory for $S, R, A$. On sparse similarity graphs with $m = |E| << n^2$ allowed exemplar links, updates restrict to edges in $E$, reducing cost to $O(T m)$ time and $O(m)$ space.

Geometrically, AP partitions space into exemplar-centered Voronoi cells under *any* similarity function (asymmetric divergences, edit distances, negative KL divergence, or non-metric graphs), selecting real observed data points as interpretable exemplars without needing vector coordinates. Its quadratic memory and iterative sweeps make dense AP slower than $K$-means and small-$n$ spectral methods, yet it automatically discovers cluster cardinality $K$ via preference $p$ (Theorem 2).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Dense / Sparse Complexity], [$O(T n^2)$ dense; $O(T m)$ on sparse similarity graph $E$.],
    [Memory Footprint], [$O(n^2)$ for $S, R, A$ (or $O(m)$ sparse).],
    [Input Generality], [Arbitrary asymmetric, negative, or non-metric $S$.],
    [Cluster Representation], [Real data exemplars (medoids); $K$ chosen automatically.],
    [Invariances], [Positive affine $S |-> alpha S + beta bb(1) bb(1)^top$ & permutation (Thm. 4).],
    [Convergence Control], [Damping $lambda in [0.5, 0.99]$ prevents oscillations (Thm. 3).],
  ),
  caption: [Operational, complexity, and theoretical profile of Affinity Propagation.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Bioinformatics Exon Discovery and Protein Sequence Summarization:* Frey and Dueck (@frey2007affinity) demonstrated AP on microarray exon clustering and protein sequence alignment matrices where similarities (BLAST / Smith-Waterman scores) are asymmetric and non-metric, and a real biological sequence (exemplar) is required to represent each family rather than an invalid averaged coordinate vector.
2. *Face Image and Keyframe Exemplar Selection:* In computer vision datasets with thousands of faces or video frames, AP selects crisp canonical prototype frames directly from pairwise perceptual or descriptor similarities, avoiding the blurry ghost prototypes produced by Euclidean averaging in pixel or latent space.
3. *Airline Routing, Facility Location, and Sensor Placement:* Because Definition 1 is mathematically identical to Uncapacitated Facility Location with opening cost $-p_k$, AP directly solves hub selection in logistics and wireless sensor network cluster-head election, naturally honoring asymmetric travel times or directed wireless link qualities $s(i, k) != s(k, i)$.
4. *Document and News-Headline Extractive Summarization:* When pairwise sentence similarities are computed via asymmetric entailment or ROUGE overlap, AP's exemplars directly constitute extractive summary sentences covering distinct subtopics without pre-specifying the summary length $K$.
5. *Non-Convex Multi-Shape Sub-Prototype Compression:* When preference $p$ is moderately set on non-convex shapes, AP covers curved manifolds with tight, well-separated local exemplars that achieve high Silhouette cohesion and serve as clean codebook compressions for downstream hierarchical merging.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

AP suffers four principal failure modes: (i) *quadratic $O(n^2)$ memory bottleneck*, precluding dense execution beyond $n approx 2 times 10^4$ without sparse graph pruning; (ii) *preference-scale sensitivity*, where default median preference often over-splits compact datasets into too many small exemplar clusters (Theorem 2); (iii) *loopy oscillation / non-convergence* when $lambda$ is too small or near degenerate symmetries (Theorem 3), requiring small random noise dithering on $S$ and $lambda in [0.5, 0.9]$; and (iv) *star-shaped exemplar bias*, which splits elongated manifolds into chains of separate star clusters unless paired with geodesic similarities.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Asymmetric / non-metric data], [Excellent: no metric or symmetry needed.], [Add tiny noise to break ties.],
    [Exemplar summarization / hubs], [Exact facility-location model (Def. 1).], [Calibrate preference $p$ carefully.],
    [Unknown moderate $K$], [High: monotone control via $p$ (Thm. 2).], [Median $p$ may over-partition.],
    [Large datasets ($n > 20,000$)], [Poor dense $O(n^2)$ memory.], [Use sparse $q$-NN AP or FasterPAM (@schubert2019fasterpam).],
    [Single elongated filaments], [Medium: breaks into local stars.], [Use Ncut (Alg. 104) or HDBSCAN (@campello2013hdbscan).],
  ),
  caption: [Domain suitability and failure-mode matrix for Affinity Propagation.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Affinity Propagation connects directly to $K$-medoids (PAM @kaufman1990finding, FasterPAM @schubert2019fasterpam) as its Lagrangian relaxation over cluster cardinality $K$, while replacing discrete swap heuristics with parallel message passing.

#proposition(title: "Lagrangian Equivalence to K-Medoids", origin: "Original Proof")[
  Let $d(i, j) >= 0$ be any dissimilarity matrix ($d(i, i) = 0$) and set $s(i, j) = -d(i, j)$ for $i != j$ with uniform preference $s(k, k) = p in bb(R)$. If $c^*$ is a global maximizer of AP's objective $cal(E)_p (c)$ using $K = K^*(p)$ exemplars, then $c^*$ is simultaneously a global minimizer of the $K$-medoids total distortion $J_"med"(c) = sum_(i=1)^n d(i, c_i)$ over all valid $K$-medoid partitions of cardinality $K$.
]

#proof[
  For any valid exemplar assignment $c$ with exactly $K$ exemplars, $c_k = k$ contributes $s(k, k) = p$ while $d(k, k) = 0$, and every non-exemplar $i$ contributes $s(i, c_i) = -d(i, c_i)$. Hence
  $ cal(E)_p (c) = K p - sum_(i=1)^n d(i, c_i) = K p - J_"med"(c). $
  If any valid $K$-medoid assignment $tilde(c)$ with $K(tilde(c)) = K$ satisfied $J_"med"(tilde(c)) < J_"med"(c^*)$, then $cal(E)_p (tilde(c)) = K p - J_"med"(tilde(c)) > K p - J_"med"(c^*) = cal(E)_p (c^*)$, contradicting global optimality of $c^*$ for $cal(E)_p$.
]

Compared to the spectral graph methods in Category V (`Algs. 103–105`), AP bypasses eigendecomposition completely and preserves exact non-metric pairwise affinities without requiring a positive-semidefinite Laplacian, trading spectral manifold unrolling for exemplar interpretability and automatic cardinality discovery.

*Empirical Benchmark Analysis:* On the 4-dataset ClusterOverview benchmark, Affinity Propagation achieves standout clustering quality on structured low-dimensional geometries while revealing its iterative runtime cost. On *Multi-Shape*, AP records the highest Silhouette of all four methods (`Sil = 0.727`, `ARI = N/A`), outperforming Ncut (`0.622`), NJW (`0.592`), and Ratio Cut (`0.390`) by placing compact, well-separated exemplars across each dense geometric lobe. On *Iris*, AP leads in external ground-truth recovery with `ARI = 0.449` (`Sil = 0.293`), surpassing Ncut (`ARI = 0.414`), NJW (`0.296`), and Ratio Cut (`0.114`). On *Customers* (`Sil = 0.261`), AP closely tracks Ncut (`0.269`), whereas on *High-D* (`Sil = 0.169`, `ARI = N/A`), Euclidean distance concentration narrows off-diagonal similarity contrasts ($s_"max" - s_"min"$), softening message margins relative to Ncut (`0.293`) and Ratio Cut (`0.219`) while still decisively beating NJW (`-0.028`). At a mean runtime of `13.438 ms`, AP is approximately $3 times$ slower than the spectral methods (`3.882–4.568 ms`), directly reflecting its damped iterative $O(T n^2)$ message-passing sweeps (Theorem 3) until fixed-point convergence.
