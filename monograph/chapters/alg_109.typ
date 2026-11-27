#import "../template/ieee.typ": *

#alg_entry(
  id: 109,
  name: "Markov Clustering (MCL)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2000",
  paradigm: "Stochastic Flow Simulation via Expansion-Inflation",
  complexity: $O(T n k^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The *Markov Clustering (MCL)* algorithm @vandongen2000mcl departs
fundamentally from cut-based and modularity objectives: instead of
optimizing a global partition functional, it simulates random-walk
flow directly on the graph and iteratively sharpens transition
probabilities until flow evaporates across sparse bottlenecks and
concentrates inside dense clusters.

Let $G = (V, E, w)$ be an undirected (or directed) weighted graph on
$n = |V|$ vertices with nonnegative adjacency matrix
$A in RR_(>= 0)^(n times n)$. To prevent bipartite periodicity and
stabilize short walks, MCL adds positive self-loops
$tilde(A) = A + alpha I$ ($alpha > 0$, commonly $alpha = 1$ or
$alpha_i = max_j A_(i j)$) and forms the column-stochastic Markov
transition matrix
$ M^((0)) = tilde(A) tilde(D)^(-1),
  quad tilde(D) = "diag"(tilde(d)_1, dots, tilde(d)_n),
  quad tilde(d)_j = sum_(i = 1)^n tilde(A)_(i j), $
so that $M_(i j)^((0)) = "Pr"(X_(t + 1) = i | X_t = j)$ and
$bold(1)^top M^((0)) = bold(1)^top$.

#definition(title: "Expansion and Hadamard Inflation Operators")[
  Let $cal(S)_n = {M in RR_(>= 0)^(n times n) : bold(1)^top M = bold(1)^top}$
  be the convex polytope of column-stochastic matrices. For integer
  *expansion order* $e >= 2$ (default $e = 2$) and real *inflation
  parameter* $r > 1$, define:
  $ cal(E)_e (M) &= M^e, \
    (Gamma_r (M))_(i j) &= (M_(i j))^r / (sum_(k = 1)^n (M_(k j))^r),
    quad i, j in {1, dots, n}. $
  The MCL dynamical system iterates
  $M^((t + 1)) = Gamma_r (cal(E)_e (M^((t))))$ until converging to a
  *doubly idempotent* limit $M^((infinity))$ satisfying both
  $(M^((infinity)))^2 = M^((infinity))$ and
  $Gamma_r (M^((infinity))) = M^((infinity))$.
]

Intuitively, expansion $M^e$ diffuses probability mass across $e$-step
paths — reinforcing flow inside clusters where many parallel short
paths exist — while inflation $Gamma_r$ boosts large column entries
and suppresses small inter-cluster leakage. At convergence, nonzero
rows of $M^((infinity))$ define *attractors*, and vertex $j$ is
assigned to the cluster of attractor $i$ whenever
$M_(i j)^((infinity)) > 0$.

#algbox(
  title: "Markov Clustering (MCL) With Sparse Column Pruning",
  inputs: [Adjacency $A in RR_(>= 0)^(n times n)$, self-loop $alpha > 0$,
    expansion $e >= 2$, inflation $r > 1$, column cap $k$, tolerance $epsilon$.],
  outputs: [Cluster partition (or overlapping attractor assignment) $cal(P)$.],
)[
  1. Form $tilde(A) <- A + alpha I$ and column-stochastic matrix
     $M^((0)) <- tilde(A) tilde(D)^(-1)$. Set $t <- 0$.
  2. *Expansion*: compute $W <- (M^((t)))^e$ via dense or sparse
     matrix multiplication.
  3. *Inflation*: compute $M_(i j)^((t + 1)) <- W_(i j)^r / sum_l W_(l j)^r$
     for all columns $j = 1, dots, n$.
  4. *Optional Pruning*: in each column $j$, retain the $k$ largest
     entries (or entries above threshold $tau$), zero the remainder,
     and renormalize column $j$ to unit $ell_1$-norm.
  5. Compute chaos metric
     $eta_t = max_j (max_i M_(i j)^((t + 1)) - sum_i (M_(i j)^((t + 1)))^2)$.
     If $eta_t < epsilon$ and $max_(i, j) |M_(i j)^((t + 1)) - M_(i j)^((t))| < epsilon$,
     halt; else set $t <- t + 1$ and return to Step 2.
  6. Extract attractor set $cal(A) = {i : exists j, M_(i j)^((t + 1)) > 0}$
     and interpret connected components of the directed attractor graph
     (nonzero entries of $M^((t + 1))$) as clusters.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish exponential-family entropy contraction and majorization
of inflation, a complete algebraic characterization of doubly
idempotent limits, superlinear local convergence rates, and exact
phase-transition thresholds on symmetric block graphs.

#subsec("A", "Entropy Contraction and Lorenz Majorization of Inflation")

#lemma(
  title: "Strict Entropy Contraction and Majorization of Inflation",
  origin: "Original Proof",
)[
  Let $p in Delta_(n - 1)$ be a probability column that is non-uniform
  on its support $S = "supp"(p)$ ($|S| >= 2$), and let
  $q = Gamma_r (p)$ with $r > 1$. Then:
  (i) $q$ strictly preserves ordering and amplifies odds ratios:
  $p_i > p_j > 0 ==> q_i / q_j = (p_i / p_j)^r > p_i / p_j$;
  (ii) $q$ strictly majorizes $p$ ($p prec q$);
  (iii) Shannon entropy strictly contracts, $H(Gamma_r (p)) < H(p)$,
  and collision probability strictly increases,
  $norm(Gamma_r (p))_2^2 > norm(p)_2^2$.
]

#proof[
  (i) For $i, j in S$,
  $q_i / q_j = (p_i^r / Z(r)) / (p_j^r / Z(r)) = (p_i / p_j)^r$,
  where $Z(r) = sum_(k in S) p_k^r$. Since $r > 1$ and $p_i / p_j > 1$,
  $(p_i / p_j)^r > p_i / p_j$.

  (ii) Sort entries as $p_((1)) >= dots >= p_((n))$. The ratio
  $q_((i)) / p_((i)) = p_((i))^(r - 1) / Z(r)$ is non-increasing on $S$.
  Hence there exists $k^*$ such that $q_((i)) >= p_((i))$ for
  $i <= k^*$ and $q_((i)) <= p_((i))$ for $i > k^*$. Partial sums
  $F_k = sum_(i = 1)^k (q_((i)) - p_((i)))$ therefore increase up to
  $k^*$ and decrease to $F_n = 0$, proving $F_k >= 0$ for all $k$ with
  strict inequality at $k^*$ (since $p$ is non-uniform on $S$).

  (iii) View $q(s)_i = exp(s ln p_i) / Z(s)$ ($i in S$) as an
  exponential family in natural parameter $s > 0$ with sufficient
  statistic $T_i = ln p_i$. Its log-partition function satisfies
  $(dif / (dif s)) ln Z(s) = EE_(q(s)) [T]$ and
  $(dif / (dif s)) EE_(q(s)) [T] = "Var"_(q(s)) (T)$. Differentiating
  the Shannon entropy $H(q(s)) = ln Z(s) - s EE_(q(s)) [T]$ gives
  $ (dif / (dif s)) H(q(s)) = EE_(q(s)) [T] - EE_(q(s)) [T]
    - s "Var"_(q(s)) (T) = - s "Var"_(q(s)) (ln p) < 0, $
  strictly negative for every $s > 0$ because $p$ is non-uniform on $S$.
  Integrating from $s = 1$ to $s = r$ gives $H(q(r)) < H(p)$.
  Similarly, $(dif q_i) / (dif s) = q_i (ln p_i - EE_(q(s)) [ln p])$,
  so
  $ (dif / (dif s)) norm(q(s))_2^2
    = 2 sum_(i in S) q_i^2 (ln p_i - EE_(q(s)) [ln p])
    = 2 "Cov"_(q(s)) (q(s), ln p) > 0 $
  by Chebyshev's sum inequality, since $q(s)_i$ is a strictly
  increasing, non-constant function of $ln p_i$ on $S$. Integrating
  from $1$ to $r$ yields $norm(q(r))_2^2 > norm(p)_2^2$.
]

#subsec("B", "Algebraic Structure of Doubly Idempotent Fixed Points")

#theorem(
  title: "Characterization of Doubly Idempotent Stochastic Matrices",
  origin: "Original Proof",
)[
  Let $M in cal(S)_n$ be doubly idempotent: $M^2 = M$ and
  $Gamma_r (M) = M$ for some $r > 1$. Denote column supports by
  $S_j = {i : M_(i j) > 0}$ and the attractor set by
  $cal(A) = {i : exists j, M_(i j) > 0}$. Then:
  (i) Every column $j$ is uniform on $S_j$: $M_(i j) = 1 / |S_j|$ for
  $i in S_j$.
  (ii) $cal(A)$ partitions into $K$ disjoint attractor equivalence
  classes $A_1, dots, A_K$ such that $S_i = A_c$ for all $i in A_c$,
  and each principal block $M_(A_c, A_c) = (1 / |A_c|) bold(1) bold(1)^top$
  is rank-$1$ uniform.
  (iii) The matrix rank equals the number of attractor classes:
  $"rank"(M) = "Tr"(M) = K$.
]

#proof[
  (i) By Lemma II.1(iii), $Gamma_r (M_(:, j)) = M_(:, j)$ forces
  $"Var" (ln M_(:, j)) = 0$ on $S_j$, so $M_(i j) = 1 / |S_j|$ for
  every $i in S_j$.

  (ii) Take any attractor $i in cal(A)$ and let
  $m_i = max_l M_(i l) > 0$, achieved at column $j^*$, so $i in S_(j^*)$.
  Because $M^2 = M$ and column $j^*$ sums to $1$,
  $ m_i = M_(i j^*) = sum_(k in S_(j^*)) M_(i k) M_(k j^*)
    <= sum_(k in S_(j^*)) m_i M_(k j^*) = m_i. $
  Since $M_(k j^*) > 0$ for all $k in S_(j^*)$, equality forces
  $M_(i k) = m_i$ for every $k in S_(j^*)$. In particular, setting
  $k = i in S_(j^*)$ yields $M_(i i) = m_i = 1 / |S_i| > 0$. Applying
  the same expansion to column $i$ gives
  $m_i = M_(i i) = sum_(k in S_i) M_(i k) M_(k i) <= m_i$, forcing
  $M_(i k) = m_i = 1 / |S_i|$ for every $k in S_i$. Hence every
  $k in S_i$ belongs to $cal(A)$, satisfies $i in S_k$, and obeys
  $1 / |S_k| = M_(i k) = 1 / |S_i|$, so $|S_k| = |S_i|$. Now for any
  $j in S_i$, $m_j = 1 / |S_j| = 1 / |S_i| = M_(j i)$, so column $i$
  achieves the row maximum $m_j$; repeating the equality step forces
  $M_(j k) = m_j = 1 / |S_i|$ for all $k in S_i$. Summing over the
  $|S_i|$ indices $j in S_i$ exhausts the unit sum of column $k$,
  proving $S_k = S_i$ for all $k in S_i$. Thus the relation
  $k ~ i <==> k in S_i$ is an equivalence relation on $cal(A)$,
  partitioning $cal(A)$ into classes $A_1, dots, A_K$ with
  $M_(A_c, A_c) = (1 / |A_c|) bold(1) bold(1)^top$.

  (iii) Because $M^2 = M$, every eigenvalue of $M$ lies in ${0, 1}$
  and $M$ is diagonalizable, so $"rank"(M) = "Tr"(M)$. Non-attractors
  satisfy $M_(i i) = 0$, whereas attractors $i in A_c$ satisfy
  $M_(i i) = 1 / |A_c|$. Thus
  $ "rank"(M) = sum_(c = 1)^K sum_(i in A_c) 1 / |A_c| = K. $
]

#subsec("C", "Superlinear Local Convergence Near Attractor Partitions")

#theorem(
  title: "Order-r Superlinear Local Convergence of MCL",
  origin: "Classical — Cited",
)[
  Let $M^* in cal(S)_n$ be a hard-partition doubly idempotent limit
  with unique attractor $a_c$ per cluster $C_c$ ($M_(a_c, j)^* = 1$
  for $j in C_c$). Define the maximum off-attractor column mass
  $delta(M) = max_(1 <= j <= n) sum_(i != a_(c(j))) M_(i j)$. Under
  MCL with $e = 2$ and $r > 1$ @vandongen2000mcl, every $M in cal(S)_n$
  with $delta(M) = delta < 1 / 6$ satisfies
  $ delta(Gamma_r (M^2)) <= ((2 delta) / (1 - delta)^2)^r < (3 delta)^r. $
  In particular, default $r = 2$ converges quadratically in a
  neighborhood of every hard-partition attractor state.
]

#proof[
  Fix any column $j in C_c$ with attractor $a = a_c$. By hypothesis,
  $M_(a, j) >= 1 - delta$ and $M_(a, a) >= 1 - delta$. Thus the
  expanded attractor entry satisfies
  $ (M^2)_(a, j) = sum_(k = 1)^n M_(a, k) M_(k, j)
    >= M_(a, a) M_(a, j) >= (1 - delta)^2. $
  Because $M^2 in cal(S)_n$, total off-attractor mass in column $j$
  of $M^2$ obeys
  $ sum_(i != a) (M^2)_(i, j) = 1 - (M^2)_(a, j)
    <= 1 - (1 - delta)^2 = 2 delta - delta^2 <= 2 delta. $
  Applying $Gamma_r$ ($r > 1$) and using the superadditivity inequality
  $sum_(i != a) x_i^r <= (sum_(i != a) x_i)^r$ for nonnegative $x_i$
  gives
  $ sum_(i != a) (Gamma_r (M^2))_(i, j)
    &= (sum_(i != a) ((M^2)_(i, j))^r)
    / (((M^2)_(a, j))^r + sum_(i != a) ((M^2)_(i, j))^r) \
    &<= ((sum_(i != a) (M^2)_(i, j))^r) / ((M^2)_(a, j))^r
    <= ((2 delta) / (1 - delta)^2)^r. $
  For $delta < 1 / 6$, $(1 - delta)^2 > 25 / 36 > 2 / 3$, so
  $2 delta / (1 - delta)^2 < 3 delta < 1 / 2$, proving order-$r$
  superlinear convergence to $M^*$.
]

#subsec("D", "Bipartite Oscillation and Symmetric Block Phase Transition")

#proposition(
  title: "Necessity of Self-Loops and Two-Block Resolution Threshold",
  origin: "Original Derivation",
)[
  (i) On any connected bipartite graph $G = (V_1 union V_2, E)$ with
  self-loop weight $alpha = 0$, $M^2$ is block-diagonal across
  $(V_1, V_2)$, permanently decoupling the two parts regardless of
  connectivity, while odd $e$ induces period-$2$ oscillation. Adding
  $alpha > 0$ breaks parity decoupling.
  (ii) On a symmetric two-block lumped chain
  $P = mat(1 - p, p; p, 1 - p)$ with inter-block escape probability
  $p in (0, 1 / 2)$ and $e = 2$, MCL separates the two blocks
  ($p^((t)) -> 0$) iff $r$ exceeds the exact critical threshold
  $ r_c (p) = (ln (p / (1 - p)))
    / (ln ((2 p (1 - p)) / (1 - 2 p + 2 p^2))), $
  and merges them ($p^((t)) -> 1 / 2$) for $1 < r < r_c (p)$.
]

#proof[
  (i) With $alpha = 0$, $M = mat(0, P_(1 2); P_(2 1), 0)$, so
  $M^2 = mat(P_(1 2) P_(2 1), 0; 0, P_(2 1) P_(1 2))$ preserves zero
  off-diagonal blocks under $Gamma_r$, trapping flow inside $V_1$ and
  $V_2$ separately; odd powers map $V_1 <-> V_2$ indefinitely.
  (ii) Squaring $P$ yields $P^2 = mat(1 - q, q; q, 1 - q)$ with
  $q = 2 p (1 - p)$ and $1 - q = 1 - 2 p + 2 p^2$. Denote the
  expanded off-to-diagonal ratio by
  $rho(p) = q / (1 - q) = (2 p (1 - p)) / (1 - 2 p + 2 p^2) in (0, 1)$.
  Applying $Gamma_r$ maps the odds ratio $p / (1 - p)$ to
  $p' / (1 - p') = (rho(p))^r$. Thus $p' < p$ (strict contraction
  toward $0$, separating the blocks) holds iff
  $(rho(p))^r < p / (1 - p)$. Taking logarithms (both sides negative
  for $p < 1 / 2$) reverses the inequality to $r > r_c (p)$.
  As $p -> 0^+$, $r_c (p) -> 1^+$; as $p -> (1 / 2)^-$ with
  $p = 1 / 2 - epsilon$, $r_c (p) = 1 / (2 epsilon) + O(1) -> infinity$.
  For instance, at $p = 0.25$, $r_c (0.25) = ln 3 / ln (5 / 3) approx 2.151$,
  so $r = 2.0$ merges the blocks whereas $r = 2.2$ separates them.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact dense MCL requires $O(n^3)$ time per expansion step and
$O(n^2)$ memory, becoming prohibitive beyond $n approx 5000$ vertices
and accounting for its higher wall-clock latency relative to linear
graph methods on small dense matrices. In practice, sparse column
pruning (Step 4) enforces at most $k$ nonzeros per column (typically
$k in [50, 1000]$), reducing sparse matrix multiplication $M^2$ to
$O(n k^2)$ time and $O(n k)$ space @vandongen2000mcl. Because local
convergence is quadratic at $r = 2$ (Theorem II.3), MCL typically
converges in $T = 10$ to $30$ iterations.

Geometrically, MCL operates purely via local random-walk transition
probabilities without any global normalization $m_w$: it is immune to
modularity's global resolution limit, scale-invariant to uniform
positive edge-weight multiplication (when $alpha$ scales with edge
weights), and permutation-equivariant. However, because expansion
explores finite-horizon walks before inflation prunes weak tails, MCL
favors compact, small-diameter clusters and tends to subdivide
high-diameter chains or elongated filaments into segments whose
diameter scales inversely with $r$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(T n k^2)$ pruned sparse; $O(T n^3)$ exact dense.],
    [Space complexity], [$O(n k)$ pruned sparse; $O(n^2)$ exact dense.],
    [Convergence rate], [Order-$r$ superlinear locally (quadratic at $r = 2$, Thm. II.3).],
    [Hyperparameters], [Inflation $r > 1$ (granularity), expansion $e >= 2$, loop $alpha$, cap $k$.],
    [Determinism], [Fully deterministic (no random initialization or sweep order).],
    [Invariances], [Positive weight scale-invariant (with scaled $alpha$); perm-equivariant.],
    [Overlap capability], [Naturally detects symmetric boundary overlaps ($|S_j| >= 2$).],
    [Structural bias], [Local-flow bias; no global resolution limit; splits long chains.],
  ),
  caption: [Operational and theoretical profile of Markov Clustering (MCL).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Protein Sequence Homology Clustering (TribeMCL & OrthoMCL):*
   MCL's flagship domain clusters millions of proteins from all-vs-all
   BLAST $E$-value similarity graphs. Multi-domain proteins act as
   promiscuous bridges between distinct protein families; MCL's
   nonlinear Hadamard inflation suppresses single-domain cross-talk
   while preserving densely connected paralog and ortholog families
   without suffering modularity's resolution limit on tiny gene families.
2. *Protein Complex Identification in PPI Networks:* Physical
   interaction assays (yeast two-hybrid, AP-MS) exhibit substantial
   false-positive edges. Flow-based random walks amplify multi-path
   co-complex cliques while inflation filters isolated spurious links.
3. *High-Dimensional k-NN Graph Clustering:* In high-dimensional
   feature spaces, hub points accumulate spurious cross-cluster k-NN
   edges that distort degree-based modularity null models. MCL's
   column normalization down-weights outgoing hub probability mass as
   $1 / tilde(d)_j$, and inflation suppresses weak cross-cluster flow.
4. *Semantic Word-Sense and Co-Occurrence Graph Clustering:* Word
   co-occurrence graphs contain polysemous hubs linking unrelated
   semantic domains; MCL isolates distinct lexical senses and naturally
   exposes symmetric overlapping attractors when a word bridges two
   equally strong contexts.
5. *Malware Call-Graph and Network-Flow Anomaly Grouping:* Clustering
   communication or API-call similarity graphs benefits from MCL's
   complete determinism, absence of a prescribed cluster count $K$,
   and tunable inflation knob $r$ matching desired granularity.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

MCL encounters four characteristic failure modes: (a) *high-diameter
filament fragmentation*, where elongated chains or 2D manifold ribbons
exceed the effective random-walk mixing length and break into multiple
short segments; (b) *computational bottleneck* at $O(n^3)$ dense or
$O(n k^2)$ pruned when graphs are dense or lack aggressive pruning;
(c) *sensitivity to $r$ and $alpha$*, as small changes near $r_c (p)$
(Proposition II.4) can trigger sharp cluster merges or splits, while
missing self-loops ($alpha = 0$) break bipartite subgraphs; and
(d) *lack of native multilevel hierarchy* in a single run, requiring
separate executions across a grid of inflation values $r$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Protein homology (BLAST)], [Ideal: hub-damped flow, tiny families intact], [Prune columns ($k approx 200$) for scale],
    [PPI complex discovery], [High: filters sparse noisy inter-links], [Tune $r in [1.4, 3.0]$ for granularity],
    [High-D k-NN graphs], [High: suppresses hub leakage via $tilde(D)^(-1)$], [Expensive dense $M^2$; use sparse pruning],
    [Elongated / chain manifolds], [Poor: finite walk horizon chops filaments], [Use HDBSCAN @campello2013hdbscan or LPA],
    [Billion-edge social graphs], [Poor: $O(n k^2)$ matrix fill-in bottleneck], [Use Leiden (Alg. 108) or LPA (Alg. 110)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for MCL.],
)

#sec("V", "Structural Relations and Empirical Behavior")

MCL belongs to the diffusion-based family alongside Walktrap
@pons2006walktrap, Infomap @rosvall2008infomap, Power Iteration
Clustering @lin2010pic, and Spectral Clustering @shi2000ncut
@vonluxburg2007tutorial. The following spectral similarity links
MCL's linear expansion step directly to the normalized graph spectrum
and quantifies how self-loop weight $alpha$ eliminates negative
bipartite eigenvalues.

#proposition(
  title: "Spectral Similarity of Expansion and Self-Loop Eigenvalue Floor",
  origin: "Original Derivation",
)[
  For undirected symmetric $A >= 0$ with $tilde(A) = A + alpha I$
  ($alpha > 0$) and $d_("max") = max_i sum_j A_(i j)$, the expansion
  operator $M^e$ is similar to the $e$-th power of the symmetric
  normalized matrix $S = tilde(D)^(-1 / 2) tilde(A) tilde(D)^(-1 / 2)$:
  $ M^e = tilde(D)^(1 / 2) S^e tilde(D)^(-1 / 2), $
  whose spectrum lies strictly inside $(-1, 1]$ with explicit lower
  bound
  $ lambda_("min") (S) >= -1 + (2 alpha) / (alpha + d_("max")). $
  In particular, choosing $alpha >= d_("max")$ guarantees
  $S succ.eq 0$ (all eigenvalues nonnegative).
]

#proof[
  Since $M = tilde(A) tilde(D)^(-1) = tilde(D)^(1 / 2) S tilde(D)^(-1 / 2)$,
  telescoping gives $M^e = tilde(D)^(1 / 2) S^e tilde(D)^(-1 / 2)$.
  For any nonzero $x in RR^n$, note that
  $x^top (D + A) x = sum_(i, j) A_(i j) (x_i + x_j)^2 / 2 >= 0$, so
  $x^top A x >= - x^top D x$, where $D = "diag"(d_1, dots, d_n)$.
  Therefore
  $ x^top tilde(A) x &= alpha norm(x)_2^2 + x^top A x
    >= sum_(i = 1)^n (alpha - d_i) x_i^2 \
    &= sum_(i = 1)^n ((alpha - d_i) / (alpha + d_i))
    (alpha + d_i) x_i^2
    >= ((alpha - d_("max")) / (alpha + d_("max"))) x^top tilde(D) x, $
  since $t |-> (alpha - t) / (alpha + t)$ is decreasing on $t >= 0$.
  Substituting $y = tilde(D)^(1 / 2) x$ gives
  $(y^top S y) / norm(y)_2^2 >= (alpha - d_("max")) / (alpha + d_("max"))
  = -1 + 2 alpha / (alpha + d_("max")) > -1$, and $S succ.eq 0$ when
  $alpha >= d_("max")$.
]

Whereas Spectral Clustering @vonluxburg2007tutorial extracts the top
eigenvectors of $S$ and clusters in Euclidean embedding space, and
Walktrap @pons2006walktrap computes fixed-$t$ transition distances
$norm(D^(-1 / 2) (P^t_i - P^t_j))_2$, MCL interleaves linear
spectral damping ($S^e$) with nonlinear Hadamard inflation ($Gamma_r$)
to collapse the matrix directly into sparse rank-$K$ attractors
(Theorem II.2) without ever forming eigenvectors or choosing $K$.

*Empirical Behavior (ClusterOverview 4-Dataset Benchmark).* On the
ClusterOverview 4-dataset benchmark (Silhouette
@rousseeuw1987silhouettes, ARI @hubert1985comparing), MCL achieves
*Iris* Sil = $0.323$, ARI = $0.448$; *Customers* Sil = $0.430$;
*Multi-Shape* Sil = $0.534$ (ARI = N/A); and *High-D* Sil = *$0.321$*
(ARI = N/A), with a mean runtime of *$11.375$ ms*. Notably, MCL posts
the *highest High-D Silhouette* ($0.321$) among all four Cat. V graph
algorithms in this group (surpassing Louvain's $0.232$, Leiden's
$0.303$, and LPA's $0.303$), substantially beats modularity methods on
*Iris* ARI ($0.448$ vs $0.333$ for Louvain and $0.340$ for Leiden),
and edges ahead on *Multi-Shape* ($0.534$ vs $0.526$). This reflects
the theoretical advantages established in Section II: column
normalization and Hadamard inflation suppress high-dimensional hub
leakage and escape modularity's resolution limit, at the predicted
cost of dense $O(n^3)$ matrix squaring ($11.375$ ms — roughly
$16 times$ to $26 times$ slower than Leiden, LPA, and Louvain).
