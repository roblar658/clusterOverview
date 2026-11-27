#import "../template/ieee.typ": *

#alg_entry(
  id: 198,
  name: "EVO-KMeans / Genetic K-Means (GKA)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1999",
  paradigm: "Hybrid Evolutionary Search with One-Step K-Means Operator",
  complexity: $cal(O)(G P n K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Given a dataset $X = {bold(x)_1, dots.h, bold(x)_n} subset bb(R)^d$ and integer
$K in {2, dots.h, n - 1}$, minimum Sum-of-Squared-Errors (SSE) clustering seeks
a hard partition encoding $bold(s) = (s_1, dots.h, s_n) in cal(S)_K = {1,
dots.h, K}^n$ minimizing the Total Within-Cluster Variation (TWCV):
$
  J(bold(s)) = sum_(k=1)^K sum_(i : s_i = k) norm(bold(x)_i - bold(c)_k (bold(s)))_2^2,
  wide
  bold(c)_k (bold(s)) = 1 / (|C_k (bold(s))|) sum_(i in C_k (bold(s))) bold(x)_i,
$
where $C_k (bold(s)) = {i : s_i = k}$ (requiring $|C_k (bold(s))| >= 1$ for
legal strings in $cal(S)_K^+$). Minimizing $J(bold(s))$ is NP-hard even for
$K = 2$ in general dimension @aloise2009np or in the Euclidean plane $d = 2$
@mahajan2012planar. Standard Lloyd's $K$-means @lloyd1982least performs greedy
coordinate descent that can take super-polynomial iterations @vattani2011kmeans
and stalls in poor local optima sensitive to initialization. Genetic $K$-Means
(GKA / EVO-KMeans) embeds a problem-specific one-step $K$-means local operator
inside an elitist genetic algorithm to guarantee global asymptotic convergence
while preserving fast local exploitation @krishna1999gka.

#definition(title: "GKA Encoding, Scaled Fitness, and Hybrid Operators")[
  A population $Z_t = (bold(s)^((1)), dots.h, bold(s)^((P))) in (cal(S)_K^+)^P$
  evolves through four stochastic operators at generation $t$:
  - *Sigma-truncated fitness*: Set raw score $g(bold(s)) = -J(bold(s))$ with
    population mean $overline(g)$ and standard deviation $sigma_g$, defining
    $
      F(bold(s)) = cases(
        max(0, g(bold(s)) - (overline(g) - c_s sigma_g)) &"if" sigma_g > 0,
        1 &"if" sigma_g = 0,
      )
    $
    for scaling constant $c_s in [1, 3]$.
  - *Elitist proportional selection*: Copy the best chromosome $bold(s)^(ast)$
    unchanged and sample remaining parents with probabilities proportional to
    $F(bold(s))$.
  - *Distance-biased mutation* ($cal(M)_(p_m)$): Independently with probability
    $p_m in (0, 1)$ per gene $i$, reassign $s_i$ to cluster $k in {1, dots.h, K}$
    with probability
    $
      p_(i, k) (bold(s)) =
      (c_m d_"max" (bold(x)_i) - d(bold(x)_i, bold(c)_k (bold(s)))) /
      (sum_(j=1)^K (c_m d_"max" (bold(x)_i) - d(bold(x)_i, bold(c)_j (bold(s))))),
    $
    where $d(bold(x)_i, bold(c)_k) = norm(bold(x)_i - bold(c)_k)_2$,
    $d_"max" (bold(x)_i) = max_j d(bold(x)_i, bold(c)_j)$, and $c_m > 1$.
  - *One-step $K$-Means Operator (KMO)* ($cal(K)$): Replace each chromosome
    $bold(s)$ by $cal(K)(bold(s))$, computing centroids $bold(c)_k (bold(s))$
    and reassigning every $i$ to $arg min_k norm(bold(x)_i - bold(c)_k (bold(s)))_2$,
    repairing any empty cluster via farthest-point seeding. (EVO-KMeans
    variants optionally pair KMO with label-aligned crossover before mutation.)
]

Notice the pivotal mathematical design insight of Krishna and Murty
@krishna1999gka: naive crossover over partition strings suffers severe
$K!$-fold label-permutation disruption, so canonical GKA replaces blind
crossover with KMO as a domain-informed search operator while retaining
strictly positive distance-biased mutation to guarantee ergodicity.

#algbox(
  title: "EVO-KMeans / Genetic K-Means Algorithm (GKA)",
  inputs: [Dataset $X subset bb(R)^d$, clusters $K$, population $P$,
    generations $G$, mutation rate $p_m in (0, 1)$, bias $c_m > 1$,
    sigma factor $c_s in [1, 3]$.],
  outputs: [Best legal partition $bold(s)^(ast) in cal(S)_K^+$ and centroids
    ${bold(c)_k^(ast)}_(k=1)^K$.],
)[
  + Initialize $P$ legal chromosomes $Z_0 = (bold(s)^((1)), dots.h,
    bold(s)^((P)))$ randomly (or via $K$-means++ seeding @arthur2007kmeanspp)
    and apply one KMO pass to ensure all clusters are non-empty.
  + Record global elite $bold(s)_"best" <- arg min_(bold(s) in Z_0) J(bold(s))$.
  + *For* generation $t = 1, dots.h, G$:
    - Evaluate $J(bold(s)^((p)))$ and $sigma$-truncated fitness
      $F(bold(s)^((p)))$ for $p = 1, dots.h, P$.
    - *Selection*: Sample $P$ strings via roulette-wheel selection with
      probabilities $F(bold(s)^((p))) / sum_q F(bold(s)^((q)))$.
    - *Optional aligned crossover (EVO-KMeans)*: Pair selected parents, align
      centroid labels via minimum-weight bipartite matching, and recombine.
    - *Biased mutation*: For each string and each point $i in {1, dots.h, n}$,
      with probability $p_m$ sample new label $s_i ~ p_(i, dot) (bold(s))$ and
      repair empty clusters if any arise.
    - *KMO local step*: Update every string via $bold(s)^((p)) <- cal(K)(bold(s)^((p)))$.
    - *Elitism*: Update $bold(s)_"best"$ if improved; replace the worst string
      in $Z_t$ with $bold(s)_"best"$.
  + Return $bold(s)_"best"$ and its exact cluster centroids.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical results: KMO Lyapunov descent and
fixed-point optimality, almost-sure global convergence via finite Markov chains,
a constructive exponential worst-case basin-escape lower bound, an algebraic
proof of $K!$-symmetry crossover disruption, and a finite-generation tail bound.

#subsec("A", "Monotonic KMO Descent and Global Fixed-Point Invariance")

Before analyzing stochastic global convergence, we prove that the deterministic
local operator $cal(K)$ never increases SSE and preserves every global optimum.

#lemma(
  title: "Strict KMO Descent and Global-Optimum Fixed Points",
  origin: "Original Proof",
)[
  For every legal partition $bold(s) in cal(S)_K^+$ whose KMO image
  $bold(s)' = cal(K)(bold(s))$ has no empty clusters:
  (i) $J(cal(K)(bold(s))) <= J(bold(s))$, with equality if and only if
  $bold(s)$ is a Voronoi-stationary fixed point of Lloyd's algorithm.
  (ii) Every global minimizer $bold(s)^(ast) in arg min_(bold(s) in cal(S)_K^+)
  J(bold(s))$ is a fixed point of KMO up to ties on Voronoi bisectors, i.e.,
  $J(cal(K)(bold(s)^(ast))) = J(bold(s)^(ast))$.
]

#proof[
  (i) Let $bold(c)_k = bold(c)_k (bold(s))$ and $bold(c)'_k = bold(c)_k (bold(s)')$
  be the centroids before and after KMO. By definition of Voronoi reassignment,
  for every $i in {1, dots.h, n}$, $norm(bold(x)_i - bold(c)_(s'_i))_2^2 <=
  norm(bold(x)_i - bold(c)_(s_i))_2^2$. Summing over $i = 1, dots.h, n$ gives
  $
    sum_(k=1)^K sum_(i in C_k (bold(s)')) norm(bold(x)_i - bold(c)_k)_2^2
    <= sum_(k=1)^K sum_(i in C_k (bold(s))) norm(bold(x)_i - bold(c)_k)_2^2
    = J(bold(s)).
  $
  Next, applying the parallel-axis (bias-variance) identity inside each new
  cluster $C'_k = C_k (bold(s)')$ with mean $bold(c)'_k$ yields
  $
    sum_(i in C'_k) norm(bold(x)_i - bold(c)_k)_2^2
    = sum_(i in C'_k) norm(bold(x)_i - bold(c)'_k)_2^2
      + |C'_k| norm(bold(c)'_k - bold(c)_k)_2^2.
  $
  Summing over $k = 1, dots.h, K$ establishes the exact descent gap:
  $
    J(bold(s)) - J(bold(s)')
    >= sum_(k=1)^K |C'_k| norm(bold(c)'_k - bold(c)_k)_2^2 >= 0.
  $
  Equality holds iff every point already sits in a nearest-centroid cell AND
  $bold(c)'_k = bold(c)_k$ for all $k$, i.e., $bold(s)$ is Voronoi-stationary.

  (ii) Suppose $bold(s)^(ast)$ is a global minimizer. By optimality,
  $J(cal(K)(bold(s)^(ast))) >= J(bold(s)^(ast))$, while part (i) yields
  $J(cal(K)(bold(s)^(ast))) <= J(bold(s)^(ast))$ (no cluster can empty under
  KMO at a global optimum when $X$ has at least $K$ distinct points, since
  reassigning a point from a multi-point cluster to an empty center strictly
  reduces SSE). Hence $J(cal(K)(bold(s)^(ast))) = J(bold(s)^(ast))$.
]

#subsec("B", "Almost-Sure Global Convergence via Markov Chains")

Unlike standard $K$-means, elitist GKA converges with probability $1$ to the
exact NP-hard global SSE optimum in the infinite-generation limit.

#theorem(
  title: "Global Asymptotic Convergence of Elitist GKA",
  origin: "Classical — Cited",
)[
  Let $Omega = (cal(S)_K^+)^P$ be the finite population state space of elitist
  GKA with $p_m in (0, 1)$ and $c_m > 1$, and let
  $Omega^(ast) = {Z in Omega : min_(bold(s) in Z) J(bold(s)) = J^(ast)}$ be the
  set of populations containing a global minimizer. Then for any initial
  distribution over $Omega$, the sequence of best objectives
  $J_t^(ast) = min_(bold(s) in Z_t) J(bold(s))$ satisfies
  $
    lim_(t -> infinity) "Pr"(Z_t in Omega^(ast)) = 1
    wide "and" wide
    J_t^(ast) -> J^(ast) "almost surely"
  $
  as established in @krishna1999gka.
]

#proof[
  Because $X$ is finite and $cal(S)_K^+$ has finite cardinality
  $|cal(S)_K^+| <= K^n$, the population sequence ${Z_t}_(t >= 0)$ is a finite
  time-homogeneous Markov chain with transition matrix $Q = S M K_"op" E$,
  corresponding to selection, mutation, KMO, and elitism.

  *Step 1 (Strict Positivity of Mutation).* For every point $bold(x)_i$ and
  cluster $k in {1, dots.h, K}$, $d(bold(x)_i, bold(c)_k) <= d_"max" (bold(x)_i)$.
  Because $c_m > 1$, the numerator satisfies
  $c_m d_"max" (bold(x)_i) - d(bold(x)_i, bold(c)_k) >= (c_m - 1) d_"max" (bold(x)_i)$.
  If $d_"max" (bold(x)_i) > 0$, then strictly
  $p_(i, k) (bold(s)) >= (c_m - 1) / (K c_m) =: p_"min" > 0$ (and if
  $d_"max" (bold(x)_i) = 0$, $p_(i, k) = 1 / K >= p_"min"$). Since each allele
  mutates independently with probability $p_m in (0, 1)$, the probability of
  mutating any legal string $bold(s)$ directly into a fixed global minimizer
  $bold(s)^(ast)$ is lower-bounded by
  $
    M_(bold(s), bold(s)^(ast))
    >= (min(1 - p_m, p_m p_"min"))^n =: delta_0 > 0,
  $
  uniformly over all $bold(s) in cal(S)_K^+$.

  *Step 2 (Reachability of $Omega^(ast)$ Through KMO).* By Lemma 1(ii),
  $cal(K)(bold(s)^(ast))$ is itself a global minimizer with objective
  $J^(ast)$. Therefore, from *any* population $Z in Omega$, after selection
  produces $P$ intermediate strings, the probability that at least one string
  mutates to $bold(s)^(ast)$ and survives KMO as a global optimum is at least
  $
    "Pr"(Z_(t+1) in Omega^(ast) | Z_t = Z) >= 1 - (1 - delta_0)^P =: delta_P > 0.
  $

  *Step 3 (Absorbing Closure and Almost-Sure Convergence).* Elitism enforces
  $J_(t+1)^(ast) <= J_t^(ast)$ deterministically, so once $Z_t in Omega^(ast)$,
  $Z_s in Omega^(ast)$ for all $s >= t$ ($Omega^(ast)$ is a closed absorbing
  class). By induction on $t$,
  $
    "Pr"(Z_t in.not Omega^(ast)) <= (1 - delta_P)^t -> 0
    quad "as" t -> infinity.
  $
  Furthermore, $sum_(t=0)^infinity "Pr"(J_t^(ast) > J^(ast)) <=
  sum_(t=0)^infinity (1 - delta_P)^t = 1 / delta_P < infinity$. By the first
  Borel-Cantelli lemma, $"Pr"(J_t^(ast) > J^(ast) "i.o.") = 0$, proving
  $J_t^(ast) -> J^(ast)$ almost surely.
]

#subsec("C", "Exponential Basin-Escape Hitting Time Lower Bound")

Theorem 2 guarantees asymptotic optimality, yet NP-hardness implies that worst-case
finite-time convergence cannot be polynomial. We exhibit a sharp lower bound
quantifying the basin-escape barrier of KMO-coupled mutation.

#theorem(
  title: "Worst-Case Basin-Transition Hitting Time Barrier",
  origin: "Original Proof",
)[
  Let $K = 3$ and $n = 4r$ for integer $r >= 4$. There exists a 1D dataset $X$
  and a suboptimal Voronoi-stationary partition $bold(s)_"loc"$ such that from
  any population initialized inside $bold(s)_"loc"$, escaping the basin of
  attraction of $bold(s)_"loc"$ requires mutating at least $r$ genes in a
  single generation, yielding an expected hitting time to the global optimum
  $tau^(ast) = inf {t >= 0 : Z_t in Omega^(ast)}$ lower-bounded by
  $
    bb(E)[tau^(ast)] >= 1 / (P binom(4r, r) p_m^r) >= 1 / (P (4 e p_m)^r).
  $
  In particular, for any mutation rate $p_m < 1 / (4 e)$, $bb(E)[tau^(ast)]$
  grows exponentially in the cluster size $r$.
]

#proof[
  Place four tight point groups on $bb(R)$ of size $r$ each: $G_1$ at $x = 0$,
  $G_2$ at $x = 1$, $G_3$ at $x = 100 r$, and $G_4$ at $x = 100 r + epsilon$
  with $epsilon = 10^(-3)$. The unique global 3-clustering separates $G_1$,
  $G_2$, and $G_3 union G_4$ with $J^(ast) = (r / 2) epsilon^2$. Consider the
  suboptimal Voronoi-stationary partition $bold(s)_"loc"$ that merges
  $C_1 = G_1 union G_2$ (centroid $0.5$) while splitting $C_2 = G_3$ (centroid
  $100 r$) and $C_3 = G_4$ (centroid $100 r + epsilon$), yielding
  $J(bold(s)_"loc") = r / 2 >> J^(ast)$.

  Consider any mutated offspring $tilde(bold(s))$ obtained from $bold(s)_"loc"$
  in which strictly fewer than $r$ genes mutate across all $n = 4r$ points.
  Let $a$ be the number of points leaving $C_2$ and $b$ the number of points
  from $G_1 union G_2$ entering $C_2$; since these are disjoint subsets of
  mutated points, $a + b <= r - 1$, so $r - a >= b + 1 > b$. Consequently, the
  mutated centroid $tilde(c)_2$ satisfies
  $
    tilde(c)_2 >= (100 r (r - a)) / (r - a + b)
    >= (100 r (b + 1)) / (2b + 1) > 50 r.
  $
  By the identical argument, $tilde(c)_3 > 50 r$. Meanwhile, let $k_1$ points
  leave $C_1$ and $k_2$ points from $G_3 union G_4$ enter $C_1$, with
  $k_1 + k_2 <= r - 1$. Then $2r - k_1 >= r + 1 + k_2 > 2 k_2$, so
  $k_2 / (2r - k_1 + k_2) < 1 / 3$ and
  $
    tilde(c)_1 <= (2r - k_1 + k_2 (100 r + epsilon)) / (2r - k_1 + k_2)
    < (100 r + epsilon) / 3 + 1 < 35 r.
  $
  Because $0 <= 1 < tilde(c)_1 < 35 r < 50 r < min(tilde(c)_2, tilde(c)_3)$,
  every point in $G_1 union G_2$ is strictly closer to $tilde(c)_1$ than to
  $tilde(c)_2$ or $tilde(c)_3$, while every point in $G_3 union G_4$ sits at
  distance $> 65 r$ from $tilde(c)_1$ and $< 50 r + epsilon$ from
  $tilde(c)_2, tilde(c)_3$. Hence applying KMO restores $C_1 = G_1 union G_2$
  and $C_2 union C_3 = G_3 union G_4$ identically, pulling the chromosome back
  into the suboptimal basin! Thus escaping requires at least $r$ simultaneous
  mutations in a single chromosome. By the binomial union bound over $P$
  individuals and Stirling's bound $binom(4r, r) <= (4 e)^r$, the per-generation
  escape probability satisfies $p_"esc" <= P binom(4r, r) p_m^r <= P (4 e p_m)^r$,
  giving $bb(E)[tau^(ast)] >= 1 / p_"esc" >= (P (4 e p_m)^r)^(-1)$.
]

#subsec("D", "Label-Permutation Degeneracy and Crossover Disruption")

Why did canonical GKA drop standard crossover in favor of KMO, and why do
EVO-KMeans variants require label alignment? We prove this algebraically.

#proposition(
  title: "Symmetry-Induced Crossover Disruption and Alignment Cure",
  origin: "Original Derivation",
)[
  Let $S_K$ be the symmetric group of $K!$ cluster-label permutations acting on
  $cal(S)_K^+$ by $(pi dot bold(s))_i = pi(s_i)$.
  (i) $J(pi dot bold(s)) = J(bold(s))$ for all $pi in S_K$ and all $bold(s)$.
  (ii) Let $bold(s)^(ast)$ be a globally optimal partition of $K$ well-separated
  clusters of radius $<= r$ and pairwise centroid separation $>= D > 4r$, and
  let $bold(s)' = pi dot bold(s)^(ast)$ for any derangement $pi in S_K$
  ($pi(k) != k$ for all $k$). Then uniform or balanced single-point crossover
  between the *identical* optimal solutions $bold(s)^(ast)$ and $bold(s)'$
  produces an offspring $bold(u)$ with strictly degraded SSE
  $J(bold(u)) >= J(bold(s)^(ast)) + Omega(n D^2)$.
  (iii) Solving the $cal(O)(K^3)$ linear assignment problem
  $hat(pi) = arg min_(pi in S_K) sum_(k=1)^K norm(bold(c)_k (bold(s)) -
  bold(c)_(pi(k)) (bold(s)'))_2^2$ prior to crossover recovers exact
  idempotence ($bold(u) = bold(s)^(ast)$) whenever $bold(s)' in S_K dot bold(s)^(ast)$.
]

#proof[
  (i) Permuting cluster indices merely reorders the outer summation over
  $k = 1, dots.h, K$ in $J(bold(s))$, leaving centroids and total SSE invariant.
  (ii) Suppose crossover selects prefix $A$ (containing half of each cluster
  $C_k^(ast)$) from $bold(s)^(ast)$ and suffix $B$ from $pi dot bold(s)^(ast)$.
  For each label $k$, offspring cluster $C_k (bold(u))$ contains half of
  $C_k^(ast)$ (near $bold(c)_k^(ast)$) merged with half of $C_(pi^(-1)(k))^(ast)$
  (near $bold(c)_(pi^(-1)(k))^(ast)$). Because $pi$ is a derangement,
  $pi^(-1)(k) != k$, so the two halves sit at distance $>= D - 2r > D/2$.
  By the parallel-axis theorem, each offspring cluster incurs extra variance at
  least $(|C_k|/4) (D/2)^2 = |C_k| D^2 / 16$, summing to
  $J(bold(u)) - J(bold(s)^(ast)) >= (n / 16) D^2$.
  (iii) When $bold(s)' = pi dot bold(s)^(ast)$ with distinct optimal centroids,
  the cost matrix $M_(j, k) = norm(bold(c)_j (bold(s)) - bold(c)_k (bold(s)'))_2^2$
  has unique zero-cost permutation $hat(pi) = pi^(-1)$, so relabeling $bold(s)'$
  by $hat(pi)$ yields $hat(pi) dot bold(s)' = bold(s)^(ast)$ identically, making
  crossover completely non-destructive on equivalent optimal parents.
]

#subsec("E", "Finite-Generation Tail Bound and Population Amplification")

#corollary(
  title: "Geometric Tail Decay Under Population Parallelism",
  origin: "Original Derivation",
)[
  For any tolerance $epsilon > 0$, let $q_epsilon > 0$ be the minimum one-step
  probability that a single mutated-and-KMO chromosome enters the
  $epsilon$-optimal level set $cal(S)_epsilon = {bold(s) : J(bold(s)) <= J^(ast)
  + epsilon}$. Then after $G$ generations with elitist population size $P$,
  $
    "Pr"(J_G^(ast) > J^(ast) + epsilon) <= (1 - q_epsilon)^(P G)
    <= exp(-P G q_epsilon).
  $
]

#proof[
  Conditioned on $Z_(t-1)$, the $P$ individuals undergo independent mutation and
  KMO transitions, so the probability that all $P$ individuals miss
  $cal(S)_epsilon$ is at most $(1 - q_epsilon)^P$. Because elitism preserves
  membership in $cal(S)_epsilon$ once reached, multiplying across $G$
  generations gives $"Pr"(J_G^(ast) > J^(ast) + epsilon) <= (1 - q_epsilon)^(P G)
  <= exp(-P G q_epsilon)$ via $1 - x <= e^(-x)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each generation of GKA evaluates $P$ chromosomes of length $n$. Computing
centroids takes $cal(O)(n d)$ per individual, evaluating point-to-centroid
distances for KMO and distance-biased mutation takes $cal(O)(n K d)$ per
individual, and selection/elitism takes $cal(O)(P)$. Across $G$ generations,
total time complexity is $cal(O)(G P n K d)$ (plus $cal(O)(G P K^3)$ if
Hungarian label alignment is enabled for EVO-KMeans crossover). Space
complexity is $cal(O)(P n + P K d + n d)$ because all $P$ full partition
strings must be retained in memory—making GKA substantially heavier per
iteration than single-pass streaming trees or continuous centroid swarms when
$n$ is large.

Geometrically, because GKA minimizes the exact Euclidean TWCV objective and
projects every candidate onto a Voronoi-stationary state via KMO, its feasible
clusters are strictly convex Voronoi polytopes (identical in inductive bias to
Lloyd's $K$-means @lloyd1982least and $K$-means++ @arthur2007kmeanspp), but
its stochastic evolutionary outer loop dramatically improves the likelihood of
locating the global Voronoi tessellation.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(G P n K d)$ across $G$ generations and
      population $P$],
    [Space Complexity], [$cal(O)(P n + P K d + n d)$ storing full length-$n$
      partition strings],
    [Objective Function], [Exact Euclidean SSE / TWCV with Voronoi-polytope
      inductive bias],
    [Convergence Guarantee], [Almost-sure global optimum as $G -> infinity$
      (Theorem 2); exponential worst-case $G$],
    [Key Hyperparameters], [Cluster count $K$, population $P$, generations
      $G$, mutation $p_m$, bias $c_m$],
    [Parallelizability], [Embarrassingly parallel across all $P$ chromosomes
      within each generation],
    [Invariances], [Translation, rotation, and label-permutation invariant;
      NOT scale invariant],
    [Local Polish], [Built-in one-step KMO guarantees every elite approaches a
      Voronoi local minimum],
  ),
  caption: [Operational and geometric profile of EVO-KMeans / GKA.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Medical MRI and Histopathology Tissue Segmentation*:
   Multi-channel MRI intensity distributions exhibit closely spaced modes with
   deep local SSE minima where single-run $K$-means frequently collapses gray
   and white matter regions. GKA's global evolutionary exploration paired with
   exact KMO refinement reliably identifies clinically consistent tissue
   partitions when moderate offline compute budgets are available @jain2010data.

2. *Gene Expression Microarray and Transcriptomic Module Discovery*:
   Biological cohort studies typically involve hundreds to a few thousand
   samples across high-dimensional expression profiles ($d ~ 50$–$500$). In
   this small-$n$, high-$d$ regime, $O(G P n K d)$ remains tractable while
   global SSE optimization materially improves biological pathway recovery over
   greedy restarts @xu2005survey.

3. *Vector Quantization Codebook Optimization*:
   In speech, audio, and image compression, Lloyd-Max / LBG vector quantizers
   @linde1980algorithm directly minimize MSE distortion $J(bold(s))$. Codebooks
   are trained once offline and deployed across millions of devices, making
   GKA's extra training runtime negligible relative to permanent SNR gains.

4. *Facility Location and Telecommunications Cell-Tower Placement*:
   Minimizing squared transmission distance across $K$ base stations is a
   combinatorial geometric optimization problem with abundant suboptimal basins.
   GKA's biased mutation naturally tests long-range cell swaps without losing
   Voronoi centroid optimality.

5. *Strategic Customer Segmentation and Actuarial Risk Tiering*:
   Offline financial and marketing cohort analyses operate on curated tables
   where solution stability and reproducibility outweigh millisecond latency;
   elitist GKA produces consistently low-variance partitions across repeated
   audits.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GKA encounters four well-defined failure regimes: (1) *non-convex or manifold
data* (rings, spirals, interlocking moons), where even the *global* SSE optimum
splits natural manifolds along linear Voronoi hyperplanes; (2) *large-$n$ or
streaming workloads* ($n > 10^5$), where storing $P$ length-$n$ strings and
computing $G P$ full distance passes becomes prohibitive; (3) *mutation-rate
mis-tuning*, where $p_m$ too small stalls in local minima (Theorem 3) while
$p_m$ too large destroys good schemata faster than KMO can repair them; and
(4) *outlier sensitivity*, as squared Euclidean norms allow a single extreme
point to hijack an entire cluster centroid.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-D Compact Clusters ($n <= 10^4$)], [Strong: escapes local SSE minima
      via KMO + biased mutation], [Normalize features; seed with $K$-means++],
    [Codebook / Facility Design], [Ideal: offline global MSE optimization pays
      lasting dividends], [Parallelize across $P$ workers],
    [Non-Convex Manifolds], [Poor: global SSE enforces convex Voronoi cuts],
      [Use D-Stream (Alg. 197), DBSCAN, or spectral],
    [Real-Time Data Streams], [Unsuitable: requires $G P$ multi-pass sweeps],
      [Use ClusTree (Alg. 196) or D-Stream (Alg. 197)],
    [Large-$n$ Continuous Centroids], [Heavy $cal(O)(P n)$ memory for string
      chromosomes], [Prefer PSO-Clustering (Alg. 199) or mini-batch],
  ),
  caption: [Domain suitability matrix for EVO-KMeans / GKA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

GKA sits at the intersection of discrete combinatorial search, local Voronoi
descent, and bio-inspired metaheuristics.

#proposition(
  title: "Unification with Lloyd's K-Means, Global K-Means, and PSO",
  origin: "Original Derivation",
)[
  EVO-KMeans / GKA reduces algebraically to classical clustering paradigms
  along three parameter axes:
  (i) Setting $p_m = 0$ and disabling selection/crossover reduces GKA
  identically to $P$ independent restarts of Lloyd's $K$-means @lloyd1982least
  executed for $G$ iterations each.
  (ii) Setting $P = 1$ and replacing roulette-wheel selection with Metropolis
  acceptance over single-gene flips recovers simulated-annealing $K$-means
  (and deterministic swap neighborhoods recover Kanungo et al.'s local-search
  $(9 + epsilon)$-approximation @kanungo2004local).
  (iii) Dualizing chromosome representation from discrete partition strings
  $bold(s) in {1, dots.h, K}^n$ to continuous centroid matrices
  $Z in bb(R)^(K times d)$ transforms GKA into PSO-Clustering (Alg. 199,
  @vandermerwe2003pso), cutting per-agent state memory from $cal(O)(n)$ to
  $cal(O)(K d)$ at the expense of losing KMO's exact stationarity projection.
]

#proof[
  (i) When $p_m = 0$ and selection/crossover are identities, each chromosome
  evolves independently by deterministic iteration $bold(s)_(t+1)^((p)) =
  cal(K)(bold(s)_t^((p)))$, which is verbatim Lloyd's alternating minimization
  from $P$ initial seeds, with elitism returning the best of the $P$ runs.
  (ii) With $P = 1$, a single incumbent string evolves by stochastic local
  perturbation followed by acceptance; constraining perturbations to $1$-center
  swaps yields the bounded-distortion local-search family of @kanungo2004local
  and the incremental seeding philosophy of Global $K$-Means @likas2003global.
  (iii) By Lemma 1, every KMO-projected partition is uniquely parameterized
  (up to ties) by its $K$ centroids in $bb(R)^(K times d)$; optimizing directly
  over centroid coordinates via continuous velocity updates yields PSO
  @vandermerwe2003pso (or pheromone-guided discrete assignment in ACO clustering
  @shelokar2004aco).
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark
(Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index
@hubert1985comparing), EVO-KMeans / GKA achieves Iris $"Sil" = 0.447$ and
$"ARI" = 0.583$, Customers $"Sil" = 0.340$, Multi-Shape $"Sil" = 0.712$
($"ARI"$ N/A), and High-D $"Sil" = 0.303$ ($"ARI"$ N/A), at a mean runtime of
*8.592 ms*. The numbers illuminate the exact optimization–runtime tradeoff
proved above:
(1) *Best global label recovery and high-D separation*: GKA attains the
*highest Iris ARI (0.583)* and *highest High-D Silhouette (0.303)* among all
four algorithms (Algs. 196–199), outperforming both streaming methods by nearly
$+0.09$ to $+0.20$ in Iris ARI and $+0.046$ to $+0.056$ in High-D Silhouette.
This empirically validates Theorem 2 and Lemma 1: global stochastic exploration
coupled with exact KMO projection escapes suboptimal high-dimensional basins
that trap single-pass trees and grids.
(2) *Voronoi geometry signature*: On Iris ($"Sil" = 0.447$) and Multi-Shape
($"Sil" = 0.712$), GKA trails D-Stream's Multi-Shape Silhouette (0.767) because
even global SSE minimization enforces convex linear boundaries.
(3) *Population compute cost*: Its 8.592 ms mean runtime—$3.0 times$ slower
than PSO-Clustering (2.857 ms) and $70 times$ slower than D-Stream (0.123
ms)—directly mirrors the $cal(O)(G P n K d)$ cost of mutating and projecting
length-$n$ partition strings across every generation.
