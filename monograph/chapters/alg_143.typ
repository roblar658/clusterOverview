#import "../template/ieee.typ": *

#alg_entry(
  id: 143,
  name: "ENCLUS",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1999",
  paradigm: "Entropy-Driven Apriori Subspace Discovery",
  complexity: $O(n d dot 2^(k_"max"))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

In high-dimensional data spaces, clusters frequently exist only within
low-dimensional coordinate subspaces while remaining obscured by irrelevant
noise attributes @beyer1999nearest, @kriegel2009subspacesurvey. Axis-parallel
subspace clustering pioneered by CLIQUE @agrawal1998clique detects dense
subspaces solely by thresholding single-cell point counts, which requires
tuning scale-sensitive density parameters and can report uninformative
Cartesian products of independent one-dimensional dense margins. *ENCLUS*
(ENtropy-based CLUStering), introduced by Cheng, Fu, and Zhang
@cheng1999enclus, replaces raw density counting with an information-theoretic
criterion grounded in *Shannon entropy* and *subspace interest*
(multi-information) to simultaneously enforce high coverage, high density,
and genuine inter-dimensional correlation.

Let $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ be a dataset over
attribute set $cal(A) = {1, dots, d}$. Discretize each coordinate axis
$j in cal(A)$ into $m >= 2$ equal-width intervals. Any $k$-dimensional
subspace $S = {j_1, dots, j_k} subset.eq cal(A)$ induces a regular grid
$cal(G)_S$ of $m^k$ cells. Letting $pi_S : bb(R)^d -> bb(R)^k$ denote
coordinate projection onto $S$, define the empirical cell probability
distribution
$
  p_S (x) = |{i in {1, dots, n} : pi_S (bold(x)_i) in x}| / n,
  wide x in cal(G)_S.
$

#definition(title: "Subspace Entropy, Interest, and Interest Gain")[
  For any non-empty subspace $S subset.eq cal(A)$ with empirical grid
  distribution $p_S$ over $cal(G)_S$, the *subspace Shannon entropy* (in bits)
  is
  $
    H(S) = - sum_(x in cal(G)_S, p_S (x) > 0) p_S (x) log_2 p_S (x).
  $
  The *interest* (Watanabe total correlation / multi-information) of $S$ is
  $
    "Int"(S) = sum_(j in S) H({j}) - H(S),
  $
  and the *interest gain* of extending $S$ by an attribute $r in cal(A) without S$
  is
  $
    "IG"(S, r) = "Int"(S union {r}) - "Int"(S) = H({r}) - H(r | S).
  $
]

#definition(title: "ENCLUS Clusterable Subspace Criterion")[
  Given entropy ceiling $omega > 0$ and interest floor $epsilon > 0$, a
  subspace $S subset.eq cal(A)$ is *good (clusterable)* if
  $
    H(S) < omega wide "and" wide "Int"(S) > epsilon.
  $
  A good subspace $S$ is *minimal significant* if no proper subset $S' subset S$
  satisfies both conditions, and *maximal non-redundant* if adding any
  attribute $r in.not S$ either violates $H(S union {r}) < omega$ or yields
  interest gain $"IG"(S, r) <= epsilon'$.
]

Crucially, $H(S) < omega$ is *downward-closed* (anti-monotone) while
$"Int"(S) > epsilon$ is *upward-closed* (monotone). ENCLUS exploits this dual
monotonicity in a level-wise Apriori search that maintains two lists at level
$k$: significant subspaces $cal(S)_k$ (already satisfying both thresholds) and
non-significant low-entropy subspaces $cal(N)_k$ ($H(S) < omega$ but
$"Int"(S) <= epsilon$) used to generate $(k+1)$-dimensional candidates. Once
the clusterable subspaces are discovered, a grid-connected or density
clusterer (such as CLIQUE @agrawal1998clique or DBSCAN @ester1996dbscan)
extracts clusters within each selected subspace.

#algbox(
  title: "ENCLUS: Entropy-Based Apriori Subspace Discovery",
  inputs: [Dataset $cal(X) subset bb(R)^d$, bins per axis $m$, entropy ceiling $omega$, interest floor $epsilon$, optional interest-gain threshold $epsilon'$],
  outputs: [Significant subspaces $cal(S)^*$ and subspace clusters ${cal(C)_S}_(S in cal(S)^*)$],
)[
  + Discretize each axis $j in {1, dots, d}$ into $m$ equal-width bins; compute 1D entropies $H({j})$.
  + Initialize $k <- 1$, $cal(S)_1 <- emptyset$, and $cal(N)_1 <- {{j} : H({j}) < omega}$.
  + *While* $cal(N)_k != emptyset$ (or extendable subspaces remain):
  + Generate $(k+1)$-subspace candidates $cal(C)_(k+1)$ via Apriori self-join on $cal(N)_k$ (and optionally $cal(S)_k$ if seeking maximal subspaces), pruning any $S$ containing a $k$-subset outside $cal(N)_k union cal(S)_k$.
  + *For each* $S in cal(C)_(k+1)$: hash all $n$ points into non-empty cells of $cal(G)_S$, compute $H(S)$ and $"Int"(S) = sum_(j in S) H({j}) - H(S)$.
  + *If* $H(S) < omega$ and $"Int"(S) > epsilon$: insert $S$ into $cal(S)_(k+1)$ (pruning redundant supersets with $"IG" <= epsilon'$ if mining maximal subspaces).
  + *Else if* $H(S) < omega$: insert $S$ into $cal(N)_(k+1)$. (Discard if $H(S) >= omega$.)
  + Increment $k <- k + 1$.
  + Set $cal(S)^* <- union.big_(j >= 2) cal(S)_j$; run dense-cell connected-component clustering in each $S in cal(S)^*$ and return clusters.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees governing ENCLUS: exact downward closure
of entropy, KL-divergence super-additivity of interest, quantitative
density-and-coverage bounds from low entropy, finite-sample Miller--Madow
bias and high-dimensional saturation, and an exact closed-form rotational
non-invariance proof.

#subsec("A", "Downward Closure of Entropy and Exact Apriori Soundness")

#theorem(title: "Entropy Monotonicity and Zero-False-Negative Pruning", origin: "Classical — Cited")[
  For any nested subspaces $S_1 subset.eq S_2 subset.eq cal(A)$ over a fixed
  per-axis discretization,
  $
    H(S_2) = H(S_1) + H(S_2 without S_1 | S_1) >= H(S_1),
  $
  with equality if and only if the cell index in $S_2 without S_1$ is
  uniquely determined by the cell index in $S_1$ for every point in $cal(X)$.
  Consequently, the predicate $H(S) < omega$ is strictly downward-closed, and
  ENCLUS's Apriori pruning never discards any clusterable subspace
  @cheng1999enclus.
]

#proof[
  Let $T = S_2 without S_1$. Each cell $z in cal(G)_(S_2)$ corresponds
  bijectively to an ordered pair $(x, y) in cal(G)_(S_1) times cal(G)_T$, with
  joint empirical mass $p_(S_2) (x, y) = p_(S_1) (x) p_(T | S_1) (y | x)$
  whenever $p_(S_1) (x) > 0$. Expanding $H(S_2)$ gives
  $
    H(S_2) &= - sum_(x, y) p_(S_1) (x) p_(T | S_1) (y | x)
              [log_2 p_(S_1) (x) + log_2 p_(T | S_1) (y | x)] \
           &= H(S_1) + sum_(x : p_(S_1) (x) > 0) p_(S_1) (x) H(T | S_1 = x).
  $
  Because discrete conditional probabilities satisfy
  $0 <= p_(T | S_1) (y | x) <= 1$, every term $-p_(T | S_1) (y | x) log_2 p_(T | S_1) (y | x) >= 0$,
  so $H(T | S_1) >= 0$ with equality iff for every occupied $x$, $p_(T | S_1) (dot | x)$
  is a point mass. Hence $H(S_2) < omega$ implies $H(S_1) <= H(S_2) < omega$
  for all $S_1 subset.eq S_2$, guaranteeing exact Apriori soundness.
]

#subsec("B", "KL-Divergence Representation and Super-Additivity of Interest")

#theorem(title: "Multi-Information Identity and Upward Closure of Interest", origin: "Original Unified Proof")[
  For any subspace $S subset.eq cal(A)$, let $q_S (x) = product_(j in S) p_j (x_j)$
  be the product distribution of its 1D marginals. Then:
  + $"Int"(S) = D_"KL" (p_S || q_S) >= 0$, with equality iff the discretized
    coordinates in $S$ are mutually empirically independent;
  + For any disjoint subspaces $S_1, T subset.eq cal(A)$, interest is
    *super-additive*:
    $
      "Int"(S_1 union T) = "Int"(S_1) + "Int"(T) + I(S_1; T) >= "Int"(S_1) + "Int"(T),
    $
    where $I(S_1; T) = H(S_1) + H(T) - H(S_1 union T) >= 0$ is the mutual
    information between $S_1$ and $T$;
  + Extending $S$ by an empirically independent noise attribute $r in.not S$
    produces zero interest gain $"IG"(S, r) = 0$ while increasing entropy by
    $H({r})$.
]

#proof[
  For (1), noting that $sum_(x : x_j = a) p_S (x) = p_j (a)$, we compute
  $
    D_"KL" (p_S || q_S) &= sum_(x in cal(G)_S) p_S (x) log_2 (p_S (x) / (product_(j in S) p_j (x_j))) \
                        &= - H(S) - sum_(j in S) sum_a p_j (a) log_2 p_j (a) \
                        &= sum_(j in S) H({j}) - H(S) = "Int"(S).
  $
  Non-negativity and the equality condition follow immediately from Gibbs'
  inequality. For (2), expand $"Int"(S_1 union T)$ directly:
  $
    "Int"(S_1 union T) &= sum_(j in S_1) H({j}) + sum_(j in T) H({j}) - H(S_1 union T) \
                       &= "Int"(S_1) + "Int"(T) + [H(S_1) + H(T) - H(S_1 union T)].
  $
  Since $H(S_1) + H(T) - H(S_1 union T) = D_"KL" (p_(S_1 union T) || p_(S_1) times.o p_T) = I(S_1; T) >= 0$,
  super-additivity holds. Setting $T = {r}$ yields $"IG"(S, r) = I(S; {r})$;
  when $r$ is independent of $S$, $I(S; {r}) = 0$ and
  $H(S union {r}) = H(S) + H({r})$, proving (3).
]

Theorem 2 rigorously justifies ENCLUS's dual criterion: low entropy alone
would admit independent skewed margins, whereas demanding $"Int"(S) > epsilon$
(and positive interest gain $"IG"$) isolates genuinely correlated subspaces.

#subsec("C", "Quantitative Density and Coverage Concentration Bounds")

#lemma(title: "Density-Concentration and Coverage Guarantees of Low Entropy", origin: "Original Proof")[
  Let $S subset.eq cal(A)$ be a $k$-dimensional subspace on $[0, 1]^k$ gridded
  into $m^k$ cells of volume $m^(-k)$, with empirical entropy $H(S) < omega$.
  Then:
  + *Peak-cell and average-point-density lower bounds*:
    $
      max_(x in cal(G)_S) p_S (x) >= 2^(-H(S)) > 2^(-omega),
      wide
      bar(rho)_S := n m^k sum_(x in cal(G)_S) p_S (x)^2 >= n m^k 2^(-omega);
    $
  + *Coverage-driven entropy upper bound*: If a dense cell subset
    $D subset.eq cal(G)_S$ of cardinality $r = |D|$ covers a fraction
    $gamma = sum_(x in D) p_S (x)$ of the dataset while the remaining $1 - gamma$
    points occupy $N_"bg"$ background cells, then
    $
      H(S) <= h_2 (gamma) + gamma log_2 r + (1 - gamma) log_2 N_"bg",
    $
    where $h_2 (gamma) = -gamma log_2 gamma - (1 - gamma) log_2 (1 - gamma)$
    is binary entropy.
]

#proof[
  Let $p_"max" = max_x p_S (x)$. For every occupied cell $x$,
  $log_2 (1 / p_S (x)) >= log_2 (1 / p_"max")$; averaging against $p_S (x)$
  gives $H(S) >= log_2 (1 / p_"max")$, i.e., $p_"max" >= 2^(-H(S)) > 2^(-omega)$.
  For the expected point density, note that $g(u) = 2^(-u)$ is strictly convex
  on $bb(R)$. Applying Jensen's inequality to the random variable
  $U = log_2 (1 / p_S (X))$ with $X tilde.op p_S$ yields
  $
    sum_(x : p_S (x) > 0) p_S (x)^2 = bb(E)[2^(-U)] >= 2^(-bb(E)[U]) = 2^(-H(S)) > 2^(-omega).
  $
  Multiplying by $n / m^(-k) = n m^k$ yields the lower bound on average
  experienced point density $bar(rho)_S$. For (2), let $Z = bb(I)(X in D)$ be
  the binary coverage indicator. By the chain rule,
  $H(S) = H(Z) + H(X | Z) = h_2 (gamma) + gamma H(X | X in D) + (1 - gamma) H(X | X in.not D)$.
  Since discrete entropy over $M$ outcomes is bounded above by $log_2 M$, we
  obtain $H(X | X in D) <= log_2 r$ and $H(X | X in.not D) <= log_2 N_"bg"$,
  completing the proof.
]

#subsec("D", "Finite-Sample Miller--Madow Bias and High-D Saturation")

#theorem(title: "Downward Entropy Bias and Spurious Interest in Sparse Grids", origin: "Original Proof")[
  Let $bold(x)_1, dots, bold(x)_n$ be i.i.d. samples generating cell counts
  from a true distribution $p$ on $cal(G)_S$ with $N_+ <= m^k$ positively
  supported cells. Let $hat(H)_n (S)$ be the plug-in empirical entropy. Then:
  + $bb(E)[hat(H)_n (S)] <= H(p)$ for all $n >= 1$, with asymptotic
    Miller--Madow downward bias
    $
      bb(E)[hat(H)_n (S)] = H(p) - (N_+ - 1) / (2 n ln 2) + O(n^(-2));
    $
  + *Sparse-grid saturation*: Regardless of the true distribution $p$ and
    subspace dimension $k$, $hat(H)_n (S) <= log_2 n$. Consequently, on pure
    i.i.d. uniform noise across $k$ dimensions with $m^k >> n$, the empirical
    interest satisfies
    $
      "Int"_n (S) >= k log_2 m - log_2 n - O((k m) / n),
    $
    which grows linearly in $k$ and triggers false-positive subspace
    detections whenever $k log_2 m > log_2 n + epsilon$.
]

#proof[
  Since $phi(u) = -u log_2 u$ is strictly concave on $[0, 1]$, Jensen's
  inequality applied termwise to $hat(p)_x$ (with $bb(E)[hat(p)_x] = p_x$)
  gives $bb(E)[phi(hat(p)_x)] <= phi(p_x)$, proving $bb(E)[hat(H)_n (S)] <= H(p)$.
  For the asymptotic bias, expand $psi(u) = -u ln u$ to second order around
  $p_x > 0$ using $psi'(u) = -1 - ln u$ and $psi''(u) = -1 / u$:
  $
    psi(hat(p)_x) = psi(p_x) + psi'(p_x)(hat(p)_x - p_x)
      - 1 / (2 p_x) (hat(p)_x - p_x)^2 + O(|hat(p)_x - p_x|^3).
  $
  Using $bb(E)[hat(p)_x - p_x] = 0$ and
  $bb(E)[(hat(p)_x - p_x)^2] = p_x (1 - p_x) / n$, taking expectations and
  summing over all $N_+$ supported cells (where $sum_x p_x = 1$) gives
  $
    bb(E)[sum_x psi(hat(p)_x)] = sum_x psi(p_x) - (N_+ - 1) / (2 n) + O(n^(-2)).
  $
  Dividing by $ln 2$ converts nats to bits, establishing (1). For (2), at
  most $n$ cells of $cal(G)_S$ can contain any of the $n$ points, so the
  empirical distribution has support size $<= n$ and thus
  $hat(H)_n (S) <= log_2 n$. Meanwhile, each 1D uniform margin has only $m$
  bins, so by (1) each $hat(H)_n ({j}) = log_2 m - O(m / n)$. Substituting
  into $"Int"_n (S) = sum_(j in S) hat(H)_n ({j}) - hat(H)_n (S)$ yields the
  stated linear-in-$k$ lower bound on spurious empirical interest.
]

#subsec("E", "Axis-Affine Invariance and Exact Rotational Breakdown")

#proposition(title: "Per-Axis Affine Invariance vs. Exact 45-Degree Spurious Interest", origin: "Original Counterexample")[
  + With equal-width binning per coordinate, ENCLUS is strictly invariant to
    coordinate permutations and independent positive per-axis affine rescalings
    $x_j |-> a_j x_j + b_j$ ($a_j > 0$).
  + Conversely, ENCLUS is *not* rotationally invariant: rotating completely
    unclustered uniform noise $(U, V) tilde.op "Unif"([0, 1]^2)$ by $pi / 4$
    yields continuous coordinates $(X, Y) = ((U + V)/sqrt(2), (U - V)/sqrt(2))$
    whose true differential mutual information (the fine-grid limit $m -> infinity$
    of $"Int"({1, 2})$) increases from $0$ to the exact positive constant
    $
      I(X; Y) = (1 - ln 2) / (ln 2) approx 0.4427 " bits" > 0.
    $
]

#proof[
  For (1), mapping $x_j |-> a_j x_j + b_j$ ($a_j > 0$) transforms the axis
  range $[x_"min", x_"max"]$ affinely to $[a_j x_"min" + b_j, a_j x_"max" + b_j]$,
  preserving the exact bin index of every point on every axis. All joint cell
  probabilities $p_S (x)$, $H(S)$, and $"Int"(S)$ remain identically unchanged.
  For (2), because mutual information is invariant under invertible scalar
  rescaling, $I(X; Y) = I(U + V; U - V)$. Set $W_1 = U + V$ and $W_2 = U - V$.
  The Jacobian determinant of $(U, V) |-> (W_1, W_2)$ has magnitude
  $|det mat(1, 1; 1, -1)| = 2$, so $(W_1, W_2)$ is uniform with density $1/2$
  on a diamond of area $2$, giving differential joint entropy
  $h_"nat" (W_1, W_2) = ln 2$ nats. The marginal $W_1 = U + V$ has symmetric
  triangular density $f_(W_1) (w) = w$ on $[0, 1]$ and $2 - w$ on $[1, 2]$, so
  $
    h_"nat" (W_1) = - 2 integral_0^1 w ln w dif w
      = - 2 [ (w^2 / 2) ln w - w^2 / 4 ]_0^1 = 1 / 2 " nat".
  $
  By symmetry, $W_2$ is a shifted triangular density with $h_"nat" (W_2) = 1/2$
  nat. Therefore $I(X; Y) = (1/2 + 1/2 - ln 2) / ln 2 = (1 - ln 2) / ln 2 approx 0.4427$
  bits, proving that a rigid $45 degree$ rotation injects $0.443$ bits of
  spurious interest into pure uniform noise.
]

#sec("III", "Complexity, Geometry and Operational Profile")

ENCLUS's computational cost is output-sensitive in the size of the Apriori
candidate lattice: evaluating each candidate $k$-subspace $S in cal(C)_k$
hashes $n$ points into at most $min(n, m^k)$ non-empty grid cells in $O(k n)$
time, giving total runtime $O(n sum_(k=1)^(k_"max") k |cal(C)_k|)$, which is
worst-case $O(n d dot 2^(k_"max"))$ if the entropy ceiling $omega$ is set too
loosely, yet remains sub-millisecond on moderate $d$ when $omega$ prunes
unclustered dimensions early. Space complexity is $O(n d + |cal(C)_k| k + n)$
since each candidate subspace requires only a hash table of at most $n$
occupied cells.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n sum_k k |cal(C)_k|)$ output-sensitive; $O(n d dot 2^(k_"max"))$ worst-case Apriori search],
    [Space Complexity], [$O(n d + |cal(C)_k| k)$; at most $n$ occupied hash buckets per candidate subspace],
    [Cluster Geometry], [Axis-parallel overlapping subspace clusters (arbitrary connected unions of grid cells per $S$)],
    [Noise Handling], [High-entropy noise axes pruned globally; sparse cells filtered during cluster extraction],
    [Hyperparameters], [Bin count $m$, entropy threshold $omega$, interest threshold $epsilon$, interest-gain $epsilon'$],
    [Invariances], [Per-axis positive affine & permutation invariant; NOT rotationally invariant (Prop. 1)],
  ),
  caption: [Operational, geometric, and complexity profile of ENCLUS.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Gene Expression and Multi-Omics Pathway Discovery*: Microarray and RNA-seq
  cohorts profile thousands of genes across hundreds of biological conditions,
  where functional modules co-vary tightly across small subsets of samples or
  genes. ENCLUS's interest criterion $"Int"(S) > epsilon$ directly captures
  nonlinear co-regulation (mutual information) missed by Pearson correlation,
  while $H(S) < omega$ suppresses constitutively flat or noisy genes
  @cheng1999enclus.
+ *Financial Portfolio Regime and Cross-Asset Co-Movement Mining*: Equity and
  derivative markets exhibit localized clusters of strongly coupled assets
  during volatility regimes. Discretized return distributions are heavy-tailed
  and non-Gaussian; Shannon entropy and total correlation detect multi-asset
  dependence structures without assuming Gaussian covariance.
+ *Targeted Marketing and Customer Attribute Selection*: High-dimensional CRM
  tables mix demographics, channel usage, and purchase histories where only
  specific 2--5-attribute combinations exhibit crisp behavioral segments.
  ENCLUS identifies exact interpretable attribute subsets whose joint
  distribution departs sharply from independence.
+ *Semiconductor Yield Root-Cause Analysis*: Manufacturing test suites record
  hundreds of parametric wafer measurements per die, with defect signatures
  concentrated in low-dimensional combinations of interacting tool parameters.
  ENCLUS's minimal significant subspaces pinpoint the smallest set of
  interacting variables driving yield excursions.
+ *Intrusion Detection and Multi-Sensor Industrial Fault Diagnosis*: Cyber and
  SCADA logs contain dozens of heterogeneous metrics whose baseline operation
  is loosely coupled, whereas coordinated faults induce low-entropy,
  high-interest signatures in specific sensor subspaces.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

ENCLUS exhibits four principled failure regimes: (i) *finite-sample grid
sparsity* (Theorem 3), where $m^k >> n$ induces severe Miller--Madow downward
bias and spurious interest on independent noise; (ii) *oblique rotated
subspaces* (Proposition 1), where axis discretization blurs diagonal clusters
or manufactures spurious $0.443$-bit diamond correlations; (iii) *exponential
candidate explosion* if $omega$ is overly permissive in high $d$; and
(iv) *bin-boundary discretization artifacts* when cluster cores straddle fixed
equal-width grid edges.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Moderate-D tabular & bioinformatics], [Excellent: detects nonlinear multi-attribute coupling via $"Int"(S)$], [Enforce $m^k << n$ or apply Miller--Madow correction],
    [Feature-subset discovery (2--6D)], [Strong: principled dual $(omega, epsilon)$ pruning of noise axes], [Use quantile bins if margins are heavy-tailed],
    [Very high-D / small-$n$ ($m^k >> n$)], [Poor: $hat(H)_n <= log_2 n$ causes spurious interest (Thm. 3)], [Use FINDIT (Alg. 144) or PROCLUS @aggarwal1999proclus],
    [Arbitrarily oriented / oblique manifolds], [Weak: axis-binning loses $O(d)$ invariance (Prop. 1)], [Use ORCLUS @aggarwal2000orclus or OptiGrid (Alg. 141)],
    [Large-scale continuous $k_"max" > 10$], [Weak: Apriori lattice grows exponentially in $k_"max"$], [Use FIRES @kriegel2005fires or SUBCLU @kailing2004subclu],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for ENCLUS.],
)

#sec("V", "Structural Relations and Empirical Behavior")

ENCLUS directly generalizes CLIQUE @agrawal1998clique and sits alongside
density-based subspace miners MAFIA @goil1999mafia, SUBCLU @kailing2004subclu,
FIRES @kriegel2005fires, and INSCY @assent2008inscy. Whereas CLIQUE enforces
only a per-cell density threshold (risking spurious high-dimensional joins of
independent skewed attributes), ENCLUS introduces an information-theoretic
filter that strictly tightens CLIQUE's candidate lattice.

#proposition(title: "Information-Theoretic Refinement of CLIQUE", origin: "Original Proof")[
  Let $cal(S)_"CLIQUE" (tau, gamma)$ be the set of subspaces containing a
  connected union $D$ of $r$ dense cells ($p_S (x) >= tau$) covering at least
  fraction $gamma$ of points, with $r <= floor(1 / tau)$ and $N_"bg"$ occupied
  background cells. If the ENCLUS ceiling satisfies
  $omega > h_2 (gamma) + gamma log_2 (1 / tau) + (1 - gamma) log_2 N_"bg"$,
  then every CLIQUE subspace passes ENCLUS's entropy filter $H(S) < omega$,
  while ENCLUS strictly refines $cal(S)_"CLIQUE"$ by eliminating all
  independent-margin subspaces with $D_"KL" (p_S || product_(j in S) p_j) <= epsilon$.
]

#proof[
  Any subspace $S in cal(S)_"CLIQUE" (tau, gamma)$ has at most
  $r <= floor(gamma / tau) <= floor(1 / tau)$ dense cells covering mass
  $gamma$. Applying Lemma 1(2) directly gives
  $H(S) <= h_2 (gamma) + gamma log_2 (1 / tau) + (1 - gamma) log_2 N_"bg" < omega$.
  Conversely, if $S = S_1 union S_2$ is formed from two independent skewed
  1D attributes each containing a dense interval, $S$ passes CLIQUE's density
  test yet has $"Int"(S) = D_"KL" (p_S || q_S) = 0 <= epsilon$ by Theorem 2,
  so ENCLUS rejects it.
]

Compared to projected-partitioning methods such as PROCLUS @aggarwal1999proclus,
DOC @procopiuc2002doc, and FINDIT (Alg. 144, @woo2004findit), ENCLUS allows
overlapping clusters across multiple distinct subspaces rather than forcing
each data point into a single global partition.

On the ClusterOverview 4-dataset benchmark (Silhouette @rousseeuw1987silhouettes
and ARI @hubert1985comparing), ENCLUS obtains: *Iris* $"Sil" = 0.406$,
$"ARI" = 0.818$; *Customers* $"Sil" = 0.263$; *Multi-Shape* $"Sil" = 0.609$
($"ARI" =$ N/A); *High-D* $"Sil" = 0.232$ ($"ARI" =$ N/A); with a mean runtime
of $0.048$ ms. The standout result is *Iris* $"ARI" = 0.818$---far exceeding
OptiGrid ($0.567$), BANG ($0.510$), and FINDIT ($0.574$)---even while its
full-space Silhouette is a modest $0.406$. This divergence is a textbook
illustration of subspace inductive bias: ENCLUS automatically identifies the
low-entropy, high-interest petal-length / petal-width subspace that cleanly
separates the three *Iris* species (driving ARI to $0.818$), whereas full-4D
Euclidean Silhouette penalizes the resulting partition along noisy sepal-width
directions. Its $0.048$ ms runtime confirms rapid Apriori pruning on small $d$,
while its moderate *High-D* Silhouette ($0.232$) mirrors Theorem 3's
finite-sample grid-saturation bias when $m^k$ approaches $n$.
