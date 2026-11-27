#import "../template/ieee.typ": *

#alg_entry(
  id: 170,
  name: "Consensus Clustering (Monti)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2003",
  paradigm: "Resampling-Based Stability Consensus & Hierarchical Partitioning",
  complexity: $O(H (T_"base" (n, d, K) + n^2) + n^2 log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = {x_1, dots, x_n} subset.eq bb(R)^d$ denote a finite dataset indexed by $V = {1, dots, n}$, and let $cal(A)_K$ be a base clustering operator mapping any nonempty subsample $D subset.eq V$ to a $K$-way partition $pi : D -> {1, dots, K}$. Partitioning methods such as $K$-means @lloyd1982least, PAM @kaufman1990finding, and spectral clustering @ng2002spectral exhibit sensitivity to random initialization, sampling fluctuations, and misspecification of $K$. Monti et al. @monti2003consensus formulate cluster discovery and model-order selection jointly as a resampling stability problem: true latent structure should induce pairwise co-membership relations that persist across repeated perturbations of the dataset, whereas spurious partitions fracture unpredictably across subsamples.

Let $D^((1)), dots, D^((H)) subset.eq V$ be $H$ independent resampled subsets, drawn either by uniform subsampling without replacement of size $m = floor(p n)$ with sampling proportion $p in (0, 1]$ (canonically $p = 0.8$) or by bootstrap resampling with replacement. For each replicate $h in {1, dots, H}$, define two symmetric binary matrices in ${0, 1}^(n times n)$: the _co-inclusion indicator_ $I^((h))$ and the _connectivity (co-clustering) matrix_ $M^((h))$, with entries
$ I^((h))_(i j) = bb(1)(i in D^((h)) and j in D^((h))), $
$ M^((h))_(i j) = I^((h))_(i j) dot bb(1)(pi^((h)) (i) = pi^((h)) (j)), $
where $pi^((h)) = cal(A)_K (D^((h)))$. Accumulating across $H$ replicates yields the pair-observation count matrix $N = sum_(h=1)^H I^((h))$ and raw co-occurrence matrix $S = sum_(h=1)^H M^((h))$.

#definition(title: "Monti Consensus Matrix, Empirical CDF, and Area Functional")[
  For a fixed cluster count $K >= 2$, the _consensus matrix_ $cal(M)(K) in [0, 1]^(n times n)$ is the entrywise ratio
  $ cal(M)_(i j) (K) = cases(S_(i j) / N_(i j) & "if " N_(i j) > 0, 0 & "if " N_(i j) = 0.) $
  Let $cal(U) = {(i, j) : 1 <= i < j <= n}$ index the $L = binom(n, 2)$ upper-triangular off-diagonal pairs, and let $u_1 <= u_2 <= dots <= u_L$ denote the sorted values of ${cal(M)_(i j) (K) : (i, j) in cal(U)}$. The _empirical consensus CDF_ is
  $ "CDF"_K (c) = 1 / L sum_(r=1)^L bb(1)(u_r <= c), quad c in [0, 1], $
  with associated area functional $A(K) = integral_0^1 "CDF"_K (c) dif c$ and relative area increment
  $ Delta(K) = cases(A(2) & "if " K = 2, (A(K) - A(K - 1)) / A(K - 1) & "if " K > 2.) $
]

An idealized, completely stable $K$-clustering produces a binary consensus matrix $cal(M)_(i j) in {0, 1}$ (a block-diagonal permutation of all-ones blocks), whose empirical CDF is a two-step function jumping at $c = 0$ and $c = 1$. Conversely, unstable partitions spread off-diagonal mass across the interior $(0, 1)$, flattening the CDF. To summarize cluster-level and object-level stability, Monti et al. @monti2003consensus define the _cluster consensus_ $m(c)$ for a final consensus cluster $C_c$ and the _item consensus_ $m_i (c)$ of item $i$ with cluster $C_c$:
$ m(c) = 1 / (abs(C_c)(abs(C_c) - 1)) sum_(i, j in C_c, i != j) cal(M)_(i j), $
$ m_i (c) = 1 / (abs(C_c) - bb(1)(i in C_c)) sum_(j in C_c, j != i) cal(M)_(i j). $
Final hard clusters are extracted by treating $D = bold(1) bold(1)^top - cal(M)(K)$ as a pairwise dissimilarity matrix and executing agglomerative hierarchical clustering (typically average linkage, UPGMA @sokal1958statistical) cut at $K$ branches.

#algbox(
  title: "Consensus Clustering (Monti et al., 2003)",
  inputs: [Dataset $X = {x_1, dots, x_n}$, candidate range ${2, dots, K_"max"}$, base algorithm $cal(A)_K$, subsample fraction $p in (0, 1]$, replicates $H$.],
  outputs: [Consensus matrices ${cal(M)(K)}$, area curve ${Delta(K)}$, selected $hat(K)$, and final partition $hat(cal(C))$.],
)[
  + *For each $K = 2, dots, K_"max"$*: initialize $S <- bold(0)_(n times n)$ and $N <- bold(0)_(n times n)$.
  + *Resampling loop ($h = 1, dots, H$)*: draw $D^((h)) subset.eq V$ of size $m = floor(p n)$ uniformly without replacement (or bootstrap); execute $pi^((h)) <- cal(A)_K (D^((h)))$.
  + *Matrix accumulation*: for all pairs $i, j in D^((h))$, increment $N_(i j) <- N_(i j) + 1$; if $pi^((h)) (i) = pi^((h)) (j)$, increment $S_(i j) <- S_(i j) + 1$.
  + *Consensus normalization*: set $cal(M)_(i j) (K) <- S_(i j) / N_(i j)$ for $N_(i j) > 0$ (and $0$ otherwise), with diagonal $cal(M)_(i i) (K) = 1$ whenever $N_(i i) > 0$.
  + *Stability diagnostics*: compute $"CDF"_K (c)$, area $A(K)$, and relative change $Delta(K)$ (optionally augmenting with Proportion of Ambiguous Clustering $"PAC"_K (u_1, u_2) = "CDF"_K (u_2) - "CDF"_K (u_1)$).
  + *Consensus extraction*: for chosen $hat(K)$ (elbow of $Delta(K)$ or minimizer of $"PAC"_K$), form dissimilarity $D = bold(1) bold(1)^top - cal(M)(hat(K))$ and cut the average-linkage dendrogram at $hat(K)$ clusters to output $hat(cal(C))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five structural theorems governing the algebraic geometry of $cal(M)$, finite-sample concentration of resampling estimates, exact integral characterization of Monti's area functional, exact recovery of planted partitions under UPGMA, and asymptotic stability limits.

#subsec("A", "Positive Semidefiniteness and Metric Breakdown Under Subsampling")

A frequent assumption in downstream kernel or spectral post-processing of $cal(M)$ is that $cal(M)$ is positive semidefinite (PSD) and $1 - cal(M)$ is a metric. The following result proves that both properties hold unconditionally at $p = 1$, yet both can fail simultaneously whenever $p < 1$.

#theorem(title: "PSD and Pseudometric Dichotomy of Consensus Matrices", origin: "Original Proof")[
  Let $cal(M) in [0, 1]^(n times n)$ be the Monti consensus matrix formed from $H$ replicates, and define $d_(i j) = 1 - cal(M)_(i j)$.
  (a) Under full-sample resampling ($p = 1$, so $N_(i j) = H$ for all $i, j$), $cal(M) succ.eq 0$ is positive semidefinite and $d$ is a pseudometric on $V$.
  (b) Under proper subsampling ($p < 1$), non-uniform pair counts $N_(i j)$ can violate both positive semidefiniteness of $cal(M)$ and the triangle inequality for $d$, even in the infinite-replicate limit $H -> infinity$.
]
#proof[
  *Part (a):* For each replicate $h$, let $z^((h))_c in {0, 1}^n$ be the indicator vector of cluster $c in {1, dots, K}$ in $pi^((h))$ (zero outside $D^((h))$). Since clusters partition $D^((h))$, we have the rank-$K$ factorization
  $ M^((h)) = sum_(c=1)^K z^((h))_c (z^((h))_c)^top succ.eq 0. $
  When $p = 1$, $N_(i j) = H$ for every pair, so $cal(M) = (1 / H) sum_(h=1)^H M^((h))$ is a convex combination of PSD matrices and hence satisfies $v^top cal(M) v = (1 / H) sum_(h=1)^H sum_(c=1)^K (v^top z^((h))_c)^2 >= 0$ for all $v in bb(R)^n$. For the dissimilarity, define $d^((h))_(i j) = 1 - M^((h))_(i j) in {0, 1}$. Because co-clustering in $pi^((h))$ is an equivalence relation on $V$, if $d^((h))_(i j) = 0$ and $d^((h))_(j l) = 0$, then $pi^((h)) (i) = pi^((h)) (j) = pi^((h)) (l)$, forcing $d^((h))_(i l) = 0 <= 0 + 0$; if either $d^((h))_(i j) = 1$ or $d^((h))_(j l) = 1$, then $d^((h))_(i l) <= 1 <= d^((h))_(i j) + d^((h))_(j l)$. Thus every $d^((h))$ is a pseudometric, and their average $d = (1 / H) sum_(h=1)^H d^((h))$ inherits symmetry, $d_(i i) = 0$, and the triangle inequality.

  *Part (b):* Construct an explicit $n = 3$, $m = 2$ ($p = 2/3$) instance with $H = 3$ replicates $D^((1)) = {1, 2}$, $D^((2)) = {2, 3}$, $D^((3)) = {1, 3}$. Let $cal(A)$ place ${1, 2}$ together on $D^((1))$, place ${2, 3}$ together on $D^((2))$, and separate ${1}, {3}$ on $D^((3))$ (a deterministic order- or geometry-sensitive base clusterer). Then $N_(i i) = 2$, $S_(i i) = 2$, while $N_(1 2) = N_(2 3) = N_(1 3) = 1$, yielding
  $ cal(M) = mat(1, 1, 0; 1, 1, 1; 0, 1, 1). $
  Evaluating the quadratic form at $v = (1, -1, 1)^top$ gives $cal(M) v = (0, 1, 0)^top$ and $v^top cal(M) v = -1 < 0$; indeed, the characteristic polynomial $(1 - lambda)((1 - lambda)^2 - 2) = 0$ yields the strictly negative eigenvalue $lambda_min (cal(M)) = 1 - sqrt(2) approx -0.4142 < 0$. Simultaneously, $d_(1 3) = 1 - 0 = 1 > 0 + 0 = d_(1 2) + d_(2 3)$, violating the triangle inequality. Repeating these three subsamples uniformly with probability $1/3$ each as $H -> infinity$ converges almost surely by the Strong Law of Large Numbers to the exact same population consensus matrix $cal(M)^infinity = cal(M)$, proving that neither pathology vanishes with large $H$.
]

Theorem 1 provides a fundamental algorithmic justification for Monti et al.'s use of agglomerative hierarchical clustering rather than unprojected spectral or Euclidean MDS methods on $cal(M)$: when $p < 1$, the Hadamard division $S_(i j) / N_(i j)$ destroys the cone structure of $S succ.eq 0$.

#subsec("B", "Finite-Sample Concentration of the Resampled Consensus Matrix")

To determine how many replicates $H$ suffice to estimate the population co-clustering probabilities reliably uniformly over all $binom(n, 2)$ pairs, we derive a finite-sample concentration bound accounting explicitly for the random denominators $N_(i j)$.

#theorem(title: "Uniform Entrywise Concentration of the Consensus Matrix", origin: "Original Derivation")[
  Suppose $D^((1)), dots, D^((H))$ are drawn i.i.d. uniformly without replacement of size $m = floor(p n) >= 2$, and conditional on $D^((h))$ the base clusterings $pi^((h))$ are conditionally independent with identical pair-wise co-clustering probability
  $ p_(i j) = "Pr"(pi^((h)) (i) = pi^((h)) (j) | i, j in D^((h))). $
  Let $q_2 = m(m - 1) / (n(n - 1))$ be the pair inclusion probability. For any accuracy $epsilon in (0, 1)$ and confidence $delta in (0, 1)$, if
  $ H >= 1 / q_2 max{ 1 / epsilon^2 ln((2 n^2) / delta), quad 8 ln(n^2 / delta) }, $
  then with probability at least $1 - delta$,
  $ max_(1 <= i < j <= n) abs(cal(M)_(i j) - p_(i j)) <= epsilon. $
]
#proof[
  Fix a pair $(i, j)$ with $i < j$. Each indicator $I^((h))_(i j)$ is i.i.d. Bernoulli($q_2$), so $N_(i j) tilde.op "Binomial"(H, q_2)$ with expectation $mu = H q_2$. Define the low-count event $E^c = {N_(i j) < mu / 2}$. By the multiplicative Chernoff bound for binomial tails with deviation parameter $gamma = 1/2$,
  $ "Pr"(N_(i j) < mu / 2) <= exp(-gamma^2 mu / 2) = exp(-H q_2 / 8). $
  Now condition on the random index set $cal(S)_(i j) = {h : I^((h))_(i j) = 1}$ with cardinality $r = abs(cal(S)_(i j)) >= ceil(mu / 2) >= 1$. Because replicates are drawn independently, conditioning on the event $i, j in D^((h))$ leaves the remaining $m - 2$ elements of $D^((h))$ independent across $h in cal(S)_(i j)$, so the variables ${M^((h))_(i j)}_(h in cal(S)_(i j))$ are conditionally i.i.d. Bernoulli($p_(i j)$) bounded in $[0, 1]$. Applying Hoeffding's inequality to their sample mean $cal(M)_(i j) = (1 / r) sum_(h in cal(S)_(i j)) M^((h))_(i j)$ yields
  $ "Pr"(abs(cal(M)_(i j) - p_(i j)) >= epsilon | N_(i j) = r) <= 2 exp(-2 r epsilon^2). $
  Since $r >= H q_2 / 2$ on $E$, we obtain $2 exp(-2 r epsilon^2) <= 2 exp(-H q_2 epsilon^2)$. Integrating over $N_(i j)$ via total probability gives
  $ "Pr"(abs(cal(M)_(i j) - p_(i j)) >= epsilon) <= exp(-H q_2 / 8) + 2 exp(-H q_2 epsilon^2). $
  Applying a union bound over all $binom(n, 2) < n^2 / 2$ unordered pairs,
  $ "Pr"(max_(i < j) abs(cal(M)_(i j) - p_(i j)) >= epsilon) & <= n^2 / 2 exp(-H q_2 / 8) \
    & quad + n^2 exp(-H q_2 epsilon^2). $
  Under the stated lower bound on $H$, we have $H q_2 / 8 >= ln(n^2 / delta)$ and $H q_2 epsilon^2 >= ln(2 n^2 / delta)$, making each term at most $delta / 2$ and their sum at most $delta$.
]

Because $q_2 approx p^2$, Theorem 2 reveals the precise sample-complexity penalty of aggressive subsampling: halving $p$ from $0.8$ to $0.4$ quadruples the required replicate budget $H = Theta(p^(-2) epsilon^(-2) log(n / delta))$ needed to maintain identical uniform precision across $cal(M)$.

#subsec("C", "Exact Integral Identity for CDF Area and Ambiguity Blindness")

Monti et al. @monti2003consensus select $K$ via the area $A(K)$ under $"CDF"_K$ and its relative change $Delta(K)$. We prove an exact closed-form reduction showing that $A(K)$ depends solely on the first moment of $cal(M)$ and is therefore blind to dispersion in the ambiguous interior $(0, 1)$.

#theorem(title: "CDF Area Identity and Mean-Preserving Ambiguity Blindness", origin: "Original Proof")[
  Let $u_1 <= dots <= u_L$ be the $L = binom(n, 2)$ off-diagonal entries of $cal(M)(K)$ in $[0, 1]$ with mean $bar(u)_K = (1 / L) sum_(r=1)^L u_r$ and variance $"Var"(u) = (1 / L) sum_(r=1)^L (u_r - bar(u)_K)^2$. Then:
  (a) The exact area under the empirical CDF satisfies
  $ A(K) = integral_0^1 "CDF"_K (c) dif c = 1 - bar(u)_K. $
  (b) The quadratic consensus ambiguity $Q(cal(M)) = (1 / L) sum_(r=1)^L u_r (1 - u_r)$ obeys the exact decomposition
  $ Q(cal(M)) = bar(u)_K (1 - bar(u)_K) - "Var"(u) >= 0, $
  with $Q(cal(M)) = 0$ if and only if $cal(M)$ is perfectly crisp ($u_r in {0, 1}$ for all $r$). Consequently, a perfectly crisp balanced $K$-clustering and a maximally ambiguous uniform matrix $u_r equiv 1 / K$ yield identical CDF areas $A(K) = 1 - 1 / K$, despite achieving the minimum ($0$) and maximum ($(K - 1) / K^2$) ambiguity respectively.
]
#proof[
  *Part (a):* Set $u_0 = 0$ and $u_(L+1) = 1$. On each $[u_r, u_(r+1))$ ($0 <= r <= L$), $"CDF"_K (c) = r / L$. Piecewise integration and summation by parts yield
  $ A(K) & = sum_(r=0)^L r / L (u_(r+1) - u_r) \
    & = 1 / L ( L u_(L+1) - sum_(r=1)^L u_r ) = 1 - bar(u)_K, $
  since $u_(L+1) = 1$.

  *Part (b):* Expanding $Q(cal(M))$ gives
  $ Q(cal(M)) = 1 / L sum_(r=1)^L u_r - 1 / L sum_(r=1)^L u_r^2 = bar(u)_K (1 - bar(u)_K) - "Var"(u). $
  Because $u_r in [0, 1]$, every term $u_r (1 - u_r) >= 0$, vanishing iff $u_r in {0, 1}$. A crisp balanced $K$-clustering ($bar(u)_K = 1 / K$, $"Var"(u) = (1 / K)(1 - 1 / K)$) and the uniform matrix $u_r equiv 1 / K$ ($"Var"(u) = 0$) both produce $A(K) = 1 - 1 / K$, while yielding $Q = 0$ versus $Q = (K - 1) / K^2$.
]

Theorem 3 explains why the Proportion of Ambiguous Clustering $"PAC"_K (u_1, u_2) = "CDF"_K (u_2) - "CDF"_K (u_1)$ ($0 < u_1 < u_2 < 1$) strictly sharpens $Delta(K)$: when $1/K in (u_1, u_2)$, $A(K)$ is blind to dispersion whereas $"PAC"_K$ separates crisp ($0$) from diffuse ($1$) consensus maximally.

#subsec("D", "Exact Recovery Under Average-Linkage Consensus Extraction")

#theorem(title: "Exact Zero-Error Recovery by UPGMA on Consensus Dissimilarity", origin: "Original Proof")[
  Let $cal(C)^* = {C_1^*, dots, C_K^*}$ be a ground-truth $K$-partition with population co-clustering margin
  $ p_"in" = min_c min_(i, j in C_c^*) p_(i j) > max_(c != c') max_(i in C_c^*, j in C_(c')^*) p_(i j) = p_"out", $
  and gap $Delta_p = p_"in" - p_"out" > 0$. If $H > q_2^(-1) max{4 Delta_p^(-2) ln(2 n^2 / delta), 8 ln(n^2 / delta)}$, then with probability at least $1 - delta$, cutting the average-linkage (UPGMA) dendrogram on $D = bold(1) bold(1)^top - cal(M)$ at $K$ clusters recovers $cal(C)^*$ with zero error ($"ARI" = 1$).
]
#proof[
  Applying Theorem 2 with $epsilon < Delta_p / 2$, the event $cal(E) = {max_(i < j) abs(cal(M)_(i j) - p_(i j)) < Delta_p / 2}$ holds with probability at least $1 - delta$. Let $tau = 1 - (p_"in" + p_"out") / 2$. On $cal(E)$, every intra-cluster pair satisfies $D_(i j) = 1 - cal(M)_(i j) < tau$, while every inter-cluster pair satisfies $D_(i j) > tau$.

  We prove by induction on UPGMA merge steps $t = 0, dots, n - K$ that every active cluster is a _pure_ subset of a unique $C_c^*$. The base case ($t = 0$, singletons) is trivial. Suppose after $t < n - K$ merges there remain $n - t > K$ pure active clusters $A_1, dots, A_(n-t)$. By the pigeonhole principle, some distinct $A_a, A_b subset.eq C_c^*$ share a true cluster, so their average dissimilarity is
  $ bar(D)(A_a, A_b) = 1 / (abs(A_a) abs(A_b)) sum_(i in A_a, j in A_b) D_(i j) < tau. $
  Conversely, any $A_u subset.eq C_(c_1)^*, A_v subset.eq C_(c_2)^*$ ($c_1 != c_2$) satisfy $bar(D)(A_u, A_v) > tau$. Hence the minimum UPGMA distance at step $t + 1$ strictly prefers a pure intra-cluster merge, preserving purity. At $t = n - K$, the $K$ remaining nonempty pure blocks must equal $cal(C)^*$ identically.
]

#subsec("E", "Asymptotic Stability Limitations Under Unimodal Null Models")

#proposition(title: "Risk-Minimization Stability Without Cluster Structure", origin: "Classical — Cited")[
  As proven by Ben-David, von Luxburg, and Pál @ben-david2006sober, asymptotic resampling stability of risk-minimizing base clusterers $cal(A)_K$ (such as $K$-means) reflects uniqueness of the population risk minimizer rather than true multimodality. For a unimodal anisotropic Gaussian $cal(N)(bold(0), "diag"(sigma_1^2, sigma_2^2))$ with $sigma_1 > sigma_2$, $cal(M)(2)$ converges to a crisp binary matrix ($Q -> 0$, $"PAC"_2 -> 0$) despite $K^* = 1$.
]
#proof[
  Any $2$-means hyperplane through the origin with normal $u = (cos theta, sin theta)^top$ yields population distortion $sigma_1^2 + sigma_2^2 - (2 / pi)(sigma_1^2 cos^2 theta + sigma_2^2 sin^2 theta)$, uniquely minimized at $theta in {0, pi}$ when $sigma_1 > sigma_2$. By consistency of empirical $K$-means minimizers @ben-david2006sober, subsample split hyperplanes converge a.s. to ${x_1 = 0}$, forcing $p_(i j) -> bb(1)(x_(i,1) x_(j,1) > 0)$ and vanishing consensus ambiguity on strictly unimodal data.
]

#sec("III", "Complexity, Geometry and Operational Profile")

For fixed $K$, $H$ base runs on size $m = floor(p n)$ cost $O(H T_"base" (p n, d, K))$; updating dense $N$ and $S$ costs $O(H p^2 n^2)$; UPGMA via nearest-neighbor chains @mullner2011modern costs $O(n^2)$ ($O(n^2 log n)$ with binary heaps); and sorting off-diagonal entries for $"CDF"_K$ takes $O(n^2 log n)$ ($O(n^2)$ via histogram binning). Sweeping $K in {2, dots, K_"max"}$ totals $O(K_"max" (H (T_"base" (p n, d) + p^2 n^2) + n^2 log n))$. Memory is dominated by $Theta(n^2)$ storage for $N, S, cal(M)$ ($approx 16 n^2$ bytes in float64, limiting exact execution to $n <= 3 times 10^4$), while replicates parallelize linearly across cores.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(H (T_"base" (p n, d, K) + p^2 n^2) + n^2 log n)$ per $K$],
    [Space Complexity], [$O(n^2 + n d)$ for dense $S$, $N$, and $cal(M)$],
    [Geometric Bias], [Inherited from $cal(A)_K$, smoothed via average linkage on $1 - cal(M)$],
    [Hyperparameters], [Base clusterer $cal(A)_K$, $K$-range, $p in [0.7, 0.9]$, $H in [50, 1000]$],
    [Determinism], [Stochastic; uniform error $O(p^(-1) sqrt(log(n)/H))$ (Thm. 2)],
    [Invariances], [Permutation invariant to base labels; inherits $cal(A)_K$ feature invariance],
    [Outputs], [Soft stability matrix $cal(M)$, $"CDF"_K$, $Delta(K)$, item/cluster consensus],
  ),
  caption: [Operational and geometric profile of Monti Consensus Clustering (Alg. 170).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Cancer Molecular Subtyping from Transcriptomics and Multi-Omics:* Monti et al. @monti2003consensus originally introduced consensus clustering for DNA microarray gene-expression cohorts ($n approx 50$–$1000$, $d approx 10^3$–$2 times 10^4$). In this $n << d$ regime, single-run $K$-means or hierarchical partitions are severely sensitive to outlier patients and probe noise. Subsampling $p = 0.8$ directly simulates cohort perturbation, while item consensus $m_i (c)$ flags borderline clinical specimens whose molecular subtype assignment lacks reproducibility.

*2. Single-Cell RNA-seq Cell-Type Discovery and Quality Auditing:* On pre-clustered or metacell-aggregated scRNA-seq cohorts, combining stochastic graph or $K$-means base runs via $cal(M)$ suppresses dropout-induced cluster splitting. Theorem 4 guarantees that well-separated cell lineages with $p_"in" > p_"out"$ are recovered cleanly by UPGMA even when individual replicates misclassify up to $(p_"in" - p_"out") / 2$ of boundary cells.

*3. Neuroimaging Functional Connectivity Parcellation:* Resting-state fMRI voxel or ROI parcellations suffer from high temporal autocorrelation and subject-level variability. Subject- or epoch-level resampling builds a cortical co-assignment consensus matrix $cal(M)$ whose cluster consensus $m(c)$ quantifies anatomical stability across scanning sessions and identifies reproducible resting-state networks.

*4. Customer and Financial Cohort Stability Validation:* In risk modeling and marketing segmentation, downstream business policies depend on cluster identities remaining stable when new customer waves arrive. Consensus clustering operationalizes temporal and sample robustness directly: low item consensus $m_i (c)$ identifies customers straddling adjacent behavioral segments who should receive hedged or multi-segment treatments.

*5. Ensemble Hyperparameter and Model-Order Screening:* Across chemoinformatics and materials informatics, consensus CDF curves and $"PAC"_K$ profiles serve as model-agnostic diagnostics to rule out unstable cluster counts $K$ before committing to costly experimental validation, pairing naturally with internal validity indices such as Silhouette @rousseeuw1987silhouettes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four mathematically distinct failure modes constrain Monti Consensus Clustering:
(i) _Anisotropic unimodal false positives_ (Prop. 5): elongated unimodal clouds produce artificially crisp $K >= 2$ consensus matrices because $K$-means splits deterministically along the principal axis @ben-david2006sober; testing $K = 1$ requires parametric or reference-null simulation (e.g., Monte Carlo Reference Consensus).
(ii) _CDF area blindness_ (Thm. 3): relying solely on $Delta(K)$ can mask diffuse interior consensus mass unless paired with $"PAC"_K$ or $Q(cal(M))$.
(iii) _Quadratic scalability barrier_: $O(n^2)$ memory prevents direct deployment on raw million-point datasets without prior coreset or mini-batch reduction.
(iv) _Base-clusterer bias inheritance_: if $cal(A)_K$ is consistently misspecified (e.g., spherical $K$-means on concentric rings), $cal(M)$ amplifies that biased partition with high apparent stability.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Genomics / Omics ($n <= 5000$)], [High: directly quantifies cohort perturbation stability and borderline samples], [Verify $K=1$ null against anisotropic unimodality],
    [Neuroimaging Parcellation], [High: aggregates noisy subject/session partitions smoothly], [Use sparse co-assignment for $> 10^5$ voxels],
    [Behavioral Segmentation], [Medium-High: item consensus identifies transitional customers], [Pair $Delta(K)$ with $"PAC"_K$ and Silhouette],
    [Non-Convex Manifolds], [Conditional: succeeds iff base $cal(A)_K$ resolves geometry], [Use spectral / density base or CSPA (Alg. 171)],
    [Web-Scale Streams ($n > 10^5$)], [Poor: $O(n^2)$ matrix $cal(M)$ exceeds memory], [Prefer MCLA (Alg. 173) or streaming coreset],
  ),
  caption: [Domain suitability and failure-mode matrix for Monti Consensus Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Monti Consensus Clustering bridges resampling-based model selection and cluster ensembles @jain2010data. When subsampling is disabled ($p = 1$) and replicates come from a fixed heterogeneous ensemble $Lambda = {pi^((1)), dots, pi^((r))}$, the consensus matrix $cal(M)$ coincides identically with the co-association similarity matrix $S = (1 / r) H H^top$ of CSPA (Alg. 171 @strehl2002cluster), differing solely in the final partitioning operator: Monti applies agglomerative average linkage (UPGMA @sokal1958statistical) to $1 - cal(M)$, whereas CSPA applies balanced $k$-way graph partitioning (METIS) to $S$.

#proposition(title: "Algebraic Reduction to Co-Association Linkage Ensembles", origin: "Original Proof")[
  When $p = 1$, Monti Consensus Clustering is exactly equivalent to average-linkage hierarchical clustering under the squared-Euclidean metric $d_"emb" (i, j)^2 = 2 (1 - cal(M)_(i j))$ induced by the normalized cluster-indicator feature map $phi(i) = (1 / sqrt(H)) H_(i, :)^top$, and differs from CSPA (Alg. 171) exclusively by replacing balanced multi-way graph cut with unconstrained agglomerative UPGMA, thereby eliminating CSPA's cluster-size balance restriction at the cost of greedy hierarchical merging.
]
#proof[
  At $p = 1$, $N_(i j) = H$ identically, so $cal(M)_(i j) = (1 / H) sum_(h=1)^H bb(1)(pi^((h)) (i) = pi^((h)) (j)) = chevron.l phi(i), phi(j) chevron.r$ with $norm(phi(i))_2 = 1$. Hence $norm(phi(i) - phi(j))_2^2 = 2 - 2 chevron.l phi(i), phi(j) chevron.r = 2 (1 - cal(M)_(i j)) = 2 D_(i j)$. Since positive scalar multiplication of a dissimilarity matrix preserves the merge order of UPGMA identically, Monti's final stage is UPGMA on squared feature distances without size constraints, whereas CSPA partitions the identical affinity matrix $S = cal(M)$ under a balanced $k$-way cut constraint @strehl2002cluster.
]

Unlike HGPA (Alg. 172), which operates on binary hyperedge cuts and suffers from single-vertex cut saturation, and MCLA (Alg. 173), which aggregates clusters via a coarse meta-graph, Monti retains full $n times n$ pairwise co-membership resolution and avoids hard cluster-size balancing constraints.

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* On the standardized ClusterOverview benchmark, Monti Consensus Clustering achieves: *Iris* $"Sil" = 0.445$, $"ARI" = 0.634$; *Customers* $"Sil" = 0.332$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.257$ ($"ARI" = "N/A"$); with mean runtime $3.809$ ms. Two empirical patterns directly reflect the theory above. First, Monti attains the *highest Iris ARI* ($0.634$) among all four consensus methods (Algs. 170–173, surpassing CSPA and MCLA at $0.563$ and HGPA at $0.492$ by $+0.071$ to $+0.142$), corroborating Theorem 4 and Proposition 6: unconstrained UPGMA on resampled co-clustering probabilities respects asymmetric boundary structure in the overlapping Iris _versicolor_ / _virginica_ region without forcing equal-size cuts. Second, its slightly lower Silhouette than CSPA on *Multi-Shape* ($0.609$ vs. $0.712$) and *High-D* ($0.257$ vs. $0.303$) alongside higher runtime ($3.809$ ms vs. $2.957$ ms for CSPA and $3.197$ ms for MCLA) reflects both the overhead of repeated resampling ($p < 1$) and the sensitivity of UPGMA to minor non-metric distortions in $1 - cal(M)$ (Theorem 1(b)) when base partitions fluctuate on non-convex or high-dimensional geometries.
