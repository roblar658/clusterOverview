#import "../template/ieee.typ": *

#alg_entry(
  id: 145,
  name: "COSA (Clustering Objects on Subsets of Attributes)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2004",
  paradigm: "Entropy-Weighted Dissimilarity / Soft Subspace Clustering",
  complexity: $O(T n^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

High-dimensional clustering frequently fails when genuine object groupings exist only within distinct, cluster-specific subsets of relevant attributes while the remaining coordinates act as unstructured noise @beyer1999nearest @kriegel2009subspacesurvey. In standard distance-based clustering, all $d$ attributes contribute equally to pairwise dissimilarities, so irrelevant dimensions swamp the signal as $d$ grows. Friedman and Meulman @friedman2004cosa introduced *COSA* (Clustering Objects on Subsets of Attributes) to learn attribute-weighted dissimilarities directly from local neighborhoods, assigning each cluster—and iteratively each individual object—its own probability simplex of attribute weights prior to applying any downstream dissimilarity-based clustering algorithm.

Let $bold(X) = [x_(i k)] in bb(R)^(n times d)$ denote $n$ objects over $d$ numerical attributes. To remove measurement-unit discrepancies, COSA defines the normalized coordinate dissimilarity between objects $i$ and $j$ on attribute $k in {1, dots, d}$ as
$ d_(i j k) = abs(x_(i k) - x_(j k)) / s_k, $
where $s_k > 0$ is a robust attribute dispersion scale (canonically the mean absolute pairwise deviation or trimmed standard deviation along coordinate $k$; categorical attributes use binary mismatch indicators $bb(1)[x_(i k) != x_(j k)]$).

#definition(title: "Entropy-Regularized COSA Criterion and Object-Specific Weights")[
  Let $Delta_d = {bold(w) in bb(R)_+^d : sum_(k=1)^d w_k = 1}$ be the $(d-1)$-dimensional probability simplex, and let $lambda > 0$ be the entropy incentive (homotopy) parameter. Given a partition $cal(C) = {C_1, dots, C_L}$ of ${1, dots, n}$ with positive cluster weights $a_l > 0$ and cluster attribute-weight vectors $bold(W) = (bold(w)_1, dots, bold(w)_L) in Delta_d^L$, the *variational COSA criterion* is
  $ Q_lambda (cal(C), bold(W)) = sum_(l=1)^L a_l sum_(i, j in C_l) [ sum_(k=1)^d w_(l k) d_(i j k) + lambda sum_(k=1)^d w_(l k) log w_(l k) ], $
  with convention $0 log 0 = 0$. In the *cluster-free object-weighted formulation*, each object $i$ maintains a local $K$-nearest-neighbor set $N_K^((t)) (i)$ under current dissimilarities $D^((t))$, local attribute dispersions
  $ S_(i k)^((t)) = 1 / K sum_(j in N_K^((t)) (i)) d_(i j k), $
  Gibbs attribute weights $w_(i k)^((t)) in Delta_d$, and symmetrized pairwise COSA dissimilarity
  $ D_(i j)^((t+1)) = max( D_(i j) (bold(w)_i^((t))), D_(i j) (bold(w)_j^((t))) ), $
  where $D_(i j) (bold(w)) = ( sum_(k=1)^d w_k^p d_(i j k)^p )^(1/p)$ ($p >= 1$, canonically $p = 1$ or $p = 2$, or the shifted soft-min form).
]

The negative Shannon entropy penalty $lambda sum_k w_k log w_k$ prevents degenerate collapse onto a single coordinate with smallest dispersion: $lambda$ smoothly interpolates between hard single-attribute selection ($lambda -> 0^+$) and uniform full-dimensional dissimilarity ($lambda -> infinity$). Because the iterative object-weighted loop only requires pairwise distances, COSA decouples *subspace metric learning* from *cluster extraction*: once the converged $n times n$ matrix $bold(D)^*$ is formed, hierarchical linkage @ward1963hierarchical, PAM @kaufman1990finding, or multidimensional scaling can extract arbitrary cluster topologies.

#algbox(
  title: "Iterative COSA Dissimilarity Learning and Partitioning",
  inputs: [Data matrix $bold(X) in bb(R)^(n times d)$; neighborhood size $K in {1, dots, n-1}$; incentive $lambda > 0$; homotopy schedule $eta_1 >= dots >= eta_T = 1$; norm order $p >= 1$; downstream clustering method $cal(A)$.],
  outputs: [Object weight matrix $bold(W)^* in Delta_d^n$; learned dissimilarity matrix $bold(D)^* in bb(R)_+^(n times n)$; final partition $cal(C)$.],
)[
  1. *Standardize coordinates:* For each $k in {1, dots, d}$, compute robust dispersion $s_k > 0$ and normalized coordinate dissimilarities $d_(i j k) = abs(x_(i k) - x_(j k)) / s_k$.
  2. *Initialize weights and distances:* Set $w_(i k)^((0)) = 1 / d$ for all $i, k$, and compute initial uniform dissimilarities $D_(i j)^((0)) = ( sum_(k=1)^d (w_(i k)^((0)) d_(i j k))^p )^(1/p)$ (or base soft-min $Phi_lambda (bold(d)_(i j))$).
  3. *Iterative local subspace refinement:* For $t = 0, dots, T - 1$ with effective temperature $lambda_t = eta_(t+1) lambda$:
     - For each object $i in {1, dots, n}$, find its $K$ nearest neighbors $N_K^((t)) (i)$ under $bold(D)^((t))$ and evaluate local dispersions $S_(i k)^((t)) = K^(-1) sum_(j in N_K^((t)) (i)) d_(i j k)$.
     - Update object-specific Gibbs weights:
       $ w_(i k)^((t+1)) = exp(-S_(i k)^((t)) / lambda_t) / (sum_(r=1)^d exp(-S_(i r)^((t)) / lambda_t)). $
     - Recompute symmetrized pairwise dissimilarities $D_(i j)^((t+1)) = max( D_(i j) (bold(w)_i^((t+1))), D_(i j) (bold(w)_j^((t+1))) )$.
     - Stop early if $max_(i, k) abs(w_(i k)^((t+1)) - w_(i k)^((t))) < tau_"tol"$ and $eta_(t+1) = 1$.
  4. *Downstream partitioning:* Apply distance-based algorithm $cal(A)$ (e.g., average-linkage hierarchical clustering or PAM) to $bold(D)^*$ to output $cal(C)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish the variational duality of COSA's Gibbs weights, exact KL-quantified monotonic convergence of alternating updates, a constructive proof that COSA dissimilarities violate the triangle inequality, and sharp finite-dimensional guarantees on coordinate invariance and exponential noise suppression.

#subsec("A", "Convex Conjugate Duality and Soft-Minimum Bounds")

#theorem(title: "Exact Variational Dual and Soft-Minimum Sandwich Bounds", origin: "Original Unified Proof")[
  Fix a cluster $C_l$ (or local neighborhood $N_K (i)$) with coordinate dispersions $bold(S)_l = (S_(l 1), dots, S_(l d)) in bb(R)_+^d$, where $S_(l k) = abs(C_l)^(-2) sum_(i, j in C_l) d_(i j k)$, and fix $lambda > 0$. Consider
  $ f_lambda (bold(w)) = sum_(k=1)^d w_k S_(l k) + lambda sum_(k=1)^d w_k log w_k, quad bold(w) in Delta_d. $
  Then: (i) $f_lambda$ is strictly convex on the relative interior $"int"(Delta_d)$ and admits the unique global minimizer
  $ w_(l k)^* = exp(-S_(l k) / lambda) / (sum_(r=1)^d exp(-S_(l r) / lambda)), quad k = 1, dots, d. $
  (ii) The shifted dual optimum $Phi_lambda (bold(S)_l) = min_(bold(w) in Delta_d) f_lambda (bold(w)) + lambda log d$ equals the exact soft-minimum
  $ Phi_lambda (bold(S)_l) = -lambda log ( 1 / d sum_(k=1)^d exp(-S_(l k) / lambda) ), $
  and satisfies the tight sandwich bounds
  $ min_k S_(l k) <= Phi_lambda (bold(S)_l) <= min_k S_(l k) + lambda log d, $
  with limits $lim_(lambda -> 0^+) Phi_lambda (bold(S)_l) = min_k S_(l k)$ and $lim_(lambda -> infinity) Phi_lambda (bold(S)_l) = d^(-1) sum_(k=1)^d S_(l k)$.
]

#proof[
  For any $bold(w) in "int"(Delta_d)$, the Hessian of $f_lambda$ is diagonal with entries $partial^2 f_lambda / (partial w_k^2) = lambda / w_k > 0$, proving strict convexity. Moreover, as $w_k -> 0^+$, $partial f_lambda / (partial w_k) = S_(l k) + lambda (1 + log w_k) -> -infinity$, so the boundary of $Delta_d$ cannot contain a minimizer; the unique global minimum lies in $"int"(Delta_d)$. Forming the Lagrangian $cal(L)(bold(w), nu) = f_lambda (bold(w)) + nu (sum_(k=1)^d w_k - 1)$ and setting $partial cal(L) / (partial w_k) = 0$ gives
  $ S_(l k) + lambda (1 + log w_k^*) + nu = 0 ==> w_k^* = exp(-(lambda + nu) / lambda) exp(-S_(l k) / lambda). $
  Enforcing $sum_(r=1)^d w_r^* = 1$ yields the partition function $Z_l = sum_(r=1)^d exp(-S_(l r) / lambda)$ and the unique Gibbs minimizer $w_(l k)^* = exp(-S_(l k) / lambda) / Z_l$. Substituting $log w_(l k)^* = -S_(l k) / lambda - log Z_l$ back into $f_lambda (bold(w)_l^*)$ cancels the linear term identically:
  $ f_lambda (bold(w)_l^*) = sum_(k=1)^d w_(l k)^* S_(l k) + lambda sum_(k=1)^d w_(l k)^* (-S_(l k) / lambda - log Z_l) = -lambda log Z_l. $
  Adding $lambda log d$ yields $Phi_lambda (bold(S)_l) = -lambda log (Z_l / d)$. Let $m = min_k S_(l k)$. Because each term satisfies $exp(-S_(l k) / lambda) <= exp(-m / lambda)$ and at least one term attains $exp(-m / lambda)$, we have
  $ 1 / d exp(-m / lambda) <= 1 / d sum_(k=1)^d exp(-S_(l k) / lambda) <= exp(-m / lambda). $
  Taking $-lambda log(dot.op)$ reverses inequalities and yields $m <= Phi_lambda (bold(S)_l) <= m + lambda log d$, immediately implying $lim_(lambda -> 0^+) Phi_lambda (bold(S)_l) = m$. For $lambda -> infinity$, Taylor expansion $exp(-S_(l k) / lambda) = 1 - S_(l k) / lambda + O(lambda^(-2))$ gives
  $ Phi_lambda (bold(S)_l) = -lambda log ( 1 - 1 / (lambda d) sum_(k=1)^d S_(l k) + O(lambda^(-2)) ) -> 1 / d sum_(k=1)^d S_(l k), $
  completing the proof.
]

#subsec("B", "KL-Quantified Lyapunov Descent and Convergence")

#theorem(title: "Monotonic Lyapunov Descent and Finite-Step Convergence of Alternating COSA", origin: "Original Proof")[
  Fix $lambda > 0$ and consider alternating minimization of $Q_lambda (cal(C), bold(W))$ over partitions $cal(C)$ of ${1, dots, n}$ into $L$ clusters and weights $bold(W) in Delta_d^L$:
  1. $bold(W)^((t+1)) = arg min_(bold(W) in Delta_d^L) Q_lambda (cal(C)^((t)), bold(W))$;
  2. $cal(C)^((t+1)) = arg min_(cal(C)) Q_lambda (cal(C), bold(W)^((t+1)))$ (or any assignment step satisfying $Q_lambda (cal(C)^((t+1)), bold(W)^((t+1))) <= Q_lambda (cal(C)^((t)), bold(W)^((t+1)))$).
  Then the weight update produces the exact KL-divergence energy drop
  $ Q_lambda (cal(C)^((t)), bold(W)^((t))) - Q_lambda (cal(C)^((t)), bold(W)^((t+1))) = lambda sum_(l=1)^L a_l abs(C_l^((t)))^2 "KL"(bold(w)_l^((t)) || bold(w)_l^((t+1))) >= 0, $
  where $"KL"(bold(p) || bold(q)) = sum_(k=1)^d p_k log (p_k / q_k)$. Consequently, $Q_lambda (cal(C)^((t)), bold(W)^((t)))$ is monotonically non-increasing and terminates in finitely many iterations at a block-coordinate stationary point.
]

#proof[
  Rewrite the cluster-$l$ inner sum as $abs(C_l^((t)))^2 f_lambda^((l)) (bold(w)_l)$, where $f_lambda^((l)) (bold(w)) = sum_(k=1)^d w_k S_(l k)^((t)) + lambda sum_(k=1)^d w_k log w_k$. By Theorem II.A, the exact minimizer $bold(w)_l^((t+1))$ satisfies $S_(l k)^((t)) = -lambda (log w_(l k)^((t+1)) + log Z_l^((t)))$ and $f_lambda^((l)) (bold(w)_l^((t+1))) = -lambda log Z_l^((t))$. For any $bold(w)_l^((t)) in Delta_d$, substituting this identity for $S_(l k)^((t))$ gives
  $ f_lambda^((l)) (bold(w)_l^((t))) - f_lambda^((l)) (bold(w)_l^((t+1))) &= sum_(k=1)^d w_(l k)^((t)) [ -lambda log w_(l k)^((t+1)) - lambda log Z_l^((t)) + lambda log w_(l k)^((t)) ] + lambda log Z_l^((t)) \
  &= lambda sum_(k=1)^d w_(l k)^((t)) log ( w_(l k)^((t)) / w_(l k)^((t+1)) ) = lambda "KL"(bold(w)_l^((t)) || bold(w)_l^((t+1))). $
  Multiplying by $a_l abs(C_l^((t)))^2$ and summing over $l = 1, dots, L$ establishes the exact energy drop identity. By Gibbs' inequality, $"KL"(bold(p) || bold(q)) >= 0$ with equality iff $bold(p) = bold(q)$. Combined with the non-increasing assignment step, $Q_lambda (cal(C)^((t+1)), bold(W)^((t+1))) <= Q_lambda (cal(C)^((t)), bold(W)^((t)))$. Because $bold(W)^((t+1))$ is uniquely determined by $cal(C)^((t))$, the mapping $cal(C) |-> min_bold(W) Q_lambda (cal(C), bold(W))$ takes values on the finite set of $L$-partitions of $n$ objects (cardinality at most $L^n / L!$). Whenever $cal(C)^((t+1)) != cal(C)^((t))$ strictly decreases $Q_lambda$, no previous partition can ever recur, forcing finite termination at a partial optimum.
]

#subsec("C", "Non-Metricity and Triangle-Inequality Distortion")

Although every coordinate dissimilarity $d_(i j k)$ is a metric on $bb(R)$, COSA's attribute-selective aggregation breaks metricity—a decisive geometric property that invalidates triangle-inequality index trees (ball trees, VP-trees, Elkan-type bounds).

#proposition(title: "Unbounded Triangle-Inequality Violation of COSA Dissimilarities", origin: "Original Counterexample")[
  (i) For any $lambda > 0$ and dimension $d >= 2$, the base soft-minimum COSA dissimilarity $D_(i j)^((0)) = Phi_lambda (bold(d)_(i j))$ violates the triangle inequality with unbounded distortion ratio: for every $M > 2 lambda log 2$, there exist three points in $bb(R)^2$ with $D_(1 3)^((0)) / (D_(1 2)^((0)) + D_(2 3)^((0))) > M / (2 lambda log 2) -> infinity$.
  (ii) In $d >= 3$, the symmetrized weighted COSA dissimilarity $D_(i j) = max( chevron.l bold(w)_i, bold(d)_(i j) chevron.r, chevron.l bold(w)_j, bold(d)_(i j) chevron.r )$ with strictly positive object weights $bold(w)_i in "int"(Delta_d)$ likewise violates the triangle inequality with arbitrarily large ratio.
]

#proof[
  (i) Take $d = 2$ with unit scales $s_1 = s_2 = 1$ and points $bold(x)_1 = (0, 0)$, $bold(x)_2 = (0, M)$, $bold(x)_3 = (M, M)$ for $M > 0$. The coordinate dissimilarity vectors are $bold(d)_(1 2) = (0, M)$, $bold(d)_(2 3) = (M, 0)$, and $bold(d)_(1 3) = (M, M)$. Evaluating $Phi_lambda$ directly gives
  $ D_(1 2)^((0)) = D_(2 3)^((0)) = -lambda log ( (1 + exp(-M / lambda)) / 2 ) < lambda log 2, $
  whereas $D_(1 3)^((0)) = -lambda log ( exp(-M / lambda) ) = M$. Hence for any $M >= 2 lambda log 2$,
  $ D_(1 3)^((0)) = M > 2 lambda log 2 > D_(1 2)^((0)) + D_(2 3)^((0)), $
  and letting $M -> infinity$ makes $D_(1 3)^((0)) / (D_(1 2)^((0)) + D_(2 3)^((0))) -> infinity$.
  (ii) Take $d = 3$, $s_1 = s_2 = s_3 = 1$, points $bold(x)_1 = (0, 0, 0)$, $bold(x)_2 = (0, M, 0)$, $bold(x)_3 = (M, M, 0)$, and strictly interior weights for $epsilon in (0, 1/3)$:
  $ bold(w)_1 = (1 - 2 epsilon, epsilon, epsilon), quad bold(w)_2 = (epsilon, epsilon, 1 - 2 epsilon), quad bold(w)_3 = (epsilon, 1 - 2 epsilon, epsilon). $
  Then $bold(d)_(1 2) = (0, M, 0)$ yields $D_(1 2) = max(epsilon M, epsilon M) = epsilon M$; $bold(d)_(2 3) = (M, 0, 0)$ yields $D_(2 3) = max(epsilon M, epsilon M) = epsilon M$; and $bold(d)_(1 3) = (M, M, 0)$ yields $D_(1 3) = (1 - epsilon) M$. Therefore $D_(1 3) / (D_(1 2) + D_(2 3)) = (1 - epsilon) / (2 epsilon) > 1$ for all $epsilon < 1/3$, diverging to $infinity$ as $epsilon -> 0^+$.
]

#subsec("D", "Equivariance and Exponential Noise Suppression")

#theorem(title: "Diagonal-Affine Invariance and Exponential Noise-Attribute Attenuation", origin: "Original Derivation")[
  (i) *Symmetries:* Under any coordinate permutation $bold(P) in S_d$ and positive diagonal affine transformation $bold(x)_i |-> bold(A) bold(x)_i + bold(b)$ ($bold(A) = "diag"(a_1, dots, a_d)$, $a_k > 0$), all normalized dissimilarities $d_(i j k)$, learned weights $bold(w)_i^((t))$, and pairwise COSA distances $D_(i j)^((t))$ are strictly invariant (with $bold(w)_i^((t))$ equivariantly permuted by $bold(P)$). Conversely, COSA is *not* invariant under general orthogonal rotations $bold(R) in O(d)$.
  (ii) *Noise Attenuation:* Suppose a cluster $C_l$ possesses an $s$-element relevant subspace $cal(R)$ with dispersions $S_(l k) <= delta$ for $k in cal(R)$, while the $d - s$ irrelevant noise coordinates $cal(N) = {1, dots, d} without cal(R)$ satisfy $S_(l k) >= Delta > delta$. Then the aggregate weight placed on all noise attributes obeys
  $ W_(cal(N)) = sum_(k in cal(N)) w_(l k)^* <= rho / (1 + rho) < rho, quad rho = ((d - s) / s) exp( - (Delta - delta) / lambda ). $
  In particular, choosing $lambda <= (Delta - delta) / log((d - s) / (epsilon s))$ guarantees $W_(cal(N)) < epsilon$ for any target tolerance $epsilon in (0, 1)$, even when $d >> s$.
]

#proof[
  (i) Because $s_k$ is positively homogeneous of degree $1$ and translation invariant, transforming $tilde(x)_(i k) = a_k x_(i k) + b_k$ rescales $tilde(s)_k = a_k s_k$, so $tilde(d)_(i j k) = a_k abs(x_(i k) - x_(j k)) / (a_k s_k) = d_(i j k)$ identically. Coordinate permutations relabel the index $k$ identically inside every sum over $k = 1, dots, d$, leaving $D_(i j)^((t))$ unchanged. To disprove $O(d)$ invariance, take a $1$-dimensional line cluster in $bb(R)^2$ along $x_2 = 0$: here $S_1 > 0, S_2 = 0$, giving $w_2^* > w_1^*$; rotating by $pi / 4$ produces identical dispersions $tilde(S)_1 = tilde(S)_2$, forcing uniform weights $tilde(w)_1^* = tilde(w)_2^* = 1/2$.
  (ii) Partition $Z_l = Z_(cal(R)) + Z_(cal(N))$, where $Z_(cal(R)) = sum_(k in cal(R)) exp(-S_(l k) / lambda) >= s exp(-delta / lambda)$ and $Z_(cal(N)) = sum_(k in cal(N)) exp(-S_(l k) / lambda) <= (d - s) exp(-Delta / lambda)$. Thus $Z_(cal(N)) / Z_(cal(R)) <= ((d - s) / s) exp(-(Delta - delta) / lambda) = rho$. Because $u |-> u / (1 + u)$ is strictly increasing on $[0, infinity)$,
  $ W_(cal(N)) = (Z_(cal(N)) / Z_(cal(R))) / (1 + Z_(cal(N)) / Z_(cal(R))) <= rho / (1 + rho) < rho. $
  Solving $rho <= epsilon$ for $lambda$ immediately yields $lambda <= (Delta - delta) / log((d - s) / (epsilon s))$.
]

#corollary(title: "High-Dimensional Contrast Restoration via Homotopy Annealing", origin: "Classical — Cited")[
  Under the planted subspace model of Theorem II.D with $d -> infinity$ at fixed signal dimension $s$ and separation $Delta - delta > 0$, standard Euclidean distance suffers distance concentration $D_"max" / D_"min" -> 1$ @beyer1999nearest, whereas COSA @friedman2004cosa with temperature scaling $lambda = O((Delta - delta) / log d)$ preserves $Theta(1)$ relative contrast between within-cluster and cross-cluster dissimilarities.
]

#proof[
  Under uniform weighting, $d - s$ noise dimensions contribute $(d - s) Delta$ to both intra-cluster and inter-cluster distances, driving relative contrast to $O(s / d) -> 0$ @beyer1999nearest. Under COSA with $lambda = (Delta - delta) / (2 log d)$ (for $d > s$), Theorem II.D gives $W_(cal(N)) <= (d / s) exp(-2 log d) = 1 / (s d) -> 0$ while $W_(cal(R)) = 1 - W_(cal(N)) -> 1$, so the weighted dissimilarity converges to a strictly convex combination supported solely on $cal(R)$, preserving $Theta(1)$ separation between clusters that differ on $cal(R)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Computing all pairwise coordinate differences on the fly within each iteration requires $O(n^2 d)$ floating-point operations per sweep, plus $O(n^2 log K)$ (or $O(n^2)$ via quickselect) to extract all $n$ local $K$-nearest-neighbor sets $N_K^((t)) (i)$. Across $T$ homotopy iterations (typically $T in [10, 30]$), the total learning cost is $O(T n^2 d)$, followed by $O(n^2 log n)$ for downstream average-linkage hierarchical clustering or $O(n^2 L)$ for PAM. Crucially, the $O(n^2 d)$ 3-tensor $d_(i j k)$ need never be materialized in RAM: streaming over coordinate blocks keeps peak memory at $O(n^2 + n d)$ for the pairwise dissimilarity matrix $bold(D)$ and object weights $bold(W)$.

Geometrically, each object $i$ induces a personalized weighted $L_p$ ball aligned with the coordinate axes; symmetrizing via $max(D_(i j)(bold(w)_i), D_(i j)(bold(w)_j))$ enforces mutual subspace consistency—two objects are judged close only if they agree on *both* objects' locally relevant feature subsets. Because Proposition II.C shows that $bold(D)^*$ violates the triangle inequality, spatial index structures cannot accelerate exact $K$-NN queries, making quadratic $O(n^2)$ scaling the chief computational constraint.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T n^2 d)$ dissimilarity learning + $O(n^2 log n)$ downstream clustering],
    [Space Complexity], [$O(n^2 + n d)$ (streaming $d_(i j k)$ avoids $O(n^2 d)$ tensor storage)],
    [Metric Structure], [Symmetric semimetric; violates triangle inequality (Prop. II.C)],
    [Invariances], [Coordinate-wise positive affine & permutation invariant; not $O(d)$ invariant],
    [Key Hyperparameters], [Neighborhood size $K$, entropy temperature $lambda$, homotopy schedule $eta_t$],
    [Cluster Geometry], [Axis-aligned soft subspaces; arbitrary within-subspace shapes via $cal(A)$],
  ),
  caption: [Operational and geometric profile of COSA (Alg. 145).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Transcriptomics and Disease Subtyping:* In microarray and single-cell RNA-seq matrices ($d approx 10^3$–$10^4$ genes, $n approx 10^2$–$10^3$ patients), distinct cancer subtypes share tight co-expression across compact 20–100 gene signaling pathways while remaining heterogeneous across the rest of the transcriptome @friedman2004cosa. COSA's exponential weight attenuation (Theorem II.D) isolates subtype-specific pathways automatically without pre-filtering genes globally.
2. *Cheminformatics and Molecular QSAR Library Grouping:* Compound databases encode hundreds of mixed-scale physicochemical, topological, and pharmacophore descriptors. Because different scaffold classes share activity along distinct subsets of descriptors (e.g., hydrogen-bond donors vs. lipophilicity), COSA's scale-invariant object weighting groups bioactive chemotypes cleanly across heterogeneous measurement scales.
3. *High-Dimensional Customer and Market Micro-Segmentation:* Retail CRM records capture spending across dozens of disjoint product categories; niche customer cohorts behave identically within 3–5 signature categories while showing noisy variance elsewhere. Object-level weights $bold(w)_i$ allow overlapping category relevance without forcing a single global feature selection.
4. *Multi-Sensor Industrial Fault Diagnostics:* Rotating machinery and semiconductor fabs monitor hundreds of vibration, acoustic, and thermal channels. Distinct fault mechanisms produce tight signatures on small channel subsets; COSA suppresses nominal background channels and enables hierarchical fault-tree taxonomy discovery.
5. *Metabolomics and Proteomics Biomarker Discovery:* Liquid-chromatography mass-spectrometry profiles contain abundant uninformative background ions and skewed intensities. Robust dispersion scaling $s_k$ paired with soft-min dissimilarity recovers phenotype-specific metabolite panels even at modest sample sizes $n$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

COSA degrades sharply in four regimes: (1) *Oblique/correlated subspaces*, where clusters lie on tilted hyperplanes not aligned with coordinate axes (where ORCLUS @aggarwal2000orclus or Sparse Subspace Clustering @elhamifar2013ssc is required); (2) *Large-$n$ datasets* ($n > 2 times 10^4$), where the $O(n^2)$ dissimilarity matrix becomes prohibitive (favoring P3C @moise2006p3c or EWKM @huang2005automated); (3) *Miscalibrated temperature $lambda$*, where $lambda -> 0^+$ latches onto single-coordinate ties or measurement artifacts while overly large $lambda$ reverts to unweighted clustering; and (4) *Extreme cluster-size imbalance*, where a global neighborhood size $K$ larger than the smallest cluster pulls foreign points into $N_K (i)$, contaminating $S_(i k)$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Genomics / Omics], [High: small $n$, huge $d$, axis-aligned gene modules], [Watch small-$K$ variance; Biclustering @cheng2000biclustering],
    [Cheminformatics QSAR], [High: mixed descriptor scales, local scaffolds], [Tune $lambda$ via stability; PROCLUS @aggarwal1999proclus],
    [Large Web / Telemetry], [Low: $O(T n^2 d)$ runtime and $O(n^2)$ RAM], [Use P3C @moise2006p3c or EWKM @huang2005automated],
    [Rotated Linear Manifolds], [Low: strictly axis-aligned coordinate weights], [Use ORCLUS @aggarwal2000orclus or SSC @elhamifar2013ssc],
    [Mixed-Shape Subspaces], [High: pairs with linkage or density clustering], [Non-metric $bold(D)^*$ disables ball trees],
  ),
  caption: [Domain suitability and failure-mode matrix for COSA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

COSA sits at the nexus between centroid-based weighted $k$-means and hard projected clustering, linked by exact limiting reductions.

#proposition(title: "Variational Reductions to EWKM and Hard Projected Clustering", origin: "Original Derivation")[
  (i) Replacing pairwise $L_1$ dispersions $sum_(i, j in C_l) d_(i j k)$ in $Q_lambda (cal(C), bold(W))$ by squared Euclidean within-cluster variances $sum_(i in C_l) (x_(i k) - mu_(l k))^2$ recovers the exact Entropy-Weighting $k$-Means (EWKM) objective @huang2005automated, preserving the identical Gibbs weight update with $S_(l k) = sum_(i in C_l) (x_(i k) - mu_(l k))^2$.
  (ii) Adding the box constraint $0 <= w_(l k) <= 1 / s_l$ to the simplex $Delta_d$ and letting $lambda -> 0^+$ converges to uniform $1 / s_l$ weights on the $s_l$ lowest-dispersion attributes of cluster $C_l$, recovering the hard axis-aligned subspace selection step of PROCLUS @aggarwal1999proclus.
]

#proof[
  (i) Setting $d_(i j k) = (x_(i k) - x_(j k))^2$ and $a_l = (2 abs(C_l))^(-1)$ satisfies the algebraic identity $(2 abs(C_l))^(-1) sum_(i, j in C_l) (x_(i k) - x_(j k))^2 = sum_(i in C_l) (x_(i k) - mu_(l k))^2$, where $mu_(l k) = abs(C_l)^(-1) sum_(i in C_l) x_(i k)$. Substituting into $Q_lambda$ yields $sum_(l=1)^L [ sum_(k=1)^d w_(l k) sum_(i in C_l) (x_(i k) - mu_(l k))^2 + tilde(lambda)_l sum_(k=1)^d w_(l k) log w_(l k) ]$, which is identically the EWKM objective @huang2005automated.
  (ii) Minimizing $sum_k w_k S_(l k) + lambda sum_k w_k log w_k$ over $Delta_d inter [0, 1 / s_l]^d$ as $lambda -> 0^+$ converges to the linear program $min { sum_k w_k S_(l k) : sum_k w_k = 1, 0 <= w_k <= 1 / s_l }$, whose extreme-point optimum assigns $w_k = 1 / s_l$ to the $s_l$ coordinates with smallest $S_(l k)$ and $0$ elsewhere—matching PROCLUS's $s_l$-dimensional attribute assignment @aggarwal1999proclus.
]

COSA also relates closely to Minkowski Weighted $k$-Means @deamorim2012minkowski (which replaces entropy regularization with power-law weights $w_(l k)^beta$) and to FIRES (Alg. 146) @kriegel2005fires, which replaces continuous soft weights with discrete 1D base-cluster merging.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized ClusterOverview suite evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing, COSA achieves *Iris* $"Sil" = 0.435$, $"ARI" = 0.727$; *Customers* $"Sil" = 0.327$; *Multi-Shape* $"Sil" = 0.609$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.259$ ($"ARI" = "N/A"$), with a mean runtime of *12.455 ms*. Several mathematical signatures emerge directly from our theory: first, pairing COSA's learned dissimilarity with flexible downstream partitioning produces strong non-convex geometry recovery (*Multi-Shape* $"Sil" = 0.609$, tying FIRES and outperforming both INSCY at $0.554$ and P3C at $0.557$) alongside competitive supervised agreement on *Iris* ($"ARI" = 0.727$). Second, on *High-D* ($"Sil" = 0.259$), soft full-space participation ($w_(i k) > 0$ on all $d$ coordinates for finite $lambda$) retains residual noise tail mass relative to INSCY's hard maximal-subspace pruning ($0.336$) and P3C's Poisson-filtered cores ($0.302$). Finally, its 12.455 ms mean runtime reflects the $O(T n^2 d)$ iterative dissimilarity matrix computation—comparable to FIRES (15.611 ms) and INSCY (14.553 ms), yet two orders of magnitude slower than P3C's histogram-based pipeline (0.107 ms).
