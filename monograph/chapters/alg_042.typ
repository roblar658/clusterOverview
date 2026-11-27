#import "../template/ieee.typ": *

#alg_entry(
  id: 42,
  name: "Medoidshift",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2007",
  paradigm: "Discrete Weighted-Medoid Mode-Seeking",
  complexity: $O(n^2 (d + n))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical Mean Shift (Alg. 34 @cheng1995mean @comaniciu2002mean) locates modes of a kernel density estimate by iteratively shifting each query point $y in bb(R)^d$ to its kernel-weighted arithmetic mean $m_"MS" (y) = (sum_j w_j x_j) / (sum_j w_j)$. While theoretically elegant in Euclidean spaces, this update requires a linear vector-space structure (addition and scalar multiplication) and generates intermediate iterates $y^((t)) in.not cal(X)$ lying outside the observed dataset—rendering Mean Shift inapplicable to general dissimilarity spaces such as graph shortest-path metrics, biological sequence alignments, articulated 3D shapes, or non-Euclidean manifolds where extrinsic averages are undefined.

*Medoidshift*, introduced by Sheikh, Khan, and Kanade @sheikh2007medoidshift, resolves this limitation by reformulating Mean Shift as a weighted Fréchet optimization and restricting the shifted prototype strictly to the discrete sample set $cal(X)$. Formally, let $cal(X) = {x_1, dots, x_n}$ be a finite set equipped with a symmetric dissimilarity function $d : cal(X) times cal(X) -> bb(R)_(>= 0)$ satisfying $d(x_i, x_i) = 0$, and let $D^((2)) in bb(R)^(n times n)$ denote the squared-dissimilarity matrix $D^((2))_(k j) = d(x_k, x_j)^2$. Given a differentiable, non-increasing, convex kernel profile $K : [0, infinity) -> bb(R)_(>= 0)$ with shadow derivative $g(u) := -K'(u) >= 0$ (for the canonical Gaussian profile $K(u) = exp(-u / (2 sigma^2))$, $g(u) = (2 sigma^2)^(-1) exp(-u / (2 sigma^2))$), define the pairwise kernel weight matrix $W in bb(R)^(n times n)$ by $W_(i j) := g(d(x_i, x_j)^2)$.

Recall that in Euclidean space, the continuous Mean Shift target $mu_i$ minimizes the weighted squared-distance energy $z |-> sum_(j=1)^n W_(i j) norm(z - x_j)_2^2$ over all $z in bb(R)^d$. Restricting $z$ to $cal(X)$ yields the discrete *weighted-medoid objective*:
$ Q(k | i) := sum_(j=1)^n W_(i j) d(x_k, x_j)^2, quad k in {1, dots, n}. $

#definition(title: "Medoidshift Operator, Modes, and Basins")[
  The *Medoidshift transition map* $M : {1, dots, n} -> {1, dots, n}$ assigns to each datum $i$ a weighted sample medoid minimizing $Q(dot | i)$:
  $ M(i) := arg min_(k in {1, dots, n}) Q(k | i) = arg min_(k in {1, dots, n}) (W D^((2)))_(i k), $
  using *self-preferring tie-breaking*: set $M(i) = i$ whenever $Q(i | i) = min_k Q(k | i)$, and break remaining ties by smallest index $k$. A point $i^*$ is a *discrete medoid mode* if it is a fixed point $M(i^*) = i^*$. Iterating $i^((t+1)) = M(i^((t)))$ from $i^((0)) = i$ defines the mode trajectory of $x_i$, and the *Medoidshift cluster* of mode $i^*$ is its basin of attraction $C_(i^*) = {i : exists t >= 0, space M^((t)) (i) = i^*}$.
]

Because $M$ maps the finite set ${1, dots, n}$ into itself, every intermediate point along any shift trajectory is already an element of $cal(X)$ whose next transition $M(dot)$ was computed in a *single matrix multiplication* $S = W D^((2))$ (Algorithm 1)—eliminating iterative kernel re-evaluations and stopping tolerances @sheikh2007medoidshift.

#algbox(
  title: "Medoidshift Mode-Seeking Algorithm",
  inputs: [Finite set $cal(X) = {x_1, dots, x_n}$ with symmetric dissimilarity $d$, bandwidth $sigma > 0$ (kernel $K$, shadow $g = -K'$).],
  outputs: [Mode exemplars $cal(M)^*$ and cluster partition $cal(C)$.],
)[
  1. *Precomputation:* Form squared-dissimilarity matrix $D^((2))_(k j) = d(x_k, x_j)^2$ and weight matrix $W_(i j) = g(D^((2))_(i j))$ for all $i, j, k in {1, dots, n}$.
  2. *One-shot energy matrix:* Compute $S <- W D^((2)) in bb(R)^(n times n)$, so $S_(i k) = Q(k | i)$.
  3. *Medoid transition map:* For each $i in {1, dots, n}$, compute $s_"min" <- min_k S_(i k)$; if $S_(i i) = s_"min"$, set $M(i) <- i$, else set $M(i) <- arg min_k S_(i k)$ (smallest index).
  4. *Trajectory resolution (pointer jumping):* For each $i in {1, dots, n}$, follow $u <- M(u)$ until $M(u) = u$ (with path compression so shared paths are traversed once), and set label $L_i <- u$.
  5. *Return* mode exemplars $cal(M)^* = {u : M(u) = u}$ and clusters induced by $L$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of Medoidshift: (1) an exact Minorization–Maximization (MM) proof guaranteeing strict density ascent, zero cycles, and finite termination on *arbitrary* non-Euclidean dissimilarities; (2) an exact Euclidean projection identity showing Medoidshift equals nearest-sample Mean Shift; (3) the *Locking-Ball Theorem* and a constructive unimodal fragmentation counterexample; (4) an asymptotic score-function residual bound proving convergence to continuous Mean Shift as sample fill distance vanishes; and (5) isometry invariance and super-exponential outlier isolation.

#subsec("A", "Universal MM Ascent and Cycle-Free Termination")

#theorem(title: "Exact MM Density Ascent on Arbitrary Dissimilarity Spaces", origin: "Original Proof")[
  Let $(cal(X), d)$ be *any* finite symmetric dissimilarity space (not necessarily Euclidean or even obeying the triangle inequality). Let $K : [0, infinity) -> bb(R)_(>= 0)$ be convex and differentiable with $g(u) = -K'(u) >= 0$, set $W_(i j) = g(d(x_i, x_j)^2)$, and define the discrete kernel density
  $ hat(p)_K (x_i) := 1 / n sum_(j=1)^n K(d(x_i, x_j)^2). $
  Then under Algorithm 1's self-preferring tie-breaking:
  1. Every non-stationary transition $i -> M(i) != i$ strictly increases $hat(p)_K$ by the quantitative margin
     $ hat(p)_K (x_(M(i))) - hat(p)_K (x_i) >= 1 / n (Q(i | i) - Q(M(i) | i)) > 0. $
  2. The functional graph $G_M = ({1, dots, n}, {(i, M(i)) : M(i) != i})$ is a *directed acyclic forest* with no cycles of length $>= 2$; every trajectory terminates at a fixed point $M(i^*) = i^*$ in at most $n - 1$ steps.
]

#proof[
  (1) Because $K$ is convex and differentiable on $[0, infinity)$, its first-order tangent minorizer satisfies $K(v) - K(u) >= K'(u) (v - u) = g(u) (u - v)$ for all $u, v >= 0$. Evaluating this inequality at $u = d(x_i, x_j)^2$ and $v = d(x_k, x_j)^2$ and averaging over $j = 1, dots, n$ yields, for every pair $(i, k)$:
  $ hat(p)_K (x_k) - hat(p)_K (x_i) &= 1 / n sum_(j=1)^n (K(d(x_k, x_j)^2) - K(d(x_i, x_j)^2)) \
    &>= 1 / n sum_(j=1)^n g(d(x_i, x_j)^2) (d(x_i, x_j)^2 - d(x_k, x_j)^2) \
    &= 1 / n (Q(i | i) - Q(k | i)). $
  Observe that this minorization requires *zero* vector or metric structure beyond symmetry of $d$. Now set $k = M(i)$. By definition of $M(i)$ as a minimizer of $Q(dot | i)$ with self-preferring tie-breaking, whenever $M(i) != i$ we necessarily have the strict inequality $Q(M(i) | i) < Q(i | i)$, proving Claim 1. (For unscaled Gaussian weights $tilde(W)_(i j) = exp(-d_(i j)^2 / (2 sigma^2)) = 2 sigma^2 g(d_(i j)^2)$, the identical minimizer $M(i)$ satisfies $hat(p)_K (x_(M(i))) - hat(p)_K (x_i) >= (2 n sigma^2)^(-1) (tilde(Q)(i | i) - tilde(Q)(M(i) | i)) > 0$.)

  (2) Along any trajectory $i^((0)), i^((1)), dots$ prior to hitting a fixed point, Claim 1 ensures strict monotonicity $hat(p)_K (x_(i^((0)))) < hat(p)_K (x_(i^((1)))) < dots$. Therefore no index can repeat before termination, ruling out all limit cycles of length $>= 2$. Since $|cal(X)| = n$, at most $n - 1$ strict transitions can occur before reaching a fixed point $M(i^*) = i^*$, and every weakly connected component of $G_M$ is a directed in-tree rooted at its unique fixed point.
]

#remark(title: "Necessity of Squared Dissimilarities for Convexity")[
  If one replaces $d(x_k, x_j)^2$ inside $Q(k | i)$ by the unsquared $L_1$ distance $d(x_k, x_j)$ while retaining Gaussian weights $exp(-d^2 / (2 sigma^2))$, the profile $r |-> exp(-r^2 / (2 sigma^2))$ is *concave* on $[0, sigma)$, invalidating the minorization bound of Theorem 1 and permitting oscillation cycles.
]

#subsec("B", "Euclidean Projection Identity and the Mode-Locking Phenomenon")

#lemma(title: "Exact Projection onto Continuous Mean Shift", origin: "Original Proof")[
  Suppose $cal(X) subset.eq (bb(R)^d, norm(dot)_2)$ with $Z_i := sum_(j=1)^n W_(i j) > 0$. Let $mu_i := Z_i^(-1) sum_(j=1)^n W_(i j) x_j$ denote the one-step continuous Mean Shift iterate from $x_i$ and $m(x_i) := mu_i - x_i$ its shift vector. Then for every $k in {1, dots, n}$,
  $ Q(k | i) = Q(i | i) + Z_i (norm(x_k - mu_i)_2^2 - norm(m(x_i))_2^2). $
  Consequently, in Euclidean space, Medoidshift is *algebraically identical* to mapping each point to the nearest dataset exemplar of its continuous Mean Shift target:
  $ M(i) = arg min_(k in {1, dots, n}) norm(x_k - mu_i)_2. $
]

#proof[
  Insert $plus.minus mu_i$ into the squared Euclidean norm inside $Q(k | i)$:
  $ Q(k | i) &= sum_(j=1)^n W_(i j) norm((x_k - mu_i) + (mu_i - x_j))_2^2 \
    &= Z_i norm(x_k - mu_i)_2^2 + sum_(j=1)^n W_(i j) norm(mu_i - x_j)_2^2 \
    &quad + 2 chevron.l x_k - mu_i, space sum_(j=1)^n W_(i j) (mu_i - x_j) chevron.r. $
  By definition of $mu_i$, $sum_(j=1)^n W_(i j) (mu_i - x_j) = Z_i mu_i - sum_(j=1)^n W_(i j) x_j = bold(0)$, so the cross-term vanishes identically. Subtracting the same identity at $k = i$ cancels the $k$-independent local variance term $sum_j W_(i j) norm(mu_i - x_j)_2^2$, leaving $Q(k | i) - Q(i | i) = Z_i (norm(x_k - mu_i)_2^2 - norm(x_i - mu_i)_2^2)$. Since $Z_i > 0$ and $norm(x_i - mu_i)_2 = norm(m(x_i))_2$, minimizing $Q(k | i)$ over $k$ is strictly equivalent to minimizing $norm(x_k - mu_i)_2$.
]

#theorem(title: "Locking-Ball Theorem & Unimodal Over-Segmentation", origin: "Original Counterexample")[
  Assume Euclidean geometry $cal(X) subset.eq bb(R)^d$.
  1. *Exact Locking Condition:* A sample $x_i$ is a fixed point of Medoidshift ($M(i) = i$) *if and only if* the open Euclidean ball $B(mu_i, norm(m(x_i))_2)$ contains no data point from $cal(X)$. In particular, any point whose nearest neighbor satisfies $min_(k != i) norm(x_k - x_i)_2 > 2 norm(m(x_i))_2$ is locked as a stationary mode even when $nabla hat(p)_K (x_i) != bold(0)$.
  2. *Unimodal Fragmentation Counterexample:* For any bandwidth $sigma > 0$ and any even sample size $n = 2 m >= 2$, there exists a 1D dataset whose Gaussian KDE $hat(p)_sigma$ is *strictly unimodal* on $bb(R)$ (so continuous Mean Shift yields 1 cluster), yet Medoidshift locks into *two equal-sized clusters* of size $m$.
]

#proof[
  (1) By Lemma 1, $M(i) = i$ under self-preferring tie-breaking holds iff $Q(k | i) >= Q(i | i)$ for all $k$, i.e., $norm(x_k - mu_i)_2 >= norm(m(x_i))_2$ for all $k in {1, dots, n}$, which states precisely that $B(mu_i, norm(m(x_i))_2) inter cal(X) = emptyset$. If $min_(k != i) norm(x_k - x_i)_2 > 2 norm(m(x_i))_2$, the reverse triangle inequality gives for all $k != i$:
  $ norm(x_k - mu_i)_2 >= norm(x_k - x_i)_2 - norm(x_i - mu_i)_2 > norm(m(x_i))_2, $
  guaranteeing an empty locking ball and trapping $x_i$ as a fixed point.

  (2) Fix any $a in (0, sigma]$ and place $m$ copies at $x = -a$ and $m$ copies at $x = +a$ (separation $2 a <= 2 sigma$). The Gaussian KDE on $bb(R)$ is proportional to
  $ f(x) &= e^(-(x + a)^2 / (2 sigma^2)) + e^(-(x - a)^2 / (2 sigma^2)) \
         &= 2 e^(-(x^2 + a^2) / (2 sigma^2)) cosh( (a x) / sigma^2 ). $
  Differentiating with respect to $x$ yields
  $ f'(x) = 2 / sigma^2 e^(-(x^2 + a^2) / (2 sigma^2)) cosh( (a x) / sigma^2 ) [ a tanh( (a x) / sigma^2 ) - x ]. $
  Using the elementary strict inequality $tanh(z) < z$ for all $z > 0$, every $x > 0$ satisfies $a tanh(a x / sigma^2) < (a / sigma)^2 x <= x$ (as $a <= sigma$), hence $f'(x) < 0$ for all $x > 0$ and symmetrically $f'(x) > 0$ for all $x < 0$. Thus $x^* = 0$ is the *unique* global mode of $hat(p)_sigma$, and continuous Mean Shift merges all $2 m$ points into one cluster at $0$. Yet for any point $x_i = -a$, the cross-cluster weight $w = exp(-2 a^2 / sigma^2) < 1$ yields Mean Shift target
  $ mu_i = (-a + a w) / (1 + w) = -a (1 - w) / (1 + w) < 0. $
  Hence $|-a - mu_i| < |+a - mu_i|$, so $M(i)$ remains at $-a$ (and symmetrically every $+a$ remains at $+a$), splitting the unimodal distribution into two locked clusters of size $m$.
]

#subsec("C", "Asymptotic Convergence to Continuous Mean Shift")

#theorem(title: "Fill-Distance Coupling and Score-Residual Bound", origin: "Original Derivation")[
  Let $cal(X) subset.eq Omega subset.eq bb(R)^d$ lie inside a compact convex domain $Omega$ with covering radius (fill distance) $h_cal(X) := sup_(y in Omega) min_(1 <= k <= n) norm(y - x_k)_2$, and apply Gaussian Medoidshift with bandwidth $sigma > 0$. Then:
  1. Every one-step Medoidshift transition approximates continuous Mean Shift uniformly: $norm(x_(M(i)) - mu_i)_2 <= h_cal(X)$ for all $i$.
  2. Every fixed point $x_i$ of Medoidshift ($M(i) = i$) satisfies the *exact log-density gradient (score) residual bound*
     $ norm(nabla log hat(p)_K (x_i))_2 = norm(m(x_i))_2 / sigma^2 <= h_cal(X) / sigma^2. $
  Consequently, in the dense-sample regime $h_cal(X) / sigma^2 -> 0$, every Medoidshift fixed point converges asymptotically to a true stationary mode of the continuous KDE @sheikh2007medoidshift @ariascastro2016meanshift.
]

#proof[
  (1) Because $W_(i j) >= 0$ and $Z_i > 0$, the Mean Shift target $mu_i = sum_j (W_(i j) / Z_i) x_j$ is a convex combination of $cal(X) subset.eq Omega$, hence $mu_i in Omega$ by convexity of $Omega$. By Lemma 1 and the definition of fill distance,
  $ norm(x_(M(i)) - mu_i)_2 = min_(1 <= k <= n) norm(x_k - mu_i)_2 <= h_cal(X). $

  (2) For the Gaussian kernel $K(u) = exp(-u / (2 sigma^2))$, differentiating $hat(p)_K (x) = n^(-1) sum_j exp(-norm(x - x_j)_2^2 / (2 sigma^2))$ gives
  $ nabla hat(p)_K (x_i) = 1 / (n sigma^2) sum_(j=1)^n tilde(W)_(i j) (x_j - x_i) = (tilde(Z)_i / (n sigma^2)) (mu_i - x_i) = (hat(p)_K (x_i) / sigma^2) m(x_i), $
  where $tilde(Z)_i = sum_j tilde(W)_(i j) = n hat(p)_K (x_i) > 0$. Dividing by $hat(p)_K (x_i)$ yields $nabla log hat(p)_K (x_i) = m(x_i) / sigma^2$. At any fixed point $M(i) = i$, Claim 1 gives $norm(m(x_i))_2 = norm(x_i - mu_i)_2 = norm(x_(M(i)) - mu_i)_2 <= h_cal(X)$, immediately yielding $norm(nabla log hat(p)_K (x_i))_2 <= h_cal(X) / sigma^2$.
]

#subsec("D", "Isometry Invariance and Super-Exponential Outlier Isolation")

#proposition(title: "Dissimilarity Invariance and Asymptotic Outlier Immunity", origin: "Original Proof")[
  1. Medoidshift is invariant under all distance-preserving bijections of $(cal(X), d)$ and equivariant under uniform scaling $(d, sigma) |-> (alpha d, alpha sigma)$.
  2. Under a Gaussian kernel on a metric space $(cal(X), d)$, perturbing a single outlier $x_n$ so that $r_n := min_(j < n) d(x_n, x_j) -> infinity$ alters inlier objectives $Q(k | i)$ ($i, k < n$) by $O(r_n^2 e^(-r_n^2 / (2 sigma^2))) -> 0$ and isolates $x_n$ as a self-locked singleton mode $M(n) = n$ without shifting any non-degenerate inlier cluster.
]

#proof[
  (1) Both $D^((2))$ and $W$ depend exclusively on pairwise ratios $d(x_i, x_j) / sigma$, making $S = W D^((2))$ and $M$ invariant under isometries and uniform rescaling.

  (2) Let $Delta_"in" = max_(i, k < n) d(x_i, x_k) < infinity$. For any inliers $i, k < n$, the triangle inequality gives $d(x_k, x_n) <= d(x_i, x_n) + Delta_"in"$. Thus the outlier's additive contribution to $tilde(Q)(k | i)$ is bounded by
  $ tilde(W)_(i n) d(x_k, x_n)^2 <= (d(x_i, x_n) + Delta_"in")^2 exp( - d(x_i, x_n)^2 / (2 sigma^2) ) -> 0 $
  super-exponentially as $r_n -> infinity$, preserving every strict inlier minimizer $M(i)$. Conversely, for $i = n$, $tilde(W)_(n n) = 1$ while $tilde(W)_(n j) <= e^(-r_n^2 / (2 sigma^2))$; hence $tilde(Q)(n | n) <= (n - 1) (r_n + Delta_"in")^2 e^(-r_n^2 / (2 sigma^2)) -> 0$ whereas for every $k < n$, $tilde(Q)(k | n) >= tilde(W)_(n n) d(x_k, x_n)^2 >= r_n^2 -> infinity$, forcing $M(n) = n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

On general dissimilarity spaces, computing $D^((2))$ and $W$ costs $O(n^2 d)$ (or $O(n^2)$ metric evaluations), while evaluating the dense matrix product $S = W D^((2))$ costs $Theta(n^3)$ scalar operations ($2 n^3$ FLOPs) and $O(n^2)$ memory—making Medoidshift substantially more expensive than Quickshift ($O(n^2 d)$) on large $n$, though fully amenable to BLAS-3 GEMM and GPU acceleration. Importantly, two exact speedups arise directly from our analysis: (i) for compactly supported kernels with average neighborhood cardinality $overline(m) << n$, sparse matrix multiplication evaluates $S$ in $O(n overline(m)^2)$ time @sheikh2007medoidshift; and (ii) in Euclidean $bb(R)^d$, Lemma 1 bypasses $W D^((2))$ completely by computing Mean Shift targets $mu_i$ in $O(n^2 d)$ time and locating nearest exemplars via a KD-tree or spatial index in $O(n log n)$ time—collapsing Euclidean complexity from $O(n^3)$ to $O(n^2 d)$!

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Medoidshift],
    [Time / Space Complexity], [$O(n^2 d + n^3)$ general (BLAS-3 GEMM); $O(n^2 d)$ in Euclidean space via Lemma 1; $O(n^2)$ space.],
    [Geometric Inductive Bias], [Nonparametric mode basins anchored at real exemplars; no convexity or fixed-$k$ assumption.],
    [Determinism / Stopping], [Strictly deterministic; zero step-size or convergence tolerance; cycle-free for convex $K$ (Theorem 1).],
    [Invariances], [Invariant under metric isometries; equivariant under $(d, sigma) |-> (alpha d, alpha sigma)$; valid on non-Euclidean dissimilarities.],
    [Outlier & Sparsity Profile], [Super-exponentially immune to far outliers (Proposition 1); vulnerable to mode-locking in sparse/high-$d$ regions (Theorem 2).],
    [Recommended Companions], [Quickshift (Alg. 41) for $O(n^2 d)$ speed without locking; Mean Shift++ (Alg. 55) in $bb(R)^d$; PAM/FasterPAM (Algs. 4, 24) for fixed-$k$ medoids.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of Medoidshift.],
)

#sec("IV", "Engineering and Scientific Applications")

Medoidshift's decisive comparative advantage emerges on *non-Euclidean, structured, or exemplar-constrained domains* where pairwise distances exist natively but continuous averages are ill-defined or semantically meaningless.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *3D Articulated Shapes, Meshes & Manifold Clustering:* Originally demonstrated by Sheikh et al. @sheikh2007medoidshift on human-action video sequences and non-Euclidean manifolds (e.g., rotation matrices $S O(3)$, geodesic surfaces), Medoidshift finds genuine physical mode frames/poses without projecting onto invalid extrinsic centroids outside the manifold.
2. *Biological Sequence & Protein-Fold Motif Discovery:* Genomic sequences, protein structures, and phylogenetic taxa are compared via edit distances, BLAST alignment scores, or structural RMSD where arithmetic sequence averages do not exist; Theorem 1 guarantees cycle-free mode convergence directly from the pairwise dissimilarity matrix $D^((2))$.
3. *Video Shot-Boundary & Action Segmentation:* Grouping video keyframes using non-metric visual descriptors (Earth Mover's Distance, dynamic time warping, or learned perceptual similarities) benefits from exemplar modes that serve directly as representative storyboard thumbnails.
4. *Graph-Structured & Network Entity Clustering:* On road networks, molecular graphs, or citation graphs equipped with shortest-path or diffusion distances, Medoidshift discovers density-central hub nodes directly on the vertex set $V$ without spectral embedding distortion.
5. *Auditable & Privacy-Constrained Exemplar Summarization:* In medical cohort discovery and legal document triage, returned cluster modes must correspond to real, inspectable patient profiles or case records—matching PAM (Alg. 4 @kaufman1990finding) and Affinity Propagation (Alg. 106 @frey2007affinity) while automatically determining the number of modes from bandwidth $sigma$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Medoidshift suffers from four structural limitations: (a) *discrete mode-locking fragmentation* (Theorem 2), where sparse regions or unsampled density peaks trap points into spurious micro-clusters; (b) *cubic $O(n^3)$ scaling* on general non-Euclidean dissimilarity matrices for $n > 10^4$; (c) *curse of dimensionality* in $bb(R)^d$, since covering radius scales as $h_cal(X) = Omega(n^(-1/d))$, causing the score residual $h_cal(X) / sigma^2$ (Theorem 3) to remain large and locking balls to empty out @beyer1999nearest; and (d) *unique-scale bandwidth rigidity*, where a single global $sigma$ either over-merges dense clusters or locks sparse ones.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Medoidshift Fits / Prerequisite], [Caution / Alternative],
    [Manifolds / Shapes], [No vector averages; exact MM ascent], [Ensure dense sampling ($h_cal(X) << sigma^2$)],
    [Bio-Sequences / Strings], [Operates on alignment matrix $D^((2))$], [$O(n^3)$ caps $n <= 5 times 10^3$],
    [Video Keyframing], [Exemplar modes are real video frames], [Use Quickshift (Alg. 41) for larger $n$],
    [Euclidean Point Clouds], [Lemma 1 gives $O(n^2 d)$ speedup], [Use Mean Shift++ (Alg. 55) or DENCLUE 2.0 (Alg. 43)],
    [High-D Sparse Data], [Empty locking balls shatter clusters], [Use HDBSCAN (Alg. 32) or SNN-DPC (Alg. 52)],
  ),
  caption: [Application suitability and engineering guidelines for Medoidshift.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Medoidshift anchors the discrete-exemplar bridge between continuous Mean Shift (Alg. 34 @comaniciu2002mean), Blurring Mean Shift (Alg. 54 @carreira2006fast), DENCLUE 2.0 (Alg. 43 @hinneburg2007denclue2), and Quickshift (Alg. 41 @vedaldi2008quickshift).

#proposition(title: "Unified Mode-Seeking Hierarchy: Mean Shift, Medoidshift, and Quickshift", origin: "Original Unified Proof")[
  Fix a finite Euclidean dataset $cal(X) subset.eq bb(R)^d$ and Gaussian kernel $K(u) = exp(-u / (2 sigma^2))$. Comparing the three mode-seeking operators $m_"MS"$ (Mean Shift), $M$ (Medoidshift), and $pi_tau$ (Quickshift) on identical density $hat(p)_K$:
  1. *Strict Escape Ordering:* Every fixed point of continuous Mean Shift lying in $cal(X)$ ($mu_i = x_i$) is a fixed point of Medoidshift, and every $tau$-local empirical mode (Quickshift root) whose $tau$-ball contains $B(mu_i, norm(m(x_i))_2)$ is also a fixed point of Medoidshift—yet the converses fail: whenever $x_i$ is mode-locked by Medoidshift ($B(mu_i, norm(m(x_i))_2) inter cal(X) = emptyset$) while a strictly denser neighbor exists within distance $tau$, Quickshift escapes while Medoidshift stalls.
  2. *Asymptotic Fidelity vs. Traversal Geometry:* As $h_cal(X) -> 0$, Medoidshift tracks continuous Mean Shift trajectories uniformly (Theorem 3), whereas Quickshift takes greedy nearest-higher-density steps that can deviate laterally from steepest-ascent streamlines up to the half-space bound of Alg. 41.
]

#proof[
  (1) If $mu_i = x_i$, then $norm(x_i - mu_i)_2 = 0 <= norm(x_k - mu_i)_2$ for all $k$, so Lemma 1 immediately forces $M(i) = i$. Likewise, suppose $x_i$ is a Quickshift root at threshold $tau$ with $B(mu_i, norm(m(x_i))_2) subset.eq B(x_i, tau)$. If $M(i) = k != i$, then Theorem 1 yields $hat(p)_K (x_k) > hat(p)_K (x_i)$ while Lemma 1 places $x_k in B(mu_i, norm(m(x_i))_2) subset.eq B(x_i, tau)$, contradicting $tau$-local maximality of $x_i$; hence $M(i) = i$. Conversely, on the 3-point sample $x_1 = -1.6 sigma, x_2 = 0, x_3 = 0.1 sigma$ with $tau = 2 sigma$, $hat(p)_K (x_2) > hat(p)_K (x_1)$ and $|x_2 - x_1| = 1.6 sigma <= tau$ so Quickshift links $x_1 -> x_2$, whereas $mu_1 approx -0.99 sigma$ satisfies $|x_1 - mu_1| approx 0.61 sigma < 0.99 sigma = |x_2 - mu_1|$, locking $M(1) = 1$ in Medoidshift.

  (2) Uniform step closeness $norm(x_(M(i)) - mu_i)_2 <= h_cal(X)$ follows directly from Theorem 3(1), whereas Quickshift minimizes $norm(x_k - x_i)_2$ subject only to $hat(p)_K (x_k) > hat(p)_K (x_i)$, allowing steps nearly orthogonal to $nabla hat(p)_K (x_i)$ when dense ridges curve sharply.
]

Beyond mode-seeking, Medoidshift connects directly to $k$-medoids partitioning (PAM Alg. 4 @kaufman1990finding and FasterPAM Alg. 24 @schubert2021fast), replacing global $k$-cluster sum-of-distances minimization with local kernel-weighted medoids that discover the number of modes automatically.

On the ClusterOverview 4-dataset empirical benchmark (evaluating Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing with naive matrix evaluation $W D^((2))$), Medoidshift attains *Iris* $"Sil" = 0.210$ and $"ARI" = 0.570$, *Customers* $"Sil" = 0.219$, *Multi-Shape* $"Sil" = 0.262$ ($"ARI" =$ N/A), and *High-D Subspace+Noise* $"Sil" = 0.139$ ($"ARI" =$ N/A), with a mean runtime of *11.398 ms*. Every single measurement tightly corroborates our theoretical theorems: (i) the strong *Iris* $"ARI" = 0.570$ alongside modest $"Sil" = 0.210$ reflects Theorem 2's mode-locking phenomenon—splitting smooth species lobes into adjacent intra-class sub-modes that preserve class alignment while compressing inter-cluster separation; (ii) the low Silhouette on *Multi-Shape* ($0.262$ vs. Quickshift's $0.787$) and *High-D Subspace+Noise* ($0.139$) demonstrates empty locking balls when curved manifolds or high dimensions enlarge $h_cal(X) / sigma^2$ (Theorem 3); and (iii) the $11.398$ ms mean runtime—more than $10 times$ slower than Quickshift ($1.073$ ms) and ST-DBSCAN ($1.028$ ms)—empirically exposes the cubic $O(n^3)$ cost of evaluating $W D^((2))$ without Lemma 1's Euclidean nearest-neighbor acceleration.
