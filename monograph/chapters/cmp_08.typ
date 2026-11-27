#import "../template/ieee.typ": *

#cmp_entry(
  roman: "VIII",
  title: "Intra-Category Comparative Analysis: Grid, Subspace, and Projection (Algs. 133-150)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 133-150",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part VI of this monograph develops eighteen grid-partitioning, axis-parallel subspace, projected-clustering, soft attribute-weighting, and algebraic linear-subspace algorithms (Algs. 133--150) tackling the dual computational and statistical pathologies of full-space Euclidean clustering: $Omega(n^2)$ pairwise scans on large spatial databases and nearest-neighbor distance concentration in high dimensions @beyer1999nearest @kriegel2009subspacesurvey @jain2010data @aggarwal2013data. Whereas prototype partitioning (Part I) and full-space density methods (Part II) evaluate global isotropic distances across all $d$ coordinates, Category VI exploits localized spatial cells, cluster-specific coordinate subsets $S_m subset.eq [d]$, or arbitrarily oriented affine subspaces $cal(A)_m = mu_m + cal(S)_m subset.eq bb(R)^d$.

Despite wide diversity across quadtrees, separable wavelets, Apriori lattices, Shannon multi-information, Monte Carlo hypercubes, entropy-weighted dissimilarities, Poisson--Binomial tests, local PCA, and convex $ell_1$ self-expressiveness, all eighteen algorithms instantiate a single *unifying projection-and-concentration functional* over $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, seeking pairs ${(C_m, cal(P)_m)}_(m=1)^M$ maximizing
$
  cal(J)({(C_m, cal(P)_m)}_(m=1)^M) &= sum_(m=1)^M Phi(C_m; cal(P)_m) \
  &quad - Omega({(C_m, cal(P)_m)}_(m=1)^M),
$
where $cal(P)_m$ encodes the spatial lattice, wavelet subband, contracting cutting plane, axis-parallel projector $Pi_(S_m)$, per-object soft weight $w_i in Delta^(d-1)$, orthogonal eigenbasis $U_m U_m^top$, or self-expressive column span $X c_i$; $Phi$ measures concentration, density contrast, statistical surprisal, or reconstruction fidelity; and $Omega$ enforces capacity control, subspace parsimony, or non-redundancy.

Four structural fault lines organize this comparative chapter into five sections: (I) pure grid and multi-resolution spatial partitioning (Algs. 133, 135, 141, 142); (II) axis-parallel bottom-up subspace monotonicity and redundancy pruning (Algs. 134, 137, 140, 143, 146, 147); (III) projected clustering, soft attribute weighting, and statistical calibration (Algs. 136, 139, 144, 145, 148, 149); (IV) arbitrarily oriented and sparse algebraic subspaces (Algs. 138, 150); and (V) a complete 18-algorithm comparison matrix, dominance lattice, and practitioner's decision guide.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Pure Grid & Multi-Resolution Spatial Partitioning")

Four algorithms in Category VI partition the ambient data domain $Omega = product_(j=1)^d [a_j, b_j] subset.eq bb(R)^d$ into structured spatial cells rather than computing $O(n^2)$ pairwise distances: Statistical Information Grid (STING, Alg. 133) @wang1997sting, WaveCluster (Alg. 135) @sheikholeslami1998wavecluster, OptiGrid (Alg. 141) @hinneburg1999optigrid, and BANG-Clustering (Alg. 142) @schikuta1998bang. While all four achieve sub-quadratic or strictly linear-in-$n$ time in low dimensions and natural immunity to data-ordering permutations, they diverge sharply along three theoretical axes: data-independent vs. data-adaptive cell counts as $d$ grows, signal-processing noise and bridge suppression across scales, and optimal low-density cutting-plane placement.

#subsec("A", "Hierarchical Moment Aggregation (STING) vs. Laminar Block Trees (BANG)")

STING (\#133, @wang1997sting) builds a static $L$-level regular tree $cal(T)$ with branching factor $r = 2^d$ (quadtree in 2D), storing in every cell $c$ the exact sufficient statistics $(n_c, m_c, s_c, "min"_c, "max"_c, cal(D)_c)$ so that spatial-attribute queries traverse $cal(T)$ top-down in $O(K)$ time *independent of $n$*, where $K = r^(L-1)$ is the leaf count. By contrast, BANG-Clustering (\#142, @schikuta1998bang) constructs a data-adaptive *BANG-file* (Balanced And Nested Grid) directory tree by cyclically bisecting only overfull regions (capacity $> b$), allowing strictly nested blocks $B_1, dots, B_q subset.eq B_0$ whose density is normalized by *effective volume* $V_"eff" (B_0) = "vol"(B_0) - sum_(r=1)^q "vol"(B_r)$ and clustered via descending-density dendrogram merges.

#theorem(title: "Moment Conservation, Laminar Effective-Volume Stability & Spatial Sparsity Dichotomy", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq [0, 1]^d$ carry scalar attribute values $y_1, dots, y_n in bb(R)$, and compare an $L$-level STING hierarchy (\#133, $r = 2^d$, $K = 2^(d (L - 1))$ leaves) against a BANG directory (\#142, bucket capacity $b >= 2$, $M_"BANG"$ non-empty blocks).
  1. *(Zero-Error Bottom-Up Aggregation & One-Sided Cantelli Pruning in STING \#133):* For any non-leaf cell $c$ with disjoint children $"ch"(c)$, the parent count, mean, and sample variance satisfy the exact algebraic identities
     $
       n_c &= sum_(u in "ch"(c)) n_u, quad m_c = 1 / n_c sum_(u in "ch"(c)) n_u m_u, \
       s_c^2 &= 1 / n_c sum_(u in "ch"(c)) n_u ( s_u^2 + (m_u - m_c)^2 ).
     $
     For any upper-tail predicate $Q = {y >= t}$ with $t > m_c$ and arbitrary unknown empirical distribution $cal(D)_c = "None"$, the Cantelli bound
     $
       p_c (Q) := (|{x_i in c : y_i >= t}|) / n_c <= (s_c^2) / (s_c^2 + (t - m_c)^2)
     $
     guarantees *zero false pruning* whenever STING retains all cells whose upper bound exceeds selectivity $p_0$.
  2. *(Exact Laminar Volume Floor & Exponential Cell-Count Separation of BANG \#142 vs. STING \#133):* Every BANG block family is *laminar* (for any two blocks $A, B$, either $A inter B = emptyset$, $A subset.eq B$, or $B subset.eq A$). If block $B_0$ encloses $q >= 1$ maximal proper sub-blocks $B_1, dots, B_q$ generated by cyclic halving at depth increments $Delta_r = "dep"(B_r) - "dep"(B_0) >= 1$ that do not exhaust $B_0$, then
     $
       V_"eff" (B_0) = "vol"(B_0) ( 1 - sum_(r=1)^q 2^(-Delta_r) ) >= 2^(-Delta_"max") "vol"(B_0) > 0,
     $
     where $Delta_"max" = max_r Delta_r$. Moreover, BANG creates at most $M_"BANG" = O(n / b)$ directory blocks independent of ambient empty volume, achieving an exponential storage and query compression over STING's static lattice:
     $
       K / M_"BANG" = Omega( (b dot 2^(d (L - 1))) / n ).
     $
  3. *(Dual Structural Pathologies — Boundary-Splitting Erasure vs. Cyclic-Split Starvation):*
     - *(STING \#133):* A single point cluster of mass $n_0$ concentrated within radius $epsilon -> 0^+$ around an interior grid vertex $v in (0, 1)^d$ can split equally across all $2^d$ incident leaf cells, attenuating its observed cell density by $2^(-d)$ and causing *complete erasure* whenever $n_0 / 2^d < tau <= n_0$.
     - *(BANG \#142):* Under cyclic coordinate splitting $j_t = (t mod d) + 1$, any dataset with $n < b dot 2^d$ terminates at maximum split depth $< d$, leaving *every* coordinate $j > floor(log_2 (n / b))$ completely unpartitioned and rendering BANG blind to clusters separated exclusively along late coordinates!
]

#proof[
  *Part 1:* Because $"ch"(c)$ partitions $cal(X)(c)$ into disjoint subsets, $n_c = sum_u n_u$ and $sum_(i in c) y_i = sum_u sum_(i in u) y_i = sum_u n_u m_u$, yielding the exact mean formula. For the variance, expanding $(y_i - m_c)^2 = ((y_i - m_u) + (m_u - m_c))^2$ and noting that $sum_(i in u) (y_i - m_u) = 0$ identically gives
  $
    n_c s_c^2 &= sum_(u in "ch"(c)) sum_(i in u) [ (y_i - m_u)^2 + (m_u - m_c)^2 ] \
    &= sum_(u in "ch"(c)) n_u ( s_u^2 + (m_u - m_c)^2 ).
  $
  For Cantelli's inequality, let $Y$ be drawn uniformly from ${y_i : x_i in c}$ with $bb(E)[Y] = m_c$, $"Var"(Y) = s_c^2$, and set $Z = Y - m_c$, $a = t - m_c > 0$. For any $lambda > 0$, the event ${Z >= a}$ implies $(Z + lambda)^2 >= (a + lambda)^2$, so Markov's inequality yields
  $
    bb(P)(Z >= a) <= (bb(E)[(Z + lambda)^2]) / ((a + lambda)^2) = (s_c^2 + lambda^2) / ((a + lambda)^2).
  $
  Minimizing over $lambda > 0$ attains its unique global minimum at $lambda^* = s_c^2 / a$, giving $p_c (Q) <= s_c^2 / (s_c^2 + a^2)$. Thus pruning $c$ only when $s_c^2 / (s_c^2 + (t - m_c)^2) < p_0$ never discards a cell with $p_c (Q) >= p_0$.

  *Part 2:* Each BANG split bisects a dyadic hyper-rectangle along a coordinate hyperplane into two disjoint dyadic halves of exact volume $2^(-1)$ times the parent. Two dyadic boxes in the same cyclic-bisection tree either are disjoint or satisfy containment, establishing laminarity. Every maximal proper sub-block $B_r subset B_0$ has volume $2^(-Delta_r) "vol"(B_0)$, and maximal sub-blocks are pairwise disjoint, so $V_"eff" (B_0) = "vol"(B_0) (1 - sum_(r=1)^q 2^(-Delta_r))$. Multiplying by $2^(Delta_"max")$ makes $2^(Delta_"max") (1 - sum_(r=1)^q 2^(-Delta_r))$ a strictly positive integer whenever $union.big_(r=1)^q B_r subset.neq B_0$, hence at least $1$, yielding $V_"eff" (B_0) >= 2^(-Delta_"max") "vol"(B_0)$. Because every BANG split is triggered only when a bucket exceeds $b$ points and buddy-merging prevents arbitrary empty buckets @schikuta1998bang, $M_"BANG" = Theta(n / b)$ while STING allocates $K = 2^(d (L - 1))$ full-lattice leaves.

  *Part 3:* Place $n_0 / 2^d$ points at $v + epsilon s$ for each sign vector $s in {-1, +1}^d$ around an interior STING grid vertex $v$. Each of the $2^d$ incident cells receives exact count $n_0 / 2^d < tau$, erasing a cluster whose true ball mass is $n_0 >= tau$. Conversely, in BANG under cyclic splitting, reaching coordinate $j$ at least once requires a binary split path of depth $>= j$, which on a balanced tree requires at least $2^(j - 1)$ non-empty buckets of size $Omega(b)$; for $n < b dot 2^(j - 1)$ on uniform data, no bucket ever splits along coordinate $j$, leaving all intervals along axis $j$ equal to $[0, 1]$.
]

Theorem 1 separates STING and BANG cleanly: STING excels on 2D--3D databases needing repeated zero-scan $O(K)$ range queries, whereas BANG eliminates static empty-cell blowup on skewed sparse data up to moderate $d$ before hitting cyclic-split starvation.

#subsec("B", "Wavelet Noise-Suppression Bounds & Valley Carving (WaveCluster)")

WaveCluster (\#135, @sheikholeslami1998wavecluster) transforms the $d$-dimensional quantized count histogram $F in bb(R)^(xi_1 times dots times xi_d)$ using $r$ levels of Mallat's separable Discrete Wavelet Transform (DWT) with 1D low-pass scaling filter $h = (h_k)$ and high-pass wavelet filter $g = (g_k)$, and extracts connected components of grid cells whose scale-$r$ approximation coefficients $W_r [bold(m)]$ exceed a threshold $tau_r$. Unlike STING's boxcar averaging, symmetric biorthogonal wavelet filters (specifically Cohen--Daubechies--Feauveau hat-shaped splines) actively *suppress isolated outliers*, *boost interior cluster SNR*, and *carve negative valleys* around dense clusters that sever thin noise bridges.

#theorem(title: "Wavelet Outlier SNR Amplification, Bridge-Disconnection Valley Carving & Shift-Variance Limit", origin: "Original Unified Proof")[
  Let $F[bold(m)] = F_"sig" [bold(m)] + F_"noise" [bold(m)]$ be a $d$-dimensional quantized grid density ($bold(m) in bb(Z)^d$), and let $W_r = (cal(L)_h^(plus.o d))^r F$ denote the level-$r$ separable low-pass DWT approximation using a 1D low-pass filter $h$ normalized so that $sum_k h_k = sqrt(2)$ (unit DC gain per dyadic scale, $norm(h)_2 = 1$ for orthonormal filters).
  1. *(Exact Exponential Outlier-Suppression & Statistical SNR Gain):* Suppose a dense cluster forms a constant plateau of height $A > 0$ over a hypercube exceeding the level-$r$ filter support, while contamination consists either of (a) an isolated point outlier of unit mass $delta_(bold(m)_0)$ or (b) i.i.d. zero-mean spatial grid noise of variance $sigma^2$. Then at scale $r >= 1$:
     - *(Deterministic Signal-to-Outlier Ratio):* On the cluster interior, $W_r^"sig" = 2^(r d / 2) A$, while an isolated outlier spreads with peak amplitude bounded by $norm(W_r^"out")_infinity <= norm(h)_infinity^(r d)$. For any orthonormal filter (where $norm(h)_infinity <= 1$, and $norm(h)_infinity = 2^(-1/2)$ for Haar), the Signal-to-Outlier Ratio amplifies exponentially in $r d$:
       $
         "SOR"_r / "SOR"_0 >= ( (sqrt(2)) / (norm(h)_infinity) )^(r d) >= 2^(r d / 2),
       $
       attaining exact factor $2^(r d)$ for the Haar filter!
     - *(Statistical White-Noise SNR):* For any orthonormal filter $h$ ($norm(h)_2 = 1$), the iterated $d$-dimensional low-pass impulse response $h^((r, d))$ preserves unit $ell_2$ norm ($norm(h^((r, d)))_2 = 1$), keeping output noise variance identically $sigma^2$ while amplifying plateau amplitude by $2^(r d / 2)$, yielding an *exact statistical SNR gain*:
       $
         ("SNR"_r) / ("SNR"_0) = 2^(r d / 2) quad ("power gain" = 2^(r d)).
       $
  2. *(Hat-Filter Valley-Carving Bridge Disconnection):* Consider a symmetric biorthogonal low-pass filter $h$ with positive central weights $h_0 > 0, h_(plus.minus 1) >= 0$ and negative outer side-lobes $h_(plus.minus 2) = -beta < 0$ (e.g., the CDF $(2, 2)$ / $(4, 2)$ hat-shaped filters of @sheikholeslami1998wavecluster), normalized to $sum_k h_k = sqrt(2)$. In 1D (or along the normal axis to a thin bridge), suppose two dense plateaus of height $A$ occupy $|m| >= 2$ and are joined by a narrow noise bridge of height $b in (0, A)$ on $|m| <= 1$. Then the filtered coefficient at the bridge center satisfies
     $
       W_1 [0] = (h_0 + 2 h_1) b - 2 beta A,
     $
     which becomes *strictly negative* (guaranteeing complete bridge disconnection at any threshold $tau_1 >= 0$!) whenever
     $
       b / A < (2 beta) / (h_0 + 2 h_1) = (2 beta) / (sqrt(2) + 2 beta).
     $
     By contrast, any nonnegative averaging kernel (such as STING \#133 or Haar WaveCluster) satisfies $W_1 [0] >= (h_0 + 2 h_1) b + 2 h_2 A > 0$ and *monotonically elevates* the bridge density by bleeding mass from the adjacent plateaus!
  3. *(Exact Haar--STING Equivalence & Dyadic Shift-Variance Counterexample):* Level-$r$ Haar WaveCluster (\#135) on a $2^(L-1)$-per-axis grid produces *identically the same connected components* as level-$(L - r)$ uniform-density STING (\#133) with threshold $tau_"STING" = 2^(-r d / 2) tau_r$. Conversely, dyadic decimation breaks translation invariance: shifting a compact cluster of width $2^(r-1)$ and mass $M$ by $2^(r-1)$ cells along every axis splits its Haar energy across $2^d$ coarse coefficients, dropping peak response from $2^(-r d / 2) M$ to $2^(-(r + 2) d / 2) M$ (a $2^(-d)$ attenuation).
]

#proof[
  *Part 1:* Because the $d$-dimensional low-pass operator is the Kronecker tensor product $cal(L)_h^(plus.o d)$ of $d$ 1D convolution--decimation operators $(cal(L)_h u)[m] = sum_k h_k u[2 m - k]$, a constant 1D array $u equiv 1$ maps to $sum_k h_k = sqrt(2)$ per stage; over $r$ stages and $d$ axes, a plateau of height $A$ scales to $(sqrt(2))^(r d) A = 2^(r d / 2) A$. For a unit impulse $delta_0$, a single 1D stage outputs at most $max_k |h_k| = norm(h)_infinity$; inductively across $r$ stages and $d$ axes, $norm(W_r^"out")_infinity <= norm(h)_infinity^(r d)$ (with exact equality $2^(-r d / 2)$ for Haar $h = (2^(-1/2), 2^(-1/2))$). Dividing the plateau response by the outlier peak yields $"SOR"_r / "SOR"_0 >= (sqrt(2) / norm(h)_infinity)^(r d)$ ($= 2^(r d)$ for Haar). For i.i.d. zero-mean grid noise $F_"noise" [bold(m)]$ with variance $sigma^2$, the rows of the orthonormal analysis operator $cal(L)_h^(plus.o d)$ are orthonormal in $ell_2 (bb(Z)^d)$:
  $
    "Var"(W_r^"noise" [bold(m)]) = sigma^2 norm(h^((r, d)))_2^2 = sigma^2 (norm(h)_2^2)^(r d) = sigma^2,
  $
  since Parseval's identity preserves unit $ell_2$ norm under cascaded orthonormal subband filters. Taking the ratio of signal amplitude $2^(r d / 2) A$ to noise standard deviation $sigma$ gives exact amplitude SNR gain $2^(r d / 2)$.

  *Part 2:* Applying the symmetric 5-tap filter $(h_(-2), h_(-1), h_0, h_1, h_2) = (-beta, h_1, h_0, h_1, -beta)$ (with $h_0 + 2 h_1 - 2 beta = sqrt(2)$) to $u[m] = b$ for $|m| <= 1$ and $u[m] = A$ for $|m| >= 2$ gives directly
  $
    W_1 [0] = h_0 b + 2 h_1 b + 2 (-beta) A = (h_0 + 2 h_1) b - 2 beta A.
  $
  Substituting $h_0 + 2 h_1 = sqrt(2) + 2 beta$ shows $W_1 [0] < 0$ if and only if $b / A < 2 beta / (sqrt(2) + 2 beta)$. For any nonnegative filter ($h_k >= 0$ for all $k$), $W_1 [0] = (h_0 + 2 h_1) b + 2 h_2 A >= sqrt(2) b > 0$, strictly raising the bridge whenever $A > b$ and $h_2 > 0$.

  *Part 3:* In Haar DWT, $r$ stages sum disjoint $2^r times dots times 2^r$ blocks of fine cells and multiply by $2^(-r d / 2)$. Since level-$(L - r)$ STING cells are precisely those same disjoint dyadic $2^r$-blocks storing exact sum $n_c$, the thresholding inequalities $W_r [bold(m)] >= tau_r$ and $n_c >= 2^(r d / 2) tau_r$ are algebraically identical! For shift variance, placing mass $M$ inside $[0, 2^(r-1) - 1]^d$ lies wholly inside coarse block $bold(0)$ (coefficient $2^(-r d / 2) M$), whereas translating by $(2^(r-1), dots, 2^(r-1))$ centers it across the dyadic boundary $2^r$, distributing $M / 2^d$ into each of $2^d$ adjacent coarse blocks and multiplying peak amplitude by $2^(-d)$.
]

Theorem 2 shows why WaveCluster (\#135) dominates STING (\#133) for low-$d$ density clustering at matching $O(n + K)$ cost via $2^(r d / 2)$ SNR gain and valley carving—though both collapse under $O(xi^d)$ grid blowup beyond $d > 8$.

#subsec("C", "Contracting Projections & Optimal Low-Density Cutting Planes (OptiGrid)")

OptiGrid (\#141, @hinneburg1999optigrid) breaks the $O(xi^d)$ grid barrier by replacing uniform regular lattices with recursive low-density hyperplanes selected along $q$ *contracting projections* $P_j : bb(R)^d -> bb(R)$ satisfying $norm(P_j x - P_j y)_2 <= norm(x - y)_2$ (i.e., $P_j x = chevron.l v_j, x chevron.r$ with $norm(v_j)_2 <= 1$). On each projection, OptiGrid estimates a 1D kernel density $hat(f)_(v_j, h) (t)$ and inserts cutting planes $H_(v_j, t^*) = {x : chevron.l v_j, x chevron.r = t^*}$ only at deep local minima.

#theorem(title: "Radon No-False-Splitting Guarantee, Exact Valley-Separation Threshold & Oblique SNR Gain", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ arise from a Gaussian mixture $p(x) = sum_(m=1)^M w_m cal(N)(x; mu_m, sigma^2 I_d)$, and let $f_(v, h) (t)$ be the population 1D Gaussian KDE with bandwidth $h > 0$ along unit vector $v in bb(S)^(d-1)$ ($norm(v)_2 = 1$, a contracting projection).
  1. *(Radon Marginalization & Zero Core-Splitting Margin):* Along any unit projection $v$, the 1D projected density is exact:
     $
       f_(v, h) (t) = sum_(m=1)^M w_m cal(N)(t; chevron.l v, mu_m chevron.r, sigma^2 + h^2).
     $
     Consequently, any cutting plane $H_(v, t^*)$ placed where $f_(v, h) (t^*) <= tau_"cut"$ is strictly excluded from the core of every mode $m$ satisfying $w_m > tau_"cut" sqrt(2 pi (sigma^2 + h^2))$ by the explicit orthogonal clearance
     $
       |t^* - chevron.l v, mu_m chevron.r| >= sqrt( 2 (sigma^2 + h^2) ln( w_m / (tau_"cut" sqrt(2 pi (sigma^2 + h^2))) ) ).
     $
  2. *(Exact Cutting-Plane Valley-Detection Threshold & $sqrt(d)$ Oblique SNR Dominance):* For two balanced modes $w_1 = w_2 = 1 / 2$ separated by $Delta = norm(mu_2 - mu_1)_2$ at angle $theta$ to projection $v$ ($chevron.l v, mu_2 - mu_1 chevron.r = Delta cos theta$), $f_(v, h) (t)$ exhibits a strict interior local minimum at the midpoint $t^* = chevron.l v, (mu_1 + mu_2) / 2 chevron.r$ (enabling an optimal separating cut) *if and only if*
     $
       Delta |cos theta| > 2 sqrt(sigma^2 + h^2).
     $
     On a diagonally separated pair $mu_2 - mu_1 = (Delta / sqrt(d)) bold(1)_d$, *every* axis-parallel projection (STING \#133, WaveCluster \#135, BANG \#142, or axis-only OptiGrid) suffers $|cos theta| = d^(-1/2)$, failing to detect any 1D valley unless $Delta > 2 sqrt(d (sigma^2 + h^2))$, whereas OptiGrid with oblique contracting projections ($v approx (mu_2 - mu_1) / Delta$) restores dimension-free separation $Delta > 2 sqrt(sigma^2 + h^2)$ — a strict $sqrt(d)$ SNR advantage!
  3. *(Polynomial $O(n d log n)$ Leaf Bound vs. Concentric-Shell Obstruction):* Recursive OptiGrid partitioning halting whenever a cell holds $< n_"min"$ points produces at most $floor(n / n_"min")$ populated leaves in $O(q n log n + n d log n)$ time (eliminating STING/WaveCluster's $O(xi^d)$ barrier). However, *no single linear cutting plane* can separate two concentric hyperspheres $norm(x)_2 = r_1 < r_2$ without cutting straight through the outer shell, whereas low-dimensional STING (\#133), WaveCluster (\#135), and BANG (\#142) separate concentric rings effortlessly via local cell adjacency.
]

#proof[
  *Part 1:* Because the convolution of a $d$-dimensional Gaussian $cal(N)(mu_m, sigma^2 I_d)$ with isotropic Gaussian kernel $cal(N)(bold(0), h^2 I_d)$ is $cal(N)(mu_m, (sigma^2 + h^2) I_d)$, projecting onto unit vector $v$ ($norm(v)_2 = 1$) integrates out the $(d-1)$-dimensional hyperplane orthogonal to $v$ and yields $f_(v, h) (t) = sum_m w_m (2 pi (sigma^2 + h^2))^(-1/2) exp(-(t - chevron.l v, mu_m chevron.r)^2 / (2(sigma^2 + h^2)))$. Retaining only the $m$-th non-negative summand gives $tau_"cut" >= f_(v, h) (t^*) >= w_m (2 pi (sigma^2 + h^2))^(-1/2) exp(-(t^* - chevron.l v, mu_m chevron.r)^2 / (2(sigma^2 + h^2)))$. Taking logarithms and solving for $|t^* - chevron.l v, mu_m chevron.r|$ yields the stated exact core-clearance margin.

  *Part 2:* Without loss of generality shift the midpoint to $t^* = 0$ and write $a := (Delta |cos theta|) / 2$ and $s^2 := sigma^2 + h^2$. Then $f_(v, h) (t) = (2 sqrt(2 pi) s)^(-1) [exp(-(t - a)^2 / (2 s^2)) + exp(-(t + a)^2 / (2 s^2))]$. By symmetry, $f'_(v, h) (0) = 0$. Differentiating twice at $t = 0$:
  $
    f''_(v, h) (0) = (1) / (sqrt(2 pi) s^3) exp( - a^2 / (2 s^2) ) ( a^2 / s^2 - 1 ).
  $
  Thus $t^* = 0$ is a strict local minimum ($f''_(v, h) (0) > 0$) if and only if $a > s$, i.e., $Delta |cos theta| > 2 sqrt(sigma^2 + h^2)$! (When $a <= s$, $f_(v, h)$ is globally unimodal with maximum at $0$.) For $mu_2 - mu_1 = (Delta / sqrt(d)) bold(1)_d$, every coordinate vector $e_j$ has $(|chevron.l e_j, mu_2 - mu_1 chevron.r|) / Delta = d^(-1/2)$, inflating the critical separation threshold by exact factor $sqrt(d)$, whereas an aligned oblique unit vector achieves $|cos theta| = 1$.

  *Part 3:* Every recursive cut partitions populated cells into disjoint subsets; discarding cells below $n_"min"$ ensures total populated leaves cannot exceed $floor(n / n_"min")$. Evaluating 1D KDEs on $q$ projections via sorting/grid binning takes $O(q n log n)$ per recursion depth ($O(log n)$ balanced depth, or $O(n d log n)$ total with $q = Theta(d)$ @hinneburg1999optigrid). For concentric spheres $norm(x)_2 = r_1 < r_2$ centered at $bold(0)$, any unit projection $v$ maps the inner sphere onto $[-r_1, r_1]$ and the outer sphere onto $[-r_2, r_2] supset [-r_1, r_1]$ with strictly positive projected density everywhere on $(-r_2, r_2)$, so every linear hyperplane separating any part of the inner sphere necessarily bisects the outer sphere!
]

#wide_figure(
  caption: [*Table I:* Head-to-Head Theoretical & Algorithmic Comparison of Pure Grid and Multi-Resolution Spatial Partitioning Methods (Section I: Algs. 133, 135, 141, 142).],
  table(
    columns: (1.02fr, 0.96fr, 0.92fr, 0.86fr, 0.94fr, 1.30fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Partition Geometry], [Time Complexity], [Effective $d$], [Noise / Bridges], [Primary Strength / Failure Mode],
    [STING (\#133) @wang1997sting], [Static $r$-ary cell tree], [$O(n + K)$ build; $O(K)$ query], [$d <= 4$], [Cantelli / Bleeds], [Zero-scan $O(K)$ queries; $2^d$ split & $r^L$ blowup],
    [WaveCluster (\#135) @sheikholeslami1998wavecluster], [Regular grid + DWT], [$O(n + K)$], [$d <= 6$], [*$2^(r d / 2)$ SNR + Valley*], [Optimal 2D--3D noise/bridge filter; $xi^d$ & shift-var.],
    [OptiGrid (\#141) @hinneburg1999optigrid], [Contracting cut planes], [$O(n d log n)$], [$d <= 50$], [KDE valley cuts], [*No empty-cell blowup*; fails concentric shells],
    [BANG (\#142) @schikuta1998bang], [Laminar nested blocks], [$O(n d log n)$], [$d <= 15$], [Density dendro.], [Adapts to skew ($O(n/b)$ blocks); cyclic starvation],
  ),
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Axis-Parallel Bottom-Up Subspace Clustering & Monotonicity")

When ambient dimensionality reaches $d >= 20$, even adaptive full-space grids (BANG \#142, OptiGrid \#141) deteriorate because irrelevant noise attributes mask cluster valleys across full-dimensional distances @beyer1999nearest. Six algorithms in Category VI resolve this via *bottom-up axis-parallel subspace discovery*, identifying clusters embedded inside coordinate subsets $S subset.eq [d]$: CLIQUE (Alg. 134) @agrawal1998clique, MAFIA (Alg. 140) @goil1999mafia, SUBCLU (Alg. 137) @kailing2004subclu, ENCLUS (Alg. 143) @cheng1999enclus, FIRES (Alg. 146) @kriegel2005fires, and INSCY (Alg. 147) @assent2008inscy. We establish first the unified downward-closure (and interest upward-closure) laws governing Apriori search, and second the exact complexity and recall trade-offs separating exhaustive lattice search from 1D filter-refinement (FIRES) and depth-first redundancy pruning (INSCY).

#subsec("A", "Unified Apriori Monotonicity of Density, Adaptive Grids & Entropy")

Exhaustively scanning all $2^d - 1$ nonempty subspaces of $[d]$ is intractable beyond $d approx 20$. CLIQUE (\#134), MAFIA (\#140), SUBCLU (\#137), and ENCLUS (\#143) prune the subspace lattice by anchoring clusterability to predicates enjoying *downward closure* (anti-monotonicity): whenever a $k$-dimensional subspace $S$ satisfies the criterion, every $(k - 1)$-dimensional projection $T subset S$ must also satisfy it. Yet the four methods instantiate fundamentally different geometric and information-theoretic functionals—with subtle soundness boundaries in MAFIA and an exact dual *upward-closure* law in ENCLUS!

#theorem(title: "Unified Downward Closure of Grid, DBSCAN & Entropy Predicates, Interest Upward Closure & Dimensionality Squeeze", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq [0, 1]^d$, let $pi_S : bb(R)^d -> bb(R)^(|S|)$ denote coordinate projection onto $S subset.eq [d]$, and consider any nested pair $emptyset != T subset.eq S subset.eq [d]$.
  1. *(CLIQUE Grid Selectivity \#134):* Under equal-width per-axis partitions into $xi$ bins, every $S$-unit $u_S = product_(j in S) I_j$ projects onto a $T$-unit $u_T = pi_T (u_S)$ satisfying $"sel"(u_T) >= "sel"(u_S)$. Hence dense $S$-units ($"sel" >= tau$) strictly imply dense $T$-projections.
  2. *(SUBCLU Density-Connectivity \#137):* For every $L_p$ metric ($p in [1, infinity]$), $norm(pi_T (x - y))_p <= norm(pi_S (x - y))_p$, so $N_epsilon^S (x) subset.eq N_epsilon^T (x)$ for all $x in cal(X)$. Consequently: (i) every $S$-core point ($|N_epsilon^S (x)| >= mu$) is a $T$-core point; (ii) every $S$-density-connected cluster $C$ is *wholly contained* in a single maximal $T$-cluster $C' supset.eq C$.
  3. *(ENCLUS Entropy Downward Closure & Exact Multi-Information Upward Closure \#143):* Let $H(S) = -sum_(u in cal(G)_S) p(u) log_2 p(u)$ be the discrete grid Shannon entropy on $S$, and define subspace *interest* (Watanabe total correlation) $I(S) := sum_(j in S) H({j}) - H(S)$. Then:
     - *(Downward-Closed Low Entropy):* $H(T) <= H(S)$, so $H(S) <= omega ==> H(T) <= omega$ with zero false negatives;
     - *(Exact Upward-Closed Interest):* For every $k in [d] without S$,
       $
         I(S union {k}) - I(S) &= H({k}) - H({k} | S) \
         &= I_"MI" ({k}; S) >= 0,
       $
       so $I(T) >= epsilon ==> I(S) >= epsilon$! Hence ENCLUS's joint criterion ($H(S) <= omega$ and $I(S) >= epsilon$) forms a *convex interval* in $2^([d])$: Apriori prunes supersets of $H(S) > omega$ while continuing upward through low-entropy, low-interest subspaces ($H(S) <= omega, I(S) < epsilon$) toward correlated supersets.
  4. *(MAFIA Adaptive-Interval Soundness Boundary \#140):* Merging adjacent equal-density 1D bins into variable-width intervals $I_j$ of normalized width $Delta_j in (0, 1]$ preserves downward closure for raw count thresholds, *but* under MAFIA's volume-normalized density test $n(u_S) >= alpha n "vol"(u_S) = alpha n product_(j in S) Delta_j$, downward closure holds universally *if and only if* the threshold is scaled by the maximum bin width (or equivalently $Delta_j <= 1$ with count floor $alpha n product_(j in S) Delta_j$ monotonically non-increasing in $|S|$). Specifically, if $alpha$ is fixed globally while $Delta_k < 1$ on $k in S without T$, then $u_T$ can *fail* $n(u_T) >= alpha n "vol"(u_T)$ even when $u_S$ passes $n(u_S) >= alpha n "vol"(u_S)$ whenever $n(u_S) / n(u_T) > Delta_k$!
  5. *(The Exponential Fixed-Threshold Dimensionality Squeeze):* In both CLIQUE (\#134) and SUBCLU (\#137), suppose a true cluster in $S^*$ ($|S^*| = k$) follows $cal(N)(mu, sigma^2 I_k)$ with mixing weight $w$. Under fixed CLIQUE cell width $Delta$ or SUBCLU Euclidean radius $epsilon$, the maximum cell/ball probability mass inside $S^*$ decays *exponentially* in $k$:
     $
       p_"CLIQUE" (k) <= w ( "erf"( Delta / (2 sqrt(2) sigma) ) )^k, \
       p_"SUBCLU" (k) <= w dot ( (epsilon / (sqrt(2) sigma))^k ) / (Gamma(k / 2 + 1)).
     $
     Thus any fixed global threshold $tau$ (or $mu / n$) either misses all clusters of dimension $k > O(log(w / tau))$ or collapses to near-zero selectivity where uniform noise floods 1D--2D projections!
]

#proof[
  *Part 1:* Because $x in u_S <==> forall j in S, x_j in I_j$, dropping coordinates $S without T$ relaxes the conjunction, so ${x in cal(X) : x in u_S} subset.eq {x in cal(X) : x in u_T}$ and $"sel"(u_T) = abs(cal(X) inter u_T) / n >= abs(cal(X) inter u_S) / n = "sel"(u_S)$.

  *Part 2:* For any $p in [1, infinity)$, $norm(pi_T (x - y))_p^p = sum_(j in T) |x_j - y_j|^p <= sum_(j in S) |x_j - y_j|^p = norm(pi_S (x - y))_p^p$ (and similarly for $p = infinity$ via $max_(j in T) <= max_(j in S)$). Hence $norm(pi_S (x - y))_p <= epsilon ==> norm(pi_T (x - y))_p <= epsilon$, giving $N_epsilon^S (x) subset.eq N_epsilon^T (x)$. Any $S$-core point therefore has $|N_epsilon^T (x)| >= |N_epsilon^S (x)| >= mu$. Likewise, if $p_0, dots, p_m$ is a direct density-reachability chain in $S$, every intermediate $p_r$ is a $T$-core point and $p_(r+1) in N_epsilon^S (p_r) subset.eq N_epsilon^T (p_r)$, so $C$ stays density-connected inside $T$ and sits wholly inside one maximal $T$-cluster $C'$.

  *Part 3:* Write $S = T union R$ with $T inter R = emptyset$. By the chain rule of discrete Shannon entropy and non-negativity of conditional entropy ($H(R | T) >= 0$),
  $
    H(S) = H(T) + H(R | T) >= H(T),
  $
  proving downward closure of ${S : H(S) <= omega}$. For interest $I(S) = sum_(j in S) H({j}) - H(S)$, adding any $k in.not S$ gives
  $
    I(S union {k}) - I(S) &= ( sum_(j in S) H({j}) + H({k}) - H(S union {k}) ) \
    &quad - ( sum_(j in S) H({j}) - H(S) ) \
    &= H({k}) - (H(S union {k}) - H(S)) \
    &= H({k}) - H({k} | S) = I_"MI" ({k}; S) >= 0,
  $
  because conditioning never increases entropy! Inductively, $T subset.eq S ==> I(T) <= I(S)$, proving exact upward closure of ${S : I(S) >= epsilon}$ and lattice convexity of the intersection ${H(S) <= omega} inter {I(S) >= epsilon}$.

  *Part 4:* Construct an explicit counterexample in $d = 2$ on $[0, 1]^2$ with $n = 100$: let axis $1$ have a single merged wide interval $I_1 = [0, 1]$ ($Delta_1 = 1$), let axis $2$ have a narrow adaptive dense interval $I_2 = [0, 0.1]$ ($Delta_2 = 0.1$) containing $30$ points, and set volume-normalized density threshold $alpha = 1.5$. On $S = {1, 2}$, $"vol"(I_1 times I_2) = 0.1$, so $n(I_1 times I_2) = 30 >= 1.5 times 100 times 0.1 = 15$ (dense!). Yet on $T = {1}$, $"vol"(I_1) = 1.0$, requiring $n(I_1) >= 1.5 times 100 times 1.0 = 150 > 100$ (impossible!). More generally, $n(u_S) >= alpha n "vol"(u_S)$ fails to imply $n(u_T) >= alpha n "vol"(u_T)$ whenever $n(u_S) / n(u_T) > "vol"(u_S) / "vol"(u_T) = Delta_k$. Thus MAFIA's volume-normalized threshold can prune dense higher-dimensional cells unless either (i) candidate generation uses a raw count floor $tau$ or (ii) per-axis adaptive bin widths are bounded below by $Delta_"min"$ and pruned with relaxed lower-dimensional factor $alpha Delta_"min"^(k - |T|)$ @goil1999mafia.

  *Part 5:* For $X tilde.op cal(N)(mu, sigma^2 I_k)$, any axis-parallel cell of side $Delta$ maximizes Gaussian mass when centered at $mu$, giving exact probability $( "erf"(Delta / (2 sqrt(2) sigma)) )^k$. For an Euclidean ball $B_k (mu, epsilon)$, the density is bounded above by $(2 pi sigma^2)^(-k/2)$ everywhere on $B_k (mu, epsilon)$ of volume $pi^(k/2) epsilon^k / Gamma(k/2 + 1)$, yielding exact upper bound $((epsilon / (sqrt(2) sigma))^k) / Gamma(k/2 + 1)$, which decays *super-exponentially* in $k$ by Stirling's formula!
]

Theorem 4 contrasts SUBCLU (grid-free DBSCAN projections) and ENCLUS (multi-information $I(S) >= epsilon$ isolating true coordinate correlation from independent 1D spikes, though plug-in $hat(H)(S)$ incurs downward sparse-grid bias $O((xi^k - 1)/(2 n ln 2))$ when $xi^k >> n$).

#subsec("B", "Exponential Lattice Bounds, FIRES Polynomial Merging & INSCY Redundancy Pruning")

Every sound Apriori bottom-up algorithm (CLIQUE, MAFIA, SUBCLU, ENCLUS) confronts an insurmountable combinatorial wall: a single cluster in a $k$-dimensional subspace $S^*$ casts $2^k - 1$ nonempty dense projections across all subsets $emptyset != T subset.eq S^*$. Two landmark algorithms dismantle this wall from opposite directions: *FIRES (\#146)* (@kriegel2005fires) replaces exponential level-wise lattice search with a polynomial 1D base-cluster similarity graph, while *INSCY (\#147)* (@assent2008inscy) preserves exact density-connected subspace discovery while performing *in-process depth-first redundancy pruning* via a compact SCY-tree index and dimensionality-normalized expected density!

#theorem(title: "Exponential Apriori Lower Bound, FIRES Polynomial Bound vs. 1D Masking, and INSCY Reverse-Inclusion Redundancy Pruning", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq [0, 1]^d$ contain $M >= 1$ disjoint true clusters, each embedded in an axis-parallel subspace of dimension $k >= 2$.
  1. *(Unconditional Exponential Lower Bound for Exhaustive Bottom-Up Search):* Every sound Apriori subspace clustering algorithm (CLIQUE \#134, MAFIA \#140, SUBCLU \#137, ENCLUS \#143) that outputs all dense projections must evaluate at least $2^k - 1$ distinct candidate subspaces per embedded $k$-dimensional cluster, incurring unavoidable worst-case runtime $Omega(n dot 2^k)$ (CLIQUE/MAFIA/ENCLUS) or $Omega(n^2 dot 2^k)$ (SUBCLU) and emitting $Omega(M (2^k - 1))$ redundant output clusters!
  2. *(FIRES \#146 Exact $O(n d^2)$ Inverted-List Filter Bound vs. Deterministic 1D Masking Failure):* In FIRES (@kriegel2005fires), let each dimension $j in [d]$ yield disjoint 1D base clusters ${B_(j, 1), dots, B_(j, m_j)}$ computed by 1D clustering (e.g., 1D DBSCAN or $k$-Means). Computing *all* pairwise base-cluster intersection sizes $|B_(j, a) inter B_(l, b)|$ across all $binom(d, 2)$ dimension pairs via point-indexed inverted lists requires exact worst-case time
     $
       T_"intersect" <= n binom(d, 2) = O(n d^2),
     $
     *completely independent of $k$ and $2^k$*, reducing total FIRES runtime to $O(d n log n + n d^2 + m_"cand" n^2 d_"avg")$. *Conversely*, there exists a clean, well-separated 2-cluster distribution in $[0, 1]^2$ whose 1D coordinate marginals are *identically uniform* ($"Unif"[0, 1]$), causing FIRES to find zero 1D base clusters and suffer $0\%$ recall, whereas SUBCLU (\#137), CLIQUE (\#134), and ENCLUS (\#143) achieve $100\%$ recovery!
  3. *(INSCY \#147 Depth-First Reverse-Inclusion Invariant, Output Compression & Non-Transitivity of $R$-Redundancy):* Define a cluster $C_T$ in $T subset S$ to be *$R$-redundant* ($R in (0, 1]$) with respect to $C_S$ if $|C_S inter C_T| >= R |C_T|$ and $C_S$ exceeds the dimensionality-normalized expected density contrast $F dot E(k)$ (where $E(k) = n dot "vol"_k (epsilon)$).
     - *(Exact Reverse-Inclusion Order):* Lexicographic depth-first post-order traversal of INSCY's SCY-tree visits every strict superset $S supset T$ *strictly before* $T$. Thus every maximal non-redundant cluster $C_S$ is already discovered before any projection $T subset S$ is evaluated, enabling exact in-process pruning of $R$-redundant projections with *output compression factor $2^k - 1$* (reducing $M(2^k - 1)$ clusters to $M$)!
     - *(Exact Dimensionality-Contrast Scaling):* Replacing fixed global $mu$ with dimensionality-normalized threshold $max(mu, F dot n "vol"_k (epsilon))$ matches the exact $k$-dimensional null volume decay $"vol"_k (epsilon) = pi^(k/2) epsilon^k / Gamma(k/2 + 1)$, neutralizing Theorem 4(5)'s fixed-threshold squeeze up to the finite-sample breakdown dimension $k_"crit" = Theta(log n / log(1 / epsilon))$ where $n "vol"_k (epsilon) < 1$.
     - *(Non-Transitivity Counterexample):* For every $R in (1 / 2, 1)$, $R$-redundancy is *not transitive*: there exist clusters $C_1, C_2, C_3$ in $S_3 subset S_2 subset S_1$ with $C_2$ $R$-redundant to $C_1$ and $C_3$ $R$-redundant to $C_2$, yet $|C_1 inter C_3| < R |C_3|$, proving that naive transitive suppression can lose valid projections unless maximal covers are tracked directly.
]

#proof[
  *Part 1:* By Theorem 4(1--3), every nonempty subset $T subset.eq S^*$ of a true dense $k$-dimensional cluster satisfies the downward-closed criterion. Because there are $sum_(j=1)^k binom(k, j) = 2^k - 1$ nonempty subsets of $S^*$, any level-wise Apriori algorithm enumerating all qualifying subspaces must inspect and verify every one of these $2^k - 1$ subsets, costing at least $Omega(n)$ per subspace in grid/entropy scans and $Omega(n^2)$ in unindexed DBSCAN (SUBCLU).

  *Part 2:* For the exact FIRES inverted-list bound, each observation $x_i in cal(X)$ belongs to at most *one* 1D base cluster per coordinate $j in [d]$ (since 1D base clusters partition or sub-partition each axis disjointly). Let $L_i subset.eq {1, dots, d}$ be the coordinates where $x_i$ belongs to a 1D cluster ($|L_i| <= d$). Iterating over each point $x_i$ and incrementing the intersection counter for every unordered pair of base clusters ${(B_(j, a(i,j)), B_(l, a(i,l))) : j < l in L_i}$ executes exactly $binom(|L_i|, 2) <= binom(d, 2)$ constant-time updates per point, totaling at most $n binom(d, 2) = O(n d^2)$ operations across the entire dataset!
  For the deterministic 1D masking counterexample, place $n / 2$ points uniformly on $Q_1 = [0, 1/2]^2$ and $n / 2$ points uniformly on $Q_2 = [1/2, 1]^2$ (density $2$ on $Q_1 union Q_2$ and $0$ on the two off-diagonal quadrants). The 1D marginal density along *both* $x_1$ and $x_2$ is identically
  $
    p_1 (t) = p_2 (t) = 1 / 2 dot 2 bb(I)(t in [0, 1/2]) + 1 / 2 dot 2 bb(I)(t in [1/2, 1]) equiv 1
  $
  for all $t in [0, 1]$! Thus neither coordinate possesses any density peak or valley, so FIRES extracts zero non-trivial 1D base clusters and misses both clusters completely, whereas 2D CLIQUE ($xi >= 4$), 2D SUBCLU, and 2D ENCLUS ($I({1, 2}) = 1$ bit $> 0$) immediately detect $Q_1$ and $Q_2$!

  *Part 3:* Represent each nonempty subspace $S subset.eq [d]$ as a strictly increasing tuple $(j_1 < dots < j_m)$. In INSCY's depth-first SCY-tree expansion, processing $T = (j_1, dots, j_m)$ recursively evaluates all extensions $(j_1, dots, j_m, j')$ ($j' > j_m$) *before* clustering $T$. Formally, for any strict superset $S supset.neq T$, let $j^* = min(S without T)$: if $j^* > max(T)$, $S$ lies in $T$'s subtree and post-order visits $S$ before $T$; if $j^* < max(T)$, $S$ precedes $T$ lexicographically at its first differing index and lies in an already-completed sibling subtree. Thus *every* strict superset $S supset T$ is finalized before $T$, enabling $2^k - 1$ output compression on $M$ disjoint $k$-dimensional clusters.
  Finally, for non-transitivity, fix $R in (1 / 2, 1)$ and choose integers $a, b$ with $a / b >= R > a^2 / b^2$. Taking nested subsets $C_1 subset C_2 subset C_3$ of sizes $|C_1| = a^2$, $|C_2| = a b$, $|C_3| = b^2$ yields $abs(C_1 inter C_2) / abs(C_2) = a / b >= R$ and $abs(C_2 inter C_3) / abs(C_3) = a / b >= R$, yet $abs(C_1 inter C_3) / abs(C_3) = a^2 / b^2 < R$!
]

Theorem 5 orders bottom-up search cleanly: MAFIA (\#140) and ENCLUS (\#143) upgrade CLIQUE (\#134) via adaptive intervals and correlation filtering, INSCY (\#147) strictly dominates SUBCLU (\#137) via in-process redundancy pruning, and FIRES (\#146) trades 1D-masking completeness for polynomial $O(n d^2)$ speed.

#wide_figure(
  caption: [*Table II:* Head-to-Head Theoretical & Algorithmic Comparison of Axis-Parallel Bottom-Up Subspace Clustering Algorithms (Section II: Algs. 134, 137, 140, 143, 146, 147).],
  table(
    columns: (1.00fr, 0.96fr, 0.98fr, 0.78fr, 0.86fr, 1.42fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Monotonicity Law], [Worst-Case Time], [Redundancy], [Dim. Threshold], [Primary Strength / Failure Mode],
    [CLIQUE (\#134) @agrawal1998clique], [Fixed-grid selectivity], [$O(c^k + m k n)$], [$2^k - 1$], [Fixed $tau$ (squeeze)], [DNF interpretable rules; $2^k$ blowup & split loss],
    [MAFIA (\#140) @goil1999mafia], [Adaptive-bin count], [$O(n(d + |cal(U)_"cand"|))$], [$2^k - 1$], [Vol.-norm. (careful)], [Adaptive bins + parallel; non-monotone if naive],
    [SUBCLU (\#137) @kailing2004subclu], [$L_p$ DBSCAN core], [$O(|cal(S)_"cand"| n^2)$], [$2^k - 1$], [Fixed $(epsilon, mu)$], [Arbitrary shapes, grid-free; $O(n^2 2^k)$ & squeeze],
    [ENCLUS (\#143) @cheng1999enclus], [$H$ down / $I$ *up*], [$O(n d dot 2^(k_"max"))$], [Maximal only], [Entropy + interest], [Detects correlation ($I > 0$); sparse-$H$ bias],
    [FIRES (\#146) @kriegel2005fires], [1D containment graph], [$O(d n log n + n d^2 + m n^2 d_"avg")$], [Low (merged)], [1D + refine], [*Breaks $2^k$ barrier*; fails 1D-masked checkerboards],
    [INSCY (\#147) @assent2008inscy], [SCY-tree DFS post], [$O(n d log n + |cal(T)_"vis"| n^2)$], [*Pruned ($R$)*], [*Vol.-normalized* $F$], [In-process non-redundant SUBCLU; $R$ non-trans.],
  ),
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Projected Clustering & Soft Attribute Weighting")

Whereas Section II's bottom-up algorithms allow every point to belong to arbitrarily many overlapping subspace clusters, *projected clustering* and *soft attribute weighting* assign each point to a *single* cluster (or background outlier set) paired with a tailored subspace or soft feature weight vector. Six algorithms in Category VI span this paradigm: PROCLUS (Alg. 136) @aggarwal1999proclus, DOC (Alg. 139) @procopiuc2002doc, FINDIT (Alg. 144) @woo2004findit, COSA (Alg. 145) @friedman2004cosa, P3C (Alg. 148) @moise2006p3c, and STATPC (Alg. 149) @moise2008statpc. We organize their mathematical trade-offs into three head-to-head contrasts: medoid $Z$-score projection vs. bounded Dimension-Oriented Distance, Monte Carlo hypercube sampling vs. entropy-regularized soft weights, and Poisson vs. exact Binomial statistical significance.

#subsec("A", "Medoid Projection (PROCLUS) vs. Bounded Dimension-Oriented Distance (FINDIT)")

PROCLUS (\#136, @aggarwal1999proclus) generalizes $k$-Medoids (PAM, @kaufman1990finding) to projected subspaces: given target cluster count $k$ and average subspace dimensionality $l$, it greedily selects a candidate pool $M$ of medoids, forms full-space Euclidean localities $cal(L)_i = B(m_i, delta_i)$ around each current medoid $m_i$, computes mean coordinate deviation $X_(i, j) = |cal(L)_i|^(-1) sum_(x in cal(L)_i) |x_j - m_(i j)|$, standardizes across dimensions via $Z_(i, j) = (X_(i, j) - Y_i) / sigma_i$ (where $Y_i, sigma_i$ are the mean and standard deviation of ${X_(i, 1), dots, X_(i, d)}$), selects subspaces ${D_i}_(i=1)^k$ minimizing $sum_(i=1)^k sum_(j in D_i) Z_(i, j)$ subject to $|D_i| >= 2$ and $sum_i |D_i| = k l$, and assigns points using the dimension-normalized *segmental Manhattan distance*
$
  d_(D_i) (x, m_i) = (1) / (|D_i|) sum_(j in D_i) |x_j - m_(i j)|.
$
FINDIT (\#144, @woo2004findit) replaces unbounded $L_1$/$L_2$ coordinate distances with the bounded *Dimension-Oriented Distance* (DOD):
$
  "DOD"_epsilon (x, y) := d - sum_(j=1)^d bb(I)( |x_j - y_j| <= epsilon ),
$
counting the exact number of coordinates along which $x$ and $y$ differ by more than $epsilon$, and determines each medoid's relevant coordinates via dimension voting across its DOD neighborhood.

#theorem(title: "PROCLUS Exact Greedy ILP Optimality vs. FINDIT Arbitrary-Noise-Variance Immunity", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ contain $k$ clusters, where cluster $C_i$ lives in relevant subspace $S_i^*$ ($|S_i^*| = r_i$) with within-cluster coordinate spread $sigma_"sig"$, inter-cluster separation $Delta_"sep"$, and $d - r_i$ irrelevant noise coordinates of variance $sigma_"noise"^2$.
  1. *(Exact Global Optimality of PROCLUS's $O(k d log(k d))$ Subspace Selection \#136):* The discrete combinatorial program
     $
       min_(D_1, dots, D_k subset.eq [d]) sum_(i=1)^k sum_(j in D_i) Z_(i, j) \
       "s.t." quad |D_i| >= 2 forall i, quad sum_(i=1)^k |D_i| = k l
     $
     is solved *globally and exactly* by (i) assigning each medoid $i$ its two smallest $Z_(i, j)$ coordinates and (ii) greedily selecting the $k (l - 2)$ smallest remaining entries across all $(i, j)$ pairs in $O(k d log(k d))$ time.
  2. *(PROCLUS Locality Pollution Under High-Variance Noise vs. FINDIT Unbounded-Variance Immunity \#144):*
     - *(PROCLUS Failure):* Because PROCLUS constructs each medoid's locality $cal(L)_i = B_2 (m_i, delta_i)$ in the *full $d$-dimensional Euclidean metric*, the squared distance between two points in the same cluster $C_i$ obeys $bb(E)[norm(x - m_i)_2^2] = 2 r_i sigma_"sig"^2 + 2 (d - r_i) sigma_"noise"^2$, while two points from distinct clusters differ only by $r_i Delta_"sep"^2$ on relevant axes. Whenever $(d - r_i) sigma_"noise"^2 >> r_i Delta_"sep"^2$, full-space Euclidean contrast vanishes ($"CV" -> 0$, @beyer1999nearest), contaminating $cal(L)_i$ with foreign clusters and collapsing the $Z_(i, j)$ separation signal! Moreover, PROCLUS's alternating medoid-swap / subspace-reselection loop is *not* guaranteed to decrease a single global Lyapunov objective monotonically.
     - *(FINDIT Exact Bernstein Contrast Preservation):* Under FINDIT's $"DOD"_epsilon$, let $p_"in" = bb(P)(|X_j - Y_j| <= epsilon | "same", j in S_i^*)$, $p_"out" = bb(P)(|X_j - Y_j| <= epsilon | "diff", j in S_i^*)$, and $q_"noise" = bb(P)(|X_j - Y_j| <= epsilon | j in.not S_i^*)$. Then for *any* noise variance $sigma_"noise"^2 in (0, +infinity)$ (even heavy-tailed Cauchy noise with infinite variance!), the expected DOD separation gap is *completely independent of $sigma_"noise"^2$*:
       $
         bb(E)["DOD"_epsilon (x, y_"diff") - "DOD"_epsilon (x, y_"same")] = r_i (p_"in" - p_"out") > 0,
       $
       and Hoeffding's inequality certifies that the probability of a foreign point beating an in-cluster neighbor satisfies
       $
         bb(P)( "DOD"_epsilon (x, y_"diff") <= "DOD"_epsilon (x, y_"same") ) \
         <= exp( - (r_i^2 (p_"in" - p_"out")^2) / (2 d) ),
       $
       which remains exponentially small whenever $r_i (p_"in" - p_"out") >> sqrt(d)$ regardless of $sigma_"noise"^2$!
  3. *(Semi-Metric Structure & Scale-Sensitivity Trade-Off of DOD):* For all $x, y, z in bb(R)^d$, $"DOD"_epsilon$ is symmetric and satisfies the exact $2 epsilon$-relaxed triangle inequality
     $
       "DOD"_(2 epsilon) (x, z) <= "DOD"_epsilon (x, y) + "DOD"_epsilon (y, z),
     $
     yet fails the standard triangle inequality at fixed $epsilon$ and loses all discrimination if coordinates are left unscaled ($epsilon << sigma_"sig"$ drives $p_"in" -> 0$; $epsilon >> Delta_"sep"$ drives $p_"out" -> 1$).
]

#proof[
  *Part 1:* Let $x_(i j) in {0, 1}$ indicate $j in D_i$. The integer program minimizes $sum_(i, j) Z_(i, j) x_(i j)$ over $sum_j x_(i j) >= 2$ and $sum_(i, j) x_(i j) = k l$. Since the constraint matrix is the incidence matrix of a partition matroid union with a uniform cardinality cut (totally unimodular!), every extreme point is integral. For any feasible solution $x$, if some medoid $i$ omits one of its two smallest entries $j_1, j_2$ while picking a larger $j_3$, swapping $j_3 -> j_1$ strictly reduces the objective while preserving all row sums. Thus every global optimum includes the two smallest entries of every row; the remaining constraint is purely $sum_((i, j) in R) x_(i j) = k(l - 2)$ over the remaining set $R$, whose exact global minimum selects the $k(l - 2)$ smallest values in $R$, found in $O(k d log(k d))$ via sorting (or $O(k d)$ via quickselect).

  *Part 2:* For PROCLUS, let $x, m_i in C_i$ and $y in C_(i')$ ($i' != i$). Under independent coordinates, $bb(E)[norm(x - m_i)_2^2] = 2 r_i sigma_"sig"^2 + 2 (d - r_i) sigma_"noise"^2$ while $bb(E)[norm(y - m_i)_2^2] = r_i (Delta_"sep"^2 + 2 sigma_"sig"^2) + 2 (d - r_i) sigma_"noise"^2$, so the relative contrast ratio obeys
  $
    (bb(E)[norm(y - m_i)_2^2] - bb(E)[norm(x - m_i)_2^2]) / (bb(E)[norm(x - m_i)_2^2]) \
    = (r_i Delta_"sep"^2) / (2 r_i sigma_"sig"^2 + 2 (d - r_i) sigma_"noise"^2) -> 0
  $
  as $(d - r_i) sigma_"noise"^2 -> infinity$, causing $B_2 (m_i, delta_i)$ to capture points from $C_(i')$ with probability approaching $1 - 1 / k$.
  For FINDIT's DOD, write $I_j (x, y) := bb(I)(|x_j - y_j| > epsilon) in {0, 1}$ so $"DOD"_epsilon (x, y) = sum_(j=1)^d I_j (x, y)$. On every irrelevant noise axis $j in.not S_i^*$, the marginal distribution of $|x_j - y_j|$ is *identical* for in-cluster and out-of-cluster pairs, so $bb(E)[I_j (x, y_"diff") - I_j (x, y_"same")] = (1 - q_"noise") - (1 - q_"noise") = 0$ identically, even when $sigma_"noise"^2 = infinity$! On each relevant axis $j in S_i^*$, $bb(E)[I_j (x, y_"diff") - I_j (x, y_"same")] = (1 - p_"out") - (1 - p_"in") = p_"in" - p_"out" > 0$. Summing over $j = 1, dots, d$ gives exact expectation gap $r_i (p_"in" - p_"out")$. Moreover, the random variable $U = "DOD"_epsilon (x, y_"diff") - "DOD"_epsilon (x, y_"same") = sum_(j=1)^d Z_j$ is a sum of $d$ independent bounded increments $Z_j in {-1, 0, +1}$ (range width $2$); Hoeffding's inequality directly gives $bb(P)(U <= 0) <= exp(-2 (bb(E)[U])^2 / (4 d)) = exp(-r_i^2 (p_"in" - p_"out")^2 / (2 d))$.

  *Part 3:* By the scalar triangle inequality, $|x_j - z_j| > 2 epsilon$ implies $|x_j - y_j| > epsilon$ or $|y_j - z_j| > epsilon$, so $bb(I)(|x_j - z_j| > 2 epsilon) <= bb(I)(|x_j - y_j| > epsilon) + bb(I)(|y_j - z_j| > epsilon)$ coordinate-wise; summing over $j = 1, dots, d$ yields $"DOD"_(2 epsilon) (x, z) <= "DOD"_epsilon (x, y) + "DOD"_epsilon (y, z)$. Setting $x = 0, y = epsilon, z = 2 epsilon$ in $d = 1$ gives $"DOD"_epsilon (x, z) = 1 > 0 + 0 = "DOD"_epsilon (x, y) + "DOD"_epsilon (y, z)$, violating the exact triangle inequality at fixed $epsilon$.
]

Thus PROCLUS (\#136) suits standardized, moderate-noise data via fast $O(k d log(k d))$ subspace allocation, while FINDIT (\#144) dominates under extreme or heavy-tailed noise variance once $epsilon$ is calibrated.

#subsec("B", "Monte Carlo Hypercube Sampling (DOC) vs. Soft Entropy Weighting (COSA)")

DOC (\#139, @procopiuc2002doc) and COSA (\#145, @friedman2004cosa) abandon fixed-$k$ medoid heuristics altogether via two radically different mathematical mechanisms. DOC frames projected clustering as a randomized combinatorial optimization problem: find an $(alpha, beta)$-*projective cluster* $(C, D)$ such that $|C| >= alpha n$, $max_(x, y in C) |x_j - y_j| <= w$ for all $j in D$, and $mu(|C|, |D|) = |C| (1 / beta)^(|D|)$ is maximized (with $beta in (0, 1 / 2)$). COSA, by contrast, avoids hard binary subspace indicators $D_i subset.eq [d]$ altogether: it solves an entropy-regularized variational program assigning each *individual object* $i$ a smooth probability simplex weight vector $w_i in Delta^(d-1)$ over coordinates, producing a learned pairwise dissimilarity matrix $D_"COSA" (i, j)$ usable inside any distance-based clustering engine (hierarchical agglomeration, PAM, or multidimensional scaling).

#theorem(title: "DOC Monte Carlo Recovery Complexity vs. COSA Exact Log-Sum-Exp Dual & Metric Violation", origin: "Original Unified Proof")[
  1. *(Exact Probabilistic Recovery Bound of DOC \#139):* Suppose $cal(X) subset.eq bb(R)^d$ contains an $(alpha, beta)$-projective cluster $(C^*, D^*)$ of size $|C^*| >= alpha n$ and width $w$ on $D^*$, such that on every irrelevant axis $j in.not D^*$ at most fraction $2 beta < 1$ of all observations fall inside any width-$2 w$ strip around $p in C^*$. In a single DOC inner trial sampling a uniform pivot $p in cal(X)$ and an independent uniform discriminating sample $X$ of size $r >= 1$, let $D(p, X) := {j in [d] : forall x in X, |x_j - p_j| <= w}$. Then exact subspace recovery occurs with probability
     $
       bb(P)(p in C^*, X subset.eq C^*, D(p, X) = D^*) \
       >= alpha^(r + 1) ( 1 - (d - |D^*|) (2 beta)^r ).
     $
     Setting $r = ceil(log_(2 beta) (1 / (2 d)))$ guarantees per-trial success probability $>= (1 / 2) alpha^(r + 1)$, so $m = ceil(2 alpha^(-(r + 1)) ln(1 / delta))$ trials recover $(C^*, D^*)$ with probability $>= 1 - delta$ in total time
     $
       T_"DOC" &= O( n d dot alpha^(-(r + 1)) ln(1 / delta) ) \
       &= O( n d^(1 + (ln(1 / alpha)) / (ln(1 / (2 beta)))) ln(1 / delta) ),
     $
     which is *linear in $n$* yet *super-polynomial in $d$* whenever $alpha$ decays with $k$ or $beta -> (1 / 2)^-$!
  2. *(Exact Variational Log-Sum-Exp Dual, Sandwich Bound & Exponential Noise Attenuation of COSA \#145):* For each object $i$ with current $k$-NN set $N_i$ and per-coordinate mean neighborhood dissimilarity $overline(d)_(i k) := |N_i|^(-1) sum_(j in N_i) d_(i j k) >= 0$, the entropy-regularized COSA subproblem
     $
       phi_lambda (i) := min_(w_i in Delta^(d-1)) sum_(k=1)^d ( w_(i k) overline(d)_(i k) + lambda w_(i k) log w_(i k) ) quad (lambda > 0)
     $
     admits the unique closed-form Boltzmann--Gibbs optimum and exact soft-minimum dual value
     $
       w_(i k)^* &= (exp( - overline(d)_(i k) / lambda )) / (sum_(l=1)^d exp( - overline(d)_(i l) / lambda )), \
       phi_lambda (i) &= -lambda log( sum_(k=1)^d exp( - overline(d)_(i k) / lambda ) ),
     $
     obeying the tight two-sided sandwich bound
     $
       min_k overline(d)_(i k) - lambda log d <= phi_lambda (i) <= min_k overline(d)_(i k)
     $
     and suppressing every irrelevant noise coordinate $k$ with excess dispersion $Delta_k = overline(d)_(i k) - min_l overline(d)_(i l) > 0$ *exponentially* as $w_(i k)^* <= exp(-Delta_k / lambda)$.
  3. *(Alternating Lyapunov Convergence vs. Unbounded Triangle-Inequality Violation of COSA \#145):* Alternating exact weight updates $W^((t))$ and exact $k$-NN neighborhood updates $N^((t))$ strictly decreases the global COSA objective at every step and terminates in finitely many iterations; however, the induced symmetrized dissimilarity $D_"COSA" (i, j) = max(D_(w_i) (i, j), D_(w_j) (i, j))$ (or average) violates the triangle inequality by an *arbitrarily large multiplicative factor* as $lambda -> 0^+$!
]

#proof[
  *Part 1:* Because $p$ and the $r$ points of $X$ are drawn independently uniformly from $cal(X)$ and $|C^*| >= alpha n$, the event $E_1 := {p in C^* "and" X subset.eq C^*}$ has probability $(abs(C^*) / n)^(r + 1) >= alpha^(r + 1)$. Conditioned on $E_1$, every $j in D^*$ satisfies $|x_j - p_j| <= w$ for all $x in X$ deterministically, so $D^* subset.eq D(p, X)$ with probability $1$. Conversely, for any irrelevant coordinate $j in.not D^*$, each point $x in X$ independently satisfies $|x_j - p_j| <= w$ with probability at most $2 beta$ @procopiuc2002doc, so $bb(P)(j in D(p, X) | E_1) <= (2 beta)^r$. Applying the union bound across all $d - |D^*|$ irrelevant coordinates yields
  $
    bb(P)(D(p, X) = D^* | E_1) >= 1 - (d - |D^*|) (2 beta)^r.
  $
  Multiplying by $bb(P)(E_1) >= alpha^(r + 1)$ gives the exact single-trial probability. Choosing $r = ceil(log_(2 beta) (1 / (2 d)))$ enforces $(d - |D^*|)(2 beta)^r <= 1 / 2$, so per-trial success is $>= p_0 := (1 / 2) alpha^(r + 1)$, requiring $m = ceil(p_0^(-1) ln(1 / delta))$ trials of $O(n d)$ cost each and exponent $alpha^(-r) approx (2 d)^(ln(1 / alpha) / ln(1 / (2 beta)))$.

  *Part 2:* Form the Lagrangian $cal(L)(w_i, nu) = sum_k (w_(i k) overline(d)_(i k) + lambda w_(i k) log w_(i k)) + nu (sum_k w_(i k) - 1)$. Setting $(partial cal(L)) / (partial w_(i k)) = overline(d)_(i k) + lambda (log w_(i k) + 1) + nu = 0$ gives $w_(i k)^* = exp(-(nu / lambda + 1)) exp(-overline(d)_(i k) / lambda) > 0$. Enforcing $sum_k w_(i k)^* = 1$ uniquely fixes the partition function $Z_i = sum_(l=1)^d exp(-overline(d)_(i l) / lambda)$ and the stated Boltzmann--Gibbs weights (unique global minimizer by strict convexity of $w log w$). Substituting $log w_(i k)^* = -overline(d)_(i k) / lambda - log Z_i$ cancels the linear term identically:
  $
    phi_lambda (i) = sum_(k=1)^d w_(i k)^* ( overline(d)_(i k) - overline(d)_(i k) - lambda log Z_i ) = -lambda log Z_i.
  $
  Let $d_"min" := min_k overline(d)_(i k)$. Since $exp(-d_"min" / lambda) <= Z_i <= d exp(-d_"min" / lambda)$, taking $-lambda log(dot)$ yields $d_"min" - lambda log d <= phi_lambda (i) <= d_"min"$ and $w_(i k)^* <= exp(-overline(d)_(i k) / lambda) / exp(-d_"min" / lambda) = exp(-Delta_k / lambda)$.

  *Part 3:* Holding neighborhoods $N$ fixed, updating $W$ globally minimizes the strictly convex separable objective; holding $W$ fixed, replacing each $N_i$ with the $k$-nearest neighbors under $D_(w_i)$ strictly minimizes $sum_(j in N_i) D_(w_i) (i, j)$ among all size-$k$ sets. Because the number of finite neighborhood graphs is finite and every strict update reduces the bounded objective, finite-step termination is guaranteed. For unbounded triangle-inequality violation, consider three points in $bb(R)^3$ ($L >> 1$): $x_1 = (0, 0, L)$, $x_2 = (0, 0, 0)$, and $x_3 = (L, 0, 0)$, paired with local cluster neighbors fixing per-object weights $w_1 = (1 - 2 eta, eta, eta)$, $w_2 = (eta, 1 - 2 eta, eta)$, and $w_3 = (eta, eta, 1 - 2 eta)$ where $eta -> 0^+$ as $lambda -> 0^+$. Evaluating coordinatewise $L_1$ differences $d_(i j) = |x_i - x_j|$ yields:
  - $d(1, 2) = (0, 0, L)$: $D_(w_1) (1, 2) = D_(w_2) (1, 2) = eta L ==> D_"COSA" (1, 2) = eta L$;
  - $d(2, 3) = (L, 0, 0)$: $D_(w_2) (2, 3) = D_(w_3) (2, 3) = eta L ==> D_"COSA" (2, 3) = eta L$;
  - $d(1, 3) = (L, 0, L)$: $D_(w_1) (1, 3) = D_(w_3) (1, 3) = (1 - eta) L ==> D_"COSA" (1, 3) = (1 - eta) L$.
  Hence the triangle-inequality violation ratio diverges without bound:
  $
    (D_"COSA" (1, 3)) / (D_"COSA" (1, 2) + D_"COSA" (2, 3)) = (1 - eta) / (2 eta) -> +infinity quad (eta -> 0^+),
  $
  holding identically for both the $max$-symmetrized and average-symmetrized COSA dissimilarities.
]

Theorem 7 highlights DOC's PAC guarantees vs. $O(d^c)$ runtime blowup at small $alpha$, contrasted against COSA's smooth object-specific weights at $O(T n^2 d)$ non-metric cost.

#subsec("C", "Statistical Significance & Non-Redundancy (P3C vs. STATPC)")

Both P3C (\#148, @moise2006p3c) and STATPC (\#149, @moise2008statpc) attack the root cause of Theorem 4(5)'s dimensionality squeeze: replacing arbitrary density thresholds ($tau$, $epsilon$, $w$) with *formal statistical hypothesis tests* against a coordinate-wise uniform null distribution on $[0, 1]^d$. P3C marks dense 1D Sturges-histogram bins via iterative chi-square uniformity rejection, intersects them into *cluster cores* $H = (S, B)$ whenever observed support $m(H)$ exceeds the Poisson expectation $lambda(H) = n p(H)$ ($p(H) = product_(j in S) p_j$) at p-value threshold $alpha_"Pois"$, and completes soft assignments via Expectation--Maximization (@dempster1977em). STATPC eliminates P3C's Poisson approximation, 1D-binning reliance, and EM post-processing altogether: it evaluates *exact Binomial upper-tail surprisal*
$
  mu(H) := -log_(10) ( sum_(r = m(H))^n binom(n, r) p(H)^r (1 - p(H))^(n - r) )
$
directly on axis-parallel hyper-rectangles $H = (S_H, B_H)$ and selects a maximal *mutually non-redundant* set of statistically significant regions.

#theorem(title: "Exact Binomial-vs.-Poisson Tail Dominance, Automatic Noise-Dimension Penalty & STATPC Non-Redundant Approximation", origin: "Original Unified Proof")[
  Let $H = (S_H, B_H)$ be an axis-parallel hyper-rectangle in $[0, 1]^d$ with volume $p = p(H) in (0, 1)$, observed point count $m = m(H) > n p$, and empirical fraction $q = m / n in (p, 1)$.
  1. *(Strict Chernoff Tail Dominance of STATPC's Exact Binomial Test \#149 Over P3C's Poisson Test \#148):* Let $M tilde.op "Bin"(n, p)$ and $Z tilde.op "Pois"(n p)$. While Le Cam's bound ensures $d_"TV" (M, Z) <= p$, the Binomial MGF strictly dominates the Poisson MGF in reverse: $bb(E)[e^(t M)] <= bb(E)[e^(t Z)]$ for all $t in bb(R)$, and their exact large-deviation upper-tail Chernoff exponents satisfy the strict positive gap
     $
       &I_"Bin" (q, p) - I_"Pois" (q, p) \
       &quad = n (1 - q) ln( (1 - q) / (1 - p) ) + n (q - p) > 0
     $
     for every $q in (p, 1)$! Consequently, P3C's Poisson p-value is systematically *conservative* (understating significance) in low-dimensional subspaces where $p$ is non-negligible, while tightly matching STATPC in high-dimensional subspaces where $p -> 0^+$.
  2. *(Constructive Non-Monotonicity of Statistical Significance & Automatic Noise-Dimension Rejection):* Neither P3C's Poisson test nor STATPC's Binomial surprisal is downward-closed (a 2D checkerboard cluster can achieve $mu(H) -> infinity$ while every 1D projection has $mu = 0$). *Crucially*, appending a uniform irrelevant noise coordinate $j' in.not S_H$ restricted to any proper sub-interval of width $w in (0, 1)$ reduces expected volume to $w p$ and expected count to $w m$, *strictly decreasing* the Binomial surprisal rate:
     $
       (dif) / (dif w) D_"KL" (w q || w p) = q ( ln x - 1 + 1 / x ) > 0, \
       x := (q (1 - w p)) / (p (1 - w q)) > 1.
     $
     Thus STATPC *automatically penalizes and rejects irrelevant noise dimensions* without requiring a target dimensionality parameter $l$ (unlike PROCLUS \#136)!
  3. *(NP-Hardness of Optimal Non-Redundant Selection & Greedy Surprisal Bound in STATPC \#149):* Selecting a maximum-total-surprisal (or maximum-cardinality) mutually non-redundant subset of significant hyper-rectangles from a candidate pool $cal(H)$ is NP-hard (by reduction from Maximum Weight Independent Set on the pairwise redundancy conflict graph $G_"red" = (cal(H), E_"red")$). Ordering candidates by decreasing surprisal $mu(H)$ and greedily admitting $H$ whenever it is non-redundant with all previously selected regions achieves a guaranteed $1 / (Delta_"red" + 1)$-approximation (where $Delta_"red"$ is the maximum conflict degree in $G_"red"$) and recovers all ground-truth clusters *exactly* on separated hyper-rectangular mixtures!
]

#proof[
  *Part 1:* For every $t in bb(R)$, applying the elementary inequality $1 + u <= e^u$ with $u = p (e^t - 1)$ yields
  $
    bb(E)[e^(t M)] = (1 - p + p e^t)^n <= exp( n p (e^t - 1) ) = bb(E)[e^(t Z)].
  $
  Optimizing the Chernoff upper bound $inf_(t > 0) (e^(-t m) bb(E)[e^(t X)])$ gives the exact large-deviation exponents $I_"Bin" (q, p) = n [q ln(q / p) + (1 - q) ln((1 - q) / (1 - p))]$ for $M$ and $I_"Pois" (q, p) = n [q ln(q / p) - q + p]$ for $Z$. Subtracting yields
  $
    g(q) := (I_"Bin" (q, p) - I_"Pois" (q, p)) / n = (1 - q) ln( (1 - q) / (1 - p) ) + (q - p).
  $
  Observe $g(p) = 0$ and differentiate with respect to $q in (p, 1)$:
  $
    g'(q) = -ln( (1 - q) / (1 - p) ) - 1 + 1 = ln( (1 - p) / (1 - q) ) > 0.
  $
  Hence $g(q) > 0$ strictly for all $q in (p, 1)$, proving that STATPC's exact Binomial tail decays strictly faster than P3C's Poisson approximation whenever $p > 0$.

  *Part 2:* For non-monotonicity, the 2D uniform-marginal checkerboard of Theorem 5(2) ($Q_1 = [0, 1/2]^2, Q_2 = [1/2, 1]^2$) yields $p(Q_1) = 1 / 4, m(Q_1) = n / 2$, giving $I_"Bin" (1/2, 1/4) = n (1/2 ln 2 + 1/2 ln(2/3)) approx 0.1438 n -> infinity$ as $n -> infinity$, whereas *both* 1D projections $[0, 1/2]$ have $m = n / 2 = n p$ with zero surprisal ($I_"Bin" = 0$)! Next, define $f(w) := D_"KL" (w q || w p) = w q ln(q / p) + (1 - w q) ln((1 - w q) / (1 - w p))$ for $w in (0, 1]$ and $0 < p < q < 1$. Differentiating with respect to $w$:
  $
    f'(w) &= q ln(q / p) - q ln( (1 - w q) / (1 - w p) ) \
    &quad + (1 - w q) ( (-q) / (1 - w q) + p / (1 - w p) ) \
    &= q ln( (q (1 - w p)) / (p (1 - w q)) ) - q + (p (1 - w q)) / (1 - w p) \
    &= q ( ln x - 1 + 1 / x ),
  $
  where $x := (q (1 - w p)) / (p (1 - w q)) > 1$ because $q > p ==> q - w p q > p - w p q$. Since $h(x) := ln x - 1 + 1 / x$ satisfies $h(1) = 0$ and $h'(x) = (x - 1) / x^2 > 0$ for all $x > 1$, we have $f'(w) > 0$ strictly on $(0, 1]$. Thus $D_"KL" (w q || w p) < D_"KL" (q || p)$ for every $w in (0, 1)$, proving that shrinking an irrelevant axis strictly diminishes Binomial surprisal!

  *Part 3:* Any finite graph $G = (V, E)$ with positive vertex weights $w_v$ embeds as a redundancy conflict graph $G_"red"$ of hyper-rectangles in $d = |E| + 1$ with surprisals proportional to $w_v$, reducing Maximum Weight Independent Set to optimal non-redundant STATPC selection @moise2008statpc. Conversely, when greedy selection admits $v^* = arg max_(v in U) mu(v)$ and deletes ${v^*} union N(v^*)$, the weight lost from the optimal set $O^*$ is at most $abs(O^* inter ({v^*} union N(v^*))) dot mu(v^*) <= (Delta_"red" + 1) mu(v^*)$. Summing across greedy steps gives $sum_(v in A_"greedy") mu(v) >= (Delta_"red" + 1)^(-1) sum_(v in O^*) mu(v)$. On separated hyper-rectangular mixtures, distinct clusters have disjoint conflict neighborhoods and each true maximal box uniquely maximizes surprisal within its projection family (Part 2), yielding exact recovery!
]

Thus STATPC (\#149) represents the statistical apex of axis-parallel projected clustering: replacing P3C's 1D-masked histograms and conservative Poisson tail with exact Binomial surprisal, built-in noise-dimension rejection, and explicit non-redundancy.

#wide_figure(
  caption: [*Table III:* Head-to-Head Theoretical & Algorithmic Comparison of Projected Clustering and Soft Attribute-Weighting Methods (Section III: Algs. 136, 139, 144, 145, 148, 149).],
  table(
    columns: (1.00fr, 0.96fr, 0.96fr, 0.80fr, 0.86fr, 1.42fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Subspace Mechanism], [Time Complexity], [Needs $k, l$?], [Noise Immunity], [Primary Strength / Failure Mode],
    [PROCLUS (\#136) @aggarwal1999proclus], [$k$-Medoids + $Z$-score], [$O(I n k d + k d log(k d))$], [Both $k, l$], [Low (full $L_2$)], [Fast exact ILP dims; locality pollution & non-descent],
    [DOC (\#139) @procopiuc2002doc], [Monte Carlo strip $w$], [$O(n d dot alpha^(-(r+1)))$], [No ($alpha, beta, w$)], [High (width $w$)], [Provable PAC recovery; $d^c$ blowup for small $alpha$],
    [FINDIT (\#144) @woo2004findit], [Bounded $"DOD"_epsilon$ vote], [$O(n k d + |S||M| d)$], [$k$ (+ $epsilon$)], [*Unbounded $sigma_"noise"^2$*], [Immune to heavy noise var.; scale-sensitive $epsilon$],
    [COSA (\#145) @friedman2004cosa], [Soft $w_i in Delta^(d-1)$], [$O(T n^2 d)$], [No ($lambda, k_"nn"$)], [Exp. $e^(-Delta_k / lambda)$], [Smooth object weights; $O(n^2 d)$ & non-metric],
    [P3C (\#148) @moise2006p3c], [1D $chi^2$ + Pois. + EM], [$O(d(n + log^2 n) + T n K d)$], [No ($alpha_"Pois"$)], [High (stat.)], [Auto-$k$ & fuzzy EM; 1D masking & Pois. gap],
    [STATPC (\#149) @moise2008statpc], [Exact Bin. surprisal], [$O(n_s d^2 n + |cal(H)|^2 n)$], [*No* ($alpha_0, alpha_K$)], [*Auto-rejects noise*], [Exact Binomial + non-redundant; axis-aligned only],
  ),
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Arbitrarily Oriented & Sparse Algebraic Subspaces")

Every algorithm analyzed in Sections II and III—from CLIQUE and INSCY to PROCLUS, COSA, and STATPC—shares one rigid geometric assumption: cluster subspaces must align with the canonical coordinate axes $e_1, dots, e_D$. In computer vision (3D rigid-body motion trajectories, face images under varying illumination), hyperspectral unmixing, and dynamical systems, clusters instead lie near a union of *arbitrarily rotated linear or affine subspaces* $cal(A)_m = mu_m + cal(S)_m subset.eq bb(R)^D$. Two flagship Category VI algorithms dismantle the axis-parallel barrier through radically contrasting mathematical engines: *ORCLUS (Alg. 138)* @aggarwal2000orclus via hierarchical local PCA projection, and *Sparse Subspace Clustering (SSC, Alg. 150)* @elhamifar2013ssc @soltanolkotabi2012geometric via convex $ell_1$ self-expressiveness and spectral clustering.

#subsec("A", "Collapse of Axis-Parallel Subspace Methods Under Orthogonal Rotation")

We first prove formally why *no* axis-parallel algorithm in Algs. 134, 136, 137, 139, 140, 143--149 can reliably recover a generically rotated low-dimensional subspace.

#proposition(title: "Exact Coordinate Delocalization & Axis-Parallel Contrast Collapse Under Rotation", origin: "Original Proof")[
  Let a cluster follow $X = U Z in bb(R)^D$ where $Z tilde.op cal(N)(bold(0), I_r)$ ($1 <= r << D$) and $U in bb(R)^(D times r)$ has orthonormal columns ($U^top U = I_r$) drawn from a balanced Hadamard / orthogonal DFT frame satisfying $|U_(j, m)| = D^(-1/2)$ for all $j in [D], m in [r]$ (or asymptotically for Haar-random $U in O(D)$). Then:
  1. *(Zero 1D Marginal Contrast):* Every coordinate marginal obeys $X_j tilde.op cal(N)(0, r / D)$ *identically* for all $j = 1, dots, D$. Consequently, PROCLUS (\#136) suffers degenerate $Z_(i, j) equiv 0$, COSA (\#145) collapses to uniform weights $w_(i j) equiv 1 / D$, and 1D-peak methods (FIRES \#146, P3C \#148) observe zero inter-coordinate discrimination!
  2. *(Irreducible Axis-Parallel Residual Energy vs. Zero Oriented Residual):* For *any* axis-parallel $r$-coordinate projector $Pi_S$ ($|S| = r$), the captured variance is at most $"Tr"(Pi_S U U^top Pi_S) = r^2 / D$, leaving fraction $1 - r / D -> 1$ of the cluster's orthogonal variance uncaptured (and inflating projected hyper-box volume exponentially relative to intrinsic $r$-volume), whereas the true oriented projector $P = U U^top$ achieves *exact zero orthogonal residual*: $bb(E)[norm((I_D - P) X)_2^2] = 0$.
]

#proof[
  For any coordinate $j in [D]$, $X_j = sum_(m=1)^r U_(j, m) Z_m$ is centered Gaussian with variance
  $
    "Var"(X_j) = sum_(m=1)^r U_(j, m)^2 = r dot (1 / D) = r / D,
  $
  identical across all $j = 1, dots, D$! Thus all $D$ coordinate dispersions $X_(i, j)$ in PROCLUS and $overline(d)_(i j)$ in COSA coincide identically in population, yielding $Z_(i, j) = 0$ and $w_(i j)^* = 1 / D$. For any size-$r$ coordinate subset $S subset.eq [D]$, the expected projected energy is $bb(E)[norm(Pi_S X)_2^2] = sum_(j in S) "Var"(X_j) = r^2 / D$, whereas total energy is $bb(E)[norm(X)_2^2] = "Tr"(U^top U) = r$, leaving orthogonal residual energy $r (1 - r / D) -> r$ as $r / D -> 0$. By contrast, $(I_D - U U^top) U Z = bold(0)$ identically!
]

#subsec("B", [Local PCA Projection (ORCLUS) vs. Exact $ell_1$ Subspace Detection (SSC)])

ORCLUS (\#138, @aggarwal2000orclus) extends PROCLUS to arbitrary orientations by maintaining $k_t$ clusters (annealed hierarchically from $k_0 > k$ down to $k$) and projected dimensionality $l_t$ (annealed from $D$ down to target $l$). At each stage, ORCLUS computes each cluster's sample covariance $hat(Sigma)_m$, extracts the orthonormal basis $E_m in bb(R)^(D times l_t)$ spanning its *smallest*-variance orthogonal directions (equivalently projecting onto its $l_t$-dimensional principal subspace $U_m$), merges the pair of clusters minimizing the projected orthogonal residual energy $E_"proj" (C_a union C_b)$, and reassigns points by projected distance.

Sparse Subspace Clustering (SSC, \#150, @elhamifar2013ssc @soltanolkotabi2012geometric) operates on an entirely different algebraic principle—*self-expressiveness*. Given data matrix $X = [x_1, dots, x_n] in bb(R)^(D times n)$ whose columns lie in a union of $K$ linear subspaces $union.big_(k=1)^K cal(S)_k$ (with $d_k = dim(cal(S)_k) << D$), noiseless SSC solves for each column $i in [n]$ the convex Basis Pursuit program
$
  min_(c_i in bb(R)^n) norm(c_i)_1 quad "s.t." quad x_i = X c_i, quad c_(i i) = 0
$
(appending $bold(1)^top c_i = 1$ for *affine* subspaces $cal(A)_k = mu_k + cal(S)_k$), forms the symmetric affinity matrix $W = |C| + |C|^top$, and recovers the $K$ clusters via Normalized Spectral Clustering (@ng2002spectral @shi2000ncut @vonluxburg2007tutorial). A solution $C^*$ satisfies the *Subspace Detection Property (SDP)* if for every $i$, $c_i^*$ is nontrivial and $c_(i j)^* = 0$ whenever $x_i$ and $x_j$ belong to different subspaces.

#theorem(title: [Exact $ell_1$ Subspace Detection (SSC \#150) vs. Local-PCA Perturbation & Intersecting-Subspace Breakdown (ORCLUS \#138)], origin: "Original Unified Proof")[
  Let $X = [X^((1)), dots, X^((K))] in bb(R)^(D times n)$ have unit-norm columns drawn from $K$ linear subspaces $cal(S)_1, dots, cal(S)_K subset.eq bb(R)^D$ of intrinsic dimensions $d_1, dots, d_K$.
  1. *(Exact SDP for Independent Subspaces — Elhamifar & Vidal @elhamifar2013ssc):* Suppose $cal(S)_1, dots, cal(S)_K$ are *independent*, meaning $dim(cal(S)_1 + dots + cal(S)_K) = sum_(k=1)^K d_k$, and every $X_(-i)^((k))$ spans $cal(S)_k$. Then *every* optimal solution $c_i^*$ of noiseless SSC satisfies exact SDP ($c_(i j)^* = 0$ for all $x_j in.not cal(S)_k$).
  2. *(Geometric Incoherence--Inradius SDP for Intersecting Subspaces — Soltanolkotabi & Candès @soltanolkotabi2012geometric):* Even when subspaces *intersect nontrivially* ($cal(S)_k inter cal(S)_l != {bold(0)}$) or satisfy $sum_(k=1)^K d_k > D$, let $cal(P)_(-i)^((k)) := "conv"(plus.minus X_(-i)^((k))) subset.eq cal(S)_k$ be the symmetrized convex hull of the remaining points in $cal(S)_k$, let $r_i := r(cal(P)_(-i)^((k)))$ be its Euclidean *inradius* inside $cal(S)_k$ (radius of the largest Euclidean ball centered at $bold(0)$ contained in $cal(P)_(-i)^((k))$), and let $P_(cal(S)_k)$ be orthogonal projection onto $cal(S)_k$. Whenever the geometric separation condition
     $
       max_(x_j in.not cal(S)_k) norm(P_(cal(S)_k) x_j)_2 < r(cal(P)_(-i)^((k)))
     $
     holds, *every* optimal solution $c_i^*$ of noiseless SSC satisfies exact SDP! (In particular, if $cos theta_(k, l) < r(cal(P)_(-i)^((k)))$ for all $l != k$, where $theta_(k, l)$ is the minimum principal angle, exact subspace separation is guaranteed despite $mu_1 = dots = mu_K = bold(0)$.)
  3. *(Davis--Kahan Local-PCA Stability vs. Intersecting-Subspace Breakdown of ORCLUS \#138):* Let cluster $C_m$ have true covariance $Sigma_m^*$ with $l$-th eigengap $gamma_m := lambda_l (Sigma_m^*) - lambda_(l+1) (Sigma_m^*) > 0$ and empirical covariance $hat(Sigma)_m$. By the Davis--Kahan $sin Theta$ theorem @davis1970rotation, ORCLUS's local PCA subspace estimate $hat(U)_m$ obeys
     $
       norm(sin Theta(hat(U)_m, U_m^*))_F <= (2 norm(hat(Sigma)_m - Sigma_m^*)_F) / gamma_m.
     $
     *Conversely*, when two balanced subspaces $cal(S)_1, cal(S)_2$ of dimension $d_1 = d_2 = l$ intersect at the origin ($mu_1 = mu_2 = bold(0)$) with principal angles $theta_1 = dots = theta_l = theta in (0, pi / 2)$, any full-space Euclidean Voronoi seed near the origin pools symmetric samples from both subspaces with mixed covariance $Sigma_"mix" = (1 / 2)(U_1 U_1^top + U_2 U_2^top)$, whose top-$l$ eigenvectors lie strictly *bisecting* $cal(S)_1$ and $cal(S)_2$ at angle $theta / 2$ to both! Thus ORCLUS's local-PCA initialization incurs an irreducible $O(sin(theta / 2))$ orientation error and fails to separate intersecting subspaces with shared centroids, whereas SSC achieves $100\%$ exact SDP recovery!
]

#proof[
  *Part 1:* Fix $x_i in cal(S)_k$ and decompose any feasible solution $c_i$ of $x_i = X c_i$ ($c_(i i) = 0$) into its in-subspace part $c_i^((k))$ (supported on $cal(S)_k$) and out-of-subspace parts $c_i^((l))$ ($l != k$). Rearranging $x_i = X c_i$ gives
  $
    underbrace(x_i - X^((k)) c_i^((k)), in cal(S)_k) = sum_(l != k) underbrace(X^((l)) c_i^((l)), in cal(S)_l).
  $
  Because $cal(S)_1, dots, cal(S)_K$ are linearly independent, $cal(S)_k inter (sum_(l != k) cal(S)_l) = {bold(0)}$! Therefore *both* sides must vanish identically: $X^((k)) c_i^((k)) = x_i$ and $sum_(l != k) X^((l)) c_i^((l)) = bold(0)$. Hence the restricted vector $tilde(c)_i$ obtained by zeroing every $c_i^((l))$ ($l != k$) remains feasible ($X tilde(c)_i = X^((k)) c_i^((k)) = x_i$), yet satisfies
  $
    norm(tilde(c)_i)_1 = norm(c_i^((k)))_1 < norm(c_i^((k)))_1 + sum_(l != k) norm(c_i^((l)))_1 = norm(c_i)_1
  $
  strictly whenever any out-of-subspace coefficient is nonzero! Thus every global $ell_1$ optimum satisfies exact SDP.

  *Part 2:* Fix $x_i in cal(S)_k$ and consider the restricted in-subspace Basis Pursuit problem:
  $
    min_(c in bb(R)^(n_k - 1)) norm(c)_1 quad "s.t." quad X_(-i)^((k)) c = x_i.
  $
  By strong Lagrange duality (feasible because $X_(-i)^((k))$ spans $cal(S)_k$), its dual program over $nu in cal(S)_k$ is
  $
    max_(nu in cal(S)_k) chevron.l x_i, nu chevron.r quad "s.t." quad norm((X_(-i)^((k)))^top nu)_infinity <= 1,
  $
  attaining optimal value $chevron.l x_i, nu_i^* chevron.r = norm(c_i^*)_1$ at an optimal dual vector $nu_i^* in cal(S)_k$. By definition of inradius $r_i = r(cal(P)_(-i)^((k)))$, the Euclidean ball $B_(cal(S)_k) (bold(0), r_i)$ is contained in $cal(P)_(-i)^((k)) = "conv"(plus.minus X_(-i)^((k)))$. Since every extreme point $plus.minus x_j$ ($x_j in cal(S)_k, j != i$) satisfies $|chevron.l x_j, nu_i^* chevron.r| <= 1$, convexity gives $|chevron.l z, nu_i^* chevron.r| <= 1$ for every $z in cal(P)_(-i)^((k))$, and in particular for $z = r_i nu_i^* / norm(nu_i^*)_2 in B_(cal(S)_k) (bold(0), r_i)$, yielding the exact dual-norm bound
  $
    r_i norm(nu_i^*)_2 <= 1 ==> norm(nu_i^*)_2 <= 1 / r_i.
  $
  Now view $nu_i^* in cal(S)_k subset.eq bb(R)^D$ as a dual certificate for the full ambient problem. For every out-of-subspace column $x_j in.not cal(S)_k$, since $nu_i^* in cal(S)_k$, Cauchy--Schwarz yields
  $
    |chevron.l x_j, nu_i^* chevron.r| = |chevron.l P_(cal(S)_k) x_j, nu_i^* chevron.r| <= norm(P_(cal(S)_k) x_j)_2 dot norm(nu_i^*)_2 \
    <= (norm(P_(cal(S)_k) x_j)_2) / r_i < 1
  $
  strictly by hypothesis! Now take *any* feasible global representation $x_i = sum_(j != i) c_j x_j$ having $c_j != 0$ for at least one $x_j in.not cal(S)_k$. Taking the inner product with $nu_i^*$ yields
  $
    norm(c_i^*)_1 &= chevron.l x_i, nu_i^* chevron.r \
    &= sum_(j in cal(S)_k, j != i) c_j chevron.l x_j, nu_i^* chevron.r + sum_(j in.not cal(S)_k) c_j chevron.l x_j, nu_i^* chevron.r \
    &< sum_(j in cal(S)_k, j != i) |c_j| dot 1 + sum_(j in.not cal(S)_k) |c_j| dot 1 = norm(c)_1.
  $
  Hence no solution with nonzero out-of-subspace support can be optimal, proving exact SDP for every optimal solution!

  *Part 3:* The Davis--Kahan $sin Theta$ inequality follows immediately from @davis1970rotation applied to symmetric matrices $Sigma_m^*$ and $hat(Sigma)_m$ with $l$-dimensional leading eigenspace gap $gamma_m = lambda_l (Sigma_m^*) - lambda_(l+1) (Sigma_m^*)$. For the intersecting-subspace construction, rotate coordinates so $U_1 = [I_l; bold(0)]$ and $U_2 = [(cos theta) I_l; (sin theta) I_l]$ in $bb(R)^(2 l)$. Then the balanced mixed covariance $Sigma_"mix" = (1 / 2)(U_1 U_1^top + U_2 U_2^top)$ decomposes into $l$ identical $2 times 2$ blocks
  $
    B(theta) = 1 / 2 mat(1 + cos^2 theta, cos theta sin theta; cos theta sin theta, sin^2 theta),
  $
  whose exact eigenvalues are $lambda_+ = (1 + cos theta) / 2 = cos^2 (theta / 2)$ and $lambda_- = (1 - cos theta) / 2 = sin^2 (theta / 2)$, with leading eigenvector $(cos(theta / 2), sin(theta / 2))^top$ strictly bisecting $cal(S)_1$ and $cal(S)_2$ at angle $theta / 2$ to both! Because both subspaces share mean $bold(0)$ and radial density, points within radius $R$ of $bold(0)$ fall identically into the bisecting PCA cluster, destroying ORCLUS's pure-seed recovery!
]

#subsec("C", "LASSO Collapse, Spectral Over-Segmentation & Scalability Inversion (SSC vs. ORCLUS)")

Theorem 9 demonstrates SSC's profound geometric supremacy on intersecting and heterogeneous-dimension subspaces. Yet a rigorous IEEE comparison must equally expose *where SSC fails and ORCLUS dominates*: noisy LASSO regularization collapse, within-subspace graph disconnection despite exact SDP, and cubic-in-$n$ vs. linear-in-$n$ computational scaling!

#proposition(title: "LASSO-SSC Collapse Threshold, Spectral Over-Segmentation Despite Exact SDP & Scalability Inversion", origin: "Original Unified Proof")[
  1. *(Exact LASSO-SSC Trivial-Zero Collapse Threshold \#150):* In noisy LASSO-SSC,
     $
       min_(c_i in bb(R)^n : c_(i i) = 0) (1 / 2) norm(x_i - X c_i)_2^2 + lambda norm(c_i)_1,
     $
     the unique global solution collapses to the trivial zero vector $c_i^* = bold(0)$ (isolating vertex $i$ completely in the affinity graph $W = |C| + |C|^top$) *if and only if*
     $
       lambda >= lambda_"max" (i) := max_(j != i) |chevron.l x_j, x_i chevron.r|.
     $
  2. *(Exact SDP Does NOT Imply Connected Clusters — Spectral Over-Segmentation Counterexample):* There exists a noiseless dataset of $n = 4$ unit-norm points lying inside a *single* 2-dimensional linear subspace $cal(S)_1 = bb(R)^2$ ($K = 1$) for which Basis Pursuit SSC satisfies exact SDP with zero cross-subspace error, *yet* the induced affinity graph $W = |C^*| + |C^*|^top$ has *two completely disconnected components* ($lambda_2 (L_"sym") = 0$), causing Normalized Spectral Clustering (@shi2000ncut @ng2002spectral) to over-segment a single subspace!
  3. *(Scalability Inversion & Dimension-Heterogeneity Trade-Off):* Solving SSC via ADMM plus truncated spectral eigendecomposition costs $O(T min(n, D) n^2 + n^3)$ time and $O(n^2)$ memory (infeasible beyond $n approx 10^4$ without landmark approximation), while handling *arbitrary heterogeneous* dimensions $d_1 != d_2 != dots != d_K$ automatically without specifying $d_k$. Conversely, ORCLUS (\#138) costs $O(k_0^3 + k_0 n D + k_0^2 D^3)$ time and $O(n D + k_0 D^2)$ memory—scaling *linearly in $n$* to $n >= 10^6$ records when $D$ is moderate—at the price of enforcing a *single uniform* target dimension $l$ across all clusters and requiring separated centroids!
]

#proof[
  *Part 1:* The LASSO objective $J_i (c) = (1 / 2) norm(x_i - X_(-i) c)_2^2 + lambda norm(c)_1$ is convex over $c in bb(R)^(n - 1)$. By Fermat's convex optimality condition, $c = bold(0)$ is a global minimizer if and only if $bold(0) in partial J_i (bold(0)) = -X_(-i)^top x_i + lambda [-1, 1]^(n - 1)$, which is algebraically equivalent to $norm(X_(-i)^top x_i)_infinity <= lambda$, i.e., $lambda >= max_(j != i) |chevron.l x_j, x_i chevron.r|$. Moreover, if $c != bold(0)$ satisfies $J_i (c) = J_i (bold(0))$ while $lambda > lambda_"max" (i)$, then
  $
    J_i (c) - J_i (bold(0)) = (1 / 2) norm(X_(-i) c)_2^2 - chevron.l X_(-i)^top x_i, c chevron.r + lambda norm(c)_1 \
    >= (lambda - lambda_"max" (i)) norm(c)_1 > 0,
  $
  confirming strict uniqueness of the trivial zero collapse for $lambda > lambda_"max" (i)$ (and at $lambda = lambda_"max" (i)$ whenever no collinear unit duplicate achieves zero residual).

  *Part 2:* Place four unit vectors in a single 2D subspace $cal(S)_1 = bb(R)^2$: $x_1 = (1, 0)^top$, $x_2 = (-1, 0)^top$, $x_3 = (0, 1)^top$, $x_4 = (0, -1)^top$. To represent $x_1 = c_2 x_2 + c_3 x_3 + c_4 x_4 = (-c_2, c_3 - c_4)^top$, feasibility forces $c_2 = -1$ and $c_3 = c_4$, so $norm(c)_1 = 1 + 2 |c_3|$, uniquely minimized at $c_2^* = -1, c_3^* = c_4^* = 0$! By symmetry, the unique global Basis Pursuit matrix $C^*$ links *only* $(1, 2)$ and $(3, 4)$:
  $
    W = |C^*| + |C^*|^top = mat(0, 2, 0, 0; 2, 0, 0, 0; 0, 0, 0, 2; 0, 0, 2, 0).
  $
  Thus $W$ consists of *two disjoint connected components* ${1, 2}$ and ${3, 4}$ with zero cut between them, even though all four points span one exact 2D subspace! (Extending to $K >= 2$ subspaces and $2 m$ points per orthogonal cone produces identical over-segmentation whenever intra-subspace data exhibit angular clustering!)

  *Part 3:* In ADMM-SSC @elhamifar2013ssc, precomputing and factoring the $n times n$ matrix $X^top X + rho I_n$ (or dual Woodbury $D times D$ update when $D < n$) costs $O(min(n, D) n^2)$ per matrix multiplication across $T$ iterations, followed by $O(n^3)$ full (or $O(K n^2)$ sparse Lanczos) spectral clustering on $W in bb(R)^(n times n)$. In ORCLUS @aggarwal2000orclus, $k_0$ initial seeds require $O(k_0 n D)$ point assignments per pass, $O(k_0 D^3)$ per eigen-update across at most $k_0$ merge rounds ($O(k_0^2 D^3)$), and $O(k_0^3)$ merge-pair priority bookkeeping, staying strictly linear in $n$!
]

#wide_figure(
  caption: [*Table IV:* Head-to-Head Theoretical & Algorithmic Comparison of Oriented and Sparse Algebraic Subspace Algorithms (Section IV: ORCLUS \#138 vs. SSC \#150, benchmarked against axis-parallel PROCLUS \#136 and STATPC \#149).],
  table(
    columns: (1.02fr, 0.94fr, 1.00fr, 0.84fr, 0.88fr, 1.32fr),
    align: (left, left, left, center, center, left),
    [Algorithm], [Subspace Geometry], [Time Complexity], [Heterog. $d_k$], [Intersecting $cal(S)_k$], [Primary Strength / Failure Mode],
    [ORCLUS (\#138) @aggarwal2000orclus], [Arbitrary affine PCA], [$O(k_0^3 + k_0 n D + k_0^2 D^3)$], [No (fixed $l$)], [Fails ($theta/2$ bias)], [*Linear in $n$*, oriented; breaks at shared centroids],
    [SSC (\#150) @elhamifar2013ssc], [Arbitrary lin./affine $ell_1$], [$O(T min(n, D) n^2 + n^3)$], [*Yes* (auto)], [*Exact SDP*], [Provable SDP on intersecting $cal(S)_k$; $O(n^3)$ & over-seg.],
    [PROCLUS (\#136) @aggarwal1999proclus], [Axis-aligned medoids], [$O(I n k D + k D log(k D))$], [Soft (avg $l$)], [Fails rot.], [Fast linear-$n$ baseline; collapses under rotation],
    [STATPC (\#149) @moise2008statpc], [Axis-aligned hyperbox], [$O(n_s D^2 n + |cal(H)|^2 n)$], [*Yes* (auto)], [Fails rot.], [Exact Binomial significance; axis-aligned only],
  ),
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 18-Algorithm Head-to-Head Matrix & Practitioner's Decision Guide")

We now synthesize all eighteen grid, subspace, and projection algorithms of Category VI (Algs. 133--150) in a unified comparative matrix (Table V), followed by exact cross-algorithm algebraic reductions, strict dominance relations, an eight-regime practitioner's selection guide, and open mathematical frontiers linking Category VI to the rest of the monograph.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 18 Grid, Subspace, and Projection Clustering Algorithms (Category VI, Algorithms 133--150). Abbreviations: *GR* = Pure Grid / Multi-Resolution, *BU* = Bottom-Up Axis-Parallel Subspace, *PR* = Projected / Soft Weighting, *OR* = Oriented / Algebraic Subspace; *Axis* = Axis-Parallel Subspaces, *Full* = Full-Space Partition, *Arb.* = Arbitrarily Oriented Subspaces; *Auto-$k$* = Automatic Cluster Count.],
  {
    set text(size: 6.7pt)
    set table(inset: (x: 2.8pt, y: 2.3pt))
    table(
      columns: (1.18fr, 0.36fr, 0.72fr, 1.12fr, 0.62fr, 0.74fr, 1.36fr),
      align: (left, center, center, left, center, center, left),
      [Algorithm], [Grp.], [Subspace], [Time Complexity], [Auto-$k$], [Redundancy / Noise], [Distinctive Strength / Primary Failure Mode],
      [\#133 STING @wang1997sting], [GR], [Full grid], [$O(n + K)$ ($O(K)$ query)], [Yes], [Cantelli / Bleeds], [Zero-scan $O(K)$ queries; $2^d$ split & $r^L$ blowup],
      [\#134 CLIQUE @agrawal1998clique], [BU], [Axis (ovlp.)], [$O(c^k + m k n)$], [Yes], [$2^k - 1$ / Squeeze], [Interpretable DNF covers; $2^k$ blowup & fixed $tau$],
      [\#135 WaveCluster @sheikholeslami1998wavecluster], [GR], [Full DWT], [$O(n + K)$], [Yes], [*$2^(r d/2)$ SNR + Valley*], [Fastest 2D--3D noise/bridge filter; $xi^d$ & shift-var.],
      [\#136 PROCLUS @aggarwal1999proclus], [PR], [Axis (disj.)], [$O(I n k d + k d log(k d))$], [No ($k, l$)], [None / $L_2$ poll.], [Fast exact ILP dims; locality pollution & non-descent],
      [\#137 SUBCLU @kailing2004subclu], [BU], [Axis (ovlp.)], [$O(|cal(S)_"cand"| n^2)$], [Yes], [$2^k - 1$ / Squeeze], [Grid-free arbitrary shapes; $O(n^2 2^k)$ & fixed $(epsilon, mu)$],
      [\#138 ORCLUS @aggarwal2000orclus], [OR], [*Arb. PCA*], [$O(k_0^3 + k_0 n d + k_0^2 d^3)$], [No ($k, l$)], [None / Centroid], [Linear-$n$ oriented subspaces; fails shared centroids],
      [\#139 DOC @procopiuc2002doc], [PR], [Axis (disj.)], [$O(n d dot alpha^(-(r + 1)))$], [Yes], [Greedy / High], [PAC Monte Carlo recovery; $d^c$ blowup at small $alpha$],
      [\#140 MAFIA @goil1999mafia], [BU], [Axis (ovlp.)], [$O(n(d + |cal(U)_"cand"|))$], [Yes], [$2^k - 1$ / Adap.], [Adaptive bins + parallel; naive vol.-norm breaks closure],
      [\#141 OptiGrid @hinneburg1999optigrid], [GR], [Full / Obl.], [$O(n d log n)$], [Yes], [None / KDE cuts], [No empty-cell blowup, $sqrt(d)$ oblique gain; fails shells],
      [\#142 BANG @schikuta1998bang], [GR], [Full nested], [$O(n d log n)$], [Yes], [None / Dendro.], [$O(n/b)$ laminar blocks; cyclic split starvation at high $d$],
      [\#143 ENCLUS @cheng1999enclus], [BU], [Axis (ovlp.)], [$O(n d dot 2^(k_"max"))$], [Yes], [Maximal / $H, I$], [Upward-closed interest $I > 0$; sparse-grid entropy bias],
      [\#144 FINDIT @woo2004findit], [PR], [Axis (disj.)], [$O(n k d + |S||M| d)$], [No ($k, epsilon$)], [None / *$sigma_"noise"^2$-free*], [DOD immune to infinite noise var.; $epsilon$ scale-sensitive],
      [\#145 COSA @friedman2004cosa], [PR], [*Soft $w_i$*], [$O(T n^2 d)$], [Via MDS], [None / $e^(-Delta/lambda)$], [Exact soft-min object weights; $O(n^2 d)$ & non-metric],
      [\#146 FIRES @kriegel2005fires], [BU], [Axis (ovlp.)], [$O(d n log n + n d^2 + m n^2 d_"avg")$], [Yes], [Low / 1D filter], [*Polynomial $O(n d^2)$* breaks $2^k$; fails 1D masking],
      [\#147 INSCY @assent2008inscy], [BU], [Axis (non-red.)], [$O(n d log n + |cal(T)_"vis"| n^2)$], [Yes], [*$R$-pruned* / Norm.], [Exact DFS reverse-inclusion pruning; $R$ non-transitive],
      [\#148 P3C @moise2006p3c], [PR], [Axis (fuzzy)], [$O(d(n + log^2 n) + T n K d)$], [Yes], [Cores / Pois.], [Auto-$k$ Poisson cores + EM; 1D masking & Pois. gap],
      [\#149 STATPC @moise2008statpc], [PR], [Axis (non-red.)], [$O(n_s d^2 n + |cal(H)|^2 n)$], [*Yes*], [*Non-red. / Bin.*], [Exact Binomial surprisal + auto noise rejection; axis-only],
      [\#150 SSC @elhamifar2013ssc], [OR], [*Arb. $ell_1$*], [$O(T min(n, D) n^2 + n^3)$], [Via spec.], [None / *Exact SDP*], [Exact recovery on intersecting $cal(S)_k$; $O(n^3)$ & over-seg.],
    )
  },
)

#subsec("A", "Exact Algebraic Reductions & Strict Dominance Lattice")

Theorems 1--9 and Propositions 1--2 yield an interlocking lattice of *exact limiting reductions* and *strict dominance pairs* across Category VI:
1. *Exact Limiting Reductions:*
   - *WaveCluster (\#135, Haar scale $r$) $equiv$ STING (\#133, level $L - r$)* up to scale factor $2^(r d / 2)$ (Theorem 2.3);
   - *CLIQUE (\#134) $subset$ MAFIA (\#140)* (equal bin widths) and *CLIQUE (\#134) $subset$ ENCLUS (\#143)* (interest $epsilon = 0$);
   - *INSCY (\#147) $->$ SUBCLU (\#137)* as $R -> 1^+, F -> 0^+$, and *FIRES (\#146) $->$ SUBCLU (\#137)* on unmasked 1D-separable clusters;
   - *PROCLUS (\#136) at $l = d$ $equiv$ PAM* (@kaufman1990finding), and axis-restricted *ORCLUS (\#138) $->$ PROCLUS (\#136)*;
   - *COSA (\#145)* bridges isotropic $k$-NN ($lambda -> infinity$), entropy-weighted W-$k$-Means (@huang2005automated), and hard 1D projection ($lambda -> 0^+$, Theorem 7.2);
   - *P3C (\#148)* core search acts as Poisson-calibrated *CLIQUE (\#134)* + EM (@dempster1977em), upgraded by *STATPC (\#149)* to exact Binomial non-redundant boxes (Theorem 8.1);
   - *SSC (\#150)* under connected SDP reduces identically to zero-cut *Normalized Spectral Clustering* (@shi2000ncut @ng2002spectral @vonluxburg2007tutorial).
2. *Strict Practical Dominance Pairs:*
   - *WaveCluster (\#135) dominates STING (\#133)* for 2D--3D density clustering via $2^(r d / 2)$ SNR gain and valley carving (Theorem 2); retain STING strictly for repeated zero-scan $O(K)$ range queries (Theorem 1.1).
   - *OptiGrid (\#141) and BANG (\#142) dominate STING/WaveCluster* for $d >= 8$ (eliminating $O(xi^d)$ cells), with OptiGrid dominating BANG beyond $d >= 15$ via $sqrt(d)$ oblique gain and no cyclic starvation (Theorems 1.3 & 3.2).
   - *INSCY (\#147) dominates SUBCLU (\#137)* via $2^k - 1$ redundancy pruning and volume-normalized density (Theorem 5.3); *MAFIA (\#140)* and *ENCLUS (\#143)* dominate bare *CLIQUE (\#134)*.
   - *STATPC (\#149) dominates P3C (\#148)* and fixed-threshold grids via exact Binomial surprisal and noise-axis rejection (Theorem 8).
   - *SSC (\#150) dominates ORCLUS (\#138)* on intersecting or heterogeneous-$d_k$ subspaces at $n <= 10^4$ (Theorem 9), while *ORCLUS* dominates raw SSC at $n >= 10^5$ on separated uniform-$l$ affine flats (Proposition 2.3).

#subsec("B", "Mathematically Grounded Practitioner's Decision Guide")

Select among Algorithms 133--150 using seven operational regimes directly grounded in Theorems 1--9 and Propositions 1--2:
1. *Low-$d$ Spatial Databases ($n >= 10^5, d <= 4$):* Use *WaveCluster (\#135)* with CDF hat filters for $2^(r d / 2)$ SNR gain and bridge disconnection (Theorem 2); use *STING (\#133)* strictly for repeated zero-scan $O(K)$ range-aggregate queries with Cantelli guarantees (Theorem 1.1).
2. *Moderate-$d$ Full-Space Partitioning ($5 <= d <= 50$):* Deploy *OptiGrid (\#141)* with oblique contracting projections to retain $O(n d log n)$ speed and $sqrt(d)$ oblique SNR gain (Theorem 3.2). Use *BANG (\#142)* on skewed $d <= 12$ data with concentric shells where cutting planes fail (Theorem 3.3), ensuring $n >> b dot 2^d$ (Theorem 1.3).
3. *Overlapping Axis-Parallel Subspaces ($d <= 30, k_"max" <= 8$):* Prefer *INSCY (\#147)* over *SUBCLU (\#137)* for grid-free, $R$-pruned clusters (Theorem 5.3); use *MAFIA (\#140)* (with sound count-floor pruning, Theorem 4.4) or *ENCLUS (\#143)* (monitoring sparse-$H$ bias, Theorem 4.3) over *CLIQUE (\#134)* when rectangular DNF rules or correlation filtering are required.
4. *High-$d$ Axis-Parallel Discovery ($d >= 50$, Breaking $2^k$):* Use *FIRES (\#146)* for $O(n d^2)$ speed when 1D peaks are unmasked (Theorem 5.2). Under 1D masking or fixed-threshold squeeze, default to *STATPC (\#149)* for exact Binomial surprisal and automatic noise rejection (Theorem 8), reserving *P3C (\#148)* for fuzzy EM Gaussian-core posteriors.
5. *Disjoint Projected Partitions & Heavy Noise ($k$ Known or PAC):* Use *PROCLUS (\#136)* on standardized low-noise data ($O(k d log(k d))$ exact dims, Theorem 6.1). Under extreme noise variance $sigma_"noise"^2 >> Delta_"sep"^2$, switch to *FINDIT (\#144)* for variance-free $"DOD"_epsilon$ separation (Theorem 6.2). Use *DOC (\#139)* when unknown $(k, l)$ and PAC guarantees justify $O(n d^c)$ Monte Carlo sampling (Theorem 7.1).
6. *Soft Per-Object Attribute Weighting ($n <= 5000$):* Deploy *COSA (\#145)* with annealed $lambda$ for Boltzmann--Gibbs soft weights $w_i in Delta^(d-1)$ feeding PAM or hierarchical clustering, avoiding triangle-inequality metric indexes (Theorem 7.2--7.3).
7. *Arbitrarily Oriented Linear / Affine Subspaces:* Never use axis-parallel methods (Proposition 1). Deploy *SSC (\#150)* for $n <= 10^4$, intersecting flats, or heterogeneous $d_k$—adding $bold(1)^top C = bold(1)^top$ for affine subspaces, setting $lambda in (lambda_"noise", 0.5 min_i lambda_"max" (i))$ against zero collapse, and checking $lambda_2 (L_"sym")$ against over-segmentation (Theorem 9, Proposition 2.1--2.2). For $n >> 10^4$ separated affine clusters of uniform $l$, deploy *ORCLUS (\#138)* for $O(k_0 n D)$ scaling (Proposition 2.3).

#subsec("C", "Open Theoretical Frontiers & Cross-Paradigm Bridges")

Four open frontiers connect Category VI across this monograph:
1. *Sub-Quadratic Provable SSC:* Breaking SSC's $O(n^3)$ wall via sketched, landmark, or streaming $ell_1$/elastic-net solvers that preserve Soltanolkotabi--Candès SDP guarantees @soltanolkotabi2012geometric while preventing intra-subspace disconnection (Proposition 2.2).
2. *Rotation-Invariant Statistically Calibrated Subspaces:* Extending STATPC's exact Binomial surprisal and noise-axis rejection (Theorem 8) from axis-aligned boxes to oriented low-rank cylinders and affine subspaces, unifying STATPC, OptiGrid, and SSC.
3. *Transitive Submodular Redundancy Pruning:* Replacing INSCY's non-transitive $R$-redundancy (Theorem 5.3) and STATPC's NP-hard conflict graph (Theorem 8.3) with monotone submodular coverage-minus-redundancy objectives admitting $(1 - 1 / e)$ polynomial-time guarantees.
4. *Deep Nonlinear Subspace Self-Expressiveness:* Lifting $ell_1$ self-expressiveness (SSC) and entropy-weighted distances (COSA) into deep autoencoder latent manifolds bridges Category VI directly to Deep Embedded Clustering (DEC @xie2016dec, IDEC @guo2017idec) in Part VIII and Comparison Chapters XII--XIII.
