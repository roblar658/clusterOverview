#import "../template/ieee.typ": *

#cmp_entry(
  roman: "IX",
  title: "Intra-Category Comparative Analysis: Fuzzy, Biclustering, Time-Series, and Domain Methods (Algs. 151-175)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 151-175",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part VII develops twenty-five specialized clustering algorithms (Algs. 151--175) that relax the three rigid assumptions of flat Euclidean partitioning: hard binary assignments, exchangeable $n times d$ tabular vectors, and single-view single-run consensus @jain2010data @aggarwal2013data @xu2005survey. Spanning fuzzy and possibilistic prototypes, two-mode genomic biclustering, elastic and convolutional time-series methods, image superpixels, topographic watersheds, 3D LiDAR clouds, trajectories, resampling and hypergraph ensembles, and multi-view subspace fusion, all twenty-five algorithms instantiate a unified *structured-domain variational template*:
$
  cal(J)(cal(M), Theta) &= sum_(omega in Omega) cal(L)_omega (cal(M), Theta; cal(D)) \
  &quad + Psi_"struct" (cal(M), Theta),
$
where $cal(M)$ encodes membership or structural correspondence (simplex grades, unconstrained typicalities, row--column masks, warping alignments, spatial zones, co-association hypergraphs, or Grassmannian subspaces), $Theta$ collects prototypes or latent factors, $cal(L)_omega$ scores domain-faithful distortion or mutual information, and $Psi_"struct"$ enforces algebraic, topological, temporal, spatial, or cross-view coupling.

Four structural fault lines organize Category VII:
1. *Soft Membership Constraints & Metric Geometry (Algs. 151--156):* Simplex coupling vs. decoupled typicality (FCM, PCM) and spherical, adaptive Mahalanobis (GK), exponential FMLE (GG), RKHS kernel (KFCM), and grid-spatial (sFCM) dissimilarities.
2. *Two-Mode Matrix Coherence & Overlap (Algs. 157--161):* Additive ANOVA residues (Cheng--Church), overlapping additive layers (PLAID), exact binary Galois closures (Bimax), thresholded bipartite singular signatures (ISA), and sparse multiplicative factor analysis (FABIA).
3. *Temporal Transformation Invariance & Smoothness (Algs. 162--165):* Shift-invariant cross-correlation (k-Shape), elastic dynamic time warping (DBA), entropically smoothed warping (Soft-DTW), and linear-time random convolutional embeddings (ROCKET-Clustering).
4. *Spatial Topology, Ensembles & Multi-View Coupling (Algs. 166--175):* Windowed 5D superpixels (SLIC), elevation-invariant immersion (Watershed), 3D $k$-d tree LiDAR clustering (DBSCAN-3D / ECE), MDL sub-trajectory mining (TRACLUS), resampling stability and hypergraph ensembles (Monti, CSPA, HGPA, MCLA), and power-weighted vs. Grassmannian multi-view fusion (MVKM, Co-regularized Spectral).

Sections I--IV prove sharp head-to-head separation theorems across these four subfamilies; Section V synthesizes all twenty-five algorithms in a master reference matrix, exact reduction lattice, and practitioner's decision guide.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Fuzzy & Possibilistic Prototype Families")

Six algorithms generalize $k$-Means @lloyd1982least via continuous memberships $u_(i j) in [0, 1]$: Fuzzy $C$-Means (FCM, \#151) @dunn1974well @bezdek1981fcm @hathaway1986fcm, Possibilistic $C$-Means (PCM, \#152) @krishnapuram1993pcm, Gustafson--Kessel (GK, \#153) @gustafson1979gk, Gath--Geva (GG, \#154) @gath1989gg, Kernel FCM (KFCM, \#155) @zhang2004kfcm, and Spatial FCM (sFCM, \#156) @chuang2006sfcm. All six optimize weighted distortions $sum_(i, j) u_(i j)^m d_j^2 (x_i, v_j)$ ($m > 1$), splitting along membership constraints and metric geometry.

#subsec("A", "Simplex Noise Trap vs. Coincident Collapse (FCM vs. PCM)")

FCM (\#151) minimizes $J_m (U, V) = sum_(i=1)^n sum_(j=1)^c u_(i j)^m norm(x_i - v_j)_A^2$ over the probabilistic fuzzy $c$-partition space $M_"fc" = {U in [0, 1]^(n times c) : sum_(j=1)^c u_(i j) = 1, 0 < sum_i u_(i j) < n}$. Because $sum_j u_(i j) = 1$ forces relative memberships even for extreme outliers, Krishnapuram and Keller @krishnapuram1993pcm proposed PCM (\#152) over typicalities $T = [t_(i j)] in [0, 1]^(n times c)$ with positive scales $bold(eta) in bb(R)_(>0)^c$:
$
  J_(m, bold(eta)) (T, V) &= sum_(j=1)^c sum_(i=1)^n t_(i j)^m d_(i j)^2 \
  &quad + sum_(j=1)^c eta_j sum_(i=1)^n (1 - t_(i j))^m.
$
Theorem 1 proves that dropping the simplex constraint induces an exact theoretical duality: PCM completely cures FCM's $1 / n$ equidistant-noise breakdown, yet destroys inter-cluster repulsion and collapses *every* global optimum onto coincident prototypes!

#theorem(title: "Simplex Noise Trap vs. Redescending Immunity & Coincident Collapse", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ ($c >= 2, m > 1, A = I_d$).
  1. *(FCM \#151 Equidistant Noise Trap & $1 / n$ Breakdown):* Place two symmetric clean clusters of $n / 2$ points at $plus.minus a e_1$ ($a > 0$) and one orthogonal outlier $x_0 = R e_2$ ($R -> infinity$). For any symmetric prototype pair $v_1 = (-v_x, v_y)$, $v_2 = (v_x, v_y)$, FCM forces $u_(0 1) = u_(0 2) = 1 / 2$ for all $R > 0$, pulling both prototypes vertically to
     $
       v_y &= (2^(-m) R) / ((n / 2)(u_"in"^m + (1 - u_"in")^m) + 2^(-m)) \
       &= Theta(R / n) -> +infinity,
     $
     giving finite-sample breakdown point $epsilon^* = 1 / n$.
  2. *(PCM \#152 Redescending M-Estimator Immunity):* For fixed $V$, PCM's unique optimal typicalities decouple across $j$:
     $
       t_(i j)^* (v_j) = (1 + (norm(x_i - v_j)_2^2 / eta_j)^(1 / (m - 1)))^(-1),
     $
     reducing $min_T J_(m, bold(eta)) (T, V)$ identically to $sum_(j=1)^c sum_(i=1)^n rho_(m, eta_j) (norm(x_i - v_j)_2)$ with bounded loss
     $
       rho_(m, eta_j) (r) = eta_j (1 + (eta_j / r^2)^(1 / (m - 1)))^(1 - m) < eta_j.
     $
     An outlier at radius $r = norm(x_0 - v_j)_2$ exerts pull $psi_j (r) := t_(0 j)^m r = O(r^(-(m + 1)/(m - 1))) -> 0$, peaking finitely at $r_j^* = sqrt(eta_j) ((m - 1)/(m + 1))^((m - 1)/2)$ with
     $
       sup_(r >= 0) psi_j (r) &= sqrt(eta_j) ((m + 1) / (2 m))^m \
       &quad dot ((m - 1) / (m + 1))^((m - 1) / 2) < infinity.
     $
  3. *(Exact Global Coincident-Cluster Collapse of PCM \#152):* If $eta_1 = dots = eta_c = eta > 0$ and $cal(V)_eta^* := arg min_v sum_(i=1)^n rho_(m, eta) (norm(x_i - v)_2)$, then *every* global minimizer $(T^*, V^*)$ satisfies $v_j^* in cal(V)_eta^*$ for all $j$. In particular, on two separated clusters of unequal sizes $n_1 > n_2$ at $plus.minus a e_1$ ($4 a^2 > eta$), all $c$ global optimal prototypes collapse identically into the denser cluster, strictly beating the separated configuration $V_"sep" = (-a e_1, a e_1)$ by margin
     $
       J_(m, eta) (T_"sep", V_"sep") - J_(m, eta) (T^*, V^*) \
       quad >= (n_1 - n_2) rho_(m, eta) (2 a) > 0.
     $
]

#proof[
  *Part 1:* By reflection $x_1 |-> -x_1$, $d_(0 1) = d_(0 2) > 0$, so FCM's KKT update $u_(i j) = [sum_ell (d_(i j) / d_(i ell))^(2 / (m - 1))]^(-1)$ enforces $u_(0 1) = u_(0 2) = 1 / 2$ for all $R$. Substituting into the vertical coordinate of $v_j = (sum_i u_(i j)^m x_i) / (sum_i u_(i j)^m)$ (clean points have $x_(i, 2) = 0$ and memberships $u_"in"$ or $1 - u_"in"$) gives $v_y = (2^(-m) R) / ((n / 2)(u_"in"^m + (1 - u_"in")^m) + 2^(-m)) >= (2^(-m) R) / (n / 2 + 2^(-m)) -> +infinity$.

  *Part 2:* Setting $partial J_(m, bold(eta)) / partial t_(i j) = m t_(i j)^(m - 1) d_(i j)^2 - m eta_j (1 - t_(i j))^(m - 1) = 0$ yields $t_(i j)^* = (1 + s)^(-1)$ where $s := (d_(i j)^2 / eta_j)^(1 / (m - 1))$. Substituting $t_(i j)^*$ gives
  $
    rho_(m, eta_j) (r) &= eta_j (s^(m - 1) + s^m) / (1 + s)^m \
    &= eta_j (1 + s^(-1))^(1 - m) < eta_j.
  $
  In PCM's centroid condition $sum_i t_(i j)^m (x_i - v_j) = bold(0)$, writing $u = (r^2 / eta_j)^(1 / (m - 1))$ gives $psi_j (r) = sqrt(eta_j) u^((m - 1) / 2) / (1 + u)^m$. Setting $(log psi_j)' = (m - 1) / (2 u) - m / (1 + u) = 0$ yields the unique maximum $u^* = (m - 1) / (m + 1)$, confirming the stated peak bound and redescending decay $psi_j (r) = O(r^(-(m + 1) / (m - 1))) -> 0$.

  *Part 3:* Minimizing out $T$ yields the completely decoupled functional
  $
    tilde(J)_(m, bold(eta)) (V) := min_T J_(m, bold(eta)) (T, V) = sum_(j=1)^c Phi_(m, eta_j) (v_j),
  $
  with $Phi_(m, eta_j) (v) := sum_(i=1)^n rho_(m, eta_j) (norm(x_i - v)_2)$. With equal $eta_j = eta$, $min_V tilde(J)_(m, eta) (V) = c dot min_v Phi_(m, eta) (v)$, attained iff every $v_j^* in cal(V)_eta^*$. On $n_1$ points at $-a e_1$ and $n_2 < n_1$ at $+a e_1$, $Phi_(m, eta) (-a e_1) = n_2 rho_(m, eta) (2 a) < n_1 rho_(m, eta) (2 a) = Phi_(m, eta) (+a e_1)$; hence for $c = 2$, $tilde(J)_(m, eta) (V^*) <= 2 n_2 rho_(m, eta) (2 a)$, whereas $tilde(J)_(m, eta) (V_"sep") = (n_1 + n_2) rho_(m, eta) (2 a)$, yielding gap at least $(n_1 - n_2) rho_(m, eta) (2 a) > 0$.
]

#remark(title: "FCM Warm-Start & PFCM Hybridization")[
  Theorem 1 explains why PCM (\#152) must be initialized from FCM (\#151) with locally estimated $eta_j$: warm-starting traps PCM inside distinct local basins of $Phi_(m, eta_j)$ without triggering coincident collapse, while Possibilistic-Fuzzy $C$-Means (PFCM) explicitly combines both simplex and typicality terms.
]

#subsec("B", "Metric Geometry: Spherical, Mahalanobis, FMLE, RKHS & Spatial")

Algorithms 153--156 transform FCM's dissimilarity geometry while retaining simplex (or filtered simplex) memberships:
- *Gustafson--Kessel (GK, \#153)* @gustafson1979gk learns per-cluster SPD matrices $A_j succ 0$ under fixed volume $det A_j = rho_j > 0$, yielding $A_j^* = [rho_j det(F_j)]^(1 / d) F_j^(-1)$ from fuzzy scatter
  $
    F_j = (sum_(i=1)^n u_(i j)^m (x_i - v_j)(x_i - v_j)^top) / (sum_(i=1)^n u_(i j)^m).
  $
- *Gath--Geva (GG, \#154)* @gath1989gg adopts the exponential FMLE dissimilarity incorporating prior weights $P_j = n^(-1) sum_i u_(i j)^m$ and unconstrained volume $det F_j$:
  $
    d_"GG"^2 (x_i, v_j) &= (sqrt(det F_j)) / P_j \
    &quad dot exp( 1 / 2 (x_i - v_j)^top F_j^(-1) (x_i - v_j) ).
  $
- *Kernel FCM (KFCM, \#155)* @zhang2004kfcm lifts data into RKHS $cal(H)$ via Mercer kernel $K(x, y) = chevron.l phi(x), phi(y) chevron.r_cal(H)$ @scholkopf1998nonlinear, either placing prototypes in $cal(H)$ (*KFCM-F*, $norm(phi(x_i) - v_j^Phi)_cal(H)^2$) or keeping input prototypes with kernel distance (*KFCM-K*, $d_K^2 = 2(1 - K(x_i, v_j))$).
- *Spatial FCM (sFCM, \#156)* @chuang2006sfcm filters memberships over pixel neighborhoods $cal(N)_i$:
  $
    u'_(i j) = (u_(i j)^p h_(i j)^q) / (sum_(ell=1)^c u_(i ell)^p h_(i ell)^q), \
    h_(i j) = sum_(r in cal(N)_i) u_(r j), quad p, q >= 0.
  $

#theorem(title: "Affine Equivariance, FMLE Singularity, RKHS Separation & Spatial Non-Monotonicity", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ ($c >= 2, m > 1$).
  1. *(Invariance & Volume Hierarchy: FCM \#151 $subset.neq$ GK \#153 $subset.neq$ GG \#154):* Under $tilde(x)_i = M x_i + b$ ($M in "GL"(d, bb(R))$): Euclidean FCM preserves $U$ iff $M = alpha Q$ ($Q in O(d)$); GK is *strictly $"GL"(d, bb(R))$-invariant* in $U$, yet assumes known volume ratios $rho_j$ and fails upon scatter rank drop ($det F_j = 0$); GG is *strictly $"GL"(d, bb(R))$-invariant* while adapting arbitrary volumes and priors, and at $m = 2$ its membership update coincides *identically* with the Gaussian Mixture Model E-step posterior $r_(i j) prop P_j cal(N)(x_i; v_j, F_j)$ @dempster1977em @gath1989gg—yet inherits GMM's unbounded singularity as $det F_j -> 0^+$.
  2. *(Concentric Nonlinear Separation & Bounded Outlier Pull of KFCM \#155):* On two concentric uniform spheres $S_1 = {norm(x)_2 = r_1}$ and $S_2 = {norm(x)_2 = r_2}$ ($0 < r_1 < r_2$), every linear/spherical prototype method (FCM, PCM, GK) admits $v_1 = v_2 = bold(0)$ as a stationary trap and cannot separate $S_1$ from $S_2$. Under quadratic kernel $K(x, y) = (chevron.l x, y chevron.r + 1)^2$, KFCM-F separates $S_1, S_2$ in $cal(H)$ with exact centroid gap $norm(v_2^Phi - v_1^Phi)_cal(H)^2 = (r_2^2 - r_1^2)^2 / d > 0$. Furthermore, Gaussian KFCM-K ($K_sigma (x, y) = exp(-norm(x - y)_2^2 / sigma^2)$) yields *strictly redescending* outlier pull $<= sigma / sqrt(2 e) -> 0$ while preserving $U in M_"fc"$!
  3. *(Impulse Denoising vs. Non-Monotonicity & Thin-Filament Erasure of sFCM \#156):* On a 2D 4-neighbor grid ($|cal(N)_i| = 4, c = 2$):
     - An isolated impulse pixel ($u_(0 1) = epsilon < 1 / 2$) inside a clean cluster-1 patch ($u_(r 1) >= 1 - delta > 1 / 2$) is restored in one step ($u'_(0 1) > 1 / 2$) iff $((1 - delta) / delta)^q > ((1 - epsilon) / epsilon)^p$.
     - Conversely, an endpoint of a clean 1-pixel-wide filament (1 neighbor in cluster 1, 3 in cluster 2, raw membership $1 - delta -> 1$) has $h_(i 1) / h_(i 2) -> 1 / 3$ and is *erased* ($u'_(i 1) < 1 / 2$) whenever $((1 - delta) / delta)^p < 3^q$, eroding the filament inward! Moreover, for $q > 0$, $"sFCM"_(1, q)$ *strictly increases* $J_m (U, V)$ across clean edges.
]

#proof[
  *Part 1:* Under $tilde(x)_i = M x_i + b$, $tilde(F)_j = M F_j M^top$ and $det tilde(F)_j = |det M|^2 det F_j$. Euclidean FCM preserves all distance ratios iff $M^top M = alpha^2 I_d$. In GK, $tilde(A)_j = [rho_j det tilde(F)_j]^(1 / d) tilde(F)_j^(-1) = |det M|^(2 / d) M^(-top) A_j M^(-1)$, so $d_"GK"^2$ scales uniformly by $|det M|^(2 / d)$, canceling identically in $u_(i j)$. In GG, Mahalanobis exponents are $"GL"(d, bb(R))$-invariant and $sqrt(det tilde(F)_j) = |det M| sqrt(det F_j)$ again cancels across $j$. At $m = 2$, letting $g_(i j) := P_j (det F_j)^(-1/2) exp(-(1/2)(x_i - v_j)^top F_j^(-1)(x_i - v_j))$,
  $
    u_(i j) = (d_"GG"^(-2) (x_i, v_j)) / (sum_(ell=1)^c d_"GG"^(-2) (x_i, v_ell)) = g_(i j) / (sum_(ell=1)^c g_(i ell)),
  $
  which matches the GMM E-step responsibility identically (after multiplying numerator and denominator by $(2 pi)^(-d / 2)$). Setting $v_1 = x_1, F_1 = sigma^2 I_d$ ($sigma -> 0^+$) sends $d_"GG"^2 (x_1, v_1) = sigma^d / P_1 -> 0$, collapsing $u_(1 1) -> 1$ onto a single atom.

  *Part 2:* Zero-mean concentric spheres satisfy $bb(E)_(S_a) [x] = bold(0)$ and $bb(E)_(S_a) [x x^top] = (r_a^2 / d) I_d$, so $v_1 = v_2 = bold(0)$ is stationary for FCM/PCM/GK. Under $phi(x) = (x x^top, sqrt(2) x, 1)$, the RKHS means are $v_a^Phi = ((r_a^2 / d) I_d, bold(0), 1)$, giving exact separation
  $
    norm(v_2^Phi - v_1^Phi)_cal(H)^2 = norm((r_2^2 - r_1^2) / d I_d)_F^2 = ((r_2^2 - r_1^2)^2) / d > 0,
  $
  and $norm(phi(x) - v_1^Phi)_cal(H)^2 < norm(phi(x) - v_2^Phi)_cal(H)^2$ on $S_1$ since $r_1^2 < (r_1^2 + r_2^2) / 2 < r_2^2$. For Gaussian KFCM-K, setting $nabla_(v_j) [2 sum_(i, j) u_(i j)^m (1 - K_sigma (x_i, v_j))] = bold(0)$ weights each point by $u_(i j)^m K_sigma (x_i, v_j)$, bounding radial outlier pull by $r exp(-r^2 / sigma^2) <= sigma / sqrt(2 e) -> 0$.

  *Part 3:* For $c = 2$, $u'_(i 1) > 1 / 2 <=> (h_(i 1) / h_(i 2))^q > (u_(i 2) / u_(i 1))^p$. At an isolated impulse pixel, $h_(0 1) / h_(0 2) >= (1 - delta) / delta$, proving one-step restoration when $((1 - delta) / delta)^q > ((1 - epsilon) / epsilon)^p$. At a 1-pixel filament tip, $h_(i 1) = 1 + 2 delta$ and $h_(i 2) = 3 - 2 delta$, so $h_(i 1) / h_(i 2) -> 1 / 3$ and $u'_(i 1) < 1 / 2$ whenever $((1 - delta) / delta)^p < 3^q$. Finally, because FCM's raw $U^*$ uniquely minimizes strictly convex $U |-> J_m (U, V)$ over $M_"fc"$, any boundary perturbation $U' != U^*$ ($q > 0$) strictly raises $J_m (U', V) > J_m (U^*, V)$.
]

#wide_figure(
  caption: [*Table I:* Head-to-Head Theoretical & Algorithmic Comparison of Fuzzy and Possibilistic Prototype Families (Section I: Algs. 151--156).],
  table(
    columns: (1.06fr, 0.90fr, 0.96fr, 0.84fr, 0.96fr, 1.28fr),
    align: (left, left, left, center, left, left),
    [Algorithm], [Membership Space], [Metric Geometry], [Invariance], [Time Complexity], [Distinctive Strength / Failure Mode],
    [FCM (\#151) @bezdek1981fcm], [Simplex $sum_j u_(i j) = 1$], [Fixed SPD $norm(dot)_A^2$], [Similarity], [$O(n c d T)$], [Simple AO descent; $1/n$ noise trap],
    [PCM (\#152) @krishnapuram1993pcm], [Cube $t_(i j) in [0, 1]$], [Fixed SPD + scale $eta_j$], [Similarity], [$O(n c d T)$], [Redescending M-est.; coincident collapse],
    [GK (\#153) @gustafson1979gk], [Simplex $M_"fc"$], [Adaptive Mahalanobis], [*$"GL"(d, bb(R))$*], [$O((n c d^2 + c d^3) T)$], [Ellipsoidal clusters; fixed $rho_j$, rank drop],
    [GG (\#154) @gath1989gg], [Simplex $M_"fc"$], [Exponential FMLE], [*$"GL"(d, bb(R))$*], [$O(T c (n d^2 + d^3))$], [Unequal vol./priors; $det F_j -> 0$ sing.],
    [KFCM (\#155) @zhang2004kfcm], [Simplex $M_"fc"$], [RKHS / Kernel $K(x, y)$], [Kernel-dep.], [$O(T c n^2)$], [Nonlinear + Gauss. robust; $O(n^2)$ mem.],
    [sFCM (\#156) @chuang2006sfcm], [Filtered $u^p h^q$], [Euclidean + grid $cal(N)_i$], [Grid shift], [$O(T c n (d + |cal(N)|))$], [1-step impulse fix; thin-edge erasure],
  ),
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Biclustering Models for Two-Mode Matrices")

Five algorithms in Category VII extract coherent row--column submatrices $(I_k, J_k) subset.eq {1, dots, n} times {1, dots, p}$ from $X in bb(R)^(n times p)$: Cheng--Church (\#157) @cheng2000biclustering, PLAID (\#158) @lazzeroni2002plaid, Bimax (\#159) @prelic2006bimax, Iterative Signature Algorithm (ISA, \#160) @bergmann2003isa, and Factor Analysis for Bicluster Acquisition (FABIA, \#161) @hochreiter2010fabia.

#subsec("A", "Five Structural Signal Models & Exact Coherence Taxonomy")

Let $P_I^perp = I_(|I|) - |I|^(-1) bold(1)_I bold(1)_I^top$ and $P_J^perp = I_(|J|) - |J|^(-1) bold(1)_J bold(1)_J^top$ denote centering projectors on $I, J$. The five methods posit incompatible algebraic structures:
1. *Cheng--Church (\#157)* @cheng2000biclustering seeks maximum-volume *additive ANOVA* blocks $x_(i j) approx mu + alpha_i + beta_j$ via Mean Squared Residue (MSR)
   $
     H(I, J) := 1 / (|I| |J|) norm(P_I^perp X_(I, J) P_J^perp)_F^2 <= delta,
   $
   extracting sequential biclusters by overwriting found entries with random noise (*destructive masking*).
2. *PLAID (\#158)* @lazzeroni2002plaid fits a simultaneous *overlapping additive layer superposition* via continuous relaxation and backfitting:
   $
     X_(i j) &= theta_0 + sum_(k=1)^K (mu_k + alpha_(i k) + beta_(j k)) rho_(i k) kappa_(j k) \
     &quad + epsilon_(i j), quad rho_(i k), kappa_(j k) in {0, 1}.
   $
3. *Bimax (\#159)* @prelic2006bimax binarizes $X$ to $B in {0, 1}^(n times p)$ and enumerates *all* $beta$ inclusion-maximal all-ones submatrices (equivalently Galois closures / maximal bicliques) via divide-and-conquer in $O(n p beta min(n, p))$ time.
4. *ISA (\#160)* @bergmann2003isa finds *multiplicative co-regulated modules* via thresholded bipartite power iteration on row-standardized $E_R$ and column-standardized $E_C$:
   $
     v^((t+1)) = cal(T)_(t_c) (E_R^top u^((t)) / norm(u^((t)))_1), \
     u^((t+1)) = cal(T)_(t_g) (E_C v^((t+1)) / norm(v^((t+1)))_1).
   $
5. *FABIA (\#161)* @hochreiter2010fabia fits a *joint sparse multiplicative factor model* $X = Lambda Z + Upsilon$ ($Upsilon_(i j) tilde.op cal(N)(0, psi_i)$) with independent Laplacian priors on $Lambda, Z$ via variational EM.

#subsec("B", "Formal Separation Across Multiplicative Signals, Overlap & Dropouts")

Theorem 3 proves unbounded structural separation across multiplicative amplitude scaling, overlapping biclusters, and binary bit dropouts.

#theorem(title: "Multiplicative Collapse, Masking Distortion & Exponential Dropout Shattering", origin: "Original Unified Proof")[
  Let $X in bb(R)^(n times p)$.
  1. *(MSR Breakdown of Cheng--Church \#157 on Multiplicative Modules vs. FABIA \#161 / ISA \#160):* On a noiseless rank-1 multiplicative bicluster $X_(I, J) = c lambda z^top$ ($c > 0$) with centered factor variances $sigma_lambda^2 = |I|^(-1) norm(P_I^perp lambda)_2^2 > 0$ and $sigma_z^2 = |J|^(-1) norm(P_J^perp z)_2^2 > 0$,
     $
       H(I, J) = c^2 sigma_lambda^2 sigma_z^2
     $
     (scaling as $c^4$ when both factors scale by $c$). Thus for any tolerance $delta > 0$, raising $c > sqrt(delta) / (sigma_lambda sigma_z)$ causes Cheng--Church to *reject the noiseless module*, whereas FABIA and rank-1 ISA attain *zero reconstruction error* for all $c > 0$.
  2. *(PLAID \#158 Embedding vs. Cheng--Church \#157 Destructive Masking Bias):* Every zero-MSR submatrix embeds losslessly into one PLAID layer. On two noiseless constant overlapping biclusters $B_1 = (I_1, J_1)$ (level $mu_1$) and $B_2 = (I_2, J_2)$ (level $mu_2$) intersecting on proper sub-rectangle $Omega = (I_1 inter I_2) times (J_1 inter J_2)$ with overlap fraction $gamma = abs(Omega) / (abs(I_2) abs(J_2))$, overwriting $B_1$ with independent noise of mean $mu_0$ and variance $sigma_"mask"^2$ inflates expected MSR on $B_2$ to
     $
       bb(E)[H(I_2, J_2)] &= (1 - 1 / abs(I_2))(1 - 1 / abs(J_2)) gamma sigma_"mask"^2 \
       &quad + gamma (1 - gamma_"row")(1 - gamma_"col") (mu_0 - mu_2)^2 > 0,
     $
     where $gamma_"row" = abs(I_1 inter I_2) / abs(I_2) < 1$ and $gamma_"col" = abs(J_1 inter J_2) / abs(J_2) < 1$, fracturing $B_2$, whereas PLAID recovers both layers with zero residual.
  3. *(Exponential Bit-Dropout Shattering of Bimax \#159 vs. Continuous Frobenius Stability):* Corrupting an all-ones $r times s$ bicluster $B$ by $b <= min(r, s)$ zero-flip dropouts on distinct rows and columns causes Bimax (\#159) to shatter the single ground-truth module into *exactly* $2^b$ incomparable inclusion-maximal submatrices ($binom(b, k)$ of size $(r - k) times (s - (b - k))$), whereas continuous models (PLAID, ISA, FABIA) incur only $O(sqrt(b / (r s)))$ Frobenius/subspace perturbation.
]

#proof[
  *Part 1:* By bilinearity, $P_I^perp (c lambda z^top) P_J^perp = c (P_I^perp lambda)(P_J^perp z)^top$, and $norm(u v^top)_F^2 = norm(u)_2^2 norm(v)_2^2$ gives
  $
    H(I, J) &= (c^2 / (|I| |J|)) norm(P_I^perp lambda)_2^2 norm(P_J^perp z)_2^2 \
    &= c^2 sigma_lambda^2 sigma_z^2.
  $
  Hence $H(I, J) > delta$ whenever $c > sqrt(delta) / (sigma_lambda sigma_z)$, whereas $c lambda z^top$ lies identically in FABIA's $Lambda Z$ and ISA's rank-1 singular model.

  *Part 2:* $P_I^perp X_(I, J) P_J^perp = bold(0)$ implies $x_(i j) = x_(I J) + (x_(i J) - x_(I J)) + (x_(I j) - x_(I J))$, matching a PLAID layer. After masking $B_1$, write $M = bb(E)[M] + Xi$ on $B_2$, where $Xi$ is zero-mean variance-$sigma_"mask"^2$ noise supported on $Omega$. By orthogonality,
  $
    bb(E)[norm(P_(I_2)^perp M P_(J_2)^perp)_F^2] &= norm(P_(I_2)^perp bb(E)[M] P_(J_2)^perp)_F^2 \
    &quad + bb(E)[norm(P_(I_2)^perp Xi P_(J_2)^perp)_F^2].
  $
  In $bb(E)[M] = mu_2 bold(1)_(I_2) bold(1)_(J_2)^top + (mu_0 - mu_2) bold(1)_(I_1 inter I_2) bold(1)_(J_1 inter J_2)^top$, the constant vanishes under centering and the indicator product has centered variances $gamma_"row" (1 - gamma_"row")$ and $gamma_"col" (1 - gamma_"col")$; by Part 1 its residue is $(mu_0 - mu_2)^2 gamma (1 - gamma_"row")(1 - gamma_"col")$. For the noise term, since $(P_(I_2)^perp)_(i i) = 1 - 1 / abs(I_2)$ and $(P_(J_2)^perp)_(j j) = 1 - 1 / abs(J_2)$,
  $
    bb(E)[norm(P_(I_2)^perp Xi P_(J_2)^perp)_F^2] &= sum_((i, j) in Omega) sigma_"mask"^2 (1 - 1 / abs(I_2))(1 - 1 / abs(J_2)) \
    &= abs(Omega) (1 - 1 / abs(I_2))(1 - 1 / abs(J_2)) sigma_"mask"^2.
  $
  Dividing by $|I_2| |J_2|$ proves the identity; PLAID fits both overlapping indicator outer products with zero residual.

  *Part 3:* Place the $b$ zeros at $(1, 1), dots, (b, b)$. An all-ones submatrix $(I, J)$ is inclusion-maximal iff it keeps all clean rows ${b + 1, dots, r}$ and columns ${b + 1, dots, s}$ and, for each $t in {1, dots, b}$, omits *exactly one* of row $t$ or column $t$ (omitting both is non-maximal; keeping both hits $(t, t) = 0$). These $b$ independent choices correspond bijectively to $S subset.eq {1, dots, b}$, yielding $2^b$ maximal blocks of sizes $(r - |S|) times (s - (b - |S|))$. Meanwhile $norm(B - bold(1)_r bold(1)_s^top)_F / norm(bold(1)_r bold(1)_s^top)_F = sqrt(b / (r s))$, so Wedin--Davis--Kahan @davis1970rotation bounds singular-vector drift by $O(sqrt(b / (r s)))$.
]

#subsec("C", "Spectral Asymmetry of ISA vs. Sparse Identifiability of FABIA")

#proposition(title: "Dual-Standardized ISA Nilpotency vs. FABIA Identifiability & MI Ranking", origin: "Original Unified Proof")[
  1. *(Unthresholded ISA $equiv$ SVD vs. Dual-Standardization Nilpotency):* On a single matrix $E_R = E_C = E$ with $t_g = t_c = 0$, ISA reduces to power iteration on symmetric PSD $E E^top$. However, under Bergmann et al.'s dual standardization ($E_R != E_C$) @bergmann2003isa, $T = E_C E_R^top$ is non-symmetric; indeed, on any pure rank-2 additive matrix $E = a bold(1)_p^top + bold(1)_n b^top$ (non-constant $a, b$), dual standardization forces $E_C E_R^top equiv bold(0)$ identically, annihilating the entire signal!
  2. *(Rotation-Free Identifiability & Exact MI Ranking of FABIA \#161):* Unlike Gaussian Factor Analysis—non-identifiable up to orthogonal rotations $(Lambda Q, Q^top Z)$—FABIA's independent Laplacian priors ($kappa_4 = 3 > 0$) guarantee identifiability up to signed permutation (for full-rank $Lambda$), estimate heteroscedastic noise variances $Psi = "diag"(psi_1, dots, psi_n)$, and score each module via exact channel mutual information
     $
       I(x; z_k) = 1 / 2 log(1 + lambda_k^top Psi^(-1) lambda_k).
     $
]

#proof[
  *Part 1:* Setting $E_R = E_C = E$ with zero thresholds and $ell_2$ normalization gives $u^((t+1)) prop E E^top u^((t))$, exact dominant-SVD power iteration. Under dual standardization on $E = a bold(1)_p^top + bold(1)_n b^top$, row-centering annihilates $a bold(1)_p^top$ yielding $E_R prop bold(1)_n tilde(b)^top$ (with $bold(1)_p^top tilde(b) = 0$), while column-centering annihilates $bold(1)_n b^top$ yielding $E_C prop tilde(a) bold(1)_p^top$; thus $E_C E_R^top prop tilde(a) (bold(1)_p^top tilde(b)) bold(1)_n^top = bold(0)$ identically!

  *Part 2:* If $Lambda Z$ and $tilde(Lambda) tilde(Z)$ share the same distribution with independent non-Gaussian Laplacian coordinates and full-column-rank $Lambda, tilde(Lambda)$, then $Z = (Lambda^+ tilde(Lambda)) tilde(Z)$; by the Darmois--Skitovich theorem, $Lambda^+ tilde(Lambda)$ must be a signed permutation matrix (scaling fixed by unit prior variance). For Gaussianized variational factor $z_k tilde.op cal(N)(0, 1)$ in $x = lambda_k z_k + upsilon$ ($upsilon tilde.op cal(N)(0, Psi)$), Sylvester's determinant identity gives
  $
    I(x; z_k) &= 1 / 2 log((det(Psi + lambda_k lambda_k^top)) / (det Psi)) \
    &= 1 / 2 log(1 + lambda_k^top Psi^(-1) lambda_k).
  $
]

#wide_figure(
  caption: [*Table II:* Head-to-Head Theoretical & Algorithmic Comparison of Two-Mode Biclustering Models (Section II: Algs. 157--161).],
  table(
    columns: (1.02fr, 0.96fr, 0.85fr, 0.92fr, 0.96fr, 1.29fr),
    align: (left, left, center, center, left, left),
    [Algorithm], [Coherence Model], [Overlap], [Noise / Dropouts], [Time Complexity], [Distinctive Strength / Failure Mode],
    [Cheng--Church (\#157) @cheng2000biclustering], [Additive ANOVA ($H <= delta$)], [Poor (masking)], [Med (mask bias)], [$O(K n p (n + p))$], [Simple MSR; breaks on mult. & overlap],
    [PLAID (\#158) @lazzeroni2002plaid], [Additive layer sum], [*Native*], [Good ($ell_2$ fit)], [$O(K T n p)$], [Clean overlapping layers; no mult.],
    [Bimax (\#159) @prelic2006bimax], [Exact binary all-1s], [*Exact all*], [*$2^b$ shatters*], [$O(n p beta min(n, p))$], [Complete exact enum.; binarization loss],
    [ISA (\#160) @bergmann2003isa], [Rank-1 thresholded], [Multi-seed], [Good (sparse)], [$O(S T n p)$], [Fast mult. modules; dual-std. nilpotency],
    [FABIA (\#161) @hochreiter2010fabia], [*Sparse mult. factor*], [*Joint native*], [*Heterosced.*], [$O(T (K^3 (n + p) + K n p))$], [Identifiable + MI rank; $O(K^3)$ in $K$],
  ),
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Time-Series Elastic, Shape & Convolutional Invariance")

Four algorithms handle ordered length-$m$ time series $x in bb(R)^m$: *k-Shape* (\#162) @paparrizos2015kshape, *DTW Barycenter Averaging (DBA)* (\#163) @petitjean2011dba @sakoe1978dtw, *Soft-DTW Clustering* (\#164) @cuturi2017softdtw, and *ROCKET-Clustering* (\#165) @dempster2020rocket.

#subsec("A", "Transformation Invariance Groups & Phase vs. Warp Separation")

Each method encodes a distinct temporal invariance:
- *k-Shape (\#162)* computes FFT-based Shape-Based Distance $"SBD"(x, y) = 1 - max_w "NCC"_c (x, y)_w$ on $z$-normalized series over all $2 m - 1$ lags in $O(m log m)$ time.
- *DBA (\#163)* minimizes squared DTW over monotone warping paths $pi in cal(A)_(m, m)$ (steps in ${(1, 0), (0, 1), (1, 1)}$) with path matrix $A_pi$ and cost matrix $Delta_(r s) = (x_r - y_s)^2$:
  $
    "DTW"^2 (x, y) = min_(pi in cal(A)_(m, m)) chevron.l A_pi, Delta(x, y) chevron.r.
  $
- *Soft-DTW (\#164)* replaces hard minimum by temperature-$gamma$ LogSumExp ($gamma > 0$):
  $
    "dtw"_gamma (x, y) = -gamma log sum_(pi in cal(A)_(m, m)) exp( - (chevron.l A_pi, Delta(x, y) chevron.r) / gamma ).
  $
- *ROCKET-Clustering (\#165)* convolves $x$ with $M$ random zero-mean dilated kernels $(w_r, b_r, d_r)$ ($ell_r in {7, 9, 11}$), extracts global $"max"$ and Proportion of Positive Values ($"PPV"$) into $bb(R)^(2 M)$, and runs Euclidean $k$-Means++.

#theorem(title: "Exact Phase vs. Elastic Warp Separation & Non-Metricity Barrier", origin: "Original Unified Proof")[
  Let $m >= 4$ be even.
  1. *(Elastic-Stretch Failure of k-Shape \#162 vs. Exact DBA \#163 Invariance):* Let $x in bb(R)^m$ be the balanced square wave ($+1$ on $[1, m / 2]$, $-1$ on $(m / 2, m]$) and $y^((p))$ its monotone warp switching at $p m$ ($p in (1 / 2, 1)$, $p m in bb(Z)$). Then $"DTW"(x, y^((p))) = 0$ identically, whereas k-Shape's SBD attains its exact minimum at lag $w^* = (p - 1 / 2) m$ with closed-form discrepancy
     $
       "SBD"(x, y^((p))) = 1 - (p + 1 / 2) sqrt((1 - p) / p) > 0,
     $
     strictly positive on $(1 / 2, 1)$ and tending to $1$ as $p -> 1^-$.
  2. *(Rigid Shift Efficiency vs. Band-Limited DBA Failure):* On shifted unit impulses $e_1, e_(1 + tau)$ ($1 <= tau <= m / 2$), k-Shape gives $"SBD" = 0$ in $O(m log m)$ time, whereas Sakoe--Chiba band-limited DBA with bandwidth $r < tau$ fails completely ($"DTW"_r^2 = 2$).
  3. *(Exact Offset Invariance of ROCKET \#165):* Zero-mean kernels ($sum_j w_(r, j) = 0$) render interior responses $(x * w_r)_t + b_r$, hence $"max"$ and $"PPV"$, identically invariant under arbitrary additive offsets $x |-> x + c bold(1)_m$.
  4. *(Universal Non-Metricity of Alignment Distances):* Both $"SBD"$ (\#162) and $"DTW"$ (\#163) violate the triangle inequality (for $x = (0, 1, 1, 2), y = (0, 1, 2, 2), z = (0, 2, 2, 2)$, $"DTW"(x, y) = 0, "DTW"(y, z) = 1$, yet $"DTW"(x, z) = sqrt(2) > 1$), while raw Soft-DTW (\#164) even violates non-negativity ($"dtw"_gamma (x, x) < 0$ for all $m >= 2, gamma > 0$). Only ROCKET's $ell_2$ embedding (\#165) defines a true metric.
]

#proof[
  *Part 1:* Matching $[1, m / 2]$ onto $[1, p m]$ and $(m / 2, m]$ onto $(p m, m]$ gives $"DTW"(x, y^((p))) = 0$. For k-Shape, $tilde(x) = x$, while $y^((p))$ has mean $2 p - 1$ and variance $4 p (1 - p)$, giving $z$-normalized entries $+a$ on $[1, p m]$ and $-b$ on $(p m, m]$ where $a := sqrt((1 - p) / p) < 1 < b := sqrt(p / (1 - p))$. Computing $C(w) = m^(-1) sum_t tilde(x)_t tilde(y)_(t + w)$ (zero-padded):
  - For $0 <= w <= (p - 1 / 2) m$, $C(w) = (1 - p + w / m) a + (1 - p) b$ increases linearly (slope $+a / m$) to its peak at $w^* = (p - 1 / 2) m$:
    $
      C(w^*) &= (1 / 2) a + (1 - p) b \
      &= (p + 1 / 2) sqrt((1 - p) / p).
    $
  - For $(p - 1 / 2) m < w <= m / 2$, $C(w) = (p - w / m) a + (p - 2 w / m) b$ decreases with slope $-(a + 2 b) / m < 0$, and decreases further for $w > m / 2$.
  - For $w = -r < 0$, $C(-r) = (1 - p - 2 r / m) a + (1 - p - r / m) b$ decreases strictly on $1 <= r <= (1 - p) m$ and is $<= 0$ for $r > (1 - p) m$.
  Hence $w^*$ is the unique global maximizer, and $1 - C(w^*)^2 = (p - 1 / 2)(p^2 + p / 2 + 1 / 2) / p > 0$ for all $p in (1 / 2, 1)$.

  *Part 2:* Lag $w = tau$ aligns $e_1$ and $e_(1 + tau)$ exactly ($"SBD" = 0$), whereas $|i - j| <= r < tau$ forbids pairing $(1, 1 + tau)$, forcing $"DTW"_r^2 >= 1 + 1 = 2$.

  *Part 3:* $((x + c bold(1)) * w_r)_t = (x * w_r)_t + c sum_j w_(r, j) = (x * w_r)_t$ identically.

  *Part 4:* Direct DP verifies $"DTW"(x, y) = 0, "DTW"(y, z) = 1, "DTW"(x, z) = sqrt(2) > 0 + 1$. For Soft-DTW ($m >= 2$), the diagonal path has cost $0$ while $|cal(A)_(m, m)| - 1 >= 2$ off-diagonal paths have finite cost $c_pi >= 0$, forcing $"dtw"_gamma (x, x) = -gamma log(1 + sum_(pi != pi_"diag") e^(-c_pi / gamma)) < 0$.
]

#subsec("B", "Fréchet Barycenter Smoothness, Soft-DTW Hessian & ROCKET Scaling")

#theorem(title: "Barycenter Smoothness, Soft-DTW Hessian & ROCKET Kernel Concentration", origin: "Original Unified Proof")[
  Let $C = {x_1, dots, x_n} subset.eq bb(R)^m$.
  1. *(Spectral Optimality of k-Shape \#162 vs. Combinatorial Traps of DBA \#163):* Holding lags fixed, k-Shape's centroid update is solved globally by the leading eigenvector of $M = Q^top (sum_i x_i^((w_i)) (x_i^((w_i)))^top) Q$ ($Q = I_m - m^(-1) bold(1) bold(1)^top$), automatically zero-mean @paparrizos2015kshape. By contrast, DBA's inertia $F_0 (c) = n^(-1) sum_i min_pi chevron.l A_pi, Delta(c, x_i) chevron.r$ is the pointwise minimum of $D(m - 1, m - 1)^n$ quadratics, where the Delannoy path count grows as
     $
       |cal(A)_(m, m)| tilde.op ((3 + 2 sqrt(2))^(m - 1)) / sqrt(pi (3 sqrt(2) - 4)(m - 1)),
     $
     creating non-differentiable kinks where DBA's MM update can stall at suboptimal local fixed points.
  2. *(Entropic Duality, Uniform Envelope & Exact Hessian of Soft-DTW \#164):* For $gamma > 0$, Soft-DTW obeys Legendre--Fenchel duality over $p in Delta(cal(A)_(m, m))$:
     $
       "dtw"_gamma (c, x) = min_p { sum_pi p_pi chevron.l A_pi, Delta(c, x) chevron.r - gamma H(p) },
     $
     yielding the exact uniform DTW envelope
     $
       "DTW"^2 (c, x) - gamma log |cal(A)_(m, m)| \
       quad <= "dtw"_gamma (c, x) <= "DTW"^2 (c, x)
     $
     (asymptotic gap $approx 1.763 gamma m$). Moreover, $c |-> "dtw"_gamma (c, x)$ is $C^infinity$-smooth with exact gradient and Hessian
     $
       nabla_c "dtw"_gamma (c, x) &= 2 ("diag"(E_gamma bold(1)) c - E_gamma x), \
       nabla_c^2 "dtw"_gamma (c, x) &= 2 "diag"(E_gamma bold(1)) \
       &quad - 1 / gamma "Cov"_(p_gamma) (g_pi (c, x)),
     $
     where $E_gamma = bb(E)_(p_gamma)[A_pi]$ ($O(m^2)$ backward DP @cuturi2017softdtw) and $g_pi := nabla_c chevron.l A_pi, Delta(c, x) chevron.r$. Since $2 "diag"(E_gamma bold(1)) succ.eq 2 I_m$, raising $gamma$ shrinks non-convex path-covariance curvature as $O(1 / gamma)$, achieving local convexity once $gamma > (1 / 2) lambda_"max" ("Cov"_(p_gamma) (g_pi))$!
  3. *(Linear-in-$m$ Scaling & Uniform Concentration of ROCKET \#165):* Extracting $M$ random $"PPV"$ features costs $O(n M m)$ (linear in $m$, vs. $O(T I n m^2)$ for DBA/Soft-DTW). Letting $hat(k)_M (x, y) = M^(-1) sum_(r=1)^M phi_r (x) phi_r (y)$ and $k_"PPV" = bb(E)[hat(k)_M]$, Hoeffding's bound guarantees with probability $>= 1 - delta$:
     $
       max_(1 <= i, j <= n) |hat(k)_M (x_i, x_j) - k_"PPV" (x_i, x_j)| \
       quad <= sqrt( (log((n (n + 1)) / delta)) / (2 M) ),
     $
     reducing $k$-Means++ on ROCKET features to $O(log K)$-approximate convolutional Kernel $k$-Means @arthur2007kmeanspp @dempster2020rocket.
]

#proof[
  *Part 1:* By Rayleigh--Ritz, maximizing $mu^top M mu$ on $bb(S)^(m - 1)$ yields the top eigenvector of $M$, orthogonal to $bold(1) in ker(Q) subset.eq ker(M)$. Delannoy asymptotics count $(1, 0), (0, 1), (1, 1)$ lattice paths from $(1, 1)$ to $(m, m)$; pointwise minima of convex quadratics have non-differentiable kinks at alignment ties.

  *Part 2:* Writing $f_pi (c) := chevron.l A_pi, Delta(c, x) chevron.r$, Gibbs' variational identity gives
  $
    -gamma log sum_pi e^(-f_pi / gamma) = min_p { sum_pi p_pi f_pi - gamma H(p) },
  $
  attained at $p_(gamma, pi) = e^(-f_pi / gamma) / Z_gamma$. Bounding $0 <= H(p) <= log |cal(A)_(m, m)|$ gives the exact two-sided envelope. Differentiating $-gamma log Z_gamma (c)$ yields $nabla_c "dtw"_gamma = bb(E)_(p_gamma)[g_pi] = 2("diag"(E_gamma bold(1)) c - E_gamma x)$. Differentiating again via $nabla_c p_(gamma, pi) = -gamma^(-1) p_(gamma, pi) (g_pi - bb(E)_(p_gamma)[g_pi])$ yields $nabla_c^2 "dtw"_gamma = 2 "diag"(E_gamma bold(1)) - gamma^(-1) "Cov"_(p_gamma)(g_pi)$; because every admissible path visits every row $r$ at least once, $(E_gamma bold(1))_r >= 1$, hence $2 "diag"(E_gamma bold(1)) succ.eq 2 I_m$.

  *Part 3:* Each $phi_r (x_i) phi_r (x_j) in [0, 1]$ is i.i.d. across $r = 1, dots, M$; Hoeffding's inequality plus a union bound over all $(n (n + 1)) / 2$ unordered pairs gives $n (n + 1) exp(-2 M t^2) = delta => t = sqrt(log((n (n + 1)) / delta) / (2 M))$. Each fixed-length convolution costs $O(m)$, yielding $O(n M m + T K n M)$ total time.
]

#remark(title: "Soft-DTW Divergence Debiasing")[
  Because Theorem 4(4) gives $"dtw"_gamma (x, x) < 0$, modern Soft-DTW Clustering (\#164) optimizes the *Soft-DTW divergence* $D_gamma (c, x) = "dtw"_gamma (c, x) - (1 / 2)("dtw"_gamma (c, c) + "dtw"_gamma (x, x)) >= 0$, eliminating self-shrinkage bias while keeping $C^infinity$ smoothness and $O(m^2)$ gradients.
]

#wide_figure(
  caption: [*Table III:* Head-to-Head Theoretical & Algorithmic Comparison of Time-Series Clustering Methods (Section III: Algs. 162--165).],
  table(
    columns: (1.02fr, 0.94fr, 0.98fr, 0.72fr, 1.06fr, 1.28fr),
    align: (left, left, left, center, left, left),
    [Algorithm], [Dissimilarity / Space], [Invariances], [Metric?], [Time Complexity], [Distinctive Strength / Failure Mode],
    [k-Shape (\#162) @paparrizos2015kshape], [SBD ($1 - max "NCC"_c$)], [Shift, scale, rigid lag], [No], [$O(T (n K m log m + n m^2 + K m^3))$], [Fast exact spectral mean; fails on warps],
    [DBA (\#163) @petitjean2011dba], [Squared DTW], [*Nonlinear time warp*], [No], [$O(T (K + I) n m^2)$], [Sharp warp alignment; non-smooth traps],
    [Soft-DTW (\#164) @cuturi2017softdtw], [Smoothed $"dtw"_gamma$], [Smooth elastic warp], [No ($< 0$)], [$O(T (K + I) n m^2)$], [$C^infinity$ smooth mean; $O(gamma m)$ bias, $O(m^2)$],
    [ROCKET (\#165) @dempster2020rocket], [Random conv. $bb(R)^(2 M)$], [Offset, multi-scale], [*Yes ($ell_2$)*], [$O(n M m + T K n M)$], [*Linear in $m$*; loses direct raw barycenter],
  ),
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Spatial Vision, 3D Clouds, Trajectories, Ensembles & Multi-View")

The final ten algorithms (Algs. 166--175) span three structured regimes: (A) spatial grids, 3D LiDAR clouds, and trajectories; (B) resampling stability and hypergraph cluster ensembles; and (C) multi-view feature and spectral fusion.

#subsec("A", "Spatial Vision, 3D Point Clouds & Trajectories (Algs. 166--169)")

*SLIC (\#166)* @achanta2012slic runs 5D CIELAB--$(x, y)$ $k$-Means ($D^2 = d_"lab"^2 + (m / S)^2 d_"xy"^2$, $S = sqrt(N / k)$) inside local $2 S times 2 S$ windows ($O(N I)$ independent of $k$). *Watershed (\#167)* @vincent1991watersheds floods gradient reliefs $g: V -> bb(R)$ from regional minima via geodesic influence zones in $O(N log N)$ time. *DBSCAN-3D / PCL ECE (\#168)* @rusu2011pcl @ester1996dbscan applies $k$-d tree $d_"th"$-ball connectivity ($"MinPts" = 1$ for ECE, $>= 2$ for DBSCAN-3D) to 3D point clouds. *TRACLUS (\#169)* @lee2007traclus splits trajectories by Minimum Description Length ($"MDL" = L(H) + L(D | H)$) and clusters line segments via DBSCAN under composite distance
$
  d(L_i, L_j) &= w_perp d_perp (L_i, L_j) + w_parallel d_parallel (L_i, L_j) \
  &quad + w_theta d_theta (L_i, L_j).
$

#theorem(title: "Spatial Edge Adherence, LiDAR Fragmentation & Sub-Trajectory Separation", origin: "Original Unified Proof")[
  1. *(SLIC \#166 Step-Edge Law vs. Watershed \#167 Monotone Invariance & Noise Breakdown):*
     - On an injective pixel relief $g$, Watershed is *identically invariant* under every strictly increasing transform $phi compose g$ and coincides with the minima-rooted Minimum Spanning Forest (MSF) of edge weights $w(u, v) = max(g(u), g(v))$. However, injecting $r$ noise minima of arbitrarily small depth $epsilon > 0$ spawns $r$ spurious basins (*zero-threshold breakdown*), requiring $h$-minima filtering ($h > 2 epsilon$).
     - In SLIC, every pixel tests $O(1)$ centers ($O(N I)$ total). Across a straight Lab step edge of contrast $Delta_c$ shifted by $delta in (0, S / 2)$ from the midpoint of two seeds spaced $S$ apart, the boundary pixel snaps to its color-matching superpixel *if and only if* $Delta_c^2 > (2 m^2 delta) / S$ (requiring $Delta_c > m$ at $delta = S / 2$).
  2. *(EMST Equivalence & LiDAR Critical Range of DBSCAN-3D / ECE \#168):* PCL ECE ($"MinPts" = 1$) outputs identically the connected components of the 3D EMST after cutting edges $> d_"th"$ (vulnerable to 1-point chaining bridges suppressed by DBSCAN-3D with $"MinPts" >= 4$). For a spinning LiDAR scanner with angular steps $(Delta phi, Delta theta)$ illuminating a surface at range $r$ and incidence angle $psi$, return density decays as $rho(r) = (cos psi) / (r^2 Delta phi Delta theta)$, inducing an *exact critical fragmentation range*
     $
       r_"crit" = d_"th" / max(Delta phi, (Delta theta) / (cos psi)),
     $
     beyond which fixed-radius extraction shatters continuous targets into disconnected rings.
  3. *(Sub-Trajectory Discovery: TRACLUS \#169 vs. Whole-Trajectory DTW \#163):* On $n >= "MinLns"$ planar trajectories sharing a central segment $[(0, 0), (L, 0)]$ while diverging radially by length $R$ at angles $theta_i = (2 pi i) / n$, whole-trajectory DTW incurs $"DTW"^2 (T_i, T_j) >= 4 R^2 sin^2 (pi / n) -> +infinity$, whereas TRACLUS's MDL partitioning isolates the shared segment with distance $0$ across all $n$ trajectories (at the price of non-metric segment distances and lost global trajectory order).
]

#proof[
  *Part 1:* Immersion processes pixels in increasing order of $g$; any strictly increasing $phi$ preserves all level sets and regional minima, and the first minimum flooding $u$ minimizes the minimax path height $max_(v in P) g(v)$, characterizing the minima-rooted MSF of $w(u, v) = max(g(u), g(v))$. Any spurious minimum of depth $epsilon > 0$ claims at least its own vertex before water reaches height $+epsilon$. For SLIC with seeds at $(plus.minus S / 2, 0)$ of colors $c_1, c_2$ ($norm(c_1 - c_2)_2 = Delta_c$) and a color-$c_1$ pixel at $(delta, 0)$, $D_1^2 < D_2^2$ holds iff
  $
    Delta_c^2 &> (m / S)^2 [ (S / 2 + delta)^2 - (S / 2 - delta)^2 ] \
    &= (2 m^2 delta) / S,
  $
  giving $Delta_c > m$ at $delta = S / 2$.

  *Part 2:* By the MST cycle property, $u, v$ are connected in the $d_"th"$-ball graph iff the EMST path between them has maximum edge $<= d_"th"$. Sparse bridges of spacing $(2 / 3) d_"th"$ keep degree $<= 2$ and fail core status when $"MinPts" >= 4$. At LiDAR range $r$, tangential azimuth spacing is $r Delta phi$ and elevation-ring surface spacing is $(r Delta theta) / (cos psi)$, yielding area per point $(r^2 Delta phi Delta theta) / (cos psi)$ and disconnecting the $d_"th"$-graph whenever $r > r_"crit"$.

  *Part 3:* Because DTW pins startpoints, $"DTW"^2 (T_i, T_j) >= norm(R e^(i theta_i) - R e^(i theta_j))_2^2 >= 4 R^2 sin^2 (pi / n) -> infinity$. Conversely, MDL splits every sharp bend at $(0, 0)$ and $(L, 0)$ (avoiding $Omega(log R)$ unsplit coding cost), yielding $n$ identical copies of $[(0, 0), (L, 0)]$ at zero pairwise distance.
]

#subsec("B", "Resampling Stability & Cluster Ensembles (Algs. 170--173)")

*Monti Consensus (\#170)* @monti2003consensus averages co-clustering over $H$ subsamples of fraction $f in (0, 1)$, forming $M_(i j) = (sum_h M_(i j)^((h))) / (sum_h I_(i j)^((h)))$, diagnosing $K$ via empirical CDF area $A(K)$, and clustering $1 - M$ by UPGMA. Given $r$ base partitions $Lambda = {lambda^((1)), dots, lambda^((r))}$ encoded in binary indicator matrix $H = [H^((1)), dots, H^((r))] in {0, 1}^(n times m)$ ($m = sum_q k^((q))$), Strehl and Ghosh @strehl2002cluster target *Average Normalized Mutual Information*:
$
  Phi^("ANMI") (Lambda, lambda) &:= 1 / r sum_(q=1)^r "NMI"(lambda, lambda^((q))), \
  "NMI"(lambda, lambda^((q))) &:= (I(lambda; lambda^((q)))) / sqrt(H(lambda) H(lambda^((q)))),
$
via three hypergraph solvers: *CSPA (\#171)* partitions object Gram matrix $S = r^(-1) H H^top$ via METIS @karypis1998metis ($O(n^2 (r + k))$); *HGPA (\#172)* directly cuts $m$ unit-weight hyperedges ($O(n k r)$); and *MCLA (\#173)* partitions the *dual* $m times m$ cluster-space Jaccard meta-graph derived from $H^top H$, collapses meta-clusters into soft columns $A_(:, c)$, and assigns objects by competition ($O(n k^2 r^2)$, linear in $n$).

#theorem(title: "Information-Theoretic Consensus Optimality, Dual Gram Symmetry & Failure Modes", origin: "Original Unified Proof")[
  Let $Lambda = {lambda^((1)), dots, lambda^((r))}$ have indicator matrix $H in {0, 1}^(n times m)$.
  1. *(ANMI Information-Theoretic Optimality & Supra-Consensus Dominance @strehl2002cluster):* $Phi^("ANMI") (Lambda, lambda) in [0, 1]$, attaining $1$ iff all $lambda^((q))$ and $lambda$ agree up to label permutation. When base partitions are conditionally independent noisy observations of a balanced ground-truth partition $lambda_"true"$ through symmetric $k$-ary label-noise channels ($p < (k - 1) / k$), maximizing average mutual information over balanced partitions is strictly consistent for $lambda_"true"$ as $r -> infinity$, and the supra-consensus rule $lambda^("supra") = arg max_(lambda in {lambda_"CSPA", lambda_"HGPA", lambda_"MCLA"}) Phi^("ANMI") (Lambda, lambda)$ provably dominates every individual solver.
  2. *(Dual Gram Spectral Equivalence & Mirkin Median Identity of CSPA \#171 vs. MCLA \#173):* CSPA's $S = r^(-1) H H^top$ is PSD of rank $<= m - r + 1$, and maximizing $sum_(c=1)^k sum_(i, j in C_c) (2 S_(i j) - 1)$ is algebraically identical to minimizing average Mirkin partition distance $r^(-1) sum_q norm(Z_lambda Z_lambda^top - H^((q)) (H^((q)))^top)_F^2$. Since $H H^top$ (CSPA) and $H^top H$ (MCLA) share *identical* nonzero eigenvalues $sigma_ell^2$, MCLA recovers the same low-rank consensus subspace in $O(n k^2 r^2)$ linear time instead of CSPA's $O(n^2)$.
  3. *(Failure Trilemma Across HGPA \#172, CSPA \#171 & Monti \#170):*
     - *(HGPA Cut Saturation):* Under unit hyperedge cut $"hcut" = sum_a bb(I)[|{lambda(i) : i in C_a}| >= 2]$, a *single* noisy vertex in base cluster $C_a$ incurs cut $1$—identical to an even 50--50 split—saturating $"hcut" = m$ with zero gradient under mild base noise.
     - *(CSPA Unanimous Imbalance Failure):* On a noise-free unanimous ensemble with unequal cluster sizes $n_1 > (1 + epsilon_"bal") n / 2$, balanced METIS partitioning of $S$ is infeasible at the unanimous truth and misassigns at least $n_1 - (1 + epsilon_"bal") n / 2$ points.
     - *(Monti Non-PSD Subsampling & CDF Blindness):* Subsampled Monti $M$ ($f < 1$) can have strictly negative eigenvalues, and its CDF area satisfies the exact identity
       $
         A(K) := integral_0^1 "CDF"(c) dif c = 1 - 1 / binom(n, 2) sum_(i < j) M_(i j),
       $
       depending *exclusively* on the mean $overline(M)$ and assigning identical area $A = 1 / 2$ to clean binary ${0, 1}$ stability and completely ambiguous $M_(i j) equiv 1 / 2$!
]

#proof[
  *Part 1:* Cauchy--Schwarz on $I(lambda; lambda^((q))) <= min(H(lambda), H(lambda^((q))))$ yields $"NMI" in [0, 1]$, equal to $1$ iff $H(lambda | lambda^((q))) = H(lambda^((q)) | lambda) = 0$. Under conditionally independent symmetric channels given $lambda_"true"$, the data-processing inequality $I(lambda; lambda^((q))) <= I(lambda_"true"; lambda^((q)))$ holds strictly for every $lambda != lambda_"true"$ up to relabeling; by the law of large numbers over $r -> infinity$, maximizing $r^(-1) sum_q I(lambda; lambda^((q)))$ (equivalent to $Phi^("ANMI")$ on balanced partitions with fixed $H(lambda) = log k$) uniquely recovers $lambda_"true"$. Supra-consensus dominance follows by taking the maximum over the three candidates.

  *Part 2:* $v^top S v = r^(-1) norm(H^top v)_2^2 >= 0$; the $r - 1$ constraints $H^((q)) bold(1) = H^((1)) bold(1)$ bound $"rank"(S) <= m - r + 1$. Since $Z = Z_lambda Z_lambda^top in {0, 1}^(n times n)$ satisfies $Z_(i j)^2 = Z_(i j)$,
  $
    1 / r sum_(q=1)^r norm(Z - H^((q)) (H^((q)))^top)_F^2 &= C_Lambda + sum_(i, j) Z_(i j) (1 - 2 S_(i j)) \
    &= C_Lambda - sum_(c=1)^k sum_(i, j in C_c) (2 S_(i j) - 1).
  $
  Via SVD $H = U Sigma V^top$, $H H^top = U Sigma^2 U^top$ and $H^top H = V Sigma^2 V^top$ share exact nonzero eigenvalues $Sigma^2$.

  *Part 3:* One outlier vertex in $C_a$ raises $|{lambda(i) : i in C_a}|$ from $1$ to $2$, matching the maximum unit cut $1$. Balanced CSPA enforces $|C_c| <= (1 + epsilon_"bal") n / 2 < n_1$, severing at least $n_1 - (1 + epsilon_"bal") n / 2$ points. For Monti with $n = 3$ and pair subsamples ${1, 2}$ (together), ${2, 3}$ (together), ${1, 3}$ (apart), $M = mat(1, 1, 0; 1, 1, 1; 0, 1, 1)$ has eigenvalue $1 - sqrt(2) < 0$. Sorting the $N_p = binom(n, 2)$ off-diagonal entries $0 =: c_0 <= c_1 <= dots <= c_(N_p) <= c_(N_p + 1) := 1$ gives
  $
    A(K) &= sum_(ell=1)^(N_p) ell / N_p (c_(ell + 1) - c_ell) \
    &= 1 - 1 / N_p sum_(ell=1)^(N_p) c_ell = 1 - overline(M),
  $
  blind to whether entries concentrate at ${0, 1}$ or at $1 / 2$.
]

#subsec("C", "Multi-View Partitional vs. Co-Regularized Spectral (Algs. 174--175)")

Given $V >= 2$ raw views $X_v in bb(R)^(n times d_v)$ ($D = sum_v d_v$), *Multi-View $k$-Means (MVKM, \#174)* @cai2013multiview couples shared partition $G in {0, 1}^(n times k)$ to view centroids ${C_v}$ with power-weighted simplex weights $alpha in Delta_(V - 1)$ ($gamma > 1$): $J_gamma = sum_(v=1)^V alpha_v^gamma E_v$, $E_v := norm(X_v - G C_v)_F^2$. *Co-regularized Multi-View Spectral Clustering (\#175)* @kumar2011coregularized builds normalized similarities $S_v = D_v^(-1/2) W_v D_v^(-1/2)$ @chung1997spectral @vonluxburg2007tutorial and maximizes Grassmannian-coupled spectral energy over $U_v in "St"(n, k)$ ($P_v = U_v U_v^top$):
$
  cal(F)_"pair" ({U_v}) &= sum_(v=1)^V "Tr"(U_v^top S_v U_v) \
  &quad + sum_(1 <= v < w <= V) lambda_(v w) "Tr"(P_v P_w).
$

#theorem(title: "Power-Mean Weighting, Grassmannian Strong-Coupling Limit & Hijacking", origin: "Original Unified Proof")[
  1. *(Power-Mean Reduction & Scale Hijacking of MVKM \#174):* For fixed $(G, {C_v})$ with $E_v > 0$, minimizing $J_gamma$ over $alpha in Delta_(V - 1)$ ($gamma > 1$) yields unique weights $alpha_v^* = (E_v^(-1 / (gamma - 1))) / (sum_w E_w^(-1 / (gamma - 1)))$ and exact power-mean objective
     $
       min_(alpha in Delta_(V - 1)) J_gamma = ( sum_(v=1)^V E_v^(-1 / (gamma - 1)) )^(-(gamma - 1)),
     $
     interpolating $min_v E_v$ ($gamma -> 1^+$) and equal averaging ($gamma -> infinity$). Scaling a pure-noise view by $s -> 0^+$ ($E_"noise" = s^2 E_0$) drives
     $
       alpha_"noise"^* = (1 + sum_(v != "noise") (s^2 E_0 / E_v)^(1 / (gamma - 1)))^(-1) -> 1,
     $
     hijacking consensus unless views are variance-normalized!
  2. *(Ky Fan Ascent, Strong-Coupling Variance Reduction & Fixed-$lambda$ Hijacking of \#175):* Alternating block updates solve $max_(U_v in "St"(n, k)) "Tr"(U_v^top tilde(S)_v U_v)$ for $tilde(S)_v = S_v + sum_(w != v) lambda_(v w) P_w$, guaranteeing monotone $O(k)$-invariant ascent in $O(T V n^3)$ time. As $lambda_(v w) = lambda -> infinity$, all subspaces collapse to $U^*$ maximizing $"Tr"(U^top (sum_v S_v) U)$; under independent perturbations $S_v = S^* + Delta_v$ with eigengap $delta_* > 0$, Davis--Kahan @davis1970rotation yields
     $
       norm(sin Theta(hat(U)^*, U^*))_F <= (2 norm(sum_(v=1)^V Delta_v)_F) / (V delta_*),
     $
     shrinking error as $O(1 / sqrt(V))$. Yet under fixed uniform $lambda > g_1$ on $V = 2$ orthogonal block views $S_1 = g_1 P_"true"$ and adversarial $S_2 = g_2 P_"adv"$ ($g_2 > g_1$), the global maximizer collapses onto the *spurious higher-gap view* $(P_"adv", P_"adv")$!
]

#proof[
  *Part 1:* Hölder's inequality with $p = gamma, q = gamma / (gamma - 1)$ gives
  $
    1 &= sum_(v=1)^V (alpha_v E_v^(1 / gamma))(E_v^(-1 / gamma)) \
    &<= (sum_(v=1)^V alpha_v^gamma E_v)^(1 / gamma) (sum_(v=1)^V E_v^(-1 / (gamma - 1)))^((gamma - 1) / gamma).
  $
  Raising to power $gamma$ yields the exact lower bound attained uniquely at $alpha_v^* prop E_v^(-1 / (gamma - 1))$. Substituting $E_"noise" = s^2 E_0 -> 0^+$ sends $alpha_"noise"^* -> 1$ while clean terms scale as $O(s^(2 gamma / (gamma - 1))) = o(s^2)$, handing the objective to noise.

  *Part 2:* Terms depending on $U_v$ equal $"Tr"(U_v^top tilde(S)_v U_v)$, maximized globally on $"St"(n, k)$ by top-$k$ eigenvectors of symmetric $tilde(S)_v$ (Ky Fan principle), invariant under $U_v |-> U_v Q_v$. Rewriting
  $
    cal(F)_"pair" &= sum_(v=1)^V "Tr"(P_v S_v) \
    &quad + (lambda / 2) (V (V - 1) k - sum_(v < w) norm(P_v - P_w)_F^2)
  $
  forces $norm(P_v - P_w)_F^2 = O(1 / lambda) -> 0$, collapsing to $max_P "Tr"(P sum_v S_v)$. Applying Davis--Kahan @davis1970rotation to $overline(S) = S^* + V^(-1) sum_v Delta_v$ gives the stated bound and $O(1 / sqrt(V))$ independent-noise contraction. For $V = 2$ with $P_"true" P_"adv" = bold(0)$, decoupled $(P_"true", P_"adv")$ scores $k (g_1 + g_2)$ whereas aligned $(P_"adv", P_"adv")$ scores $k (g_2 + lambda) > k (g_1 + g_2)$ whenever $lambda > g_1$, strictly hijacking both views!
]

#wide_figure(
  caption: [*Table IV:* Head-to-Head Theoretical & Algorithmic Comparison of Spatial Vision, 3D, Trajectory, Ensemble & Multi-View Methods (Section IV: Algs. 166--175).],
  {
    set text(size: 7.3pt)
    set table(inset: (x: 3.2pt, y: 2.3pt))
    table(
      columns: (1.06fr, 0.88fr, 0.96fr, 0.92fr, 0.94fr, 1.24fr),
      align: (left, left, left, left, left, left),
      [Algorithm], [Domain / Input], [Core Mechanism], [Invariance / Opt.], [Time Complexity], [Distinctive Strength / Failure Mode],
      [SLIC (\#166) @achanta2012slic], [2D RGB/Lab image], [$2 S times 2 S$ 5D $k$-Means], [Grid shift; $O(N)$], [$O(N I)$ (indep. $k$)], [Fast compact superpixels; $m$ edge trade-off],
      [Watershed (\#167) @vincent1991watersheds], [Scalar relief grid], [Immersion / MSF], [*All mono. $phi(g)$*], [$O(N log N)$], [Exact topology; $0$-threshold over-seg.],
      [DBSCAN-3D (\#168) @rusu2011pcl], [3D LiDAR cloud], [$k$-d tree EMST / core], [$E(3)$ rigid], [$O(N log N)$], [Fast 3D objects; $r_"crit"$ range split],
      [TRACLUS (\#169) @lee2007traclus], [2D/3D polylines], [MDL split + DBSCAN], [Sub-traj. match], [$O(M^2)$], [Shared corridors; non-metric, no order],
      [Monti (\#170) @monti2003consensus], [Subsampled runs], [Co-assoc. + UPGMA], [Resample stab.], [$O(H (T_"base" + n^2))$], [Stability K-choice; non-PSD, CDF blind],
      [CSPA (\#171) @strehl2002cluster], [Cluster ensemble], [$S = r^(-1) H H^top$ cut], [Mirkin median], [$O(n^2 (r + k))$], [Exact Gram fusion; $O(n^2)$ + balance fail],
      [HGPA (\#172) @strehl2002cluster], [Cluster ensemble], [Direct hypergraph cut], [Label perm.], [$O(n k r)$], [Fastest linear pass; $0$-gradient cut sat.],
      [MCLA (\#173) @strehl2002cluster], [Cluster ensemble], [Dual $H^top H$ Jaccard], [Dual Gram + soft], [$O(n k^2 r^2)$], [*Linear in $n$* + soft; small-$k$ meta-bias],
      [MVKM (\#174) @cai2013multiview], [Multi-view $bb(R)^(d_v)$], [$alpha_v^gamma$-weighted $k$-Means], [Power-mean opt.], [$O(T n k D)$], [Linear adaptive views; scale hijacking],
      [Co-Reg Spec (\#175) @kumar2011coregularized], [Multi-view graphs], [Grassmannian Ky Fan], [$O(k)$ subspace], [$O(T V n^3)$], [Non-convex manifolds; $O(n^3)$ + fixed-$lambda$ hijack],
    )
  },
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 25-Algorithm Matrix & Practitioner's Decision Guide")

Table V synthesizes *all twenty-five* algorithms of Category VII (Algs. 151--175), followed by exact cross-algorithm reductions, an operational decision guide grounded in Theorems 1--8, and open theoretical frontiers.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 25 Fuzzy, Biclustering, Time-Series, and Domain Clustering Algorithms (Category VII, Algorithms 151--175). Subfamily codes: *FZ* = Fuzzy / Possibilistic, *BC* = Biclustering, *TS* = Time-Series, *SP* = Spatial / 3D / Trajectory, *EN* = Ensemble, *MV* = Multi-View.],
  {
    set text(size: 6.4pt)
    set table(inset: (x: 2.3pt, y: 1.85pt))
    table(
      columns: (1.22fr, 0.34fr, 0.86fr, 0.98fr, 1.04fr, 0.66fr, 1.40fr),
      align: (left, center, left, left, left, center, left),
      [Algorithm], [Subf.], [Data Domain], [Core Geometry / Model], [Time Complexity], [Robust.], [Distinctive Strength / Critical Weakness],
      [\#151 FCM @bezdek1981fcm], [FZ], [Cont. $bb(R)^d$], [Simplex $u_(i j)^m norm(dot)_A^2$], [$O(n c d T)$], [Low ($1/n$)], [Smooth soft partition; equidistant noise trap],
      [\#152 PCM @krishnapuram1993pcm], [FZ], [Cont. $bb(R)^d$], [Decoupled typicality $t_(i j)$], [$O(n c d T)$], [*Redesc.*], [Immune to outliers; coincident global collapse],
      [\#153 GK @gustafson1979gk], [FZ], [Cont. $bb(R)^d$], [Adaptive $A_j$ ($det = rho_j$)], [$O((n c d^2 + c d^3) T)$], [Low], [Full $"GL"(d, bb(R))$ ellipsoids; fixed $rho_j$, rank drop],
      [\#154 GG @gath1989gg], [FZ], [Cont. $bb(R)^d$], [Exponential FMLE $d_"GG"^2$], [$O(T c (n d^2 + d^3))$], [Low], [Unequal vol./priors (= GMM at $m=2$); singular],
      [\#155 KFCM @zhang2004kfcm], [FZ], [Cont. / Kernel], [RKHS / $2(1 - K(x, v))$], [$O(T c n^2)$], [*Gauss. red.*], [Separates shells + robust; $O(n^2)$ kernel cost],
      [\#156 sFCM @chuang2006sfcm], [FZ], [2D/3D Grid], [Filtered $u_(i j)^p h_(i j)^q$], [$O(T c n (d + |cal(N)|))$], [*Impulse*], [1-step pixel denoising; non-monotone, erases 1px],
      [\#157 Cheng--Church @cheng2000biclustering], [BC], [Matrix $bb(R)^(n times p)$], [Additive MSR $H <= delta$], [$O(K n p (n + p))$], [Low], [Fast greedy ANOVA blocks; fails mult. & overlap],
      [\#158 PLAID @lazzeroni2002plaid], [BC], [Matrix $bb(R)^(n times p)$], [Additive layer sum], [$O(K T n p)$], [Med], [Clean overlapping additive blocks; no mult.],
      [\#159 Bimax @prelic2006bimax], [BC], [Binary ${0, 1}^(n times p)$], [Maximal all-1s Galois], [$O(n p beta min(n, p))$], [*Fragile*], [Exact exhaustive enumeration; $2^b$ dropout shatter],
      [\#160 ISA @bergmann2003isa], [BC], [Matrix $bb(R)^(n times p)$], [Thresholded bipartite SVD], [$O(S T n p)$], [Med--High], [Fast mult. modules; dual-std. nilpotency],
      [\#161 FABIA @hochreiter2010fabia], [BC], [Matrix $bb(R)^(n times p)$], [Sparse Laplace $Lambda Z + Upsilon$], [$O(T (K^3 (n + p) + K n p))$], [*Heterosc.*], [Identifiable mult. + MI rank; $O(K^3)$ in factors],
      [\#162 k-Shape @paparrizos2015kshape], [TS], [Series $bb(R)^m$], [FFT SBD + Rayleigh mean], [$O(T (n K m log m + n m^2 + K m^3))$], [Med], [Exact shift/scale invariance; fails elastic warp],
      [\#163 DBA @petitjean2011dba], [TS], [Series $bb(R)^m$], [DTW alignment MM mean], [$O(T (K + I) n m^2)$], [Low], [Exact nonlinear warp alignment; non-smooth traps],
      [\#164 Soft-DTW @cuturi2017softdtw], [TS], [Series $bb(R)^m$], [LogSumExp $"dtw"_gamma$], [$O(T (K + I) n m^2)$], [Med], [$C^infinity$ differentiable mean; $O(gamma m)$ bias, $O(m^2)$],
      [\#165 ROCKET @dempster2020rocket], [TS], [Series $bb(R)^m$], [Random conv. max + PPV], [$O(n M m + T K n M)$], [*High*], [Linear in $m$, true $ell_2$ metric; no raw barycenter],
      [\#166 SLIC @achanta2012slic], [SP], [2D/3D Image], [Windowed 5D Lab--$x y$], [$O(N I)$], [Med], [$O(N)$ indep. of $k$; $Delta_c > m$ edge trade-off],
      [\#167 Watershed @vincent1991watersheds], [SP], [Relief Grid], [Immersion / MSF basins], [$O(N log N)$], [*0-thresh.*], [Monotone relief invariant; severe over-seg. w/o $h$],
      [\#168 DBSCAN-3D @rusu2011pcl], [SP], [3D Point Cloud], [$k$-d tree EMST / core], [$O(N log N)$], [Core-dep.], [Fast LiDAR objects; $r_"crit"$ range fragmentation],
      [\#169 TRACLUS @lee2007traclus], [SP], [Trajectories], [MDL segments + DBSCAN], [$O(M^2)$], [Med], [Finds local shared corridors; non-metric, no order],
      [\#170 Monti @monti2003consensus], [EN], [Any + Base Alg.], [Subsampled co-assoc. CDF], [$O(H (T_"base" + n^2) + n^2 log n)$], [High], [Resampling stab. K-choice; non-PSD, CDF blind],
      [\#171 CSPA @strehl2002cluster], [EN], [Base Partitions], [Object Gram $r^(-1) H H^top$], [$O(n^2 (r + k))$], [Med], [Mirkin median link; $O(n^2)$ & imbalance failure],
      [\#172 HGPA @strehl2002cluster], [EN], [Base Partitions], [Unit hypergraph $k$-cut], [$O(n k r)$], [*Low*], [Fastest linear ensemble; $0$-gradient cut saturation],
      [\#173 MCLA @strehl2002cluster], [EN], [Base Partitions], [Dual Jaccard meta-graph], [$O(n k^2 r^2)$], [*High*], [Linear in $n$ + soft labels; meta-balance bias],
      [\#174 MVKM @cai2013multiview], [MV], [Multi-view $bb(R)^(d_v)$], [$alpha_v^gamma$-weighted $k$-Means], [$O(T n k D)$], [Scale-sens.], [Linear adaptive power-mean; scale hijacking],
      [\#175 Co-Reg Spec @kumar2011coregularized], [MV], [Multi-view Graph], [Grassmannian $P_v P_w$], [$O(T V n^3)$], [Gap-sens.], [Non-convex manifolds; $O(n^3)$ & fixed-$lambda$ hijack],
    )
  },
)

#subsec("A", "Exact Parametric Reductions & Cross-Paradigm Equivalence Lattice")

Theorems 1--8 interlock all twenty-five algorithms into exact limiting bridges across the monograph:
- *Fuzzy & Kernel Limits:* FCM (\#151) $->_(m -> 1^+)$ Lloyd $k$-Means (Alg. 1, @lloyd1982least); GK (\#153, $A_j = I_d$) $equiv$ FCM $equiv$ $"sFCM"_(1, 0)$ (\#156); GG (\#154, $m = 2$) $equiv$ GMM E-step (Alg. 81, @dempster1977em); KFCM-F (\#155) $->_(m -> 1^+)$ Kernel $k$-Means / Ncut @dhillon2004kernel.
- *Biclustering Embeddings:* Zero-MSR Cheng--Church (\#157) $arrow.hook$ 1-layer PLAID (\#158); Bimax (\#159) $equiv$ binary zero-MSR maximal bicliques; unthresholded symmetric ISA (\#160) $equiv$ rank-1 SVD; Gaussian-prior FABIA (\#161) $equiv$ classical Factor Analysis.
- *Temporal Limits:* Zero-lag k-Shape (\#162) $equiv$ Spherical $k$-Means (Alg. 12, @dhillon2001concept); zero-band DBA (\#163) $equiv$ Lloyd $k$-Means; Soft-DTW (\#164) $->_(gamma -> 0^+)$ DBA; ROCKET (\#165) $->_(M -> infinity)$ multi-scale convolutional Kernel $k$-Means.
- *Spatial, Ensemble & Multi-View Bridges:* Unwindowed SLIC (\#166, $m -> 0$) $equiv$ 3D Lab $k$-Means; injective Watershed (\#167) $equiv$ minima-rooted MSF; ECE (\#168, $"MinPts" = 1$) $equiv$ 3D EMST cut, while DBSCAN-3D and TRACLUS (\#169) instantiate GDBSCAN @sander1998gdbscan; full-sample Monti (\#170, $f = 1$) $equiv$ CSPA Gram $S = r^(-1) H H^top$ (\#171), whose dual $H^top H$ drives MCLA (\#173); MVKM (\#174, $V = 1$) $equiv$ Lloyd $k$-Means; Co-regularized Spectral (\#175) reduces at $lambda = 0$ to NJW @ng2002spectral and at $lambda -> infinity$ to normalized graph-sum spectral clustering.

#subsec("B", "Mathematically Grounded Practitioner's Decision Guide")

1. *Soft Prototypes (Clean vs. Noisy Euclidean):* Default to *FCM (\#151)* ($m in [1.5, 2.0]$) on clean spherical data. Under gross outliers (Theorem 1.1), either run *PCM (\#152)* *warm-started from FCM* with locally calibrated $eta_j$ (avoiding coincident collapse, Theorem 1.3) or deploy *Gaussian KFCM-K (\#155)* to retain simplex memberships alongside redescending outlier immunity (Theorem 2.2).
2. *Ellipsoidal, Manifold & Raster Fuzzy Data:* Choose *GK (\#153)* (with ridge regularization) for equal-volume skewed ellipsoids; upgrade to warm-started *GG (\#154)* (with eigenvalue floors) only when unequal volumes/priors are essential (Theorem 2.1). Use *KFCM-F (\#155)* on non-convex shells ($n <= 10^4$), and *sFCM (\#156)* ($q <= 1$) on impulse-noisy raster imagery while guarding against Theorem 2.3's 1-pixel filament erosion.
3. *Two-Mode Matrix Biclustering:* For *additive* possibly overlapping blocks, prefer *PLAID (\#158)* over *Cheng--Church (\#157)* (avoiding Theorem 3.2's masking distortion). For *multiplicative* modules ($lambda_k z_k^top$), reject Cheng--Church (Theorem 3.1): choose *FABIA (\#161)* for joint identifiable factors, heteroscedastic noise calibration, and MI ranking (Proposition 1.2), or *ISA (\#160)* on very large matrices. Restrict *Bimax (\#159)* strictly to noise-free binary tables (Theorem 3.3).
4. *Time-Series Clustering:* Use *k-Shape (\#162)* for rigid phase/scale shifts ($m <= 2000$). Under nonlinear time warping (Theorem 4.1), select *DBA (\#163)* for short series ($m <= 500$) or *Soft-DTW divergence (\#164)* for smooth gradient barycenters escaping DBA kinks (Theorem 5.2). For long series ($m >= 10^3$) or large $n$, *ROCKET-Clustering (\#165)* dominates in $O(n M m)$ linear time with exact metricity (Theorem 5.3).
5. *Superpixels & Relief Segmentation:* Choose *SLIC (\#166)* for $O(N)$ compact superpixels with $m < Delta_c$ (Theorem 6.1). Choose *Watershed (\#167)* for contrast-invariant closed contours—*always* preceded by $h$-minima persistence filtering or markers to prevent zero-threshold over-segmentation!
6. *3D LiDAR & Moving Trajectories:* Use PCL *ECE (\#168)* on clean dense scans, switching to *DBSCAN-3D* ($"MinPts" >= 4$) against thin bridges and adopting range-adaptive $d_"th" (r) prop r$ (or HDBSCAN @campello2013hdbscan) beyond Theorem 6.2's critical fragmentation range $r_"crit"$. Use *TRACLUS (\#169)* whenever trajectories share local sub-corridors despite global divergence (Theorem 6.3).
7. *Stability & Cluster Ensembles:* Pair *Monti (\#170)* with Proportion of Ambiguous Clustering (PAC) rather than raw CDF area $A(K)$ (Theorem 7.3). For ensemble fusion, default to *MCLA (\#173)* ($O(n)$ with soft confidences), include *CSPA (\#171)* for moderate balanced $n <= 5000$, avoid standalone *HGPA (\#172)* under base noise, and run Strehl--Ghosh's *ANMI supra-consensus* selector whenever feasible (Theorem 7.1).
8. *Multi-View Fusion:* Use *MVKM (\#174)* ($O(T n k D)$) for large convex feature tables *strictly after Frobenius-normalizing every view* (eliminating Theorem 8.1's scale hijacking). Use *Co-regularized Spectral (\#175)* for non-convex manifolds/graphs ($n <= 5000$), auditing view eigengaps against Theorem 8.2's fixed-$lambda$ adversarial hijacking.

#subsec("C", "Open Theoretical Frontiers & Cross-Paradigm Bridges")

Four open frontiers bridge Category VII to deep clustering (Part VIII) and Chapters X--XIII:
1. *Coincidence-Free Affine-Invariant Robust Fuzzy Clustering:* Uniting $"GL"(d, bb(R))$ unequal-volume adaptivity (GG), redescending breakdown (PCM / Gaussian KFCM), and inter-cluster repulsion in one $O(n c d^2)$ variational objective.
2. *Hybrid Additive--Multiplicative Triclustering:* Establishing minimax recovery and polynomial-time solvers for 3-way multi-omics tensors mixing PLAID offsets, FABIA factors, and missing entries.
3. *Sub-Quadratic Differentiable Elastic Divergences:* Designing $O(m "poly"(log m))$ smooth Soft-DTW approximations with exact gradients, uniting explicit barycenters with ROCKET's linear scaling.
4. *Scale-Invariant Sketched Multi-View & ANMI Consensus:* Combining adaptive view weights (curing MVKM scale hijacking and fixed-$lambda$ spectral vulnerability) with $O(n)$ Grassmannian sketching and ANMI optimality.
