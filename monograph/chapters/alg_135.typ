#import "../template/ieee.typ": *

#alg_entry(
  id: 135,
  name: "WaveCluster (Wavelet-Based Grid Clustering)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1998",
  paradigm: "Multiresolution Discrete Wavelet Grid-Density Clustering",
  complexity: $cal(O)(n + K)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

WaveCluster @sheikholeslami1998wavecluster reformulates spatial clustering as a *multiresolution signal-processing problem*. Given a finite dataset $cal(X) = {bold(x)_1, dots.h, bold(x)_n} subset Omega = product_(j=1)^d [a_j, b_j] subset RR^d$, WaveCluster quantizes the point cloud onto a regular $d$-dimensional lattice, treats the cell counts as a discrete $d$-dimensional density signal, and applies a separable Discrete Wavelet Transform (DWT). Low-pass wavelet filtering suppresses isolated outlier spikes while reinforcing coherent dense interiors, negative-sidelobe biorthogonal filters sharpen inter-cluster valleys, and iterating the transform across scales $r = 1, dots.h, R$ produces a nested coarse-to-fine clustering hierarchy in time strictly linear in $n$.

Quantize each coordinate axis $j in {1, dots.h, d}$ into $m$ uniform intervals of width $Delta_j = (b_j - a_j) / m$ (taking $m = 2^p$ for dyadic convenience), defining the fine lattice $cal(G)_0 = {0, dots.h, m - 1}^d$ of $K = m^d$ cells. Each observation $bold(x)_i$ maps to grid index $bold(u)(bold(x)_i) in cal(G)_0$ via $u_j = min(m - 1, floor((x_(i,j) - a_j) / Delta_j))$, yielding the discrete grid-density signal $f in ell^2 (cal(G)_0)$:
$ f[bold(u)] = |{bold(x)_i in cal(X) : bold(u)(bold(x)_i) = bold(u)}|. $

#definition(title: "Separable Multiresolution DWT and Scale-r Wavelet Clusters")[
  Let $(h, g)$ be a compactly supported 1D low-pass / high-pass wavelet filter pair (such as Haar, Daubechies, or biorthogonal Cohen-Daubechies-Feauveau CDF(2,2) hat-shaped filters @sheikholeslami1998wavecluster), with unitary DC gain $sum_k h_k = sqrt(2)$, balanced polyphase sums $sum_m h_(2m) = sum_m h_(2m+1) = 1 / sqrt(2)$, and vanishing-moment high-pass sum $sum_k g_k = 0$. Define the 1D analysis operators along coordinate axis $j$ by convolution followed by dyadic downsampling:
  $ (cal(L)_j A)[dots.h, v_j, dots.h] = sum_k h_k A[dots.h, 2 v_j - k, dots.h], $
  $ (cal(H)_j A)[dots.h, v_j, dots.h] = sum_k g_k A[dots.h, 2 v_j - k, dots.h]. $
  Initializing $A_0 := f$, the *scale-$r$ approximation subband* ($r in {1, dots.h, R}$) on the coarsened lattice $cal(G)_r = {0, dots.h, m / 2^r - 1}^d$ is generated recursively via separable application across all $d$ axes:
  $ A_r = (cal(L)_d compose dots.c compose cal(L)_1) A_(r-1). $
  Given a density threshold $tau_r > 0$, the *active wavelet cells* are $cal(V)_r = {bold(v) in cal(G)_r : A_r [bold(v)] > tau_r}$. The *scale-$r$ WaveCluster partition* consists of the maximal connected components $C_1^((r)), dots.h, C_(M_r)^((r))$ of $cal(V)_r$ under $2d$- or $(3^d - 1)$-adjacency, mapped back to fine cells $bold(u) in cal(G)_0$ and observations $bold(x)_i in cal(X)$ via the dyadic lookup map $bold(v) = floor(bold(u) / 2^r)$.
]

Because $A_r$ aggregates spatial neighborhoods through localized multiresolution scaling functions, isolated outliers falling inside sparse grid regions yield sub-threshold responses and are automatically discarded, whereas dense cluster cores receive coherent amplitude reinforcement across all scales $r = 1, dots.h, R$.

#algbox(
  title: "WaveCluster Multiresolution Grid-Density Clustering",
  inputs: [Dataset $cal(X) subset Omega subset RR^d$, grid size $m$ per axis ($K = m^d$), wavelet filter $h$, target scale $r >= 1$, threshold $tau_r$.],
  outputs: [Scale-$r$ cluster labels $ell(bold(x)_i) in {0, 1, dots.h, M_r}$ (with $0$ denoting noise).],
)[
  + *Spatial Quantization:* Initialize zero array $A_0$ of size $m^d$. In one pass over $cal(X)$, hash each $bold(x)_i$ to $bold(u)_i = bold(u)(bold(x)_i)$ and increment $A_0 [bold(u)_i] <- A_0 [bold(u)_i] + 1$.
  + *Separable Pyramidal DWT:* For scale $s = 1$ to $r$:
    - Set $B <- A_(s-1)$.
    - For axis $j = 1$ to $d$: apply 1D convolution with filter $h$ along axis $j$ and downsample by $2$: $B <- cal(L)_j B$.
    - Set $A_s <- B$ (optionally retaining detail subbands $cal(H)_j$ for multi-scale edge inspection).
  + *Thresholding and Connected Components:* Extract active coarse cells $cal(V)_r = {bold(v) in cal(G)_r : A_r [bold(v)] > tau_r}$ and run linear-time BFS/DFS connected-component labeling over $cal(V)_r$ to assign component IDs $"cid"[bold(v)] in {1, dots.h, M_r}$.
  + *Dyadic Inverse Lookup:* Build lookup table $T[bold(u)] = "cid"[floor(bold(u) / 2^r)]$ if $floor(bold(u) / 2^r) in cal(V)_r$, else $0$.
  + *Point Assignment:* In one final pass over $cal(X)$, assign $ell(bold(x)_i) <- T[bold(u)(bold(x)_i)]$ and return labels.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four formal theorems characterizing WaveCluster's mass conservation and outlier SNR amplification, the spectral overshoot and valley-carving bridge separation of biorthogonal hat-shaped filters, exact linear pyramidal complexity, and shift-variance / Nyquist resolution barriers.

#subsec("A", "Mass Conservation and Impulse-Noise SNR Amplification")

#theorem(title: "Exact Mass Conservation and Outlier SNR Gain", origin: "Original Proof")[
  Let $h$ be a compactly supported 1D low-pass wavelet scaling filter with $sum_k h_k = sqrt(2)$ and polyphase balance $sum_m h_(2m) = sum_m h_(2m+1) = 1 / sqrt(2)$, and let $h^((r))$ denote the equivalent $r$-scale 1D cascade filter satisfying $(cal(L)^r a)[v] = sum_u h^((r))_(2^r v - u) a[u]$. On a periodic (or zero-padded infinite) lattice:
  + *Mass Conservation:* At every scale $r >= 1$, total approximation sum scales exactly as
    $ sum_(bold(v) in cal(G)_r) A_r [bold(v)] = 2^(-r d / 2) sum_(bold(u) in cal(G)_0) f[bold(u)] = 2^(-r d / 2) n. $
  + *Impulse-Noise SNR Amplification:* Suppose $f$ contains a flat cluster plateau of density $rho_c > 0$ wider than $2^r |"supp"(h)|$ per axis and an isolated outlier impulse of mass $eta > 0$ at $bold(u)_0$. Then inside the plateau interior $A_r [bold(v)] = 2^(r d / 2) rho_c$, while the outlier peak obeys $max_bold(v) |A_r^"out" [bold(v)]| <= eta norm(h^((r)))_infinity^d$, yielding exact Signal-to-Noise Ratio ("SNR") gain
    $ "SNR"_r / "SNR"_0 >= (2^(r / 2) / norm(h^((r)))_infinity)^d. $
    In particular, every orthogonal wavelet ($norm(h^((r)))_2 = 1$) guarantees $"SNR"_r / "SNR"_0 >= 2^(r d / 2)$, and Haar achieves exact gain $2^(r d)$.
]

#proof[
  For (1), evaluate the sum of a single 1D analysis pass $b[v] = sum_k h_k a[2 v - k]$ by substituting $u = 2 v - k$ and splitting by parity of $u$:
  $ sum_v b[v] &= sum_u a[u] sum_(v) h_(2 v - u) \
    &= sum_(u "even") a[u] sum_m h_(2 m) + sum_(u "odd") a[u] sum_m h_(2 m + 1) \
    &= 2^(-1 / 2) sum_u a[u]. $
  Applying this identity separably across all $d$ coordinate axes multiplies the total grid sum by $2^(-d / 2)$ per scale; iterating over $r$ scales gives $sum_(bold(v)) A_r [bold(v)] = 2^(-r d / 2) n$.

  For (2), on a constant plateau $f[bold(u)] equiv rho_c$ over the filter support, each 1D pass multiplies the amplitude by $sum_k h_k = sqrt(2)$, so $d$ axes across $r$ scales yield $A_r [bold(v)] = (sqrt(2))^(r d) rho_c = 2^(r d / 2) rho_c$. By separability, the response to a single impulse $f^"out" [bold(u)] = eta delta_(bold(u), bold(u)_0)$ is
  $ A_r^"out" [bold(v)] = eta product_(j=1)^d h^((r))_(2^r v_j - u_(0,j)), $
  whose maximum magnitude across all $bold(v) in cal(G)_r$ is bounded by $eta norm(h^((r)))_infinity^d$. Taking the ratio against the input $"SNR"_0 = rho_c / eta$ yields
  $ "SNR"_r / "SNR"_0 = (2^(r d / 2) rho_c / max_bold(v) |A_r^"out" [bold(v)]|) / (rho_c / eta) >= (2^(r / 2) / norm(h^((r)))_infinity)^d. $
  For any orthogonal wavelet, the $r$-scale cascade filter $h^((r))$ forms an orthonormal basis shift vector with $norm(h^((r)))_2 = 1$, hence $norm(h^((r)))_infinity <= norm(h^((r)))_2 = 1$, yielding $"SNR"_r / "SNR"_0 >= 2^(r d / 2)$. For the Haar filter, $h^((r))_k = 2^(-r / 2)$ on ${0, dots.h, 2^r - 1}$, so $norm(h^((r)))_infinity = 2^(-r / 2)$ and $"SNR"_r / "SNR"_0 = 2^(r d)$.
]

#subsec("B", "Spectral Overshoot and Valley-Carving Bridge Separation of CDF(2,2)")

A signature design choice of WaveCluster @sheikholeslami1998wavecluster is preferring symmetric biorthogonal hat-shaped filters — specifically Cohen-Daubechies-Feauveau CDF(2,2) with taps $h = sqrt(2)(-1/8, 1/4, 3/4, 1/4, -1/8)$ — over simple non-negative box/Haar averaging. We prove analytically why CDF(2,2) simultaneously enhances cluster cores and severs low-density bridges.

#theorem(title: "Zero-Phase Passband Boost and Valley-Carving Bridge Disconnection", origin: "Original Derivation")[
  Let $h = sqrt(2)(-1/8, 1/4, 3/4, 1/4, -1/8)$ (indexed on $k in {-2, -1, 0, 1, 2}$) be the 1D biorthogonal CDF(2,2) scaling filter.
  + *Spectral Characterization:* Its discrete-time frequency response $H(omega) = sum_(k=-2)^2 h_k e^(-i k omega)$ is purely real (zero phase distortion) and satisfies the exact trigonometric factorization
    $ H(omega) = sqrt(2) cos^2 (omega / 2) (1 + 2 sin^2 (omega / 2)), $
    which has a double zero at Nyquist ($H(pi) = H'(pi) = 0$) yet exhibits a mid-band cluster-enhancing gain peak $H(pi / 3) = (9 / 8) sqrt(2) > H(0)$.
  + *Bridge Disconnection:* Consider two 1D cluster plateaus of density $rho_c > 0$ on $u <= -2$ and $u >= 2$ joined by a 3-cell low-density bridge of height $rho_b in (0, rho_c)$ on $u in {-1, 0, 1\}$. At coarse cell $v = 0$, CDF(2,2) produces
    $ A_1^"CDF" [0] = sqrt(2) (5 / 4 rho_b - 1 / 4 rho_c) < sqrt(2) rho_b, $
    dropping strictly below zero whenever $rho_b <= rho_c / 5$. Thus for every threshold $tau_1 in (A_1^"CDF" [0], sqrt(2) rho_b)$, CDF(2,2) severs the bridge and recovers both clusters separately, whereas Haar and every non-negative low-pass filter of equal DC gain merge them into one cluster.
]

#proof[
  For (1), symmetry $h_(-k) = h_k$ eliminates imaginary sine terms and gives
  $ H(omega) / sqrt(2) &= 3 / 4 + 1 / 2 cos(omega) - 1 / 4 cos(2 omega) \
    &= 3 / 4 + 1 / 2 cos(omega) - 1 / 4 (2 cos^2 (omega) - 1) \
    &= 1 + 1 / 2 cos(omega) - 1 / 2 cos^2 (omega) \
    &= 1 / 2 (1 + cos(omega)) (2 - cos(omega)). $
  Substituting the half-angle identities $1 + cos(omega) = 2 cos^2 (omega / 2)$ and $2 - cos(omega) = 1 + 2 sin^2 (omega / 2)$ yields $H(omega) = sqrt(2) cos^2 (omega / 2) (1 + 2 sin^2 (omega / 2))$. Setting $s = sin^2 (omega / 2) in [0, 1]$ expresses $H(omega) / sqrt(2) = (1 - s)(1 + 2 s) = 1 + s - 2 s^2$. At $omega = pi$ ($s = 1$), the factor $cos^2 (omega / 2)$ produces a double zero $H(pi) = H'(pi) = 0$, strongly suppressing high-frequency checkerboard noise. Differentiating $1 + s - 2 s^2$ gives $1 - 4 s = 0 => s^* = 1 / 4$ ($omega^* = pi / 3$), where $H(pi / 3) = sqrt(2)(1 + 1 / 4 - 1 / 8) = (9 / 8) sqrt(2)$, proving mid-band amplification.

  For (2), evaluating $(cal(L) f)[0] = sum_(k=-2)^2 h_k f[-k]$ on the bridge profile ($f[-2] = f[2] = rho_c$, $f[-1] = f[0] = f[1] = rho_b$) yields
  $ A_1^"CDF" [0] &= sqrt(2) (-1 / 8 rho_c + (1 / 4 + 3 / 4 + 1 / 4) rho_b - 1 / 8 rho_c) \
    &= sqrt(2) (5 / 4 rho_b - 1 / 4 rho_c) = sqrt(2) rho_b - sqrt(2) / 4 (rho_c - rho_b). $
  Because $rho_c > rho_b$, $A_1^"CDF" [0] < sqrt(2) rho_b$ strictly, and $A_1^"CDF" [0] <= 0$ whenever $rho_b <= rho_c / 5$. In contrast, any non-negative filter $tilde(h) >= 0$ with $sum_k tilde(h)_k = sqrt(2)$ satisfies $(tilde(cal(L)) f)[0] = sum_k tilde(h)_k f[-k] >= sum_k tilde(h)_k rho_b = sqrt(2) rho_b$ (since $f[u] >= rho_b$ everywhere on the support). Hence any threshold $tau_1 in (A_1^"CDF" [0], sqrt(2) rho_b)$ strictly excludes $v = 0$ under CDF(2,2) while retaining $v = 0$ under every non-negative filter, disconnecting the two plateaus cleanly.
]

#subsec("C", "Linear Pyramidal Complexity and Order Invariance")

#theorem(title: "Linear-in-n Complexity via Mallat's Separable Pyramid", origin: "Classical — Cited")[
  For a $d$-dimensional grid of $K = m^d$ cells and a compactly supported wavelet filter of length $M = |"supp"(h)|$, executing WaveCluster up to any scale $R >= 1$ takes exact total time @sheikholeslami1998wavecluster
  $ T_"Wave" = cal(O)(n d + d M K) = cal(O)(n + K) quad ("for fixed" d, M), $
  with space $cal(O)(K)$, independent of the number or shape of clusters and strictly invariant to input point ordering.
]

#proof[
  Step 1 hashes each of the $n$ points across $d$ axes in $cal(O)(n d)$ time to form $A_0$, accumulating commutative integer counts invariant to any permutation $pi in S_n$. In Step 2, at scale $s in {1, dots.h, R}$, the input tensor $A_(s-1)$ contains $K / 2^((s-1) d)$ entries. Applying the separable 1D length-$M$ convolution along axis $j = 1$ reads $K / 2^((s-1) d)$ cells and writes $K / (2 dot.op 2^((s-1) d))$ downsampled cells; across all $d$ axes at scale $s$, arithmetic operations are bounded by $d M K / 2^((s-1) d)$. Summing the geometric series over $s = 1, dots.h, R$ yields
  $ T_"DWT" <= d M K sum_(s=1)^R (2^(-d))^(s - 1) < (2^d) / (2^d - 1) d M K <= 2 d M K. $
  Step 3 runs linear-time BFS/DFS over $cal(G)_r$ ($K / 2^(r d)$ cells, degree $d_"adj" <= 3^d - 1$) in $cal(O)(d_"adj" K / 2^(r d))$ time; Step 4 builds $T$ in $cal(O)(K)$ time (or computes $floor(bold(u) / 2^r)$ on the fly in $cal(O)(d)$); Step 5 assigns all $n$ points in $cal(O)(n d)$ time. Summing all stages yields $cal(O)(n d + d M K)$ time and $cal(O)(K)$ space.
]

#subsec("D", "Shift-Variance Aliasing and Dimensionality Barriers")

#theorem(title: "Dyadic Shift-Variance Counterexample and Nyquist Resolution Bound", origin: "Original Counterexample")[
  Because the critically sampled DWT applies dyadic downsampling ($arrow.b 2$), WaveCluster is *not translation-invariant* at sub-coarse scales and faces both a Nyquist separation limit and an exponential dimensionality wall:
  + *Single-Cell Shift Peak Halving:* Under scale-1 Haar WaveCluster, translating a 2-cell 1D cluster of density $rho$ by a single fine cell ($Delta$) halves its peak approximation coefficient from $sqrt(2) rho$ to $rho / sqrt(2)$, toggling cluster existence for every threshold $tau_1 in (rho / sqrt(2), sqrt(2) rho]$.
  + *Nyquist Gap Aliasing:* Two clusters separated by a gap of width $gamma < 2^r Delta$ along an axis can map into the same coarse cell at scale $r$, guaranteeing false merging regardless of $tau_r$.
  + *Curse of Dimensionality:* Supporting $R$ dyadic scales with at least $2$ coarse cells per axis requires $m >= 2^(R + 1)$, imposing a strict memory floor $K = m^d >= 2^((R + 1) d)$.
]

#proof[
  For (1), define two 1D grid signals differing by a single-cell shift: $f[u] = rho$ on $u in {2 t, 2 t + 1}$ (zero elsewhere), and $tilde(f)[u] = f[u - 1] = rho$ on $u in {2 t + 1, 2 t + 2}$. Under Haar analysis $A_1 [v] = (f[2 v] + f[2 v + 1]) / sqrt(2)$:
  $ A_1 [t] = (rho + rho) / sqrt(2) = sqrt(2) rho, quad A_1 [v] = 0 quad (v != t), $
  whereas the 1-cell shifted signal splits equally across two coarse cells:
  $ tilde(A)_1 [t] = (0 + rho) / sqrt(2) = rho / sqrt(2), quad tilde(A)_1 [t + 1] = (rho + 0) / sqrt(2) = rho / sqrt(2). $
  Thus the peak amplitude drops by an exact factor of $2$, erasing the cluster completely whenever $tau_1 in (rho / sqrt(2), sqrt(2) rho]$.

  For (2), each coarse cell $v in cal(G)_r$ aggregates a fine block of width $2^r Delta$. If two disjoint clusters occupy cells $u_1 = 2^r v$ and $u_2 = 2^r v + 2^r - 1$ separated by an empty gap of $2^r - 2$ fine cells (width $(2^r - 2) Delta < 2^r Delta$), both clusters map to the identical coarse index $v = floor(u_1 / 2^r) = floor(u_2 / 2^r)$, making separation impossible at scale $r$.

  For (3), resolving even two non-trivial cells per axis after $r = R$ dyadic halvings requires $m / 2^R >= 2 => m >= 2^(R + 1)$, so $K = m^d >= 2^((R + 1) d)$ (exceeding $10^9$ cells by $d = 10, R = 2$), proving WaveCluster is viable only in low-to-moderate ambient dimensions.
]

#sec("III", "Complexity, Geometry and Operational Profile")

WaveCluster achieves one of the lowest asymptotic costs in all of clustering: $cal(O)(n d + d M K)$ time and $cal(O)(K)$ space, strictly linear in dataset cardinality $n$ for fixed grid resolution $K = m^d$. Because the separable DWT sweeps contiguously across memory axes (analogous to 2D/3D image convolutions), it achieves exceptional SIMD/GPU throughput and cache bandwidth. Crucially, multi-scale exploration across scales $r = 1, dots.h, R$ requires *zero re-scanning of $cal(X)$* — each coarser scale shrinks the active array by factor $2^d$, costing a tiny geometric fraction of scale 1.

Geometrically, WaveCluster discovers *arbitrary non-convex, nested, and elongated manifolds* in 2D/3D without specifying the cluster count $k$, while simultaneously performing automatic noise suppression (Theorem II.A) and valley carving (Theorem II.B). Its main operational constraints are tuning the grid size $m$, scale $r$, wavelet family $h$, and threshold $tau_r$, along with the $m^d$ memory ceiling for $d >= 5$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(n d + d M K)$: strictly linear in $n$, independent of cluster count $k$],
    [Space Complexity], [$cal(O)(K) = cal(O)(m^d)$ dense lattice tensor + $cal(O)(K)$ lookup map],
    [Multiresolution], [Native nested coarse-to-fine hierarchy across scales $r = 1, dots.h, R$],
    [Cluster Geometry], [Arbitrary non-convex, curved, ring, and multi-connected low-D shapes],
    [Outlier Rejection], [Inherent low-pass SNR gain $>= 2^(r d / 2)$ (orthogonal) or $2^(r d)$ (Haar)],
    [Order / Streaming], [Strictly permutation-invariant; single-pass additive grid accumulation],
    [Primary Hyperparams], [Grid resolution $m$, wavelet filter $h$, scale $r$, density threshold $tau_r$],
    [Key Limitations], [$m^d$ memory barrier for $d >= 5$; dyadic shift variance; Nyquist gap limit],
  ),
  caption: [Computational, geometric, and operational profile of WaveCluster @sheikholeslami1998wavecluster.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Medical and Microscopy Image Segmentation:* Fluorescence microscopy, histopathology slides, and MRI/CT volumetric scans natively inhabit 2D/3D Cartesian grids corrupted by shot noise. WaveCluster's CDF(2,2) hat-shaped filter suppresses single-pixel speckle noise (Theorem II.A) while carving valleys between touching cellular nuclei (Theorem II.B) in real time @sheikholeslami1998wavecluster.

*2. Geospatial Hot-Spot and Urban Mobility Mapping:* Billions of GPS trajectories, rideshare pickups, and seismic epicenters project onto 2D latitude-longitude lattices. Multiresolution wavelet scales $r = 1, 2, 3$ simultaneously reveal neighborhood-level hotspots, city-district hubs, and metropolitan corridors from a single linear scan of $cal(X)$ without quadratic spatial joins.

*3. LiDAR and 3D Point-Cloud Object Extraction:* Autonomous-vehicle LiDAR frames produce $10^5$–$10^6$ $(x, y, z)$ returns per sweep with stray atmospheric/rain returns. Voxelizing into a 3D lattice and applying a 1-scale DWT attenuates isolated rain speckles by up to $8 times$ (Haar) while extracting arbitrary non-convex vehicle, pedestrian, and vegetation contours.

*4. Astronomical Multi-Scale Structure and Filament Detection:* Cosmological N-body simulations and galaxy surveys exhibit hierarchical clustering across orders of magnitude in spatial scale (clusters, superclusters, filaments). Wavelet multiresolution naturally matches astrophysical structure scales, separating compact halos at fine $r$ from large-scale superclusters at coarse $r$.

*5. High-Throughput Flow Cytometry and Bi-Axial Gating:* Immunology assays record millions of single-cell fluorescence intensities projected onto 2D/3D marker scatter plots containing skewed, non-ellipsoidal cell populations and sparse debris. WaveCluster's $cal(O)(n)$ speed and inherent outlier suppression automate multi-population gating far faster than iterative KDE or DBSCAN @ester1996dbscan.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

WaveCluster fails in four mathematical regimes: (i) *high-dimensional data ($d >= 5$)* where $K = m^d$ exhausts memory or turns sparse (Theorem II.D(3)); (ii) *tightly adjacent micro-clusters* separated by gaps $< 2^r Delta$ that alias together under downsampling (Theorem II.D(2)); (iii) *extreme multi-density scenes*, since a single global threshold $tau_r$ cannot isolate a dense cluster and a $100 times$ sparser cluster simultaneously; and (iv) *boundary-critical tasks* sensitive to dyadic grid alignment (mitigated by using the undecimated / stationary wavelet transform at $cal(O)(R d M K)$ cost).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [2D/3D Spatial & GIS Data], [Linear $cal(O)(n + K)$; non-convex shapes + multi-scale hierarchy], [Shift aliasing; use stationary DWT or OPTICS @ankerst1999optics],
    [Noisy Microscopy / LiDAR], [Automatic impulse-noise SNR gain + CDF(2,2) valley carving], [Tight gaps $< 2^r Delta$ merge; reduce scale $r$],
    [Multi-Density Point Clouds], [Moderate: global $tau_r$ struggles with wide density ratios], [Use OPTICS @ankerst1999optics or DENCLUE @hinneburg1998denclue],
    [High-D Tabular / Omics ($d >= 6$)], [Poor: $K = m^d$ tensor explosion and zero occupancy], [Use PROCLUS @aggarwal1999proclus or CLIQUE @agrawal1998clique],
    [Evolving Data Streams], [Single-pass additive $A_0$; fast periodic DWT refresh], [For temporal decay use D-Stream @chen2007dstream],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for WaveCluster.],
)

#sec("V", "Structural Relations and Empirical Behavior")

WaveCluster unifies grid-based summarization (STING @wang1997sting) with kernel density level-set clustering (DENCLUE @hinneburg1998denclue, Mean Shift @ariascastro2016meanshift, and DBSCAN @ester1996dbscan).

#proposition(title: "Exact Equivalence to STING (Haar) and B-Spline KDE Level-Set Limit", origin: "Original Proof")[
  Let $A_r$ be WaveCluster's scale-$r$ approximation subband on $cal(G)_r$ with mesh $Delta$ and filter $h$.
  + *Exact STING Reduction (Haar Filter):* When $h = (1 / sqrt(2), 1 / sqrt(2))$ (Haar), $A_r [bold(v)] = 2^(-r d / 2) n_(c(bold(v)))$, where $n_(c(bold(v)))$ is the exact point count of the level-$(L - r)$ dyadic STING cell @wang1997sting covering block $bold(v)$. Thus Haar WaveCluster with threshold $tau_r$ produces the *exact same* cluster partition as STING at tree level $L - r$ with count threshold $2^(r d / 2) tau_r$.
  + *Continuous KDE Level-Set Limit:* As fine mesh $Delta -> 0$ with physical scale $lambda = 2^r Delta$ held fixed, the normalized approximation density $hat(p)_lambda (bold(x)) = (n 2^(r d / 2) Delta^d)^(-1) A_r [floor(bold(x) / lambda)]$ converges to the product-kernel density estimator
    $ hat(f)_lambda (bold(x)) = 1 / (n lambda^d) sum_(i=1)^n Phi((bold(x) - bold(x)_i) / lambda), $
    where $Phi(bold(z)) = product_(j=1)^d phi(z_j)$ is the continuous scaling function generated by the dilation equation $phi(t) = sqrt(2) sum_k h_k phi(2 t - k)$. Hence WaveCluster computes the connected superlevel sets ${bold(x) : hat(f)_lambda (bold(x)) > tau}$ of a fast $cal(O)(n + K)$ discretized KDE @hinneburg1998denclue @ariascastro2016meanshift.
]

#proof[
  For (1), iterating the 1D Haar operator $r$ times sums $2^r$ consecutive fine bins and multiplies by $2^(-r / 2)$. Across $d$ axes, $A_r [bold(v)] = 2^(-r d / 2) sum_(bold(u) : floor(bold(u) / 2^r) = bold(v)) f[bold(u)] = 2^(-r d / 2) n_(c(bold(v)))$. Inequality $A_r [bold(v)] > tau_r$ is therefore algebraically identical to $n_(c(bold(v))) > 2^(r d / 2) tau_r$, yielding identical active vertices and BFS components.

  For (2), by the cascade algorithm of Mallat's multiresolution analysis, the scaled equivalent impulse response $2^(r / 2) h^((r))_(floor(2^r t))$ converges in $L^2 (RR)$ as $r -> infinity$ to the continuous scaling function $phi(t)$ satisfying $phi(t) = sqrt(2) sum_k h_k phi(2 t - k)$ with $integral phi(t) dif t = 1$. Substituting $A_r [bold(v)] = sum_(i=1)^n product_(j=1)^d h^((r))_(2^r v_j - u_(i,j))$ into $hat(p)_lambda (bold(x))$ and letting $Delta -> 0$ with $lambda = 2^r Delta$ fixed yields
  $ hat(p)_lambda (bold(x)) = 1 / (n lambda^d) sum_(i=1)^n product_(j=1)^d (2^(r / 2) h^((r))_(2^r v_j - u_(i,j))) -> 1 / (n lambda^d) sum_(i=1)^n Phi((bold(x) - bold(x)_i) / lambda), $
  identifying WaveCluster's active components directly with the connected superlevel sets of a compactly supported KDE (converging to a Gaussian kernel as B-spline order increases), bridging grid clustering with DENCLUE @hinneburg1998denclue, Grid-DBSCAN @mahran2008gridbscan, BANG @schikuta1998bang, and OPTIGRID @hinneburg1999optigrid @xu2005survey @aggarwal2013data @jain2010data.
]

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, WaveCluster records *Iris* $"Sil" = 0.517$, $"ARI" = 0.527$; *Customers* $"Sil" = 0.353$; *Multi-Shape* $"Sil" = 0.593$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.257$ ($"ARI" = "N/A"$), with a mean runtime of *$0.058$ ms*. These empirical metrics provide a striking validation of Proposition V and Section II:
- *Exact STING Equivalence on Three Benchmarks:* WaveCluster matches STING to three decimal places on *Iris* ($0.517 / 0.527$), *Multi-Shape* ($0.593$, tied best across Algs. 133–136), and *High-D* ($0.257$), directly reflecting Proposition V(1) whenever low-scale wavelet approximation preserves the same active grid superlevel components as hierarchical cell aggregation.
- *Filtering Divergence on Customers:* On *Customers*, WaveCluster scores $"Sil" = 0.353$ (trailing STING's $0.449$ yet doubling CLIQUE's $0.176$ and beating PROCLUS's $0.273$), illustrating how wavelet spatial smoothing alters borderline cell thresholding on skewed continuous tabular densities.
- *Near-Zero DWT Overhead:* Its $0.058$ ms runtime adds only $0.006$ ms over STING ($0.052$ ms), confirming Theorem II.C's $cal(O)(n + K)$ linear efficiency while outrunning PROCLUS ($4.477$ ms) by $77 times$.
