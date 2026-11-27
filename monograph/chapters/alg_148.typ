#import "../template/ieee.typ": *

#alg_entry(
  id: 148,
  name: "P3C (Projected Clustering via Cluster Cores)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2006",
  paradigm: "Poisson-Significant 1D Projections & Fuzzy EM Refinement",
  complexity: $O(d (n + log^2 n) + abs(cal(K)_"cand") n + T_"EM" n K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Most subspace and projected clustering methods suffer from brittle, unintuitive parameterization: PROCLUS @aggarwal1999proclus requires pre-specifying both the exact cluster count $K$ and the average subspace dimensionality $l$, while density-based methods like SUBCLU @kailing2004subclu and CLIQUE @agrawal1998clique demand arbitrary absolute density thresholds that fail to distinguish genuine clusters from random chance fluctuations. Moise, Sander, and Ester @moise2006p3c introduced *P3C* (Projected Clustering via Cluster Cores) to replace heuristic density cutoffs with formal *statistical significance testing*: P3C extracts statistically non-uniform 1D intervals via iterative chi-square histogram marking, combines them into multi-dimensional *cluster cores* whose observed point support exceeds conditional independence under a *Poisson tail test*, and refines the resulting cores into a full projected clustering with automatic outlier rejection via fuzzy Expectation-Maximization (EM) @dempster1977em.

Let $cal(D) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ with coordinate bounds $[a_k, b_k]$ along each attribute $k in {1, dots, d}$. Following Sturges' histogram rule, P3C partitions each range $[a_k, b_k]$ into $beta = ceil(1 + log_2 n)$ equal-width bins $B_(k, 1), dots, B_(k, beta)$ with observed bin counts $c_(k, 1), dots, c_(k, beta)$ summing to $n$.

#definition(title: "Sturges-Marked 1D Intervals and Poisson-Significant Cluster Cores")[
  1. *1D Iterative Chi-Square Marking:* Initialize the unmarked bin set $U_k <- {1, dots, beta}$. At each step, let $m = abs(U_k)$, $N_"rem" = sum_(b in U_k) c_(k, b)$, expected uniform count $E_k = N_"rem" / m$, and Pearson statistic
     $ chi_k^2 = sum_(b in U_k) (c_(k, b) - E_k)^2 / E_k. $
     While $m >= 2$ and $"Pr"(chi^2_(m - 1) >= chi_k^2) < alpha_"mark"$, mark $b^* = arg max_(b in U_k) c_(k, b)$ as dense and set $U_k <- U_k without {b^*}$. Adjacent marked bins merge into maximal contiguous 1D dense intervals $cal(I)_k = {I_(k, 1), dots, I_(k, m_k)}$ with supports $S(I) = {bold(x) in cal(D) : x_k in I}$.
  2. *Poisson-Significant $q$-Core:* A $q$-dimensional rectangular box $R = I_(k_1) times dots times I_(k_q)$ ($k_1 < dots < k_q$) has support $C(R) = inter_(r=1)^q S(I_(k_r))$ and observed count $O(R) = abs(C(R))$. For any $(q - 1)$-dimensional parent $R' = R without {k_r}$ and extending 1D interval $I_(k_r)$, the conditional-independence expected count is
     $ E(R', I_(k_r)) = abs(C(R')) dot.op ( abs(S(I_(k_r))) / n ). $
     Modeling $O(R) tilde.op "Pois"(lambda_R)$ with $lambda_R = E(R', I_(k_r))$, the core $R$ is *$alpha_P$-significant* if its upper Poisson tail satisfies
     $ p_"Pois" (O(R); lambda_R) = sum_(j = O(R))^infinity e^(-lambda_R) lambda_R^j / j! < alpha_P $
     (requiring $O(R) > lambda_R$ across parent splits). A *cluster core* is a maximal significant box $R$ having no significant $(q + 1)$-dimensional extension.
]

The set of maximal cluster cores $cal(K) = {R_1, dots, R_K}$ automatically determines *both* the number of clusters $K$ and initial subspace prototypes without user guesswork. In Phase 3, P3C initializes a $K$-component Gaussian mixture from $cal(K)$, runs fuzzy EM @dempster1977em to compute posterior memberships $bold(U) in [0, 1]^(n times K)$, flags points whose squared Mahalanobis distance exceeds a $chi^2$ critical value as outliers, and assigns each cluster $C_c$ its relevant attributes $S_c subset.eq {1, dots, d}$ via 1D statistical departure from uniformity.

#algbox(
  title: "P3C: Projected Clustering via Poisson-Significant Cluster Cores",
  inputs: [Dataset $cal(D) subset bb(R)^d$; 1D marking level $alpha_"mark"$; Poisson significance threshold $alpha_P$ (e.g., $10^(-5)$); outlier Mahalanobis cutoff $alpha_"out"$.],
  outputs: [Projected clusters ${(C_1, S_1), dots, (C_K, S_K)}$ and outlier set $O$.],
)[
  1. *Phase 1 — 1D Dense Interval Discovery:* Partition each coordinate $k = 1, dots, d$ into $beta = ceil(1 + log_2 n)$ equal-width bins. Iteratively apply the $chi^2_(m-1)$ uniformity test at level $alpha_"mark"$, removing the highest-count bin until uniformity holds on $U_k$, and merge adjacent marked bins into intervals $cal(I)_k$.
  2. *Phase 2 — Apriori Poisson Core Assembly:* Initialize level-1 cores $cal(L)_1 = union_(k=1)^d cal(I)_k$. For $q = 2, 3, dots, d$ while $cal(L)_(q-1) != emptyset$:
     - Generate candidate $q$-boxes $R = R' times I_j$ from $R' in cal(L)_(q-1)$ and $I_j in cal(I)_j$ ($j > max(S_(R'))$).
     - Evaluate observed support $O(R) = abs(C(R))$ and expected count $lambda_R = abs(C(R')) abs(S(I_j)) / n$; retain $R in cal(L)_q$ iff $O(R) > lambda_R$ and $p_"Pois" (O(R); lambda_R) < alpha_P$.
     - Retain maximal cores $cal(K) = {R_1, dots, R_K}$ not strictly contained in any higher-dimensional significant core (merging heavily overlapping cores).
  3. *Phase 3 — Fuzzy EM Refinement & Subspace Selection:*
     - Initialize fuzzy memberships $u_(i c)^((0))$ from core memberships $bb(1)[bold(x)_i in C(R_c)]$ and run EM updates for $(pi_c, bold(mu)_c, bold(Sigma)_c)$ until log-likelihood convergence.
     - Assign object $bold(x)_i$ to $c^* = arg max_c u_(i c)$ unless its Mahalanobis distance $D_M^2 (bold(x)_i, c^*) > chi^2_(d, 1 - alpha_"out")$, in which case place $bold(x)_i in O$.
     - For each cluster $C_c$, select relevant dimensions $S_c$ where $C_c$'s 1D marginal rejects uniformity (or spans marked intervals of $R_c$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish exact Le Cam approximation and Chernoff calibration bounds for P3C's Poisson core test, prove deterministic $O(log n)$ termination and ANOVA power bounds for Phase-1 chi-square marking, construct a 3D parity counterexample proving non-monotonicity of Poisson significance, and verify EM monotonicity and coordinate-scale invariance.

#subsec("A", "Le Cam Approximation and Poisson Chernoff Calibration")

#theorem(title: "Total-Variation Accuracy and Chernoff Inversion of the Poisson Core Test", origin: "Original Unified Proof")[
  Under conditional independence between a parent core $R'$ of size $m = abs(C(R'))$ and a 1D interval $I_j$ of marginal fraction $p = abs(S(I_j)) / n$, let $Y tilde.op "Bin"(m, p)$ be the exact with-replacement null count and $Z tilde.op "Pois"(lambda)$ with $lambda = m p$ be P3C's Poisson model. Then:
  (i) *Exact Total-Variation Error:*
  $ d_"TV" ("Bin"(m, p), "Pois"(lambda)) = 1 / 2 sum_(r=0)^infinity abs("Pr"(Y = r) - "Pr"(Z = r)) <= min(p, lambda p). $
  (ii) *Explicit Core-Detection Threshold:* For any relative excess factor $c > 1$, the upper Poisson tail obeys
  $ p_"Pois" (ceil(c lambda); lambda) <= exp( - lambda ( c log c - c + 1 ) ). $
  Consequently, any candidate core with observed count $O(R) >= c lambda$ is guaranteed to pass $p_"Pois" (O(R); lambda) < alpha_P$ whenever
  $ lambda >= log(1 / alpha_P) / (c log c - c + 1). $
]

#proof[
  (i) Express $Y = sum_(i=1)^m B_i$ with $B_i tilde.op "Bern"(p)$ i.i.d. and $Z = sum_(i=1)^m Z_i$ with $Z_i tilde.op "Pois"(p)$ i.i.d. By maximal coupling of a single coordinate $(B_i, Z_i)$, $"Pr"(B_i != Z_i) = (1 - p) - e^(-p) + (1 - e^(-p) - p e^(-p))$; specifically, since $"Pr"(B_i = 0) = 1 - p <= e^(-p) = "Pr"(Z_i = 0)$ and $"Pr"(B_i = 1) = p >= p e^(-p) = "Pr"(Z_i = 1)$, the single-term total variation is
  $ d_"TV" ("Bern"(p), "Pois"(p)) = p (1 - e^(-p)) <= p^2. $
  By subadditivity of total variation for sums of independent variables, $d_"TV" (Y, Z) <= sum_(i=1)^m d_"TV" (B_i, Z_i) <= m p^2 = lambda p$ (Le Cam's bound); the sharper Chen–Stein coupling refines the prefactor to $(1 - e^(-lambda)) p <= p$, giving $min(p, lambda p)$.
  (ii) For any $t > 0$, the moment generating function of $Z tilde.op "Pois"(lambda)$ is $bb(E)[e^(t Z)] = exp(lambda (e^t - 1))$. By Markov's inequality for any $k > lambda$,
  $ "Pr"(Z >= k) <= inf_(t > 0) exp( -t k + lambda (e^t - 1) ). $
  Differentiating $-t k + lambda (e^t - 1)$ with respect to $t$ yields the unique stationary minimizer $e^(t^*) = k / lambda > 1$, i.e., $t^* = log(k / lambda) > 0$. Substituting $t^*$ gives
  $ "Pr"(Z >= k) <= exp( -k log(k / lambda) + k - lambda ). $
  Setting $k = ceil(c lambda) >= c lambda$ and noting that $u |-> u log(u / lambda) - u + lambda$ is strictly increasing on $u > lambda$ (derivative $log(u / lambda) > 0$) yields $p_"Pois" (ceil(c lambda); lambda) <= exp(-lambda (c log c - c + 1))$. Since $h(c) = c log c - c + 1 > 0$ for all $c > 1$, requiring $exp(-lambda h(c)) <= alpha_P$ solves directly to $lambda >= log(1 / alpha_P) / h(c)$.
]

#subsec("B", "Termination and ANOVA Detection Power of Phase-1 Marking")

#theorem(title: "Finite O(log n) Termination and Exact ANOVA Power Bound of Chi-Square Marking", origin: "Original Proof")[
  Let attribute $k$ be partitioned into $beta = ceil(1 + log_2 n)$ Sturges bins.
  (i) *Termination & Complexity:* Phase-1 marking terminates in at most $beta - 1 = O(log n)$ iterations on coordinate $k$, executing in $O(n + log^2 n)$ time per dimension ($O(d (n + log^2 n))$ globally).
  (ii) *ANOVA Lower Bound for Planted Modes:* At any iteration with $m = abs(U_k)$ unmarked bins and total remaining count $N = sum_(b in U_k) c_(k, b)$, suppose $r in {1, dots, m - 1}$ elevated cluster bins have mean count $bar(c)_A = bar(c)_B + Delta$ ($Delta > 0$) above the $m - r$ background bins' mean $bar(c)_B$. Then the Pearson statistic satisfies the exact lower bound
  $ chi_k^2 >= ( r (m - r) / N ) Delta^2, $
  with equality iff bin counts are constant within both the elevated and background groups.
  (iii) *Null False-Discovery Control:* Under pure uniform noise across all $d$ coordinates, the probability of marking any bin on the initial sweep is at most $1 - (1 - alpha_"mark")^d <= d alpha_"mark"$ (asymptotically exact as $n -> infinity$).
]

#proof[
  (i) Constructing the $beta$-bin histogram takes $O(n)$ time. Each iteration evaluates $chi_k^2$ over $abs(U_k) <= beta$ bins in $O(beta)$ time and strictly removes one bin $b^*$ whenever $H_0$ is rejected, halting once $abs(U_k) <= 1$. Thus at most $beta - 1$ sweeps occur, totaling $O(beta^2) = O(log^2 n)$ arithmetic steps per attribute.
  (ii) Let $U_k = A union B$ with $abs(A) = r$, $abs(B) = m - r$, overall mean $bar(c) = N / m = (r / m) bar(c)_A + ((m - r) / m) bar(c)_B$. Then $bar(c)_A - bar(c) = ((m - r) / m) Delta$ and $bar(c)_B - bar(c) = -(r / m) Delta$. By the exact orthogonal ANOVA decomposition of squared deviations,
  $ sum_(b in U_k) (c_(k, b) - bar(c))^2 &= sum_(b in A) (c_(k, b) - bar(c)_A)^2 + sum_(b in B) (c_(k, b) - bar(c)_B)^2 \
  &quad + r (bar(c)_A - bar(c))^2 + (m - r) (bar(c)_B - bar(c))^2 \
  &>= [ r (m - r)^2 / m^2 + (m - r) r^2 / m^2 ] Delta^2 = (r (m - r) / m) Delta^2. $
  Dividing both sides by $bar(c) = N / m$ cancels $m$ and gives $chi_k^2 >= (r (m - r) / N) Delta^2$, with equality precisely when both within-group variance sums vanish.
  (iii) Under $H_0$, each coordinate's initial $chi_k^2$ converges in distribution to $chi^2_(beta - 1)$, rejecting with probability $alpha_"mark"$. Independence across $d$ noise coordinates (or Boole's union bound generally) bounds the family-wise initial false-marking probability by $1 - (1 - alpha_"mark")^d <= d alpha_"mark"$.
]

#subsec("C", "Failure of Anti-Monotonicity for Poisson Core Significance")

Unlike raw support in CLIQUE @agrawal1998clique, statistical significance is *not* downward-closed—exposing an intrinsic completeness gap in P3C's Apriori Phase-2 assembly.

#proposition(title: "Constructive Non-Monotonicity of Poisson Core Significance", origin: "Original Counterexample")[
  The Poisson significance predicate $p_"Pois" (O(R); E(R)) < alpha_P$ is neither monotone nor anti-monotone on the subspace lattice. Specifically, for every threshold $alpha_P in (10^(-20), 1/2)$, there exists a $3$-dimensional dataset ($n = 1000$) and dense 1D intervals $I_1, I_2, I_3$ such that *every* $2$-dimensional parent box $I_a times I_b$ ($1 <= a < b <= 3$) fails significance with $p_"Pois" > 1/2$, while the $3$-dimensional child box $R = I_1 times I_2 times I_3$ is overwhelmingly significant with $p_"Pois" < 10^(-20) << alpha_P$.
]

#proof[
  Partition $[0, 1]^3$ into $8$ octants defined by $I_k = [0, 1/2]$ and $I_k^c = (1/2, 1]$ for $k = 1, 2, 3$. Place exactly $250$ points in each of the $4$ even-parity octants $(I_1, I_2, I_3)$, $(I_1, I_2^c, I_3^c)$, $(I_1^c, I_2, I_3^c)$, and $(I_1^c, I_2^c, I_3)$, and $0$ points in the $4$ odd-parity octants ($n = 1000$; assume $I_1, I_2, I_3$ are marked via tight sub-peaks inside $[0, 1/2]$). Each 1D marginal has $abs(S(I_k)) = 250 + 250 = 500$ ($p_k = 1/2$). For any pair $a < b$, only one even-parity octant lies in $I_a times I_b$, giving observed count $O(I_a times I_b) = 250$ and expected independence count $E(I_a times I_b) = 500 dot.op 500 / 1000 = 250$. Hence $O(I_a times I_b) = E(I_a times I_b)$, yielding $p_"Pois"(250; 250) > 1/2 > alpha_P$, so *zero* 2D cores survive level $q = 2$. Yet for the 3D core $R = I_1 times I_2 times I_3$, $O(R) = 250$ while every parent split yields conditional expectation $lambda_R = O(I_a times I_b) dot.op (abs(S(I_c)) / n) = 250 dot.op (1/2) = 125$. Applying Theorem II.A(ii) with $lambda = 125$ and $c = 2$ gives
  $ p_"Pois"(250; 125) <= exp( -125 (2 log 2 - 1) ) < exp(-48.2) < 10^(-20) << alpha_P. $
  Thus level-wise Apriori pruning on Poisson significance misses $R$ completely.
]

#subsec("D", "EM Monotonicity and Coordinate-Scale Invariance")

#theorem(title: "Monotonic EM Likelihood Ascent and Diagonal-Affine Invariance", origin: "Classical — Cited")[
  (i) *Likelihood Monotonicity:* Phase-3 fuzzy EM refinement @dempster1977em @moise2006p3c satisfies $cal(L)(Theta^((t+1))) >= cal(L)(Theta^((t)))$ at every iteration with exact gap
  $ cal(L)(Theta^((t+1))) - cal(L)(Theta^((t))) >= sum_(i=1)^n "KL"( u_i^((t)) || u_i^((t+1)) ) >= 0, $
  converging to a stationary point of the mixture log-likelihood under bounded-variance regularization @wu1983convergence @mclachlan2000finite.
  (ii) *Coordinate-Affine Invariance:* All three phases of P3C (Sturges binning, Poisson core testing, and Mahalanobis EM refinement) are strictly invariant under arbitrary positive diagonal affine rescalings $bold(x)_i |-> bold(A) bold(x)_i + bold(t)$ ($bold(A) = "diag"(s_1, dots, s_d)$, $s_k > 0$) and coordinate permutations.
]

#proof[
  (i) Let $cal(L)(Theta) = sum_(i=1)^n log sum_(c=1)^K pi_c p(bold(x)_i | theta_c)$ and $u_(i c)^((t)) = "Pr"(z_i = c | bold(x)_i, Theta^((t)))$. Writing $log p(bold(x)_i | Theta) = log p(bold(x)_i, z_i | Theta) - log p(z_i | bold(x)_i, Theta)$ and taking expectation over $u_i^((t))$ yields the exact identity
  $ cal(L)(Theta) = cal(Q)(Theta | Theta^((t))) + H(u^((t))) + sum_(i=1)^n "KL"(u_i^((t)) || u_i (Theta)), $
  where $cal(Q)(Theta | Theta^((t))) = sum_(i, c) u_(i c)^((t)) log (pi_c p(bold(x)_i | theta_c))$ and $H(u^((t)))$ is constant in $Theta$. Because the M-step chooses $Theta^((t+1)) = arg max_Theta cal(Q)(Theta | Theta^((t)))$, we have $cal(Q)(Theta^((t+1)) | Theta^((t))) >= cal(Q)(Theta^((t)) | Theta^((t)))$; subtracting the identity at $Theta^((t))$ (where the KL term vanishes) yields $cal(L)(Theta^((t+1))) - cal(L)(Theta^((t))) >= sum_(i=1)^n "KL"(u_i^((t)) || u_i^((t+1))) >= 0$ @dempster1977em @wu1983convergence.
  (ii) Under $tilde(x)_(i k) = s_k x_(i k) + t_k$ ($s_k > 0$), the coordinate range $[a_k, b_k]$ maps to $[s_k a_k + t_k, s_k b_k + t_k]$, so every equal-width Sturges bin boundary maps affinely, leaving all bin counts $c_(k, b)$, marked intervals $cal(I)_k$, core supports $O(R)$, and Poisson $p$-values identical. Finally, Mahalanobis distance $(bold(x) - bold(mu))^(top) bold(Sigma)^(-1) (bold(x) - bold(mu))$ is invariant under all invertible affine transforms, preserving Phase-3 memberships and outlier flags identically.
]

#sec("III", "Complexity, Geometry and Operational Profile")

P3C achieves exceptional practical speed because Phase 1 scans each coordinate once in $O(n)$ time and tests at most $beta - 1 = O(log n)$ bin-removal steps ($O(d (n + log^2 n))$ total), completely avoiding $O(n^2)$ pairwise distance matrices. Phase 2 evaluates $abs(cal(K)_"cand")$ candidate cores in $O(abs(cal(K)_"cand") n)$ time using bitset or inverted-list intersections—small in typical tabular data due to aggressive Poisson pruning, though worst-case exponential $O(2^d n)$ if many high-dimensional cores nest. Phase 3 runs $T_"EM"$ iterations of diagonal-covariance EM in $O(T_"EM" n K d)$ time (or $O(T_"EM" n K d^2)$ for full covariance). Memory is $O(n d + n K + abs(cal(K)_"cand") d)$.

Geometrically, Phase 2 identifies axis-aligned hyper-rectangular cluster cores, while Phase 3 expands each core into a smooth elliptical Gaussian basin with explicit $chi^2$ Mahalanobis outlier boundaries.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(d (n + log^2 n) + abs(cal(K)_"cand") n + T_"EM" n K d)$],
    [Space Complexity], [$O(n d + n K + abs(cal(K)_"cand") d)$ (no $O(n^2)$ matrix)],
    [Statistical Model], [Sturges $chi^2$ 1D marking + Poisson core test + Gaussian EM],
    [Automatic $K$ & Outliers], [Determines $K = abs(cal(K))$ automatically; $chi^2$ Mahalanobis outlier filter],
    [Invariances], [Coordinate-wise positive affine & permutation invariant (Thm. II.D)],
    [Cluster Geometry], [Axis-aligned cores refined into convex elliptical clusters],
  ),
  caption: [Operational and geometric profile of P3C (Alg. 148).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *High-Throughput Tabular Fraud and Telemetry Screening:* Financial transaction logs and server telemetry feature large $n$, mixed attribute scales, unknown cluster counts $K$, and heavy outlier contamination. P3C's $O(n d)$ histogramming, automatic $K$ selection, and built-in Mahalanobis outlier filter fit these requirements natively @moise2006p3c.
2. *Automated Biomedical Subtyping Without Pre-Specifying $K$:* Clinical and biomarker cohorts frequently lack prior knowledge of how many patient subpopulations exist. P3C discovers statistically significant core populations via Poisson $p$-values while quantifying posterior membership uncertainty $u_(i c)$ through EM.
3. *Industrial Sensor Regime Discovery with Background Noise:* Manufacturing test benches log dozens of noisy continuous gauges where operational modes appear as elevated modes above near-uniform sensor ranges. Theorem II.B's ANOVA sensitivity isolates these active channels rapidly.
4. *Retail Customer Segmentation on Large CRM Tables:* When customer tables exceed tens of thousands of rows, $O(n^2)$ dissimilarity methods (COSA, FIRES) stall, whereas P3C scales linearly in $n$ across both Phase 1 and Phase 3 while delivering compact, interpretable cluster prototypes.
5. *Materials and Chemical Descriptor Screening:* High-throughput computational screening produces multi-property descriptor catalogs with mixed units; P3C's strict coordinate-affine invariance (Theorem II.D) eliminates manual feature-rescaling heuristics.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

P3C exhibits four distinct failure modes: (1) *Synergistic higher-order clusters* (Proposition II.C), where 3D+ interactions lack pairwise excess and fail Apriori Poisson pruning (later addressed by STATPC @moise2008statpc); (2) *Non-uniform unimodal noise*, where a heavy-tailed or Gaussian background attribute triggers false 1D Sturges marking because P3C's null hypothesis $H_0$ assumes uniformity; (3) *Non-convex manifold geometries*, where Phase-3 Gaussian EM over-splits curved spirals or moons; and (4) *Very small sample sizes*, where Sturges' rule $beta = ceil(1 + log_2 n)$ creates overly coarse histograms that blur adjacent modes.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Large Tabular / Telemetry], [High: $O(n d)$ histograms, sub-ms speed, outlier filter], [Watch non-uniform noise; STATPC @moise2008statpc],
    [Unknown-$K$ Convex Subspaces], [High: Poisson cores self-determine $K$ & subspaces], [Verify $H_0$ uniformity assumption],
    [Non-Convex Shapes], [Moderate: Phase-3 Gaussian EM favors ellipsoids], [Use FIRES @kriegel2005fires or COSA @friedman2004cosa],
    [XOR / High-Order Synergy], [Low: Apriori Poisson pruning misses zero-2D excess], [Use STATPC @moise2008statpc or INSCY @assent2008inscy],
    [Oblique Linear Subspaces], [Low: axis-aligned rectangular cluster cores], [Use ORCLUS @aggarwal2000orclus or SSC @elhamifar2013ssc],
  ),
  caption: [Domain suitability and failure-mode matrix for P3C.],
)

#sec("V", "Structural Relations and Empirical Behavior")

P3C bridges statistical grid-core discovery (CLIQUE @agrawal1998clique) and disjoint projected clustering (PROCLUS @aggarwal1999proclus), while directly inspiring STATPC @moise2008statpc.

#proposition(title: "Structural Reduction to CLIQUE Grid Cores and Elimination of PROCLUS Priors", origin: "Original Derivation")[
  (i) Replacing P3C's adaptive Sturges $chi^2$ marking and Poisson test $p_"Pois"(O(R); E(R)) < alpha_P$ in Phases 1–2 with a fixed uniform $xi$-grid and absolute support threshold $O(R) >= tau n$ reproduces the exact axis-aligned dense-unit discovery phase of CLIQUE @agrawal1998clique.
  (ii) Unlike PROCLUS @aggarwal1999proclus—which constrains output to exactly $K$ medoid clusters of mean dimensionality $l$—P3C generates an unrestricted $K = abs(cal(K))$ from maximal significant cores and assigns each cluster an individualized dimensionality $abs(S_c)$ while replacing hard medoids with fuzzy Gaussian EM posteriors.
]

#proof[
  (i) Under a fixed $xi$-bin partition of each coordinate, setting every bin with $c_(k, b) >= tau n$ as marked and replacing the Poisson tail test by $O(R) >= tau n$ turns Phase 2 into exact level-wise Apriori enumeration of all hyper-rectangular grid cells of support $>= tau n$, which is verbatim CLIQUE's dense-unit generation @agrawal1998clique.
  (ii) In PROCLUS, the $L_1$ medoid assignment minimizes within-cluster average dispersion subject to fixed $K$ and global cardinality budget $sum_(c=1)^K abs(S_c) = K l$. P3C removes both constraints: $K$ emerges directly as the cardinality of maximal statistically significant cores $cal(K)$, and each $S_c$ is tested independently per cluster against the background distribution, generalizing PROCLUS's hard projected partition to soft EM posteriors with explicit outlier class $O$.
]

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), P3C achieves *Iris* $"Sil" = 0.439$, $"ARI" = 0.587$; *Customers* $"Sil" = 0.347$; *Multi-Shape* $"Sil" = 0.557$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.302$ ($"ARI" = "N/A"$), with a remarkable mean runtime of *0.107 ms*. Every metric mirrors our formal analysis: first, P3C attains the *highest Silhouette on both Iris ($0.439$) and Customers ($0.347$)* among Chapters 145–148—and the second-highest on *High-D* ($0.302$, trailing only INSCY's $0.336$)—because Phase-3 Gaussian EM refinement produces maximally compact convex clusters around Poisson-validated cores. Second, its slightly lower *Iris* ARI ($0.587$ vs. FIRES $0.745$ and COSA $0.727$) reflects coarse Sturges binning at $n = 150$ ($beta = 9$), which can merge the overlapping *versicolor*–*virginica* modes into fewer cores. Third, on *Multi-Shape* ($"Sil" = 0.557$), elliptical EM refinement trails non-convex COSA and FIRES ($0.609$). Finally, its *0.107 ms* runtime is *116$times$ to 146$times$ faster* than COSA (12.455 ms), INSCY (14.553 ms), and FIRES (15.611 ms), vividly confirming the $O(d(n + log^2 n))$ histogram complexity advantage over $O(n^2)$ pairwise methods.
