#import "../template/ieee.typ": *

#alg_entry(
  id: 153,
  name: "Gustafson-Kessel (GK)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "1979",
  paradigm: "Adaptive Mahalanobis Fuzzy Clustering",
  complexity: $O((n c d^2 + c d^3) T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Fuzzy C-Means (Alg. 151 @bezdek1981fcm) enforces a single global inner-product matrix $A succ 0$ across all clusters, implicitly assuming that every cluster shares identical spherical or hyper-ellipsoidal geometry. As proven in Alg. 151, fixed-metric FCM fails whenever clusters exhibit distinct orientations, feature correlations, or elongated aspect ratios. Gustafson and Kessel @gustafson1979gk overcame this geometric rigidity by endowing each cluster $j in {1, dots, c}$ with its own adaptive SPD norm-inducing matrix $A_j in bb(S)_(++)^d$, allowing every cluster to deform independently to match local ellipsoidal covariance structure.

Given $X = {bold(x)_1, dots, bold(x)_n} subset RR^d$, clusters $c in {2, dots, n - 1}$, fuzzifier $m > 1$, and positive cluster volume constants $bold(rho) = (rho_1, dots, rho_c)^top in RR_(>0)^c$ (default $rho_j = 1$), GK minimizes
$
  J_m (U, V, cal(A)) = sum_(i=1)^n sum_(j=1)^c u_(i j)^m (bold(x)_i - bold(v)_j)^top A_j (bold(x)_i - bold(v)_j),
$
subject to $U in M_(f c)$, $V in RR^(d times c)$, and the per-cluster determinant constraints
$
  A_j in bb(S)_(++)^d, quad det(A_j) = rho_j, quad j = 1, dots, c.
$
Without the determinant constraint, the linear functional $A_j mapsto J_m$ would trivially collapse to $A_j -> 0$. Constraining $det(A_j) = rho_j$ fixes the hyper-volume of each unit Mahalanobis ellipsoid $\{bold(z) : bold(z)^top A_j bold(z) <= 1\}$ while permitting arbitrary principal-axis rotations and aspect-ratio stretching.

#definition(title: "Fuzzy Cluster Scatter and Covariance")[
  For $U in M_(f c)$ and $V in RR^(d times c)$, the _fuzzy scatter matrix_ $S_j$ and _fuzzy covariance matrix_ $F_j$ of cluster $j$ are
  $
    S_j &= sum_(i=1)^n u_(i j)^m (bold(x)_i - bold(v)_j) (bold(x)_i - bold(v)_j)^top, \
    F_j &= (S_j) / (sum_(i=1)^n u_(i j)^m).
  $
]

Because finite samples or collinear clusters can render $F_j$ singular or severely ill-conditioned, practical GK implementations apply shrinkage regularization $tilde(F)_j = (1 - gamma) F_j + gamma tau I_d$ (with $gamma in (0, 1)$, typically $tau = det(F_0)^(1 / d)$ from global covariance $F_0$) and optional eigenvalue-ratio clamping ($lambda_("max") / lambda_("min") <= beta$) prior to inversion.

#algbox(
  title: "Gustafson-Kessel Adaptive Fuzzy Clustering",
  inputs: [Data $X = {bold(x)_i}_(i=1)^n subset RR^d$, clusters $c$, fuzzifier $m > 1$, volumes $rho_j > 0$, shrinkage $gamma in [0, 1)$, tolerance $epsilon > 0$, max iterations $T$.],
  outputs: [Fuzzy partition $U^* in M_(f c)$, prototypes $V^*$, adaptive metrics $cal(A)^* = {A_j^*}_(j=1)^c$.],
)[
  1. *Initialize* $U^((0)) in M_(f c)$ (e.g., via FCM warm-start) and set $t <- 0$.
  2. *Repeat* until $norm(U^((t)) - U^((t-1)))_F < epsilon$ or $t = T$:
  3. #h(1.0em) *Prototype update:* For $j = 1, dots, c$,
     $ bold(v)_j^((t+1)) = (sum_(i=1)^n (u_(i j)^((t)))^m bold(x)_i) / (sum_(i=1)^n (u_(i j)^((t)))^m). $
  4. #h(1.0em) *Fuzzy covariance & metric update:* For $j = 1, dots, c$, compute $F_j^((t+1))$, apply shrinkage $tilde(F)_j = (1 - gamma) F_j^((t+1)) + gamma tau I_d$ (if $gamma > 0$), and form
     $ A_j^((t+1)) = (rho_j det(tilde(F)_j))^(1 / d) tilde(F)_j^(-1). $
  5. #h(1.0em) *Mahalanobis distance:* Compute $d_(i j)^2 = (bold(x)_i - bold(v)_j^((t+1)))^top A_j^((t+1)) (bold(x)_i - bold(v)_j^((t+1)))$.
  6. #h(1.0em) *Membership update:* Update $U^((t+1))$ via the FCM rule using $d_(i j)^2$ (handling $I_i = {j : d_(i j) = 0}$ identically).
  7. #h(1.0em) Set $t <- t + 1$.
  8. *Return* $(U^((t)), V^((t)), cal(A)^((t)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish global optimality of the GK metric update via an eigenvalue inequality, prove exact general-linear affine invariance, derive the reduced determinant-volume objective, and expose the algorithm's collinearity singularity and equal-volume blind spot.

#subsec("A", "Global Optimality of the Adaptive Metric Matrix")

#theorem(title: "Unique Global Minimizer on the Volume-Constrained SPD Cone", origin: "Original Proof")[
  Fix $U in M_(f c)$ and $V in RR^(d times c)$ such that the fuzzy scatter matrix $S_j$ is strictly positive definite ($S_j succ 0$). Then the constrained matrix minimization problem
  $ min_(A_j succ 0, det(A_j) = rho_j) sum_(i=1)^n u_(i j)^m (bold(x)_i - bold(v)_j)^top A_j (bold(x)_i - bold(v)_j) $
  admits the *unique global minimizer*
  $ A_j^* = (rho_j det(S_j))^(1 / d) S_j^(-1) = (rho_j det(F_j))^(1 / d) F_j^(-1), $
  attaining the exact minimum value $d dot rho_j^(1 / d) det(S_j)^(1 / d)$.
]

#proof[
  Using trace cyclicity, rewrite the cluster-$j$ objective as
  $ sum_(i=1)^n u_(i j)^m (bold(x)_i - bold(v)_j)^top A_j (bold(x)_i - bold(v)_j) = "Tr"(A_j S_j). $
  Because $S_j succ 0$, its unique SPD square root $S_j^(1 / 2)$ exists and is invertible. Define the congruence transform $B_j = S_j^(1 / 2) A_j S_j^(1 / 2) succ 0$, and let $mu_1, dots, mu_d > 0$ be the eigenvalues of $B_j$. By trace cyclicity and multiplicativity of the determinant,
  $
    "Tr"(A_j S_j) &= "Tr"(B_j) = sum_(r=1)^d mu_r, \
    product_(r=1)^d mu_r &= det(B_j) = det(A_j) det(S_j) = rho_j det(S_j).
  $
  Applying the arithmetic-geometric mean (AM-GM) inequality to $mu_1, dots, mu_d > 0$ immediately yields
  $
    "Tr"(A_j S_j) = sum_(r=1)^d mu_r >= d (product_(r=1)^d mu_r)^(1 / d) = d dot rho_j^(1 / d) det(S_j)^(1 / d).
  $
  Strictness of AM-GM guarantees that equality holds *if and only if* $mu_1 = dots = mu_d = c_j > 0$, i.e., $B_j = c_j I_d$, which is equivalent to $A_j = c_j S_j^(-1)$. Enforcing the constraint $det(A_j) = c_j^d det(S_j)^(-1) = rho_j$ uniquely determines $c_j = (rho_j det(S_j))^(1 / d)$. Finally, since $S_j = w_j F_j$ with $w_j = sum_(i=1)^n u_(i j)^m > 0$,
  $ det(S_j)^(1 / d) S_j^(-1) = (w_j^d det(F_j))^(1 / d) (w_j^(-1) F_j^(-1)) = det(F_j)^(1 / d) F_j^(-1), $
  proving invariance to scatter normalization and global optimality over the entire non-convex constraint manifold $det(A_j) = rho_j$.
]

Whereas the classical derivation @gustafson1979gk applies Lagrange multipliers to obtain a first-order stationary condition, the AM-GM argument above establishes *global* optimality directly without requiring second-order manifold Hessians.

#subsec("B", "Exact General Linear Affine Invariance")

#theorem(title: "Full GL(d, R) Affine Invariance of GK", origin: "Original Proof")[
  Let $T(bold(x)) = M bold(x) + bold(b)$ be any invertible affine transformation ($M in "GL"(d, RR)$, $bold(b) in RR^d$), and run unregularized GK ($gamma = 0$) on the transformed dataset $tilde(X) = {T(bold(x)_i)}_(i=1)^n$ with identical initial memberships $tilde(U)^((0)) = U^((0))$ and volume constants $rho_j$. Assuming non-singular scatter matrices at every step, for all $t >= 0$:
  1. $tilde(bold(v))_j^((t)) = M bold(v)_j^((t)) + bold(b)$ and $tilde(F)_j^((t)) = M F_j^((t)) M^top$;
  2. $tilde(A)_j^((t)) = abs(det(M))^(2 / d) M^(-top) A_j^((t)) M^(-1)$;
  3. Every squared distance scales uniformly by $abs(det(M))^(2 / d)$, leaving memberships identically invariant: $tilde(U)^((t)) = U^((t))$.
]

#proof[
  We proceed by induction on $t$. At $t = 0$, $tilde(U)^((0)) = U^((0))$ by hypothesis. Assuming $tilde(U)^((t)) = U^((t))$, linearity of weighted averages gives
  $ tilde(bold(v))_j = w_j^(-1) sum_(i=1)^n u_(i j)^m (M bold(x)_i + bold(b)) = M bold(v)_j + bold(b). $
  Consequently $tilde(bold(x))_i - tilde(bold(v))_j = M (bold(x)_i - bold(v)_j)$, so the transformed fuzzy covariance satisfies
  $ tilde(F)_j = w_j^(-1) sum_(i=1)^n u_(i j)^m M (bold(x)_i - bold(v)_j)(bold(x)_i - bold(v)_j)^top M^top = M F_j M^top. $
  Taking determinants yields $det(tilde(F)_j) = det(M)^2 det(F_j)$, hence $det(tilde(F)_j)^(1 / d) = abs(det(M))^(2 / d) det(F_j)^(1 / d)$. Substituting into Theorem 1 produces
  $
    tilde(A)_j &= (rho_j det(tilde(F)_j))^(1 / d) tilde(F)_j^(-1) \
               &= abs(det(M))^(2 / d) (rho_j det(F_j))^(1 / d) M^(-top) F_j^(-1) M^(-1) \
               &= abs(det(M))^(2 / d) M^(-top) A_j M^(-1).
  $
  The transformed Mahalanobis distance therefore becomes
  $
    tilde(d)_(i j)^2 &= (M (bold(x)_i - bold(v)_j))^top tilde(A)_j (M (bold(x)_i - bold(v)_j)) \
                     &= abs(det(M))^(2 / d) (bold(x)_i - bold(v)_j)^top A_j (bold(x)_i - bold(v)_j) \
                     &= abs(det(M))^(2 / d) d_(i j)^2.
  $
  Crucially, the positive scalar $abs(det(M))^(2 / d)$ is completely independent of cluster index $j$. It cancels identically in every ratio $tilde(d)_(i j) / tilde(d)_(i l) = d_(i j) / d_(i l)$ and preserves the zero set $I_i$, yielding $tilde(U)^((t+1)) = U^((t+1))$ and completing the induction.
]

#subsec("C", "Three-Block Descent and Determinant-Product Reformulation")

#proposition(title: "Three-Block Monotone Descent and Determinant Criterion", origin: "Classical — Cited")[
  Assuming $S_j^((t)) succ 0$ along the trajectory, GK three-block alternating optimization ($V$, $cal(A)$, $U$) monotonically decreases $J_m (U, V, cal(A))$ @gustafson1979gk @bezdek1981fcm. Eliminating $cal(A)$ reduces GK exactly to the fuzzy determinant-volume minimization problem
  $ min_(U in M_(f c), V in RR^(d times c)) d sum_(j=1)^c rho_j^(1 / d) det(S_j (U, V))^(1 / d). $
  When eigenvalues of $A_j$ are restricted to a compact box $[lambda_("min"), lambda_("max")] subset (0, infinity)$, every accumulation point of the iterate sequence is a KKT stationary point of the regularized objective.
]

#proof[
  Each iteration updates three coordinate blocks in sequence: $V^((t+1))$ uniquely minimizes $J_m (U^((t)), dot, cal(A)^((t)))$ because $A_j^((t)) succ 0$ cancels out of the first-order optimality condition $sum_i u_(i j)^m A_j (bold(x)_i - bold(v)_j) = bold(0)$; $cal(A)^((t+1))$ globally minimizes $J_m (U^((t)), V^((t+1)), dot)$ by Theorem 1; and $U^((t+1))$ uniquely minimizes $J_m (dot, V^((t+1)), cal(A)^((t+1)))$ by Lemma 1 of Alg. 151. Hence
  $
    J_m (U^((t+1)), V^((t+1)), cal(A)^((t+1))) &<= J_m (U^((t)), V^((t+1)), cal(A)^((t+1))) \
                                               &<= J_m (U^((t)), V^((t+1)), cal(A)^((t))) \
                                               &<= J_m (U^((t)), V^((t)), cal(A)^((t))).
  $
  Substituting Theorem 1's minimum value $d rho_j^(1 / d) det(S_j)^(1 / d)$ for each $j$ yields the exact reduced determinant criterion. Under compact spectral bounds $lambda_("min") I_d prec.eq A_j prec.eq lambda_("max") I_d$, the feasible set is compact and bounded away from the singularity $det(S_j) = 0$, so Zangwill's theorem guarantees subsequential convergence to stationary points.
]

#subsec("D", "Collinearity Singularity and Equal-Volume Blind Spot")

#theorem(title: "Rank-Deficiency Collapse and Isotropic Volume Cancellation", origin: "Original Counterexample")[
  Consider unregularized GK ($gamma = 0$) with equal volume parameters $rho_1 = dots = rho_c = 1$.
  1. _Subspace singularity_: If any $d$ observations lie in an affine hyperplane $H$ of dimension $d - 1$, then $inf J_m (U, V, cal(A)) = 0$ is approached degenerately by collapsing one cluster onto $H$ with $norm(A_j)_2 -> infinity$.
  2. _Equal-volume blind spot_: For any two isotropic clusters with covariances $F_1 = sigma_1^2 I_d$ and $F_2 = sigma_2^2 I_d$ ($sigma_2 / sigma_1 = kappa > 1$), GK produces *identical* metric matrices $A_1 = A_2 = I_d$, completely canceling the variance ratio $kappa^2$ and reverting to the exact Euclidean Voronoi bisector of FCM.
]

#proof[
  For (1), let $bold(n) in RR^d$ ($norm(bold(n))_2 = 1$) be normal to the affine hyperplane $H$ containing a prototype $bold(v)_1$ and points $bold(x)_i in H$. Consider the sequence $A_1 (delta) = delta^((1 - d) / d) bold(n) bold(n)^top + delta^(1 / d) (I_d - bold(n) bold(n)^top)$ for $delta > 0$. Its eigenvalues are $delta^((1 - d) / d)$ (multiplicity $1$) and $delta^(1 / d)$ (multiplicity $d - 1$), so $det(A_1 (delta)) = delta^((1 - d) / d) (delta^(1 / d))^(d - 1) = 1$ for all $delta > 0$. Because every point $bold(x)_i in H$ satisfies $chevron.l bold(x)_i - bold(v)_1, bold(n) chevron.r = 0$, its squared distance equals $delta^(1 / d) norm(bold(x)_i - bold(v)_1)_2^2 -> 0$ as $delta -> 0^+$, while $norm(A_1 (delta))_2 = delta^((1 - d) / d) -> infinity$. If a cluster attracts only points near $H$, $det(S_1) -> 0$ drives the reduced objective to zero degenerately.

  For (2), substitute $F_j = sigma_j^2 I_d$ into Theorem 1 with $rho_j = 1$:
  $
    A_j = det(sigma_j^2 I_d)^(1 / d) (sigma_j^2 I_d)^(-1) = (sigma_j^(2 d))^(1 / d) sigma_j^(-2) I_d = I_d,
  $
  independent of $sigma_j > 0$. Hence $d_(i 1)^2 = norm(bold(x)_i - bold(v)_1)_2^2$ and $d_(i 2)^2 = norm(bold(x)_i - bold(v)_2)_2^2$ identically, erasing all cluster-volume information.
]

#subsec("E", "Conditioning Guarantee Under Shrinkage Regularization")

#proposition(title: "Spectral Condition-Number Bound for Shrinkage GK", origin: "Original Derivation")[
  Let $tilde(F)_j = (1 - gamma) F_j + gamma tau I_d$ with $gamma in (0, 1)$ and $tau > 0$. Then $tilde(F)_j succ 0$ even when $"rank"(F_j) < d$, the regularized metric $tilde(A)_j = (rho_j det(tilde(F)_j))^(1 / d) tilde(F)_j^(-1)$ preserves orthogonal equivariance, and its spectral condition number satisfies
  $ kappa(tilde(A)_j) = (lambda_("max") (tilde(A)_j)) / (lambda_("min") (tilde(A)_j)) <= 1 + ((1 - gamma) / (gamma tau)) "Tr"(F_j). $
]

#proof[
  Since $F_j succ.eq 0$, every eigenvalue $lambda_r (F_j) >= 0$ satisfies $lambda_r (F_j) <= "Tr"(F_j)$. Therefore the eigenvalues of $tilde(F)_j$ lie in the strictly positive interval
  $ gamma tau <= lambda_r (tilde(F)_j) <= (1 - gamma) "Tr"(F_j) + gamma tau. $
  Because $tilde(A)_j$ is a positive scalar multiple of $tilde(F)_j^(-1)$, its condition number equals that of $tilde(F)_j$:
  $ kappa(tilde(A)_j) = (lambda_("max") (tilde(F)_j)) / (lambda_("min") (tilde(F)_j)) <= ((1 - gamma) "Tr"(F_j) + gamma tau) / (gamma tau) = 1 + ((1 - gamma) / (gamma tau)) "Tr"(F_j). $
  For any $Q in O(d)$, $Q tilde(F)_j Q^top = (1 - gamma) Q F_j Q^top + gamma tau I_d$, preserving orthogonal equivariance while eliminating rank-deficiency blow-up.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Relative to FCM's $O(n c d T)$ runtime, GK introduces substantial quadratic and cubic dimensional overhead. In each sweep, computing the $c$ fuzzy covariance matrices $F_j = w_j^(-1) sum_i u_(i j)^m (bold(x)_i - bold(v)_j)(bold(x)_i - bold(v)_j)^top$ takes $O(n c d^2)$ operations; inverting $F_j$ and computing $det(F_j)$ via Cholesky factorization $F_j = L_j L_j^top$ requires $O(c d^3)$ operations (with $log det(F_j) = 2 sum_(r=1)^d ln (L_j)_(r r)$ avoiding determinant overflow/underflow in large $d$); and evaluating the $n c$ Mahalanobis distances via triangular solves $L_j bold(z)_(i j) = bold(x)_i - bold(v)_j$ costs $O(n c d^2)$. Thus total time complexity is $O((n c d^2 + c d^3) T)$ and space complexity is $O(n c + n d + c d^2)$.

Geometrically, the decision boundary $d_(i j)^2 = d_(i l)^2$ between clusters $j$ and $l$ is a general quadric hypersurface $(bold(x) - bold(v)_j)^top A_j (bold(x) - bold(v)_j) = (bold(x) - bold(v)_l)^top A_l (bold(x) - bold(v)_l)$, allowing curved, obliquely tilted separating boundaries that match correlated ellipsoidal populations of comparable volume.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O((n c d^2 + c d^3) T)$ via Cholesky solves],
    [Space Complexity], [$O(n c + n d + c d^2)$ ($c$ SPD covariance matrices)],
    [Inductive Geometry], [Arbitrarily oriented equal-volume hyper-ellipsoids],
    [Hyperparameters], [$c$, fuzzifier $m > 1$, volumes $rho_j$, shrinkage $gamma$],
    [Invariance Group], [Full general linear affine group $"GL"(d, RR) times.l RR^d$ ($gamma = 0$)],
    [Breakdown Point], [$epsilon_n^* = 1 / n$ (quadratic scatter sensitive to outliers)],
    [Convergence Type], [Monotone 3-block descent (regularized KKT stationarity)],
  ),
  caption: [Operational and theoretical profile of Gustafson-Kessel (Alg. 153).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Takagi-Sugeno Nonlinear System Identification:* Input-output samples of piecewise-linear dynamical systems concentrate along thin, obliquely oriented hyperplanes in product state space. GK's adaptive covariance ellipsoids align their smallest eigenvectors orthogonal to each local hyperplane, recovering local linear models with high precision.
2. *Biological Morphometrics and Taxonomic Classification:* Anatomical measurements (such as sepal and petal dimensions in botany or cranial landmarks in zoology) exhibit strong positive allometric correlations. GK's affine invariance (Theorem 2) eliminates arbitrary unit-scaling effects and isolates tilted taxonomic ellipsoids without manual feature whitening.
3. *Industrial Process Monitoring and Fault Isolation:* Multivariate sensor arrays in chemical reactors and power plants produce tightly correlated variables whose fault states drift along oblique covariance directions rather than spherical halos; GK models each operating mode's correlation geometry directly.
4. *Computer Vision for Linear and Elliptical Structures:* Detecting elongated structures—road lanes, cracks, blood vessels, and elliptical cells—requires adaptive aspect ratios that spherical clustering fragments into multiple false sub-clusters.
5. *Quantitative Finance Regime Clustering:* Market volatility regimes differ dramatically in cross-asset correlation structure. GK groups joint return vectors by their local covariance geometry while assigning smooth regime-transition weights.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

GK faces four well-defined failure regimes: (i) _small-sample covariance singularity_ when effective cluster mass satisfies $n_("eff", j) <= d$ (Theorem 3(1)), necessitating shrinkage regularization (Proposition 2); (ii) _unequal cluster volumes_ with default $rho_j = 1$ (Theorem 3(2)), where dense cores steal boundary points from adjacent wide clusters; (iii) _extreme outlier vulnerability_, since quadratic outer products $(bold(x)_i - bold(v)_j)(bold(x)_i - bold(v)_j)^top$ inflate scatter determinants quadratically; and (iv) _computational scaling_ in high dimension $d >> 100$ due to $O(c d^3)$ inversions and $O(c d^2)$ covariance parameters.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Correlated Morphometrics], [Excellent: full $"GL"(d)$ invariance], [Apply shrinkage if $n_j$ small],
    [Piecewise-Linear Modeling], [Strong: aligns with hyperplanes], [Clamp eigen-ratio $lambda_("max") / lambda_("min")$],
    [Variable-Volume Clusters], [Poor under fixed $rho_j = 1$], [Use Gath-Geva (Alg. 154) or GMM],
    [High-$d$ Small-$n$ Data], [Rank-deficient $F_j$, $O(c d^3)$ cost], [Use PCA or diagonal covariances],
    [Heavy-Tailed Noise], [Scatter $F_j$ corrupted by outliers], [Combine with PCM or Student-$t$],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for GK.],
)

#sec("V", "Structural Relations and Empirical Behavior")

GK occupies the crucial structural bridge between spherical fuzzy clustering (FCM, Alg. 151) and full-covariance probabilistic mixture modeling (Gath-Geva, Alg. 154 @gath1989gg, and EM for Gaussian mixtures @dempster1977em @fraley2002mclust). Fixing $A_j = I_d$ recovers FCM; freeing cluster volumes and prior probabilities via Gaussian likelihood densities yields GG. Within the model-based clustering taxonomy @fraley2002mclust, GK's constant-volume free-shape free-orientation geometry corresponds directly to the equal-volume ellipsoidal Gaussian family (`EVV`), while its reduced objective $sum_j rho_j^(1 / d) det(S_j)^(1 / d)$ is the smooth fuzzy analogue of the classical Friedman-Rubin determinant criterion.

#proposition(title: "Volume-Invariant Equivalence with Whitened FCM", origin: "Original Proof")[
  If a dataset $X$ consists of clusters sharing a common covariance $Sigma succ 0$ ($F_1 = dots = F_c = Sigma$), then GK converges to the exact same partition as FCM applied to the globally whitened data $tilde(bold(x))_i = Sigma^(-1 / 2) bold(x)_i$.
]

#proof[
  When $F_j = Sigma$ for all $j$ and $rho_j = 1$, Theorem 1 gives $A_j = det(Sigma)^(1 / d) Sigma^(-1)$ identically across all clusters. Hence $d_(i j)^2 = det(Sigma)^(1 / d) norm(Sigma^(-1 / 2) (bold(x)_i - bold(v)_j))_2^2$. The shared constant $det(Sigma)^(1 / d)$ cancels in the membership ratios, matching Euclidean FCM on $tilde(bold(x))_i = Sigma^(-1 / 2) bold(x)_i$ identically.
]

Empirical results on the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing highlight GK's distinctive inductive bias vividly. On *Iris*, GK attains an extraordinary *ARI $= 0.859$* (Silhouette $= 0.393$), dramatically outperforming FCM ($0.594$), PCM ($0.510$), and GG ($0.559$). This leap occurs because Iris's three species form obliquely oriented, strongly correlated 4D ellipsoids of comparable volume—matching Theorem 2's affine-invariant equal-volume bias while avoiding GG's extra free-volume degrees of freedom on $n = 150$ points. Tellingly, GK's Euclidean Silhouette ($0.393$) is *lower* than FCM's ($0.447$) precisely because Euclidean Silhouette penalizes elliptical Mahalanobis boundaries that cut diagonally across Euclidean Voronoi midplanes, underscoring why Silhouette alone can undervalue anisotropic models. On *Customers* (Silhouette $= 0.318$) and *Multi-Shape* (Silhouette $= 0.709$, ARI $=$ N/A), GK matches or edges FCM, while on *High-D* (Silhouette $= 0.239$, ARI $=$ N/A) regularized covariance adaptation clearly improves over FCM ($0.188$) and avoids PCM's single-mode collapse. The *11.652 ms* mean runtime—roughly $14 times$ FCM's $0.813$ ms—directly reflects the $O((n c d^2 + c d^3) T)$ covariance assembly and inversion overhead.
