#import "../template/ieee.typ": *

#cmp_entry(
  roman: "III",
  title: "Intra-Category Comparative Analysis: Centroid, Medoid, and Partitioning (Algs. 1-30)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 1-30",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part I of this monograph develops thirty prototype-based partitioning algorithms (Algs. 1--30) spanning Euclidean vector quantization, adaptive randomized seeding, exact metric pruning, immediate-update combinatorial transfer, discrete medoid local search, impartial trimming, power-mean homotopy smoothing, automatic model selection, directional and kernelized geometries, feature-weighted subspaces, and instance-constrained clustering @jain2010data @xu2005survey. Despite their surface diversity, all thirty methods are instantiations of a single *unified prototype-partitioning variational template* minimizing an empirical distortion over an active subset $A subset.eq {1, dots, n}$, a $k$-partition $cal(C) = {C_1, dots, C_k}$ of $A$, prototype parameters $bold(theta) = (theta_1, dots, theta_k) in Theta^k$, and optional feature weights $W$:
$ cal(J)_Phi (A, cal(C), bold(theta), W) = sum_(j=1)^k sum_(i in C_j) Phi_j (x_i, theta_j; W) + Omega(k, W, cal(C)), $
where $Phi_j$ encodes the local dissimilarity (squared Euclidean, $ell_1$, $ell_p^p$, cosine, RKHS squared norm, Hamming, mixed Euclidean--Hamming, or arbitrary metric $d$), $Theta$ restricts prototypes to continuous means, coordinate medians, unit-sphere directions, discrete exemplars $X$, or RKHS spans, and $Omega$ incorporates BIC penalties, constraint feasibility, or simplex regularization.

This chapter conducts a rigorous head-to-head mathematical and algorithmic comparison of *all thirty* algorithms across five structural axes: (I) seeding dynamics, homotopy continuation, and approximation guarantees; (II) exact geometric acceleration and update-schedule dominance; (III) metric generalizations, discrete exemplars, and the exact finite-sample breakdown hierarchy; (IV) automatic model selection, divisive splitting, feature weighting, and constrained feasibility; and (V) a complete 30-algorithm comparative reference matrix and engineering decision architecture.

#sec("I", "Seeding Dynamics, Homotopy Smoothing, and Approximation Guarantees")

Minimizing Euclidean WCSS $J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2$ is NP-hard both for $k = 2$ in arbitrary dimension $d$ @aloise2009np and in the plane $d = 2$ for arbitrary $k$ @mahajan2012planar, while Lloyd's local refinement (Alg. 1) requires up to $2^(Omega(sqrt(n)))$ iterations in worst-case planar gadgets @vattani2011kmeans (albeit polynomial smoothed runtime @arthur2011smoothed). Consequently, global solution quality across Euclidean partitioning is governed primarily by *how* the initial prototype configuration or continuation landscape is constructed. We directly contrast seven foundational mechanisms: Lloyd (Alg. 1) @lloyd1982least, Forgy (Alg. 22) @forgy1965cluster, MacQueen (Alg. 21) @macqueen1967some, K-Means++ (Alg. 2) @arthur2007kmeanspp, K-Means|| (Alg. 23) @bahmani2012scalable, Global K-Means (Alg. 17) @likas2003global, and Power K-Means (Alg. 27) @xu2019power.

#subsec("A", "Initial Spread, Variance Gain, and One-Step Reduction")

A subtle yet consequential divide separates *random partition* seeding (assigning each point uniformly to ${1, dots, k}$ and taking initial means) from *Forgy data seeding* (Alg. 22, drawing $k$ distinct points uniformly without replacement), *$D^2$-weighted seeding* (Alg. 2), and *deterministic greedy insertion* (Alg. 17). Let $overline(x) = n^(-1) sum_(i=1)^n x_i$ and $sigma_cal(X)^2 = n^(-1) sum_(i=1)^n norm(x_i - overline(x))_2^2$ denote the empirical mean and total variance of $cal(X)$.

#theorem(title: "Exact Initial Centroid Spread and One-Step Gain Ordering", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ have $n >= 2 k$ and non-zero total variance $sigma_cal(X)^2 > 0$.
  1. *Initial Spread Collapse (Random Partition vs. Forgy):* If $bold(mu)^("RP")$ is formed by a balanced random $k$-partition (each cluster of size $m = n / k$) and $bold(mu)^("Forgy")$ by Forgy seeding (Alg. 22), then for every distinct pair $j != ell$,
     $ bb(E) [norm(mu_j^("RP") - mu_ell^("RP"))_2^2] &= (2 k) / (n - 1) sigma_cal(X)^2, \
       bb(E) [norm(mu_j^("Forgy") - mu_ell^("Forgy"))_2^2] &= (2 n) / (n - 1) sigma_cal(X)^2, $
     so Forgy seeds start strictly $n / k$ times farther apart in mean-squared separation, whereas random-partition centroids collapse onto $overline(x)$ as $n -> infinity$.
  2. *One-Step Seeding Gain Ordering ($k = 1 -> 2$):* Fix any first center $mu_1 in bb(R)^d$ with pointwise costs $w_i := norm(x_i - mu_1)_2^2$, average cost $overline(w) = n^(-1) sum_(i=1)^n w_i > 0$, and empirical variance $"Var"(w) = n^(-1) sum_(i=1)^n (w_i - overline(w))^2$. Let $Delta(y) := Phi({mu_1}) - Phi({mu_1, y})$ denote the exact distortion drop upon inserting candidate $y in cal(X)$, where $Phi(S) = sum_(i=1)^n min_(c in S) norm(x_i - c)_2^2$. Then the guaranteed self-elimination lower bounds $w(Y) <= Delta(Y)$ satisfy
     $ &bb(E)_("Forgy") [w(Y)] = overline(w) <= overline(w) + ("Var"(w)) / overline(w) \
       &quad = bb(E)_(D^2) [w(Y)] <= max_(y in cal(X)) w(y) <= max_(y in cal(X)) Delta(y), $
     with strict inequality between Forgy and $D^2$ (K-Means++, Alg. 2) whenever the $w_i$ are not all identical, and with the rightmost bound achieved deterministically by Global K-Means (Alg. 17) prior to Lloyd refinement.
]

#proof[
  (1) For Forgy seeding, $(mu_j^("Forgy"), mu_ell^("Forgy"))$ is an ordered pair of distinct points drawn uniformly without replacement from $cal(X)$. Hence
  $ bb(E) [norm(mu_j^("Forgy") - mu_ell^("Forgy"))_2^2] &= 1 / (n (n - 1)) sum_(i != r) norm(x_i - x_r)_2^2 \
    &= (2 n^2 sigma_cal(X)^2) / (n (n - 1)) = (2 n) / (n - 1) sigma_cal(X)^2, $
  using the Lagrange identity $sum_(i, r) norm(x_i - x_r)_2^2 = 2 n sum_i norm(x_i - overline(x))_2^2 = 2 n^2 sigma_cal(X)^2$. For a balanced random partition with $|C_j| = |C_ell| = m = n / k$, write $Z_i = x_i - overline(x)$. Sampling $m$ elements without replacement from a finite population with mean $bold(0)$ and covariance trace $sigma_cal(X)^2$ gives
  $ bb(E) [norm(mu_j^("RP") - overline(x))_2^2] = (n - m) / (m (n - 1)) sigma_cal(X)^2. $
  Moreover, for $i in C_j$ and $r in C_ell$ ($j != ell$), $bb(E)[chevron.l Z_i, Z_r chevron.r] = (n(n-1))^(-1) sum_(a != b) chevron.l Z_a, Z_b chevron.r = -(n - 1)^(-1) sigma_cal(X)^2$, since $sum_(a, b) chevron.l Z_a, Z_b chevron.r = norm(sum_a Z_a)_2^2 = 0$. Therefore,
  $ bb(E) [chevron.l mu_j^("RP") - overline(x), mu_ell^("RP") - overline(x) chevron.r] = - 1 / (n - 1) sigma_cal(X)^2. $
  Expanding $norm(mu_j^("RP") - mu_ell^("RP"))_2^2$ yields
  $ &bb(E) [norm(mu_j^("RP") - mu_ell^("RP"))_2^2] \
    &quad = 2 (n - m) / (m (n - 1)) sigma_cal(X)^2 + 2 / (n - 1) sigma_cal(X)^2 \
    &quad = (2 n) / (m (n - 1)) sigma_cal(X)^2 = (2 k) / (n - 1) sigma_cal(X)^2. $

  (2) Inserting $y = x_r in cal(X)$ reduces the cost of $x_r$ itself from $w_r$ to $0$ and never increases the cost of any other point, so $Delta(x_r) >= w_r$ for every $r$. Under Forgy, $bb(P)[Y = x_r] = 1 / n$, giving $bb(E)_("Forgy")[w(Y)] = overline(w)$. Under $D^2$-sampling, $bb(P)[Y = x_r] = w_r / (n overline(w))$, yielding
  $ bb(E)_(D^2) [w(Y)] &= (sum_(r=1)^n w_r^2) / (sum_(r=1)^n w_r) = (n^(-1) sum_(r=1)^n w_r^2) / overline(w) \
    &= (overline(w)^2 + "Var"(w)) / overline(w) = overline(w) + ("Var"(w)) / overline(w). $
  Since $"Var"(w) >= 0$ with equality iff all $w_r$ coincide, the second inequality is strict for every non-spherical shell around $mu_1$. Finally, every weighted average of ${w_r}$ is bounded above by $max_r w_r <= max_r Delta(x_r)$, and Global K-Means evaluates *every* candidate $x_r in cal(X)$ followed by monotone Lloyd refinement (Theorem 1 of Alg. 1), guaranteeing a post-insertion drop of at least $max_r Delta(x_r)$.
]

#subsec("B", "Strict Approximation-Ratio Hierarchy Across Seven Paradigms")

Theorem 1 quantifies local one-step progress, but global clustering accuracy hinges on worst-case approximation ratios $rho(cal(A)) := (bb(E)[J(cal(A))]) / J^*$ across multiscale geometries. The next unified theorem separates all seven seeding and smoothing paradigms sharply on explicit analytical instances.

#theorem(title: "Unified Seeding & Homotopy Approximation Hierarchy", origin: "Original Unified Proof")[
  Let $J^* := min_(cal(C), bold(mu)) J(cal(C), bold(mu))$ denote the global $k$-means optimum on $cal(X) subset.eq bb(R)^d$. Across Algs. 1, 2, 17, 21, 22, 23, and 27:
  1. *Exponential Miss Failure of Forgy (Alg. 22) and Unseeded Lloyd (Alg. 1):* For any $k >= 3$ and any target factor $M > 0$, there exists a balanced dataset of $k$ separated clusters ($n = k m$ with $m >= 4 k$) such that Forgy misses at least one cluster with probability $p_"miss" >= 1 - sqrt(2 pi k e) e^(-k) e^(1/(12 k))$, forcing $(bb(E)[J_"Forgy"]) / J^* > M$.
  2. *Adversarial Permutation Drift of MacQueen (Alg. 21):* Even when Forgy succeeds with positive probability, there exists an arrival ordering $pi$ of a 4-point 1D instance ($k = 3$) such that sequential MacQueen initialization freezes at a CVT with $(J_"MacQueen" (pi)) / J^* > M$.
  3. *Guaranteed $O(log k)$ Ceiling of K-Means++ (Alg. 2) vs. Fixed-Round Collapse of K-Means|| (Alg. 23):* K-Means++ unconditionally guarantees $bb(E)[J_(++)] <= 8 (ln k + 2) J^*$ @arthur2007kmeanspp across *every* dataset in $k$ sequential passes. K-Means|| with oversampling $ell = Theta(k)$ achieves $O(alpha)$-approximation (yielding $O(log k)$ when reclustered by K-Means++) in $r = O(log psi)$ parallel rounds (where $psi = Phi_1 / J^*$) @bahmani2012scalable; *however*, for any *fixed* round budget $r >= 1$ independent of the data aspect ratio, there exists a multiscale cascade on $bb(R)$ ($k >= r + 3$) on which K-Means|| with $r$ rounds incurs unbounded expected ratio $(bb(E)[J_(||, r)]) / J^* -> infinity$.
  4. *Strict Constant-Factor Suboptimality of Deterministic Global K-Means (Alg. 17):* Although Global K-Means deterministically solves the 4-point Forgy trap ${0, 1, L, L + D}$ ($k = 3$) to exact global optimality, on the symmetric 4-point instance $cal(X)_a = {-a, -1, 1, a}$ with $a in (3, 5)$ and $k = 3$, Global K-Means terminates deterministically at a suboptimal partition with exact ratio
     $ (J_"Global" (cal(X)_a)) / (J^* (cal(X)_a)) = ((a - 1)^2) / 4 -> 4^- quad (a -> 5^-), $
     whereas K-Means++ finds the exact global optimum with probability strictly bounded below by $1 / 2$.
  5. *Homotopy Trap Escape by Power K-Means (Alg. 27):* On the adversarial 4-point Lloyd/Forgy trap $cal(X) = {0, 1, L, L + D}$ ($L > 2 D + 2$) initialized directly *at* the trapping CVT $bold(mu)^((0)) = (0, 1, L + D / 2)$ where Lloyd and Forgy freeze with ratio $D^2$, Power K-Means at any finite negative exponent $s < 0$ breaks stationarity on iteration 1 with strictly positive escape gradient $nabla_(mu_2) F_(s, epsilon) (bold(mu)^((0))) < 0$, while satisfying uniform objective approximation $F_(-infinity) (bold(mu)) <= F_(s, 0) (bold(mu)) <= k^(-1 / s) F_(-infinity) (bold(mu))$ @xu2019power.
]

#proof[
  (1) Place $k$ clusters of equal size $m >= 4 k$ inside disjoint balls of unit diameter centered at pairwise distances $L > 2 sqrt(M k m)$. Any $k$-sample drawn without replacement covers all $k$ clusters with probability
  $ p_"hit" = m^k / binom(k m, k) = (k!) / (k^k) product_(t=0)^(k-1) (1 - t / (k m))^(-1). $
  Because $t <= k - 1 < k$ implies $(1 - t / (k m))^(-1) <= (1 - 1 / m)^(-1) <= e^(2 / m)$ for $m >= 2$, choosing $m >= 4 k$ gives $p_"hit" <= (k! / k^k) e^(2 k / m) <= sqrt(2 pi k e) e^(-k) e^(1 / (12 k))$ by Stirling's inequality. Thus $p_"miss" = 1 - p_"hit"$ approaches $1$ exponentially fast in $k$ (already $p_"miss" > 0.63$ for $k >= 3$). Whenever a cluster is missed, two distant clusters share a center, paying at least $m (L / 2)^2$ while $J^* <= k m$, so $(bb(E)[J_"Forgy"]) / J^* >= p_"miss" L^2 / (4 k) > M$ for $L > 2 sqrt(M k / p_"miss")$.

  (2) On $cal(X) = {0, 1, L, L + D}$ ($k = 3$, $L > 2 D + 2$, $D = sqrt(M) + 1$), present the stream in order $(x_1, x_2, x_3, x_4) = (0, 1, L, L + D)$. MacQueen sets initial centers $(0, 1, L)$ from the first $k = 3$ observations, assigns $x_4 = L + D$ to cluster $3$, updates $mu_3 <- L + D / 2$, and halts at the strict CVT ${{0}, {1}, {L, L + D}}$ with $(J_"MacQueen") / J^* = D^2 > M$.

  (3) Arthur and Vassilvitskii @arthur2007kmeanspp prove $bb(E)[J_(++)] <= 8 (ln k + 2) J^*$ by induction on uncovered optimal clusters using the factor-2 uniform first-center identity and factor-8 $D^2$-hit lemma. For K-Means|| (Alg. 23), Bahmani et al. @bahmani2012scalable prove that each oversampling round contracts expected cost by a constant factor $beta < 1$ plus $O(J^*)$, requiring $r = O(log(Phi_1 / J^*))$ rounds. Now fix any constant $r >= 1$, set $k = r + 3$, and construct a 1D geometric cascade of $r + 3$ clusters $C_0, C_1, dots, C_(r+2)$ at scales $R_0 = 0 < R_1 << R_2 << dots << R_(r+2)$ with $R_(t+1) / R_t = Lambda -> infinity$: cluster $C_0$ holds $n_0$ points at $0$, while each cluster $C_t$ ($t >= 1$) holds two points at $R_t$ and $R_t + 1$ (so $J^* = (r + 2) / 2$). Starting from a uniform seed in $C_0$ (which occurs with probability $1 - O(1 / n_0)$), the initial potential is dominated by the *outermost* scale $Phi_0 = Theta(R_(r+2)^2)$. In round $1$, every point in $C_t$ ($t <= r + 1$) has sampling probability $p_i <= ell R_t^2 / R_(r+2)^2 = O(ell Lambda^(-2)) -> 0$ as $Lambda -> infinity$. Thus round 1 covers $C_(r+2)$ with high probability, reducing the residual potential to $Theta(R_(r+1)^2)$, while missing all inner clusters $C_1, dots, C_(r+1)$ with probability $1 - O(ell Lambda^(-2))$. Inductively, each round peels off *at most one* scale level from the outside in. After $r$ rounds, clusters $C_1$ and $C_2$ remain completely unseeded with probability $1 - O(r ell Lambda^(-2))$, forcing any subset of $k$ centers chosen from the coreset to leave at least one scale $R_1$ uncovered, incurring cost $Omega(R_1^2) / J^* -> infinity$.

  (4) Consider $cal(X)_a = {-a, -1, 1, a} subset.eq bb(R)$ with $a in (3, 5)$ and $k = 3$.
  - *Step $m = 1$:* The global mean is $0$.
  - *Step $m = 2$:* Evaluate all possible 2-partitions. The symmetric split ${{-a, -1}, {1, a}}$ has centroids $plus.minus (a + 1) / 2$ and cost $2 [ ((a - 1) / 2)^2 + ((a - 1) / 2)^2 ] = (a - 1)^2$. Any 1-vs.-3 split such as ${{-a}, {-1, 1, a}}$ has centroids $-a$ and $a / 3$ and cost $(-1 - a/3)^2 + (1 - a/3)^2 + (2 a / 3)^2 = 2 a^2 / 3 + 2$. Subtracting yields
    $ (2 / 3 a^2 + 2) - (a - 1)^2 = - 1 / 3 (a^2 - 6 a - 3) > 0 $
    for all $a in (3, 3 + 2 sqrt(3)) approx (3, 6.464)$, which strictly contains $(3, 5)$. Hence Global K-Means uniquely selects the 2-means solution $bold(mu)^((2)) = (-(a+1)/2, (a+1)/2)$.
  - *Step $m = 3$:* Global K-Means initializes Lloyd from $(-(a+1)/2, (a+1)/2, x_r)$ for each $x_r in cal(X)_a$.
    - Trying $x_r = -a$: initial Voronoi assignments map $-a -> -a$, $-1 -> -(a+1)/2$ (since $(a-1)/2 < a - 1$), and $1, a -> (a+1)/2$. Updating means produces $(-a, -1, (a+1)/2)$. Check Voronoi stability at $x = 1$: its distance to $-1$ is $2$, whereas its distance to $(a+1)/2$ is $(a - 1) / 2 < 2$ because $a < 5$! Thus $1$ stays assigned to $(a+1)/2$, so ${{-a}, {-1}, {1, a}}$ is an immediate *strict* Lloyd fixed point of cost $((a-1)/2)^2 dot 2 = (a - 1)^2 / 2$.
    - Trying $x_r = -1$: initial assignments map $-a -> -(a+1)/2$ ($(a-1)/2 < a-1$), $-1 -> -1$, and $1 -> (a+1)/2$ ($(a-1)/2 < 2$ since $a < 5$), $a -> (a+1)/2$. One mean update yields the exact same fixed point $(-a, -1, (a+1)/2)$ with cost $(a - 1)^2 / 2$.
    - By mirror symmetry $x |-> -x$, trying $x_r = 1$ or $x_r = a$ yields ${{-a, -1}, {1}, {a}}$ with identical cost $(a - 1)^2 / 2$.
  - Meanwhile, the optimal 3-partition ${{-a}, {-1, 1}, {a}}$ has centroids $(-a, 0, a)$ and cost $J^* = 1^2 + 1^2 = 2$. Hence $J_"Global" / J^* = (a - 1)^2 / 4 -> 4^-$ as $a -> 5^-$. Under K-Means++ on $cal(X)_a$, if the first seed is $c_1 in {-a, a}$ (prob $1 / 2$), then $c_2 = -c_1$ has squared distance $4 a^2$ vs. $(a plus.minus 1)^2$, so $c_2 = -c_1$ is picked with prob $(4 a^2) / (4 a^2 + (a+1)^2 + (a-1)^2) = (2 a^2) / (3 a^2 + 1) > 9 / 14$; then the third seed is either $-1$ or $1$, whose Voronoi partition is ${{-a}, {-1, 1}, {a}}$ immediately! If instead $c_1 in {-1, 1}$ (prob $1 / 2$), $c_2 in {-a, a}$ is selected with prob $(a^2 + 1) / (a^2 + 3) > 5 / 6$ (since $(a-1)^2 + (a+1)^2 = 2 a^2 + 2$ vs. $4$ for $-c_1$), and $c_3 = -c_2$ with prob $> (a-1)^2 / ((a-1)^2 + 4) > 1/2$ (for $a > 3$), also yielding the exact global optimum; summing both branches exceeds $1 / 2$ uniformly on $a in (3, 5)$.

  (5) As established via the power-mean gradient expansion in Alg. 27, at $bold(mu)^((0)) = (0, 1, L + D / 2)$ with $s < 0$ and $epsilon -> 0^+$, the non-nearest responsibility weights $w_(3 2), w_(4 2) = Theta((D / (2 L))^(2 (1 - s))) > 0$ remain strictly positive while $w_(1 2) = O(epsilon^(1 - s)) -> 0$, yielding $nabla_(mu_2) F_(s, epsilon) (bold(mu)^((0))) = -2 [w_(3 2)(L - 1) + w_(4 2)(L + D - 1) - w_(1 2)] < 0$. Hence the first MM iterate pulls $mu_2$ strictly rightward toward ${x_3, x_4}$, breaking the trapping barrier while Lemma 1(3) of Alg. 27 guarantees $F_(-infinity) <= F_(s, 0) <= k^(-1/s) F_(-infinity)$ uniformly.
]

#remark(title: "Outlier Center Starvation Trade-Off")[
  The very property that gives K-Means++ (Alg. 2) and Global K-Means (Alg. 17) their superior cluster-discovery power — prioritizing points with maximal squared distance to existing seeds — makes both *more vulnerable than Forgy* to isolated gross outliers! Specifically, a single outlier $x_"out"$ at distance $R$ from a clean dataset of WCSS $J_"clean"$ is chosen on the second step of K-Means++ with probability $R^2 / (R^2 + J_"clean") -> 1$ and deterministically by Global K-Means, starving clean clusters of prototypes, whereas Forgy selects $x_"out"$ with probability only $k / n$.
]

#subsec("C", "Parametric Unification of Seeding and Homotopy Families")

Rather than isolated heuristics, the seeding and smoothing algorithms of Category I form three continuous parametric bridges:

#proposition(title: "Three Continuous Unification Bridges in Category I", origin: "Original Unified Proof")[
  1. *Inverse-Temperature Seeding Family ($D^(2 gamma)$):* Sampling each new seed with probability $p_i (gamma) prop d(x_i, S)^(2 gamma)$ ($gamma >= 0$) recovers Forgy (Alg. 22) at $gamma = 0$, K-Means++ (Alg. 2) at $gamma = 1$, and deterministic farthest-first traversal (the greedy upper-bound surrogate of Global K-Means, Alg. 17) as $gamma -> +infinity$.
  2. *Power-Mean Homotopy Family ($M_s$):* Power K-Means (Alg. 27) with exponent $s < 0$ satisfies $lim_(s -> -infinity) F_(s, 0) (bold(mu)) = F_(-infinity) (bold(mu))$ (Lloyd, Alg. 1) at uniform multiplicative rate $k^(-1/s) - 1 = O(|s|^(-1) log k)$, while at $s = -1$ it minimizes the harmonic $k$-means potential whose smooth basins dissolve Voronoi walls.
  3. *Mini-Batch Interpolation ($b in {1, dots, n}$):* Mini-Batch K-Means (Alg. 3) with batch size $b$ and per-center step size $eta_(j, t) = 1 / v_j$ recovers sequential MacQueen (Alg. 21) at $b = 1$ and exact batch Lloyd (Alg. 1) at $b = n$ when cumulative counts $v_j$ are reset each epoch @sculley2010webscale.
]

#proof[
  (1) At $gamma = 0$, $d(x_i, S)^0 = 1$ (on $cal(X) backslash S$), yielding uniform without-replacement sampling (Forgy). At $gamma = 1$, $p_i (1) = (d(x_i, S)^2) / (sum_r d(x_r, S)^2)$ is exact $D^2$-sampling (K-Means++). As $gamma -> infinity$, let $M = max_r d(x_r, S)$ and $A = {r : d(x_r, S) = M}$; dividing numerator and denominator by $M^(2 gamma)$ gives $p_i (gamma) -> |A|^(-1) bb(1)_(i in A)$, concentrating identically on farthest-point selection.
  (2) By Lemma 1(3) of Alg. 27, $min_j y_(i j) <= M_s (y_i) <= k^(-1/s) min_j y_(i j)$ for every $i$ and $s < 0$; summing over $i = 1, dots, n$ yields $0 <= F_(s, 0) - F_(-infinity) <= (k^(-1/s) - 1) F_(-infinity)$, and $k^(-1/s) - 1 = exp((ln k) / (|s|)) - 1 <= (2 ln k) / (|s|)$ for $|s| >= ln k$.
  (3) For $b = 1$, a single point $x_t$ increments $v_(c_t) <- v_(c_t) + 1$ and sets $mu_(c_t) <- (1 - v_(c_t)^(-1)) mu_(c_t) + v_(c_t)^(-1) x_t$, which is algebraically identical to MacQueen's running arithmetic mean update. For $b = n$ with reset $v_j <- 0$ at each epoch start, aggregating all points $i in C_j$ sequentially gives $(|C_j|)^(-1) sum_(i in C_j) x_i$, recovering Lloyd's exact batch centroid update.
]

#wide_figure(
  caption: [*Table I:* Head-to-Head Comparison of Seeding & Homotopy Smoothing Algorithms (Algs. 1, 2, 17, 21, 22, 23, 27).],
  table(
    columns: (1.1fr, 1.25fr, 0.7fr, 1.35fr, 1.6fr),
    align: (left, left, center, center, left),
    [Algorithm], [Seeding / Landscape], [Passes], [Worst-Case Ratio], [Primary Failure Mode],
    [Alg. 1 Lloyd], [External / Rand. Part.], [$T$], [Unbounded ($infinity$)], [Centroid collapse & local CVTs],
    [Alg. 22 Forgy], [Uniform $k$-subset], [$1 + T$], [Unbounded ($infinity$)], [Exp. cluster-miss prob. in $k$],
    [Alg. 21 MacQueen], [First $k$ online pts.], [$1$--$T$], [Unbounded ($infinity$)], [Input-order & primacy drift],
    [Alg. 2 K-Means++], [Sequential $D^2$], [$k + T$], [$8(ln k + 2)$], [Outlier center starvation],
    [Alg. 23 K-Means||], [$r$-round $ell D^2$], [$r + T$], [$O(log k)$ if $r = Omega(log psi)$], [Unbounded if $r = O(1)$ fixed],
    [Alg. 17 Global], [Exact $m -> m + 1$], [$O(n k^2 T)$], [Deterministic; $-> 4^-$ on $n=4$], [Quadratic $O(n^2 k^2)$ scaling],
    [Alg. 27 Power], [$M_s$ MM homotopy], [$T_"hom"$], [Escapes Lloyd traps], [Slow annealing cost; $1/n$ brk.],
  ),
)

#sec("II", "Exact Geometric Acceleration and Update-Schedule Dominance")

Once prototypes are initialized, Euclidean $k$-means variants diverge along two orthogonal axes: *exact geometric pruning* of Lloyd's batch iteration (Elkan Alg. 13 @elkan2003using and Hamerly Alg. 14 @hamerly2010making) versus *changing the optimization trajectory itself* via immediate single-point transfers (Hartigan--Wong Alg. 19 @hartigan1979algorithm), sequential online updates (MacQueen Alg. 21), or stochastic mini-batches (Mini-Batch K-Means Alg. 3 @sculley2010webscale).

#subsec("A", "Exact Trajectory Equivalence & Pruning Dominance: Elkan vs. Hamerly vs. Lloyd")

Let $delta_j^((t)) := norm(mu_j^((t+1)) - mu_j^((t)))_2$ denote the drift of centroid $j$ after iteration $t$, and let $s^((t)) (j) := (1 / 2) min_(ell != j) norm(mu_j^((t)) - mu_ell^((t)))_2$ be the half-margin to the nearest rival centroid. Elkan (Alg. 13) maintains one upper bound $u_E (i) >= norm(x_i - mu_(c(i)))_2$ and $k$ lower bounds $l_E (i, j) <= norm(x_i - mu_j)_2$ per point, whereas Hamerly (Alg. 14) compresses the entire lower-bound row into a single scalar $l_H (i) <= min_(j != c(i)) norm(x_i - mu_j)_2$.

#theorem(title: "Exact Lloyd Equivalence and Strict Per-Point Pruning Dominance", origin: "Original Unified Proof")[
  Fix any dataset $cal(X) subset.eq bb(R)^d$, initial centers $bold(mu)^((0))$, and deterministic smallest-index tie-breaking. Let $cal(N)_E^((t)) (i)$ and $cal(N)_H^((t)) (i)$ denote the number of point-to-centroid Euclidean distance evaluations performed for point $x_i$ at iteration $t$ by Elkan (Alg. 13) and Hamerly (Alg. 14), respectively.
  1. *Bit-Exact Equivalence:* At every iteration $t >= 0$, Lloyd (Alg. 1), Elkan (Alg. 13), and Hamerly (Alg. 14) produce *identical* assignments $c^((t))$ and centroids $bold(mu)^((t))$.
  2. *Strict Pruning Dominance of Elkan over Hamerly:* Under synchronized drift updates $u_E (i) <- u_E (i) + delta_(c(i))$, $l_E (i, j) <- max(0, l_E (i, j) - delta_j)$, $u_H (i) <- u_H (i) + delta_(c(i))$, and $l_H (i) <- max(0, l_H (i) - delta_"max" (c(i)))$ (where $delta_"max" (a) := max_(ell != a) delta_ell$), the maintained bounds obey $l_H (i) <= min_(j != c(i)) l_E (i, j)$ and $u_E (i) <= u_H (i)$ whenever $c(i)$ is preserved across steps. Consequently, whenever Hamerly's outer gate $u_H (i) <= max(l_H (i), s(c(i)))$ prunes point $x_i$ with zero distance evaluations, Elkan *also* incurs zero distance evaluations ($cal(N)_E^((t)) (i) = 0$); moreover, for every $t$ where both enter inner loops, $cal(N)_E^((t)) (i) <= cal(N)_H^((t)) (i)$ up to tight-bound refreshes.
  3. *Drift-Spillover Separation ($Theta(k)$ Gap):* For every $k >= 3$, there exists a dataset and iteration $t$ in which moving *only one* distant centroid $mu_k$ by $delta_k$ forces Hamerly to execute $k$ distance evaluations for every point in stationary clusters $C_1, dots, C_(k-2)$ ($cal(N)_H^((t)) (i) = k$), whereas Elkan executes at most $1$ distance evaluation per point (testing only $j = k$, or $0$ when $l_E (i, k) - delta_k >= u_E (i)$).
]

#proof[
  (1) By induction on $t$. At $t = 0$, all three algorithms compute exact distances to $bold(mu)^((0))$ and break ties by smallest index, giving identical $c^((0))$, $bold(mu)^((1))$, and exact initial bounds $u_E (i) = u_H (i) = norm(x_i - mu_(c(i))^((0)))_2$, $l_E (i, j) = norm(x_i - mu_j^((0)))_2$, and $l_H (i) = min_(j != c(i)) norm(x_i - mu_j^((0)))_2$. Assuming admissibility at step $t$, the triangle inequality guarantees after moving $mu_j^((t)) -> mu_j^((t+1))$:
  $ norm(x_i - mu_(c(i))^((t+1)))_2 &<= norm(x_i - mu_(c(i))^((t)))_2 + delta_(c(i))^((t)) <= u(i) + delta_(c(i))^((t)), \
    norm(x_i - mu_j^((t+1)))_2 &>= norm(x_i - mu_j^((t)))_2 - delta_j^((t)) >= l_E (i, j) - delta_j^((t)), $
  and for every $j != c(i)$, $norm(x_i - mu_j^((t+1)))_2 >= min_(ell != c(i)) norm(x_i - mu_ell^((t)))_2 - max_(ell != c(i)) delta_ell^((t)) >= l_H (i) - delta_"max" (c(i))$. Furthermore, if $u(i) <= s(c(i))$, then for all $j != c(i)$,
  $ norm(x_i - mu_j)_2 >= norm(mu_(c(i)) - mu_j)_2 - norm(x_i - mu_(c(i)))_2 >= 2 s(c(i)) - u(i) >= u(i). $
  Strict pruning gates ($u < l$ or $u < s$, using $<$ on tie indices $j < c(i)$ and $<=$ for $j > c(i)$) therefore never skip a strictly closer or lower-indexed tied centroid. Hence every pruned point retains the exact Lloyd assignment, preserving $c^((t+1))$ and $bold(mu)^((t+2))$ identically across all three algorithms.

  (2) Assume $u_E (i) <= u_H (i)$ and $l_H (i) <= min_(j != c(i)) l_E (i, j)$ prior to drift propagation. Since $delta_j <= delta_"max" (c(i))$ for every $j != c(i)$, subtracting drifts preserves
  $ l_H' (i) &= max(0, l_H (i) - delta_"max" (c(i))) \
    &<= min_(j != c(i)) max(0, l_E (i, j) - delta_j) = min_(j != c(i)) l_E' (i, j), $
  while $u_E' (i) = u_E (i) + delta_(c(i)) <= u_H (i) + delta_(c(i)) = u_H' (i)$. If Hamerly's outer gate passes ($u_H' (i) <= max(l_H' (i), s(c(i)))$), then either $u_E' (i) <= u_H' (i) <= s(c(i))$ (Elkan's Lemma 1 gate triggers immediately, zero distances) or $u_E' (i) <= l_H' (i) <= l_E' (i, j)$ for *every* $j != c(i)$ (Elkan skips every candidate $j$ without even tightening $u_E$, zero distances). Whenever Hamerly fails both gates, it scans *all* $k$ centroids ($k$ distance evaluations), whereas Elkan evaluates at most $1 + |{j != c(i) : u_E (i) > max(l_E (i, j), (1/2) norm(mu_(c(i)) - mu_j)_2)}| <= k$ distances.

  (3) Place $k - 1$ stationary clusters $C_1, dots, C_(k-1)$ in $bb(R)^2$ with centroids $mu_1 = (0, 0)$, $mu_2 = (4, 0)$, and $mu_3, dots, mu_(k-1)$ at mutual pairwise distance $>= 4$ (so $s(1) = 2.0$). Consider a point $x_i in C_1$ placed at $(1.95, 0.64)$ so that $u_E (i) = u_H (i) = norm(x_i - mu_1)_2 = sqrt(1.95^2 + 0.64^2) approx 2.052 > s(1) = 2.0$ while $l_H (i) = l_E (i, 2) = norm(x_i - mu_2)_2 = sqrt(2.05^2 + 0.64^2) approx 2.148$ and $l_E (i, ell) >= 2.148$ for all $ell in {2, dots, k - 1}$. Place a single distant cluster $C_k$ centered at $mu_k = (100, 0)$ whose centroid shifts by $delta_k = 5$ on iteration $t$ while $mu_1, dots, mu_(k-1)$ remain completely stationary ($delta_1 = dots = delta_(k-1) = 0$). Then:
  - In Hamerly, $delta_"max" (1) = delta_k = 5$, so $l_H' (i) = max(0, 2.148 - 5) = 0$. Since $u_H' (i) = 2.052 > max(0, 2.0)$, Hamerly's outer gate fails, tightening $u_H$ keeps $2.052 > 2.0$, and Hamerly computes *all $k$ distances* for $x_i$!
  - In Elkan, $l_E' (i, j) = l_E (i, j) >= 2.148 > 2.052 = u_E' (i)$ for all $j in {2, dots, k - 1}$, while $l_E' (i, k) = norm(x_i - mu_k)_2 - 5 >= 93 > 2.052$. Every single candidate $j != 1$ is pruned immediately, so Elkan executes *$0$ distance evaluations* (and even if $delta_k = 98$, at most $1$ evaluation for $j = k$), yielding an exact $Theta(k)$ gap!
]

Despite Theorem 3(2)--(3), Elkan does *not* unconditionally dominate Hamerly in wall-clock runtime. Elkan stores $n k$ lower bounds ($8 n k$ bytes in double precision), whereas Hamerly stores only $2 n$ scalars ($16 n$ bytes). Every Elkan pass touches $O(n k)$ memory entries for drift updates even when zero distances are computed; whenever memory-bandwidth latency per $k$-entry bound scan exceeds $O(k d)$ cache-resident fused distance arithmetic (typically $d <= 16$--$32$ with $k >= 50$, or $n k > 10^8$ overflowing L3 cache/RAM), Hamerly runs strictly faster @hamerly2010making! Furthermore, in ultra-high dimension ($d -> infinity$) with isotropic noise, distance concentration @beyer1999nearest shrinks $s(c(i)) / u(i)$ toward $0$, causing *both* Elkan and Hamerly to degenerate toward Lloyd's $O(n k d)$ inner work.

#subsec("B", "Equilibrium Hierarchy & Update Schedules: Hartigan–Wong vs. Lloyd vs. MacQueen vs. Mini-Batch")

Unlike Elkan and Hamerly, Hartigan--Wong (Alg. 19) @hartigan1979algorithm alters the optimization trajectory itself by testing whether transferring a single point $x_i$ from $C_a$ ($n_a >= 2$) to $C_b$ strictly decreases exact combinatorial WCSS $J^* (cal(C))$ *after* moving both centroids simultaneously. By Lemma 1 of Alg. 19, the exact net energy reduction is
$ Delta J^* (i : a -> b) = n_a / (n_a - 1) norm(x_i - mu_a)_2^2 - n_b / (n_b + 1) norm(x_i - mu_b)_2^2. $
We establish the exact mathematical relationship among all four update trajectories — Hartigan--Wong, Lloyd, MacQueen, and Mini-Batch K-Means — including the asymptotic rate at which Hartigan--Wong's equilibrium advantage vanishes in massive clusters.

#theorem(title: "Update-Schedule Equilibrium Hierarchy and Asymptotic Boundary-Shell Decay", origin: "Original Unified Proof")[
  Let $cal(E)_"HW" (cal(X), k)$ and $cal(E)_"Lloyd" (cal(X), k)$ denote the sets of non-empty Hartigan--Wong equilibria and discrete Lloyd CVTs on $cal(X) subset.eq bb(R)^d$.
  1. *Strict Inclusion & Local Escape Advantage:* $cal(E)_"HW" (cal(X), k) subset.eq cal(E)_"Lloyd" (cal(X), k)$ universally, and every $cal(C) in cal(E)_"HW"$ obeys the strict interior Voronoi separation $norm(x_i - mu_a)_2 < norm(x_i - mu_b)_2$ whenever $|C_a| >= 2$ and $x_i != mu_a$. On $cal(X)_alpha = {0, 2, 2 + alpha}$ ($alpha in (1, 2)$, $k = 2$), ${{0, 2}, {2 + alpha}} in cal(E)_"Lloyd" backslash cal(E)_"HW"$, and Hartigan--Wong reduces WCSS by factor $4 / alpha^2 -> 4^-$ in one single transfer.
  2. *Exact $O(n_"min"^(-1))$ Decay of the Hartigan--Wong Advantage Shell:* For any non-empty partition $cal(C) in cal(E)_"Lloyd"$ with minimum cluster size $n_"min" := min_j |C_j| >= 2$, a point $x_i in C_a$ can trigger an improving Hartigan--Wong transfer to $C_b$ *only if* its squared distance ratio falls inside the thin Voronoi-boundary shell
     $ 1 <= (norm(x_i - mu_b)_2^2) / (norm(x_i - mu_a)_2^2) < (n_a (n_b + 1)) / ((n_a - 1) n_b) <= 1 + 2 / (n_"min" - 1). $
     Geometrically, letting $h_(a b) = norm(mu_a - mu_b)_2$ and letting $t_i = chevron.l x_i - (mu_a + mu_b)/2, (mu_b - mu_a)/h_(a b) chevron.r <= 0$ be the signed distance from $x_i$ to the Voronoi bisector $H_(a b)$, any HW-transferable Lloyd point satisfies
     $ 0 <= -t_i < (norm(x_i - mu_a)_2^2) / (h_(a b) (n_"min" - 1)). $
     Thus for continuous data densities with bounded support, the fraction of points in the escape shell decays as $O(n_"min"^(-1))$ as $n_"min" -> infinity$.
  3. *Monotonicity vs. Stochastic Variance Floor (Mini-Batch Alg. 3 vs. Batch Methods):* Lloyd, Elkan, Hamerly, and Hartigan--Wong are *strictly monotone* Lyapunov descent methods with finite termination. By contrast, Mini-Batch K-Means (Alg. 3) with batch size $b < n$ is *strictly non-monotone* on $J$ (with probability $1$ on non-degenerate clusters of diameter $> 0$), achieving an $O(n / b)$ per-iteration speedup while sustaining a local per-center MSE around an interior CVT $bold(mu)^*$ scaling as $bb(E)[norm(mu_j^((t)) - mu_j^*)_2^2] = Theta(sigma_j^2 / (b p_j t))$ under decaying counts $v_j^((t)) approx b p_j t$ @sculley2010webscale (where $p_j = (|C_j^*|) / n$ and $sigma_j^2$ is cluster variance).
]

#proof[
  (1) If $cal(C) in cal(E)_"HW"$ and $n_a >= 2$, equilibrium states $(n_a / (n_a - 1)) norm(x_i - mu_a)_2^2 <= (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2$. Because $n_a / (n_a - 1) > 1 > n_b / (n_b + 1)$, this implies $norm(x_i - mu_a)_2^2 <= norm(x_i - mu_b)_2^2$, strictly whenever $x_i != mu_a$, proving $cal(E)_"HW" subset.eq cal(E)_"Lloyd"$. On ${0, 2, 2 + alpha}$ with $C_1 = {0, 2}$ ($mu_1 = 1$, $W = 2$) and $C_2 = {2 + alpha}$ ($mu_2 = 2 + alpha$), $|2 - 1| = 1 < alpha = |2 - (2 + alpha)|$ confirms strict Lloyd stationarity, whereas $R_1 (2) = 2 (1)^2 = 2 > alpha^2 / 2 = I_2 (2)$ triggers an immediate transfer to ${{0}, {2, 2 + alpha}}$ of WCSS $alpha^2 / 2$.

  (2) Suppose $cal(C) in cal(E)_"Lloyd"$ so $norm(x_i - mu_b)_2^2 >= norm(x_i - mu_a)_2^2 > 0$, yet $x_i$ improves $J^*$ under a Hartigan--Wong transfer $a -> b$. Then $I_b (x_i) < R_a (x_i)$, which rearranges to
  $ (norm(x_i - mu_b)_2^2) / (norm(x_i - mu_a)_2^2) < (n_a / (n_a - 1)) / (n_b / (n_b + 1)) = (n_a (n_b + 1)) / ((n_a - 1) n_b). $
  Expand the right-hand side:
  $ (n_a (n_b + 1)) / ((n_a - 1) n_b) &= (1 + 1 / (n_a - 1)) (1 + 1 / n_b) \
    &= 1 + 1 / (n_a - 1) + 1 / n_b + 1 / ((n_a - 1) n_b). $
  Because each term on the right-hand side is strictly decreasing in $n_a$ and $n_b$, substituting $n_a, n_b >= n_"min" >= 2$ yields the exact bound $(n_a - 1)^(-1) + n_b^(-1) + ((n_a - 1) n_b)^(-1) <= (n_"min" + (n_"min" - 1) + 1) / (n_"min" (n_"min" - 1)) = 2 / (n_"min" - 1)$, attained with equality at $n_a = n_b = n_"min"$. Next, apply the Euclidean bisector identity $norm(x_i - mu_b)_2^2 - norm(x_i - mu_a)_2^2 = 2 h_(a b) (-t_i)$: substituting $norm(x_i - mu_b)_2^2 - norm(x_i - mu_a)_2^2 < (2 / (n_"min" - 1)) norm(x_i - mu_a)_2^2$ and dividing by $2 h_(a b)$ yields $0 <= -t_i < norm(x_i - mu_a)_2^2 / (h_(a b) (n_"min" - 1))$. Because every $x_i$ lies within bounded distance $R$ of $mu_a$, the escape region is contained inside a boundary slab of thickness $O(R^2 / (h_(a b) n_"min")) = O(n_"min"^(-1))$.

  (3) For Mini-Batch K-Means (Alg. 3), initialize at the exact global optimum $bold(mu)^*$ of any dataset containing a cluster $C_j^*$ with two distinct points $x_1 != x_2$. Drawing a batch $B$ containing $x_1$ and not $x_2$ shifts $mu_j$ away from $mu_j^* = overline(x)_(C_j^*)$, strictly *increasing* $J$ by Lemma 1 of Alg. 1 (parallel-axis identity). Inside an interior Voronoi neighborhood where assignments remain constant ($c_i equiv c_i^*$), the cumulative sequential update over $t$ independent batches of size $b$ with $N_j (t) = sum_(s=1)^t |B_s inter C_j^*|$ produces the exact sample mean of $N_j (t)$ i.i.d. draws from $C_j^*$; conditioning on $N_j (t) >= 1$ gives $bb(E)[norm(mu_j^((t)) - mu_j^*)_2^2 | N_j (t)] = sigma_j^2 / N_j (t)$, and by concentration of $N_j (t) tilde.op "Binomial"(b t, p_j)$ around $b p_j t$, the expected MSE scales as $Theta(sigma_j^2 / (b p_j t))$ while costing $O(b k d)$ per step instead of $O(n k d)$.
]

Theorem 4 resolves a classical empirical debate: Hartigan--Wong (Alg. 19) is strictly superior to Lloyd (Alg. 1) on small-to-medium datasets or partitions containing small clusters ($n_"min" <= 50$), where the $2 / (n_"min" - 1)$ margin is substantial, whereas on million-point balanced datasets ($n_"min" >= 10^4$) the escape shell collapses to $< 0.02\%$ thickness and Elkan/Hamerly-accelerated Lloyd (Algs. 13--14) or Mini-Batch K-Means (Alg. 3) wins decisively on wall-clock efficiency. Note also that Hartigan--Wong itself still has unbounded global worst-case approximation ratio on $k >= 3$ (the 4-point singleton-locked trap of Alg. 19), so pairing K-Means++ seeding (Alg. 2) with Hartigan--Wong or Hamerly refinement combines $O(log k)$ global protection with crisp local convergence.

#sec("III", "Metric Generalizations, Discrete Exemplars, and Breakdown-Point Hierarchy")

Half of Category I moves beyond unweighted squared Euclidean centroids to handle arbitrary dissimilarities, discrete exemplars, directional text embeddings, nonlinear manifolds, categorical/mixed tables, and heavy-tailed contamination. We compare twelve algorithms directly: K-Means (Alg. 1), PAM (Alg. 4) @kaufman1990finding, CLARA (Alg. 5) @srikantha2011clara, CLARANS (Alg. 6) @ng2002clarans, K-Modes (Alg. 7) @huang1998extensions, K-Prototypes (Alg. 8) @huang1998extensions, Spherical K-Means (Alg. 12) @dhillon2001concept, K-Medians (Alg. 15) @bradley1997clustering, Kernel K-Means (Alg. 16) @dhillon2004kernel, FastPAM / FasterPAM (Alg. 24) @schubert2019fasterpam @schubert2021fast, CLARA-MED (Alg. 25), and Trimmed K-Means (Alg. 26) @cuesta1997trimmed.

#subsec("A", "Complete Finite-Sample Breakdown-Point Hierarchy")

A persistent folklore misconception holds that replacing squared Euclidean distance by $ell_1$ distance (K-Medians) or discrete medoids (PAM/FasterPAM/CLARA) yields a $50\%$ global breakdown point for all $k >= 2$. We disprove that claim uniformly and establish the exact four-tier robustness hierarchy governing all of Category I.

#definition(title: "Local vs. Global Finite-Sample Breakdown Points")[
  Given a clean finite dataset $X_"clean"$ of size $n - b$ inside a bounded ball of diameter $D_0 < infinity$, let $X_b (R)$ denote the dataset obtained by replacing $b$ arbitrary points with adversarial outliers at distance $>= R$ from $X_"clean"$. For a $k$-prototype estimator returning $bold(theta)(X_b (R)) = (theta_1, dots, theta_k)$:
  - The *global breakdown point* $epsilon_k^*$ is the smallest fraction $b / n$ such that $sup_R max_(1 <= j <= k) "dist"(theta_j, X_"clean") = infinity$ (or, for scoring-based estimators such as CLARA-MED, $sup_R J_"score" (M_"clean") = infinity$).
  - The *intra-cluster (1-prototype) breakdown point* $epsilon_1^*$ is the breakdown point at $k = 1$ (equivalently, holding a cluster's membership fixed).
]

#theorem(title: "Unified Four-Tier Breakdown Hierarchy of Category I", origin: "Original Unified Proof")[
  Let $X_"clean"$ have diameter $D_0 < infinity$ and $k >= 2$ (with clean cluster sizes $n_1^*, dots, n_k^*$). Under arbitrary unbounded contamination $X_b (R)$ ($R -> infinity$):
  1. *Tier 0 — Total Breakdown ($epsilon_1^* = epsilon_k^* = 1 / n$):* For K-Means (Alg. 1) and every untrimmed squared-Euclidean method (Algs. 2, 3, 9, 13, 14, 17, 19, 21, 22, 23, 27, 30), a single outlier ($b = 1$) drives its assigned centroid to infinity even at $k = 1$: $norm(mu_j)_2 >= (R - D_0) / n -> infinity$.
  2. *Tier 1 — Bounded-Domain / Directional Shielding (Zero Radial Divergence, Bounded Deflection):*
     - In *Spherical K-Means* (Alg. 12), normalizing $x |-> x / norm(x)_2 in bb(S)^(d-1)$ forces $norm(mu_j)_2 = 1$ identically; a single adversarial direction deflects a clean cluster of resultant norm $R_j = norm(sum_(i in C_j) x_i)_2 > 1$ by angle at most $arcsin(1 / R_j)$.
     - In *Kernel K-Means* (Alg. 16) with any bounded continuous kernel $K(x, x) <= M < infinity$ vanishing at infinity ($K(x, y) -> 0$ as $norm(y)_2 -> infinity$, e.g., Gaussian RBF), every RKHS centroid stays inside the compact ball $norm(mu_j)_cal(H) <= sqrt(M)$; an outlier at $R -> infinity$ decouples orthogonally from $X_"clean"$ (contribution $O(1 / (|C_j|))$ to clean inner products), though at $k >= 2$ it can still hijack a cluster slot if the clean clusters have internal RKHS dispersion exceeding $M$.
     - In *K-Modes* (Alg. 7) on finite discrete alphabets, prototypes remain inside the finite product alphabet; a mode along attribute $r$ resists up to $floor((f_((1), r) - f_((2), r) - 1) / 2)$ adversarial substitutions (where $f_((1), r), f_((2), r)$ are the top two category counts).
  3. *Tier 2 — Local 50% Robustness with Global $1 / n$ Center-Hijacking ($epsilon_1^* = ceil(n / 2) / n -> 1 / 2$, yet $epsilon_k^* = 1 / n$ for $k >= 2$):* Let $J(bold(theta)) = sum_(i=1)^n min_(1 <= j <= k) d(x_i, theta_j)$ be *any* untrimmed additive metric distortion — encompassing *K-Medians* (Alg. 15, $ell_1$), *PAM / FastPAM / FasterPAM* (Algs. 4, 24), *CLARA* (Alg. 5), and *CLARANS* (Alg. 6). Then:
     - At $k = 1$, $epsilon_1^* = ceil(n / 2) / n -> 1 / 2$.
     - For *every* $k >= 2$, placing a *single* outlier $x_"out"$ at metric distance $R > n D_0$ forces *every* global minimizer $bold(theta)^*$ to allocate a prototype to $x_"out"$ ($"dist"(theta_(j^*), X_"clean") -> infinity$), yielding $epsilon_k^* = 1 / n$ — yet unlike Tier 0, all remaining $k - 1$ prototypes remain strictly inside $"conv"(X_"clean")$ (or $X_"clean"$ for medoids).
     - *K-Prototypes* (Alg. 8) exhibits *asymmetric breakdown*: $epsilon_1^* = 1 / n$ on numeric attributes and bounded finite-alphabet resistance on categorical attributes.
  4. *Tier 3 — True Positive Global Breakdown ($epsilon_k^* > 0$):*
     - *Trimmed K-Means* (Alg. 26) with trimming fraction $alpha in (0, 1)$ ($h = n - floor(alpha n)$) achieves exact positive global breakdown against $b <= floor(alpha n)$ outliers whenever clean inter-cluster separation prevents collapsing two clean clusters into one while sacrificing clean points; specifically, if $b <= floor(alpha n)$ and $b$ is smaller than the minimum clean cluster margin, all $k$ centroids remain uniformly bounded inside $"conv"(X_"clean")$ as $R -> infinity$, whereas $b >= floor(alpha n) + 1$ drives at least one active residual — and hence its centroid — to $infinity$ (giving sharp upper ceiling $epsilon_k^* <= (floor(alpha n) + 1) / n -> alpha$).
     - *CLARA-MED* (Alg. 25) achieves exact $50\%$ *scoring* breakdown $epsilon_"CLARA-MED"^* = ceil(n / 2) / n -> 1 / 2$ for *every* $k >= 2$ ($J_"med" (M_"clean") <= D_0$ uniformly in $R$ iff $b <= ceil(n / 2) - 1$), at the exact dual cost of *minority-cluster masking* for any true cluster of relative mass $< 1 / 2$.
]

#proof[
  (1) For $k = 1$, shifting $x_1$ to $R u$ ($norm(u)_2 = 1$) moves $mu = n^(-1)(R u + sum_(i=2)^n x_i)$ so $norm(mu - x_2)_2 >= (R - n D_0) / n -> infinity$. For $k >= 2$, the outlier belongs to some cluster $C_j$ of size $1 <= |C_j| <= n$, driving $mu_j$ at least $(R - |C_j| D_0) / (|C_j|) -> infinity$.

  (2) For Spherical K-Means, let $S = sum_(i in C_j) x_i$ with $R_j = norm(S)_2 > 1$ and unit direction $u = S / R_j$. Adding one unit outlier $v in bb(S)^(d-1)$ gives new resultant $S' = S + v$. The sine of the angle $theta$ between $S$ and $S'$ is maximized when $S'$ is tangent to the unit sphere centered at $S$, i.e., $v perp S'$, yielding $sin theta <= norm(v)_2 / norm(S)_2 = 1 / R_j$, hence $theta <= arcsin(1 / R_j)$. For Kernel K-Means with $K(x, x) <= M$, every centroid $mu_j = |C_j|^(-1) sum_(i in C_j) phi(x_i)$ satisfies $norm(mu_j)_cal(H) <= |C_j|^(-1) sum_(i in C_j) sqrt(K(x_i, x_i)) <= sqrt(M)$ by Minkowski's inequality. For K-Modes, flipping $b$ entries on coordinate $r$ decreases $f_((1), r)$ by at most $b$ and increases $f_((2), r)$ by at most $b$; hence the original mode remains unique whenever $f_((1), r) - b > f_((2), r) + b$, i.e., $b < (f_((1), r) - f_((2), r)) / 2$.

  (3) For $k = 1$ under any metric $d$ (or $ell_1$) with $b < n / 2$ outliers, fix any clean point $theta_0 in X_"clean"$ and consider any candidate prototype $theta_"far"$ at distance $L := min_(x in X_"clean") d(theta_"far", x)$ from $X_"clean"$. By the triangle inequality, every outlier $x in X_"out"$ satisfies $d(x, theta_0) <= d(x, theta_"far") + d(theta_"far", theta_0) <= d(x, theta_"far") + L + D_0$. Summing over $X_"out"$ (size $b$) and using $d(x, theta_0) <= D_0$, $d(x, theta_"far") >= L$ on $X_"clean"$ (size $n - b$) gives
  $ J(theta_"far") - J(theta_0) &>= (n - b) L - (n - b) D_0 - b (L + D_0) \
    &= (n - 2 b) L - n D_0. $
  Whenever $b < n / 2$, $n - 2 b >= 1$, so any global minimizer must obey $L <= (n D_0) / (n - 2 b) <= n D_0$ *uniformly for all $R$ and all outlier configurations*! Conversely, placing $b = ceil(n / 2)$ coincident outliers at distance $R -> infinity$ forces $J(theta_"out") = (n - b) R <= b R - epsilon < J(theta_0)$ (for odd $n$, or perturbed by $epsilon$ for even $n$), proving $epsilon_1^* = ceil(n / 2) / n$ simultaneously for K-Medians and all $k$-medoids methods!

  Now let $k >= 2$ and move a single point $x_"out"$ ($b = 1$) to distance $R > n D_0$ from $X_"clean"$. Any $k$-prototype set $bold(theta)$ with $min_j d(x_"out", theta_j) >= R / 2$ incurs $J(bold(theta)) >= R / 2$. Conversely, placing $theta_1 = x_"out"$ and $theta_2, dots, theta_k in X_"clean"$ pays $0$ for $x_"out"$ and at most $D_0$ for each of the $n - 1$ clean points, yielding $J(bold(theta)') <= (n - 1) D_0 < R / 2$ whenever $R > 2 n D_0$. Thus every global minimizer must place a prototype within $R / 2$ of $x_"out"$ (distance $>= R / 2 -> infinity$ from $X_"clean"$), proving $epsilon_k^* = 1 / n$ uniformly across K-Medians, PAM, FastPAM/FasterPAM, CLARA, and CLARANS! Yet every non-hijacked cluster $C_j subset.eq X_"clean"$ has its prototype minimizing $sum_(i in C_j) d(x_i, theta_j)$, keeping all $k - 1$ remaining prototypes inside $"conv"(X_"clean")$ (or $X_"clean"$).

  (4) For Trimmed K-Means (Alg. 26) with $h = n - floor(alpha n)$, if $b >= floor(alpha n) + 1$, then by the pigeonhole principle *every* $h$-element active subset $A in cal(A)_h$ contains at least one outlier $x_"out" in X_"out"$. If all $k$ centroids remained within distance $B$ of $X_"clean"$, the active objective would exceed $(R - B)^2 / h -> infinity$, whereas placing one centroid on $x_"out"$ keeps $J_alpha <= D_0^2$; hence at least one centroid escapes to $infinity$, proving $epsilon_k^* <= (floor(alpha n) + 1) / n$. Conversely, when $b <= floor(alpha n)$, the clean subset $X_"clean"$ has cardinality $n - b >= h$; choosing $A subset.eq X_"clean"$ and $k$ centroids in $X_"clean"$ achieves $J_alpha <= D_0^2$ uniformly in $R$. If $b$ is strictly smaller than the smallest clean cluster margin so that sacrificing clean points to capture an outlier cannot lower the $h$-point WCSS below the clean $k$-cluster optimum (Theorem 2 of Alg. 26), all $k$ centroids stay inside $"conv"(X_"clean")$ for all $R$. Finally, for CLARA-MED (Alg. 25), let $u = floor(n / 2) + 1$: if $b <= n - u = ceil(n / 2) - 1$, then at least $u$ points lie in $X_"clean"$ with $d(x_i, M_"clean") <= D_0$, bounding the median $J_"med" (M_"clean") <= D_0$ independently of $R$; if $b >= ceil(n / 2)$, then $d_((u)) (M_"clean") >= R$, driving $J_"med" (M_"clean") >= R / 2 -> infinity$. (Meanwhile, by Theorem 3 of Alg. 25, on a dataset of $3 m$ majority points and $m - 1$ distant minority points, $J_"med"$ assigns both medoids to the majority cluster with $J_"med" = 0$, completely masking the minority cluster!)
]

#wide_figure(
  caption: [*Table II:* Four-Tier Finite-Sample Breakdown & Robustness Comparison Across Prototype Paradigms (Algs. 1, 4--8, 12, 15, 16, 24--26).],
  table(
    columns: (0.65fr, 1.85fr, 0.75fr, 0.75fr, 2.0fr),
    align: (left, left, center, center, left),
    [Tier], [Algorithms], [$epsilon_1^*$ (Local)], [$epsilon_k^*$ ($k >= 2$)], [Outlier Behavior & Trade-Off],
    [Tier 0], [Algs. 1--3, 9, 13--14, 17, 19, 21--23, 27, 30], [$1 / n$], [$1 / n$], [Continuous centroid drag to $infinity$],
    [Tier 1], [Algs. 7, 12, 16], [Bounded], [Bounded\*], [Unit-sphere / RKHS / finite alphabet shield],
    [Tier 2], [Algs. 4, 5, 6, 15, 24], [$-> 1 / 2$], [$1 / n$], [1 outlier hijacks 1 slot; other $k - 1$ safe],
    [Tier 2\*], [Alg. 8 (K-Prototypes)], [$1 / n$ num.], [$1 / n$], [Asymmetric: numeric fragile, cat. safe],
    [Tier 3], [Alg. 26 (Trimmed)], [$-> alpha$], [$<= alpha$], [Positive global brk.; over-trim erasure risk],
    [Tier 3\*], [Alg. 25 (CLARA-MED)], [$-> 1 / 2$], [$-> 1 / 2$], [50% scoring brk.; masks $< 50\%$ clusters],
  ),
)

#subsec("B", "Discrete Exemplar Gaps, Johnson-Graph Scaling & Non-Euclidean Isometries")

Beyond robustness, Category I's non-Euclidean and discrete-exemplar algorithms obey tight quantitative distortion bridges that explain when simpler Euclidean surrogates suffice versus when specialized objectives are indispensable.

#theorem(title: "Tight Factor-of-2 Discrete/Categorical Bridges and Johnson-Graph Hierarchy", origin: "Original Unified Proof")[
  1. *Discrete-to-Continuous Factor-2 Integrality Gap (Medoids Algs. 4, 24 vs. K-Means Alg. 1 & K-Medians Alg. 15):*
     - In any metric space $(cal(M), d)$ and finite $X subset.eq cal(M)$, the optimal $k$-medoids cost $"TD"_"medoid"^*$ and continuous $k$-median cost $"TD"_"median"^* = min_(bold(theta) in cal(M)^k) sum_(x in X) min_j d(x, theta_j)$ satisfy
       $ "TD"_"median"^* <= "TD"_"medoid"^* <= 2 "TD"_"median"^*, $
       tight as $n -> infinity$ on star metrics.
     - Under squared Euclidean distance on $X subset.eq bb(R)^d$, restricting prototypes to $M subset.eq X$ satisfies the identical tight factor-of-2 bound:
       $ J_"means"^* <= J_"sq-medoid"^* <= 2 J_"means"^*, $
       tight on the vertices of a regular $(n - 1)$-simplex ($k = 1$).
  2. *One-Hot Factor-of-2 Sandwich for K-Modes (Alg. 7) & K-Prototypes (Alg. 8):* Map each mixed observation $x_i = (x_i^((r)), x_i^((c))) in bb(R)^(d_r) times product_(q=1)^(d_c) cal(A)_q$ via scaled one-hot encoding $psi(x_i) := (x_i^((r)), sqrt(gamma / 2) e_(x_(i 1)^((c))), dots, sqrt(gamma / 2) e_(x_(i d_c)^((c))))$. Then for *every* partition $cal(C)$, the exact K-Prototypes objective $J_"KP" (cal(C))$ (specializing to K-Modes at $d_r = 0, gamma = 1$) and the continuous Euclidean K-Means WCSS $J_"euc" (cal(C))$ on $psi(cal(X))$ obey
     $ J_"euc" (cal(C)) <= J_"KP" (cal(C)) <= 2 J_"euc" (cal(C)), $
     with the upper constant $2$ attained exactly by uniform binary attributes.
  3. *Spherical (Alg. 12) vs. Euclidean (Alg. 1) vs. Kernel K-Means (Alg. 16):* On unit-normalized data $cal(X) subset.eq bb(S)^(d-1)$, maximizing Spherical K-Means cosine similarity $sum_j sum_(i in C_j) chevron.l x_i, hat(mu)_j chevron.r$ ($norm(hat(mu)_j)_2 = 1$) equals maximizing $sum_(j=1)^k norm(sum_(i in C_j) x_i)_2$, whereas Euclidean K-Means maximizes $sum_(j=1)^k |C_j|^(-1) norm(sum_(i in C_j) x_i)_2^2$ (strictly distinct partitions whenever cluster sizes or concentrations differ!). Both are exact specializations of Kernel K-Means (Alg. 16) in RKHS $cal(H)$ via $K_"lin" (x, y) = chevron.l x, y chevron.r$ and $K_"cos" (x, y) = (chevron.l x, y chevron.r) / (norm(x)_2 norm(y)_2)$, and for any weighted graph $G = (V, W)$ with degree matrix $D$, Normalized Cut @shi2000ncut is algebraically equivalent to weighted Kernel K-Means with $K = sigma D^(-1) + D^(-1) W D^(-1)$ @dhillon2004kernel.
  4. *Exact Johnson-Graph Comparison of Medoid Solvers (Algs. 4, 5, 6, 24, 25):* Define the Johnson graph $J(n, k)$ whose vertices are $k$-subsets $M subset.eq X$ linked whenever $|M inter M'| = k - 1$ (degree $k(n - k)$).
     - *PAM (Alg. 4) and FastPAM1 (Alg. 24)* trace the *exact same* steepest-descent trajectory on $J(n, k)$ to a 1-swap local optimum guaranteed to satisfy $"TD"(M) <= 5 "TD"(M^*)$ in metric spaces @kanungo2004local @schubert2021fast, while FastPAM1 slashes per-pass cost from $O(k (n - k)^2)$ to $O(n(n - k))$ via simultaneous delta decomposition $Delta_"TD" (m_i, x_c) = Delta_i^- + Delta^+ (x_c) + gamma_i (x_c)$; *FasterPAM* traverses eager first-improvement edges on $J(n, k)$ with $O(n)$ balanced cache updates, converging to the same certified 5-approximation class in $O(1)$ passes.
     - *CLARANS (Alg. 6)* performs randomized bounded-degree walks on the full $J(n, k)$ ($m_"max"$ neighbors sampled per step, cost $O(r m_"max" n)$), escaping subsample blindness but forfeiting the 5-approximation certificate whenever $m_"max" < k(n - k)$.
     - *CLARA (Alg. 5) and CLARA-MED (Alg. 25)* restrict local search to $q$ random induced subgraphs $J(s, k)$ ($O(q(s^2 T + n k))$), where classical $s = 40 + 2 k$ fails with coverage probability $< e^(s / (k - 2)) / k -> 0$ and $s = Omega(k log k)$ is both necessary and sufficient for constant-probability cluster coverage.
]

#proof[
  (1) Let $cal(C)^* = {C_1^*, dots, C_k^*}$ and $bold(theta)^*$ attain $"TD"_"median"^*$. For each cluster $C_j^*$ ($n_j = |C_j^*|$), choosing a random exemplar $M_j in C_j^*$ uniformly and applying the triangle inequality gives
  $ min_(m in C_j^*) sum_(x in C_j^*) d(x, m) &<= 1 / n_j sum_(y in C_j^*) sum_(x in C_j^*) d(x, y) \
    &<= 1 / n_j sum_(y in C_j^*) sum_(x in C_j^*) (d(x, theta_j^*) + d(y, theta_j^*)) \
    &= 2 sum_(x in C_j^*) d(x, theta_j^*). $
  Summing over $j = 1, dots, k$ proves $"TD"_"medoid"^* <= 2 "TD"_"median"^*$. Tightness holds on a star metric with $n$ leaves at distance $1$ from a virtual hub ($k = 1$): $"TD"_"median"^* = n$ while every leaf pays $2(n - 1)$, ratio $2(1 - 1/n) -> 2$. For squared Euclidean distance, the parallel-axis identity yields for each $C_j^*$ with mean $mu_j^*$:
  $ 1 / n_j sum_(y in C_j^*) sum_(x in C_j^*) norm(x - y)_2^2 &= 1 / n_j sum_(y in C_j^*) ( sum_(x in C_j^*) norm(x - mu_j^*)_2^2 + n_j norm(y - mu_j^*)_2^2 ) \
    &= 2 sum_(x in C_j^*) norm(x - mu_j^*)_2^2. $
  Since the minimum over $y in C_j^*$ cannot exceed the average, $J_"sq-medoid"^* <= 2 J_"means"^*$. On the $n$ vertices of a regular simplex inscribed in the unit sphere ($norm(x_i)_2 = 1$, $chevron.l x_i, x_r chevron.r = -1/(n - 1)$), $J_"means"^* = n$ while every vertex has $J_"sq-medoid" = (n - 1) (2 + 2/(n-1)) = 2 n$, attaining factor $2$ identically!

  (2) Fix any cluster $C_j$ of size $n_j$ and any single categorical attribute $q$ with empirical category proportions $p_v = n_j^(-1) |{i in C_j : x_(i q)^((c)) = v}|$ and maximum frequency $p_"max" = max_v p_v$. Under K-Prototypes with weight $gamma$, the optimal mode incurs exact categorical cost $n_j gamma (1 - p_"max")$. Under Euclidean K-Means on the scaled one-hot block $sqrt(gamma / 2) e_(x_(i q)^((c)))$, the continuous centroid is $sqrt(gamma / 2) p$, whose exact WCSS is
  $ sum_(i in C_j) gamma / 2 norm(e_(x_(i q)^((c))) - p)_2^2 = n_j gamma / 2 (1 - sum_v p_v^2). $
  Observe the elementary two-sided Gini--mode inequality: because $p_"max"^2 <= sum_v p_v^2 <= p_"max" sum_v p_v = p_"max"$, we have
  $ 1 / 2 (1 - p_"max") &<= 1 / 2 (1 - sum_v p_v^2) <= 1 / 2 (1 - p_"max"^2) \
    &= (1 + p_"max") / 2 (1 - p_"max") <= 1 - p_"max". $
  Multiplying by $n_j gamma$, summing over all categorical attributes $q = 1, dots, d_c$, and adding the identical numeric WCSS $sum_(i in C_j) norm(x_i^((r)) - mu_j^((r)))_2^2$ (for which $1 <= 1 <= 2$ trivially holds) proves $J_"euc" (cal(C)) <= J_"KP" (cal(C)) <= 2 J_"euc" (cal(C))$ for every partition $cal(C)$. Equality $J_"KP" = 2 J_"euc"$ holds when $d_r = 0$ and $p = (1/2, 1/2)$ ($1 - p_"max" = 1/2 = 2 dot 1/4$).

  (3) For $x_i in bb(S)^(d-1)$, $sum_(i in C_j) chevron.l x_i, hat(mu)_j chevron.r = chevron.l R_j, hat(mu)_j chevron.r <= norm(R_j)_2$ (where $R_j = sum_(i in C_j) x_i$), with equality uniquely at $hat(mu)_j = R_j / norm(R_j)_2$. Meanwhile, Euclidean K-Means satisfies $sum_j sum_(i in C_j) norm(x_i - overline(x)_(C_j))_2^2 = n - sum_(j=1)^k |C_j|^(-1) norm(R_j)_2^2$. Because maximizing $sum_j norm(R_j)_2$ vs. $sum_j |C_j|^(-1) norm(R_j)_2^2$ weights resultant norms linearly vs. quadratically divided by $|C_j|$, they disagree whenever cluster densities and sizes are imbalanced (for instance, on $bb(S)^1$, moving a boundary point between a large diffuse cluster and a tight small cluster changes $norm(R_j)_2$ and $|C_j|^(-1) norm(R_j)_2^2$ with opposite signs). The RKHS identities follow directly from $norm(phi(x) - mu_j)_cal(H)^2 = K(x, x) - 2 |C_j|^(-1) sum_(r in C_j) K(x, x_r) + |C_j|^(-2) sum_(r, s in C_j) K(x_r, x_s)$ and Dhillon--Guan--Kulis @dhillon2004kernel.

  (4) The exact algebraic identity $Delta_"TD" (m_i, x_c) = Delta_i^- + Delta^+ (x_c) + gamma_i (x_c)$ (proved across all five cases in Alg. 24) establishes bit-for-bit trajectory identity of FastPAM1 with PAM on $J(n, k)$ at $O(n(n - k))$ cost per pass; FasterPAM executes eager swaps along strictly descending edges of $J(n, k)$ and halts only when a complete pass verifies $Delta_"TD" (m_i, x_c) >= 0$ for all $k(n - k)$ neighbors, certifying the identical metric 5-approximation @kanungo2004local. The CLARA / CLARA-MED without-replacement coverage ceiling $bb(P)[forall j, S_r inter C_j^* != emptyset] <= e^(s / (k - 2)) / k$ follows from negative covariance of cluster-miss indicators and Chebyshev's inequality (Theorem 2 of Alg. 25), proving necessity of $s = Omega(k log k)$.
]

#sec("IV", "Model Selection, Divisive Splitting, Feature Weighting, and Constraints")

The remaining eight algorithms of Category I address four fundamental structural extensions of classical $k$-means: *automatic estimation of $k$* (X-Means Alg. 10 @pelleg2000xmeans, G-Means Alg. 11 @hamerly2003gmeans, ISODATA Alg. 18 @ball1965isodata), *hierarchical divisive quantization* (Bisecting K-Means Alg. 9 @steinbach2000comparison and LBG Alg. 30 @linde1980algorithm), *automated feature/subspace weighting* (W-K-Means Alg. 28 @huang2005automated and Minkowski Weighted K-Means Alg. 29 @deamorim2012minkowski), and *semi-supervised instance constraints* (COP-KMeans Alg. 20 @wagstaff2001constrained).

#subsec("A", "Model-Selection & Divisive Splitting Consistency Trichotomy")

X-Means, G-Means, and ISODATA all promise to eliminate manual specification of $k$, yet their statistical decision criteria behave in sharply opposite ways across anisotropic, non-Gaussian, and spherical geometries. We prove an exact analytical trichotomy separating all five dynamic/divisive methods.

#theorem(title: "Exact Model-Selection Trichotomy: X-Means vs. G-Means vs. ISODATA vs. Bisecting/LBG", origin: "Original Unified Proof")[
  1. *Sharp Anisotropic Over-Splitting Threshold of X-Means (Alg. 10):* Let $X_1, dots, X_n$ be i.i.d. draws from a *single* $d$-dimensional Gaussian cluster $cal(N)(bold(0), Sigma)$ ($d >= 2$) with $Sigma = "diag"(lambda_1, lambda_2, dots, lambda_2)$ and aspect ratio $kappa := lambda_1 / lambda_2 >= 1$. Under Pelleg--Moore's spherical BIC criterion @pelleg2000xmeans @schwarz1978bic:
     - When $kappa = 1$ (isotropic), X-Means *never* splits as $n -> infinity$ ($(Delta "BIC") / n -> - (d / 2) ln(1 - 2 / (pi d)) - ln 2 < 0$ for all $d >= 1$).
     - For every $d >= 2$, splitting along the principal axis asymptotically increases spherical BIC with probability $1$ whenever
       $ kappa > kappa_"crit" (d) := (d - 1) (1 - 4^(-1 / d)) / (2 / pi - (1 - 4^(-1 / d))). $
       In particular, $kappa_"crit" (2) approx 3.66$ and $lim_(d -> infinity) kappa_"crit" (d) = pi ln 2 approx 2.178$; thus even mild $2.2 times$ elongation in high dimensions forces X-Means to over-split a single Gaussian cluster!
  2. *Exact Reversal for G-Means (Alg. 11) — Anisotropic Immunity vs. Platykurtic Over-Splitting:* Under the identical elongated single Gaussian $cal(N)(bold(0), Sigma)$ for *arbitrary* $kappa in [1, infinity)$, every linear projection $v^T X$ is exact univariate Gaussian, so G-Means' Anderson--Darling test @anderson1952asymptotic @hamerly2003gmeans retains exact asymptotic false-split probability $alpha$ independently of $kappa$. Conversely, on a *single* uniform hypercube cluster $X tilde.op "Unif"([0, 1]^d)$ ($d >= 2$), X-Means asymptotic log-likelihood change is $-(d / 2) ln(1 - 3 / (4 d)) - ln 2 <= 3 / 8 - ln 2 + O(d^(-1)) < 0$ (never splitting for $d >= 2$), whereas G-Means rejects normality with probability $-> 1$ as $n -> infinity$ and recursively shatters the uniform cluster!
  3. *Split--Merge Cycling of ISODATA (Alg. 18):* Unlike X-Means and G-Means (whose divisive splits strictly decrease WCSS and terminate finitely), ISODATA alternates threshold-driven splits and merges without a global Lyapunov energy; whenever its merge distance threshold $theta_C$ and split standard-deviation threshold $theta_S$ violate the hysteresis deadband $theta_C < 2 sqrt(1 - 2 / pi) theta_S$ (or in general $theta_C < 2 delta_"split"$ where $delta_"split"$ is half the post-split centroid separation), ISODATA admits *infinite period-2 limit cycles*!
  4. *Global Voronoi Optimality Gap: Bisecting K-Means (Alg. 9) vs. LBG (Alg. 30):* Both construct $k$ clusters by binary splitting, but Bisecting K-Means freezes ancestor boundaries (cost $O(n d R T log k)$ for balanced trees), yielding non-Voronoi leaf partitions whose global WCSS can remain strictly suboptimal relative to one global Lloyd pass, whereas LBG executes full global GLA after each codebook doubling $m |-> 2 m$, guaranteeing global Lloyd--Max CVT optimality at cost $O(n k d T_"max")$ (with total splitting overhead bounded by $sum_(r=0)^(ceil(log_2 k) - 1) 2^r < 2 k$ passes).
]

#proof[
  (1) Under spherical covariance $sigma^2 I_d$, the single-cluster MLE variance is $hat(sigma)_1^2 -> d^(-1) "Tr"(Sigma) = d^(-1) (lambda_1 + (d - 1) lambda_2)$. Splitting $cal(N)(bold(0), Sigma)$ symmetrically at $x_1 = 0$ replaces coordinate 1's variance $lambda_1$ by the half-Gaussian conditional variance $"Var"(Z | Z > 0) lambda_1 = (1 - 2 / pi) lambda_1$ (since $bb(E)[Z | Z > 0] = sqrt(2 / pi)$ for $Z tilde.op cal(N)(0, 1)$), while coordinates $2, dots, d$ are independent of $x_1$ and retain variance $lambda_2$. Hence the pooled 2-cluster MLE variance converges almost surely to
  $ hat(sigma)_2^2 -> 1 / d ((1 - 2 / pi) lambda_1 + (d - 1) lambda_2) = hat(sigma)_1^2 (1 - 2 / pi dot kappa / (kappa + d - 1)). $
  For $k$ equal-weight spherical clusters, Pelleg--Moore's log-likelihood is $ell_k = -n ln k - (n d / 2) ln(2 pi hat(sigma)_k^2) - n d / 2$, and the BIC parameter penalty increases by $((p_2 - p_1) / 2) ln n = ((d + 1) / 2) ln n = o(n)$. Therefore,
  $ ("BIC"_2 - "BIC"_1) / n -> g(d, kappa) := - d / 2 ln(1 - 2 / pi dot kappa / (kappa + d - 1)) - ln 2. $
  At $kappa = 1$, letting $u = 2 / (pi d) in (0, 2 / pi]$, note that $-ln(1 - u) / u$ is strictly decreasing in $d$: at $d = 1$, $g(1, 1) = -1/2 ln(1 - 2/pi) - ln 2 approx -0.187 < 0$, and as $d -> infinity$, $g(d, 1) -> 1 / pi - ln 2 approx -0.375 < 0$, proving isotropic non-splitting for all $d >= 1$. Solving $g(d, kappa) > 0$ for $d >= 2$ gives
  $ 1 - 2 / pi dot kappa / (kappa + d - 1) < 4^(-1 / d) <==> 2 / pi dot kappa / (kappa + d - 1) > 1 - 4^(-1 / d). $
  Because $1 - 4^(-1 / d) <= 1 / 2 < 2 / pi$ for every $d >= 2$, the denominator $2 / pi - (1 - 4^(-1 / d))$ is strictly positive, isolating $kappa > kappa_"crit" (d)$ explicitly. Using $1 - 4^(-1 / d) = (ln 4) / d + O(d^(-2))$ yields $lim_(d -> infinity) kappa_"crit" (d) = (ln 4) / (2 / pi) = pi ln 2 approx 2.178$.

  (2) If $X tilde.op cal(N)(mu, Sigma)$ with arbitrary PSD $Sigma$, then for any non-zero direction $v in bb(R)^d$, $Z = v^T X tilde.op cal(N)(v^T mu, v^T Sigma v)$, whose standardized Anderson--Darling statistic $A_*^2$ converges to the universal null distribution of @anderson1952asymptotic independently of $kappa$, keeping G-Means' split probability at nominal level $alpha$ (up to finite-sample data-snooping of $v = mu_+ - mu_-$ analyzed in Alg. 11). Conversely, for $X tilde.op "Unif"([0, 1]^d)$ ($d >= 2$), splitting coordinate 1 at $1 / 2$ reduces its variance from $1 / 12$ to $1 / 48$ (a drop of $3 / 4$ on one coordinate), giving $hat(sigma)_2^2 / hat(sigma)_1^2 -> 1 - 3 / (4 d)$. At $d = 2$, $-(2 / 2) ln(1 - 3 / 8) - ln 2 = -ln(5 / 8) - ln 2 = ln(4 / 5) < 0$, and for $d >= 2$, $- (d / 2) ln(1 - 3 / (4 d))$ decreases monotonically toward $3 / 8 < ln 2$, so X-Means *never* splits $[0, 1]^d$ for $d >= 2$! Meanwhile, along any coordinate projection, the uniform CDF $F_"unif"$ differs from the Gaussian CDF $Phi$ with strictly positive Anderson--Darling distance $integral (F_"unif" (t) - Phi(t))^2 / (Phi(t)(1 - Phi(t))) dif Phi(t) > 0$, so $A^2 = Theta(n) -> infinity$ and G-Means splits $[0, 1]^d$ with probability $-> 1$ at every level!

  (3) Consider two balanced symmetric 1D sub-populations centered at $plus.minus c$ with tiny internal radius $r << c$, initially merged into one cluster of mean $0$ and standard deviation $sigma approx c$. If $sigma > theta_S$, ISODATA splits the cluster into two clusters at centroids $plus.minus c$ (separation $2 c approx 2 sigma$, or $2 delta_"split"$ for general distributions). If on the next step $2 c < theta_C$ (violating the hysteresis deadband), ISODATA's merge step immediately recombines them into the exact initial single cluster, generating an infinite period-2 cycle (Theorem 2 of Alg. 18).

  (4) In Bisecting K-Means, once a root 2-means hyperplane separates $cal(X)$ into $L$ and $R$, subsequent bisections operate strictly inside $L$ and $R$ independently; a point in $L$ near the original root cut can end up strictly closer to a newly created child centroid inside $R$ than to its own leaf centroid in $L$, violating the Voronoi condition. By contrast, LBG's inner GLA runs globally over all $m$ active codewords after every doubling step, enforcing both NNC and CC globally while summing geometric stage costs $sum_(r=1)^(ceil(log_2 k)) n dot 2^r d T_"max" <= 2 n k d T_"max" = O(n k d T_"max")$.
]

#subsec("B", "Feature Weighting & Subspace Expressiveness: W-K-Means vs. MWK-Means vs. K-Means")

High-dimensional tabular data often bury cluster structure in a low-dimensional signal subspace surrounded by noise coordinates, where unweighted K-Means (Alg. 1) suffers distance concentration @beyer1999nearest. Both W-K-Means (Alg. 28) @huang2005automated and Minkowski Weighted K-Means (Alg. 29) @deamorim2012minkowski combat noise via simplex-constrained coordinate weights, yet differ fundamentally in *global* versus *cluster-local* subspace expressiveness.

#theorem(title: "Exact Noise Suppression and Strict Subspace Separation: MWK-Means vs. W-K-Means vs. K-Means", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ and fix weight-regularization exponent $gamma > 1$ ($gamma = beta$ for W-K-Means, $gamma = p$ for MWK-Means with $epsilon_0 = 0$).
  1. *Unified Noise-Attenuation Law:* Suppose $d = d_S + d_N$ coordinates consist of $d_S$ signal features each with within-cluster $gamma$-dispersion $D_S > 0$ and $d_N$ noise features each with dispersion $D_N > D_S$. Under optimal Hölder weights $w_r^* prop D_r^(-1 / (gamma - 1))$, the effective *metric-weight* ratio and *weighted-dispersion* ratio of a single noise coordinate to a signal coordinate are identically
     $ ((w_N^*)^gamma) / ((w_S^*)^gamma) = (D_S / D_N)^(gamma / (gamma - 1)), quad ((w_N^*)^gamma D_N) / ((w_S^*)^gamma D_S) = (D_S / D_N)^(1 / (gamma - 1)). $
     Consequently, total noise energy is suppressed below signal energy whenever $d_N / d_S < (D_N / D_S)^(1 / (gamma - 1))$, whereas isotropic K-Means ($gamma -> infinity$) fails as soon as $d_N / d_S > D_S / D_N$!
  2. *Strict Cluster-Local Subspace Dominance of MWK-Means over W-K-Means:* There exists a standardized 2-cluster orthogonal-subspace dataset $cal(X) subset.eq bb(R)^2$ on which:
     - Global W-K-Means (Alg. 28, any $beta > 1$) is forced to assign trivial uniform weights $w_1^* = w_2^* = 1 / 2$, degenerating identically to isotropic K-Means (Alg. 1) with reduced distortion $(1 / 2)^(beta - 1) (D_"tight" + D_"wide")$;
     - Cluster-local MWK-Means (Alg. 29, $p = beta > 1$) assigns strictly asymmetric weights $w_(1 1)^* = w_(2 2)^* > 1 / 2 > w_(1 2)^* = w_(2 1)^*$ achieving strictly lower reduced objective
       $ J_"MWK"^* / J_"WK"^* = (2 (D_"tight"^(-1/(p-1)) + D_"wide"^(-1/(p-1)))^(1 - p)) / (2^(1 - p) (D_"tight" + D_"wide")) < 1, $
       with ratio vanishing as $(D_"tight" / D_"wide") -> 0$!
  3. *Shared Structural Fragility — Scale Collapse & Rotation Non-Invariance:* Both W-K-Means and MWK-Means lose $O(d)$ orthogonal invariance and suffer *inverse scale collapse*: rescaling a single coordinate $x_(i 1) |-> s x_(i 1)$ with $s -> 0^+$ drives $D_1 -> 0^+$ and concentrates $w_1^* -> 1$ onto the shrunk coordinate regardless of its clustering informativeness.
]

#proof[
  (1) By Lemma 1 of Alg. 28 and Lemma 1 of Alg. 29, minimizing $sum_(r=1)^d w_r^gamma D_r$ over $w in Delta_(d-1)$ for $gamma > 1$ via Hölder's inequality uniquely sets $w_r^* = D_r^(-1/(gamma - 1)) / Z$ ($Z = sum_q D_q^(-1/(gamma - 1))$). Dividing $w_N^* / w_S^* = (D_S / D_N)^(1 / (gamma - 1))$, raising to power $gamma$ gives the metric-weight ratio $(D_S / D_N)^(gamma / (gamma - 1))$, and multiplying by $D_N / D_S$ gives $(D_S / D_N)^(1 / (gamma - 1))$. Total noise contribution equals $d_N (w_N^*)^gamma D_N$, which falls strictly below total signal contribution $d_S (w_S^*)^gamma D_S$ iff $d_N / d_S < (D_N / D_S)^(1 / (gamma - 1))$.

  (2) Construct $cal(X) = C_1 union C_2 subset.eq bb(R)^2$ with equal cluster sizes $m$: let $C_1$ lie tightly along coordinate 1 with $p$-dispersion $D_(1 1) = D_"tight"$ and spread widely along coordinate 2 with $D_(1 2) = D_"wide" > D_"tight"$, while $C_2$ mirrors $C_1$ across the diagonal $x_1 = x_2$ so $D_(2 1) = D_"wide"$ and $D_(2 2) = D_"tight"$ (choosing cluster centers $(0, 0)$ and $(c, c)$ or symmetric locations so both coordinates have identical marginal mean and variance). For W-K-Means, the global coordinate dispersions are equal:
  $ D_1 = D_(1 1) + D_(2 1) = D_"tight" + D_"wide" = D_(1 2) + D_(2 2) = D_2. $
  Hence $w_1^* = w_2^* = 1 / 2$ identically for every $beta > 1$, yielding minimum objective $J_"WK"^* = 2 (1 / 2)^beta (D_"tight" + D_"wide") = 2^(1 - beta) (D_"tight" + D_"wide")$. Setting $p = beta > 1$, MWK-Means optimizes $w_1$ and $w_2$ independently per cluster, attaining
  $ J_"MWK"^* = 2 (D_"tight"^(-1 / (p - 1)) + D_"wide"^(-1 / (p - 1)))^(1 - p). $
  Let $r := D_"tight" / D_"wide" in (0, 1)$. Then
  $ J_"MWK"^* / J_"WK"^* = 2^p / (1 + r) dot (r^(-1 / (p - 1)) + 1)^(1 - p) = (2^p r) / ((1 + r) (1 + r^(1 / (p - 1)))^(p - 1)). $
  By strict power-mean ordering $M_(-1/(p-1)) (D_"tight", D_"wide") < M_1 (D_"tight", D_"wide")$ (valid whenever $D_"tight" != D_"wide"$ and $-1/(p - 1) < 0 < 1$):
  $ ( 1 / 2 (D_"tight"^(-1/(p-1)) + D_"wide"^(-1/(p-1))) )^(1 - p) < 1 / 2 (D_"tight" + D_"wide"). $
  Multiplying both sides by $2 dot 2^(1 - p)$ gives directly
  $ 2 (D_"tight"^(-1/(p-1)) + D_"wide"^(-1/(p-1)))^(1 - p) < 2^(1 - p) (D_"tight" + D_"wide")! $
  Moreover, as $r = D_"tight" / D_"wide" -> 0^+$, $J_"MWK"^* / J_"WK"^* = 2^p r (1 + o(1)) -> 0$ linearly in $r$, proving an arbitrarily large expressiveness gap on axis-aligned multi-subspace clusters!

  (3) Rescaling coordinate $1$ by $s > 0$ multiplies $D_1$ by $s^gamma$, replacing $D_1^(-1/(gamma - 1))$ by $s^(-gamma / (gamma - 1)) D_1^(-1/(gamma - 1)) -> infinity$ as $s -> 0^+$; hence $w_1^* -> 1$ and $w_r^* -> 0$ ($r >= 2$). Likewise, rotating an axis-aligned signal subspace by $45 degree$ distributes signal dispersion equally across rotated coordinates, destroying the weight separation achieved before rotation.
]

#subsec("C", "Instance Constraints & Feasibility Phase Transitions (COP-KMeans)")

While all other 29 algorithms are strictly unsupervised, COP-KMeans (Alg. 20) @wagstaff2001constrained integrates domain-expert Must-Link ($E_"ML"$) and Cannot-Link ($E_"CL"$) graphs into Lloyd's assignment loop. Comparing COP-KMeans directly with unconstrained Lloyd (Alg. 1) and Hartigan--Wong (Alg. 19) exposes a sharp structural asymmetry between positive (`ML`) and negative (`CL`) constraints.

#proposition(title: "ML/CL Structural Asymmetry & Exact Penalty Bridge for COP-KMeans", origin: "Original Unified Proof")[
  Let $G_"ML" = (V, E_"ML")$ and $G_"CL" = (V, E_"CL")$ be Must-Link and Cannot-Link graphs on $cal(X) subset.eq bb(R)^d$.
  1. *Exact `ML` Reduction to Weighted Unconstrained K-Means:* Collapsing each connected component $S_m$ of $G_"ML"$ into a single super-point $overline(x)_(S_m) = |S_m|^(-1) sum_(i in S_m) x_i$ with integer weight $omega_m = |S_m|$ transforms any `ML`-only COP-KMeans problem into an *exact* weighted Euclidean K-Means problem (shifting WCSS by the partition-independent constant $sum_m sum_(i in S_m) norm(x_i - overline(x)_(S_m))_2^2$), preserving Voronoi convexity, strict Lyapunov descent, and polynomial feasibility.
  2. *`CL` Complexity Barrier & Greedy Failure:* With `CL` constraints and $k >= 3$, deciding whether *any* feasible $k$-partition exists is NP-complete (by direct reduction from Graph $k$-Coloring @wagstaff2001constrained). Moreover, even for $k = 2$ on a bipartite 4-cycle $C_4$ of `CL` edges, greedy point-by-point assignment in COP-KMeans can terminate in an infeasible dead-end whenever two opposite vertices of $C_4$ are first assigned to distinct clusters.
  3. *Exact Soft-Penalty Bridge:* Consider the soft-constrained objective $J_lambda (cal(C), bold(mu)) = J(cal(C), bold(mu)) + lambda V_"CL" (cal(C))$ penalizing each violated `CL` edge by $lambda > 0$. Whenever a feasible partition $cal(C)_"feas"$ exists, setting $lambda > J(cal(C)_"feas")$ guarantees that every global minimizer of $J_lambda$ has $V_"CL" = 0$ and coincides identically with the exact global optimum of hard COP-KMeans!
]

#proof[
  (1) By the parallel-axis identity (Lemma 1 of Alg. 1), for every `ML` component $S_m$ forced into a single cluster $C_j$,
  $ sum_(i in S_m) norm(x_i - mu_j)_2^2 = sum_(i in S_m) norm(x_i - overline(x)_(S_m))_2^2 + |S_m| norm(overline(x)_(S_m) - mu_j)_2^2. $
  The first sum is constant over all feasible partitions, so optimizing $(cal(C), bold(mu))$ over `ML`-feasible partitions is algebraically identical to weighted unconstrained K-Means on ${(overline(x)_(S_m), |S_m|)}$.

  (2) Given any graph $G = (V, E)$, setting $E_"CL" = E$ and $E_"ML" = emptyset$ makes a $k$-partition `CL`-feasible iff it is a proper $k$-coloring of $G$, proving NP-completeness for $k >= 3$. For greedy dead-end on $k = 2$, let $C_4$ have `CL` edges $(1, 2), (2, 3), (3, 4), (4, 1)$ and place $x_1$ closest to $mu_1$, $x_3$ closest to $mu_2$. Scanning $(1, 3, 2, 4)$ assigns $c_1 = 1$, $c_3 = 2$; then vertex $2$ has `CL` edges to *both* $1$ and $3$, leaving zero feasible clusters despite the valid 2-coloring $c_1 = c_3 = 1$, $c_2 = c_4 = 2$!

  (3) If $V_"CL" (cal(C)) >= 1$, then $J_lambda (cal(C), bold(mu)) >= 0 + lambda > J(cal(C)_"feas") = J_lambda (cal(C)_"feas")$, so no infeasible partition can minimize $J_lambda$. Among feasible partitions ($V_"CL" = 0$), $J_lambda$ equals $J$ identically.
]

#sec("V", "Complete 30-Algorithm Head-to-Head Matrix & Practitioner's Selection Guide")

We now synthesize all thirty algorithms of Category I (Algs. 1--30) in a unified comparative reference matrix (Table III), followed by an actionable engineering selection architecture linking dataset scale, geometry, contamination level, and structural constraints directly to the optimal algorithm combination.

#set page(columns: 1)

#block(
  width: 100%,
  above: 2pt,
  below: 4pt,
  {
    set align(center)
    text(size: 8.3pt, weight: "bold", fill: ieee_navy)[TABLE III: COMPLETE 30-ALGORITHM HEAD-TO-HEAD COMPARATIVE MATRIX (CATEGORY I: ALGS. 1--30)]
    v(3pt)
    set text(size: 6.45pt)
    table(
      columns: (1.42fr, 1.15fr, 1.15fr, 0.88fr, 1.28fr, 0.76fr, 1.02fr, 1.54fr),
      align: (left, left, left, left, left, center, left, left),
      inset: (x: 2.6pt, y: 1.95pt),
      [Alg. / Name], [Prototype & Metric], [Time Complexity], [Space], [Guarantee / Converg.], [Brk. ($epsilon_1^* slash epsilon_k^*$)], [Auto-$k$ / Key Params], [Strengths & Failure Modes],
      [*\#1* K-Means (Lloyd)], [Mean / $ell_2^2$], [$O(n k d T)$], [$O(n d)$], [Finite CVT; $2^(Omega(sqrt(n)))$ worst], [$1/n$ / $1/n$], [No ($k$)], [Fast baseline; unbounded ratio on bad seeds],
      [*\#2* K-Means++], [Mean / $ell_2^2$ ($D^2$ seed)], [$O(n k d (1 + T))$], [$O(n d)$], [$bb(E)[J] <= 8(ln k + 2) J^*$], [$1/n$ / $1/n$], [No ($k$)], [Gold-standard seeding; outlier starvation],
      [*\#3* Mini-Batch], [Mean / $ell_2^2$ (SGD)], [$O(b k d T)$], [$O(b d + k d)$], [Non-monotone; $O(1/(b t))$ MSE], [$1/n$ / $1/n$], [No ($k, b$)], [Web-scale $O(n/b)$ speed; variance floor],
      [*\#4* K-Medoids (PAM)], [Medoid $in X$ / gen. $d$], [$O(k (n - k)^2 T)$], [$O(n^2)$], [1-swap local opt.; 5-approx.], [$1/2$ / $1/n$], [No ($k$)], [Any dissimilarity; quadratic-times-$k$ slow],
      [*\#5* CLARA], [Medoid $in S_r$ / gen. $d$], [$O(R(k s^2 T + n k))$], [$O(s^2 + n k)$], [Subsample excess risk], [$1/2$ / $1/n$], [No ($k, R, s$)], [Linear in $n$; needs $s = Omega(k log k)$],
      [*\#6* CLARANS], [Medoid $in X$ / $J(n,k)$], [$O(r m_"max" n T)$], [$O(n d)$], [Stochastic Johnson walk], [$1/2$ / $1/n$], [No ($k, r, m_"max"$)], [No subsample lock; stalls if $m_"max"$ small],
      [*\#7* K-Modes], [Mode / Hamming], [$O(n k d T)$], [$O(n d)$], [Finite int. descent; $2 times$ euc.], [Bnd. / Bnd.], [No ($k$)], [Pure categorical; non-unique tie modes],
      [*\#8* K-Prototypes], [Mean+Mode / Mixed], [$O(n k (d_r + d_c) T)$], [$O(n d)$], [Finite CVT; $2 times$ one-hot], [$1/n$ / $1/n$], [No ($k, gamma$)], [Mixed tables; sensitive to $gamma$ scaling],
      [*\#9* Bisecting KM], [Divisive mean / $ell_2^2$], [$O(n d R T log k)$], [$O(n d)$], [Monotone splits; non-Voronoi], [$1/n$ / $1/n$], [Tree ($k, R$)], [Fast hierarchy; frozen ancestor cuts],
      [*\#10* X-Means], [Mean / $ell_2^2$ + BIC], [$O(n k_"max" d T)$], [$O(n d)$], [Consistent on isotropic Gauss.], [$1/n$ / $1/n$], [Yes ($k_"min"..k_"max"$)], [Auto-$k$ BIC; over-splits $kappa > kappa_"crit" (d)$],
      [*\#11* G-Means], [Mean / $ell_2^2$ + AD test], [$O(n k d T + n k log n)$], [$O(n d)$], [Anisotropic-Gaussian immune], [$1/n$ / $1/n$], [Yes ($alpha_"AD"$)], [Handles stretched Gauss.; shatters Uniform],
      [*\#12* Spherical KM], [Unit mean / Cosine], [$O("nnz"(X) k T)$], [$O("nnz"(X))$], [Monotone ascent; vMF limit], [Bnd. / Bnd.], [No ($k$)], [Sparse text/embeddings; ignores norms],
      [*\#13* Elkan's KM], [Mean / $ell_2^2$ (pruned)], [$O(n k d T)$ pruned], [$O(n k + n d)$], [Exact Lloyd trajectory equiv.], [$1/n$ / $1/n$], [No ($k$)], [Max pruning; $O(n k)$ RAM bottleneck],
      [*\#14* Hamerly's KM], [Mean / $ell_2^2$ (2-gate)], [$O(n k d T)$ pruned], [$O(n d)$], [Exact Lloyd trajectory equiv.], [$1/n$ / $1/n$], [No ($k$)], [Tiny $2n$ overhead; drift-spillover drop],
      [*\#15* K-Medians], [Coord. median / $ell_1$], [$O(n k d T)$], [$O(n d)$], [Finite $ell_1$ descent], [$1/2$ / $1/n$], [No ($k$)], [Robust local median; not $O(d)$-invariant],
      [*\#16* Kernel KM], [RKHS mean / $K(x,y)$], [$O(n^2 (d + T))$], [$O(n^2)$], [PSD monotone; Ncut dual], [Bnd.\* / Bnd.\*], [No ($k, K$)], [Non-convex manifolds; $O(n^2)$ Gram cost],
      [*\#17* Global KM], [Mean / $ell_2^2$ (greedy)], [$O(n^2 k^2 d T)$], [$O(n d)$], [Deterministic; $-> 4^-$ on $n=4$], [$1/n$ / $1/n$], [Nested $1..k$], [Seed-free deterministic; $O(n^2 k^2)$ slow],
      [*\#18* ISODATA], [Mean / Split--Merge], [$O((n k + k^2) d T)$], [$O(n d)$], [Can cycle w/o hysteresis], [$1/n$ / $1/n$], [Yes ($theta_S, theta_C, theta_N$)], [Domain-tunable; 6+ brittle thresholds],
      [*\#19* Hartigan--Wong], [Mean / 1-pt transfer], [$O(n k d T)$], [$O(n d)$], [$cal(E)_"HW" subset.neq cal(E)_"Lloyd"$], [$1/n$ / $1/n$], [No ($k$)], [Escapes Lloyd traps; $O(n_"min"^(-1))$ gain],
      [*\#20* COP-KMeans], [Mean / `ML`+`CL`], [$O(T k(n d + |E|))$], [$O(n d + |E|)$], [NP-complete `CL` ($k >= 3$)], [$1/n$ / $1/n$], [No ($k, "ML", "CL"$)], [Semi-supervised; greedy dead-ends],
      [*\#21* MacQueen's KM], [Online mean / $ell_2^2$], [$O(n k d T)$], [$O(k d)$], [Robbins--Monro a.s. conv.], [$1/n$ / $1/n$], [No ($k$)], [Single-pass streaming; order-dependent],
      [*\#22* Forgy's KM], [Mean / Data seeds], [$O(n k d T)$], [$O(n d)$], [Finite CVT; $O(n/k)$ spread], [$1/n$ / $1/n$], [No ($k$)], [Beats rand. part. spread; misses small clusters],
      [*\#23* K-Means||], [Parallel $ell D^2$ seed], [$O(n d r ell + r ell k d)$], [$O(n d)$], [$O(log k)$ in $O(log psi)$ rds.], [$1/n$ / $1/n$], [No ($k, r, ell$)], [MapReduce/GPU parallel; fails fixed $r$],
      [*\#24* Fast/FasterPAM], [Medoid $in X$ / gen. $d$], [$O((n - k)^2 T)$], [$O(n^2)$], [Exact PAM / 5-approx.], [$1/2$ / $1/n$], [No ($k$)], [$Theta(k)$ faster than PAM; strictly dominates \#4],
      [*\#25* CLARA-MED], [Medoid / Median score], [$O(q(s^2 T + n k))$], [$O(s^2 + n)$], [Moment-free $O(sqrt((log q)/n))$], [$1/2$ / $1/2$\*], [No ($k, q, s$)], [50% scoring brk.; masks $< 50\%$ clusters],
      [*\#26* Trimmed KM], [Mean / Impartial $alpha$], [$O(n(k d + log n) T)$], [$O(n d)$], [Moment-free strong cons.], [$alpha$ / $<= alpha$], [No ($k, alpha$)], [True global brk.; over-trim erases small clusters],
      [*\#27* Power KM], [MM $M_s$ homotopy], [$O(n k d T)$], [$O(n d + n k)$], [$k^(-1/s)$ unif.; escapes traps], [$1/n$ / $1/n$], [No ($k, s_0, eta$)], [Dissolves Voronoi walls; extra $s$-annealing],
      [*\#28* W-K-Means], [Global wtd. $ell_2^2$], [$O(n k d T)$], [$O(n d)$], [3-block finite; noise suppr.], [$1/n$ / $1/n$], [No ($k, beta > 1$)], [Global feature selection; scale collapse],
      [*\#29* MWK-Means], [Local wtd. $ell_p^p$], [$O(n d(k + T_"in") T)$], [$O(n d + k d)$], [Cluster-local Hölder opt.], [$1/n$ / $1/n$], [No ($k, p > 1, epsilon_0$)], [Local subspaces; non-convex cells, scale-sens.],
      [*\#30* LBG (VQ)], [Codebook split / $rho$], [$O(n k d T_"max")$], [$O(n d)$], [Lloyd--Max CVT; $2^B$ VQ], [$1/n$ / $1/n$], [No ($k = 2^B, epsilon, tau$)], [Deterministic VQ/Bregman; greedy split traps],
    )
    v(2pt)
    set align(left)
    set par(first-line-indent: 0pt, justify: true)
    text(size: 6.9pt, fill: ieee_slate)[*Notes:* Bnd. $=$ bounded prototype domain (unit sphere, finite alphabet, or bounded RKHS kernel); Brk. ($epsilon_1^* slash epsilon_k^*$) lists local ($k = 1$) vs. global ($k >= 2$) finite-sample breakdown ceilings (Theorem 5); CLARA-MED's $1/2$\* refers to its outer median-scoring breakdown; Kernel KM's Bnd.\* assumes bounded kernel $K(x, x) <= M$ vanishing at infinity.]
  },
)

#set page(columns: 2)

#subsec("A", "Strict Dominance Relations & Deprecated Baselines")

Our comparative theorems identify three *strict dominance* pairs within Category I where one algorithm matches or outperforms another mathematically across essentially all practical regimes:
1. *FasterPAM (Alg. 24) strictly dominates classical PAM (Alg. 4):* FastPAM1 is provably bit-identical to PAM while saving an exact $Theta(k)$ factor per swap pass, and FasterPAM preserves the exact metric 5-approximation certificate while slashing total iterations by an additional order of magnitude @schubert2021fast. Classical PAM should be retained only for pedagogical reference.
2. *K-Means++ (Alg. 2) / K-Means|| (Alg. 23) dominate bare Forgy (Alg. 22) and Random Partition (Alg. 1) on clean data:* Whenever gross outlier contamination is absent or pre-filtered, $D^2$-seeding replaces exponential cluster-miss failure with an unconditional $8(ln k + 2)$ guarantee at negligible $O(n k d)$ seeding overhead (or $O(log psi)$ parallel passes via K-Means||).
3. *Elkan (Alg. 13) and Hamerly (Alg. 14) strictly dominate naive Lloyd (Alg. 1) for exact batch refinement:* Both return bit-identical partitions to Lloyd while skipping $60\%$--$99\%$ of distance evaluations on clustered data in low-to-moderate dimensions ($d <= 100$); bare Lloyd only remains competitive on GPU matrix-multiplication kernels (GEMM) where dense batched inner products outperform branchy bound checks.

#subsec("B", "Practitioner's Engineering Selection Guide")

To operationalize Table III, practitioners should navigate Category I along seven quantitative decision regimes:

1. *Standard Clean Euclidean Data ($n <= 10^6$, fixed $k$):* Initialize with *K-Means++* (Alg. 2; multi-start $R = 10$--$20$). Refine using *Hamerly* (Alg. 14) when $d <= 50$ or $n k > 10^7$ (low memory footprint), *Elkan* (Alg. 13) when $50 < d <= 200$ and $n k$ fits in L3 cache/RAM, or *Hartigan--Wong* (Alg. 19) when cluster sizes are small or imbalanced ($n_"min" <= 100$, where its $O(n_"min"^(-1))$ escape shell yields measurable WCSS improvements). If severe non-convex local traps persist even across multi-starts, replace hard Lloyd refinement with *Power K-Means* (Alg. 27, geometric annealing $s_0 = -1$, $eta in [1.05, 1.2]$). Use *Global K-Means* (Alg. 17) only on small benchmark datasets ($n <= 2000$) requiring strict seed-free determinism.
2. *Web-Scale, Streaming, and Distributed Workloads ($n > 10^6$):* Seed with *K-Means||* (Alg. 23, $ell = 2 k$, $r = 5$--$8$ rounds—never $r = 1$ on multiscale data) and refine via *Mini-Batch K-Means* (Alg. 3, batch size $b = max(1024, 10 k)$ to keep per-center MSE $O(sigma_j^2 / (b p_j t))$ small). Reserve *MacQueen* (Alg. 21) exclusively for strict single-pass $O(k d)$-memory online streams where revisiting historical observations is physically impossible.
3. *Arbitrary Dissimilarities, Non-Euclidean Metrics & Exemplar Prototypes:* Run *FasterPAM* (Alg. 24) whenever $n <= 2 times 10^4$ (full $O(n^2)$ distance matrix fits in memory). For $n > 2 times 10^4$, deploy *CLARA* paired with FasterPAM (Alg. 5 + Alg. 24) using principled subsample size $s = max(200, ceil(5 k ln k))$ rather than classical $40 + 2 k$ (Theorem 6(4)), or *CLARANS* (Alg. 6) when subsampling risks missing rare spatial clusters. Under heavy metric contamination with balanced clusters (no critical minority cluster below $50\%$), switch outer scoring to *CLARA-MED* (Alg. 25).
4. *Heavy-Tailed Noise & Adversarial Outliers in $bb(R)^d$:* When an upper bound $alpha$ on the outlier fraction is known and satisfies $alpha < min_j (|C_j^*|) / n$, *Trimmed K-Means* (Alg. 26) is the sole Euclidean centroid method offering positive global breakdown ($epsilon_k^* > 0$) and moment-free strong consistency @cuesta1997trimmed; diagnose over-trimming via trimming-curve elbows across $alpha in [0.01, 0.20]$. When coordinate-wise $ell_1$ geometry is native and outliers are moderate (or $k$ can be slightly over-provisioned to absorb hijacked slots), use *K-Medians* (Alg. 15).
5. *Directional Text/Embeddings, Nonlinear Manifolds, and Categorical/Mixed Tables:* For TF--IDF vectors or normalized neural embeddings on $bb(S)^(d-1)$, use *Spherical K-Means* (Alg. 12) exploiting sparse $O("nnz"(X) k T)$ dot products. For concentric or nonlinearly separable clusters with $n <= 10^4$, apply *Kernel K-Means* (Alg. 16) with RBF or graph-cut kernels @dhillon2004kernel (scaling beyond $10^4$ via Nyström approximation). For pure categorical schemas use *K-Modes* (Alg. 7), and for mixed numeric--categorical tables use *K-Prototypes* (Alg. 8) with variance-calibrated $gamma$ (matching mean numeric coordinate variance to prevent numeric or categorical dominance).
6. *Unknown $k$, Divisive Hierarchies, and Signal Compression:* If clusters are approximately Gaussian with unknown aspect ratios, prefer *G-Means* (Alg. 11) over *X-Means* (Alg. 10) to avoid Theorem 7(1)'s $kappa > kappa_"crit" (d)$ over-splitting collapse; conversely, if clusters are near-spherical or bounded/uniform (platykurtic), prefer *X-Means* (Alg. 10) over G-Means to avoid uniform shattering. Use *Bisecting K-Means* (Alg. 9) for ultra-fast divisive document taxonomies (optionally followed by one global Lloyd polish to restore Voronoi optimality), *LBG* (Alg. 30) for $2^B$-level vector quantization codebooks and Bregman speech/signal distortions, and *ISODATA* (Alg. 18) only in domain-calibrated remote-sensing pipelines with explicit hysteresis deadbands.
7. *High-Dimensional Noise Subspaces & Semi-Supervised Constraints:* Always standardize features ($z$-score or range normalization) before weighting to prevent Theorem 8(3)'s inverse scale collapse; then apply *W-K-Means* (Alg. 28, $beta in [1.5, 4]$) when all clusters share one relevant global subspace, or *MWK-Means* (Alg. 29, $p in [1.5, 3]$ with floor $epsilon_0 > 0$) when different clusters live in distinct local feature subspaces or exhibit non-Gaussian $ell_p$ tails. When domain knowledge supplies pairwise `ML`/`CL` supervision, apply *COP-KMeans* (Alg. 20) after collapsing `ML` connected components into weighted super-nodes, and fall back to soft-penalty optimization ($lambda > J_"feas"$) whenever `CL` graph density induces greedy infeasibility.

#subsec("C", "Unresolved Boundary with External Paradigms")

Finally, our intra-category analysis exposes the intrinsic ceiling shared by *all thirty* algorithms in Category I (with the partial exception of Kernel K-Means Alg. 16 and non-convex local-metric MWK-Means Alg. 29): prototype minimization partitions space into convex (or generalized Voronoi) cells anchored by $k$ discrete representatives. Whenever clusters form elongated density filaments, arbitrary non-convex manifolds with background clutter, or nested multiscale hierarchies, prototype partitioning must yield to the density-connectivity (Part II), hierarchical/ultrametric (Part III), probabilistic mixture (Part IV), and spectral-graph (Part V) paradigms analyzed in Comparison Chapters IV--VII.
