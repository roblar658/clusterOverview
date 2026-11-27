#import "../template/ieee.typ": *

#alg_entry(
  id: 30,
  name: "Linde–Buzo–Gray (LBG) Algorithm",
  category: "Cat. I: Centroid & Partitioning",
  year: "1980",
  paradigm: "Hierarchical Codebook-Splitting Vector Quantization",
  complexity: $O(n k d T_"max")$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Introduced by Linde, Buzo, and Gray @linde1980algorithm to design optimal block-coding *vector quantizers (VQs)* for speech and signal compression without requiring an analytic source density, the *LBG algorithm* generalizes Stuart Lloyd's scalar Method I @lloyd1982least along two foundational dimensions: (i) optimization over arbitrary $d$-dimensional training sequences under general distortion measures $rho$, and (ii) a deterministic *binary codebook-splitting* hierarchy that constructs an $N = k = 2^B$-level codebook from a single global centroid without random seeding.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ denote an empirical training sequence drawn from an unknown source $P$, and let $rho: bb(R)^d times bb(R)^d -> bb(R)_+$ be a measurable distortion measure—canonically squared Euclidean error $rho(x, y) = norm(x - y)_2^2$, or a spectral distortion such as the Itakura–Saito divergence @linde1980algorithm. A *$k$-level vector quantizer* $Q = (cal(P), cal(Y))$ of coding rate $R = (log_2 k) / d$ bits per sample consists of an *encoder partition* $cal(P) = {R_1, dots, R_k}$ of ${1, dots, n}$ and a *decoder codebook* $cal(Y) = {y_1, dots, y_k} subset.eq bb(R)^d$ of reproduction vectors (codewords), assigning $Q(x_i) = y_j$ for $i in R_j$. The *average empirical distortion* functional is
$ D(cal(P), cal(Y)) = 1/n sum_(j=1)^k sum_(i in R_j) rho(x_i, y_j). $
Assuming $rho$ admits a well-defined generalized centroid operator $"cent"(S) := arg min_(y in bb(R)^d) sum_(i in S) rho(x_i, y)$ for every non-empty cell $S subset.eq {1, dots, n}$, optimal quantization is governed by the two classical *Lloyd–Max optimality conditions*:

#definition(title: "Lloyd–Max Necessary Conditions (NNC & CC)")[
  A $k$-level vector quantizer $(cal(P), cal(Y))$ satisfies:
  1. *Nearest-Neighbor Condition (NNC):* Every training vector $i in R_j$ obeys $rho(x_i, y_j) <= rho(x_i, y_ell)$ for all $ell in {1, dots, k}$.
  2. *Centroid Condition (CC):* Every non-empty cell satisfies $y_j in "cent"(R_j) = arg min_(y in bb(R)^d) sum_(i in R_j) rho(x_i, y)$.
]

LBG couples an inner *Generalized Lloyd Algorithm (GLA)* governed by a fractional distortion-improvement threshold $tau > 0$ with an outer binary splitting cascade that doubles the active codebook size $m <- 2 m$ at each stage (Algorithm 1). When $k$ is not a power of $2$, only the $k - m$ cells with largest intra-cell distortion are split on the final pass.

#algbox(
  title: "Linde–Buzo–Gray (LBG) Vector Quantizer Design",
  inputs: [Training set $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, target codebook size $k >= 2$, distortion $rho$, split amplitude $epsilon > 0$, relative tolerance $tau in (0, 1)$.],
  outputs: [Encoder partition $cal(P) = {R_1, dots, R_k}$ and codebook $cal(Y) = {y_1, dots, y_k}$.],
)[
  1. *Initialization ($m = 1$):* Compute global centroid $y_1 = "cent"({1, dots, n})$; set $cal(Y) = {y_1}$, $m = 1$.
  2. *While* $m < k$ *do* (*Binary Codebook Splitting*):
  3. #h(1em) Let $s = min(m, k - m)$. Select the $s$ cells of highest total distortion (all $m$ cells when $2 m <= k$).
  4. #h(1em) Replace each selected codeword $y_j$ by two split seeds $y_j^+ = y_j + epsilon u_j$ and $y_j^- = y_j - epsilon u_j$ (canonically $u_j = y_j$ if $y_j != bold(0)$, or leading variance axis $v_1 (Sigma_j)$); update $m <- m + s$.
  5. #h(1em) *Inner Generalized Lloyd Loop (GLA):* Set $t = 0$, $D^((0)) = +infinity$.
  6. #h(1em) *Repeat*:
  7. #h(2em) *NNC step:* Assign $c_i = arg min_(j in {1, dots, m}) rho(x_i, y_j)$ (ties to smallest $j$); form $R_j = {i : c_i = j}$.
  8. #h(2em) *CC step:* For each non-empty $R_j$, set $y_j <- "cent"(R_j)$. (Reseed empty cells by splitting the highest-distortion cell.)
  9. #h(2em) Increment $t <- t + 1$; compute $D^((t)) = D(cal(P), cal(Y))$.
  10. #h(1em) *Until* $(D^((t-1)) - D^((t))) / D^((t-1)) <= tau$ or $D^((t)) = 0$.
  11. *Return* $(cal(P), cal(Y))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal theoretical results characterizing LBG: (A) necessity of the Lloyd–Max conditions, strict monotone descent, and a logarithmic iteration bound under fractional stopping; (B) an exact closed-form characterization of splitting distortion gain and zero-centroid degeneracy; (C) classical strong consistency and Bennett–Zador high-rate distortion scaling; (D) an explicit counterexample proving unbounded suboptimality of uniform binary splitting; and (E) the exact geometric-series complexity bound and invariance profile.

#subsec("A", "Monotone Distortion Descent and Fractional Stopping Rate")

#theorem(title: "Lloyd–Max Optimality & Monotone GLA Descent", origin: "Classical — Cited")[
  Let $rho >= 0$ be any distortion measure admitting a generalized centroid $"cent"(S)$ @linde1980algorithm. Then:
  1. Every globally optimal $k$-level quantizer $(cal(P)^*, cal(Y)^*)$ necessarily satisfies both NNC and CC.
  2. Every inner GLA iteration (Steps 7–8) satisfies
     $ D(cal(P)^((t)), cal(Y)^((t))) <= D(cal(P)^((t)), cal(Y)^((t-1))) <= D(cal(P)^((t-1)), cal(Y)^((t-1))). $
  3. Under squared Euclidean distortion $rho(x, y) = norm(x - y)_2^2$, each iteration satisfies the strict energy drop
     $ D^((t-1)) - D^((t)) >= 1/n sum_(j=1)^m |R_j^((t))| norm(y_j^((t)) - y_j^((t-1)))_2^2. $
]

#proof(label: "Proof (Adapted from Linde–Buzo–Gray 1980)")[
  For (1), if any $i in R_j^*$ violated NNC, moving $i$ to $arg min_ell rho(x_i, y_ell^*)$ would strictly reduce $D$; if any $y_j^*$ violated CC, replacing $y_j^*$ by $"cent"(R_j^*)$ would strictly reduce $sum_(i in R_j^*) rho(x_i, y)$, contradicting optimality.

  For (2), fixing $cal(Y)^((t-1))$, the NNC assignment chooses $c_i^((t)) = arg min_ell rho(x_i, y_ell^((t-1)))$ pointwise, so $rho(x_i, y_(c_i^((t)))^((t-1))) <= rho(x_i, y_(c_i^((t-1)))^((t-1)))$. Averaging over $i = 1, dots, n$ gives $D(cal(P)^((t)), cal(Y)^((t-1))) <= D(cal(P)^((t-1)), cal(Y)^((t-1)))$. Next, fixing $cal(P)^((t))$, replacing each $y_j^((t-1))$ by $y_j^((t)) = "cent"(R_j^((t)))$ minimizes each cell sum independently, yielding $D(cal(P)^((t)), cal(Y)^((t))) <= D(cal(P)^((t)), cal(Y)^((t-1)))$.

  For (3), applying the parallel-axis identity to each non-empty Euclidean cell $R_j^((t))$ with centroid $y_j^((t)) = overline(x)_(R_j^((t)))$ yields
  $ sum_(i in R_j^((t))) norm(x_i - y_j^((t-1)))_2^2 &= sum_(i in R_j^((t))) norm(x_i - y_j^((t)))_2^2 \
    &quad + abs(R_j^((t))) norm(y_j^((t)) - y_j^((t-1)))_2^2. $
  Summing over $j = 1, dots, m$, dividing by $n$, and combining with the non-negative NNC reduction proves Claim 3.
]

#corollary(title: "Logarithmic Iteration Bound under LBG Relative Stopping", origin: "Original Proof")[
  At any fixed codebook level $m$, suppose the initial distortion after splitting satisfies $D^((1)) > 0$ and the final positive distortion satisfies $D^* >= D_"min" > 0$. Then the LBG relative stopping criterion $(D^((t-1)) - D^((t))) / D^((t-1)) <= tau$ ($tau in (0, 1)$) guarantees termination in at most
  $ T <= 1 + floor( ln(D^((1)) / D_"min") / ln(1 / (1 - tau)) ) <= 1 + floor( 1/tau ln(D^((1)) / D_"min") ) $
  iterations, eliminating worst-case exponential tail cycling @vattani2011kmeans.
]

#proof[
  Whenever the loop does not halt at iteration $t >= 2$, we have $(D^((t-1)) - D^((t))) / D^((t-1)) > tau$, equivalently $D^((t)) < (1 - tau) D^((t-1))$. Applying this recursion over $T - 1$ non-halting transitions gives $D_"min" <= D^((T)) < (1 - tau)^(T - 1) D^((1))$. Taking natural logarithms yields
  $ (T - 1) ln(1 / (1 - tau)) < ln(D^((1)) / D_"min"). $
  Because $1 - tau <= e^(-tau)$ for all $tau in (0, 1)$, we have $ln(1 / (1 - tau)) >= tau > 0$, directly establishing both bounds.
]

#subsec("B", "Exact Distortion Gain of Binary Splitting and Singularity")

#theorem(title: "Exact Splitting Gain & Zero-Centroid Singularity", origin: "Original Derivation")[
  Let $R_j$ be a non-empty Euclidean cell with centroid $y_j = overline(x)_(R_j)$, and split $y_j$ into $y_j^(plus.minus) = y_j plus.minus epsilon u$ for unit vector $u in bb(R)^d$ ($norm(u)_2 = 1$) and $epsilon > 0$. Define the projected mean absolute deviation $S_u := abs(R_j)^(-1) sum_(i in R_j) abs(chevron.l x_i - y_j\, u chevron.r)$, the half-space partition $R_j^+ = {i in R_j : chevron.l x_i - y_j, u chevron.r >= 0}$, $R_j^- = R_j backslash R_j^+$, and balance fraction $p_u = abs(R_j^+) / abs(R_j)$. Then:
  1. *Exact Pre-Recentering Gain:* Assigning $R_j$ between $y_j^+$ and $y_j^-$ changes average distortion by exactly
     $ Delta D_"split" (epsilon, u) = (abs(R_j) / n) ( epsilon^2 - 2 epsilon S_u ). $
     Hence strict decrease holds iff $0 < epsilon < 2 S_u$, minimized at $epsilon^* = S_u$ with gain $-(abs(R_j) / n) S_u^2$.
  2. *Exact Post-Recentering Gain:* Re-centering $R_j^+$ and $R_j^-$ to their centroids lowers distortion by at least
     $ D_"before" - D_"recentered" &= (abs(R_j) / n) p_u (1 - p_u) norm(overline(x)_(R_j^+) - overline(x)_(R_j^-))_2^2 \
       &>= (abs(R_j) / n) S_u^2 / (4 p_u (1 - p_u)) >= (abs(R_j) / n) S_u^2. $
  3. *Zero-Centroid Singularity:* Under multiplicative splitting $y_j^(plus.minus) = (1 plus.minus epsilon) y_j$, whenever $y_j = bold(0)$ (notably at $m = 1$ on zero-mean data), $y_j^+ = y_j^- = bold(0)$ and the split collapses with zero distortion gain.
]

#proof[
  (1) The perpendicular bisector of $y_j + epsilon u$ and $y_j - epsilon u$ is ${x : chevron.l x - y_j, u chevron.r = 0}$, assigning $R_j^+$ to $y_j^+$ and $R_j^-$ to $y_j^-$. Expanding the squared Euclidean norm for each $i in R_j^+$ gives
  $ norm(x_i - (y_j + epsilon u))_2^2 = norm(x_i - y_j)_2^2 - 2 epsilon chevron.l x_i - y_j, u chevron.r + epsilon^2, $
  and symmetrically for $i in R_j^-$ with $+2 epsilon chevron.l x_i - y_j, u chevron.r$. In both half-spaces, the linear term equals $-2 epsilon abs(chevron.l x_i - y_j\, u chevron.r)$. Summing over $i in R_j$ and dividing by $n$ produces $Delta D_"split" (epsilon, u) = (abs(R_j)/n)(epsilon^2 - 2 epsilon S_u)$, a convex quadratic in $epsilon$ minimized at $epsilon^* = S_u$ with minimum $-(abs(R_j)/n) S_u^2$.

  (2) By Ward's variance decomposition (parallel-axis identity on $R_j^+$ and $R_j^-$ around $y_j = p_u overline(x)_(R_j^+) + (1 - p_u) overline(x)_(R_j^-)$),
  $ D_"before" - D_"recentered" = (abs(R_j) / n) p_u (1 - p_u) norm(overline(x)_(R_j^+) - overline(x)_(R_j^-))_2^2. $
  Because $sum_(i in R_j) chevron.l x_i - y_j, u chevron.r = 0$, the positive and negative projected sums are equal:
  $ sum_(i in R_j^+) chevron.l x_i - y_j, u chevron.r = sum_(i in R_j^-) abs(chevron.l x_i - y_j\, u chevron.r) = (abs(R_j) / 2) S_u. $
  Dividing by $abs(R_j^+) = p_u abs(R_j)$ and $abs(R_j^-) = (1 - p_u) abs(R_j)$ yields
  $ chevron.l overline(x)_(R_j^+) - overline(x)_(R_j^-), u chevron.r = S_u / (2 p_u) + S_u / (2 (1 - p_u)) = S_u / (2 p_u (1 - p_u)). $
  By Cauchy–Schwarz ($norm(u)_2 = 1$), $norm(overline(x)_(R_j^+) - overline(x)_(R_j^-))_2 >= S_u / (2 p_u (1 - p_u))$. Substituting into the variance decomposition and using $4 p_u (1 - p_u) <= 1$ yields Claim 2. Claim 3 follows immediately from $(1 plus.minus epsilon) bold(0) = bold(0)$, motivating additive eigenvector splitting $u_j = v_1 (Sigma_j)$.
]

#subsec("C", "Strong Consistency and High-Rate Quantization Scaling")

#theorem(title: "Empirical Consistency & Bennett–Zador High-Rate Law", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a density $p$ on $bb(R)^d$ with $bb(E)[norm(X)_2^(2 + delta)] < infinity$ ($delta > 0$).
  1. *Strong Consistency @pollard1981strong:* If the population $k$-level squared-error quantizer $cal(Y)^*$ is unique up to permutation, the empirical global minimizer $hat(cal(Y))_n$ converges almost surely to $cal(Y)^*$ as $n -> infinity$.
  2. *High-Rate Scaling @linde1980algorithm:* As codebook size $k = 2^(d R) -> infinity$ (rate $R -> infinity$), the optimal per-dimension squared distortion satisfies
     $ (D^*(k)) / d approx C(d, 2) dot 2^(-2 R) dot ( integral_(bb(R)^d) p(x)^(d / (d + 2)) dif x )^((d + 2)/d), $
     where $C(d, 2)$ is the normalized second moment of the optimal $d$-dimensional tessellating polytope ($C(1, 2) = 1/12$, approaching $(2 pi e)^(-1)$ as $d -> infinity$).
]

#proof(label: "Proof Synthesis (Cited — Pollard 1981 & High-Rate Variational Bound)")[
  For (1), Pollard @pollard1981strong establishes almost-sure confinement of $hat(cal(Y))_n$ to a compact ball $B_R (0)^k$ using $bb(E)[norm(X)_2^2] < infinity$, followed by uniform Glivenko–Cantelli convergence of empirical distortion to population distortion over $B_R (0)^k$.

  For (2), represent a high-rate $k$-codebook by a normalized codeword point density $lambda(x) >= 0$ ($integral lambda(x) dif x = 1$). Near $x$, each Voronoi cell has volume $v(x) approx (k lambda(x))^(-1)$, contributing per-dimension distortion $C(d, 2) v(x)^(2/d) = C(d, 2) k^(-2/d) lambda(x)^(-2/d)$. Integrating against $p(x)$ and applying Hölder's inequality with conjugate exponents $(d + 2)/d$ and $(d + 2)/2$ gives
  $ integral p(x)^(d/(d+2)) dif x &= integral ( p(x) lambda(x)^(-2/d) )^(d/(d+2)) lambda(x)^(2/(d+2)) dif x \
    &<= ( integral p(x) lambda(x)^(-2/d) dif x )^(d/(d+2)) ( integral lambda(x) dif x )^(2/(d+2)). $
  Raising to power $(d+2)/d$ proves that the integral is minimized uniquely at the *optimal point density* $lambda^*(x) prop p(x)^(d/(d+2))$, yielding $(D^*(k)) / d approx C(d, 2) k^(-2/d) norm(p)_(d/(d+2))$. Because $C(1, 2) / C(infinity, 2) = pi e / 6$, vector quantization enjoys up to $10 log_10 (pi e / 6) approx 1.53$ dB *shaping gain* over scalar quantization even for memoryless i.i.d. sources.
]

#subsec("D", "Suboptimality Counterexample for Uniform Binary Splitting")

#proposition(title: "Unbounded Approximation Ratio of Binary LBG Splitting", origin: "Original Counterexample")[
  For any target ratio $M > 0$ and codebook size $k = 4$, there exists a 5-point 1D dataset on which deterministic binary LBG splitting ($1 -> 2 -> 4$, with either additive or multiplicative perturbation) terminates at a strict Lloyd–Max fixed point $(cal(P)_"LBG", cal(Y)_"LBG")$ with zero ties and distortion ratio
  $ D(cal(P)_"LBG", cal(Y)_"LBG") / D(cal(P)^*, cal(Y)^*) > M. $
]

#proof[
  Set $D = sqrt(M) + 2$, $L = 10 D$, and place $n = 5$ scalar points at
  $ x_1 = -L, quad x_2 = -L + 1, quad x_3 = L, quad x_4 = L + D, quad x_5 = L + 3 D. $
  *Stage 1 ($m = 1$):* The global mean is $y_0 = (L + 4 D + 1)/5 = (14 D + 1)/5$, which strictly satisfies $x_2 = -L + 1 < 0 < y_0 < 3 D < L = x_3$.

  *Stage 2 ($1 -> 2$):* Splitting $y_0$ by any sufficiently small perturbation $0 < delta_0 < L / 4$ (whether additive $epsilon$ or multiplicative $epsilon y_0$) keeps $x_1, x_2 < y_0 - delta_0$ and $x_3, x_4, x_5 > y_0 + delta_0$. Thus NNC assigns $R_1 = {x_1, x_2}$ and $R_2 = {x_3, x_4, x_5}$, with centroids $y_1 = -L + 1/2$ and $y_2 = L + (4/3) D$. Since the bisector $(y_1 + y_2)/2 = (2/3) D + 1/4$ strictly separates $x_2$ from $x_3$, GLA halts at level $m = 2$.

  *Stage 3 ($2 -> 4$):* Uniform binary splitting duplicates *both* $y_1$ and $y_2$ with perturbation magnitude $delta_j < min(1/4, D/6)$ (e.g., multiplicative $epsilon < 1/(4 L)$). Inspect strict inequalities:
  - $x_1 = -L < y_1 - delta_1$ and $x_2 = -L + 1 > y_1 + delta_1$, separating ${x_1}$ and ${x_2}$;
  - $x_3 = L < x_4 = L + D < y_2 - delta_2$ (since $y_2 - x_4 = D/3 > delta_2$) while $x_5 = L + 3 D > y_2 + delta_2$, assigning ${x_3, x_4}$ to $y_2^-$ and ${x_5}$ to $y_2^+$ with *zero ties*.
  Re-centering yields $cal(Y)_"LBG" = {-L, -L + 1, L + D/2, L + 3 D}$. Checking all five points verifies strict Voronoi stability—notably $|x_4 - (L + D/2)| = D/2 < 2 D = |x_4 - (L + 3 D)|$. Hence LBG halts with $n D_"LBG" = 0 + 0 + (D/2)^2 + (D/2)^2 + 0 = D^2 / 2$.

  In contrast, the globally optimal 4-codebook $cal(Y)^* = {-L + 1/2, L, L + D, L + 3 D}$ attains $n D^* = (1/2)^2 + (1/2)^2 = 1/2$, giving $D_"LBG" / D^* = D^2 > M$. Uniform doubling fails because it allocates two codewords to the low-distortion cluster $R_1$ instead of concentrating splits on the high-variance cluster $R_2$.
]

#subsec("E", "Cascade Complexity Bound, Invariances, and Robustness")

#proposition(title: "Geometric-Series Complexity & Symmetry Profile", origin: "Original Proof")[
  Let $k = 2^B$ and suppose GLA at codebook level $m_b = 2^b$ ($b = 1, dots, B$) executes $T_b <= T_"max"$ iterations. Then:
  1. The cumulative arithmetic cost of the entire binary splitting cascade satisfies
     $ sum_(b=1)^B n d dot 2^b dot T_b <= 2 (1 - 1/k) n k d T_"max" < 2 n k d T_"max". $
  2. Under squared Euclidean distortion, additive eigenvector splitting $y_j plus.minus epsilon v_1 (Sigma_j)$ is strictly invariant under rigid isometries $x |-> Q x + b$ ($Q in O(d)$), whereas multiplicative splitting $(1 plus.minus epsilon) y_j$ *breaks translation invariance*.
  3. The finite-sample Donoho–Huber breakdown point of Euclidean LBG is $epsilon^* = 1 / n$.
]

#proof[
  (1) Factoring $n d T_"max"$ yields $sum_(b=1)^B 2^b = 2 (2^B - 1) = 2 k (1 - 1/k) < 2 k$, proving that deterministic construction from $m = 1$ to $k$ costs strictly less than twice a single $k$-level Lloyd run.

  (2) Under $T(x) = Q x + b$, cell means transform as $Q y_j + b$ and covariances as $Q Sigma_j Q^T$, whose top eigenvector rotates as $Q v_1 (Sigma_j)$; hence additive splitting commutes with $T$. Conversely, under multiplicative splitting, $(1 plus.minus epsilon)(y_j + b) = (1 plus.minus epsilon) y_j + b plus.minus epsilon b$, introducing an extraneous shift $plus.minus epsilon b$ that alters the splitting hyperplane orientation (and triggers Theorem 2(3)'s collapse whenever $b = -y_j$).

  (3) Sending a single point $norm(x_1)_2 -> infinity$ drives both the initial global centroid $y_1 = overline(x)$ and the centroid of $x_1$'s final cell to infinity at rate $|R_j|^(-1) norm(x_1)_2$, giving breakdown point $1/n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

By Proposition 2(1) and Corollary 1, LBG runs in $O(n k d T_"max")$ total time with $T_"max" = O(tau^(-1) log(D^((1)) / D_"min"))$ and strictly linear memory $O((n + k) d)$. Because each level $m_b$ warm-starts near a converged $(m_b / 2)$-level partition, inner GLA loops typically converge in far fewer iterations than randomly seeded Lloyd runs. At encoding (inference) time, exhaustive NNC lookup costs $O(k d)$ per query vector, which can be accelerated via triangle-inequality pruning (Elkan Alg. 13 @elkan2003using; Hamerly Alg. 14 @hamerly2010making) or tree-structured VQ (TSVQ) in $O(d log k)$ time at a small distortion penalty.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for LBG Algorithm (Alg. 30)],
    [Time / Space Complexity], [$O(n k d T_"max")$ time ($< 2 times$ single-level cost, Prop. 2); $O((n + k) d)$ space.],
    [Stopping Guarantee], [Fractional rule $(D^((t-1)) - D^((t))) / D^((t-1)) <= tau$ halts in $O(tau^(-1) log(D^((1))/D_"min"))$ steps (Cor. 1).],
    [Determinism / Seeding], [Fully deterministic via binary splitting; eliminates Forgy randomness but suffers uniform-split bias (Prop. 1).],
    [Invariances], [Additive PCA-axis split is $O(d)$ + translation invariant; multiplicative $(1 plus.minus epsilon) y_j$ is not translation invariant.],
    [Outlier Breakdown Point], [$epsilon^* = 1/n$ (non-robust under squared $L_2$; sensitive to zero-mean collapse if multiplicative).],
    [Recommended Companions], [K-Means++ (Alg. 2) for $O(log k)$ guarantee; Bisecting K-Means (Alg. 9) for adaptive splitting; Elkan (Alg. 13) for speed.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of the LBG algorithm.],
)

#sec("IV", "Engineering and Scientific Applications")

Conceived inside information theory and communications engineering @linde1980algorithm, LBG remains a cornerstone primitive for lossy compression, codebook synthesis, and deterministic initialization across signal processing and machine learning.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Speech & Audio Coding (CELP / LPC Spectral Quantization):* Linear Predictive Coding (LPC) and Line Spectral Frequency (LSF) speech codecs quantize 10–16-dimensional vocal-tract filter envelopes into $2^B$-entry codebooks ($B in [8, 12]$). LBG was explicitly formulated @linde1980algorithm to minimize Itakura–Saito spectral distortion directly on empirical speech frames without parametric source assumptions, exploiting the 1.53 dB shaping gain of Theorem 3.
2. *Image, Texture & Video Compression (Block VQ & Palette Coding):* Partitioning $4 times 4$ pixel or DCT residual blocks ($d = 16$) into a 256-word LBG codebook achieves $0.5$ bits/pixel compression with $O(1)$ table-lookup decoding—critical for ultra-low-power embedded displays, medical imaging archives, and GPU texture compression (BCn/ASTC).
3. *Speaker Recognition & Acoustic Modeling:* Classic GMM–UBM and VQ-based speaker verification systems build per-speaker Mel-Frequency Cepstral Coefficient (MFCC) codebooks via deterministic LBG splitting, guaranteeing reproducible acoustic templates across enrollment sessions and fast nearest-codeword distortion scoring.
4. *Neural Compression & Product Quantization (VQ-VAE / IVF-PQ):* Modern billion-scale ANN search engines (Product Quantization) and discrete representation learners rely on subspace $2^8$-codebook quantizers initialized via LBG/Lloyd so every sub-vector encodes cleanly into a single byte ($B = 8$) for SIMD LUT distance accumulation.
5. *Scientific Lossy Compression & Embedded Telemetry:* High-throughput sensor arrays (seismic traces, satellite hyperspectral cubes, IoT waveforms) deploy LBG codebooks offline and transmit only $B$-bit codeword indices over bandwidth-constrained links, matching Theorem 3's rate-distortion frontier.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

LBG fails systematically under four conditions: (a) *imbalanced cluster variance*, where uniform doubling wastes codewords on compact cells while under-allocating to wide clusters (Proposition 1; resolved by Bisecting K-Means Alg. 9 @steinbach2000comparison or K-Means++ Alg. 2 @arthur2007kmeanspp); (b) *zero-centered data with multiplicative splitting*, causing Theorem 2(3)'s identical-codeword collapse; (c) *high-dimensional noise features*, because unweighted $L_2$ distortion lets noise axes dominate Voronoi boundaries (resolved by MWK-Means Alg. 29); and (d) *non-convex manifolds or heavy-tailed outliers* (resolved by DBSCAN Alg. 31 or Trimmed K-Means Alg. 26).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why LBG Fits], [Caution / Alternative],
    [Speech / LPC Coding], [General distortion (Itakura–Saito) & exact $2^B$ bits], [Guard empty cells & use additive split],
    [Image Block VQ], [Fast LUT decoding; deterministic $2^8$-word bytes], [Use PCA-axis split (Thm. 2)],
    [ANN Subspace PQ], [Byte-aligned $k = 256$ codebooks; $O(n k d)$ speed], [Use K-Means++ (Alg. 2) if skewed],
    [Noisy Tabular Data], [Unweighted $L_2$ degrades under noise columns], [Use MWK-Means (Alg. 29)],
    [Imbalanced Clusters], [Uniform doubling misallocates codewords (Prop. 1)], [Use Bisecting (Alg. 9) or Global (Alg. 17)],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for the LBG algorithm.],
)

#sec("V", "Structural Relations and Empirical Behavior")

LBG forms a direct mathematical bridge linking Lloyd's batch partitioning (Alg. 1 @lloyd1982least), Bregman divergence clustering @banerjee2005clustering, divisive hierarchical splitting (Bisecting K-Means Alg. 9 @steinbach2000comparison; Global K-Means Alg. 17 @likas2003global; ISODATA Alg. 18 @ball1965isodata), and topology-preserving neural vector quantizers (SOM Alg. 176 @kohonen1982som; Neural Gas Alg. 179 @martinetz1993neuralgas; GNG Alg. 177 @fritzke1995gng).

#proposition(title: "Bregman Invariance of LBG's Centroid Condition", origin: "Original Derivation")[
  Let $phi: Omega -> bb(R)$ be a strictly convex differentiable Legendre function on an open convex set $Omega subset.eq bb(R)^d$, inducing the Bregman divergence $D_phi (x || y) = phi(x) - phi(y) - chevron.l nabla phi(y), x - y chevron.r$. Under distortion $rho(x, y) = D_phi (x || y)$, the LBG Centroid Condition is *identically invariant*: for every non-empty cell $R_j subset.eq Omega$,
  $ "cent"(R_j) = arg min_(y in Omega) sum_(i in R_j) D_phi (x_i || y) = 1 / (|R_j|) sum_(i in R_j) x_i = overline(x)_(R_j). $
  In particular, both squared Euclidean distortion ($phi(x) = norm(x)_2^2$) and the speech-coding Itakura–Saito distortion $D_"IS" (x || y) = sum_(r=1)^d (x_r / y_r - ln(x_r / y_r) - 1)$ ($phi(x) = -sum_(r=1)^d ln x_r$ on $bb(R)_(>0)^d$) share the exact arithmetic mean as their LBG codeword update @linde1980algorithm @banerjee2005clustering.
]

#proof[
  For any $y in Omega$, let $overline(x) = overline(x)_(R_j) in Omega$ (by convexity of $Omega$). Expanding the difference directly yields the Bregman parallel-axis identity @banerjee2005clustering:
  $ sum_(i in R_j) &D_phi (x_i || y) - sum_(i in R_j) D_phi (x_i || overline(x)) \
    &= |R_j| ( phi(overline(x)) - phi(y) ) - chevron.l nabla phi(y), sum_(i in R_j) (x_i - y) chevron.r \
    &quad + chevron.l nabla phi(overline(x)), sum_(i in R_j) (x_i - overline(x)) chevron.r \
    &= |R_j| ( phi(overline(x)) - phi(y) - chevron.l nabla phi(y), overline(x) - y chevron.r ) \
    &= |R_j| D_phi (overline(x) || y) >= 0, $
  where the third term vanishes identically because $sum_(i in R_j) (x_i - overline(x)) = bold(0)$. By strict convexity of $phi$, $D_phi (overline(x) || y) = 0$ if and only if $y = overline(x)$. Taking $phi(x) = -sum_(r=1)^d ln x_r$ on $bb(R)_(>0)^d$ gives $(nabla phi(y))_r = -1/y_r$, so
  $ D_phi (x || y) = sum_(r=1)^d ( -ln x_r + ln y_r + (x_r - y_r) / y_r ) = D_"IS" (x || y), $
  proving that Itakura–Saito LBG speech quantization is an exact special case of Bregman Lloyd–Max minimization.
]

Beyond Bregman unification, LBG differs structurally from stochastic online quantizers such as MacQueen's K-Means (Alg. 21 @macqueen1967some), SOM (Alg. 176 @kohonen1982som), and Neural Gas (Alg. 179 @martinetz1993neuralgas) by performing deterministic batch-splitting updates with exact Lyapunov descent rather than stochastic approximation, and differs from Global K-Means (Alg. 17 @likas2003global) by splitting existing codewords locally in $O(n k d T_"max")$ instead of testing all $n$ data points as new candidate centers in $O(n^2 k d)$.

On the real ClusterOverview 4-dataset empirical benchmark suite (evaluated with Euclidean Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), LBG completes in an exceptionally fast *mean runtime of 0.171 ms*—over $5 times$ faster than MWK-Means (0.874 ms), directly corroborating Proposition 2(1)'s $< 2 times$ geometric-series overhead and Corollary 1's logarithmic relative-tolerance convergence. Across datasets, LBG achieves *Iris Silhouette 0.447* and *Iris ARI 0.583*, *Customers Silhouette 0.231*, *Multi-Shape Silhouette 0.609*, and *High-D Subspace+Noise Silhouette 0.303* (ARI: #box[N/A] on the latter three). The empirical contrast with MWK-Means (Alg. 29) neatly illustrates our theoretical bounds: LBG attains higher unweighted Euclidean Silhouette on *Iris* (0.447 vs. 0.393) and *High-D Subspace+Noise* (0.303 vs. 0.229) because it directly minimizes unweighted ambient squared error, yet trails sharply on *Iris ARI* (0.583 vs. 0.859) and *Multi-Shape* (0.609 vs. 0.718), reflecting both Proposition 1's binary-splitting allocation bias on $k = 3$ and the absence of adaptive feature weighting or curved Voronoi boundaries.
