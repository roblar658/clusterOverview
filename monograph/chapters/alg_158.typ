#import "../template/ieee.typ": *

#alg_entry(
  id: 158,
  name: "PLAID",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2002",
  paradigm: "Additive Overlapping Layer Matrix Decomposition",
  complexity: $O(K T n p)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A fundamental limitation of heuristic submatrix-extraction methods such
as Cheng-Church (Alg. 157) @cheng2000biclustering is destructive masking:
once a bicluster is discovered, overwriting its entries with random noise
corrupts any overlapping module sharing those rows and columns. In gene
expression biology, however, a single gene can be simultaneously
up-regulated by two distinct transcription factors under overlapping
experimental conditions, so its observed log-expression is naturally the
_additive superposition_ of multiple active biological processes. To
model overlapping biclusters generatively without destructive masking,
Lazzeroni and Owen @lazzeroni2002plaid introduced the *PLAID* model,
representing a real data matrix $bold(Y) = [Y_(i j)] in bb(R)^(n times p)$
as a sum of a background layer and $K$ overlapping two-way additive
bicluster layers plus noise:
$
  Y_(i j) = theta_(i j 0) + sum_(k=1)^K theta_(i j k) rho_(i k) kappa_(j k)
  + epsilon_(i j),
$
where $rho_(i k) in {0, 1}$ and $kappa_(j k) in {0, 1}$ are binary row
and column membership indicators for layer $k$ (with $rho_(i 0) = kappa_(j 0) = 1$
for the global background layer $k = 0$), and each layer effect admits an
additive two-way ANOVA parameterization:
$
  theta_(i j k) = mu_k + alpha_(i k) + beta_(j k),
$
subject to weighted sum-to-zero identifiability constraints
$
  sum_(i=1)^n rho_(i k)^2 alpha_(i k) = 0,
  quad sum_(j=1)^p kappa_(j k)^2 beta_(j k) = 0.
$
Setting $alpha_(i k) = beta_(j k) = 0$ yields constant-shift layers
($theta_(i j k) = mu_k$); retaining only $alpha_(i k)$ or $beta_(j k)$
yields row-specific or column-specific additive layers.

#definition(title: "PLAID Additive Layer Model and Least-Squares Objective")[
  A _$K$-layer PLAID decomposition_ of $bold(Y) in bb(R)^(n times p)$
  consists of layer parameters
  $Theta = {(mu_k, bold(alpha)_k, bold(beta)_k)}_(k=0)^K$ and membership
  matrices $bold(R) = [rho_(i k)] in {0, 1}^(n times K)$,
  $bold(C) = [kappa_(j k)] in {0, 1}^(p times K)$ minimizing the global
  Frobenius residual energy
  $
    Q(Theta, bold(R), bold(C)) = 1 / 2 sum_(i=1)^n sum_(j=1)^p
    (Y_(i j) - theta_(i j 0) - sum_(k=1)^K theta_(i j k) rho_(i k) kappa_(j k))^2.
  $
  The support of layer $k$ is the bicluster rectangle
  $I_k times J_k$ where $I_k = {i : rho_(i k) = 1}$ and
  $J_k = {j : kappa_(j k) = 1}$.
]

Direct global minimization over binary $(bold(R), bold(C))$ is
combinatorially intractable. Lazzeroni and Owen @lazzeroni2002plaid fit
layers sequentially via forward stagewise residual extraction with
continuous relaxation $rho_(i k), kappa_(j k) in [0, 1]$. At stage $k$,
given the current residual matrix
$Z_(i j)^((k)) = Y_(i j) - theta_(i j 0) - sum_(l=1)^(k-1) theta_(i j l) rho_(i l) kappa_(j l)$,
the algorithm minimizes the single-layer objective
$
  Q_k (mu_k, bold(alpha)_k, bold(beta)_k, bold(rho)_k, bold(kappa)_k)
  = 1 / 2 sum_(i=1)^n sum_(j=1)^p
  (Z_(i j)^((k)) - (mu_k + alpha_(i k) + beta_(j k)) rho_(i k) kappa_(j k))^2
$
by alternating closed-form parameter updates with projected
box-constrained membership updates over $[0, 1]^n$ and $[0, 1]^p$,
gradually rounding (releasing) memberships to ${0, 1}$, refitting
$Theta_k$ on the binary support, and retaining layer $k$ only if its
explained sum of squares $"SS"_k$ exceeds a permutation-test significance
threshold. Optional Gauss–Seidel backfitting re-optimizes previously
extracted layers against updated residuals.

#algbox(
  title: "PLAID Sequential Layer Extraction with Backfitting",
  inputs: [Data $bold(Y) in bb(R)^(n times p)$, max layers $K$, inner cycles $T$, release threshold $tau in [0.5, 1)$, permutations $B$, significance level $gamma$.],
  outputs: [Background $hat(bold(L))_0$ and accepted additive layers ${(hat(Theta)_k, bold(rho)_k, bold(kappa)_k)}_(k=1)^(K')$.],
)[
  1. Fit background layer $hat(theta)_(i j 0) = hat(mu)_0 + hat(alpha)_(i 0) + hat(beta)_(j 0)$ on $bold(Y)$ and initialize residual $bold(Z)^((1)) <- bold(Y) - hat(bold(L))_0$.
  2. For $k = 1, dots, K$:
  3. *SVD Initialization:* Extract leading singular pair $(bold(u)_1, bold(v)_1)$ of $bold(Z)^((k))$ and initialize $bold(rho)_k, bold(kappa)_k in [0, 1]^n times [0, 1]^p$ (e.g., centered sign-aligned thresholds).
  4. *Alternating Continuous Relaxation:* For $t = 1, dots, T$:
     - Let $R_2 = sum_i rho_(i k)^2$, $C_2 = sum_j kappa_(j k)^2$; set $hat(mu)_k = (sum_(i,j) rho_(i k) kappa_(j k) Z_(i j)^((k))) / (R_2 C_2)$, $hat(alpha)_(i k) = (sum_j kappa_(j k) Z_(i j)^((k))) / (rho_(i k) C_2) - hat(mu)_k$ ($rho_(i k) > 0$), and $hat(beta)_(j k) = (sum_i rho_(i k) Z_(i j)^((k))) / (kappa_(j k) R_2) - hat(mu)_k$ ($kappa_(j k) > 0$).
     - Update $rho_(i k) <- Pi_([0,1])((sum_j hat(theta)_(i j k) kappa_(j k) Z_(i j)^((k))) / (sum_j hat(theta)_(i j k)^2 kappa_(j k)^2))$ and $kappa_(j k) <- Pi_([0,1])((sum_i hat(theta)_(i j k) rho_(i k) Z_(i j)^((k))) / (sum_i hat(theta)_(i j k)^2 rho_(i k)^2))$.
  5. *Binary Release & Refit:* Threshold $rho_(i k) <- bb(I)(rho_(i k) >= tau)$, $kappa_(j k) <- bb(I)(kappa_(j k) >= tau)$, refit $(hat(mu)_k, hat(bold(alpha))_k, hat(bold(beta))_k)$ on $I_k times J_k$, and compute $"SS"_k = sum_(i,j) (hat(theta)_(i j k) rho_(i k) kappa_(j k))^2$.
  6. *Permutation Pruning:* Compare $"SS"_k$ against $B$ permuted residuals of $bold(Z)^((k))$; stop if $p$-value $> gamma$. Otherwise set $bold(Z)^((k+1)) <- bold(Z)^((k)) - hat(bold(L))_k$ (and optionally run backfitting over $l = 0, dots, k$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees for PLAID: the exact orthogonal
ANOVA projection structure of each continuous layer, monotone Lyapunov
descent of block-coordinate updates and backfitting, algebraic
identifiability conditions alongside an explicit non-uniqueness
counterexample, sharp rank-2 Eckart–Young spectral bounds, and exact
finite-sample permutation Type-I error control.

#subsec("A", "Weighted Orthogonal ANOVA Projection of a PLAID Layer")

#lemma(title: "Weighted Orthogonal ANOVA Decomposition and Closed-Form Optimality", origin: "Original Derivation")[
  Fix any nonzero membership vectors $bold(rho) in [0, 1]^n$,
  $bold(kappa) in [0, 1]^p$ with $R_2 = norm(bold(rho))_2^2 > 0$ and
  $C_2 = norm(bold(kappa))_2^2 > 0$. Define the weighted additive layer
  subspace $cal(V)(bold(rho), bold(kappa)) = {[(mu + alpha_i + beta_j) rho_i kappa_j]_(i,j)}$.
  Then:
  1. Under the constraints $sum_i rho_i^2 alpha_i = 0$ and
     $sum_j kappa_j^2 beta_j = 0$, the three component matrices
     $mu bold(rho) bold(kappa)^top$, $"diag"(bold(rho)) bold(alpha) bold(kappa)^top$,
     and $bold(rho) bold(beta)^top "diag"(bold(kappa))$ are mutually
     Frobenius-orthogonal.
  2. The unique least-squares estimates $(hat(mu), hat(bold(alpha)), hat(bold(beta)))$
     minimizing $Q_k$ take the closed form of Step 4, satisfy both
     constraints identically, and yield exact Pythagorean energy decay
     $
       norm(bold(Z)^((k)) - hat(bold(L))_k)_F^2
       = norm(bold(Z)^((k)))_F^2 - "SS"_k,
     $
     where $"SS"_k = hat(mu)^2 R_2 C_2 + C_2 sum_i rho_i^2 hat(alpha)_i^2 + R_2 sum_j kappa_j^2 hat(beta)_j^2$.
]

#proof[
  For (1), compute the three pairwise Frobenius inner products:
  $
    chevron.l mu bold(rho) bold(kappa)^top,
    "diag"(bold(rho)) bold(alpha) bold(kappa)^top chevron.r_F
    &= mu (sum_(i=1)^n rho_i^2 alpha_i) (sum_(j=1)^p kappa_j^2) = 0, \
    chevron.l "diag"(bold(rho)) bold(alpha) bold(kappa)^top,
    bold(rho) bold(beta)^top "diag"(bold(kappa)) chevron.r_F
    &= (sum_(i=1)^n rho_i^2 alpha_i) (sum_(j=1)^p kappa_j^2 beta_j) = 0,
  $
  and symmetrically $chevron.l mu bold(rho) bold(kappa)^top, bold(rho) bold(beta)^top "diag"(bold(kappa)) chevron.r_F = 0$.
  Hence $cal(V)(bold(rho), bold(kappa))$ splits as an orthogonal direct
  sum of three subspaces.

  For (2), because the map $(mu, bold(alpha), bold(beta)) |-> bold(L) in cal(V)(bold(rho), bold(kappa))$
  is linear, $Q_k$ is convex quadratic and its stationarity equations
  decouple across the three orthogonal subspaces. Setting
  $partial Q_k / partial mu = 0$ yields
  $sum_(i,j) rho_i kappa_j Z_(i j)^((k)) = mu R_2 C_2 + C_2 (sum_i rho_i^2 alpha_i) + R_2 (sum_j kappa_j^2 beta_j) = mu R_2 C_2$,
  giving $hat(mu) = (R_2 C_2)^(-1) sum_(i,j) rho_i kappa_j Z_(i j)^((k))$.
  Setting $partial Q_k / partial alpha_i = 0$ for any $i$ with $rho_i > 0$
  gives $rho_i sum_j kappa_j Z_(i j)^((k)) = rho_i^2 (hat(mu) + alpha_i) C_2$,
  so $hat(alpha)_i = (rho_i C_2)^(-1) sum_j kappa_j Z_(i j)^((k)) - hat(mu)$
  (with $hat(alpha)_i = 0$ when $rho_i = 0$). Summing $rho_i^2 hat(alpha)_i$
  gives $C_2^(-1) sum_(i,j) rho_i kappa_j Z_(i j)^((k)) - R_2 hat(mu) = 0$,
  verifying the constraint identically; the derivation for $hat(beta)_j$
  is symmetric. Because $hat(bold(L))_k$ is the orthogonal projection of
  $bold(Z)^((k))$ onto $cal(V)(bold(rho), bold(kappa))$,
  $chevron.l bold(Z)^((k)) - hat(bold(L))_k, hat(bold(L))_k chevron.r_F = 0$,
  yielding $norm(bold(Z)^((k)) - hat(bold(L))_k)_F^2 = norm(bold(Z)^((k)))_F^2 - norm(hat(bold(L))_k)_F^2$.
  Applying pairwise orthogonality from (1) expands
  $"SS"_k = norm(hat(bold(L))_k)_F^2$ into the stated three-term ANOVA sum.
]

#subsec("B", "Block-Coordinate Monotonicity and Backfitting Convergence")

#theorem(title: "Lyapunov Monotonicity of Continuous Relaxation and Multi-Layer Backfitting", origin: "Original Proof")[
  1. *(Single-Layer Continuous Relaxation)* During Step 4 of PLAID, the
     alternating block updates $Theta_k -> bold(rho)_k -> bold(kappa)_k$
     satisfy $Q_k^((t+1)) <= Q_k^((t))$ monotonically for every $t >= 0$.
     Every accumulation point whose membership denominators
     $w_(i, rho) = sum_j hat(theta)_(i j k)^2 kappa_(j k)^2$ and
     $w_(j, kappa) = sum_i hat(theta)_(i j k)^2 rho_(i k)^2$ are positive
     is a block-coordinate KKT stationary point of $Q_k$.
  2. *(Multi-Layer Backfitting Descent)* Given fixed binary supports
     ${(I_k, J_k)}_(k=0)^K$, cyclic Gauss–Seidel backfitting updating
     $Theta_k$ against $bold(Y) - sum_(l != k) bold(L)_l$ decreases the
     global objective $Q(Theta, bold(R), bold(C))$ monotonically and
     converges to a global minimizer of $Q(dot, bold(R), bold(C))$.
]

#proof[
  For (1), by Lemma 1 the parameter update $Theta_k^((t+1))$ globally
  minimizes $Q_k (dot, bold(rho)_k^((t)), bold(kappa)_k^((t)))$. Fixing
  $(Theta_k^((t+1)), bold(kappa)_k^((t)))$, $Q_k$ separates across rows
  as $1/2 sum_(i=1)^n (w_(i, rho) rho_(i k)^2 - 2 b_(i, rho) rho_(i k) + c_i)$
  with $w_(i, rho) = sum_j hat(theta)_(i j k)^2 kappa_(j k)^2 >= 0$ and
  $b_(i, rho) = sum_j hat(theta)_(i j k) kappa_(j k) Z_(i j)^((k))$.
  When $w_(i, rho) > 0$, this univariate quadratic is strictly convex on
  $[0, 1]$ and uniquely minimized by the Euclidean projection
  $Pi_([0,1])(b_(i, rho) / w_(i, rho)) = min(1, max(0, b_(i, rho) / w_(i, rho)))$.
  The column update for $kappa_(j k) in [0, 1]$ is symmetric. Hence each
  substep decreases $Q_k >= 0$, so ${Q_k^((t))}$ converges. Because each
  block optimization is uniquely attained whenever $w_(i, rho), w_(j, kappa) > 0$
  over compact convex domains, standard block-coordinate optimality
  guarantees KKT stationarity at every limit point.

  For (2), with binary $(bold(R), bold(C))$ fixed, the global residual
  $Q(Theta, bold(R), bold(C))$ is a convex quadratic function of the
  stacked parameter vector $Theta$, and each layer subspace
  $cal(V)(bold(rho)_k, bold(kappa)_k)$ is closed and linear. Updating
  $Theta_k$ via Lemma 1 executes an exact orthogonal projection of the
  partial residual onto $cal(V)(bold(rho)_k, bold(kappa)_k)$, decreasing
  $Q$ at every step unless already orthogonal to every layer subspace.
  By cyclic alternating projection onto a finite collection of linear
  subspaces, the residual converges to the orthogonal projection of
  $bold(Y)$ onto $sum_(k=0)^K cal(V)(bold(rho)_k, bold(kappa)_k)$, i.e.
  a global minimizer of $Q(dot, bold(R), bold(C))$.
]

#subsec("C", "Algebraic Identifiability and Overlap Non-Uniqueness")

#theorem(title: "Identifiability Criterion and Non-Uniqueness Counterexample for Overlapping Layers", origin: "Original Proof")[
  1. *(Constructive Non-Uniqueness Counterexample)* Even in the noiseless
     constant-layer regime ($theta_(i j k) = mu_k$, $bold(epsilon) = bold(0)$)
     with $K = 2$, the PLAID representation is *not* unique without
     structural regularization: any nested pair $I_2 subset I_1$,
     $J_2 = J_1$ admits at least three distinct exact ($Q = 0$) layer
     decompositions with identical matrix $bold(Y)$.
  2. *(Exact Amplitude Identifiability Given Supports)* For fixed
     constant-layer supports ${I_k times J_k}_(k=1)^K$ with indicator
     matrices $bold(M)_k = bold(1)_(I_k) bold(1)_(J_k)^top$ and
     $bold(M)_0 = bold(1)_n bold(1)_p^top$, the amplitude vector
     $(mu_0, mu_1, dots, mu_K)$ is uniquely recoverable from $bold(Y)$
     if and only if ${bold(M)_0, bold(M)_1, dots, bold(M)_K}$ is linearly
     independent in $bb(R)^(n times p)$. In particular, identifiability is
     guaranteed whenever there exists a background-only cell
     $(i_0, j_0) in.not union_(k=1)^K (I_k times J_k)$ and each layer $k$
     contains an exclusive cell
     $(i_k, j_k) in (I_k times J_k) without union_(l != k) (I_l times J_l)$.
]

#proof[
  For (1), partition $I_1 = A union B$ with nonempty disjoint $A, B$ and
  fix nonempty $J subset.eq {1, dots, p}$. Let $bold(Y)$ take value
  $mu_1$ on $A times J$, $mu_1 + mu_2$ on $B times J$, and $0$ elsewhere
  ($mu_1, mu_2 != 0$). Then $bold(Y)$ is represented with $Q = 0$ by:
  (a) nested overlapping layers $(I_1, J; mu_1)$ and $(B, J; mu_2)$;
  (b) disjoint layers $(A, J; mu_1)$ and $(B, J; mu_1 + mu_2)$; and
  (c) complementary nested layers $(I_1, J; mu_1 + mu_2)$ and
  $(A, J; -mu_2)$. Thus overlapping and disjoint layer structures are
  observationally equivalent under pure least squares.

  For (2), the noiseless constant-layer model is the linear system
  $bold(Y) = sum_(k=0)^K mu_k bold(M)_k$, whose solution
  $(mu_0, dots, mu_K)$ is unique iff $sum_(k=0)^K c_k bold(M)_k = bold(0) => c_0 = dots = c_K = 0$,
  i.e. linear independence of ${bold(M)_k}_(k=0)^K$. Under the
  exclusive-cell condition, evaluating $sum_(k=0)^K c_k bold(M)_k = bold(0)$
  at $(i_0, j_0)$ forces $c_0 = 0$; subsequently evaluating at each
  exclusive cell $(i_k, j_k)$ yields $c_0 + c_k = 0 => c_k = 0$ for all
  $k = 1, dots, K$, proving linear independence.
]

#subsec("D", "Spectral Rank-2 Ceiling and Exact Permutation Calibration")

#proposition(title: "Rank-2 Spectral Energy Ceiling and Permutation Type-I Error Control", origin: "Classical — Cited")[
  Let $sigma_1 >= sigma_2 >= dots$ denote the singular values of the
  current residual matrix $bold(Z)^((k)) in bb(R)^(n times p)$.
  1. *(Spectral Ceiling)* Every fitted PLAID layer
     $hat(bold(L))_k = [(hat(mu)_k + hat(alpha)_(i k) + hat(beta)_(j k)) rho_(i k) kappa_(j k)]$
     has matrix rank at most $2$ (rank at most $1$ when $hat(bold(alpha))_k = bold(0)$
     or $hat(bold(beta))_k = bold(0)$), and therefore its explained sum
     of squares obeys @lazzeroni2002plaid @hastie2009elements
     $
       "SS"_k <= sigma_1^2 (bold(Z)^((k))) + sigma_2^2 (bold(Z)^((k))),
     $
     tightening to $"SS"_k <= sigma_1^2 (bold(Z)^((k)))$ for
     constant-shift, row-only, or column-only layers.
  2. *(Exact Permutation Type-I Control)* Under the null hypothesis $H_0$
     that entries of $bold(Z)^((k))$ are exchangeable under a finite
     permutation group $Pi$ (e.g., independent row/column permutations),
     accepting layer $k$ only when $"SS"_k (bold(Z)^((k))) > max_(1 <= b <= B) "SS"_k (pi_b bold(Z)^((k)))$
     for $B$ i.i.d. uniform draws $pi_b in Pi$ bounds the false-discovery
     probability by $"Pr"_(H_0) ("accept") <= 1 / (B + 1)$ @lazzeroni2002plaid.
]

#proof[
  For (1), rewrite $hat(bold(L))_k$ explicitly as the sum of two rank-1
  outer products:
  $
    hat(bold(L))_k = bold(rho)_k (hat(mu)_k bold(kappa)_k + "diag"(bold(kappa)_k) hat(bold(beta))_k)^top
    + ("diag"(bold(rho)_k) hat(bold(alpha))_k) bold(kappa)_k^top.
  $
  Hence $"rank"(hat(bold(L))_k) <= 2$ (and $<= 1$ whenever either
  $hat(bold(alpha))_k = bold(0)$ or $hat(bold(beta))_k = bold(0)$). By
  the Eckart–Young–Mirsky theorem @hastie2009elements, any rank-$r$
  matrix $bold(L)$ obeys
  $norm(bold(Z)^((k)) - bold(L))_F^2 >= sum_(s > r) sigma_s^2 (bold(Z)^((k)))$.
  Combining this lower bound with the Pythagorean identity of Lemma 1
  ($"SS"_k = norm(bold(Z)^((k)))_F^2 - norm(bold(Z)^((k)) - hat(bold(L))_k)_F^2$)
  gives $"SS"_k <= sum_(s=1)^r sigma_s^2 (bold(Z)^((k)))$ with $r = 2$
  (or $r = 1$ in the one-way/constant submodels).

  For (2), under $H_0$ the $B + 1$ random variables
  $(S_0, S_1, dots, S_B) = ("SS"_k (bold(Z)^((k))), "SS"_k (pi_1 bold(Z)^((k))), dots, "SS"_k (pi_B bold(Z)^((k))))$
  are jointly exchangeable. At most one index in ${0, 1, dots, B}$ can be
  strictly larger than all others, so by symmetry across the $B + 1$
  positions, $"Pr"_(H_0) (S_0 > max_(1 <= b <= B) S_b) <= 1 / (B + 1)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each inner iteration of Step 4 evaluates row and column weighted sums of
$bold(Z)^((k))$ in $O(n p)$ flops, followed by $O(n + p)$ scalar updates
for $hat(mu)_k, hat(bold(alpha))_k, hat(bold(beta))_k, bold(rho)_k, bold(kappa)_k$.
Extracting $K$ layers with $T$ inner relaxation iterations each (plus
truncated rank-1 power-method SVD initialization in $O(n p)$ per layer)
therefore costs $O(K T n p)$ flops and $O(n p)$ working memory for the
residual matrix—strictly linear in the matrix entry count $n p$ and
dramatically faster than Cheng-Church's $O(K n p (n + p))$ worst-case
deletion loop! Adding $B$ permutation tests scales runtime by factor
$(B + 1)$ to $O(K (B + 1) T n p)$, which is embarrassingly parallel
across permutations. Optional $R_"bf"$ full Gauss–Seidel backfitting
sweeps add $O(R_"bf" K n p)$ flops.

Geometrically, PLAID constructs a sparse union of axis-aligned rank-2
(or rank-1) rectangular blocks superposed over a global additive
background. Unlike disjoint clustering, rows and columns may belong to
zero, one, or multiple layers simultaneously, making PLAID naturally
suited to overlapping combinatorial latent factors. It is equivariant to
global matrix rescaling and invariant to row/column permutations, but
sensitive to column scale heterogeneity and feature rotations.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(K T n p)$ without permutations; $O(K (B + 1) T n p)$ with $B$ permutation tests],
    [Space Complexity], [$O(n p + K (n + p))$ (residual matrix plus sparse layer factors)],
    [Generative Model], [Additive superposition $bold(Y) = bold(L)_0 + sum_(k=1)^K bold(L)_k + bold(epsilon)$ with $"rank"(bold(L)_k) <= 2$],
    [Overlap Handling], [Native additive superposition (no destructive random masking)],
    [Hyperparameters], [Max layers $K$, relaxation cycles $T$, release threshold $tau$, permutations $B$],
    [Statistical Pruning], [Exact exchangeability permutation test with Type-I error $<= 1 / (B + 1)$],
  ),
  caption: [Complexity, generative structure, and operational profile of PLAID.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Combinatorial Gene Regulation in Microarray & Bulk RNA-Seq:* Genes
   controlled by multiple transcription factors exhibit additive
   log-expression shifts whenever two regulatory programs activate
   concurrently under shared perturbations. PLAID @lazzeroni2002plaid
   disentangles overlapping regulons directly via residual subtraction
   while absorbing global housekeeping baselines into layer $k = 0$.
2. *Multi-Factor Consumer Basket and Retail Demand Decomposition:* Store
   sales matrices (SKU $times$ week) combine baseline seasonal demand
   (background $k = 0$) with overlapping promotional campaigns, holiday
   spikes, and regional shocks affecting specific product subsets over
   specific weeks—mirroring PLAID's overlapping rectangular additive
   layers.
3. *Neuroscience Multi-Task fMRI Activation Mapping:* Cortical regions
   (rows) participate in overlapping cognitive networks across behavioral
   tasks (columns). Additive layer decomposition separates shared sensory
   background responses from task-specific overlapping functional
   subnetworks without forcing mutually exclusive parcellation.
4. *Industrial Sensor-Array Fault Superposition:* In chemical plants and
   power grids, simultaneous localized faults induce additive anomaly
   signatures over distinct subsets of sensors and time intervals.
   Sequential residual extraction isolates concurrent fault signatures
   even when sensor subsets overlap.
5. *Epidemiological and Environmental Spatio-Temporal Anomalies:*
   Region-by-time incidence matrices reflect endemic background rates
   plus overlapping localized outbreaks or pollution episodes; PLAID
   quantifies each outbreak's spatial footprint ($I_k$), temporal window
   ($J_k$), and excess magnitude ($theta_(i j k)$) with permutation-tested
   statistical significance.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure modes govern PLAID in practice. First, *overlap
non-identifiability* (Theorem 2) allows nested or collinear rectangular
layers to trade amplitude arbitrarily, causing greedy forward extraction
to split or merge true blocks depending on initialization order. Second,
*greedy stagewise residual capture*: if an early layer partially absorbs
variance from a strongly overlapping secondary layer before backfitting,
binary thresholding ($tau$) locks in suboptimal supports. Third,
*binary-release discontinuity*: hard thresholding continuous
$rho_(i k), kappa_(j k) in [0, 1]$ to ${0, 1}$ can abruptly increase $Q_k$
or produce empty layers if $tau$ is too strict. Fourth, *non-additive
multiplicative interactions* ($Y_(i j) = u_i v_j$) exceed PLAID's
two-way additive parameterization unless log-transformed. Fifth,
*global Euclidean partition mismatch*: rows unassigned to any layer fall
into a heterogeneous background pool ($k = 0$), breaking geometric
compactness metrics such as Silhouette.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Overlapping Log-Expression Modules], [Strong: native additive superposition + $p$-values], [Enable backfitting to mitigate greedy order bias],
    [Multiplicative Factor Data], [Moderate–Weak without log transform], [Use FABIA @hochreiter2010fabia, NMF @lee1999nmf, or Spectral Biclustering @kluger2003spectral],
    [Large Sparse Matrices ($n, p > 10^4$)], [Strong base $O(K T n p)$; permutation test costly], [Parallelize permutations or use asymptotic SVD screening],
    [Disjoint Geometric Point Clouds], [Poor: background layer merges unrelated regions], [Use KFCM (Alg. 155), sFCM (Alg. 156), or k-means @lloyd1982least],
    [Small-Column Tabular Benchmarks], [Poor: limited column subset structure ($p <= 5$)], [Use standard partitioning or density-based clustering],
  ),
  caption: [Domain suitability and failure-mode matrix for PLAID.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PLAID bridges residue-based biclustering (Cheng-Church, Alg. 157) and
low-rank matrix factorization (SVD, NMF @lee1999nmf @lee2001algorithms,
and sparse factor analysis FABIA @hochreiter2010fabia). The following
original proposition establishes the exact algebraic embedding between
Cheng-Church and PLAID.

#proposition(title: "Exact Embedding of Cheng-Church Biclusters into PLAID Layers", origin: "Original Derivation")[
  Let $(I, J)$ be any submatrix of residual matrix $bold(Z)^((k))$ with
  Cheng-Church Mean Squared Residue $H(I, J)$ @cheng2000biclustering.
  Fitting a single binary PLAID layer supported on $bold(rho) = bold(1)_I$,
  $bold(kappa) = bold(1)_J$ achieves:
  1. Exact in-block residual energy
     $sum_(i in I) sum_(j in J) (Z_(i j)^((k)) - hat(theta)_(i j k))^2 = |I| |J| H(I, J)$,
     so zero-MSR Cheng-Church blocks are reproduced with zero error.
  2. Exact global variance reduction
     $"SS"_k = norm(bold(Z)_(I J)^((k)))_F^2 - |I| |J| H(I, J)$, showing
     that PLAID balances *high signal energy* $norm(bold(Z)_(I J)^((k)))_F^2$
     against *low internal residue* $|I| |J| H(I, J)$, automatically
     rejecting flat zero-variance blocks that satisfy $H = 0$ trivially.
]

#proof[
  On binary support $bold(rho) = bold(1)_I$, $bold(kappa) = bold(1)_J$,
  Lemma 1 yields $R_2 = |I|$, $C_2 = |J|$, $hat(mu)_k = z_(I J)^((k))$,
  $hat(alpha)_(i k) = z_(i J)^((k)) - z_(I J)^((k))$ ($i in I$), and
  $hat(beta)_(j k) = z_(I j)^((k)) - z_(I J)^((k))$ ($j in J$). Hence
  $hat(theta)_(i j k) = z_(i J)^((k)) + z_(I j)^((k)) - z_(I J)^((k))$ on
  $I times J$ and $0$ outside $I times J$. Therefore the in-block error
  is $Z_(i j)^((k)) - hat(theta)_(i j k) = r_(i j) (I, J)$, whose sum of
  squares equals $|I| |J| H(I, J)$ identically. Because $hat(bold(L))_k$
  vanishes outside $I times J$, the global residual energy decomposes as
  $norm(bold(Z)^((k)) - hat(bold(L))_k)_F^2 = norm(bold(Z)^((k)))_F^2 - norm(bold(Z)_(I J)^((k)))_F^2 + |I| |J| H(I, J)$.
  Subtracting from $norm(bold(Z)^((k)))_F^2$ via Lemma 1 gives
  $"SS"_k = norm(bold(Z)_(I J)^((k)))_F^2 - |I| |J| H(I, J)$. In
  particular, a constant zero block has $norm(bold(Z)_(I J)^((k)))_F^2 = 0$
  and $"SS"_k = 0$, preventing trivial flat-submatrix extraction.
]

Unlike NMF @lee1999nmf (which requires non-negative entries and
approximates $bold(Y) approx bold(W) bold(H)$ via continuous factors) and
FABIA @hochreiter2010fabia (which fits multiplicative sparse Bayesian
factors), PLAID handles signed matrices with explicit binary supports
and additive two-way effects.

*Empirical Benchmark Behavior (ClusterOverview Suite).* Across the
four-dataset ClusterOverview benchmark (Silhouette
@rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing),
PLAID is the fastest method among Algorithms 155–158 at a mean runtime
of $0.250 "ms"$—nearly $3 times$ faster than Cheng-Church ($0.707 "ms"$)
and $10 times$ faster than KFCM ($2.723 "ms"$) and sFCM ($2.409 "ms"$),
directly reflecting its lightweight $O(K T n p)$ matrix-vector loop.
When evaluated as a hard global row partition, however, PLAID records
$"Sil" = -0.177$, $"ARI" = 0.095$ on *Iris*, $"Sil" = 0.039$ on
*Customers*, $"Sil" = 0.508$ ($"ARI" = "N/A"$) on *Multi-Shape*, and
$"Sil" = 0.072$ ($"ARI" = "N/A"$) on *High-D*. The negative Silhouette on
*Iris* ($-0.177$) is an illuminating consequence of PLAID's generative
bias: rows excluded from active layers collapse into a single catch-all
background class ($k = 0$) spanning opposite sides of the feature space,
while overlapping rows forced into one hard label group points sharing
additive shifts rather than Euclidean proximity. Conversely, on
*Multi-Shape* ($"Sil" = 0.508$, outperforming Cheng-Church's $0.404$),
strong localized mean contrasts between coordinate blocks yield cohesive
dominant layers, and on *High-D* ($0.072$ vs. Cheng-Church's $0.043$),
additive residual subtraction avoids random-masking corruption,
confirming PLAID's comparative advantage for overlapping submatrix
discovery @prelic2006bimax.
