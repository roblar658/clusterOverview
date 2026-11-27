#import "../template/ieee.typ": *

#alg_entry(
  id: 130,
  name: "SBM (Stochastic Block Model)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "1983",
  paradigm: "Generative Probabilistic Graph Clustering",
  complexity: $cal(O)(T (m K + n K^2))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Combinatorial graph clustering criteria such as modularity @newman2004finding
and normalized cuts @shi2000ncut assume *assortative* communities whose
internal edge density strictly exceeds external density. Real networks,
however, frequently exhibit *disassortative* (bipartite or ecological
predator--prey), *core--periphery*, and *hierarchical role* architectures.
Introduced by Holland, Laskey, and Leinhardt @holland1983sbm, the Stochastic
Block Model (SBM) places community discovery within an exact generative
probabilistic framework capable of learning arbitrary block-interaction
topologies while supporting rigorous statistical inference and model selection
@abbe2018sbmreview.

Let $G = (V, E)$ be an undirected simple graph on $n$ vertices with symmetric
binary adjacency matrix $A in {0, 1}^(n times n)$ ($A_(i i) = 0$) and fixed
latent block count $K >= 2$.

#definition(title: "Stochastic Block Model and Mean-Field Variational ELBO")[
  Under the *$K$-block Stochastic Block Model* $"SBM"(n, pi, B)$
  @holland1983sbm, each vertex $i in {1, dots, n}$ independently draws a latent
  block membership $Z_i tilde.op "Cat"(pi)$ with mixing vector
  $pi in Delta^(K - 1)$ ($pi_k > 0$, $sum_(k=1)^K pi_k = 1$). Conditioned on
  $Z = (Z_1, dots, Z_n)$, upper-triangular edges ${A_(i j)}_(1 <= i < j <= n)$
  are mutually independent Bernoulli variables governed by a symmetric
  connectivity matrix $B = B^top in (0, 1)^(K times K)$:
  $ A_(i j) bar.v (Z_i = k, Z_j = l) tilde.op "Bern"(B_(k l)). $
  The complete-data log-likelihood is
  $ log P(A, Z bar.v pi, B) &= sum_(i=1)^n log pi_(Z_i) + sum_(1 <= i < j <= n) [ A_(i j) log B_(Z_i Z_j) \
                            & quad + (1 - A_(i j)) log(1 - B_(Z_i Z_j)) ]. $
  Approximating the intractable posterior $P(Z bar.v A, pi, B)$ by a factorized
  mean-field family $q(Z) = product_(i=1)^n "Cat"(Z_i bar.v tau_i)$
  ($tau_i in Delta^(K - 1)$) @blei2006variational defines the *Evidence Lower
  Bound* (ELBO):
  $ cal(L)(tau, pi, B) &= bb(E)_q [log P(A, Z bar.v pi, B)] + H(q) \
                       &= log P(A bar.v pi, B) - D_"KL" (q(Z) || P(Z bar.v A, pi, B)). $
]

Because $B$ is unconstrained beyond symmetry, SBM unifies assortative
($B_(k k) > B_(k l)$), disassortative ($B_(k k) < B_(k l)$), and
core--periphery ($B_(1 1) > B_(1 2) > B_(2 2)$) structures. Exact EM
@dempster1977em requires summing over $K^n$ latent assignments; Variational EM
(VEM) replaces the intractable E-step by maximizing $cal(L)$ over $tau$.

#algbox(
  title: "Sparse Mean-Field Variational EM for the Stochastic Block Model",
  inputs: [Adjacency $A in {0, 1}^(n times n)$, block count $K$, tolerance $delta$, iterations $T$],
  outputs: [Responsibilities $tau in Delta_K^n$, parameters $(hat(pi), hat(B))$, MAP labels $hat(z)$],
)[
  + *Spectral Initialization:* Compute top-$K$ eigenvectors of $A$ (or
    normalized Laplacian @ng2002spectral @vonluxburg2007tutorial), run
    $k$-means to form initial soft responsibilities $tau^((0))$, and compute
    initial $(pi^((0)), B^((0)))$.
  + *Sparse Variational E-Step (Gauss--Seidel):* Precompute
    $N_l <- sum_(j=1)^n tau_(j l)$, logit matrix
    $W_(k l) = log(B_(k l) / (1 - B_(k l)))$, and $U_(k l) = log(1 - B_(k l))$.
    Sweep sequentially over $i = 1, dots, n$:
    - Accumulate neighbor sums $S_(i l) = sum_(j in N(i)) tau_(j l)$ for all $l$.
    - Update for $k = 1, dots, K$:
      $tau_(i k) prop pi_k exp(sum_(l=1)^K [S_(i l) W_(k l) + (N_l - tau_(i l)) U_(k l)])$,
      normalize $sum_k tau_(i k) = 1$, and update $N$ incrementally in $cal(O)(K)$.
  + *Exact M-Step:* Form $S = A tau$, $E_(k l) = (tau^top S)_(k l)$, and
    $N_(k l) = N_k N_l - sum_(i=1)^n tau_(i k) tau_(i l)$. Update
    $pi_k <- N_k / n$ and $B_(k l) <- E_(k l) / N_(k l)$ (clipped to
    $[epsilon_0, 1 - epsilon_0]$).
  + Repeat Steps 2--3 until $abs(cal(L)^((t)) - cal(L)^((t-1))) < delta$ or
    $t = T$. Output $hat(z)_i = arg max_k tau_(i k)$, $(hat(pi), hat(B))$, and
    BIC / ICL score @schwarz1978bic.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four central theoretical pillars of SBM: closed-form variational
ELBO concavity and Gauss--Seidel monotonicity, the fundamental information-theoretic
detection and exact-recovery phase transitions, an original closed-form proof
of vanilla SBM's degree-sorting collapse, and a Davis--Kahan spectral
consistency bound justifying spectral VEM initialization.

#subsec("A", "Closed-Form ELBO and Monotonic Gauss--Seidel VEM Ascent")

Because self-loops are excluded ($A_(i i) = 0$), the mean-field ELBO is
bilinear across distinct vertices and strictly concave in every individual
vertex distribution $tau_i$.

#theorem(
  title: "Exact Mean-Field ELBO and Monotonic Coordinate VEM Convergence",
  origin: "Original Proof",
)[
  Define expected block-pair edge and pair counts
  $E_(k l) (tau) = sum_(i != j) tau_(i k) tau_(j l) A_(i j)$ and
  $N_(k l) (tau) = sum_(i != j) tau_(i k) tau_(j l)$. Under
  $q(Z) = product_(i=1)^n "Cat"(tau_i)$:
  + The ELBO evaluates in closed form to
    $ cal(L)(tau, pi, B) &= sum_(i=1)^n sum_(k=1)^K tau_(i k) log (pi_k / tau_(i k)) \
                         & quad + 1 / 2 sum_(k, l = 1)^K [ E_(k l) (tau) log B_(k l) \
                         & quad + (N_(k l) (tau) - E_(k l) (tau)) log(1 - B_(k l)) ]. $
  + Holding $(tau_(-i), pi, B)$ fixed, $cal(L)$ is strictly concave in
    $tau_i in Delta^(K - 1)$ with unique maximizer given by Step 2 of the
    algorithm. Holding $tau$ fixed, $cal(L)$ is concave in $(pi, B)$ with
    unique maximizer $pi_k^* = N_k / n$ and
    $B_(k l)^* = E_(k l) (tau) / N_(k l) (tau)$ ($N_(k l) > 0$). Hence
    Gauss--Seidel VEM satisfies $cal(L)^((t+1)) >= cal(L)^((t))$ at every step
    and converges.
]

#proof[
  *(1)* Because $q$ factorizes across vertices, for $i != j$ we have
  $q(Z_i = k, Z_j = l) = tau_(i k) tau_(j l)$. Taking the expectation of
  $log P(A, Z bar.v pi, B) - log q(Z)$, using the symmetry of $A$ and $B$ to
  rewrite $sum_(i < j) = (1 / 2) sum_(i != j)$, and collecting coefficients of
  $log B_(k l)$ and $log(1 - B_(k l))$ directly yields the closed-form ELBO.

  *(2)* Fix vertex $i$ and hold $(tau_(-i), pi, B)$ constant. Because the
  double sum excludes $i = j$, $tau_(i k)$ appears linearly in the expected
  edge log-likelihood: by symmetry of $(i, k) <-> (j, l)$, both $(r = i, s != i)$
  and $(s = i, r != i)$ contribute identically, cancelling the factor $1 / 2$.
  Differentiating with respect to $tau_(i k)$ gives
  $ (partial cal(L)) / (partial tau_(i k)) &= log pi_k - log tau_(i k) - 1 \
                                           & quad + sum_(j != i) sum_(l=1)^K tau_(j l) [ A_(i j) W_(k l) + U_(k l) ], $
  where $W_(k l) = log(B_(k l) / (1 - B_(k l)))$ and $U_(k l) = log(1 - B_(k l))$.
  The Hessian with respect to $tau_i$ is
  $nabla_(tau_i)^2 cal(L) = "diag"(-1 / tau_(i 1), dots, -1 / tau_(i K)) prec 0$,
  which is strictly negative definite on the interior of $Delta^(K - 1)$.
  Imposing the simplex multiplier $sum_k tau_(i k) = 1$ yields the unique
  global maximizer
  $ tau_(i k)^* prop pi_k exp(sum_(l=1)^K [S_(i l) W_(k l) + (N_l - tau_(i l)) U_(k l)]), $
  using $sum_(j != i) tau_(j l) = N_l - tau_(i l)$ and
  $sum_(j != i) tau_(j l) A_(i j) = S_(i l)$ ($A_(i i) = 0$).

  For the M-step at fixed $tau$, $sum_k N_k log pi_k$ is uniquely maximized on
  $Delta^(K - 1)$ by $pi_k^* = N_k / n$ via Gibbs' inequality. For each block
  pair $(k, l)$,
  $ (partial cal(L)) / (partial B_(k l)) = (E_(k l)) / B_(k l) - (N_(k l) - E_(k l)) / (1 - B_(k l)) = 0 <=> B_(k l)^* = E_(k l) / N_(k l), $
  with strictly negative second derivative
  $-E_(k l) / B_(k l)^2 - (N_(k l) - E_(k l)) / (1 - B_(k l))^2 < 0$. Because
  every Gauss--Seidel vertex update and M-step maximizes $cal(L)$ over its
  coordinate block, $cal(L)^((t+1)) >= cal(L)^((t))$; since
  $cal(L) <= log P(A bar.v pi, B) <= 0$, the monotone bounded sequence
  converges.
]

Note that *synchronous* (Jacobi) updates of all $tau_i$ simultaneously can
oscillate due to bilinear coupling between adjacent vertices; Gauss--Seidel
coordinate sweeps are essential for guaranteed monotonicity.

#subsec("B", "Fundamental Information-Theoretic Phase Transitions")

A landmark achievement of modern network theory is the exact characterization
of when SBM community recovery is statistically possible @abbe2018sbmreview.

#theorem(
  title: "Kesten--Stigum Detection and Chernoff--Hellinger Exact Recovery",
  origin: "Classical — Cited",
)[
  Let $G tilde.op "SBM"(n, p, q)$ be a symmetric balanced $2$-block SBM
  ($abs(V_1) = abs(V_2) = n / 2$) @abbe2018sbmreview:
  + *(Weak recovery, constant-degree regime)* For $p = a / n$, $q = b / n$
    with fixed $a, b > 0$, non-trivial correlation with $Z$ is achievable iff
    the Kesten--Stigum signal-to-noise ratio satisfies
    $ "SNR" = (a - b)^2 / (2 (a + b)) > 1. $
  + *(Exact recovery, logarithmic-degree regime)* For
    $p = alpha (ln n) / n$, $q = beta (ln n) / n$ with $alpha > beta > 0$,
    exact recovery $bb(P)(hat(Z) = Z "up to swap") -> 1$ is achievable iff
    $ D_+ (alpha, beta) := 1 / 2 (sqrt(alpha) - sqrt(beta))^2 > 1 <==> abs(sqrt(alpha) - sqrt(beta)) > sqrt(2). $
]

#proof[
  We derive the sharp Chernoff--Hellinger exact-recovery threshold directly
  via a genie-aided hypothesis test @abbe2018sbmreview. Given $Z_(-i)$ for a
  vertex $i in V_1$, let $X tilde.op "Bin"(n / 2 - 1, p)$ and
  $Y tilde.op "Bin"(n / 2, q)$ be $i$'s degree into $V_1 without {i}$ and
  $V_2$, respectively. Because $alpha > beta$, the MAP log-likelihood ratio
  misclassifies $i$ whenever $Y >= X$. For any Chernoff parameter $t > 0$,
  $ bb(P)(Y - X >= 0) &<= bb(E)[e^(t (Y - X))] \
                      &= (1 - q + q e^t)^(n / 2) (1 - p + p e^(-t))^(n / 2 - 1). $
  Applying $1 + x <= e^x$ bounds the log-probability by
  $ (n / 2) [q (e^t - 1) + p (e^(-t) - 1)] + o(1) = - ((ln n) / 2) (alpha + beta - alpha e^(-t) - beta e^t) + o(1). $
  Minimizing $alpha e^(-t) + beta e^t$ over $t > 0$ sets
  $e^t = sqrt(alpha / beta)$, yielding minimum value $2 sqrt(alpha beta)$ and
  $ bb(P)(Y >= X) <= exp(- 1 / 2 (sqrt(alpha) - sqrt(beta))^2 ln n + o(ln n)) = n^(-D_+ (alpha, beta) + o(1)). $
  By the union bound over all $n$ vertices, the expected number of genie errors
  is at most $n^(1 - D_+ (alpha, beta) + o(1)) -> 0$ whenever
  $D_+ (alpha, beta) > 1$, enabling exact recovery via spectral initialization
  plus one local majority cleanup round @abbe2018sbmreview. Conversely, when
  $D_+ (alpha, beta) < 1$, a matching Poisson second-moment bound shows that
  with probability $1 - o(1)$ both blocks contain vertices with $Y > X$,
  making even the MAP estimator fail. For weak recovery (1), the expected
  branching matrix $M = (1 / 2) mat(a, b; b, a)$ has eigenvalues
  $lambda_1 = (a + b) / 2$ and $lambda_2 = (a - b) / 2$; the Kesten--Stigum
  condition $lambda_2^2 > lambda_1$ gives $(a - b)^2 / (2 (a + b)) > 1$
  @abbe2018sbmreview.
]

#subsec("C", "Degree-Sorting Collapse of Vanilla SBM")

Because vanilla SBM assumes every vertex inside block $k$ shares identical
expected degree, strong degree heterogeneity induces catastrophic
misclassification unless degree correction @karrer2011dcsbm is introduced.

#theorem(
  title: "Exact Log-Likelihood Collapse of Vanilla SBM under Core--Periphery Degrees",
  origin: "Original Counterexample",
)[
  Consider a $2$-community graph $G$ on $n$ vertices ($4 divides n$) with zero
  inter-community edges ($E(V_1, V_2) = emptyset$, $abs(V_1) = abs(V_2) = n / 2$),
  where each community $V_c = H_c union.sq P_c$ splits equally into a dense
  core $H_c$ and sparse periphery $P_c$ ($abs(H_c) = abs(P_c) = n / 4$) with
  edge densities $p_(H H) = 1$ and $p_(H P) = p_(P P) = epsilon in (0, 1)$. For
  all sufficiently small $epsilon > 0$ and large $n$, the maximum-likelihood
  $2$-block vanilla SBM strictly prefers the *degree-sorted* partition
  $cal(P)_"deg" = {H_1 union H_2, P_1 union P_2}$ over the *true* zero-cut
  partition $cal(P)_"true" = {V_1, V_2}$ by an asymptotic log-likelihood margin
  $ ell(cal(P)_"deg") - ell(cal(P)_"true") >= binom(n, 2) (3 / 8 ln(4 / 3) - o_epsilon (1)) = Omega(n^2), $
  despite $cal(P)_"deg"$ achieving $"ARI"(cal(P)_"deg", cal(P)_"true") = 0$.
]

#proof[
  For any balanced $2$-partition $cal(P) = {C_1, C_2}$, maximizing the
  Bernoulli log-likelihood over $B_(k l)$ sets $hat(B)_(k l) = e_(k l) / N_(k l)$
  and evaluates the profile log-likelihood to
  $ ell(cal(P)) = -binom(n, 2) [ 1 / 4 h(hat(B)_(1 1)) + 1 / 4 h(hat(B)_(2 2)) + 1 / 2 h(hat(B)_(1 2)) ] + cal(O)(n), $
  where $h(x) = -x ln x - (1 - x) ln(1 - x)$ is binary entropy (with
  $h(0) = 0$) and the weights $1 / 4, 1 / 4, 1 / 2$ reflect the fractions of
  unordered pairs inside $C_1$, inside $C_2$, and across $(C_1, C_2)$.

  Evaluate the empirical block densities for both candidate partitions (up to
  $cal(O)(1 / n)$ diagonal corrections):
  - For $cal(P)_"true" = {V_1, V_2}$: inside $V_c$, one quarter of pairs are
    $H_c$--$H_c$ (density $1$) and three quarters involve $P_c$ (density
    $epsilon$), giving $hat(B)_(1 1)^"true" = hat(B)_(2 2)^"true" = 1 / 4 + (3 / 4) epsilon$,
    while $hat(B)_(1 2)^"true" = 0$. Thus the normalized negative
    log-likelihood is
    $ J(cal(P)_"true") = 1 / 2 h(1 / 4 + 3 / 4 epsilon) ->_(epsilon -> 0^+) 1 / 2 h(1 / 4). $
  - For $cal(P)_"deg" = {H_1 union H_2, P_1 union P_2}$: inside $H_1 union H_2$,
    half of pairs are intra-$H_c$ (density $1$) and half cross $H_1$--$H_2$
    (density $0$), giving $hat(B)_(1 1)^"deg" = 1 / 2$. Similarly
    $hat(B)_(2 2)^"deg" = epsilon / 2$ and $hat(B)_(1 2)^"deg" = epsilon / 2$,
    yielding
    $ J(cal(P)_"deg") = 1 / 4 h(1 / 2) + 3 / 4 h(epsilon / 2) ->_(epsilon -> 0^+) 1 / 4 h(1 / 2). $

  Compute the exact analytical entropy gap at $epsilon -> 0^+$:
  $ 1 / 2 h(1 / 4) - 1 / 4 h(1 / 2) &= 1 / 2 (2 ln 2 - 3 / 4 ln 3) - 1 / 4 ln 2 \
                                    &= 3 / 4 ln 2 - 3 / 8 ln 3 = 3 / 8 ln(4 / 3) approx 0.1079 > 0. $
  By continuity of $h$, for all sufficiently small $epsilon > 0$ (e.g.,
  $epsilon = 0.04$, where each $V_c$ is connected w.h.p. yet
  $J(cal(P)_"true") - J(cal(P)_"deg") > 0.05$), $ell(cal(P)_"deg")$ strictly
  exceeds $ell(cal(P)_"true")$ by $Omega(n^2)$, while each block of
  $cal(P)_"deg"$ splits $V_1$ and $V_2$ evenly, giving $"ARI" = 0$.
]

Replacing $B_(k l)$ with $theta_i theta_j omega_(k l)$ in the Degree-Corrected
SBM (DC-SBM) @karrer2011dcsbm absorbs the core--periphery degree disparity into
the vertex propensities $theta_i$, restoring $cal(P)_"true"$ as the unique MLE.

#subsec("D", "Davis--Kahan Spectral Consistency of SBM Embedding")

Because the VEM ELBO is non-convex jointly in $(tau, pi, B)$, theoretical
guarantees rely on spectral initialization @vonluxburg2007tutorial
@abbe2018sbmreview.

#theorem(
  title: "Davis--Kahan Eigenspace Concentration for Balanced SBM",
  origin: "Original Unified Proof",
)[
  Let $tilde(P) = Theta B Theta^top$ be the population adjacency matrix (with
  diagonal included) of a balanced $K$-block SBM, where
  $Theta in {0, 1}^(n times K)$ satisfies $Theta^top Theta = (n / K) I_K$ and
  $B$ has smallest non-zero eigenvalue magnitude $mu_"min" = min_r abs(lambda_r (B)) > 0$.
  Let $U, hat(U) in bb(R)^(n times K)$ contain orthonormal eigenvectors for the
  $K$ largest-magnitude eigenvalues of $tilde(P)$ and $A$, respectively. If
  $norm(A - tilde(P))_2 <= (n mu_"min") / (2 K)$, then:
  + Rows of $U$ satisfy $U_(i:) = U_(j:)$ iff $Z_i = Z_j$, and
    $norm(U_(i:) - U_(j:))_2 = sqrt(2 K / n)$ whenever $Z_i != Z_j$.
  + There exists an orthogonal matrix $O in cal(O)(K)$ such that
    $ norm(hat(U) O - U)_F <= (2 sqrt(2) K^(3 / 2) norm(A - tilde(P))_2) / (n mu_"min"), $
    and the fraction of vertices whose rotated row $hat(U)_(i:) O$ lies closer
    to a wrong block centroid than to $U_(i:)$ is at most
    $(16 K^2 norm(A - tilde(P))_2^2) / (n^2 mu_"min"^2)$.
]

#proof[
  *(1)* Set $W = sqrt(K / n) Theta$, which has orthonormal columns
  ($W^top W = I_K$), so $tilde(P) = (n / K) W B W^top$. Diagonalizing symmetric
  $B = Q Lambda Q^top$ ($Q in cal(O)(K)$) yields
  $tilde(P) = (W Q) ((n / K) Lambda) (W Q)^top$. Hence $U = W Q$ spans the
  non-zero eigenspace of $tilde(P)$ with smallest non-zero eigenvalue magnitude
  $delta = (n / K) mu_"min"$. Since $Q$ is orthogonal,
  $ norm(U_(i:) - U_(j:))_2 = norm((W_(i:) - W_(j:)) Q)_2 = norm(W_(i:) - W_(j:))_2, $
  which equals $0$ when $Z_i = Z_j$ and $sqrt(K / n + K / n) = sqrt(2 K / n)$
  when $Z_i != Z_j$.

  *(2)* The non-zero eigenvalues of $tilde(P)$ are separated from $0$ by gap
  $delta = (n / K) mu_"min"$. Under the Weyl margin
  $norm(A - tilde(P))_2 <= delta / 2$, the $K$ perturbed eigenvalues of $A$
  remain at distance at least $delta - norm(A - tilde(P))_2 >= delta / 2$ from
  $0$. Applying the Davis--Kahan $sin Theta$ theorem @davis1970rotation gives
  $ norm(sin Theta(hat(U), U))_F <= (norm((A - tilde(P)) U)_F) / (delta / 2) <= (2 sqrt(K) norm(A - tilde(P))_2) / delta. $
  Using the standard subspace alignment inequality
  $inf_(O in cal(O)(K)) norm(hat(U) O - U)_F <= sqrt(2) norm(sin Theta(hat(U), U))_F$
  and substituting $delta = (n / K) mu_"min"$ gives
  $ norm(hat(U) O - U)_F <= (2 sqrt(2) K^(3 / 2) norm(A - tilde(P))_2) / (n mu_"min"). $
  Finally, a vertex $i$ can be closer to a wrong centroid $U_(j:)$ ($Z_j != Z_i$)
  only if $norm(hat(U)_(i:) O - U_(i:))_2 >= (1 / 2) sqrt(2 K / n) = sqrt(K / (2 n))$.
  By Markov's inequality on squared row norms, the number of such vertices $M_"err"$
  satisfies
  $ M_"err" <= (norm(hat(U) O - U)_F^2) / (K / (2 n)) <= (16 K^2 norm(A - tilde(P))_2^2) / (n mu_"min"^2), $
  yielding error fraction $M_"err" / n <= (16 K^2 norm(A - tilde(P))_2^2) / (n^2 mu_"min"^2)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

A naive mean-field VEM implementation sums over all $binom(n, 2)$ vertex pairs
per iteration, costing $cal(O)(T n^2 K^2)$ time and becoming prohibitive beyond
a few thousand vertices. Exploiting sparsity via the global responsibility
totals $N_l = sum_(j=1)^n tau_(j l)$ (Step 2 of the algorithm) separates the
edge term from the non-edge term: computing $S = A tau$ requires $2 m K$ flops
(sparse matrix--dense matrix multiply), and evaluating the logit and non-edge
combinations $S W^top + (bold(1) N^top - tau) U^top$ requires $2 n K^2$ flops.
The M-step evaluates $E = tau^top S$ and $N_(k l)$ in $cal(O)(n K^2)$ flops.
Thus sparse VEM runs in exact $cal(O)(T (m K + n K^2))$ time per initialization
and $cal(O)(n K + m + K^2)$ memory (excluding spectral initialization, which
costs $cal(O)(T_"eig" m K + n K^2)$ via Lanczos / power iteration).

Geometrically, SBM maps vertices onto the $(K - 1)$-simplex $Delta^(K - 1)$
via soft posterior probabilities $tau_i$ rather than Euclidean coordinates,
supporting uncertainty quantification and link prediction
$hat(P)(A_(i j) = 1) = tau_i^top hat(B) tau_j$. Selecting $K$ requires
evaluating penalized likelihood criteria such as BIC @schwarz1978bic or
Integrated Completed Likelihood (ICL), multiplying total wall-clock time across
candidate values $K in {2, dots, K_"max"}$ and multiple random or spectral
restarts.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Sparse VEM time], [$cal(O)(T (m K + n K^2))$ (dense naive: $cal(O)(T n^2 K^2)$)],
    [Space complexity], [$cal(O)(n K + m + K^2)$ (CSR + responsibility matrix)],
    [Partition type], [Soft probabilistic block posteriors $tau_i in Delta^(K-1)$],
    [Structural scope], [Assortative, disassortative, bipartite, core--periphery],
    [Determinism], [Non-convex ELBO; depends on initialization / restarts],
    [Degree sensitivity], [High in vanilla SBM (Section II-C); solved by DC-SBM],
    [Hyperparameters], [Block count $K$ (selected via BIC / ICL / MDL)],
    [Link prediction], [Native generative score $hat(p)_(i j) = tau_i^top hat(B) tau_j$],
  ),
  caption: [Operational and complexity profile of SBM (Alg. 130).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

SBM is the premier graph clustering model whenever block roles go beyond simple
assortative cliques or when rigorous statistical inference and link prediction
are required:

+ *Ecological Food Webs and Trophic Role Discovery:* Predator--prey networks
  are strongly *disassortative*: species occupying the same trophic level
  rarely eat one another ($B_(k k) approx 0$) yet share predators and prey
  ($B_(k l) >> 0$). Modularity and CPM fail completely here, whereas SBM's
  general $B$ matrix recovers exact trophic guilds.
+ *Political Polarization and Bipartite Interaction Networks @karrer2011dcsbm:*
  In legislative co-voting, blog citation, and user--product graphs, DC-SBM
  separates ideological blocs while absorbing celebrity/hub degree disparities
  via vertex propensities $theta_i$.
+ *Structural Brain Connectomics and Functional Roles:* Anatomical cortical
  networks combine densely interconnected rich-club cores with sparsely linked
  sensory peripheries. SBM's core--periphery block matrices capture both
  assortative modules and hub--periphery hierarchy simultaneously.
+ *Knowledge-Graph Link Prediction and Anomaly Detection:* Because SBM yields
  an explicit generative edge probability $hat(p)_(i j) = tau_i^top hat(B) tau_j$,
  missing edges with high $hat(p)_(i j)$ serve as predicted interactions and
  observed edges with tiny $hat(p)_(i j)$ flag anomalous or spurious links.
+ *Social Network Role Equivalence @holland1983sbm:* Originally formulated in
  mathematical sociology to operationalize *stochastic structural equivalence*,
  SBM groups actors who exhibit statistically identical interaction patterns
  toward all other groups (e.g., managers vs. engineers vs. suppliers).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four operational pitfalls merit care: (i) *heavy-tailed degree heterogeneity*
in vanilla SBM triggers Section II-C's core--periphery collapse—always use DC-SBM
@karrer2011dcsbm on scale-free graphs; (ii) *below-Kesten--Stigum sparsity*
($"SNR" <= 1$, Section II-B), where no algorithm can beat random guessing;
(iii) *computational cost* of iterative VEM plus BIC/ICL sweeps over $K$ on
very large graphs; and (iv) *local ELBO traps* under random initialization,
requiring spectral warm starts (Section II-D).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Ecological food webs], [High: recovers disassortative trophic roles], [Select $K$ via ICL; spectral warm start],
    [Core--periphery / web], [High with DC-SBM @karrer2011dcsbm], [Vanilla SBM collapses on hubs (Sec. II-C)],
    [Connectomics / bio], [High: generative uncertainty + link prediction], [Slow for $n > 10^5$; try fast spectral @chung1997spectral],
    [Ultra-large assort.], [Moderate: VEM slower than linear sweeps], [Use Fast Unfolding (Alg. 129) @blondel2008louvain],
    [Dense clique overlap], [Moderate: hard-block prior (needs MMSBM)], [Use CPM (Alg. 127) or DEMON (Alg. 128)],
  ),
  caption: [Domain suitability and failure-mode matrix for SBM (Alg. 130).],
)

#sec("V", "Structural Relations and Empirical Behavior")

SBM provides the unifying statistical umbrella beneath both modularity
optimization (Fast Unfolding, Alg. 129 @blondel2008louvain) and spectral
clustering @vonluxburg2007tutorial @chung1997spectral @shi2000ncut.

#proposition(
  title: "Exact Reduction of Modularity to Planted-Partition DC-SBM",
  origin: "Original Derivation",
)[
  Consider the Poisson Degree-Corrected SBM @karrer2011dcsbm with fixed vertex
  propensities $theta_i = k_i$ (so $sum_i theta_i = 2 m$) and symmetric
  planted-partition affinity matrix
  $omega_(a b) = omega_"out" + (omega_"in" - omega_"out") delta(a, b)$
  satisfying $omega_"in" > omega_"out" > 0$, where
  $A_(i j) tilde.op "Pois"(k_i k_j omega_(c_i c_j) / (2 m))$. Then the DC-SBM
  log-likelihood reduces identically to resolution-parameterized modularity:
  $ log P(A bar.v c, omega) = "const" + m log(omega_"in" / omega_"out") Q_(gamma^*) (c; G), $
  with exact optimal resolution parameter
  $ gamma^* = (omega_"in" - omega_"out") / (log(omega_"in" / omega_"out")) in (omega_"out", omega_"in"). $
]

#proof[
  Up to terms independent of the partition $c$, the Poisson DC-SBM
  log-likelihood is
  $ log P(A bar.v c, omega) = "const" + 1 / 2 sum_(i, j = 1)^n [ A_(i j) log omega_(c_i c_j) - (k_i k_j) / (2 m) omega_(c_i c_j) ]. $
  Substituting
  $log omega_(c_i c_j) = log omega_"out" + log(omega_"in" / omega_"out") delta(c_i, c_j)$
  and $omega_(c_i c_j) = omega_"out" + (omega_"in" - omega_"out") delta(c_i, c_j)$,
  the $delta$-independent terms fold into $"const"$, leaving
  $ log P(A bar.v c, omega) &= "const" + 1 / 2 sum_(i, j) [ A_(i j) log(omega_"in" / omega_"out") \
                            & quad - (omega_"in" - omega_"out") (k_i k_j) / (2 m) ] delta(c_i, c_j) \
                            &= "const" + m log(omega_"in" / omega_"out") Q_(gamma^*) (c; G), $
  where $gamma^* = (omega_"in" - omega_"out") / log(omega_"in" / omega_"out")$
  lies strictly in $(omega_"out", omega_"in")$ by the logarithmic mean
  inequality. Thus Newman--Girvan modularity (and Fast Unfolding, Alg. 129) is
  a strict assortative, equal-parameter special case of DC-SBM, explaining why
  modularity suffers the resolution limit @fortunato2007resolution whenever
  $gamma = 1$ disagrees with $gamma^*$.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite).* Evaluated on
$k$-NN proximity graphs across the four ClusterOverview datasets (Silhouette
@rousseeuw1987silhouettes and ARI @hubert1985comparing):
- *Iris (`Sil = 0.445`, `ARI = 0.531`):* Achieves competitive ARI (`0.531`,
  trailing CPM / Louvain `0.567` and beating DEMON `0.498`) with moderate
  Silhouette (`0.445`). Because vanilla SBM penalizes degree variance inside
  blocks (Section II-C), boundary-density gradients within Versicolor--Virginica
  induce slight core--periphery splitting compared to purely assortative cuts.
- *Customers (`Sil = 0.356`):* Highest Silhouette of all four graph methods
  (decisively outperforming CPM `0.262`, DEMON `0.250`, and Fast Unfolding
  `0.190`), directly demonstrating SBM's ability to fit asymmetric
  core--periphery and multi-density customer tiers via its full $K times K$
  matrix $B$.
- *Multi-Shape (`Sil = 0.713`, `ARI = N/A`):* Ties Fast Unfolding (Alg. 129)
  at the top Silhouette (`0.713`), cleanly separating the disconnected or
  weakly bridged geometric manifold components into distinct blocks.
- *High-D (`Sil = 0.303`, `ARI = N/A`):* Highest Silhouette by a clear margin
  (beating CPM and Louvain `0.188` and DEMON `0.163`): global likelihood
  pooling across both observed edges and non-edges resists high-dimensional
  $k$-NN distance concentration @beyer1999nearest and local triangle loss.
- *Runtime (`Mean = 130.533 ms`):* Substantially the slowest of the four
  algorithms ($191.7 times$ slower than Fast Unfolding `0.681 ms`, $46.6 times$
  slower than CPM `2.800 ms`, and $29.5 times$ slower than DEMON `4.424 ms`),
  reflecting the iterative $cal(O)(T (m K + n K^2))$ VEM optimization and
  spectral initialization overhead traded for superior probabilistic modeling
  on heterogeneous and high-dimensional graphs.
