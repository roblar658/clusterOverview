#import "../template/ieee.typ": *

#alg_entry(
  id: 100,
  name: "Correlated Topic Model (CTM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2007",
  paradigm: "Logistic-Normal Variational Topic Model",
  complexity: $cal(O)(T n (K bar(L) + K^2) + T K^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Latent Dirichlet Allocation (LDA) @blei2003lda represents each document as a
probability mixture over $K$ latent topics drawn from a Dirichlet prior
$bold(theta)_d ~ "Dir"_K (bold(alpha))$. Despite its analytical conjugacy, the
Dirichlet distribution enforces a severe structural restriction: apart from the
simplex sum-to-one constraint, components of $bold(theta)_d$ are nearly
independent and strictly negatively correlated (Theorem 1), making LDA unable
to express the pervasive empirical fact that certain themes—such as genetics
and evolutionary biology, or monetary policy and inflation—co-occur positively
within documents. The *Correlated Topic Model (CTM)* of Blei and Lafferty
@blei2007ctm resolves this limitation by replacing the Dirichlet prior with a
*logistic-normal* distribution on the simplex, modeling arbitrary pairwise
topic correlations through a full $K times K$ covariance matrix $bold(Sigma)$.

Let $cal(D) = {d_1, ..., d_n}$ be a corpus over vocabulary $cal(V) = {1, ..., V}$,
where document $d$ contains $N_d$ word tokens $(w_(d,1), ..., w_(d, N_d))$
($U_d$ unique terms). Let $bold(B) in bb(R)^(K times V)$ collect the $K$
topic-word distributions $bold(beta)_k in Delta^(V-1)$, and let
$(bold(mu), bold(Sigma)) in bb(R)^K times bb(S)_(++)^K$ parameterize the latent
Gaussian prior. The CTM generative process for each document $d = 1, ..., n$ is:
$
  bold(eta)_d &~ cal(N)_K (bold(mu), bold(Sigma)), \
  theta_(d,k) &= f_k (bold(eta)_d) :=
  exp(eta_(d,k)) / (sum_(j=1)^K exp(eta_(d,j))), quad k = 1, ..., K, \
  z_(d,m) | bold(theta)_d &~ "Cat"(bold(theta)_d), quad m = 1, ..., N_d, \
  w_(d,m) | z_(d,m), bold(B) &~ "Cat"(bold(beta)_(z_(d,m))).
$
Document clustering is obtained either softly via posterior topic proportions
$hat(bold(theta))_d = f(bold(lambda)_d)$ or as hard partitions via dominant
topic assignment $hat(c)_d = arg max_k lambda_(d,k)$.

#definition(title: "Mean-Field Variational CTM Family and Tangent ELBO")[
  Because the logistic-normal prior is non-conjugate to the multinomial
  likelihood, CTM introduces a factorized per-document variational distribution
  $
    q_d (bold(eta)_d, bold(z)_d) =
    ( product_(k=1)^K cal(N)(eta_(d,k) | lambda_(d,k), nu_(d,k)^2) )
    ( product_(m=1)^(N_d) "Cat"(z_(d,m) | bold(phi)_(d,m)) ),
  $
  with variational parameters $bold(lambda)_d in bb(R)^K$,
  $bold(nu)_d^2 in bb(R)_(> 0)^K$, $bold(phi)_(d,m) in Delta^(K-1)$, plus an
  auxiliary scalar tangent parameter $zeta_d > 0$ that linearizes the intractable
  log-sum-exp expectation $bb(E)_(q_d) [log sum_(j=1)^K exp(eta_(d,j))]$.
]

Inference proceeds via Variational Expectation-Maximization (VEM) @blei2007ctm
@blei2006variational: the variational E-step tightens the lower bound
$cal(L)_d$ per document via coordinate ascent, and the M-step updates
$(bold(mu), bold(Sigma), bold(B))$ in exact closed form.

#algbox(
  title: "Variational EM for the Correlated Topic Model (CTM)",
  inputs: [Corpus $cal(D) = {d_1, ..., d_n}$, topics $K$, outer iterations $T$,
    inner tolerance $epsilon_"in"$],
  outputs: [Model parameters $(bold(mu), bold(Sigma), bold(B))$, topic
    coordinates ${bold(lambda)_d}_(d=1)^n$, hard clusters ${hat(c)_d}$],
)[
  + *Initialize:* Set $bold(mu) = bold(0)$, $bold(Sigma) = bold(I)_K$, and
    initialize $bold(B)$ via smoothed random documents or LDA warm-start.
  + *For* outer iteration $t = 1, ..., T$ *do:*
    + *Variational E-step:* Precompute $bold(Sigma)^(-1)$; for each document
      $d = 1, ..., n$, iterate until $cal(L)_d$ changes by $< epsilon_"in"$:
      + Update tangent parameter:
        $zeta_d <- sum_(k=1)^K exp(lambda_(d,k) + nu_(d,k)^2 / 2)$.
      + Update token-topic responsibilities for each unique token $v in d$:
        $phi_(d,v,k) prop beta_(k,v) exp(lambda_(d,k))$.
      + Maximize $cal(L)_d$ over $bold(lambda)_d in bb(R)^K$ via conjugate
        gradient / Newton (strictly concave; Theorem 2).
      + Maximize $cal(L)_d$ over each $nu_(d,k)^2 > 0$ via 1D Newton with
        positivity safeguard (strictly concave; Theorem 2).
    + *Variational M-step:* Update corpus parameters in closed form:
      $
        bold(mu) &<- n^(-1) sum_(d=1)^n bold(lambda)_d, \
        bold(Sigma) &<- n^(-1) sum_(d=1)^n [ "diag"(bold(nu)_d^2)
        + (bold(lambda)_d - bold(mu))(bold(lambda)_d - bold(mu))^top ], \
        beta_(k,v) &<- (sum_(d=1)^n w_(d,v) phi_(d,v,k))
        / (sum_(v'=1)^V sum_(d=1)^n w_(d,v') phi_(d,v',k)).
      $
  + *Return:* $(bold(mu), bold(Sigma), bold(B))$, posterior topic profiles
    $f(bold(lambda)_d)$, and hard clusters $hat(c)_d = arg max_k lambda_(d,k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Universal Negative Correlation of Dirichlet vs. CTM Expressiveness")

We rigorously quantify the correlation barrier of LDA and construct an explicit
positive-definite CTM prior achieving arbitrarily strong positive topic
correlation.

#theorem(
  title: "Dirichlet Negative-Correlation Barrier vs. CTM Positive Correlation",
  origin: "Original Proof",
)[
  Let $bold(theta) ~ "Dir"_K (bold(alpha))$ with $alpha_k > 0$ and
  $alpha_0 = sum_(k=1)^K alpha_k$. For every distinct pair $j != k$,
  $
    "Cov"(theta_j, theta_k) = -(alpha_j alpha_k) / (alpha_0^2 (alpha_0 + 1)) < 0,
  $
  and $"Corr"(theta_j, theta_k) = -sqrt((alpha_j alpha_k) / ((alpha_0 - alpha_j)(alpha_0 - alpha_k)))$
  is strictly negative and rigidly fixed by the topic means. Conversely, in a
  $K >= 3$ CTM, for every $rho_* in (0, 1)$ there exists a strictly
  positive-definite covariance $bold(Sigma) in bb(S)_(++)^K$ such that
  $"Corr"(theta_1, theta_2) > rho_*$.
]

#proof[
  For $bold(theta) ~ "Dir"_K (bold(alpha))$, integrating the polynomial moments
  via the multivariate Beta function gives
  $bb(E)[theta_j] = alpha_j / alpha_0$,
  $bb(E)[theta_j^2] = (alpha_j (alpha_j + 1)) / (alpha_0 (alpha_0 + 1))$, and
  $bb(E)[theta_j theta_k] = (alpha_j alpha_k) / (alpha_0 (alpha_0 + 1))$ for
  $j != k$. Subtracting $bb(E)[theta_j] bb(E)[theta_k]$ gives
  $
    "Cov"(theta_j, theta_k) = (alpha_j alpha_k) / (alpha_0 (alpha_0 + 1))
    - (alpha_j alpha_k) / alpha_0^2
    = -(alpha_j alpha_k) / (alpha_0^2 (alpha_0 + 1)) < 0.
  $
  Dividing by $sqrt("Var"(theta_j) "Var"(theta_k))$ where
  $"Var"(theta_j) = alpha_j (alpha_0 - alpha_j) / (alpha_0^2 (alpha_0 + 1))$
  yields the exact rigid correlation formula.

  For the CTM construction with $K >= 3$, set $bold(mu) = bold(0)$ and define
  $bold(eta) = Z bold(u) + bold(epsilon)$, where $Z ~ cal(N)(0, sigma^2)$ with
  $sigma > 0$, $bold(u) = (1, 1, 0, ..., 0)^top in bb(R)^K$, and
  $bold(epsilon) ~ cal(N)_K (bold(0), tau^2 bold(I)_K)$ is independent of $Z$.
  Then $bold(Sigma) = tau^2 bold(I)_K + sigma^2 bold(u) bold(u)^top in bb(S)_(++)^K$
  is strictly positive definite for every $tau > 0$. By the law of total
  covariance conditioning on $Z$,
  $
    "Cov"(theta_1, theta_2) = "Cov"(bb(E)[theta_1 | Z], bb(E)[theta_2 | Z])
    + bb(E)["Cov"(theta_1, theta_2 | Z)].
  $
  Because $(epsilon_1, epsilon_2)$ are identically distributed and enter
  $f(bold(eta))$ symmetrically, $bb(E)[theta_1 | Z] = bb(E)[theta_2 | Z] =: g_tau (Z)$
  and $"Var"(theta_1) = "Var"(theta_2) = "Var"(g_tau (Z)) + bb(E)["Var"(theta_1 | Z)]$.
  As $tau -> 0^+$, $theta_1, theta_2 -> g_0 (Z) := e^Z / (2 e^Z + K - 2)$
  pointwise almost surely. Since $|theta_k| <= 1$, the dominated convergence
  theorem gives $bb(E)["Cov"(theta_1, theta_2 | Z)] -> 0$,
  $bb(E)["Var"(theta_1 | Z)] -> 0$, and
  $"Var"(g_tau (Z)) -> "Var"(g_0 (Z)) > 0$ (as $g_0$ is strictly increasing and
  non-constant on $bb(R)$). Therefore
  $
    lim_(tau -> 0^+) "Corr"(theta_1, theta_2)
    = lim_(tau -> 0^+) ("Var"(g_tau (Z)) + o(1)) / ("Var"(g_tau (Z)) + o(1))
    = 1.
  $
  Choosing $tau > 0$ sufficiently small yields $"Corr"(theta_1, theta_2) > rho_*$
  while maintaining $bold(Sigma) succ 0$.
]

#subsec("B", "Tangent Variational ELBO and Strict Block Concavity")

#theorem(
  title: "Validity and Strict Coordinate Concavity of the CTM ELBO",
  origin: "Classical — Cited",
)[
  For every document $d$ and any $zeta_d > 0$, the evidence lower bound
  $
    cal(L)_d &= -1/2 log det(bold(Sigma))
    - 1/2 "Tr"(bold(Sigma)^(-1) ("diag"(bold(nu)_d^2) + (bold(lambda)_d - bold(mu))
    (bold(lambda)_d - bold(mu))^top)) \
    &quad + sum_(m=1)^(N_d) sum_(k=1)^K phi_(d,m,k) (lambda_(d,k) + log beta_(k, w_(d,m))
    - log phi_(d,m,k)) \
    &quad - N_d ( zeta_d^(-1) sum_(k=1)^K exp(lambda_(d,k) + nu_(d,k)^2 / 2)
    + log zeta_d - 1 )
    + 1/2 sum_(k=1)^K (log(nu_(d,k)^2) + 1)
  $
  satisfies $cal(L)_d <= log p(d | bold(mu), bold(Sigma), bold(B))$
  @blei2007ctm. Moreover, holding other blocks fixed: (i) $zeta_d^* = sum_k exp(lambda_(d,k) + nu_(d,k)^2 / 2)$
  uniquely maximizes $cal(L)_d$; (ii) $cal(L)_d$ is strictly concave in
  $bold(lambda)_d in bb(R)^K$; and (iii) $cal(L)_d$ is strictly concave in each
  variance coordinate $v_k = nu_(d,k)^2 > 0$.
]

#proof[
  Applying Jensen's inequality to $log p(d | bold(mu), bold(Sigma), bold(B))$
  with respect to $q_d$ yields the exact variational lower bound
  $bb(E)_(q_d) [log p(bold(eta)_d, bold(z)_d, d)] + H(q_d)$. Every expectation
  evaluates directly in closed form except
  $-N_d bb(E)_(q_d) [log sum_(k=1)^K exp(eta_(d,k))]$. Since $x |-> log x$ is
  concave on $(0, infinity)$, its tangent line at any $zeta_d > 0$ satisfies
  $log x <= x / zeta_d + log zeta_d - 1$ with equality iff $x = zeta_d$.
  Applying this inequality to $x = sum_k exp(eta_(d,k))$ and using the exact
  Gaussian moment generating function
  $bb(E)_(q_d) [exp(eta_(d,k))] = exp(lambda_(d,k) + nu_(d,k)^2 / 2)$ produces
  $cal(L)_d <= log p(d | bold(mu), bold(Sigma), bold(B))$.

  Differentiating $cal(L)_d$ with respect to $zeta_d$ gives
  $partial cal(L)_d / partial zeta_d = N_d (zeta_d^(-2) sum_k exp(lambda_(d,k) + nu_(d,k)^2 / 2) - zeta_d^(-1)) = 0$,
  yielding the unique positive maximizer $zeta_d^*$ (since
  $partial^2 cal(L)_d / partial zeta_d^2 |_(zeta_d^*) = -N_d / (zeta_d^*)^2 < 0$).
  For $bold(lambda)_d$, the exact Hessian is
  $
    nabla_(bold(lambda)_d)^2 cal(L)_d = -bold(Sigma)^(-1)
    - (N_d / zeta_d) "diag"(exp(lambda_(d,k) + nu_(d,k)^2 / 2)) prec 0,
  $
  which is strictly negative definite because $bold(Sigma)^(-1) succ 0$ and all
  diagonal entries are positive. Finally, writing $v_k = nu_(d,k)^2 > 0$,
  $
    (partial^2 cal(L)_d) / (partial v_k^2) =
    -1 / (2 v_k^2) - (N_d / (4 zeta_d)) exp(lambda_(d,k) + v_k / 2) < 0
  $
  strictly for all $v_k > 0$, proving strict concavity in each variance coordinate.
]

#subsec("C", "Exact Variational M-Step and Automatic Non-Singularity")

Unlike maximum-likelihood GMMs, whose sample covariance matrix becomes singular
whenever $n <= K$, we prove that CTM's variational covariance update is
*unconditionally positive definite* for every corpus size $n >= 1$.

#theorem(
  title: "Exact M-Step and Unconditional Positive Definiteness of Sigma",
  origin: "Original Proof",
)[
  At every variational EM iteration, the global maximizer of the corpus ELBO
  $cal(L) = sum_(d=1)^n cal(L)_d$ over $(bold(mu), bold(Sigma), bold(B))$ is
  given uniquely by the closed-form M-step of the algorithm box. Furthermore,
  even when $n < K$, the updated covariance matrix $bold(Sigma)^((t+1))$ is
  strictly positive definite, with minimum eigenvalue bounded below by
  $
    lambda_"min" (bold(Sigma)^((t+1))) >=
    min_(k in {1, ..., K}) ( n^(-1) sum_(d=1)^n nu_(d,k)^2 ) > 0.
  $
]

#proof[
  Isolating terms in $cal(L) = sum_(d=1)^n cal(L)_d$ that depend on
  $(bold(mu), bold(Sigma))$ gives
  $
    cal(L)_(bold(mu), bold(Sigma)) =
    -n/2 log det(bold(Sigma)) - 1/2 "Tr"(bold(Sigma)^(-1) bold(M)(bold(mu))),
  $
  where
  $bold(M)(bold(mu)) = sum_(d=1)^n [ "diag"(bold(nu)_d^2) + (bold(lambda)_d - bold(mu))(bold(lambda)_d - bold(mu))^top ]$.
  Setting $nabla_(bold(mu)) cal(L) = bold(Sigma)^(-1) sum_d (bold(lambda)_d - bold(mu)) = bold(0)$
  yields $bold(mu)^((t+1)) = n^(-1) sum_d bold(lambda)_d$. For fixed $bold(mu)$,
  standard matrix calculus over $bold(Sigma)^(-1) in bb(S)_(++)^K$ gives the
  unique maximizer $bold(Sigma)^((t+1)) = n^(-1) bold(M)(bold(mu)^((t+1)))$.
  Similarly, maximizing $sum_(d,m,k) phi_(d,m,k) log beta_(k, w_(d,m))$
  subject to $sum_v beta_(k,v) = 1$ via Lagrange multipliers yields
  $beta_(k,v)^((t+1)) prop sum_d w_(d,v) phi_(d,v,k)$.

  To prove strict positive definiteness, examine the E-step first-order
  optimality condition for $v_(d,k) = nu_(d,k)^2 > 0$:
  $
    (partial cal(L)_d) / (partial v_(d,k)) =
    -1/2 (bold(Sigma)^(-1))_(k,k) + 1 / (2 v_(d,k))
    - (N_d / (2 zeta_d)) exp(lambda_(d,k) + v_(d,k) / 2) = 0.
  $
  As $v_(d,k) -> 0^+$, the logarithmic entropy barrier $1 / (2 v_(d,k))$
  diverges to $+infinity$ while the remaining two terms remain finite, and as
  $v_(d,k) -> +infinity$, the derivative tends to $-infinity$. By Theorem 2,
  this guarantees a unique finite, strictly positive solution $nu_(d,k)^2 > 0$
  for every $d, k$. Since the outer-product sum
  $sum_d (bold(lambda)_d - bold(mu))(bold(lambda)_d - bold(mu))^top$ is
  positive semidefinite, Weyl's inequality gives
  $
    bold(Sigma)^((t+1)) succ.eq n^(-1) sum_(d=1)^n "diag"(bold(nu)_d^2) succ 0,
  $
  yielding the explicit strictly positive lower bound on $lambda_"min" (bold(Sigma)^((t+1)))$.
]

#subsec("D", "Shift-Invariance Gauge Symmetry and Identifiable Projection")

#proposition(
  title: "Additive Shift Gauge Invariance and Rank-(K-1) Identifiability",
  origin: "Original Proof",
)[
  The softmax transformation satisfies $f(bold(eta) + c bold(1)_K) = f(bold(eta))$
  for every $c in bb(R)$, inducing a non-identifiable flat shift direction in
  unconstrained $K$-dimensional CTM coordinates: replacing $bold(eta)_d$ by
  $bold(eta)_d + Z_d bold(1)_K$ with $Z_d ~ cal(N)(c_0, sigma_0^2)$ maps
  $(bold(mu), bold(Sigma)) |-> (bold(mu) + c_0 bold(1)_K, bold(Sigma) + sigma_0^2 bold(1)_K bold(1)_K^top)$
  without altering the likelihood $p(cal(D))$. Applying the additive log-ratio
  contrast matrix $bold(J) = [bold(I)_(K-1), -bold(1)_(K-1)] in bb(R)^((K-1) times K)$
  yields a bijective, strictly identifiable parameterization
  $(tilde(bold(mu)), tilde(bold(Sigma))) = (bold(J) bold(mu), bold(J) bold(Sigma) bold(J)^top)$
  on $bb(R)^(K-1) times bb(S)_(++)^(K-1)$.
]

#proof[
  Adding $c bold(1)_K$ multiplies every numerator and denominator term in
  $f_k (bold(eta))$ by $e^c$, which cancels identically; hence adding an
  independent Gaussian shift $Z_d bold(1)_K$ changes $(bold(mu), bold(Sigma))$
  by $(c_0 bold(1)_K, sigma_0^2 bold(1)_K bold(1)_K^top)$ while leaving the
  law of $bold(theta)_d = f(bold(eta)_d)$ invariant. Under the contrast map
  $tilde(bold(eta)) = bold(J) bold(eta) = (eta_1 - eta_K, ..., eta_(K-1) - eta_K)^top$,
  note that $bold(J) bold(1)_K = bold(0)$, annihilating the entire shift kernel.
  Moreover, the inverse softmax map $tilde(eta)_k = log(theta_k / theta_K)$
  ($k = 1, ..., K - 1$) is a smooth diffeomorphism between the open simplex
  $"int"(Delta^(K-1))$ and $bb(R)^(K-1)$. Since a multivariate Gaussian on
  $bb(R)^(K-1)$ is uniquely identified by its mean and covariance,
  $(tilde(bold(mu)), tilde(bold(Sigma)))$ uniquely identifies the logistic-normal
  topic distribution on $Delta^(K-1)$.
]

#subsec("E", "Monotone Variational Ascent and Precision-Matrix Topic Graphs")

#proposition(
  title: "Monotone ELBO Convergence and Conditional Topic-Graph Structure",
  origin: "Classical — Cited",
)[
  Variational EM generates a monotonically non-decreasing sequence
  $cal(L)^((0)) <= cal(L)^((1)) <= ... <= 0$ that converges to a finite limit
  @blei2007ctm @dempster1977em. Furthermore, under the Gaussian prior
  $bold(eta)_d ~ cal(N)_K (bold(mu), bold(Sigma))$ with precision matrix
  $bold(Omega) = bold(Sigma)^(-1)$, two latent topic log-odds satisfy
  conditional independence $eta_(d,j) perp eta_(d,k) | bold(eta)_(d, -{j,k})$
  if and only if $Omega_(j,k) = 0$.
]

#proof[
  By Theorems 2 and 3, each coordinate block update in the E-step and each
  closed-form update in the M-step maximizes $cal(L)$ with respect to its block
  holding all other blocks fixed, ensuring $cal(L)^((t+1)) >= cal(L)^((t))$.
  Because every document likelihood $p(d | bold(mu), bold(Sigma), bold(B))$ is
  a discrete probability in $(0, 1]$, we have
  $cal(L)^((t)) <= sum_(d=1)^n log p(d | Theta^((t))) <= 0$. A bounded
  monotone real sequence converges. For the Gaussian graphical-model property,
  the joint density of $bold(eta)_d$ factors as
  $exp(-(1/2) sum_(j,k) Omega_(j,k) (eta_(d,j) - mu_j)(eta_(d,k) - mu_k))$; by
  the Hammersley-Clifford factorization theorem, the cross-term between
  $eta_(d,j)$ and $eta_(d,k)$ vanishes iff $Omega_(j,k) = 0$, which is
  equivalent to conditional independence given all remaining $K - 2$ topics
  @bishop2006pattern @blei2007ctm.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each outer VEM iteration inverts $bold(Sigma)$ once in $cal(O)(K^3)$ time.
Inside the E-step, updating $phi_(d,v,k)$ costs $cal(O)(K U_d)$ per unique
token, updating $bold(lambda)_d$ via conjugate gradient costs $cal(O)(K^2)$ per
inner step (avoiding the $cal(O)(K^3)$ per-document dense Newton solve), and
updating $bold(nu)_d^2$ costs $cal(O)(K)$. Across $n$ documents of average
unique length $bar(L)$ and $T$ outer iterations (suppressing the small inner
iteration factor), total complexity is $cal(O)(T n (K bar(L) + K^2) + T K^3)$.
Space complexity is $cal(O)(K V + K^2 + n K + n bar(U))$. Geometrically, CTM
embeds documents into the Aitchison simplex geometry via the additive log-ratio
diffeomorphism (Proposition 1), capturing compositional ratios and cross-topic
ellipsoids rather than Euclidean distances.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T n (K bar(L) + K^2) + T K^3)$ with CG E-step; $cal(O)(T n K^3)$ with full Newton],
    [Space Complexity], [$cal(O)(K V + K^2 + n K + n bar(U))$],
    [Cluster Geometry], [Soft mixed-membership simplex cones + Gaussian log-odds ellipsoids],
    [Topic Correlation], [Arbitrary positive/negative pairwise correlations (Theorem 1)],
    [Covariance Safety], [$bold(Sigma)^((t+1)) succ 0$ unconditionally even for $n < K$ (Theorem 3)],
    [Identifiability], [Shift-invariant along $bold(1)_K$; identifiable on $bb(R)^(K-1)$ (Prop. 1)],
  ),
  caption: [Operational and complexity profile of the Correlated Topic Model.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Interdisciplinary Scientific Literature Mapping:* As originally demonstrated
  on the *Science* archive @blei2007ctm, research articles span correlated
  fields (e.g., genomics and ecology, physics and astrophysics). CTM recovers
  both document topic profiles and a sparse topic-correlation graph via
  $ell_1$-regularized estimation of $bold(Omega) = bold(Sigma)^(-1)$
  (Proposition 2).
+ *Clinical EHR Comorbidity Phenotyping:* Patient diagnostic code histories
  (ICD-10 bags of codes) exhibit strong positive correlations across comorbid
  chronic conditions (e.g., diabetes, renal failure, cardiovascular disease),
  directly contradicting LDA's negative-correlation constraint (Theorem 1).
+ *Financial 10-K Risk-Factor Correlation Analysis:* Corporate regulatory
  disclosures combine co-occurring macroeconomic, supply-chain, and liquidity
  risk narratives whose joint correlation structure informs portfolio systemic
  risk clustering.
+ *Multi-Theme News Streams and Digital Humanities:* Long-form journalism
  systematically couples geopolitical, energy, and market themes; CTM models
  positive topic co-occurrence while supporting soft multi-label retrieval.
+ *Bag-of-Visual-Words Scene Understanding:* Images represented as quantized
  local patch histograms contain strongly correlated visual topics (e.g., sky,
  ocean, sand) that CTM captures far more accurately than independent Dirichlet
  priors.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CTM exhibits four primary failure modes: (i) *low-dimensional continuous
Euclidean tabular data*, where normalizing features onto the simplex discards
Euclidean magnitude and collapses metric separation; (ii) *large-$K$
over-parameterization*, because $bold(Sigma)$ introduces $K (K + 1) / 2$
parameters that overfit when $K^2 >= n$ unless regularized via graphical
lasso @blei2007ctm; (iii) *short texts*, where sparse $N_d$ leaves
$bold(eta)_d$ dominated by the prior rather than lexical evidence (prefer DMM
@yin2014dmm or BTM @yan2013btm); and (iv) *variational posterior covariance
underestimation*, because the diagonal factorization
$product_k cal(N)(lambda_(d,k), nu_(d,k)^2)$ ignores posterior correlation
among components of $bold(eta)_d$ for individual documents.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Long multi-topic corpora], [Strong; captures positive topic graphs & mixed membership], [Regularize $bold(Sigma)^(-1)$ when $K > 50$],
    [EHR comorbidity codes], [Ideal for positively correlated disease themes], [Check document-length imbalance],
    [Short texts / tweets], [Weak; sparse $N_d$ under-determines $bold(eta)_d$], [Use DMM @yin2014dmm or BTM @yan2013btm],
    [Continuous Euclidean data], [Severe mismatch; discards scale & Euclidean balls], [Use GMM / IGMM @rasmussen2000igmm or Skew-Normal @lin2009skewnormal],
    [Very large $K$ ($> 200$)], [$cal(O)(K^3)$ inversion & $cal(O)(K^2)$ overfitting], [Use sparse precision / factor CTM or LDA @blei2003lda],
  ),
  caption: [Domain suitability matrix for the Correlated Topic Model.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CTM directly extends LDA @blei2003lda by lifting the prior on $Delta^(K-1)$
from the $K$-parameter Dirichlet family to the $(K + K (K+1)/2)$-parameter
logistic-normal family (Theorem 1), sitting alongside pLSA @hofmann1999plsa and
Non-Negative Matrix Factorization @lee1999nmf @ding2008nmfplsa in the
mixed-membership topic-factorization paradigm. Within Category IV, CTM forms an
informative structural polarity with DMM (Alg. 97 @yin2014dmm): DMM enforces
hard single-topic membership with exact collapsed conjugacy, whereas CTM allows
correlated multi-topic membership at the price of non-conjugate tangent
variational approximation (Theorem 2). Unlike IGMM (Alg. 99 @rasmussen2000igmm)
and HDP @teh2006hdp, CTM requires fixing $K$ a priori, though its covariance
matrix $bold(Sigma)$ can absorb redundant topic correlations gracefully.

On the ClusterOverview 4-dataset continuous tabular benchmark, CTM produces a
textbook diagnostic signature of inductive-bias mismatch: Iris Silhouette
$= -0.189$ and ARI $= 0.049$ @rousseeuw1987silhouettes @hubert1985comparing,
Customers Silhouette $= 0.251$, Multi-Shape Silhouette $= 0.509$ (ARI N/A),
High-D Silhouette $= 0.206$ (ARI N/A), and the fastest mean runtime of all four
algorithms at $2.256$ ms. The negative Iris Silhouette ($-0.189$) and near-zero
ARI ($0.049$) have a clean mathematical explanation: Iris comprises only $d = 4$
dense, collinear continuous measurements; casting them into a 4-word
compositional simplex via softmax discards absolute floral scale (total petal
and sepal length) and mixes distinct species that share similar feature ratios
across soft overlapping topics, causing hard $arg max_k lambda_(d,k)$ labels to
cut straight across Euclidean clusters. By contrast, Customers ($0.251$) and
Multi-Shape ($0.509$) contain partially separable positive compositional cones,
allowing modest clustering recovery, while the small feature vocabulary $V$ and
low topic count $K$ let deterministic variational EM converge in just $2.256$ ms.
