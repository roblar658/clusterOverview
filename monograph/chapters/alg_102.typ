#import "../template/ieee.typ": *

#alg_entry(
  id: 102,
  name: "SNOB (Minimum Message Length Mixture Clustering)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1968",
  paradigm: "Minimum Message Length Mixture Modeling",
  complexity: $O(T N K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Introduced by Wallace and Boulton in 1968 @wallace1968snob, SNOB is the foundational statistical mixture clustering algorithm governed by the _Minimum Message Length_ (MML) principle of inductive inference. Predating Schwarz's Bayesian Information Criterion (BIC) @schwarz1978bic by a decade and AutoClass @cheeseman1996autoclass by nearly three decades, SNOB simultaneously solves parameter estimation and automatic model-order selection ($K$) by framing clustering as lossless two-part data compression over a noiseless channel.

Let $cal(X) = \{bold(x)_1, dots, bold(x)_N\}$ be $N$ observations over $d$ attributes, where $bold(x)_i = (x_(i,1), dots, x_(i,d))$ may combine continuous measurements (recorded to finite precision quantum $epsilon_j > 0$) and discrete categorical attributes with $S_j$ states. SNOB models $cal(X)$ as an i.i.d. sample from a $K$-class finite mixture @mclachlan2000finite satisfying _intra-class attribute independence_ (generalizing Latent Class Analysis @lazarsfeld1968latent):
$ P(bold(x)_i | bold(Theta)) = sum_(k=1)^K pi_k product_(j=1)^d P_j (x_(i,j) | bold(theta)_(k,j)), $
where $bold(pi) = (pi_1, dots, pi_K) in Delta_(K-1)$ are class proportions, $bold(theta)_(k,j) in bb(R)^(nu_j)$ parameterizes the $j$-th attribute in class $k$ (e.g., Gaussian $(mu_(k,j), sigma_(k,j))$ with $nu_j = 2$ and discrete probability $P_j (x | mu, sigma) approx epsilon_j cal(N)(x; mu, sigma^2)$, or categorical $bold(p)_(k,j) in Delta_(S_j - 1)$ with $nu_j = S_j - 1$), and $bold(Theta) = (bold(pi), \{bold(theta)_(k,j)\})$ has total free parameter dimension $D(K) := (K - 1) + K sum_(j=1)^d nu_j$.

#definition(title: "Two-Part MML Objective for SNOB")[
  Given a prior $h(K)$ on cluster count, a smooth parameter prior $h(bold(Theta))$ on $bold(Theta)$, the expected Fisher information matrix $F(bold(Theta)) := bb(E)[nabla^2 (-log P(cal(X) | bold(Theta)))]$, and the optimal $D$-dimensional lattice quantizer constant $kappa_D$ ($kappa_1 = 1 / 12$, $lim_(D -> oo) 2 pi e kappa_D = 1$), the _Wallace-Boulton / MML87 two-part message length_ (in nats) is
  $ I_"MML" (K, bold(Theta); cal(X)) &:= underbrace(-log h(K) - log(K!) - log h(bold(Theta)) \
    &quad + 1 / 2 log det F(bold(Theta)) + D / 2 log kappa_D + D / 2, "Part 1: Hypothesis Code Length " I(cal(H))) \
    &quad + underbrace(-sum_(i=1)^N log P(bold(x)_i | bold(Theta)), "Part 2: Data Code Length " I(cal(X) | cal(H))). $
  Here $-log(K!)$ discounts the $K!$ equivalent label permutations of the mixture components, and SNOB seeks $(hat(K), hat(bold(Theta))) = arg min_(K, bold(Theta)) I_"MML" (K, bold(Theta); cal(X))$.
]

In the factorized intra-class model with effective component counts $N_k := sum_(i=1)^N r_(i k)$, the complete-data Fisher determinant separates cleanly across mixing weights and per-class attributes:
$ det F(bold(Theta)) approx (N^(K-1) / (product_(k=1)^K pi_k)) product_(k=1)^K product_(j=1)^d det F_(k,j) (bold(theta)_(k,j)), $
where for a Gaussian attribute $det F_(k,j) (mu_(k,j), sigma_(k,j)) = 2 N_k^2 / sigma_(k,j)^4$, and for an $S_j$-state categorical attribute $det F_(k,j) (bold(p)_(k,j)) = N_k^(S_j - 1) / product_(s=1)^(S_j) p_(k,j,s)$.

#algbox(
  title: "SNOB (MML Mixture Estimation and Structural Search)",
  inputs: [Data $cal(X) = \{bold(x)_i\} _ (i=1)^N$, attribute types and measurement quanta $\{epsilon_j\}$, priors $h(K), h(bold(Theta))$, max sweeps $T$.],
  outputs: [Optimal class count $hat(K)$, parameters $hat(bold(Theta))$, soft posteriors $\{r_(i k)\}$, hard assignments $\{c_i\}$.],
)[
  + *Initialization:* Start with $K = 1$ (or seed $K_0$ classes) and compute initial parameters $bold(Theta)^((0))$ and message length $I^((0)) = I_"MML" (K, bold(Theta)^((0)); cal(X))$.
  + *Repeat* for $t = 1, dots, T$ until no move decreases $I_"MML"$:
    - *E-Step (Partial Assignment):* For each $i in [N]$, $k in [K]$, compute
      $ r_(i k) = (pi_k product_(j=1)^d P_j (x_(i,j) | bold(theta)_(k,j))) / (sum_(l=1)^K pi_l product_(j=1)^d P_j (x_(i,j) | bold(theta)_(l,j))), quad N_k = sum_(i=1)^N r_(i k). $
    - *MML M-Step:* Update mixing weights $hat(pi)_k = (N_k + 1 / 2) / (N + K / 2)$; for Gaussian attributes set $hat(mu)_(k,j) = N_k^(-1) sum_i r_(i k) x_(i,j)$ and $hat(sigma)_(k,j)^2 = (N_k - 1)^(-1) sum_i r_(i k) (x_(i,j) - hat(mu)_(k,j))^2$ (with $N_k > 1$); for categorical attributes set $hat(p)_(k,j,s) = (N_(k,j,s) + 1 / 2) / (N_k + S_j / 2)$.
    - *Structural Search:* Propose (a) splitting a class $k$ along its highest-variance attribute into two sub-classes, (b) merging the closest pair $(k, l)$, and (c) reallocating low-probability classes; run inner EM and accept any move with $Delta I_"MML" < 0$.
  + *Output:* Return $hat(K)$, $hat(bold(Theta))$, $\{r_(i k)\}$, and $c_i = arg max_k r_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Optimal Fisher-Metric Quantization and Two-Part MML")

Why does SNOB penalize model complexity via $1 / 2 log det F(bold(Theta)) + (D / 2)(log kappa_D + 1)$ rather than merely counting parameters? We derive the exact MML87 formula from rate-distortion geometry in Fisher metric.

#theorem(title: "Exact Optimal Quantization Volume in Fisher Metric", origin: "Original Derivation")[
  Let $L(bold(Theta)) := -log P(cal(X) | bold(Theta))$ be twice continuously differentiable near a local minimizer $hat(bold(Theta)) in bb(R)^D$ with positive-definite Hessian $H approx F(hat(bold(Theta))) succ 0$, and let $h(bold(Theta))$ be locally constant across a quantization cell $cal(Q) subset bb(R)^D$ of Lebesgue volume $V > 0$ centered at $hat(bold(Theta))$. If $cal(Q)$ is an optimal lattice quantizer isotropized in the Fisher metric (normalized second moment $kappa_D$), then the expected two-part message length
  $ overline(I)(V) := -log(h(hat(bold(Theta))) V) + bb(E)_(bold(delta) in cal(Q))[L(hat(bold(Theta)) + bold(delta))] $
  is uniquely minimized at the Fisher-metric volume
  $ V^* (hat(bold(Theta))) = (kappa_D^(D / 2) sqrt(det F(hat(bold(Theta)))))^(-1), $
  yielding the exact minimum expected message length
  $ overline(I)(V^*) = -log h(hat(bold(Theta))) + 1 / 2 log det F(hat(bold(Theta))) + D / 2 log kappa_D + D / 2 + L(hat(bold(Theta))). $
]

#proof[
  Perform a second-order Taylor expansion of $L$ around the stationary point $hat(bold(Theta))$ (where $nabla L(hat(bold(Theta))) = bold(0)$):
  $ L(hat(bold(Theta)) + bold(delta)) approx L(hat(bold(Theta))) + 1 / 2 bold(delta)^top F(hat(bold(Theta))) bold(delta). $
  Apply the whitening transformation $bold(u) := F(hat(bold(Theta)))^(1 / 2) bold(delta)$, which maps $cal(Q)$ to a cell $tilde(cal(Q))$ of volume $tilde(V) := V sqrt(det F(hat(bold(Theta))))$. By definition of the optimal $D$-dimensional lattice quantizer constant $kappa_D$, a uniform random vector $bold(u)$ over $tilde(cal(Q))$ satisfies $bb(E)[||bold(u)||_2^2] = D kappa_D tilde(V)^(2 / D)$. Therefore,
  $ bb(E)_(bold(delta) in cal(Q))[1 / 2 bold(delta)^top F(hat(bold(Theta))) bold(delta)] = D / 2 kappa_D (V sqrt(det F(hat(bold(Theta)))))^(2 / D). $
  Substituting $tilde(V) = V sqrt(det F(hat(bold(Theta)))) > 0$ into $overline(I)(V)$ gives
  $ g(tilde(V)) &= -log h(hat(bold(Theta))) + 1 / 2 log det F(hat(bold(Theta))) - log tilde(V) \
    &quad + D / 2 kappa_D tilde(V)^(2 / D) + L(hat(bold(Theta))). $
  Differentiating with respect to $tilde(V)$ yields
  $ g'(tilde(V)) = -1 / tilde(V) + kappa_D tilde(V)^(2 / D - 1) = 0 ==> tilde(V)^* = kappa_D^(-D / 2). $
  At $tilde(V)^*$, the second derivative is $g''(tilde(V)^*) = 2 / (D (tilde(V)^*)^2) > 0$, confirming a strict global minimum. Finally, $V^* = tilde(V)^* / sqrt(det F(hat(bold(Theta))))$, and plugging $tilde(V)^* = kappa_D^(-D / 2)$ back into $g$ gives $-log tilde(V)^* = (D / 2) log kappa_D$ and $(D / 2) kappa_D (tilde(V)^*)^(2 / D) = D / 2$, completing the derivation.
]

#subsec("B", "Bessel-Corrected Variance and Singularity Immunity")

Unconstrained Gaussian mixture MLE suffers from catastrophic likelihood singularities whenever a component collapses onto a single point ($sigma -> 0^+$) @mclachlan2000finite. We prove that SNOB's MML criterion simultaneously yields Bessel's unbiased variance denominator $N_k - 1$ and eliminates variance-collapse singularities.

#lemma(title: "MML Variance Estimator and Singularity Immunity", origin: "Original Proof")[
  Consider a Gaussian attribute $j$ in class $k$ with effective weight $N_k = sum_(i=1)^N r_(i k) > 0$, measurement quantum $epsilon_j > 0$, Jeffreys scale prior $h(mu, sigma) prop sigma^(-1)$, and sufficient statistic $S_(k,j) := sum_(i=1)^N r_(i k) (x_(i,j) - mu)^2$.
  + In the fine-quantization regime ($sigma >> epsilon_j$) with $N_k > 1$, the MML message length in $(mu, sigma)$ is uniquely minimized by
    $ hat(mu)_(k,j) = 1 / N_k sum_(i=1)^N r_(i k) x_(i,j), quad hat(sigma)_(k,j)^2 = S_(k,j) / (N_k - 1). $
  + Under the exact discrete-quantum likelihood $P_j (x | mu, sigma) := integral_(x - epsilon_j / 2)^(x + epsilon_j / 2) cal(N)(t; mu, sigma^2) dif t$, for any $N_k > 0$ (including singleton clusters $N_k = 1$ with $S_(k,j) = 0$),
    $ lim_(sigma -> 0^+) I_"MML" (mu, sigma) = +oo, $
    strictly ruling out variance-collapse singularities.
]

#proof[
  For $N_k$ i.i.d. observations from $cal(N)(mu, sigma^2)$, the negative log-likelihood is $L(mu, sigma) = N_k log(sqrt(2 pi) sigma / epsilon_j) + S_(k,j) / (2 sigma^2)$. Taking expected second derivatives gives
  $ bb(E)[(partial^2 L) / (partial mu^2)] = N_k / sigma^2, quad bb(E)[(partial^2 L) / (partial sigma^2)] = (2 N_k) / sigma^2, quad bb(E)[(partial^2 L) / (partial mu partial sigma)] = 0, $
  so $det F_(k,j) (mu, sigma) = 2 N_k^2 / sigma^4$ and $1 / 2 log det F_(k,j) = 1 / 2 log 2 + log N_k - 2 log sigma$. Combining this with $-log h(mu, sigma) = log sigma + c_0$ and $L(mu, sigma)$ produces the exact $(mu, sigma)$-dependent message length:
  $ I_(k,j) (mu, sigma) = (N_k - 1) log sigma + S_(k,j) / (2 sigma^2) + C(N_k, epsilon_j). $
  For fixed $sigma > 0$, minimizing $S_(k,j)$ over $mu$ yields the weighted mean $hat(mu)_(k,j)$. Differentiating with respect to $sigma$ gives $(N_k - 1) / sigma - S_(k,j) / sigma^3 = 0$, whose unique positive root for $N_k > 1$ is $hat(sigma)_(k,j)^2 = S_(k,j) / (N_k - 1)$, proving Claim 1.

  For Claim 2, the exact quantized probability satisfies $P_j (x | mu, sigma) <= integral_(-oo)^(oo) cal(N)(t; mu, sigma^2) dif t = 1$, hence the exact data term obeys $-sum_i r_(i k) log P_j (x_(i,j) | mu, sigma) >= 0$ for all $sigma > 0$. Meanwhile, the hypothesis coding cost retains $-log h(mu, sigma) + 1 / 2 log det F_(k,j) (mu, sigma) = -log sigma + c_1$, which diverges to $+oo$ as $sigma -> 0^+$. Thus $I_"MML" (mu, sigma) >= -log sigma + c_1 -> +oo$ as $sigma -> 0^+$, preventing any singleton spike.
]

#subsec("C", "Bits-Back Coding Identity and Lyapunov Monotonicity")

Wallace and Boulton @wallace1968snob originally formulated hard assignments and subsequently introduced _partial assignments_ (soft responsibilities). The exact information-theoretic link is the _bits-back_ coding identity.

#theorem(title: "Bits-Back Identity and Monotone MML Descent", origin: "Original Proof")[
  For any responsibility matrix $R = (r_(i k)) in [0,1]^(N times K)$ with $sum_(k=1)^K r_(i k) = 1$, define the variational partial-assignment message length by subtracting the bits-back entropy $H(r_(i, dot.op)) := -sum_(k=1)^K r_(i k) log r_(i k)$ from the expected complete-data code length:
  $ cal(F)(bold(Theta), R) &:= I(cal(H)) + sum_(i=1)^N sum_(k=1)^K r_(i k) (-log pi_k - log P(bold(x)_i | bold(Theta)_k)) \
    &quad - sum_(i=1)^N H(r_(i, dot.op)). $
  Then for every $(bold(Theta), R)$,
  $ cal(F)(bold(Theta), R) = I_"MML" (K, bold(Theta); cal(X)) + sum_(i=1)^N D_"KL" (r_(i, dot.op) || p(dot.op | bold(x)_i, bold(Theta))), $
  where $p(k | bold(x)_i, bold(Theta)) := pi_k P(bold(x)_i | bold(Theta)_k) / P(bold(x)_i | bold(Theta))$. Consequently, alternating E-steps ($R^((t+1)) = arg min_R cal(F)(bold(Theta)^((t)), R)$), MML M-steps, and strict-descent structural moves generates a non-increasing sequence $I_"MML"^((t+1)) <= I_"MML"^((t))$ bounded below by $0$.
]

#proof[
  Expand the per-datum terms of $cal(F)(bold(Theta), R)$:
  $ &sum_(k=1)^K r_(i k) (-log pi_k - log P(bold(x)_i | bold(Theta)_k) + log r_(i k)) \
    &= sum_(k=1)^K r_(i k) log ( r_(i k) / (pi_k P(bold(x)_i | bold(Theta)_k)) ) \
    &= sum_(k=1)^K r_(i k) log ( r_(i k) / (P(bold(x)_i | bold(Theta)) p(k | bold(x)_i, bold(Theta))) ) \
    &= -log P(bold(x)_i | bold(Theta)) + sum_(k=1)^K r_(i k) log ( r_(i k) / (p(k | bold(x)_i, bold(Theta))) ). $
  Summing over $i = 1, dots, N$ and adding $I(cal(H))$ gives the exact identity $cal(F)(bold(Theta), R) = I_"MML" (K, bold(Theta); cal(X)) + sum_(i=1)^N D_"KL" (r_(i, dot.op) || p(dot.op | bold(x)_i, bold(Theta)))$. By Gibbs' inequality, $D_"KL" >= 0$ with equality if and only if $r_(i k) = p(k | bold(x)_i, bold(Theta))$. Therefore, the E-step sets $cal(F)(bold(Theta)^((t)), R^((t+1))) = I_"MML" (bold(Theta)^((t)))$, and any M-step reducing $cal(F)$ satisfies
  $ I_"MML" (bold(Theta)^((t+1))) <= cal(F)(bold(Theta)^((t+1)), R^((t+1))) <= cal(F)(bold(Theta)^((t)), R^((t+1))) = I_"MML" (bold(Theta)^((t))). $
  Structural split/merge moves are accepted only when $Delta I_"MML" < 0$. Because exact discrete-quantum data probabilities satisfy $P(bold(x)_i | bold(Theta)) <= 1$ and $I(cal(H))$ is bounded below (Lemma 1), the monotone sequence converges to a finite limit @dempster1977em @wu1983convergence.
]

#subsec("D", "Asymptotic Equivalence to BIC and Model Selection")

#theorem(title: "Asymptotic Expansion of SNOB MML to BIC", origin: "Classical — Cited")[
  For a regular $K$-component mixture of dimension $D(K)$ with positive-definite per-datum Fisher information $overline(F)(bold(Theta)) := N^(-1) F(bold(Theta)) succ 0$ and continuous strictly positive prior $h(bold(Theta)) > 0$ in a neighborhood of the MLE $hat(bold(Theta))$, as $N -> oo$ @wallace1968snob @schwarz1978bic @fraley2002mclust:
  $ I_"MML" (K, hat(bold(Theta)); cal(X)) = -log P(cal(X) | hat(bold(Theta))) + (D(K)) / 2 log N + R_N (K) = 1 / 2 "BIC"(K) + R_N (K), $
  where the finite-sample geometric residual $R_N (K) = O_p (1)$ equals
  $ R_N (K) &= -log h(K) - log(K!) - log h(hat(bold(Theta))) \
    &quad + 1 / 2 log det overline(F)(hat(bold(Theta))) + (D(K)) / 2 (log kappa_D + 1). $
]

#proof[
  Because $F(hat(bold(Theta))) = N overline(F)(hat(bold(Theta)))$ is a $D(K) times D(K)$ matrix, determinant homogeneity gives $det F(hat(bold(Theta))) = N^(D(K)) det overline(F)(hat(bold(Theta)))$, hence
  $ 1 / 2 log det F(hat(bold(Theta))) = (D(K)) / 2 log N + 1 / 2 log det overline(F)(hat(bold(Theta))). $
  Substituting into Definition 1 isolates $-log P(cal(X) | hat(bold(Theta))) + (D(K) / 2) log N = 1 / 2 "BIC"(K)$ @schwarz1978bic plus the exact residual $R_N (K)$. By consistency of $hat(bold(Theta))$ and continuity of $h$ and $det overline(F) > 0$, $R_N (K) = O_p (1)$ as $N -> oo$.
]

Critically, in finite samples, $1 / 2 log det F$ replaces the uniform $N$ inside BIC with component-specific effective sizes $N_k = pi_k N$ (contributing $(nu / 2) log pi_k$ inside $R_N (K)$), calibrating complexity for small clusters far more accurately than generic BIC @fraley2002mclust.

#subsec("E", "Correlated-Attribute Over-Splitting Counterexample")

SNOB's classical conditional-independence model assumes axis-aligned (diagonal) covariances inside each class. We prove constructively that any nonzero within-cluster correlation forces SNOB to over-split a single Gaussian cluster as $N -> oo$.

#proposition(title: "Asymptotic Over-Splitting of a Single Correlated Gaussian", origin: "Original Counterexample")[
  Let $bold(x)_1, dots, bold(x)_N tilde.op cal(N)(bold(0), bold(Sigma)_rho)$ in $bb(R)^2$ with $bold(Sigma)_rho = mat(1, rho; rho, 1)$ and $0 < |rho| < 1$. Consider diagonal-covariance SNOB comparing the true single-component model ($K = 1$, $D(1) = 4$) against a symmetric two-component diagonal mixture ($K = 2$, $D(2) = 9$). There exists a diagonal two-component mixture $p_2$ satisfying
  $ bb(E)_(bold(x) tilde.op cal(N)(bold(0), bold(Sigma)_rho)) [ log ( p_2 (bold(x)) / p_1 (bold(x)) ) ] >= (3 rho^2) / 128 > 0, $
  where $p_1 (bold(x)) = cal(N)(bold(x); bold(0), I_2)$ is the unique $K = 1$ diagonal MLE. Consequently, $"Pr"(I_"MML" (2) < I_"MML" (1)) -> 1$ as $N -> oo$, so diagonal SNOB systematically over-splits any correlated Gaussian cluster once $N = Omega(rho^(-2) log(rho^(-1)))$.
]

#proof[
  Set $s := "sign"(rho) in \{-1, +1\}$, choose $m := sqrt(|rho| / 32)$ (so $m^2 = |rho| / 32 <= 1 / 32$), and set $sigma^2 := 1 - m^2$. Define $p_2 (bold(x)) := 1 / 2 cal(N)(bold(x); bold(mu), sigma^2 I_2) + 1 / 2 cal(N)(bold(x); -bold(mu), sigma^2 I_2)$ with $bold(mu) := (m, s m)^top$. Using $1 / 2 (e^u + e^(-u)) = cosh(u)$ with $u := m(x_1 + s x_2) / sigma^2$, exact cancellation gives
  $ log ( p_2 (bold(x)) / p_1 (bold(x)) ) &= -log(1 - m^2) + 1 / 2 (1 - 1 / (1 - m^2))(x_1^2 + x_2^2) \
    &quad - m^2 / (1 - m^2) + log cosh(u). $
  First, observe the global lower bound $log cosh(u) >= u^2 / 2 - u^4 / 12$ for all $u in bb(R)$: letting $psi(u) := log cosh(u) - u^2 / 2 + u^4 / 12$, we have $psi(0) = psi'(0) = 0$ and $psi''(u) = 1 / cosh(u)^2 - 1 + u^2 = u^2 - tanh(u)^2 >= 0$ since $|tanh(u)| <= |u|$, so $psi$ is convex with global minimum $0$.

  Under $bold(x) tilde.op cal(N)(bold(0), bold(Sigma)_rho)$, $bb(E)[x_1^2 + x_2^2] = 2$ and $Z := x_1 + s x_2 tilde.op cal(N)(0, 2(1 + |rho|))$, yielding $bb(E)[Z^2] = 2(1 + |rho|)$ and $bb(E)[Z^4] = 12(1 + |rho|)^2$. Taking expectations and applying $-log(1 - m^2) >= m^2$:
  $ bb(E)[ log ( p_2 (bold(x)) / p_1 (bold(x)) ) ] &>= m^2 - (2 m^2) / (1 - m^2) + (m^2 (1 + |rho|)) / (1 - m^2)^2 - (m^4 (1 + |rho|)^2) / (1 - m^2)^4 \
    &= (m^2 (|rho| + m^4)) / (1 - m^2)^2 - (m^4 (1 + |rho|)^2) / (1 - m^2)^4. $
  Since $(1 - m^2)^2 < 1$, $(1 + |rho|)^2 <= 4$, and $(1 - m^2)^4 >= (31 / 32)^4 > 1 / 2$, the expectation is bounded below by
  $ |rho| m^2 - 8 m^4 = |rho| m^2 - 8 (|rho| / 32) m^2 = 3 / 4 |rho| m^2 = (3 rho^2) / 128 > 0. $
  By the Weak Law of Large Numbers, the $K = 2$ negative log-likelihood improves over $K = 1$ by at least $(3 rho^2 / 128) N - o_p (N)$, whereas by Theorem 3 the MML complexity penalty increases by only $(D(2) - D(1)) / 2 log N + O_p (1) = (5 / 2) log N + O_p (1)$. For $N = Omega(rho^(-2) log(rho^(-1)))$, the linear likelihood gain overwhelms the logarithmic penalty, proving $"Pr"(I_"MML" (2) < I_"MML" (1)) -> 1$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because intra-class attributes are conditionally independent, each E-step evaluates $d$ univariate densities per class and datum in $O(N K d)$ time, and each M-step updates $K d$ scalar sufficient statistics in $O(N K d)$ time. Over $T$ EM and structural-search iterations (where split proposals evaluate candidate binary partitions on class subsets of size $N_k$ with $sum_k N_k = N$, and merge checks examine $O(K^2)$ pairs in $O(K^2 d)$ time), total time complexity is $O(T N K d)$ for $K <= N$. Crucially, factorizing covariance avoids the $O(N K d^2 + K d^3)$ Cholesky and matrix-inversion bottlenecks of full-covariance GMMs @fraley2002mclust.

Space complexity is $O(N d + N K + K d)$: $O(N d)$ stores the data, $O(N K)$ stores soft responsibilities $r_(i k)$, and $O(K d)$ stores per-class attribute parameters. Because $I_"MML"$ evaluates parameter Fisher volumes separately per attribute, SNOB seamlessly handles _mixed-type_ datasets (Gaussian, categorical, Poisson, von Mises @banerjee2005vmf) and missing values (by omitting unobserved factors $j$ from both the likelihood and effective counts $N_(k,j)$).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T N K d)$ across EM and split/merge search ($O(K^2 d)$ merge screen).],
    [Space Complexity], [$O(N d + N K + K d)$ ($O(K d)$ parameter footprint, linear in $d$).],
    [Cluster Geometry], [Axis-aligned hyper-ellipsoids (continuous) and product multinomials (discrete).],
    [Invariances], [Invariant to sample/class permutations and coordinate rescaling (with scaled $epsilon_j$); _not_ rotation-invariant.],
    [Hyperparameters], [Measurement quanta $\{epsilon_j\}$ and priors $h(K), h(bold(Theta))$; $K$ selected automatically.],
    [Regularity], [Immune to singleton variance collapse (Lemma 1); monotone MML descent (Theorem 2).],
  ),
  caption: [Operational and geometric profile of SNOB.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

SNOB's combination of mixed-type factorized likelihoods, automatic model-order selection, and singularity immunity makes it well-suited to scientific taxonomy discovery:

+ *Biological and Medical Taxonomy Discovery:* Wallace and Boulton originally applied SNOB to seal morphology and clinical records @wallace1968snob. Patient cohorts routinely mix continuous biomarkers, binary symptoms, and multi-state genomic markers with missing entries; SNOB's factorized likelihood handles mixed types natively while MML prevents spurious subtyping on small cohorts.
+ *Astronomical Source Catalog Classification:* Multi-band sky surveys measure continuous fluxes, colors, and discrete morphology flags across thousands of objects. Because stellar, galactic, and quasar populations have unknown multiplicity $K$ and widely varying class sizes $N_k$, MML's effective-count penalty $log N_k$ identifies rare classes without collapsing onto individual outliers (Lemma 1).
+ *Customer Segmentation and Survey Profiling:* Marketing databases combine continuous spending totals, purchase counts, and categorical demographics. SNOB's axis-aligned classes produce directly interpretable per-attribute profiles, while the $-log(K!)$ symmetry discount avoids redundant micro-segments.
+ *Ecological and Soil-Survey Zonification:* Environmental transects record heterogeneous soil chemistry, categorical substrate types, and angular slope orientations (modeled via circular von Mises factors @banerjee2005vmf). SNOB automatically infers the parsimonious number of ecological zones justified by measurement precision $epsilon_j$.
+ *High-Dimensional Gene-Expression and Proteomics Screening:* When $d$ approaches or exceeds $N_k$, full covariance matrices become singular. SNOB's $O(K d)$ diagonal parameterization remains well-conditioned and computationally tractable where full-covariance MCLUST @fraley2002mclust cannot invert $d times d$ matrices.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SNOB's primary failure modes stem directly from its inductive assumptions: (i) *Strong intra-cluster attribute correlations* induce systematic over-splitting into chains of axis-aligned blobs (Proposition 1); (ii) *Heavy-tailed outliers* distort Gaussian means and variances unless replaced by Student-$t$ mixtures @peel2000robust; (iii) *Non-convex manifolds* (interlocking moons, rings) are approximated by multiple ellipsoids; and (iv) *Mis-specified quantization accuracy $epsilon_j$* biases the balance between data precision and class count $K$ when $epsilon_j$ approaches within-cluster scale $sigma_(k,j)$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Mixed-Type Clinical / Survey], [Optimal: native continuous + categorical + missing data.], [Verify measurement quanta $epsilon_j$ scaling.],
    [Correlated Continuous ($d <= 20$)], [Biased: diagonal assumption over-splits (Prop. 1).], [Use full-covariance MCLUST @fraley2002mclust or PCA.],
    [High-D Independent ($d >> 20$)], [Strong: $O(K d)$ parameters avoid singular $bold(Sigma)_k$.], [Check feature redundancy; compare IGMM @rasmussen2000igmm.],
    [Heavy-Tailed / Noisy Data], [Sensitive: Gaussian quadratic tails inflate $S_(k,j)$.], [Use Student-$t$ mixtures @peel2000robust.],
    [Non-Convex Manifolds], [Approximate: splits curved structures into sub-blobs.], [Use density-based or spectral clustering.],
  ),
  caption: [Domain suitability and failure-mode matrix for SNOB.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SNOB sits at the historical root of probabilistic model-based clustering, unifying Latent Class Analysis @lazarsfeld1968latent, Bayesian AutoClass @cheeseman1996autoclass, conceptual clustering (COBWEB @fisher1987cobweb, CLASSIT @gennari1989classit), and BIC-guided mixtures (X-means @pelleg2000xmeans, MCLUST @fraley2002mclust, MIXMOD @biernacki2006mixmod).

#proposition(title: "Unified Structural Bridge to LCA, AutoClass, and MCLUST", origin: "Original Derivation")[
  Under the $K$-class factorized mixture model with parameter dimension $D = D(K)$:
  + *Latent Class Analysis:* Restricting all attributes to discrete categorical variables and fixing $K$ reduces SNOB's Part 2 likelihood identically to Lazarsfeld's Latent Class Analysis @lazarsfeld1968latent.
  + *AutoClass Marginal Evidence:* Let $I_"Laplace" (K) := -log integral P(cal(X) | bold(Theta)) h(bold(Theta)) dif bold(Theta)$ denote the Laplace-approximated negative log marginal likelihood optimized by AutoClass @cheeseman1996autoclass (including $-log(h(K) K!)$). Then the exact discrepancy between SNOB's MML criterion and AutoClass's Laplace criterion is
    $ I_"MML" (K, hat(bold(Theta))) - I_"Laplace" (K) = D / 2 log(2 pi e kappa_D) + O(N^(-1)), $
    whose per-parameter gap $1 / 2 log(2 pi e kappa_D)$ vanishes monotonically as $D -> oo$ (and equals at most $0.176$ nats at $D = 1$).
  + *MCLUST Reduction:* Relaxing diagonal independence to general eigen-decomposed covariances $bold(Sigma)_k$ and dropping the $O_p (1)$ Fisher residual $R_N (K)$ of Theorem 3 yields MCLUST's BIC criterion @fraley2002mclust.
]

#proof[
  Claim 1 is immediate since $P(bold(x)_i | bold(Theta)) = sum_k pi_k product_j p_(k,j, x_(i,j))$ with fixed $K$ removes Part 1 model-order comparison. For Claim 2, expanding $L(bold(Theta)) = -log P(cal(X) | bold(Theta))$ to second order around $hat(bold(Theta))$ inside the marginal likelihood integral gives the Gaussian integral
  $ integral_(bb(R)^D) h(bold(Theta)) e^(-L(bold(Theta))) dif bold(Theta) approx h(hat(bold(Theta))) e^(-L(hat(bold(Theta)))) (2 pi)^(D / 2) (det F(hat(bold(Theta))))^(-1 / 2). $
  Taking negative logarithms and adding $-log(h(K) K!)$ yields
  $ I_"Laplace" (K) &= -log h(K) - log(K!) - log h(hat(bold(Theta))) \
    &quad + 1 / 2 log det F(hat(bold(Theta))) - D / 2 log(2 pi) + L(hat(bold(Theta))) + O(N^(-1)). $
  Subtracting this from Definition 1 cancels all data, prior, permutation, and Fisher-determinant terms, leaving
  $ I_"MML" - I_"Laplace" = D / 2 (log kappa_D + 1 + log(2 pi)) = D / 2 log(2 pi e kappa_D). $
  Because optimal lattice quantizer constants satisfy $1 / (2 pi e) <= kappa_D <= 1 / 12$ with $lim_(D -> oo) 2 pi e kappa_D = 1$, the per-parameter gap obeys $0 <= 1 / 2 log(2 pi e kappa_D) <= 1 / 2 log(2 pi e / 12) approx 0.176$ nats and vanishes as $D -> oo$. Claim 3 follows directly from Theorem 3 by substituting full/structured covariance parameterizations into $D(K)$.
]

*Empirical Behavior on Tabular Benchmarks.* On the ClusterOverview 4-dataset benchmark (measured via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), SNOB achieves:
- *Iris:* $"Sil" = 0.306$, $"ARI" = 0.456$.
- *Customers:* $"Sil" = 0.397$.
- *Multi-Shape:* $"Sil" = 0.750$, $"ARI" = "N/A"$.
- *High-D:* $"Sil" = 0.321$, $"ARI" = "N/A"$.
- *Mean Runtime:* $30.551$ ms.

These empirical measurements closely validate our theoretical profile. On Multi-Shape (`Sil = 0.750`) and Customers (`Sil = 0.397`), MML's automatic model-order selection and singularity immunity (Lemma 1) isolate well-separated structures with high cluster cohesion. On Iris (`Sil = 0.306`, `ARI = 0.456`), strong positive correlations among petal length, petal width, and sepal length trigger the diagonal-covariance over-splitting bias established in Proposition 1, reducing label recovery relative to full-covariance Gaussian mixtures @fraley2002mclust. Conversely, on High-D (`Sil = 0.321`), SNOB's $O(K d)$ conditional-independence factorization avoids high-dimensional covariance ill-conditioning and markedly outperforms BTM (`Sil = 0.162`). The $30.551$ ms mean runtime reflects the modest overhead of iterative MML split/merge evaluations alongside EM updates.
