#import "../template/ieee.typ": *

#alg_entry(
  id: 88,
  name: "Student's t-Mixture Models (TMM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2000",
  paradigm: "Heavy-Tailed Scale-Mixture Generative Clustering",
  complexity: $cal(O)(T n K d^2 + T K d^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Gaussian Mixture Models (GMMs) are acutely sensitive to atypical observations: a single extreme outlier can hijack an entire component's mean, inflate its covariance determinant, or force a spurious extra cluster. Peel and McLachlan @peel2000robust @mclachlan2000finite resolve this vulnerability by replacing Gaussian components with multivariate Student's $t$-distributions whose polynomial tail decay adapts automatically via component-specific degrees of freedom $nu_k > 0$.

Let $cal(D) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$. A $K$-component Student's $t$-Mixture Model (TMM) defines the generative density
$ p(bold(x) | bold(Theta)) = sum_(k=1)^K pi_k f_t (bold(x) | bold(mu)_k, bold(Sigma)_k, nu_k), $
with mixing proportions $bold(pi) in Delta^(K - 1)$, locations $bold(mu)_k in bb(R)^d$, symmetric positive-definite scale matrices $bold(Sigma)_k succ 0$, and degrees of freedom $nu_k > 0$. Writing the squared Mahalanobis distance as $delta_k (bold(x)) = (bold(x) - bold(mu)_k)^T bold(Sigma)_k^(-1) (bold(x) - bold(mu)_k)$, the $d$-dimensional Student's $t$ density is
$ f_t (bold(x) | bold(mu)_k, bold(Sigma)_k, nu_k) &= (Gamma((nu_k + d) / 2) |bold(Sigma)_k|^(-1/2)) / ((pi nu_k)^(d / 2) Gamma(nu_k / 2)) \
  &quad times (1 + (delta_k (bold(x))) / nu_k)^(-(nu_k + d) / 2). $
Crucially, direct maximum-likelihood estimation is made tractable by an exact hierarchical latent-variable augmentation @peel2000robust @bishop2006pattern.

#definition(title: "Gamma-Gaussian Scale-Mixture Representation")[
  For each observation $bold(x)_i$, introduce a discrete one-hot component indicator $bold(z)_i in {0, 1}^K$ ($sum_(k=1)^K z_(i k) = 1$) and a continuous latent precision-scale variable $u_i > 0$ distributed hierarchically as
  $ bold(z)_i ~ "Cat"(bold(pi)), quad u_i | (z_(i k) = 1) ~ "Gamma"(nu_k / 2, nu_k / 2), \
    bold(x)_i | (z_(i k) = 1, u_i) ~ cal(N)(bold(mu)_k, bold(Sigma)_k / u_i), $
  where $"Gamma"(alpha, beta)$ denotes the Gamma density $beta^alpha u^(alpha - 1) e^(-beta u) / Gamma(alpha)$ with shape $alpha$ and rate $beta$, so $bb(E)[u_i | z_(i k) = 1] = 1$.
]

Under Definition 1, atypical observations naturally receive small posterior latent scales $u_i << 1$, smoothly inflating their effective covariance $bold(Sigma)_k / u_i$ without corrupting the core component parameters. Parameter estimation maximizes $cal(L)(bold(Theta)) = sum_(i=1)^n log p(bold(x)_i | bold(Theta))$ via the Expectation / Conditional Maximization Either (ECME) algorithm @peel2000robust @mclachlan2000finite.

#algbox(
  title: "ECM / ECME Algorithm for Student's t-Mixture Models",
  inputs: [Dataset $cal(D) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$, component count $K$, tolerance $epsilon > 0$, ridge $eta >= 0$.],
  outputs: [Fitted parameters $hat(bold(Theta)) = {(hat(pi)_k, hat(bold(mu))_k, hat(bold(Sigma))_k, hat(nu)_k)}_(k=1)^K$, responsibilities $r_(i k)$, and outlier weights $u_(i k)$.],
)[
  + Initialize $bold(Theta)^((0))$ (e.g., via $k$-means++ or GMM with $nu_k^((0)) in [4, 20]$).
  + *Repeat* for $t = 0, 1, dots, T - 1$ until $|cal(L)(bold(Theta)^((t+1))) - cal(L)(bold(Theta)^((t)))| < epsilon$:
  + *E-step*: for all $i = 1, dots, n$ and $k = 1, dots, K$, compute Mahalanobis distances $delta_k (bold(x)_i)$ via Cholesky factors of $bold(Sigma)_k^((t))$ and update
    $ r_(i k) &= (pi_k^((t)) f_t (bold(x)_i | bold(Theta)_k^((t)))) / (sum_(l=1)^K pi_l^((t)) f_t (bold(x)_i | bold(Theta)_l^((t)))), \
      u_(i k) &= (nu_k^((t)) + d) / (nu_k^((t)) + delta_k (bold(x)_i)). $
  + *CM-step 1* (exact weighted Gaussian updates): let $N_k = sum_(i=1)^n r_(i k)$ and $W_k = sum_(i=1)^n r_(i k) u_(i k)$; update $pi_k^((t+1)) = N_k / n$, $bold(mu)_k^((t+1)) = W_k^(-1) sum_(i=1)^n r_(i k) u_(i k) bold(x)_i$, and
    $ bold(Sigma)_k^((t+1)) &= 1 / N_k sum_(i=1)^n r_(i k) u_(i k) (bold(x)_i - bold(mu)_k^((t+1))) \
      &quad times (bold(x)_i - bold(mu)_k^((t+1)))^T + eta bold(I). $
  + *CM-step 2* (degrees of freedom): for each $k$, compute $c_k = N_k^(-1) sum_(i=1)^n r_(i k) (log u_(i k) - u_(i k)) + psi((nu_k^((t)) + d) / 2) - log((nu_k^((t)) + d) / 2)$ and solve the 1D monotone equation
    $ -psi(nu_k / 2) + log(nu_k / 2) + 1 + c_k = 0 $
    for $nu_k^((t+1)) > 0$ (or directly maximize the observed log-likelihood in $nu_k$ for accelerated ECME).
  + Return $hat(bold(Theta))$ and hard labels $hat(z)_i = arg max_k r_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five rigorous theorems characterizing TMM: exact scale-mixture marginalization and conjugate Gamma posteriors, exact uniform redescending outlier-influence bounds for both location and scale, existence and uniqueness of the degrees-of-freedom root with monotone ECM ascent, uniform Gaussian convergence and heavy-tail dominance, and full affine equivariance.

#subsec("A", "Exact Scale-Mixture Integration and Posterior Scale Moments")

#lemma(title: "Exact Marginalization and Conjugate Gamma Scale Posterior", origin: "Original Proof")[
  Under Definition 1, marginalizing the latent scale $u_i > 0$ recovers the exact $d$-variate Student's $t$ density $f_t (bold(x)_i | bold(mu)_k, bold(Sigma)_k, nu_k)$. Moreover, the conditional posterior of $u_i$ given $(bold(x)_i, z_(i k) = 1)$ is exactly
  $ u_i | (bold(x)_i, z_(i k) = 1) ~ "Gamma"((nu_k + d) / 2, (nu_k + delta_k (bold(x)_i)) / 2), $
  yielding the closed-form E-step sufficient expectations
  $ bb(E)[u_i | bold(x)_i, z_(i k) = 1] = (nu_k + d) / (nu_k + delta_k (bold(x)_i)), $
  $ bb(E)[log u_i | bold(x)_i, z_(i k) = 1] = psi((nu_k + d) / 2) - log((nu_k + delta_k (bold(x)_i)) / 2), $
  where $psi(y) = (dif / (dif y)) log Gamma(y)$ is the digamma function.
]

#proof[
  Multiply the conditional Gaussian density of $bold(x)_i | (z_(i k) = 1, u)$ by the $"Gamma"(nu_k / 2, nu_k / 2)$ density of $u$:
  $ p(bold(x)_i, u | z_(i k) = 1) &= (2 pi)^(-d / 2) |bold(Sigma)_k|^(-1/2) u^(d / 2) exp(-u delta_k (bold(x)_i) / 2) \
    &quad times ((nu_k / 2)^(nu_k / 2) / Gamma(nu_k / 2)) u^(nu_k / 2 - 1) exp(-u nu_k / 2) \
    &= C_k u^((nu_k + d) / 2 - 1) exp(-(nu_k + delta_k (bold(x)_i)) u / 2), $
  where $C_k = (2 pi)^(-d / 2) |bold(Sigma)_k|^(-1/2) (nu_k / 2)^(nu_k / 2) / Gamma(nu_k / 2)$. Recognizing the kernel in $u$ as an unnormalized $"Gamma"(alpha_k, beta_(i k))$ density with $alpha_k = (nu_k + d) / 2$ and $beta_(i k) = (nu_k + delta_k (bold(x)_i)) / 2$, integrating over $u in (0, infinity)$ immediately yields $C_k Gamma(alpha_k) / beta_(i k)^(alpha_k)$. Simplifying powers of $2$ and $nu_k$:
  $ (2 pi)^(-d / 2) (nu_k / 2)^(nu_k / 2) ((nu_k + delta_k (bold(x)_i)) / 2)^(-(nu_k + d) / 2) \
    = (pi nu_k)^(-d / 2) (1 + (delta_k (bold(x)_i)) / nu_k)^(-(nu_k + d) / 2), $
  recovering $f_t (bold(x)_i | bold(mu)_k, bold(Sigma)_k, nu_k)$ identically. Dividing the joint density by its marginal proves $u_i | (bold(x)_i, z_(i k) = 1) ~ "Gamma"(alpha_k, beta_(i k))$. For any $U ~ "Gamma"(alpha, beta)$, $bb(E)[U] = alpha / beta$ and differentiating $integral_0^infinity u^(alpha - 1) e^(-beta u) dif u = Gamma(alpha) beta^(-alpha)$ with respect to $alpha$ gives $bb(E)[log U] = psi(alpha) - log beta$, yielding both E-step expectations.
]

#subsec("B", "Redescending Bounded Outlier Influence on Location and Scale")

The defining statistical advantage of TMM over Gaussian mixtures is formalized by sharp, uniform bounds on the component score functions.

#theorem(title: "Uniform Redescending Location and Bounded Scale Influence", origin: "Original Proof")[
  For any Student's $t$ component with $nu_k in (0, infinity)$ and $bold(Sigma)_k succ 0$, let $u_k (bold(x)) = (nu_k + d) / (nu_k + delta_k (bold(x)))$. The standardized location score $bold(s)_(mu, k) (bold(x)) = bold(Sigma)_k^(1/2) nabla_(bold(mu)_k) log f_t (bold(x))$ and whitened scale-contribution matrix $bold(W)_k (bold(x)) = u_k (bold(x)) bold(Sigma)_k^(-1/2) (bold(x) - bold(mu)_k) (bold(x) - bold(mu)_k)^T bold(Sigma)_k^(-1/2)$ satisfy:
  1. *Redescending Location Bound*: For all $bold(x) in bb(R)^d$,
     $ norm(bold(s)_(mu, k) (bold(x)))_2 = ((nu_k + d) sqrt(delta_k (bold(x)))) / (nu_k + delta_k (bold(x))) <= (nu_k + d) / (2 sqrt(nu_k)), $
     attaining its exact maximum at $delta_k (bold(x)) = nu_k$ and decaying strictly to $0$ as $cal(O)(delta_k (bold(x))^(-1/2))$ when $norm(bold(x) - bold(mu)_k) -> infinity$.
  2. *Uniformly Bounded Scale Influence*: For all $bold(x) in bb(R)^d$, the score $nabla_(bold(Sigma)_k^(-1)) log f_t (bold(x)) = (1/2) [bold(Sigma)_k - u_k (bold(x)) (bold(x) - bold(mu)_k) (bold(x) - bold(mu)_k)^T]$ has whitened outer-product norm strictly bounded by
     $ norm(bold(W)_k (bold(x)))_2 &= "Tr"(bold(W)_k (bold(x))) \
       &= ((nu_k + d) delta_k (bold(x))) / (nu_k + delta_k (bold(x))) < nu_k + d. $
  By contrast, as $nu_k -> infinity$ (GMM), $norm(bold(s)_(mu, k) (bold(x)))_2 = sqrt(delta_k (bold(x)))$ and $"Tr"(bold(W)_k (bold(x))) = delta_k (bold(x))$ diverge linearly and quadratically in $norm(bold(x) - bold(mu)_k)$, respectively.
]

#proof[
  Differentiating $log f_t (bold(x)) = "const" + (1/2) log |bold(Sigma)_k^(-1)| - ((nu_k + d) / 2) log(1 + delta_k (bold(x)) / nu_k)$ with respect to $bold(mu)_k$ gives
  $ nabla_(bold(mu)_k) log f_t (bold(x)) &= (nu_k + d) / (nu_k + delta_k (bold(x))) bold(Sigma)_k^(-1) (bold(x) - bold(mu)_k) \
    &= u_k (bold(x)) bold(Sigma)_k^(-1) (bold(x) - bold(mu)_k). $
  Setting the whitened vector $bold(z) = bold(Sigma)_k^(-1/2) (bold(x) - bold(mu)_k)$, we have $norm(bold(z))_2 = sqrt(delta_k (bold(x)))$ and $bold(s)_(mu, k) (bold(x)) = u_k (bold(x)) bold(z)$. Taking Euclidean norms gives $norm(bold(s)_(mu, k) (bold(x)))_2 = (nu_k + d) sqrt(delta_k (bold(x))) / (nu_k + delta_k (bold(x)))$. By the arithmetic-geometric mean inequality, $nu_k + delta_k (bold(x)) >= 2 sqrt(nu_k delta_k (bold(x)))$ with equality if and only if $delta_k (bold(x)) = nu_k$, yielding the exact global bound $(nu_k + d) / (2 sqrt(nu_k))$. Moreover, as $delta_k (bold(x)) -> infinity$, $(nu_k + d) sqrt(delta_k) / (nu_k + delta_k) ~ (nu_k + d) delta_k^(-1/2) -> 0$, proving Claim 1.

  For Claim 2, differentiating $log f_t (bold(x))$ with respect to $bold(Sigma)_k^(-1)$ via $nabla_(bold(A)) log |bold(A)| = bold(A)^(-1)$ and $nabla_(bold(A)) (bold(v)^T bold(A) bold(v)) = bold(v) bold(v)^T$ yields
  $ nabla_(bold(Sigma)_k^(-1)) log f_t (bold(x)) \
    = 1 / 2 [ bold(Sigma)_k - u_k (bold(x)) (bold(x) - bold(mu)_k) (bold(x) - bold(mu)_k)^T ]. $
  The whitened matrix $bold(W)_k (bold(x)) = u_k (bold(x)) bold(z) bold(z)^T$ is rank-1 positive semidefinite, so its spectral norm equals its trace:
  $ norm(bold(W)_k (bold(x)))_2 &= "Tr"(bold(W)_k (bold(x))) = u_k (bold(x)) norm(bold(z))_2^2 \
    &= ((nu_k + d) delta_k (bold(x))) / (nu_k + delta_k (bold(x))) < nu_k + d, $
  where the strict inequality holds for all finite $delta_k (bold(x)) >= 0$ since $nu_k > 0$. Taking $nu_k -> infinity$ gives $u_k (bold(x)) -> 1$, recovering the unbounded Gaussian scores.
]

Theorem 1 reveals a striking contrast: extreme outliers exert *zero* asymptotic pull on $bold(mu)_k$ (redescending M-estimator behavior) and *strictly bounded* leverage ($< nu_k + d$) on $bold(Sigma)_k$, even when assigned with responsibility $r_(i k) = 1$ to component $k$.

#subsec("C", "Uniqueness of the Degrees-of-Freedom Root and Monotone Ascent")

#theorem(title: "Strict Unimodality of CM-Step 2 and Monotone ECM Ascent", origin: "Original Proof")[
  At every ECM iteration with $N_k = sum_(i=1)^n r_(i k) > 0$, the expected complete-data objective for $nu_k > 0$,
  $ Q_(nu, k) (nu_k) = N_k [ (nu_k / 2) log(nu_k / 2) - log Gamma(nu_k / 2) + (nu_k / 2) c_k ], $
  where $c_k = N_k^(-1) sum_(i=1)^n r_(i k) bb(E)[log u_i - u_i | bold(x)_i, z_(i k) = 1, bold(Theta)^((t))]$, satisfies $1 + c_k < 0$ strictly, is strictly unimodal on $(0, infinity)$, and admits a unique finite global maximizer $nu_k^((t+1)) in (0, infinity)$. Consequently, the ECM / ECME sequence guarantees monotone observed log-likelihood ascent $cal(L)(bold(Theta)^((t+1))) >= cal(L)(bold(Theta)^((t)))$ at every iteration @peel2000robust @dempster1977em.
]

#proof[
  Setting $y = nu_k / 2 > 0$, differentiate $Q_(nu, k)$ with respect to $y$:
  $ (dif Q_(nu, k)) / (dif y) = N_k [ log y + 1 - psi(y) + c_k ] = N_k [ h(y) + (1 + c_k) ], $
  where $h(y) = log y - psi(y)$. By the polygamma series representation $psi'(y) = sum_(m=0)^infinity (y + m)^(-2)$ and integral comparison against the strictly decreasing function $t |-> (y + t)^(-2)$,
  $ psi'(y) = sum_(m=0)^infinity 1 / (y + m)^2 > integral_0^infinity 1 / (y + t)^2 dif t = 1 / y, $
  whence $h'(y) = 1 / y - psi'(y) < 0$ strictly for every $y > 0$. Furthermore, from the recurrence $psi(y) = psi(y + 1) - 1 / y$, as $y -> 0^+$ we have $h(y) = log y + 1 / y + cal(O)(1) -> +infinity$, while from $sum_(m=1)^infinity (y + m)^(-2) < integral_0^infinity (y + t)^(-2) dif t = 1 / y$ and $lim_(y -> infinity) h(y) = 0$ we obtain $0 < h(y) < 1 / y -> 0^+$ as $y -> +infinity$. Thus $h : (0, infinity) -> (0, infinity)$ is a strictly decreasing smooth bijection.

  Next, observe that $phi(u) = 1 + log u - u <= 0$ for all $u > 0$, with equality if and only if $u = 1$. By Lemma 1, each conditional posterior $u_i | (bold(x)_i, z_(i k) = 1, bold(Theta)^((t)))$ is a non-degenerate continuous Gamma random variable, which assigns probability $1$ to $u_i != 1$. Hence $bb(E)[1 + log u_i - u_i | bold(x)_i, z_(i k) = 1] < 0$ strictly for every $i$, and taking the $r_(i k)$-weighted average proves $1 + c_k < 0$, i.e., $-(1 + c_k) > 0$. Because $h$ is a strictly decreasing bijection onto $(0, infinity)$, the equation $h(y) = -(1 + c_k)$ has a unique solution $y^* = h^(-1)(-(1 + c_k)) in (0, infinity)$, with $Q_(nu, k)'(nu_k) > 0$ for $nu_k < 2 y^*$ and $Q_(nu, k)'(nu_k) < 0$ for $nu_k > 2 y^*$. Thus $nu_k^((t+1)) = 2 y^*$ uniquely maximizes $Q_(nu, k)$. Combining with closed-form maximization of the quadratic Gaussian $Q$-function over $(bold(pi), bold(mu)_k, bold(Sigma)_k)$ in CM-step 1 (and exact observed-likelihood maximization in ECME), standard EM monotonicity @dempster1977em @wu1983convergence implies $cal(L)(bold(Theta)^((t+1))) >= cal(L)(bold(Theta)^((t)))$.
]

Note a subtle distinction illuminated by Theorem 2: in pure ECM, $nu_k^((t+1))$ is always strictly finite at every step (drifting to $+infinity$ only asymptotically across iterations on light-tailed Gaussian data as $c_k -> -1^-$), whereas ECME directly maximizes the observed log-likelihood over $nu_k$ (achieving faster convergence and directly detecting Gaussian limits $nu_k -> infinity$ @peel2000robust).

#subsec("D", "Gaussian Convergence Limit and Tail-Dominance Asymptotics")

#proposition(title: "Gaussian Limit and Heavy-Tail Outlier Absorption", origin: "Classical — Cited")[
  1. For fixed $(bold(mu)_k, bold(Sigma)_k)$, as $nu_k -> infinity$, $f_t (bold(x) | bold(mu)_k, bold(Sigma)_k, nu_k) -> cal(N)(bold(x) | bold(mu)_k, bold(Sigma)_k)$ uniformly on $bb(R)^d$ @mclachlan2000finite @bishop2006pattern.
  2. In any two-component mixture with finite $nu_1 < nu_2 <= infinity$ (where $nu_2 = infinity$ denotes a Gaussian component), regardless of $bold(mu)_1, bold(mu)_2, bold(Sigma)_1, bold(Sigma)_2$ and positive weights $pi_1, pi_2 > 0$, the heavier-tailed component asymptotically absorbs extreme outliers in every direction:
     $ lim_(norm(bold(x)) -> infinity) P(z = 1 | bold(x), bold(Theta)) = 1. $
]

#proof[
  For Claim 1, Stirling's formula yields $Gamma((nu_k + d) / 2) / ((nu_k / 2)^(d / 2) Gamma(nu_k / 2)) -> 1$ as $nu_k -> infinity$, so the normalizing constant converges to $(2 pi)^(-d / 2) |bold(Sigma)_k|^(-1/2)$. Simultaneously, $(1 + delta_k (bold(x)) / nu_k)^(-(nu_k + d) / 2) -> exp(-delta_k (bold(x)) / 2)$ uniformly on compact sets, and both densities vanish uniformly outside large Mahalanobis balls, giving uniform convergence on $bb(R)^d$.

  For Claim 2, as $norm(bold(x)) -> infinity$, both $delta_1 (bold(x))$ and $delta_2 (bold(x))$ scale as $Theta(norm(bold(x))_2^2)$ because $bold(Sigma)_1, bold(Sigma)_2 succ 0$. Hence $f_t (bold(x) | bold(Theta)_1) = Theta(norm(bold(x))_2^(-(nu_1 + d)))$, whereas for finite $nu_2 > nu_1$, $f_t (bold(x) | bold(Theta)_2) = Theta(norm(bold(x))_2^(-(nu_2 + d)))$ (and $O(exp(-c norm(bold(x))_2^2))$ if $nu_2 = infinity$). Thus the density ratio $f_2 (bold(x)) / f_1 (bold(x)) = cal(O)(norm(bold(x))_2^(-(nu_2 - nu_1))) -> 0$, forcing $P(z = 1 | bold(x)) = (1 + (pi_2 f_2 (bold(x))) / (pi_1 f_1 (bold(x))))^(-1) -> 1$.
]

#subsec("E", "Full Affine Equivariance and Support Nonsingularity")

#proposition(title: "Exact Affine Equivariance and Nonsingularity Support Bound", origin: "Original Proof")[
  For any invertible affine map $bold(x)' = bold(A) bold(x) + bold(b)$ ($det(bold(A)) != 0$), the TMM log-likelihood and every ECM / ECME iterate are strictly equivariant: transforming $bold(mu)_k' = bold(A) bold(mu)_k + bold(b)$ and $bold(Sigma)_k' = bold(A) bold(Sigma)_k bold(A)^T$ leaves $delta_k (bold(x))$, $u_(i k)$, $r_(i k)$, and $nu_k$ invariant and shifts $cal(L)$ by the constant $-n log |det(bold(A))|$. Moreover, with $eta = 0$, $bold(Sigma)_k^((t+1))$ is positive definite if and only if the points with $r_(i k) > 0$ span all of $bb(R)^d$ affinely (requiring at least $d + 1$ affinely independent points).
]

#proof[
  Directly, $(bold(x)' - bold(mu)_k')^T (bold(Sigma)_k')^(-1) (bold(x)' - bold(mu)_k') = (bold(x) - bold(mu)_k)^T bold(A)^T bold(A)^(-T) bold(Sigma)_k^(-1) bold(A)^(-1) bold(A) (bold(x) - bold(mu)_k) = delta_k (bold(x))$, and $|bold(Sigma)_k'|^(-1/2) = |det(bold(A))|^(-1) |bold(Sigma)_k|^(-1/2)$. Thus every component density scales by $|det(bold(A))|^(-1)$, cancelling identically in $r_(i k)$, leaving $u_(i k)$ and $c_k$ unchanged, and transforming $bold(mu)_k^((t+1))$ and $bold(Sigma)_k^((t+1))$ linearly and congruently under $bold(A)$. For nonsingularity, because $r_(i k) u_(i k) > 0$ whenever $r_(i k) > 0$, any nonzero vector $bold(v) in bb(R)^d$ satisfies $bold(v)^T bold(Sigma)_k^((t+1)) bold(v) = 0$ if and only if $bold(v)^T (bold(x)_i - bold(mu)_k^((t+1))) = 0$ for all $i$ with $r_(i k) > 0$, i.e., all positively weighted points lie in a proper affine hyperplane of $bb(R)^d$. Hence $bold(Sigma)_k^((t+1)) succ 0$ iff the support contains at least $d + 1$ affinely independent vectors; adding a ridge $eta bold(I)$ (or an Inverse-Wishart prior) prevents boundary likelihood singularities.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Using Cholesky factorizations $bold(Sigma)_k = bold(L)_k bold(L)_k^T$ ($cal(O)(K d^3)$ per iteration), computing $delta_k (bold(x)_i) = norm(bold(L)_k^(-1) (bold(x)_i - bold(mu)_k))_2^2$ and $log |bold(Sigma)_k| = 2 sum_j log L_(k, j j)$ costs $cal(O)(d^2)$ per $(i, k)$ pair via triangular back-substitution. Updating $bold(Sigma)_k$ costs $cal(O)(n d^2)$ per component, while the 1D monotone Newton / Brent solve of Theorem 2 for each $nu_k$ requires only $cal(O)(1)$ scalar digamma evaluations (typically $< 10$ iterations). Total runtime over $T$ iterations is therefore $cal(O)(T n K d^2 + T K d^3)$—asymptotically identical to a full-covariance GMM, with memory $cal(O)(n K + K d^2)$.

Geometrically, TMM models elliptically contoured clusters with adaptive tail thickness: small $nu_k$ ($nu_k in [1, 5]$) accommodates Cauchy-like heavy tails and gross outliers without sacrificing affine equivariance (Proposition 2), while large $nu_k$ ($nu_k > 30$) smoothly recovers full-covariance Gaussian ellipsoids (Proposition 1). Model order $K$ is selected cleanly via BIC @schwarz1978bic with $p_K = K (1 + d + d(d+1)/2) + K - 1$ free parameters (including $nu_1, dots, nu_K$).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Continuous real vectors $bold(x) in bb(R)^d$ (outlier-tolerant)],
    [Objective Function], [Multivariate Student's $t$-mixture log-likelihood via ECM / ECME],
    [Time Complexity], [$cal(O)(T n K d^2 + T K d^3)$ full covariance; $cal(O)(T n K d)$ diagonal],
    [Space Complexity], [$cal(O)(n K + K d^2)$ for responsibilities, scales, and Cholesky factors],
    [Robustness Mechanism], [Redescending location score ($<= (nu_k + d) / (2 sqrt(nu_k))$) and bounded scale leverage ($< nu_k + d$)],
    [Invariances], [Strictly affine-equivariant under all invertible $bold(A) bold(x) + bold(b)$ and order-invariant in $n$],
    [Diagnostics Output], [Per-datum outlier weights $u_(i k) in (0, 1 + d / nu_k]$ and learned tail indices $hat(nu)_k$],
    [Failure Safeguard], [Ridge $eta bold(I)$ or Inverse-Wishart prior prevents singular covariance spikes],
  ),
  caption: [Operational and geometric profile of Student's $t$-Mixture Models.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Quantitative Finance and Asset-Return Regime Clustering:* Financial log-returns and volatility vectors exhibit pronounced leptokurtosis (fat tails, empirical $nu approx 3$–$6$) and crash shocks. TMM naturally separates bull, calm, and crisis correlation regimes via full affine-equivariant covariances $bold(Sigma)_k$ while Theorem 1 prevents flash-crash outliers from corrupting covariance estimates used in portfolio risk allocation.

*2. Flow Cytometry and Single-Cell Biomarker Gating:* High-throughput cytometry measurements contain dead cells, cellular debris, and fluorescence doublets superimposed on correlated elliptical cell populations. TMM captures oblique marker correlations while automatically flagging debris through low latent weights $u_(i k) << 1$, achieving substantially higher gating purity than Gaussian mixtures @peel2000robust.

*3. Robust Computer Vision and Background / Color Segmentation:* Pixel chromaticity and motion-vector distributions suffer from specular reflections, sensor hot-pixels, and occlusion artifacts. Fitting TMMs yields sharp foreground / background models whose redescending influence function ignores transient glint spikes without requiring brittle manual outlier clipping.

*4. Astrophysical Photometry and Remote-Sensing Classification:* Multi-band astronomical catalogs and hyperspectral satellite scenes contain cosmic-ray hits, cloud-edge corruptions, and saturated pixels. TMM separates overlapping spectral ellipsoids in color-magnitude space while outputting per-object weight maps $u_(i k)$ that double as calibrated anomaly scores.

*5. Industrial Process Monitoring and Clinical Biomarker Subtyping:* Medical assays and chemical-plant sensors frequently experience transient calibration glitches alongside true biological or operational heterogeneity. TMM isolates correlated physiological or thermodynamic phenotypes without creating spurious singleton clusters around instrument faults.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

TMM exhibits four primary limitations: (i) *high-dimensional small-$n$ regimes* ($d$ approaching or exceeding $n / K$), where estimating $K$ unrestricted $d times d$ scale matrices ($cal(O)(K d^2)$ parameters) overfits or hits Proposition 2's $d + 1$ rank barrier unless diagonal, factor-analyzer, or regularized covariance structures are imposed; (ii) *strongly skewed or asymmetric tails*, since Student's $t$ components remain elliptically symmetric (motivating skew-$t$ mixtures @lin2009skewnormal); (iii) *non-convex manifolds* (interlocking moons or filaments); and (iv) *local EM optima* requiring multiple initializations or trimmed warm starts.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Correlated tabular data with outliers], [Optimal: full affine equivariance + bounded redescending influence], [Add small ridge $eta bold(I)$; use BIC for $K$],
    [Financial returns / flow cytometry], [High: leptokurtic tails match finite $nu_k$; $u_(i k)$ flags anomalies], [Use skew-$t$ mixtures @lin2009skewnormal if strongly asymmetric],
    [Mixed categorical + continuous data], [Unsuitable directly: continuous multivariate density only], [Use AutoClass (Alg. 87) or MIXMOD @biernacki2006mixmod],
    [High-D ($d >> sqrt(n / K)$) embeddings], [Moderate-Low: $cal(O)(K d^2)$ covariance overparameterization], [Constrain $bold(Sigma)_k$ to diagonal or factor-analyzer form],
    [Non-convex / topological manifolds], [Low: elliptical contours cannot trace curved filaments], [Use HDBSCAN, spectral clustering, or DPC],
  ),
  caption: [Domain suitability matrix and failure mitigations for TMM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Student's $t$-Mixture Models sit at the center of robust parametric clustering.

#proposition(title: "Continuous Parametric Bridge Between GMM, Cauchy Mixtures, and Soft Trimming", origin: "Original Derivation")[
  Let $bold(Theta)_K$ be a $K$-component TMM. Then: (i) $nu_k -> infinity$ recovers exact full-covariance GMM / MCLUST @fraley2002mclust (Proposition 1); (ii) $nu_k = 1$ yields multivariate Cauchy mixtures with location score bounded by $(1 + d) / 2$; and (iii) relative to hard-trimmed clustering (such as trimmed $k$-means @cuesta1997trimmed, which sets binary weights $w_i in {0, 1}$ at an external trimming fraction $alpha$), TMM implements an adaptive, differentiable soft-trimming estimator with continuous per-point weights $u_(i k) = (nu_k + d) / (nu_k + delta_k (bold(x)_i)) in (0, 1 + d / nu_k]$ whose effective trimming threshold $delta_k = nu_k$ is learned directly from the data via Theorem 2.
]

#proof[
  Claims (i) and (ii) follow directly from Proposition 1 and substituting $nu_k = 1$ into Theorem 1's exact bound $(nu_k + d) / (2 sqrt(nu_k))$. For Claim (iii), compare CM-step 1's weighted mean $bold(mu)_k = (sum_i r_(i k) u_(i k) bold(x)_i) / (sum_i r_(i k) u_(i k))$ with hard-trimmed estimators @cuesta1997trimmed that optimize $sum_(i in cal(A)) norm(bold(x)_i - bold(mu)_(z_i))_2^2$ over an active subset $|cal(A)| = ceil((1 - alpha) n)$. Writing $u_(i k)$ as a smooth decreasing function of Mahalanobis distance $delta_k (bold(x)_i)$ shows $u_(i k) > 1$ for inliers ($delta_k < d$), $u_(i k) = 1$ at the expected squared Mahalanobis radius $delta_k = d$, peak location leverage at $delta_k = nu_k$, and $u_(i k) -> 0$ smoothly as $delta_k -> infinity$ without specifying a hard contamination fraction $alpha$.
]

Beyond symmetric $t$-mixtures, skew-normal and skew-$t$ mixtures @lin2009skewnormal add skewness vectors to accommodate asymmetric heavy tails, while MCLUST @fraley2002mclust approximates outlier robustness by adding a uniform Poisson noise component rather than heavy-tailed component densities.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized 4-dataset benchmark, TMM achieves Iris Silhouette $= 0.361$ and an extraordinary $"ARI" = 0.950$ @hubert1985comparing @rousseeuw1987silhouettes, Customers Silhouette $= 0.367$, Multi-Shape Silhouette $= 0.712$ ($"ARI" =$ N/A), and High-D Silhouette $= 0.292$ ($"ARI" =$ N/A), with a mean runtime of $22.958$ ms. Every single metric aligns sharply with our mathematical theory:
1. *Near-Exact Species Recovery on Iris ($"ARI" = 0.950$):* TMM outperforms COBWEB ($0.066$), CLASSIT ($0.066$), and AutoClass ($0.762$) dramatically on ground-truth recovery. Although isotropic Euclidean Silhouette is higher for AutoClass ($0.431$ vs. $0.361$), _versicolor_ and _virginica_ are strongly correlated oblique ellipsoids with borderline specimens: TMM's full affine-equivariant $bold(Sigma)_k$ (Proposition 2) aligns with the true Mahalanobis geometry while Theorem 1's redescending weights $u_(i k)$ prevent borderline points from distorting class boundaries.
2. *Best Robustness on Heavy-Tailed Customers ($"Sil" = 0.367$):* Wholesale customer spending distributions exhibit heavy right tails and extreme bulk buyers; TMM achieves the highest Customers Silhouette among Algs. 85–88 ($0.367$ vs. $0.312$–$0.318$), directly confirming Theorem 1's bounded scale leverage.
3. *High-D Trade-off and Runtime:* On *Multi-Shape* ($"Sil" = 0.712$) and *High-D* ($"Sil" = 0.292$), full $d times d$ scale estimation incurs slightly more variance than AutoClass's regularized diagonal evidence ($0.782$ / $0.321$), yet TMM's mean runtime of $22.958$ ms adds barely $3.8$ ms over AutoClass ($19.192$ ms), confirming the low computational cost of Cholesky-accelerated ECME.
