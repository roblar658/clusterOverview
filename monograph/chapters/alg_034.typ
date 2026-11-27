#import "../template/ieee.typ": *

#alg_entry(
  id: 34,
  name: "Mean Shift",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1975 / 1995 / 2002",
  paradigm: "Nonparametric Kernel Mode-Seeking",
  complexity: $O(n^2 d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset in Euclidean space equipped with norm $norm(dot)_2$, and let $h > 0$ be a scalar bandwidth parameter. Rather than imposing a fixed cluster count $k$ or spherical Voronoi cells, *Mean Shift* @cheng1995mean @comaniciu2002mean models clusters nonparametrically as the *basins of attraction* of local maxima (modes) of a smooth multivariate kernel density estimator (KDE). Given a radially symmetric kernel $K(u) = c_(k, d) k(norm(u)_2^2)$ with normalization constant $c_(k, d) > 0$ and non-increasing, piecewise smooth profile $k : [0, infinity) -> [0, infinity)$ satisfying $integral_0^infinity k(r) r^(d/2 - 1) dif r < infinity$, the empirical density estimate at $x in bb(R)^d$ is
$ hat(f)_(h, K) (x) = c_(k, d) / (n h^d) sum_(i=1)^n k( norm((x - x_i) / h)_2^2 ). $
Define the *shadow profile* $g(r) := -k'(r) >= 0$ for all $r >= 0$ where $k'$ exists, with associated shadow kernel $G(u) = c_(g, d) g(norm(u)_2^2)$ @cheng1995mean @comaniciu2002mean. Writing $g_i (x) := g(norm((x - x_i)/h)_2^2)$, differentiating $hat(f)_(h, K) (x)$ yields the fundamental gradient identity:
$ nabla hat(f)_(h, K) (x) &= (2 c_(k, d)) / (n h^(d+2)) sum_(i=1)^n (x_i - x) g_i (x) \
  &= (2 c_(k, d)) / (n h^(d+2)) [ sum_(i=1)^n g_i (x) ] m_(h, G) (x), $
where the *mean shift vector* $m_(h, G) (x) in bb(R)^d$ is defined wherever $sum_(i=1)^n g_i (x) > 0$ by
$ m_(h, G) (x) := ( sum_(i=1)^n x_i g_i (x) ) / ( sum_(i=1)^n g_i (x) ) - x. $
Rearranging reveals that $m_(h, G) (x)$ is positively collinear with $nabla hat(f)_(h, K) (x)$ with a density-adaptive step size:
$ m_(h, G) (x) = (h^2 c_(g, d)) / (2 c_(k, d)) dot (nabla hat(f)_(h, K) (x)) / (hat(f)_(h, G) (x)). $
Thus the fixed-point iteration $y^((t+1)) = T_h (y^((t))) := y^((t)) + m_(h, G) (y^((t)))$ performs *automatic variable-step gradient ascent*, taking large steps across low-density valleys where $hat(f)_(h, G) (y)$ is small and damped steps near high-density modes without requiring an external learning-rate schedule.

#definition(title: "Mode-Seeking Basins and Induced Partition")[
  Let $T_h : bb(R)^d -> bb(R)^d$ be the mean shift operator $T_h (y) = (sum_(i=1)^n x_i g_i (y)) / (sum_(i=1)^n g_i (y))$. A point $mu^* in bb(R)^d$ is an *attracting density mode* if $T_h (mu^*) = mu^*$ and $nabla^2 hat(f)_(h, K) (mu^*) prec 0$. The *basin of attraction* of $mu^*$ is $cal(B)_h (mu^*) := {x_i in cal(X) : lim_(t -> infinity) T_h^t (x_i) = mu^*}$. The *Mean Shift clustering* partitions $cal(X)$ by merging modes within numerical tolerance $epsilon_"merge" > 0$ and grouping points sharing the same merged attractor.
]

Canonical profiles include the Gaussian profile $k(r) = exp(-r / 2)$ (self-shadowing, $g(r) = (1/2) k(r)$, $C^infinity$-smooth) and the Epanechnikov profile $k(r) = max(0, 1 - r)$ (compact support on $B_h$, whose shadow $g(r) = bold(1)_([0, 1]) (r)$ is the uniform flat kernel). Algorithm 1 details the exact batch mode-seeking procedure.

#algbox(
  title: "Mean Shift Mode-Seeking Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, bandwidth $h > 0$, shadow profile $g$, shift tolerance $tau > 0$, merge radius $epsilon_"merge" > 0$, max steps $T$.],
  outputs: [Endogenous mode set $cal(M) = {mu_1, dots, mu_hat(k)}$ and partition $cal(C)$.],
)[
  1. *For* each datum $i = 1, dots, n$ (in parallel) *do*:
  2. #h(1em) Initialize trajectory $y_i^((0)) <- x_i$.
  3. #h(1em) *For* $t = 0, 1, dots, T - 1$ *do*:
  4. #h(2em) Compute kernel weights $w_(i j)^((t)) <- g(norm((y_i^((t)) - x_j)/h)_2^2)$ for all $j in {1, dots, n}$.
  5. #h(2em) Update $y_i^((t+1)) <- (sum_(j=1)^n w_(i j)^((t)) x_j) / (sum_(j=1)^n w_(i j)^((t)))$.
  6. #h(2em) *If* $norm(y_i^((t+1)) - y_i^((t)))_2 <= tau$, *set* $hat(y)_i <- y_i^((t+1))$ and *break*.
  7. *Mode consolidation:* Form the graph on ${hat(y)_1, dots, hat(y)_n}$ linking $(i, j)$ whenever $norm(hat(y)_i - hat(y)_j)_2 <= epsilon_"merge"$; extract connected components $C_1, dots, C_hat(k)$ and set mode $mu_r <- |C_r|^(-1) sum_(i in C_r) hat(y)_i$.
  8. *Return* $cal(M) = {mu_1, dots, mu_hat(k)}$ and partition $cal(C) = {C_1, dots, C_hat(k)}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of Mean Shift: (1) an exact Majorization–Minimization (MM) quadratic minorizer proving strict density ascent and square-summable shift convergence; (2) exact EM equivalence and the closed-form local Jacobian $J_(T_h) (y^*) = h^(-2) Sigma_(w^*)$ governing non-oscillatory linear convergence; (3) an exact pitchfork phase transition characterizing bandwidth-induced mode collapse; (4) rigid-motion invariance and outlier localization; and (5) asymptotic gradient-line consistency under smooth Morse densities.

#subsec("A", "Majorization–Minimization Ascent and Stationarity")

#lemma(title: "Exact Quadratic Minorizer for Convex Profiles", origin: "Original Proof")[
  Suppose the kernel profile $k : [0, infinity) -> [0, infinity)$ is convex and differentiable on $[0, infinity)$ with shadow $g(r) = -k'(r) >= 0$. For any $y in bb(R)^d$ with $W(y) := sum_(i=1)^n g_i (y) > 0$, define the quadratic surrogate
  $ Q(z | y) &:= hat(f)_(h, K) (y) \
    &quad - c_(k, d) / (n h^(d+2)) sum_(i=1)^n g_i (y) (norm(z - x_i)_2^2 - norm(y - x_i)_2^2). $
  Then (i) $Q(y | y) = hat(f)_(h, K) (y)$, (ii) $Q(z | y) <= hat(f)_(h, K) (z)$ for all $z in bb(R)^d$, and (iii) $T_h (y)$ is the unique global maximizer of $z |-> Q(z | y)$, satisfying the quantitative ascent bound
  $ hat(f)_(h, K) (T_h (y)) - hat(f)_(h, K) (y) >= (c_(k, d) W(y)) / (n h^(d+2)) norm(T_h (y) - y)_2^2. $
]

#proof[
  Identity (i) is immediate upon setting $z = y$. Because $k$ is convex and differentiable on $[0, infinity)$, its tangent line globally minorizes $k$: for every $r, r' >= 0$, $k(r') >= k(r) + k'(r)(r' - r) = k(r) - g(r)(r' - r)$. Substituting $r' = norm((z - x_i)/h)_2^2$ and $r = norm((y - x_i)/h)_2^2$, summing over $i = 1, dots, n$, and scaling by $c_(k, d) / (n h^d) > 0$ yields
  $ hat(f)_(h, K) (z) &>= hat(f)_(h, K) (y) \
    &quad - c_(k, d) / (n h^(d+2)) sum_(i=1)^n g_i (y) (norm(z - x_i)_2^2 - norm(y - x_i)_2^2) \
    &= Q(z | y), $
  proving (ii). To establish (iii), observe that $Q(z | y)$ is strictly concave in $z$ with Hessian $nabla_z^2 Q(z | y) = -(2 c_(k, d) W(y) / (n h^(d+2))) I_d prec 0$. Applying the weighted parallel-axis identity around the barycenter $y^+ := T_h (y) = W(y)^(-1) sum_(i=1)^n g_i (y) x_i$ gives
  $ sum_(i=1)^n g_i (y) norm(y - x_i)_2^2 &= sum_(i=1)^n g_i (y) norm(y^+ - x_i)_2^2 \
    &quad + W(y) norm(y - y^+)_2^2. $
  Thus $z = y^+$ uniquely maximizes $Q(z | y)$, and combining (ii) with this exact identity gives
  $ hat(f)_(h, K) (y^+) - hat(f)_(h, K) (y) &>= Q(y^+ | y) - Q(y | y) \
    &= (c_(k, d) W(y)) / (n h^(d+2)) norm(y^+ - y)_2^2, $
  completing the proof.
]

#theorem(title: "Monotonic Ascent, Convex-Hull Confinement & Rate", origin: "Original Proof")[
  Let $k$ be convex, continuously differentiable, and strictly decreasing on its support with $k(0) < infinity$. For any initialization $y^((0)) in bb(R)^d$ with $W(y^((0))) > 0$, the iterates $y^((t+1)) = T_h (y^((t)))$ satisfy:
  1. $y^((t)) in "conv"(cal(X))$ for every $t >= 1$, hence ${y^((t))}$ is bounded.
  2. The density sequence ${hat(f)_(h, K) (y^((t)))}$ is monotonically non-decreasing and converges to a finite limit $f^* <= c_(k, d) k(0) / h^d$.
  3. Suppose $W_"min" := inf_(t >= 0) W(y^((t))) > 0$ (automatically guaranteed for the Gaussian profile since $"conv"(cal(X))$ is compact). Then $sum_(t=0)^infinity norm(m_(h, G) (y^((t))))_2^2 < infinity$,
     $ min_(0 <= s < t) norm(m_(h, G) (y^((s))))_2^2 <= (n h^(d+2) (f^* - hat(f)_(h, K) (y^((0))))) / (t c_(k, d) W_"min"), $
     and every accumulation point $y^*$ of ${y^((t))}$ is a stationary point ($nabla hat(f)_(h, K) (y^*) = bold(0)$).
]

#proof[
  (1) Because $g_i (y^((t))) >= 0$ and $W(y^((t))) > 0$, the normalized weights $alpha_i^((t)) = g_i (y^((t))) / W(y^((t)))$ satisfy $alpha_i^((t)) >= 0$ and $sum_i alpha_i^((t)) = 1$, so $y^((t+1)) = sum_i alpha_i^((t)) x_i in "conv"(cal(X))$.

  (2) By Lemma 1, $hat(f)_(h, K) (y^((t+1))) >= hat(f)_(h, K) (y^((t)))$ with strict inequality whenever $y^((t+1)) != y^((t))$. Since $k$ is non-increasing, $hat(f)_(h, K) (y) <= c_(k, d) k(0) / h^d$ globally; by the monotone convergence theorem, $hat(f)_(h, K) (y^((t))) -> f^* < infinity$.

  (3) Summing Lemma 1's ascent inequality from $s = 0$ to $t - 1$ telescopes:
  $ sum_(s=0)^(t-1) norm(y^((s+1)) - y^((s)))_2^2 &<= (n h^(d+2)) / (c_(k, d) W_"min") \
    &quad times (hat(f)_(h, K) (y^((t))) - hat(f)_(h, K) (y^((0)))). $
  Taking $t -> infinity$ proves $sum_(s=0)^infinity norm(m_(h, G) (y^((s))))_2^2 < infinity$, and dividing by $t$ yields the $O(1/t)$ minimum squared shift bound. By compactness of $"conv"(cal(X))$, ${y^((t))}$ admits at least one convergent subsequence $y^((t_r)) -> y^*$. Since $m_(h, G) (y^((t_r))) -> bold(0)$ and $g$ is continuous, $m_(h, G) (y^*) = bold(0)$, whence $nabla hat(f)_(h, K) (y^*) = bold(0)$.
]

#remark(title: "Point Convergence vs. Limit-Point Stationarity")[
  Classical analyses @comaniciu2002mean inferred single-point convergence $y^((t)) -> y^*$ directly from $norm(y^((t+1)) - y^((t)))_2 -> 0$, which strictly requires isolated stationary points or real-analyticity @ariascastro2016meanshift. For the Gaussian kernel, $hat(f)_(h, K)$ is real-analytic, so the Łojasiewicz gradient inequality guarantees finite trajectory length $sum_t norm(y^((t+1)) - y^((t)))_2 < infinity$ and genuine point convergence.
]

#subsec("B", "Gaussian EM Identity and Exact Local Contraction Spectrum")

#proposition(title: "EM Equivalence & Exact Spectral Radius at Modes", origin: "Original Derivation")[
  For the Gaussian profile $k(r) = exp(-r / 2)$, the Mean Shift update $y^((t+1)) = T_h (y^((t)))$ is algebraically identical to an Expectation–Maximization (EM) step @dempster1977em @carreira2006fast maximizing $log hat(f)_(h, K) (y)$. Furthermore, at any strict local maximum $y^*$ ($nabla hat(f)_(h, K) (y^*) = bold(0)$, $H^* := nabla^2 hat(f)_(h, K) (y^*) prec 0$), the Jacobian $J_(T_h) (y^*)$ satisfies the exact covariance identity
  $ J_(T_h) (y^*) = 1 / h^2 Sigma_(w^*), quad Sigma_(w^*) := sum_(i=1)^n w_i^* (x_i - y^*)(x_i - y^*)^T, $
  where $w_i^* = exp(-norm(x_i - y^*)_2^2 / (2 h^2)) / sum_(j=1)^n exp(-norm(x_j - y^*)_2^2 / (2 h^2))$. Consequently, every eigenvalue of $J_(T_h) (y^*)$ lies strictly in $[0, 1)$, guaranteeing non-oscillatory local linear convergence at exact rate $rho = lambda_"max" (Sigma_(w^*)) / h^2 < 1$.
]

#proof[
  Viewing $hat(f)_(h, K) (y) = sum_(i=1)^n (1/n) cal(N)(y | x_i, h^2 I_d)$ as a uniform Gaussian mixture evaluated at a single variable $y$, the E-step posterior responsibilities are precisely $w_i (y) = exp(-norm(y - x_i)_2^2 / (2 h^2)) / sum_j exp(-norm(y - x_j)_2^2 / (2 h^2))$, and the M-step maximizes $-1/(2 h^2) sum_i w_i (y) norm(z - x_i)_2^2$, uniquely yielding $z = sum_i w_i (y) x_i = T_h (y)$.

  For the local Jacobian, denote $phi_i (y) = exp(-norm(y - x_i)_2^2 / (2 h^2))$. Direct differentiation gives $nabla hat(f)_h (y) = (hat(f)_h (y) / h^2) (T_h (y) - y)$, hence $T_h (y) = y + h^2 nabla hat(f)_h (y) / hat(f)_h (y)$. Differentiating once more:
  $ J_(T_h) (y) = I_d + h^2 ( (nabla^2 hat(f)_h (y)) / (hat(f)_h (y)) - (nabla hat(f)_h (y) nabla hat(f)_h (y)^T) / (hat(f)_h (y)^2) ). $
  Meanwhile, differentiating $nabla phi_i (y) = h^(-2) phi_i (y) (x_i - y)$ yields
  $ nabla^2 hat(f)_h (y) = (hat(f)_h (y)) / h^4 [ sum_(i=1)^n w_i (y) (x_i - y)(x_i - y)^T - h^2 I_d ]. $
  Evaluating at $y^*$ where $nabla hat(f)_h (y^*) = bold(0)$ cancels the outer-product term, substituting $nabla^2 hat(f)_h (y^*)$ yields
  $ J_(T_h) (y^*) = I_d + 1 / h^2 (Sigma_(w^*) - h^2 I_d) = 1 / h^2 Sigma_(w^*). $
  Since $Sigma_(w^*) succ.eq 0$, all eigenvalues of $J_(T_h) (y^*)$ are non-negative (ruling out sign-flipping oscillations). Because $y^*$ is a strict local maximum, $H^* = (hat(f)_h (y^*) / h^4)(Sigma_(w^*) - h^2 I_d) prec 0$, which holds if and only if $lambda_"max" (Sigma_(w^*)) < h^2$, hence $rho(J_(T_h) (y^*)) = lambda_"max" (Sigma_(w^*)) / h^2 in [0, 1)$.
]

#subsec("C", "Bandwidth Phase Transition and Mode Collapse")

Bandwidth selection $h$ is the single most critical hyperparameter of Mean Shift, governing a sharp topological bifurcation in the number of modes.

#theorem(title: "Exact Mode-Merging Bifurcation in Bimodal Data", origin: "Original Counterexample")[
  Let $cal(X) subset.eq bb(R)^d$ consist of two balanced point clusters of equal size $n / 2$ concentrated at $-a e_1$ and $+a e_1$ ($a > 0$), estimated with the Gaussian kernel of bandwidth $h > 0$. Then:
  1. For all $h >= a$, the midpoint $y = bold(0)$ is the *unique* critical point and global mode of $hat(f)_(h, K)$; Mean Shift collapses both clusters into one.
  2. For all $0 < h < a$, $y = bold(0)$ becomes a strict saddle (local minimum along $e_1$, maximum transversely for $d >= 2$), and two symmetric strict modes $plus.minus y^*(h) e_1$ emerge via a supercritical pitchfork bifurcation satisfying $y^* = a tanh(a y^* / h^2)$.
]

#proof[
  Write $y = t e_1 + y_perp$ with $y_perp perp e_1$. Because $norm(y plus.minus a e_1)_2^2 = (t plus.minus a)^2 + norm(y_perp)_2^2$, the KDE factorizes as $hat(f)_h (t e_1 + y_perp) = exp(-norm(y_perp)_2^2 / (2 h^2)) tilde(f)_h (t)$, where
  $ tilde(f)_h (t) = c_(k, d) / h^d exp(-(t^2 + a^2) / (2 h^2)) cosh(a t / h^2). $
  Any critical point must satisfy $y_perp = bold(0)$ and $tilde(f)'_h (t) = 0$. Computing the longitudinal derivative yields
  $ tilde(f)'_h (t) = (tilde(f)_h (t)) / h^2 [ a tanh(a t / h^2) - t ] =: (tilde(f)_h (t)) / h^2 psi(t). $
  Note that $psi(0) = 0$ and
  $ psi'(t) = a^2 / (h^2 cosh^2 (a t / h^2)) - 1. $
  Since $cosh^2 (u) > 1$ for all $u != 0$:
  - If $h >= a$, then $a^2 / h^2 <= 1$, so $psi'(t) < 0$ for every $t != 0$. Hence $psi(t) > 0$ for $t < 0$ and $psi(t) < 0$ for $t > 0$, making $t = 0$ the unique global maximum of $tilde(f)_h$ and collapsing all trajectories to $bold(0)$.
  - If $0 < h < a$, then $psi'(0) = a^2 / h^2 - 1 > 0$, whereas $lim_(t -> +infinity) psi(t) = -infinity$. Moreover, for $t > 0$, $psi''(t) = -2 (a^3 / h^4) tanh(a t / h^2) / cosh^2 (a t / h^2) < 0$, so $psi$ is strictly concave on $(0, infinity)$ and admits a unique positive root $y^* in (0, a)$ solving $y^* = a tanh(a y^* / h^2)$ (and symmetrically $-y^*$). At $t = 0$, $tilde(f)''_h (0) = (tilde(f)_h (0) / h^2)(a^2 / h^2 - 1) > 0$, confirming a strict local minimum along $e_1$.
]

#subsec("D", "Isometry Invariance and Outlier Localization")

#proposition(title: "Rigid Invariance and Super-Exponential Outlier Immunity", origin: "Original Proof")[
  Let $T(x) = Q x + b$ be any rigid isometry ($Q in O(d)$, $b in bb(R)^d$). Then:
  1. Mean Shift trajectories commute identically with $T$ and isotropic rescaling $(x, h) |-> (alpha x, alpha h)$, but are *not* invariant under anisotropic diagonal scaling $S(x) = "diag"(s_1, dots, s_d) x$.
  2. For any compactly supported kernel ($g(r) = 0$ for $r > 1$), moving an adversarial outlier $x_1$ to distance $R > h$ from $cal(X) backslash {x_1}$ induces *exactly zero* displacement on all remaining modes. For the Gaussian kernel, perturbing $x_1$ to distance $R = norm(x_1 - y^*)_2$ from a mode $y^*$ with unnormalized clean weight $W_0 = sum_(i=2)^n exp(-norm(x_i - y^*)_2^2 / (2 h^2))$ perturbs the shift at $y^*$ by at most
     $ norm(T_h^"corrupt" (y^*) - y^*)_2 <= R / W_0 exp(-R^2 / (2 h^2)) -> 0 quad (R -> infinity). $
]

#proof[
  (1) Because $norm(Q y + b - (Q x_i + b))_2 = norm(y - x_i)_2$, all weights $g_i$ are invariant under $T$, and convex combinations commute with affine maps: $sum_i alpha_i (Q x_i + b) = Q (sum_i alpha_i x_i) + b$. Anisotropic scaling distorts spherical kernel level sets into ellipsoids, altering $g_i$ and mode locations.

  (2) If $g(r) = 0$ for $r > 1$ and $"dist"(x_1, cal(X) backslash {x_1}) > h$, then $g_1 (y) = 0$ on every $h$-ball intersecting $cal(X) backslash {x_1}$, leaving all non-outlier trajectories identical (only $x_1$ forms a singleton mode). For the Gaussian kernel, since $sum_(i=2)^n g_i (y^*) (x_i - y^*) = bold(0)$ at the clean mode $y^*$, adding $x_1$ yields
  $ T_h^"corrupt" (y^*) - y^* = (g_1 (y^*) (x_1 - y^*)) / (W_0 + g_1 (y^*)), $
  whose norm is bounded above by $(g_1 (y^*) / W_0) R = (R / W_0) exp(-R^2 / (2 h^2))$, decaying super-exponentially to $0$ as $R -> infinity$ (unlike K-Means' breakdown point $1/n$ @lloyd1982least).
]

#subsec("E", "Statistical Consistency of Gradient Lines and Basins")

#theorem(title: "Uniform Consistency of Empirical Mean Shift Flows", origin: "Classical — Cited")[
  Let $X_1, dots, X_n$ be i.i.d. draws from a compactly supported, twice continuously differentiable Morse density $f$ on $bb(R)^d$ with bounded third derivatives, and let $K$ be a smooth radially symmetric kernel (e.g., Gaussian). Under a bandwidth schedule $h_n -> 0$ satisfying $n h_n^(d+4) / log n -> infinity$, the empirical density modes converge to the true modes of $f$ at rate $O_p (h_n^2 + (n h_n^(d+2))^(-1/2))$, and the discrete Mean Shift trajectories initialized on any compact set disjoint from the stable manifolds of saddle points converge uniformly to the true integral curves of $nabla f$ @ariascastro2016meanshift @cheng1995mean.
]

#proof(label: "Proof Synthesis (Cited)")[
  Arias-Castro, Mason, and Pelletier @ariascastro2016meanshift establish via empirical-process concentration inequalities that $n h_n^(d+4) / log n -> infinity$ guarantees uniform almost-sure convergence of $hat(f)_h$, $nabla hat(f)_h$, and $nabla^2 hat(f)_h$ to $(f, nabla f, nabla^2 f)$. Non-degeneracy of the Hessian at Morse critical points ensures via the implicit function theorem that each empirical mode lies within $O_p (h_n^2 + (n h_n^(d+2))^(-1/2))$ of a unique population mode. Combining Lipschitz stability of normalized gradient flows (Grönwall's inequality) with a step-size discretization bound showing $T_h$ tracks the Euler polygon of $nabla hat(f)_h / hat(f)_h$ proves uniform convergence of discrete Mean Shift paths to population gradient ascent lines. Notably, the exponent $d + 4$ exposes the severe nonparametric curse of dimensionality.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact batch Mean Shift evaluates $O(n^2)$ pairwise kernel interactions per iteration across $d$ dimensions for $T$ iterations, costing $O(n^2 d T)$ arithmetic operations and only $O(n d)$ working memory (storing data and active trajectories without materializing an $n times n$ distance matrix). For compactly supported kernels (Epanechnikov, uniform) in low dimension ($d <= 8$), spatial trees (KD-trees, ball trees) prune range queries to $O(n d T log n)$; for Gaussian kernels, Fast Gauss Transforms, dual-tree methods, and grid-binned Mean Shift++ (Alg. 55) @jang2021meanshiftpp accelerate iteration cost toward $O(n d T)$. Trajectories that enter an already visited $epsilon$-neighborhood can also be short-circuited, substantially lowering empirical step counts.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Mean Shift],
    [Time / Space Complexity], [$O(n^2 d T)$ exact; $O(n d T log n)$ with compact-kernel spatial trees; $O(n d)$ memory.],
    [Inductive Bias], [Arbitrary non-convex density basins around smooth KDE modes; $hat(k)$ chosen automatically.],
    [Determinism / Seeding], [Fully deterministic given $(h, k, tau)$; zero random seeding required.],
    [Invariances], [Strictly invariant under $Q in O(d)$, translation, and $(x, h)$ scaling; sensitive to anisotropic scaling.],
    [Outlier Robustness], [High local immunity: exact zero shift for compact kernels, $O(R exp(-R^2 / (2 h^2)))$ for Gaussian.],
    [Recommended Upgrades], [Blurring Mean Shift (Alg. 54) for cubic local speed; Mean Shift++ (Alg. 55) for grid-linear scale; Quickshift (Alg. 41) for arbitrary metrics.],
  ),
  caption: [Theoretical guarantees and operational profile of Mean Shift (Alg. 34).],
)

#sec("IV", "Engineering and Scientific Applications")

Mean Shift excels in low-to-moderate dimensional continuous spaces ($d <= 10$) where cluster count is unknown, shapes are non-convex, and localized outliers must not distort cluster centers.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Discontinuity-Preserving Image Segmentation & Filtering:* Introduced in joint spatial-range space $(x, y, L^*, u^*, v^*) in bb(R)^5$ with product kernel $K_(h_s, h_r)$ @comaniciu2002mean, Mean Shift smooths homogeneous texture regions while preserving sharp object boundaries because gradient ascent never crosses low-density color discontinuities—unlike isotropic Gaussian blur or convex K-Means/SLIC (Alg. 166) @achanta2012slic.
2. *Real-Time Visual Object Tracking:* Modeling target appearance via $m$-bin spatial-weighted color histograms $q$, maximizing the Bhattacharyya similarity $rho(p(y), q) = sum_u sqrt(p_u (y) q_u)$ via Taylor linearization reduces directly to a Mean Shift iteration in image coordinates, achieving millisecond tracking with guaranteed monotonic similarity ascent @comaniciu2002mean.
3. *LiDAR Point-Cloud & Forestry Canopy Segmentation:* Individual tree crowns in airborne 3D LiDAR point clouds form local height-density maxima over irregular ground topologies. Mean Shift's automatic mode discovery extracts variable crown counts without forcing spherical partitions, while isolated bird/sensor returns remain harmless singleton modes (Proposition 2).
4. *Flow Cytometry & Low-Dimensional Cellular Gating:* Fluorescent marker panels ($d approx 3$–$8$) exhibit curved, asymmetric cell subpopulations with unknown subtype counts and debris outliers. Nonparametric KDE mode-seeking delineates skewed immunological populations cleanly without Gaussian covariance assumptions.
5. *Geospatial Hotspot & Epidemiological Peak Detection:* On 2D planar event coordinates (urban incidents, seismic swarms, disease outbreaks), bandwidth $h$ directly encodes physical spatial resolution (meters/kilometers), identifying statistically coherent density peaks invariant to arbitrary map rotation.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Mean Shift degrades sharply under: (a) *high dimensionality ($d > 15$)*, where KDE convergence deteriorates as $n^(-1/(d+4))$ (Theorem 3) and pairwise distance concentration @beyer1999nearest forces either total fragmentation ($hat(k) approx n$) or single-mode collapse (Theorem 2); (b) *strongly varying cluster densities*, since a single global bandwidth $h$ simultaneously over-smooths tight modes and fragments sparse clusters; (c) *flat plateaus and saddle ridges*, where slow gradient decay inflates iteration count $T$; and (d) *non-Euclidean / categorical spaces* where weighted Euclidean means are undefined (resolved by Medoidshift Alg. 42 @sheikh2007medoidshift or Quickshift Alg. 41 @vedaldi2008quickshift).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit / Prerequisite], [Caution / Alternative],
    [Image Segmentation], [Joint spatial-range edge preservation], [Mean Shift++ (Alg. 55) for HD video],
    [Visual Tracking], [Monotonic Bhattacharyya mode climb], [Kalman/particle filter for occlusion],
    [3D LiDAR / GIS], [Unknown $hat(k)$; physical bandwidth $h$], [HDBSCAN (Alg. 32) for varied density],
    [Flow Cytometry], [Skewed non-Gaussian basins; debris robust], [DPC (Alg. 51) if $d > 10$],
    [High-D Embeddings], [Poor fit: $n h^(d+4)$ KDE breakdown], [SNN (Alg. 36) or Spherical (Alg. 12)],
  ),
  caption: [Engineering domain suitability and failure-mode guide for Mean Shift.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Mean Shift anchors the continuous mode-seeking branch of Category II and connects algebraically to density-attractor and graph-based algorithms across the catalog.

#proposition(title: "Exact Equivalence with DENCLUE 2.0 and Blurring Limit", origin: "Original Unified Proof")[
  Under the Gaussian kernel $K_h (u) prop exp(-norm(u)_2^2 / (2 h^2))$, the adaptive-step hill-climbing update of DENCLUE 2.0 (Alg. 43) @hinneburg2007denclue2 is algebraically identical to Mean Shift (Alg. 34). Conversely, replacing static data anchors $x_j$ in Algorithm 1 with evolving iterates $y_j^((t))$ yields Gaussian Blurring Mean Shift (Alg. 54) @carreira2006fast, which accelerates local convergence from linear (Proposition 1) to cubic while asymptotically collapsing all points to a single centroid unless stopped early.
]

#proof[
  In DENCLUE 2.0 @hinneburg2007denclue2, setting the partial derivative of the auxiliary EM lower bound $Q(z | y) = sum_i w_i (y) log cal(N)(z | x_i, sigma^2 I_d)$ to zero gives $z = sum_i w_i (y) x_i$, which coincides term-for-term with $T_sigma (y)$ in Proposition 1. In Blurring Mean Shift @carreira2006fast, replacing $x_j$ by $y_j^((t))$ replaces the fixed local covariance $Sigma_(w^*)$ by a shrinking covariance $Sigma^((t+1)) approx h^(-4) (Sigma^((t)))^3$ for a Gaussian cloud, yielding cubic contraction of intra-cluster spread while still shifting separated clusters toward the global mean unless truncated by an explicit stopping criterion.
]

Beyond DENCLUE (Alg. 35) @hinneburg1998denclue and its adaptive successor (Alg. 43), restricting mode-seeking shifts to observed data points yields Medoidshift (Alg. 42) @sheikh2007medoidshift and Quickshift (Alg. 41) @vedaldi2008quickshift, whereas replacing continuous ascent with one-hop links to nearest higher-density points yields Density Peaks Clustering (DPC, Alg. 51) @rodriguez2014dpc. Unlike DBSCAN (Alg. 31) @ester1996dbscan, which extracts connected superlevel sets at a fixed density cut, Mean Shift partitions space by gradient-flow basins—splitting multimodal ridges even when connected above a high density threshold.

In empirical evaluation on the ClusterOverview 4-dataset benchmark, Mean Shift attains a mean runtime of $7.649$ ms and excels on low-dimensional non-convex geometry, achieving $"Sil" = 0.787$ on Multi-Shape ($"ARI"$ N/A) and solid cluster recovery on Iris ($"Sil" = 0.586$, $"ARI" = 0.567$). Conversely, on Customers ($"Sil"$ N/A) and High-D Subspace+Noise ($"Sil"$ N/A, $"ARI"$ N/A), the fixed global bandwidth collapses points into a degenerate single-basin regime where Silhouette is undefined—directly corroborating Theorem 2's bandwidth bifurcation and Theorem 3's $d + 4$ dimensionality barrier, motivating multi-scale density or rank-based SNN methods (Alg. 36) in high-dimensional settings.
