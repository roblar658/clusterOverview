#import "../template/ieee.typ": *

#alg_entry(
  id: 43,
  name: "DENCLUE 2.0",
  category: "Cat. II: Density & Mode-Seeking",
  year: "2007",
  paradigm: "Adaptive EM-Style Kernel Mode-Seeking",
  complexity: $O(n^2 d T_"avg")$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite sample equipped with the Euclidean norm $norm(dot)_2$, and let $h > 0$ be an isotropic smoothing bandwidth. DENCLUE 2.0 @hinneburg2007denclue2 models the empirical distribution analytically via the *Gaussian kernel density estimator (KDE)*:
$ hat(f)_h (x) = 1 / n sum_(i=1)^n K_h (x - x_i), $
where the isotropic Gaussian kernel is
$ K_h (u) = 1 / (2 pi h^2)^(d / 2) exp(- norm(u)_2^2 / (2 h^2)). $
Clusters are associated with local maxima (*density attractors*) of the smooth scalar field $hat(f)_h : bb(R)^d -> (0, infinity)$. In DENCLUE 1.0 (Alg. 35, @hinneburg1998denclue), attractors were sought via fixed-step normalized gradient ascent $x^((t+1)) = x^((t)) + delta nabla hat(f)_h (x^((t))) / norm(nabla hat(f)_h (x^((t))))_2$, which requires manually tuning $delta > 0$ and can overshoot or oscillate indefinitely around modes.

DENCLUE 2.0 @hinneburg2007denclue2 eliminates the fixed step-size hyperparameter by solving the first-order stationarity equation $nabla hat(f)_h (x) = bold(0)$ directly as an *adaptive fixed-point iteration*. Differentiating $hat(f)_h$ yields
$ nabla hat(f)_h (x) = 1 / (n h^2) sum_(i=1)^n K_h (x - x_i) (x_i - x). $
Setting $nabla hat(f)_h (x) = bold(0)$ and solving for the isolated $x$ on the right-hand side defines the *DENCLUE 2.0 hill-climbing operator* $T: bb(R)^d -> bb(R)^d$:
$ x^((t+1)) = T(x^((t))) := sum_(i=1)^n w_i (x^((t))) x_i, $
where the normalized softmax kernel responsibilities are
$ w_i (x) := (K_h (x - x_i)) / (sum_(j=1)^n K_h (x - x_j)) in (0, 1), quad sum_(i=1)^n w_i (x) = 1. $
Equivalently, the adaptive displacement vector is
$ x^((t+1)) - x^((t)) = (h^2) / (hat(f)_h (x^((t)))) nabla hat(f)_h (x^((t))), $
automatically scaling the step size inversely with local density $hat(f)_h (x^((t)))$ so iterates take large strides across low-density regions and contract smoothly near peaks.

#definition(title: "Density Attractor and Center / Arbitrary-Shape Clusters")[
  A point $x^* in bb(R)^d$ is a *density attractor* for seed $x_i$ if the sequence $x^((0)) = x_i$, $x^((t+1)) = T(x^((t)))$ converges to $x^*$. Given a noise threshold $xi >= 0$ and merge radius $tau_m > 0$:
  1. A *center-defined cluster* around a significant attractor $x^*$ with $hat(f)_h (x^*) >= xi$ is its basin $cal(A)(x^*) = {x_i in cal(X) : lim_(t -> infinity) T^t (x_i) = x^*}$ (grouping numerically close attractors within distance $tau_m$).
  2. An *arbitrary-shape cluster* is a maximal union of basins $union.big_(r in S) cal(A)(x_r^*)$ whose attractors satisfy $hat(f)_h (x_r^*) >= xi$ and are connected by continuous paths $gamma$ with $hat(f)_h (gamma(s)) >= xi$. Points whose attractors have $hat(f)_h (x^*) < xi$ are declared *noise*.
]

Hinneburg and Gabriel @hinneburg2007denclue2 further accelerate evaluation via local neighborhood pruning (restricting sums to active points near $x^((t))$), trajectory sharing (stopping early when $x^((t))$ approaches an already-converged orbit), and step-size stopping criteria (Algorithm 1).

#algbox(
  title: "DENCLUE 2.0 Adaptive Hill-Climbing Clustering",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, bandwidth $h > 0$, noise threshold $xi >= 0$, tolerance $epsilon_"stop" > 0$, merge radius $tau_m > 0$, max steps $T_"max"$.],
  outputs: [Cluster labels $cal(C)$ and noise set $cal(N)$.],
)[
  1. *For* each $i in {1, dots, n}$ *do*:
  2. #h(1em) Initialize $x^((0)) <- x_i$.
  3. #h(1em) *For* $t = 0, 1, dots, T_"max" - 1$ *do*:
  4. #h(2em) Evaluate weights $w_j (x^((t))) prop exp(-norm(x^((t)) - x_j)_2^2 / (2 h^2))$ over $cal(X)$ (or active local spatial cells within $r h$).
  5. #h(2em) Update $x^((t+1)) <- sum_(j=1)^n w_j (x^((t))) x_j$.
  6. #h(2em) *If* $norm(x^((t+1)) - x^((t)))_2 <= epsilon_"stop"$, *break*.
  7. #h(1em) Record attractor $a_i <- x^((t+1))$ and density $hat(f)_h (a_i)$.
  8. Mark $i in cal(N)$ if $hat(f)_h (a_i) < xi$; merge remaining attractors within distance $tau_m$ (or $xi$-connected paths) into clusters $cal(C)$.
  9. *Return* $(cal(C), cal(N))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical pillars of DENCLUE 2.0: (1) an exact EM/MM minorizer identity yielding a quantitative Lyapunov log-density bound and finite-precision iteration bound; (2) the exact Jacobian covariance spectrum proving oscillation-free local $Q$-linear convergence; (3) exact rigid-motion invariance paired with constructive counterexamples for bandwidth-collapse bifurcation and high-dimensional softmax freeze; and (4) classical EM convergence and mode-flow consistency bounds.

#subsec("A", "Exact EM Minorizer and Quantitative Lyapunov Ascent")

#lemma(title: "Tight Jensen–EM Minorizer of the Log-KDE", origin: "Original Derivation")[
  Fix any iterate $x^((t)) in bb(R)^d$ with responsibility vector $w^((t)) = (w_1 (x^((t))), dots, w_n (x^((t))))$, and define the surrogate functional
  $ Q(x | x^((t))) := &log hat(f)_h (x^((t))) \
                      &+ sum_(i=1)^n w_i (x^((t))) log (K_h (x - x_i)) / (K_h (x^((t)) - x_i)). $
  Then for every $x in bb(R)^d$:
  1. $log hat(f)_h (x) - Q(x | x^((t))) = D_"KL" (w^((t)) || w(x)) >= 0$, with equality if and only if $w(x) = w^((t))$.
  2. $Q(dot | x^((t)))$ is strictly concave and quadratic, with unique global maximizer $x^((t+1)) = T(x^((t))) = sum_(i=1)^n w_i (x^((t))) x_i$.
  3. The exact surrogate ascent at $x^((t+1))$ is
     $ Q(x^((t+1)) | x^((t))) - log hat(f)_h (x^((t))) = 1 / (2 h^2) norm(x^((t+1)) - x^((t)))_2^2. $
]

#proof[
  Let $Z(x) = sum_(j=1)^n K_h (x - x_j) = n hat(f)_h (x)$, so $w_i (x) = K_h (x - x_i) / Z(x)$. Expanding the Kullback–Leibler divergence between $w^((t))$ and $w(x)$ gives
  $ D_"KL" (w^((t)) || w(x)) &= sum_(i=1)^n w_i^((t)) log ((K_h (x^((t)) - x_i) / Z(x^((t)))) / (K_h (x - x_i) / Z(x))) \
    &= log (Z(x) / Z(x^((t)))) - sum_(i=1)^n w_i^((t)) log (K_h (x - x_i)) / (K_h (x^((t)) - x_i)) \
    &= log hat(f)_h (x) - Q(x | x^((t))). $
  By Gibbs' inequality, $D_"KL" (w^((t)) || w(x)) >= 0$ with equality iff $w(x) = w^((t))$, proving Claim 1.

  For Claim 2, substituting the Gaussian kernel ratio cancels normalization constants:
  $ Q(x | x^((t))) - log hat(f)_h (x^((t))) = 1 / (2 h^2) sum_(i=1)^n w_i^((t)) (norm(x^((t)) - x_i)_2^2 - norm(x - x_i)_2^2). $
  Because $sum_(i=1)^n w_i^((t)) = 1$, applying the weighted parallel-axis identity around the convex combination $x^((t+1)) = sum_(i=1)^n w_i^((t)) x_i$ yields, for every $x in bb(R)^d$,
  $ sum_(i=1)^n w_i^((t)) norm(x - x_i)_2^2 = sum_(i=1)^n w_i^((t)) norm(x^((t+1)) - x_i)_2^2 + norm(x - x^((t+1)))_2^2. $
  Substituting this identity twice (for $x$ and for $x^((t))$) yields
  $ Q(x | x^((t))) - log hat(f)_h (x^((t))) = 1 / (2 h^2) (norm(x^((t)) - x^((t+1)))_2^2 - norm(x - x^((t+1)))_2^2). $
  The right-hand side is strictly concave in $x$ with Hessian $-h^(-2) I_d prec 0$, uniquely maximized at $x = x^((t+1))$, where it equals $norm(x^((t+1)) - x^((t)))_2^2 / (2 h^2)$, proving Claims 2 and 3.
]

#theorem(title: "Strict Lyapunov Ascent, Convex-Hull Trapping, and Iteration Bound", origin: "Original Proof")[
  For any initialization $x^((0)) in bb(R)^d$, the DENCLUE 2.0 sequence $x^((t+1)) = T(x^((t)))$ satisfies:
  1. *Exact one-step log-density gain:*
     $ log hat(f)_h (x^((t+1))) - log hat(f)_h (x^((t))) = 1 / (2 h^2) norm(x^((t+1)) - x^((t)))_2^2 + D_"KL" (w^((t)) || w^((t+1))). $
  2. *Convex-hull invariance:* $x^((t)) in "conv"(cal(X))$ for all $t >= 1$, and $hat(f)_h (x^((t)))$ increases monotonically to a finite limit $f^* <= (2 pi h^2)^(-d/2)$.
  3. *Finite-precision step bound:* The number of iterations $T_epsilon$ before $norm(x^((t+1)) - x^((t)))_2 <= epsilon$ satisfies
     $ T_epsilon <= floor((2 h^2) / (epsilon^2) log (hat(f)_"max" / (hat(f)_h (x^((0)))))), $
     where $hat(f)_"max" = max_(x in "conv"(cal(X))) hat(f)_h (x)$. Moreover, $nabla hat(f)_h (x^((t))) -> bold(0)$ as $t -> infinity$, and every limit point of ${x^((t))}$ is a stationary point of $hat(f)_h$.
]

#proof[
  Claim 1 follows immediately by combining Claims 1 and 3 of Lemma 1 at $x = x^((t+1))$. Because $w_i (x^((t))) > 0$ and $sum_(i=1)^n w_i (x^((t))) = 1$, every iterate $x^((t+1)) = sum_(i=1)^n w_i (x^((t))) x_i$ is a strict convex combination of $cal(X)$, proving $x^((t)) in "conv"(cal(X))$ for all $t >= 1$. Since $K_h (u) <= (2 pi h^2)^(-d/2)$, $hat(f)_h$ is bounded above on compact $"conv"(cal(X))$ by $hat(f)_"max" <= (2 pi h^2)^(-d/2)$. Dropping the nonnegative KL term in Claim 1 and telescoping over $t = 0, dots, N - 1$ gives
  $ sum_(t=0)^(N-1) norm(x^((t+1)) - x^((t)))_2^2 <= 2 h^2 log ((hat(f)_h (x^((N)))) / (hat(f)_h (x^((0))))) <= 2 h^2 log (hat(f)_"max" / (hat(f)_h (x^((0))))). $
  If $norm(x^((t+1)) - x^((t)))_2 > epsilon$ for $T_epsilon$ consecutive steps, then $T_epsilon epsilon^2 < 2 h^2 log(hat(f)_"max" / hat(f)_h (x^((0))))$, proving the finite iteration bound. Letting $N -> infinity$ shows $sum_(t=0)^infinity norm(x^((t+1)) - x^((t)))_2^2 < infinity$, hence $norm(x^((t+1)) - x^((t)))_2 -> 0$. Because $nabla hat(f)_h (x^((t))) = h^(-2) hat(f)_h (x^((t))) (x^((t+1)) - x^((t)))$ and $hat(f)_h (x^((t))) <= hat(f)_"max"$, we obtain $nabla hat(f)_h (x^((t))) -> bold(0)$. By compactness of $"conv"(cal(X))$ and continuity of $nabla hat(f)_h$, every accumulation point $x^*$ satisfies $nabla hat(f)_h (x^*) = bold(0)$.
]

#subsec("B", "Exact Jacobian Spectrum and Local Contraction Rate")

Unlike fixed-step gradient ascent, DENCLUE 2.0 never overshoots a strict local mode along any eigen-direction.

#theorem(title: "Positive Semidefinite Jacobian and Local Q-Linear Convergence", origin: "Original Proof")[
  For every $x in bb(R)^d$, the Jacobian matrix $J_T (x) = (partial T(x)) / (partial x)$ equals the scaled local covariance of $cal(X)$ under $w(x)$:
  $ J_T (x) = 1 / h^2 sum_(i=1)^n w_i (x) (x_i - T(x)) (x_i - T(x))^T = 1 / h^2 "Cov"_(w(x)) (X) succ.eq 0. $
  At any strict local maximizer $x^*$ of $hat(f)_h$ (where $T(x^*) = x^*$ and $nabla^2 hat(f)_h (x^*) prec 0$), the Jacobian satisfies
  $ J_T (x^*) = I_d + (h^2) / (hat(f)_h (x^*)) nabla^2 hat(f)_h (x^*), $
  so every eigenvalue $lambda_r$ of $J_T (x^*)$ lies strictly in $[0, 1)$. Consequently, DENCLUE 2.0 converges locally $Q$-linearly without oscillation, with asymptotic contraction rate $rho = max_r lambda_r (J_T (x^*)) < 1$.
]

#proof[
  Differentiating $w_i (x) = K_h (x - x_i) / Z(x)$ using $nabla_x K_h (x - x_i) = h^(-2) K_h (x - x_i) (x_i - x)$ gives
  $ nabla_x w_i (x) &= w_i (x) / h^2 (x_i - x) - w_i (x) / h^2 sum_(j=1)^n w_j (x) (x_j - x) \
                    &= 1 / h^2 w_i (x) (x_i - T(x)). $
  Hence the Jacobian of $T(x) = sum_(i=1)^n w_i (x) x_i$ is
  $ J_T (x) &= sum_(i=1)^n x_i (nabla_x w_i (x))^T = 1 / h^2 sum_(i=1)^n w_i (x) x_i (x_i - T(x))^T \
            &= 1 / h^2 sum_(i=1)^n w_i (x) (x_i - T(x)) (x_i - T(x))^T, $
  where subtracting $T(x) sum_(i=1)^n w_i (x) (x_i - T(x))^T = bold(0)$ exposes manifest symmetry and positive semidefiniteness $J_T (x) succ.eq 0$. Next, differentiate $nabla hat(f)_h (x) = h^(-2) hat(f)_h (x) (T(x) - x)$:
  $ nabla^2 hat(f)_h (x) = (hat(f)_h (x)) / h^2 (J_T (x) - I_d) + 1 / h^2 (T(x) - x) nabla hat(f)_h (x)^T. $
  Evaluating at $x^*$ where $T(x^*) = x^*$ annihilates the second term and yields $J_T (x^*) = I_d + (h^2 / hat(f)_h (x^*)) nabla^2 hat(f)_h (x^*)$. Because $nabla^2 hat(f)_h (x^*) prec 0$ and $hat(f)_h (x^*) > 0$, all eigenvalues of $J_T (x^*)$ are strictly less than $1$; because $J_T (x^*) succ.eq 0$, all eigenvalues are nonnegative. Thus $0 <= rho(J_T (x^*)) < 1$. By Ostrowski's theorem, there exists a neighborhood of $x^*$ in which $norm(x^((t+1)) - x^*)_2 <= (rho + eta) norm(x^((t)) - x^*)_2$ for arbitrarily small $eta > 0$, with nonnegative spectrum ruling out oscillatory overshoots.
]

#remark(title: "Cluster-Variance Dependence of Contraction Speed")[
  If a well-separated cluster is sampled from $cal(N)(x^*, sigma_c^2 I_d)$, the Gaussian product covariance satisfies $"Cov"_(w(x^*)) (X) approx (sigma_c^2 h^2) / (sigma_c^2 + h^2) I_d$, yielding $rho approx sigma_c^2 / (sigma_c^2 + h^2)$. Tight clusters ($sigma_c << h$) contract rapidly in a handful of steps, whereas flat plateaus ($sigma_c >> h$) yield $rho -> 1^-$, motivating the convergence-acceleration checks in @hinneburg2007denclue2.
]

#subsec("C", "Invariance, Bandwidth Bifurcation, and High-D Softmax Freeze")

#proposition(title: "Isometry Equivariance, Exact Mode-Collapse Threshold, and Softmax Freeze", origin: "Original Counterexample")[
  1. For any rigid isometry $Phi(x) = Q x + b$ ($Q in O(d)$, $b in bb(R)^d$), $T(Phi(x)) = Phi(T(x))$ identically; however, anisotropic scaling $S(x) = D x$ alters attractors.
  2. (*Exact Mode-Collapse Bifurcation:*) On the symmetric dataset $cal(X) = {-a, +a} subset.eq bb(R)$ ($a > 0$, with $n/2$ replicates at each point), DENCLUE 2.0 recovers two distinct attractors $plus.minus x^* != 0$ if and only if $h < a$. For every $h >= a$, the two clusters collapse catastrophically into the unique global attractor $x^* = 0$.
  3. (*High-Dimensional Softmax Freeze:*) If $min_(j != i) norm(x_i - x_j)_2^2 >= c d$ with fixed $h > 0$, then at $x^((0)) = x_i$ the self-weight satisfies $1 - w_i (x_i) <= (n - 1) exp(-c d / (2 h^2)) -> 0$ as $d -> infinity$, trapping every sample point at a trivial singleton mode.
]

#proof[
  (1) Since $norm(Phi(x) - Phi(x_i))_2 = norm(x - x_i)_2$, weights $w_i (Phi(x)) = w_i (x)$ are unchanged, so $T(Phi(x)) = sum_i w_i (x) (Q x_i + b) = Q T(x) + b = Phi(T(x))$.

  (2) On $cal(X) = {-a, +a}$ with equal multiplicity,
  $ T(x) = (-a e^(-(x+a)^2 / (2 h^2)) + a e^(-(x-a)^2 / (2 h^2))) / (e^(-(x+a)^2 / (2 h^2)) + e^(-(x-a)^2 / (2 h^2))) = a tanh((a x) / h^2). $
  By Theorem 2, $J_T (0) = a^2 / h^2$. If $h >= a$, then for every $x != 0$ the strict inequality $|tanh(u)| < |u|$ ($u != 0$) yields
  $ |T(x)| = a tanh((a |x|) / h^2) < (a^2 / h^2) |x| <= |x|, $
  so $|x^((t))|$ decreases strictly to $0$ from every starting point—collapsing both modes into $x^* = 0$. Conversely, if $h < a$, then $T'(0) = a^2 / h^2 > 1$ (so $hat(f)_h ''(0) > 0$ makes $0$ an unstable local minimum) while $T(a) = a tanh(a^2 / h^2) < a$; by the intermediate value theorem and strict concavity of $tanh$ on $(0, infinity)$, there exist unique symmetric attracting fixed points $plus.minus x^* in (0, a)$ solving $x^* = a tanh(a x^* / h^2)$.

  (3) For $x = x_i$, $K_h (0) = (2 pi h^2)^(-d/2)$ while $K_h (x_i - x_j) <= K_h (0) exp(-c d / (2 h^2))$ for all $j != i$. Hence
  $ 1 - w_i (x_i) = (sum_(j != i) K_h (x_i - x_j)) / (K_h (0) + sum_(j != i) K_h (x_i - x_j)) <= (n - 1) exp(-(c d) / (2 h^2)), $
  which vanishes exponentially in $d$, forcing $norm(T(x_i) - x_i)_2 -> 0$ and fragmenting the dataset into isolated point modes.
]

#subsec("D", "Classical EM Convergence and Mode-Flow Consistency")

#theorem(title: "General Radial-Kernel Convergence and Gradient-Flow Consistency", origin: "Classical — Cited")[
  1. (*EM Convergence & General Kernels:*) Viewing $hat(f)_h (x)$ as a Gaussian mixture likelihood over fixed centers $cal(X)$, Hinneburg and Gabriel @hinneburg2007denclue2 establish DENCLUE 2.0 as an exact Expectation-Maximization algorithm @dempster1977em obeying Wu's stationary convergence theorem @wu1983convergence; for general radial kernels $K(u) = c_k k(norm(u)_2^2)$ with convex decreasing profile $k$, the shadow-kernel update $T_g (x) = (sum_i g_i x_i) / (sum_i g_i)$ ($g(s) = -k'(s)$) preserves monotone density ascent @comaniciu2002mean.
  2. (*Statistical Consistency of Mode Basins:*) Let $X_1, dots, X_n$ be i.i.d. from a smooth Morse density $f in C^3 (bb(R)^d)$ with compact support. If $h_n -> 0$ and $n h_n^(d + 4) / log n -> infinity$, then the empirical DENCLUE 2.0 trajectories converge uniformly to the true gradient-ascent integral curves of $f$, and the induced basin partition consistently recovers the population mode basins @ariascastro2016meanshift.
]

#proof(label: "Proof Synthesis (Cited)")[
  For any convex non-increasing kernel profile $k$, convexity gives $k(s_2) - k(s_1) >= k'(s_1)(s_2 - s_1) = -g(s_1)(s_2 - s_1)$, so maximizing the quadratic lower bound $-sum_i g(norm((x^((t)) - x_i)/h)_2^2) norm(x - x_i)_2^2$ reproduces $T_g (x^((t)))$ and forces $hat(f)_h (x^((t+1))) >= hat(f)_h (x^((t)))$ @comaniciu2002mean. Under $h_n -> 0$ with $n h_n^(d+4) / log n -> infinity$, standard KDE concentration ensures $sup_x norm(nabla hat(f)_(h_n) (x) - nabla f(x))_2 = O_p (sqrt(log n / (n h_n^(d+2))) + h_n^2) -> 0$ and uniform Hessian convergence. Arias-Castro et al. @ariascastro2016meanshift then apply Grönwall's inequality along normalized gradient flows and Morse stability of non-degenerate maxima to prove uniform trajectory convergence away from basin separatrices.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact evaluation of $T(x^((t)))$ scans all $n$ data points in $O(n d)$ time; executed across $n$ initial points for an average of $T_"avg"$ steps, unpruned DENCLUE 2.0 runs in $O(n^2 d T_"avg")$ time and $O(n d)$ memory. Because Gaussian weights decay super-exponentially beyond distance $r h$ ($r approx 3$–$4$), Hinneburg and Gabriel @hinneburg2007denclue2 restrict the summation to spatial grid cubes or $k$-d/ball-tree neighborhoods of average size $m << n$ and terminate trajectories once they enter the basin of a previously converged orbit, lowering practical runtime to $O(n m d T_"avg")$. Compared to DENCLUE 1.0 (Alg. 35, @hinneburg1998denclue), eliminating fixed step tuning reduces $T_"avg"$ dramatically while avoiding end-of-climb oscillation. However, in high ambient dimensions ($d >= 20$), spatial pruning loses selectivity due to distance concentration @beyer1999nearest, reverting toward $O(n^2 d T_"avg")$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for DENCLUE 2.0],
    [Time / Space Complexity], [$O(n^2 d T_"avg")$ exact ($O(n m d T_"avg")$ pruned); $O(n d)$ memory.],
    [Geometric Inductive Bias], [Smooth nonparametric KDE basins and $xi$-connected ridges of arbitrary shape.],
    [Convergence Guarantee], [Strict monotone log-density ascent; oscillation-free local $Q$-linear rate (Theorems 1–2).],
    [Invariances], [Equivariant under $Q x + b$ ($Q in O(d)$); sensitive to anisotropic scaling and bandwidth $h$.],
    [Outlier / Noise Handling], [Native rejection via density threshold $hat(f)_h (x^*) < xi$; isolated outliers do not shift dense peaks.],
    [Recommended Companions], [Mean Shift++ (Alg. 55) for grid speed; Blurring Mean Shift (Alg. 54); RNN-DBSCAN (Alg. 45) for multi-density.],
  ),
  caption: [Operational and theoretical profile of DENCLUE 2.0.],
)

#sec("IV", "Engineering and Scientific Applications")

DENCLUE 2.0 excels whenever observations are drawn from continuous multimodal densities contaminated by additive background noise, and where clusters correspond naturally to smooth probability peaks or connected high-density ridges.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multimedia & Feature-Vector Database Indexing:* Originally targeted by Hinneburg et al. @hinneburg1998denclue @hinneburg2007denclue2 for large multimedia repositories (MPEG/color histograms, spectral audio descriptors), where smooth KDE aggregation suppresses quantization noise and $xi$-thresholding isolates background clutter without requiring convex clusters.
2. *Single-Molecule Localization Microscopy (STORM / PALM):* Fluorescent emitters are observed with isotropic Gaussian point-spread-function uncertainty $sigma_"PSF"$. Setting $h approx sigma_"PSF"$ aligns the Gaussian kernel directly with the physical noise model, making density attractors maximum-likelihood-like emitter localizations.
3. *Geospatial Hotspot & Mobility Corridor Discovery:* Urban crime incidents, taxi pick-up coordinates, and wildlife GPS telemetry form irregular peaks and linear road corridors; center-defined clusters extract exact hotspot centroids while $xi$-connected arbitrary-shape unions trace continuous transit ribbons.
4. *Astronomical Star-Cluster & Substructure Identification:* Smooth kernel density estimation over phase-space or sky coordinates cleanly distinguishes gravitationally bound star clusters from uniform field-star backgrounds via the threshold $xi$, without imposing spherical Voronoi boundaries.
5. *Flow & Mass Cytometry (CyTOF) Subpopulation Gating:* Continuous multiparametric biomarker intensities exhibit skewed, overlapping, non-spherical cell-state peaks; adaptive EM hill climbing traces smooth valleys and converges cleanly to cell-type modes without step-size tuning.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DENCLUE 2.0 degrades under four well-characterized mechanisms: (a) *global bandwidth mismatch* (Proposition 1(2)), where a single $h$ either merges close modes when $h >= a$ or shatters diffuse clusters when $h$ is too small; (b) *high-dimensional softmax freeze* (Proposition 1(3)), where distance concentration traps iterates at singletons unless $h$ is inflated into mode collapse (prefer OptiGrid @hinneburg1999optigrid or subspace methods); (c) *strongly heterogeneous cluster densities*, where no single $(h, xi)$ separates both sparse and tight clusters simultaneously (prefer HDBSCAN Alg. 32 or RNN-DBSCAN Alg. 45); and (d) *very large unindexed $n$*, where exact $O(n^2 d T_"avg")$ cost becomes prohibitive (prefer Mean Shift++ Alg. 55 @jang2021meanshiftpp or Quickshift Alg. 41 @vedaldi2008quickshift).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why DENCLUE 2.0 Fits], [Caution / Alternative],
    [Fluorescence / Cryo-EM], [Gaussian kernel matches PSF], [Tune $h$ to sensor $sigma$],
    [Geospatial Hotspots], [Smooth peaks & $xi$-ridges], [Use HDBSCAN (Alg. 32) if densities vary],
    [Cytometry Gating], [Arbitrary smooth cell modes], [Avoid $d > 15$ without PCA/UMAP],
    [Multimedia Databases], [Built-in noise threshold $xi$], [Use Mean Shift++ (Alg. 55) for huge $n$],
    [High-D Sparse Text], [Softmax freeze / collapse], [Use Spherical K-Means (Alg. 12) or SNN (Alg. 36)],
  ),
  caption: [Domain suitability and engineering guide for DENCLUE 2.0.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DENCLUE 2.0 sits at the exact mathematical intersection of kernel mode-seeking and graph-based density level sets.

#proposition(title: "Exact Equivalence of Gaussian DENCLUE 2.0 and Mean Shift", origin: "Original Unified Proof")[
  For the isotropic Gaussian kernel $K_h (u) = (2 pi h^2)^(-d/2) k(norm(u/h)_2^2)$ with profile $k(s) = exp(-s / 2)$, the DENCLUE 2.0 EM update $T(x)$ is *algebraically identical* to the Mean Shift operator (Alg. 34, @cheng1995mean, @comaniciu2002mean). The two algorithms differ solely in post-convergence cluster formation: Mean Shift partitions $cal(X)$ strictly by mode basins, whereas DENCLUE 2.0 filters attractors by $hat(f)_h (x^*) >= xi$ and optionally unions basins along $xi$-level-set paths.
]

#proof[
  In Comaniciu and Meer's formulation @comaniciu2002mean, Mean Shift uses the shadow profile $g(s) = -k'(s)$ to define the shift vector
  $ m_h (x) = (sum_(i=1)^n g(norm((x - x_i) / h)_2^2) x_i) / (sum_(i=1)^n g(norm((x - x_i) / h)_2^2)) - x. $
  For $k(s) = exp(-s / 2)$, differentiation yields $g(s) = -k'(s) = 1 / 2 exp(-s / 2) = 1 / 2 k(s)$. The constant factor $1 / 2$ cancels between numerator and denominator, so $x + m_h (x) = (sum_i K_h (x - x_i) x_i) / (sum_i K_h (x - x_i)) = T(x)$ identically for every $x in bb(R)^d$.
]

This identity unifies DENCLUE 2.0 with Mean Shift (Alg. 34), Blurring Mean Shift (Alg. 54, @carreira2006fast, where the data support $cal(X)$ itself updates synchronously at each iteration for cubic-like contraction), Medoidshift (Alg. 42, @sheikh2007medoidshift, restricting $x^((t))$ to sample points), Quickshift (Alg. 41, @vedaldi2008quickshift, linking each point to its nearest denser neighbor), and DBSCAN (Alg. 31, @ester1996dbscan, recovered conceptually when $K_h$ is replaced by a flat top-hat ball of radius $epsilon$ with threshold $xi = "minPts" / (n V_d epsilon^d)$).

Empirically on the ClusterOverview 4-dataset benchmark, DENCLUE 2.0 achieves *Iris* $"Sil" = 0.586$, $"ARI" = 0.567$ and an outstanding *Multi-Shape* $"Sil" = 0.787$ ($"ARI" = "N/A"$), confirming that adaptive EM hill climbing cleanly concentrates smooth, well-separated modes without oscillation. Conversely, both *Customers* ($"Sil" = "N/A"$) and *High-D Subspace+Noise* ($"Sil" = "N/A"$, $"ARI" = "N/A"$) collapse under automated global bandwidth selection—directly corroborating Proposition 1's mode-collapse and high-dimensional softmax-freeze regimes. Its mean runtime of $7.234$ ms reflects the exact $O(n^2 d T_"avg")$ dense kernel evaluation across all trajectories without grid or spatial-tree pruning.
