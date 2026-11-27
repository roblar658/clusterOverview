#import "../template/ieee.typ": *

#alg_entry(
  id: 54,
  name: "Blurring Mean Shift",
  category: "Cat. II: Density & Mode-Seeking",
  year: "1995",
  paradigm: "Iterative Dataset-Updating Kernel Mode-Seeking",
  complexity: $cal(O)(T n^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Nonparametric mode-seeking clustering identifies clusters with basins of attraction of local maxima of a kernel density estimate (KDE). Given a dataset $X^((0)) = {x_1^((0)), dots, x_n^((0))} subset bb(R)^d$, bandwidth $h > 0$, and a radially symmetric kernel $K(u) = c_(k, d) k(norm(u)_2^2)$ with differentiable, nonincreasing profile $k : [0, infinity) -> [0, infinity)$ and shadow profile $g(r) = -k'(r) >= 0$, the KDE is
$
  hat(p)_h (x) = 1 / (n h^d) sum_(j=1)^n K((x - x_j^((0))) / h).
$
In *standard (non-blurring) Mean Shift* @cheng1995mean @comaniciu2002mean, each query trajectory $y_i^((t))$ ascends $hat(p)_h$ while the data points $X^((0))$ defining the density surface remain frozen. *Blurring Mean Shift* (BMS), formulated by Cheng @cheng1995mean and analyzed in depth for Gaussian kernels (GBMS) by Carreira-Perpiñán @carreira2006fast, instead updates the dataset itself synchronously at every step: the KDE is continuously recomputed on the contracting point cloud $X^((t))$, accelerating local cluster collapse from linear to cubic convergence.

#definition(title: "Blurring Mean Shift Operator and Matrix Dynamics")[
  Let $bold(X)^((t)) in bb(R)^(n times d)$ denote the matrix whose $i$-th row is $(x_i^((t)))^top$ at iteration $t >= 0$, and let $g : [0, infinity) -> [0, infinity)$ be a nonincreasing shadow profile with $g(0) > 0$. Define the symmetric affinity matrix $bold(W)^((t)) in bb(R)^(n times n)$ and degree matrix $bold(D)^((t)) = "diag"(d_1^((t)), dots, d_n^((t)))$ by
  $
    W_(i j)^((t)) &= g((norm(x_i^((t)) - x_j^((t)))_2^2) / h^2), \
    d_i^((t)) &= sum_(j=1)^n W_(i j)^((t)) > 0.
  $
  The *Blurring Mean Shift* dynamical system is the synchronous nonlinear Markov recursion
  $
    bold(X)^((t+1)) = bold(P)^((t)) bold(X)^((t)),
    quad bold(P)^((t)) = (bold(D)^((t)))^(-1) bold(W)^((t)),
  $
  equivalently written point-wise as
  $
    x_i^((t+1)) &= sum_(j=1)^n P_(i j)^((t)) x_j^((t)) \
                &= (sum_(j=1)^n W_(i j)^((t)) x_j^((t))) / (sum_(j=1)^n W_(i j)^((t))).
  $
  For *Gaussian Blurring Mean Shift* (GBMS), $k(r) = exp(-r / 2)$ implies $g(r) = 1 / 2 exp(-r / 2)$, and the constant $1 / 2$ cancels so $W_(i j)^((t)) = exp(-norm(x_i^((t)) - x_j^((t)))_2^2 / (2 h^2))$ @carreira2006fast.
]

Because a strictly positive kernel (such as Gaussian) eventually collapses the entire dataset to a single point as $t -> infinity$, GBMS requires a stopping criterion that halts execution during its *metastable plateau*—after local clusters have contracted to sharp point-masses but before exponentially slow inter-cluster attraction merges them @carreira2006fast. Following Carreira-Perpiñán @carreira2006fast, one monitors the mean step magnitude $e_t = n^(-1) sum_(i=1)^n norm(x_i^((t+1)) - x_i^((t)))_2$ or the normalized Shannon entropy $H_t$ of the step-length histogram, halting when $e_t <= epsilon$ (or when all shifts equalize and inter-step variance collapses), and then extracting clusters as connected components of the $delta$-proximity graph on $X^((t))$.

#algbox(
  title: "Gaussian Blurring Mean Shift (GBMS)",
  inputs: [Dataset $bold(X)^((0)) in bb(R)^(n times d)$, bandwidth $h > 0$, tolerances $epsilon, delta > 0$, max steps $T_"max"$],
  outputs: [Cluster labels $c_1, dots, c_n$ and converged mode representatives ${m_1, dots, m_K}$],
)[
  1. Set $t <- 0$ and initialize $bold(X)^((0)) <- [x_1^((0)), dots, x_n^((0))]^top$.
  2. *Repeat* for $t = 0, 1, dots, T_"max" - 1$:
     (a) Compute pairwise squared distances $Delta_(i j)^((t)) <- norm(x_i^((t)) - x_j^((t)))_2^2$ and Gaussian affinities $W_(i j)^((t)) <- exp(-Delta_(i j)^((t)) / (2 h^2))$.
     (b) Compute row sums $d_i^((t)) <- sum_(j=1)^n W_(i j)^((t))$ and row-stochastic transition matrix $bold(P)^((t)) <- (bold(D)^((t)))^(-1) bold(W)^((t))$.
     (c) Update the entire dataset synchronously: $bold(X)^((t+1)) <- bold(P)^((t)) bold(X)^((t))$.
     (d) Compute step metric $e_t <- n^(-1) sum_(i=1)^n norm(x_i^((t+1)) - x_i^((t)))_2$ and step-length variance $v_t$.
     (e) *If* $e_t <= epsilon$ or ($t >= 2$ and $v_t <= epsilon^2$ and $e_t$ stabilizes), *break*.
  3. *Cluster Extraction*: Form connected components of the graph connecting $i, j$ whenever $norm(x_i^((t+1)) - x_j^((t+1)))_2 <= delta$; assign cluster label $c_i$ by component index and set mode $m_a$ to the component mean.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five formal results establishing convex-hull and spectral energy dissipation, an exact finite-sample cubic contraction bound for arbitrary point clouds via the Dobrushin ergodicity coefficient, Carreira-Perpiñán's continuous Gaussian variance recursion, exact multi-cluster separation for compact kernels, and a constructive proof of mode-displacement bias relative to standard Mean Shift.

#subsec("A", "Convex Hull Nesting and Spectral Energy Dissipation")

#lemma(title: "Convex Hull Monotonicity and Weighted Scatter Dissipation", origin: "Original Proof")[
  For any nonnegative kernel profile $g$ with $g(0) > 0$, the Blurring Mean Shift iterates satisfy:
  1. *(Convex hull nesting)* $"conv"(X^((t+1))) subset.eq "conv"(X^((t)))$ for every $t >= 0$.
  2. *(Positive-definite spectral dissipation)* If the kernel $K(x, y) = g(norm(x - y)_2^2 / h^2)$ is positive semidefinite (in particular, Gaussian), then $bold(P)^((t))$ is self-adjoint in the $bold(D)^((t))$-weighted inner product $chevron.l u, v chevron.r_(bold(D)^((t))) = u^top bold(D)^((t)) v$ with all eigenvalues in $[0, 1]$, and for any fixed center $c in bb(R)^d$,
     $
       sum_(i=1)^n d_i^((t)) norm(x_i^((t+1)) - c)_2^2
       &+ sum_(i=1)^n d_i^((t)) norm(x_i^((t+1)) - x_i^((t)))_2^2 \
       &<= sum_(i=1)^n d_i^((t)) norm(x_i^((t)) - c)_2^2.
     $
]

#proof[
  *Part 1.* Because $W_(i j)^((t)) >= 0$ and $W_(i i)^((t)) = g(0) > 0$, each row of $bold(P)^((t))$ consists of nonnegative weights summing to $1$. Thus every $x_i^((t+1)) = sum_(j=1)^n P_(i j)^((t)) x_j^((t))$ is a convex combination of $X^((t))$, implying $X^((t+1)) subset.eq "conv"(X^((t)))$ and hence $"conv"(X^((t+1))) subset.eq "conv"(X^((t)))$.

  *Part 2.* Suppressing the superscript $(t)$, observe that $bold(D) bold(P) = bold(W) = bold(W)^top = bold(P)^top bold(D)$, so $bold(P)$ is self-adjoint with respect to $chevron.l dot, dot chevron.r_bold(D)$. Moreover, $bold(P)$ is similar to the symmetric matrix $bold(A) = bold(D)^(-1/2) bold(W) bold(D)^(-1/2)$. By positive semidefiniteness of the kernel (for Gaussian kernels, by Bochner's theorem), $bold(W) succ.eq 0$, whence $bold(A) succ.eq 0$. Since $bold(P)$ is row-stochastic, its spectral radius is $rho(bold(A)) = rho(bold(P)) <= norm(bold(P))_infinity = 1$; thus every eigenvalue $lambda_r$ of $bold(A)$ satisfies $0 <= lambda_r <= 1$.

  Fix a coordinate $s in {1, dots, d}$ and let $u in bb(R)^n$ have entries $u_i = x_(i, s)^((t)) - c_s$. Because $bold(P) bold(1) = bold(1)$, we have $(bold(P) u)_i = x_(i, s)^((t+1)) - c_s$. Setting $v = bold(D)^(1/2) u$, we obtain
  $
    norm(bold(P) u)_bold(D)^2 + norm(bold(P) u - u)_bold(D)^2
    &= norm(bold(A) v)_2^2 + norm((bold(A) - bold(I)) v)_2^2 \
    &= v^top (bold(I) - 2 bold(A)(bold(I) - bold(A))) v.
  $
  Since $lambda_r (1 - lambda_r) >= 0$ for every eigenvalue $lambda_r in [0, 1]$ of $bold(A)$, the symmetric matrix $bold(A)(bold(I) - bold(A))$ is positive semidefinite, giving $norm(bold(P) u)_bold(D)^2 + norm(bold(P) u - u)_bold(D)^2 <= norm(v)_2^2 = norm(u)_bold(D)^2$. Summing over coordinates $s = 1, dots, d$ proves the inequality.
]

#subsec("B", "Universal Finite-Sample Cubic Convergence via Dobrushin Contraction")

Standard analyses of GBMS prove cubic contraction only for continuous Gaussian distributions. We now establish a *distribution-free*, finite-sample cubic contraction theorem holding for *every* finite dataset in $bb(R)^d$.

#theorem(title: "Universal Dobrushin Diameter Contraction and Cubic Convergence", origin: "Original Proof")[
  Let $D_t = "diam"(X^((t))) = max_(1 <= i, l <= n) norm(x_i^((t)) - x_l^((t)))_2$ denote the Euclidean diameter of the $t$-th Blurring Mean Shift iterate under any strictly positive nonincreasing profile $g$ with $g(0) = 1$. Then for every $t >= 0$:
  $
    D_(t+1) <= (1 - g(D_t^2 / h^2)) D_t <= (1 - g(D_0^2 / h^2))^(t+1) D_0.
  $
  Consequently, $X^((t))$ collapses globally to a single point $x^* in "conv"(X^((0)))$ as $t -> infinity$. Moreover, for Gaussian Blurring Mean Shift ($g(r) = exp(-r / 2)$), for every $t >= 0$ and any finite dataset $X^((0)) subset bb(R)^d$,
  $
    D_(t+1) <= (1 - exp(-D_t^2 / (2 h^2))) D_t <= (D_t^3) / (2 h^2),
  $
  so whenever $D_(t_0) < sqrt(2) h$, the diameter shrinks *cubically* as $D_(t_0 + s) <= sqrt(2) h (D_(t_0) / (sqrt(2) h))^(3^s)$ for all $s >= 0$.
]

#proof[
  Fix $t >= 0$ and any pair of indices $i, l in {1, dots, n}$. Let $eta_t = min_(r, j) P_(r j)^((t))$. Because $0 <= norm(x_r^((t)) - x_j^((t)))_2 <= D_t$ and $g$ is nonincreasing with $g(0) = 1$, we have $g(D_t^2 / h^2) <= W_(r j)^((t)) <= 1$, whence $d_r^((t)) = sum_(j=1)^n W_(r j)^((t)) <= n$ and
  $
    eta_t = min_(r, j) (W_(r j)^((t))) / (d_r^((t))) >= (g(D_t^2 / h^2)) / n.
  $
  Decompose the signed row difference $P_(i j)^((t)) - P_(l j)^((t)) = a_j - b_j$ into nonnegative parts $a_j = max(P_(i j)^((t)) - P_(l j)^((t)), 0)$ and $b_j = max(P_(l j)^((t)) - P_(i j)^((t)), 0)$. Since $sum_(j=1)^n (P_(i j)^((t)) - P_(l j)^((t))) = 1 - 1 = 0$, we have
  $
    sum_(j=1)^n a_j = sum_(j=1)^n b_j
    &= 1 / 2 sum_(j=1)^n abs(P_(i j)^((t)) - P_(l j)^((t))) \
    &=: theta_(i l).
  $
  Using $min(p, q) = 1 / 2 (p + q - abs(p - q))$ and summing over $j = 1, dots, n$ yields the Dobrushin ergodicity bound:
  $
    1 - theta_(i l) = sum_(j=1)^n min(P_(i j)^((t)), P_(l j)^((t))) >= n eta_t >= g(D_t^2 / h^2),
  $
  so $theta_(i l) <= 1 - g(D_t^2 / h^2)$. If $theta_(i l) = 0$, then $x_i^((t+1)) = x_l^((t+1))$. Otherwise, define the two convex combinations $mu = theta_(i l)^(-1) sum_(j=1)^n a_j x_j^((t))$ and $nu = theta_(i l)^(-1) sum_(j=1)^n b_j x_j^((t))$, both lying in $"conv"(X^((t)))$. Because the Euclidean diameter of a compact set equals that of its convex hull,
  $
    norm(x_i^((t+1)) - x_l^((t+1)))_2
    = theta_(i l) norm(mu - nu)_2
    <= (1 - g(D_t^2 / h^2)) D_t.
  $
  Maximizing over $(i, l)$ gives $D_(t+1) <= (1 - g(D_t^2 / h^2)) D_t$. Since $D_t <= D_0$ inductively and $g$ is nonincreasing, $1 - g(D_t^2 / h^2) <= 1 - g(D_0^2 / h^2) < 1$, yielding $D_t <= (1 - g(D_0^2 / h^2))^t D_0 -> 0$. By Cantor's intersection theorem applied to the nested compact sets $"conv"(X^((t)))$ (Lemma 1), $inter.big_(t=0)^infinity "conv"(X^((t))) = {x^*}$ is a single point.

  Finally, for GBMS ($g(r) = exp(-r / 2)$), the elementary inequality $1 - exp(-u) <= u$ ($u >= 0$) with $u = D_t^2 / (2 h^2)$ immediately yields
  $
    D_(t+1) <= (1 - exp(-D_t^2 / (2 h^2))) D_t <= (D_t^3) / (2 h^2).
  $
  Rewriting as $D_(t+1) / (sqrt(2) h) <= (D_t / (sqrt(2) h))^3$ and iterating $s$ times yields $D_(t_0 + s) <= sqrt(2) h (D_(t_0) / (sqrt(2) h))^(3^s)$.
]

#subsec("C", "Exact Gaussian Covariance Recursion and Compact-Kernel Decoupling")

#theorem(title: "Exact Cubic Covariance Collapse for Gaussian Clusters", origin: "Classical — Cited")[
  Consider the population limit of GBMS with bandwidth $h > 0$ initialized at a $d$-dimensional Gaussian distribution $X^((0)) tilde cal(N)(mu, Sigma_0)$ with $Sigma_0 succ 0$ @carreira2006fast. Then for every $t >= 0$, the blurred distribution remains exactly Gaussian $X^((t)) tilde cal(N)(mu, Sigma_t)$ with covariance recursion
  $
    Sigma_(t+1) = Sigma_t^3 (Sigma_t + h^2 bold(I)_d)^(-2).
  $
  Along each principal eigenvector of $Sigma_0$, the standard deviation $sigma_r^((t)) = sqrt(lambda_r (Sigma_t))$ satisfies $sigma_r^((t+1)) = (sigma_r^((t)))^3 / ((sigma_r^((t)))^2 + h^2)$, achieving cubic convergence, whereas non-blurring Mean Shift converges linearly at rate $(sigma_r^((0)))^2 / ((sigma_r^((0)))^2 + h^2)$.
]

#proof[
  Without loss of generality translate so $mu = bold(0)$ and assume inductively $X^((t)) tilde cal(N)(bold(0), Sigma_t)$. The population GBMS update map at $x in bb(R)^d$ is the conditional expectation
  $
    f_t (x) = (integral_(bb(R)^d) y exp(-norm(x - y)_2^2 / (2 h^2)) p_t (y) dif y) / (integral_(bb(R)^d) exp(-norm(x - y)_2^2 / (2 h^2)) p_t (y) dif y),
  $
  where $p_t (y) prop exp(-1 / 2 y^top Sigma_t^(-1) y)$. Combining the quadratic forms in the exponent gives precision matrix $Lambda_t = Sigma_t^(-1) + h^(-2) bold(I)_d$ and linear term $h^(-2) x^top y$, so the normalized integrand in $y$ is Gaussian with mean
  $
    f_t (x) = (Sigma_t^(-1) + h^(-2) bold(I)_d)^(-1) (h^(-2) x)
            = Sigma_t (Sigma_t + h^2 bold(I)_d)^(-1) x.
  $
  Because $f_t$ is linear, $X^((t+1)) = f_t (X^((t)))$ is zero-mean Gaussian with covariance
  $
    Sigma_(t+1) = bold(M)_t Sigma_t bold(M)_t^top
                = Sigma_t^3 (Sigma_t + h^2 bold(I)_d)^(-2),
  $
  where $bold(M)_t = Sigma_t (Sigma_t + h^2 bold(I)_d)^(-1)$ commutes with $Sigma_t$. Diagonalizing $Sigma_0 = bold(U) "diag"(lambda_1^((0)), dots, lambda_d^((0))) bold(U)^top$ preserves the eigenbasis $bold(U)$ for all $t$, giving $lambda_r^((t+1)) = (lambda_r^((t)))^3 / (lambda_r^((t)) + h^2)^2$ and thus $sigma_r^((t+1)) = (sigma_r^((t)))^3 / ((sigma_r^((t)))^2 + h^2) <= (sigma_r^((t)))^3 / h^2$. By contrast, non-blurring Mean Shift applies the fixed matrix $bold(M)_0$ repeatedly, yielding linear decay $y^((t)) = bold(M)_0^t y^((0))$ @carreira2006fast.
]

#proposition(title: "Block-Diagonal Decoupling Under Compactly Supported Kernels", origin: "Original Proof")[
  Suppose the shadow profile $g$ has compact support on $[0, 1]$ ($g(r) = 0$ for $r > 1$, e.g., Epanechnikov or flat kernel). If at some iteration $t_0 >= 0$, $X^((t_0))$ partitions into disjoint subsets $C_1^((t_0)), dots, C_K^((t_0))$ whose convex hulls are separated by strictly more than $h$,
  $
    min_(a != b) "dist"("conv"(C_a^((t_0))), "conv"(C_b^((t_0)))) > h,
  $
  then for every $t >= t_0$, $bold(P)^((t))$ is strictly block-diagonal across $C_1, dots, C_K$, the $K$ components never merge, and any component satisfying $"diam"(C_a^((t_0))) <= h$ with $g$ strictly positive on $[0, 1]$ collapses to its own distinct limit mode $m_a^* in "conv"(C_a^((t_0)))$ with $norm(m_a^* - m_b^*)_2 > h$ for $a != b$.
]

#proof[
  We induct on $t >= t_0$ for the invariant $"conv"(C_a^((t))) subset.eq "conv"(C_a^((t_0)))$ for all $a in {1, dots, K}$. The base case $t = t_0$ is trivial. Assuming the invariant holds at step $t$, any $i in C_a$ and $j in C_b$ with $a != b$ satisfy
  $
    norm(x_i^((t)) - x_j^((t)))_2
    >= "dist"("conv"(C_a^((t_0))), "conv"(C_b^((t_0)))) > h.
  $
  Hence $W_(i j)^((t)) = g(norm(x_i^((t)) - x_j^((t)))_2^2 / h^2) = 0$, so $bold(P)^((t))$ is block-diagonal with zero off-block entries. Consequently, every $x_i^((t+1))$ ($i in C_a$) is a convex combination of $C_a^((t))$ alone, implying $"conv"(C_a^((t+1))) subset.eq "conv"(C_a^((t))) subset.eq "conv"(C_a^((t_0)))$, which closes the induction. Finally, if $"diam"(C_a^((t_0))) <= h$ and $g > 0$ on $[0, 1]$, then restricting Theorem 1 to block $C_a$ yields $"diam"(C_a^((t))) -> 0$, so $C_a^((t))$ converges to a unique $m_a^* in "conv"(C_a^((t_0)))$, and $norm(m_a^* - m_b^*)_2 > h$ follows from the convex-hull separation.
]

#subsec("D", "Constructive Mode-Bias Discrepancy vs. Standard Mean Shift")

#proposition(title: "Systematic Mode Displacement of Blurring Mean Shift", origin: "Original Counterexample")[
  Consider the three-point dataset $X^((0)) = {0, 0, 1} subset bb(R)$ with Gaussian kernel bandwidth $h = 1$. The unique global mode $m_"MS"^*$ of the original KDE $hat(p)_1 (x)$ and the unique GBMS limit point $m_"BMS"^*$ satisfy the strict separation
  $
    m_"MS"^* < 3 / 10 < 668 / 2223 < m_"BMS"^* < 1 / 2.
  $
  Hence Blurring Mean Shift does not in general converge to the stationary modes of the initial KDE, exhibiting a systematic drift toward the local sample barycenter.
]

#proof[
  *Upper bound on $m_"MS"^*$.* Up to a positive constant, $hat(p)_1 (x) = 2 exp(-x^2 / 2) + exp(-(x - 1)^2 / 2)$, whose derivative is $hat(p)'_1 (x) = -2 x exp(-x^2 / 2) + (1 - x) exp(-(1 - x)^2 / 2)$. For $x in (0, 1)$, $hat(p)'_1 (x) < 0$ is equivalent to $2 x / (1 - x) > exp(x - 1 / 2)$, or $psi(x) := ln(2 x) - ln(1 - x) - x + 1 / 2 > 0$. Differentiating yields
  $
    psi'(x) = 1 / x + 1 / (1 - x) - 1 = 1 / (x (1 - x)) - 1 >= 4 - 1 = 3 > 0,
  $
  so $psi$ is strictly increasing on $(0, 1)$. At $x = 3 / 10$, we have $2 x / (1 - x) = 6 / 7$, whereas $exp(-1 / 5) < 1 / (1 + 1 / 5) = 5 / 6 < 6 / 7$; hence $psi(3 / 10) > 0$ and $hat(p)'_1 (x) < 0$ for all $x in [3 / 10, 1]$. Since $hat(p)'_1 (0) > 0$, the unique global mode obeys $m_"MS"^* < 3 / 10$.

  *Lower bound on $m_"BMS"^*$.* Let $w = exp(-1 / 2)$. Since $2.6 < e < 2.75$, squaring verifies $9 / 25 < 1 / e < 25 / 64$, so $3 / 5 < w < 5 / 8$. After one GBMS step from $x_1^((0)) = x_2^((0)) = 0$, $x_3^((0)) = 1$:
  $
    x_1^((1)) = x_2^((1)) &= w / (2 + w) > (3 / 5) / (13 / 5) = 3 / 13, \
    x_3^((1)) &= 1 / (1 + 2 w) > 1 / (1 + 5 / 4) = 4 / 9.
  $
  Moreover, $D_1 = x_3^((1)) - x_1^((1)) < 1 / (1 + 6 / 5) - 3 / 13 = 5 / 11 - 3 / 13 = 32 / 143 < 1 / 4$, so $w_1 = exp(-D_1^2 / 2) > 1 - D_1^2 / 2 > 31 / 32$. At step $t = 2$, $x_1^((2)) = (2 x_1^((1)) + w_1 x_3^((1))) / (2 + w_1)$ is strictly increasing in $x_1^((1))$, $x_3^((1))$, and $w_1$ (since $x_3^((1)) > x_1^((1))$), yielding
  $
    x_1^((2)) &> (2 (3 / 13) + (31 / 32)(4 / 9)) / (2 + 31 / 32) \
              &= (6 / 13 + 31 / 72) / (95 / 32) = 668 / 2223 > 3 / 10.
  $
  Since $x_3^((1)) = 1 / (1 + 2 w) < 1 / (1 + 6 / 5) = 5 / 11 < 1 / 2$, convex-hull nesting (Lemma 1) guarantees $m_"BMS"^* in [x_1^((2)), x_3^((1))] subset (668 / 2223, 1 / 2)$, completing the proof.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Blurring Mean Shift evaluates $n (n - 1) / 2$ pairwise distances and performs a dense matrix product $bold(P)^((t)) bold(X)^((t))$, costing $cal(O)(n^2 d)$ arithmetic operations per step and $cal(O)(T n^2 d)$ overall. Crucially, however, cubic local contraction (Theorems 1--2) reduces the iteration count $T$ to typically $3$--$8$ steps—an order-of-magnitude reduction over the $T = cal(O)(log(1 / epsilon))$ steps required by non-blurring Mean Shift @carreira2006fast @comaniciu2002mean. Furthermore, because the entire update is expressed as BLAS-3 dense matrix multiplication (distance Gram matrix plus $bold(P)^((t)) bold(X)^((t))$), GBMS exhibits exceptional vectorization, cache locality, and GPU parallelism for moderate $n$ ($n <= 10^4$). Memory consumption is $cal(O)(n^2 + n d)$ for explicit $bold(W)^((t))$ storage, or $cal(O)(n d + b n)$ via tiled block-row evaluation of block size $b$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(T n^2 d)$ with very small $T$ (typically $3$--$8$ steps via cubic contraction)],
    [Space complexity], [$cal(O)(n^2 + n d)$ dense; $cal(O)(n d + b n)$ with row-blocked matrix streaming],
    [Hyperparameters], [Bandwidth $h > 0$, stopping tolerance $epsilon$, merge radius $delta$],
    [Convergence rate], [Cubic ($D_(t+1) <= D_t^3 / (2 h^2)$) locally; globally collapses without stopping],
    [Symmetries], [Exact Euclidean isometry, permutation, and joint $(x, h)$ scale invariance],
    [Cluster geometry], [Non-convex manifolds & arbitrary cluster counts $K$ discovered automatically],
    [Main trade-off], [Dramatic speedup vs. stopping sensitivity and slight barycentric mode bias],
  ),
  caption: [Operational and computational profile of Blurring Mean Shift.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Fast Color-Space Image Segmentation and Quantization*: In 3D $(L^*, a^*, b^*)$ or 5D spatial-color spaces, pixel clouds form dense, non-convex color manifolds. GBMS contracts homogeneous color regions in under six iterations @carreira2006fast, achieving segmentation quality comparable to standard Mean Shift @comaniciu2002mean at $10 times$--$50 times$ lower wall-clock latency when combined with spatial-color subsampling or superpixels.
2. *Point-Cloud Denoising and Manifold Sharpening*: Running 1--3 iterations of BMS acts as an anisotropic, structure-preserving bilateral filter that contracts Gaussian sensor noise orthogonally onto underlying 1D curves or 2D surfaces while preserving sharp gaps wider than $h$ (Proposition 1).
3. *Flow Cytometry and Single-Cell Phenotyping*: Biological cell populations in 4D--10D fluorescence marker spaces are skewed and non-spherical, with unknown cluster count $K$. GBMS discovers irregular sub-populations automatically without requiring $K$ upfront, and cubic contraction cleanly separates closely packed immune-cell lineages.
4. *Tabular Customer and Market Segmentation*: Mid-sized customer-behavior datasets ($n tilde 10^3$--$10^4$, $d <= 10$) often exhibit overlapping spending modes and skewed tails. Because BMS preserves Euclidean invariance and rapidly consolidates dense behavioral cohorts before global collapse, it yields higher compactness than single-scale density cuts.
5. *Multi-Target Radar and LiDAR Detection Consolidation*: Raw sensor detections produce noisy point clouds around each physical target. A handful of GBMS updates collapses each target cloud onto a crisp centroid estimate while keeping distinct vehicles separated whenever inter-target spacing exceeds the bandwidth $h$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Blurring Mean Shift encounters four structural failure modes: (i) *over-blurring collapse* (Theorem 1), where an overly loose stopping tolerance $epsilon$ lets metastable clusters drift together and merge into a single trivial mode; (ii) *high-dimensional Gaussian flatness* @beyer1999nearest, where pairwise squared distances concentrate around $2 d sigma^2$, making $W_(i j) = exp(-Delta_(i j) / (2 h^2))$ either near-uniform (instant collapse) or near-identity (zero movement); (iii) *quadratic $cal(O)(n^2)$ scaling* for $n > 10^5$, where dense affinity matrices become memory- and compute-bound; and (iv) *asymmetric mode bias* (Proposition 2), where sharp density peaks adjacent to heavy shoulders drift toward the local barycenter rather than the exact KDE mode.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Image & color segmentation ($d <= 5$)], [High: cubic collapse on non-convex color modes], [For $n > 10^5$ pixels, use Mean Shift++ (Alg. 55)],
    [Point-cloud denoising & sharpening], [High: 1--3 iterations denoise manifolds smoothly], [Halt early before global collapse (Theorem 1)],
    [Cytometry & tabular cohorts ($d <= 10$)], [High: automatic $K$, strong non-convex fidelity], [Standardize feature scales; tune $h$ via $k$-NN],
    [High-D embeddings ($d > 20$)], [Poor: Gaussian kernel degenerates to $bold(I)$ or $bold(1) bold(1)^top$], [Apply PCA / UMAP first or cosine clustering],
    [Exact KDE mode estimation], [Moderate: barycentric shift bias (Proposition 2)], [Use non-blurring Mean Shift or DENCLUE 2.0],
  ),
  caption: [Domain suitability matrix for Blurring Mean Shift.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Blurring Mean Shift bridges nonparametric mode-seeking (@cheng1995mean, @comaniciu2002mean, @ariascastro2016meanshift) and spectral graph diffusion (@vonluxburg2007tutorial). The following algebraic identity makes this structural link exact.

#proposition(title: "Equivalence to Time-Varying Graph-Laplacian Diffusion", origin: "Original Derivation")[
  Let $bold(L)_"rw"^((t)) = bold(I)_n - (bold(D)^((t)))^(-1) bold(W)^((t))$ denote the random-walk normalized graph Laplacian of the $t$-th affinity graph $(X^((t)), bold(W)^((t)))$. Then one Blurring Mean Shift iteration is an exact unit-step explicit Euler integration of the nonlinear heat equation on the evolving point cloud:
  $
    bold(X)^((t+1)) - bold(X)^((t)) = -bold(L)_"rw"^((t)) bold(X)^((t)).
  $
  In particular, the shift vector $m^((t)) (x_i^((t))) = x_i^((t+1)) - x_i^((t))$ is the negative random-walk Laplacian applied to the coordinate functions, directly linking BMS to spectral clustering's diffusion operator @vonluxburg2007tutorial while updating the graph weights adaptively at every step.
]

#proof[
  By Definition 1, $bold(X)^((t+1)) = bold(P)^((t)) bold(X)^((t)) = (bold(D)^((t)))^(-1) bold(W)^((t)) bold(X)^((t))$. Subtracting $bold(X)^((t))$ from both sides gives $bold(X)^((t+1)) - bold(X)^((t)) = -(bold(I)_n - (bold(D)^((t)))^(-1) bold(W)^((t))) bold(X)^((t)) = -bold(L)_"rw"^((t)) bold(X)^((t))$, which is the explicit Euler step with step size $Delta tau = 1$ for the coordinate diffusion $dif bold(X) / dif tau = -bold(L)_"rw" (bold(X)) bold(X)$.
]

Within Category II, Blurring Mean Shift sits directly between standard Mean Shift @comaniciu2002mean (which freezes $X^((0))$, converging linearly to exact KDE modes), DENCLUE 2.0 @hinneburg1998denclue @hinneburg2007denclue2 (which accelerates standard Mean Shift via hill-climbing steps), Quick Shift @vedaldi2008quickshift and Medoid Shift @sheikh2007medoidshift (which restrict updates to observed data points), and Mean Shift++ @jang2021meanshiftpp (Alg. 55, which accelerates blurring-style updates to $cal(O)(n)$ via a Cartesian grid hash table). Unlike spectral clustering @vonluxburg2007tutorial, which extracts the leading eigenvectors of a static $bold(L)_"rw"^((0))$ for a pre-specified $K$, BMS repeatedly applies $bold(I) - bold(L)_"rw"^((t))$ directly to coordinates and discovers $K$ automatically @hastie2009elements @jain2010data.

Empirically, on the ClusterOverview 4-dataset benchmark evaluated by Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, Blurring Mean Shift exhibits standout geometric fidelity with a mean runtime of `4.568 ms`. On Iris (`Sil = 0.586`, `ARI = 0.567`), BMS improves ARI substantially over OPTICS-OF (`0.481`) while maintaining strong Silhouette, confirming that cubic Gaussian contraction sharpens the _Setosa_ and _Versicolor_/_Virginica_ modes effectively. On Customers (`Sil = 0.382`), BMS achieves the strongest tabular compactness among Algorithms 53--55 (beating OPTICS-OF's `0.281` and producing a valid multi-cluster partition where grid-based Mean Shift++ yields `N/A`), illustrating the value of smooth rotationally invariant Gaussian weights on heterogeneous feature spreads. On Multi-Shape (`Sil = 0.787`, `ARI = N/A`), BMS attains an exceptionally high Silhouette of `0.787`—far surpassing OPTICS-OF (`0.446`) and Mean Shift++ (`0.567`)—because its metastable plateau contracts each curved shape cleanly before inter-shape drift occurs. Finally, on High-D (`Sil = N/A`, `ARI = N/A`), BMS fails to form a non-degenerate partition, matching our analysis of high-dimensional Gaussian kernel flatness and distance concentration @beyer1999nearest, where $exp(-Delta_(i j) / (2 h^2))$ either freezes all points or collapses them uniformly into one cluster.
