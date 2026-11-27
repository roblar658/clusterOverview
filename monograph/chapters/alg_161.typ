#import "../template/ieee.typ": *

#alg_entry(
  id: 161,
  name: "FABIA (Factor Analysis for Bicluster Acquisition)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2010",
  paradigm: "Sparse Multiplicative Variational Factor Analysis",
  complexity: $O(T (p^3 (n + m) + p n m))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

*FABIA* (Factor Analysis for Bicluster Acquisition) @hochreiter2010fabia formulates biclustering as a sparse, multiplicative latent-variable generative model. Unlike hard combinatorial extraction in Bimax @prelic2006bimax or heuristic iterative thresholding in ISA @bergmann2003isa, FABIA models an observed data matrix $X in bb(R)^(n times m)$ ($n$ features/genes by $m$ samples/conditions) as a superposition of $p$ sparse rank-$1$ outer products corrupted by heteroscedastic row-specific Gaussian noise:
$ X = sum_(k=1)^p bold(lambda)_k bold(z)_k^top + bold(Upsilon) = bold(Lambda) bold(Z) + bold(Upsilon), $
where $bold(Lambda) = [bold(lambda)_1, dots, bold(lambda)_p] in bb(R)^(n times p)$ collects the sparse prototype (loading) vectors, $bold(Z) in bb(R)^(p times m)$ collects the sparse factor weights with $j$-th sample column $bold(z)_j in bb(R)^p$ (so row $k$ of $bold(Z)$ is $bold(z)_k^top$ selecting samples belonging to bicluster $k$), and each noise column $bold(epsilon)_j in bb(R)^n$ of $bold(Upsilon)$ is drawn independently from $cal(N)(bold(0), bold(Psi))$ with diagonal covariance $bold(Psi) = "diag"(psi_1, dots, psi_n)$, $psi_i > 0$.

A bicluster is encoded jointly by the non-zero support of $bold(lambda)_k$ (rows) and the $k$-th row of $bold(Z)$ (columns). To enforce simultaneous sparsity on both rows and columns — and to break the orthogonal rotational ambiguity of classical Gaussian factor analysis — FABIA places independent zero-mean *Laplace priors* on both the latent factors $bold(z)_j$ (standardized to unit variance, rate $sqrt(2)$) and the loadings $bold(Lambda)$ (rate $alpha > 0$) @hochreiter2010fabia:
$ p(bold(z)_j) &= product_(k=1)^p (1 / sqrt(2)) exp(-sqrt(2) |z_(k j)|), \
  p(bold(Lambda)) &= product_(i=1)^n product_(k=1)^p (alpha / 2) exp(-alpha |Lambda_(i k)|). $

#definition(title: "FABIA Variational MAP Biclustering Objective")[
  Given $X = [bold(x)_1, dots, bold(x)_m] in bb(R)^(n times m)$, factor count $p >= 1$, sparsity rate $alpha > 0$, and variance floor $psi_"min" > 0$, FABIA maximizes a variational lower bound $cal(L)(bold(Lambda), bold(Psi), {q_j, bold(Xi)_j}_(j=1)^m)$ on the marginalized log-posterior
  $ log p(bold(Lambda) | X, bold(Psi)) + "const" = log p(bold(Lambda)) + sum_(j=1)^m log integral_(bb(R)^p) p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) p(bold(z)_j) dif bold(z)_j, $
  where each variational posterior $q_j (bold(z)_j) = cal(N)(bold(mu)_(z_j), bold(Sigma)_(z_j))$ is Gaussian and each sample carries a positive diagonal variational scale matrix $bold(Xi)_j = "diag"(xi_(1 j), dots, xi_(p j)) succ 0$ arising from the quadratic lower bound on the Laplace prior.
]

After convergence, biclusters are ranked by their variational mutual information content $I_k$, and discrete row/column memberships $(I_k, J_k)$ are extracted by thresholding $|Lambda_(i k)|$ and $|(bold(mu)_(z_j))_k|$ against scale-normalized cutoffs @hochreiter2010fabia.

#algbox(
  title: "FABIA Variational EM Algorithm",
  inputs: [Centered matrix $X in bb(R)^(n times m)$; number of biclusters $p$; sparsity $alpha > 0$; variance floor $psi_"min" > 0$; max iterations $T$.],
  outputs: [Sparse loadings $bold(Lambda) in bb(R)^(n times p)$, factor means $M_Z = [bold(mu)_(z_1), dots, bold(mu)_(z_m)]$, noise covariance $bold(Psi)$, and information scores ${I_k}_(k=1)^p$.],
)[
  1. *Initialization:* Initialize $bold(Lambda)^((0))$ randomly (small variance), $bold(Psi)^((0)) <- "diag"("Var"(X_(1,:)), dots, "Var"(X_(n,:)))$, and $bold(Xi)_j^((0)) <- I_p$ for $j = 1, dots, m$.
  2. *Repeat* for $t = 1, dots, T$ until $cal(L)^((t)) - cal(L)^((t-1)) < epsilon$:
  3. *Variational E-Step:* Compute $G <- bold(Lambda)^top bold(Psi)^(-1) bold(Lambda)$. For each sample $j = 1, dots, m$:
     $ bold(Sigma)_(z_j) &<- (bold(Xi)_j^(-1) + G)^(-1), \
       bold(mu)_(z_j) &<- bold(Sigma)_(z_j) bold(Lambda)^top bold(Psi)^(-1) bold(x)_j, \
       bold(S)_j &<- bold(mu)_(z_j) bold(mu)_(z_j)^top + bold(Sigma)_(z_j). $
  4. *Variational Scale Update:* For all $k = 1, dots, p$ and $j = 1, dots, m$:
     $ xi_(k j) <- sqrt((bold(S)_j)_(k k) / 2) = sqrt(((bold(mu)_(z_j))_k^2 + (bold(Sigma)_(z_j))_(k k)) / 2). $
  5. *Sparse M-Step for $bold(Lambda)$:* Form $A <- sum_(j=1)^m bold(S)_j$ and $B <- X M_Z^top$. For each row $i = 1, dots, n$, update $bold(Lambda)_(i,:)$ by solving the strictly convex $ell_1$-penalized quadratic program
     $ min_(bold(lambda) in bb(R)^p) (1 / (2 psi_i)) (bold(lambda)^top A bold(lambda) - 2 bold(b)_i^top bold(lambda)) + alpha norm(bold(lambda))_1. $
  6. *Heteroscedastic M-Step for $bold(Psi)$:* For each row $i = 1, dots, n$:
     $ psi_i <- max(psi_"min", (1 / m) (norm(X_(i,:))_2^2 - 2 bold(Lambda)_(i,:) bold(b)_i + bold(Lambda)_(i,:) A bold(Lambda)_(i,:)^top)). $
  7. *Bicluster Ranking:* Compute $I_k = (1 / 2) sum_(j=1)^m log(1 + xi_(k j) bold(lambda)_k^top bold(Psi)^(-1) bold(lambda)_k)$ for $k = 1, dots, p$ and threshold $(bold(Lambda), M_Z)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Quadratic Lower Bound for the Multivariate Laplace Prior")

Exact Bayesian inference under a Laplace factor prior is intractable because the posterior over $bold(z)_j$ is a piecewise-Gaussian density over $2^p$ orthants. We prove the exactness and global minorization property of FABIA's convex-duality Gaussian envelope @hochreiter2010fabia @blei2006variational.

#lemma(title: "Tight Laplace-Gaussian Variational Envelope", origin: "Original Derivation")[
  For every $bold(z) in bb(R)^p$ and every diagonal matrix $bold(Xi) = "diag"(xi_1, dots, xi_p) succ 0$, define the Gaussian-form surrogate
  $ tilde(p)(bold(z); bold(Xi)) = 2^(-p / 2) exp(- (1 / 2) bold(z)^top bold(Xi)^(-1) bold(z) - "Tr"(bold(Xi))). $
  Then $p(bold(z)) >= tilde(p)(bold(z); bold(Xi))$ pointwise for all $bold(z) in bb(R)^p$, with equality if and only if $xi_k = |z_k| / sqrt(2)$ for all $k$ (for $z_k != 0$, and in the limit $xi_k -> 0^+$ when $z_k = 0$). Consequently, the variational objective
  $ cal(L)(bold(Lambda), bold(Psi), {q_j, bold(Xi)_j}) &= log p(bold(Lambda)) + sum_(j=1)^m (H(q_j) \
    &+ bb(E)_(q_j) [log p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) + log tilde(p)(bold(z)_j; bold(Xi)_j)]) $
  globally lower-bounds $log p(bold(Lambda)) + sum_(j=1)^m log p(bold(x)_j | bold(Lambda), bold(Psi))$.
]

#proof[
  For any scalar $z in bb(R)$ and any $xi > 0$, completing the square yields
  $ z^2 / (2 xi) + xi - sqrt(2) |z| = (|z| - sqrt(2) xi)^2 / (2 xi) >= 0, $
  with equality if and only if $|z| = sqrt(2) xi$, i.e., $xi = |z| / sqrt(2)$. Negating, exponentiating, and multiplying across $k = 1, dots, p$ establishes
  $ p(bold(z)) = 2^(-p / 2) exp(-sqrt(2) sum_(k=1)^p |z_k|) >= 2^(-p / 2) exp(-sum_(k=1)^p (z_k^2 / (2 xi_k) + xi_k)) = tilde(p)(bold(z); bold(Xi)). $
  For each sample $j$ and any density $q_j$ with support $bb(R)^p$, Jensen's inequality applied to the concave logarithm gives
  $ log p(bold(x)_j | bold(Lambda), bold(Psi)) &= log integral q_j (bold(z)_j) (p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) p(bold(z)_j)) / (q_j (bold(z)_j)) dif bold(z)_j \
    &>= bb(E)_(q_j) [log p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) + log tilde(p)(bold(z)_j; bold(Xi)_j)] + H(q_j), $
  where the inequality combines Jensen's bound with $log p(bold(z)_j) >= log tilde(p)(bold(z)_j; bold(Xi)_j)$. Summing over $j = 1, dots, m$ and adding $log p(bold(Lambda))$ completes the proof.
]

#subsec("B", "Monotone Variational EM Ascent and Boundedness")

We establish rigorous ascent across all four coordinate blocks of FABIA's variational EM loop and prove that the variance floor $psi_"min" > 0$ prevents singular Heywood-case likelihood divergence.

#theorem(title: "Monotonicity and Bounded Convergence of FABIA Variational EM", origin: "Classical — Cited")[
  Under the exact convex $bold(Lambda)$-step and variance floor $psi_i >= psi_"min" > 0$, every block update in Algorithm 161 monotonically non-decreases the variational ELBO $cal(L)$ @hochreiter2010fabia @dempster1977em @wu1983convergence:
  $ cal(L)^((t-1)) <= cal(L)(q^((t))) <= cal(L)(bold(Xi)^((t))) <= cal(L)(bold(Lambda)^((t))) <= cal(L)(bold(Psi)^((t))) = cal(L)^((t)). $
  Moreover, $cal(L)^((t))$ is uniformly bounded above by $n p log(alpha / 2) - (n m / 2) log(2 pi psi_"min") < infinity$, guaranteeing convergence of $cal(L)^((t))$ to a finite stationary value $cal(L)^star$.
]

#proof[
  *1. E-step ($q_j$):* Fixing $(bold(Lambda), bold(Psi), bold(Xi)_j)$, the integrand exponent $log p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) + log tilde(p)(bold(z)_j; bold(Xi)_j)$ is a strictly concave quadratic in $bold(z)_j$ with precision $bold(Sigma)_(z_j)^(-1) = bold(Xi)_j^(-1) + bold(Lambda)^top bold(Psi)^(-1) bold(Lambda) succ 0$ and linear coefficient $bold(Lambda)^top bold(Psi)^(-1) bold(x)_j$. Absorbing $bold(z)_j$-independent constants, maximizing $bb(E)_(q_j) [dots] + H(q_j)$ over all densities $q_j$ equals minimizing $"KL"(q_j || cal(N)(bold(mu)_(z_j), bold(Sigma)_(z_j)))$, uniquely maximized at $q_j^star = cal(N)(bold(mu)_(z_j), bold(Sigma)_(z_j))$ with $bold(mu)_(z_j) = bold(Sigma)_(z_j) bold(Lambda)^top bold(Psi)^(-1) bold(x)_j$.

  *2. Variational scale step ($bold(Xi)_j$):* Isolating terms in $cal(L)$ depending on $xi_(k j) > 0$ yields $phi(xi_(k j)) = - (bold(S)_j)_(k k) / (2 xi_(k j)) - xi_(k j)$, where $(bold(S)_j)_(k k) = (bold(mu)_(z_j))_k^2 + (bold(Sigma)_(z_j))_(k k) > 0$ (since $bold(Sigma)_(z_j) succ 0$). Differentiating gives $phi'(xi_(k j)) = (bold(S)_j)_(k k) / (2 xi_(k j)^2) - 1 = 0 <==> xi_(k j) = sqrt((bold(S)_j)_(k k) / 2)$, with $phi''(xi_(k j)) = - (bold(S)_j)_(k k) / xi_(k j)^3 < 0$, confirming the unique global maximum.

  *3. Loading step ($bold(Lambda)$):* Collecting $bold(Lambda)$-dependent terms in $cal(L)$ gives
  $ sum_(i=1)^n (- (1 / (2 psi_i)) (bold(Lambda)_(i,:) A bold(Lambda)_(i,:)^top - 2 bold(Lambda)_(i,:) bold(b)_i) - alpha norm(bold(Lambda)_(i,:))_1), $
  which decouples across rows $i = 1, dots, n$. Because $A = sum_(j=1)^m (bold(mu)_(z_j) bold(mu)_(z_j)^top + bold(Sigma)_(z_j)) succ 0$ is strictly positive definite, each row subproblem is strictly convex and admits a unique global maximizer, non-decreasing $cal(L)$.

  *4. Noise variance step ($bold(Psi)$):* Let $r_i = norm(X_(i,:))_2^2 - 2 bold(Lambda)_(i,:) bold(b)_i + bold(Lambda)_(i,:) A bold(Lambda)_(i,:)^top = sum_(j=1)^m bb(E)_(q_j) [(x_(i j) - bold(Lambda)_(i,:) bold(z)_j)^2] >= 0$. The $psi_i$-dependent term in $cal(L)$ is $-(m / 2) log psi_i - r_i / (2 psi_i)$, which is strictly increasing on $(0, r_i / m)$ and strictly decreasing on $(r_i / m, infinity)$. Maximizing over $psi_i >= psi_"min"$ uniquely yields $psi_i = max(psi_"min", r_i / m)$.

  *5. Upper Bound:* By Lemma 1, $cal(L)^((t)) <= log p(bold(Lambda)^((t))) + sum_(j=1)^m log p(bold(x)_j | bold(Lambda)^((t)), bold(Psi)^((t)))$. For every $bold(z)_j$, $p(bold(x)_j | bold(z)_j, bold(Lambda), bold(Psi)) <= (2 pi)^(-n / 2) (det bold(Psi))^(-1 / 2) <= (2 pi psi_"min")^(-n / 2)$; integrating against the probability density $p(bold(z)_j)$ preserves this bound, while $log p(bold(Lambda)) <= n p log(alpha / 2)$. Hence $cal(L)^((t))$ is monotone and bounded above, so it converges.
]

#subsec("C", "Breaking Rotational Ambiguity and Information Content Ranking")

In classical Gaussian factor analysis ($bold(z)_j tilde.op cal(N)(bold(0), I_p)$), any orthogonal matrix $Q in O(p)$ maps $(bold(Lambda), bold(Z)) |-> (bold(Lambda) Q, Q^top bold(Z))$ while leaving both the likelihood and prior invariant — destroying bicluster identity. We prove that FABIA's Laplace priors collapse both rotational and scaling ambiguities, and derive its exact mutual-information ranking criterion.

#theorem(title: "Identifiability Up to Signed Permutations and Exact Channel Mutual Information", origin: "Original Proof")[
  1. *(Rotational Identifiability):* An orthogonal matrix $Q in O(p)$ preserves the independent Laplace factor prior $p(Q^top bold(z)) = p(bold(z))$ for all $bold(z) in bb(R)^p$ if and only if $Q$ is a _signed permutation matrix_ ($Q in cal(P)_(plus.minus) (p)$).
  2. *(Scale Fixing):* For any positive diagonal rescaling $D = "diag"(c_1, dots, c_p)$ ($c_k > 0$) preserving the product $bold(Lambda) bold(Z) = (bold(Lambda) D)(D^(-1) bold(Z))$, the joint negative log-prior $alpha norm(bold(Lambda) D)_1 + sqrt(2) norm(D^(-1) bold(Z))_1$ is strictly convex in $(c_1, dots, c_p)$ on $(0, infinity)^p$ with unique minimizer $c_k^star = (sqrt(2) norm(bold(z)_k)_1 / (alpha norm(bold(lambda)_k)_1))^(1 / 2)$ (whenever $bold(lambda)_k, bold(z)_k != bold(0)$).
  3. *(Exact Gaussianized Mutual Information):* Conditioning on all other factors $bold(z)_(-k, j)$, the residual channel $bold(r)_j = bold(lambda)_k z_(k j) + bold(epsilon)_j$ under $z_(k j) tilde.op cal(N)(0, xi_(k j))$ and $bold(epsilon)_j tilde.op cal(N)(bold(0), bold(Psi))$ has exact Shannon mutual information
     $ I(bold(r)_j; z_(k j)) = (1 / 2) log(1 + xi_(k j) bold(lambda)_k^top bold(Psi)^(-1) bold(lambda)_k). $
]

#proof[
  *Part 1:* Sufficiency is immediate since signed permutations preserve $ell_1$ norm. Conversely, suppose $Q in O(p)$ satisfies $norm(Q^top bold(z))_1 = norm(bold(z))_1$ for all $bold(z) in bb(R)^p$. Evaluating at each standard basis vector $bold(z) = bold(e)_r$ gives $sum_(k=1)^p |Q_(r k)| = 1$. Orthogonality of $Q$ simultaneously gives $sum_(k=1)^p Q_(r k)^2 = 1$, implying $|Q_(r k)| <= 1$ for all $r, k$. Subtracting yields
  $ sum_(k=1)^p |Q_(r k)| (1 - |Q_(r k)|) = 0. $
  Because every summand is non-negative, $|Q_(r k)| in {0, 1}$ for all $r, k$. Combined with $sum_k Q_(r k)^2 = sum_r Q_(r k)^2 = 1$, each row and column contains exactly one non-zero entry equal to $plus.minus 1$, i.e., $Q in cal(P)_(plus.minus) (p)$.

  *Part 2:* Rescaling factor $k$ maps $(bold(lambda)_k, bold(z)_k) |-> (c_k bold(lambda)_k, c_k^(-1) bold(z)_k)$, leaving $bold(Lambda) bold(Z)$ unchanged while modifying the $k$-th joint prior penalty to $h_k (c_k) = alpha norm(bold(lambda)_k)_1 c_k + sqrt(2) norm(bold(z)_k)_1 c_k^(-1)$. Differentiating twice gives $h''_k (c_k) = 2 sqrt(2) norm(bold(z)_k)_1 c_k^(-3) > 0$ on $(0, infinity)$, so $h_k$ is strictly convex with unique global minimum at $h'_k (c_k^star) = 0 <==> c_k^star = (sqrt(2) norm(bold(z)_k)_1 / (alpha norm(bold(lambda)_k)_1))^(1 / 2)$, uniquely balancing $alpha norm(c_k^star bold(lambda)_k)_1 = sqrt(2) norm((c_k^star)^(-1) bold(z)_k)_1$.

  *Part 3:* Under $z_(k j) tilde.op cal(N)(0, xi_(k j))$ independent of $bold(epsilon)_j tilde.op cal(N)(bold(0), bold(Psi))$, the marginal covariance of $bold(r)_j = bold(lambda)_k z_(k j) + bold(epsilon)_j$ is $bold(Sigma)_(r_j) = bold(Psi) + xi_(k j) bold(lambda)_k bold(lambda)_k^top$. By the matrix determinant lemma,
  $ det(bold(Psi) + xi_(k j) bold(lambda)_k bold(lambda)_k^top) = det(bold(Psi)) (1 + xi_(k j) bold(lambda)_k^top bold(Psi)^(-1) bold(lambda)_k). $
  Since differential entropy of $cal(N)(bold(0), bold(Sigma))$ is $(1 / 2) log((2 pi e)^n det bold(Sigma))$, the exact mutual information is $I(bold(r)_j; z_(k j)) = H(bold(r)_j) - H(bold(r)_j | z_(k j)) = (1 / 2) log(det bold(Sigma)_(r_j) / det bold(Psi)) = (1 / 2) log(1 + xi_(k j) bold(lambda)_k^top bold(Psi)^(-1) bold(lambda)_k)$. Summing over $j = 1, dots, m$ gives FABIA's factor-ranking score $I_k$ @hochreiter2010fabia.
]

#subsec("D", "Exact KKT Sparsity Threshold and Total Factor Collapse")

FABIA's $ell_1$ penalty induces exact zeros in $bold(Lambda)$, but over-penalization or small sample count $m$ can trigger complete collapse to the trivial zero solution.

#proposition(title: "Exact Zero-Loading Criterion and Absorbing Factor Collapse", origin: "Original Derivation")[
  At any M-step of Algorithm 161 with $A = sum_(j=1)^m bold(S)_j succ 0$ and $bold(b)_i = sum_(j=1)^m x_(i j) bold(mu)_(z_j)$:
  1. Loading entry $Lambda_(i k)^star$ is set to _exact zero_ whenever the partial residual satisfies
     $ | b_(i k) - sum_(ell != k) A_(k ell) Lambda_(i ell)^star | <= alpha psi_i, $
     demonstrating that noisier rows (larger $psi_i$) automatically enforce stronger sparsity thresholds.
  2. *(Absorbing Collapse Threshold):* If $alpha >= max_(1 <= i <= n) (norm(bold(b)_i)_infinity / psi_i)$ at any iteration $t$, then $bold(Lambda)^((t)) = bold(0)_(n times p)$ is the unique global optimum of the M-step, and $(bold(Lambda), M_Z) = (bold(0), bold(0))$ is an exact absorbing fixed point for all subsequent iterations $t' > t$.
]

#proof[
  The row-$i$ M-step minimizes $g_i (bold(lambda)) = (1 / (2 psi_i))(bold(lambda)^top A bold(lambda) - 2 bold(b)_i^top bold(lambda)) + alpha norm(bold(lambda))_1$. By convex subdifferential optimality, $bold(lambda)^star$ is the unique global minimizer iff $bold(0) in (1 / psi_i)(A bold(lambda)^star - bold(b)_i) + alpha partial norm(bold(lambda)^star)_1$. For coordinate $k$, since $partial |lambda_k^star| = [-1, 1]$ at $lambda_k^star = 0$, setting $lambda_k^star = 0$ satisfies coordinate-$k$ KKT optimality iff $|(A bold(lambda)^star - bold(b)_i)_k| <= alpha psi_i$, which reduces to $|b_(i k) - sum_(ell != k) A_(k ell) lambda_ell^star| <= alpha psi_i$.

  For Part 2, substituting $bold(lambda)^star = bold(0)$ into the KKT condition yields $norm(bold(b)_i)_infinity <= alpha psi_i$. If this inequality holds for every row $i = 1, dots, n$, strict convexity guarantees $bold(Lambda)^((t)) = bold(0)_(n times p)$ uniquely. In iteration $t + 1$, the E-step computes $bold(mu)_(z_j)^((t+1)) = bold(Sigma)_(z_j) bold(0)^top bold(Psi)^(-1) bold(x)_j = bold(0)$ for all $j$, yielding $bold(b)_i^((t+1)) = bold(0)$ and preserving $bold(Lambda)^((t')) = bold(0)$ inductively for all $t' > t$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each variational EM iteration requires three principal computational stages. First, forming $G = bold(Lambda)^top bold(Psi)^(-1) bold(Lambda)$ costs $O(n p^2)$, and inverting the $p times p$ positive definite precision matrix $bold(Xi)_j^(-1) + G$ for each of the $m$ samples via Cholesky factorization takes $O(m p^3)$, while computing all posterior means $M_Z$ via $bold(Lambda)^top bold(Psi)^(-1) X$ costs $O(p n m + m p^2)$. Second, updating $bold(Xi)_j$ takes $O(p m)$. Third, forming $A = sum_(j=1)^m bold(S)_j$ and $B = X M_Z^top$ costs $O(m p^2 + p n m)$, and solving the $n$ decoupled $p$-dimensional Lasso-type M-steps takes $O(n p^3)$ (or $O(n p^2)$ per coordinate-descent sweep, matching Hochreiter et al.'s sign-fixed linear solve @hochreiter2010fabia). Across $T$ iterations, total runtime is $O(T(p^3 (n + m) + p n m))$, scaling linearly in $n m$ whenever $p << min(n, m)$ (typical $p in [5, 20]$).

Working memory is $O(n m + p(n + m) + p^2 m)$ — or $O(n m + p(n + m) + p^2)$ when accumulating $A = sum_j bold(S)_j$ on the fly without storing all $m$ matrices $bold(Sigma)_(z_j)$ simultaneously — enabling cache-resident Cholesky solves for small $p$. Geometrically, FABIA models data as a union of overlapping sparse rank-$1$ rays in $bb(R)^n$ modulated by heavy-tailed sample activations, capturing both positive and anti-correlated (sign-flipped) co-regulation inside a single bicluster — a structural capability absent from binary Bimax @prelic2006bimax and non-negative NMF @lee1999nmf. However, because FABIA assumes zero-mean latent factors ($bb(E)[bold(z)_j] = bold(0)$), an uncentered global additive offset $mu_0 bold(1)_n bold(1)_m^top$ consumes a dense factor or inflates $bold(Psi)$ to $O(mu_0^2)$, making feature centering mandatory prior to fitting.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Input Domain], [Centered continuous real matrices $X in bb(R)^(n times m)$],
    [Time Complexity], [$O(T(p^3 (n + m) + p n m))$ ($O(T p n m)$ for small $p$)],
    [Space Complexity], [$O(n m + p(n + m) + p^2)$ with streaming $A$-accumulation],
    [Optimality], [Monotone variational ELBO ascent to local MAP stationarity],
    [Identifiability], [Unique up to signed permutations $cal(P)_(plus.minus) (p)$ and balanced scale],
    [Noise Model], [Heteroscedastic diagonal Gaussian $bold(Psi) = "diag"(psi_1, dots, psi_n)$],
    [Bicluster Polarity], [Multiplicative signed co-regulation (up- and down-regulated)],
    [Key Hyperparameters], [Factor count $p$; sparsity rate $alpha$; variance floor $psi_"min"$],
  ),
  caption: [Operational and complexity profile of FABIA @hochreiter2010fabia.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Toxicogenomics and Drug-Response Transcriptomics:* In drug-perturbation gene-expression studies (the flagship benchmark of @hochreiter2010fabia), compounds frequently up-regulate one gene subset while simultaneously down-regulating an antagonistic subset in proportion to dose. FABIA's multiplicative outer-product model $bold(lambda)_k bold(z)_k^top$ captures both opposing signs and continuous dose-response amplitudes while absorbing measurement noise per gene via $psi_i$.
2. *Cancer Subtype Stratification from High-Dimensional Omics:* Across heterogeneous tumor cohorts ($n >> m$), distinct molecular subtypes activate sparse overlapping gene signatures. FABIA's Laplace priors isolate subtype-specific biomarker modules while ranking genuine biological signals above technical artifacts via channel mutual information $I_k$.
3. *Haplotype and Rare-Variant IBD Segment Detection:* In population genetics (FABIA's HapFABIA extension), rare shared identity-by-descent (IBD) alleles form extremely sparse, high-dimensional biclusters across hundreds of thousands of SNPs and thousands of individuals, matching the sparse heavy-tailed prior structure.
4. *Financial Multi-Factor Alpha and Regime Attribution:* In equity return panels, heteroscedastic idiosyncratic risk ($psi_i$) varies drastically across stocks. FABIA separates sparse, transient sector-specific factors from stock-specific volatility far more accurately than homoscedastic PCA.
5. *Hyperspectral and Chemometric Spectral Unmixing:* When chemical mixtures exhibit sparse localized absorption bands with varying component concentrations across samples, FABIA recovers signed sparse component profiles and sample concentrations under band-dependent sensor noise $psi_i$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

FABIA exhibits four main failure regimes: (i) *low-dimensional / low-$m$ factor collapse* (Proposition 1): when $m$ is very small, inner-product magnitudes $norm(bold(b)_i)_infinity$ fail to exceed $alpha psi_i$, collapsing $bold(Lambda) -> bold(0)$; (ii) *uncentered additive offsets*, where non-zero global means waste a sparse factor on a dense constant vector; (iii) *nonlinear manifold geometry*, where curved non-linear clusters (e.g., concentric rings or spirals) violate the linear multiplicative subspace assumption $X approx bold(Lambda) bold(Z)$; and (iv) *high-$p$ cubic overhead*, where extracting $p > 100$ biclusters incurs $O(m p^3)$ Cholesky bottlenecks per EM iteration.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [High-D Signed Omics ($n, m >> p$)], [Exact multiplicative signed model + heteroscedastic $psi_i$], [Center data; tune $alpha$ below collapse],
    [Low-D Continuous Tables ($m <= 6$)], [Poor; small $norm(bold(b)_i)_infinity$ triggers $bold(Lambda) -> bold(0)$], [Use GMM, $k$-Means, or hierarchical],
    [Curved / Non-Convex Manifolds], [Poor; strictly linear rank-$p$ generative subspace], [Use DBSCAN, spectral, or HDBSCAN],
    [Exact Binary Incidence Mining], [Suboptimal; continuous Gaussian likelihood blurs $0/1$], [Use Bimax @prelic2006bimax],
    [Very Large $p$ ($p > 100$)], [Slow; $O(m p^3)$ per-sample precision inversion], [Prefer ISA @bergmann2003isa or sparse NMF],
  ),
  caption: [Domain suitability and failure-mode matrix for FABIA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

FABIA unifies classical factor analysis @bishop2006pattern, sparse matrix factorization @lee1999nmf @lee2001algorithms, and probabilistic biclustering models @lazzeroni2002plaid within a single parametric family.

#proposition(title: "Parametric Reductions to Factor Analysis, Sparse NMF, and Plaid Models", origin: "Original Proof")[
  Let $cal(L)_"FABIA"$ denote the variational objective of Algorithm 161.
  1. *(Gaussian Factor Analysis & PPCA):* Replacing the Laplace factor prior with $bold(z)_j tilde.op cal(N)(bold(0), I_p)$ (fixing $bold(Xi)_j equiv I_p$) and setting $alpha = 0$ eliminates variational relaxation error identically, reducing Algorithm 161 to exact EM for classical Gaussian Factor Analysis @bishop2006pattern; additionally constraining $bold(Psi) = sigma^2 I_n$ recovers Probabilistic PCA.
  2. *(Penalized NMF Equivalence):* Fixing homoscedastic noise $bold(Psi) = sigma^2 I_n$, taking the MAP point-estimate limit ($bold(Sigma)_(z_j) -> bold(0)$, $bold(mu)_(z_j) = bold(z)_j$), and restricting $bold(Lambda) >= 0, bold(Z) >= 0$ converts FABIA's negative log-posterior directly into the $ell_1$-regularized NMF objective @lee1999nmf @lee2001algorithms
     $ min_(bold(Lambda) >= 0, bold(Z) >= 0) (1 / (2 sigma^2)) norm(X - bold(Lambda) bold(Z))_F^2 + alpha norm(bold(Lambda))_1 + sqrt(2) norm(bold(Z))_1. $
  3. *(Plaid Model Subclass):* Restricting $Z_(k j) in {0, 1}$ and $Lambda_(i k) in {0, mu_k + a_(i k)}$ recovers the multiplicative layer structure of the Plaid biclustering model @lazzeroni2002plaid without global row/column background terms.
]

#proof[
  For (1), when $p(bold(z)_j) = cal(N)(bold(0), I_p)$, the exact posterior $p(bold(z)_j | bold(x)_j, bold(Lambda), bold(Psi))$ is already Gaussian with covariance $(I_p + bold(Lambda)^top bold(Psi)^(-1) bold(Lambda))^(-1)$ and mean $bold(Sigma)_(z_j) bold(Lambda)^top bold(Psi)^(-1) bold(x)_j$, matching Step 3 with $bold(Xi)_j = I_p$. With $alpha = 0$, Step 5 solves $bold(Lambda) = B A^(-1)$ in closed form, which is verbatim the classical Rubin-Thayer EM M-step for Factor Analysis @dempster1977em @bishop2006pattern; setting $bold(Psi) = sigma^2 I_n$ reduces the model to isotropic PPCA. For (2), substituting $bold(Psi) = sigma^2 I_n$ and point-mass posteriors $q_j (bold(z)_j) = delta(bold(z)_j - bold(mu)_(z_j))$ into $-log p(X, bold(Lambda), bold(Z) | bold(Psi))$ yields $(2 sigma^2)^(-1) sum_(j=1)^m norm(bold(x)_j - bold(Lambda) bold(z)_j)_2^2 + alpha norm(bold(Lambda))_1 + sqrt(2) norm(bold(Z))_1 + "const"$, which on the non-negative orthant is exact sparse NMF. Statement (3) follows directly by restricting the continuous support of $(bold(Lambda), bold(Z))$ to binary sample indicators and layer-specific row effects.
]

*Empirical Behavior on ClusterOverview Benchmarks:* On the ClusterOverview four-dataset suite evaluated with Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, FABIA exhibits a mean runtime of *6.321 ms*, noticeably heavier than Bimax (`0.215 ms`) and ISA (`0.401 ms`) due to iterative variational EM updates and $p times p$ matrix inversions per sample, yet faster than $k$-Shape (`14.685 ms`). Crucially, the benchmark metrics dramatically validate Proposition 1 (Absorbing Factor Collapse): on low-dimensional *Iris* (`Sil = N/A`, `ARI = 0.000`), *Customers* (`Sil = N/A`), and *Multi-Shape* (`Sil = N/A`, `ARI = N/A`), small feature dimensionality $m$ causes insufficient cross-sample inner products $norm(bold(b)_i)_infinity$ relative to $alpha psi_i$ (and non-linear manifold structure in Multi-Shape), collapsing FABIA's hard-converted factor assignments into a single trivial label (`ARI = 0.000`, `Sil = N/A`). In stark contrast, on *High-D* (`Sil = 0.154`, `ARI = N/A`) — where high feature dimensionality supplies abundant signal energy to overcome the KKT sparsity threshold $alpha psi_i$ and heteroscedastic noise variances $psi_i$ separate signal dimensions from irrelevant noise dimensions — FABIA achieves the *highest High-D Silhouette (`0.154`) of all four algorithms in this group* (surpassing Bimax `0.095`, ISA `0.075`, and $k$-Shape `0.055`), confirming that FABIA's inductive bias is specifically tailored to high-dimensional sparse factor regimes.
