#import "../template/ieee.typ": *

#alg_entry(
  id: 123,
  name: "Spectral Biclustering (Kluger)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2003",
  paradigm: "Bistochastic Spectral Co-Clustering",
  complexity: $cal(O)(n m min(n, m))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical one-way clustering partitions observations across all features simultaneously, failing whenever distinct sample subsets exhibit coherent structure only within specific feature subsets. Biclustering @cheng2000biclustering addresses this limitation by jointly grouping rows and columns of a data matrix $A in bb(R)_(>= 0)^(n times m)$ into coherent submatrices. Kluger, Basri, Chang, and Gerstein @kluger2003spectral introduced *Spectral Biclustering*, formulating joint row-column partitioning as the discovery of a latent *checkerboard* block structure hidden beneath arbitrary multiplicative row and column nuisance scalings.

Specifically, suppose $A in bb(R)_(> 0)^(n times m)$ represents nonnegative measurements (originally expression levels of $n$ genes across $m$ clinical conditions @kluger2003spectral) governed by a row partition $cal(P)_r = {I_1, ..., I_(k_r)}$ and a column partition $cal(P)_c = {J_1, ..., J_(k_c)}$ with indicator matrices $Phi_r in {0, 1}^(n times k_r)$ and $Phi_c in {0, 1}^(m times k_c)$. In empirical matrices, raw row totals (e.g., baseline gene abundance) and column totals (e.g., array illumination or sequencing depth) distort the singular vectors of $A$. Kluger et al. @kluger2003spectral propose three spectral normalizations—*independent rescaling*, *bistochastic normalization*, and *log-interaction double centering*—to isolate piecewise-constant singular vectors aligned with $(cal(P)_r, cal(P)_c)$.

#definition(title: "Checkerboard Model and Kluger Normalizations")[
  A strictly positive matrix $A in bb(R)_(> 0)^(n times m)$ follows a *multiplicative checkerboard model* with partitions $(cal(P)_r, cal(P)_c)$ if
  $
    A = D_alpha (Phi_r M Phi_c^top) D_beta,
  $
  where $M in bb(R)_(> 0)^(k_r times k_c)$ is the latent block-intensity matrix and $D_alpha = "diag"(alpha_1, ..., alpha_n) succ 0$, $D_beta = "diag"(beta_1, ..., beta_m) succ 0$ are positive diagonal nuisance scalings. The three canonical Kluger normalizations are:
  1. *Independent rescaling*: $tilde(A)_"ind" = D_R^(-1/2) A D_C^(-1/2)$, where $D_R = "diag"(A bold(1)_m)$ and $D_C = "diag"(A^top bold(1)_n)$.
  2. *Bistochastic scaling*: $tilde(A)_"bs" = D_1 A D_2$, where positive diagonals $D_1, D_2$ satisfy $tilde(A)_"bs" bold(1)_m = c_r bold(1)_n$ and $tilde(A)_"bs"^top bold(1)_n = c_c bold(1)_m$ for constants $c_r, c_c > 0$ with $n c_r = m c_c$.
  3. *Log-interaction centering*: $tilde(A)_"log" = H_n L H_m$, where $L_(i j) = log A_(i j)$ and $H_n = I_n - (1 / n) bold(1)_n bold(1)_n^top$.
]

Under bistochastic scaling, the trivial top singular pair is constant, $u_1 = n^(-1/2) bold(1)_n$ and $v_1 = m^(-1/2) bold(1)_m$, with largest singular value $sigma_1 = sqrt(c_r c_c)$. The nontrivial singular vectors $(u_ell, v_ell)$ for $ell >= 2$ are orthogonal to the constant vectors and encode the latent row and column checkerboard partitions. Under independent rescaling, the rescaled vectors $D_R^(-1/2) u_ell$ and $D_C^(-1/2) v_ell$ solve the coupled generalized eigenproblems $A D_C^(-1) A^top x = lambda^2 D_R x$ and $A^top D_R^(-1) A y = lambda^2 D_C y$.

#algbox(
  title: "Spectral Biclustering (Kluger)",
  inputs: [Data matrix $A in bb(R)_(> 0)^(n times m)$, normalization mode in `{"ind", "bs", "log"}`, cluster counts $(k_r, k_c)$, spectral rank $p >= max(k_r, k_c)$, Sinkhorn tolerance $epsilon_"sk" > 0$.],
  outputs: [Row partition $cal(P)_r = {I_1, ..., I_(k_r)}$, column partition $cal(P)_c = {J_1, ..., J_(k_c)}$, and checkerboard blocks ${A_(I_a, J_b)}$.],
)[
  1. *Normalize matrix*:
     - If mode is `"ind"`, set $D_R = "diag"(A bold(1)_m)$, $D_C = "diag"(A^top bold(1)_n)$, and $tilde(A) = D_R^(-1/2) A D_C^(-1/2)$.
     - If mode is `"bs"`, iterate $B <- "diag"(B bold(1)_m)^(-1) B$ and $B <- B "diag"(B^top bold(1)_n)^(-1)$ via Sinkhorn-Knopp scaling until marginal error $< epsilon_"sk"$; set $tilde(A) = B$.
     - If mode is `"log"`, compute $L_(i j) = log A_(i j)$ and double-center $tilde(A) = H_n L H_m$.
  2. *Compute truncated SVD*: Obtain the top $p$ singular triplets ${(sigma_ell, u_ell, v_ell)}_(ell = 1)^p$ of $tilde(A)$ with $sigma_1 >= ... >= sigma_p >= 0$.
  3. *Remove trivial component*: If mode is `"ind"` or `"bs"`, discard $(u_1, v_1)$ and retain $U = [u_2, ..., u_p]$ and $V = [v_2, ..., v_p]$; if `"ind"`, rescale $tilde(U) = D_R^(-1/2) U$ and $tilde(V) = D_C^(-1/2) V$, else set $tilde(U) = U$, $tilde(V) = V$.
  4. *Project and cluster*: Score singular vectors by best piecewise-constant fit (or retain all $p - 1$ components), project rows $Z_r = tilde(A) tilde(V)$ and columns $Z_c = tilde(A)^top tilde(U)$ (or directly cluster rows of $tilde(U)$ and $tilde(V)$), and run $k$-means @lloyd1982least with $k_r$ row centroids and $k_c$ column centroids.
  5. *Extract checkerboard*: Return $(cal(P)_r, cal(P)_c)$ and reordered block submatrices $A_(I_a, J_b)$ for $a in {1, ..., k_r}$, $b in {1, ..., k_c}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational results for Kluger's Spectral Biclustering: exact invariance and piecewise-constant eigenvector recovery under multiplicative scaling, geometric contraction of Sinkhorn-Knopp preprocessing, subspace stability under additive noise, and constructive identifiability failures on rank-deficient or overlapping blocks.

#subsec("A", "Exact Checkerboard Recovery under Bistochastic Scaling")

A central mathematical advantage of bistochastic normalization over unnormalized SVD is its exact invariance to arbitrary positive row and column rescalings, which preserves the piecewise-constant structure of all nontrivial singular vectors.

#theorem(title: "Exact Piecewise-Constant Recovery under Multiplicative Nuisance Scaling", origin: "Original Proof")[
  Let $A = D_alpha (Phi_r M Phi_c^top) D_beta in bb(R)_(> 0)^(n times m)$, where $Phi_r in {0, 1}^(n times k_r)$ and $Phi_c in {0, 1}^(m times k_c)$ are partition indicator matrices with block sizes $n_a = abs(I_a) > 0$ and $m_b = abs(J_b) > 0$, $D_alpha, D_beta succ 0$ are arbitrary positive diagonal matrices, and $M in bb(R)_(> 0)^(k_r times k_c)$ has rank $r <= min(k_r, k_c)$. Let $B = D_1 A D_2$ be the bistochastic scaling satisfying $B bold(1)_m = (1 / n) bold(1)_n$ and $B^top bold(1)_n = (1 / m) bold(1)_m$. Then:
  1. $B$ is uniquely determined, strictly invariant to $(D_alpha, D_beta)$, and admits the exact block representation $B = Phi_r tilde(M) Phi_c^top$, where $tilde(M)_(a b) = e^(xi_a) M_(a b) e^(zeta_b)$ for unique (up to additive shift) block potentials $(xi, zeta) in bb(R)^(k_r) times bb(R)^(k_c)$.
  2. Every singular triplet $(sigma_ell, u_ell, v_ell)$ of $B$ with $sigma_ell > 0$ satisfies $u_ell in "span"(Phi_r)$ and $v_ell in "span"(Phi_c)$; that is, $u_ell$ is strictly piecewise constant on $cal(P)_r$ and $v_ell$ is strictly piecewise constant on $cal(P)_c$.
  3. If $M$ has distinct rows and distinct columns after block scaling and $r = max(k_r, k_c)$, the rows of $U_(2:r)$ and $V_(2:r)$ collapse to exactly $k_r$ and $k_c$ distinct points, respectively, yielding zero-error $k$-means recovery.
]

#proof[
  Consider the strictly convex dual entropic potential $Psi_A : bb(R)^n times bb(R)^m -> bb(R)$ defined by
  $
    Psi_A (x, y) = sum_(i = 1)^n sum_(j = 1)^m A_(i j) e^(x_i + y_j) - (1 / n) sum_(i = 1)^n x_i - (1 / m) sum_(j = 1)^m y_j .
  $
  Setting $nabla Psi_A (x, y) = 0$ yields the exact bistochastic marginal equations $sum_j A_(i j) e^(x_i + y_j) = 1 / n$ for all $i$ and $sum_i A_(i j) e^(x_i + y_j) = 1 / m$ for all $j$. Since $A_(i j) > 0$ everywhere, the Hessian of $Psi_A$ on the quotient space orthogonal to $(bold(1)_n, -bold(1)_m)$ is strictly positive definite: for any $(u, v)$ not proportional to $(bold(1)_n, -bold(1)_m)$,
  $
    [u^top, v^top] nabla^2 Psi_A (x, y) mat(u; v) = sum_(i, j) B_(i j) (u_i + v_j)^2 > 0 .
  $
  Hence the minimizer $(x^*, y^*)$ is unique up to $(x^* + c bold(1)_n, y^* - c bold(1)_m)$, and $B_(i j) = A_(i j) e^(x_i^* + y_j^*)$ is unique.

  Now let $A_0 = Phi_r M Phi_c^top$, so $A_(i j) = alpha_i (A_0)_(i j) beta_j$. Substituting $tilde(x)_i = x_i + log alpha_i$ and $tilde(y)_j = y_j + log beta_j$ gives
  $
    Psi_A (x, y) = Psi_(A_0) (tilde(x), tilde(y)) + (1 / n) sum_i log alpha_i + (1 / m) sum_j log beta_j .
  $
  Because the last two terms are independent of $(x, y)$, minimizing $Psi_A$ produces the identical matrix $B_(i j) = (A_0)_(i j) e^(tilde(x)_i^* + tilde(y)_j^*)$, proving exact invariance to $(D_alpha, D_beta)$.

  Next, let $Pi_r$ and $Pi_c$ be any permutations preserving blocks of $cal(P)_r$ and $cal(P)_c$. Since $Pi_r A_0 Pi_c^top = A_0$ and the marginal targets $(1 / n) bold(1)_n$, $(1 / m) bold(1)_m$ are permutation-invariant, uniqueness implies $Pi_r B Pi_c^top = B$. Consequently, $tilde(x)_i^* = xi_a$ for all $i in I_a$ and $tilde(y)_j^* = zeta_b$ for all $j in J_b$, giving $B = Phi_r tilde(M) Phi_c^top$ with $tilde(M)_(a b) = e^(xi_a) M_(a b) e^(zeta_b)$.

  For part (2), any singular vector pair $(u_ell, v_ell)$ with $sigma_ell > 0$ satisfies $u_ell = sigma_ell^(-1) B v_ell = Phi_r (sigma_ell^(-1) tilde(M) Phi_c^top v_ell) in "span"(Phi_r)$ and $v_ell = sigma_ell^(-1) B^top u_ell = Phi_c (sigma_ell^(-1) tilde(M)^top Phi_r^top u_ell) in "span"(Phi_c)$. Thus $u_ell$ is constant on each $I_a$ and $v_ell$ is constant on each $J_b$. Finally, for part (3), write $B = (Phi_r N_r^(-1/2)) (N_r^(1/2) tilde(M) N_c^(1/2)) (Phi_c N_c^(-1/2))^top$, where $N_r = "diag"(n_1, ..., n_(k_r))$ and $N_c = "diag"(m_1, ..., m_(k_c))$. Because $Phi_r N_r^(-1/2)$ and $Phi_c N_c^(-1/2)$ have orthonormal columns, the rows of $U_(1:r)$ take at most $k_r$ distinct values; since $u_1$ is constant, distinct rows of $tilde(M)$ map to $k_r$ distinct points in $bb(R)^(r - 1)$, and $k$-means recovers $cal(P)_r$ and $cal(P)_c$ with zero error.
]

#subsec("B", "Entropic Descent and Hilbert Contraction of Bistochastic Scaling")

Because Kluger's bistochastic normalization relies on iterative Sinkhorn-Knopp scaling @cuturi2013sinkhorn, establishing its exact convergence rate directly governs preprocessing complexity.

#lemma(title: "KL Monotonicity and Geometric Hilbert Contraction of Sinkhorn Scaling", origin: "Original Derivation")[
  Let $A in bb(R)_(> 0)^(n times m)$ have dynamic range $kappa = (max_(i, j) A_(i j)) / (min_(i, j) A_(i j)) < infinity$. Consider alternating Sinkhorn-Knopp iterates $x^((t+1)), y^((t+1))$ updating $e^(-x_i^((t+1))) = n sum_j A_(i j) e^(y_j^((t)))$ and $e^(-y_j^((t+1))) = m sum_i A_(i j) e^(x_i^((t+1)))$, with scaled matrix $B^((t))_(i j) = A_(i j) e^(x_i^((t)) + y_j^((t)))$. Then:
  1. Every row or column scaling half-step monotonically decreases the dual potential $Psi_A (x, y)$, with exact drop at each row update equal to the marginal KL divergence $"KL"(p || r^((t))) >= 0$, where $p = (1 / n) bold(1)_n$ and $r^((t)) = B^((t)) bold(1)_m$.
  2. In the Hilbert projective metric $d_H (u, v) = log(max_(i, i') (u_i v_(i')) / (u_(i') v_i))$ on $bb(R)_(> 0)^n$, the positive scaling vectors $w^((t)) = exp(x^((t)))$ contract geometrically at every full iteration:
     $
       d_H (w^((t+1)), w^*) <= ((kappa - 1) / (kappa + 1))^2 d_H (w^((t)), w^*) .
     $
]

#proof[
  For part (1), hold $y = y^((t))$ fixed and let $r_i = sum_j A_(i j) e^(x_i + y_j)$ before the row update. Replacing $x_i$ with $x_i^+ = x_i - log(n r_i)$ sets the new row sums to $1 / n$. Evaluating the change in $Psi_A$ gives
  $
    Psi_A (x, y) - Psi_A (x^+, y) &= sum_(i = 1)^n (r_i - 1 / n) + (1 / n) sum_(i = 1)^n log(n r_i) \
                                  &= (sum_(i = 1)^n r_i - 1) - (1 / n) sum_(i = 1)^n log(1 / (n r_i)) .
  $
  After the preceding column step, $sum_i r_i = sum_j (1 / m) = 1$, so the first term vanishes and the second equals
  $
    Psi_A (x, y) - Psi_A (x^+, y) = sum_(i = 1)^n (1 / n) log((1 / n) / r_i) = "KL"(p || r) >= 0 ,
  $
  with equality if and only if $r = p$. Symmetrically, each column update drops $Psi_A$ by $"KL"(q || c) >= 0$ with $q = (1 / m) bold(1)_m$.

  For part (2), observe that $w^((t+1)) = (1 / n) (A z^((t)))^(-1)$ and $z^((t)) = (1 / m) (A^top w^((t)))^(-1)$ elementwise, where inversion $u |-> u^(-1)$ is an exact isometry of $d_H$ ($d_H (u^(-1), v^(-1)) = d_H (u, v)$). By Birkhoff's contraction theorem for positive linear maps, any strictly positive matrix $A$ satisfies $d_H (A u, A v) <= tanh(Delta(A) / 4) d_H (u, v)$, where the projective diameter is
  $
    Delta(A) = max_(i, i', j, j') log((A_(i j) A_(i' j')) / (A_(i j') A_(i' j))) <= log(kappa^2) = 2 log kappa .
  $
  Evaluating the contraction constant yields
  $
    tanh(Delta(A) / 4) <= tanh((log kappa) / 2) = (sqrt(kappa) - 1 / sqrt(kappa)) / (sqrt(kappa) + 1 / sqrt(kappa)) = (kappa - 1) / (kappa + 1) .
  $
  Since $Delta(A^top) = Delta(A)$ and $(w^*, z^*)$ is a fixed point, composing the column step $A^top$ and row step $A$ gives
  $
    d_H (w^((t+1)), w^*) <= tanh(Delta(A) / 4) d_H (z^((t)), z^*) <= ((kappa - 1) / (kappa + 1))^2 d_H (w^((t)), w^*) ,
  $
  establishing $cal(O)(log(1 / epsilon_"sk") / log((kappa + 1) / (kappa - 1)))$ iteration complexity.
]

#subsec("C", "Spectral Perturbation Stability under Additive Noise")

Real matrices deviate from exact checkerboard structure due to measurement noise $B = B_0 + E$. Using the Hermitian dilation trick together with the Davis-Kahan $\sin Theta$ bound @davis1970rotation @vonluxburg2007tutorial, we obtain an explicit finite-sample stability bound for both row and column subspaces.

#theorem(title: "Wedin Subspace Perturbation Bound for Noisy Checkerboard Matrices", origin: "Classical — Cited")[
  Let $B_0 = Phi_r tilde(M) Phi_c^top in bb(R)^(n times m)$ have exact rank $r$ with smallest nonzero singular value $sigma_r (B_0) > 0$ and singular subspaces $(U_0, V_0) in bb(R)^(n times r) times bb(R)^(m times r)$. Let $B = B_0 + E$ be an observed bistochastic matrix with perturbation spectral norm $norm(E)_2 < sigma_r (B_0) / 2$, and let $(hat(U), hat(V))$ be the top-$r$ left and right singular subspaces of $B$. Then the canonical angle matrices $Theta_r, Theta_c$ between $(U_0, hat(U))$ and $(V_0, hat(V))$ satisfy (@davis1970rotation, @vonluxburg2007tutorial):
  $
    sqrt(norm(sin Theta_r)_F^2 + norm(sin Theta_c)_F^2) <= (2 sqrt(r) norm(E)_2) / (sigma_r (B_0) - norm(E)_2) .
  $
  Furthermore, if $hat(U)$ is clustered by a $(1 + epsilon)$-approximate $k$-means solver and the minimum separation between distinct row clusters in $U_0$ is $delta_r > 0$, then the number of misclassified rows $N_"err"$ is bounded by
  $
    N_"err" <= (16 (2 + epsilon) r norm(E)_2^2) / (delta_r^2 (sigma_r (B_0) - norm(E)_2)^2) .
  $
]

#proof[
  Form the symmetric $(n + m) times (n + m)$ Jordan-Wielandt Hermitian dilation
  $
    cal(S)(B_0) = mat(0, B_0; B_0^top, 0), quad cal(S)(E) = mat(0, E; E^top, 0) .
  $
  The eigenvalues of $cal(S)(B_0)$ are $plus.minus sigma_1 (B_0), ..., plus.minus sigma_r (B_0)$ together with $n + m - 2r$ zeros, and the positive invariant subspace $W_0 in bb(R)^((n + m) times r)$ has orthonormal columns $w_ell = 2^(-1/2) [u_(0, ell)^top, v_(0, ell)^top]^top$. Note that $norm(cal(S)(E))_2 = norm(E)_2$. By Weyl's inequality, every singular value of $B$ in the top-$r$ block satisfies $hat(sigma)_ell >= sigma_r (B_0) - norm(E)_2 > sigma_r (B_0) / 2$, whereas all remaining eigenvalues of $cal(S)(B_0)$ lie in $(-infinity, 0]$. Thus the spectral gap separating $hat(W)$ from the complementary spectrum of $cal(S)(B_0)$ is at least $gamma = sigma_r (B_0) - norm(E)_2 > 0$.

  Applying the Davis-Kahan $\sin Theta$ theorem @davis1970rotation to $cal(S)(B_0)$ and $cal(S)(B) = cal(S)(B_0) + cal(S)(E)$ yields
  $
    norm(sin Theta(W_0, hat(W)))_F <= norm(cal(S)(E) hat(W))_F / gamma <= (sqrt(r) norm(E)_2) / (sigma_r (B_0) - norm(E)_2) .
  $
  Now examine the orthogonal projector difference $P_(hat(W)) - P_(W_0)$: its upper-left $n times n$ block is $(1 / 2)(hat(U) hat(U)^top - U_0 U_0^top)$ and its lower-right $m times m$ block is $(1 / 2)(hat(V) hat(V)^top - V_0 V_0^top)$. Since $norm(P_(hat(W)) - P_(W_0))_F^2 = 2 norm(sin Theta(W_0, hat(W)))_F^2$, bounding the sum of squared Frobenius norms of the two diagonal blocks by the total squared Frobenius norm gives
  $
    (1 / 4) (2 norm(sin Theta_r)_F^2 + 2 norm(sin Theta_c)_F^2) <= 2 norm(sin Theta(W_0, hat(W)))_F^2 ,
  $
  which immediately simplifies to the stated $\sin Theta$ bound.

  For the misclassification count, there exists an orthogonal rotation $Q in bb(O)(r)$ minimizing $norm(hat(U) Q - U_0)_F <= sqrt(2) norm(sin Theta_r)_F$. Since $k$-means is orthogonally invariant, the optimal $k$-means cost on $hat(U)$ is bounded by $norm(hat(U) Q - U_0)_F^2$, and any row whose assigned centroid lies closer to a wrong population centroid than to its own contributes at least $(delta_r / 2)^2$ to the squared deviation $norm(hat(C) - U_0)_F^2 <= (4 + 2 epsilon) norm(hat(U) Q - U_0)_F^2$ via the triangle inequality, yielding the stated bound on $N_"err"$.
]

#subsec("D", "Constructive Identifiability Failures")

Despite its elegance on full-rank disjoint checkerboards, Kluger's spectral biclustering suffers structural blind spots that no sample size can overcome.

#proposition(title: "Rank-1 Multiplicative Collapse and Overlapping-Block Mixing", origin: "Original Counterexample")[
  There exist strictly positive block-structured matrices on which Kluger's bistochastic spectral biclustering provably fails:
  1. *Rank-1 block annihilation*: If $A = Phi_r M Phi_c^top$ has a rank-1 block matrix $M_(a b) = p_a q_b$ with distinct intensities $p_1 != p_2$ and $q_1 != q_2$, then bistochastic normalization collapses $A$ identically to the rank-1 uniform matrix $B = (1 / (n m)) bold(1)_n bold(1)_m^top$ with $sigma_2 (B) = 0$, erasing all bicluster structure.
  2. *Overlapping-bicluster subspace mixing*: Let $A in bb(R)_(> 0)^(4 times 4)$ be formed by two overlapping $3 times 3$ elevated biclusters $B_1 = {1, 2, 3} times {1, 2, 3}$ and $B_2 = {2, 3, 4} times {2, 3, 4}$ on background $1$, namely $A_(i j) = 1 + gamma (bb(I)[(i, j) in B_1] + bb(I)[(i, j) in B_2])$ with $gamma > 0$. Then $B$ has $sigma_2 (B) > 0$, yet its second singular vector assigns value $0$ to the overlapping rows ${2, 3}$ while separating only ${1}$ from ${4}$; consequently, any 2-cluster partition of $u_2$ misassigns at least one overlapping vertex from its true multi-cluster membership.
]

#proof[
  For part (1), set $D_1 = "diag"(Phi_r (N_r^(-1) p^(-1)))$ scaled by $1 / (k_r k_c)$? More directly, notice $A = (Phi_r p) (Phi_c q)^top = D_alpha (bold(1)_n bold(1)_m^top) D_beta$ where $alpha_i = p_a$ for $i in I_a$ and $beta_j = q_b$ for $j in J_b$. By Theorem 1, bistochastic scaling is invariant to $(D_alpha, D_beta)$, so $B$ equals the bistochastic scaling of $bold(1)_n bold(1)_m^top$, which is $(1 / (n m)) bold(1)_n bold(1)_m^top$. All singular values $sigma_ell (B)$ for $ell >= 2$ vanish identically, so no nontrivial spectral embedding exists.

  For part (2), observe that $A$ is symmetric and centrosymmetric under the reversal involution $J$ swapping $1 <-> 4$ and $2 <-> 3$: $J A J = A$. By uniqueness of bistochastic scaling, $B = D_1 A D_1$ is also symmetric and centrosymmetric ($J B J = B$). Explicitly, writing $A$ as a $4 times 4$ matrix with entries $A_(1, 4) = A_(4, 1) = 1$, $A_(1, 1..3) = A_(2..3, 1) = A_(2..3, 4) = A_(4, 2..4) = 1 + gamma$, and $A_(2..3, 2..3) = 1 + 2 gamma$, rows $2$ and $3$ are identical, so $B$ has rank at most $3$. The antisymmetric subspace ${x : J x = -x}$ is spanned by $e_1 - e_4$ and $e_2 - e_3$; because rows $2$ and $3$ are identical, $B (e_2 - e_3) = 0$. Therefore the unique nonzero antisymmetric eigenvector is $u_2 = 2^(-1/2) (1, 0, 0, -1)^top$ with strictly positive eigenvalue $B_(1, 1) - B_(1, 4) > 0$, while the symmetric subspace is spanned by $u_1 = (1 / 2) bold(1)_4$ and $u_3 = (1 / 2) (1, -1, -1, 1)^top$. Thus $u_2$ places overlapping rows ${2, 3}$ strictly at $0$, forcing any disjoint 2-way $k$-means clustering of $u_2$ to merge ${2, 3}$ into a single side or split identical rows arbitrarily, failing to recover the true overlapping cover ${B_1, B_2}$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Let $N_"nz" = "nnz"(A)$ denote the number of nonzero entries (with $N_"nz" = n m$ for dense strictly positive matrices). Each Sinkhorn-Knopp iteration requires one row-sum pass and one column-sum pass in $cal(O)(N_"nz")$ time, totaling $cal(O)(T_"sk" N_"nz")$, where by Lemma 1 $T_"sk" = cal(O)(kappa log(1 / epsilon_"sk"))$ in the worst case and typically $T_"sk" <= 50$ on well-conditioned matrices. Computing a full SVD takes $cal(O)(n m min(n, m))$ time and $cal(O)(n m)$ space; when $p << min(n, m)$, truncated Lanczos or randomized SVD reduces the spectral step to $cal(O)(p N_"nz" + p^2 (n + m))$ operations. Post-hoc $k$-means over $U in bb(R)^(n times (p - 1))$ and $V in bb(R)^(m times (p - 1))$ requires $cal(O)(T_"km" p (n k_r + m k_c))$ time, which is negligible compared to dense SVD.

*Geometric and Cache Behavior.* Geometrically, bistochastic scaling projects $A$ onto the Birkhoff-type transportation polytope of fixed marginals via iterative Bregman projections in Kullback-Leibler divergence @cuturi2013sinkhorn, whereas log-interaction normalization projects $log A$ orthogonally onto the subspace orthogonal to additive row and column effects. Dense SVD leverages BLAS-3 blocked Householder bidiagonalization with high arithmetic intensity and predictable L2/L3 cache locality, though memory scales quadratically ($cal(O)(n m)$) for dense arrays.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [Dense SVD: $cal(O)(n m min(n, m) + T_"sk" n m)$; truncated rank-$p$: $cal(O)((T_"sk" + p) "nnz"(A))$.],
    [Space complexity], [$cal(O)(n m)$ dense; $cal(O)("nnz"(A) + p (n + m))$ sparse.],
    [Cluster geometry], [Global disjoint checkerboard $cal(P)_r times cal(P)_c$ ($k_r$ row blocks, $k_c$ column blocks).],
    [Overlap / noise], [Strictly disjoint partition of rows and columns; no native background/outlier exclusion.],
    [Scale invariance], [Exact invariance to positive row and column scalings $D_alpha A D_beta$ under `"bs"` mode.],
    [Hyperparameters], [Row/column counts $(k_r, k_c)$, normalization mode, spectral rank $p$, Sinkhorn tolerance $epsilon_"sk"$.],
    [Determinism], [Spectral step deterministic (up to subspace sign/rotation); post-hoc $k$-means requires seeding.],
  ),
  caption: [Operational and theoretical profile of Kluger Spectral Biclustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Microarray and Bulk Transcriptomics Co-Expression.* Originally designed for cancer microarray cohorts @kluger2003spectral, gene-by-patient expression matrices exhibit massive technical nuisance factors (probe affinity per gene, RNA amplification per patient) layered atop subtype-specific gene modules. Bistochastic scaling eliminates both multiplicative factors simultaneously (Theorem 1), isolating disease-subtype checkerboards.
2. *Document-Term Matrix Co-Clustering.* In information retrieval and text mining @dhillon2001concept, term-frequency matrices vary dramatically in document length and inverse term frequency. Independent and bistochastic spectral normalizations map co-occurring word groups and topical document clusters into shared piecewise-constant singular vectors, outperforming unnormalized LSA when $k_r != k_c$.
3. *Collaborative Filtering and User-Item Interaction Blocks.* E-commerce click, purchase, and rating-count matrices exhibit heavy-tailed user activity and item popularity power laws ($D_alpha, D_beta$). Factoring out user/item volume biases via Sinkhorn scaling uncovers latent niche preference blocks matching customer segments to product categories.
4. *Metabolomics and Proteomics Cohort Stratification.* Mass-spectrometry intensity matrices contain strong ionization efficiency differences across analytes and sample-dilution offsets across biological specimens. Log-interaction double centering removes multiplicative ionization and dilution factors linearly in log space, recovering phenotype-specific biochemical pathways.
5. *Multi-Asset Financial Factor Regime Discovery.* Nonnegative volatility-spillover or absolute-correlation matrices between assets and macroeconomic regimes often admit rectangular block structures ($k_r$ asset sectors across $k_c$ market regimes); spectral biclustering jointly segments sectors and regimes while normalizing away baseline asset volatility.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Kluger's spectral biclustering fails in four well-characterized regimes: (i) *non-positive or signed data* without shifting, where Sinkhorn-Knopp scaling is undefined or diverges if support lacks total support; (ii) *sparse isolated or overlapping biclusters* embedded in a noisy background (Proposition 1), where disjoint global checkerboards force background rows into spurious blocks—here Plaid @lazzeroni2002plaid, ISA @bergmann2003isa, or FABIA @hochreiter2010fabia are preferable; (iii) *metric Euclidean clustering* where cluster separation lives in overall feature magnitude or nonlinear manifolds rather than row-column interactions, because double scaling destroys radial distance cues; and (iv) *rank-1 proportional block matrices* where bistochastic scaling annihilates all signal.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Gene expression / omics], [Strong: removes multiplicative probe and sample bias.], [Use FABIA @hochreiter2010fabia or Plaid @lazzeroni2002plaid for overlapping pathways.],
    [Document-term corpora], [Strong: normalized bipartite word-topic blocks.], [Use Dhillon @dhillon2001concept when $k_r = k_c$ bijection holds.],
    [User-item interaction], [Good: eliminates popularity and activity skew.], [Sparse zeros break Sinkhorn; add small pseudocount.],
    [Euclidean tabular data], [Weak: scaling erases isotropic/radial separation.], [Prefer $k$-means @lloyd1982least, GMM, or spectral clustering @ng2002spectral.],
    [Non-convex manifolds], [Weak: linear SVD misses geodesic geometry.], [Prefer HDBSCAN @campello2013hdbscan or graph-cut methods.],
  ),
  caption: [Domain suitability matrix for Kluger Spectral Biclustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Kluger's Spectral Biclustering sits at the intersection of bipartite spectral graph theory @chung1997spectral @vonluxburg2007tutorial, normalized graph cuts @shi2000ncut, and matrix-factorization biclustering @cheng2000biclustering @lee1999nmf. The following unified result formalizes its exact algebraic bridge to Dhillon's bipartite spectral co-clustering @dhillon2001concept.

#proposition(title: "Algebraic Equivalence with Bipartite Normalized Cut and Dhillon Co-Clustering", origin: "Original Unified Proof")[
  Let $G_"bip" = (V_r union V_c, E, A)$ be the weighted bipartite graph with adjacency matrix $W = mat(0, A; A^top, 0)$ and degree matrix $D = "diag"(D_R, D_C)$, where $D_R = "diag"(A bold(1)_m)$ and $D_C = "diag"(A^top bold(1)_n)$. Then $(sigma_ell, u_ell, v_ell)$ is a singular triplet of Kluger's independently rescaled matrix $tilde(A)_"ind" = D_R^(-1/2) A D_C^(-1/2)$ if and only if
  $
    z_ell^(plus.minus) = 2^(-1/2) mat(D_R^(-1/2) u_ell; plus.minus D_C^(-1/2) v_ell)
  $
  are generalized eigenvectors of the bipartite Normalized Cut Laplacian $(D - W) z = lambda D z$ @shi2000ncut with paired eigenvalues $lambda_ell^(plus.minus) = 1 minus.plus sigma_ell$. Consequently, Dhillon's spectral co-clustering @dhillon2001concept uses the single-cut bijective embedding ($k_r = k_c$) derived from $z_ell^(+)$, whereas Kluger's method @kluger2003spectral decouples row and column clustering ($k_r$ need not equal $k_c$) and replaces single-step degree normalization with full bistochastic projection.
]

#proof[
  Consider the symmetric normalized Laplacian $L_"sym" = I_(n + m) - D^(-1/2) W D^(-1/2) = mat(I_n, -tilde(A)_"ind"; -tilde(A)_"ind"^top, I_m)$. Applying $L_"sym"$ to $w_ell^(plus.minus) = 2^(-1/2) [u_ell^top, plus.minus v_ell^top]^top$ and using $tilde(A)_"ind" v_ell = sigma_ell u_ell$ and $tilde(A)_"ind"^top u_ell = sigma_ell v_ell$ yields
  $
    L_"sym" w_ell^(plus.minus) = 2^(-1/2) mat(u_ell minus.plus sigma_ell u_ell; plus.minus v_ell - sigma_ell v_ell) = (1 minus.plus sigma_ell) w_ell^(plus.minus) .
  $
  Left-multiplying by $D^(1/2)$ and substituting $z_ell^(plus.minus) = D^(-1/2) w_ell^(plus.minus)$ gives $(D - W) z_ell^(plus.minus) = (1 minus.plus sigma_ell) D z_ell^(plus.minus)$. Conversely, any generalized eigenvector splits into $(u, v)$ blocks satisfying $tilde(A)_"ind" v = (1 - lambda) u$ and $tilde(A)_"ind"^top u = (1 - lambda) v$, proving exact equivalence. Minimizing the Shi-Malik Normalized Cut @shi2000ncut relaxes to taking the smallest nontrivial eigenvalues $1 - sigma_ell$, which are identically the largest nontrivial singular values $sigma_ell$ of $tilde(A)_"ind"$.
]

Beyond Dhillon @dhillon2001concept, Kluger's checkerboard model contrasts sharply with greedy residue-based biclustering (Cheng-Church @cheng2000biclustering), binary inclusion-maximal biclustering (BiMax @prelic2006bimax), additive overlapping Plaid models @lazzeroni2002plaid, and sparse factor analysis (FABIA @hochreiter2010fabia): Kluger enforces a *global, mutually exclusive checkerboard* solved in closed spectral form without local-search initialization traps during the embedding phase.

*Empirical Behavior on ClusterOverview Benchmarks.* On the 4-dataset ClusterOverview benchmark suite—evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing—Spectral Biclustering (Kluger) records *Iris* Sil $= -0.029$, ARI $= -0.012$; *Customers* Sil $= 0.108$; *Multi-Shape* Sil $= 0.709$ (ARI N/A); and *High-D* Sil $= -0.022$ (ARI N/A), with a fast mean runtime of $2.809$ ms. This strongly bimodal empirical profile aligns directly with our theoretical analysis. On *Iris* (four positively correlated morphological features) and *High-D* (isotropic Gaussian clusters in Euclidean space), cluster differences reside primarily in global row-scale and mean shifts rather than multiplicative feature-block interactions; Theorem 1 and Proposition 1 explain why bistochastic/log-interaction normalization strips out these global scale cues, yielding near-zero Silhouette ($-0.029$, $-0.022$) and chance-level ARI ($-0.012$). On *Customers*, modest co-variation between customer groups and spending categories produces a positive Silhouette ($0.108$). Conversely, on *Multi-Shape*, where feature coordinates exhibit sharp block-aligned spatial structure after normalization, the decoupled spectral projection achieves a high Silhouette of $0.709$ in under $3$ ms, illustrating that Kluger's algorithm excels specifically when the target structure aligns with matrix row-column interactions rather than pure Euclidean ball geometry.
