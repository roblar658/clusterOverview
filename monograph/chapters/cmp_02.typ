#import "../template/ieee.typ": *

#cmp_entry(
  roman: "II",
  title: "Formal Mathematical Equivalences, Asymptotic Limits, and Exact Cross-Family Reductions",
  scope_label: "16 Complete Cross-Family Equivalence & Reduction Proofs Across All 200 Algorithms",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

#sec("I", "Taxonomy of Cross-Family Mathematical Reductions")

A central illusion of the clustering literature—spanning six decades across statistics, numerical linear algebra, graph theory, computational geometry, information theory, and deep learning @jain2010data, @xu2005survey—is that its eight algorithmic paradigms operate on mutually incommensurable mathematical principles. In Parts I–VIII of this monograph, we analyzed 200 distinct algorithms under their native variational, probabilistic, topological, or algorithmic formulations. In this chapter, we dismantle the apparent boundaries separating those eight families by establishing *sixteen formal theorems with complete, self-contained mathematical proofs* demonstrating that flagship algorithms from disparate families are exact algebraic reparameterizations, zero-temperature thermodynamic limits, metric-space functorial images, or multiplicative distortion-preserving reductions of one another.

To prevent imprecise analogies from masquerading as mathematical identities, we first codify the exact hierarchy of four formal reduction relations proved throughout Sections II–VI.

#definition(title: "Formal Cross-Family Reduction Relations")[
  Let $cal(A)$ and $cal(B)$ denote two clustering procedures defined on input spaces $cal(I)_cal(A)$ and $cal(I)_cal(B)$, with associated feasible partition/parameter spaces $Omega_cal(A), Omega_cal(B)$ and objective functionals $J_cal(A), J_cal(B)$ (or fixed-point update operators $T_cal(A), T_cal(B)$):
  1. *Exact Objective Isomorphism ($cal(A) arrow.l.r.double cal(B)$):* There exists an explicit bijection (or quotient bijection modulo scale/permutation symmetries) $Psi: Omega_cal(A) -> Omega_cal(B)$ and constants $alpha > 0$, $beta in bb(R)$ independent of the clustering variables such that $J_cal(A)(theta) = alpha J_cal(B)(Psi(theta)) + beta$ for every $theta in Omega_cal(A)$. Hence global optima and KKT stationary points coincide identically.
  2. *Asymptotic / Zero-Temperature Limit ($cal(A)_(tau) ->^(tau -> tau_0) cal(B)$):* For a family of parameterized algorithms indexed by temperature, dispersion, or fuzzifier $tau$, there exist explicit normalization constants $a(tau) > 0, b(tau) in bb(R)$ such that $a(tau) J_(cal(A)_tau)(theta) + b(tau) -> J_cal(B)(theta)$ with quantitative convergence rates, and the one-step iterative map satisfies $norm(T_(cal(A)_tau)(theta) - T_cal(B)(theta)) -> 0$ away from measure-zero tie hyperplanes.
  3. *Constructive Functorial Reduction ($cal(A) = cal(B) compose Phi$):* There exists a deterministic transformation $Phi: cal(I)_cal(A) -> cal(I)_cal(B)$ (such as a positive semi-definite kernel shift, a mutual-reachability metric transform, a spectral filter, or an $ell_1$ self-expressive affinity graph) such that executing $cal(B)$ on $Phi(cal(X))$ produces the exact hierarchy, partition, or embedding of $cal(A)$ on $cal(X)$.
  4. *Multiplicative Coreset Reduction ($cal(A) <=_(1 plus.minus epsilon) cal(B)$):* A streaming or hierarchical compression operator produces a weighted summary $(S, w)$ of size $|S| << n$ satisfying $(1 - epsilon) J_cal(B)(cal(X), Q) <= J_cal(B)((S, w), Q) <= (1 + epsilon) J_cal(B)(cal(X), Q)$ uniformly over all candidate prototype sets $Q$, transferring any $c$-approximation algorithm on $(S, w)$ into a $((1 + epsilon)/(1 - epsilon)) c$-approximation on $cal(X)$.
]

Table 1 maps the six cross-family bridges developed across Sections II–VI, showing how every Part of the monograph connects into a single transitive mathematical network.

#figure(
  table(
    columns: (0.7fr, 1.35fr, 0.95fr),
    align: (left, left, left),
    [Sec. / Thms], [Connected Parts & Algorithms], [Reduction Type],
    [II (Thms. 1–2)], [Part I (Alg. 16) $arrow.l.r.double$ Part V (Algs. 103–105, 114, 132)], [Isomorphism & Spectral Filter],
    [III (Thms. 3–6)], [Part IV (Algs. 81, 89, 95–97) & Part VII (Alg. 151) $->$ Part I (Algs. 1, 7, 12)], [Zero-Temp Limit & Legendre Dual],
    [IV (Thms. 7–10)], [Part II (Algs. 31, 32, 34, 43) $arrow.l.r.double$ Part III (Algs. 56, 59) $arrow.l.r.double$ Part I (Alg. 1)], [Metric Functor, EM & Greedy Delta],
    [V (Thms. 11–12)], [Part IV (Algs. 84, 92) & Part V (Algs. 107, 108, 117, 131)], [Exact Objective Isomorphism],
    [VI (Thms. 13–14)], [Part VIII (Algs. 190, 195) & Part III (Alg. 78) $->$ Part I (Algs. 1, 12)], [OT Dual & $(1 plus.minus epsilon)$ Coreset],
    [VI (Thms. 15–16)], [Part VI (Alg. 150) & Part VIII (Algs. 181–184) $->$ Parts I, IV, V], [Dual Certificate & Variational EM],
  ),
  caption: [Architectural map of the 16 cross-family equivalence, limit, and reduction theorems proved in Comparison Chapter II.],
)

#sec("II", "Spectral, Kernel, and Random-Walk Equivalences")

Graph partitioning, kernel prototype clustering, and Markov random-walk community detection were developed in separate literatures—combinatorial optimization @hagen1992ratiocut, @shi2000ncut, reproducing kernel Hilbert spaces (RKHS) @scholkopf1998nonlinear, and discrete-time Markov chains @chung1997spectral, @pons2006walktrap. The first two theorems prove that these three families are governed by the exact same matrix-trace bilinear forms and spectral filters of the normalized graph Laplacian.

#subsec("A", "Weighted Kernel K-Means and Spectral Graph Cuts")

Let $G = (V, E, W)$ be an undirected weighted similarity graph on $n = |V|$ vertices with symmetric non-negative adjacency matrix $W = W^T in bb(R)_(>= 0)^(n times n)$, vertex degrees $d_i = sum_(j=1)^n W_(i j) > 0$, degree matrix $D = "diag"(d_1, dots, d_n)$, unnormalized Laplacian $L = D - W$, and symmetric normalized Laplacian
$ L_"sym" = I_n - D^(-1/2) W D^(-1/2) = D^(-1/2) L D^(-1/2). $
For a disjoint $k$-partition $cal(C) = {C_1, dots, C_k}$ of $V$ with non-empty clusters, recall the *cut*, *volume*, *Ratio Cut* (Alg. 105, @hagen1992ratiocut), and *Normalized Cut* (Alg. 104, @shi2000ncut) functionals:
$ &"cut"(C_c, overline(C)_c) := sum_(i in C_c) sum_(j in.not C_c) W_(i j), quad "vol"(C_c) := sum_(i in C_c) d_i, \
  &"Rcut"(cal(C)) := sum_(c=1)^k ("cut"(C_c, overline(C)_c)) / (|C_c|), \
  &"Ncut"(cal(C)) := sum_(c=1)^k ("cut"(C_c, overline(C)_c)) / ("vol"(C_c)). $
Separately, let $phi: cal(X) -> cal(H)$ be a feature map into a Hilbert space $cal(H)$ with Gram matrix $K_(i j) = chevron.l phi(x_i), phi(x_j) chevron.r_cal(H)$ ($K succ.eq 0$), and let $w = (w_1, dots, w_n)^T in bb(R)_(> 0)^n$ be positive observation weights with $D_w = "diag"(w_1, dots, w_n)$. *Weighted Kernel $k$-Means* (Alg. 16) minimizes
$ J_w (cal(C); K) := sum_(c=1)^k sum_(i in C_c) w_i norm(phi(x_i) - m_c)_cal(H)^2, $
where the weighted RKHS centroid of cluster $C_c$ is $m_c := s_c^(-1) sum_(j in C_c) w_j phi(x_j)$ and $s_c := sum_(j in C_c) w_j$.

#theorem(title: "Weighted Kernel K-Means <=> Normalized Cut & Ratio Cut", origin: "Classical — Cited")[
  Let $G = (V, E, W)$ be a weighted graph with $d_i > 0$ for all $i in {1, dots, n}$, and let $cal(C) = {C_1, dots, C_k}$ range over all partitions of $V$ into $k$ non-empty clusters @dhillon2004kernel.
  1. *Unified Trace Identity:* For any positive weight vector $w > 0$ and symmetric kernel matrix $K in bb(R)^(n times n)$, let $Z(cal(C)) in bb(R)^(n times k)$ be the weighted orthonormal indicator matrix with entries $Z_(i c) = sqrt(w_i / s_c) bb(1)_(i in C_c)$ ($Z^T Z = I_k$). Then
     $ J_w (cal(C); K) = "Tr"(D_w^(1/2) K D_w^(1/2)) - "Tr"(Z^T D_w^(1/2) K D_w^(1/2) Z). $
  2. *Exact Ratio Cut Isomorphism (Alg. 16 $arrow.l.r.double$ Alg. 105):* Set uniform weights $D_w = I_n$ and shifted Laplacian kernel $K_"RC" := sigma I_n - L$ with any shift $sigma >= lambda_"max" (L)$. Then $K_"RC" succ.eq 0$ is a valid Mercer kernel and
     $ J_1 (cal(C); K_"RC") = "Rcut"(cal(C)) + (n - k) sigma - "Tr"(L). $
  3. *Exact Normalized Cut Isomorphism (Alg. 16 $arrow.l.r.double$ Alg. 104):* Set degree weights $D_w = D$ and shifted normalized kernel $K_"NC" := sigma D^(-1) + D^(-1) W D^(-1)$ with any shift $sigma >= 1$. Then $K_"NC" succ.eq 0$ is a valid Mercer kernel and
     $ J_d (cal(C); K_"NC") = "Ncut"(cal(C)) + (n - k)(sigma + 1) - "Tr"(L_"sym"). $
  In both cases, the additive offsets are strictly independent of $cal(C)$, so minimizing Weighted Kernel $k$-Means is identically equivalent to minimizing Ratio Cut and Normalized Cut, respectively.
]

#proof[
  *Part 1 (Matrix Trace Reformulation):* Let $H in {0, 1}^(n times k)$ be the hard partition indicator matrix ($H_(i c) = 1$ iff $i in C_c$), so $H^T D_w H = "diag"(s_1, dots, s_k)$ and $Z = D_w^(1/2) H (H^T D_w H)^(-1/2)$. Observe that
  $ Z^T Z = (H^T D_w H)^(-1/2) (H^T D_w H) (H^T D_w H)^(-1/2) = I_k. $
  Expanding the squared RKHS norm in $J_w (cal(C); K)$ for each cluster $C_c$:
  $ &sum_(i in C_c) w_i norm(phi(x_i) - m_c)_cal(H)^2 \
    &= sum_(i in C_c) w_i K_(i i) - 2 sum_(i in C_c) w_i chevron.l phi(x_i), m_c chevron.r_cal(H) + (sum_(i in C_c) w_i) norm(m_c)_cal(H)^2 \
    &= sum_(i in C_c) w_i K_(i i) - (2 / s_c) sum_(i, j in C_c) w_i w_j K_(i j) + (s_c / s_c^2) sum_(i, j in C_c) w_i w_j K_(i j) \
    &= sum_(i in C_c) w_i K_(i i) - (1 / s_c) sum_(i, j in C_c) w_i w_j K_(i j). $
  Denote the $c$-th column of $H$ by $h_c$ and the $c$-th column of $Z$ by $z_c = s_c^(-1/2) D_w^(1/2) h_c$. Then
  $ (1 / s_c) sum_(i, j in C_c) w_i w_j K_(i j) &= (1 / s_c) h_c^T D_w K D_w h_c \
    &= (s_c^(-1/2) h_c^T D_w^(1/2)) (D_w^(1/2) K D_w^(1/2)) (s_c^(-1/2) D_w^(1/2) h_c) \
    &= z_c^T D_w^(1/2) K D_w^(1/2) z_c. $
  Summing over $c = 1, dots, k$ and using $sum_(c=1)^k sum_(i in C_c) w_i K_(i i) = sum_(i=1)^n w_i K_(i i) = "Tr"(D_w^(1/2) K D_w^(1/2))$ gives
  $ J_w (cal(C); K) = "Tr"(D_w^(1/2) K D_w^(1/2)) - "Tr"(Z^T D_w^(1/2) K D_w^(1/2) Z), $
  establishing Part 1. Note that this algebraic identity holds for any symmetric $K$, while $K succ.eq 0$ ensures an explicit Hilbert space representation $K = Phi^T Phi$ and guarantees monotone descent of kernel Lloyd iterations.

  *Part 2 (Ratio Cut Equivalence):* For every cluster indicator vector $h_c in {0, 1}^n$, the quadratic form of the unnormalized Laplacian satisfies
  $ h_c^T L h_c &= 1/2 sum_(i, j = 1)^n W_(i j) (h_(i c) - h_(j c))^2 \
    &= sum_(i in C_c) sum_(j in.not C_c) W_(i j) = "cut"(C_c, overline(C)_c). $
  Taking $D_w = I_n$ gives $s_c = |C_c|$ and $z_c = |C_c|^(-1/2) h_c$, whence
  $ "Rcut"(cal(C)) = sum_(c=1)^k (h_c^T L h_c) / (|C_c|) = sum_(c=1)^k z_c^T L z_c = "Tr"(Z^T L Z). $
  Because $L = D - W$ is symmetric PSD with spectrum contained in $[0, lambda_"max" (L)]$ (and $lambda_"max" (L) <= 2 max_i d_i$ by Gershgorin's circle theorem), any shift $sigma >= lambda_"max" (L)$ yields $K_"RC" = sigma I_n - L succ.eq 0$. Substituting $D_w = I_n$ and $K = K_"RC"$ into Part 1 and utilizing $"Tr"(Z^T Z) = "Tr"(I_k) = k$:
  $ J_1 (cal(C); K_"RC") &= "Tr"(sigma I_n - L) - "Tr"(Z^T (sigma I_n - L) Z) \
    &= n sigma - "Tr"(L) - k sigma + "Tr"(Z^T L Z) \
    &= "Rcut"(cal(C)) + (n - k) sigma - "Tr"(L). $

  *Part 3 (Normalized Cut Equivalence):* Take $D_w = D$, so $s_c = sum_(i in C_c) d_i = "vol"(C_c)$ and $z_c = "vol"(C_c)^(-1/2) D^(1/2) h_c$. Then
  $ z_c^T L_"sym" z_c &= (1 / "vol"(C_c)) h_c^T D^(1/2) (D^(-1/2) L D^(-1/2)) D^(1/2) h_c \
    &= (h_c^T L h_c) / ("vol"(C_c)) = ("cut"(C_c, overline(C)_c)) / ("vol"(C_c)). $
  Summing over $c = 1, dots, k$ produces $"Ncut"(cal(C)) = "Tr"(Z^T L_"sym" Z)$. Now form $K_"NC" = sigma D^(-1) + D^(-1) W D^(-1)$ and compute its degree-weighted similarity matrix:
  $ D^(1/2) K_"NC" D^(1/2) &= sigma I_n + D^(-1/2) W D^(-1/2) \
    &= (sigma + 1) I_n - (I_n - D^(-1/2) W D^(-1/2)) \
    &= (sigma + 1) I_n - L_"sym". $
  Since every eigenvalue of $L_"sym"$ lies in $[0, 2]$ (as $v^T L_"sym" v = 1/2 sum_(i, j) W_(i j) (v_i / sqrt(d_i) - v_j / sqrt(d_j))^2 <= sum_(i, j) W_(i j) (v_i^2 / d_i + v_j^2 / d_j) = 2 norm(v)_2^2$), choosing any $sigma >= 1$ implies $(sigma + 1) I_n - L_"sym" succ.eq (sigma - 1) I_n succ.eq 0$. By Sylvester's law of inertia with $D^(-1/2) succ 0$, $K_"NC" succ.eq 0$ is strictly positive semi-definite. Substituting into Part 1:
  $ J_d (cal(C); K_"NC") &= "Tr"((sigma + 1) I_n - L_"sym") \
    &quad - "Tr"(Z^T ((sigma + 1) I_n - L_"sym") Z) \
    &= n (sigma + 1) - "Tr"(L_"sym") - k (sigma + 1) + "Tr"(Z^T L_"sym" Z) \
    &= "Ncut"(cal(C)) + (n - k)(sigma + 1) - "Tr"(L_"sym"), $
  completing the proof.
]

#corollary(title: "Spectral Stiefel Relaxation & Exact Monotone Refinement", origin: "Original Proof")[
  Relaxing the discrete constraint $Z = D_w^(1/2) H (H^T D_w H)^(-1/2)$ in Theorem 1 to the continuous Stiefel manifold $cal(S)_(n, k) := {Z in bb(R)^(n times k) : Z^T Z = I_k}$ gives, by the Ky Fan variational principle,
  $ min_(Z in cal(S)_(n, k)) "Tr"(Z^T L_"sym" Z) = sum_(r=1)^k lambda_r (L_"sym") <= min_(cal(C)) "Ncut"(cal(C)), $
  attained by any orthogonal rotation $Z^* = U_k R$ ($R in O(k)$) of the bottom-$k$ eigenvectors $U_k = [u_1, dots, u_k]$ of $L_"sym"$. Consequently:
  1. Spectral Clustering NJW (Alg. 103, @ng2002spectral) solves the continuous relaxation of Theorem 1 and rounds $U_k$ back to a discrete partition.
  2. Executing Weighted Kernel $k$-Means (Alg. 16) with $K_"NC"$ ($sigma >= 1$) initialized from any spectral cut $cal(C)^((0))$ monotonically decreases the *exact discrete* $"Ncut"(cal(C))$ at every iteration without computing any eigenvectors.
]

#proof[
  Because $L_"sym" = sum_(r=1)^n lambda_r u_r u_r^T$ with $0 = lambda_1 <= lambda_2 <= dots <= lambda_n$, for any $Z in cal(S)_(n, k)$ let $p_r = norm(Z^T u_r)_2^2 in [0, 1]$. Then $sum_(r=1)^n p_r = "Tr"(Z^T Z) = k$, so $"Tr"(Z^T L_"sym" Z) = sum_(r=1)^n lambda_r p_r >= sum_(r=1)^k lambda_r$, with equality when $"span"(Z) = "span"(u_1, dots, u_k)$. Since every discrete indicator $Z(cal(C))$ lies in $cal(S)_(n, k)$, the Ky Fan minimum lower-bounds $min_cal(C) "Ncut"(cal(C))$. Finally, because $K_"NC" succ.eq 0$ for $sigma >= 1$, $J_d (cal(C); K_"NC")$ is a genuine weighted squared-norm objective in $cal(H)$, so each alternating Voronoi-and-centroid step of Weighted Kernel $k$-Means strictly reduces $J_d$, and hence strictly reduces $"Ncut"(cal(C))$ by Theorem 1(3).
]

#remark(title: "Role of the Diagonal Shift Sigma")[
  Although the objective identity in Theorem 1 holds for *every* $sigma in bb(R)$, choosing $sigma < 1$ in $K_"NC"$ can render $K_"NC"$ indefinite, breaking the Parallel-Axis Lemma in $cal(H)$ and causing Kernel $k$-Means iterations to oscillate or increase $"Ncut"$. Conversely, choosing excessively large $sigma >> 1$ adds $(sigma / d_i) (1 - 2 d_i / "vol"(C_c))$ to point-to-cluster squared distances, biasing assignments toward high-volume clusters and slowing local movement; the optimal operational choice is the minimal PSD shift $sigma = 1$ (or $sigma = max(0, -lambda_"min" (D^(-1/2) W D^(-1/2)))$).
]

#subsec("B", "Spectral NJW, Commute-Time Embedding, Walktrap, and PIC")

We next prove that four seemingly distinct graph algorithms—Ng–Jordan–Weiss Spectral Clustering (Alg. 103, @ng2002spectral), Random-Walk Commute-Time embedding @vonluxburg2007tutorial, Walktrap agglomerative diffusion clustering (Alg. 114, @pons2006walktrap), and Power Iteration Clustering (PIC, Alg. 132, @lin2010pic)—all operate on the *exact same* right eigenfunctions of the Markov transition operator $P = D^(-1) W$, differing exclusively in the spectral filter $g(lambda_r)$ applied to the Laplacian eigenvalues $lambda_r (L_"sym")$.

Assume $G = (V, E, W)$ is connected with volume $"vol"(G) = sum_(i=1)^n d_i$. Let $S := D^(-1/2) W D^(-1/2) = I_n - L_"sym"$ have orthonormal eigendecomposition $S = sum_(r=1)^n mu_r u_r u_r^T$ with $1 = mu_1 > mu_2 >= dots >= mu_n >= -1$, so $lambda_r := 1 - mu_r$ are the eigenvalues of $L_"sym"$ ($0 = lambda_1 < lambda_2 <= dots <= lambda_n <= 2$) and $u_1 = "vol"(G)^(-1/2) D^(1/2) bold(1)$. Define the $D$-orthonormal right eigenvectors of $P = D^(-1) W$ and $L_"rw" = I_n - P$:
$ v_r := D^(-1/2) u_r, quad v_r^T D v_s = delta_(r s), quad P v_r = mu_r v_r. $
We first prove a subtle pseudoinverse lemma resolving a common algebraic error between the Moore–Penrose pseudoinverse $L^+$ of the unnormalized Laplacian and the similarity-transformed pseudoinverse $D^(-1/2) L_"sym"^+ D^(-1/2)$.

#lemma(title: "Exact Resistance Equivalence of Normalized and Unnormalized Pseudoinverses", origin: "Original Proof")[
  Let $G = (V, E, W)$ be a connected weighted graph, let $L^+$ be the Moore–Penrose pseudoinverse of $L = D - W$, and define $M := D^(-1/2) L_"sym"^+ D^(-1/2) = sum_(r=2)^n lambda_r^(-1) v_r v_r^T$. Although $M != L^+$ whenever $G$ is non-regular, for *every* pair of vertices $i, j in V$ the effective resistance satisfies
  $ R_"eff" (i, j) := (e_i - e_j)^T L^+ (e_i - e_j) = (e_i - e_j)^T M (e_i - e_j). $
]

#proof[
  Since $G$ is connected, $"ker"(L) = "span"(bold(1))$ and $"im"(L) = bold(1)^perp$. Because $e_i - e_j in bold(1)^perp$, the vector $y := L^+ (e_i - e_j)$ satisfies $L y = (I_n - n^(-1) bold(1) bold(1)^T)(e_i - e_j) = e_i - e_j$. Next compute $L M$ using $L = D^(1/2) L_"sym" D^(1/2)$ and $L_"sym" L_"sym"^+ = I_n - u_1 u_1^T$:
  $ L M &= D^(1/2) L_"sym" L_"sym"^+ D^(-1/2) \
    &= D^(1/2) (I_n - "vol"(G)^(-1) D^(1/2) bold(1) bold(1)^T D^(1/2)) D^(-1/2) \
    &= I_n - "vol"(G)^(-1) (D bold(1)) bold(1)^T. $
  Applying $L M$ to $e_i - e_j$ annihilates the rank-1 term because $bold(1)^T (e_i - e_j) = 0$, yielding $L (M (e_i - e_j)) = e_i - e_j$. Setting $x := M (e_i - e_j)$, we obtain $L (x - y) = bold(0)$, hence $x - y = c bold(1)$ for some scalar $c in bb(R)$. Taking the inner product with $e_i - e_j$ gives
  $ (e_i - e_j)^T x - (e_i - e_j)^T y = c (e_i - e_j)^T bold(1) = 0, $
  proving $(e_i - e_j)^T L^+ (e_i - e_j) = (e_i - e_j)^T M (e_i - e_j)$.
]

#theorem(title: "Spectral NJW <=> Commute-Time, Walktrap, and PIC via Unified Spectral Filtering", origin: "Original Unified Proof")[
  Let $G = (V, E, W)$ be a connected weighted graph with Markov transition matrix $P = D^(-1) W$ and eigenpairs $(lambda_r, v_r)$ ($v_r = D^(-1/2) u_r$, $mu_r = 1 - lambda_r$). Every spectral/random-walk distance in the family takes the universal weighted Euclidean form
  $ d_g^2 (i, j) = sum_(r=2)^n g(lambda_r) (v_(r i) - v_(r j))^2 $
  under a method-specific non-negative spectral filter $g: (0, 2] -> bb(R)_(>= 0)$:
  1. *Commute-Time Distance:* The expected round-trip hitting time $C(i, j) = H(i, j) + H(j, i)$ satisfies $C(i, j) = d_(g_"CT")^2 (i, j)$ with exact inverse-eigenvalue filter $g_"CT" (lambda_r) = "vol"(G) / lambda_r$.
  2. *Walktrap $t$-Step Diffusion Distance (Alg. 114, @pons2006walktrap):* The Walktrap vertex dissimilarity $r_(i j)^2 (t) := sum_(ell=1)^n d_ell^(-1) (P^t_(i ell) - P^t_(j ell))^2$ satisfies $r_(i j)^2 (t) = d_(g_"WT")^2 (i, j)$ with exact exponential-decay filter $g_"WT" (lambda_r) = (1 - lambda_r)^(2 t) = mu_r^(2 t)$.
  3. *Spectral Clustering (Shi–Malik & NJW, Alg. 103, @ng2002spectral):* Truncated $L_"rw"$ embedding uses the ideal low-pass boxcar filter $g_"SM" (lambda_r) = bb(1)_(2 <= r <= k)$. Moreover, on any $k$-component block graph $G_0$ (or perturbed graph with $norm(L_"sym" - L_"sym"^(0))_2 = delta << lambda_(k+1)$), the NJW row-normalized embedding $tilde(U)_(i, :) = U_(i, :) / norm(U_(i, :))_2$ is an orthogonal transformation of the $L_"rw"$ embedding scaled by cluster volumes, mapping every cluster $C_c$ to a single orthogonal unit vector on $bb(S)^(k-1)$ up to perturbation $O(delta / lambda_(k+1))$.
  4. *Power Iteration Clustering (PIC, Alg. 132, @lin2010pic):* For any non-degenerate starting vector $x^((0)) in bb(R)^n$ with $c_r := v_r^T D x^((0))$, the $t$-step power iterate $x^((t)) = P^t x^((0))$ satisfies
     $ x^((t)) = c_1 bold(1) + sum_(r=2)^k c_r (1 - lambda_r)^t v_r + O((max_(r > k) |1 - lambda_r|)^t). $
     Hence in the spectral-gap plateau where $max_(r > k) |mu_r / mu_k|^t << 1$ and $|mu_k|^t$ remains $Omega(1)$, PIC embeds vertices into a one-dimensional random linear projection of the exact top-$k$ $L_"rw"$ eigenspace $"span"(v_1, dots, v_k)$.
]

#proof[
  *Part 1 (Commute-Time Filter):* By the classical electrical-network theorem (Chandra et al.; see @vonluxburg2007tutorial, @chung1997spectral), the commute time of the reversible Markov chain $P = D^(-1) W$ equals $C(i, j) = "vol"(G) R_"eff" (i, j) = "vol"(G) (e_i - e_j)^T L^+ (e_i - e_j)$. Applying Lemma 1 immediately gives
  $ C(i, j) &= "vol"(G) (e_i - e_j)^T (sum_(r=2)^n lambda_r^(-1) v_r v_r^T) (e_i - e_j) \
    &= sum_(r=2)^n ("vol"(G) / lambda_r) (v_(r i) - v_(r j))^2. $

  *Part 2 (Walktrap Diffusion Filter):* Using $P = D^(-1/2) S D^(1/2)$, the $t$-step transition matrix decomposes spectrally as
  $ P^t = D^(-1/2) S^t D^(1/2) = sum_(r=1)^n mu_r^t (D^(-1/2) u_r) (D^(1/2) u_r)^T = sum_(r=1)^n mu_r^t v_r (D v_r)^T. $
  Entrywise, $P^t_(i ell) = sum_(r=1)^n mu_r^t v_(r i) u_(r ell) sqrt(d_ell)$. For $r = 1$, $mu_1 = 1$, $v_(1 i) = "vol"(G)^(-1/2)$ is constant across all vertices $i$, and $u_(1 ell) sqrt(d_ell) = d_ell / sqrt("vol"(G))$, giving the stationary distribution $pi_ell = d_ell / "vol"(G)$. Hence the $r = 1$ term cancels identically in the difference:
  $ (P^t_(i ell) - P^t_(j ell)) / sqrt(d_ell) = sum_(r=2)^n mu_r^t (v_(r i) - v_(r j)) u_(r ell). $
  Squaring and summing over $ell = 1, dots, n$, then invoking orthonormality $sum_(ell=1)^n u_(r ell) u_(s ell) = delta_(r s)$:
  $ r_(i j)^2 (t) &= sum_(ell=1)^n (sum_(r=2)^n mu_r^t (v_(r i) - v_(r j)) u_(r ell))^2 \
    &= sum_(r, s = 2)^n mu_r^t mu_s^t (v_(r i) - v_(r j)) (v_(s i) - v_(s j)) (sum_(ell=1)^n u_(r ell) u_(s ell)) \
    &= sum_(r=2)^n mu_r^(2 t) (v_(r i) - v_(r j))^2 = sum_(r=2)^n (1 - lambda_r)^(2 t) (v_(r i) - v_(r j))^2. $
  Moreover, because Walktrap subsequently runs Ward's agglomerative method (Alg. 59) on the Euclidean coordinates $Phi_"WT" (i) = (mu_2^t v_(2 i), dots, mu_n^t v_(n i))^T$, Theorem 10 below shows that Walktrap performs greedy WCSS minimization directly in the $t$-step spectrally filtered $L_"rw"$ embedding!

  *Part 3 (NJW Alignment & Davis–Kahan Stability):* On an ideal $k$-component graph $G_0$ with connected components $C_1, dots, C_k$, the nullspace of $L_"sym"^(0)$ is spanned by the orthonormal indicator vectors $w_c = "vol"(C_c)^(-1/2) D^(1/2) h_c$ ($c = 1, dots, k$). Any orthonormal eigenbasis $U_k^(0) = [u_1^(0), dots, u_k^(0)]$ therefore satisfies $U_k^(0) = W_0 R$ for some orthogonal $R in O(k)$, where $W_0 = [w_1, dots, w_k]$. For any vertex $i in C_c$, the $i$-th row of $U_k^(0)$ is
  $ (U_k^(0))_(i, :) = sqrt(d_i / "vol"(C_c)) e_c^T R, $
  whose Euclidean norm is strictly $norm((U_k^(0))_(i, :))_2 = sqrt(d_i / "vol"(C_c))$. Dividing row $i$ by its norm cancels the degree factor $sqrt(d_i)$ identically—precisely mimicking the left multiplication $D^(-1/2) U_k$ of $L_"rw"$—and yields
  $ (tilde(U)^(0))_(i, :) = e_c^T R = sqrt("vol"(C_c)) dot (D^(-1/2) U_k^(0))_(i, :). $
  Thus every vertex in $C_c$ maps to the *exact* same unit vector $R^T e_c in bb(S)^(k-1)$, with pairwise inter-cluster Euclidean distance $sqrt(2)$. Under an inter-cluster perturbation $E = L_"sym" - L_"sym"^(0)$ with $norm(E)_2 = delta < lambda_(k+1)(L_"sym"^(0)) / 2$, the Davis–Kahan $sin Theta$ theorem @davis1970rotation guarantees an orthogonal matrix $Q in O(k)$ such that
  $ norm(U_k Q - U_k^(0))_F <= (sqrt(2 k) delta) / (lambda_(k+1)(L_"sym"^(0)) - delta), $
  preserving cluster separability under row normalization whenever $d_"min" / "vol"(C_"max")$ stays bounded away from zero.

  *Part 4 (PIC Spectral-Gap Truncation):* Because ${v_1, dots, v_n}$ forms a $D$-orthonormal basis of $bb(R)^n$, expand the initial vector as $x^((0)) = sum_(r=1)^n c_r v_r$ with $c_r = v_r^T D x^((0))$. Applying $P^t$ and using $P^t v_r = mu_r^t v_r = (1 - lambda_r)^t v_r$ together with $v_1 = "vol"(G)^(-1/2) bold(1)$ yields
  $ x^((t)) = P^t x^((0)) &= c_1 "vol"(G)^(-1/2) bold(1) + sum_(r=2)^k c_r (1 - lambda_r)^t v_r \
    &quad + sum_(r=k+1)^n c_r (1 - lambda_r)^t v_r. $
  Let $rho := max_(r > k) |1 - lambda_r| < |1 - lambda_k|$. Since the residual tail satisfies $norm(sum_(r > k) c_r mu_r^t v_r)_D <= rho^t norm(x^((0)))_D$, whenever $t$ is chosen in the window
  $ (log(1 / epsilon)) / (log(|mu_k| / rho)) <= t <= (log(1 / eta)) / (log(1 / |mu_k|)), $
  the high-frequency components are attenuated by a factor $epsilon$ relative to the $k$-th cluster mode while the $k$-th mode retains amplitude at least $eta$. On a nearly block-diagonal graph where each $v_r$ ($r <= k$) is approximately piecewise constant on $C_1, dots, C_k$, any linear combination $sum_(r=2)^k c_r mu_r^t v_r$ is likewise piecewise constant on the $k$ clusters (with distinct cluster values almost surely for random $x^((0))$), allowing 1D $k$-means on $x^((t)) / norm(x^((t)))_1$ to recover the spectral partition in $O(t |E|)$ time.
]

#remark(title: "Why Low-Pass Filtering Beats Raw Commute Time in Large Graphs")[
  Theorem 2 exposes a striking asymptotic pathology noted by von Luxburg et al. @vonluxburg2007tutorial: in large graphs with bounded spectral gap ($lambda_2 = Omega(1)$) and growing $n$, the tail $sum_(r > k) lambda_r^(-1) (v_(r i) - v_(r j))^2$ dominates $C(i, j) / "vol"(G)$, forcing $C(i, j) / "vol"(G) -> d_i^(-1) + d_j^(-1)$—which depends *solely* on local vertex degrees and loses all global cluster structure! By replacing $g_"CT" (lambda_r) = "vol"(G) / lambda_r$ with either the sharp cutoff $g_"SM" (lambda_r) = bb(1)_(r <= k)$ (NJW) or the exponential damping $g_"WT" (lambda_r) = (1 - lambda_r)^(2 t)$ (Walktrap/PIC), spectral and diffusion methods explicitly annihilate the degree-dominated high-frequency tail.
]


#sec("III", "Zero-Temperature Thermodynamic and Variational Limits of Partitioning")

Whereas Section II revealed exact algebraic isomorphisms among graph and kernel objectives, the connection between Part I (hard prototype partitioning), Part IV (probabilistic statistical mixtures), and Part VII (fuzzy soft partitioning) is governed by *statistical-mechanical temperature limits*. In each case, a probabilistic or fuzzy algorithm introduces an entropic or power-law smoothing parameter ($sigma^2$, $T$, $m - 1$, $kappa^(-1)$, or dispersion $tau$), and taking the zero-temperature limit collapses soft posterior responsibilities onto deterministic Voronoi indicators while converging the free-energy functional onto hard distortion minimization.

#subsec("A", "Isotropic GMM EM to Lloyd's K-Means as Sigma -> 0+")

Consider a $k$-component homoscedastic isotropic Gaussian Mixture Model (Alg. 81, @dempster1977em, @mclachlan2000finite, @bishop2006pattern) on $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ with shared variance parameter $sigma^2 > 0$, mixing weights $bold(pi) = (pi_1, dots, pi_k) in "int"(Delta_k)$ ($pi_j > 0$, $sum_(j=1)^k pi_j = 1$), and component means $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$:
$ p(x | bold(mu), bold(pi); sigma) = sum_(j=1)^k pi_j (2 pi sigma^2)^(-d / 2) exp(- (norm(x - mu_j)_2^2) / (2 sigma^2)). $
Its empirical log-likelihood is $ell(bold(mu), bold(pi); sigma) := sum_(i=1)^n log p(x_i | bold(mu), bold(pi); sigma)$. Given current parameters $(bold(mu), bold(pi))$, one EM iteration @dempster1977em computes the posterior responsibility matrix $Gamma(sigma) = (gamma_(i j)(sigma))$ and updated centroids $bold(mu)^"EM" (sigma)$:
$ &gamma_(i j) (sigma) = (pi_j exp(- norm(x_i - mu_j)_2^2 / (2 sigma^2))) / (sum_(ell=1)^k pi_ell exp(- norm(x_i - mu_ell)_2^2 / (2 sigma^2))), \
  &mu_j^"EM" (sigma) = (sum_(i=1)^n gamma_(i j) (sigma) x_i) / (sum_(i=1)^n gamma_(i j) (sigma)), quad pi_j^"EM" (sigma) = 1 / n sum_(i=1)^n gamma_(i j) (sigma). $
For any $bold(mu)$, define the hard Voronoi minimizer set $cal(A)_i (bold(mu)) := arg min_(j in {1, dots, k}) norm(x_i - mu_j)_2^2$. When $|cal(A)_i (bold(mu))| = 1$ with unique nearest centroid $c_i$, define the strict Voronoi separation margin
$ Delta_i (bold(mu)) := min_(ell != c_i) norm(x_i - mu_ell)_2^2 - norm(x_i - mu_(c_i))_2^2 > 0, $
and the dataset-wide margin $Delta(bold(mu)) := min_(1 <= i <= n) Delta_i (bold(mu))$. Let $C_j (bold(mu)) := {i : c_i = j}$ and $mu_j^"Lloyd" := |C_j (bold(mu))|^(-1) sum_(i in C_j (bold(mu))) x_i$ be the hard Lloyd partition and centroid update (Alg. 1, @lloyd1982least), with optimal WCSS $J^*(bold(mu)) := sum_(i=1)^n min_j norm(x_i - mu_j)_2^2$.

#theorem(title: "Laplace Zero-Temperature Limit: Isotropic GMM EM -> Lloyd's K-Means", origin: "Original Proof")[
  Let $cal(X) subset.eq bb(R)^d$ satisfy $R_"max" := max_i norm(x_i)_2 < infinity$, let $pi_"min" := min_j pi_j > 0$, and fix $bold(mu) in (bb(R)^d)^k$. Then as $sigma -> 0^+$:
  1. *Exponential Posterior Concentration (E-Step Limit):* For every point $x_i$ with unique nearest centroid $c_i$ ($Delta_i (bold(mu)) > 0$) and every $ell != c_i$,
     $ gamma_(i ell) (sigma) <= (pi_ell / pi_(c_i)) exp(- (Delta_i (bold(mu))) / (2 sigma^2)), $
     $ 1 - gamma_(i c_i) (sigma) <= ((1 - pi_"min") / pi_"min") exp(- (Delta_i (bold(mu))) / (2 sigma^2)). $
     On any tie hyperplane ($|cal(A)_i (bold(mu))| >= 2$), $gamma_(i j) (sigma)$ converges to the prior-weighted split $pi_j / (sum_(r in cal(A)_i) pi_r)$ for $j in cal(A)_i$ and $0$ for $j in.not cal(A)_i$.
  2. *Quantitative Centroid Convergence (M-Step Limit):* If $Delta(bold(mu)) > 0$ and every Voronoi cell is non-empty ($|C_j (bold(mu))| >= 1$), then for all sufficiently small $sigma > 0$ satisfying $n pi_"min"^(-1) exp(- Delta(bold(mu)) / (2 sigma^2)) <= 1 / 2$,
     $ norm(mu_j^"EM" (sigma) - mu_j^"Lloyd")_2 <= (4 n R_"max") / (|C_j (bold(mu))| pi_"min") exp(- (Delta(bold(mu))) / (2 sigma^2)). $
  3. *Exact Two-Sided Laplace Free-Energy Sandwich:* For *every* $bold(mu) in (bb(R)^d)^k$ (including tie configurations) and every $sigma > 0$,
     $ &J^*(bold(mu)) <= -2 sigma^2 ell(bold(mu), bold(pi); sigma) - n d sigma^2 log(2 pi sigma^2) \
       &quad <= J^*(bold(mu)) - 2 sigma^2 sum_(i=1)^n log pi_(c_i), $
     where $c_i in cal(A)_i (bold(mu))$. Hence $-2 sigma^2 ell(bold(mu), bold(pi); sigma) - n d sigma^2 log(2 pi sigma^2) -> J^*(bold(mu))$ uniformly on any compact subset with $pi_j >= epsilon_0 > 0$ at exact rate $O(sigma^2)$.
  4. *Finite-Horizon Trajectory Equivalence:* If Lloyd's algorithm initialized at $bold(mu)^((0))$ generates $T$ tie-free iterations with non-empty clusters ($min_(0 <= t < T) Delta(bold(mu)^((t))) > 0$), then there exists $sigma_0 > 0$ such that for all $sigma in (0, sigma_0)$, the EM iterates satisfy $norm(bold(mu)_"EM"^((t)) (sigma) - bold(mu)_"Lloyd"^((t)))_2 = O(exp(- c / (2 sigma^2)))$ and induce the *identical* hard MAP partition at every step $t = 1, dots, T$.
]

#proof[
  *Part 1 (E-Step Concentration):* Write $d_(i j)^2 := norm(x_i - mu_j)_2^2$ and $d_(i *)^2 := min_j d_(i j)^2$. Multiply numerator and denominator of $gamma_(i j) (sigma)$ by $exp(d_(i *)^2 / (2 sigma^2))$:
  $ gamma_(i j) (sigma) = (pi_j exp(- (d_(i j)^2 - d_(i *)^2) / (2 sigma^2))) / (sum_(ell=1)^k pi_ell exp(- (d_(i ell)^2 - d_(i *)^2) / (2 sigma^2))). $
  For every $r in cal(A)_i (bold(mu))$, $d_(i r)^2 - d_(i *)^2 = 0$, so the denominator is strictly bounded below by $sum_(r in cal(A)_i) pi_r$. When $cal(A)_i (bold(mu)) = {c_i}$, for any $ell != c_i$ we have $d_(i ell)^2 - d_(i *)^2 >= Delta_i (bold(mu))$, whence
  $ gamma_(i ell) (sigma) <= (pi_ell exp(- Delta_i (bold(mu)) / (2 sigma^2))) / (pi_(c_i)). $
  Summing over all $ell != c_i$ yields
  $ 1 - gamma_(i c_i) (sigma) = sum_(ell != c_i) gamma_(i ell) (sigma) <= ((1 - pi_(c_i)) / pi_(c_i)) exp(- (Delta_i (bold(mu))) / (2 sigma^2)) $
  which is at most $((1 - pi_"min") / pi_"min") exp(- Delta_i (bold(mu)) / (2 sigma^2))$. When $|cal(A)_i (bold(mu))| >= 2$, every term with $ell in.not cal(A)_i$ vanishes exponentially as $sigma -> 0^+$, while every term with $r in cal(A)_i$ remains constant at $pi_r$, giving the exact prior-weighted limit $pi_j / (sum_(r in cal(A)_i) pi_r)$ for $j in cal(A)_i$.

  *Part 2 (M-Step Centroid Bound):* Let $epsilon(sigma) := pi_"min"^(-1) exp(- Delta(bold(mu)) / (2 sigma^2))$. By Part 1, for each cluster $j$:
  - If $i in C_j (bold(mu))$, then $|gamma_(i j)(sigma) - 1| = 1 - gamma_(i j)(sigma) <= epsilon(sigma)$.
  - If $i in.not C_j (bold(mu))$, then $0 <= gamma_(i j)(sigma) <= epsilon(sigma)$.
  Denote $N_j (sigma) := sum_(i=1)^n gamma_(i j)(sigma)$ and $n_j := |C_j (bold(mu))| >= 1$. Then $|N_j (sigma) - n_j| <= sum_(i=1)^n |gamma_(i j)(sigma) - bb(1)_(i in C_j)| <= n epsilon(sigma)$. Whenever $n epsilon(sigma) <= 1 / 2 <= n_j / 2$, we have $N_j (sigma) >= n_j / 2$. Now decompose the centroid discrepancy:
  $ &mu_j^"EM" (sigma) - mu_j^"Lloyd" = 1 / (N_j (sigma)) sum_(i=1)^n (gamma_(i j)(sigma) - bb(1)_(i in C_j)) x_i \
    &quad + (1 / (N_j (sigma)) - 1 / n_j) sum_(i in C_j) x_i. $
  Taking Euclidean norms and applying the triangle inequality with $norm(x_i)_2 <= R_"max"$ and $norm(sum_(i in C_j) x_i)_2 <= n_j R_"max"$:
  $ norm(mu_j^"EM" (sigma) - mu_j^"Lloyd")_2 &<= (R_"max" / (N_j (sigma))) sum_(i=1)^n |gamma_(i j)(sigma) - bb(1)_(i in C_j)| \
    &quad + (|n_j - N_j (sigma)|) / (n_j N_j (sigma)) (n_j R_"max") \
    &<= (2 n epsilon(sigma) R_"max") / (N_j (sigma)) <= (4 n R_"max") / (n_j pi_"min") exp(- (Delta(bold(mu))) / (2 sigma^2)). $

  *Part 3 (Laplace Free-Energy Sandwich):* For each datum $x_i$, isolate the factor $(2 pi sigma^2)^(-d / 2) exp(- d_(i *)^2 / (2 sigma^2))$ from $p(x_i | bold(mu), bold(pi); sigma)$:
  $ p(x_i | bold(mu), bold(pi); sigma) = (2 pi sigma^2)^(-d / 2) exp(- d_(i *)^2 / (2 sigma^2)) S_i (sigma), $
  where $S_i (sigma) := sum_(j=1)^k pi_j exp(- (d_(i j)^2 - d_(i *)^2) / (2 sigma^2))$. Because every exponent satisfies $-(d_(i j)^2 - d_(i *)^2) / (2 sigma^2) <= 0$ and the $c_i$-th exponent equals $0$, we obtain the exact non-asymptotic bounds
  $ pi_(c_i) <= S_i (sigma) <= sum_(j=1)^k pi_j = 1. $
  Taking logarithms, multiplying by $-2 sigma^2$, and summing over $i = 1, dots, n$:
  $ &sum_(i=1)^n d_(i *)^2 + n d sigma^2 log(2 pi sigma^2) <= -2 sigma^2 ell(bold(mu), bold(pi); sigma) \
    &quad <= sum_(i=1)^n d_(i *)^2 + n d sigma^2 log(2 pi sigma^2) - 2 sigma^2 sum_(i=1)^n log pi_(c_i). $
  Subtracting $n d sigma^2 log(2 pi sigma^2)$ from all three sides and noting $sum_(i=1)^n d_(i *)^2 = J^*(bold(mu))$ yields Part 3. Since $-2 sigma^2 sum_(i=1)^n log pi_(c_i) <= 2 n sigma^2 log(1 / epsilon_0) -> 0$, convergence is uniform on ${bold(pi) : min_j pi_j >= epsilon_0}$.

  *Part 4 (Finite-Horizon Induction):* Let $Delta_"min" := min_(0 <= t < T) Delta(bold(mu)_"Lloyd"^((t))) > 0$. Because $bold(mu) |-> Delta(bold(mu))$ is continuous (as a finite minimum of continuous functions), there exists $rho > 0$ such that whenever $max_j norm(mu_j - mu_(j, "Lloyd")^((t)))_2 < rho$, the hard Voronoi partition $cal(C)(bold(mu))$ equals $cal(C)(bold(mu)_"Lloyd"^((t)))$ identically and $Delta(bold(mu)) >= Delta_"min" / 2$. Proceeding by induction on $t = 0, dots, T - 1$: at $t = 0$, $bold(mu)_"EM"^((0)) = bold(mu)_"Lloyd"^((0))$; if $bold(mu)_"EM"^((t))$ lies within $rho$ of $bold(mu)_"Lloyd"^((t))$ with mixing weights bounded below by $pi_"min"^((t)) >= 1 / (2 n)$ (since $|C_j| >= 1$ and $pi_j^"EM" = N_j / n -> |C_j| / n >= 1 / n$), Part 2 guarantees that $norm(mu_(j, "EM")^((t+1)) - mu_(j, "Lloyd")^((t+1)))_2 <= C_0 exp(- Delta_"min" / (4 sigma^2)) < rho$ for all sufficiently small $sigma < sigma_0$, preserving both the Voronoi partition and the exponential error bound across all $T$ steps.
]

#remark(title: "Fixed-Sigma Homoscedasticity vs. Free-Variance MLE Singularity")[
  Theorem 3 requires *fixed* (or common lower-bounded) variance $sigma^2$ sent parametrically to $0^+$. In an *unconstrained* heteroscedastic GMM (Alg. 81 with free covariance matrices $Sigma_j succ 0$), placing a single centroid $mu_1 = x_1$ and letting $Sigma_1 -> 0$ drives $ell -> +infinity$ while leaving all other $n - 1$ points arbitrarily assigned @bishop2006pattern. Thus $k$-means is the zero-temperature limit of *isotropic fixed-scale* GMM, whereas unconstrained GMM requires covariance regularization $Sigma_j succ.eq epsilon I_d$ or Bayesian Wishart priors (Alg. 82) to stay well-posed.
]

#subsec("B", "Fuzzy C-Means as m -> 1+ and Entropy-Regularized K-Means <=> Softmax EM")

Soft partitioning splits into two major regularization mechanisms in Parts IV and VII: Bezdek's power-law *Fuzzy C-Means (FCM)* (Alg. 151, @bezdek1981fcm, @hathaway1986fcm), which raises memberships to a fuzzifier exponent $m > 1$:
$ J_m (U, bold(mu)) := sum_(i=1)^n sum_(j=1)^k u_(i j)^m norm(x_i - mu_j)_2^2, quad U in Delta_k^n, $
and *Entropy-Regularized $k$-Means* (deterministic annealing / Kullback–Leibler regularized distortion) at temperature $T > 0$ with prior $bold(pi) in "int"(Delta_k)$:
$ cal(F)_T (U, bold(mu), bold(pi)) := sum_(i=1)^n sum_(j=1)^k u_(i j) norm(x_i - mu_j)_2^2 + T sum_(i=1)^n "KL"(u_i || bold(pi)), $
where $"KL"(u_i || bold(pi)) = sum_(j=1)^k u_(i j) log(u_(i j) / pi_j)$ (which equals $-H(u_i) + log k$ when $pi_j = 1 / k$). Theorem 4 proves both the exact algebraic zero-fuzzifier reduction of FCM to Lloyd's $k$-means (and its degenerate $m -> infinity$ collapse) and the exact variational isomorphism between Entropy-Regularized $k$-Means and Isotropic Softmax GMM EM.

#theorem(title: "FCM Zero-Fuzzifier Limit (m -> 1+) and Entropy-Regularized K-Means <=> Softmax EM", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and fix prototypes $bold(mu) = (mu_1, dots, mu_k)$ with squared distances $d_(i j)^2 = norm(x_i - mu_j)_2^2$.
  1. *FCM Hard-Partition Limit ($m -> 1^+$):* For any $m > 1$, the unique minimizer $U^*(m) = arg min_(U in Delta_k^n) J_m (U, bold(mu))$ satisfies, whenever $d_(i j) > 0$ for all $j$,
     $ u_(i j)^* (m) = (sum_(ell=1)^k (d_(i j)^2 / d_(i ell)^2)^(1 / (m - 1)))^(-1). $
     If $x_i$ has unique nearest centroid $c_i$ with distance ratio $rho_i := max_(ell != c_i) (d_(i c_i) / d_(i ell)) in [0, 1)$, then for every $ell != c_i$,
     $ u_(i ell)^* (m) <= rho_i^(2 / (m - 1)), quad 1 - u_(i c_i)^* (m) <= (k - 1) rho_i^(2 / (m - 1)). $
     Consequently, as $m -> 1^+$, $U^*(m) -> H(bold(mu))$ (the hard Voronoi matrix), the fuzzy centroids $mu_j^"FCM" (m) = (sum_i (u_(i j)^*)^m x_i) / (sum_i (u_(i j)^*)^m)$ converge to Lloyd's centroids $mu_j^"Lloyd"$, and $min_(U in Delta_k^n) J_m (U, bold(mu)) -> J^*(bold(mu))$.
  2. *FCM Fuzzifier Collapse ($m -> infinity$):* For any configuration with $d_(i j) > 0$ for all $i, j$, as $m -> infinity$, $u_(i j)^* (m) -> 1 / k$ for all $(i, j)$, while a single centroid update approaches the geometric-ratio average $mu_j^"FCM" (m) -> (sum_i (tilde(d)_i / d_(i j))^2 x_i) / (sum_i (tilde(d)_i / d_(i j))^2)$ where $tilde(d)_i := (product_(ell=1)^k d_(i ell))^(1 / k)$. Moreover, letting $S_"dir" := n^(-1) sum_(i=1)^n (x_i - overline(x))(x_i - overline(x))^T / norm(x_i - overline(x))_2^2$ (for $x_i != overline(x)$), whenever $lambda_"max" (S_"dir") < 1 / 2$ and $m > m_"crit" := (1 - 2 lambda_"max" (S_"dir"))^(-1)$, the degenerate global-mean configuration $mu_1 = dots = mu_k = overline(x)$ with $u_(i j) = 1 / k$ becomes a *strict local minimizer* of the reduced FCM objective $R_m (bold(mu))$, trapping FCM in complete collapse.
  3. *Exact Isomorphism of Entropy-Regularized K-Means and Softmax EM:* Under the exact temperature-variance dictionary $T = 2 sigma^2$:
     (a) Minimizing $cal(F)_T (U, bold(mu), bold(pi))$ over $U in Delta_k^n$ uniquely attains the exact isotropic GMM E-step responsibility matrix $u_(i j)^* = gamma_(i j) (sqrt(T / 2))$.
     (b) Minimizing $cal(F)_T (U, bold(mu), bold(pi))$ over $(bold(mu), bold(pi))$ reproduces the exact GMM M-step updates $mu_j^* = (sum_i u_(i j) x_i) / (sum_i u_(i j))$ and $pi_j^* = n^(-1) sum_i u_(i j)$.
     (c) Profiling out $U$ yields the exact marginal free-energy identity
     $ min_(U in Delta_k^n) cal(F)_T (U, bold(mu), bold(pi)) &= -T sum_(i=1)^n log(sum_(j=1)^k pi_j exp(- (norm(x_i - mu_j)_2^2) / T)) \
       &= -2 sigma^2 ell(bold(mu), bold(pi); sigma) - n d sigma^2 log(2 pi sigma^2). $
     In particular, alternating coordinate descent on $cal(F)_T$ is *identically* isotropic GMM EM (Alg. 81), and converges to Lloyd's $k$-means (Alg. 1) as $T -> 0^+$.
]

#proof[
  *Part 1 (FCM KKT Derivation and $m -> 1^+$ Rate):* Fix $i in {1, dots, n}$. Because $m > 1$, the map $u_i |-> sum_(j=1)^k u_(i j)^m d_(i j)^2$ is strictly convex on $Delta_k$ whenever $d_(i j) > 0$ for all $j$. Forming the Lagrangian $L_i (u_i, lambda) = sum_(j=1)^k u_(i j)^m d_(i j)^2 - lambda (sum_(j=1)^k u_(i j) - 1)$ and setting $partial L_i / partial u_(i j) = m u_(i j)^(m - 1) d_(i j)^2 - lambda = 0$ yields $u_(i j) = (lambda / m)^(1 / (m - 1)) d_(i j)^(-2 / (m - 1))$. Enforcing $sum_(ell=1)^k u_(i ell) = 1$ eliminates $(lambda / m)^(1 / (m - 1))$ and gives
  $ u_(i j)^* (m) = (d_(i j)^(-2 / (m - 1))) / (sum_(ell=1)^k d_(i ell)^(-2 / (m - 1))) = (sum_(ell=1)^k (d_(i j)^2 / d_(i ell)^2)^(1 / (m - 1)))^(-1). $
  (If the singularity set $I_i := {j : d_(i j) = 0}$ is non-empty, then setting $u_(i j) = 0$ for $j in.not I_i$ and arbitrary $sum_(j in I_i) u_(i j) = 1$ achieves the global minimum $0$; when $I_i = {c_i}$, $u_(i c_i)^* (m) = 1$ and $u_(i ell)^* (m) = 0$ for $ell != c_i$, so the bound holds trivially with $rho_i = 0$.)

  Now assume $d_(i j) > 0$ for all $j$ and $c_i$ is unique, so $rho_i = max_(ell != c_i) (d_(i c_i) / d_(i ell)) < 1$. For any $ell != c_i$, retain only the $c_i$-th term in the denominator:
  $ u_(i ell)^* (m) = (sum_(r=1)^k (d_(i ell) / d_(i r))^(2 / (m - 1)))^(-1) <= (d_(i c_i) / d_(i ell))^(2 / (m - 1)) <= rho_i^(2 / (m - 1)). $
  Summing over all $k - 1$ non-winning clusters gives $1 - u_(i c_i)^* (m) <= (k - 1) rho_i^(2 / (m - 1))$. Since $rho_i < 1$ and $2 / (m - 1) -> +infinity$ as $m -> 1^+$, we obtain $u_(i j)^* (m) -> bb(1)_(j = c_i)$ and $(u_(i j)^* (m))^m -> bb(1)_(j = c_i)$ exponentially in $(m - 1)^(-1)$. Substituting into the FCM centroid update $mu_j^"FCM" (m) = (sum_i (u_(i j)^*)^m x_i) / (sum_i (u_(i j)^*)^m)$ (for every non-empty Voronoi cell $|C_j| >= 1$) yields $mu_j^"FCM" (m) -> |C_j|^(-1) sum_(i in C_j) x_i = mu_j^"Lloyd"$. Furthermore, evaluating $J_m (U^*(m), bold(mu))$ directly gives the harmonic-power representation
  $ J_m (U^*(m), bold(mu)) = sum_(i=1)^n (sum_(j=1)^k (d_(i j)^2)^(-1 / (m - 1)))^(-(m - 1)). $
  Letting $p := 1 / (m - 1) -> +infinity$, the generalized $ell_(-p)$ mean satisfies
  $ k^(-1 / p) d_(i c_i)^2 <= (sum_(j=1)^k d_(i j)^(-2 p))^(-1 / p) <= d_(i c_i)^2, $
  whose left-hand side converges monotonically to $d_(i c_i)^2$ as $p -> +infinity$ (since $k^(-1 / p) = k^(-(m - 1)) -> 1$), proving $min_U J_m (U, bold(mu)) -> sum_(i=1)^n d_(i c_i)^2 = J^*(bold(mu))$ even when ties occur!

  *Part 2 ($m -> infinity$ Limit & Hessian Collapse Criterion):* When $d_(i j) > 0$ for all $i, j$, as $m -> infinity$ we have $1 / (m - 1) -> 0^+$, so $(d_(i j)^2 / d_(i ell)^2)^(1 / (m - 1)) -> 1$ for every pair $(j, ell)$, whence $u_(i j)^* (m) -> 1 / k$. For the single-step centroid update, divide numerator and denominator by $k^(-m)$:
  $ mu_j^"FCM" (m) = (sum_(i=1)^n (k u_(i j)^* (m))^m x_i) / (sum_(i=1)^n (k u_(i j)^* (m))^m). $
  Expanding $log(k u_(i j)^* (m)) = - log(k^(-1) sum_(ell=1)^k exp((2 / (m - 1)) log(d_(i j) / d_(i ell)))) = -(2 / (m - 1)) (log d_(i j) - k^(-1) sum_(ell=1)^k log d_(i ell)) + O(m^(-2))$ and multiplying by $m$ gives $(k u_(i j)^* (m))^m -> (tilde(d)_i / d_(i j))^2$ as $m -> infinity$.

  To establish global-mean collapse of the *joint* optimization (@bezdek1981fcm, @hathaway1986fcm), examine the reduced objective $R_m (bold(mu)) := min_(U in Delta_k^n) J_m (U, bold(mu)) = sum_(i=1)^n S_i^(-1 / p)$, where $p = 1 / (m - 1)$ and $S_i := sum_(j=1)^k (d_(i j)^2)^(-p)$. Writing $v_i := x_i - overline(x)$ and $r_i := norm(v_i)_2 > 0$, differentiate $R_m$ twice at $bold(mu)_0 = (overline(x), dots, overline(x))$ where $S_i = k r_i^(-2 p)$:
  $ nabla_(mu_j) R_m (bold(mu)_0) = 2 k^(-m) sum_(i=1)^n (overline(x) - x_i) = bold(0), $
  confirming stationarity, while the $(j, ell)$ $d times d$ Hessian block evaluates via the product rule to
  $ H_(j ell) (bold(mu)_0) = 2 n k^(-m) [delta_(j ell) I_d - (2 m) / (m - 1) (delta_(j ell) - 1 / k) S_"dir"]. $
  Decompose any perturbation $delta bold(mu) in (bb(R)^d)^k$ orthogonally into a uniform translation $delta mu_j = w$ (where $delta_(j ell) - 1 / k$ vanishes, yielding strictly positive curvature $2 n k^(1 - m) norm(w)_2^2 > 0$) and a zero-sum split perturbation $sum_(j=1)^k delta mu_j = bold(0)$ (where $delta_(j ell) - 1 / k$ acts as $delta_(j ell)$). On the zero-sum subspace, the quadratic form is
  $ sum_(j, ell = 1)^k delta mu_j^T H_(j ell) delta mu_ell = 2 n k^(-m) sum_(j=1)^k delta mu_j^T (I_d - (2 m) / (m - 1) S_"dir") delta mu_j. $
  This operator is strictly positive definite if and only if $(2 m / (m - 1)) lambda_"max" (S_"dir") < 1$, which rearranges identically to $lambda_"max" (S_"dir") < 1 / 2$ and $m > (1 - 2 lambda_"max" (S_"dir"))^(-1) = m_"crit"$. Note that since $"Tr"(S_"dir") = 1$, isotropic high-dimensional data satisfy $lambda_"max" (S_"dir") approx 1 / d$, driving $m_"crit" approx d / (d - 2) -> 1^+$ and explaining why standard $m = 2$ FCM suffers catastrophic collapse to $overline(x)$ in high dimensions!

  *Part 3 (Entropy-Regularized K-Means <=> Softmax EM):* Fix $(bold(mu), bold(pi))$ and minimize $cal(F)_T (U, bold(mu), bold(pi))$ over $U in Delta_k^n$. Since $u_i |-> "KL"(u_i || bold(pi))$ is strictly convex on $Delta_k$ with gradient blowing up at the boundary $u_(i j) -> 0^+$, the unique minimizer lies in $"int"(Delta_k)$. Form the Lagrangian with multiplier $nu_i$ for $sum_j u_(i j) = 1$:
  $ partial / (partial u_(i j)) [sum_(ell=1)^k u_(i ell) d_(i ell)^2 + T sum_(ell=1)^k u_(i ell) log(u_(i ell) / pi_ell) - nu_i (sum_(ell=1)^k u_(i ell) - 1)] = 0. $
  Differentiating gives $d_(i j)^2 + T (log(u_(i j) / pi_j) + 1) - nu_i = 0$, hence
  $ u_(i j)^* = pi_j exp((nu_i - T) / T) exp(- d_(i j)^2 / T). $
  Enforcing $sum_(ell=1)^k u_(i ell)^* = 1$ determines the normalizing constant $Z_i := sum_(ell=1)^k pi_ell exp(- d_(i ell)^2 / T)$, yielding
  $ u_(i j)^* = (pi_j exp(- d_(i j)^2 / T)) / (sum_(ell=1)^k pi_ell exp(- d_(i ell)^2 / T)). $
  Substituting $T = 2 sigma^2$ gives $u_(i j)^* equiv gamma_(i j) (sigma)$, proving 3(a). Next, holding $U$ fixed: minimizing $sum_(i, j) u_(i j) norm(x_i - mu_j)_2^2$ over $mu_j$ gives $mu_j^* = (sum_i u_(i j) x_i) / (sum_i u_(i j))$ by the weighted Parallel-Axis identity; minimizing $T sum_i "KL"(u_i || bold(pi)) = "const" - T sum_(i, j) u_(i j) log pi_j$ over $bold(pi) in Delta_k$ via Gibbs' inequality yields $pi_j^* = n^(-1) sum_i u_(i j)$, proving 3(b). Finally, substituting $u_(i j)^* = pi_j exp(- d_(i j)^2 / T) / Z_i$ back into $cal(F)_T$:
  $ cal(F)_T (U^*, bold(mu), bold(pi)) &= sum_(i=1)^n sum_(j=1)^k u_(i j)^* [d_(i j)^2 + T (- d_(i j)^2 / T - log Z_i)] \
    &= - T sum_(i=1)^n (sum_(j=1)^k u_(i j)^*) log Z_i = - T sum_(i=1)^n log Z_i. $
  Since $p(x_i | bold(mu), bold(pi); sigma) = (2 pi sigma^2)^(-d / 2) Z_i$ when $T = 2 sigma^2$, we have
  $ -2 sigma^2 ell(bold(mu), bold(pi); sigma) = n d sigma^2 log(2 pi sigma^2) - 2 sigma^2 sum_(i=1)^n log Z_i, $
  so $min_U cal(F)_(2 sigma^2) (U, bold(mu), bold(pi)) = -2 sigma^2 ell(bold(mu), bold(pi); sigma) - n d sigma^2 log(2 pi sigma^2)$, completing the proof.
]

#remark(title: "Algebraic vs. Exponential Tail Decay and Scale Invariance")[
  Theorem 4 reveals a deep structural contrast between Bezdek's FCM (Alg. 151) and Softmax / Entropy-Regularized EM (Alg. 81):
  - *Scale Invariance:* Rescaling all data $x_i |-> c x_i$ ($c > 0$) leaves FCM memberships $u_(i j)^* (m)$ *strictly invariant* (since $d_(i j) / d_(i ell)$ is homogeneous of degree zero), whereas Entropy-Regularized $k$-Means requires quadratic temperature rescaling $T |-> c^2 T$ to preserve assignments.
  - *Outlier Tail Behavior:* For a distant outlier with $d_(i 1) approx d_(i 2) >> 1$, FCM decays algebraically ($u_(i 2) / u_(i 1) = (d_(i 1) / d_(i 2))^(2 / (m - 1)) -> 1$), pulling all centroids symmetrically, whereas Softmax EM decays exponentially in squared-distance difference ($exp(-(d_(i 2)^2 - d_(i 1)^2) / T)$), snapping far points onto a single winning cluster unless heavy-tailed Student-$t$ mixtures (Alg. 88, @peel2000robust) are employed. Notice also that replacing Euclidean distance by RKHS distance immediately extends Part 1 to *Kernel FCM* (Alg. 155, @zhang2004kfcm), while replacing the inner power mean by $M_s$ recovers *Power $k$-Means* (Alg. 27, @xu2019power).
]

#subsec("C", "Spherical K-Means as the Zero-Temperature Limit of von Mises–Fisher Mixtures")

For high-dimensional directional data (normalized TF–IDF corpora, unit-norm deep embeddings), observations reside on the unit hypersphere $bb(S)^(d - 1) := {x in bb(R)^d : norm(x)_2 = 1}$ ($d >= 2$). *Spherical $k$-Means* (Alg. 12, @dhillon2001concept) maximizes total cosine similarity over unit prototypes $bold(mu) = (mu_1, dots, mu_k) in (bb(S)^(d - 1))^k$:
$ J_"cos" (cal(C), bold(mu)) := sum_(j=1)^k sum_(i in C_j) chevron.l x_i, mu_j chevron.r = n - 1 / 2 sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2, $
iterating cosine assignment $c_i = arg max_j chevron.l x_i, mu_j chevron.r$ and normalized resultant update $mu_j^"Sph" = R_j / norm(R_j)_2$ where $R_j := sum_(i in C_j) x_i$. Its probabilistic counterpart in Part IV is the *von Mises–Fisher (vMF) Mixture Model* (Alg. 89, @banerjee2005vmf) with shared concentration $kappa > 0$:
$ f_d (x | mu_j, kappa) := c_d (kappa) exp(kappa chevron.l x, mu_j chevron.r), quad c_d (kappa) := (kappa^(d / 2 - 1)) / ((2 pi)^(d / 2) I_(d / 2 - 1) (kappa)), $
where $I_nu (kappa)$ is the modified Bessel function of the first kind of order $nu$. We first establish a sharp elementary Lipschitz bound for radial projection onto $bb(S)^(d - 1)$.

#lemma(title: "Spherical Normalization Lipschitz Bound", origin: "Original Proof")[
  For any two nonzero vectors $u, v in bb(R)^d backslash {bold(0)}$,
  $ norm(u / norm(u)_2 - v / norm(v)_2)_2 <= (2 norm(u - v)_2) / (norm(v)_2). $
]

#proof[
  Add and subtract $u / norm(v)_2$ inside the norm and apply the triangle and reverse triangle inequalities:
  $ norm(u / norm(u)_2 - v / norm(v)_2)_2 &= norm((u - v) / norm(v)_2 + u (1 / norm(u)_2 - 1 / norm(v)_2))_2 \
    &<= (norm(u - v)_2) / (norm(v)_2) + norm(u)_2 (|norm(v)_2 - norm(u)_2|) / (norm(u)_2 norm(v)_2) \
    &<= (norm(u - v)_2) / (norm(v)_2) + (norm(u - v)_2) / (norm(v)_2) = (2 norm(u - v)_2) / (norm(v)_2). $
]

#theorem(title: "Spherical K-Means <=> Zero-Temperature Limit (Kappa -> Infinity) of vMF Mixture", origin: "Original Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(S)^(d - 1)$ ($d >= 2$), $bold(mu) in (bb(S)^(d - 1))^k$, $bold(pi) in "int"(Delta_k)$ with $pi_"min" = min_j pi_j > 0$, and let $ell_"vMF" (bold(mu), bold(pi); kappa) := sum_(i=1)^n log(sum_(j=1)^k pi_j f_d (x_i | mu_j, kappa))$ denote the homoscedastic vMF mixture log-likelihood.
  1. *Exact vMF EM Equations:* For any $kappa > 0$, the E-step posterior and M-step mean-direction update are
     $ &gamma_(i j) (kappa) = (pi_j exp(kappa chevron.l x_i, mu_j chevron.r)) / (sum_(ell=1)^k pi_ell exp(kappa chevron.l x_i, mu_ell chevron.r)), \
       &r_j (kappa) := sum_(i=1)^n gamma_(i j) (kappa) x_i, quad mu_j^"vMF" (kappa) = (r_j (kappa)) / (norm(r_j (kappa))_2), $
     whenever $r_j (kappa) != bold(0)$.
  2. *Exponential Convergence to Spherical K-Means Updates:* Suppose every $x_i$ has unique cosine maximizer $c_i$ with strict angular margin $eta_i (bold(mu)) := chevron.l x_i, mu_(c_i) chevron.r - max_(ell != c_i) chevron.l x_i, mu_ell chevron.r > 0$, let $eta(bold(mu)) := min_i eta_i (bold(mu)) > 0$, and assume each hard cluster resultant satisfies $norm(R_j)_2 > 0$. Then for every $ell != c_i$ and every $kappa > 0$,
     $ gamma_(i ell) (kappa) <= (pi_ell / pi_(c_i)) exp(- kappa eta_i (bold(mu))), $
     and for all sufficiently large $kappa$ ensuring $n pi_"min"^(-1) exp(- kappa eta(bold(mu))) < norm(R_j)_2$:
     $ norm(mu_j^"vMF" (kappa) - mu_j^"Sph")_2 <= (2 n) / (pi_"min" norm(R_j)_2) exp(- kappa eta(bold(mu))). $
  3. *Asymptotic Bessel Free-Energy Limit:* As $kappa -> +infinity$, the concentration-rescaled vMF log-likelihood converges uniformly on ${bold(pi) : min_j pi_j >= epsilon_0}$ to both the negative spherical cosine distortion and (equivalently) the negative half-WCSS on $bb(S)^(d - 1)$:
     $ lim_(kappa -> +infinity) 1 / kappa ell_"vMF" (bold(mu), bold(pi); kappa) &= sum_(i=1)^n (max_(1 <= j <= k) chevron.l x_i, mu_j chevron.r - 1) \
       &= J_"cos"^* (bold(mu)) - n = - 1 / 2 J_"WCSS"^* (bold(mu)), $
     with explicit asymptotic error $O((log kappa) / kappa)$.
]

#proof[
  *Part 1 (Exact EM Updates):* The E-step computes $P(Z_i = j | x_i) prop pi_j c_d (kappa) exp(kappa chevron.l x_i, mu_j chevron.r)$, where the shared normalization factor $c_d (kappa)$ cancels identically between numerator and denominator, yielding $gamma_(i j) (kappa)$. In the M-step, maximizing the expected complete-data log-likelihood over $mu_j in bb(S)^(d - 1)$ reduces to
  $ max_(norm(mu_j)_2 = 1) kappa sum_(i=1)^n gamma_(i j) (kappa) chevron.l x_i, mu_j chevron.r = max_(norm(mu_j)_2 = 1) kappa chevron.l r_j (kappa), mu_j chevron.r. $
  By the Cauchy–Schwarz inequality, $chevron.l r_j (kappa), mu_j chevron.r <= norm(r_j (kappa))_2 norm(mu_j)_2 = norm(r_j (kappa))_2$, with unique equality at $mu_j^"vMF" (kappa) = r_j (kappa) / norm(r_j (kappa))_2$ whenever $r_j (kappa) != bold(0)$ @banerjee2005vmf.

  *Part 2 (Exponential Concentration & Spherical Centroid Bound):* Let $s_(i *) := max_j chevron.l x_i, mu_j chevron.r = chevron.l x_i, mu_(c_i) chevron.r$. Multiplying numerator and denominator of $gamma_(i j) (kappa)$ by $exp(- kappa s_(i *))$ yields
  $ gamma_(i ell) (kappa) &= (pi_ell exp(- kappa (s_(i *) - chevron.l x_i, mu_ell chevron.r))) / (pi_(c_i) + sum_(r != c_i) pi_r exp(- kappa (s_(i *) - chevron.l x_i, mu_r chevron.r))) \
    &<= (pi_ell / pi_(c_i)) exp(- kappa eta_i (bold(mu))) $
  for every $ell != c_i$, and $1 - gamma_(i c_i) (kappa) <= pi_"min"^(-1) exp(- kappa eta_i (bold(mu)))$. Consequently, comparing the soft resultant $r_j (kappa)$ to the hard resultant $R_j = sum_(i in C_j) x_i$ using $norm(x_i)_2 = 1$:
  $ norm(r_j (kappa) - R_j)_2 &= norm(sum_(i=1)^n (gamma_(i j) (kappa) - bb(1)_(i in C_j)) x_i)_2 \
    &<= sum_(i=1)^n |gamma_(i j) (kappa) - bb(1)_(i in C_j)| <= (n / pi_"min") exp(- kappa eta(bold(mu))). $
  When this bound is strictly smaller than $norm(R_j)_2$, we have $r_j (kappa) != bold(0)$, and applying Lemma 2 with $u = r_j (kappa)$ and $v = R_j$ gives
  $ norm(mu_j^"vMF" (kappa) - mu_j^"Sph")_2 &<= (2 norm(r_j (kappa) - R_j)_2) / (norm(R_j)_2) \
    &<= (2 n) / (pi_"min" norm(R_j)_2) exp(- kappa eta(bold(mu))). $

  *Part 3 (Hankel Bessel Asymptotics & Laplace Limit):* Recall the classical Hankel large-argument asymptotic expansion of the modified Bessel function $I_nu (kappa)$ for fixed order $nu = d / 2 - 1$:
  $ I_nu (kappa) = (e^kappa) / sqrt(2 pi kappa) (1 - (4 nu^2 - 1) / (8 kappa) + O(kappa^(-2))) quad (kappa -> +infinity). $
  Taking logarithms of $c_d (kappa) = kappa^nu / ((2 pi)^(nu + 1) I_nu (kappa))$ gives
  $ log c_d (kappa) &= nu log kappa - (nu + 1) log(2 pi) \
    &quad - [kappa - 1 / 2 log(2 pi kappa) + O(kappa^(-1))] \
    &= - kappa + ((d - 1) / 2) log(kappa / (2 pi)) + O(kappa^(-1)). $
  Hence $(n / kappa) log c_d (kappa) = - n + ((n (d - 1)) / (2 kappa)) log(kappa / (2 pi)) + O(kappa^(-2)) -> - n$ as $kappa -> +infinity$. Meanwhile, factoring $exp(kappa s_(i *))$ out of the mixture sum $M_i (kappa) := sum_(j=1)^k pi_j exp(kappa chevron.l x_i, mu_j chevron.r)$ yields the exact two-sided bound
  $ pi_(c_i) exp(kappa s_(i *)) <= M_i (kappa) <= exp(kappa s_(i *)). $
  Taking logarithms, summing over $i = 1, dots, n$, adding $n log c_d (kappa)$, and dividing by $kappa$:
  $ &sum_(i=1)^n s_(i *) + (1 / kappa) sum_(i=1)^n log pi_(c_i) + (n / kappa) log c_d (kappa) \
    &quad <= 1 / kappa ell_"vMF" (bold(mu), bold(pi); kappa) <= sum_(i=1)^n s_(i *) + (n / kappa) log c_d (kappa). $
  Sending $kappa -> +infinity$ squeezes $kappa^(-1) ell_"vMF" (bold(mu), bold(pi); kappa)$ to $sum_(i=1)^n (s_(i *) - 1) = J_"cos"^* (bold(mu)) - n$ at rate $O((log kappa) / kappa)$. Finally, because $norm(x_i)_2 = norm(mu_j)_2 = 1$, we have $1 - chevron.l x_i, mu_j chevron.r = 1 / 2 norm(x_i - mu_j)_2^2$, whence $J_"cos"^* (bold(mu)) - n = - 1 / 2 J_"WCSS"^* (bold(mu))$.
]

#subsec("D", "Universal Exponential-Family Mixture <=> Bregman Partitioning Duality")

Theorems 3 and 5 are not isolated coincidences: they are particular instances of a grand duality theorem due to Banerjee et al. @banerjee2005clustering linking *all* regular exponential-family mixture models in Part IV—Gaussian (Alg. 81), von Mises–Fisher (Alg. 89), Poisson (Alg. 95), Bernoulli (Alg. 96), and Dirichlet–Multinomial / Multinomial mixtures (Alg. 97)—to hard Bregman divergence partitioning in Part I (Algs. 1, 7, 12, 16).

Let $p_psi (x | theta) = exp(chevron.l x, theta chevron.r - psi(theta)) p_0 (x)$ be a minimal regular exponential family on $cal(X) subset.eq bb(R)^d$ with strictly convex, smooth log-partition function $psi: Theta -> bb(R)$ on open natural domain $Theta subset.eq bb(R)^d$. Let $phi(mu) := psi^* (mu) = sup_(theta in Theta) {chevron.l mu, theta chevron.r - psi(theta)}$ be its Legendre–Fenchel conjugate on the expectation-parameter domain $cal(M) = nabla psi(Theta)$, so the bijection $mu = nabla psi(theta) arrow.l.r.double theta = nabla phi(mu)$ holds globally and generates the Bregman divergence
$ D_phi (x || mu) := phi(x) - phi(mu) - chevron.l x - mu, nabla phi(mu) chevron.r >= 0. $
To formalize a temperature scale across arbitrary exponential families, introduce the *exponential dispersion family* with dispersion parameter $tau > 0$:
$ p_(psi, tau) (x | theta) := exp(1 / tau [chevron.l x, theta chevron.r - psi(theta)]) p_(0, tau) (x), $
where $log p_(0, tau) (x) = -(1 / tau) phi(x) + h_tau (x)$ (or equivalently, $tau = 1 / m$ represents averaging $m$ i.i.d. observations from $p_psi$).

#theorem(title: "Universal Exponential-Family Mixture <=> Bregman K-Means Duality", origin: "Classical — Cited")[
  Let $p_(psi, tau) (x | theta)$ be a regular exponential dispersion family with conjugate pair $(psi, phi = psi^*)$ @banerjee2005clustering.
  1. *Exact Density–Bregman Bijection:* For every $x in "dom"(phi)$ and $mu = nabla psi(theta) in cal(M)$,
     $ p_(psi, tau) (x | theta) = exp(- (D_phi (x || mu)) / tau) b_(phi, tau) (x), $
     where $b_(phi, tau) (x) := exp(phi(x) / tau) p_(0, tau) (x)$ is strictly independent of the cluster parameter $mu$.
  2. *Exact Soft-EM / Bregman Identity:* For a $k$-component mixture $sum_(j=1)^k pi_j p_(psi, tau) (x | theta_j)$, the E-step responsibilities and M-step expectation parameters are given identically in Bregman form by
     $ &gamma_(i j) (tau) = (pi_j exp(- D_phi (x_i || mu_j) / tau)) / (sum_(ell=1)^k pi_ell exp(- D_phi (x_i || mu_ell) / tau)), \
       &mu_j^"EM" (tau) = (sum_(i=1)^n gamma_(i j) (tau) x_i) / (sum_(i=1)^n gamma_(i j) (tau)), $
     because the arithmetic centroid uniquely minimizes expected Bregman distortion $mu |-> sum_i gamma_(i j) D_phi (x_i || mu)$.
  3. *Zero-Dispersion Reduction ($tau -> 0^+$):* Away from Bregman Voronoi bisectors ($D_phi (x_i || mu_(c_i)) < D_phi (x_i || mu_ell)$ for all $ell != c_i$), as $tau -> 0^+$ the mixture EM operator converges exponentially to hard *Bregman $k$-Means*, and the rescaled log-likelihood satisfies
     $ &lim_(tau -> 0^+) [- tau ell(bold(mu), bold(pi); tau) + tau sum_(i=1)^n log b_(phi, tau) (x_i)] \
       &quad = sum_(i=1)^n min_(1 <= j <= k) D_phi (x_i || mu_j). $
]

#proof[
  *Part 1:* By Legendre–Fenchel duality on the open domain $Theta$, $theta = nabla phi(mu)$ and $psi(theta) + phi(mu) = chevron.l mu, theta chevron.r = chevron.l mu, nabla phi(mu) chevron.r$. Therefore
  $ chevron.l x, theta chevron.r - psi(theta) &= chevron.l x, nabla phi(mu) chevron.r - (chevron.l mu, nabla phi(mu) chevron.r - phi(mu)) \
    &= phi(x) - [phi(x) - phi(mu) - chevron.l x - mu, nabla phi(mu) chevron.r] \
    &= phi(x) - D_phi (x || mu). $
  Dividing by $tau$, exponentiating, and multiplying by $p_(0, tau) (x)$ yields $p_(psi, tau) (x | theta) = exp(- D_phi (x || mu) / tau) exp(phi(x) / tau) p_(0, tau) (x) = exp(- D_phi (x || mu) / tau) b_(phi, tau) (x)$.

  *Part 2:* In the E-step, the data-only factor $b_(phi, tau) (x_i)$ cancels between numerator and denominator, giving $gamma_(i j) (tau)$ directly in terms of $D_phi (x_i || mu_j)$. In the M-step, maximizing $sum_(i=1)^n gamma_(i j) log p_(psi, tau) (x_i | theta_j)$ over $theta_j$ is equivalent to minimizing $F_j (mu_j) := sum_(i=1)^n gamma_(i j) D_phi (x_i || mu_j)$ over $mu_j in cal(M)$. Let $overline(mu)_j := (sum_i gamma_(i j) x_i) / (sum_i gamma_(i j)) in cal(M)$ (which lies in $cal(M)$ by convexity of $cal(M)$). Apply the *Bregman three-point identity*:
  $ &D_phi (x_i || mu_j) - D_phi (x_i || overline(mu)_j) - D_phi (overline(mu)_j || mu_j) \
    &quad = chevron.l x_i - overline(mu)_j, space nabla phi(overline(mu)_j) - nabla phi(mu_j) chevron.r. $
  Multiplying by $gamma_(i j)$ and summing over $i = 1, dots, n$ makes the right-hand side vanish identically since $sum_i gamma_(i j) (x_i - overline(mu)_j) = bold(0)$, yielding the generalized Parallel-Axis identity:
  $ sum_(i=1)^n gamma_(i j) D_phi (x_i || mu_j) &= sum_(i=1)^n gamma_(i j) D_phi (x_i || overline(mu)_j) \
    &quad + (sum_(i=1)^n gamma_(i j)) D_phi (overline(mu)_j || mu_j). $
  Since $phi$ is strictly convex, $D_phi (overline(mu)_j || mu_j) >= 0$ with equality if and only if $mu_j = overline(mu)_j$, proving uniqueness of the weighted arithmetic centroid across *every* regular exponential family!

  *Part 3:* Define $delta_i^* := min_j D_phi (x_i || mu_j)$ achieved uniquely at $c_i$ with gap $Delta_i^phi := min_(ell != c_i) D_phi (x_i || mu_ell) - delta_i^* > 0$. Repeating the exact Laplace sandwich of Theorem 3 with $d_(i j)^2 / (2 sigma^2)$ replaced by $D_phi (x_i || mu_j) / tau$:
  $ &pi_(c_i) b_(phi, tau) (x_i) exp(- delta_i^* / tau) <= sum_(j=1)^k pi_j p_(psi, tau) (x_i | theta_j) \
    &quad <= b_(phi, tau) (x_i) exp(- delta_i^* / tau), $
  whence $gamma_(i ell) (tau) <= (pi_ell / pi_(c_i)) exp(- Delta_i^phi / tau) -> 0$ and
  $ sum_(i=1)^n delta_i^* &<= - tau ell(bold(mu), bold(pi); tau) + tau sum_(i=1)^n log b_(phi, tau) (x_i) \
    &<= sum_(i=1)^n delta_i^* - tau sum_(i=1)^n log pi_(c_i), $
  which converges uniformly to $sum_(i=1)^n min_j D_phi (x_i || mu_j)$ as $tau -> 0^+$.
]

Table 2 summarizes the concrete instances of Theorem 6 connecting Part IV mixture algorithms to Part I/IV Bregman partitioners, including the exact reduction of binary Bernoulli Mixture Models (Alg. 96) with coordinate-wise thresholding to *Binary $k$-Modes* (Alg. 7, @huang1998extensions).

#figure(
  table(
    columns: (0.95fr, 0.95fr, 1.1fr),
    align: (left, left, left),
    [Mixture Model (Part IV)], [Conjugate $phi(x)$ & $D_phi (x || mu)$], [Zero-Temp Limit ($tau -> 0^+$)],
    [Isotropic GMM (Alg. 81)], [$1/2 norm(x)_2^2$; $1/2 norm(x - mu)_2^2$], [Lloyd's $k$-Means (Alg. 1)],
    [vMF Mixture (Alg. 89)], [On $bb(S)^(d-1)$: $1 - chevron.l x, mu chevron.r$], [Spherical $k$-Means (Alg. 12)],
    [Poisson Mixture (Alg. 95)], [$sum x_r log x_r - x_r$; I-div.], [KL-Bregman / Hard KL-NMF (Alg. 92)],
    [Bernoulli Mixture (Alg. 96)], [Binary neg. entropy; Logistic KL], [Logistic $k$-Means / Binary $k$-Modes (Alg. 7)],
    [Multinomial / DMM (Alg. 97)], [Simplex neg. entropy; $"KL"(x || mu)$], [Information-Theoretic $k$-Means / pLSA (Alg. 84)],
  ),
  caption: [Exact specializations of Theorem 6 bridging Part IV exponential-family mixtures and Part I/IV Bregman partitioners @banerjee2005clustering.],
)


#sec("IV", "Graph-Theoretic, Ultrametric, Density, and Agglomerative Equivalences")

We now bridge Parts II (Density and Mode-Seeking), III (Hierarchical and Tree-Based), and I (Partitioning). Four exact theorems prove that Single Linkage, Kruskal's Minimum Spanning Tree, DBSCAN, and HDBSCAN are manifestations of a single ultrametric level-set functor; that Mean Shift and DENCLUE 2.0 are exact Gaussian-KDE EM algorithms; and that Ward's agglomerative hierarchy is exact greedy step-by-step minimization of the $k$-means WCSS functional.

#subsec("A", "Single Linkage, Kruskal MST Thresholding, and Core-Point DBSCAN")

Let $(cal(X), delta)$ be a finite metric space with $|cal(X)| = n$, and let $K_n = (cal(X), E, delta)$ be the complete undirected graph with edge weights $delta(x_i, x_j)$. For any scale $epsilon >= 0$, define the *$epsilon$-neighborhood graph* $G_epsilon := (cal(X), E_epsilon)$ where $E_epsilon := {{x_i, x_j} : i != j, space delta(x_i, x_j) <= epsilon}$, and for any spanning tree $T$ of $K_n$, let $T_epsilon := (cal(X), {e in E(T) : delta(e) <= epsilon})$ denote its $epsilon$-thresholded spanning forest. Recall that an *ultrametric* $u: cal(X) times cal(X) -> bb(R)_(>= 0)$ satisfies the strong triangle inequality $u(x, z) <= max{u(x, y), u(y, z)}$ for all $x, y, z in cal(X)$.

In DBSCAN (Alg. 31, @ester1996dbscan, @schubert2017dbscanrevisited) with parameters $(epsilon, "MinPts")$, the closed $epsilon$-ball is $B_epsilon (x) := {y in cal(X) : delta(x, y) <= epsilon}$, and the *core-point set* is
$ cal(X)_"core" (epsilon, "MinPts") := {x in cal(X) : |B_epsilon (x)| >= "MinPts"}. $
A point $y$ is *directly density-reachable* from $x$ if $x in cal(X)_"core"$ and $y in B_epsilon (x)$; *density-reachable* ($x ->^* y$) if linked by a directed chain of direct density-reachability steps; and two points $x, y$ are *density-connected* if some $o in cal(X)_"core"$ satisfies $o ->^* x$ and $o ->^* y$.

#theorem(title: "Single Linkage <=> Kruskal MST Thresholding <=> Core-Point DBSCAN", origin: "Original Proof")[
  Let $(cal(X), delta)$ be a finite metric space and let $T^*$ be *any* Minimum Spanning Tree (MST) of $(K_n, delta)$ @gower1969minimum, @sibson1973slink, @carlsson2010characterization.
  1. *Single Linkage <=> MST Thresholding:* For every threshold $epsilon >= 0$, the following three partitions of $cal(X)$ are *identically equal*:
     (a) the connected components of the $epsilon$-neighborhood graph $G_epsilon$;
     (b) the connected components of the $epsilon$-thresholded MST forest $T^*_epsilon$;
     (c) the partition obtained by cutting the Single Linkage (Alg. 56) dendrogram at height $epsilon$.
  2. *Subdominant Ultrametric Identity:* The Single Linkage cophenetic distance $u_"SL" (x, y)$ equals both the minimax path bottleneck across $K_n$ and the unique maximum edge weight along the unique path $P_(T^*) (x, y)$ in $T^*$:
     $ u_"SL" (x, y) = min_(P: x -> y) max_(e in P) delta(e) = max_(e in P_(T^*) (x, y)) delta(e), $
     and $u_"SL"$ is the *unique maximal (subdominant) ultrametric* dominated by $delta$ ($u <= delta$ pointwise).
  3. *Core-Point DBSCAN Level-Set Equivalence (Alg. 31 $arrow.l.r.double$ Alg. 56):* Fix $(epsilon, "MinPts")$ and let $cal(X)_"core" = cal(X)_"core" (epsilon, "MinPts")$. Restricted to $cal(X)_"core"$, density-reachability is an equivalence relation whose equivalence classes (the DBSCAN core clusters) coincide *identically* with:
     (i) the connected components of the induced graph $G_epsilon [cal(X)_"core"]$;
     (ii) the connected components of $T^*_epsilon (cal(X)_"core")$, where $T^*(cal(X)_"core")$ is any MST of $(cal(X)_"core", delta)$; and
     (iii) the Single Linkage (Alg. 56) dendrogram of $cal(X)_"core"$ cut at height $epsilon$.
  4. *Exact Border-Point Dichotomy:*
     - When $"MinPts" = 1$, $cal(X)_"core" = cal(X)$ and DBSCAN equals Single Linkage cut at $epsilon$ identically on all of $cal(X)$.
     - When $"MinPts" = 2$, $cal(X)_"core"$ is precisely the set of non-isolated vertices of $G_epsilon$, so DBSCAN clusters are identically the non-singleton connected components of $T^*_epsilon$, with singletons classified as noise.
     - When $"MinPts" >= 3$, DBSCAN differs from Single Linkage on $cal(X)_"core"$ *solely* on border points $x_b in cal(X) backslash cal(X)_"core"$ adjacent in $G_epsilon$ to two or more distinct core components (which depend on vertex scan order); under DBSCAN\* (classifying non-core points as noise), the equivalence holds globally across all $cal(X)$ without order dependence.
]

#proof[
  *Part 1 (Equivalence of $G_epsilon$, $T^*_epsilon$, and Single Linkage):* Fix $epsilon >= 0$. Because $E(T^*_epsilon) subset.eq E_epsilon$, every connected component of $T^*_epsilon$ is contained in a connected component of $G_epsilon$. For the reverse inclusion, suppose $x, y in cal(X)$ are connected in $G_epsilon$; then there exists a path $P = (v_0, v_1, dots, v_r)$ with $v_0 = x, v_r = y$ and $delta(v_(s - 1), v_s) <= epsilon$ for every $s = 1, dots, r$. Consider any edge $e = {v_(s - 1), v_s}$ on $P$. If $e in E(T^*)$, then $e in E(T^*_epsilon)$. If $e in.not E(T^*)$, adding $e$ to $T^*$ creates a unique fundamental cycle $C_e = P_(T^*) (v_(s - 1), v_s) union {e}$. By the MST Cycle Property, $e$ must have weight at least as large as every other edge on $C_e$ (otherwise replacing a strictly heavier edge $e' in P_(T^*) (v_(s - 1), v_s)$ with $e$ would yield a strictly lighter spanning tree, contradicting minimality of $T^*$). Thus every edge $e' in P_(T^*) (v_(s - 1), v_s)$ satisfies $delta(e') <= delta(e) <= epsilon$, so $v_(s - 1)$ and $v_s$ are connected in $T^*_epsilon$. Concatenating these $T^*_epsilon$-paths across $s = 1, dots, r$ connects $x$ to $y$ in $T^*_epsilon$. Hence $G_epsilon$ and $T^*_epsilon$ have the *identical* connected components.

  Next compare with Single Linkage (Alg. 56): Single Linkage merges at step $t$ the pair of current clusters $(A, B)$ minimizing $d_"SL" (A, B) := min_(a in A, b in B) delta(a, b)$. By induction on $t$, this is *identically* Kruskal's MST algorithm scanning edges of $K_n$ in non-decreasing weight order and uniting distinct components whenever an edge connects two different clusters—each accepted Kruskal edge $e_t in E(T^*)$ triggers a Single Linkage merge at exact dendrogram height $h_t = delta(e_t)$ (with ties handled identically because adding all edges of weight $<= epsilon$ yields the same connected components regardless of tie order!). Therefore cutting the Single Linkage dendrogram at height $epsilon$ retains precisely the merges with $delta(e_t) <= epsilon$, yielding the exact connected components of $T^*_epsilon$.

  *Part 2 (Subdominant Ultrametric):* By Part 1, $x$ and $y$ lie in the same Single Linkage cluster at height $epsilon$ if and only if they are connected in $T^*_epsilon$ (iff $max_(e in P_(T^*)(x, y)) delta(e) <= epsilon$) and if and only if they are connected in $G_epsilon$ (iff $min_(P: x -> y) max_(e in P) delta(e) <= epsilon$). Taking the infimum over all such $epsilon$ establishes both equalities for $u_"SL" (x, y)$. Clearly $u_"SL" (x, y) <= delta(x, y)$ (taking the 1-edge path $P = {x, y}$), and concatenating bottleneck-optimal paths from $x$ to $z$ and $z$ to $y$ proves $u_"SL" (x, y) <= max{u_"SL" (x, z), u_"SL" (z, y)}$, so $u_"SL"$ is an ultrametric. Finally, let $u$ be *any* ultrametric satisfying $u <= delta$. For any path $P = (v_0, dots, v_r)$ from $x$ to $y$, repeated application of the strong triangle inequality gives
  $ u(x, y) <= max_(1 <= s <= r) u(v_(s - 1), v_s) <= max_(1 <= s <= r) delta(v_(s - 1), v_s). $
  Taking the minimum over all paths $P$ yields $u(x, y) <= u_"SL" (x, y)$, proving $u_"SL"$ is the unique maximal ultrametric dominated by $delta$.

  *Part 3 (Core-Point DBSCAN Equivalence):* Restrict attention to $cal(X)_"core"$. For $x, y in cal(X)_"core"$, direct density-reachability $x -> y$ holds iff $delta(x, y) <= epsilon$, which is symmetric because $delta$ is a metric and both $x, y in cal(X)_"core"$. Furthermore, in any density-reachability chain $x = v_0 -> v_1 -> dots -> v_r = y$ between core points, every intermediate vertex $v_0, dots, v_(r - 1)$ must be a core point by definition of direct density-reachability, and $v_r = y in cal(X)_"core"$ by hypothesis; hence the entire chain lies inside $cal(X)_"core"$ and corresponds bijectively to a path in $G_epsilon [cal(X)_"core"]$. Therefore two core points are density-connected iff they lie in the same connected component of $G_epsilon [cal(X)_"core"]$. Applying Part 1 to the finite metric subspace $(cal(X)_"core", delta)$ immediately identifies these DBSCAN core clusters with both the connected components of $T^*_epsilon (cal(X)_"core")$ and the Single Linkage cut of $cal(X)_"core"$ at height $epsilon$.

  *Part 4 (Border-Point Dichotomy):* Since $x in B_epsilon (x)$, $|B_epsilon (x)| >= 1$ always, so $"MinPts" = 1$ forces $cal(X)_"core" = cal(X)$ identically. For $"MinPts" = 2$, $|B_epsilon (x)| >= 2$ holds iff $x$ has at least one neighbor $y != x$ with $delta(x, y) <= epsilon$, i.e., iff $x$ is non-isolated in $G_epsilon$; no non-core point can have any neighbor in $G_epsilon$, so no border points exist, and all non-singleton components of $G_epsilon$ (hence $T^*_epsilon$) are exact DBSCAN clusters. For $"MinPts" >= 3$, every non-core point $x_b$ that has a core neighbor $x_c in B_epsilon (x_b) inter cal(X)_"core"$ is density-reachable from $x_c$'s core component; if $x_b$ touches only one core component, its assignment is unique, whereas if $x_b$ lies within $epsilon$ of two distinct core components $C_1, C_2$, DBSCAN assigns $x_b$ to whichever component is expanded first @schubert2017dbscanrevisited. Dropping border points into noise (DBSCAN\*) eliminates this non-uniqueness and makes the partition canonical.
]

#subsec("B", "HDBSCAN as Single Linkage on the Mutual Reachability Metric")

A classic limitation of DBSCAN (Alg. 31) is its reliance on a single global density scale $epsilon$. HDBSCAN (Alg. 32, @campello2013hdbscan, @campello2015hierarchical) resolves this by lifting the entire continuous filtration of DBSCAN\* partitions across *all* $epsilon in [0, infinity)$ into a single ultrametric tree via a density-adaptive metric transformation. Given smoothing parameter $m_"pts" in {1, dots, n}$, define the *core distance* of each $x_i in cal(X)$ as the minimal radius enclosing at least $m_"pts"$ points:
$ d_"core" (x_i) := inf {r >= 0 : |B_r (x_i)| >= m_"pts"}, $
so $x_i in cal(X)_"core" (epsilon, m_"pts")$ if and only if $d_"core" (x_i) <= epsilon$. Define the *mutual reachability distance* $d_"mreach" : cal(X) times cal(X) -> bb(R)_(>= 0)$ by $d_"mreach" (x_i, x_i) := 0$ and, for $i != j$,
$ d_"mreach" (x_i, x_j) := max{d_"core" (x_i), space d_"core" (x_j), space delta(x_i, x_j)}. $

#theorem(title: "HDBSCAN <=> Single Linkage on the Mutual Reachability Metric", origin: "Constructive Proof")[
  Let $(cal(X), delta)$ be a finite metric space and fix $m_"pts" in {1, dots, n}$ @campello2013hdbscan, @campello2015hierarchical.
  1. *Metric-Space Structure:* The pair $(cal(X), d_"mreach")$ is a genuine *metric space* (satisfying non-negativity, identity of indiscernibles, symmetry, and the triangle inequality), with $d_"mreach" (x, y) >= delta(x, y)$ everywhere.
  2. *Exact Continuous Filtration Isomorphism:* For *every* density threshold $epsilon >= 0$, the non-noise clusters of DBSCAN\* at parameters $(epsilon, m_"pts")$ coincide *identically* with the connected components of the $epsilon$-thresholded graph $G_epsilon^"mr" := (cal(X), {{x_i, x_j} : d_"mreach" (x_i, x_j) <= epsilon})$ restricted to active vertices ${x_i : d_"core" (x_i) <= epsilon}$ (equivalently, every vertex with $d_"core" (x_i) > epsilon$ is isolated in $G_epsilon^"mr"$ and labeled noise).
  3. *Single Linkage / MST Equivalence & Exact EOM Dynamic Program:* The complete HDBSCAN cluster hierarchy is *isomorphic* to the Single Linkage (Alg. 56) dendrogram—and Kruskal MST—of the metric space $(cal(X), d_"mreach")$. Furthermore, HDBSCAN's Excess-of-Mass (EOM) flat clustering solves the exact global maximum-weight antichain problem on the condensed tree $cal(T)_"cond"$:
     $ max_(cal(F) subset.eq cal(T)_"cond" " antichain") sum_(C in cal(F)) S(C), quad S(C) := sum_(x in C) (lambda_"max" (x, C) - lambda_"birth" (C)), $
     where $lambda = 1 / epsilon$, in exact $O(n)$ bottom-up dynamic programming time.
]

#proof[
  *Part 1 (Metric Axioms):* Non-negativity ($d_"mreach" (x, y) >= delta(x, y) >= 0$), symmetry (since $max{a, b, c}$ is symmetric in $a, b$ and $delta(x, y) = delta(y, x)$), and identity of indiscernibles ($d_"mreach" (x, y) = 0 ==> delta(x, y) = 0 ==> x = y$) are immediate. For the triangle inequality, take any $x, y, z in cal(X)$. If $x = z$, $d_"mreach" (x, z) = 0 <= d_"mreach" (x, y) + d_"mreach" (y, z)$; if $x = y$ or $y = z$, equality holds trivially. Assume $x, y, z$ are pairwise distinct. We must verify that each of the three arguments inside $d_"mreach" (x, z) = max{d_"core" (x), d_"core" (z), delta(x, z)}$ is bounded above by $R := d_"mreach" (x, y) + d_"mreach" (y, z)$:
  - $d_"core" (x) <= max{d_"core" (x), d_"core" (y), delta(x, y)} = d_"mreach" (x, y) <= R$, since $d_"mreach" (y, z) >= 0$.
  - $d_"core" (z) <= max{d_"core" (y), d_"core" (z), delta(y, z)} = d_"mreach" (y, z) <= R$, since $d_"mreach" (x, y) >= 0$.
  - By the triangle inequality of the underlying metric $delta$:
    $ delta(x, z) <= delta(x, y) + delta(y, z) <= d_"mreach" (x, y) + d_"mreach" (y, z) = R. $
  Taking the maximum over the three inequalities proves $d_"mreach" (x, z) <= d_"mreach" (x, y) + d_"mreach" (y, z)$. Hence $(cal(X), d_"mreach")$ is a metric space!

  *Part 2 (Filtration Equivalence):* Fix any $epsilon >= 0$. For any distinct pair $x_i != x_j$, by definition of the maximum,
  $ d_"mreach" (x_i, x_j) <= epsilon &arrow.l.r.double (d_"core" (x_i) <= epsilon) and (d_"core" (x_j) <= epsilon) \
    &quad and (delta(x_i, x_j) <= epsilon) \
    &arrow.l.r.double (x_i, x_j in cal(X)_"core" (epsilon, m_"pts")) and ({x_i, x_j} in E_epsilon). $
  Therefore, two consequences follow simultaneously for every $epsilon >= 0$:
  1. Every non-core vertex ($d_"core" (x_i) > epsilon$) has *zero* incident edges in $G_epsilon^"mr"$, so it is strictly isolated and cleanly identified as noise (even without recomputing neighborhoods!).
  2. On the core vertex subset $cal(X)_"core" (epsilon, m_"pts")$, the edge set of $G_epsilon^"mr"$ coincides *identically* with the edge set of $G_epsilon [cal(X)_"core" (epsilon, m_"pts")]$.
  By Theorem 7(3), the connected components of $G_epsilon^"mr"$ on $cal(X)_"core" (epsilon, m_"pts")$ are *identically* the DBSCAN\* clusters at $(epsilon, m_"pts")$ across the entire continuum $epsilon in [0, infinity)$. (Note that adding self-loops of weight $d_"core" (x_i)$ cleanly distinguishes singleton core clusters with $m_"pts" = 1$ or isolated cores from noise.)

  *Part 3 (Single Linkage Equivalence & Exact EOM DP):* Because $(cal(X), d_"mreach")$ is a finite metric space by Part 1, Theorem 7(1) applies verbatim: for every $epsilon >= 0$, the connected components of $G_epsilon^"mr"$ equal the connected components of the $epsilon$-thresholded MST $T^*_epsilon (cal(X), d_"mreach")$ and the height-$epsilon$ cut of Single Linkage (Alg. 56) run on $d_"mreach"$. In HDBSCAN's condensed tree $cal(T)_"cond"$ (pruning child splits smaller than $m_"clSize"$ into falling-out points), define for each node $C in cal(T)_"cond"$ the optimal subtree stability value $hat(S)(C)$ recursively from leaves upward:
  $ hat(S)(C) = cases(S(C) & "if " C " is a leaf,", max{S(C), space sum_(C' in "children"(C)) hat(S)(C')} & "otherwise.") $
  Because any valid flat clustering corresponds to an antichain (at most one selected node on any root-to-leaf path) and the children's subtrees are vertex-disjoint, induction on subtree height proves that $hat(S)("root")$ attains the *exact global maximum* of $sum_(C in cal(F)) S(C)$ over all antichains $cal(F)$ in $O(|cal(T)_"cond"|) = O(n)$ time!
]

#remark(title: "Link to OPTICS and Statistical Density-Tree Consistency")[
  Theorem 8 also embeds OPTICS (Alg. 33, @ankerst1999optics) into the same hierarchy: OPTICS executes Prim's MST traversal directly on the asymmetric reachability distance $max{d_"core" (x_i), delta(x_i, x_j)}$, whose symmetrized version is $d_"mreach"$. Furthermore, because $d_"core" (x_i)^(-d)$ is proportional to the $k$-NN density estimator $hat(f)_(m_"pts") (x_i)$, Chaudhuri et al. @chaudhuri2014consistent proved that single-linkage clustering on $d_"mreach"$ (with $m_"pts" = Theta(log n)$) achieves Hartigan's strong *cluster-tree consistency* @hartigan1981consistency for recovering the population density level-set tree of smooth densities on $bb(R)^d$, whereas raw Euclidean Single Linkage ($m_"pts" = 1$) fails consistency in $d >= 2$ due to noise-bridge chaining!
]

#subsec("C", "Mean Shift, Gaussian KDE EM, and DENCLUE 2.0")

Mode-seeking algorithms partition $bb(R)^d$ into basins of attraction of local maxima of a radial kernel density estimator (KDE) with bandwidth $h > 0$:
$ hat(f)_h (x) = (c_(k, d)) / (n h^d) sum_(i=1)^n k(norm((x - x_i) / h)_2^2), $
where $k: [0, infinity) -> [0, infinity)$ is a continuously differentiable, strictly decreasing, convex profile with shadow profile $g(r) := -k'(r) > 0$. Setting $nabla hat(f)_h (x) = bold(0)$ yields the fixed-point Mean Shift operator (Alg. 34, @cheng1995mean, @comaniciu2002mean):
$ T_"MS" (x) := (sum_(i=1)^n x_i g(norm((x - x_i) / h)_2^2)) / (sum_(i=1)^n g(norm((x - x_i) / h)_2^2)), quad m(x) := T_"MS" (x) - x. $
For the Gaussian kernel $k(r) = exp(-r / 2)$, we have $g(r) = 1 / 2 k(r)$ and $hat(f)_h (x) = n^(-1) sum_(i=1)^n cal(N)(x | x_i, h^2 I_d)$. Independently, Hinneburg and Gabriel @hinneburg2007denclue2 replaced the step-size-dependent gradient ascent of DENCLUE 1.0 (Alg. 35, @hinneburg1998denclue) with the step-size-free DENCLUE 2.0 hill-climbing update (Alg. 43). Theorem 9 proves the exact triple equivalence among Mean Shift, EM for Gaussian KDE mode-seeking @carreira2006fast, and DENCLUE 2.0, armed with a razor-sharp Lyapunov ascent identity.

#theorem(title: "Mean Shift <=> EM Monotone Ascent for Gaussian KDE <=> DENCLUE 2.0", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and $h > 0$.
  1. *Exact EM Isomorphism for Gaussian KDE (Alg. 34 $arrow.l.r.double$ Alg. 43):* Viewing $hat(f)_h (x) = sum_(i=1)^n n^(-1) cal(N)(x | x_i, h^2 I_d)$ as a fixed-parameter mixture over latent kernel index $Z in {1, dots, n}$, maximizing $log hat(f)_h (x)$ over $x in bb(R)^d$ by EM produces the exact iterate
     $ &w_i^((t)) = (exp(- norm(x^((t)) - x_i)_2^2 / (2 h^2))) / (sum_(j=1)^n exp(- norm(x^((t)) - x_j)_2^2 / (2 h^2))), \
       &x^((t+1)) = sum_(i=1)^n w_i^((t)) x_i = T_"MS" (x^((t))), $
     which coincides *identically* with Gaussian Mean Shift (Alg. 34) and DENCLUE 2.0 (Alg. 43) @carreira2006fast, @hinneburg2007denclue2.
  2. *Exact Quantitative Lyapunov Identity (Gaussian Case):* At every iteration $t >= 0$,
     $ &log hat(f)_h (x^((t+1))) - log hat(f)_h (x^((t))) \
       &quad = 1 / (2 h^2) norm(x^((t+1)) - x^((t)))_2^2 + "KL"(w^((t)) || w^((t+1))) \
       &quad >= 1 / (2 h^2) norm(x^((t+1)) - x^((t)))_2^2. $
     Consequently, $sum_(t=0)^infinity norm(x^((t+1)) - x^((t)))_2^2 <= 2 h^2 log(hat(f)_"max" / hat(f)_h (x^((0)))) < infinity$, and every limit point $x^*$ lies in $"conv"(cal(X))$ and satisfies $nabla hat(f)_h (x^*) = bold(0)$.
  3. *General Convex-Profile MM Ascent:* For *any* convex, strictly decreasing profile $k$ with $g = -k' > 0$,
     $ hat(f)_h (x^((t+1))) - hat(f)_h (x^((t))) >= (c_(k, d)) / (n h^(d + 2)) (sum_(i=1)^n g_i^((t))) norm(x^((t+1)) - x^((t)))_2^2 >= 0, $
     where $g_i^((t)) := g(norm((x^((t)) - x_i) / h)_2^2)$.
  4. *Exact Fixed-Point Set Identity:* For every smooth radial kernel with $g > 0$,
     $ x^* = T_"MS" (x^*) arrow.l.r.double nabla hat(f)_h (x^*) = bold(0), $
     so the fixed points of Mean Shift, Gaussian KDE EM, and DENCLUE 1.0/2.0 coincide identically.
]

#proof[
  *Part 1 (EM Derivation):* Define the joint density $p(x, Z = i) := n^(-1) (2 pi h^2)^(-d / 2) exp(- norm(x - x_i)_2^2 / (2 h^2))$, whose marginal over $Z$ is $p(x) = hat(f)_h (x)$. At current iterate $x^((t))$, the E-step computes the posterior $w_i^((t)) = p(Z = i | x^((t))) = p(x^((t)), i) / hat(f)_h (x^((t)))$, yielding the stated softmax weights $w_i^((t))$. The EM auxiliary function is
  $ Q(x | x^((t))) &:= sum_(i=1)^n w_i^((t)) log p(x, i) \
    &= C_h - 1 / (2 h^2) sum_(i=1)^n w_i^((t)) norm(x - x_i)_2^2, $
  where $C_h := -log n - (d / 2) log(2 pi h^2)$. Since $Q(x | x^((t)))$ is a strictly concave isotropic quadratic in $x$ with $sum_(i=1)^n w_i^((t)) = 1$, its unique global maximizer is given by setting $nabla_x Q(x | x^((t))) = -h^(-2) (x - sum_i w_i^((t)) x_i) = bold(0)$, i.e., $x^((t+1)) = sum_(i=1)^n w_i^((t)) x_i$, which is identically both Gaussian Mean Shift and the DENCLUE 2.0 update!

  *Part 2 (Exact Lyapunov Ascent Identity):* For any distribution $q = (q_1, dots, q_n)$ on ${1, dots, n}$ and any $x in bb(R)^d$, write the exact variational log-likelihood identity:
  $ log hat(f)_h (x) &= sum_(i=1)^n q_i log(p(x, i) / q_i) + sum_(i=1)^n q_i log(q_i / p(i | x)) \
    &= Q_q (x) + H(q) + "KL"(q || p(dot | x)), $
  where $Q_q (x) = sum_i q_i log p(x, i)$. Choose $q = w^((t)) = p(dot | x^((t)))$. Evaluating at $x = x^((t))$ gives $log hat(f)_h (x^((t))) = Q(x^((t)) | x^((t))) + H(w^((t)))$ (since $"KL"(w^((t)) || w^((t))) = 0$), and evaluating at $x = x^((t+1))$ gives $log hat(f)_h (x^((t+1))) = Q(x^((t+1)) | x^((t))) + H(w^((t))) + "KL"(w^((t)) || w^((t+1)))$. Subtracting:
  $ &log hat(f)_h (x^((t+1))) - log hat(f)_h (x^((t))) \
    &quad = [Q(x^((t+1)) | x^((t))) - Q(x^((t)) | x^((t)))] + "KL"(w^((t)) || w^((t+1))). $
  By the weighted Parallel-Axis identity around the weighted mean $x^((t+1)) = sum_(i=1)^n w_i^((t)) x_i$:
  $ sum_(i=1)^n w_i^((t)) norm(x^((t)) - x_i)_2^2 = sum_(i=1)^n w_i^((t)) norm(x^((t+1)) - x_i)_2^2 + norm(x^((t)) - x^((t+1)))_2^2. $
  Dividing by $2 h^2$ yields $Q(x^((t+1)) | x^((t))) - Q(x^((t)) | x^((t))) = (2 h^2)^(-1) norm(x^((t+1)) - x^((t)))_2^2$ with *exact equality*, proving the stated identity! Since $"KL"(dot || dot) >= 0$ and $hat(f)_h (x) <= hat(f)_"max" := (2 pi h^2)^(-d / 2)$, summing from $t = 0$ to $T - 1$ telescopes to
  $ sum_(t=0)^(T - 1) norm(x^((t+1)) - x^((t)))_2^2 &<= 2 h^2 [log hat(f)_h (x^((T))) - log hat(f)_h (x^((0)))] \
    &<= 2 h^2 log(hat(f)_"max" / hat(f)_h (x^((0)))). $
  Sending $T -> infinity$ forces $norm(x^((t+1)) - x^((t)))_2 -> 0$. Because $x^((t+1)) = sum_i w_i^((t)) x_i in "conv"(cal(X))$ (which is compact), the sequence admits limit points, and by continuity of $T_"MS"$, every limit point $x^*$ satisfies $x^* = T_"MS" (x^*)$ (and even stronger, finite-length convergence holds when modes are isolated @ariascastro2016meanshift).

  *Part 3 (General Convex Profile MM Bound):* For any convex $k$, the supporting tangent inequality $k(r_2) - k(r_1) >= k'(r_1)(r_2 - r_1) = -g(r_1)(r_2 - r_1)$ applied at $r_1 = norm((x^((t)) - x_i) / h)_2^2$ and $r_2 = norm((x^((t+1)) - x_i) / h)_2^2$ gives
  $ &hat(f)_h (x^((t+1))) - hat(f)_h (x^((t))) \
    &quad >= (c_(k, d)) / (n h^(d + 2)) sum_(i=1)^n g_i^((t)) [norm(x^((t)) - x_i)_2^2 - norm(x^((t+1)) - x_i)_2^2]. $
  Since $x^((t+1)) = (sum_i g_i^((t)) x_i) / (sum_i g_i^((t)))$, applying the weighted Parallel-Axis identity with weights $g_i^((t))$ gives
  $ &sum_(i=1)^n g_i^((t)) [norm(x^((t)) - x_i)_2^2 - norm(x^((t+1)) - x_i)_2^2] \
    &quad = (sum_(i=1)^n g_i^((t))) norm(x^((t+1)) - x^((t)))_2^2, $
  establishing Part 3 directly!

  *Part 4 (Fixed-Point Set Identity):* Differentiating $hat(f)_h (x)$ yields
  $ nabla hat(f)_h (x) &= - (2 c_(k, d)) / (n h^(d + 2)) sum_(i=1)^n g(norm((x - x_i) / h)_2^2) (x - x_i) \
    &= (2 c_(k, d)) / (n h^(d + 2)) (sum_(i=1)^n g(norm((x - x_i) / h)_2^2)) (T_"MS" (x) - x). $
  Because $c_(k, d) > 0$ and $sum_(i=1)^n g(norm((x - x_i) / h)_2^2) > 0$, we have $nabla hat(f)_h (x^*) = bold(0)$ if and only if $T_"MS" (x^*) = x^*$.
]

#subsec("D", "Ward's Method as Greedy Agglomerative K-Means WCSS Minimization")

Ward's hierarchical clustering method (Alg. 59, @ward1963hierarchical) appears structurally distant from Lloyd's $k$-means (Alg. 1): one constructs a nested sequence of $n$ partitions bottom-up via the Lance–Williams recursion (Alg. 72, @lance1967general), while the other refines a fixed $k$-partition by Voronoi reassignment. Theorem 10 proves that Ward's linkage distance is the *exact* one-step WCSS increment of merging two clusters, derives the exact Lance–Williams coefficients from first principles, and proves inversion-freedom and exact WCSS telescoping.

For any finite non-empty cluster $A subset.eq cal(X) subset.eq bb(R)^d$ of cardinality $n_A = |A|$ and centroid $mu_A = n_A^(-1) sum_(x in A) x$, denote its within-cluster sum of squares by $E(A) := sum_(x in A) norm(x - mu_A)_2^2$. For a partition $cal(C)$, the $k$-means objective is $J(cal(C)) = sum_(A in cal(C)) E(A)$.

#theorem(title: "Ward's Method <=> Greedy Agglomerative K-Means WCSS Minimization", origin: "Original Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and define the Ward merge-cost functional $Delta(A, B) := E(A union B) - E(A) - E(B)$ for disjoint non-empty subsets $A, B subset.eq cal(X)$ @ward1963hierarchical, @lance1967general, @mullner2011modern.
  1. *Exact Huygens–Steiner Variance Delta:* For all disjoint non-empty $A, B subset.eq cal(X)$,
     $ Delta(A, B) = (n_A n_B) / (n_A + n_B) norm(mu_A - mu_B)_2^2. $
     Consequently, merging the pair $(A, B)$ minimizing $Delta(A, B)$ at every agglomerative step is *identically* greedy step-by-step minimization of the $k$-means WCSS objective $J(cal(C))$.
  2. *Exact Lance–Williams Recursion:* For any three pairwise disjoint non-empty clusters $A, B, C subset.eq cal(X)$ with $N := n_A + n_B + n_C$, the Ward dissimilarity $D(A, B) := Delta(A, B)$ (and equivalently $2 Delta(A, B)$) satisfies the *exact* Lance–Williams linear update:
     $ D(A union B, C) &= (n_A + n_C) / N D(A, C) + (n_B + n_C) / N D(B, C) \
       &quad - n_C / N D(A, B). $
  3. *Reducibility and Zero Dendrogram Inversions:* Whenever $D(A, B) <= min{D(A, C), D(B, C)}$, we have $D(A union B, C) >= D(A, B)$. Hence the Ward merge heights $0 <= Delta_1 <= Delta_2 <= dots <= Delta_(n - 1)$ are *strictly non-decreasing* (ultrametric, zero inversions), and the Nearest-Neighbor Chain algorithm computes the exact Ward dendrogram in $O(n^2)$ time.
  4. *Exact Total-Inertia Telescoping:* Let $cal(C)_k^"Ward"$ be the $k$-cluster cut of Ward's dendrogram ($1 <= k <= n$). Then its exact $k$-means WCSS equals the partial sum of the first $n - k$ merge heights, and the full sum equals total dataset inertia:
     $ J(cal(C)_k^"Ward") = sum_(t=1)^(n - k) Delta_t, quad sum_(t=1)^(n - 1) Delta_t = sum_(i=1)^n norm(x_i - overline(x))_2^2. $
]

#proof[
  *Part 1 (Huygens–Steiner Merge Delta):* The merged centroid is $mu_(A union B) = (n_A mu_A + n_B mu_B) / (n_A + n_B)$, so
  $ &mu_A - mu_(A union B) = n_B / (n_A + n_B) (mu_A - mu_B), \
    &mu_B - mu_(A union B) = - n_A / (n_A + n_B) (mu_A - mu_B). $
  Applying the Parallel-Axis identity separately to $A$ and $B$ around $c = mu_(A union B)$:
  $ E(A union B) &= sum_(x in A) norm(x - mu_(A union B))_2^2 + sum_(x in B) norm(x - mu_(A union B))_2^2 \
    &= E(A) + n_A norm(mu_A - mu_(A union B))_2^2 \
    &quad + E(B) + n_B norm(mu_B - mu_(A union B))_2^2. $
  Substituting the two centroid displacements:
  $ Delta(A, B) &= [n_A (n_B / (n_A + n_B))^2 + n_B (n_A / (n_A + n_B))^2] norm(mu_A - mu_B)_2^2 \
    &= (n_A n_B (n_B + n_A)) / ((n_A + n_B)^2) norm(mu_A - mu_B)_2^2 = (n_A n_B) / (n_A + n_B) norm(mu_A - mu_B)_2^2. $
  When merging $A$ and $B$ in a current partition $cal(C)$, all other clusters remain unchanged, so $J(cal(C)') - J(cal(C)) = Delta(A, B)$. Minimizing $J(cal(C)')$ greedily over all pairs $(A, B)$ in $cal(C)$ therefore coincides identically with Ward's merge rule!

  *Part 2 (Exact Lance–Williams Proof):* Define displacement vectors $u := mu_A - mu_C$ and $v := mu_B - mu_C$, noting that $u - v = mu_A - mu_B$ and $mu_(A union B) - mu_C = (n_A u + n_B v) / (n_A + n_B)$. First verify the elementary algebraic vector identity
  $ &norm(n_A u + n_B v)_2^2 = n_A^2 norm(u)_2^2 + n_B^2 norm(v)_2^2 + 2 n_A n_B chevron.l u, v chevron.r \
    &quad = n_A (n_A + n_B) norm(u)_2^2 + n_B (n_A + n_B) norm(v)_2^2 - n_A n_B norm(u - v)_2^2, $
  which holds identically because $-n_A n_B norm(u - v)_2^2 = -n_A n_B norm(u)_2^2 - n_A n_B norm(v)_2^2 + 2 n_A n_B chevron.l u, v chevron.r$. Now evaluate $D(A union B, C)$ using Part 1 with $n_(A union B) = n_A + n_B$ and $N = n_A + n_B + n_C$:
  $ D(A union B, C) &= ((n_A + n_B) n_C) / N norm((n_A u + n_B v) / (n_A + n_B))_2^2 \
    &= n_C / ((n_A + n_B) N) norm(n_A u + n_B v)_2^2 \
    &= (n_A n_C) / N norm(u)_2^2 + (n_B n_C) / N norm(v)_2^2 \
    &quad - (n_A n_B n_C) / ((n_A + n_B) N) norm(u - v)_2^2. $
  Rewriting each term in terms of $D(A, C) = (n_A n_C / (n_A + n_C)) norm(u)_2^2$, $D(B, C) = (n_B n_C / (n_B + n_C)) norm(v)_2^2$, and $D(A, B) = (n_A n_B / (n_A + n_B)) norm(u - v)_2^2$ immediately yields
  $ D(A union B, C) &= (n_A + n_C) / N D(A, C) + (n_B + n_C) / N D(B, C) \
    &quad - n_C / N D(A, B), $
  establishing the exact Lance–Williams coefficients $alpha_A = (n_A + n_C) / N$, $alpha_B = (n_B + n_C) / N$, $beta = -n_C / N$, $gamma = 0$.

  *Part 3 (Reducibility & Inversion-Freedom):* Suppose $D(A, B) <= min{D(A, C), D(B, C)}$. Because $alpha_A > 0$, $alpha_B > 0$, and $alpha_A + alpha_B + beta = (n_A + n_C + n_B + n_C - n_C) / N = 1$, substituting $D(A, C) >= D(A, B)$ and $D(B, C) >= D(A, B)$ into Part 2 gives
  $ D(A union B, C) >= (alpha_A + alpha_B + beta) D(A, B) = D(A, B). $
  Thus every newly formed dissimilarity after merging a mutually closest pair $(A, B)$ is at least $D(A, B)$, proving both Bruynooghe's reducibility property (the exact sufficient condition for $O(n^2)$ Nearest-Neighbor Chain correctness @mullner2011modern) and monotonicity $Delta_1 <= Delta_2 <= dots <= Delta_(n - 1)$!

  *Part 4 (Telescoping):* At initialization $cal(C)_n = {{x_1}, dots, {x_n}}$, every singleton has $E({x_i}) = 0$, so $J(cal(C)_n) = 0$. Each merge step $t$ increases $J$ by exactly $Delta_t$, whence $J(cal(C)_(n - r)^"Ward") = sum_(t=1)^r Delta_t$. Setting $r = n - k$ gives $J(cal(C)_k^"Ward") = sum_(t=1)^(n - k) Delta_t$, and setting $k = 1$ gives $J({cal(X)}) = sum_(i=1)^n norm(x_i - overline(x))_2^2 = sum_(t=1)^(n - 1) Delta_t$.
]

#remark(title: "Ward Initialization for Lloyd's K-Means")[
  Because Ward's Method commits permanently to hierarchical merges without point reassignment, $J(cal(C)_k^"Ward") >= J(cal(C)_k^*)$ strictly in general; however, seeding Lloyd's $k$-means (Alg. 1) with Ward's $k$-cluster centroids initializes Lloyd directly at a structured low-WCSS nested partition and guarantees monotone WCSS improvement $J(cal(C)_k^"Ward+Lloyd") <= J(cal(C)_k^"Ward")$ by Theorem 1 of Alg. 1.
]

#sec("V", "Matrix Factorization, Topic Models, and Community Detection")

We next prove two landmark exact isomorphisms connecting latent-variable topic models, non-negative matrix factorization, statistical physics, and generative random-graph models across Parts IV and V.

#subsec("A", "KL-NMF and Probabilistic Latent Semantic Analysis (pLSA)")

Let $X in bb(R)_(>= 0)^(m times n)$ be a non-negative term–document co-occurrence matrix with strictly positive total count $N := sum_(i=1)^m sum_(j=1)^n X_(i j) > 0$.
- *KL-NMF (Alg. 92, @lee1999nmf, @lee2001algorithms)* seeks non-negative matrices $W in bb(R)_(>= 0)^(m times k)$ and $H in bb(R)_(>= 0)^(k times n)$ minimizing the generalized Kullback–Leibler (I-) divergence
  $ D_"KL" (X || W H) &:= sum_(i=1)^m sum_(j=1)^n (X_(i j) log((X_(i j)) / ((W H)_(i j))) \
    &quad - X_(i j) + (W H)_(i j)), $
  with convention $0 log 0 = 0$.
- *pLSA (Alg. 84, @hofmann1999plsa)* models the joint discrete probability $hat(P)_(i j) = P(w_i, d_j) = sum_(r=1)^k P(z_r) P(w_i | z_r) P(d_j | z_r) = (U S V^T)_(i j)$, where $U in bb(R)_(>= 0)^(m times k)$ has unit column sums ($bold(1)_m^T U = bold(1)_k^T$), $V in bb(R)_(>= 0)^(n times k)$ has unit column sums ($bold(1)_n^T V = bold(1)_k^T$), and $S = "diag"(s_1, dots, s_k)$ with $s in Delta_k$ ($sum_(r=1)^k s_r = 1$), maximizing the multinomial log-likelihood
  $ cal(L)_"pLSA" (U, S, V) := sum_(i=1)^m sum_(j=1)^n X_(i j) log((U S V^T)_(i j)). $

#theorem(title: "KL-NMF <=> pLSA Maximum Likelihood Equivalence", origin: "Classical — Cited")[
  Let $X in bb(R)_(>= 0)^(m times n)$ with $N = sum_(i, j) X_(i j) > 0$ @ding2008nmfplsa.
  1. *Exact Mass Conservation of KL-NMF:* Every KKT stationary point $(W^*, H^*)$ of $D_"KL" (X || W H)$—as well as every iterate after a single Lee–Seung multiplicative update of $W$ or $H$—satisfies exact total mass conservation:
     $ sum_(i=1)^m sum_(j=1)^n (W H)_(i j) = sum_(i=1)^m sum_(j=1)^n X_(i j) = N. $
  2. *Exact Objective Isomorphism (Alg. 92 $arrow.l.r.double$ Alg. 84):* On the mass-conserving manifold $cal(M)_N := {(W, H) >= 0 : sum_(i, j) (W H)_(i j) = N}$, define the scale-canonicalization map $(W, H) |-> (U, S, V)$ by setting $d_r^W := sum_(i=1)^m W_(i r)$, $d_r^H := sum_(j=1)^n H_(r j)$, $U_(i r) := W_(i r) / d_r^W$, $V_(j r) := H_(r j) / d_r^H$, and $s_r := (d_r^W d_r^H) / N$ (for active factors $d_r^W d_r^H > 0$, setting $s_r = 0$ otherwise). Then $W H = N U S V^T$, $(U, S, V)$ is a valid pLSA parameterization, and
     $ D_"KL" (X || W H) = - cal(L)_"pLSA" (U, S, V) + sum_(i=1)^m sum_(j=1)^n X_(i j) log((X_(i j)) / N). $
     Consequently, KL-NMF and pLSA have the *identical* set of global optima and KKT stationary points (modulo the positive diagonal scale gauge $(W, H) |-> (W D, D^(-1) H)$).
  3. *Alternating vs. Simultaneous MM Dynamics:* Both Lee–Seung KL-NMF and pLSA EM optimize the identical auxiliary surrogate induced by Jensen's inequality with posterior $q_(i j r) = (W_(i r) H_(r j)) / (W H)_(i j)$, but pLSA EM updates $(U, S, V)$ *simultaneously* from one frozen $q$, whereas Lee–Seung KL-NMF updates $H$ and $W$ *alternatingly* (recomputing $q$ between half-steps), yielding a strictly tighter monotonic descent step per full pass.
]

#proof[
  *Part 1 (Mass Conservation):* Let $(W^*, H^*)$ be any KKT stationary point of $D_"KL" (X || W H)$ on the non-negative orthant. Consider the one-dimensional radial perturbation $phi(alpha) := D_"KL" (X || alpha W^* H^*)$ for $alpha > 0$, which stays strictly feasible since $alpha W^* >= 0$. Differentiating at $alpha = 1$:
  $ 0 = phi'(1) &= sum_(i=1)^m sum_(j=1)^n (- (X_(i j)) / alpha + (W^* H^*)_(i j)) |_(alpha = 1) \
    &= - N + sum_(i, j) (W^* H^*)_(i j), $
  proving $sum_(i, j) (W^* H^*)_(i j) = N$. Moreover, examine a single Lee–Seung multiplicative update of $W$ holding $H$ fixed @lee2001algorithms:
  $ W'_(i r) = W_(i r) (sum_(j=1)^n H_(r j) X_(i j) / (W H)_(i j)) / (sum_(j=1)^n H_(r j)). $
  Summing $(W' H)_(i j)$ over all $(i, j)$ cancels the denominator $sum_j H_(r j)$ identically:
  $ sum_(i=1)^m sum_(j=1)^n (W' H)_(i j) &= sum_(i=1)^m sum_(r=1)^k W'_(i r) (sum_(j=1)^n H_(r j)) \
    &= sum_(i=1)^m sum_(j=1)^n (X_(i j)) / ((W H)_(i j)) (sum_(r=1)^k W_(i r) H_(r j)) \
    &= sum_(i, j) X_(i j) = N! $
  By exact symmetry, a single $H$-update likewise enforces $sum_(i, j) (W H')_(i j) = N$ identically.

  *Part 2 (Objective & Stationary-Point Bijection):* Given $(W, H) in cal(M)_N$, verify first that $U, V$ have unit column sums by definition, $s_r >= 0$, and
  $ sum_(r=1)^k s_r = 1 / N sum_(r=1)^k (sum_(i=1)^m W_(i r)) (sum_(j=1)^n H_(r j)) = 1 / N sum_(i, j) (W H)_(i j) = 1. $
  Furthermore, for every $(i, j)$,
  $ N (U S V^T)_(i j) &= N sum_(r=1)^k (W_(i r) / d_r^W) ((d_r^W d_r^H) / N) (H_(r j) / d_r^H) \
    &= sum_(r=1)^k W_(i r) H_(r j) = (W H)_(i j). $
  Conversely, any valid pLSA triplet $(U, S, V)$ maps into $cal(M)_N$ via $W := N U S^(1/2)$, $H := S^(1/2) V^T$ (or $W = N U S, H = V^T$), giving a bijection between pLSA parameters and equivalence classes of $cal(M)_N$ modulo $(W D, D^(-1) H)$ for diagonal $D succ 0$. Substituting $(W H)_(i j) = N (U S V^T)_(i j)$ and $sum_(i, j) (W H)_(i j) = N = sum_(i, j) X_(i j)$ into $D_"KL" (X || W H)$:
  $ D_"KL" (X || W H) &= sum_(i, j) [X_(i j) log(X_(i j) / N) - X_(i j) log((U S V^T)_(i j)) \
    &quad - X_(i j) + (W H)_(i j)] \
    &= - cal(L)_"pLSA" (U, S, V) + sum_(i, j) X_(i j) log(X_(i j) / N), $
  since the linear terms $-sum X_(i j) + sum (W H)_(i j) = -N + N = 0$ cancel identically! Because Part 1 guarantees that all global and local KKT stationary points of unconstrained KL-NMF lie inside $cal(M)_N$, global optima and KKT stationary points of KL-NMF and pLSA coincide identically.

  *Part 3 (Surrogate / EM Comparison):* Introduce auxiliary latent probabilities $q_(i j r) := (W_(i r) H_(r j)) / (W H)_(i j)$ ($sum_r q_(i j r) = 1$). By concavity of $log$,
  $ - log((tilde(W) tilde(H))_(i j)) <= - sum_(r=1)^k q_(i j r) log((tilde(W)_(i r) tilde(H)_(r j)) / q_(i j r)), $
  yielding the tight separable majorizer $G(tilde(W), tilde(H) | W, H) >= D_"KL" (X || tilde(W) tilde(H))$ with equality at $(tilde(W), tilde(H)) = (W, H)$. Differentiating $G$ with respect to $tilde(W)_(i r)$ and $tilde(H)_(r j)$ gives:
  - In *pLSA EM*, minimizing $G(N tilde(U) tilde(S) tilde(V)^T | W, H)$ *jointly* over $(tilde(U), tilde(S), tilde(V))$ for fixed $q$ produces the closed-form simultaneous M-step $tilde(s)_r = N^(-1) sum_(i, j) X_(i j) q_(i j r)$, $tilde(U)_(i r) prop sum_j X_(i j) q_(i j r)$, $tilde(V)_(j r) prop sum_i X_(i j) q_(i j r)$.
  - In *Lee–Seung KL-NMF*, minimizing $G(W, tilde(H) | W, H)$ over $tilde(H)$ alone yields the exact $H$-update, then *re-evaluating* $q'$ using $(W, H')$ and minimizing $G(tilde(W), H' | W, H')$ yields the $W$-update. Because re-tightening the surrogate at $(W, H')$ strictly reduces the upper bound before the second half-step, alternating KL-NMF achieves strictly steeper per-epoch auxiliary reduction while converging to the exact same stationary manifold @ding2008nmfplsa.
]

#remark(title: "Dirichlet-Regularized Extension to LDA")[
  Adding symmetric Dirichlet log-priors $(alpha - 1) sum_(r, j) log H_(r j) + (beta - 1) sum_(i, r) log W_(i r)$ to the pLSA / KL-NMF objective directly yields the MAP objective of *Latent Dirichlet Allocation* (LDA, Alg. 83, @blei2003lda), placing KL-NMF, pLSA, and LDA along one continuous probabilistic-factorization spectrum.
]

#subsec("B", "Modularity, Potts Spinglass Hamiltonian, and Uniform-Planting DC-SBM MLE")

In Part V, community detection appears divided between three separate scientific traditions: (i) combinatorial *Newman–Girvan modularity* maximization—optimized by Louvain (Alg. 107, @blondel2008louvain), Leiden (Alg. 108, @traag2019leiden), Clauset–Newman–Moore greedy modularity (Alg. 115, @clauset2004cnm), and Newman's leading eigenvector method (Alg. 116, @newman2006leading); (ii) statistical-mechanics *Potts Spinglass* ground-state minimization (Alg. 117, @reichardt2006spinglass); and (iii) probabilistic generative inference under the *Degree-Corrected Stochastic Block Model* (DC-SBM, Alg. 131, @karrer2011dcsbm, @abbe2018sbmreview). Theorem 12 proves that all three are mathematically isomorphic!

Let $G = (V, E, A)$ be an undirected graph (possibly weighted or multigraph) with symmetric adjacency $A_(i j) >= 0$, degrees $d_i = sum_(j=1)^n A_(i j)$, total volume $2 m := sum_(i=1)^n d_i > 0$, and community assignment vector $g = (g_1, dots, g_n) in {1, dots, k}^n$. Define the *resolution-$gamma$ Newman–Girvan modularity* (@newman2004finding, @reichardt2006spinglass):
$ Q_gamma (g) := 1 / (2 m) sum_(i, j = 1)^n (A_(i j) - gamma (d_i d_j) / (2 m)) delta(g_i, g_j), $
where $delta(g_i, g_j) = bb(1)_(g_i = g_j)$ is the Kronecker delta and $gamma = 1$ recovers standard Louvain/Leiden modularity.

#theorem(title: "Modularity <=> Potts Spinglass Hamiltonian <=> Uniform-Planting DC-SBM MLE", origin: "Original Unified Proof")[
  Let $G = (V, E, A)$ be an undirected graph with total edge volume $2 m = sum_i d_i > 0$ and partition $g in {1, dots, k}^n$.
  1. *Modularity <=> Reichardt–Bornholdt Potts Spinglass (Algs. 107–108 $arrow.l.r.double$ Alg. 117):* Define the $k$-state Potts Spinglass Hamiltonian @reichardt2006spinglass with ferromagnetic edge rewards $A_(i j)$ and configuration-model antiferromagnetic couplings $gamma p_(i j) = gamma d_i d_j / (2 m)$:
     $ cal(H)_"Potts" (g) := - sum_(1 <= i < j <= n) (A_(i j) - gamma (d_i d_j) / (2 m)) delta(g_i, g_j). $
     For any simple graph ($A_(i i) = 0$),
     $ cal(H)_"Potts" (g) = - m Q_gamma (g) - gamma / (4 m) sum_(i=1)^n d_i^2, $
     where the second term is strictly independent of $g$ (and vanishes if diagonal self-terms $i = j$ are included in $cal(H)_"Potts"$). Hence ground states of the Potts Spinglass Hamiltonian coincide *identically* with maximizers of resolution-$gamma$ modularity.
  2. *Modularity <=> Uniform-Planting DC-SBM MLE (Algs. 107–108 $arrow.l.r.double$ Alg. 131):* Consider the Poisson Degree-Corrected Stochastic Block Model @karrer2011dcsbm with independent edge counts $A_(i j) tilde.op "Poisson"(theta_i theta_j omega_(g_i g_j))$ ($i < j$) and $A_(i i) / 2 tilde.op "Poisson"(1 / 2 theta_i^2 omega_(g_i g_i))$, parameterized by a uniformly planted assortative affinity matrix
     $ omega_(r s) = cases(omega_"in" & "if " r = s,, omega_"out" & "if " r != s,) quad (omega_"in" > omega_"out" > 0), $
     and degree propensities $theta_i = d_i / sqrt(2 m)$. Then the exact Poisson DC-SBM log-likelihood $ell_"DC" (g)$ satisfies
     $ ell_"DC" (g) = m log(omega_"in" / omega_"out") Q_(gamma^*) (g) + C(A, omega_"in", omega_"out"), $
     where $C(A, omega_"in", omega_"out")$ is strictly independent of $g$ and the *exact critical resolution parameter* is the logarithmic mean of $(omega_"in", omega_"out")$:
     $ gamma^* = (omega_"in" - omega_"out") / (log omega_"in" - log omega_"out") in (sqrt(omega_"in" omega_"out"), space (omega_"in" + omega_"out") / 2). $
     Consequently, maximizing resolution-$gamma^*$ modularity is *identically equivalent* to Maximum Likelihood Estimation of the uniformly planted Poisson DC-SBM!
]

#proof[
  *Part 1 (Potts Hamiltonian Identity):* Because $A_(i j) = A_(j i)$ and $delta(g_i, g_j) = delta(g_j, g_i)$, splitting the full double sum in $Q_gamma (g)$ into off-diagonal ($i != j$) and diagonal ($i = j$) terms gives
  $ 2 m Q_gamma (g) &= sum_(i != j) (A_(i j) - gamma (d_i d_j) / (2 m)) delta(g_i, g_j) \
    &quad + sum_(i=1)^n (A_(i i) - gamma (d_i^2) / (2 m)) delta(g_i, g_i) \
    &= 2 sum_(1 <= i < j <= n) (A_(i j) - gamma (d_i d_j) / (2 m)) delta(g_i, g_j) \
    &quad + sum_(i=1)^n A_(i i) - gamma / (2 m) sum_(i=1)^n d_i^2. $
  Dividing by $-2$ and substituting $A_(i i) = 0$ for simple graphs yields
  $ cal(H)_"Potts" (g) = - m Q_gamma (g) - gamma / (4 m) sum_(i=1)^n d_i^2. $
  Since $m > 0$ and $-gamma (4 m)^(-1) sum_i d_i^2$ depends solely on the degree sequence and not on $g$, minimizing $cal(H)_"Potts" (g)$ is identically equivalent to maximizing $Q_gamma (g)$.

  *Part 2 (Uniform-Planting DC-SBM MLE):* Under the independent Poisson multigraph DC-SBM @karrer2011dcsbm, letting $lambda_(i j) := theta_i theta_j omega_(g_i g_j)$, the log-likelihood of observing $A$ given $g$ is
  $ ell_"DC" (g) &= sum_(1 <= i < j <= n) [A_(i j) log lambda_(i j) - lambda_(i j) - log(A_(i j)!)] \
    &quad + sum_(i=1)^n [(A_(i i) / 2) log(lambda_(i i) / 2) - lambda_(i i) / 2 - log((A_(i i) / 2)!)] \
    &= 1 / 2 sum_(i, j = 1)^n [A_(i j) log(theta_i theta_j omega_(g_i g_j)) - theta_i theta_j omega_(g_i g_j)] + C_0 (A), $
  where $C_0 (A)$ collects all factorial and factor-of-$2$ constants depending only on $A$. Next, rewrite both $omega_(g_i g_j)$ and $log omega_(g_i g_j)$ linearly in $delta(g_i, g_j)$:
  $ &omega_(g_i g_j) = omega_"out" + (omega_"in" - omega_"out") delta(g_i, g_j), \
    &log omega_(g_i g_j) = log omega_"out" + log(omega_"in" / omega_"out") delta(g_i, g_j). $
  Substitute $theta_i theta_j = d_i d_j / (2 m)$ together with these two linear identities into $ell_"DC" (g)$:
  $ ell_"DC" (g) &= C_0 (A) + 1 / 2 sum_(i, j = 1)^n [A_(i j) log(theta_i theta_j omega_"out") - (d_i d_j) / (2 m) omega_"out"] \
    &quad + 1 / 2 sum_(i, j = 1)^n [A_(i j) log(omega_"in" / omega_"out") - (omega_"in" - omega_"out") (d_i d_j) / (2 m)] delta(g_i, g_j). $
  Observe that every term in the first line depends exclusively on $(A, d, omega_"out")$ and is *completely independent* of the partition $g$! Calling that constant $C(A, omega_"in", omega_"out")$ and factoring out $m log(omega_"in" / omega_"out") > 0$ (since $omega_"in" > omega_"out" > 0$) from the second line:
  $ ell_"DC" (g) &= C(A, omega_"in", omega_"out") + m log(omega_"in" / omega_"out") \
    &quad dot [1 / (2 m) sum_(i, j = 1)^n (A_(i j) - (omega_"in" - omega_"out") / (log omega_"in" - log omega_"out") (d_i d_j) / (2 m)) delta(g_i, g_j)] \
    &= m log(omega_"in" / omega_"out") Q_(gamma^*) (g) + C(A, omega_"in", omega_"out"), $
  where $gamma^* = (omega_"in" - omega_"out") / (log omega_"in" - log omega_"out")$. Finally, for any $a > b > 0$, let $t = sqrt(a / b) > 1$; then $log(a / b) = 2 log t$, while the standard inequalities $2 (t - 1) / (t + 1) < log t < (t^2 - 1) / (2 t)$ rearrange identically to $sqrt(a b) < (a - b) / log(a / b) < (a + b) / 2$, establishing the exact logarithmic-mean bounds on $gamma^*$.
]

#corollary(title: "Generative Origin of the Fortunato–Barthélemy Resolution Limit", origin: "Original Proof")[
  Let $G$ contain two disjoint cliques $K_r, K_r'$ of size $r >= 3$ connected by a single bridge edge, inside a larger connected graph of $m$ total edges @fortunato2007resolution. Under resolution-$gamma$ modularity, merging $K_r$ and $K_r'$ strictly increases $Q_gamma$ whenever
  $ l_"in" := binom(r, 2) < sqrt((m) / (2 gamma)) - 1 / 2 quad ("equivalently, " d_"tot" (K_r) < sqrt(2 m / gamma)). $
  Via Theorem 12(2), fixing $gamma = 1$ (standard Louvain/Leiden) implicitly imposes the generative constraint $omega_"in" - omega_"out" = log(omega_"in" / omega_"out")$. In a growing graph with $k = Theta(m)$ constant-size communities, the true within-community affinity scales as $omega_"in" = Omega(m)$ while $omega_"out" = O(1)$, forcing the true Bayes-optimal resolution to diverge as $gamma^* = Omega(m / log m) -> infinity$; thus $gamma = 1$ fails statistically by underestimating $gamma^*$!
]

#proof[
  Denote the volume of each clique by $D_1 = D_2 = r (r - 1) + 1 = 2 l_"in" + 1$ (counting the single bridge edge incident to each clique). Separating $K_r$ and $K_r'$ contributes $2 [(2 l_"in") / (2 m) - gamma (D_1 / (2 m))^2]$ to $Q_gamma$, whereas merging them into $K_r union K_r'$ adds the bridge weight $2$ to the internal volume and changes the null penalty from $2 D_1^2$ to $(2 D_1)^2 = 4 D_1^2$. Hence the exact modularity gain of merging is
  $ Delta Q_gamma = 2 / (2 m) - gamma (2 D_1^2) / (4 m^2) = 1 / m (1 - gamma (D_1^2) / (2 m)). $
  Thus $Delta Q_gamma > 0$ if and only if $D_1 < sqrt(2 m / gamma)$, i.e., $l_"in" < sqrt(m / (2 gamma)) - 1 / 2$ (at $gamma = 1$, every clique of fewer than $sqrt(m / 2)$ internal edges is erroneously merged with its neighbor @fortunato2007resolution). By Theorem 12(2), in a planted DC-SBM with $k$ balanced cliques of constant volume $D_1 = O(1)$ and total volume $2 m = k D_1$, the expected internal edge count per block is $omega_"in" D_1^2 / (2 m) = Theta(1)$, requiring $omega_"in" = Theta(m)$, while sparse cross-edges give $omega_"out" = O(1)$; substituting into $gamma^* = (omega_"in" - omega_"out") / log(omega_"in" / omega_"out")$ yields $gamma^* = Theta(m / log m) >> 1$, pinpointing the exact generative source of the resolution limit.
]

#remark(title: "Louvain vs. Leiden on the Identical Modularity Objective")[
  While Theorem 12 proves that Louvain (Alg. 107), Leiden (Alg. 108), CNM (Alg. 115), and Potts Spinglass (Alg. 117) all target the exact same objective $Q_gamma (g)$, Traag et al. @traag2019leiden demonstrated that Louvain's greedy node-moving phase can orphan internal bridge vertices and yield *internally disconnected* communities; Leiden inserts a refinement phase guaranteeing that every returned community is internally $gamma$-connected while monotonically improving $Q_gamma$, illustrating why identical variational objectives can still exhibit dramatically different discrete-optimization geometry (see also @brandes2008modularity on NP-hardness of exact $Q_gamma$ maximization).
]


#sec("VI", "Optimal Transport, Coreset Streaming, Subspace Recovery, and Deep Variational Clustering")

The final four theorems extend our reduction calculus to modern large-scale, high-dimensional, and neural clustering paradigms in Parts III, VI, and VIII: entropic optimal transport in self-supervised deep clustering, $(1 plus.minus epsilon)$ coreset compression for data streams, $ell_1$ self-expressive subspace recovery, and amortized variational deep embedding architectures.

#subsec("A", "SwAV Sinkhorn–Knopp as Entropic Optimal Transport Balanced K-Means")

Self-supervised visual clustering methods such as DeepCluster (Alg. 191, @caron2018deepcluster) face a chronic failure mode when jointly updating an encoder and cluster prototypes: trivial collapse, in which all embeddings map onto a single prototype. *SwAV* (*Swapping Assignments between Views*, Alg. 190, @caron2020swav) prevents collapse online by solving an equipartitioned entropic optimal transport problem @cuturi2013sinkhorn on every mini-batch.

Let $Z = [z_1, dots, z_B] in bb(R)^(d times B)$ be a mini-batch of $B$ unit-normalized neural embeddings ($norm(z_i)_2 = 1$), and let $C = [c_1, dots, c_k] in bb(R)^(d times k)$ be $k$ unit-normalized cluster prototypes ($norm(c_j)_2 = 1$). Define the spherical Euclidean cost matrix $M in bb(R)^(k times B)$ with entries $M_(j i) := 1 / 2 norm(c_j - z_i)_2^2 = 1 - chevron.l c_j, z_i chevron.r in [0, 2]$, uniform marginals $a := k^(-1) bold(1)_k in Delta_k$ and $b := B^(-1) bold(1)_B in Delta_B$, and the *balanced transportation polytope*
$ cal(U)(a, b) := {Q in bb(R)_(>= 0)^(k times B) : Q bold(1)_B = a, space Q^T bold(1)_k = b}. $
For regularization parameter $epsilon > 0$, SwAV computes the online code matrix $Q^* in cal(U)(a, b)$ solving
$ max_(Q in cal(U)(a, b)) "Tr"(Q^T C^T Z) + epsilon H(Q), quad H(Q) := - sum_(j=1)^k sum_(i=1)^B Q_(j i) log Q_(j i). $

#theorem(title: "SwAV Sinkhorn–Knopp <=> Entropic Optimal Transport Balanced K-Means", origin: "Classical — Cited")[
  Let $Z in bb(R)^(d times B)$, $C in bb(R)^(d times k)$ have unit-norm columns, $M_(j i) = 1 - chevron.l c_j, z_i chevron.r$, and $epsilon > 0$ @cuturi2013sinkhorn, @caron2020swav.
  1. *Exact Entropic OT & Balanced Spherical K-Means Equivalence:* On $cal(U)(a, b)$,
     $ &max_(Q in cal(U)(a, b)) ["Tr"(Q^T C^T Z) + epsilon H(Q)] \
       &quad = 1 - min_(Q in cal(U)(a, b)) ["Tr"(Q^T M) - epsilon H(Q)]. $
     Thus SwAV code computation is *identically* Cuturi's Entropic Optimal Transport problem $W_epsilon (a, b)$ @cuturi2013sinkhorn. In the zero-entropy limit $epsilon -> 0^+$ (when $k$ divides $B$), every extreme point of $cal(U)(a, b)$ scaled by $B$ is an integral $0$-$1$ assignment with *exact* cluster cardinalities $|C_j| = B / k$, so SwAV's code step converges to exact *cardinality-balanced Spherical $k$-Means* (Alg. 12).
  2. *Fenchel–Lagrange Dual & Sinkhorn–Knopp Matrix Scaling:* The primal problem has a unique strictly positive optimizer of diagonal-scaling form
     $ Q^*_(j i) = u_j^* K_(j i) v_i^*, quad K_(j i) := exp(- M_(j i) / epsilon) = e^(-1 / epsilon) exp(chevron.l c_j, z_i chevron.r / epsilon), $
     where $u_j^* = exp(alpha_j^* / epsilon) > 0$ and $v_i^* = exp(beta_i^* / epsilon) > 0$ maximize the smooth concave unconstrained dual
     $ g(alpha, beta) = alpha^T a + beta^T b - epsilon sum_(j=1)^k sum_(i=1)^B exp((alpha_j + beta_i - M_(j i)) / epsilon) + epsilon. $
     Exact alternating block coordinate ascent on $(alpha, beta)$ is *identically* the Sinkhorn–Knopp matrix-scaling iteration:
     $ u^((t + 1)) = a / (K v^((t))), quad v^((t + 1)) = b / (K^T u^((t + 1))) quad ("elementwise"). $
  3. *Global Linear Contraction in the Hilbert Projective Metric:* Defining $d_H (u, u') := log(max_j (u_j / u'_j) / min_j (u_j / u'_j))$ on $bb(R)_(> 0)^k$, Sinkhorn–Knopp converges globally from *any* $v^((0)) > 0$ at geometric rate
     $ d_H (u^((t + 1)), u^*) <= tanh^2 (Delta(K) / 4) d_H (u^((t)), u^*), $
     where the Birkhoff projective diameter satisfies the explicit sphere bound $Delta(K) := log max_(j, j', i, i') ((K_(j i) K_(j' i')) / (K_(j i') K_(j' i))) <= 4 / epsilon$.
]

#proof[
  *Part 1 (Primal Reformulation & Integrality as $epsilon -> 0^+$):* Because every $Q in cal(U)(a, b)$ satisfies $sum_(j, i) Q_(j i) = bold(1)_k^T Q bold(1)_B = bold(1)_k^T a = 1$, we have
  $ "Tr"(Q^T M) = sum_(j=1)^k sum_(i=1)^B Q_(j i) (1 - chevron.l c_j, z_i chevron.r) = 1 - "Tr"(Q^T C^T Z). $
  Negating and adding $epsilon H(Q)$ yields the exact identity of Part 1. As $epsilon -> 0^+$, $min_(Q in cal(U)(a, b)) "Tr"(Q^T M)$ is a standard bipartite transportation linear program with marginal vectors $B a = (B / k) bold(1)_k in bb(Z)^k$ and $B b = bold(1)_B in bb(Z)^B$ (when $k | B$). Since the node-edge incidence matrix of a bipartite graph is *totally unimodular*, every vertex of the rescaled polytope $B cal(U)(a, b)$ is integral in ${0, 1}^(k times B)$ with every column summing to $1$ (hard assignment) and every row summing to $B / k$ (exact equipartition), minimizing $sum_(j=1)^k sum_(i in C_j) (1 - chevron.l c_j, z_i chevron.r)$—precisely balanced Spherical $k$-Means!

  *Part 2 (Strong Fenchel–Lagrange Duality):* Because $x |-> x log x$ is strictly convex on $bb(R)_(>= 0)$ with vertical slope $-infinity$ at $0^+$ and $cal(U)(a, b)$ is compact with strictly positive interior point $Q_0 = a b^T > 0$, Slater's condition holds and the primal minimizer $Q^*$ is unique and strictly positive. Form the Lagrangian with dual potentials $alpha in bb(R)^k$, $beta in bb(R)^B$ (adding the constant $epsilon (1 - sum Q_(j i)) = 0$ to simplify derivatives):
  $ cal(L)(Q, alpha, beta) &= sum_(j, i) [M_(j i) Q_(j i) + epsilon Q_(j i) (log Q_(j i) - 1)] + epsilon \
    &quad - alpha^T (Q bold(1)_B - a) - beta^T (Q^T bold(1)_k - b). $
  Setting $partial cal(L) / partial Q_(j i) = M_(j i) + epsilon log Q_(j i) - alpha_j - beta_i = 0$ uniquely solves the inner minimization over $Q > 0$:
  $ Q_(j i) (alpha, beta) = exp((alpha_j + beta_i - M_(j i)) / epsilon) = e^(alpha_j / epsilon) K_(j i) e^(beta_i / epsilon). $
  Substituting $Q(alpha, beta)$ back into $cal(L)$ cancels $(M_(j i) - alpha_j - beta_i + epsilon log Q_(j i)) Q_(j i) = 0$ identically, yielding the exact smooth concave dual function
  $ g(alpha, beta) = alpha^T a + beta^T b - epsilon sum_(j, i) exp((alpha_j + beta_i - M_(j i)) / epsilon) + epsilon. $
  Now maximize $g(alpha, beta)$ over $alpha$ for fixed $beta$: setting $partial g / partial alpha_j = a_j - e^(alpha_j / epsilon) sum_(i=1)^B K_(j i) e^(beta_i / epsilon) = 0$ yields
  $ u_j := e^(alpha_j / epsilon) = (a_j) / ((K v)_j), quad v_i := e^(beta_i / epsilon). $
  Symmetrically, setting $partial g / partial beta_i = b_i - e^(beta_i / epsilon) sum_(j=1)^k K_(j i) u_j = 0$ yields $v_i = b_i / (K^T u)_i$. Thus alternating exact dual block maximization is *identically* the Sinkhorn–Knopp scaling iteration! (Notice that the constant factor $e^(-1/epsilon)$ in $K$ scales out identically between $u$ and $v$, so implementations work directly with $tilde(K)_(j i) = exp(chevron.l c_j, z_i chevron.r / epsilon)$.)

  *Part 3 (Birkhoff–Hopf Contraction):* Because $K_(j i) > 0$ strictly everywhere, the positive linear operator $x |-> K x$ maps the positive orthant into its interior. By the Birkhoff–Hopf theorem, any strictly positive matrix $K$ is a strict contraction in the Hilbert projective metric $d_H$ with Lipschitz constant $tanh(Delta(K) / 4) < 1$. Because inversion $x |-> a / x$ is an exact $d_H$-isometry ($d_H (a / x, a / y) = d_H (x, y)$), one full Sinkhorn cycle $u^((t)) |-> v^((t+1)) = b / (K^T u^((t))) |-> u^((t+1)) = a / (K v^((t+1)))$ composes two projective contractions ($K^T$ and $K$) and two isometries, contracting $d_H$ by at least $tanh^2 (Delta(K) / 4)$ per iteration. Finally, since $M_(j i) in [0, 2]$ on $bb(S)^(d - 1)$,
  $ (K_(j i) K_(j' i')) / (K_(j i') K_(j' i)) = exp((M_(j i') + M_(j' i) - M_(j i) - M_(j' i')) / epsilon) <= exp(4 / epsilon), $
  giving $Delta(K) <= 4 / epsilon$ and guaranteeing rapid geometric convergence (typically $3$–$5$ iterations in SwAV @caron2020swav).
]

#remark(title: "Why SwAV Eliminates Collapse in Amortized Deep Clustering")[
  In SwAV's outer loop, two augmented views $(z_t, z_s)$ yield Sinkhorn balanced codes $(q_t, q_s)$ (where $q_i := B Q^*_(:, i) in Delta_k$), and the encoder/prototypes minimize the swapped cross-entropy $ell(z_t, q_s) + ell(z_s, q_t)$ with softmax logits $C^T z / tau$. If all prototypes collapsed ($c_1 = dots = c_k$), then $K$ would be rank-1, forcing $q_i = (1/k) bold(1)_k$ uniformly across views—a saddle point repelled by cross-view semantic variance—while the row marginal $Q bold(1)_B = k^(-1) bold(1)_k$ forbids assigning the entire batch to a single prototype!
]

#subsec("B", "BICO and StreamKM++ Coreset Distortion Preservation for K-Means")

When $n$ reaches tens of millions in single-pass streaming regimes, running Lloyd's $k$-means (Alg. 1) or $k$-means++ (Alg. 2, @arthur2007kmeanspp) directly on $cal(X)$ exceeds memory budgets. Both *BICO* (Alg. 78, @fichtenberger2013bico—combining BIRCH's CF-tree, Alg. 64 @zhang1996birch, with geometric coreset theory @harpeled2004coresets) and *StreamKM++* (Alg. 195, @ackermann2012streamkmpp) solve this by constructing a compact weighted point set $(S, w)$ ($|S| << n$, $w: S -> bb(R)_(> 0)$) that preserves the $k$-means objective uniformly across *all* possible center sets $Q subset.eq bb(R)^d$ of cardinality $|Q| = k$.

Denote the unweighted and weighted $k$-means costs for any $Q subset.eq bb(R)^d$ ($|Q| = k$) by
$ &phi(cal(X), Q) := sum_(x in cal(X)) min_(q in Q) norm(x - q)_2^2, \
  &phi_w (S, Q) := sum_(p in S) w(p) min_(q in Q) norm(p - q)_2^2, $
and let $"OPT"_k (cal(X)) := min_(|Q| = k) phi(cal(X), Q)$. A weighted set $(S, w)$ is a *$(k, epsilon)$-strong coreset* for $cal(X)$ if, for *every* $Q subset.eq bb(R)^d$ with $|Q| = k$,
$ (1 - epsilon) phi(cal(X), Q) <= phi_w (S, Q) <= (1 + epsilon) phi(cal(X), Q). $

#theorem(title: "BICO & StreamKM++ Coreset Distortion Preservation for K-Means", origin: "Original Proof")[
  Let $cal(X) subset.eq bb(R)^d$ be a finite dataset with $"OPT"_k (cal(X)) > 0$ and fix $epsilon in (0, 1)$.
  1. *Approximation Transfer Reduction:* If $(S, w)$ is a $(k, epsilon)$-strong coreset of $cal(X)$ and $hat(Q)$ is any $c$-approximate $k$-means solution on $(S, w)$ (e.g., $c = 8(log k + 2)$ in expectation via weighted $k$-means++ @arthur2007kmeanspp, or $c = 9 + delta$ via local search @kanungo2004local), then
     $ phi(cal(X), hat(Q)) <= ((1 + epsilon) / (1 - epsilon)) c dot "OPT"_k (cal(X)). $
  2. *Merge-and-Reduce Streaming Composition (StreamKM++, Alg. 195):* Strong coresets satisfy two exact algebraic closure laws @harpeled2004coresets, @ackermann2012streamkmpp:
     (a) *Union Closure:* If $(S_1, w_1)$ and $(S_2, w_2)$ are $(k, epsilon)$-coresets of disjoint sets $cal(X)_1, cal(X)_2$, then $(S_1 union S_2, w_1 + w_2)$ is a $(k, epsilon)$-coreset of $cal(X)_1 union cal(X)_2$.
     (b) *Composition Closure:* If $(S_2, w_2)$ is a $(k, epsilon_2)$-coreset of $(S_1, w_1)$ and $(S_1, w_1)$ is a $(k, epsilon_1)$-coreset of $cal(X)$, then $(S_2, w_2)$ is a $(k, epsilon_1 + epsilon_2 + epsilon_1 epsilon_2)$-coreset of $cal(X)$.
     Consequently, a binary streaming merge-and-reduce tree of height $h = ceil(log_2 (n / m))$ invoking a $(k, epsilon')$-coreset reducer with $epsilon' = epsilon / (2 h)$ produces a guaranteed $(k, epsilon)$-strong coreset of the entire stream using only $O(m log(n / m))$ memory.
  3. *Exact BICO / CF-Tree Centroid Compression Bound (Alg. 78):* Let ${P_v}_(v in cal(V))$ be any partition of $cal(X)$ into compressed micro-clusters summarized by their exact Clustering Feature centroids $S = {mu_v : v in cal(V)}$ with weights $w(mu_v) = n_v := |P_v|$ and total internal compression inertia
     $ E_"tot" := sum_(v in cal(V)) E_v, quad E_v := sum_(x in P_v) norm(x - mu_v)_2^2. $
     Then for *every* center set $Q subset.eq bb(R)^d$ ($|Q| = k$) and *every* scalar $alpha > 0$,
     $ phi(cal(X), Q) - E_"tot" <= phi_w (S, Q) <= (1 + alpha) phi(cal(X), Q) + (1 + 1 / alpha) E_"tot". $
     In particular, whenever the BICO threshold-doubling tree (Alg. 78, @fichtenberger2013bico) or the StreamKM++ $D^2$-splitting coreset tree (Alg. 195, @ackermann2012streamkmpp) ensures $E_"tot" <= (epsilon^2 / 8) "OPT"_k (cal(X))$, the weighted centroid set $(S, w)$ is *deterministically* a $(k, epsilon)$-strong coreset of $cal(X)$!
]

#proof[
  *Part 1 (Approximation Transfer):* Let $Q^* = arg min_(|Q| = k) phi(cal(X), Q)$ achieve $"OPT"_k (cal(X))$, and let $Q_S^* = arg min_(|Q| = k) phi_w (S, Q)$. Applying the coreset lower bound at $hat(Q)$, the $c$-approximation property on $(S, w)$, optimality of $Q_S^*$, and the coreset upper bound at $Q^*$:
  $ phi(cal(X), hat(Q)) &<= (1 / (1 - epsilon)) phi_w (S, hat(Q)) <= (c / (1 - epsilon)) phi_w (S, Q_S^*) \
    &<= (c / (1 - epsilon)) phi_w (S, Q^*) <= ((1 + epsilon) / (1 - epsilon)) c dot phi(cal(X), Q^*). $

  *Part 2 (Merge-and-Reduce Closure):* For (a), by additivity of $phi$ over disjoint unions,
  $ &(1 - epsilon) [phi(cal(X)_1, Q) + phi(cal(X)_2, Q)] <= phi_(w_1) (S_1, Q) + phi_(w_2) (S_2, Q) \
    &quad <= (1 + epsilon) [phi(cal(X)_1, Q) + phi(cal(X)_2, Q)]. $
  For (b), multiplying the two-sided inequalities yields
  $ &(1 - epsilon_1)(1 - epsilon_2) phi(cal(X), Q) <= phi_(w_2) (S_2, Q) \
    &quad <= (1 + epsilon_1)(1 + epsilon_2) phi(cal(X), Q), $
  and $(1 + epsilon_1)(1 + epsilon_2) = 1 + epsilon_1 + epsilon_2 + epsilon_1 epsilon_2$, while $(1 - epsilon_1)(1 - epsilon_2) = 1 - (epsilon_1 + epsilon_2 - epsilon_1 epsilon_2) >= 1 - (epsilon_1 + epsilon_2 + epsilon_1 epsilon_2)$. In a binary merge tree of height $h$, every input block undergoes at most $h$ reductions (unions do not degrade accuracy), so with $epsilon' = epsilon / (2 h)$ and $epsilon in (0, 1]$:
  $ &(1 + epsilon / (2 h))^h <= exp(epsilon / 2) <= 1 + epsilon, \
    &(1 - epsilon / (2 h))^h >= 1 - h (epsilon / (2 h)) = 1 - epsilon / 2 >= 1 - epsilon, $
  where we used $1 + u <= e^u$, $e^u <= 1 + 2 u$ for $u in [0, 1 / 2]$, and Bernoulli's inequality $(1 - u)^h >= 1 - h u$!

  *Part 3 (Exact CF-Centroid Perturbation Sandwich):* Fix any center set $Q$ of size $k$ and any micro-cluster $P_v$ with centroid $mu_v$, size $n_v$, and internal inertia $E_v = sum_(x in P_v) norm(x - mu_v)_2^2$. Let $q_v := arg min_(q in Q) norm(mu_v - q)_2^2$ and $q(x) := arg min_(q in Q) norm(x - q)_2^2$ for each $x in P_v$.

  *Lower bound on $phi_w (S, Q)$:* Because $q(x)$ minimizes squared distance from $x$ to $Q$, $norm(x - q(x))_2^2 <= norm(x - q_v)_2^2$. Summing over $x in P_v$ and invoking the exact Parallel-Axis identity around $mu_v$:
  $ phi(P_v, Q) &= sum_(x in P_v) norm(x - q(x))_2^2 <= sum_(x in P_v) norm(x - q_v)_2^2 \
    &= E_v + n_v norm(mu_v - q_v)_2^2 = E_v + phi_w ({mu_v}, Q). $
  Summing over all $v in cal(V)$ gives $phi(cal(X), Q) - E_"tot" <= phi_w (S, Q)$.

  *Upper bound on $phi_w (S, Q)$:* For every $x in P_v$, optimality of $q_v$ gives $norm(mu_v - q_v)_2^2 <= norm(mu_v - q(x))_2^2$. Apply Young's inequality $norm(a + b)_2^2 <= (1 + alpha) norm(a)_2^2 + (1 + 1 / alpha) norm(b)_2^2$ (valid for all $alpha > 0$ since $alpha norm(a)_2^2 + alpha^(-1) norm(b)_2^2 >= 2 norm(a)_2 norm(b)_2 >= 2 chevron.l a, b chevron.r$) with $a = x - q(x)$ and $b = mu_v - x$:
  $ norm(mu_v - q_v)_2^2 &<= norm((x - q(x)) + (mu_v - x))_2^2 \
    &<= (1 + alpha) norm(x - q(x))_2^2 + (1 + 1 / alpha) norm(x - mu_v)_2^2. $
  Summing over all $x in P_v$ yields
  $ phi_w ({mu_v}, Q) = n_v norm(mu_v - q_v)_2^2 <= (1 + alpha) phi(P_v, Q) + (1 + 1 / alpha) E_v. $
  Summing over all $v in cal(V)$ proves the exact two-sided inequality of Part 3!

  Finally, suppose $E_"tot" <= (epsilon^2 / 8) "OPT"_k (cal(X))$. Since $"OPT"_k (cal(X)) <= phi(cal(X), Q)$ for *every* $Q$ of size $k$, setting $alpha := epsilon / 2$ yields:
  - Upper bound: $(1 + epsilon / 2) phi(cal(X), Q) + (1 + 2 / epsilon)(epsilon^2 / 8) "OPT"_k (cal(X)) <= [1 + epsilon / 2 + epsilon / 4 + epsilon^2 / 8] phi(cal(X), Q) <= (1 + epsilon) phi(cal(X), Q)$ for all $epsilon in (0, 1]$ (since $epsilon / 2 + epsilon / 4 + epsilon^2 / 8 <= 7 epsilon / 8 < epsilon$).
  - Lower bound: $phi_w (S, Q) >= phi(cal(X), Q) - (epsilon^2 / 8) phi(cal(X), Q) >= (1 - epsilon) phi(cal(X), Q)$.
  Thus $(S, w)$ is strictly a $(k, epsilon)$-strong coreset of $cal(X)$! In BICO @fichtenberger2013bico, storing exact additive CF triples $(n_v, sum_(x in P_v) x, sum_(x in P_v) norm(x)_2^2)$ inside a multi-level projection tree with adaptive radius doubling directly monitors and bounds $E_"tot"$ in $O(1)$ per insertion; in StreamKM++ @ackermann2012streamkmpp, building an $m$-leaf coreset tree by recursively splitting the leaf sampled with probability $E_v / E_"tot"$ around a $D^2$-sampled pivot drives $E_"tot" <= (epsilon^2 / 8) "OPT"_k (cal(X))$ with probability $>= 1 - delta$ using $m = O(k d epsilon^(-2) log(1 / delta))$ leaves.
]

#subsec("C", "Sparse Subspace Clustering to Subspace-Preserving Spectral Cut")

In high-dimensional Euclidean space $bb(R)^D$, data often reside not near isotropic centroids but along a union of $k$ low-dimensional linear subspaces $union_(ell=1)^k cal(S)_ell$ ($dim(cal(S)_ell) = d_ell < D$). *Sparse Subspace Clustering (SSC)* (Alg. 150, @elhamifar2013ssc) bridges Part VI (Subspace Clustering) directly to Part V (Spectral Clustering NJW, Alg. 103) via an $ell_1$ self-expressive basis-pursuit program.

Let $X = [x_1, dots, x_n] in bb(R)^(D times n)$ have unit-norm columns ($norm(x_i)_2 = 1$) partitioned into true subspace blocks $X = [X^((1)), dots, X^((k))]$ (up to unknown column permutation). For each datum $x_i$, SSC solves the convex linear program
$ c_i^* in arg min_(c in bb(R)^n) norm(c)_1 quad "subject to" quad x_i = X c, space c_i = 0, $
assembles the coefficient matrix $C^* = [c_1^*, dots, c_n^*] in bb(R)^(n times n)$, forms the symmetric non-negative affinity graph $W := |C^*| + |C^*|^T$, and executes Spectral Clustering NJW (Alg. 103, @ng2002spectral) on $W$. A matrix $C^*$ is *subspace-preserving* if $C^*_(j i) = 0$ whenever $x_i$ and $x_j$ belong to different subspaces.

#theorem(title: "Sparse Subspace Clustering -> Exact Subspace-Preserving Spectral Cut", origin: "Classical — Cited")[
  Let $X = [x_1, dots, x_n] in bb(R)^(D times n)$ ($norm(x_i)_2 = 1$) lie on a union of $k$ subspaces $cal(S)_1, dots, cal(S)_k$ such that every $x_i in cal(S)_ell$ lies in $"span"(X_(-i)^((ell)))$ @elhamifar2013ssc, @soltanolkotabi2012geometric.
  1. *Independent Subspaces (Elhamifar–Vidal):* If $cal(S)_1, dots, cal(S)_k$ are linearly independent—i.e., $dim(cal(S)_1 + dots + cal(S)_k) = sum_(ell=1)^k dim(cal(S)_ell)$—then *every* optimal solution $c_i^*$ of the SSC $ell_1$ program is strictly subspace-preserving.
  2. *Disjoint / Intersecting Subspaces via Geometric Dual Certificate (Soltanolkotabi–Candès):* Without assuming linear independence, let $X_(-i)^((ell))$ denote the columns in $cal(S)_ell$ excluding $x_i$, and let $cal(P)_(-i)^((ell)) := "conv"(plus.minus X_(-i)^((ell)))$ be their symmetrized convex hull with Euclidean inradius $r(cal(P)_(-i)^((ell)))$ inside $cal(S)_ell$. If there exists an intra-subspace dual optimal vector $nu_i in cal(S)_ell$ satisfying $(X_(-i)^((ell)))^T nu_i in [-1, 1]^(n_ell - 1)$, $chevron.l nu_i, x_i chevron.r = min {norm(c)_1 : X_(-i)^((ell)) c = x_i}$, and the strict dual separation bound
     $ max_(x_j in.not cal(S)_ell) |chevron.l nu_i, x_j chevron.r| < 1, $
     then *every* optimal SSC solution $c_i^*$ is strictly subspace-preserving. In particular, the deterministic geometric condition
     $ max_(x_j in.not cal(S)_ell) norm(Pi_(cal(S)_ell) x_j)_2 < r(cal(P)_(-i)^((ell))) $
     is strictly sufficient for subspace preservation (where $Pi_(cal(S)_ell)$ is orthogonal projection onto $cal(S)_ell$).
  3. *Exact Zero-Cut Spectral Recovery (Alg. 150 $->$ Alg. 103):* Whenever $C^*$ is subspace-preserving and each intra-subspace affinity subgraph $G_ell$ induced by $W = |C^*| + |C^*|^T$ on $cal(S)_ell$ is connected, $W$ is block-diagonal with zero inter-subspace cut ($"cut"(C_ell, overline(C)_ell) = 0$), $"dim"("ker"(L_"sym" (W))) = k$, and NJW Spectral Clustering (Alg. 103) maps every point in $cal(S)_ell$ onto the exact orthogonal vertex $R^T e_ell in bb(S)^(k - 1)$, recovering the true union-of-subspaces partition with *zero* clustering error!
]

#proof[
  *Part 1 (Independent Subspaces):* Fix $x_i in cal(S)_ell$ and let $c_i^*$ be any optimal solution of the SSC program, decomposed by subspace blocks as $c_i^* = (c_i^((1)), dots, c_i^((k)))$. Feasibility $x_i = X c_i^*$ rearranges to
  $ x_i - X^((ell)) c_i^((ell)) = sum_(r != ell) X^((r)) c_i^((r)). $
  The left-hand side lies in $cal(S)_ell$, whereas the right-hand side lies in $plus.o_(r != ell) cal(S)_r$. By linear independence of $cal(S)_1, dots, cal(S)_k$, we have $cal(S)_ell inter (plus.o_(r != ell) cal(S)_r) = {bold(0)}$, forcing both sides to vanish identically: $x_i = X^((ell)) c_i^((ell))$ and $sum_(r != ell) X^((r)) c_i^((r)) = bold(0)$. Consequently, the truncated vector $tilde(c)_i$ defined by $tilde(c)_i^((ell)) := c_i^((ell))$ and $tilde(c)_i^((r)) := bold(0)$ for all $r != ell$ remains strictly feasible ($X tilde(c)_i = X^((ell)) c_i^((ell)) = x_i$, $tilde(c)_(i i) = 0$) while satisfying
  $ norm(tilde(c)_i)_1 = norm(c_i^((ell)))_1 <= norm(c_i^((ell)))_1 + sum_(r != ell) norm(c_i^((r)))_1 = norm(c_i^*)_1. $
  If any cross-subspace block $c_i^((r))$ ($r != ell$) were nonzero, the inequality would be *strict*, contradicting optimality of $c_i^*$. Hence $c_i^((r)) = bold(0)$ for all $r != ell$.

  *Part 2 (Dual Certificate & Inradius Bound):* Let $c_i^((ell)*)$ solve the *restricted* intra-subspace $ell_1$ problem $min norm(c)_1$ s.t. $X_(-i)^((ell)) c = x_i$ (feasible by hypothesis). Its Lagrange dual is
  $ max_(nu in cal(S)_ell) chevron.l nu, x_i chevron.r quad "subject to" quad norm((X_(-i)^((ell)))^T nu)_infinity <= 1. $
  By strong LP duality, any dual optimal $nu_i in cal(S)_ell$ achieves $chevron.l nu_i, x_i chevron.r = norm(c_i^((ell)*))_1$. Now let $c = (c^((ell)), c^((-ell)))$ be *any* feasible solution of the full SSC problem ($X_(-i)^((ell)) c^((ell)) + X^((-ell)) c^((-ell)) = x_i$) with $c^((-ell)) != bold(0)$. Taking the inner product of the feasibility equation with $nu_i$ and applying Hölder's inequality ($chevron.l a, b chevron.r <= norm(a)_infinity norm(b)_1$) separately to both blocks:
  $ norm(c_i^((ell)*))_1 &= chevron.l nu_i, x_i chevron.r = chevron.l (X_(-i)^((ell)))^T nu_i, c^((ell)) chevron.r + chevron.l (X^((-ell)))^T nu_i, c^((-ell)) chevron.r \
    &<= norm((X_(-i)^((ell)))^T nu_i)_infinity norm(c^((ell)))_1 + (max_(x_j in.not cal(S)_ell) |chevron.l nu_i, x_j chevron.r|) norm(c^((-ell)))_1 \
    &< 1 dot norm(c^((ell)))_1 + 1 dot norm(c^((-ell)))_1 = norm(c)_1, $
  where the strict inequality follows from $max_(x_j in.not cal(S)_ell) |chevron.l nu_i, x_j chevron.r| < 1$ and $norm(c^((-ell)))_1 > 0$! Because padding $c_i^((ell)*)$ with zeros is feasible for the full problem with objective $norm(c_i^((ell)*))_1 < norm(c)_1$, *no* solution with nonzero cross-subspace support can be optimal!

  To obtain the geometric inradius bound @soltanolkotabi2012geometric, observe that the polar body of $cal(P)_(-i)^((ell)) = "conv"(plus.minus X_(-i)^((ell)))$ inside $cal(S)_ell$ is precisely the dual feasible set $cal(K)_(-i)^((ell)) := {nu in cal(S)_ell : norm((X_(-i)^((ell)))^T nu)_infinity <= 1}$. By definition of inradius, the Euclidean ball $B_(cal(S)_ell) (bold(0), r(cal(P)_(-i)^((ell)))) subset.eq cal(P)_(-i)^((ell))$, so taking polars reverses inclusion and yields $cal(K)_(-i)^((ell)) subset.eq B_(cal(S)_ell) (bold(0), 1 / r(cal(P)_(-i)^((ell))))$, i.e., $norm(nu_i)_2 <= 1 / r(cal(P)_(-i)^((ell)))$ for *every* dual feasible $nu_i in cal(S)_ell$. Since $nu_i in cal(S)_ell$, Cauchy–Schwarz gives
  $ |chevron.l nu_i, x_j chevron.r| = |chevron.l nu_i, Pi_(cal(S)_ell) x_j chevron.r| <= norm(nu_i)_2 norm(Pi_(cal(S)_ell) x_j)_2 <= (norm(Pi_(cal(S)_ell) x_j)_2) / (r(cal(P)_(-i)^((ell)))) < 1, $
  verifying the dual certificate deterministically!

  *Part 3 (Exact NJW Spectral Recovery):* Because $C^*$ is subspace-preserving, $W = |C^*| + |C^*|^T$ decomposes (after block permutation) into $W = "diag"(W^((1)), dots, W^((k)))$ with zero off-block weights, hence $L_"sym" (W) = "diag"(L_"sym" (W^((1))), dots, L_"sym" (W^((k))))$. By hypothesis each block subgraph $G_ell$ is connected (with positive degrees since $x_i = X^((ell)) c_i^((ell))$ with $norm(x_i)_2 = 1$ forces $c_i^((ell)) != bold(0)$), so each block $L_"sym" (W^((ell)))$ has simple eigenvalue $0$ with unique unit eigenvector $w_ell = "vol"(C_ell)^(-1 / 2) D_ell^(1 / 2) bold(1)_(n_ell)$ @vonluxburg2007tutorial. Therefore $0$ has exact multiplicity $k$ in $L_"sym" (W)$ with $lambda_(k + 1) > 0$, and by Theorem 2(3), NJW row-normalization maps every point $i in C_ell$ identically to the unit vertex $R^T e_ell in bb(S)^(k - 1)$ with pairwise inter-cluster distance $sqrt(2)$ and zero intra-cluster variance, guaranteeing exact recovery!
]

#remark(title: "Connectivity Caveat in High-Dimensional Subspaces")[
  Notice the essential connectivity hypothesis in Theorem 15(3): while $ell_1$ minimization guarantees zero *false* cross-subspace edges, Nasihatkon and Hartley proved that when $d_ell >= 4$, sparse $ell_1$ graphs within a single subspace can split into disconnected components (over-segmentation), explaining why practical SSC variants enforce connectivity via Elastic-Net ($ell_1 + lambda ell_2$) or post-hoc subspace merging.
]

#subsec("D", "Deep Clustering (DEC/IDEC/DCN/VaDE) as Amortized Variational EM / Manifold K-Means")

Deep clustering architectures in Part VIII—*DEC* (Alg. 181, @xie2016dec), *IDEC* (Alg. 182, @guo2017idec), *DCN* (Alg. 183, @yang2017dcn), and *VaDE* (Alg. 184, @jiang2017vade)—jointly learn a nonlinear neural representation $z = f_phi (x)$ (or stochastic variational encoder $q_phi (z | x) = cal(N)(z | tilde(mu)_phi (x), "diag"(tilde(sigma)_phi^2 (x)))$ @kingma2014vae) alongside cluster parameters $Theta = {(pi_j, mu_j, Sigma_j)}_(j=1)^k$. Rather than four ad-hoc neural heuristics, Theorem 16 proves that all four algorithms are exact specializations of ONE *Master Amortized Latent-Cluster Free Energy*:
$ cal(J)(phi, psi, Theta, P) &:= lambda_"rec" sum_(i=1)^n bb(E)_(q_phi (z | x_i)) [ell_"rec" (x_i, g_psi (z))] \
  &quad + cal(D)_"clus" (q_phi, P, Theta), $
where $g_psi : cal(Z) -> cal(X)$ is a neural decoder, $lambda_"rec" >= 0$ weights manifold reconstruction fidelity, and $cal(D)_"clus"$ couples the latent distribution $q_phi$ to cluster parameters $Theta$ and auxiliary soft/hard assignments $P = (p_(i j))$.

#theorem(title: "Deep Clustering (Algs. 181–184) <=> Amortized Variational EM & Manifold K-Means", origin: "Original Unified Proof")[
  Under the Master Amortized Latent-Cluster Free Energy $cal(J)(phi, psi, Theta, P)$:
  1. *DCN (Alg. 183, @yang2017dcn) <=> Regularized Alternating Manifold K-Means:* Setting $q_phi (z | x_i) = delta_(f_phi (x_i)) (z)$, $ell_"rec" (x, hat(x)) = norm(x - hat(x))_2^2$, hard one-hot assignments $P = H in {0, 1}^(n times k)$, and $cal(D)_"clus" = (lambda / 2) sum_(i=1)^n norm(f_phi (x_i) - sum_(j=1)^k h_(i j) mu_j)_2^2$ ($lambda_"rec" = 1$, $lambda > 0$):
     (a) For fixed $(phi, psi)$, minimizing $cal(J)$ over $(H, bold(mu))$ executes *exact Lloyd's $k$-means* (Alg. 1) in latent space ${f_phi (x_i)}$;
     (b) As $lambda / lambda_"rec" -> 0^+$ around a pretrained autoencoder, DCN converges to $k$-means on fixed latent codes;
     (c) Conversely, if $lambda_"rec" = 0$ and the encoder class contains constant maps, the unconstrained infimum collapses trivially to $cal(J) = 0$ via $f_phi (x) equiv mu_1$ (trivial representation collapse).
  2. *VaDE (Alg. 184, @jiang2017vade) <=> Exact Amortized Variational EM for Latent GMM:* Let $p_theta (x, z, c) = pi_c cal(N)(z | mu_c, "diag"(sigma_c^2)) p_psi (x | z)$ and mean-field posterior $q(z, c | x_i) = q_phi (z | x_i) q_i (c)$. Choosing $lambda_"rec" = 1$, $ell_"rec" (x_i, g_psi (z)) = - log p_psi (x_i | z)$, and $cal(D)_"clus" = sum_(i=1)^n "KL"(q_phi (z | x_i) q_i (c) || p(z, c))$ gives $cal(J) = - sum_(i=1)^n "ELBO"(x_i)$. Furthermore:
     (a) Holding $q_phi (z | x_i)$ and $Theta$ fixed, the unique globally optimal variational cluster posterior is
     $ q_i^* (c = j) = (exp(bb(E)_(q_phi (z | x_i)) [log p(z, c = j)])) / (sum_(ell=1)^k exp(bb(E)_(q_phi (z | x_i)) [log p(z, c = ell)])), $
     where $p(z, c = j) = pi_j cal(N)(z | mu_j, "diag"(sigma_j^2))$.
     (b) In the deterministic encoder limit $tilde(sigma)_phi (x_i) -> 0^+$ ($q_phi (z | x_i) -> delta_(z_i)$), $q_i^* (c = j)$ reduces *identically* to the exact GMM E-step responsibility $gamma_(i j)$ (Alg. 81) on $z_i = tilde(mu)_phi (x_i)$, and maximizing the ELBO over $(pi_j, mu_j, sigma_j^2)$ yields the *exact diagonal-covariance GMM M-step* in latent space!
  3. *DEC (Alg. 181, @xie2016dec) & IDEC (Alg. 182, @guo2017idec) <=> Self-Sharpening Student-$t$ Variational EM:* Setting $q_phi = delta_(z_i)$, $cal(D)_"clus" = "KL"(P || Q)$ with Student-$t$ ($nu = 1$) soft assignments $q_(i j) = (1 + norm(z_i - mu_j)_2^2)^(-1) / sum_ell (1 + norm(z_i - mu_ell)_2^2)^(-1)$, and alternating target $p_(i j) = (q_(i j)^2 / f_j) / sum_ell (q_(i ell)^2 / f_ell)$ ($f_j = sum_i q_(i j)$) recovers DEC ($lambda_"rec" = 0$) and IDEC ($lambda_"rec" > 0$). Moreover:
     (a) Holding $P$ fixed, the exact centroid gradient is
     $ &nabla_(mu_j) "KL"(P || Q) \
       &quad = -2 sum_(i=1)^n (1 + norm(z_i - mu_j)_2^2)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j), $
     which performs Cauchy-reweighted robust M-step attraction toward high-confidence targets;
     (b) Every fixed point of the self-sharpening map $P(Q) = Q$ satisfies $q_(i j) in {0, c_i}$ with $c_i$ constant on its active support, meaning only *hard partitions* or *degenerate uniform-collapse equilibria* ($q_(i j) = 1 / k$) are fixed points;
     (c) Without $lambda_"rec" > 0$ (DEC), minimizing $"KL"(P || Q)$ over $f_phi$ can map arbitrary non-adjacent inputs onto identical prototypes, whereas IDEC's reconstruction penalty $lambda_"rec" > 0$ enforces approximate encoder injectivity on the data manifold.
]

#proof[
  *Part 1 (DCN & Collapse Theorem):* Substituting $q_phi (z | x_i) = delta_(f_phi (x_i)) (z)$ collapses the expectation to deterministic evaluation:
  $ cal(J)_"DCN" (phi, psi, H, bold(mu)) &= lambda_"rec" sum_(i=1)^n norm(x_i - g_psi (f_phi (x_i)))_2^2 \
    &quad + lambda / 2 sum_(i=1)^n sum_(j=1)^k h_(i j) norm(f_phi (x_i) - mu_j)_2^2. $
  When $(phi, psi)$ are frozen, the first sum is constant, so alternating minimization over $H$ (nearest latent centroid) and $bold(mu)$ (mean of assigned latent codes $z_i = f_phi (x_i)$) is *identically* Lloyd's $k$-means (Alg. 1). When $lambda_"rec" = 0$, setting $f_phi (x) equiv mu_1$ for all $x in cal(X)$ and $h_(i 1) = 1$ achieves $cal(J)_"DCN" = 0$, the absolute global lower bound of a sum of squares, destroying all data-discriminative structure. Conversely, when $lambda_"rec" > 0$, if $f_phi$ collapsed two distinct points $x_a != x_b$ to $f_phi (x_a) = f_phi (x_b) = z_0$, then $g_psi (z_0)$ can match at most their midpoint $(x_a + x_b) / 2$, incurring irreducible reconstruction penalty $>= (lambda_"rec" / 2) norm(x_a - x_b)_2^2 > 0$; thus $lambda_"rec" > 0$ acts as a topological homeomorphism regularizer!

  *Part 2 (VaDE Exact Variational EM Proof):* For each observation $x_i$, apply Jensen's inequality with mean-field variational distribution $q(z, c | x_i) = q_phi (z | x_i) q_i (c)$ @blei2006variational, @jiang2017vade:
  $ log p_theta (x_i) &>= bb(E)_(q_phi (z | x_i) q_i (c)) [log((p_psi (x_i | z) p(z, c)) / (q_phi (z | x_i) q_i (c)))] \
    &= bb(E)_(q_phi (z | x_i)) [log p_psi (x_i | z)] \
    &quad - "KL"(q_phi (z | x_i) q_i (c) || p(z, c)) =: "ELBO"(x_i). $
  Isolating the terms in $-"ELBO"(x_i)$ that depend on the discrete cluster posterior $q_i (c) in Delta_k$:
  $ &"KL"(q_phi (z | x_i) q_i (c) || p(z, c)) + H(q_phi (z | x_i)) \
    &quad = sum_(j=1)^k q_i (j) [log q_i (j) - bb(E)_(q_phi (z | x_i)) [log p(z, c = j)]]. $
  Define the normalized Gibbs distribution $tilde(q)_i (j) := Z_i^(-1) exp(bb(E)_(q_phi (z | x_i)) [log p(z, c = j)])$, where $Z_i := sum_(ell=1)^k exp(bb(E)_(q_phi (z | x_i)) [log p(z, c = ell)])$. Then the bracketed sum equals $"KL"(q_i || tilde(q)_i) - log Z_i$, which is *uniquely minimized* over $q_i in Delta_k$ by Gibbs' inequality at $q_i^* = tilde(q)_i$, proving 2(a)! (Notice that in the standard VaDE implementation @jiang2017vade, substituting $q_i (c) = p(c | z)$ inside the single-sample Monte Carlo expectation yields the exact same identity in the small-variance regime, while $q_i^* = tilde(q)_i$ is the *strictly tighter* coordinate-ascent variational E-step because by Jensen's inequality $bb(E)_(q_phi) [log p(z, c)] <= log bb(E)_(q_phi) [p(z, c)]$ and $log Z_i$ minimizes exact mean-field KL.)

  Now send the encoder variance $tilde(sigma)_phi (x_i) -> 0^+$, so $q_phi (z | x_i) -> delta_(z_i)$ with $z_i = tilde(mu)_phi (x_i)$. Evaluating the expectation at $z = z_i$ gives
  $ q_i^* (c = j) -> (pi_j cal(N)(z_i | mu_j, "diag"(sigma_j^2))) / (sum_(ell=1)^k pi_ell cal(N)(z_i | mu_ell, "diag"(sigma_ell^2))) = gamma_(i j) (z_i), $
  the *exact* GMM E-step responsibility in latent space! Similarly, differentiating $sum_(i=1)^n "ELBO"(x_i)$ with respect to the Gaussian prior parameters $(pi_j, mu_j, sigma_(j r)^2)$ while holding $(phi, psi, q^*)$ fixed gives
  $ &pi_j^* = 1 / n sum_(i=1)^n q_i^* (j), quad mu_j^* = (sum_(i=1)^n q_i^* (j) tilde(mu)_phi (x_i)) / (sum_(i=1)^n q_i^* (j)), \
    &(sigma_(j r)^2)^* = (sum_(i=1)^n q_i^* (j) [tilde(sigma)_(phi, r)^2 (x_i) + (tilde(mu)_(phi, r) (x_i) - mu_(j r)^*)^2]) / (sum_(i=1)^n q_i^* (j)). $
  As $tilde(sigma)_phi (x_i) -> 0^+$, the encoder-variance floor $tilde(sigma)_(phi, r)^2 (x_i)$ vanishes and these three equations become *identically* the exact diagonal-covariance GMM M-step (Alg. 81)!

  *Part 3 (DEC / IDEC Student-$t$ Gradient & Fixed-Point Analysis):* Write $a_(i j) := (1 + norm(z_i - mu_j)_2^2)^(-1)$ and $S_i := sum_(ell=1)^k a_(i ell)$, so $q_(i j) = a_(i j) / S_i$ and $log q_(i j) = -log(1 + norm(z_i - mu_j)_2^2) - log S_i$. Differentiating $L_"DEC" := sum_(i, j) p_(i j) log(p_(i j) / q_(i j)) = "const" - sum_(i, j) p_(i j) log q_(i j)$ with respect to $mu_j$ (for fixed $P$ with $sum_ell p_(i ell) = 1$):
  $ nabla_(mu_j) L_"DEC" &= - sum_(i=1)^n [p_(i j) nabla_(mu_j) log a_(i j) - (sum_(ell=1)^k p_(i ell)) nabla_(mu_j) log S_i] \
    &= - sum_(i=1)^n (p_(i j) - (a_(i j)) / S_i) (2 (z_i - mu_j)) / (1 + norm(z_i - mu_j)_2^2) \
    &= - 2 sum_(i=1)^n (1 + norm(z_i - mu_j)_2^2)^(-1) (p_(i j) - q_(i j)) (z_i - mu_j), $
  proving 3(a)! Observe that the prefactor $(1 + norm(z_i - mu_j)_2^2)^(-1)$ downweights far outliers exactly as in the Student-$t$ mixture M-step (Alg. 88, @peel2000robust), while $(p_(i j) - q_(i j))$ pulls $mu_j$ toward points whose sharpened target exceeds current probability.

  To prove 3(b), suppose $P(Q) = Q$ is a fixed point, i.e., for every $(i, j)$,
  $ q_(i j) = (q_(i j)^2 / f_j) / (sum_(ell=1)^k q_(i ell)^2 / f_ell) =: (q_(i j)^2 / f_j) / Z_i. $
  Thus for every $(i, j)$ with $q_(i j) > 0$, dividing both sides by $q_(i j)$ gives $q_(i j) = f_j Z_i$. Summing over $j in "supp"(q_i)$ yields $1 = Z_i sum_(j in "supp"(q_i)) f_j$, hence $q_(i j) = f_j / (sum_(ell in "supp"(q_i)) f_ell)$ for all $j in "supp"(q_i)$. In particular, if $q_(i j) > 0$ for all $(i, j)$ (as strictly enforced by the strictly positive Student-$t$ kernel $a_(i j) > 0$ for finite embeddings!), then $"supp"(q_i) = {1, dots, k}$ for every $i$, forcing $q_(i j) = f_j / n$ to be *strictly identical across all $i$* (complete collapse, with symmetric point $q_(i j) = 1 / k$)! Thus the *only* non-degenerate boundary attractors of DEC's self-sharpening dynamics are hard $0$-$1$ partitions ($norm(z_i - mu_(c_i))_2 -> 0$ or $norm(z_i - mu_ell)_2 -> infinity$), and without IDEC's reconstruction penalty $lambda_"rec" > 0$ @guo2017idec, nothing prevents the encoder $f_phi$ from distorting local neighborhood distances to force $q_(i j) -> bb(1)_(j = c_i)$.
]

#page(columns: 1)[
  #sec("VII", "Master Cross-Family Reduction Summary and Architectural Synthesis")

  #subsec("A", "Complete 16-Theorem Cross-Family Master Reduction Table")

  Table 3 consolidates all sixteen formal equivalence, asymptotic limit, functorial, and coreset reduction theorems proved in Sections II–VI into a single unified reference matrix spanning Parts I–VIII of the 200-algorithm monograph.

  #show figure: set block(breakable: true)
  #show table: set text(size: 7.4pt)
  #set table(inset: (x: 4pt, y: 2.8pt))
  #figure(
    table(
      columns: (0.38fr, 1.38fr, 0.72fr, 1.42fr, 1.15fr, 1.35fr),
      align: (center, left, left, left, left, left),
      table.header(
        [Thm], [Algorithms Connected (Catalog IDs)], [Reduction Class], [Exact Condition / Limit / Transform], [Rate / Error / Exact Offset], [Algorithmic & Engineering Implication],
      ),
      [1], [Weighted Kernel $k$-Means (Alg. 16) $arrow.l.r.double$ Ratio Cut (Alg. 105) & Ncut (Alg. 104)], [Exact Objective Isomorphism], [Rcut: $D_w = I_n$, $K = sigma I_n - L$ ($sigma >= lambda_"max" (L)$). Ncut: $D_w = D$, $K = sigma D^(-1) + D^(-1) W D^(-1)$ ($sigma >= 1$).], [Exact affine offset $(n - k)(sigma + 1) - "Tr"(L_"sym")$; zero gap.], [Eliminates $O(n^3)$ eigendecomposition; kernel Lloyd monotonically refines discrete Ncut @dhillon2004kernel.],
      [2], [Spectral NJW (Alg. 103) $arrow.l.r.double$ Commute-Time, Walktrap (Alg. 114), PIC (Alg. 132)], [Spectral Filter Unification], [Universal eigenbasis $v_r = D^(-1/2) u_r$ of $L_"rw"$ with filter $g(lambda_r)$: $"vol"(G)/lambda_r$ (CT), $(1 - lambda_r)^(2 t)$ (WT), $bb(1)_(r <= k)$ (NJW).], [Davis–Kahan $O(delta / lambda_(k+1))$; PIC tail decay $O(max_(r > k) |1 - lambda_r|^t)$.], [Low-pass filtering defeats commute-time degree degeneracy; PIC cuts $O(n^3)$ to $O(t |E|)$ @lin2010pic.],
      [3], [Isotropic GMM EM (Alg. 81) $->$ Lloyd's $k$-Means (Alg. 1)], [Zero-Temp Laplace Limit], [Fixed homoscedastic $sigma -> 0^+$ with $pi_"min" > 0$ and tie-free Voronoi margin $Delta(bold(mu)) > 0$.], [E/M-step error $O(exp(-Delta / (2 sigma^2)))$; free energy $O(sigma^2)$.], [Explains why $k$-means fails on anisotropic/unequal covariances; avoids free-$Sigma$ MLE singularity.],
      [4], [Fuzzy C-Means (Alg. 151) $->$ $k$-Means (Alg. 1); Entropy-$k$-Means $arrow.l.r.double$ Softmax EM (Alg. 81)], [Fuzzifier Limit & Isomorphism], [FCM: $m -> 1^+$ (hard $k$-means); $m > (1 - 2 lambda_"max" (S_"dir"))^(-1)$ (mean collapse). Entropy: $T = 2 sigma^2$.], [FCM rate $O(rho_i^(2 / (m - 1)))$; exact identity $min_U cal(F)_(2 sigma^2) equiv -2 sigma^2 ell + "const"$.], [Standard $m = 2$ FCM collapses in high $d$ ($m_"crit" approx d / (d - 2)$); Softmax EM trades scale-invariance for exponential tails.],
      [5], [vMF Mixture (Alg. 89) $->$ Spherical $k$-Means (Alg. 12)], [Zero-Temp Limit ($kappa -> infinity$)], [Unit hypersphere $bb(S)^(d - 1)$, shared concentration $kappa -> +infinity$, cosine margin $eta(bold(mu)) > 0$.], [Centroid error $O(exp(-kappa eta))$; rescaled $kappa^(-1) ell_"vMF"$ error $O((log kappa) / kappa)$.], [Justifies Spherical $k$-Means for high-$d$ text/deep embeddings without Bessel evaluation @banerjee2005vmf.],
      [6], [Exp.-Family Mixtures (Algs. 81, 89, 95–97) $->$ Bregman $k$-Means (Algs. 1, 7, 12, 92)], [Legendre Dual & $tau -> 0^+$], [Conjugate pair $phi = psi^*$, $D_phi (x || mu) = phi(x) - phi(mu) - chevron.l x - mu, nabla phi(mu) chevron.r$, dispersion $tau -> 0^+$.], [Exact density bijection; $O(exp(-Delta^phi / tau))$ hard limit.], [Arithmetic mean is universal M-step across all regular exponential families @banerjee2005clustering.],
      [7], [Single Linkage (Alg. 56) $arrow.l.r.double$ Kruskal MST $arrow.l.r.double$ Core-Point DBSCAN (Alg. 31)], [Graph-Theoretic Isomorphism], [Cut $T^*$ at height $epsilon$ on $cal(X)_"core" (epsilon, "MinPts")$; exact on all $cal(X)$ for $"MinPts" <= 2$ or DBSCAN\*.], [Zero error on core sets; subdominant ultrametric $u_"SL" <= delta$.], [One $O(n^2)$ MST computes *all* DBSCAN\* core partitions across every $epsilon$ simultaneously.],
      [8], [HDBSCAN (Alg. 32) $arrow.l.r.double$ Single Linkage (Alg. 56) on $d_"mreach"$], [Constructive Metric Functor], [Transform $delta |-> d_"mreach" (x, y) :=$ $max{d_"core" (x), d_"core" (y), delta(x, y)}$; apply SLINK + $O(n)$ EOM DP.], [Exact metric space & exact continuum DBSCAN\* filtration isomorphism.], [Unifies multi-scale density & hierarchical trees; achieves Hartigan cluster-tree consistency @campello2015hierarchical.],
      [9], [Mean Shift (Alg. 34) $arrow.l.r.double$ Gaussian KDE EM $arrow.l.r.double$ DENCLUE 2.0 (Alg. 43)], [Exact EM / MM Fixed-Point], [Gaussian KDE $hat(f)_h (x)$ treated as fixed mixture; convex profile $k(r)$ with shadow $g = -k'$.], [Exact ascent $Delta log hat(f)_h = (2 h^2)^(-1) norm(Delta x)_2^2 + "KL"(w^((t)) || w^((t+1)))$.], [Guarantees step-size-free monotone mode ascent, eliminating DENCLUE 1.0 learning-rate tuning @carreira2006fast.],
      [10], [Ward's Method (Alg. 59) $arrow.l.r.double$ Greedy Agglomerative $k$-Means WCSS (Alg. 1)], [Exact Merge Delta & Recursion], [Merge cost $Delta(A, B) = (n_A n_B / (n_A + n_B)) norm(mu_A - mu_B)_2^2$; Lance–Williams $alpha_A, alpha_B, beta$.], [Exact WCSS telescope $J(cal(C)_k^"Ward") = sum_(t=1)^(n - k) Delta_t$; zero inversions.], [Guarantees $O(n^2)$ Nearest-Neighbor Chain exactness and ideal structured initialization for Lloyd's $k$-means @mullner2011modern.],
      [11], [KL-NMF (Alg. 92) $arrow.l.r.double$ pLSA Maximum Likelihood (Alg. 84)], [Exact Objective Isomorphism], [Mass-conserving manifold $sum (W H)_(i j) = N$; canonical map $W H = N U S V^T$ up to diagonal $D$.], [Exact identity $D_"KL" (X || W H) = -cal(L)_"pLSA" + C(X)$.], [Identical stationary points, yet alternating Lee–Seung KL-NMF tightens MM surrogate faster than simultaneous pLSA EM @ding2008nmfplsa.],
      [12], [Modularity (Louvain/Leiden, Algs. 107–108) $arrow.l.r.double$ Potts Spinglass (Alg. 117) $arrow.l.r.double$ DC-SBM (Alg. 131)], [Exact Objective Isomorphism], [Uniformly planted Poisson DC-SBM ($omega_"in" > omega_"out"$, $theta_i = d_i / sqrt(2 m)$) with $gamma^* = (omega_"in" - omega_"out") / log(omega_"in" / omega_"out")$.], [Exact affine scaling $ell_"DC" = m log(omega_"in" / omega_"out") Q_(gamma^*) + C$.], [Derives Fortunato–Barthélemy resolution limit ($l_"in" < sqrt(m / (2 gamma))$) from fixed $gamma = 1$ misspecification.],
      [13], [SwAV Sinkhorn–Knopp (Alg. 190) $arrow.l.r.double$ Entropic OT Balanced $k$-Means (Alg. 12)], [Fenchel–Lagrange OT Dual], [Transportation polytope $cal(U)(k^(-1) bold(1)_k, B^(-1) bold(1)_B)$, spherical cost $M_(j i) = 1 - chevron.l c_j, z_i chevron.r$, $epsilon -> 0^+$.], [Hilbert projective contraction $tanh^2 (Delta(K) / 4)$, $Delta(K) <= 4 / epsilon$.], [3–5 GPU matrix multiplications enforce exact equipartition and eliminate trivial deep-clustering collapse @caron2020swav.],
      [14], [BICO (Alg. 78) & StreamKM++ (Alg. 195) $<=_(1 plus.minus epsilon)$ Euclidean $k$-Means (Alg. 1)], [Multiplicative Coreset Reduction], [CF centroid inertia $E_"tot" <= (epsilon^2 / 8) "OPT"_k (cal(X))$; binary merge tree height $h$ with $epsilon' = epsilon / (2 h)$.], [Exact $(1 plus.minus epsilon)$ uniform distortion; $((1 + epsilon) / (1 - epsilon)) c$ transfer.], [Compresses massive streams into $O(m log(n / m))$ memory with provable offline $k$-means++ guarantees @fichtenberger2013bico.],
      [15], [Sparse Subspace Clustering (Alg. 150) $->$ Subspace-Preserving NJW (Alg. 103)], [Dual-Certificate Functor], [$ell_1$ basis pursuit $x_i = X c_i$ ($c_(i i) = 0$) under independent subspaces or $max_(j in.not cal(S)_ell) norm(Pi_(cal(S)_ell) x_j)_2 < r(cal(P)_(-i)^((ell)))$.], [Zero inter-subspace cut; exact $k$-dim nullspace of $L_"sym" (W)$ if connected.], [Reduces nonlinear union-of-subspaces discovery directly to zero-cut spectral clustering @elhamifar2013ssc.],
      [16], [Deep Clustering: DEC/IDEC/DCN/VaDE (Algs. 181–184) $arrow.l.r.double$ Amortized Variational EM / Manifold $k$-Means], [Unified Amortized Free Energy], [Master functional $lambda_"rec" bb(E)[ell_"rec"] + cal(D)_"clus"$; VaDE $tilde(sigma)_phi -> 0^+$; DCN $q_phi = delta_(f_phi)$; DEC/IDEC Student-$t$ KL.], [Exact latent GMM E/M-step for VaDE; exact $cal(J) = 0$ collapse proof when $lambda_"rec" = 0$.], [Proves reconstruction loss ($lambda_"rec" > 0$) or OT balancing (Thm. 13) is mathematically mandatory to prevent manifold collapse.],
    ),
    caption: [Master Cross-Family Mathematical Reduction Table synthesizing Theorems 1–16 across all eight paradigms of the 200-algorithm monograph.],
  )
]

#subsec("B", "Transitive Reduction Chains and Practical Algorithmic Architecture")

Because the reduction relations of Definition 1 compose transitively—exact isomorphisms compose with asymptotic limits, metric functors compose with spectral filters, and coresets compose with any target partitioner—Theorems 1–16 collapse the 200-algorithm catalog into four *master transitive highways* that guide both theoretical research and production system design.

*1. The Variational–Thermodynamic Prototype Highway:*
$ &"VaDE (184)" ->^(tilde(sigma)_phi -> 0) "Latent GMM (81)" arrow.l.r.double "Entropy-"k"-Means" \
  &quad ->^(sigma, T -> 0) "DCN / Lloyd "k"-Means (183, 1)" <-^(m -> 1^+) "FCM (151)" \
  &quad <-_(1 plus.minus epsilon) "BICO / StreamKM++ (78, 195)" <-^"greedy" "Ward (59)". $
Every method on this chain implicitly minimizes an isotropic quadratic (or, via Theorem 6, Bregman) distortion functional in either raw feature space or an autoencoder latent manifold. Two concrete engineering laws follow immediately: (i) whenever clusters are convex and nearly spherical, streaming coresets (Theorem 14) and Ward seeding (Theorem 10) provide provable speed and initialization gains for Lloyd's $k$-means without changing the underlying variational objective; (ii) whenever a practitioner replaces hard $k$-means with FCM (Alg. 151) in dimension $d >= 5$, Theorem 4(2) warns that the default fuzzifier $m = 2$ exceeds $m_"crit" approx d / (d - 2)$ and triggers global-mean collapse, making temperature-controlled Softmax EM / GMM (Alg. 81) or calibrated $m in (1, 1 + 2 / d)$ strictly preferable.

*2. The Spectral–Kernel–Subspace Highway:*
$ &"SSC (150)" ->^(ell_1 " dual") "Block-"W" NJW (103)" \
  &quad arrow.l.r.double^("filter " g(lambda)) "Walktrap (114) / PIC (132)" <-^"Stiefel" "Ncut / Rcut (104–105)" \
  &quad arrow.l.r.double^("shift " sigma) "Weighted Kernel "k"-Means (16)". $
This chain proves that nonlinear manifold separation, random-walk diffusion, and sparse subspace recovery all converge onto the *same* bottom-$k$ eigenspace of $L_"sym"$. Practically, Theorem 1 and Theorem 2 unlock an $O(n^3) -> O(t |E|)$ computational shortcut: on million-node sparse graphs where explicit NJW eigendecomposition is prohibitive, initializing Weighted Kernel $k$-Means (Alg. 16 with $K_"NC"$, $sigma = 1$) from a $O(t |E|)$ Power Iteration embedding (Alg. 132) monotonically optimizes the *exact discrete* Normalized Cut objective without ever forming a dense $n times n$ matrix!

*3. The Density–Ultrametric–Mode Highway:*
$ &"HDBSCAN (32)" arrow.l.r.double "SLINK"(d_"mreach") " (56)" arrow.l.r.double "MST"(d_"mreach") \
  &quad arrow.l.r.double_(forall epsilon) "DBSCAN* (31)"; \
  &"Mean Shift (34)" arrow.l.r.double "KDE-EM" arrow.l.r.double "DENCLUE 2.0 (43)". $
Theorems 7–9 separate non-parametric density clustering into two exact mathematical engines: *level-set connectivity* (all unified under Kruskal MST on raw or mutual-reachability metrics) and *gradient-flow mode seeking* (unified under step-size-free Gaussian KDE EM). Computing a single mutual-reachability MST (Theorem 8) simultaneously encodes every possible DBSCAN\* cut across all $epsilon > 0$, rendering grid search over $epsilon$ obsolete whenever variable-density clusters coexist.

*4. The Latent-Factor, Network-Physics, and Balanced-OT Highway:*
$ &"pLSA (84)" arrow.l.r.double "KL-NMF (92)" subset.eq "LDA (83)"; \
  &"DC-SBM MLE (131)" arrow.l.r.double_(gamma^*) "Modularity (107–108)" \
  &quad arrow.l.r.double "Potts Spinglass (117)"; \
  &"SwAV (190)" arrow.l.r.double "Entropic OT" \
  &quad ->^(epsilon -> 0^+) "Balanced Spherical "k"-Means (12)". $
Theorems 11–13 reveal how seemingly heuristic algorithmic choices encode exact probabilistic priors: KL-NMF encodes multinomial topic likelihoods (and alternating Lee–Seung updates strictly outperform simultaneous pLSA EM); standard modularity ($gamma = 1$) silently encodes the rigid DC-SBM constraint $omega_"in" - omega_"out" = log(omega_"in" / omega_"out")$, whose violation causes the Fortunato–Barthélemy resolution limit (Corollary 2); and SwAV's Sinkhorn code step imposes exact equipartition constraints via a strongly contractive OT dual, solving the trivial collapse pathology established for unregularized deep clustering in Theorem 16.

Finally, these reductions foreshadow *Comparison Chapter XI*: because every isomorphism and continuous limit preserves or transfers scale-invariance, richness, and consistency properties across families, *no* reduction chain can escape Kleinberg's axiomatic impossibility boundary @kleinberg2002impossibility—every exact equivalence transfers both the statistical strengths and the structural inductive biases of its target functional across the entire 200-algorithm landscape.
