#import "../template/ieee.typ": *

#alg_entry(
  id: 103,
  name: "Spectral Clustering (Ng-Jordan-Weiss)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2002",
  paradigm: "Normalized Laplacian Spectral Embedding",
  complexity: $O(n^2 d + n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Given a dataset $X = {x_1, dots, x_n} subset.eq bb(R)^d$ and a target cluster count $k in {2, dots, n}$, Euclidean prototype methods such as $K$-means (@lloyd1982least) restrict clusters to convex Voronoi polytopes and fail on intertwined manifolds, rings, or non-convex filaments. Ng, Jordan, and Weiss (@ng2002spectral) formulate clustering as a spectral embedding problem on a weighted similarity graph $G = (V, E, W)$ with vertex set $V = {1, dots, n}$, symmetric nonnegative affinity matrix $W in bb(R)_(>= 0)^(n times n)$ having zero diagonal ($W_(i i) = 0$), and diagonal degree matrix $D = "diag"(d_1, dots, d_n)$ where $d_i = sum_(j=1)^n W_(i j) > 0$.

For Euclidean data, the canonical radial basis affinity is $W_(i j) = exp(-norm(x_i - x_j)_2^2 / (2 sigma^2))$ for $i != j$ with scale parameter $sigma > 0$, or a symmetric $q$-nearest-neighbor sparsification. Rather than clustering directly in $bb(R)^d$, NJW forms the normalized affinity matrix $A = D^(-1/2) W D^(-1/2)$---equivalently the symmetric normalized graph Laplacian $cal(L)_"sym" = I_n - A$---extracts its top-$k$ orthonormal eigenvectors, projects each row onto the unit hypersphere $bb(S)^(k-1)$, and partitions the resulting spherical coordinates via $K$-means (@ng2002spectral, @vonluxburg2007tutorial).

#definition(title: "Symmetric Normalized Laplacian and NJW Spherical Map")[
  Let $G = (V, E, W)$ be a weighted undirected graph with strictly positive degrees $d_i = sum_(j=1)^n W_(i j) > 0$ and $D = "diag"(d_1, dots, d_n)$. The *symmetric normalized Laplacian* is
  $ cal(L)_"sym" = I_n - D^(-1/2) W D^(-1/2) = D^(-1/2) (D - W) D^(-1/2). $
  Let $U = [u_1, dots, u_k] in bb(R)^(n times k)$ hold orthonormal eigenvectors corresponding to the $k$ smallest eigenvalues $0 = lambda_1 <= lambda_2 <= dots <= lambda_k$ of $cal(L)_"sym"$ (equivalently the $k$ largest eigenvalues $mu_m = 1 - lambda_m$ of $A = D^(-1/2) W D^(-1/2)$). Assuming $norm(U_(i,:))_2 > 0$ for all $i in V$, the *NJW spherical embedding* $Phi_"NJW" : V -> bb(S)^(k-1)$ maps vertex $i$ to
  $ t_i = (U_(i,:)^top) / norm(U_(i,:))_2, quad T_(i j) = U_(i j) / (sum_(ell=1)^k U_(i ell)^2)^(1/2). $
]

Geometrically, row normalization removes degree-induced radial dispersion so that vertices belonging to a weakly coupled subgraph concentrate around a single canonical ray on $bb(S)^(k-1)$, transforming non-convex connectivity in $bb(R)^d$ into angular separation of tight spherical clusters.

#algbox(
  title: "Ng-Jordan-Weiss (NJW) Spectral Clustering",
  inputs: [Data $X = {x_i}_(i=1)^n subset.eq bb(R)^d$ (or affinity $W$), clusters $k$, scale $sigma > 0$.],
  outputs: [Disjoint partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$.],
)[
  1. *Affinity Construction:* Form $W in bb(R)_(>= 0)^(n times n)$ with $W_(i i) = 0$ and $W_(i j) = exp(-norm(x_i - x_j)_2^2 / (2 sigma^2))$ for $i != j$ (or symmetric $q$-NN affinity).
  2. *Symmetric Normalization:* Compute $d_i = sum_(j=1)^n W_(i j)$, set $D = "diag"(d_1, dots, d_n)$, and form $A = D^(-1/2) W D^(-1/2)$ (or $cal(L)_"sym" = I_n - A$).
  3. *Spectral Decomposition:* Compute $k$ orthonormal eigenvectors $U = [u_1, dots, u_k] in bb(R)^(n times k)$ for the $k$ largest eigenvalues of $A$ (smallest of $cal(L)_"sym"$).
  4. *Spherical Row Normalization:* Form $T in bb(R)^(n times k)$ by normalizing each row $t_i = U_(i,:)^top / norm(U_(i,:))_2$ onto $bb(S)^(k-1)$ (setting $t_i = 0$ if $norm(U_(i,:))_2 = 0$).
  5. *Spherical Partitioning:* Run $K$-means (with $K$-means++ seeding @arthur2007kmeanspp) on ${t_1, dots, t_n}$ to obtain clusters $C_1, dots, C_k$, assigning $x_i$ to $C_m$ iff $t_i$ is assigned to cluster $m$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Spectral Envelope and Connected-Component Nullity")

#lemma(title: "Spectral Bounds and Kernel Characterization", origin: "Original Proof")[
  For any weighted undirected graph $G = (V, E, W)$ with $W_(i j) >= 0$ and $d_i > 0$ for all $i in V$, the symmetric normalized Laplacian $cal(L)_"sym"$ satisfies:
  (i) $"spec"(cal(L)_"sym") subset.eq [0, 2]$, with $lambda_n (cal(L)_"sym") = 2$ if and only if $G$ contains a nontrivial bipartite connected component;
  (ii) the multiplicity of eigenvalue $0$ equals the number $c$ of connected components $C_1, dots, C_c$ of $G$, and $"ker"(cal(L)_"sym") = "span"{D^(1/2) bb(1)_(C_m) : m = 1, dots, c}$.
]

#proof[
  Let $f in bb(R)^n$ be arbitrary and set $g = D^(-1/2) f$, so $f_i = sqrt(d_i) g_i$. Expanding the quadratic form of $cal(L)_"sym"$ yields
  $ f^top cal(L)_"sym" f
    &= sum_(i=1)^n d_i g_i^2 - sum_(i,j=1)^n W_(i j) g_i g_j \
    &= 1/2 sum_(i,j=1)^n W_(i j) (g_i - g_j)^2 >= 0, $
  establishing positive semidefiniteness ($lambda_1 >= 0$). For the upper bound, note that $2 I_n - cal(L)_"sym" = D^(-1/2) (D + W) D^(-1/2)$, whose quadratic form is
  $ f^top (2 I_n - cal(L)_"sym") f = 1/2 sum_(i,j=1)^n W_(i j) (g_i + g_j)^2 >= 0. $
  Hence $0 <= f^top cal(L)_"sym" f <= 2 norm(f)_2^2$ for all $f in bb(R)^n$, proving $"spec"(cal(L)_"sym") subset.eq [0, 2]$ by the Rayleigh-Ritz theorem. Equality $f^top (2 I_n - cal(L)_"sym") f = 0$ holds for $f != 0$ iff $g_i = -g_j$ across every edge $W_(i j) > 0$. On any connected component where $g$ is not identically zero, this forces $|g_i|$ to be a positive constant and partitions its vertices by the sign of $g_i$ with no intra-sign edges, i.e., a nontrivial bipartite component.

  For (ii), $cal(L)_"sym" f = 0$ iff $f^top cal(L)_"sym" f = 0$ (since $cal(L)_"sym" succ.eq 0$), which holds iff $g_i = g_j$ whenever $W_(i j) > 0$. By path induction, $g$ is constant on each connected component $C_m$, so $g = sum_(m=1)^c alpha_m bb(1)_(C_m)$ and $f = D^(1/2) g = sum_(m=1)^c alpha_m D^(1/2) bb(1)_(C_m)$. Because the supports $C_1, dots, C_c$ are disjoint and $d_i > 0$, the vectors ${D^(1/2) bb(1)_(C_m)}_(m=1)^c$ are mutually orthogonal and nonzero, so $dim("ker"(cal(L)_"sym")) = c$.
]

#subsec("B", "Exact Orthogonal Collapse Under Row Normalization")

A defining innovation of NJW over unnormalized spectral embeddings is the row-normalization step $T_(i j) = U_(i j) / norm(U_(i,:))_2$. The following theorem proves why this projection is both necessary and sufficient for exact degree-invariant recovery in the ideal block-graph regime.

#theorem(title: "Rotation- and Degree-Invariant Simplex Collapse", origin: "Original Derivation")[
  Suppose $G$ consists of $k$ disconnected, non-bipartite connected components $C_1, dots, C_k$ with volumes $"vol"(C_m) = sum_(i in C_m) d_i > 0$. Let $U in bb(R)^(n times k)$ be *any* orthonormal basis for $"ker"(cal(L)_"sym")$. Then:
  (i) there exists a unique orthogonal matrix $R in O(k)$ with rows $r_1^top, dots, r_k^top$ such that every row of the NJW matrix $T$ satisfies $t_i = r_(m(i))$ whenever $i in C_(m(i))$;
  (ii) the pairwise distances in $T$ satisfy $norm(t_i - t_j)_2 = 0$ if $m(i) = m(j)$ and $norm(t_i - t_j)_2 = sqrt(2)$ if $m(i) != m(j)$, with zero intra-cluster $K$-means distortion;
  (iii) conversely, without row normalization, the raw rows $u_i = U_(i,:)^top$ satisfy $norm(u_i)_2 = sqrt(d_i / "vol"(C_(m(i))))$, so whenever two nodes $i in C_1, j in C_2$ satisfy $d_i / "vol"(C_1) + d_j / "vol"(C_2) < d_i / "vol"(C_1) + d_(i') / "vol"(C_1) - 2 sqrt(d_i d_(i')) / "vol"(C_1)$ for a high-degree node $i' in C_1$, $i$ is strictly closer in $bb(R)^k$ to $j in C_2$ than to $i' in C_1$.
]

#proof[
  Define the canonical orthonormal basis $V in bb(R)^(n times k)$ for $"ker"(cal(L)_"sym")$ via Lemma 1 by setting $V_(i m) = sqrt(d_i / "vol"(C_m)) bb(I)(i in C_m)$. Clearly $V^top V = I_k$ because disjoint supports vanish off-diagonally and $sum_(i in C_m) d_i / "vol"(C_m) = 1$. Since $U$ and $V$ span the same $k$-dimensional eigenspace $"ker"(cal(L)_"sym")$, there exists a unique orthogonal matrix $R in O(k)$ such that $U = V R$. For any $i in C_m$, row $i$ of $V$ is $sqrt(d_i / "vol"(C_m)) e_m^top$, so row $i$ of $U$ is
  $ u_i = U_(i,:)^top = sqrt(d_i / "vol"(C_m)) R^top e_m = sqrt(d_i / "vol"(C_m)) r_m, $
  where $r_m = R^top e_m$ is the $m$-th row of $R$ written as a column vector. Since $R in O(k)$, $norm(r_m)_2 = 1$ and $chevron.l r_m, r_(m') chevron.r = delta_(m m')$. Consequently $norm(u_i)_2 = sqrt(d_i / "vol"(C_m)) > 0$, and normalizing row $i$ cancels the degree factor identically:
  $ t_i = u_i / norm(u_i)_2 = r_m. $
  For any $i in C_m$ and $j in C_(m')$,
  $ norm(t_i - t_j)_2^2 = norm(r_m)_2^2 + norm(r_(m'))_2^2 - 2 chevron.l r_m, r_(m') chevron.r = 2 (1 - delta_(m m')). $
  Thus $T$ collapses all $n$ points onto $k$ mutually orthogonal unit vertices on $bb(S)^(k-1)$ at mutual distance $sqrt(2)$, giving global minimum $K$-means SSE $= 0$ exclusively for the true partition ${C_1, dots, C_k}$.

  For (iii), without row normalization, the squared distance between $u_i$ ($i in C_1$) and $u_j$ ($j in C_2$) is $norm(u_i - u_j)_2^2 = norm(u_i)_2^2 + norm(u_j)_2^2 = p_i + p_j$ where $p_i = d_i / "vol"(C_(m(i)))$, while between $u_i$ and $u_(i')$ ($i' in C_1$) it is $(sqrt(p_i) - sqrt(p_(i')))^2$. Whenever $p_i, p_j -> 0$ while $p_(i') = Theta(1)$, we have $norm(u_i - u_j)_2^2 = p_i + p_j < (sqrt(p_(i')) - sqrt(p_i))^2 = norm(u_i - u_(i'))_2^2$, forcing $K$-means on raw $U$ to merge low-degree nodes across distinct components around the origin.
]

#subsec("C", "Perturbation Stability and Subspace Distortion")

In realistic graphs, inter-cluster edges are nonzero rather than identically zero. Let $A = A^((0)) + E$, where $A^((0)) = "diag"(A_1^((0)), dots, A_k^((0)))$ is an ideal $k$-block normalized affinity with internal spectral gaps, and $E$ captures inter-cluster affinities and degree perturbation.

#theorem(title: "Davis-Kahan Perturbation Bound for NJW", origin: "Classical — Cited")[
  Let $A = D^(-1/2) W D^(-1/2)$ and ideal block-diagonal $A^((0))$ have eigenvalues $1 = mu_1^((0)) = dots = mu_k^((0)) > mu_(k+1)^((0)) >= dots >= mu_n^((0))$ with eigengap $gamma = 1 - mu_(k+1)^((0)) > 0$ (@ng2002spectral, @davis1970rotation). Suppose $norm(E)_F = norm(A - A^((0)))_F < gamma / 2$. Let $U, V in bb(R)^(n times k)$ contain the top-$k$ orthonormal eigenvectors of $A$ and $A^((0))$, respectively. Then there exists an orthogonal matrix $R in O(k)$ such that
  $ norm(U - V R)_F <= (sqrt(2) norm(E)_F) / (gamma - norm(E)_2) <= (2 sqrt(2) norm(E)_F) / gamma. $
  Moreover, if $rho_"min" = min_i sqrt(d_i^((0)) / "vol"^((0))(C_(m(i)))) > 0$ and the normalized perturbation satisfies $eta_i = norm(u_i - (V R)_(i,:)^top)_2 < rho_"min" / 2$, then the spherical embedding error obeys
  $ norm(t_i - r_(m(i)))_2 <= (2 eta_i) / rho_"min", $
  and the fraction of vertices with spherical displacement exceeding $tau in (0, sqrt(2))$ is bounded by $8 norm(E)_F^2 / (n gamma^2 rho_"min"^2 tau^2)$.
]

#proof[
  By Weyl's inequality, $|mu_(k+1)(A) - mu_(k+1)^((0))| <= norm(E)_2 <= norm(E)_F < gamma / 2$, so the spectral separation between the top-$k$ spectrum of $A^((0))$ (all equal to $1$) and the bottom-$(n-k)$ spectrum of $A$ is at least $tilde(gamma) = 1 - mu_(k+1)(A) >= gamma - norm(E)_2 > gamma / 2$. Applying the Davis-Kahan $sin Theta$ theorem (@davis1970rotation) to the invariant subspaces spanned by $U$ and $V$ yields $norm(sin Theta(U, V))_F <= norm(E)_F / tilde(gamma)$. Since the Procrustes distance satisfies $inf_(R in O(k)) norm(U - V R)_F <= sqrt(2) norm(sin Theta(U, V))_F$, choosing the minimizing $R in O(k)$ yields
  $ norm(U - V R)_F <= (sqrt(2) norm(E)_F) / (gamma - norm(E)_2) <= (2 sqrt(2) norm(E)_F) / gamma. $
  For any vertex $i$, let $v_i = (V R)_(i,:)^top = rho_i r_(m(i))$ with $rho_i >= rho_"min"$. Using the elementary unit-sphere Lipschitz identity $norm(a / norm(a)_2 - b / norm(b)_2)_2 <= 2 norm(a - b)_2 / norm(b)_2$ (valid for all $a, b != 0$, since $norm(a / norm(a)_2 - b / norm(b)_2)_2 <= norm(a - b)_2 / norm(b)_2 + |norm(b)_2 - norm(a)_2| / norm(b)_2 <= 2 norm(a - b)_2 / norm(b)_2$), we obtain $norm(t_i - r_(m(i)))_2 <= 2 eta_i / rho_i <= 2 eta_i / rho_"min"$ whenever $eta_i < rho_"min" / 2$ (guaranteeing $u_i != 0$). Finally, by Markov's inequality applied to $sum_(i=1)^n norm(t_i - r_(m(i)))_2^2 <= (4 / rho_"min"^2) sum_(i=1)^n eta_i^2 = (4 / rho_"min"^2) norm(U - V R)_F^2$, the number $n_"bad"$ of vertices with $norm(t_i - r_(m(i)))_2 >= tau$ satisfies
  $ n_"bad" / n <= (4 norm(U - V R)_F^2) / (n rho_"min"^2 tau^2) <= (32 norm(E)_F^2) / (n gamma^2 rho_"min"^2 tau^2) $
  (with tighter factor $8$ when $norm(E)_2 << gamma$, as $gamma - norm(E)_2 -> gamma$).
]

#subsec("D", "Spherical K-Means Equivalence, Invariances, and Scale Failure")

#proposition(title: "Spherical Cosine Equivalence and Invariance Group", origin: "Original Proof")[
  For any NJW embedding ${t_i}_(i=1)^n subset.eq bb(S)^(k-1)$:
  (i) minimizing the Euclidean $K$-means objective $J_"SSE"(cal(C)) = sum_(m=1)^k sum_(i in C_m) norm(t_i - mu_m)_2^2$ over partitions $cal(C)$ is strictly equivalent to maximizing the within-cluster mean-resultant energy $sum_(m=1)^k |C_m| norm(mu_m)_2^2$, where $mu_m = |C_m|^(-1) sum_(i in C_m) t_i$;
  (ii) the NJW embedding (up to global orthogonal rotation in $O(k)$) is strictly invariant under positive affinity scaling $W |-> alpha W$ ($alpha > 0$), vertex relabeling $W |-> P W P^top$, and rigid Euclidean motions $x_i |-> Q x_i + b$ ($Q in O(d)$) for isotropic affinities.
]

#proof[
  Because $norm(t_i)_2^2 = 1$ for every $i in V$, expanding the $K$-means SSE gives
  $ J_"SSE"(cal(C))
    &= sum_(m=1)^k sum_(i in C_m) (norm(t_i)_2^2 - 2 chevron.l t_i, mu_m chevron.r + norm(mu_m)_2^2) \
    &= n - sum_(m=1)^k |C_m| norm(mu_m)_2^2 \
    &= n - sum_(m=1)^k 1 / |C_m| sum_(i,j in C_m) chevron.l t_i, t_j chevron.r. $
  Since $n$ is constant, minimizing $J_"SSE"$ is strictly equivalent to maximizing $sum_(m=1)^k |C_m| norm(mu_m)_2^2$, i.e., maximizing average intra-cluster cosine similarity on $bb(S)^(k-1)$.

  For (ii), scaling $W |-> alpha W$ scales $D |-> alpha D$, leaving $A = (alpha D)^(-1/2) (alpha W) (alpha D)^(-1/2) = D^(-1/2) W D^(-1/2)$ identically unchanged. Under vertex permutation $P$, $A' = P A P^top$ has eigenvectors $U' = P U$ (up to right-multiplication by $R in O(k)$ when eigenvalues are degenerate), so $T' = P T R$, preserving all pairwise distances on $bb(S)^(k-1)$. Rigid motions leave $norm(x_i - x_j)_2$ invariant, hence preserve $W$ and $T$.
]

#proposition(title: "Fixed-Bandwidth Collapse Under Multi-Scale Density", origin: "Original Counterexample")[
  Consider $k = 2$ one-dimensional clusters: a dense cluster $C_1$ of $m$ points spaced by $epsilon$ in $[0, (m-1) epsilon]$ and a sparse cluster $C_2$ of $m$ points spaced by $Delta$ in $[L, L + (m-1) Delta]$, where $epsilon << L << Delta$. For *every* global Gaussian bandwidth $sigma > 0$, the eigengap $lambda_3 (cal(L)_"sym") - lambda_2 (cal(L)_"sym")$ cannot remain bounded away from zero as $Delta / L -> infinity$ and $L / epsilon -> infinity$: either $sigma <= L / sqrt(2 ln m)$ disconnects $C_2$ into $m$ singletons ($lambda_3 -> 0$), or $sigma >= Delta / 2$ merges $C_1$ and the nearest point of $C_2$ with affinity $exp(-L^2 / (2 sigma^2)) -> 1$.
]

#proof[
  Let $x_a$ be the rightmost point of $C_1$ and $x_b$ the leftmost point of $C_2$, separated by $L$, while consecutive points in $C_2$ are separated by $Delta >> L$.
  Case 1: $sigma <= L / sqrt(2 ln m)$. Every edge incident to an interior point $x_r in C_2$ ($r >= 2$) has length at least $Delta$, hence weight at most $exp(-Delta^2 / (2 sigma^2)) <= exp(-((Delta / L)^2) ln m) = m^(-(Delta/L)^2)$. Removing the $m - 1$ points of $C_2$ beyond $x_b$ cuts at most $O(m^2)$ edges of weight $m^(-(Delta/L)^2)$, while the isolated subgraph has $m - 1$ components up to $O(m^(2 - (Delta/L)^2))$ perturbation, forcing $lambda_3 (cal(L)_"sym") <= dots <= lambda_m (cal(L)_"sym") -> 0$ as $Delta / L -> infinity$.
  Case 2: $sigma >= Delta / 2$. Then the inter-cluster affinity between $C_1$ and $x_b$ satisfies $W_(a b) = exp(-L^2 / (2 sigma^2)) >= exp(-2 (L / Delta)^2) -> 1$, while every intra-$C_2$ edge has weight at most $exp(-Delta^2 / (2 sigma^2))$. Hence $x_b$ is strictly more strongly coupled to $C_1$ than to $C_2$, breaking the $k = 2$ cluster separation.
  For intermediate $sigma in (L / sqrt(2 ln m), Delta / 2)$, the ratio $W_(i, i+1)^(C_2) / W_(a b) = exp(-(Delta^2 - L^2) / (2 sigma^2)) <= exp(-2 (1 - (L/Delta)^2)) < 1$ still makes the inter-cluster bridge $W_(a b)$ exponentially stronger than every internal link of $C_2$, proving impossibility for *all* $sigma > 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact dense NJW requires three sequential phases: pairwise squared-distance and affinity construction in $O(n^2 d)$ time and $O(n^2)$ memory; dense symmetric eigendecomposition (Householder tridiagonalization plus implicit QR or divide-and-conquer) in $O(n^3)$ worst-case time (or $O(n^2 k)$ via truncated Lanczos / Arnoldi when $k << n$ and $A$ is well-gapped); and $K$-means on $T in bb(R)^(n times k)$ in $O(I n k^2)$ time for $I$ Lloyd iterations. When $W$ is sparsified to a symmetric $q$-nearest-neighbor graph with $m = O(q n)$ nonzeros, matrix-vector products drop to $O(q n)$, reducing ARPACK / Lanczos spectral extraction to $O(T_"iter" q n k + n k^2)$ time and $O(q n + n k)$ memory, though convergence slows as $O(gamma^(-1/2))$ when the eigengap $gamma = lambda_(k+1) - lambda_k$ shrinks.

Geometrically, NJW maps arbitrary non-convex manifolds into tight clusters near $k$ orthogonal vertices on $bb(S)^(k-1)$, rendering final $K$-means largely immune to elongated shapes in the ambient space. However, NJW lacks a closed-form out-of-sample extension (requiring Nystr\"om interpolation for unseen queries), depends critically on the affinity graph bandwidth $sigma$ or neighbor count $q$, and degrades in high ambient dimension $d$ when Euclidean distances concentrate and flatten the Gaussian kernel.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Exact Dense Complexity], [$O(n^2 d + n^3)$ time; $O(n^2)$ memory for $W$ and $A$.],
    [Sparse Lanczos Profile], [$O(T_"iter" q n k + I n k^2)$ time; $O(q n + n k)$ memory.],
    [Cluster Geometry], [Non-convex manifolds, rings, filaments via $bb(S)^(k-1)$ rays.],
    [Invariances], [Positive $W$-scale, permutation, rigid ambient motion.],
    [Hyperparameters], [Cluster count $k$; kernel scale $sigma$ or $q$-NN graph.],
    [Key Sensitivity], [Multi-scale density gaps; eigengap collapse; $d$-concentration.],
  ),
  caption: [Operational and geometric profile of Ng-Jordan-Weiss spectral clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell Transcriptomics and Cell-Type Subpopulation Discovery:* On $q$-NN gene-expression graphs built over PCA/diffusion coordinates, cell lineages form curved developmental manifolds with heterogeneous local sampling densities. Symmetric degree normalization $D^(-1/2) W D^(-1/2)$ dampens hyper-connected housekeeping hubs, while row projection onto $bb(S)^(k-1)$ equalizes rare progenitor cells and abundant mature types before final $K$-means partitioning (@vonluxburg2007tutorial).
2. *Perceptual Image and Superpixel Grouping:* When pixel or superpixel affinities combine spatial proximity and color/texture similarity, semantic regions correspond to weakly coupled subgraphs with irregular silhouettes. NJW embeds multi-region affinity graphs directly into $k$ orthogonal spherical modes, avoiding the greedy recursive bias of sequential bipartitioning (@ng2002spectral, @shi2000ncut).
3. *Subspace and Motion Segmentation in Computer Vision:* In sparse subspace clustering (@elhamifar2013ssc) and multi-body structure-from-motion, self-expressive affinity matrices $W = |C| + |C|^top$ are block-diagonal up to noise even when trajectories lie on intersecting low-dimensional subspaces. Theorem 1 and Theorem 2 guarantee that NJW collapses each subspace's affinity block onto a distinct orthogonal vertex of $bb(S)^(k-1)$ regardless of point density variation within each rigid body.
4. *Multi-View and Kernel Graph Fusion in Bioinformatics:* Co-regularized multi-view spectral clustering (@kumar2011coregularized) aligns normalized Laplacian eigenspaces across heterogeneous omics modalities (expression, methylation, PPI networks). Because NJW normalizes every node's spectral embedding to unit norm, views with vastly different degree distributions combine additively on $bb(S)^(k-1)$ without high-degree modalities dominating the consensus partition.
5. *Biomolecular Conformational Metastability Analysis:* In molecular dynamics trajectory graphs where transition weights reflect kinetic overlap between microstates, the leading $k$ eigenvectors of $A = D^(-1/2) W D^(-1/2)$ span the slow metastable macrostates. Row normalization isolates transition-region microstates cleanly onto geodesic arcs between orthogonal macrostate vertices.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

NJW exhibits four characteristic failure modes: (i) *bandwidth / multi-scale mismatch* (Proposition 2), where a single global $sigma$ disconnects sparse clusters or fuses dense ones; (ii) *high-dimensional kernel flattening*, where distance concentration drives off-diagonal $W_(i j)$ toward a constant or identity matrix, collapsing the eigengap $gamma$; (iii) *extreme degree outliers on disconnected or near-disconnected graphs*, where row normalization $u_i / norm(u_i)_2$ amplifies floating-point noise on near-zero rows; and (iv) *cubic scaling* $O(n^3)$ for dense graphs beyond $n approx 10^4$, necessitating Nystr\"om sampling or Power Iteration Clustering (@lin2010pic).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Non-convex 2D/3D manifolds], [High: $bb(S)^(k-1)$ collapses curved components.], [Tune $sigma$ or use $q$-NN graph.],
    [Subspace / affinity graphs], [High: Theorem 2 gives exact block recovery.], [Check eigengap $gamma$ for $k$.],
    [Multi-scale density data], [Low: global $sigma$ fails (Prop. 2).], [Use HDBSCAN (@campello2013hdbscan) or self-tuning $sigma_i$.],
    [High-D raw Euclidean data], [Low: RBF kernel concentrates.], [Use Ncut (Alg. 104) or PCA + cosine $q$-NN.],
    [Web-scale graphs ($n > 10^5$)], [Poor dense $O(n^3)$ cost.], [Use Leiden (@traag2019leiden) or PIC (@lin2010pic).],
  ),
  caption: [Domain suitability and failure-mode matrix for NJW spectral clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

NJW occupies a central vertex in the spectral-clustering triangle alongside Normalized Cuts (Alg. 104, @shi2000ncut) and Ratio Cut (Alg. 105, @hagen1992ratiocut), and connects algebraically to weighted Kernel $K$-means (@dhillon2004kernel).

#proposition(title: "Spectral Coupling Between NJW and Shi-Malik Ncut", origin: "Original Proof")[
  Let $U = [u_1, dots, u_k]$ be the top-$k$ eigenvectors of $A = D^(-1/2) W D^(-1/2)$ (NJW) and $Y = [y_1, dots, y_k]$ the bottom-$k$ eigenvectors of the random-walk Laplacian $cal(L)_"rw" = I_n - D^(-1) W$ (Ncut). Then $Y = D^(-1/2) U$ (up to column scaling) has identical eigenvalues $lambda_m (cal(L)_"sym") = lambda_m (cal(L)_"rw")$, and the NJW row-normalized vectors $t_i = U_(i,:)^top / norm(U_(i,:))_2$ equal the row-normalized rows of $Y$: $t_i = Y_(i,:)^top / norm(Y_(i,:))_2$ for all $i in V$.
]

#proof[
  If $cal(L)_"sym" u = lambda u$, multiplying on the left by $D^(-1/2)$ gives $D^(-1/2) (I_n - D^(-1/2) W D^(-1/2)) u = (I_n - D^(-1) W) (D^(-1/2) u) = lambda (D^(-1/2) u)$, so $y = D^(-1/2) u$ is an eigenvector of $cal(L)_"rw"$ with identical eigenvalue $lambda$. Thus row $i$ of $Y = D^(-1/2) U$ satisfies $Y_(i,:) = d_i^(-1/2) U_(i,:)$. Since $d_i^(-1/2) > 0$ is a positive scalar common to all $k$ coordinates of row $i$, dividing by its Euclidean norm cancels $d_i^(-1/2)$ identically: $Y_(i,:)^top / norm(Y_(i,:))_2 = U_(i,:)^top / norm(U_(i,:))_2 = t_i$.
]

The practical divergence between Alg. 103 and Alg. 104 therefore arises entirely from *row normalization*: whereas Shi-Malik Ncut preserves the degree-rescaled coordinates $Y = D^(-1/2) U$ (corresponding directly to weighted Kernel $K$-means with node weights $d_i$, @dhillon2004kernel), NJW projects every row to unit length before unweighted $K$-means, giving low-degree boundary points equal geometric leverage on $bb(S)^(k-1)$.

*Empirical Benchmark Analysis:* On the 4-dataset ClusterOverview benchmark (evaluated with Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), NJW achieves its strongest performance on *Multi-Shape* (`Sil = 0.592`, `ARI = N/A`), confirming its theoretical ability (Theorem 1) to untangle non-convex geometric shapes via normalized Laplacian eigenvectors. However, its sensitivity to global bandwidth calibration and row-norm amplification of border noise surfaces sharply on the remaining benchmarks: on *Iris* (`Sil = 0.132`, `ARI = 0.296`), overlapping Versicolor/Virginica border samples suffer from spherical amplification compared to Ncut (`Sil = 0.382`, `ARI = 0.414`); on *Customers* (`Sil = 0.008`), multi-scale feature dispersion triggers the exact bandwidth collapse predicted by Proposition 2; and on *High-D* (`Sil = -0.028`, `ARI = N/A`), isotropic Gaussian kernel concentration flattens the affinity spectrum into negative Silhouette territory. With a mean runtime of `4.226 ms`, NJW matches its spectral peers (Ncut `3.882 ms`, Ratio Cut `4.568 ms`) while running roughly $3 times$ faster than iterative message-passing Affinity Propagation (`13.438 ms`), underscoring that NJW excels specifically when affinity graphs are locally adaptive or $q$-NN sparsified on low-dimensional manifolds.
