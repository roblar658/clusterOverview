#import "../template/ieee.typ": *

#alg_entry(
  id: 105,
  name: "Ratio Cut",
  category: "Cat. V: Graph, Spectral & Network",
  year: "1992",
  paradigm: "Unnormalized Laplacian Spectral Partitioning",
  complexity: $O(n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Early graph partitioning algorithms such as Kernighan-Lin (@kernighan1970efficient) imposed rigid hard cardinality constraints ($|A| = |B| = n / 2$) to prevent minimum-cut partitioning from peeling off isolated vertices. In VLSI circuit layout, scientific mesh decomposition, and data clustering, however, natural modules rarely have exact equal size, yet grossly unbalanced cuts remain undesirable. Hagen and Kahng (@hagen1992ratiocut) resolved this dilemma by introducing the *Ratio Cut* objective, which penalizes small cluster cardinality directly in the denominator of the cut cost and admits an exact spectral relaxation via the *unnormalized* (combinatorial) graph Laplacian $L = D - W$.

Let $G = (V, E, W)$ be a connected, weighted undirected graph on $V = {1, dots, n}$ with symmetric nonnegative weight matrix $W in bb(R)_(>= 0)^(n times n)$, vertex degrees $d_i = sum_(j=1)^n W_(i j)$, degree matrix $D = "diag"(d_1, dots, d_n)$, and combinatorial Laplacian $L = D - W$. Unlike Normalized Cuts (Alg. 104, @shi2000ncut), which normalizes cut weight by edge volume $"vol"(C_m)$, Ratio Cut normalizes by vertex cardinality $|C_m|$.

#definition(title: "Hagen-Kahng Ratio Cut and k-Way Ratio Cut")[
  For a bipartition $(A, B)$ of $V$ with $|A|, |B| >= 1$, the *Hagen-Kahng two-way Ratio Cut* (@hagen1992ratiocut) is
  $ "Rcut"_"HK"(A, B) = ("cut"(A, B)) / (|A| |B|). $
  For a general $k$-partition $cal(C) = (C_1, dots, C_k)$ with $|C_m| >= 1$ for all $m$, the standard *$k$-way Ratio Cut* objective (@vonluxburg2007tutorial) is
  $ "Rcut"(C_1, dots, C_k) = sum_(m=1)^k ("cut"(C_m, overline(C)_m)) / (|C_m|), $
  where $overline(C)_m = V without C_m$. Note that for $k = 2$, $"Rcut"(A, B) = "cut"(A, B) (1 / |A| + 1 / |B|) = n dot "Rcut"_"HK"(A, B)$, so both two-way formulations share identical optimal partitions.
]

Relaxing the cardinality-normalized cluster indicator vectors to arbitrary real orthonormal vectors produces the standard eigenproblem $L u = lambda u$. In two-way mode, the Fiedler vector $u_2$ (second smallest eigenvector of $L$) orders vertices for a linear-time exact Ratio Cut sweep; in $k$-way mode, the bottom-$k$ eigenvectors $U = [u_1, dots, u_k]$ of $L$ provide the spectral coordinates clustered by $K$-means (@vonluxburg2007tutorial).

#algbox(
  title: "Spectral Ratio Cut Partitioning",
  inputs: [Graph $G = (V, E, W)$ (or data $X$ + affinity), clusters $k$, mode (recursive 2-way or $k$-way).],
  outputs: [Partition $cal(C) = (C_1, dots, C_k)$ approximating $min "Rcut"(cal(C))$.],
)[
  1. *Unnormalized Laplacian:* Compute $d_i = sum_(j=1)^n W_(i j)$, set $D = "diag"(d_1, dots, d_n)$, and form $L = D - W$.
  2. *Eigendecomposition:* Compute the $k$ smallest eigenvalues $0 = lambda_1 <= lambda_2 <= dots <= lambda_k$ and corresponding orthonormal eigenvectors $U = [u_1, dots, u_k] in bb(R)^(n times k)$ of $L$.
  3. *2-Way Fiedler Sweep (Recursive Mode):* Sort vertices by $u_(2,i)$, compute exact $"Rcut"(A_r, overline(A)_r) = n dot "cut"(A_r, overline(A)_r) / (r (n - r))$ over all $n - 1$ prefix cuts $A_r$ in $O(m + n log n)$ time, split at $r^* = arg min_r "Rcut"(A_r, overline(A)_r)$, and recurse until $k$ blocks exist.
  4. *Simultaneous $k$-Way Rounding:* Cluster the $n$ rows of $U in bb(R)^(n times k)$ using $K$-means (with $K$-means++ initialization @arthur2007kmeanspp) to produce $C_1, dots, C_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Combinatorial Laplacian Trace Identity and Ky Fan Bound")

#theorem(title: "Exact Ratio Cut Trace Identity and Spectral Lower Bound", origin: "Original Proof")[
  Let $cal(C) = (C_1, dots, C_k)$ be any $k$-partition of $V$ with $|C_m| >= 1$. Define $H in bb(R)^(n times k)$ by $H_(i m) = |C_m|^(-1/2) bb(I)(i in C_m)$. Then:
  (i) $H^top H = I_k$ and $"Tr"(H^top L H) = "Rcut"(C_1, dots, C_k)$;
  (ii) for $k = 2$ with $cal(C) = (A, B)$, defining $f_i = sqrt(|B| / |A|)$ for $i in A$ and $f_i = -sqrt(|A| / |B|)$ for $i in B$ yields $f^top bb(1) = 0$, $norm(f)_2^2 = n$, and $f^top L f / norm(f)_2^2 = "Rcut"(A, B) = n dot "Rcut"_"HK"(A, B)$;
  (iii) every $k$-partition obeys the spectral lower bound
  $ "Rcut"(C_1, dots, C_k) >= sum_(m=1)^k lambda_m (L), quad "and" quad min_(A, B) "Rcut"_"HK"(A, B) >= (lambda_2 (L)) / n. $
]

#proof[
  For (i), the inner product of columns $h_m, h_(m')$ of $H$ is
  $ chevron.l h_m, h_(m') chevron.r = sum_(i=1)^n H_(i m) H_(i m') = delta_(m m') / |C_m| sum_(i in C_m) 1 = delta_(m m'), $
  proving $H^top H = I_k$. Using the combinatorial Laplacian quadratic form $v^top L v = 1/2 sum_(i,j) W_(i j) (v_i - v_j)^2$, each column contributes
  $ h_m^top L h_m
    &= 1/2 sum_(i,j=1)^n W_(i j) (H_(i m) - H_(j m))^2 \
    &= 1 / |C_m| sum_(i in C_m, j in overline(C)_m) W_(i j) = ("cut"(C_m, overline(C)_m)) / |C_m|. $
  Summing over $m = 1, dots, k$ establishes $"Tr"(H^top L H) = "Rcut"(C_1, dots, C_k)$.

  For (ii), direct evaluation yields $f^top bb(1) = |A| sqrt(|B| / |A|) - |B| sqrt(|A| / |B|) = 0$ and $norm(f)_2^2 = |A| (|B| / |A|) + |B| (|A| / |B|) = |B| + |A| = n$. Across every cut edge $(i, j) in A times B$, $(f_i - f_j)^2 = (sqrt(|B|/|A|) + sqrt(|A|/|B|))^2 = (|A| + |B|)^2 / (|A||B|) = n^2 / (|A||B|)$. Thus
  $ (f^top L f) / norm(f)_2^2
    &= 1 / n dot (n^2 "cut"(A, B)) / (|A| |B|) = n dot "Rcut"_"HK"(A, B) = "Rcut"(A, B). $
  For (iii), relaxing $H$ from discrete indicators to the full Stiefel manifold ${H in bb(R)^(n times k) : H^top H = I_k}$ applies Ky Fan's minimum trace principle:
  $ min_(cal(C)) "Rcut"(cal(C)) >= min_(H^top H = I_k) "Tr"(H^top L H) = sum_(m=1)^k lambda_m (L). $
  For $k = 2$, $lambda_1 (L) = 0$ (with eigenvector $bb(1)$), so $min_(A, B) "Rcut"(A, B) >= lambda_2 (L)$, and dividing by $n$ yields Hagen and Kahng's lower bound $min_(A, B) "Rcut"_"HK"(A, B) >= lambda_2 (L) / n$.
]

#subsec("B", "Unnormalized Cheeger Inequality and Maximum-Degree Slack")

Define the *edge expansion* (combinatorial isoperimetric number) of $G$ as $h_"iso"(G) = min_(1 <= |A| <= n/2) "cut"(A, overline(A)) / |A|$. Notice that for $|A| <= n/2$, $h_"iso"(A) <= "Rcut"(A, overline(A)) <= 2 h_"iso"(A)$. Unlike Normalized Cut's degree-free Cheeger constant (Alg. 104), the unnormalized sweep guarantee inherits an explicit dependence on the maximum degree $d_"max" = max_i d_i$.

#theorem(title: "Unnormalized Cheeger Bound with Degree Factor", origin: "Classical — Cited")[
  Let $G = (V, E, W)$ be a connected weighted graph with maximum weighted degree $d_"max" = max_(i in V) d_i$ and Fiedler value $lambda_2 = lambda_2 (L)$ (@hagen1992ratiocut, @chung1997spectral). The optimal prefix sweep cut $A^*$ along the Fiedler vector $u_2$ satisfies
  $ lambda_2 / 2 <= h_"iso"(G) <= h_"iso"(A^*) <= sqrt(2 d_"max" lambda_2), $
  and hence
  $ lambda_2 <= "Rcut"(A^*, overline(A)^*) <= 2 sqrt(2 d_"max" lambda_2). $
]

#proof[
  The lower bound $lambda_2 / 2 <= h_"iso"(G)$ follows immediately from Theorem 1(iii) since $lambda_2 <= "Rcut"(A, overline(A)) <= 2 h_"iso"(A)$ for every $1 <= |A| <= n / 2$.

  For the upper bound, shift $u_2$ by its unweighted median $m_"med"$ so that both $f^+ = max(u_2 - m_"med" bb(1), 0)$ and $f^- = max(m_"med" bb(1) - u_2, 0)$ have support size at most $n / 2$. Choosing $g in {f^+, f^-}$ with $(g^top L g) / norm(g)_2^2 <= lambda_2$ and sorting vertices so $g_1 >= g_2 >= dots >= g_n >= 0$, set $S_r = {1, dots, r}$ and $h^* = min_(r <= n/2) "cut"(S_r, overline(S)_r) / |S_r| >= h_"iso"(A^*)$. By the co-area expansion,
  $ sum_(i < j) W_(i j) |g_i^2 - g_j^2|
    &= sum_(r=1)^(n-1) (g_r^2 - g_(r+1)^2) "cut"(S_r, overline(S)_r) \
    &>= h^* sum_(r=1)^(n-1) (g_r^2 - g_(r+1)^2) |S_r| = h^* norm(g)_2^2. $
  Applying Cauchy-Schwarz to the left-hand side gives
  $ sum_(i < j) W_(i j) |g_i - g_j| (g_i + g_j)
    <= (g^top L g)^(1/2) (2 sum_(i=1)^n d_i g_i^2)^(1/2). $
  Crucially, because the inner product in $L u = lambda u$ is unweighted ($norm(g)_2^2$ rather than $g^top D g$), we must bound $sum_i d_i g_i^2 <= d_"max" norm(g)_2^2$. Substituting $g^top L g <= lambda_2 norm(g)_2^2$ gives
  $ h^* norm(g)_2^2 <= (lambda_2 norm(g)_2^2)^(1/2) (2 d_"max" norm(g)_2^2)^(1/2) = sqrt(2 d_"max" lambda_2) norm(g)_2^2. $
  Canceling $norm(g)_2^2 > 0$ yields $h_"iso"(A^*) <= sqrt(2 d_"max" lambda_2)$, and multiplying by $2$ gives the Ratio Cut bound.
]

#subsec("C", "Degree-Outlier Localization and Asymptotic Inconsistency")

The factor $sqrt(d_"max")$ in Theorem 2 is not a proof artifact; it reflects a fundamental pathology of unnormalized $L = D - W$: low-degree peripheral vertices hijack the bottom spectrum of $L$ before global cluster structure appears.

#theorem(title: "Outlier Eigenvector Localization Counterexample", origin: "Original Counterexample")[
  Consider a graph $G_n$ on $n = 2 m + 1$ vertices formed by two cliques $K_m^((1)), K_m^((2))$ of unit-weight edges joined by $b >= 1$ unit-weight bridge edges (with $2 b < m^2 / 4$), plus a single peripheral outlier vertex $v_0$ connected to $K_m^((1))$ by one edge of weight $epsilon > 0$. Whenever $epsilon < (2 b) / (m + 1)$, the Fiedler value of $L(G_n)$ obeys $lambda_2 (L) <= epsilon dot n / (n - 1)$, whereas the balanced inter-clique cut has $"Rcut"(K_m^((1)) union {v_0}, K_m^((2))) = (2 m + 1) b / (m(m + 1)) > epsilon dot n / (n - 1)$. Moreover, as $epsilon -> 0^+$, the Fiedler vector $u_2$ converges in $ell_2$ to the localized spike $(1 - 1 / n)^(-1/2) (e_(v_0) - n^(-1) bb(1))$, completely ignoring the two-clique community structure.
]

#proof[
  Consider the zero-mean test vector $v = e_(v_0) - n^(-1) bb(1)$, which satisfies $v^top bb(1) = 0$ and $norm(v)_2^2 = (1 - 1 / n)^2 + (n - 1) (1 / n)^2 = (n - 1) / n$. Because $v$ is constant equal to $-1 / n$ on all vertices in $V without {v_0}$, every edge internal to $K_m^((1)) union K_m^((2))$ contributes zero to $v^top L v$. Only the single edge of weight $epsilon$ incident to $v_0$ contributes, with squared difference $(1 - 1 / n - (-1 / n))^2 = 1$. Hence
  $ lambda_2 (L) <= (v^top L v) / norm(v)_2^2 = epsilon / ((n - 1) / n) = (n / (n - 1)) epsilon. $
  By contrast, the balanced partition separating $K_m^((1)) union {v_0}$ from $K_m^((2))$ cuts $b$ edges and achieves $"Rcut" = b (1 / (m + 1) + 1 / m) = (2 m + 1) b / (m (m + 1))$. Meanwhile, the combinatorial Laplacian of the core $G_"core" = K_m^((1)) union K_m^((2))$ has strictly positive Fiedler value $lambda_2 (L_"core") >= 2 b / m$ (for $b << m^2$). Writing $L(G_n) = L_0 + epsilon E_0$ where $L_0$ disconnects $v_0$ from $G_"core"$, $L_0$ has a two-dimensional kernel spanned by $bb(1)$ and $v$, separated from $lambda_3 (L_0) = lambda_2 (L_"core") > 0$. By analytic perturbation of the isolated eigenspace as $epsilon -> 0^+$, $u_2 (epsilon) = v / norm(v)_2 + O(epsilon / lambda_2 (L_"core"))$, proving complete localization onto the single outlier $v_0$.
]

This finite-sample mechanism directly mirrors von Luxburg's asymptotic spectral consistency theorem (@vonluxburg2007tutorial): as $n -> infinity$, eigenvectors of unnormalized $L / n$ converge to smooth cluster eigenfunctions *only if* the corresponding eigenvalues lie strictly below the minimum degree function $inf_x d(x)$; whenever sparse tails depress $d_"min"$ below the inter-cluster cut value, unnormalized spectral clustering becomes asymptotically inconsistent and outputs localized outlier spikes.

#subsec("D", "Regular-Graph Equivalence and Distortion Envelope")

#proposition(title: "Degree-Homogeneity Envelope Relating Rcut and Ncut", origin: "Original Proof")[
  Let $G = (V, E, W)$ have vertex degrees satisfying $0 < d_"min" <= d_i <= d_"max"$ for all $i in V$. Then for *every* $k$-partition $cal(C) = (C_1, dots, C_k)$:
  $ d_"min" dot "Ncut"(cal(C)) <= "Rcut"(cal(C)) <= d_"max" dot "Ncut"(cal(C)). $
  In particular, if $G$ is $d$-regular ($d_"min" = d_"max" = d$), then $"Rcut"(cal(C)) = d dot "Ncut"(cal(C))$ identically, $L = d cal(L)_"rw" = d cal(L)_"sym"$, and Ratio Cut, Normalized Cuts (Alg. 104), and un-row-normalized NJW share the exact same eigenspaces and optimal partitions.
]

#proof[
  For any nonempty cluster $C_m subset.eq V$, summing $d_"min" <= d_i <= d_"max"$ over $i in C_m$ gives $d_"min" |C_m| <= "vol"(C_m) <= d_"max" |C_m|$. Taking reciprocals and multiplying by $"cut"(C_m, overline(C)_m) >= 0$ yields
  $ d_"min" ("cut"(C_m, overline(C)_m)) / ("vol"(C_m)) <= ("cut"(C_m, overline(C)_m)) / |C_m| <= d_"max" ("cut"(C_m, overline(C)_m)) / ("vol"(C_m)). $
  Summing over $m = 1, dots, k$ establishes the two-sided inequality for every partition $cal(C)$. When $d_i = d$ for all $i$, $D = d I_n$, so $L = D - W = d (I_n - D^(-1) W) = d cal(L)_"rw" = d cal(L)_"sym"$, making both discrete objectives and continuous spectral relaxations identical up to the positive scalar $d$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because Ratio Cut operates directly on the symmetric positive semidefinite matrix $L = D - W$ without degree-root inversions, exact dense factorization takes $O(n^3)$ time and $O(n^2)$ space (slightly simpler constant factors than Ncut's $D^(-1/2)$ conjugations, though vulnerable to wider condition numbers $lambda_n (L) / lambda_2 (L) <= 2 d_"max" / lambda_2 (L)$ on heavy-tailed graphs). On sparse graphs with $m = |E|$ edges, Lanczos or LOBPCG extracts the bottom-$k$ eigenpairs in $O(T_"iter" m k + n k^2)$ time and $O(m + n k)$ memory, followed by $O(m + n log n)$ exact 2-way sweep or $O(I n k^2)$ $K$-means rounding.

Geometrically, Ratio Cut favors equal-cardinality partitions on regular or near-regular grids, meshes, and $q$-regular nearest-neighbor graphs where $d_"max" / d_"min" approx 1$. Conversely, whenever degree heterogeneity $d_"max" / d_"min"$ is large, Proposition 1's distortion envelope widens and Theorem 3's outlier localization severely degrades cluster purity.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Dense / Sparse Complexity], [$O(n^3)$ dense; $O(T_"iter" m k + m + n log n)$ sparse Lanczos + sweep.],
    [Memory Footprint], [$O(n^2)$ dense $L$; $O(m + n k)$ sparse CSR + $k$ eigenvectors.],
    [Balancing Criterion], [Vertex cardinality $|C_m|$ (equal-size cluster bias).],
    [Cheeger Guarantee], [$lambda_2 <= "Rcut" <= 2 sqrt(2 d_"max" lambda_2)$ (Thm. 2).],
    [Regular-Graph Limit], [Exact equivalence to Ncut when $d_"min" = d_"max"$ (Prop. 1).],
    [Primary Vulnerability], [Low-degree outlier spike localization (Thm. 3).],
  ),
  caption: [Operational, complexity, and theoretical profile of Ratio Cut.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *VLSI Circuit Netlist Partitioning and Floorplanning:* Hagen and Kahng (@hagen1992ratiocut) introduced Ratio Cut specifically for VLSI netlist decomposition, where vertices represent logic gates of comparable silicon area and edges represent interconnect wires. Minimizing $"cut"(A, B) / (|A| |B|)$ balances silicon chip area ($|A|, |B|$) directly while minimizing inter-block wire crossings, avoiding artificial 50/50 cuts that break cohesive functional blocks.
2. *Parallel Scientific Computing and Finite-Element Mesh Decomposition:* In domain decomposition for PDE solvers (Chaco @hendrickson1995chaco, METIS @karypis1998metis), mesh elements have uniform computational cost and bounded, nearly regular stencil degree ($d_"max" / d_"min" approx 1$). Ratio Cut directly balances per-processor element counts $|C_m|$ while minimizing MPI halo communication $"cut"(C_m, overline(C)_m)$, achieving optimal parallel load balance without Proposition 1's degree-skew penalty.
3. *Load-Balanced Microservice and Distributed Graph Sharding:* When partitioning homogeneous state tables or uniform-cost actor graphs across $k$ compute shards, memory footprint scales with vertex count $|C_m|$ while network RPC traffic scales with cross-shard edge weight $"cut"(C_m, overline(C)_m)$. Ratio Cut directly optimizes the communication-to-memory ratio per shard.
4. *Regular $k$-Nearest-Neighbor Manifold Partitioning:* When data are converted to a mutual or symmetrized $q$-NN graph on an approximately uniformly sampled manifold, degrees concentrate tightly around $q$. By Proposition 1, Ratio Cut preserves non-convex manifold separation while directly discouraging degenerate tiny clusters.
5. *Retail Customer Base Segmentation with Balanced Cohort Targets:* In commercial marketing operations where campaign teams require roughly comparable customer cohort sizes $|C_m|$ for A/B capacity allocation, Ratio Cut's cardinality denominator naturally yields balanced customer segments without distorting dense commercial similarity graphs.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Ratio Cut fails predictably in four regimes: (i) *heavy-tailed / scale-free networks* ($d_"max" >> d_"min"$), where Theorem 2's $sqrt(d_"max")$ slack and Theorem 3's low-degree spike localization cause Fiedler vectors to isolate peripheral pendants; (ii) *variable-density Euclidean point clouds* (heteroscedastic clusters), where sparse cluster fringes have minimal degree $d_i$ and violate von Luxburg's consistency condition (@vonluxburg2007tutorial); (iii) *severely imbalanced true cluster sizes* ($|C_1| << |C_2|$), where the $1 / |C_m|$ denominator heavily penalizes genuine small clusters; and (iv) *large-scale dense graphs* ($n > 10^4$), where multilevel coarsening (METIS @karypis1998metis) outperforms exact spectral $O(n^3)$ solvers.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [VLSI / FE mesh partitioning], [Excellent: regular degrees + area balance.], [Use METIS (@karypis1998metis) for $n > 10^5$.],
    [Balanced-size customer cohorts], [High: $|C_m|$ denominator enforces balance.], [Check degree outliers first.],
    [Scale-free / social networks], [Poor: leaf spikes hijack $u_2$ (Thm. 3).], [Use Ncut (Alg. 104) or Leiden (@traag2019leiden).],
    [Heterogeneous-density data], [Low: fails $lambda_k < d_"min"$ consistency.], [Use Ncut (Alg. 104) or HDBSCAN (@campello2013hdbscan).],
    [Imbalanced rare-class discovery], [Low: $1 / |C_m|$ suppresses small classes.], [Use density or AP (Alg. 106).],
  ),
  caption: [Domain suitability and failure-mode matrix for Ratio Cut.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Mirroring Ncut's connection to *weighted* Kernel $K$-means (Alg. 104), Ratio Cut reduces exactly to *unweighted* Kernel $K$-means via a simple spectral shift of the unnormalized Laplacian (@dhillon2004kernel).

#proposition(title: "Equivalence to Unweighted Kernel K-Means", origin: "Classical — Cited")[
  For any constant $sigma > 2 d_"max"$, let $K = sigma I_n - L = sigma I_n - D + W$. Then $K$ is symmetric positive definite, and minimizing the $k$-way Ratio Cut objective $"Rcut"(C_1, dots, C_k)$ over partitions of $G$ is strictly equivalent to minimizing the unweighted Kernel $K$-means objective $J_"KKM"(cal(C)) = sum_(m=1)^k sum_(i in C_m) norm(phi(i) - mu_m)_2^2$ in the RKHS induced by $K$ (@dhillon2004kernel).
]

#proof[
  By Gershgorin's circle theorem, every eigenvalue of $L = D - W$ satisfies $0 <= lambda_i (L) <= 2 d_"max"$, so $K = sigma I_n - L succ 0$ whenever $sigma > 2 d_"max"$. With cardinality-normalized indicators $H_(i m) = |C_m|^(-1/2) bb(I)(i in C_m)$ ($H^top H = I_k$), the unweighted Kernel $K$-means objective satisfies
  $ J_"KKM"(cal(C))
    &= "Tr"(K) - "Tr"(H^top K H) \
    &= "Tr"(K) - sigma "Tr"(H^top H) + "Tr"(H^top L H) \
    &= "Tr"(K) - sigma k + "Rcut"(cal(C)), $
  using Theorem 1(i). Because $"Tr"(K) - sigma k$ is independent of the partition $cal(C)$, the discrete minimizers of $J_"KKM"$ and $"Rcut"$ coincide identically.
]

Consequently, the progression across Algorithms 103--105 forms a unified hierarchy: Ratio Cut (`Alg. 105`) solves unweighted Kernel $K$-means on $sigma I_n - L$; Normalized Cuts (`Alg. 104`) upgrades node weights from $1$ to $d_i$ to eliminate degree-outlier localization; and NJW (`Alg. 103`) additionally normalizes spectral rows onto $bb(S)^(k-1)$ before Euclidean $K$-means.

*Empirical Benchmark Analysis:* The 4-dataset ClusterOverview benchmark starkly confirms both the strengths and theoretical vulnerabilities of Ratio Cut. On *Customers* (`Sil = 0.356`), Ratio Cut achieves the *highest* Silhouette score among all four algorithms in this group (surpassing Ncut `0.269`, Affinity Propagation `0.261`, and NJW `0.008`), illustrating its effectiveness when balanced cohort cardinality aligns with moderate degree variation. On *High-D* (`Sil = 0.219`, `ARI = N/A`), it remains competitive behind Ncut (`0.293`) and well ahead of NJW (`-0.028`). Conversely, on *Iris* (`Sil = 0.118`, `ARI = 0.114`) and *Multi-Shape* (`Sil = 0.390`, `ARI = N/A`), Ratio Cut trails Ncut (`Iris ARI = 0.414`, `Multi-Shape Sil = 0.622`) and Affinity Propagation (`Iris ARI = 0.449`, `Multi-Shape Sil = 0.727`) substantially---a direct empirical manifestation of Theorem 3, where peripheral low-degree boundary points in non-uniformly sampled shapes distort the unnormalized Fiedler subspace. Its mean runtime of `4.568 ms` closely tracks NJW (`4.226 ms`) and Ncut (`3.882 ms`), reflecting shared $O(n^3)$ dense eigendecomposition cost.
