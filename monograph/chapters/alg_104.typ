#import "../template/ieee.typ": *

#alg_entry(
  id: 104,
  name: "Normalized Cuts (Ncut)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2000",
  paradigm: "Normalized Graph Cut Relaxation",
  complexity: $O(n^3)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical minimum-cut graph partitioning minimizes $"cut"(A, B) = sum_(i in A, j in B) W_(i j)$, which exhibits a notorious structural bias toward severing isolated singleton vertices or tiny low-degree cliques whose boundary weight is trivially small. To eliminate this degenerate solution while capturing both *between-cluster dissimilarity* and *within-cluster cohesion*, Shi and Malik (@shi2000ncut) introduced the *Normalized Cut* (Ncut) criterion, which normalizes the cut weight by the total edge volume of each induced subgraph.

Let $G = (V, E, W)$ be a connected, weighted undirected graph on $V = {1, dots, n}$ with symmetric nonnegative affinity matrix $W in bb(R)_(>= 0)^(n times n)$, vertex degrees $d_i = sum_(j=1)^n W_(i j) > 0$, degree matrix $D = "diag"(d_1, dots, d_n)$, and combinatorial Laplacian $L = D - W$. For any subsets $A, B subset.eq V$, define the edge association $"assoc"(A, B) = sum_(i in A, j in B) W_(i j)$, the cut $"cut"(A, B) = "assoc"(A, B)$ for disjoint $A, B$, and the graph volume $"vol"(A) = "assoc"(A, V) = sum_(i in A) d_i$.

#definition(title: "k-Way Normalized Cut and Normalized Association")[
  For a disjoint $k$-partition $cal(C) = (C_1, dots, C_k)$ of $V$ with $"vol"(C_m) > 0$ for all $m$, the *$k$-way Normalized Cut* and *Normalized Association* are
  $ "Ncut"(C_1, dots, C_k) &= sum_(m=1)^k ("cut"(C_m, overline(C)_m)) / ("vol"(C_m)), \
    "Nassoc"(C_1, dots, C_k) &= sum_(m=1)^k ("assoc"(C_m, C_m)) / ("vol"(C_m)), $
  where $overline(C)_m = V without C_m$. Minimizing $"Ncut"(cal(C))$ over all $k$-partitions is equivalent to solving a discrete constrained trace program whose continuous spectral relaxation yields the generalized eigenproblem $(D - W) y = lambda D y$, equivalently $cal(L)_"rw" y = lambda y$ for the random-walk Laplacian $cal(L)_"rw" = D^(-1) (D - W)$.
]

In the two-way case ($k = 2$, $cal(C) = (A, B)$), Shi and Malik compute the Fiedler eigenvector $y_2$ corresponding to the second smallest generalized eigenvalue $lambda_2$ of $(D - W) y = lambda D y$, sweep across threshold splits $A(tau) = {i in V : y_(2,i) > tau}$ to minimize exact discrete $"Ncut"(A(tau), overline(A)(tau))$, and recurse on subgraphs until a target $k$ or cut threshold is reached. For simultaneous $k$-way partitioning, the bottom-$k$ generalized eigenvectors $Y = [y_1, dots, y_k]$ are discretized via $K$-means or optimal orthogonal rotation (@shi2000ncut, @vonluxburg2007tutorial).

#algbox(
  title: "Normalized Cuts (Shi-Malik Spectral Partitioning)",
  inputs: [Graph $G = (V, E, W)$ (or data $X$ + affinity), target $k$, mode (recursive 2-way or simultaneous $k$-way).],
  outputs: [Partition $cal(C) = (C_1, dots, C_k)$ approximating $min "Ncut"(cal(C))$.],
)[
  1. *Laplacian Assembly:* Compute $d_i = sum_(j=1)^n W_(i j)$, $D = "diag"(d_1, dots, d_n)$, and $L = D - W$.
  2. *Generalized Spectral Solve:* Transform $(D - W) y = lambda D y$ into the symmetric problem $D^(-1/2) L D^(-1/2) z = lambda z$, compute the bottom-$k$ eigenvectors $z_1, dots, z_k$, and back-transform $y_m = D^(-1/2) z_m$.
  3. *2-Way Sweep (Recursive Mode):* Sort vertices by $y_(2,i)$, evaluate exact discrete $"Ncut"(A_r, overline(A)_r)$ across all $n - 1$ prefix cuts $A_r$ in $O(m + n log n)$ time, split at $r^* = arg min_r "Ncut"(A_r, overline(A)_r)$, and recurse on induced subgraphs until $k$ clusters are formed.
  4. *Simultaneous $k$-Way Rounding:* Assemble $Y = [y_1, dots, y_k] in bb(R)^(n times k)$ and cluster its rows via $K$-means (or iterative orthogonal discretization) to obtain $C_1, dots, C_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Rayleigh-Quotient Formulation and Duality")

#theorem(title: "Exact Trace Identity and Ncut-Nassoc Duality", origin: "Original Proof")[
  Let $cal(C) = (C_1, dots, C_k)$ be any partition of $V$ with $"vol"(C_m) > 0$. Define the volume-normalized indicator matrix $H in bb(R)^(n times k)$ by $H_(i m) = "vol"(C_m)^(-1/2) bb(I)(i in C_m)$. Then:
  (i) $H^top D H = I_k$ and $"Tr"(H^top L H) = "Ncut"(C_1, dots, C_k) = k - "Nassoc"(C_1, dots, C_k)$;
  (ii) for $k = 2$ with $cal(C) = (A, B)$ and $b = "vol"(A) / "vol"(B)$, the centered indicator $y_i = 1$ ($i in A$) and $y_i = -b$ ($i in B$) satisfies $y^top D bb(1) = 0$ and the exact Rayleigh quotient identity
  $ (y^top (D - W) y) / (y^top D y) = "Ncut"(A, B). $
]

#proof[
  For (i), the $(m, m')$ entry of $H^top D H$ is
  $ (H^top D H)_(m m') = sum_(i=1)^n d_i H_(i m) H_(i m') = delta_(m m') / ("vol"(C_m)) sum_(i in C_m) d_i = delta_(m m'), $
  so $H^top D H = I_k$. Letting $h_m$ be the $m$-th column of $H$, its Laplacian quadratic form is
  $ h_m^top L h_m
    &= 1/2 sum_(i,j=1)^n W_(i j) (H_(i m) - H_(j m))^2 \
    &= 1 / ("vol"(C_m)) sum_(i in C_m, j in overline(C)_m) W_(i j) = ("cut"(C_m, overline(C)_m)) / ("vol"(C_m)), $
  because $(H_(i m) - H_(j m))^2 = "vol"(C_m)^(-1)$ iff exactly one of $i, j$ lies in $C_m$, and each undirected cut edge is counted twice in the double sum. Summing over $m = 1, dots, k$ yields $"Tr"(H^top L H) = "Ncut"(C_1, dots, C_k)$. Moreover, since $"vol"(C_m) = "assoc"(C_m, C_m) + "cut"(C_m, overline(C)_m)$, dividing by $"vol"(C_m)$ and summing over $m$ gives $"Ncut"(cal(C)) + "Nassoc"(cal(C)) = k$, proving that minimizing $"Ncut"$ simultaneously maximizes within-cluster normalized association.

  For (ii), direct summation gives $y^top D bb(1) = "vol"(A) - b "vol"(B) = 0$. The numerator is
  $ y^top L y
    &= 1/2 sum_(i,j=1)^n W_(i j) (y_i - y_j)^2 = (1 + b)^2 "cut"(A, B), $
  while the denominator is
  $ y^top D y
    &= "vol"(A) + b^2 "vol"(B) = "vol"(A) (1 + b). $
  Dividing numerator by denominator yields
  $ (y^top L y) / (y^top D y)
    &= (1 + b) ("cut"(A, B)) / ("vol"(A)) \
    &= ("cut"(A, B)) / ("vol"(A)) + ("cut"(A, B)) / ("vol"(B)) = "Ncut"(A, B). $
]

#subsec("B", "NP-Hardness and Spectral Relaxation Lower Bound")

#theorem(title: "NP-Hardness and Ky Fan Spectral Lower Bound", origin: "Classical — Cited")[
  Exact minimization of $"Ncut"(A, B)$ over partitions of a weighted graph $G$ is NP-complete, even on planar graphs (@shi2000ncut). Relaxing the discrete matrix $H$ in Theorem 1 to any real matrix $H in bb(R)^(n times k)$ satisfying $H^top D H = I_k$ yields the tight spectral lower bound
  $ min_(cal(C)) "Ncut"(C_1, dots, C_k) >= sum_(m=1)^k lambda_m (cal(L)_"rw"), $
  where $0 = lambda_1 <= lambda_2 <= dots <= lambda_n$ are the eigenvalues of $cal(L)_"rw" = D^(-1) L$. For $k = 2$, $min_(A,B) "Ncut"(A, B) >= lambda_2 (cal(L)_"rw")$, with equality if and only if $cal(L)_"rw"$ admits a two-valued Fiedler eigenvector.
]

#proof[
  Substitute $Z = D^(1/2) H in bb(R)^(n times k)$. The constraint $H^top D H = I_k$ becomes exact orthonormality $Z^top Z = I_k$, and the trace objective becomes
  $ "Tr"(H^top L H) = "Tr"(Z^top D^(-1/2) L D^(-1/2) Z) = "Tr"(Z^top cal(L)_"sym" Z). $
  By Ky Fan's extremal trace theorem, the minimum of $"Tr"(Z^top cal(L)_"sym" Z)$ over all $Z in bb(R)^(n times k)$ with $Z^top Z = I_k$ is $sum_(m=1)^k lambda_m (cal(L)_"sym")$, attained when the columns of $Z$ span the bottom-$k$ eigenspace of $cal(L)_"sym"$. Since $cal(L)_"rw" = D^(-1/2) cal(L)_"sym" D^(1/2)$ is similar to $cal(L)_"sym"$, their spectra coincide identically ($lambda_m (cal(L)_"rw") = lambda_m (cal(L)_"sym")$). Because discrete volume-normalized indicators form a strict subset of ${H : H^top D H = I_k}$, every $k$-partition satisfies $"Ncut"(cal(C)) >= sum_(m=1)^k lambda_m (cal(L)_"rw")$. For $k = 2$, $lambda_1 = 0$ corresponds to $y_1 = bb(1)$, so orthogonality $y^top D bb(1) = 0$ in Theorem 1(ii) gives $min_(A,B) "Ncut"(A, B) >= min_(y perp_D bb(1)) (y^top L y) / (y^top D y) = lambda_2 (cal(L)_"rw")$, with equality iff the discrete vector $y in {1, -b}^n$ achieves the Rayleigh minimum, i.e., is an eigenvector for $lambda_2$.
]

#subsec("C", "Discrete Cheeger Guarantee for Sweep Rounding")

The continuous relaxation is only useful if rounding the Fiedler vector $y_2$ yields a discrete cut provably close to optimal. Define the *Cheeger conductance* of $G$ as $h_G = min_(0 < "vol"(A) <= "vol"(V)/2) "cut"(A, overline(A)) / "vol"(A)$. Notice that for every bipartition $(A, overline(A))$ with $"vol"(A) <= "vol"(V)/2$, we have $"vol"(overline(A)) >= "vol"(A)$ and therefore
$ h_G (A) <= "Ncut"(A, overline(A)) = "cut"(A, overline(A)) (1 / ("vol"(A)) + 1 / ("vol"(overline(A)))) <= 2 h_G (A). $

#theorem(title: "Cheeger Conductance Bound for Fiedler Sweep", origin: "Classical — Cited")[
  Let $G$ be a connected weighted graph with Fiedler value $lambda_2 = lambda_2 (cal(L)_"rw")$ and eigenvector $y_2$ satisfying $cal(L)_"rw" y_2 = lambda_2 y_2$. Let $A^*$ be the optimal prefix sweep cut over vertices sorted by $y_2$ (@chung1997spectral, @shi2000ncut). Then
  $ lambda_2 / 2 <= h_G <= h_G (A^*) <= sqrt(2 lambda_2), $
  and consequently the sweep-rounded partition satisfies the constructive approximation guarantee
  $ lambda_2 <= "Ncut"(A^*, overline(A)^*) <= 2 sqrt(2 lambda_2). $
]

#proof[
  For the lower bound, Theorem 2 gives $lambda_2 <= "Ncut"(A, overline(A)) <= 2 h_G (A)$ for every $A$ with $"vol"(A) <= "vol"(V)/2$; taking the minimum over $A$ yields $lambda_2 / 2 <= h_G$.

  For the upper bound, shift $y_2$ by its volume-weighted median $m_"med"$ (chosen so that both ${i : y_(2,i) > m_"med"}$ and ${i : y_(2,i) < m_"med"}$ have volume at most $"vol"(V)/2$), and define the positive and negative parts $f^+ = max(y_2 - m_"med" bb(1), 0)$ and $f^- = max(m_"med" bb(1) - y_2, 0)$. Since $y_2 - m_"med" bb(1)$ is supported on disjoint positive/negative sets and $f^+, f^-$ have opposite signs across edges crossing $m_"med"$, at least one of $g in {f^+, f^-}$ satisfies $(g^top L g) / (g^top D g) <= lambda_2$ while $"vol"("supp"(g)) <= "vol"(V)/2$. Order vertices so $g_1 >= g_2 >= dots >= g_n >= 0$ and scale so $g_1 = 1$. Let $S_r = {1, dots, r}$ and $h^* = min_(r : "vol"(S_r) <= "vol"(V)/2) "cut"(S_r, overline(S)_r) / "vol"(S_r) >= h_G (A^*)$. Applying Cauchy-Schwarz to the telescoping sum of squared differences yields
  $ sum_(i < j) W_(i j) |g_i^2 - g_j^2|
    &= sum_(i < j) W_(i j) |g_i - g_j| (g_i + g_j) \
    &<= (sum_(i < j) W_(i j) (g_i - g_j)^2)^(1/2) (sum_(i < j) W_(i j) (g_i + g_j)^2)^(1/2). $
  The first factor equals $(g^top L g)^(1/2) <= (lambda_2 g^top D g)^(1/2)$, while the second obeys $(g_i + g_j)^2 <= 2(g_i^2 + g_j^2)$, so $sum_(i < j) W_(i j) (g_i + g_j)^2 <= 2 sum_i d_i g_i^2 = 2 g^top D g$. Hence the left-hand side is at most $sqrt(2 lambda_2) g^top D g$. On the other hand, expanding $|g_i^2 - g_j^2| = sum_(r=i)^(j-1) (g_r^2 - g_(r+1)^2)$ for $i < j$ gives the co-area identity
  $ sum_(i < j) W_(i j) |g_i^2 - g_j^2|
    &= sum_(r=1)^(n-1) (g_r^2 - g_(r+1)^2) "cut"(S_r, overline(S)_r) \
    &>= h^* sum_(r=1)^(n-1) (g_r^2 - g_(r+1)^2) "vol"(S_r) = h^* g^top D g. $
  Dividing by $g^top D g > 0$ gives $h_G (A^*) <= h^* <= sqrt(2 lambda_2)$, and $"Ncut"(A^*, overline(A)^*) <= 2 h_G (A^*) <= 2 sqrt(2 lambda_2)$.
]

#subsec("D", "Markov-Chain Metastability and Geometric Path Failure")

#theorem(title: "Stationary Random-Walk Escape Probability", origin: "Original Derivation")[
  Consider the reversible Markov chain on $V$ with transition matrix $P = D^(-1) W$ and stationary distribution $pi_i = d_i / "vol"(V)$. For any partition $cal(C) = (C_1, dots, C_k)$ with $"vol"(C_m) > 0$,
  $ "Ncut"(C_1, dots, C_k) = sum_(m=1)^k "Pr"_(X_0 ~ pi) (X_1 in.not C_m | X_0 in C_m). $
  Thus minimizing $k$-way Normalized Cut is mathematically identical to identifying the $k$ most metastable subsets minimizing total single-step stationary escape probability.
]

#proof[
  Under $X_0 ~ pi$, the probability of starting in $C_m$ is $"Pr"(X_0 in C_m) = sum_(i in C_m) d_i / "vol"(V) = "vol"(C_m) / "vol"(V)$. The joint probability of transitioning out of $C_m$ in one step is
  $ "Pr"(X_0 in C_m, X_1 in.not C_m)
    &= sum_(i in C_m) sum_(j in overline(C)_m) pi_i P_(i j) \
    &= sum_(i in C_m) sum_(j in overline(C)_m) (d_i / ("vol"(V))) (W_(i j) / d_i) = ("cut"(C_m, overline(C)_m)) / ("vol"(V)). $
  Dividing by $"Pr"(X_0 in C_m) = "vol"(C_m) / "vol"(V)$ cancels $"vol"(V)$ and yields $"Pr"(X_1 in.not C_m | X_0 in C_m) = "cut"(C_m, overline(C)_m) / "vol"(C_m)$. Summing over $m = 1, dots, k$ recovers $"Ncut"(C_1, dots, C_k)$ exactly.
]

#proposition(title: "Elongated-Region Over-Segmentation Bias", origin: "Original Counterexample")[
  Let $P_n$ be an unweighted path graph on $n >= 4$ vertices (modeling a homogeneous elongated image strip or filament) with unit edge weights. Then the optimal Ncut bipartition bisects the homogeneous path cleanly in half with $"Ncut"(A^*, B^*) = 2 / (n - 1) = Theta(1 / n)$, whereas isolating a degree-$1$ endpoint yields $"Ncut"({1}, V without {1}) = 1 + 1 / (2 n - 3) > 1$. Thus Ncut rigorously suppresses singleton outliers yet systematically fragments long homogeneous structures.
]

#proof[
  In $P_n$, total volume is $"vol"(V) = 2(n - 1)$, with $d_1 = d_n = 1$ and $d_i = 2$ for $2 <= i <= n - 1$. Cutting the single edge $(r, r + 1)$ yields $A_r = {1, dots, r}$ with $"cut"(A_r, overline(A)_r) = 1$, $"vol"(A_r) = 2 r - 1$, and $"vol"(overline(A)_r) = 2(n - r) - 1$. Thus
  $ "Ncut"(A_r, overline(A)_r) = 1 / (2 r - 1) + 1 / (2(n - r) - 1) = (2(n - 1)) / ((2 r - 1)(2(n - r) - 1)). $
  For $r = 1$, this equals $1 + 1 / (2 n - 3) > 1$. The denominator is strictly concave in $r$ and maximized at $r^* = floor(n / 2)$, yielding $"Ncut"(A_(r^*), overline(A)_(r^*)) = 2 / (n - 1)$ (for even $n$), which vanishes as $n -> infinity$ despite uniform internal connectivity.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact dense Ncut transforms $(D - W) y = lambda D y$ via $z = D^(1/2) y$ into the symmetric eigensystem $cal(L)_"sym" z = lambda z$, requiring $O(n^3)$ time and $O(n^2)$ memory for full dense factorization. In image segmentation and sparse network graphs with $m = |E|$ edges (e.g., local pixel neighborhoods of radius $r$ where $m = O(r^2 n)$), Lanczos iteration (@shi2000ncut) requires only sparse matrix-vector products $O(m)$ per step, reducing extraction of the top $k$ generalized eigenvectors to $O(T_"iter" m k + n k^2)$ time and $O(m + n k)$ memory, where $T_"iter" = O(gamma^(-1/2) log(1 / epsilon))$ depends on the spectral gap $gamma = lambda_(k+1) - lambda_k$. The exact 2-way sweep over $n - 1$ sorted thresholds maintains running cut and volume incrementally in $O(m + n log n)$ time by updating each vertex's incident edge counts once as it crosses the threshold.

Because the denominator uses graph volume $"vol"(C_m)$ rather than vertex count $|C_m|$, Ncut balances total connectivity mass rather than raw cardinality, remaining remarkably stable under degree heterogeneity and high-dimensional noise provided the underlying affinity graph preserves local neighborhood structure.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Dense / Sparse Complexity], [$O(n^3)$ dense; $O(T_"iter" m k + m + n log n)$ sparse Lanczos + sweep.],
    [Memory Footprint], [$O(n^2)$ dense affinity; $O(m + n k)$ sparse CSR + $k$ vectors.],
    [Optimization Objective], [Minimizes sum of stationary Markov escape probabilities (Thm. 4).],
    [Approximation Guarantee], [Fiedler sweep yields $lambda_2 <= "Ncut" <= 2 sqrt(2 lambda_2)$ (Thm. 3).],
    [Outlier Resilience], [High: singleton cut incurs $"Ncut" > 1$ (Prop. 1).],
    [Structural Bias], [Balanced-volume cuts; splits long thin regions (Prop. 1).],
  ),
  caption: [Operational, algorithmic, and theoretical profile of Normalized Cuts.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Perceptual Image Segmentation and Superpixel Grouping:* Shi and Malik (@shi2000ncut) originally designed Ncut for vision graphs where pixel affinities combine brightness/color/texture cue distance and spatial proximity $W_(i j) = exp(-norm(F_i - F_j)_2^2 / sigma_F^2 - norm(X_i - X_j)_2^2 / sigma_X^2)$. Volume normalization prevents noisy isolated pixels from forming trivial cuts while enforcing global perceptual closure across weak local edges; modern pipelines pair SLIC superpixels (@achanta2012slic) with Ncut to cut runtime by $10^3 times$.
2. *Connectomics and Functional Brain Parcellation:* Resting-state fMRI voxel correlation graphs exhibit strong spatial autocorrelation and heterogeneous hub degrees across cortical regions. Minimizing Ncut identifies functionally coherent parcels whose internal correlation volume dominates cross-parcel boundary flux, directly modeling neural information containment via Theorem 4.
3. *Gene Co-Expression and Protein Interaction Modules:* Biological interaction networks have heavy-tailed degree distributions with essential hub proteins and sparsely connected peripheral genes. Unlike Ratio Cut (Alg. 105), which is easily hijacked by low-degree peripheral outliers, Ncut weights vertices by degree $d_i$, grouping hub-anchored pathways cleanly without fragmenting leaf vertices.
4. *Spatial Mobility and Urban Traffic Partitioning:* In origin-destination flow graphs where edge weight $W_(i j)$ measures commuter trips between zones $i$ and $j$, Theorem 4 gives Ncut an exact operational meaning: it partitions a metropolitan area into $k$ self-contained mobility basins minimizing the total conditional fraction of trips crossing basin boundaries.
5. *Document-Term Bipartite Co-Clustering:* Spectral co-clustering (@kluger2003spectral) applies Ncut to the bipartite document-word matrix $B$ by forming $W = mat(0, B; B^top, 0)$. The normalized cut simultaneously partitions documents and words into matched topical blocks whose SVD singular vectors solve the generalized Ncut eigensystem directly.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four practical limitations govern Ncut deployment: (i) *elongated-strip over-segmentation* (Proposition 1), where long thin objects (roads, blood vessels, filaments) yield small $"Ncut"$ bisections across their waist; (ii) *volume-imbalance sensitivity*, where a small genuine cluster with tiny volume surrounded by a massive background can have worse $"Ncut"$ than splitting the background in half; (iii) *greedy recursion myopia*, where recursive 2-way Ncut commits early cuts that cannot be undone globally (mitigated by simultaneous $k$-way rounding); and (iv) *spectral convergence slowdown* on poorly separated graphs with small $lambda_(k+1) - lambda_k$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Image / superpixel vision], [Excellent: global perceptual closure.], [Pre-segment via SLIC (@achanta2012slic).],
    [Heterogeneous-degree graphs], [High: volume normalization resists leaves.], [Watch tiny minority communities.],
    [Origin-destination flow maps], [Exact Markov escape-rate objective (Thm. 4).], [Direct asymmetries need directed cut.],
    [Long filaments / vessels], [Low: bisects homogeneous paths (Prop. 1).], [Use DBSCAN (@ester1996dbscan) or HDBSCAN.],
    [Large sparse graphs ($n > 10^6$)], [Slow Lanczos if $gamma$ is tiny.], [Use METIS (@karypis1998metis) or Leiden (@traag2019leiden).],
  ),
  caption: [Domain suitability and failure-mode matrix for Normalized Cuts.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Normalized Cuts forms an exact mathematical bridge between graph partitioning and kernel prototype clustering, established rigorously by Dhillon, Guan, and Kulis (@dhillon2004kernel).

#proposition(title: "Algebraic Equivalence to Weighted Kernel K-Means", origin: "Classical — Cited")[
  Minimizing the $k$-way Normalized Cut objective $"Ncut"(C_1, dots, C_k)$ over partitions of $G = (V, E, W)$ is strictly equivalent to minimizing the weighted Kernel $K$-means objective $sum_(m=1)^k sum_(i in C_m) w_i norm(phi(i) - mu_m)_2^2$ with node weights $w_i = d_i$ and positive-definite kernel matrix $K = sigma D^(-1) + D^(-1) W D^(-1)$ (for any shift $sigma > 0$ making $K succ 0$) (@dhillon2004kernel).
]

#proof[
  For weights $w_i = d_i$ and cluster weight $W_m = sum_(i in C_m) d_i = "vol"(C_m)$, the weighted Kernel $K$-means objective with $mu_m = W_m^(-1) sum_(i in C_m) d_i phi(i)$ expands as
  $ J_"WKKM"(cal(C))
    &= sum_(i=1)^n d_i K_(i i) - sum_(m=1)^k 1 / ("vol"(C_m)) sum_(i,j in C_m) d_i d_j K_(i j) \
    &= "Tr"(D K) - "Tr"(H^top D K D H), $
  where $H_(i m) = "vol"(C_m)^(-1/2) bb(I)(i in C_m)$ satisfies $H^top D H = I_k$. Substituting $K = sigma D^(-1) + D^(-1) W D^(-1)$ gives $D K D = sigma D + W$, so
  $ "Tr"(H^top D K D H) = sigma "Tr"(H^top D H) + "Tr"(H^top W H) = sigma k + "Nassoc"(cal(C)). $
  Since $"Tr"(D K)$ and $sigma k$ are partition-independent constants and $"Nassoc"(cal(C)) = k - "Ncut"(cal(C))$ (Theorem 1), minimizing $J_"WKKM"$ is strictly equivalent to minimizing $"Ncut"(cal(C))$.
]

Relative to its immediate catalog neighbors, Ncut (`Alg. 104`) replaces Ratio Cut's (`Alg. 105`) cardinality denominator $|C_m|$ with volume $"vol"(C_m)$, eliminating the unnormalized $sqrt(d_"max")$ Cheeger degradation, and retains $D^(-1/2)$-weighted coordinates rather than applying NJW's (`Alg. 103`) unit-sphere row normalization.

*Empirical Benchmark Analysis:* On the 4-dataset ClusterOverview benchmark, Normalized Cuts demonstrates the strongest all-around robustness among the three spectral algorithms (`Algs. 103–105`). It attains `Sil = 0.622` (`ARI = N/A`) on *Multi-Shape*, surpassing both NJW (`0.592`) and Ratio Cut (`0.390`) by accurately tracing non-convex connectivity while balancing subgraph volume. On *Iris* (`Sil = 0.382`, `ARI = 0.414`), Ncut outperforms NJW (`Sil = 0.132`, `ARI = 0.296`) and Ratio Cut (`Sil = 0.118`, `ARI = 0.114`) by a wide margin, confirming that keeping $D^(-1/2)$-weighted coordinates prevents low-degree overlap points from destabilizing the partition. Crucially, Ncut also sustains strong separation on *Customers* (`Sil = 0.269`) and leads all four algorithms in this group on *High-D* (`Sil = 0.293`, `ARI = N/A`, versus NJW `-0.028`, Ratio Cut `0.219`, and Affinity Propagation `0.169`). At a mean runtime of `3.882 ms`, Ncut is also the fastest of the four methods, reflecting compact generalized eigenproblem execution and clean convergence on well-conditioned graph spectra.
