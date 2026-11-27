#import "../template/ieee.typ": *

#alg_entry(
  id: 116,
  name: "Leading Eigenvector (Newman)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2006",
  paradigm: "Spectral Divisive Modularity Maximization",
  complexity: $cal(O)(n (m + n))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $G = (V, E, w)$ be an undirected weighted graph with $n = |V|$ vertices, $m = |E|$ edges, symmetric nonnegative adjacency matrix $A in bb(R)^(n times n)$, vertex degrees $k_i = sum_(j=1)^n A_(i j)$, and total weight $m = 1 / 2 sum_(i=1)^n k_i$ (coinciding with $|E|$ for unweighted graphs). Newman @newman2006leading reformulated modularity maximization @newman2004finding as an indefinite spectral quadratic optimization problem by introducing the symmetric *modularity matrix*
$ B = A - 1 / (2 m) k k^T in bb(R)^(n times n), $
whose entries $B_(i j) = A_(i j) - k_i k_j / (2 m)$ measure excess edge weight over the configuration-model expectation. Because $sum_(j=1)^n A_(i j) = k_i$ and $sum_(j=1)^n k_j = 2 m$, every row and column of $B$ sums identically to zero:
$ B bold(1) = k - (k^T bold(1)) / (2 m) k = bold(0). $

Consider first partitioning $V$ into two communities encoded by an Ising spin vector $s in {-1, +1}^n$. Since the Kronecker delta satisfies $delta(c_i, c_j) = 1 / 2 (s_i s_j + 1)$ and $bold(1)^T B bold(1) = 0$, the modularity of the bipartition reduces to the pure quadratic form
$ Q(s) = 1 / (2 m) sum_(i, j) B_(i j) (s_i s_j + 1) / 2 = 1 / (4 m) s^T B s. $
Unlike the graph Laplacian $L = D - A$ (which is positive semidefinite and minimizes cuts via its second-smallest eigenvector @hagen1992ratiocut @shi2000ncut), $B$ is generally indefinite, and modularity maximization seeks its *largest positive* eigenvalues.

To partition $G$ into more than two communities, Newman @newman2006leading recursively bisects existing communities $g subset.eq V$ ($n_g = |g|$). Crucially, deleting inter-community edges and recomputing modularity on $G[g]$ would alter the null model $2 m$ and violate global modularity additivity. Instead, one evaluates the exact change $Delta Q$ in *global* modularity caused by splitting $g$ into $g_+ = {i in g : s_i = +1}$ and $g_- = {i in g : s_i = -1}$.

#definition(title: "Generalized Subgraph Modularity Matrix")[
  For any non-empty vertex subset $g subset.eq V$ of size $n_g = |g|$, the *generalized modularity matrix* $B^((g)) in bb(R)^(n_g times n_g)$ has entries for $i, j in g$ defined by
  $ B^((g))_(i j) = B_(i j) - delta_(i j) sum_(l in g) B_(i l). $
  Its symmetric spectral decomposition is $B^((g)) = sum_(r=1)^(n_g) lambda_r u_r u_r^T$ with eigenvalues ordered $lambda_1 >= lambda_2 >= dots >= lambda_(n_g)$ and orthonormal eigenvectors $u_1, dots, u_(n_g) in bb(R)^(n_g)$.
]

Relaxing the discrete constraint $s in {-1, +1}^(n_g)$ to the Euclidean sphere $||s||_2^2 = n_g$ maximizes $s^T B^((g)) s$ along the leading eigenvector $u_1$ associated with $lambda_1 (B^((g)))$. Rounding $s_i = "sign"(u_(1, i))$ maps the continuous relaxation back to a bipartition, optionally followed by Kernighan-Lin (KL) greedy vertex-swap fine-tuning @kernighan1970efficient.

#algbox(
  title: "Newman Leading Eigenvector Divisive Algorithm",
  inputs: [Undirected graph $G = (V, E, w)$ with modularity operator $B = A - (2 m)^(-1) k k^T$; tolerance $epsilon > 0$.],
  outputs: [Community partition $hat(cal(C)) = {C_1, dots, C_K}$ maximizing $Q$ via recursive spectral bisection.],
)[
  1. *Initialize queue:* Set active split queue $cal(Q) <- {V}$ and final partition $hat(cal(C)) <- emptyset$.
  2. *Divisive loop:* While $cal(Q) != emptyset$, pop subgraph $g$ ($n_g = |g|$):
     - If $n_g <= 1$, append $g$ to $hat(cal(C))$ and continue.
     - *Implicit matvec:* Compute the leading eigenpair $(lambda_1, u_1)$ of $B^((g))$ via Lanczos or shifted power iteration using the $cal(O)(|E_g| + n_g)$ matrix-vector action
       $ (B^((g)) v)_i = (A_g v)_i - (k_g^T v) / (2 m) k_i - d^((g))_i v_i, $
       where $d^((g))_i = sum_(l in g) B_(i l) = k_i^((g)) - k_i "vol"(g) / (2 m)$.
     - *Spectral indivisibility check:* If $lambda_1 <= epsilon$, append $g$ to $hat(cal(C))$ and continue.
     - *Sign rounding:* Set $s_i <- +1$ if $u_(1, i) >= 0$ and $s_i <- -1$ otherwise.
     - *Optional KL refinement:* Repeatedly select the unmoved vertex maximizing $Delta Q$ upon flipping $s_i -> -s_i$, keeping the best prefix state $s^*$.
     - *Accept or halt:* Compute $Delta Q = (4 m)^(-1) (s^*)^T B^((g)) s^*$. If $Delta Q > epsilon$ and both $g_+, g_-$ are non-empty, push $g_+$ and $g_-$ onto $cal(Q)$; otherwise append $g$ to $hat(cal(C))$.
  3. *Return:* $hat(cal(C))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Universal K-Way Spectral Bound and Indivisibility Certificate")

Newman @newman2006leading established that $lambda_1 (B^((g))) <= 0$ prevents any modularity-improving *bipartition* of $g$. We prove a strictly stronger universal result: the exact same operator $B^((g))$ governs arbitrary $K$-way refinements of $g$, certifying global indivisibility across *all* $K >= 2$.

#theorem(title: "Universal K-Way Spectral Modularity Bound and Indivisibility Certificate", origin: "Original Proof")[
  For any subgraph $g subset.eq V$ ($n_g = |g|$) and any $K$-way partition $g = g_1 union dots union g_K$ ($K >= 2$) with indicator vectors $z_r = bold(1)_(g_r) in {0, 1}^(n_g)$, the exact change in global modularity is
  $ Delta Q(g -> {g_r}_(r=1)^K) = 1 / (2 m) sum_(r=1)^K z_r^T B^((g)) z_r. $
  Moreover, $B^((g)) bold(1)_g = bold(0)$, and the modularity gain satisfies the universal upper bound
  $ Delta Q(g -> {g_r}_(r=1)^K) &<= (n_g lambda_1 (B^((g)))) / (2 m) (1 - sum_(r=1)^K (|g_r|^2) / n_g^2) \
                                &<= (n_g (K - 1)) / (2 m K) lambda_1 (B^((g))). $
  In particular, if $lambda_1 (B^((g))) <= 0$, then $g$ is universally indivisible: no $K$-way split for any $K >= 2$ can strictly increase global modularity.
]

#proof[
  Splitting $g$ into disjoint parts $g_1, dots, g_K$ alters only pairs $(i, j) in g times g$, yielding
  $ Delta Q = 1 / (2 m) (sum_(r=1)^K sum_(i, j in g_r) B_(i j) - sum_(i, j in g) B_(i j)). $
  Expanding the quadratic sum $sum_(r=1)^K z_r^T B^((g)) z_r$ via the definition $B^((g))_(i j) = B_(i j) - delta_(i j) sum_(l in g) B_(i l)$ gives
  $ sum_(r=1)^K z_r^T B^((g)) z_r &= sum_(r=1)^K sum_(i, j in g_r) B_(i j) - sum_(r=1)^K sum_(i in g_r) sum_(l in g) B_(i l) \
                                  &= sum_(r=1)^K sum_(i, j in g_r) B_(i j) - sum_(i, l in g) B_(i l) \
                                  &= 2 m Delta Q, $
  since ${g_r}$ partitions $g$. Checking the row sum of $B^((g))$ gives
  $ (B^((g)) bold(1)_g)_i = sum_(j in g) B_(i j) - sum_(l in g) B_(i l) = 0, $
  so $bold(1)_g in ker(B^((g)))$. For each $r$, orthogonally decompose $z_r = (|g_r| / n_g) bold(1)_g + z_(r, perp)$ with $z_(r, perp) perp bold(1)_g$. By Pythagoras,
  $ ||z_(r, perp)||_2^2 = ||z_r||_2^2 - n_g (|g_r| / n_g)^2 = |g_r| (1 - |g_r| / n_g). $
  Because $B^((g))$ is symmetric and $B^((g)) bold(1)_g = bold(0)$, the Rayleigh-Ritz theorem yields
  $ z_r^T B^((g)) z_r = z_(r, perp)^T B^((g)) z_(r, perp) <= lambda_1 (B^((g))) |g_r| (1 - |g_r| / n_g). $
  Summing over $r = 1, dots, K$ using $sum_r |g_r| = n_g$ yields the first inequality. Applying Cauchy-Schwarz, $sum_(r=1)^K |g_r|^2 >= n_g^2 / K$, gives $1 - sum_r |g_r|^2 / n_g^2 <= (K - 1) / K$, establishing the second inequality (which reduces to $(n_g / (4 m)) lambda_1$ at $K = 2$). Finally, if $lambda_1 (B^((g))) <= 0$, then $B^((g)) prec.eq 0$, so $Delta Q <= 0$ for every $K >= 2$ and every refinement of $g$.
]

#subsec("B", "Sign-Rounding Spectral Lower Bound and Localization")

While $lambda_1 (B^((g))) > 0$ is necessary for a positive split, rounding the continuous eigenvector $u_1$ to $s = "sign"(u_1)$ projects energy onto negative-eigenvalue modes. We derive an exact Parseval lower bound characterizing when sign rounding is guaranteed to succeed—and when eigenvector localization breaks it.

#lemma(title: "Sign-Rounding Lower Bound and Delocalization Threshold", origin: "Original Derivation")[
  Let $B^((g))$ have leading eigenpair $(lambda_1, u_1)$ with $lambda_1 > 0$, $||u_1||_2 = 1$, and smallest eigenvalue $lambda_("min") = lambda_(n_g) (B^((g))) <= 0$. For the sign-rounded vector $s = "sign"(u_1) in {-1, +1}^(n_g)$,
  $ Delta Q(s) >= 1 / (4 m) [ lambda_1 ||u_1||_1^2 - |lambda_("min")| (n_g - ||u_1||_1^2) ]. $
  Consequently, sign rounding is guaranteed to achieve $Delta Q(s) > 0$ without KL refinement whenever the inverse participation ratio satisfies
  $ (||u_1||_1^2) / n_g > (|lambda_("min")|) / (lambda_1 + |lambda_("min")|). $
]

#proof[
  Expand $s$ in the orthonormal eigenbasis ${u_1, dots, u_(n_g)}$ of $B^((g))$ as $s = sum_(r=1)^(n_g) c_r u_r$, where $c_r = chevron.l u_r, s chevron.r$. Since $s_i = "sign"(u_(1, i))$ (with $"sign"(0) = +1$), the leading coefficient is exact:
  $ c_1 = sum_(i in g) s_i u_(1, i) = sum_(i in g) |u_(1, i)| = ||u_1||_1. $
  By Parseval's identity, $sum_(r=1)^(n_g) c_r^2 = ||s||_2^2 = n_g$, whence $sum_(r=2)^(n_g) c_r^2 = n_g - ||u_1||_1^2 >= 0$ (by Cauchy-Schwarz, $||u_1||_1 <= sqrt(n_g) ||u_1||_2 = sqrt(n_g)$). Evaluating the quadratic form and lower-bounding $lambda_r >= lambda_("min") = -|lambda_("min")|$ (noting $lambda_("min") <= 0$ because $0$ is an eigenvalue of $B^((g))$) gives
  $ s^T B^((g)) s &= lambda_1 c_1^2 + sum_(r=2)^(n_g) lambda_r c_r^2 \
                  &>= lambda_1 ||u_1||_1^2 - |lambda_("min")| (n_g - ||u_1||_1^2). $
  Dividing by $4 m$ proves the first claim. Requiring the bracketed expression to be strictly positive rearranges directly to $||u_1||_1^2 / n_g > |lambda_("min")| / (lambda_1 + |lambda_("min")|)$. Conversely, when $u_1$ is localized on a tiny dense subgraph ($||u_1||_1^2 << n_g$) while the rest of $g$ is strongly bipartite ($|lambda_("min")| >> lambda_1$), the negative spectral tail can overwhelm $lambda_1 ||u_1||_1^2$, rendering unrefined sign rounding suboptimal or even negative-gain.
]

#subsec("C", "Exact Spectrum and Multi-Way Degeneracy on Block Graphs")

When a graph contains $K >= 3$ symmetric communities, Newman's single-vector bisection encounters eigenspace degeneracy, as we prove in closed form on the Cartesian product graph $K_K square K_q$.

#theorem(title: "Tightness of K-Way Bound and Leading Eigenspace Degeneracy", origin: "Original Proof")[
  For integers $K >= 3$ and $q > K$, consider the symmetric regular $K$-community graph $G = K_K square K_q$ on $n = K q$ vertices (comprising $K$ cliques $C_1, dots, C_K tilde.op K_q$ linked by $q$ vertex-disjoint inter-clique $K_K$ matchings). Then:
  (i) The exact spectrum of $B$ is $lambda_1 = q - 2$ (multiplicity $K - 1$, eigenspace $E_1 = bold(1)_K^perp times.o bold(1)_q$), $lambda = K - 2$ (multiplicity $q - 1$), $lambda = 0$ (multiplicity $1$, vector $bold(1)$), and $lambda_("min") = -2$ (multiplicity $(K - 1)(q - 1)$);
  (ii) The true $K$-partition ${C_1, dots, C_K}$ achieves the exact universal spectral bound of Theorem II.A:
  $ Q({C_1, dots, C_K}) = ((K - 1)(q - 2)) / (K (q + K - 2)) = (n (K - 1)) / (2 m K) lambda_1 (B); $
  (iii) Because $dim(E_1) = K - 1 >= 2$, valid leading eigenvectors such as $u_1 prop (1, -1, 0, dots, 0)^T times.o bold(1)_q$ vanish identically on $K - 2$ full cliques, causing arbitrary zero-sign tie-breaking to bisect intact cliques irreversibly.
]

#proof[
  Write the adjacency matrix of $K_K square K_q$ via Kronecker products as
  $ A = I_K times.o (J_q - I_q) + (J_K - I_K) times.o I_q, $
  where $J_r = bold(1)_r bold(1)_r^T$. Every vertex has exact degree $k = (q - 1) + (K - 1) = q + K - 2$, so $2 m = K q (q + K - 2)$ and $B = A - (q + K - 2) / (K q) J_(K q)$. The matrices $J_K - I_K$ and $J_q - I_q$ commute with $J_K$ and $J_q$ respectively: $J_K - I_K$ has eigenvalue $K - 1$ on $bold(1)_K$ and $-1$ on $bold(1)_K^perp$; $J_q - I_q$ has eigenvalue $q - 1$ on $bold(1)_q$ and $-1$ on $bold(1)_q^perp$. Since $J_(K q) = J_K times.o J_q$ acts as $K q$ on $bold(1)_K times.o bold(1)_q$ and annihilates the three orthogonal subspaces, the eigenvalues of $B$ are:
  - on $bold(1)_K times.o bold(1)_q$: $(q - 1) + (K - 1) - (q + K - 2) = 0$ (mult. $1$);
  - on $bold(1)_K^perp times.o bold(1)_q$: $(q - 1) + (-1) = q - 2$ (mult. $K - 1$);
  - on $bold(1)_K times.o bold(1)_q^perp$: $(-1) + (K - 1) = K - 2$ (mult. $q - 1$);
  - on $bold(1)_K^perp times.o bold(1)_q^perp$: $(-1) + (-1) = -2$ (mult. $(K - 1)(q - 1)$).
  Because $q > K$, $lambda_1 (B) = q - 2$ is strictly the largest eigenvalue with multiplicity $K - 1$.

  For (ii), each clique $C_r$ contains $binom(q, 2)$ internal edges and volume fraction $a_r = 1 / K$, while $m = K q (q + K - 2) / 2$. Thus
  $ Q({C_r}) &= (K binom(q, 2)) / m - K (1 / K)^2 = (q - 1) / (q + K - 2) - 1 / K \
             &= ((K - 1)(q - 2)) / (K (q + K - 2)). $
  Evaluating Theorem II.A's upper bound with $n = K q$, $2 m = K q (q + K - 2)$, and $lambda_1 (B) = q - 2$ yields the identical value.

  For (iii), since $K >= 3$, the vector $w = (1, -1, 0, dots, 0)^T in bold(1)_K^perp$ is non-zero, so $u_1 = (2 q)^(-1 / 2) (w times.o bold(1)_q)$ is a unit leading eigenvector of $B$. On cliques $C_3, dots, C_K$, every entry of $u_1$ is identically zero (or perturbed to arbitrary tiny signs under floating-point roundoff), splitting intact cliques across $g_+$ and $g_-$ where divisive bisection without KL refinement can never reunite them.
]

#subsec("D", "Davis-Kahan Recovery in Balanced Block Models")

#theorem(title: "Spectral Consistency in the Balanced 2-Block SBM", origin: "Classical — Cited")[
  Let $G$ be drawn from a symmetric 2-block Stochastic Block Model (SBM) @holland1983sbm @abbe2018sbmreview on even $n$ vertices with balanced ground-truth spins $s^* in {-1, +1}^n$ ($bold(1)^T s^* = 0$), intra-block edge probability $p$, and inter-block edge probability $q < p$. Define the expected modularity matrix $bar(B) = bb(E)[A] - (bar(d) / n) J_n$ where $bar(d) = n (p + q) / 2 - p$, and spectral gap $gamma = n (p - q) / 2 - p > 0$. If $||B - bar(B)||_2 < gamma / 2$, then by the Davis-Kahan $sin Theta$ theorem @davis1970rotation, the fraction $eta$ of vertices misclassified by sign rounding $u_1 (B)$ satisfies
  $ eta <= (8 ||B - bar(B)||_2^2) / gamma^2. $
]

#proof[
  In the balanced 2-block SBM, $bb(E)[A] = ((p + q) / 2) J_n + ((p - q) / 2) s^* (s^*)^T - p I_n$. Subtracting $(bar(d) / n) J_n = ((p + q) / 2 - p / n) J_n$ yields
  $ bar(B) = ((p - q) / 2) s^* (s^*)^T + (p / n) J_n - p I_n. $
  Check its eigenstructure: $bar(B) bold(1) = bold(0)$ (since $bold(1)^T s^* = 0$); $bar(B) s^* = (n (p - q) / 2 - p) s^* = gamma s^*$; and $bar(B) v = -p v$ for all $v perp "span"{bold(1), s^*}$. Thus $u_1 (bar(B)) = s^* / sqrt(n)$ is the unique leading unit eigenvector with spectral gap $lambda_1 (bar(B)) - lambda_2 (bar(B)) = gamma - 0 = gamma$. By the Davis-Kahan $sin Theta$ theorem @davis1970rotation, when $||B - bar(B)||_2 < gamma / 2$,
  $ min_(sigma in {-1, +1}) ||u_1 (B) - sigma s^* / sqrt(n)||_2 <= (2 sqrt(2) ||B - bar(B)||_2) / gamma. $
  Every misclassified vertex $i$ has opposite sign to $sigma s^*_i$, hence $|u_(1, i) (B) - sigma s^*_i / sqrt(n)| >= 1 / sqrt(n)$. Summing squared errors over the $M = eta n$ misclassified vertices gives $eta = M / n <= ||u_1 (B) - sigma s^* / sqrt(n)||_2^2 <= 8 ||B - bar(B)||_2^2 / gamma^2$. In sparse regimes with $n(p - q)^2 / (p + q) -> infinity$, standard matrix concentration gives $||B - bar(B)||_2 = cal(O)_("Pr")(sqrt(n(p + q)))$, implying $eta -> 0$ (weak consistency) @abbe2018sbmreview.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Implicit Matrix-Vector Complexity.* Explicitly materializing $B$ destroys sparsity and requires $Theta(n^2)$ memory and $cal(O)(n^3)$ eigensolver work. Instead, Newman's implementation exploits the rank-1 structure $B = A - (2 m)^(-1) k k^T$: each matrix-vector product $B^((g)) v$ costs only $cal(O)(|E_g| + n_g)$ time and $cal(O)(|E_g| + n_g)$ memory using CSR sparse adjacency storage plus one inner product $k_g^T v$. Using Lanczos or shifted power iteration (shifting $B^((g)) + beta I$ with $beta >= ||B^((g))||_2$ so the algebraically largest eigenvalue becomes dominant in magnitude), convergence to tolerance $epsilon$ takes $T_("eig") = cal(O)(delta_("gap")^(-1) log(n / epsilon))$ iterations, or $cal(O)(n)$ worst-case iterations. Thus a single split costs $cal(O)(n (m + n))$ in the worst case @newman2006leading, and full recursive splitting of depth $d_("tree") = cal(O)(log n)$ on sparse graphs runs in $cal(O)(n^2 log n)$ time (and often near-linear $cal(O)((m + n) log^2 n)$ when spectral gaps remain bounded away from zero). Optional KL fine-tuning adds $cal(O)(n_g (|E_g| + n_g))$ per split if run naively, or $cal(O)((|E_g| + n_g) log n_g)$ with priority queues.

*Spectral Geometry.* Leading Eigenvector embeds each subgraph $g$ onto a 1D spectral axis $u_1 (B^((g)))$ and cuts at the origin $0$. This global embedding resists single-edge local noise (unlike CNM), but inherits two geometric liabilities: (i) recursive bipartitioning imposes a binary hierarchy that struggles with symmetric multi-way star/ring geometries (Theorem II.C), and (ii) modularities still obey the $sqrt(m / 2)$ resolution limit @fortunato2007resolution.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(n (m + n))$ per split worst-case; $cal(O)(n^2 log n)$ sparse tree],
    [Space Complexity], [$cal(O)(m + n)$ via implicit rank-1 operator $A v - (2 m)^(-1) (k^T v) k$],
    [Cluster Count $K$], [Automatic via exact certificate $lambda_1 (B^((g))) <= 0$ (Thm. II.A)],
    [Determinism], [Deterministic up to sign flip and degenerate eigenspaces],
    [Connectedness], [Not guaranteed without post-processing (zero-crossings can disconnect)],
    [Invariances], [Positive weight-scale invariant; permutation equivariant],
    [Key Limitation], [Binary split bias (Thm. II.C); localization trap (Lem. II.B)],
  ),
  caption: [Operational and complexity profile of Newman's Leading Eigenvector algorithm.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Dense Weighted Correlation and Financial Networks:* Stock return correlation matrices and functional-MRI connectomes are dense ($m = Theta(n^2)$) and naturally centered around a market or global-signal mode. Because $B$ explicitly subtracts the rank-1 degree outer product $k k^T / (2 m)$ and operates via Lanczos spectral projections, Leading Eigenvector isolates sector blocks cleanly where sparse greedy tree structures (CNM) lose their advantage.
2. *Social and Political Polarization Analysis:* Networks exhibiting strong ideological bipartition (e.g., congressional co-voting, blogosphere link graphs @newman2006leading) have a dominant spectral gap $lambda_1 >> lambda_2$, placing them directly in the Davis-Kahan high-SNR regime (Theorem II.D) where sign rounding achieves near-exact recovery.
3. *Modularity Indivisibility Auditing:* Because Theorem II.A certifies that $lambda_1 (B^((g))) <= 0$ precludes *any* $K$-way modularity-improving split, practitioners use the leading eigenvalue of $B^((g))$ as a principled stopping diagnostic for over-partitioning in hierarchical clustering pipelines.
4. *Medium-Scale Scientific Collaboration Networks:* Co-authorship graphs with $10^3$--$10^5$ nodes and well-separated disciplinary clusters benefit from the global view of $u_1$, avoiding early greedy bridges between junior cross-disciplinary authors.
5. *Biophysical Contact-Map Domain Decomposition:* Residue-residue interaction graphs of multi-domain proteins have moderate size ($n <= 2000$) and dense intra-domain contacts; recursive spectral bisection identifies rigid folding domains and hinge boundaries cleanly.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

The method degrades in four settings: (i) *symmetric multi-way clusters or multi-shape manifolds* where eigenspace degeneracy (Theorem II.C) and 1D zero-crossing cuts shatter curved or multi-component $k$NN shapes; (ii) *sparse heavy-tailed hubs* where $u_1$ localizes on a tiny clique while negative bipartite modes dominate (Lemma II.B); (iii) *disconnected or tiny sub-resolution communities* below $sqrt(m / 2)$ @fortunato2007resolution; and (iv) *massive billion-edge graphs* where repeated $cal(O)(n)$ Lanczos sweeps are far slower than linear-time Louvain @blondel2008louvain or Leiden @traag2019leiden.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Polarized 2-Block Networks], [Very high: Davis-Kahan optimal], [Use KL refinement @kernighan1970efficient],
    [Dense Correlation Graphs], [High: implicit rank-1 Lanczos], [Check multi-way $K >= 3$ degeneracy],
    [Multi-Shape $k$NN Manifolds], [Low: 1D spectral zero-crossing cuts], [Use CNM (Alg. 115) or HDBSCAN @campello2013hdbscan],
    [Degree-Skewed Mega-Graphs], [Low: localization; $cal(O)(n^2 log n)$], [Use Leiden @traag2019leiden or PIC @lin2010pic],
    [Multi-Scale Micro-Modules], [Low: modularity resolution limit], [Use Spinglass @reichardt2006spinglass],
  ),
  caption: [Domain suitability matrix for Newman's Leading Eigenvector method.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Leading Eigenvector bridges modularity maximization (CNM Alg. 115, Spinglass Alg. 117) and classical spectral graph partitioning (Ratio Cut @hagen1992ratiocut, Normalized Cut @shi2000ncut, Ng-Jordan-Weiss @ng2002spectral, and spectral clustering tutorials @vonluxburg2007tutorial). On regular graphs, this connection becomes an exact eigenbasis identity.

#proposition(title: "Exact Spectral Equivalence to Fiedler Partitioning on Regular Graphs", origin: "Original Proof")[
  Let $G$ be a $d$-regular connected graph ($k_i = d$ for all $i in V$) with unnormalized Laplacian $L = d I - A$ having eigenvalues $0 = mu_1 < mu_2 <= dots <= mu_n$. Then the modularity matrix satisfies $B = d I - L - (d / n) J_n$, shares the exact same eigenvectors as $L$ on $bold(1)^perp$, and maps eigenvalues via $lambda_r (B) = d - mu_(r + 1) (L)$ for $r = 1, dots, n - 1$ (with $lambda = 0$ on $bold(1)$). Consequently, Newman's leading eigenvector $u_1 (B)$ is identically the Fiedler vector $v_2 (L)$ @chung1997spectral @hagen1992ratiocut.
]

#proof[
  For a $d$-regular graph, $k = d bold(1)$ and $2 m = n d$, so
  $ B = A - (d^2 / (n d)) J_n = d I - L - (d / n) J_n. $
  On $bold(1)$, $L bold(1) = bold(0)$ and $(d I - (d / n) J_n) bold(1) = bold(0)$, giving $B bold(1) = bold(0)$. For every Laplacian eigenvector $v_j perp bold(1)$ ($j = 2, dots, n$) with $L v_j = mu_j v_j$, we have $J_n v_j = bold(0)$, hence
  $ B v_j = (d I - L) v_j = (d - mu_j) v_j. $
  Thus maximizing $lambda(B) = d - mu_j$ over $j >= 2$ is strictly equivalent to minimizing $mu_j$ over $j >= 2$, selecting the Fiedler vector $v_2 (L)$ (with $lambda_1 (B) = d - mu_2 (L)$). On non-regular graphs, the equivalence breaks: $L_("sym") = I - D^(-1 / 2) A D^(-1 / 2)$ normalizes degrees multiplicatively, whereas $B = A - (2 m)^(-1) k k^T$ subtracts an additive rank-1 hub penalty.
]

*Empirical Behavior on ClusterOverview Benchmarks.* Across the 4-dataset ClusterOverview suite, Leading Eigenvector achieves *Iris* $"Sil" = 0.309$, $"ARI" = 0.495$; *Customers* $"Sil" = 0.370$; *Multi-Shape* $"Sil" = 0.267$; and *High-D* $"Sil" = 0.302$, with a mean runtime of $1.774$ ms. The empirical profile aligns cleanly with our theorems: (1) on *Iris*, its global spectral view outperforms CNM's greedy ARI ($0.495$ vs. $0.373$) by avoiding local low-degree bridge merges between Versicolor and Virginica, consistent with Davis-Kahan stability (Theorem II.D); (2) its implicit rank-1 Lanczos iteration runs nearly twice as fast as CNM ($1.774$ ms vs. $3.417$ ms) on these moderate-$n$ graphs; yet (3) on *Multi-Shape*, Silhouette drops sharply to $0.267$ (compared with $0.782$ for CNM and $0.543$ for HCS) because multi-component non-convex $k$NN manifolds exhibit near-degenerate positive eigenspaces (Theorem II.C) and zero-crossing sign cuts that slice elongated curves prematurely. On *Customers* ($0.370$) and *High-D* ($0.302$), it delivers balanced, reproducible spectral baselines at sub-2 ms runtime.
