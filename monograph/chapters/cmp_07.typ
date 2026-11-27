#import "../template/ieee.typ": *

#cmp_entry(
  roman: "VII",
  title: "Intra-Category Comparative Analysis: Graph, Spectral, and Network Algorithms (Algs. 103-132)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 103-132",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part V of this monograph develops thirty graph, spectral, and network clustering algorithms (Algs. 103--132) spanning Laplacian relaxations, bipartite co-clustering, power-iteration embeddings, greedy/spectral/Potts modularity optimizers, random-walk flow compressors, matrix-inflation dynamics, linear-time label propagation, factor-graph message passing, divisive combinatorial cuts, multilevel balanced partitioners, overlapping and local statistical covers, and generative block models @vonluxburg2007tutorial @fortunato2007resolution @abbe2018sbmreview @aggarwal2013data @jain2010data. Whereas prototype methods (Part I) partition metric coordinates and statistical mixtures (Part IV) fit feature densities, Category V clusters *relational topology*: a weighted graph $G = (V, E, W)$ with $n = |V|$ vertices, $m = |E|$ edges, adjacency $W in bb(R)_(>= 0)^(n times n)$, vertex strengths $d_i = sum_j W_(i j)$, $D = "diag"(d_1, dots, d_n)$, and volume $"vol"(V) = sum_i d_i = 2 m$.

Across all thirty methods, partitions or covers $cal(C) = {C_1, dots, C_K}$ optimize a unified *relational variational functional*:
$
  cal(J)(cal(C); G) &= sum_(r=1)^K Phi_r (W_(C_r C_r), W_(C_r overline(C)_r), \
  &quad P_(C_r C_r), |C_r|, "vol"(C_r)),
$
balancing intra-module affinity $W_(C_r C_r)$ against boundary cut $W_(C_r overline(C)_r) = "cut"(C_r, overline(C)_r)$, a null expectation $P_(C_r C_r)$ (configuration model $d_i d_j / (2 m)$ or constant baseline $gamma$), and cardinality or volume normalization. Four structural fault lines organize this head-to-head study: (I) *Spectral Cuts & Laplacian Families* (Algs. 103--105, 123, 132); (II) *Modularity, Resolution Limits & Connectedness* (Algs. 107--108, 115--117, 129); (III) *Random-Walk, Flow, Label & Message-Passing Dynamics* (Algs. 106, 109--111, 113--114, 119); and (IV) *Divisive Cuts, Multilevel Partitioning, Overlapping & Generative Block Models* (Algs. 112, 118, 120--122, 124--128, 130--131), followed by (V) a complete 30-algorithm matrix, reduction lattice, decision guide, and open frontiers.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Spectral Cuts & Laplacian Families")

Five Category V algorithms relax NP-hard graph cuts via linear spectral operators: *Ng--Jordan--Weiss Spectral Clustering* (NJW, \#103) @ng2002spectral, *Normalized Cuts* (Ncut, \#104) @shi2000ncut, *Ratio Cut* (\#105) @hagen1992ratiocut, *Spectral Biclustering* (Kluger, \#123) @kluger2003spectral, and *Power Iteration Clustering* (PIC, \#132) @lin2010pic. On a connected weighted graph $G = (V, E, W)$ ($d_i > 0$), they operate on three canonical operators:
$
  L &:= D - W, \
  L_"rw" &:= D^(-1) L = I_n - P, quad P := D^(-1) W, \
  L_"sym" &:= D^(-1/2) L D^(-1/2) = I_n - D^(-1/2) W D^(-1/2).
$
Symmetric PSD $L_"sym"$ is similar to $L_"rw" = D^(-1/2) L_"sym" D^(1/2)$, sharing real eigenvalues $0 = lambda_1 < lambda_2 <= dots <= lambda_n <= 2$. Yet selecting $L$ vs. $L_"rw"$ or $L_"sym"$, bipartite scaling, or 1D power truncation creates decisive separations in cut tightness, asymptotic consistency, and scalability @vonluxburg2007tutorial @chung1997spectral.

#subsec("A", "Cut Relaxations, Cheeger Bounds & NJW–Ncut Equivalence")

For a $K$-partition $cal(C) = {C_1, dots, C_K}$, *Ratio Cut* (\#105) @hagen1992ratiocut normalizes cuts by cardinality $|C_r|$, whereas *Normalized Cut* (\#104) @shi2000ncut normalizes by degree volume $"vol"(C_r) = sum_(i in C_r) d_i$:
$
  "Rcut"(cal(C)) &:= sum_(r=1)^K ("cut"(C_r, overline(C)_r)) / (|C_r|), \
  "Ncut"(cal(C)) &:= sum_(r=1)^K ("cut"(C_r, overline(C)_r)) / ("vol"(C_r)).
$
For bisections $(S, overline(S))$, conductance is $phi(S) := "cut"(S, overline(S)) / min("vol"(S), "vol"(overline(S)))$ and $phi(G) := min_(0 < "vol"(S) <= "vol"(V)/2) phi(S)$ @chung1997spectral.

#theorem(title: "Cut Relaxation Equivalence, Cheeger Sandwich & NJW Row Normalization", origin: "Original Unified Proof")[
  Let $G = (V, E, W)$ be connected ($d_i > 0$) with $L_"sym"$ eigenvalues $0 = lambda_1 < lambda_2 <= dots <= lambda_n$.
  1. *(Exact Trace Identities):* Defining $H_(i r) = |C_r|^(-1/2) bold(1)[i in C_r]$ and $T_(i r) = sqrt(d_i / "vol"(C_r)) bold(1)[i in C_r]$, both obey $H^top H = T^top T = I_K$ and satisfy identically
     $
       "Rcut"(cal(C)) = "Tr"(H^top L H), \
       "Ncut"(cal(C)) = "Tr"(T^top L_"sym" T),
     $
     yielding Ky Fan spectral lower bounds $"Rcut"^* >= sum_(r=1)^K lambda_r (L)$ and $"Ncut"^* >= sum_(r=1)^K lambda_r (L_"sym")$.
  2. *(Two-Sided Cheeger Sandwich for Shi--Malik Sweep Cut \#104 @chung1997spectral @shi2000ncut):* Let $u_2 = D^(-1/2) v_2$ be the Fiedler vector of $L_"rw"$ ($L u_2 = lambda_2 D u_2$, $u_2^top D bold(1) = 0$). The best threshold sweep cut $S_"sweep"$ along $u_2$ satisfies
     $
       (lambda_2 (L_"sym")) / 2 &<= phi(G) <= phi(S_"sweep") <= sqrt(2 lambda_2 (L_"sym")), \
       lambda_2 (L_"sym") &<= "Ncut"^* (G) <= "Ncut"(S_"sweep") <= 2 sqrt(2 lambda_2 (L_"sym")).
     $
  3. *(Exact Degree-Cancellation of NJW Row Normalization \#103 vs. $L_"rw"$ \#104):* On an ideal $K$-component graph, every orthonormal basis $V in bb(R)^(n times K)$ of $ker(L_"sym")$ satisfies $V = T Q$ ($Q in O(K)$), so $V_(i :) = sqrt(d_i / "vol"(C_r)) Q_(r :)$ for $i in C_r$. NJW row normalization $Y_(i :) := V_(i :) / norm(V_(i :))_2$ cancels $sqrt(d_i)$ identically, collapsing every vertex of $C_r$ onto the unit simplex vertex $Q_(r :) in bb(S)^(K-1)$ (pairwise distance $sqrt(2)$), matching piecewise-constant $L_"rw"$ eigenvectors up to orthogonal rotation.
]

#proof[
  *Part 1:* Since $x^top L x = (1 / 2) sum_(i, j) W_(i j) (x_i - x_j)^2$, column $r$ of $H$ gives
  $
    h_r^top L h_r &= 1 / 2 sum_(i in C_r, j in.not C_r) W_(i j) (2 / (|C_r|)) = ("cut"(C_r, overline(C)_r)) / (|C_r|).
  $
  Similarly, $u_r = D^(-1/2) t_r$ has $(u_r)_i = "vol"(C_r)^(-1/2) bold(1)[i in C_r]$, giving $t_r^top L_"sym" t_r = u_r^top L u_r = "cut"(C_r, overline(C)_r) / "vol"(C_r)$. Summing over $r = 1, dots, K$ and applying Ky Fan's variational principle over $X^top X = I_K$ yields the bounds.

  *Part 2 (Lower Bound):* For $S$ attaining $phi(G)$ ($"vol"(S) <= "vol"(V) / 2$), set $x_i = "vol"(overline(S))$ ($i in S$) and $x_i = -"vol"(S)$ ($i in overline(S)$). Then $x^top D bold(1) = 0$, so by Courant--Fischer:
  $
    lambda_2 (L_"sym") &<= (x^top L x) / (x^top D x) \
    &= ("cut"(S, overline(S)) ("vol"(S) + "vol"(overline(S)))^2) / ("vol"(S) "vol"(overline(S)) "vol"(V)) \
    &= "cut"(S, overline(S)) ( 1 / ("vol"(S)) + 1 / ("vol"(overline(S))) ) \
    &= "Ncut"(S, overline(S)) <= 2 phi(G).
  $

  *Part 2 (Sweep-Cut Upper Bound):* Order vertices $u_1 <= dots <= u_n$ ($u = u_2$) and pick weighted median $m^*$ with $"vol"({1, dots, m^* - 1}), "vol"({m^* + 1, dots, n}) <= "vol"(V) / 2$. Flipping sign if needed, assume $u_(m^*) <= 0$. Since $u != bold(0)$ and $u^top D bold(1) = 0$, $max_i u_i > 0 >= u_(m^*)$, so $y^+ := max(u - u_(m^*) bold(1), 0) != bold(0)$ with $"vol"("supp"(y^+)) <= "vol"(V) / 2$. On every edge $(i, j)$, $(y_i^+ - y_j^+)^2 <= (y_i^+ - y_j^+)(u_i - u_j)$ (checking all sign cases), hence
  $
    (y^+)^top L y^+ &<= (y^+)^top L u = lambda_2 (y^+)^top D u \
    &= lambda_2 (y^+)^top D y^+ + lambda_2 u_(m^*) (y^+)^top D bold(1) \
    &<= lambda_2 (y^+)^top D y^+,
  $
  as $u_(m^*) <= 0$. Rescale $z := y^+ / max_i y_i^+ in [0, 1]^n$ and let $S_t := {i : z_i^2 >= t}$ ($t in (0, 1]$), each obeying $"vol"(S_t) <= "vol"(V) / 2$. Integrating $"cut"(S_t, overline(S)_t) >= phi(S_"sweep") "vol"(S_t)$ via the layer-cake identity yields
  $
    phi(S_"sweep") z^top D z &<= integral_0^1 "cut"(S_t, overline(S)_t) dif t \
    &= sum_(i < j) W_(i j) |z_i^2 - z_j^2| \
    &= sum_(i < j) W_(i j) |z_i - z_j| (z_i + z_j).
  $
  By Cauchy--Schwarz and $(z_i + z_j)^2 <= 2(z_i^2 + z_j^2)$,
  $
    ( phi(S_"sweep") z^top D z )^2 &<= (z^top L z) ( 2 sum_(i < j) W_(i j) (z_i^2 + z_j^2) ) \
    &= 2 (z^top L z) (z^top D z) <= 2 lambda_2 (z^top D z)^2,
  $
  giving $phi(S_"sweep") <= sqrt(2 lambda_2)$ and $"Ncut"(S_"sweep") <= 2 sqrt(2 lambda_2)$.

  *Part 3:* Since $T$ spans $ker(L_"sym")$, $V = T Q$ gives $norm(V_(i :))_2 = sqrt(d_i / "vol"(C_r))$ ($i in C_r$), so $Y_(i :) = Q_(r :)$, independent of $d_i$ with $norm(Q_(r :) - Q_(s :))_2 = sqrt(2)$ ($r != s$). Without row normalization, vertices of $C_r$ scatter radially by $sqrt(d_"max" / d_"min")$, causing $k$-Means @lloyd1982least to split clusters by degree.
]

Moreover, $K$-way Ncut (\#104) / NJW (\#103) is algebraically equivalent to weighted kernel $k$-Means with kernel $sigma D^(-1) + D^(-1) W D^(-1)$ @dhillon2004kernel, while subspace robustness under noise scales inversely with the eigengap $delta_K = lambda_(K+1) - lambda_K$ by the Davis--Kahan $sin Theta$ theorem @davis1970rotation.

#subsec("B", "Operator Consistency: Unnormalized L vs. Normalized Laplacians")

In the large-sample limit $n -> infinity$, von Luxburg, Belkin, and Bousquet @vonluxburg2007tutorial exposed a sharp consistency dichotomy separating Ratio Cut (\#105) from Ncut/NJW (\#104/\#103): unnormalized $L$ can localize onto single peripheral vertices even as $n -> infinity$.

#theorem(title: "Unnormalized vs. Normalized Laplacian Consistency & Leaf-Localization Counterexample", origin: "Original Unified Proof")[
  1. *(Constructive Leaf-Snip Breakdown of Ratio Cut \#105 vs. Ncut \#104 / NJW \#103):* Take two $p$-cliques $K_p^((1)), K_p^((2))$ ($p >= 4$) linked by $b = ceil(c p)$ bridge edges ($c in (1/2, 1)$), plus one degree-1 leaf $v_0$ attached to $K_p^((1))$. Comparing the trivial leaf cut $S_"leaf" = {v_0}$ against the true bisection $S_"true" = K_p^((1)) union {v_0}$ yields:
     $
       "Rcut"(S_"leaf") &= 1 + 1 / (2 p) < 2 c - O(p^(-1)) = "Rcut"(S_"true"), \
       "Ncut"(S_"true") &= (2 c) / p + O(p^(-2)) << 1 + O(p^(-2)) = "Ncut"(S_"leaf").
     $
     For every $p >= 6$ (e.g., $p = 10, b = 6$: $"Rcut"(S_"leaf") = 1.05 < 1.145 = "Rcut"(S_"true")$, whereas $"Ncut"(S_"true") = 0.124 << 1.005 = "Ncut"(S_"leaf")$), Ratio Cut strictly prefers cutting off the single leaf, while Ncut/NJW recover the true bisection by an $Omega(p)$ margin!
  2. *(Continuous-Spectrum Localization vs. Compact Operator Consistency @vonluxburg2007tutorial):* On i.i.d. samples from continuous density $p$ on compact $cal(X)$ with continuous kernel $k(x, y) >= k_"min" > 0$, degree function $d(x) = integral k(x, y) p(y) dif y in [d_"min", d_"max"]$, and integral operator $(K f)(x) = integral k(x, y) f(y) p(y) dif y$:
     - $n^(-1) L_n$ converges to $U = M_d - K$, whose essential spectrum is the entire interval $[d_"min", d_"max"]$; Ratio Cut (\#105) is consistent *only* for isolated eigenvalues $lambda_k < d_"min"$, failing via spike localization whenever $lambda_k >= d_"min"$.
     - $L_"rw"$ converges to $I - T$ with *compact* $T = M_(1 / d) K$, collapsing the essential spectrum to ${1}$; every $lambda_k < 1$ is isolated of finite multiplicity, guaranteeing almost-sure consistency of Ncut (\#104) and NJW (\#103).
]

#proof[
  *Part 1:* Cutting $S_"leaf" = {v_0}$ ($|S_"leaf"| = 1$, $|overline(S)_"leaf"| = 2 p$) gives $"Rcut"(S_"leaf") = 1 + 1 / (2 p)$, whereas cutting $S_"true"$ ($|S_"true"| = p + 1, |overline(S)_"true"| = p$) gives
  $
    "Rcut"(S_"true") = b ( 1 / (p + 1) + 1 / p ) = (b (2 p + 1)) / (p (p + 1)) >= 2 c - O(p^(-1)),
  $
  exceeding $"Rcut"(S_"leaf")$ for $p > 1 / (2(2 c - 1))$. Conversely, with $"vol"({v_0}) = 1$, $"vol"(S_"true") = p(p - 1) + b + 2$, and $"vol"(overline(S)_"true") = p(p - 1) + b$:
  $
    "Ncut"(S_"leaf") &= 1 + 1 / (2 p(p - 1) + 2 b + 1) = 1 + O(p^(-2)), \
    "Ncut"(S_"true") &= b ( 1 / (p(p - 1) + b + 2) + 1 / (p(p - 1) + b) ) \
    &= (2 c) / p + O(p^(-2)).
  $

  *Part 2:* On compact $cal(X)$, continuity of $k$ and $d(x) >= k_"min" > 0$ render both $K$ and $T = M_(1 / d) K$ compact on $C(cal(X))$ by Arzelà--Ascoli. By Weyl's theorem, compact perturbations preserve essential spectra: $sigma_"ess" (M_d - K) = sigma_"ess" (M_d) = [d_"min", d_"max"]$ (trapping eigenvalues $>= d_"min"$ inside continuous spectrum), whereas $sigma_"ess" (I - T) = {1}$. By Riesz--Schauder and collective compact convergence $T_n -> T$ a.s., all spectral projections of $L_"rw"$ (and $D^(1/2)$-rescaled $L_"sym"$) with $lambda_k < 1$ converge almost surely @vonluxburg2007tutorial.
]

#subsec("C", "Bipartite Co-Clustering (Kluger) vs. Truncated Power Iteration (PIC)")

*Spectral Biclustering* (Kluger, \#123) @kluger2003spectral lifts normalized cuts to rectangular matrices $A in bb(R)_(>= 0)^(n_r times n_c)$, while *Power Iteration Clustering* (PIC, \#132) @lin2010pic avoids $O(n^3)$ eigendecomposition via a single $O(T m)$ truncated power vector.

#theorem(title: "Kluger #123 Bipartite Ncut Duality & PIC #132 Window / 1D Collision Bounds", origin: "Original Unified Proof")[
  1. *(Bipartite Ncut Equivalence & Bistochastic Scaling of Kluger \#123 @kluger2003spectral):* For $A in bb(R)_(>= 0)^(n_r times n_c)$ with $D_r = "diag"(A bold(1))$, $D_c = "diag"(A^top bold(1))$, and bipartite graph $W_"bip" = mat(0, A; A^top, 0)$, the eigenvalues of $L_"sym" (W_"bip")$ occur in symmetric pairs ${1 - sigma_j, 1 + sigma_j}$, where $(sigma_j, u_j, v_j)$ are singular triples of $tilde(A) := D_r^(-1/2) A D_c^(-1/2)$, and $(D_r^(-1/2) u_j, plus.minus D_c^(-1/2) v_j)$ are simultaneous row--column eigenvectors of $L_"rw" (W_"bip")$. Sinkhorn--Knopp bistochastic scaling $B = R A C$ cancels multiplicative main effects $A_(i j) = alpha_i beta_j M_(c(i), d(j))$, turning ideal checkerboards into exact piecewise-constant singular vectors.
  2. *(PIC \#132 Spectral Filtering Window & $Theta(K^(-2))$ 1D Collision Failure @lin2010pic):* Iterating $v^((t+1)) = (P v^((t))) / norm(P v^((t)))_1$ from random $v^((0)) = sum_j c_j u_j$ yields
     $
       v^((t)) prop c_1 bold(1) + sum_(j=2)^K c_j lambda_j^t u_j + sum_(j=K+1)^n c_j lambda_j^t u_j.
     $
     - Stopping inside the window
       $
         (ln(1 / epsilon)) / (ln(lambda_K / abs(lambda_(K+1)))) <= t <= (ln(1 / eta)) / (ln(1 / lambda_K))
       $
       suppresses noise modes by $<= epsilon$ while preserving cluster modes $>= eta$ in $O(t m)$ time.
     - *1D Collision Pathology:* Even on a noiseless $K$-component graph mapping cluster $r$ to scalar $Z_r tilde.op "Unif"[0, 1]$ i.i.d., the tightest inter-cluster gap collapses quadratically in $K$:
       $
         bb(E)[ min_(1 <= a < b <= K) |Z_a - Z_b| ] = 1 / (K^2 - 1) = Theta(K^(-2)),
       $
       whereas $K$-dimensional NJW (\#103) maintains constant $sqrt(2)$ separation for every $K$.
]

#proof[
  *Part 1:* With $D_"bip" = "diag"(D_r, D_c)$, $I - L_"sym" (W_"bip") = mat(0, tilde(A); tilde(A)^top, 0)$. Applying this block matrix to $(u_j, plus.minus v_j)^top$ yields $plus.minus sigma_j (u_j, plus.minus v_j)^top$, giving eigenvalues $1 minus.plus sigma_j$ and $L_"rw"$ eigenvectors $(D_r^(-1/2) u_j, plus.minus D_c^(-1/2) v_j)$. When $A_(i j) = alpha_i beta_j M_(c(i), d(j))$, unique positive Sinkhorn--Knopp scalings $R, C$ absorb $alpha_i^(-1), beta_j^(-1)$ identically, restoring uniform-marginal checkerboard singular vectors.

  *Part 2:* Comparing mode ratios gives $(|lambda_(K+1)| / lambda_K)^t <= epsilon$ and $lambda_K^t >= eta$, producing the exact stopping window. For $Z_1, dots, Z_K tilde.op "Unif"[0, 1]$ with order statistics $Z_((1)) <= dots <= Z_((K))$, the $K + 1$ spacings $G_0, dots, G_K$ are uniform on $Delta^K$. Shifting each interior gap $G_r >= s$ ($r = 1, dots, K - 1$) shrinks the residual simplex side length to $1 - (K - 1) s$, giving $bb(P)(min_(1 <= r <= K-1) G_r >= s) = (1 - (K - 1) s)^K$ for $s in [0, 1 / (K - 1)]$. Integrating yields
  $
    bb(E)[ min_(a < b) |Z_a - Z_b| ] &= integral_0^(1 / (K - 1)) (1 - (K - 1) s)^K dif s \
    &= 1 / ((K - 1)(K + 1)) = 1 / (K^2 - 1).
  $
  At $K = 20$, the expected tightest gap is $0.0025$, easily swamped by small perturbation noise, limiting single-vector PIC to small $K <= 8$ unless multiple deflated vectors are used.
]

#wide_figure(
  caption: [*Table I:* Head-to-Head Comparison of Spectral Cuts & Laplacian Families (Section I: Algs. 103, 104, 105, 123, 132).],
  table(
    columns: (1.06fr, 0.96fr, 0.92fr, 0.88fr, 0.96fr, 1.22fr),
    align: (left, left, center, center, left, left),
    [Algorithm], [Operator / Target], [Consistency], [Time Cost], [Embedding Dim.], [Primary Strength / Failure Mode],
    [NJW (\#103) @ng2002spectral], [$L_"sym"$ + row norm.], [*Strong* (compact)], [$O(n^2 d + n^3)$], [$K$-dim. on $bb(S)^(K-1)$], [Degree-invariant $K$-way; dense $O(n^3)$],
    [Ncut (\#104) @shi2000ncut], [$L_"rw"$ / gen. $(L, D)$], [*Strong* (compact)], [$O(n^3)$ / Lanczos], [$1$-dim. rec. or $K$-dim.], [Cheeger $sqrt(2 lambda_2)$ bound; slow dense],
    [Ratio Cut (\#105) @hagen1992ratiocut], [Unnorm. $L = D - W$], [*Fails* ($lambda >= d_"min"$)], [$O(n^3)$ / Lanczos], [$K$-dim. Euclidean], [Cardinality balance; leaf-snip failure],
    [Kluger (\#123) @kluger2003spectral], [Bipartite $D_r^(-1/2) A D_c^(-1/2)$], [*Strong* (bipartite)], [$O(n_r n_c min(n_r, n_c))$], [Coupled row+col SVD], [Exact checkerboard co-clusters; nonneg. $A$],
    [PIC (\#132) @lin2010pic], [Truncated $P^t v^((0))$], [Gap-window dep.], [$O(T m + I n K)$], [*1D scalar* $v^((t))$], [Ultra-fast sparse; $Theta(K^(-2))$ 1D collisions],
  ),
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Modularity, Resolution Limits & Connectedness")

Six Category V algorithms discover communities with unknown $K$ by optimizing Newman--Girvan modularity or Potts spin-glass energy: *Fast Greedy Modularity* (CNM, \#115) @clauset2004cnm, *Leading Eigenvector* (Newman, \#116) @newman2006leading, *Spinglass* (Reichardt--Bornholdt, \#117) @reichardt2006spinglass, *Fast Unfolding* (Blondel, \#129) @blondel2008louvain, *Louvain* (\#107) @blondel2008louvain, and *Leiden* (\#108) @traag2019leiden. Here *Fast Unfolding (\#129)* denotes Blondel et al.'s original two-phase $Q_1$ formalization, while *Louvain (\#107)* denotes the broader modern framework adding tunable resolution $gamma$ and randomized move orders—yet inheriting the exact same disconnected-community defect resolved by *Leiden (\#108)* @traag2019leiden.

#subsec("A", "Generalized Modularity, Potts Equivalence & Spectral Bisection")

For adjacency $A_(i j) = W_(i j)$, degrees $k_i = sum_j A_(i j)$, $m = (1 / 2) sum_i k_i$, and resolution $gamma > 0$, generalized modularity is @newman2004finding @reichardt2006spinglass:
$
  Q_gamma (cal(C)) &:= 1 / (2 m) sum_(i, j = 1)^n ( A_(i j) - gamma (k_i k_j) / (2 m) ) delta(c_i, c_j) \
  &= sum_(r=1)^K [ (e_r) / m - gamma ( (K_r) / (2 m) )^2 ],
$
where $e_r = (1 / 2) sum_(i, j in C_r) A_(i j)$ and $K_r = sum_(i in C_r) k_i$. Exact $Q_1$ maximization is NP-complete @brandes2008modularity.

#proposition(title: "Exact Potts–Modularity Equivalence & Subgraph Modularity Matrix", origin: "Original Unified Proof")[
  1. *(Spinglass \#117 $equiv$ Generalized Modularity @reichardt2006spinglass):* Setting $p_(i j) = k_i k_j / (2 m)$ in the $q$-state Potts Hamiltonian $cal(H)_"RB" (bold(sigma)) := -sum_(i != j) (A_(i j) - gamma p_(i j)) delta(sigma_i, sigma_j)$ gives the exact affine identity
     $
       cal(H)_"RB" (bold(sigma)) = - 2 m Q_gamma (bold(sigma)) - gamma / (2 m) sum_(i=1)^n k_i^2.
     $
     Hence zero-temperature simulated annealing (\#117) targets the exact global maximum of $Q_gamma$ while supporting multi-scale $gamma$ and signed edges natively.
  2. *(Exact Recursive Subgraph Matrix of Leading Eigenvector \#116 @newman2006leading):* With $B_(i j) := A_(i j) - k_i k_j / (2 m)$, bisecting community $g subset.eq V$ via $s in {-1, +1}^(|g|)$ alters global modularity by
     $
       Delta Q &= 1 / (4 m) s^top B^((g)) s, \
       B^((g))_(i j) &:= B_(i j) - delta_(i j) sum_(l in g) B_(i l) quad (i, j in g).
     $
     Community $g$ is spectrally indivisible whenever $lambda_"max" (B^((g))) <= 0$.
]

#proof[
  *Part 1:* Because $A_(i i) = 0$ and $delta(sigma_i, sigma_i) = 1$, adding and subtracting the $i = j$ diagonal terms yields
  $
    cal(H)_"RB" (bold(sigma)) &= - sum_(i, j = 1)^n ( A_(i j) - gamma (k_i k_j) / (2 m) ) delta(sigma_i, sigma_j) \
    &quad - gamma / (2 m) sum_(i=1)^n k_i^2 = - 2 m Q_gamma (bold(sigma)) - "const".
  $
  *Part 2:* Pair $(i, j) in g times g$ remains co-assigned iff $(1 + s_i s_j) / 2 = 1$, so subtracting the pre-split contribution gives
  $
    Delta Q &= 1 / (2 m) sum_(i, j in g) B_(i j) ( (1 + s_i s_j) / 2 - 1 ) \
    &= 1 / (4 m) sum_(i, j in g) B_(i j) s_i s_j - 1 / (4 m) sum_(i in g) s_i^2 ( sum_(l in g) B_(i l) ) \
    &= 1 / (4 m) s^top B^((g)) s,
  $
  using $s_i^2 = 1$. Since $B^((g)) bold(1)_g = bold(0)$, $lambda_"max" (B^((g))) <= 0$ forces $Delta Q <= 0$ for every split $s$.
]

Unlike Section I's conductance cuts, Newman (\#116) spectrally maximizes an indefinite null-centered operator $B^((g))$ (typically refined by Kernighan--Lin \#120), automatically choosing $K$ yet remaining vulnerable to rank-1 bisection traps on symmetric $K >= 3$ graphs and modularity's resolution limit.

#subsec("B", "The Fortunato–Barthélemy Resolution Limit & Resolution-Free CPM")

Every faithful optimizer of $Q_1$—CNM (\#115), Leading Eigenvector (\#116), Louvain/Fast Unfolding (\#107/\#129), Spinglass at $gamma = 1$ (\#117), and Leiden-Modularity (\#108)—suffers the *Fortunato--Barthélemy resolution limit* @fortunato2007resolution, which Leiden's Constant Potts Model (CPM) eliminates completely.

#theorem(title: "Fortunato–Barthélemy Resolution Limit vs. Sub-Partition Invariance of Leiden-CPM", origin: "Original Unified Proof")[
  1. *(Generalized Resolution Limit @fortunato2007resolution):* Merging two communities $C_1, C_2$ (degree volumes $K_1, K_2$, $l_(1 2) >= 1$ cross-edges) strictly increases $Q_gamma$ iff
     $
       Delta Q_gamma (C_1, C_2) &= (l_(1 2)) / m - gamma (K_1 K_2) / (2 m^2) > 0 \
       <==> K_1 K_2 &< (2 m l_(1 2)) / gamma.
     $
     On a ring of $r >= 3$ identical $p$-cliques $K_p$ ($p >= 3$) joined by single bridges ($l_(1 2) = 1$, $K_c = p(p - 1) + 2$, $2 m = r K_c$), merging adjacent cliques strictly increases $Q_gamma$ iff
     $
       r > gamma (p(p - 1) + 2) quad ("at" gamma = 1, p = 5: r > 22).
     $
     Hence every fixed-$gamma$ modularity optimizer fails to resolve modules of internal edge count $l < sqrt(m / (2 gamma))$.
  2. *(Resolution-Free Sub-Partition Invariance of Leiden-CPM \#108 @traag2019leiden):* Under the Constant Potts Model
     $
       cal(H)_"CPM" (cal(C)) := sum_(r=1)^K ( e_r - gamma binom(|C_r|, 2) ),
     $
     merging $C_1, C_2$ improves $cal(H)_"CPM"$ iff $e(C_1, C_2) / (|C_1| |C_2|) > gamma$—strictly independent of $n, m$. Any globally optimal partition $cal(C)^*$ satisfies *sub-partition invariance*: for any subset of communities $S = union.big_(r in I) C_r^*$, ${C_r^* : r in I}$ remains globally optimal on $G[S]$, resolving the $r$-clique ring for *all* $r >= 3$ whenever $1 / p^2 < gamma < 1$.
]

#proof[
  *Part 1:* Subtracting pre-merge from post-merge contributions yields
  $
    Delta Q_gamma (C_1, C_2) &= (l_(1 2)) / m - gamma ((K_1 + K_2)^2 - K_1^2 - K_2^2) / (4 m^2) \
    &= (l_(1 2)) / m - gamma (K_1 K_2) / (2 m^2).
  $
  On the $r$-clique ring, $K_1 = K_2 = K_c = p(p - 1) + 2$, $l_(1 2) = 1$, and $2 m = r K_c$, so
  $
    Delta Q_gamma = 1 / m ( 1 - gamma (K_c^2) / (2 m) ) = 1 / m ( 1 - gamma (K_c) / r ) > 0
  $
  iff $r > gamma K_c$ ($r > 22$ for $K_5$ at $gamma = 1$). For balanced modules with $K_1 approx K_2 approx 2 l$ and $l_(1 2) = 1$, $4 l^2 < 2 m / gamma <==> l < sqrt(m / (2 gamma))$.

  *Part 2:* Merging $C_1, C_2$ alters $cal(H)_"CPM"$ by
  $
    Delta cal(H)_"CPM" &= e(C_1, C_2) - gamma [ binom(|C_1| + |C_2|, 2) - binom(|C_1|, 2) - binom(|C_2|, 2) ] \
    &= e(C_1, C_2) - gamma |C_1| |C_2|,
  $
  which equals $1 - gamma p^2 < 0$ on adjacent $p$-cliques for every $gamma in (1 / p^2, 1)$ and *every* $r >= 3$. Because $cal(H)_"CPM"$ is strictly additive across disjoint unions with zero global $2 m$ coupling, replacing ${C_r^* : r in I}$ on $G[S]$ by any higher-scoring sub-partition would strictly increase global $cal(H)_"CPM"$ on $G$, proving sub-partition invariance.
]

Sweeping $gamma$ in Spinglass (\#117) or Louvain (\#107) shifts the $sqrt(m / (2 gamma))$ threshold but cannot resolve multi-scale graphs with tiny dense cliques and large sparse clusters simultaneously; only local-density CPM (Leiden \#108) eliminates the resolution limit outright.

#subsec("C", "Louvain's Disconnected-Community Defect vs. Leiden Refinement & CNM Skew")

Traag, Waltman, and van Eck @traag2019leiden revealed that iterative Louvain can leave up to $16\%$ of communities *internally disconnected*. We establish this defect via a 14-vertex constructive counterexample, prove Leiden's exact connectedness invariant, and explain CNM's super-hub scaling blowup.

#theorem(title: "Louvain/Fast Unfolding Disconnected-Community Counterexample vs. Leiden Connectedness Guarantee", origin: "Original Counterexample")[
  1. *(14-Vertex Counterexample Against Louvain \#107 & Fast Unfolding \#129):* Let $G$ consist of three disjoint cliques $A = K_4$, $B = K_4$, $D = K_5$ plus one bridge vertex $v_0$ connected by $1$ edge to $A$, $1$ edge to $B$, and $3$ edges to $D$ ($n = 14, m = 27$). From the connected state $cal(P)_0 = {C_0, D}$ ($C_0 = A union B union {v_0}$), Phase 1 moves $v_0$ from $C_0$ to $D$ with strictly positive gain
     $
       Delta Q(v_0: C_0 -> D) = 69 / 1458 approx +0.0473 > 0,
     $
     halting at the *strict Phase-1 local optimum* $cal(P)_1 = {A union B, D union {v_0}}$ where $A union B$ has *zero edges between $A$ and $B$*! Phase 2 collapses $A union B$ into one super-node, permanently trapping the disconnected community.
  2. *(Exact Connectedness & $gamma$-Separation of Leiden \#108 @traag2019leiden):* Inserting a *Refinement Phase* between Local Moving and Aggregation—merging only $gamma$-well-connected sub-communities inside each $C in cal(P)$—guarantees every converged Leiden community $C$ is *strictly internally connected* and obeys for every proper non-empty $S subset C$:
     $
       e(S, C without S) >= gamma (K_S K_(C without S)) / (2 m) \
       ("or" e(S, C without S) >= gamma |S| |C without S| "for CPM"),
     $
     yielding asymptotic sub-partition optimality and faster convergence via a neighbor-restricted fast local-moving queue.
  3. *(CNM \#115 Dendrogram-Skew Blowup @clauset2004cnm):* Greedy max-$Delta Q$ agglomeration repeatedly absorbs leaves into one growing super-community on heavy-tailed graphs, forcing unbalanced dendrogram depth $d = Theta(n)$ and degrading CNM's $O(m d log n)$ time to $O(m n log n)$, versus $O(L m)$ for Louvain/Leiden.
]

#proof[
  *Part 1:* Here $m = 6 + 6 + 10 + 5 = 27$ ($2 m = 54, 2 m^2 = 1458$), $C_0' := A union B$ has $K_(C_0') = 26$, $K_D = 23$, $k_(v_0) = 5$, $k_(v_0, C_0') = 2$, and $k_(v_0, D) = 3$. Moving $v_0$ to $D$ changes $Q_1$ by
  $
    Delta Q(v_0: C_0 -> D) &= (k_(v_0, D) - k_(v_0, C_0')) / m \
    &quad - (k_(v_0) (K_D - K_(C_0'))) / (2 m^2) \
    &= 1 / 27 - (5 (-3)) / 1458 = 69 / 1458 > 0.
  $
  In $cal(P)_1 = {C_0', D'}$ ($K_(D') = 28$), no single-vertex move improves $Q_1$: for $a_1 in A$ linked to $v_0$ ($k_(a_1) = 4, k_(a_1, C_0' without {a_1}) = 3, k_(a_1, D') = 1$), moving to $D'$ gives $-2/27 - (4(28 - 22))/1458 = -132/1458 < 0$ and isolating gives $-3/27 + (4 dot 22)/1458 = -74/1458 < 0$; every other $a in A$ (or $b in B$) has $k_(a, D') = 0$ and $3$ internal edges, giving $-3/27 + (3 dot 23)/1458 = -93/1458 < 0$ on isolation; and every vertex in $D'$ has $>= 3$ internal edges vs. $<= 2$ external edges. Thus $cal(P)_1$ is a strict Phase-1 fixed point with $e(A, B) = 0$, which Phase 2 seals permanently into a single super-node!

  *Part 2:* Leiden's refinement initializes singletons within each $C in cal(P)$ and merges sub-communities $S, T subset C$ *only if* $e(S, T) > 0$, $Delta Q_gamma (S, T) > 0$, and both are $gamma$-well-connected to $C$. On $C_0' = A union B$, $e(A, B) = 0$ forbids merging $A$ and $B$, splitting them cleanly into separate super-nodes! Induction on merges with $e(S, T) > 0$ guarantees strict internal connectedness, and convergence ($cal(P)_"agg" = cal(P)_"ref"$) ensures no subset $S subset C$ can be split off with $Delta Q_gamma > 0$, rearranging directly to the stated cut bound @traag2019leiden.

  *Part 3:* Absorbing a degree-$k_v$ leaf into a growing hub community $C$ yields $Delta Q = 1 / m - k_v K_C / (2 m^2) > 0$ until $K_C >= 2 m / k_v$, producing a caterpillar merge tree of depth $d = Theta(n)$ and $O(m n log n)$ heap maintenance @clauset2004cnm, whereas multilevel local moving contracts vertices geometrically in $L = O(log n)$ passes.
]

Theorem 5 proves strict algorithmic dominance inside Section II: *Leiden (\#108) strictly dominates Louvain (\#107), Fast Unfolding (\#129), and CNM (\#115)* across connectedness, quality, speed, and CPM resolution-limit immunity.

#wide_figure(
  caption: [*Table II:* Head-to-Head Comparison of Modularity, Potts & Multilevel Optimizers (Section II: Algs. 107, 108, 115, 116, 117, 129).],
  table(
    columns: (1.06fr, 0.96fr, 0.92fr, 0.86fr, 0.88fr, 1.32fr),
    align: (left, left, center, center, center, left),
    [Algorithm], [Objective & Strategy], [Time Cost], [Res. Limit?], [Connected?], [Primary Strength / Failure Mode],
    [CNM (\#115) @clauset2004cnm], [Agglom. greedy $Q_1$], [$O(m d log n)$], [Yes ($sqrt(m/2)$)], [Yes (merges)], [Simple heap; $d = Theta(n)$ hub snowball],
    [Lead. Eig. (\#116) @newman2006leading], [Spectral div. $B^((g))$], [$O(n(m + n))$], [Yes ($sqrt(m/2)$)], [Not guar.], [Principled $B^((g))$ stop; rank-1 trap],
    [Spinglass (\#117) @reichardt2006spinglass], [Potts $cal(H)_"RB"$ SA], [$O(N_"sw" (m + n q))$], [Tunable $gamma$], [Not guar.], [Signed/multi-scale; slow SA, needs $q$],
    [Fast Unf. (\#129) @blondel2008louvain], [2-phase aggr. $Q_1$], [$O(L T m)$], [Yes ($sqrt(m/2)$)], [*No* (Thm. 5)], [Original Blondel baseline; disconnected],
    [Louvain (\#107) @blondel2008louvain], [2-phase aggr. $Q_gamma$], [$O(L m)$], [Tunable $gamma$], [*No* (Thm. 5)], [Fast massive graphs; bridge disconnect],
    [Leiden (\#108) @traag2019leiden], [*3-phase ref.* $Q_gamma$ / CPM], [$O(L m)$ (queue)], [*Free* (CPM)], [*Provably Yes*], [Strictly dominates Louvain; gold std.],
  ),
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Random-Walk, Flow, Label & Message-Passing Dynamics")

Seven Category V algorithms replace static null contrasts with dynamical processes on the graph: *Infomap* (\#113) @rosvall2008infomap, *Walktrap* (\#114) @pons2006walktrap, *Markov Clustering* (MCL, \#109) @vandongen2000mcl, *Label Propagation* (LPA, \#110) @raghavan2007lpa, *Chinese Whispers* (\#111) @biemann2006chinese, *Fluid Communities* (FluidC, \#119) @pares2018fluidc, and *Affinity Propagation* (AP, \#106) @frey2007affinity.

#subsec("A", "Flow Compression (Infomap), Diffusion Ward (Walktrap) & Inflation (MCL)")

Infomap (\#113), Walktrap (\#114), and MCL (\#109) all begin from the Markov walk $P = D^(-1) W$ (or PageRank teleportation on directed graphs in Infomap), yet extract structure via MDL source coding, $t$-step spectral diffusion distances, and nonlinear column-entropy inflation respectively.

#theorem(title: "Map Equation Field-of-View Bound, Walktrap Spectral Identity & MCL Entropy Contraction", origin: "Original Unified Proof")[
  Let $G = (V, E, W)$ be connected undirected with $P = D^(-1) W$, $pi_u = d_u / (2 m)$, and $D$-orthonormal right eigenpairs $(lambda_alpha, u_alpha)$ ($1 = lambda_1 > lambda_2 >= dots >= lambda_n$).
  1. *(Infomap \#113 Map Equation & Exponential Resolution Advantage @rosvall2008infomap):* With module exit flows $q_i^"out" = sum_(u in M_i, v in.not M_i) pi_u P_(u v)$, $q^"out" = sum_i q_i^"out"$, and $p_i^"use" = q_i^"out" + sum_(u in M_i) pi_u$, the two-level random-walk code length is
     $
       L(cal(M)) = q^"out" H(cal(Q)) + sum_(i=1)^K p_i^"use" H(cal(P)^i).
     $
     On the $r$-clique ring of Theorem 4 ($K_c = p(p - 1) + 2, 2 m = r K_c$), merging two adjacent cliques changes $L$ by
     $
       Delta L_"pair" &= 2 / (r K_c) [ (K_c + 1) log_2 (1 + K_c / (K_c + 2)) \
       &quad - log_2 (r / 2) + O(1) ].
     $
     Thus Infomap defeats modularity's $O(p^2)$ resolution limit by an *exponential* factor—resolving cliques up to $r = 2^(Omega(K_c))$—yet exhibits its own *field-of-view limit* on sparse low-$K_c$ chains/rings once index overhead $log_2 r$ dominates!
  2. *(Exact Spectral Diffusion Identity of Walktrap \#114 @pons2006walktrap):* Walktrap's $t$-step distance $r_(i j) (t) := sqrt(sum_k (P_(i k)^t - P_(j k)^t)^2 / d_k)$ obeys identically
     $
       r_(i j)^2 (t) = sum_(alpha=2)^n lambda_alpha^(2 t) ( u_alpha (i) - u_alpha (j) )^2.
     $
     Hence Walktrap is *algebraically identical to adjacency-constrained Ward clustering* @ward1963hierarchical in the diffusion-map embedding $Phi_t (i) = (lambda_2^t u_2 (i), dots, lambda_n^t u_n (i))$, requiring $tau_"in" << t << (1 - lambda_2)^(-1)$.
  3. *(Strict Column-Entropy Contraction of MCL \#109 @vandongen2000mcl):* For any non-uniform stochastic column $p in Delta^(n-1)$ and inflation $r > 1$, $Gamma_r (p)_i := p_i^r / (sum_l p_l^r)$ satisfies
     $
       (dif) / (dif r) H(Gamma_r (p)) = - r dot "Var"_(Gamma_r (p)) (ln p) < 0,
     $
     amplifying ratios as $(p_i / p_j)^r$ and driving $M <- Gamma_r (M^2)$ toward sparse doubly idempotent cluster attractors.
]

#proof[
  *Part 1:* Setting $h(x) := x log_2 x$, expanding $L(cal(M))$ gives
  $
    L(cal(M)) &= h(q^"out") - 2 sum_(i=1)^K h(q_i^"out") \
    &quad + sum_(i=1)^K h(p_i^"use") - sum_(u=1)^n h(pi_u),
  $
  whose last term is partition-independent. On the $r$-clique ring, $q_i^"out" = 2 / (r K_c)$, $q^"out" = 2 / K_c$, and $p_i^"use" = (K_c + 2) / (r K_c)$. Merging an adjacent pair reduces $q^"out"$ by $2 / (r K_c)$, eliminates one $q_i^"out"$ term (keeping $2$ external bridges), and replaces two $p_i^"use"$ terms by $(2 K_c + 2) / (r K_c)$. Expanding $h(q^"out" - 2 / (r K_c)) - h(q^"out") = -(2 / (r K_c))(log_2 (2 / K_c) + 1 / ln 2) + O(r^(-2))$ and combining with $+2 h(2 / (r K_c))$ and $h((2 K_c + 2) / (r K_c)) - 2 h((K_c + 2) / (r K_c))$ cancels $log_2 (r K_c)$ against intra-module expansion, leaving $(2 / (r K_c)) [ (K_c + 1) log_2 ((2 K_c + 2) / (K_c + 2)) - log_2 (r / 2) + O(1) ]$, which stays positive until $r = 2^(Omega(K_c))$.

  *Part 2:* From $P = D^(-1/2) (sum_alpha lambda_alpha v_alpha v_alpha^top) D^(1/2)$ and $u_alpha = D^(-1/2) v_alpha$, $P_(i k)^t = sum_(alpha=1)^n lambda_alpha^t u_alpha (i) u_alpha (k) d_k$. The constant $alpha = 1$ term cancels in $P_(i k)^t - P_(j k)^t$, so
  $
    r_(i j)^2 (t) &= sum_(k=1)^n 1 / d_k ( sum_(alpha=2)^n lambda_alpha^t (u_alpha (i) - u_alpha (j)) u_alpha (k) d_k )^2 \
    &= sum_(alpha, beta = 2)^n lambda_alpha^t lambda_beta^t (u_alpha (i) - u_alpha (j)) \
    &quad times (u_beta (i) - u_beta (j)) ( u_alpha^top D u_beta ) \
    &= sum_(alpha=2)^n lambda_alpha^(2 t) ( u_alpha (i) - u_alpha (j) )^2.
  $

  *Part 3:* With cumulant generator $psi(r) := ln(sum_(p_i > 0) p_i^r)$, $Gamma_r (p)_i = exp(r ln p_i - psi(r))$ has $psi'(r) = bb(E)_(Gamma_r)[ln p]$ and $psi''(r) = "Var"_(Gamma_r)(ln p) > 0$ (non-uniform $p$). Since $H(Gamma_r (p)) = psi(r) - r psi'(r)$, differentiating gives
  $
    (dif) / (dif r) H(Gamma_r (p)) = - r psi''(r) = - r dot "Var"_(Gamma_r (p)) (ln p) < 0,
  $
  driving columns strictly toward zero-entropy indicator attractors @vandongen2000mcl.
]

#subsec("B", "Label Propagation (LPA, CW, FluidC) vs. Factor-Graph BP (AP)")

Three $O(T m)$ local voting algorithms—*LPA* (\#110) @raghavan2007lpa, *Chinese Whispers* (\#111) @biemann2006chinese, and *FluidC* (\#119) @pares2018fluidc—contrast starkly with $O(T n^2)$ factor-graph *Affinity Propagation* (\#106) @frey2007affinity.

#theorem(title: "LPA/CW Oscillation & Monster Collapse vs. FluidC Repulsion & AP Monotonicity", origin: "Original Unified Proof")[
  1. *(Oscillation, Asynchronous Lyapunov Convergence & Monster Optimum of LPA \#110 & CW \#111):*
     - On any connected bipartite graph initialized with label $1$ on $U$ and $2$ on $V$, *synchronous* LPA cycles permanently with period $2$ @raghavan2007lpa.
     - *Asynchronous* single-vertex LPA (\#110) and `top` Chinese Whispers (\#111) perform monotone ascent on ferromagnetic agreement energy $E(bold(c)) := (1 / 2) sum_(i, j) W_(i j) bold(1)[c_i = c_j]$, terminating finitely; *however*, the unique global maximum of $E$ on connected $G$ is the *trivial single monster community* $c_1 = dots = c_n$ ($E = m$), causing percolation collapse without degree damping (`lin`/`log` in CW \#111 @biemann2006chinese).
  2. *(Exact Monster Instability of FluidC \#119 @pares2018fluidc):* FluidC updates $K$ communities asynchronously with inverse-size density $rho(c) = 1 / abs(C_c)$ via
     $
       c_i^((t+1)) in arg max_(c in {1, dots, K}) rho(c) sum_(j in N(i) union {i}) bold(1)[c_j = c].
     $
     On any connected $d$-regular graph with $n > d + K - 1$, the monster state $|C_1| = n - K + 1$, $|C_2| = dots = |C_K| = 1$ is *strictly unstable*, forcing boundary nodes to leave $C_1$ and driving balanced equilibrium $|C_c| = Theta(n / K)$.
  3. *(Affinity Propagation \#106 Exact Facility-Location Program & Monotonicity @frey2007affinity):* AP's responsibility--availability recursions
     $
       r(i, k) &<- s(i, k) - max_(k' != k) { a(i, k') + s(i, k') }, \
       a(i, k) &<- min( 0, r(k, k) + sum_(i' in.not {i, k}) max(0, r(i', k)) ), \
       a(k, k) &<- sum_(i' != k) max(0, r(i', k))
     $
     execute exact max-sum belief propagation on the binary factor graph for exemplar facility location:
     $
       max_(c in {0, 1}^(n times n)) sum_(i, k) s(i, k) c_(i k) \
       "s.t." quad sum_(k=1)^n c_(i k) = 1, quad c_(i k) <= c_(k k).
     $
     Under uniform diagonal preference $s(k, k) = p$, optimal exemplar count $K(c^*(p)) = sum_k c_(k k)^*$ is *monotonically non-decreasing* in $p$: $p_1 < p_2 ==> K(c^*(p_1)) <= K(c^*(p_2))$.
]

#proof[
  *Part 1:* Synchronous bipartite updating swaps $U$ and $V$ labels every step forever. Under asynchronous update of $i$ from $c_i$ to $c'_i$,
  $
    E(bold(c)') - E(bold(c)) &= sum_(j in N(i)) W_(i j) ( bold(1)[c_j = c'_i] - bold(1)[c_j = c_i] ) >= 0,
  $
  forcing finite termination since $E <= m$ takes finitely many values; $E = m$ uniquely when $c$ is constant on connected $G$.

  *Part 2:* In the monster state $|C_1| = n - K + 1$, connectivity guarantees an edge between some $i in C_1$ and a singleton $j in C_r$ ($|C_r| = 1$). Inside $N(i) union {i}$, at most $d$ vertices belong to $C_1$ (since $j in.not C_1$), giving $C_1$ score $<= d / (n - K + 1) < 1 = rho(r) dot 1$ whenever $n > d + K - 1$, immediately breaking the monster state.

  *Part 3:* Write $F(c; p) = S_"off" (c) + p K(c)$ where $S_"off" (c) = sum_(i != k) s(i, k) c_(i k)$. Optimality of $c_1, c_2$ at $p_1 < p_2$ gives
  $
    S_"off" (c_1) + p_1 K(c_1) &>= S_"off" (c_2) + p_1 K(c_2), \
    S_"off" (c_2) + p_2 K(c_2) &>= S_"off" (c_1) + p_2 K(c_1).
  $
  Adding cancels $S_"off"$, yielding $(p_2 - p_1)(K(c_2) - K(c_1)) >= 0$, hence $K(c_1) <= K(c_2)$. Factor-graph max-sum message reduction follows @frey2007affinity, with damping $lambda in [0.5, 0.9]$ preventing loopy oscillations.
]

AP (\#106) uniquely supports arbitrary signed, asymmetric, non-metric similarities and returns real data-point exemplars at $O(T n^2)$ cost, whereas LPA/CW (\#110/\#111) deliver unmatched $O(T m)$ speed for billion-edge graphs and FluidC (\#119) prevents monster collapse when $K$ is fixed.

#wide_figure(
  caption: [*Table III:* Head-to-Head Comparison of Random-Walk, Flow, Label & Message-Passing Dynamics (Section III: Algs. 106, 109--111, 113--114, 119).],
  table(
    columns: (1.04fr, 0.98fr, 0.92fr, 0.72fr, 0.86fr, 1.48fr),
    align: (left, left, center, center, center, left),
    [Algorithm], [Dynamical Mechanism], [Time Cost], [Auto-$K$], [Directed / Non-Met.], [Primary Strength / Failure Mode],
    [AP (\#106) @frey2007affinity], [Max-sum BP exemplars], [$O(T n^2)$], [Via pref. $p$], [*Yes / Yes*], [Non-metric exemplars; $O(n^2)$ & cycles],
    [MCL (\#109) @vandongen2000mcl], [Expand $M^2$, inflate $Gamma_r$], [$O(T n k^2)$], [Via $r > 1$], [Yes / No], [Bio/PPI gold std.; self-loop & $r$ tuning],
    [LPA (\#110) @raghavan2007lpa], [Async. neighbor mode], [$O(T m)$], [*Yes*], [Weak / No], [Billion-edge speed; monster collapse],
    [CW (\#111) @biemann2006chinese], [Weighted async. mode], [$O(T m)$], [*Yes*], [Weak / No], [Hub-safe (`lin`/`log`) NLP; non-determ.],
    [Infomap (\#113) @rosvall2008infomap], [2-level MDL Map Eq.], [$O(m log n)$], [*Yes*], [*Native* / No], [Beats $Q$ res. limit; sparse field-of-view],
    [Walktrap (\#114) @pons2006walktrap], [Ward on $t$-step $r_(i j)$], [$O(m n H)$], [Dendro.], [No / No], [Exact diffusion Ward; high-diam. fail],
    [FluidC (\#119) @pares2018fluidc], [$1 / abs(C_c)$ fluid voting], [$O(T(n + m))$], [No (fixed $K$)], [No / No], [Monster-proof linear; needs $K$ & conn.],
  ),
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Divisive Cuts, Multilevel Partitioning, Overlapping & Generative Block Models")

The remaining twelve Category V algorithms span divisive combinatorial cuts (*Girvan--Newman* \#112 @girvan2002community, *HCS* \#118 @hartuv2000hcs), balanced multilevel partitioning (*Kernighan--Lin* \#120 @kernighan1970efficient, *METIS* \#121 @karypis1998metis, *Chaco* \#122 @hendrickson1995chaco), overlapping and local statistical covers (*OSLOM* \#124 @lancichinetti2011oslom, *BigCLAM* \#125 @yang2013bigclam, *SLPA* \#126 @xie2011slpa, *CPM* \#127 @palla2005cpm, *DEMON* \#128 @coscia2012demon), and generative stochastic block models (*SBM* \#130 @holland1983sbm @abbe2018sbmreview, *DC-SBM* \#131 @karrer2011dcsbm).

#subsec("A", "Divisive Cuts (GN, HCS) vs. Multilevel Partitioning (KL, METIS, Chaco)")

#theorem(title: "HCS Diameter-2 Invariant, GN Bottleneck & Exact Multilevel Cut-Preservation", origin: "Original Unified Proof")[
  1. *(HCS \#118 Exact Diameter-2 Bound vs. Girvan--Newman \#112):* Every highly connected cluster $H$ accepted by *HCS (\#118)* ($kappa'(H) > n_H / 2$, @hartuv2000hcs) satisfies simultaneously
     $
       "diam"(H) <= 2, quad delta(H) > n_H / 2, quad |E(H)| > (n_H^2) / 4.
     $
     Hence HCS guarantees ultra-compact diameter-2 clusters in $O(N_"cut" (m n + n^2 log n))$ Stoer--Wagner min-cut time, yet provably shatters every sparse community with $overline(d)_H <= n_H / 2$; conversely, *Girvan--Newman (\#112)* @girvan2002community removes maximum edge-betweenness bridges in $O(m^2 n)$ time, cleanly separating sparse high-diameter modules at cubic-to-quartic cost.
  2. *(Exact Multilevel Cut-Invariance of METIS \#121 & Chaco \#122 over Kernighan--Lin \#120):* Under any vertex-contraction coarsening hierarchy $G_0 -> dots -> G_L$ with summed vertex and parallel-edge weights (Heavy-Edge Matching in METIS \#121 @karypis1998metis and multilevel Chaco \#122 @hendrickson1995chaco), projecting any balanced $K$-partition $cal(P)_L$ back to $G_0$ preserves edge cut and part weights identically:
     $
       "cut"_(G_0) (cal(P)_0) &= "cut"_(G_L) (cal(P)_L), \
       w_(G_0) (V_r^((0))) &= w_(G_L) (V_r^((L))).
     $
     Uncoarsening with lookahead boundary refinement (*Kernighan--Lin \#120* @kernighan1970efficient or linear-time Fiduccia--Mattheyses) monotonically reduces cut cost while accelerating KL's flat $O(P n^2 log n)$ search to $O((n + m) log K)$ (METIS) or $O((m + n log n) log K)$ (Chaco multilevel spectral bisection/quadrisection).
]

#proof[
  *Part 1:* Since $delta(H) >= kappa'(H) > n_H / 2$, $|E(H)| = (1 / 2) sum_u "deg"_H (u) > n_H^2 / 4$. For any non-adjacent $u, v in V(H)$, $u, v in.not N_H (u) union N_H (v)$ implies $|N_H (u) union N_H (v)| <= n_H - 2$, so by inclusion--exclusion:
  $
    |N_H (u) inter N_H (v)| &= |N_H (u)| + |N_H (v)| - |N_H (u) union N_H (v)| \
    &> n_H / 2 + n_H / 2 - (n_H - 2) = 2,
  $
  guaranteeing at least $3$ common neighbors and $"diam"(H) <= 2$. Any sparse module with $|E(H)| <= n_H^2 / 4$ has $kappa'(H) <= delta(H) <= n_H / 2$ and is forcibly split.

  *Part 2:* For one coarsening surjection $pi: V_(ell) -> V_(ell + 1)$ and projected labels $c_(ell)(u) = c_(ell+1)(pi(u))$,
  $
    w_(G_(ell)) (V_r^((ell))) &= sum_(U : c_(ell+1)(U) = r) sum_(u in U) w(u) = w_(G_(ell+1)) (V_r^((ell+1))).
  $
  Edges inside contracted super-nodes are never cut, while cross-edges sum identically to super-edge weights:
  $
    "cut"_(G_(ell)) (cal(P)_(ell)) &= sum_(c_(ell+1)(U) != c_(ell+1)(U')) sum_(u in U, v in U') W(u, v) \
    &= sum_(c_(ell+1)(U) != c_(ell+1)(U')) W(U, U') = "cut"_(G_(ell+1)) (cal(P)_(ell+1)).
  $
  Induction from $L - 1$ down to $0$ yields cut/balance invariance; KL/FM lookahead prefix swaps apply only positive cumulative gains $G_k = max_t sum_(i=1)^t g_i > 0$, strictly never worsening the cut across $L = O(log n)$ geometric HEM levels @karypis1998metis.
]

#subsec("B", "Overlapping & Local Statistical Communities (OSLOM, BigCLAM, SLPA, CPM, DEMON)")

Empirical networks frequently exhibit *dense overlaps* where vertices sharing multiple affiliations connect *more* densely than non-overlapping cores @yang2013bigclam. Five Category V methods address overlapping covers through orthogonal principles.

#theorem(title: "BigCLAM Dense-Overlap Monotonicity, CPM Percolation Threshold & OSLOM Order-Statistic Shield", origin: "Original Unified Proof")[
  1. *(Dense-Overlap Inversion & Linear-Time Cached Gradient of BigCLAM \#125 @yang2013bigclam):* Under nonnegative affiliations $F in bb(R)_(>= 0)^(n times K)$ with $p_(u v) = 1 - exp(-chevron.l F_u, F_v chevron.r)$, vertices sharing $s >= 1$ communities (strength $alpha > 0$) obey strictly increasing edge probability $p_(u v) (s) = 1 - e^(-s alpha^2)$ (reversing sparse-boundary cut assumptions). Caching $S := sum_(w=1)^n F_w$ factors the log-likelihood gradient as
     $
       nabla_(F_u) ell(F) &= sum_(v in N(u)) (exp(-chevron.l F_u, F_v chevron.r)) / (1 - exp(-chevron.l F_u, F_v chevron.r)) F_v \
       &quad - ( S - F_u - sum_(v in N(u)) F_v ),
     $
     slashing per-epoch coordinate ascent from $O(K n^2)$ to *strictly linear* $O(K (m + n))$!
  2. *(Exact Critical Percolation Threshold of CPM \#127 @palla2005cpm):* Rolling $k$-cliques across shared $(k - 1)$-faces ($k >= 3$) on $G(n, p)$ triggers a giant percolation web at the exact threshold
     $
       p_c (k) = [ (k - 1) n ]^(-1 / (k - 1)) quad (n -> infinity).
     $
     Above $p_c (k)$, CPM merges the graph into a single giant component; on sparse or triangle-free communities ($p < (2 n)^(-1/2)$), CPM detects nothing while risking $O(d n 3^(d/3))$ clique-enumeration cost.
  3. *(OSLOM \#124 Order-Statistic Null Control vs. SLPA \#126 & DEMON \#128):*
     - In *OSLOM (\#124)* @lancichinetti2011oslom, scoring the minimum configuration-model $p$-value $r_((1))$ among $N_"ext"$ external neighbors via the exact extreme-order-statistic CDF
       $
         Omega(r_((1))) = 1 - (1 - r_((1)))^(N_"ext")
       $
       (or Beta CDF $I_(r_((j))) (j, N_"ext" - j + 1)$) prevents multiple-testing hallucination, outputting *zero* false clusters on Erdős--Rényi noise (unlike modularity, LPA, and BigCLAM).
     - *SLPA (\#126)* @xie2011slpa accumulates $T$-step label memory in $O(T m)$ time, thresholding posterior frequency $p(c | i) >= r$ ($r <= 0.5$ overlapping; $r > 0.5$ disjoint), whereas *DEMON (\#128)* @coscia2012demon runs local LPA inside each ego-minus-ego subgraph $G[N(v)]$ and merges modules above overlap $epsilon$, scaling as $O(T tau(G) + R M overline(s))$ in total triangles $tau(G)$ (locally resolution-free on social graphs, yet blind to edges outside triangles).
]

#proof[
  *Part 1:* Shared support of size $s$ gives $chevron.l F_u, F_v chevron.r = s alpha^2$ and $p_(u v) (s) = 1 - e^(-s alpha^2)$, strictly increasing in $s$. Differentiating
  $
    ell(F) &= sum_((u, v) in E) log(1 - e^(-chevron.l F_u, F_v chevron.r)) \
    &quad - sum_((u, v) in.not E, u != v) chevron.l F_u, F_v chevron.r
  $
  with respect to $F_u$ gives a non-neighbor sum $-sum_(v in.not N(u) union {u}) F_v = -(S - F_u - sum_(v in N(u)) F_v)$, computable in $O(K(d_u + 1))$ per vertex and $O(K(m + n))$ per epoch @yang2013bigclam.

  *Part 2:* On $G(n, p)$, entering a $k$-clique via one $(k - 1)$-face leaves $k - 1$ unused $(k - 1)$-faces, each extended by any of $n - k$ external vertices with probability $p^(k - 1)$. Setting expected Galton--Watson offspring $mu_k = (k - 1)(n - k) p^(k - 1) = 1$ yields $p_c (k) tilde.op [(k - 1) n]^(-1 / (k - 1))$ @palla2005cpm.

  *Part 3:* Under the null hypothesis with i.i.d. uniform mid-$p$ scores $R_1, dots, R_(N_"ext")$, $bb(P)(min_j R_j > r) = (1 - r)^(N_"ext")$, yielding exact CDF $Omega(r) = 1 - (1 - r)^(N_"ext")$ @lancichinetti2011oslom. In DEMON (\#128), edges in $G[N(v)]$ biject to triangles containing $v$, so $sum_v |E(G[N(v)])| = 3 tau(G)$, proving $O(T tau(G))$ local LPA cost and zero sensitivity to triangle-free bridges @coscia2012demon.
]

#subsec("C", "Generative Stochastic Block Models: Exact Recovery & DC-SBM Hub Invariance")

Unlike assortative cuts and modularity, *SBM* (\#130) @holland1983sbm @abbe2018sbmreview and *Degree-Corrected SBM* (DC-SBM, \#131) @karrer2011dcsbm infer arbitrary $K times K$ block affinities—capturing assortative, bipartite/disassortative, and core--periphery structure—while exhibiting sharp information-theoretic phase transitions and unifying Section II's modularity from first principles.

#theorem(title: "SBM Exact Recovery Threshold, Hub-Induced SBM Collapse & DC-SBM–Modularity Equivalence", origin: "Original Unified Proof")[
  1. *(Sharp Exact-Recovery & Weak-Detection Thresholds of Balanced 2-Block SBM \#130 @abbe2018sbmreview):*
     - In the logarithmic regime $p = a (ln n) / n$, $q = b (ln n) / n$ ($a > b > 0$), exact recovery (zero errors w.h.p.) is achievable iff the Chernoff--Hellinger condition holds:
       $
         | sqrt(a) - sqrt(b) | > sqrt(2).
       $
     - In the sparse constant-degree regime $p = a / n$, $q = b / n$, efficient weak detection succeeds above the Kesten--Stigum threshold $(a - b)^2 > 2 (a + b)$.
  2. *(Exact Hub-Induced Breakdown of SBM \#130 vs. Universal Invariance of DC-SBM \#131 @karrer2011dcsbm):* Consider a balanced 2-community graph (assortative ratio $c = omega_"in" / omega_"out" > 1$) where each true community $C_1, C_2$ contains half hubs ($theta_H = alpha > 1$) and half spokes ($theta_L = 1$). Let $h(x) := (x / (x + 1)) ln(2 x / (x + 1)) + (1 / (x + 1)) ln(2 / (x + 1)) in (0, ln 2)$ ($x > 1$). Comparing the true partition $g_"true" = (C_1, C_2)$ against the *wrong degree-sorted split* $g_"deg" = (H_1 union H_2, L_1 union L_2)$ yields the exact Poisson profile log-likelihood margins:
     $
       (ell_"SBM" (g_"deg") - ell_"SBM" (g_"true")) / (2 m) &= 2 h(alpha) - h(c), \
       (ell_"DC" (g_"true") - ell_"DC" (g_"deg")) / (2 m) &= h(c) > 0.
     $
     Hence for *every* hub ratio $alpha >= 9$ (where $2 h(alpha) > ln 2$), *vanilla SBM (\#130) strictly prefers splitting hubs from spokes over the true communities for ANY finite assortativity $c > 1$* (even $c = 10^6$!), whereas *DC-SBM (\#131)* recovers $g_"true"$ by margin $2 m h(c) > 0$ for *every* $alpha > 0$!
  3. *(Exact Equivalence of Planted DC-SBM \#131 and Modularity $Q_(gamma^*)$):* For a symmetric planted $K$-block Poisson DC-SBM with $p_(i j) = omega_(g_i g_j) k_i k_j / (2 m)$ ($omega_(r r) = omega_"in" > omega_"out" = omega_(r s)$),
     $
       ln P(A | g, omega, k) &= m ln(omega_"in" / omega_"out") Q_(gamma^*) (g) + C, \
       gamma^* &= (omega_"in" - omega_"out") / (ln omega_"in" - ln omega_"out"),
     $
     establishing generalized modularity as the exact MLE of symmetric DC-SBM at the logarithmic-mean resolution $gamma^*$!
]

#proof[
  *Part 1:* Exact MAP recovery fails if a vertex $v in C_1$ has more external neighbors $Y_v tilde.op "Bin"(n/2, q)$ than internal neighbors $X_v tilde.op "Bin"(n/2 - 1, p)$. Approximating by independent $X tilde.op "Pois"((a/2) ln n)$, $Y tilde.op "Pois"((b/2) ln n)$ and optimizing the Chernoff bound over $t > 0$:
  $
    bb(P)(Y >= X) &<= inf_(t > 0) exp( (ln n) / 2 [ b(e^t - 1) + a(e^(-t) - 1) ] ).
  $
  The unique minimizer $e^(t^*) = sqrt(a / b)$ attains $b(sqrt(a / b) - 1) + a(sqrt(b / a) - 1) = -(sqrt(a) - sqrt(b))^2$, giving $bb(P)(Y_v >= X_v) = n^(-(sqrt(a) - sqrt(b))^2 / 2 + o(1))$. Union bounding over $n$ vertices yields error $n^(1 - (sqrt(a) - sqrt(b))^2 / 2 + o(1)) -> 0$ iff $|sqrt(a) - sqrt(b)| > sqrt(2)$ (tight by second-moment necessity @abbe2018sbmreview).

  *Part 2:* Standard SBM maximizes $ell_"SBM" (g) = sum_(r, s) m_(r s) ln(m_(r s) / (n_r n_s))$ while DC-SBM maximizes $ell_"DC" (g) = sum_(r, s) m_(r s) ln(m_(r s) / (kappa_r kappa_s))$ @karrer2011dcsbm. Rescaling each equal subpopulation $H_1, L_1, H_2, L_2$ to size $1$ ($n_r = 2$ for both partitions):
  - For $g_"true"$: $m_(1 1) = m_(2 2) = c(alpha + 1)^2$, $m_(1 2) = m_(2 1) = (alpha + 1)^2$, $2 m = 2(c + 1)(alpha + 1)^2$, $kappa_1 = kappa_2 = m$. Thus
    $
      (ell_"SBM" (g_"true")) / (2 m) &= c / (c + 1) ln c + ln((alpha + 1)^2 / 4), \
      (ell_"DC" (g_"true")) / (2 m) &= h(c) - ln(2 m).
    $
  - For $g_"deg"$: $m_(H H) = 2(c + 1) alpha^2$, $m_(H L) = m_(L H) = 2(c + 1) alpha$, $m_(L L) = 2(c + 1)$, with $kappa_H = 2(c + 1) alpha(alpha + 1)$ and $kappa_L = 2(c + 1)(alpha + 1)$. Hence
    $
      (ell_"SBM" (g_"deg")) / (2 m) &= ln((c + 1) / 2) + (2 alpha) / (alpha + 1) ln alpha.
    $
    Subtracting gives $(ell_"SBM" (g_"deg") - ell_"SBM" (g_"true")) / (2 m) = 2 h(alpha) - h(c)$. Because $h(c) < ln 2$ for all $c > 1$ while $2 h(9) approx 0.736 > 0.693 = ln 2$, every $alpha >= 9$ causes SBM to prefer $g_"deg"$! Conversely, under $g_"deg"$, $m_(r s) / (kappa_r kappa_s) = (2 m)^(-1)$ for all $(r, s)$, so $(ell_"DC" (g_"deg")) / (2 m) = -ln(2 m)$, yielding exact DC-SBM margin $h(c) > 0$ independent of $alpha$.

  *Part 3:* Substituting $omega_(g_i g_j) = omega_"out" + (omega_"in" - omega_"out") delta(g_i, g_j)$ into the Poisson log-likelihood and collecting $g$-independent baselines into $C$:
  $
    ln P(A | g) &= 1 / 2 sum_(i, j) [ A_(i j) ln(omega_"in" / omega_"out") \
    &quad - (omega_"in" - omega_"out") (k_i k_j) / (2 m) ] delta(g_i, g_j) + C \
    &= m ln(omega_"in" / omega_"out") Q_(gamma^*) (g) + C,
  $
  with exact resolution $gamma^* = (omega_"in" - omega_"out") / (ln omega_"in" - ln omega_"out")$.
]

#wide_figure(
  caption: [*Table IV:* Head-to-Head Comparison of Divisive Cuts, Multilevel Partitioners, Overlapping & Generative Block Models (Section IV: Algs. 112, 118, 120--122, 124--128, 130--131).],
  table(
    columns: (1.06fr, 0.96fr, 0.94fr, 0.72fr, 0.86fr, 1.46fr),
    align: (left, left, center, center, center, left),
    [Algorithm], [Core Principle], [Time Cost], [Overlap?], [Non-Assort.?], [Primary Strength / Failure Mode],
    [Girvan--N. (\#112) @girvan2002community], [Divisive edge betw.], [$O(m^2 n)$], [No], [No], [Sparse bridges; cubic/quartic bottleneck],
    [HCS (\#118) @hartuv2000hcs], [Min-cut $kappa' > n/2$], [$O(N_"cut" (m n + n^2 log n))$], [No], [No], [Exact $"diam" <= 2$; shatters sparse nets],
    [Kernighan--L. (\#120) @kernighan1970efficient], [Lookahead pair swaps], [$O(P n^2 log n)$], [No], [No], [Escapes 1-swap traps; fails flat large $n$],
    [METIS (\#121) @karypis1998metis], [HEM multilevel $k$-way], [$O((n + m) log K)$], [No], [No], [Balanced HPC gold std.; needs $K$],
    [Chaco (\#122) @hendrickson1995chaco], [Multilevel spectral], [$O((m + n log n) log K)$], [No], [No], [Best FEM mesh cuts; balanced $K$ only],
    [OSLOM (\#124) @lancichinetti2011oslom], [Order-stat $p$-values], [$O(n^2)$], [*Yes*], [No], [Zero false clusters on ER; slower $O(n^2)$],
    [BigCLAM (\#125) @yang2013bigclam], [Nonneg. affil. $F F^top$], [$O(T K (m + n))$], [*Yes*], [No], [Dense overlaps at scale; needs $K$],
    [SLPA (\#126) @xie2011slpa], [Memory-based LPA], [$O(T m)$], [*Yes*], [No], [Linear-time overlap; thresh. $r$ & var.],
    [CPM (\#127) @palla2005cpm], [$k$-clique rolling], [$O(d n 3^(d/3) + c_k^2 k)$], [*Yes*], [No], [Exact clique motifs; blind w/o triangles],
    [DEMON (\#128) @coscia2012demon], [Ego-net LPA + merge], [$O(T tau(G) + R M overline(s))$], [*Yes*], [No], [Local res.-free social; blind w/o triangles],
    [SBM (\#130) @holland1983sbm], [Block Poisson MLE], [$O(T(m K + n K^2))$], [No], [*Yes*], [Exact phase theory; hub split ($alpha >= 9$)],
    [DC-SBM (\#131) @karrer2011dcsbm], [Deg.-corr. block MLE], [$O(T K(m + n))$], [No], [*Yes*], [Hub-invariant + non-assort.; needs $K$/MDL],
  ),
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 30-Algorithm Head-to-Head Matrix & Practitioner's Decision Guide")

Table V synthesizes *all thirty* graph, spectral, and network clustering algorithms of Category V (Algs. 103--132), followed by exact cross-algorithm reductions, dominance laws, a quantitative selection guide, and open frontiers.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Comparison Matrix of All 30 Graph, Spectral, and Network Algorithms (Category V, Algorithms 103--132). Group codes: *SP* = Spectral & Laplacian, *MD* = Modularity & Potts, *DY* = Dynamical / Walk / Label / BP, *CM* = Combinatorial & Multilevel Cuts, *OV* = Overlapping / Local Statistical, *GM* = Generative Block Models.],
  {
    set text(size: 6.35pt)
    set table(inset: (x: 2.3pt, y: 1.8pt))
    table(
      columns: (1.24fr, 0.34fr, 0.94fr, 1.04fr, 0.56fr, 0.58fr, 0.62fr, 1.38fr),
      align: (left, center, left, left, center, center, center, left),
      [Algorithm], [Grp.], [Core Objective / Op.], [Time Complexity], [Auto-$K$], [Overlap], [Hub-Safe], [Distinctive Strength / Failure Mode],
      [\#103 NJW @ng2002spectral], [SP], [$L_"sym"$ + row norm.], [$O(n^2 d + n^3)$], [No], [No], [*Yes*], [Consistent $K$-way cuts; dense $O(n^3)$],
      [\#104 Ncut @shi2000ncut], [SP], [$L_"rw"$ / gen. $(L, D)$], [$O(n^3)$ / Lanczos], [No], [No], [*Yes*], [Cheeger $sqrt(2 lambda_2)$ bound; heavy eigen],
      [\#105 Ratio Cut @hagen1992ratiocut], [SP], [Unnormalized $L$], [$O(n^3)$ / Lanczos], [No], [No], [*No*], [Cardinality balance; leaf-snip failure],
      [\#106 AP @frey2007affinity], [DY], [Max-sum BP exemplars], [$O(T n^2)$], [Via $p$], [No], [Yes], [Non-metric asymmetric; $O(n^2)$ & cycles],
      [\#107 Louvain @blondel2008louvain], [MD], [2-phase aggr. $Q_gamma$], [$O(L m)$], [*Yes*], [No], [Yes], [Fast massive nets; *disconnected defect*],
      [\#108 Leiden @traag2019leiden], [MD], [3-phase ref. $Q_gamma$ / CPM], [$O(L m)$], [*Yes*], [No], [*Yes*], [*Connected + res.-free CPM*; top default],
      [\#109 MCL @vandongen2000mcl], [DY], [Expand $M^2$, inflate $Gamma_r$], [$O(T n k^2)$], [Via $r$], [No], [Yes], [Bio/PPI clusters; $r$ & self-loop sensitive],
      [\#110 LPA @raghavan2007lpa], [DY], [Async. neighbor mode], [$O(T m)$], [*Yes*], [No], [Weak], [Fastest linear time; monster collapse],
      [\#111 Chinese Wh. @biemann2006chinese], [DY], [Deg.-wt. async. mode], [$O(T m)$], [*Yes*], [No], [Med], [NLP hub-damped (`lin`/`log`); stochastic],
      [\#112 Girvan--N. @girvan2002community], [CM], [Div. edge betweenness], [$O(m^2 n)$], [Via $Q$], [No], [Yes], [Interpretable bridges; $O(m^2 n)$ slow],
      [\#113 Infomap @rosvall2008infomap], [DY], [2-level MDL Map Eq.], [$O(m log n)$], [*Yes*], [No], [*Yes*], [Directed flows + beats $Q$ res. limit],
      [\#114 Walktrap @pons2006walktrap], [DY], [Ward on $t$-step $r_(i j)$], [$O(m n H)$], [Dendro.], [No], [Yes], [Exact diffusion Ward; high-diam. fail],
      [\#115 CNM @clauset2004cnm], [MD], [Agglom. heap $Q_1$], [$O(m d log n)$], [*Yes*], [No], [Weak], [Exact merge tree; $d = Theta(n)$ hub blowup],
      [\#116 Lead. Eig. @newman2006leading], [MD], [Spectral div. $B^((g))$], [$O(n(m + n))$], [*Yes*], [No], [Yes], [Principled $B^((g))$ stop; rank-1 trap],
      [\#117 Spinglass @reichardt2006spinglass], [MD], [Potts $cal(H)_"RB"$ SA], [$O(N_"sw" (m + n q))$], [Up to $q$], [No], [Yes], [Signed & multi-scale $gamma$; slow SA],
      [\#118 HCS @hartuv2000hcs], [CM], [Rec. min-cut $kappa' > n/2$], [$O(N_"cut" (m n + n^2 log n))$], [*Yes*], [No], [No], [Guaranteed $"diam" <= 2$; shatters sparse],
      [\#119 FluidC @pares2018fluidc], [DY], [$1 / abs(C_c)$ fluid voting], [$O(T(n + m))$], [No], [No], [Yes], [Monster-proof linear; fixed $K$, conn. $G$],
      [\#120 Kernighan--L. @kernighan1970efficient], [CM], [Lookahead pair swaps], [$O(P n^2 log n)$], [No], [No], [N/A], [Escapes local traps; needs multilevel],
      [\#121 METIS @karypis1998metis], [CM], [HEM multilevel $k$-way], [$O((n + m) log K)$], [No], [No], [*Yes*], [Balanced HPC gold std.; fixed $K$],
      [\#122 Chaco @hendrickson1995chaco], [CM], [Multilevel spectral], [$O((m + n log n) log K)$], [No], [No], [*Yes*], [Top scientific FEM mesh cuts; fixed $K$],
      [\#123 Spec. Biclust. @kluger2003spectral], [SP], [Bistochastic SVD], [$O(n_r n_c min(n_r, n_c))$], [No], [No], [*Yes*], [Row+col checkerboard co-clustering],
      [\#124 OSLOM @lancichinetti2011oslom], [OV], [Order-stat $p$-values], [$O(n^2)$], [*Yes*], [*Yes*], [*Yes*], [Zero false ER clusters + outliers; slow],
      [\#125 BigCLAM @yang2013bigclam], [OV], [Nonneg. $1 - e^(-F_u^top F_v)$], [$O(T K (m + n))$], [No], [*Yes*], [*Yes*], [Dense overlap at scale; needs $K$],
      [\#126 SLPA @xie2011slpa], [OV], [Memory-based LPA], [$O(T m)$], [*Yes*], [*Yes*], [Med], [Linear-time overlap; thresh. $r$ sensitive],
      [\#127 CPM @palla2005cpm], [OV], [$k$-clique rolling], [$O(d n 3^(d/3) + c_k^2 k)$], [*Yes*], [*Yes*], [No], [Exact clique motifs; blind w/o triangles],
      [\#128 DEMON @coscia2012demon], [OV], [Ego-net LPA + merge], [$O(T tau(G) + R M overline(s))$], [*Yes*], [*Yes*], [Yes], [Local social overlap; blind w/o triangles],
      [\#129 Fast Unf. @blondel2008louvain], [MD], [Orig. 2-phase $Q_1$], [$O(L T m)$], [*Yes*], [No], [Yes], [Blondel baseline; disconnected defect],
      [\#130 SBM @holland1983sbm], [GM], [Block Poisson MLE], [$O(T(m K + n K^2))$], [No], [No], [*No*], [Exact phase theory; hub split ($alpha >= 9$)],
      [\#131 DC-SBM @karrer2011dcsbm], [GM], [Deg.-corr. block MLE], [$O(T K (m + n))$], [No], [No], [*Yes*], [Hub-safe + non-assort.; exact $Q_(gamma^*)$],
      [\#132 PIC @lin2010pic], [SP], [Truncated $P^t v^((0))$], [$O(T m + I n K)$], [No], [No], [*Yes*], [Ultra-fast sparse; $Theta(K^(-2))$ 1D fail],
    )
  },
)

#subsec("A", "Exact Cross-Algorithm Reductions & Strict Dominance Lattice")

Theorems 1--10 and Proposition 1 establish an interlocking algebraic lattice across Category V:
1. *Exact Reductions:* *Ncut (\#104) $equiv$ NJW (\#103)* via $u = D^(-1/2) v$ and weighted kernel $k$-Means @dhillon2004kernel (Theorem 1.3), coinciding with *Ratio Cut (\#105)* on $d$-regular graphs; *Kluger Biclustering (\#123)* equals bipartite *Ncut* (Theorem 3.1); *Spinglass (\#117)* equals generalized modularity $Q_gamma$ (Proposition 1.1) and symmetric planted *DC-SBM (\#131)* equals $Q_(gamma^*)$ at logarithmic mean $gamma^*$ (Theorem 10.3); *Walktrap (\#114)* equals diffusion-map *Ward clustering* (Theorem 6.2); *SLPA (\#126)* ($T = 1, r > 0.5$) reduces to *LPA (\#110)* $equiv$ `top` *Chinese Whispers (\#111)*; and *Fast Unfolding (\#129)* is the $gamma = 1$ core of *Louvain (\#107)*.
2. *Strict Practical Dominance:* *Leiden (\#108)* strictly dominates *Louvain (\#107)*, *Fast Unfolding (\#129)*, and *CNM (\#115)* (Theorem 5); *NJW (\#103) / Ncut (\#104)* strictly dominate *Ratio Cut (\#105)* on heterogeneous graphs (Theorem 2); *DC-SBM (\#131)* strictly dominates *SBM (\#130)* on skewed-degree networks (Theorem 10.2); and *METIS (\#121) / Chaco (\#122)* strictly dominate standalone *Kernighan--Lin (\#120)* on large sparse graphs (Theorem 8.2).

#subsec("B", "Mathematically Grounded Practitioner's Decision Guide")

Select among Algorithms 103--132 by matching network regime to proved guarantees:
1. *Default Large Undirected Networks (Unknown $K$):* Default to *Leiden (\#108)*—using $Q_gamma$ for degree-normalized exploration and *Leiden-CPM* whenever resolving small dense cliques inside massive graphs is critical (Theorems 4--5). Retire *Louvain (\#107)*, *Fast Unfolding (\#129)*, and *CNM (\#115)* in new production code.
2. *Directed Flow Networks & Citation/Transport Graphs:* Choose *Infomap (\#113)* to exploit directional PageRank asymmetry and exponentially push back the resolution limit (Theorem 6.1), auditing sparse chains for field-of-view over-merging.
3. *Geometric Affinity Graphs, Manifolds & Rectangular Matrices (Fixed $K$):* Use *NJW (\#103)* or *Ncut (\#104)* on sparse symmetric $k$NN affinities with Lanczos solvers (avoiding *Ratio Cut \#105*, Theorem 2). Substitute *PIC (\#132)* for millisecond $O(T m)$ embeddings only when $K <= 6$--$8$ and the eigengap is crisp (avoiding Theorem 3.2's $Theta(K^(-2))$ collisions). For nonnegative bipartite matrices, use *Kluger (\#123)* with Sinkhorn bistochastic scaling.
4. *Balanced $K$-Way HPC & Mesh Partitioning:* Deploy *METIS (\#121)* for fastest irregular sparse $k$-way partitioning and *Chaco (\#122)* when multilevel spectral bisection/quadrisection improves geometric PDE mesh cuts (Theorem 8.2), reserving standalone *Kernighan--Lin (\#120)* for small graphs or post-processing (e.g., *Leading Eigenvector \#116*).
5. *Non-Assortative & Statistically Principled Block Inference:* Deploy *DC-SBM (\#131)* (selecting $K$ via MDL/BIC) for bipartite, disassortative, or core--periphery structure without hub collapse (Theorem 10.2); use *SBM (\#130)* only when homogeneous degrees within blocks are physically known.
6. *Overlapping Communities & Null-Tested Significance:* Choose *BigCLAM (\#125)* for massive dense-overlap graphs (Theorem 9.1), *SLPA (\#126)* for fastest $O(T m)$ unknown-$K$ overlapping covers, *OSLOM (\#124)* when strict Erdős--Rényi noise rejection and outlier identification are mandatory (Theorem 9.3), *DEMON (\#128)* on triangle-rich ego-centric social graphs, and *CPM (\#127)* strictly for exact rolling $k$-clique motifs.
7. *Bio/PPI Graphs, Signed Networks, Non-Metric Exemplars & Small Dendrograms:* Choose *MCL (\#109)* for weighted protein--homology graphs ($r in [1.4, 3.0]$ with self-loops, Theorem 6.3), *Spinglass (\#117)* for small signed graphs ($n <= 10^4$), *Affinity Propagation (\#106)* for signed/asymmetric non-metric exemplar clustering (Theorem 7.3), *HCS (\#118)* strictly for dense diameter-2 cliques, and *Girvan--Newman (\#112)* / *Walktrap (\#114)* for small-to-medium exploratory hierarchies ($n <= 5 times 10^3$).
8. *Billion-Edge Ultra-Fast Regimes:* Deploy degree-damped asynchronous *Chinese Whispers (\#111)* (`lin`/`log`) or multi-run *LPA (\#110)* when $K$ is unknown, and *FluidC (\#119)* when $K$ is known and monster-community immunity is mandatory (Theorem 7.2).

#subsec("C", "Open Theoretical Frontiers & Cross-Paradigm Bridges")

1. *Optimal Sparse Spectral Operators:* Standard Laplacians suffer hub localization near the sparse Kesten--Stigum threshold $(a - b)^2 > 2(a + b)$ (Theorem 10.1); scalable Bethe--Hessian and non-backtracking operators achieving optimal weak recovery under heavy-tailed DC-SBM degrees bridge Sections I and IV.
2. *Unified Connected, Resolution-Free, Overlapping Statistical Detection:* Combining Leiden-CPM's connectedness and sub-partition invariance, OSLOM's order-statistic noise rejection, BigCLAM's dense-overlap geometry, and DC-SBM's non-assortative identifiability inside one $O(m log n)$ algorithm remains open.
3. *Higher-Order Hypergraph & Simplicial Clustering:* Extending Cheeger bounds, Leiden refinement, Infomap coding, and Chernoff--Hellinger thresholds beyond pairwise projections to nonlinear hypergraph Laplacians and simplicial block models avoids clique-expansion distortion.
4. *Neural Spectral Pooling & Differentiable Block Models:* End-to-end Graph Neural Networks optimizing differentiable normalized-cut, modularity, and DC-SBM objectives bridge Category V directly to Part VIII's deep clustering architectures (Comparison Chapters I, II, XII--XIII).
