#import "../template/ieee.typ": *

#cmp_entry(
  roman: "X",
  title: "Intra-Category Comparative Analysis: Neural, Deep, Streaming, and Evolutionary (Algs. 176-200)",
  scope_label: "Head-to-Head Theoretical & Practical Comparison of Algorithms 176-200",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Part VIII develops the final twenty-five clustering algorithms of this monograph (Algs. 176--200), dismantling four foundational restrictions of classical batch clustering: fixed ambient feature coordinates, static finite datasets, fixed topological grids, and myopic single-trajectory local descent @jain2010data @aggarwal2013data @xu2005survey. Spanning competitive topological neural networks, deep representation and self-supervised clustering, single-pass and anytime data-stream engines, and population-based stochastic metaheuristics, all twenty-five methods instantiate a unified *coupled dynamical clustering operator*:
$
  (Theta_(t+1), cal(G)_(t+1), cal(S)_(t+1)) &= cal(T)_t (Theta_t, cal(G)_t, cal(S)_t; cal(D)_t),
$
where $Theta_t$ parameterizes prototypes, neural encoder--decoder weights, or candidate centroid populations; $cal(G)_t$ tracks explicit or latent neighborhood geometry (fixed grids, pure $k$-simplices, Delaunay graphs, resonance hyperboxes, pretext $k$-NN graphs, micro-cluster trees, or density lattices); $cal(S)_t$ maintains fading sufficient statistics, coreset summaries, or stigmergic pheromone and velocity memories; and $cal(D)_t$ arrives as sequential samples, augmented multi-view pairs, mini-batches, or unbounded non-stationary streams.

Four deep theoretical fault lines divide Category VIII:
1. *Topological Fidelity vs. Stability--Plasticity (Algs. 176--180):* Fixed lattices and $k$-simplices (SOM, GCS) vs. rank-order potentials and dynamic Delaunay graphs (Neural Gas, GNG) vs. vigilance-bounded resonance hyperboxes (ART / Fuzzy ART).
2. *Representational Collapse vs. Non-Trivial Self-Supervision (Algs. 181--192):* How unregularized deep clustering (DEC, DeepCluster, unpenalized CC) collapses onto degenerate solutions, and how Lipschitz reconstruction (IDEC, DCN), generative ELBOs (VaDE), balanced denoising targets (DEPICT), hyperspherical/one-hot geometry (JULE, DAC), mutual information (IIC, CC), Fano-calibrated neighbors (SCAN), and entropic Sinkhorn equipartition (SwAV) mathematically prevent collapse.
3. *Temporal Forgetting vs. Provable Coreset Approximation (Algs. 193--197):* Pyramidal snapshots (CluStream) and exponential decay semigroups (DenStream, ClusTree, D-Stream) with stream-length-independent memory vs. Merge-and-Reduce $D^2$-coresets (StreamKM++) guaranteeing end-to-end $O(log k)$ approximation yet suffering linear concept-drift inertia.
4. *Asymptotic Markov Optimality vs. Finite-Time Premature Lock-In (Algs. 198--200):* Discrete elitist Markov chains (GKA, bounded ACO) achieving almost-sure NP-hard global SSE convergence vs. continuous second-order swarm dynamics (PSO) admitting non-stationary collapse, alongside symmetry washout and exponential basin-escape barriers.

Sections I--IV establish eight head-to-head unified theorems with full constructive proofs across these four subfamilies; Section V integrates all twenty-five algorithms in a master comparison matrix, exact reduction lattice, practitioner's decision guide, and open-problem roadmap.

// =============================================================================
// SECTION I
// =============================================================================
#sec("I", "Classical Competitive & Topological Neural Networks")

Five biologically inspired neural architectures learn quantization prototypes $W = {w_1, dots, w_M} subset.eq bb(R)^d$ via competitive adaptation: Self-Organizing Maps (SOM, \#176) @kohonen1982som, Growing Neural Gas (GNG, \#177) @fritzke1995gng, Adaptive Resonance Theory (ART / Fuzzy ART, \#178) @carpenter1987art @carpenter1991fuzzyart, Neural Gas (NG, \#179) @martinetz1993neuralgas, and Growing Cell Structures (GCS, \#180) @fritzke1994gcs. While SOM, NG, GCS, and GNG share neighborhood-coupled prototype attraction, they split sharply over potential-function existence and manifold topology; ART breaks from competitive drift altogether via vigilance-gated hyperbox resonance.

#subsec("A", "Potential Existence, Topological Folding & Delaunay Fidelity")

Given input $x tilde.op p(x)$ on $bb(R)^d$, competitive vector quantizers update prototypes $w_j$ ($j = 1, dots, M$) via
$
  Delta w_j = epsilon_t h_t (j, x; W) (x - w_j),
$
differing fundamentally in the coupling kernel $h_t$ and graph evolution $G_t = (V_t, E_t)$:
- *SOM (\#176)* @kohonen1982som fixes a discrete lattice (e.g., 1D chain or 2D grid) with metric $d_"grid" (j, i)$, assigns the Best Matching Unit (BMU) $i^*(x) = arg min_i norm(x - w_i)_2^2$, and sets $h_t (j, x; W) = h_(sigma_t) (d_"grid" (j, i^*(x)))$ (typically Gaussian $exp(-d_"grid"^2 / (2 sigma_t^2))$). Heskes' batch variant replaces $i^*(x)$ by the neighborhood-averaged winner $i_"H"^*(x) = arg min_i sum_j h_(i j) norm(x - w_j)_2^2$.
- *Neural Gas (\#179)* @martinetz1993neuralgas discards any fixed lattice: it ranks all prototypes by Euclidean distance to $x$, letting $k_j (x, W) in {0, dots, M - 1}$ count prototypes strictly closer to $x$ than $w_j$, and sets $h_lambda (k_j (x, W)) = exp(-(k_j (x, W)) / lambda)$.
- *Growing Cell Structures (GCS, \#180)* @fritzke1994gcs enforces a pure $k$-dimensional simplicial complex $cal(K)$ (every edge belongs to a $k$-simplex), updates the BMU $s_1$ ($epsilon_b$) and its direct graph neighbors ($epsilon_n$), accumulates local quantization error $tau_(s_1) <- tau_(s_1) + norm(x - w_(s_1))_2^2$, and every $lambda_"ins"$ steps splits the longest edge incident to the maximum-error vertex while preserving pure $k$-simpliciality.
- *Growing Neural Gas (GNG, \#177)* @fritzke1995gng drops fixed dimensionality: upon each $x$, it finds the 1st and 2nd BMUs $(s_1, s_2)$, inserts or resets edge $(s_1, s_2)$ to age $0$ (Competitive Hebbian Learning), increments ages of all other edges incident to $s_1$, prunes edges exceeding $a_"max"$ (removing isolated nodes), and periodically bisects the edge $(q, f)$ between the highest-error vertex $q$ and its highest-error neighbor $f$.

Theorem 1 establishes a sharp four-way theoretical separation across integrability, topological ordering, simplicial rigidity, and induced Delaunay recovery.

#theorem(title: "Potential Existence, Topological Folding & Induced Delaunay Fidelity", origin: "Original Unified Proof")[
  Let $p(x)$ be a continuous density supported on compact $cal(M) subset.eq bb(R)^d$ with $M >= 2$ distinct prototypes $W = (w_1, dots, w_M)$.
  1. *(Non-Existence of Online SOM \#176 Potential vs. Exact Gradient of Neural Gas \#179 & Batch MM Descent):* For Kohonen's online SOM with hard Euclidean BMU $i^*(x) = arg min_i norm(x - w_i)_2^2$ and non-diagonal symmetric kernel $h_(j i) > 0$ ($j != i$), the expected drift field
     $
       F_j^"SOM" (W) := integral h_(j, i^*(x)) (x - w_j) p(x) dif x
     $
     is *not* conservative: across every shared Voronoi facet $V_(a b)$ with $p > 0$, mixed partial derivatives violate Clairaut symmetry ($partial_(w_b) F_a^"SOM" != (partial_(w_a) F_b^"SOM")^top$), so *no* $C^2$ Lyapunov potential $E(W)$ has $-nabla E = F^"SOM"$. Conversely, online Neural Gas (\#179) is the *exact negative gradient* of the continuous rank-distortion energy
     $
       E_"NG" (W; lambda) := 1 / (2 C(lambda)) sum_(j=1)^M integral h_lambda (k_j (x, W)) norm(x - w_j)_2^2 p(x) dif x
     $
     (with $C(lambda) = sum_(r=0)^(M-1) h_lambda (r)$) because every bisector jump cancels identically; likewise, Heskes Batch SOM and Batch NG strictly decrease their Majorize--Minimize (MM) Lyapunov functionals at every iteration.
  2. *(1D SOM Ordered Invariance vs. Fixed-Dimension Folding of SOM \#176 & GCS \#180):* In 1D SOM ($cal(M) = [0, 1]$) with strictly decreasing $h(|i - j|)$ and $epsilon in (0, 1]$, the strictly ordered cone $cal(C)_+ = {w_1 < w_2 < dots < w_M}$ is positively invariant. However, whenever $cal(M)$ has either intrinsic dimension $d_"int" > k$ or $c >= 2$ disconnected components separated by gap $delta > 0$:
     - Every connected $k$-dimensional SOM lattice ($sigma > 0$) places interpolating *dead prototypes* in zero-density inter-component voids or exhibits topological foldovers ($d_"int" > k$).
     - Every pure $k$-simplicial GCS (\#180) maintains an invariant Euler characteristic $chi(cal(K)) = 1$ under interior edge bisection (until explicit pruning), locking local neighborhood dimension strictly to $k$ and preventing topological adaptation across mixed-dimension manifolds $cal(M)_1 union cal(M)_2$ with $dim cal(M)_1 != dim cal(M)_2$.
  3. *(Exact Induced Delaunay Recovery & Exponential Edge-Age Phase Separation of GNG \#177):* Fix prototypes $W$ with non-degenerate Voronoi diagram ${V_j}$ on $cal(M)$. An edge $(u, v)$ is generated by GNG's Competitive Hebbian rule with positive probability *if and only if* $(u, v)$ lies in the *induced Delaunay triangulation* $cal(D)_cal(M) (W) := {(u, v) : V_u inter V_v inter cal(M) != emptyset}$. Moreover, if a spurious cross-void edge has co-activation probability $p_(u v) = 0$ while a true Delaunay edge satisfies $p_(u v) := bb(P)[{s_1, s_2} = {u, v} | s_1 in {u, v}] > 0$, then:
     - Every spurious edge is deterministically pruned within at most $a_"max"$ visits to ${u, v}$;
     - A true edge survives any given age-reset cycle with exact probability $1 - (1 - p_(u v))^(a_"max")$, yielding exponential phase separation as $a_"max" -> infinity$.
]

#proof[
  *Part 1:* Partition $bb(R)^d$ into Voronoi cells $V_i (W) = {x : norm(x - w_i)_2 <= norm(x - w_ell)_2 forall ell}$ with shared facet $S_(a b) = V_a inter V_b$ and outward unit normal $n_(a b) = (w_b - w_a) / norm(w_b - w_a)_2$. By Reynolds' transport theorem, differentiating $F_a^"SOM" (W) = sum_i h_(a i) integral_(V_i) (x - w_a) p(x) dif x$ with respect to $w_b$ ($b != a$) picks up only the boundary flux across facets touching $V_b$. On $S_(a b)$, the normal velocity of the bisector $(x - (w_a + w_b) / 2)^top (w_b - w_a) = 0$ with respect to $w_b$ is $v_n (x) = (x - w_b)^top / norm(w_b - w_a)_2$ (entering $V_a$ as $w_b$ moves away). Thus the cross-Jacobian contribution along $S_(a b)$ is
  $
    partial_(w_b) F_a^"SOM" |_(S_(a b)) &= (h_(a a) - h_(a b)) / norm(w_b - w_a)_2 \
    &quad dot integral_(S_(a b)) (x - w_a)(x - w_b)^top p(x) dif cal(H)^(d-1) (x).
  $
  Swapping $a$ and $b$ and transposing yields the identical surface integral multiplied by $(h_(b b) - h_(b a)) / norm(w_b - w_a)_2$. However, $F_a^"SOM"$ also receives flux across every *third* facet $S_(b c)$ ($c != a, b$) weighted by $(h_(a c) - h_(a b))$, whereas $(partial_(w_a) F_b^"SOM")^top$ integrates over $S_(a c)$ weighted by $(h_(b c) - h_(b a))$. Choosing $M = 3$ prototypes on a line ($w_1 < w_2 < w_3$) with $S_(1 3) = emptyset$ and $p > 0$ on $S_(1 2), S_(2 3)$ gives nonzero $S_(2 3)$ contribution to $partial_(w_3) F_1^"SOM"$ (proportional to $h_(1 2) - h_(1 3) != 0$) while $S_(1 2)$ gives a distinct point-valued contribution to $partial_(w_1) F_3^"SOM"$ evaluated at $(w_1 + w_2)/2 != (w_2 + w_3)/2$; hence $partial_(w_3) F_1^"SOM" != (partial_(w_1) F_3^"SOM")^top$, refuting existence of any $C^2$ potential @kohonen1982som.

  In Neural Gas (\#179), $k_j (x, W)$ is piecewise constant away from bisectors $H_(a b) = {x : norm(x - w_a)_2 = norm(x - w_b)_2}$. Differentiating $E_"NG"$ with respect to $w_a$ produces the interior term $-C(lambda)^(-1) integral h_lambda (k_a (x, W)) (x - w_a) p(x) dif x$ plus a boundary integral over each bisector $H_(a b)$ where ranks $k_a$ and $k_b$ swap values ($r <-> r + 1$). Across $H_(a b)$, the integrand jump is
  $
    Delta_((a b)) (x) &= [h_lambda (r) norm(x - w_a)_2^2 + h_lambda (r + 1) norm(x - w_b)_2^2] \
    &quad - [h_lambda (r + 1) norm(x - w_a)_2^2 + h_lambda (r) norm(x - w_b)_2^2] \
    &= (h_lambda (r) - h_lambda (r + 1)) (norm(x - w_a)_2^2 - norm(x - w_b)_2^2) equiv 0
  $
  identically on $H_(a b)$! Every boundary term vanishes identically, proving $-nabla_(w_a) E_"NG" = C(lambda)^(-1) bb(E)[h_lambda (k_a (x, W)) (x - w_a)]$ @martinetz1993neuralgas. For Heskes Batch SOM, $E_"H" (W, c) = sum_i sum_j h_(c_i, j) norm(x_i - w_j)_2^2$ decreases strictly under alternating minimization of assignments $c_i <- arg min_c sum_j h_(c, j) norm(x_i - w_j)_2^2$ and weighted centroids $w_j <- (sum_i h_(c_i, j) x_i) / (sum_i h_(c_i, j))$ (identical MM descent holds for Batch NG).

  *Part 2:* In 1D SOM with $w_1 < dots < w_M$ and BMU $i^* = i^*(x)$, for any adjacent pair $j < j + 1$ the updated gap is
  $
    w'_(j+1) - w'_j &= (1 - epsilon h_(j+1, i^*)) w_(j+1) - (1 - epsilon h_(j, i^*)) w_j \
    &quad + epsilon (h_(j+1, i^*) - h_(j, i^*)) x.
  $
  If $i^* <= j$, then $x <= (w_j + w_(j+1)) / 2 < w_(j+1)$ and $h_(j, i^*) > h_(j+1, i^*)$, so rewriting the right side as $(1 - epsilon h_(j, i^*))(w_(j+1) - w_j) + epsilon (h_(j, i^*) - h_(j+1, i^*))(w_(j+1) - x) > 0$ proves $w'_(j+1) > w'_j$; the symmetric rearrangement holds for $i^* >= j + 1$. Thus $cal(C)_+$ is positively invariant. Conversely, if $cal(M) = cal(M)_1 union cal(M)_2$ has two convex components separated by distance $delta > 0$ and a connected SOM lattice spans both, adjacent lattice nodes $(u, v)$ straddle the void; for $sigma > 0$, a node $j$ coupled symmetrically to both components equilibrates at a convex combination strictly inside the zero-density gap. For GCS (\#180), inserting a vertex $r$ at the midpoint of edge $(q, f)$ subdivides every $k$-simplex $sigma supset {q, f}$ into two $k$-simplices sharing a new $(k - 1)$-face; by Euler's formula, adding $1$ vertex, $m_1$ edges, ..., $m_k$ $k$-simplices across the link of $(q, f)$ preserves $chi(cal(K)) = sum_(j=0)^k (-1)^j f_j$ identically (since the subdivided star is a PL $k$-ball whose Euler characteristic equals $1$ before and after stellar subdivision) and keeps every simplex dimension fixed at $k$ @fritzke1994gcs.

  *Part 3:* By definition, ${s_1 (x), s_2 (x)} = {u, v}$ holds iff $norm(x - w_u)_2, norm(x - w_v)_2 <= norm(x - w_ell)_2$ for all $ell != u, v$, i.e., $x$ lies in the second-order Voronoi cell $V_(u v)$ whose interior bisector is $V_u inter V_v$. Because $p > 0$ continuously on $cal(M)$, $bb(P)[{s_1, s_2} = {u, v}] > 0 <=> V_u inter V_v inter cal(M) != emptyset <=> (u, v) in cal(D)_cal(M) (W)$ @fritzke1995gng. Whenever $s_1 in {u, v}$ occurs (a visit to endpoints ${u, v}$), edge $(u, v)$ resets to age $0$ with probability $p_(u v)$ and increments by $+1$ with probability $1 - p_(u v)$. If $p_(u v) = 0$, the edge ages by $+1$ on every visit and is deterministically deleted upon reaching $a_"max"$. If $p_(u v) > 0$, deletion from age $0$ requires $a_"max"$ consecutive non-co-activations, occurring with probability $(1 - p_(u v))^(a_"max")$, so survival probability across a cycle is exactly $1 - (1 - p_(u v))^(a_"max")$.
]

#subsec("B", "Stability--Plasticity Dilemma & Vigilance Resonance (ART vs. VQ)")

Online competitive vector quantizers (SOM, NG, GCS, GNG) face Grossberg's *stability--plasticity dilemma*: keeping learning rate $epsilon > 0$ enables continuous adaptation to novel clusters yet erases previously learned memories under task shift; annealing $epsilon_t -> 0$ freezes memory yet kills plasticity. Adaptive Resonance Theory (ART1 / Fuzzy ART, \#178) @carpenter1987art @carpenter1991fuzzyart resolves this dilemma geometrically via *complement-coded hyperboxes* and *vigilance gating*. Each input $a in [0, 1]^d$ is embedded as $I = (a, bold(1) - a) in [0, 1]^(2 d)$ ($norm(I)_1 = d$). Category $j$ maintains weight $w_j = (u_j, bold(1) - v_j) in [0, 1]^(2 d)$ encoding hyperbox $R_j = [u_j, v_j] subset.eq [0, 1]^d$. With fuzzy componentwise minimum $(x and y)_r = min(x_r, y_r)$, choice scores $T_j (I) = norm(I and w_j)_1 / (alpha + norm(w_j)_1)$ are searched in descending order until the *vigilance test* passes:
$
  norm(I and w_j)_1 / norm(I)_1 >= rho in [0, 1],
$
updating $w_j^"new" = beta (I and w_j) + (1 - beta) w_j$ ($beta in (0, 1]$), or instantiating a new category $w_(K+1) = I$ if all existing categories fail vigilance.

#theorem(title: "Catastrophic VQ Drift vs. Vigilance-Bounded Monotone Resonance & High-d Concentration", origin: "Original Unified Proof")[
  1. *(Catastrophic Forgetting of Constant-Rate Competitive VQ vs. Plasticity Lock):* Consider any prototype $w_j$ of SOM (\#176), Neural Gas (\#179), GNG (\#177), or GCS (\#180) converged at old-task cluster mean $mu_1 in bb(R)^d$, followed by $T$ i.i.d. samples from a shifted task distribution $cal(D)_2$ with mean $mu_2$ at distance $Delta = norm(mu_2 - mu_1)_2 > 0$. If $w_j$ receives expected learning-rate coupling $epsilon_t c_0$ ($c_0 > 0$) on $cal(D)_2$:
     - With constant rate $epsilon_t equiv epsilon in (0, 1 / c_0)$, $w_j$ drifts exponentially toward $mu_2$:
       $
         norm(bb(E)[w_j^((T))] - mu_1)_2 = (1 - (1 - epsilon c_0)^T) Delta -> Delta,
       $
       destroying old-task quantization ($"SSE"_1$ inflates by $+n_1 Delta^2$).
     - With summable annealing $sum_(t=1)^infinity epsilon_t = S < infinity$ after $T_0$ steps, $norm(w_j^((T)) - w_j^((T_0)))_2 <= R_cal(X) sum_(t > T_0) epsilon_t -> 0$, making novel clusters unreachable.
  2. *(Exact $L_1$-Perimeter Identity, Universal Stability Bound & Complement-Coding Necessity of Fuzzy ART \#178):* Under complement coding ($norm(I)_1 = d$):
     - Every category weight $w_j = (u_j, bold(1) - v_j)$ satisfies the exact $L_1$-perimeter identity
       $
         norm(w_j)_1 = d - norm(v_j - u_j)_1 = d - |R_j|_1,
       $
       where $|R_j|_1 := sum_(r=1)^d (v_(j r) - u_(j r))$ is the Manhattan hyperbox size.
     - Resonance $norm(I and w_j)_1 / d >= rho$ holds *if and only if* the expanded box satisfies $|R_j^"new"|_1 <= (1 - rho) d$. Hence *no* novel task at $L_1$-distance $> (1 - rho) d$ from $R_j$ can *ever* alter $w_j$ (zero catastrophic drift while retaining full plasticity via new category creation!).
     - Under fast learning ($beta = 1$) on any finite dataset $cal(X) subset.eq [0, 1]^d$, weights are monotonically non-increasing and converge in finitely many passes (strictly bounded by $2 d$ distinct coordinate extrema per category). Conversely, *without* complement coding, $w_j^"new" = I and w_j <= w_j$ erodes $norm(w_j)_1 -> 0$ on noisy/varied data, causing unbounded category proliferation.
  3. *(High-Dimensional Vigilance Concentration of Fuzzy ART \#178):* On i.i.d. uniform vectors $a_1, a_2 tilde.op "Unif"([0, 1]^d)$ with complement coding $I_1, I_2$,
     $
       bb(E)[norm(I_1 and I_2)_1 / d] = 2 / 3, wide
       "Var"(norm(I_1 and I_2)_1 / d) = 1 / (18 d),
     $
     and by Hoeffding's inequality, for every $eta > 0$,
     $
       "Pr"(abs(norm(I_1 and I_2)_1 / d - 2 / 3) >= eta) <= 2 exp(-2 d eta^2).
     $
     Consequently, as $d -> infinity$, any fixed $rho > 2 / 3$ forces $Omega(n)$ singleton proliferation, while any $rho < 2 / 3$ collapses unrelated points into giant hyperboxes!
]

#proof[
  *Part 1:* Taking expectations in $w_j^((t+1)) = w_j^((t)) + epsilon_t h_t (x_t - w_j^((t)))$ with $bb(E)[h_t] = c_0$ and $bb(E)[x_t] = mu_2$ yields the linear contraction $bb(E)[w_j^((t+1)) - mu_2] = (1 - epsilon c_0) bb(E)[w_j^((t)) - mu_2]$. Iterating from $w_j^((0)) = mu_1$ gives $bb(E)[w_j^((T)) - mu_1] = (1 - (1 - epsilon c_0)^T)(mu_2 - mu_1)$, sending drift to $Delta$ and old-cluster squared error on $n_1$ points centered at $mu_1$ up by $n_1 norm(bb(E)[w_j^((T))] - mu_1)_2^2 -> n_1 Delta^2$. Under summable $epsilon_t$, triangle inequality yields $norm(w_j^((T)) - w_j^((T_0)))_2 <= sum_(t=T_0+1)^T epsilon_t norm(x_t - w_j^((t-1)))_2 <= R_cal(X) sum_(t > T_0) epsilon_t -> 0$.

  *Part 2:* Since $0 <= u_(j r) <= v_(j r) <= 1$,
  $
    norm(w_j)_1 &= sum_(r=1)^d (u_(j r) + 1 - v_(j r)) = d - sum_(r=1)^d (v_(j r) - u_(j r)) \
    &= d - |R_j|_1.
  $
  For input $I = (a, bold(1) - a)$, the fuzzy intersection $I and w_j$ has coordinates $min(a_r, u_(j r))$ and $min(1 - a_r, 1 - v_(j r)) = 1 - max(a_r, v_(j r))$, encoding the minimal enclosing hyperbox $R_j^"new" = [u_j and a, v_j or a]$. Applying the perimeter identity to $I and w_j$ gives
  $
    norm(I and w_j)_1 / d = 1 - |R_j^"new"|_1 / d >= rho <=> |R_j^"new"|_1 <= (1 - rho) d.
  $
  Thus any sample $a$ whose inclusion would stretch $|R_j^"new"|_1 > (1 - rho) d$ fails vigilance, leaving $w_j$ untouched and spawning a fresh category. Under $beta = 1$, each coordinate of $u_j$ only decreases to elements of ${a_(i r) : i <= n}$ and each coordinate of $v_j$ only increases to elements of ${a_(i r) : i <= n}$ (finitely many values); because every weight change strictly reduces $norm(w_j)_1$ within a finite lattice, all weights stabilize in finitely many epochs @carpenter1991fuzzyart. Without complement coding, $I = a$ lacks the upper-bound half $bold(1) - a$, so $w_j^"new" = a and w_j$ tracks only coordinatewise minima; presenting points near opposite faces of $[0, 1]^d$ drives $w_j -> bold(0)$, causing $norm(a and w_j)_1 / norm(a)_1 -> 0 < rho$ and triggering unbounded category creation.

  *Part 3:* For independent $U, V tilde.op "Unif"([0, 1])$, coordinate $r$ of $I_1 and I_2$ contributes $Z_r := min(U, V) + min(1 - U, 1 - V) = 1 - |U - V| in [0, 1]$. Since $bb(E)[|U - V|] = 1 / 3$ and $bb(E)[(U - V)^2] = 1 / 6$, we obtain $bb(E)[Z_r] = 2 / 3$ and $"Var"(Z_r) = 1 / 6 - 1 / 9 = 1 / 18$. Averaging $d$ i.i.d. $[0, 1]$-valued variables yields mean $2 / 3$, variance $1 / (18 d)$, and Hoeffding's tail bound $"Pr"(|overline(Z) - 2 / 3| >= eta) <= 2 exp(-2 d eta^2)$ @beyer1999nearest.
]

#remark(title: "Order Sensitivity & GCS vs. GNG Trade-off")[
  Theorem 2 explains both ART's hard stability and its caveats: box growth is order-dependent and $rho$ requires $1 - O(1 / sqrt(d))$ tuning in high $d$. Meanwhile, Theorem 1 shows GNG (\#177) subsumes SOM (\#176) and GCS (\#180) on arbitrary manifolds, while SOM excels at 2D visual projection and Neural Gas (\#179) yields the sharpest smooth VQ potential.
]

#wide_figure(
  caption: [*Table I:* Head-to-Head Theoretical & Algorithmic Comparison of Classical Competitive & Topological Neural Networks (Section I: Algs. 176--180).],
  {
    set text(size: 7.4pt)
    set table(inset: (x: 3.4pt, y: 2.5pt))
    table(
      columns: (1.06fr, 0.96fr, 0.88fr, 0.88fr, 0.94fr, 1.28fr),
      align: (left, left, center, center, left, left),
      [Algorithm], [Topology / Structure], [Exact Potential], [Adaptive $M, d_"int"$], [Time Complexity], [Distinctive Strength / Failure Mode],
      [SOM (\#176) @kohonen1982som], [Fixed $k$-D lattice], [No online (Batch MM)], [No (fixed grid)], [$O(E n M d)$], [2D topo. visualization; folding & void dead units],
      [GNG (\#177) @fritzke1995gng], [Dynamic Delaunay graph], [Local error heuristic], [*Yes (both)*], [$O(T M_"max" d)$], [Exact induced Delaunay; age-$a_"max"$ sensitivity],
      [ART (\#178) @carpenter1991fuzzyart], [Complement hyperboxes], [Monotone $L_1$ lattice], [Yes ($K$ via $rho$)], [$O(T K d)$], [*Zero catastrophic drift*; order & $rho -> 2/3$ conc.],
      [Neural Gas (\#179) @martinetz1993neuralgas], [Rank-order $k_j (x, W)$], [*Exact $C^1$ $E_"NG"$*], [Fixed $M$, any $d_"int"$], [$O(T M (d + log M))$], [Exact gradient + no dead units; no explicit graph],
      [GCS (\#180) @fritzke1994gcs], [Pure $k$-simplicial], [Local error heuristic], [Adaptive $M$, fixed $k$], [$O(T (M_"max" d + k^2 / lambda))$], [Euler-invariant mesh; fails mixed-dim. manifolds],
    )
  },
)

// =============================================================================
// SECTION II
// =============================================================================
#sec("II", "Deep Representation & Self-Supervised Clustering")

Twelve algorithms in Category VIII (Algs. 181--192) jointly learn a nonlinear neural feature map $z_i = f_theta (x_i) in bb(R)^p$ and cluster assignments, eliminating the fixed-metric bottleneck of classical partitioning. Yet joint optimization over $(theta, C)$ introduces an existential pathology absent from fixed-feature clustering: *representational collapse*, wherein the network minimizes the clustering objective trivially by mapping all inputs to a constant vector, a rank-1 subspace, or arbitrary random noise. We classify all twelve methods into three architectural families and prove exact collapse and non-collapse criteria for every single one.

#subsec("A", "Feature-Collapse Modes vs. Autoencoder & Variational Safeguards (Algs. 181--185)")

Five methods build upon neural autoencoders $x_i |-> z_i = f_theta (x_i) |-> hat(x)_i = g_psi (z_i)$:
- *DEC (\#181)* @xie2016dec pretrains an autoencoder, *discards the decoder*, and fine-tunes $(theta, {mu_j})$ by minimizing $"KL"(P || Q)$ between Student-$t$ soft assignments $q_(i j) = (1 + norm(z_i - mu_j)_2^2)^(-1) / (sum_ell (1 + norm(z_i - mu_ell)_2^2)^(-1))$ and sharpened targets $p_(i j) = (q_(i j)^2 / f_j) / (sum_ell q_(i ell)^2 / f_ell)$ ($f_j = sum_i q_(i j)$).
- *IDEC (\#182)* @guo2017idec retains the decoder throughout clustering:
  $
    cal(L)_"IDEC" (theta, psi, M) = 1 / N sum_(i=1)^N norm(x_i - g_psi (f_theta (x_i)))_2^2 + gamma "KL"(P || Q).
  $
- *DCN (\#183)* @yang2017dcn couples reconstruction directly to hard $k$-Means ($s_i in {e_1, dots, e_K}$):
  $
    cal(L)_"DCN" (theta, psi, M, S) &= 1 / N sum_(i=1)^N norm(x_i - g_psi (f_theta (x_i)))_2^2 \
    &quad + lambda / (2 N) sum_(i=1)^N norm(f_theta (x_i) - M s_i)_2^2.
  $
- *VaDE (\#184)* @jiang2017vade embeds a Gaussian Mixture prior $p(z, c) = pi_c cal(N)(z; mu_c, "diag"(sigma_c^2))$ inside a Variational Autoencoder @kingma2014vae @blei2006variational, maximizing the Evidence Lower Bound over mean-field $q(z, c | x) = q_phi (z | x) q(c | x)$:
  $
    "ELBO"(x) = bb(E)_(q_phi (z | x)) [log p_psi (x | z)] - "KL"(q(z, c | x) || p(z, c)).
  $
- *DEPICT (\#185)* @ghasedi2017depict trains a deep denoising autoencoder with layerwise noise-robust reconstruction plus a *balanced* KL target solving $min_P "KL"(P || Q) + "KL"(overline(p) || u)$ with uniform marginal prior $u_j = 1 / K$.

Theorem 3 rigorously contrasts DEC's zero-loss collapse against Lipschitz preservation (IDEC/DCN), bottleneck scale-bypass, VaDE's generative ELBO geometry, and DEPICT's balancing distortion.

#theorem(title: "Zero-Loss DEC Collapse, Lipschitz Non-Collapse, Scale-Bypass & VaDE/DEPICT Variational Bounds", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_N} subset.eq bb(R)^d$ with $K >= 2$.
  1. *(Zero-Loss Representational Collapse of Unregularized DEC \#181):* For *any* arbitrary non-empty $K$-partition $Pi = (C_1, dots, C_K)$ of $cal(X)$ (even pure random noise) and any expressive encoder family $f_theta$, DEC achieves $inf_(theta, M) "KL"(P || Q) = 0$ by mapping every $x_i in C_k$ to $z_i -> mu_k$ while separating $min_(a != b) norm(mu_a - mu_b)_2 -> infinity$. Thus DEC has zero structural barrier against distorting the latent manifold to memorize arbitrary initial noise!
  2. *(Lipschitz Non-Collapse & Bottleneck Scale-Bypass of IDEC \#182 and DCN \#183):* Suppose decoder $g_psi$ is $L_g$-Lipschitz and achieves pointwise residuals $r_i := norm(x_i - g_psi (z_i))_2$ with mean squared error $epsilon_"rec"^2 := N^(-1) sum_i r_i^2$. Then:
     - *(Pairwise Separation Guarantee):* Every pair $(i, i')$ obeys
       $
         norm(z_i - z_(i'))_2 >= 1 / L_g (norm(x_i - x_(i'))_2 - r_i - r_(i')),
       $
       and latent variance obeys $"Var"(Z) >= L_g^(-2) (sqrt("Var"(X)) - epsilon_"rec")_+^2 > 0$ whenever $epsilon_"rec" < sqrt("Var"(X))$, strictly precluding constant collapse!
     - *(Bottleneck Scale-Bypass Counterexample):* Conversely, if $g_psi (z) = tilde(g)(W_1 z + b_1)$ has an *unregularized* linear first layer $W_1$, then rescaling $(f_theta, M, W_1) |-> (alpha f_theta, alpha M, alpha^(-1) W_1)$ as $alpha -> 0^+$ leaves $cal(L)_"rec"$ identically constant while driving DCN's $k$-Means penalty $(lambda / (2 N)) sum_i norm(f_theta (x_i) - M s_i)_2^2 = O(alpha^2) -> 0$, neutralizing clustering unless $norm(W_1)_2$ is bounded (via weight decay or spectral normalization).
  3. *(Exact VaDE \#184 Optimal Posterior, Mean-Field Gap & Posterior Collapse vs. DEPICT \#185 Balanced Distortion):*
     - In VaDE (\#184), for fixed $q_phi (z | x)$, the unique ELBO-maximizing cluster responsibility is
       $
         q^*(c | x) = (pi_c exp(-D_c (x))) / (sum_(c'=1)^K pi_(c') exp(-D_(c') (x))),
       $
       where $D_c (x) := bb(E)_(q_phi (z | x)) [log ((q_phi (z | x)) / (p(z | c)))]$, with exact non-negative mean-field gap
       $
         &"ELBO"_"struct" - "ELBO"_"MF" \
         &quad = bb(E)_(q(c | x)) ["KL"(q_phi (z | x) || p(z | x, c))] >= 0.
       $
       Moreover, if decoder $p_psi (x | z)$ attains log-likelihood within $Delta_"dec"$ of the marginal $log p(x)$ independently of $z$, then $q_phi (z | x) = p(z)$ (trivial *posterior collapse*, $I_q (x; z) = 0$) beats any informative encoder whose KL exceeds $Delta_"dec"$.
     - In DEPICT (\#185), minimizing $"KL"(P || Q) + "KL"(overline(p) || u)$ with fixed $Q$ yields the exact closed-form square-root balancing update
       $
         p_(i j)^* = (q_(i j) / sqrt(f_j)) / (sum_(ell=1)^K q_(i ell) / sqrt(f_ell)),
       $
       preventing empty-cluster collapse unconditionally—yet on ground-truth imbalanced priors $pi != u$, forcing $overline(p) -> u$ incurs an irreducible total-variation misassignment floor $>= 1 / 2 norm(pi - u)_1$ whenever $K$-way hard labels track $P$!
]

#proof[
  *Part 1:* Given arbitrary $Pi = (C_1, dots, C_K)$, place prototypes $mu_k = R e_k$ ($R > 0$) and map $f_theta (x_i) = mu_(c(i))$ for $x_i in C_(c(i))$. Then $norm(z_i - mu_(c(i)))_2 = 0$ while $norm(z_i - mu_ell)_2^2 = 2 R^2$ for all $ell != c(i)$, so
  $
    q_(i, c(i)) = (1 + (K - 1) / (1 + 2 R^2))^(-1) -> 1, wide
    q_(i, ell) -> 0 quad (ell != c(i))
  $
  as $R -> infinity$. Since $f_j = sum_i q_(i j) -> |C_j| >= 1$, $p_(i, c(i)) -> 1$ and $p_(i, ell) -> 0$ simultaneously, yielding $"KL"(P || Q) = sum_(i, j) p_(i j) log(p_(i j) / q_(i j)) -> 0$ identically @xie2016dec.

  *Part 2:* By triangle inequality and $L_g$-Lipschitz continuity of $g_psi$,
  $
    norm(x_i - x_(i'))_2 &<= norm(x_i - g_psi (z_i))_2 + norm(g_psi (z_i) - g_psi (z_(i')))_2 \
    &quad + norm(g_psi (z_(i')) - x_(i'))_2 \
    &<= r_i + L_g norm(z_i - z_(i'))_2 + r_(i'),
  $
  rearranging immediately to the pairwise lower bound @guo2017idec. Averaging $norm(x_i - overline(x))_2 <= r_i + norm(g_psi (z_i) - g_psi (overline(z)))_2 + norm(g_psi (overline(z)) - overline(x))_2$ via Minkowski's inequality gives $sqrt("Var"(X)) <= epsilon_"rec" + L_g sqrt("Var"(Z))$, proving $"Var"(Z) >= L_g^(-2) (sqrt("Var"(X)) - epsilon_"rec")_+^2 > 0$. For the scale-bypass counterexample in DCN (\#183), replacing $(f_theta, M, W_1)$ by $(alpha f_theta, alpha M, alpha^(-1) W_1)$ preserves $W_1 (alpha f_theta (x_i)) + b_1 = W_1 f_theta (x_i) + b_1$ identically (so $cal(L)_"rec"$ is invariant), while scaling the $k$-Means term by $alpha^2 -> 0^+$ @yang2017dcn; constraining $norm(W_1)_2 <= L_1$ eliminates this gauge orbit.

  *Part 3:* Expanding VaDE's KL term for $q(z, c | x) = q_phi (z | x) q(c | x)$ gives
  $
    "KL"(q(z, c | x) || p(z, c)) &= sum_(c=1)^K q(c | x) log (q(c | x) / pi_c) \
    &quad + sum_(c=1)^K q(c | x) D_c (x).
  $
  Minimizing this strictly convex functional over $q(dot | x) in Delta_(K - 1)$ via Lagrange multiplier yields $log q^*(c | x) = log pi_c - D_c (x) - log Z(x)$, proving the stated formula @jiang2017vade. Comparing with the exact posterior $p(z, c | x) = p(c | x) p(z | x, c)$ gives
  $
    "KL"(q(z, c | x) || p(z, c | x)) &= "KL"(q(c | x) || p(c | x)) \
    &quad + bb(E)_(q(c | x)) ["KL"(q_phi (z | x) || p(z | x, c))],
  $
  identifying the exact nonnegative mean-field penalty. If $log p_psi (x | z) <= log p(x)$ while an uninformative decoder attains $log p(x) - Delta_"dec"$ at $q_phi (z | x) = p(z)$ (zero KL), any latent code with $"KL"(q_phi (z | x) || p(z)) > Delta_"dec"$ strictly loses in ELBO, proving posterior collapse. For DEPICT (\#185), consider the Lagrangian of $sum_(i, j) p_(i j) log(p_(i j) / q_(i j)) + sum_j overline(p)_j log(K overline(p)_j)$ (under the first-order surrogate fixing $overline(p)_j approx f_j / N$ in the denominator balance factor, or exact block-alternating target @ghasedi2017depict): setting $partial / (partial p_(i j)) [p_(i j) log(p_(i j) / q_(i j)) + (1 / 2) p_(i j) log f_j] = 0$ yields $p_(i j)^* prop q_(i j) / sqrt(f_j)$ (attaining exact uniform marginal fixed point when $f_j = N / K$). If true class frequencies $pi$ differ from $u = (1 / K) bold(1)$, any partition whose empirical class frequencies match $u$ disagrees with the ground-truth counts by at least $sum_(c : pi_c > 1 / K) (pi_c - 1 / K) = (1 / 2) norm(pi - u)_1$ fraction of samples!
]

#subsec("B", "Recurrent, Pairwise & Pseudo-Label Discriminative Models (Algs. 186, 187, 191)")

Three deep clustering methods discard decoders entirely, coupling CNN representations directly to discrete clustering operators:
- *JULE (\#186)* @yang2016jule unrolls agglomerative clustering as a recurrent graph-merging network: at step $t$, it merges the highest-affinity cluster pair under Graph Degree Linkage (GDL) affinity $A(C_a, C_b) = abs(C_a)^(-1) sum_(i in C_a) (W_(i, C_b) W_(C_b, i)) / abs(C_b) + (a <-> b)$ (built on directed $k$-NN weights $W$), and backpropagates a weighted hyperspherical triplet-affinity loss into $z_i = (f_theta (x_i)) / norm(f_theta (x_i))_2 in bb(S)^(p - 1)$.
- *DAC (\#187)* @chang2017dac outputs normalized non-negative label features $l_i in bb(S)_+^(K - 1) := {v in bb(R)_(>= 0)^K : norm(v)_2 = 1}$, defines pairwise cosine similarities $s_(i j) = chevron.l l_i, l_j chevron.r in [0, 1]$, and trains binary cross-entropy on high-confidence pairs $r_(i j) = 1$ ($s_(i j) >= u(lambda)$) and $r_(i j) = 0$ ($s_(i j) < l(lambda)$) while annealing $u(lambda) arrow.b l(lambda)$.
- *DeepCluster (\#191)* @caron2018deepcluster alternates offline $k$-Means @lloyd1982least on PCA-whitened CNN features $z_i = f_theta (x_i)$ with supervised cross-entropy updates of $(theta, W_"cls", b_"cls")$ against pseudo-labels $y_i^((t))$, employing *uniform class resampling* and periodic head resets to combat empty-cluster collapse.

Theorem 4 proves exact algebraic identifiability and failure modes separating these three discriminative paradigms.

#theorem(title: "One-Hot DAC Identifiability, JULE Merge Amplification & DeepCluster Resampling Threshold", origin: "Original Unified Proof")[
  1. *(Exact One-Hot Identifiability & Concentrated Collapse of DAC \#187):* Let $l_1, dots, l_N in bb(S)_+^(K - 1)$ satisfy $chevron.l l_i, l_j chevron.r in {0, 1}$ for all $i, j$ and span $K$ distinct equivalence classes. Then the set of distinct vectors ${l_1, dots, l_N}$ equals the standard basis ${e_1, dots, e_K}$ *uniquely* (exact one-hot cluster indicators)! However, if initial pretrained features satisfy $min_(i, j) chevron.l l_i, l_j chevron.r >= u(lambda_0)$, *every* pair is labeled positive ($r_(i j) equiv 1$), turning the binary cross-entropy minimizer into single-cluster collapse $l_1 = dots = l_N = l_0$.
  2. *(Hyperspherical Non-Collapse & Irreversible Early-Merge Amplification of JULE \#186):* Constraining $z_i in bb(S)^(p - 1)$ prevents norm collapse, and GDL's product $W_(i, C_b) W_(C_b, i)$ vanishes identically on one-way outlier bridges. Nevertheless, because JULE's recurrent unrolling is strictly monotone ($|cal(C)_(t+1)| < |cal(C)_t|$ with no split operator), once an early step $t_0$ merges impure pair $(C_a, C_b)$ (containing true classes $1$ and $2$ with fraction $eta_"err"$), every subsequent affinity loss step $t > t_0$ enforces intra-cluster attraction across $C_a union C_b$, so the fraction of misclassified pairs grows monotonically and can never decrease below the step-$t_0$ impurity floor, while affinity recomputation scales as $O(N^2 d)$!
  3. *(Bias-Domination Collapse vs. Uniform-Resampling Equilibrium & Permutation Shock of DeepCluster \#191):*
     - Without uniform resampling, if $k$-Means assigns fraction $rho > 1 / K$ to a majority cluster $c^*$, a collapsed constant encoder $f_theta (x) equiv bold(0)$ with optimal bias $b_c = log hat(pi)_c$ attains cross-entropy $H(hat(pi)) < log K$; in the extreme $rho -> 1$, loss vanishes to $0$ at zero feature variance!
     - Under DeepCluster's *uniform class resampling* ($tilde(pi)_c = 1 / K$ for all $c$), *every* collapsed constant encoder achieves cross-entropy *at least* $log K$, strictly inferior to any representation with nonzero mutual information $I(Z; Y) > 0$ (whose Bayes cross-entropy is $log K - I(Z; Y) < log K$).
     - Yet because $k$-Means cluster indices are defined only up to arbitrary permutation $pi in S_K$ across epochs, keeping an unreset linear head $(W_"cls", b_"cls")$ across a cyclic relabeling $pi$ injects an immediate cross-entropy jump $Delta cal(L) >= 2 norm(W_"cls")_F^2 / K dot Omega(1)$ (gradient shock) unless the head is reinitialized or aligned via optimal permutation matching.
]

#proof[
  *Part 1:* Since $l_i, l_j in bb(S)_+^(K - 1)$ have unit $ell_2$ norm, $chevron.l l_i, l_j chevron.r = 1 <=> l_i = l_j$, and $chevron.l l_i, l_j chevron.r = 0 <=> "supp"(l_i) inter "supp"(l_j) = emptyset$ (because all coordinates are non-negative!). Thus $K$ distinct equivalence classes yield $K$ nonempty pairwise disjoint support subsets $S_1, dots, S_K subset.eq {1, dots, K}$. By the pigeonhole principle, $|S_k| = 1$ for every $k = 1, dots, K$. Unit norm plus non-negativity forces each $l_i$ to equal a standard basis vector $e_(sigma(k))$ @chang2017dac. Conversely, if $min_(i, j) chevron.l l_i, l_j chevron.r >= u(lambda_0)$, then $r_(i j) = 1$ for all $(i, j)$, so binary cross-entropy $-sum_(i, j) log chevron.l l_i, l_j chevron.r$ achieves its exact global minimum $0$ iff $l_i = l_j$ for all $i, j$ (complete single-cluster collapse).

  *Part 2:* On a directed $k$-NN graph where an outlier bridge links $C_a -> C_b$ without reciprocal reverse edges ($W_(C_b, i) = 0$), $W_(i, C_b) W_(C_b, i) = 0$, neutralizing one-way bridges @yang2016jule. Once $C_a union C_b$ merges at step $t_0$ without a split move, the partition hierarchy $cal(C)_t$ is strictly coarser than $cal(C)_(t_0)$ for all $t >= t_0$; hence every false same-cluster pair $(i, j) in C_a times C_b$ remains co-clustered forever, and JULE's positive-pair affinity gradient actively pulls $z_i$ and $z_j$ closer on $bb(S)^(p - 1)$, reducing inter-class margin monotonically.

  *Part 3:* When $f_theta (x) equiv bold(0)$, softmax logits reduce to constant vector $b in bb(R)^K$ with predicted probabilities $q_c = exp(b_c) / (sum_ell exp(b_ell))$. Under empirical class distribution $hat(pi)$, cross-entropy is $H(hat(pi)) + "KL"(hat(pi) || q) >= H(hat(pi))$, attained at $q = hat(pi)$; if $hat(pi)_(c^*) -> 1$, $H(hat(pi)) -> 0$ even though $f_theta$ carries zero information. Under uniform resampling ($tilde(pi)_c = 1 / K$), the identical identity forces constant-encoder cross-entropy to equal $log K + "KL"(u || q) >= log K$, whereas any non-degenerate classifier predicting $p(y | z)$ achieves $H(Y | Z) = H(Y) - I(Z; Y) = log K - I(Z; Y) < log K$ @caron2018deepcluster. Finally, if epoch $t$ aligns class prototype $w_c$ with mean representation $mu_c$ ($chevron.l w_c, mu_c chevron.r = gamma > 0$, $chevron.l w_c, mu_(c') chevron.r = -gamma / (K - 1)$ for $c' != c$), a derangement $pi in S_K$ ($pi(c) != c$) flips every sample's target logit from $+gamma$ to $-gamma / (K - 1)$, increasing cross-entropy instantaneously by $(K / (K - 1)) gamma > 0$!
]

#subsec("C", "Mutual Information, Fano Neighbors, Sinkhorn Equipartition & Contrastive Guarantees (Algs. 188--190, 192)")

Four self-supervised frameworks achieve collapse-resistant end-to-end clustering without decoders or offline $k$-Means:
- *IIC (\#188)* @ji2019iic forms the symmetrized empirical joint cluster-assignment matrix $P = (2 B)^(-1) sum_(i=1)^B (Phi(x_i) Phi(x'_i)^top + Phi(x'_i) Phi(x_i)^top)$ across augmented pairs $(x_i, x'_i)$ and maximizes mutual information
  $
    I(Z; Z') = sum_(c=1)^K sum_(c'=1)^K P_(c c') log (P_(c c') / (P_c P_(c'))).
  $
- *SCAN (\#189)* @vangansbeke2020scan decouples learning into Stage 1 (self-supervised pretext $k$-NN mining producing neighbor sets $cal(N)_i$), Stage 2 (minimizing neighbor-consistency loss with marginal entropy regularization):
  $
    cal(L)_"SCAN" (Phi) &= -1 / (N k) sum_(i=1)^N sum_(j in cal(N)_i) log chevron.l Phi(x_i), Phi(x_j) chevron.r \
    &quad - lambda H(overline(Phi)),
  $
  and Stage 3 (confidence-thresholded self-labeling).
- *SwAV (\#190)* @caron2020swav computes online *swapped* cross-entropy between augmented views $z_i^((1)), z_i^((2))$ and unit prototypes $C = [c_1, dots, c_K]$ using equipartitioned target codes $Q^* in cal(Q) := {Q in bb(R)_+^(K times B) : Q bold(1)_B = (1 / K) bold(1)_K, Q^top bold(1)_K = (1 / B) bold(1)_B}$ via entropic Sinkhorn--Knopp optimal transport @cuturi2013sinkhorn:
  $
    Q^* = arg max_(Q in cal(Q)) "Tr"(Q^top C^top Z) + epsilon H(Q).
  $
- *Contrastive Clustering (CC, \#192)* @li2021contrastive maps views $(x_i^a, x_i^b)$ to an instance matrix $Z^a, Z^b in bb(R)^(B times M)$ and a cluster probability matrix $Y^a, Y^b in [0, 1]^(B times K)$, simultaneously minimizing *row-wise* instance InfoNCE $cal(L)_"ins"$ @oord2018infonce, *column-wise* cluster InfoNCE $cal(L)_"clu"$ over cluster vectors $y_(:, c)^a, y_(:, c)^b in bb(R)^B$, and negative marginal entropy $-H(overline(y))$.

Theorem 5 delivers exact, unified mathematical guarantees explaining *why* each self-supervised mechanism prevents collapse—and its exact blind spots!

#theorem(title: "Mutual-Information Optimality, Exact Fano Threshold, Sinkhorn Contraction & Dual InfoNCE Bounds", origin: "Original Unified Proof")[
  Let $K >= 2$ and let $(x, x')$ denote positive semantic view pairs.
  1. *(IIC \#188 Automatic Collapse Rejection & Majority-Splitting Theorem):* Every degenerate constant assignment ($H(Z) = 0$) yields $I(Z; Z') = 0$, whereas the unique global maximum $I(Z; Z') = log K$ over symmetric joint distributions $P$ is attained *if and only if* $P = (1 / K) I_K$ (exact uniform deterministic clustering). However, on a ground-truth imbalanced distribution with class prior $pi != u$ ($H(pi) < log K$) and noiseless view invariance, the ground-truth partition scores only $I_"true" = H(pi) < log K$: splitting a majority class $c^*$ of mass $pi_(c^*) > 1 / K$ along an augmentation-invariant nuisance attribute raises $I(Z; Z')$ strictly above $H(pi)$, causing IIC to *over-split* majority classes!
  2. *(Exact Fano Collapse Threshold & Optimal Soft Calibration of SCAN \#189):* Suppose a balanced $K$-class pretext graph has uniform neighbor purity $alpha = bb(P)[y(x_j) = y(x_i) | j in cal(N)_i] in (1 / K, 1]$. Among symmetric calibrated semantic classifiers $Phi_gamma (x_i) = (1 - gamma) e_(y_i) + (gamma / K) bold(1)_K$ ($gamma in [0, 1]$), the Stage-2 loss $cal(L)_"SCAN" (Phi_gamma)$ is uniquely minimized at calibration smoothing
     $
       gamma^* = 1 - sqrt((K alpha - 1) / (K - 1)),
     $
     achieving exact Fano-entropy loss $cal(L)_"SCAN" (Phi_(gamma^*)) = H_"Fano" (alpha, K) - lambda log K$, where $H_"Fano" (alpha, K) := H_2 (alpha) + (1 - alpha) log(K - 1)$. Since trivial one-cluster collapse ($Phi equiv e_1$) scores $cal(L)_"SCAN" (e_1) = 0$, semantic clustering strictly defeats trivial collapse *if and only if*
     $
       lambda > lambda_"crit" (alpha, K) := (H_2 (alpha) + (1 - alpha) log(K - 1)) / (log K) < 1.
     $
     In particular, setting $lambda >= 1$ strictly defeats collapse for *every* better-than-random pretext graph ($alpha > 1 / K$).
  3. *(Sinkhorn Hilbert Contraction, Rank-1 Immunity & Zero-Entropy Limit of SwAV \#190):* Setting $K_epsilon = exp(C^top Z / epsilon) in bb(R)_(>0)^(K times B)$ with oscillation $Delta_epsilon := max_(i, j, k, ell) log ((K_(epsilon, i j) K_(epsilon, k ell)) / (K_(epsilon, i ell) K_(epsilon, k j)))$, Sinkhorn--Knopp scaling converges geometrically in Birkhoff's projective Hilbert metric at rate $tanh^2 (Delta_epsilon / 4) < 1$ to the *unique* optimal code $Q^* = "diag"(u^*) K_epsilon "diag"(v^*) in cal(Q)$. Under any rank-1 collapsed encoder ($z_i equiv z_0$), $Q^*$ collapses to the uniform matrix $(K B)^(-1) bold(1)_K bold(1)_B^top$ with swapped cross-entropy $log K$, strictly beaten by any clusterable representation; as $epsilon -> 0^+$, SwAV reduces *identically* to balanced spherical $k$-Means!
  4. *(Dual InfoNCE Mutual-Information Bound & Entropy-Penalty Necessity of CC \#192):* In Contrastive Clustering (\#192), the instance-level loss satisfies the exact mutual-information lower bound $I(X^a; X^b) >= log(2 B - 1) - cal(L)_"ins"$. For the cluster-level loss $cal(L)_"clu"$ at temperature $tau_C > 0$, orthogonal one-hot column assignments ($Y^a = Y^b in {0, 1}^(B times K)$ with nonempty columns) achieve exact cosine similarity $+1$ on positives and $0$ on negatives, attaining the unique minimum $log(1 + (2 K - 2) exp(-1 / tau_C))$. Yet *without* the negative entropy term $-H(overline(y))$, an assignment placing $B - (K - 1)$ samples in cluster $1$ and $1$ sample each in clusters $2, dots, K$ achieves the *identical global minimum* of $cal(L)_"clu"$, proving that $-H(overline(y))$ is mathematically indispensable to prevent imbalanced collapse!
]

#proof[
  *Part 1:* Since $I(Z; Z') = H(Z) - H(Z | Z') <= H(Z) <= log K$ with $H(Z | Z') >= 0$, constant collapse ($H(Z) = 0$) yields $I = 0$. Equality $I = log K$ requires $H(Z) = log K$ (uniform marginals $P_c = 1 / K$) and $H(Z | Z') = 0$ (each row/column has one nonzero entry $1 / K$); symmetry $P = P^top$ plus positive diagonal preference forces $P = (1 / K) I_K$ @ji2019iic. Under imbalanced truth $pi$ ($H(pi) < log K$), partitioning a majority class $c^*$ of mass $pi_(c^*)$ into two view-invariant sub-halves raises marginal entropy by $+pi_(c^*) log 2 > 0$ while keeping $H(Z | Z') = 0$, strictly raising $I(Z; Z')$ above $H(pi)$.

  *Part 2:* Let $s := (1 - gamma)^2 in [0, 1]$. For $Phi_gamma (x_i) = (1 - gamma) e_(y_i) + (gamma / K) bold(1)_K$, the marginal mean is $overline(Phi)_gamma = (1 / K) bold(1)_K$, so $H(overline(Phi)_gamma) = log K$. Direct inner-product expansion gives
  $
    chevron.l Phi_gamma (x_i), Phi_gamma (x_j) chevron.r = cases(
      s + (1 - s) / K & "if" y_i = y_j,
      (1 - s) / K & "if" y_i != y_j.
    )
  $
  Thus the expected neighbor term equals $-f(s)$ where $f(s) := alpha log(s + (1 - s) / K) + (1 - alpha) log((1 - s) / K)$. Differentiating:
  $
    f'(s) = (alpha (1 - 1 / K)) / (s + (1 - s) / K) - (1 - alpha) / (1 - s) = 0 \
    <=> s^* = (K alpha - 1) / (K - 1) in (0, 1],
  $
  valid iff $alpha > 1 / K$, giving $gamma^* = 1 - sqrt(s^*)$. Substituting $s^*$ back yields $s^* + (1 - s^*) / K = alpha$ and $(1 - s^*) / K = (1 - alpha) / (K - 1)$, whence
  $
    -f(s^*) &= -alpha log alpha - (1 - alpha) log((1 - alpha) / (K - 1)) \
    &= H_2 (alpha) + (1 - alpha) log(K - 1) = H_"Fano" (alpha, K).
  $
  Since $cal(L)_"SCAN" (e_1) = -log 1 - lambda dot 0 = 0$, $cal(L)_"SCAN" (Phi_(gamma^*)) < 0$ holds iff $lambda log K > H_"Fano" (alpha, K)$, proving the exact Fano threshold @vangansbeke2020scan. Because $H_"Fano" (alpha, K)$ is strictly decreasing on $[1 / K, 1]$ with $H_"Fano" (1 / K, K) = log K$, $lambda_"crit" (alpha, K) < 1$ for all $alpha > 1 / K$.

  *Part 3:* Because $Q |-> "Tr"(Q^top C^top Z) + epsilon H(Q)$ is strictly concave on compact polytope $cal(Q)$, KKT optimality yields unique $Q^* = "diag"(u^*) K_epsilon "diag"(v^*)$ with $u^* > 0, v^* > 0$ @cuturi2013sinkhorn. By Birkhoff's theorem on positive linear maps in projective Hilbert metric $d_H (u, tilde(u)) = log max_(i, j) ((u_i tilde(u)_j) / (u_j tilde(u)_i))$, applying $K_epsilon$ contracts $d_H$ by factor $tanh(Delta_epsilon / 4)$, so one full Sinkhorn round $(u -> v -> u)$ contracts by $tanh^2 (Delta_epsilon / 4) < 1$. If $z_i equiv z_0$, every column of $K_epsilon$ is identical ($K_(epsilon, k i) = a_k$), so $u_k^* = (K a_k)^(-1)$ and $v_i^* = B^(-1)$ satisfy all marginal constraints uniquely, forcing $Q^* = (K B)^(-1) bold(1)_K bold(1)_B^top$ and swapped cross-entropy $-sum_k (1 / K) log p_k >= log K$ (equality at uniform $p_k = 1 / K$). As $epsilon -> 0^+$, $epsilon H(Q) -> 0$ and the linear program $max_(Q in cal(Q)) "Tr"(Q^top C^top Z)$ selects extreme-point balanced assignments (when $B / K in bb(N)$), maximizing cosine similarity $chevron.l c_k, z_i chevron.r$ subject to equal cluster sizes—exact balanced spherical $k$-Means @caron2020swav @dhillon2001concept!

  *Part 4:* For instance InfoNCE with $1$ positive and $2 B - 2$ negatives ($2 B - 1$ total candidates), Oord et al.'s variational bound @oord2018infonce gives $I(X^a; X^b) >= log(2 B - 1) - cal(L)_"ins"$. For cluster-level InfoNCE on column vectors $y_(:, c)^a, y_(:, c)^b in [0, 1]^B$, cosine similarity obeys $s(u, v) in [0, 1]$ (since entries are non-negative!). Thus for every column $c$, the numerator $exp(s(y_(:, c)^a, y_(:, c)^b) / tau_C) <= exp(1 / tau_C)$ and each of the $2 K - 2$ negative terms obeys $exp(s(y_(:, c)^a, y_(:, c')^b) / tau_C) >= 1$, yielding the sharp lower bound
  $
    cal(L)_"clu" &>= -log (exp(1 / tau_C) / (exp(1 / tau_C) + (2 K - 2))) \
    &= log(1 + (2 K - 2) exp(-1 / tau_C)),
  $
  attained iff $y_(:, c)^a = y_(:, c)^b$ and $"supp"(y_(:, c)) inter "supp"(y_(:, c')) = emptyset$ for all $c != c'$. Crucially, cosine similarity is strictly *scale-invariant* in column $ell_2$-norm ($norm(y_(:, c))_2$ cancels identically!), so placing $B - K + 1$ points in column $1$ and $1$ point in each of columns $2, dots, K$ maintains orthogonal unit-normalized columns and achieves the *exact same* minimum of $cal(L)_"clu"$! Adding $-H(overline(y))$ strictly penalizes unequal column masses, forcing balanced clusters @li2021contrastive.
]

#remark(title: "The Deep Clustering Imbalance Dilemma")[
  Theorems 3--5 expose a universal trade-off: nearly every collapse-proof objective—DEPICT's $q_(i j) / sqrt(f_j)$, DeepCluster's uniform resampling, IIC's $log K$ peak, SCAN's entropy penalty, SwAV's Sinkhorn equipartition, and CC's $-H(overline(y))$—biases sizes toward $u = (1 / K) bold(1)_K$. Under severe imbalance ($H(pi) << log K$), only *IDEC (\#182)*, *DCN (\#183)*, *VaDE (\#184)* (learned $pi_c$), and *weighted Sinkhorn / SCAN* prevent majority splitting.
]

#wide_figure(
  caption: [*Table II:* Head-to-Head Theoretical & Algorithmic Comparison of All 12 Deep Representation & Self-Supervised Clustering Algorithms (Section II: Algs. 181--192).],
  {
    set text(size: 6.7pt)
    set table(inset: (x: 2.6pt, y: 2.0pt))
    table(
      columns: (1.06fr, 0.96fr, 0.96fr, 0.96fr, 0.90fr, 1.32fr),
      align: (left, left, left, left, left, left),
      [Algorithm], [Architecture / Signal], [Clustering Objective], [Collapse Safeguard], [Imbalance], [Distinctive Strength / Critical Weakness],
      [DEC (\#181) @xie2016dec], [Pretrained AE (no dec.)], [Student-$t$ $"KL"(P || Q)$], [*None* (0-loss collapse)], [Mod.], [Simple fine-tune; arbitrary manifold distortion],
      [IDEC (\#182) @guo2017idec], [Joint AE ($L_g$ decoder)], [$"MSE" + gamma "KL"(P || Q)$], [Lipschitz $g_psi$ + WD], [*Good*], [Preserves local geometry; needs $L_g$ / WD],
      [DCN (\#183) @yang2017dcn], [Joint AE + $k$-Means], [$"MSE" + (lambda / 2) "SSE"$], [Lipschitz $g_psi$ + WD], [*Good*], [Lyapunov alternating descent; scale-bypass risk],
      [VaDE (\#184) @jiang2017vade], [VAE + GMM prior], [Variational $"ELBO"$], [Generative $p(x | z)$], [*Native $pi_c$*], [Full generative density; posterior collapse SNR],
      [DEPICT (\#185) @ghasedi2017depict], [Denoising AE], [Balanced $"KL" +$ recon.], [$q_(i j) / sqrt(f_j)$ prior], [Poor ($u$ bias)], [Noise-robust + balanced; distorts $pi != u$],
      [JULE (\#186) @yang2016jule], [Recurrent CNN + GDL], [Agglom. triplet affinity], [Unit sphere $bb(S)^(p - 1)$], [*Good*], [Bridge-immune GDL; irreversible early merges, $O(N^2)$],
      [DAC (\#187) @chang2017dac], [CNN on $bb(S)_+^(K - 1)$], [Self-paced binary pairs], [Positive-orthant 1-hot], [Mod.], [Exact 1-hot basis theorem; dense init collapse],
      [IIC (\#188) @ji2019iic], [Siamese augmented CNN], [Max $I(Z; Z')$], [*Automatic* ($I = 0$)], [Poor (splits)], [Zero hyperparameters; splits majority classes],
      [SCAN (\#189) @vangansbeke2020scan], [3-stage pretext $k$-NN], [Neighbor agree $- lambda H$], [Fano $lambda > lambda_"crit"$], [Mod. (tunable)], [State-of-art semantics; pretext misalignment],
      [SwAV (\#190) @caron2020swav], [Online multi-crop CNN], [Swapped CE + Sinkhorn], [*Exact equipartition*], [Poor ($1/K$)], [Online scalable, no neg. pairs; uniform $Q^*$ bias],
      [DeepCluster (\#191) @caron2018deepcluster], [CNN + offline $k$-Means], [Pseudo-label CE], [Uniform resampling], [Mod.--Poor], [Simple strong baseline; epoch perm. shock],
      [CC (\#192) @li2021contrastive], [Dual-head Siamese CNN], [Row + col InfoNCE $- H$], [Dual MI + $-H(overline(y))$], [Poor ($-H$)], [End-to-end one-stage; needs $-H$ & large $B$],
    )
  },
)

// =============================================================================
// SECTION III
// =============================================================================
#sec("III", "Single-Pass & Anytime Data-Stream Engines")

Five algorithms in Category VIII (Algs. 193--197) process unbounded, non-stationary data streams under $O(1)$ or $O(log n)$ working memory and strict single-pass latency: CluStream (\#193) @aggarwal2003clustream, DenStream (\#194) @cao2006denstream, StreamKM++ (\#195) @ackermann2012streamkmpp, BIRCH-Stream / ClusTree (\#196) @kranen2011clustree, and D-Stream (\#197) @chen2007dstream. All five decouple fast online summarization from on-demand offline macro-clustering, yet diverge fundamentally across two axes: *temporal forgetting algebra* (Pyramidal snapshots vs. exponential decay vs. infinite window) and *spatial approximation guarantees* (fixed-radius / grid heuristics vs. provable coresets).

#subsec("A", "Temporal Decay Semigroups, Micro-Cluster Additivity & Density Caps (Algs. 193, 194, 196, 197)")

Let stream records $(x_tau, tau)$ arrive over discrete ticks $tau = 0, 1, dots, t$. Four engines maintain additive sufficient statistics or cell densities:
- *CluStream (\#193)* @aggarwal2003clustream maintains $q$ micro-clusters $"CF" = (overline("CF2"^x), overline("CF1"^x), "CF2"^t, "CF1"^t, n)$ stored across a *Pyramidal Time Frame* (PTF) retaining up to $alpha^l + 1$ snapshots of order $i$ (taken at multiples of $alpha^i$), answering arbitrary horizon queries $[t_c - h, t_c]$ via CF subtraction and offline weighted $k$-Means.
- *DenStream (\#194)* @cao2006denstream replaces fixed $q$ with exponentially damped micro-clusters ($f(Delta t) = 2^(-lambda Delta t)$, $lambda > 0$), partitioning summaries into *potential* micro-clusters $cal(P)$ ($w >= beta mu$) and *outlier* micro-clusters $cal(O)$ ($w < beta mu$) of bounded RMS radius $r_"rms" <= epsilon$, pruned every $T_p = ceil(lambda^(-1) log_2 ((beta mu) / (beta mu - 1)))$ ticks against threshold $xi(t_c, t_0) = (2^(-lambda (t_c - t_0 + T_p)) - 1) / (2^(-lambda T_p) - 1)$, with offline DBSCAN @ester1996dbscan on $cal(P)$.
- *BIRCH-Stream / ClusTree (\#196)* @kranen2011clustree organizes exponentially faded CFs $(n, bold("LS"), "SS")$ inside a height-balanced R-tree/CF-tree @zhang1996birch equipped with *internal hitchhiker buffers* $"CF"_e^("buf")$ that pause descent immediately upon stream interruption and carry buffered mass downward on future visits (*anytime* adaptation to arbitrary stream velocity).
- *D-Stream (\#197)* @chen2007dstream quantizes $bb(R)^d$ into $N = product_(j=1)^d m_j$ regular grid cells $g in cal(G)$ with discrete decay $D(g, t) = sum_(x in g, t_x <= t) lambda^(t - t_x)$ ($lambda in (0, 1)$), classifying cells as *dense* ($D >= D_m = C_m / (N (1 - lambda))$), *sparse* ($D <= D_l = C_l / (N (1 - lambda))$), or *transitional*, and extracting connected dense components every $t_g$ ticks.

Theorem 6 unifies the exact algebraic invariants, memory ceilings, and failure pathologies of all four temporal-window engines.

#theorem(title: "CF Semigroup Additivity, Exponential Weight Caps, Pruning Soundness & Grid Phase-Transition Gap", origin: "Original Unified Proof")[
  1. *(CluStream \#193 Exact Horizon Error, Pythagorean Macro-Loss & Unweighted Merge Breakdown):*
     - CluStream's CF tuples satisfy exact additivity/subtractivity across disjoint time intervals, and a PTF of order $(alpha, l)$ stores at most $(alpha^l + 1)(floor(log_alpha T) + 1) = O(alpha^l log_alpha T)$ snapshots while answering any horizon query $h <= T$ with relative horizon error at most $alpha^(-(l - 1))$.
     - Replacing each micro-cluster $C_j$ (mass $n_j$, centroid $mu_j$, internal SSE $V_j$) by its weighted centroid yields the exact Pythagorean decomposition
       $
         "SSE"_"full" (cal(M)) = sum_(j=1)^q V_j + sum_(j=1)^q n_j norm(mu_j - c_(pi(j)))_2^2.
       $
       Yet because CluStream merges the *closest unweighted* pair of micro-clusters when exceeding budget $q$, an adversarial arrival sequence forces arbitrarily large weighted-merge SSE blowup ($"SSE"_"CluStream" / "OPT"_k -> infinity$).
  2. *(Linear Decay Semigroup of ClusTree \#196 & DenStream \#194, Exact Buffer Ceilings & 4-$epsilon$ Non-Bridging):* Define $Phi_Delta ("CF") := 2^(-lambda Delta) "CF"$ ($Delta >= 0$). Then ${Phi_Delta}$ is a commutative linear semigroup preserving centroids $mu_(Phi_Delta ("CF")) = mu_("CF")$ and RMS radii $R_(Phi_Delta ("CF")) = R_("CF")$ identically, rendering *lazy decay* and *ClusTree hitchhiker buffering* algebraically exact at every visited node. At arrival rate $v >= 1$ points/tick:
     - Total decayed stream mass obeys $W(t) < v / (1 - 2^(-lambda))$ for all $t >= 0$, bounding DenStream's potential buffer deterministically independently of stream length:
       $
         |cal(P)| <= floor(v / (beta mu (1 - 2^(-lambda)))).
       $
     - DenStream's pruning curve $xi(t_c, t_0)$ satisfies the exact affine recurrence $a^m xi(t_c, t_0) + sum_(j=0)^(m-1) a^j = xi(t_c + m T_p, t_0) < 1 / (1 - a) <= beta mu$ ($a = 2^(-lambda T_p)$), proving that any pruned $o$-MC receiving $<= 1$ point per $T_p$ ticks can *never* reach $beta mu$.
     - Every $p$-MC of RMS radius $r_"rms" <= epsilon$ contains an actual witness $x_p in c_p$ with $norm(x_p - c_p)_2 <= epsilon$; hence two subsets $A, B subset.eq bb(R)^d$ separated by $"dist"(A, B) > 4 epsilon$ can *never* be falsely bridged by DenStream's offline DBSCAN!
  3. *(D-Stream \#197 Conservation, Exact Minimum Phase-Transition Gap & Dimensional Curse):* Under one point/tick and decay $lambda in (0, 1)$, total grid density satisfies $D_"tot" (t) = (1 - lambda^(t + 1)) / (1 - lambda) < 1 / (1 - lambda)$, hard-capping active cells independently of $t$:
     $
       |cal(G)_"dense" (t)| <= floor(N / C_m), wide
       |cal(G)_"ns" (t)| <= floor(N / C_l).
     $
     Setting offline inspection gap
     $
       t_g = floor(min(log_lambda (C_l / C_m), log_lambda ((N - C_m) / (N - C_l))))
     $
     is *necessary and sufficient* to guarantee no cell transitions between dense and sparse without passing through an offline pass @chen2007dstream. However, $N = product_(j=1)^d m_j$ grows exponentially in $d$, and axis-aligned grids break rotational invariance: rotating a thin slab by $45 degree$ inflates its covering cell count by $Omega(2^(d / 2))$!
]

#proof[
  *Part 1:* Additivity follows from coordinatewise linearity of sums $(n, sum x, sum x^2, sum t, sum t^2)$. In a PTF of order $(alpha, l)$, level $i$ keeps at most $alpha^l + 1$ snapshots spaced by $alpha^i$ up to level $floor(log_alpha T)$, bounding storage by $(alpha^l + 1)(floor(log_alpha T) + 1)$. If $t_c - h$ falls outside level $i - 1$ coverage, then $h > alpha^l dot alpha^(i - 1) = alpha^(l + i - 1)$, while level $i$ provides a snapshot within distance $<= alpha^i$, yielding relative error $<= alpha^i / alpha^(l + i - 1) = alpha^(-(l - 1))$ @aggarwal2003clustream. Expanding $norm(x - c_(pi(j)))_2^2 = norm(x - mu_j)_2^2 + 2 chevron.l x - mu_j, mu_j - c_(pi(j)) chevron.r + norm(mu_j - c_(pi(j)))_2^2$ and summing over $x in C_j$ cancels the cross term identically, yielding the Pythagorean decomposition. For the unweighted-merge breakdown ($k = 2, q = 3$ on $bb(R)$), place $n_1 = M$ points at $0$, $n_2 = M$ points at $1$, $n_3 = 1$ point at $4$, and then arrive $n_4 = 1$ point at $6$: the closest centroid pair is $(0, 1)$ (distance $1 < 2$), so CluStream merges the two massive clusters into $(2 M, 1 / 2)$ with irreducible internal SSE $V = M / 2 -> infinity$, whereas optimal 2-Means places centers at $0$ and $1$ with cost $<= (4 - 1)^2 + (6 - 1)^2 = 34$, giving approximation ratio $>= M / 68 -> infinity$.

  *Part 2:* Scalar linearity of $2^(-lambda Delta) > 0$ gives $Phi_(Delta_1) compose Phi_(Delta_2) = Phi_(Delta_1 + Delta_2)$ and exact cancellation of $2^(-lambda Delta)$ in $mu = bold("LS") / n$ and $R^2 = "SS" / n - norm(mu)_2^2$, proving inductive exactness of lazy updates and ClusTree hitchhiker pushes @kranen2011clustree. Summing $v dot 2^(-lambda (t - tau))$ over $tau = 0, dots, t$ gives $W(t) = v (1 - 2^(-lambda (t + 1))) / (1 - 2^(-lambda)) < v / (1 - 2^(-lambda))$; since each $c_p in cal(P)$ has weight $>= beta mu$ and points partition across $cal(P)$, $|cal(P)| <= floor(v / (beta mu (1 - 2^(-lambda))))$ @cao2006denstream. Writing $b = 2^(-lambda (t_c - t_0))$ gives $xi(t_c, t_0) = (1 - a b) / (1 - a)$, so
  $
    a^m xi(t_c, t_0) + sum_(j=0)^(m-1) a^j &= (a^m (1 - a b) + 1 - a^m) / (1 - a) \
    &= (1 - a^(m+1) b) / (1 - a) = xi(t_c + m T_p, t_0) < 1 / (1 - a) <= beta mu
  $
  because $T_p >= lambda^(-1) log_2 ((beta mu) / (beta mu - 1)) <=> 1 / (1 - a) <= beta mu$. Thus any $o$-MC starting below $xi$ with at most one arrival per $T_p$ ticks remains strictly below $beta mu$ forever! Finally, since $r_"rms"^2 (c_p) = sum_j (w_j / w) norm(x_j - c_p)_2^2 <= epsilon^2$, at least one $x_p in c_p$ satisfies $norm(x_p - c_p)_2 <= epsilon$. If $norm(c_p - c_q)_2 <= 2 epsilon$ (DBSCAN direct density reachability), the triangle inequality gives $norm(x_p - x_q)_2 <= epsilon + 2 epsilon + epsilon = 4 epsilon$, proving exact $4 epsilon$ non-bridging by contraposition.

  *Part 3:* Summing $lambda^(t - s)$ over $s = 0, dots, t$ yields $D_"tot" (t) = (1 - lambda^(t + 1)) / (1 - lambda) < (1 - lambda)^(-1)$. Markov counting against $D_m = C_m / (N (1 - lambda))$ gives $|cal(G)_"dense"| dot D_m < (1 - lambda)^(-1) => |cal(G)_"dense"| <= floor(N / C_m)$ (and similarly $floor(N / C_l)$ for non-sparse cells). For phase transitions: starting at $D >= D_m$ with zero arrivals, decaying to $<= D_l$ requires $lambda^Delta D_m <= D_l <=> Delta >= log_lambda (C_l / C_m) =: delta_0$; conversely, starting at $D <= D_l$ with *every* arrival landing in $g$ yields $D(g, t + Delta) = lambda^Delta D_l + (1 - lambda^Delta) / (1 - lambda) >= D_m <=> lambda^Delta <= (N - C_m) / (N - C_l) <=> Delta >= log_lambda ((N - C_m) / (N - C_l)) =: delta_1$. Hence $t_g = floor(min(delta_0, delta_1))$ catches every transition @chen2007dstream. Rotating a $d$-dimensional hyperplane slab of width $h$ onto the main diagonal $d^(-1 / 2) (1, dots, 1)$ intersects $Omega(2^(d / 2))$ unit hypercubes per unit cross-section vs. $1$ axis-aligned layer.
]

#subsec("B", "Coreset Approximation Guarantees vs. Concept-Drift Inertia (StreamKM++ vs. Faded Engines)")

While CluStream, DenStream, ClusTree, and D-Stream provide strong memory and drift control, *none* of the four offers a worst-case approximation ratio for global $k$-Means SSE. *StreamKM++ (\#195)* @ackermann2012streamkmpp fills this theoretical void via a Bentley--Saxe *Merge-and-Reduce* binary bucket cascade over chunks of size $m$, compressing $2 m |-> m$ weighted points using a randomized $D^2$-importance Coreset Tree @harpeled2004coresets @arthur2007kmeanspp. Theorem 7 proves the exact precision-vs.-adaptivity duality separating StreamKM++ from all four faded engines.

#theorem(title: "End-to-End Streaming Coreset Approximation vs. Concept-Drift Stale-Mass Inertia", origin: "Original Unified Proof")[
  Let $P = {x_1, dots, x_n} subset.eq bb(R)^d$ be a stream ingested in chunks of size $m$ ($M = floor(n / m) >= 1$), with binary popcount $s_2 (M)$ and tree height $h <= floor(log_2 M)$.
  1. *(Exact Carry Count, Coreset Composability & End-to-End $O(log k)$ Guarantee of StreamKM++ \#195):*
     - StreamKM++ holds at most $floor(log_2 M) + 2$ buckets ($O(m d log(n / m))$ space) and executes `Reduce` *exactly* $N_"red" (n) = M - s_2 (M) < n / m$ times across the entire stream.
     - If each `Reduce` call outputs a strong $(k, epsilon')$-coreset with $epsilon' = epsilon / (2 h)$ ($epsilon in (0, 1 / 2]$) and failure probability $<= delta m / n$, then with probability at least $1 - delta$ the root summary $S$ is a strong $(k, epsilon)$-coreset of $P$, and running weighted $k$-Means++ plus weighted Lloyd on $S$ yields end-to-end expected cost
       $
         bb(E)[Phi(P, hat(C))] <= (1 + epsilon) / (1 - epsilon) dot 8 (log k + 2) "OPT"_k (P).
       $
  2. *(Quantitative Concept-Drift Separation: Linear Stale-Mass Inertia of StreamKM++ \#195 vs. $O(lambda^(-1))$ Forgetting of Faded Engines):* Suppose the stream emits $n_1$ points from distribution $cal(D)_1$ ($k$ clusters at $mu_1, dots, mu_k$) and then abruptly switches to $n_2$ points from a disjoint shifted distribution $cal(D)_2$ separated by distance $Delta$.
     - In StreamKM++ (\#195), every coreset preserves total uniform weight $W_1 = n_1$ and $W_2 = n_2$ within $(1 plus.minus epsilon)$, so stale data retains fractional weight $>= ((1 - epsilon) n_1) / ((1 - epsilon) n_1 + (1 + epsilon) n_2)$: driving stale influence below $eta in (0, 1 / 2)$ requires $n_2 = Omega(n_1 / eta)$ post-drift arrivals (unboundedly slow as $n_1 -> infinity$)!
     - By exact contrast, in DenStream (\#194), ClusTree (\#196), and D-Stream (\#197) with decay rate $lambda > 0$ (or $beta = 2^(-lambda)$), stale fractional mass drops below $eta$ in at most
       $
         n_2 >= 1 / lambda log_2 (1 / eta)
       $
       post-drift ticks—*strictly independent* of pre-drift stream length $n_1$!
]

#proof[
  *Part 1:* After $M$ chunks, bucket $B_i$ ($i >= 1$) is occupied iff bit $i - 1$ of $M$ is $1$, requiring at most $floor(log_2 M) + 2$ buckets (including input buffer $B_0$). Each chunk adds $+1$ bucket and each `Reduce` merges $2 |-> 1$ ($-1$), so final occupancy $s_2 (M) = M - N_"red" (n) => N_"red" (n) = M - s_2 (M) < n / m$. By union bound over $< n / m$ calls, all reductions succeed with probability $>= 1 - delta$. Disjoint union preserves $(k, epsilon')$-coresets, and sequential composition across tree height $h$ gives $(1 - epsilon')^h Phi(P, C) <= Phi(S, C) <= (1 + epsilon')^h Phi(P, C)$. Since $epsilon' = epsilon / (2 h)$ with $epsilon in (0, 1 / 2]$, Bernoulli's inequality yields $(1 - epsilon / (2 h))^h >= 1 - epsilon / 2 > 1 - epsilon$, and $1 + u <= e^u$ yields $(1 + epsilon / (2 h))^h <= e^(epsilon / 2) <= 1 + epsilon$ (as $e^(epsilon / 2) - 1 - epsilon$ is convex with values $0$ at $0$ and $sqrt(e) - 2 < 0$ at $1$), proving $S$ is a strong $(k, epsilon)$-coreset @harpeled2004coresets @ackermann2012streamkmpp. Chaining Arthur--Vassilvitskii's $8 (log k + 2)$ bound @arthur2007kmeanspp on $S$, weighted Lloyd monotonicity $Phi(S, hat(C)) <= Phi(S, C^((0)))$, and the two-sided coreset inequalities gives
  $
    bb(E)[Phi(P, hat(C))] &<= 1 / (1 - epsilon) bb(E)[Phi(S, hat(C))] \
    &<= (8 (log k + 2)) / (1 - epsilon) "OPT"_k (S) \
    &<= (1 + epsilon) / (1 - epsilon) dot 8 (log k + 2) "OPT"_k (P).
  $

  *Part 2:* In StreamKM++, testing $C$ placing $k$ centers far from $cal(D)_1$ shows the coreset retains mass $(1 plus.minus epsilon) n_1$ on $cal(D)_1$ and $(1 plus.minus epsilon) n_2$ on $cal(D)_2$, so the stale fraction exceeds $((1 - epsilon) n_1) / ((1 - epsilon) n_1 + (1 + epsilon) n_2)$, requiring $n_2 >= ((1 - epsilon) / (1 + epsilon)) ((1 - eta) / eta) n_1 = Omega(n_1 / eta)$ to fall below $eta$. In exponentially decayed engines (rate $v = 1$, base $2^(-lambda)$), at $n_2$ ticks after drift, stale mass is $W_"old" = 2^(-lambda n_2) (1 - 2^(-lambda n_1)) / (1 - 2^(-lambda))$ while total mass is $W_"tot" = (1 - 2^(-lambda (n_1 + n_2))) / (1 - 2^(-lambda))$. Hence
  $
    W_"old" / W_"tot" = (2^(-lambda n_2) (1 - 2^(-lambda n_1))) / (1 - 2^(-lambda (n_1 + n_2))) < 2^(-lambda n_2) <= eta
  $
  for every $n_2 >= lambda^(-1) log_2 (1 / eta)$, strictly independent of $n_1$!
]

#remark(title: "High-Dimensional Radius Collapse & Anytime Velocity Adaptation")[
  Theorems 6 and 7 pinpoint every streaming trade-off: DenStream (\#194) and D-Stream (\#197) track non-convex drifting manifolds yet degrade in high $d$ (fixed-$epsilon$ concentration @beyer1999nearest and $N = product_j m_j$ cell blowup). ClusTree (\#196) adapts insertion depth to bursty stream rates via exact lazy hitchhiker CFs, whereas #box[StreamKM++] (\#195) alone guarantees end-to-end $O(log k)$ SSE accuracy (wrap in sliding-window resets under drift).
]

#wide_figure(
  caption: [*Table III:* Head-to-Head Theoretical & Algorithmic Comparison of Single-Pass & Anytime Data-Stream Engines (Section III: Algs. 193--197).],
  {
    set text(size: 7.2pt)
    set table(inset: (x: 3.1pt, y: 2.3pt))
    table(
      columns: (1.00fr, 0.94fr, 0.92fr, 0.88fr, 0.96fr, 1.30fr),
      align: (left, left, left, center, left, left),
      [Algorithm], [Temporal Window], [Summary / Geometry], [Drift / Anytime], [Online / Space], [Distinctive Strength / Critical Weakness],
      [CluStream (\#193) @aggarwal2003clustream], [Pyramidal PTF $[t_1, t_2]$], [$q$ MCs / spherical], [Horizon / No], [$O(q d)$ / $O(q d alpha^l log T)$], [Exact arbitrary horizon queries; unweighted merge blowup],
      [DenStream (\#194) @cao2006denstream], [Exp. decay $2^(-lambda Delta t)$], [$cal(P), cal(O)$ MCs / DBSCAN], [*Fast* / No], [$O((|cal(P)| + |cal(O)|) d)$], [Non-convex + $4 epsilon$ non-bridging; high-$d$ $epsilon$ collapse],
      [StreamKM++ (\#195) @ackermann2012streamkmpp], [Infinite uniform], [Merge-Reduce coreset], [None / No], [$O(m d)$ / $O(m d log(n/m))$], [*Exact $O(log k)$ SSE guarantee*; $Omega(n_1)$ drift inertia],
      [ClusTree (\#196) @kranen2011clustree], [Exp. decay $2^(-lambda Delta t)$], [Hitchhiker CF-tree], [*Fast / Native*], [$O(d M log_m K)$], [Self-adapts to bursty stream rate; greedy leaf splits],
      [D-Stream (\#197) @chen2007dstream], [Exp. decay $lambda^(Delta t)$], [Grid cells / density], [*Fast* / No], [$O(d)$ amort. / $O(|cal(G)_"act"|)$], [Fast low-$d$ non-convex + exact $t_g$; $O(2^d)$ & rot. break],
    )
  },
)

// =============================================================================
// SECTION IV
// =============================================================================
#sec("IV", "Global Stochastic Metaheuristics")

Because global Euclidean $k$-Means SSE minimization is NP-hard even for $k = 2$ or $d = 2$ @aloise2009np, Lloyd's algorithm @lloyd1982least and Hartigan--Wong relocation @hartigan1979algorithm routinely stall in poor local basins. The final three algorithms of Category VIII deploy population-based stochastic search to escape local minima: *EVO-KMeans / Genetic $K$-Means (GKA, \#198)* @krishna1999gka, *PSO-Clustering (\#199)* @vandermerwe2003pso, and *ACO-Clustering (\#200)* @shelokar2004aco. Crucially, they inhabit two distinct state spaces—discrete partition strings $s in cal(S)_K^+ subset.eq {1, dots, K}^n$ (GKA, ACO) vs. continuous centroid matrices $Z in bb(R)^(K times d)$ (PSO)—yielding sharply different global convergence and premature-stagnation guarantees.

#subsec("A", "Markov Global Convergence vs. Swarm Stability & Premature Stagnation")

- *Genetic $K$-Means (GKA, \#198)* @krishna1999gka evolves a population of $P$ partition strings $s in cal(S)_K^+$ via fitness-proportional selection, *distance-biased mutation* (independently mutating each allele $s_i$ with probability $p_m in (0, 1)$ to cluster $k$ with probability $p_(i, k) (s) prop c_m d_"max" (x_i) - d(x_i, c_k)$, $c_m > 1$), the deterministic one-step *$K$-Means Operator (KMO)* $cal(K)(s)$ replacing disruptive crossover, and strict elitism.
- *PSO-Clustering (\#199)* @vandermerwe2003pso evolves $S$ continuous particles $Z_i (t) in bb(R)^(K times d)$ with velocities $V_i (t)$, personal bests $P_i (t)$, and global best $G(t)$:
  $
    V_i (t + 1) &= w V_i (t) + c_1 r_(1, t) (P_i (t) - Z_i (t)) \
    &quad + c_2 r_(2, t) (G(t) - Z_i (t)), \
    Z_i (t + 1) &= Z_i (t) + V_i (t + 1),
  $
  optimizing either SSE or average quantization error (QE), optionally seeded by Lloyd's output (*Hybrid $K$-Means--PSO*).
- *ACO-Clustering (\#200)* @shelokar2004aco maintains an $n times K$ pheromone matrix $T^((t)) = [tau_(i j)^((t))]$ clipped to $[tau_"min", tau_"max"]$ (with evaporation $rho in (0, 1)$): $R$ ants construct partitions via pseudo-random proportional choice (greedy with probability $q_0 in [0, 1)$, proportional to $tau_(i j)^alpha eta_(i j)^beta$ otherwise), top-$L$ elites undergo deterministic local refinement $cal(L)$ (Lloyd or Hartigan--Wong), and elites deposit $Delta tau = Q / (J(tilde(s)^((ell))) + epsilon)$.

Theorem 8 delivers an airtight head-to-head mathematical verdict separating all three metaheuristics.

#theorem(title: "Markov Global Convergence of GKA/ACO vs. PSO Stability, Non-Stationary Collapse & Finite-Time Barriers", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ ($n >= K >= 2$) with exact global SSE optimum $J^* = min_(s in cal(S)_K^+) J(s)$ and optimal set $cal(S)^*$.
  1. *(Almost-Sure Markov-Chain Global Optimality of Elitist GKA \#198 & Bounded ACO \#200):*
     - In elitist GKA (\#198) with $p_m in (0, 1)$ and $c_m > 1$, KMO satisfies exact Lyapunov descent $J(cal(K)(s)) <= J(s)$ leaving every $s^* in cal(S)^*$ invariant, and uniform mutation probability obeys $p_(i, k) (s) >= (c_m - 1) / (K c_m) =: p_"min"^"GKA" > 0$. Setting $delta_0 := (min(1 - p_m, p_m p_"min"^"GKA"))^n > 0$ and $delta_P := 1 - (1 - delta_0)^P$, the finite Markov chain ${Z_t}$ obeys
       $
         "Pr"(J_t^* > J^*) <= (1 - delta_P)^t -> 0, wide
         J_t^* -> J^* "a.s."
       $
     - In ACO-Clustering (\#200) with $0 < tau_"min" <= tau_"max" < infinity$, $q_0 in [0, 1)$, and monotone local search $cal(L)$, setting $p_"min"^"ACO" := (1 - q_0) (tau_"min"^alpha eta_"min"^beta) / (tau_"min"^alpha eta_"min"^beta + (K - 1) tau_"max"^alpha eta_"max"^beta) > 0$, the $K!$ optimal labelings guarantee
       $
         "Pr"(J(s^(ast (T))) > J^*) <= exp(-R T K! (p_"min"^"ACO")^n),
       $
       yielding $J(s^(ast (T))) -> J^*$ almost surely with expected hitting time $bb(E)[tau_"hit"] <= (K! (p_"min"^"ACO")^n)^(-1)$.
  2. *(Exact Order-2 Stability Region & Non-Stationary Stagnation of Pure PSO \#199):* Under fixed attractors $p = g = z^*$ and symmetric acceleration $c_1 = c_2 = c > 0$ with inertia $w in [0, 1)$, canonical PSO-Clustering is mean-square contractive ($bb(E)[(z_t - z^*)^2] -> 0$) *if and only if*
     $
       0 < c < (12 (1 - w^2)) / (7 - 5 w).
     $
     Yet for *any* $(w, c_1, c_2) > 0$, canonical PSO admits absorbing stagnation points $Z_i (t) = P_i (t) = G(t) = Z_0$ ($V_i (t) = bold(0)$) where $J_"SSE"$ is smooth yet $nabla J_"SSE" (Z_0) != bold(0)$—failing even *first-order local stationarity* without hybrid Lloyd refinement!
  3. *(Finite-Time Failure Trilemma: Unclipped ACO Lock-In, $S_K$-Symmetry Washout & Exponential Escape Barrier):*
     - *(Unclipped ACO Lock-In):* Setting $tau_"min" = 0, tau_"max" = infinity$ in pure-pheromone ACO ($alpha = 1, beta = 0$) on the 4-point rectangle ${(0, 0), (0, delta), (1, 0), (1, delta)}$ with $delta in (1 / sqrt(2), 1)$ and $K = 2$ traps the colony *permanently* at suboptimal $s_"sub" = (1, 2, 1, 2)$ ($J = 1 > delta^2 = J^*$) with strictly positive probability $"Pr"(forall t >= 0 : J(s^(ast (t))) = 1) > 0$.
     - *(Label-Permutation Interference):* Because SSE is $S_K$-invariant, naive crossover of identical-partition parents $(1, 1, 2, 2)$ and $(2, 2, 1, 1)$ produces maximally corrupt offspring $(1, 1, 1, 1)$ (explaining why GKA replaces crossover with KMO!), while depositing $m$ full copies of $S_K$ in ACO yields identical entry deposits $(m (K - 1)! Q) / (J^* + epsilon)$, washing out pheromone contrast completely!
     - *(Exponential Basin-Escape Cost):* Escaping a strict local minimum separated by Hamming distance $r = Theta(n)$ from $cal(S)^*$ through independent per-coordinate perturbations requires expected generations $Omega((p_m p_"min")^(-r))$ in GKA and ACO unless local-search basins $cal(B)(cal(S)^*)$ overlap across population restarts.
]

#proof[
  *Part 1:* For GKA (\#198), Voronoi reassignment obeys $sum_i norm(x_i - c_(s'_i))_2^2 <= sum_i norm(x_i - c_(s_i))_2^2$, and the parallel-axis theorem inside each new cluster $C'_k$ yields exact descent gap $J(s) - J(cal(K)(s)) >= sum_k |C'_k| norm(c'_k - c_k)_2^2 >= 0$; hence global minimizers $s^*$ are fixed points of $cal(K)$. Since $c_m d_"max" (x_i) - d(x_i, c_k) >= (c_m - 1) d_"max" (x_i)$ while each denominator term is $<= c_m d_"max" (x_i)$, $p_(i, k) (s) >= (c_m - 1) / (K c_m) = p_"min"^"GKA" > 0$. Independent coordinate mutation transitions *any* string $s$ to $s^*$ with probability $>= delta_0 = (min(1 - p_m, p_m p_"min"^"GKA"))^n > 0$; across $P$ population members, hitting the absorbing optimal class $Omega^*$ in one generation has probability $>= delta_P = 1 - (1 - delta_0)^P > 0$. By elitism, $"Pr"(J_t^* > J^*) <= (1 - delta_P)^t$, whose summability $sum_t (1 - delta_P)^t = delta_P^(-1) < infinity$ yields $J_t^* -> J^*$ a.s. via the first Borel--Cantelli lemma @krishna1999gka.

  For bounded ACO (\#200), clipping $tau_(i j) in [tau_"min", tau_"max"]$ ensures every ant selects any target label $j$ at coordinate $i$ with probability $>= p_"min"^"ACO" > 0$ when exploring ($u > q_0$). Summing over all $K!$ distinct permutations in $cal(S)^*$ gives single-ant hit probability $>= K! (p_"min"^"ACO")^n$. Across $R$ independent ants and $T$ elitist iterations with monotone $cal(L)$, $"Pr"(J(s^(ast (T))) > J^*) <= (1 - K! (p_"min"^"ACO")^n)^(R T) <= exp(-R T K! (p_"min"^"ACO")^n)$, giving Borel--Cantelli almost-sure convergence and geometric hitting-time expectation $<= (K! (p_"min"^"ACO")^n)^(-1)$ @shelokar2004aco.

  *Part 2:* Eliminating velocity from the scalar PSO recurrence at $p = g = z^*$ ($y_t = z_t - z^*$) yields $y_(t+1) = alpha_t y_t - w y_(t-1)$ where $alpha_t = 1 + w - phi_t$, $phi_t = c (r_(1, t) + r_(2, t))$ with mean $overline(phi) = c$, variance $sigma_phi^2 = c^2 / 6$, and $a := bb(E)[alpha_t] = 1 + w - c$. Setting $u_t = bb(E)[y_t^2]$ and $v_t = bb(E)[y_t y_(t-1)]$ yields the exact linear second-moment system
  $
    u_(t+1) &= (a^2 + sigma_phi^2) u_t + w^2 u_(t-1) - 2 w a v_t, \
    v_(t+1) &= a u_t - w v_t,
  $
  whose $3 times 3$ transition matrix has characteristic polynomial $p(lambda) = lambda^3 + (w - a^2 - sigma_phi^2) lambda^2 + (w a^2 - w sigma_phi^2 - w^2) lambda - w^3$. By Jury's Schur-stability criterion ($|w| < 1$), $p(-1) < 0$ holds identically while
  $
    p(1) = (1 - w)((1 + w)^2 - a^2) - (1 + w) sigma_phi^2 > 0 \
    <=> c^2 / 6 < (1 - w) / (1 + w) c (2 + 2 w - c) \
    <=> c < (12 (1 - w^2)) / (7 - 5 w).
  $
  For non-stationary stagnation, take $d = 1, K = 1, cal(X) = {-1, 1\}$ (unique global minimum at $z^* = 0$) and initialize $Z_i (0) = P_i (0) = G(0) = 1$, $V_i (0) = 0$ for all $i$: then $P_i - Z_i = G - Z_i = 0$ identically, so $V_i (t) equiv 0$ and $Z_i (t) equiv 1$ for all $t >= 0$ even though $nabla J_"SSE" (1) = 4 != 0$ @vandermerwe2003pso!

  *Part 3:* On the 4-point rectangle with $delta in (1 / sqrt(2), 1)$, $s^* = (1, 1, 2, 2)$ gives $J^* = delta^2 < 1$ while $s_"sub" = (1, 2, 1, 2)$ gives $J = 1$ and is strictly stationary under both Lloyd and Hartigan--Wong (relocating any vertex yields $J' = (2 / 3)(1 + delta^2) > 1$). Let $E_t = {forall r <= R, s^((r, t)) = s_"sub"}$; $"Pr"(E_0) >= ((1 - q_0) / 2)^(4 R) > 0$. Conditioned on $E_0 inter dots inter E_(t-1)$, unclipped pheromones satisfy $tau_(i, s_("sub", i))^((t)) = (1 - rho)^t tau_0 + (1 - (1 - rho)^t) c$ ($c := L Q / (rho (1 + epsilon))$) versus $tau_(i, "other")^((t)) = (1 - rho)^t tau_0$. Thus per-coordinate deviation probability decays geometrically: $q_t <= ((1 - q_0) tau_0) / (rho c) (1 - rho)^t$. Since $1 - u >= e^(-2 u)$ on $[0, 1 / 2]$,
  $
    "Pr"(inter.big_(t=1)^infinity E_t mid(|) E_0) &= product_(t=1)^infinity (1 - q_t)^(4 R) \
    &>= exp(-(8 R (1 - q_0) tau_0 (1 - rho)) / (rho^2 c)) > 0,
  $
  proving permanent lock-in with strictly positive probability! For $S_K$ symmetry washout, exactly $(K - 1)!$ permutations assign item $i$ to label $j$, giving identical deposit $(m (K - 1)! Q) / (J^* + epsilon)$ on every entry $(i, j)$ and leaving uniform $T^((t))$ strictly uniform. Finally, if a local optimum requires flipping $r$ coordinates simultaneously to exit its KMO/local-search attraction basin, the single-trial escape probability via independent per-allele perturbation is bounded by $(p_m)^r$ (modulo combinatorial neighborhood factor), explaining the exponential worst-case hitting barrier consistent with NP-hardness @aloise2009np.
]

#remark(title: "Hybridization Imperative for Metaheuristic Clustering")[
  Theorem 8 yields a clear rule: metaheuristics must *wrap* fast local operators rather than replace them. GKA's KMO, ACO's Hartigan--Wong refinement + clipping ($tau_"min" > 0$) + elite canonicalization, and Hybrid $k$-Means++--PSO convert slow random walks into basin-hopping over Voronoi-stationary points.
]

#wide_figure(
  caption: [*Table IV:* Head-to-Head Theoretical & Algorithmic Comparison of Global Stochastic Metaheuristics (Section IV: Algs. 198--200).],
  {
    set text(size: 7.3pt)
    set table(inset: (x: 3.3pt, y: 2.4pt))
    table(
      columns: (1.02fr, 0.90fr, 0.94fr, 0.94fr, 0.96fr, 1.24fr),
      align: (left, left, left, center, left, left),
      [Algorithm], [State Space], [Escape Mechanism], [Global Conv.], [Time Complexity], [Distinctive Strength / Failure Mode],
      [GKA (\#198) @krishna1999gka], [Discrete $s in cal(S)_K^+$], [Dist.-biased mut. + KMO], [*A.s. (Markov)*], [$O(G P n K d)$], [Exact Markov global opt.; $O(n)$ string / exp. hit time],
      [PSO (\#199) @vandermerwe2003pso], [Cont. $Z in bb(R)^(K times d)$], [Momentum + $P_i, G$ pull], [*No* (stagnates)], [$O(T S n K d)$], [Compact $K d$ state + hybrid opt.; $nabla J != 0$ collapse],
      [ACO (\#200) @shelokar2004aco], [Pheromone $T in bb(R)_+^(n times K)$], [Stigmergic $T$ + local $cal(L)$], [*A.s. ($tau_"min" > 0$)*], [$O(T n (R (K + d) + L I K d))$], [Black-box fitness + parallel; $S_K$ washout & $tau_"min" = 0$ lock],
    )
  },
)

// =============================================================================
// SECTION V
// =============================================================================
#sec("V", "Complete 25-Algorithm Matrix & Practitioner's Decision Guide")

Table V synthesizes *all twenty-five* algorithms of Category VIII (Algs. 176--200), followed by exact parametric reductions linking Part VIII across the entire monograph, an operational decision guide grounded directly in Theorems 1--8, and four open mathematical frontiers.

#wide_figure(
  caption: [*Table V:* Master Head-to-Head Theoretical and Practical Comparison Matrix of All 25 Neural, Deep, Streaming, and Evolutionary Clustering Algorithms (Category VIII, Algorithms 176--200). Subfamily codes: *CN* = Competitive Neural, *DC* = Deep Clustering, *ST* = Streaming, *EV* = Evolutionary / Swarm.],
  {
    set text(size: 6.25pt)
    set table(inset: (x: 2.1pt, y: 1.50pt))
    table(
      columns: (1.18fr, 0.32fr, 0.90fr, 0.98fr, 0.98fr, 0.56fr, 1.48fr),
      align: (left, center, left, left, left, center, left),
      [Algorithm], [Subf.], [State / Geometry], [Core Mechanism], [Time Complexity], [Auto $K$], [Distinctive Strength / Critical Weakness],
      [\#176 SOM @kohonen1982som], [CN], [Fixed $k$-D grid], [Lattice kernel $h_(sigma_t)$], [$O(E n M d)$], [No], [2D topo. map; no online potential, folding & voids],
      [\#177 GNG @fritzke1995gng], [CN], [Induced Delaunay], [Hebbian edge + age $a_"max"$], [$O(T M_"max" d)$], [*Yes*], [Adapts arbitrary $d_"int"$ & comp.; drift if $epsilon > 0$],
      [\#178 ART @carpenter1991fuzzyart], [CN], [Comp. hyperboxes], [Vigilance $|R_j|_1 <= (1 - rho) d$], [$O(T K d)$], [*Yes*], [*Zero catastrophic drift*; order & $rho -> 2/3$ conc.],
      [\#179 Neural Gas @martinetz1993neuralgas], [CN], [Rank-order VQ], [Exact gradient $E_"NG" (W; lambda)$], [$O(T M (d + log M))$], [No], [Exact $C^1$ potential + no dead units; no graph],
      [\#180 GCS @fritzke1994gcs], [CN], [Pure $k$-simplices], [Error edge-split ($chi = 1$)], [$O(T (M_"max" d + k^2 / lambda))$], [No ($M$)], [Euler-invariant mesh; fixed intrinsic dim. $k$],
      [\#181 DEC @xie2016dec], [DC], [Latent Student-$t$], [Self-training $"KL"(P || Q)$], [$O(E N (W_"enc" + K p))$], [No], [Simple AE fine-tuning; 0-loss collapse w/o dec.],
      [\#182 IDEC @guo2017idec], [DC], [Joint AE latent], [$"MSE" + gamma "KL"(P || Q)$], [$O(E N (W_"AE" + K p))$], [No], [Lipschitz non-collapse; needs WD against scale-bypass],
      [\#183 DCN @yang2017dcn], [DC], [Joint AE + Voronoi], [$"MSE" + (lambda / 2) "SSE"$], [$O(E N (W_"AE" + K p))$], [No], [Lyapunov alt. descent; scale-bypass & $1/N$ outliers],
      [\#184 VaDE @jiang2017vade], [DC], [VAE + GMM prior], [Variational $"ELBO"$], [$O(E N (P + K d))$], [No], [Full generative density + $pi_c$; posterior collapse],
      [\#185 DEPICT @ghasedi2017depict], [DC], [Denoising AE], [Balanced $p_(i j) prop q_(i j) / sqrt(f_j)$], [$O(E N (P + K d))$], [No], [Noise-robust + balanced; distorts $pi != u$],
      [\#186 JULE @yang2016jule], [DC], [Hypersphere + GDL], [Recurrent agglom. CNN], [$O(T (N P + N^2 d))$], [No], [Bridge-immune GDL; irreversible early merges, $O(N^2)$],
      [\#187 DAC @chang2017dac], [DC], [$bb(S)_+^(K - 1)$ label cone], [Self-paced binary pairs], [$O(E (N P + N B K))$], [No], [Exact 1-hot identifiability; dense init collapse],
      [\#188 IIC @ji2019iic], [DC], [Augmented views], [Max mutual info $I(Z; Z')$], [$O(E N (C_"net" + K^2))$], [No], [Automatic 0-MI collapse immunity; splits majority],
      [\#189 SCAN @vangansbeke2020scan], [DC], [Pretext $k$-NN graph], [Neighbor agree $- lambda H$], [$O(N k log N + E N C_"net")$], [No], [Exact Fano $lambda > lambda_"crit"$ + top vision acc.; pretext risk],
      [\#190 SwAV @caron2020swav], [DC], [Multi-crop sphere], [Swapped CE + Sinkhorn], [$O(E N (C_"net" + M K L_"sk"))$], [No], [Online equipartition, no neg. pairs; uniform $Q^*$ bias],
      [\#191 DeepCluster @caron2018deepcluster], [DC], [CNN + $k$-Means], [Pseudo-label CE + resamp.], [$O(E (N C_"net" + I N K p))$], [No], [Uniform resampling immunity; epoch perm. shock],
      [\#192 CC @li2021contrastive], [DC], [Dual Siamese heads], [Row + col InfoNCE $- H$], [$O(E N (C_"enc" + N m + K^2))$], [No], [End-to-end 1-stage MI; needs $-H$ & large batch $B$],
      [\#193 CluStream @aggarwal2003clustream], [ST], [$q$ MCs + PTF], [Horizon CF subtraction], [$O(n q d + I_"off" q K d)$], [No], [Exact arbitrary horizon queries; unweighted merge fail],
      [\#194 DenStream @cao2006denstream], [ST], [Faded $cal(P), cal(O)$ MCs], [Decay + $xi$ prune + DBSCAN], [$O(n (|cal(P)| + |cal(O)|) d)$], [*Yes*], [Non-convex + $4 epsilon$ non-bridging; high-$d$ $epsilon$ collapse],
      [\#195 StreamKM++ @ackermann2012streamkmpp], [ST], [Binary coreset tree], [Merge-Reduce $D^2$ coreset], [$O(n m d + I_"off" m K d)$], [No], [*End-to-end $O(log k)$ SSE*; $Omega(n_1)$ drift inertia],
      [\#196 ClusTree @kranen2011clustree], [ST], [Faded CF-tree], [Lazy hitchhiker anytime], [$O(d M log_m K)$ / pt], [No], [Adapts to bursty stream speed; greedy leaf splits],
      [\#197 D-Stream @chen2007dstream], [ST], [Faded dense grid], [Grid decay + gap $t_g$], [$O(d + (d abs(cal(G)_"act")) / t_g)$], [*Yes*], [Fast low-$d$ non-convex + exact $t_g$; $O(2^d)$ & rot. fail],
      [\#198 GKA @krishna1999gka], [EV], [Discrete $s in cal(S)_K^+$], [Dist. mutation + KMO], [$O(G P n K d)$], [No], [*Almost-sure Markov global opt.*; exp. worst-case time],
      [\#199 PSO @vandermerwe2003pso], [EV], [Cont. $Z in bb(R)^(K times d)$], [Velocity + $P_i, G$ swarm], [$O(T S n K d)$], [No], [Compact $K d$ state + hybrid opt.; $nabla J != 0$ stagnation],
      [\#200 ACO @shelokar2004aco], [EV], [Pheromone $T_(n times K)$], [Stigmergic $T$ + local $cal(L)$], [$O(T n (R (K + d) + L I K d))$], [No], [*A.s. opt. ($tau_"min" > 0$)* + black-box; $S_K$ washout],
    )
  },
)

#subsec("A", "Exact Parametric Reductions & Cross-Paradigm Equivalence Lattice")

Theorems 1--8 interlock all twenty-five Category VIII methods into exact limiting bridges across Parts I--VII of this monograph:
- *Competitive Neural Limits:* SOM (\#176, $sigma_t -> 0^+$) and Neural Gas (\#179, $lambda -> 0^+$) reduce identically to MacQueen's online $k$-Means @macqueen1967some (batch variants to Lloyd's $k$-Means @lloyd1982least, Alg. 1); Fuzzy ART (\#178) on binary inputs ${0, 1}^d$ reduces identically to binary ART1 @carpenter1987art; fixed-vertex GNG (\#177, $lambda_"ins" = infinity$) reduces to Competitive Hebbian Delaunay triangulation; 1-simplex GCS (\#180, $k = 1$) reduces to a growing topological chain.
- *Deep Clustering Reductions:* IDEC (\#182, $gamma -> infinity$) recovers DEC (\#181), while $gamma -> 0^+$ recovers a standard autoencoder; DCN (\#183) with frozen identity encoder $f_theta = I$ reduces to Lloyd's $k$-Means (and with orthonormal linear encoder to PCA + $k$-Means); VaDE (\#184) with linear isotropic decoder collapses to probabilistic PCA + GMM EM @dempster1977em (Alg. 81); SwAV (\#190, $epsilon -> 0^+$) reduces identically to balanced spherical $k$-Means @dhillon2001concept (Alg. 12), and removing Sinkhorn recovers online DeepCluster (\#191); SCAN (\#189) relaxes continuously to soft normalized-cut spectral clustering @shi2000ncut @vonluxburg2007tutorial on the pretext $k$-NN graph; JULE (\#186) without backpropagation recovers GDL agglomerative clustering (sharing hierarchical neighbor-merge affinity with FINCH @sarfraz2019finch); CC (\#192) in the high-temperature limit recovers normalized-cut Gram alignment.
- *Streaming & Metaheuristic Bridges:* CluStream (\#193, infinite horizon) and ClusTree (\#196, $lambda = 0$, no interruption) reduce algebraically to leaf-level BIRCH @zhang1996birch (Alg. 60, and with coreset leaf radii to BICO @fichtenberger2013bico); DenStream (\#194, $lambda -> 0^+$) recovers static DBSCAN @ester1996dbscan (Alg. 31) on compressed $p$-MCs; D-Stream (\#197, $lambda -> 1^-$) recovers grid-connected CLIQUE / density-grid clustering @agrawal1998clique; StreamKM++ (\#195, $m >= n$) recovers exact $k$-Means++ @arthur2007kmeanspp; GKA (\#198, $p_m = 0, P = 1$) and ACO (\#200, $R = 1, q_0 = 1$) collapse to deterministic Lloyd / Hartigan--Wong local descent, while ACO's exponentially decayed elite pheromone matrix acts as an online discounted cluster ensemble!

#subsec("B", "Mathematically Grounded Practitioner's Decision Guide")

1. *Topological Visualization vs. VQ:* Choose *SOM (\#176)* (Heskes Batch for Lyapunov descent, Thm. 1.1) only for 2D interpretable maps. Prefer *Neural Gas (\#179)* for distortion-optimal VQ without folding, *GNG (\#177)* when intrinsic dimension or component count is unknown (Thm. 1.3), and *GCS (\#180)* only for pure $k$-simplicial meshes.
2. *Continual Learning Without Catastrophic Drift:* Reject constant-rate VQ (Thm. 2.1) and select complement-coded *Fuzzy ART (\#178)* for strict memory stability and automatic spawning—verifying $norm(I)_1 = d$ and tuning $rho$ (after PCA when $d > 50$) against Thm. 2.3's $rho = 2 / 3$ concentration cliff.
3. *Tabular & Imbalanced Deep Clustering:* Never run unregularized *DEC (\#181)* without collapse checks (Thm. 3.1). Default to *IDEC (\#182)* or *DCN (\#183)* *with decoder weight decay or spectral norm* (Thm. 3.2 scale-bypass cure). Use *VaDE (\#184)* (with KL warm-up against Thm. 3.3 posterior collapse) for generative sampling or learned imbalanced priors $pi_c$, and *DEPICT (\#185)* on noisy balanced data.
4. *Large-Scale Vision & Self-Supervised Semantics:* Avoid *JULE (\#186)* ($O(N^2)$ and irreversible merges, Thm. 4.2) and *DAC (\#187)* at scale. Default to *SCAN (\#189)* with $lambda >= 1$ (above Thm. 5.2's Fano threshold for all $alpha > 1 / K$) on strong pretext backbones. Choose *SwAV (\#190)* for online multi-crop learning without negative pairs (Thm. 5.3), or *Contrastive Clustering (\#192)* (retaining $-H(overline(y))$, Thm. 5.4) / *IIC (\#188)* for compact 1-stage training on near-balanced classes (Thm. 5.1).
5. *Historical Horizons vs. Non-Convex Drifting Streams:* Use *CluStream (\#193)* for arbitrary past horizons $[t_1, t_2]$ (Thm. 6.1), upgrading online merges to mass-weighted Ward costs. For non-convex drifting streams with noise, choose *DenStream (\#194)* ($d <= 20$, Thm. 6.2's $4 epsilon$ non-bridging) or *D-Stream (\#197)* ($d <= 5$, Thm. 6.3's exact gap $t_g$).
6. *Bursty Real-Time vs. Provably Accurate Streams:* Use *ClusTree (\#196)* under bursty anytime arrival budgets (Thm. 6.2 lazy hitchhiker exactness), and *#box[StreamKM++] (\#195)* for provable $O(log k)$ SSE on stationary big streams (Thm. 7.1)—wrapping in sliding resets under concept drift (Thm. 7.2).
7. *Global Non-Convex SSE Metaheuristics ($n <= 10^4$):* Never run pure *PSO (\#199)* alone (Thm. 8.2 stagnation); deploy *Hybrid $k$-Means++--PSO* inside $c < (12 (1 - w^2)) / (7 - 5 w)$. Choose *Elitist GKA (\#198)* for discrete Markov optimality with KMO, and *Bounded ACO (\#200)* ($tau_"min" > 0, q_0 < 1$, canonicalized elites against Thm. 8.3 washout) for black-box side constraints.
8. *Cross-Subfamily Pipelines:* Compose primitives directly: pretext *SwAV (\#190)* into *SCAN (\#189)*; *ClusTree (\#196)* anytime summaries into *#box[StreamKM++] (\#195)* coresets; or *IDEC (\#182)* latent spaces refined by *Hybrid PSO (\#199)*.

#subsec("C", "Open Theoretical Frontiers & Bridges to Part IX Synthesis")

Four open problems bridge Category VIII to Chapters XI--XIII:
1. *Imbalance-Aware Collapse-Free Self-Supervision:* Replacing uniform Sinkhorn (SwAV), $log K$ MI (IIC), and $-H(overline(y))$ (SCAN/CC) with collapse-free semi-relaxed OT or Dirichlet bounds recovering arbitrary long-tailed priors $pi$ without majority splitting.
2. *Drift-Adaptive Streaming Coresets:* Uniting #box[StreamKM++]'s $(1 plus.minus epsilon)$ distortion (Thm. 7.1) with $O(lambda^(-1))$ exponential forgetting (Thm. 7.2) in $O(m d log(1 / lambda))$ space for provable $O(log k)$ faded-window clustering in high $d$.
3. *Smoothed Basin-Transition Theory for Metaheuristics:* Replacing worst-case $O(p_"min"^(-n))$ Markov bounds (Thm. 8.1) with polynomial perturbation-resilient hitting times for KMO-GKA, canonicalized ACO, and Hybrid PSO.
4. *Continual Self-Supervised Topological Networks:* Unifying Fuzzy ART's vigilance-bounded zero forgetting (Thm. 2.2), GNG's induced Delaunay graph (Thm. 1.3), and deep encoders into continual clustering immune to both collapse and drift.
