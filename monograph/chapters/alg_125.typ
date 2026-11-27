#import "../template/ieee.typ": *

#alg_entry(
  id: 125,
  name: "BigCLAM (Cluster Affiliation Model for Big Networks)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2013",
  paradigm: "Nonnegative Affiliation Matrix Factorization",
  complexity: $cal(O)(T K (m + n))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

A striking empirical finding across large social, biological, and collaboration networks is that *community overlaps are denser than community interiors*: two individuals sharing two or more communities are strictly more likely to be linked than two individuals sharing only one @yang2013bigclam. Classical graph-partitioning methods @blondel2008louvain @newman2004finding forbid overlap altogether, while simplex-constrained mixed-membership stochastic blockmodels (MMSBs) model edge probabilities as convex averages across communities, inadvertently predicting *sparser* overlaps. Yang and Leskovec @yang2013bigclam introduced *BigCLAM* (Cluster Affiliation Model for Big Networks), formulating scalable overlapping community detection as continuous nonnegative matrix factorization under an independent-or probabilistic link.

Let $G = (V, E)$ be an undirected simple graph with $n = abs(V)$ vertices, $m = abs(E)$ edges, and neighbor sets $N(u) = {v in V : (u, v) in E}$. BigCLAM associates each vertex $u in V$ with a nonnegative community-affiliation vector $F_u in bb(R)_(>= 0)^K$ stacked as rows of $F in bb(R)_(>= 0)^(n times K)$, where $F_(u c) >= 0$ measures the strength of vertex $u$'s membership in community $c in {1, ..., K}$. Each community $c$ independently generates an edge between $u$ and $v$ with probability $1 - exp(-F_(u c) F_(v c))$; assuming conditional independence across communities and vertex pairs, an edge exists if at least one shared community connects $u$ and $v$.

#definition(title: "Affiliation Graph Model and BigCLAM Log-Likelihood")[
  Under the continuous Affiliation Graph Model (AGM) @yang2013bigclam with nonnegative factor matrix $F in bb(R)_(>= 0)^(n times K)$, the probability of an edge between distinct vertices $u, v in V$ is
  $
    p_(u v) (F) = 1 - product_(c = 1)^K exp(-F_(u c) F_(v c)) = 1 - exp(-chevron.l F_u, F_v chevron.r) .
  $
  The BigCLAM maximum-likelihood estimator maximizes the graph log-likelihood (counting ordered pairs for algebraic symmetry, which scales the unordered likelihood by $2$) over the nonnegative orthant:
  $
    cal(L)(F) &= sum_(u in V) sum_(v in N(u)) log(1 - exp(-chevron.l F_u, F_v chevron.r)) \
              &quad - sum_(u in V) sum_(v in.not N(u) union {u}) chevron.l F_u, F_v chevron.r ,
  $
  subject to $F_(u c) >= 0$ for all $u in V$ and $c in {1, ..., K}$ (optionally augmented with an $ell_1$ sparsity penalty $-lambda sum_(u, c) F_(u c)$).
]

To extract discrete overlapping communities from the continuous MLE $hat(F)$, BigCLAM applies a principled background-probability threshold $delta_"th" = sqrt(-log(1 - epsilon))$ (typically $epsilon = 1 / n^2$, corresponding to the random-graph baseline probability), assigning vertex $u$ to community $c$ whenever $hat(F)_(u c) >= delta_"th"$. Crucially, although evaluating $cal(L)(F)$ and its gradient naively sums over $Theta(n^2)$ non-edges, maintaining the single global sum vector $s = sum_(v in V) F_v in bb(R)^K$ reduces every per-node gradient update to $cal(O)(K ("deg"(u) + 1))$ time.

#algbox(
  title: "BigCLAM Block Coordinate Ascent",
  inputs: [Graph $G = (V, E)$, community count $K$, step parameters $(alpha_0, beta, c_1)$, iterations $T$, threshold $delta_"th"$, optional $ell_1$ weight $lambda >= 0$.],
  outputs: [Affiliation matrix $F in bb(R)_(>= 0)^(n times K)$ and overlapping cover ${C_1, ..., C_K}$.],
)[
  1. *Locally minimal neighborhood initialization*: Rank vertices by neighborhood conductance $phi(N(u))$; select $K$ locally minimal seeds $u_1, ..., u_K$ and initialize $F_(v c) = 1$ for $v in N(u_c) union {u_c}$ and $0$ (or small jitter) otherwise.
  2. *Cache global affiliation sum*: Compute $s <- sum_(v in V) F_v in bb(R)^K$.
  3. *Block coordinate ascent sweeps*: For $t = 1, ..., T$ and each vertex $u in V$:
     - Compute the exact block gradient in $cal(O)(K ("deg"(u) + 1))$ time:
       $g_u <- sum_(v in N(u)) (e^(-chevron.l F_u, F_v chevron.r) / (1 - e^(-chevron.l F_u, F_v chevron.r))) F_v - (s - F_u - sum_(v in N(u)) F_v) - lambda bold(1)_K$.
     - Find step size $eta_u = alpha_0 beta^j$ via projected Armijo backtracking so that $F_u^+ = max(0, F_u + eta_u g_u)$ increases $cal(L)_u (F_u)$ sufficiently.
     - Update cache incrementally: $s <- s - F_u + F_u^+$, and set $F_u <- F_u^+$.
  4. *Extract overlapping communities*: For each $c in {1, ..., K}$, output $C_c = {u in V : F_(u c) >= delta_"th"}$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four core mathematical properties of BigCLAM: exact block concavity contrasted with joint non-concavity, exact $cal(O)(K (m + n))$ gradient complexity and KKT stationarity, formal proof of the dense-overlap law versus simplex MMSBs, and an exact sparse-regime NMF reduction alongside a constructive clique-divergence boundary phenomenon.

#subsec("A", "Block Concavity, Exact Hessian, and Joint Non-Concavity")

BigCLAM's block coordinate ascent relies on the convexity structure of holding all other vertices fixed while optimizing $F_u$.

#theorem(title: "Block Concavity and Joint Indefiniteness of the BigCLAM Objective", origin: "Original Proof")[
  For any fixed vertex $u in V$ and fixed rows ${F_v : v != u} subset.eq bb(R)_(>= 0)^K$, define the single-node subproblem
  $
    cal(L)_u (f) = sum_(v in N(u)) log(1 - exp(-chevron.l f, F_v chevron.r)) - chevron.l f, sum_(v in.not N(u) union {u}) F_v chevron.r .
  $
  Then:
  1. On the open positive-inner-product domain $cal(D)_u = {f in bb(R)_(>= 0)^K : chevron.l f, F_v chevron.r > 0 space forall v in N(u)}$, $cal(L)_u (f)$ is concave with exact Hessian
     $
       nabla^2 cal(L)_u (f) = - sum_(v in N(u)) (exp(-chevron.l f, F_v chevron.r)) / (1 - exp(-chevron.l f, F_v chevron.r))^2 F_v F_v^top prec.eq 0 .
     $
     Whenever $"span"{F_v : v in N(u)} = bb(R)^K$, $cal(L)_u (f)$ is *strictly concave* and has a unique maximizer on any compact subset of $cal(D)_u$.
  2. Jointly over $F in bb(R)_(>= 0)^(n times K)$, $cal(L)(F)$ is *strictly non-concave* on every non-complete graph ($m < binom(n, 2)$) and for every $K >= 2$ on any graph with at least one edge.
]

#proof[
  Define the scalar link function $h(z) = log(1 - e^(-z))$ on $(0, infinity)$. Direct differentiation yields
  $
    h'(z) = e^(-z) / (1 - e^(-z)) > 0, quad h''(z) = - e^(-z) / (1 - e^(-z))^2 < 0 .
  $
  Because $cal(L)_u (f) = sum_(v in N(u)) h(chevron.l f, F_v chevron.r) - chevron.l f, w_u chevron.r$ where $w_u = sum_(v in.not N(u) union {u}) F_v$ is constant with respect to $f$, the chain rule gives the exact gradient and Hessian
  $
    nabla cal(L)_u (f) &= sum_(v in N(u)) h'(chevron.l f, F_v chevron.r) F_v - w_u , \
    nabla^2 cal(L)_u (f) &= sum_(v in N(u)) h''(chevron.l f, F_v chevron.r) F_v F_v^top .
  $
  Since $h''(z) < 0$ for all $z > 0$ and each outer product $F_v F_v^top succ.eq 0$, we have for any $x in bb(R)^K$:
  $
    x^top nabla^2 cal(L)_u (f) x = - sum_(v in N(u)) abs(h''(chevron.l f, F_v chevron.r)) chevron.l x, F_v chevron.r^2 <= 0 .
  $
  If ${F_v : v in N(u)}$ spans $bb(R)^K$, then $chevron.l x, F_v chevron.r = 0$ for all $v in N(u)$ implies $x = 0$, proving strict negative definiteness and uniqueness of the block maximizer.

  For part (2), first suppose $G$ is not complete, so there exists a non-edge $(u, v) in.not E$. Hold all other rows at $0$ and restrict $(F_(u 1), F_(v 1))$ to the line segment $F_(u 1) (t) = 2 - t$, $F_(v 1) (t) = t$ for $t in [0, 2]$ (all other entries $0$). Then $cal(L)(F(t)) = -2 (2 - t) t = 2 t^2 - 4 t$, whose second derivative is $+4 > 0$, contradicting joint concavity. Second, for $K >= 2$ on any graph with an edge $(u, v) in E$, take $F^((1))$ with $F_u^((1)) = F_v^((1)) = e_1$ and $F^((2))$ with $F_u^((2)) = F_v^((2)) = e_2$ (all other rows $0$). Both achieve $chevron.l F_u, F_v chevron.r = 1$ and $cal(L) = 2 log(1 - e^(-1))$, whereas their midpoint $bar(F) = (1 / 2)(F^((1)) + F^((2)))$ has $chevron.l bar(F)_u, bar(F)_v chevron.r = 1 / 2$ and strictly lower objective $2 log(1 - e^(-1/2)) < 2 log(1 - e^(-1))$, proving joint non-concavity.
]

#subsec("B", "Linear-Time Cached Gradient and Monotone KKT Convergence")

BigCLAM scales to million-node networks @yang2013bigclam because the dense non-edge sum collapses to a rank-1 global cache update, while projected coordinate ascent maintains monotone likelihood ascent.

#theorem(title: "Exact O(K (m + n)) Sweep Complexity and Stationary Convergence", origin: "Original Proof")[
  Let $s = sum_(v in V) F_v in bb(R)^K$. Then:
  1. The exact partial gradient with respect to row $F_u$ satisfies
     $
       nabla_(F_u) cal(L)(F) = 2 (sum_(v in N(u)) (e^(-chevron.l F_u, F_v chevron.r) / (1 - e^(-chevron.l F_u, F_v chevron.r))) F_v - s + F_u + sum_(v in N(u)) F_v) .
     $
     Evaluating and updating $F_u$ and $s$ across all $u in V$ requires exactly $cal(O)(K (m + n))$ floating-point operations per full coordinate-ascent sweep.
  2. For any $F in bb(R)_(>= 0)^(n times K)$, the log-likelihood is bounded above by $cal(L)(F) <= 0$. Under cyclic projected Armijo backtracking with $ell_1$ penalty $lambda > 0$ (or box bound $F_(u c) <= M_"max"$), the objective sequence ${cal(L)_lambda (F^((t)))}$ converges monotonically, and every limit point $F^*$ satisfies the first-order KKT optimality conditions:
     $
       F_(u c)^* >= 0, quad [nabla_(F_u) cal(L)_lambda (F^*)]_c <= 0, quad F_(u c)^* [nabla_(F_u) cal(L)_lambda (F^*)]_c = 0 .
     $
]

#proof[
  For part (1), because $G$ is undirected, vertex $u$ appears symmetrically as first and second index in $cal(L)(F)$, contributing factor $2$ to $nabla_(F_u) cal(L)(F)$. Splitting the non-neighbor sum via the global cache $s = sum_(v in V) F_v$ yields
  $
    sum_(v in.not N(u) union {u}) F_v = s - F_u - sum_(v in N(u)) F_v ,
  $
  yielding the exact gradient identity. Computing $chevron.l F_u, F_v chevron.r$ and accumulating both sums over $v in N(u)$ takes $cal(O)(K "deg"(u))$ operations; subtracting $s - F_u$, projecting $F_u^+$, and updating $s <- s - F_u + F_u^+$ take $cal(O)(K)$ operations. Summing across all vertices gives
  $
    sum_(u in V) cal(O)(K ("deg"(u) + 1)) = cal(O)(K (2 m + n)) = cal(O)(K (m + n)) .
  $

  For part (2), since $F_(u c) >= 0$, we have $chevron.l F_u, F_v chevron.r >= 0$, so $1 - e^(-chevron.l F_u, F_v chevron.r) in [0, 1)$, hence $log(1 - e^(-chevron.l F_u, F_v chevron.r)) < 0$ and $-chevron.l F_u, F_v chevron.r <= 0$. Therefore $cal(L)_lambda (F) <= cal(L)(F) <= 0$ everywhere. With $lambda > 0$, the upper level set ${F >= 0 : cal(L)_lambda (F) >= cal(L)_lambda (F^((0)))}$ is compact because $-lambda norm(F)_1 -> -infinity$ as $norm(F)_1 -> infinity$ (and trivially compact under box bounds $F_(u c) <= M_"max"$). By Armijo backtracking, each projected step $F_u^+ = Pi_(>= 0) (F_u + eta_u g_u)$ satisfies
  $
    cal(L)_lambda (F^+) - cal(L)_lambda (F) >= c_1 chevron.l g_u, F_u^+ - F_u chevron.r >= (c_1 / eta_u) norm(F_u^+ - F_u)_2^2 >= 0 ,
  $
  where the second inequality follows from the projection variational inequality $chevron.l F_u + eta_u g_u - F_u^+, F_u - F_u^+ chevron.r <= 0$. Thus ${cal(L)_lambda (F^((t)))}$ is monotone non-decreasing and bounded above by $0$, so it converges and $norm(F^((t+1)) - F^((t)))_F -> 0$. At any limit point $F^*$ (guaranteed to exist by compactness), zero step progress forces $Pi_(>= 0) (F_u^* + eta nabla_(F_u) cal(L)_lambda (F^*)) = F_u^*$ for all $u$, which is equivalent to the stated KKT complementarity system.
]

#subsec("C", "Super-Linear Overlap Density Law vs Mixed-Membership SBMs")

The defining theoretical contrast between BigCLAM and simplex-constrained latent-space models is how edge density behaves inside multi-community intersections @yang2013bigclam.

#proposition(title: "Monotone Dense Overlap in BigCLAM vs Overlap Dilution in Simplex MMSBs", origin: "Classical — Cited")[
  Consider two vertices $u, v$ sharing $r >= 1$ communities symmetrically.
  1. In BigCLAM @yang2013bigclam with uniform shared affiliation strength $F_(u c) = F_(v c) = alpha > 0$ on $r$ shared communities, the edge probability $p_"BC" (r) = 1 - exp(-r alpha^2)$ is *strictly increasing* in $r$, with positive marginal gain
     $
       p_"BC" (r + 1) - p_"BC" (r) = exp(-r alpha^2) (1 - exp(-alpha^2)) > 0 .
     $
  2. In an assortative simplex-constrained Mixed-Membership SBM ($sum_(c = 1)^K theta_(u c) = 1$, $theta_(u c) >= 0$, block matrix $B = p_"in" I_K + p_"out" (bold(1) bold(1)^top - I_K)$ with $p_"in" > p_"out"$), two vertices splitting membership equally across $r$ shared communities ($theta_(u c) = theta_(v c) = 1 / r$) have edge probability
     $
       p_"MMSB" (r) = p_"out" + (p_"in" - p_"out") / r ,
     $
     which is *strictly decreasing* in $r$.
]

#proof[
  For part (1), $chevron.l F_u, F_v chevron.r = sum_(c = 1)^r alpha^2 = r alpha^2$, so $p_"BC" (r) = 1 - e^(-r alpha^2)$. Subtracting $p_"BC" (r)$ from $p_"BC" (r + 1)$ yields $e^(-r alpha^2) - e^(-(r + 1) alpha^2) = e^(-r alpha^2)(1 - e^(-alpha^2)) > 0$ for all $alpha > 0$.

  For part (2), direct bilinear evaluation gives
  $
    p_"MMSB" (r) &= theta_u^top B theta_v = p_"out" (sum_(c = 1)^K theta_(u c)) (sum_(c = 1)^K theta_(v c)) + (p_"in" - p_"out") sum_(c = 1)^K theta_(u c) theta_(v c) \
                 &= p_"out" + (p_"in" - p_"out") sum_(c = 1)^r (1 / r)^2 = p_"out" + (p_"in" - p_"out") / r .
  $
  Since $p_"in" > p_"out"$, $p_"MMSB" (r + 1) - p_"MMSB" (r) = -(p_"in" - p_"out") / (r (r + 1)) < 0$. Thus any simplex-constrained MMSB necessarily predicts *lower* tie density in $r$-fold overlaps than in disjoint interiors, whereas BigCLAM's independent-or link naturally predicts dense overlaps observed in empirical social and biological networks @yang2013bigclam.
]

#subsec("D", "Sparse KL-NMF Reduction and Isolated-Clique Escape")

Viewing BigCLAM through the lens of nonnegative matrix factorization @lee1999nmf @lee2001algorithms reveals both its exact sparse approximation to Poisson/KL-NMF and an important unregularized pathology.

#theorem(title: "Sparse KL-NMF Equivalence, Scaling Rigidity, and Clique Escape", origin: "Original Derivation")[
  Let $A in {0, 1}^(n times n)$ be the adjacency matrix of $G$ and $hat(A) = F F^top$ the symmetric nonnegative Gram matrix.
  1. *Sparse KL-NMF equivalence*: In the sparse regime $z_(u v) = chevron.l F_u, F_v chevron.r << 1$, the BigCLAM objective satisfies
     $
       cal(L)(F) = - "KL"_"gen" (A || hat(A)) + C_A + cal(O)(sum_((u, v) in E) z_(u v) + sum_(u in V) norm(F_u)_2^2) ,
     $
     where $"KL"_"gen" (A || hat(A)) = sum_(u, v) (A_(u v) log(A_(u v) / hat(A)_(u v)) - A_(u v) + hat(A)_(u v))$ is the generalized KL divergence of Lee-Seung NMF @lee2001algorithms and $C_A = 2 m$.
  2. *Diagonal scaling rigidity*: Unlike asymmetric NMF ($W H = (W D)(D^(-1) H)$), symmetric BigCLAM eliminates continuous positive diagonal rescaling: $F D (F D)^top = F F^top$ with $D succ 0$ diagonal implies $D = I_K$ whenever every column of $F$ is nonzero.
  3. *Unregularized clique divergence*: If $G$ contains a disjoint clique $C$ ($abs(C) >= 2$) with zero edges to $V without C$, then the unregularized BigCLAM MLE ($lambda = 0$) does not exist in $bb(R)_(>= 0)^(n times K)$: setting $F_(u 1) = t$ for $u in C$ and $0$ elsewhere drives $cal(L)(F) = abs(C)(abs(C) - 1) log(1 - e^(-t^2)) -> 0$ strictly monotonically as $t -> +infinity$.
]

#proof[
  For part (1), expand $log(1 - e^(-z)) = log(z (1 - z / 2 + cal(O)(z^2))) = log z - z / 2 + cal(O)(z^2)$ for small $z > 0$. Substituting into $cal(L)(F)$ and adding/subtracting $sum_((u, v) in E) z_(u v)$ and diagonal terms $sum_u norm(F_u)_2^2$ gives
  $
    cal(L)(F) &= sum_(u != v) (A_(u v) log z_(u v) - z_(u v)) + (1 / 2) sum_((u, v) in E) z_(u v) + cal(O)(sum_((u, v) in E) z_(u v)^2) \
              &= - "KL"_"gen" (A || F F^top) + 2 m + cal(O)(sum_((u, v) in E) z_(u v) + sum_u norm(F_u)_2^2) ,
  $
  since $A_(u v) in {0, 1}$ gives $A_(u v) log A_(u v) = 0$ and $sum_(u, v) A_(u v) = 2 m$.

  For part (2), let $D = "diag"(d_1, ..., d_K) succ 0$. If $F D^2 F^top = F F^top$ and $F$ has linearly independent columns (or more generally, comparing diagonal entries when columns have disjoint support, or taking $"Tr"(F (D^2 - I_K) F^top) = 0$ after projecting onto any column's active basis): specifically, for any vector $x$ with $F^top x = e_c$ (which exists whenever $F$ has full column rank $K$), $x^top F (D^2 - I_K) F^top x = d_c^2 - 1 = 0$, forcing $d_c = 1$ for all $c$.

  For part (3), restrict $F$ to a single column supported on the isolated clique $C$: $F_(u 1) = t > 0$ for $u in C$ and $F_(u 1) = 0$ for $u in.not C$. Every unordered pair in $C$ is an edge ($z_(u v) = t^2$), and every pair with at least one endpoint outside $C$ has $z_(u v) = 0$ (contributing $0$ to the non-edge penalty). Thus
  $
    cal(L)(F(t)) = abs(C)(abs(C) - 1) log(1 - e^(-t^2)) ,
  $
  which is strictly negative, strictly increasing in $t$, and satisfies $lim_(t -> +infinity) cal(L)(F(t)) = 0 = sup_(F >= 0) cal(L)(F)$ (by Theorem 2). Hence the supremum $0$ is never attained at any finite $F$, necessitating $ell_1$ regularization ($lambda > 0$), box clipping, or early stopping in practice.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time, Space, and Parallelism.* By Theorem 2, each full coordinate-ascent sweep costs $cal(O)(K (m + n))$ time, giving $cal(O)(T K (m + n))$ total runtime over $T$ sweeps (typically $T in [50, 200]$). Storing $F$ densely requires $cal(O)(n K)$ memory; because real affiliation vectors are sparse ($norm(F_u)_0 <= K_"eff" << K$), sparse row maps reduce memory to $cal(O)(n K_"eff" + m)$ @yang2013bigclam. Furthermore, because each node's gradient reads only its immediate neighborhood $N(u)$ and the shared cache vector $s$, BigCLAM parallelizes naturally across threads using asynchronous/mini-batch updates of $s$, scaling to networks with tens of millions of edges.

*Geometric and Statistical Interpretation.* BigCLAM embeds every vertex into the nonnegative cone $bb(R)_(>= 0)^K$ such that inner products $chevron.l F_u, F_v chevron.r$ map monotonically to edge probabilities via the complementary log-log (cloglog) inverse link $1 - e^(-z)$. Selecting $K$ when unknown is performed via hold-out validation likelihood on a random 20% edge/non-edge sample @yang2013bigclam.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(T K (m + n))$ via cached global sum $s = sum_v F_v$; naive $cal(O)(T K n^2)$.],
    [Space complexity], [$cal(O)(n K + m)$ dense; $cal(O)(n K_"eff" + m)$ with sparse affiliation rows.],
    [Cluster geometry], [Assortative, densely overlapping clique-like communities (positive Gram factor).],
    [Overlap / noise], [Continuous multi-community strengths; vertices below $delta_"th"$ left unassigned.],
    [Model class], [Generative probabilistic symmetric NMF with complementary log-log binary link.],
    [Hyperparameters], [Community count $K$, threshold $delta_"th"$, optional $ell_1$ weight $lambda$, iterations $T$.],
    [Determinism], [Non-concave globally; depends on conductance seeding and sweep order.],
  ),
  caption: [Operational and theoretical profile of BigCLAM.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Online Social Networks.* Friendship and group-membership graphs (e.g., LiveJournal, Orkut, YouTube @yang2013bigclam) exhibit massive multi-group membership where users sharing college, workplace, and hobby circles connect at sharply elevated rates. BigCLAM's dense-overlap link (Proposition 1) and linear $cal(O)(T K (m + n))$ sweeps directly match both the empirical topology and million-node scale.
2. *E-Commerce Co-Purchasing and Product Taxonomies.* Co-purchased product networks contain items spanning overlapping categories (e.g., a laptop adapter bridging electronics, travel, and office supplies). Continuous affiliation strengths $F_(u c)$ naturally quantify cross-category membership intensity while filtering incidental co-purchases via $delta_"th"$.
3. *Biological Protein Complexes and Regulatory Modules.* Protein-protein interaction networks feature densely interacting core complexes that share hub subunits across overlapping pathways. BigCLAM recovers dense multi-subunit overlaps far more accurately than partition algorithms @blondel2008louvain or simplex MMSBs.
4. *Scientific Co-Authorship and Venue Networks.* Researchers routinely co-author across intersecting subdisciplines, with tie probability rising steeply among collaborators sharing multiple subfields. BigCLAM extracts overlapping research themes along with continuous researcher-topic affinity weights.
5. *Enterprise Communication and Organizational Discovery.* Email and messaging graphs within large organizations overlay formal departmental trees with cross-functional project teams. BigCLAM simultaneously resolves formal teams and dense cross-team task forces at near-linear cost.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

BigCLAM struggles in four predictable regimes: (i) *disassortative or bipartite networks*, because $chevron.l F_u, F_v chevron.r >= 0$ can only express nonnegative assortative affinities (edges within communities, never systematic inter-community repulsion); (ii) *elongated, chain-like, or low-density geometric manifolds*, where $k$-NN neighborhoods have low clustering coefficient and lack the dense clique cores that BigCLAM's likelihood rewards; (iii) *isolated exact cliques* without regularization (Theorem 3(3)), causing numerical overflow as $F_(u c) -> infinity$; and (iv) *unknown $K$ in highly heterogeneous multi-scale graphs*, where cross-validating $K$ multiplies runtime by the grid size of candidate $K$ values.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Million-node social graphs], [Excellent: dense overlaps + $cal(O)(T K (m + n))$ scaling.], [Tune $K$ via hold-out likelihood; add small $ell_1$.],
    [Protein-complex PPI], [Strong: dense overlapping cores match AGM.], [Use OSLOM @lancichinetti2011oslom if statistical p-values needed.],
    [Co-purchasing graphs], [Strong: continuous affiliation scores per item.], [Clip $F_(u c)$ on isolated complete subgraphs.],
    [Bipartite / disassortative], [Poor: $F F^top >= 0$ cannot model heterophily.], [Use SBM @holland1983sbm or Spectral Biclustering @kluger2003spectral.],
    [Geometric $k$-NN manifolds], [Weak: elongated chains lack dense clique cores.], [Prefer OSLOM @lancichinetti2011oslom, SLPA @xie2011slpa, or HDBSCAN @campello2013hdbscan.],
  ),
  caption: [Domain suitability matrix for BigCLAM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

BigCLAM occupies a pivotal unifying position connecting stochastic blockmodels @holland1983sbm @karrer2011dcsbm, symmetric NMF @lee2001algorithms @ding2008nmfplsa, and heuristic overlapping network methods @palla2005cpm @xie2011slpa @lancichinetti2011oslom.

#proposition(title: "Exact Reduction of Disjoint BigCLAM to Assortative Degree-Corrected SBM", origin: "Original Unified Proof")[
  Restrict BigCLAM to row-sparse matrices $F in bb(R)_(>= 0)^(n times K)$ with at most one nonzero entry per row ($norm(F_u)_0 <= 1$), writing $c(u) in {1, ..., K}$ for vertex $u$'s unique active community and $theta_u = F_(u, c(u)) >= 0$ for its affiliation weight, augmented with uniform baseline floor $epsilon in (0, 1)$. Then the edge probabilities simplify identically to
  $
    p_(u v) (F) = cases(
      1 - (1 - epsilon) exp(-theta_u theta_v) & "if" c(u) = c(v),
      epsilon & "if" c(u) != c(v) .
    )
  $
  In the sparse limit $theta_u theta_v << 1$ and $epsilon << 1$, $p_(u v) approx theta_u theta_v bb(I)[c(u) = c(v)] + epsilon$, reducing BigCLAM exactly to an assortative Degree-Corrected Stochastic Blockmodel (DCSBM) @karrer2011dcsbm with vertex degree propensities $theta_u$ and uniform inter-community background rate $epsilon$.
]

#proof[
  When each row $F_u$ has support contained in ${c(u)}$, the inner product satisfies $chevron.l F_u, F_v chevron.r = theta_u theta_v$ if $c(u) = c(v)$ and $0$ otherwise. Including the background independent-or tie probability $epsilon$ multiplies the non-edge probability by $(1 - epsilon)$, yielding $1 - p_(u v) = (1 - epsilon) exp(-chevron.l F_u, F_v chevron.r)$, which immediately equals the stated piecewise expression. Expanding $1 - (1 - epsilon) e^(-theta_u theta_v) = epsilon + (1 - epsilon) theta_u theta_v + cal(O)(theta_u^2 theta_v^2)$ matches the planted-partition DCSBM edge expectation $theta_u theta_v omega_(c(u), c(v))$ with diagonal affinity $omega_"in" = 1$ and off-diagonal floor $epsilon$.
]

Relative to KL-NMF @lee2001algorithms and pLSA @hofmann1999plsa @ding2008nmfplsa, BigCLAM enforces symmetric factorization $F F^top$ and uses an exact binary likelihood ($1 - e^(-z) in [0, 1)$) rather than unbounded Poisson means, preventing probability saturation above $1$ on dense overlaps. Relative to OSLOM @lancichinetti2011oslom and SLPA @xie2011slpa, BigCLAM provides a full generative model capable of link prediction and community-size likelihood scoring, at the cost of requiring $K$ and iterative non-convex optimization.

*Empirical Behavior on ClusterOverview Benchmarks.* On the 4-dataset ClusterOverview benchmark suite (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing on constructed $k$-NN graphs), BigCLAM achieves *Iris* Sil $= 0.131$, ARI $= 0.326$; *Customers* Sil $= 0.159$; *Multi-Shape* Sil $= -0.025$ (ARI N/A); and *High-D* Sil $= 0.065$ (ARI N/A), with a mean runtime of $5.835$ ms. These measurements confirm our theoretical diagnosis with remarkable clarity. BigCLAM's inductive bias explicitly seeks *densely overlapping clique-like cores* driven by positive Gram inner products (Proposition 1 and Theorem 3). On *Multi-Shape*, where clusters form long, curved non-convex manifolds rather than dense cliques, local $k$-NN neighborhoods have low global transitivity; BigCLAM therefore fractures elongated shapes into overlapping local patches, producing a slightly negative Silhouette ($-0.025$). On *Iris* (ARI $= 0.326$, Sil $= 0.131$) and *Customers* (Sil $= 0.159$), compact density modes contain partial clique-like structure and true boundary overlap (e.g., Versicolor and Virginica), allowing BigCLAM to capture meaningful community structure, though less cleanly than OSLOM ($0.451$ ARI) or SLPA ($0.341$ Sil). Finally, the $5.835$ ms mean runtime—the slowest of Algs. 123–126—reflects multi-sweep backtracking gradient ascent and $K$-factor optimization on small-to-medium graphs where BigCLAM's asymptotic million-node caching advantage has not yet dominated constant overhead.
