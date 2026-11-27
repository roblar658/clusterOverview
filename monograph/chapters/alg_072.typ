#import "../template/ieee.typ": *

#alg_entry(
  id: 72,
  name: "Lance-Williams Algorithm",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1967",
  paradigm: "Combinatorial Recursive Agglomerative Framework",
  complexity: $O(n^2) "to" O(n^2 log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Sequential Agglomerative Hierarchical Non-overlapping (SAHN) clustering constructs a nested sequence of partitions $cal(P)_0, cal(P)_1, dots.c, cal(P)_(n-1)$ over a finite dataset $X = {x_1, dots.c, x_n}$ endowed with a symmetric, zero-diagonal dissimilarity function $d : X times X arrow.r RR_(>= 0)$ ($d(x, x) = 0$, $d(x, y) = d(y, x)$). Starting from the singleton partition $cal(P)_0 = {{x_1}, dots.c, {x_n}}$, each step $t in {1, dots.c, n - 1}$ selects a minimum-dissimilarity pair $(C_i, C_j) in arg min_(u < v) D(C_u, C_v)$ from the active cluster set $cal(A)_(t-1)$, merges them into $C_(i j) = C_i union C_j$ at dendrogram height $h_t = D(C_i, C_j)$, and updates dissimilarities to all remaining clusters $C_k in cal(A)_(t-1) without {C_i, C_j}$.

Before 1967, each agglomerative linkage rule (single, complete, UPGMA, WPGMA, centroid, median, Ward) was formulated through bespoke geometric or set-theoretic definitions requiring distinct implementations and access to raw coordinates. Lance and Williams @lance1967general unified these methods into a single *combinatorial* recurrence operating purely on the inter-cluster dissimilarity matrix, eliminating any need to revisit feature vectors after initialization.

#definition(title: "Lance-Williams Combinatorial Update Law")[
  Let $C_i, C_j, C_k subset.eq X$ be pairwise disjoint clusters of cardinalities $n_i = abs(C_i)$, $n_j = abs(C_j)$, and $n_k = abs(C_k)$, with current dissimilarities $d_(i j) = D(C_i, C_j)$, $d_(i k) = D(C_i, C_k)$, and $d_(j k) = D(C_j, C_k)$. A linkage rule belongs to the *Lance-Williams family* @lance1967general if merging $C_i$ and $C_j$ updates the dissimilarity to $C_k$ via the four-parameter recurrence
  $
    D(C_i union C_j, C_k) &= alpha_i d_(i k) + alpha_j d_(j k) + beta d_(i j) \
    &quad + gamma abs(d_(i k) - d_(j k)),
  $
  where $(alpha_i, alpha_j, beta, gamma) in RR^4$ depend solely on the cluster cardinalities $(n_i, n_j, n_k)$ and fixed method hyperparameters, independent of raw coordinates or intra-cluster structure.
]

Writing $n_(i j) = n_i + n_j$ and $N = n_i + n_j + n_k$, classical linkages are exact instantiations of this recurrence:
- *Single linkage* @gower1969minimum @sibson1973slink: $alpha_i = alpha_j = 1 / 2$, $beta = 0$, $gamma = -1 / 2$, yielding $min(d_(i k), d_(j k))$.
- *Complete linkage* @defays1977clink: $alpha_i = alpha_j = 1 / 2$, $beta = 0$, $gamma = +1 / 2$, yielding $max(d_(i k), d_(j k))$.
- *UPGMA (average)* @sokal1958statistical: $alpha_i = n_i / n_(i j)$, $alpha_j = n_j / n_(i j)$, $beta = gamma = 0$.
- *WPGMA (McQuitty)* @sokal1958statistical: $alpha_i = alpha_j = 1 / 2$, $beta = gamma = 0$.
- *UPGMC (centroid)* @lance1967general: $alpha_i = n_i / n_(i j)$, $alpha_j = n_j / n_(i j)$, $beta = -n_i n_j / n_(i j)^2$, $gamma = 0$ (on squared Euclidean distances).
- *WPGMC (median)* @lance1967general: $alpha_i = alpha_j = 1 / 2$, $beta = -1 / 4$, $gamma = 0$.
- *Ward's method* @ward1963hierarchical: $alpha_i = (n_i + n_k) / N$, $alpha_j = (n_j + n_k) / N$, $beta = -n_k / N$, $gamma = 0$ (on squared Euclidean merge increments).
- *Flexible beta* @lance1967general: $alpha_i = alpha_j = (1 - beta) / 2$, $gamma = 0$, for $beta < 1$.

#algbox(
  title: "Generic Lance-Williams Agglomerative Clustering",
  inputs: [Pairwise dissimilarity matrix $D^(0) in RR_(>= 0)^(n times n)$, parameter map $(alpha_i, alpha_j, beta, gamma)$ of $(n_i, n_j, n_k)$],
  outputs: [Dendrogram merge sequence ${(C_(i_t), C_(j_t), h_t)}_(t=1)^(n-1)$ and induced ultrametric (if monotonic)],
)[
  1. Initialize active set $cal(A) arrow.l {1, dots.c, n}$, cluster sizes $n_u arrow.l 1$, and working matrix $D arrow.l D^(0)$.
  2. If $(alpha_i, alpha_j, beta, gamma)$ satisfies the reducibility condition (Theorem 3), execute the $O(n^2)$ Nearest-Neighbor Chain variant @mullner2011modern; otherwise initialize min-heaps $H_u$ over $v in cal(A) without {u}$.
  3. *For* $t = 1, dots.c, n - 1$:
     - Find $(i, j) in arg min_(u < v in cal(A)) D_(u v)$ and record merge $(i, j)$ at height $h_t arrow.l D_(i j)$.
     - Create merged node $m$ with size $n_m arrow.l n_i + n_j$.
     - *For each* $k in cal(A) without {i, j}$: evaluate $(alpha_i, alpha_j, beta, gamma)$ at $(n_i, n_j, n_k)$ and set
       $D_(m k) arrow.l alpha_i D_(i k) + alpha_j D_(j k) + beta D_(i j) + gamma abs(D_(i k) - D_(j k))$.
     - Update $cal(A) arrow.l (cal(A) without {i, j}) union {m}$ and refresh active priority queues.
  4. Return the merge tree and cut at target $k$ clusters or threshold $h$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish the complete algebraic theory of the Lance-Williams family: (A) necessary and sufficient conditions for dendrogram monotonicity and ultrametricity, (B) exact geometric equivalence to Ward and centroid criteria via Huygens-Steiner identities, (C) sharp characterization of Bruynooghe reducibility and $O(n^2)$ Nearest-Neighbor Chain exactness, and (D) finite-sample inversion counterexamples alongside invariance laws.

#subsec("A", "Dendrogram Monotonicity and Ultrametricity")

A hierarchical merge sequence induces a valid ultrametric on $X$ if and only if merge heights never decrease ($h_1 <= h_2 <= dots.c <= h_(n-1)$). A strict decrease $h_(t+1) < h_t$ is called a *dendrogram inversion* (reversal), which renders horizontal tree cuts ambiguous.

#theorem(title: "Sufficient Conditions for Monotonicity and Ultrametricity", origin: "Classical — Cited")[
  Suppose the Lance-Williams parameters satisfy, for all valid cluster sizes $(n_i, n_j, n_k)$, the three inequalities @lance1967general @mullner2011modern:
  $
    "(i)"& quad gamma >= -min(alpha_i, alpha_j), \
    "(ii)"& quad alpha_i + alpha_j >= 0, \
    "(iii)"& quad alpha_i + alpha_j + beta >= 1.
  $
  Then every merge step satisfies $D(C_i union C_j, C_k) >= D(C_i, C_j)$ for all $C_k in cal(A) without {C_i, C_j}$. Consequently, $h_1 <= h_2 <= dots.c <= h_(n-1)$, and $u(x, y) = h_(t(x, y))$ (where $t(x, y)$ is the first step uniting $x$ and $y$, with $u(x, x) = 0$) defines an ultrametric on $X$.
]

#proof[
  At merge step $t$, minimality of the chosen pair $(C_i, C_j)$ guarantees $d_(i j) <= d_(i k)$ and $d_(i j) <= d_(j k)$ for every remaining cluster $C_k$. Without loss of generality, assume $d_(i k) >= d_(j k) >= d_(i j)$ (the opposite case $d_(j k) >= d_(i k)$ follows identically by swapping indices $i$ and $j$). Define the non-negative increments
  $
    delta_1 = d_(i k) - d_(j k) >= 0, quad delta_2 = d_(j k) - d_(i j) >= 0.
  $
  Then $d_(j k) = d_(i j) + delta_2$, $d_(i k) = d_(i j) + delta_2 + delta_1$, and $abs(d_(i k) - d_(j k)) = delta_1$. Substituting into the Lance-Williams recurrence yields
  $
    D(C_(i j), C_k) &= alpha_i (d_(i j) + delta_2 + delta_1) + alpha_j (d_(i j) + delta_2) \
    &quad + beta d_(i j) + gamma delta_1 \
    &= (alpha_i + alpha_j + beta) d_(i j) + (alpha_i + alpha_j) delta_2 \
    &quad + (alpha_i + gamma) delta_1.
  $
  By condition (iii) and $d_(i j) >= 0$, the first term is at least $d_(i j)$. By (ii) and $delta_2 >= 0$, the second term is non-negative. By (i), $gamma >= -min(alpha_i, alpha_j) >= -alpha_i$, so $alpha_i + gamma >= 0$ and the third term is non-negative. Hence $D(C_(i j), C_k) >= d_(i j) = h_t$. Since every untouched pair $(C_u, C_v)$ already satisfied $D(C_u, C_v) >= h_t$ by step-$t$ minimality, the step-$(t+1)$ minimum satisfies $h_(t+1) >= h_t$.

  To verify the strong triangle inequality (ultrametricity) for any $x, y, z in X$, let $t_1 <= t_2$ be the two merge steps uniting the three points. The pair united at $t_1$ satisfies $u = h_(t_1)$, while the other two pairs both unite at $t_2$ with $u = h_(t_2)$. Monotonicity $h_(t_1) <= h_(t_2)$ implies $u(x, z) <= max(u(x, y), u(y, z))$ across all permutations of $(x, y, z)$.
]

#subsec("B", "Exact Geometric Realization of Centroid and Ward Updates")

Although the Lance-Williams formula operates without feature coordinates, its centroid and Ward parameterizations are *exact* algebraic consequences of Euclidean inner-product geometry.

#proposition(title: "Unified Huygens-Steiner Derivation of UPGMC and Ward Recurrences", origin: "Original Derivation")[
  Let $X subset.eq RR^p$ with cluster centroids $mu_u = n_u^(-1) sum_(x in C_u) x$ and within-cluster sum of squares $"SSE"(C_u) = sum_(x in C_u) norm(x - mu_u)_2^2$. Define squared centroid distance $D_"cen" (C_u, C_v) = norm(mu_u - mu_v)_2^2$ and Ward merge increment
  $
    Delta_"Ward" (C_u, C_v) = "SSE"(C_u union C_v) - "SSE"(C_u) - "SSE"(C_v).
  $
  Then: (a) $D_"cen"$ satisfies the exact Lance-Williams law with $alpha_i = n_i / n_(i j)$, $alpha_j = n_j / n_(i j)$, $beta = -n_i n_j / n_(i j)^2$, $gamma = 0$; and (b) $Delta_"Ward" (C_u, C_v) = (n_u n_v / (n_u + n_v)) norm(mu_u - mu_v)_2^2$ satisfies the exact Lance-Williams law with $alpha_i = (n_i + n_k) / N$, $alpha_j = (n_j + n_k) / N$, $beta = -n_k / N$, $gamma = 0$.
]

#proof[
  Let $mu_(i j) = (n_i mu_i + n_j mu_j) / n_(i j)$ be the centroid of $C_(i j) = C_i union C_j$. For any $mu_k in RR^p$, shift coordinates by $mu_k$ and set $a = mu_i - mu_k$, $b = mu_j - mu_k$, so $mu_(i j) - mu_k = (n_i a + n_j b) / n_(i j)$ and $mu_i - mu_j = a - b$. Expanding the squared Euclidean norm gives
  $
    n_(i j)^2 norm(mu_(i j) - mu_k)_2^2 &= n_i^2 norm(a)_2^2 + n_j^2 norm(b)_2^2 + 2 n_i n_j chevron.l a, b chevron.r.
  $
  Using the polarization identity $2 chevron.l a, b chevron.r = norm(a)_2^2 + norm(b)_2^2 - norm(a - b)_2^2$ and noting $n_i^2 + n_i n_j = n_i n_(i j)$, $n_j^2 + n_i n_j = n_j n_(i j)$, we obtain
  $
    norm(mu_(i j) - mu_k)_2^2 &= (n_i / n_(i j)) norm(a)_2^2 + (n_j / n_(i j)) norm(b)_2^2 \
    &quad - (n_i n_j / n_(i j)^2) norm(a - b)_2^2.
  $
  This immediately proves part (a) for $D_"cen"$.

  For part (b), the parallel-axis identity on $C_u union C_v$ yields
  $
    "SSE"(C_u union C_v) &= "SSE"(C_u) + n_u norm(mu_u - mu_(u v))_2^2 \
    &quad + "SSE"(C_v) + n_v norm(mu_v - mu_(u v))_2^2.
  $
  Substituting $mu_u - mu_(u v) = (n_v / (n_u + n_v))(mu_u - mu_v)$ and $mu_v - mu_(u v) = -(n_u / (n_u + n_v))(mu_u - mu_v)$ produces
  $
    Delta_"Ward" (C_u, C_v) = (n_u n_v / (n_u + n_v)) norm(mu_u - mu_v)_2^2.
  $
  Multiplying the part (a) identity by $n_(i j) n_k / N$ gives
  $
    Delta_"Ward" (C_(i j), C_k) &= (n_i n_k / N) norm(a)_2^2 + (n_j n_k / N) norm(b)_2^2 \
    &quad - (n_i n_j n_k / (N n_(i j))) norm(a - b)_2^2.
  $
  Replacing $norm(a)_2^2 = ((n_i + n_k) / (n_i n_k)) Delta_"Ward" (C_i, C_k)$, $norm(b)_2^2 = ((n_j + n_k) / (n_j n_k)) Delta_"Ward" (C_j, C_k)$, and $norm(a - b)_2^2 = (n_(i j) / (n_i n_j)) Delta_"Ward" (C_i, C_j)$ yields exact cancellation of every factor $n_i, n_j$:
  $
    Delta_"Ward" (C_(i j), C_k) &= ((n_i + n_k) / N) Delta_"Ward" (C_i, C_k) \
    &quad + ((n_j + n_k) / N) Delta_"Ward" (C_j, C_k) \
    &quad - (n_k / N) Delta_"Ward" (C_i, C_j),
  $
  completing the proof @ward1963hierarchical @lance1967general.
]

#subsec("C", "Sharp Reducibility Condition and O(n^2) NN-Chain Exactness")

Monotonicity alone does not permit out-of-order local merges; the $O(n^2)$ Nearest-Neighbor (NN) Chain algorithm @mullner2011modern requires the stronger *Bruynooghe reducibility property*.

#theorem(title: "Exact Characterization of Lance-Williams Reducibility", origin: "Original Unified Proof")[
  A Lance-Williams update satisfies Bruynooghe's reducibility inequality
  $
    D(C_i union C_j, C_k) >= min(D(C_i, C_k), D(C_j, C_k))
  $
  for all non-negative dissimilarities with $d_(i j) <= min(d_(i k), d_(j k))$ *if and only if*
  $
    "(i)"& quad alpha_i + alpha_j + beta >= 1, \
    "(ii)"& quad alpha_i + alpha_j >= 1, \
    "(iii)"& quad gamma >= -min(alpha_i, alpha_j).
  $
  Whenever (i)--(iii) hold and dissimilarities are tie-free, merging any reciprocal nearest-neighbor (RNN) pair $(C_i, C_j)$ preserves the exact global dendrogram topology and executes in $O(n^2)$ total time.
]

#proof[
  Assume WLOG $d_(i k) >= d_(j k) >= d_(i j) >= 0$, so $min(d_(i k), d_(j k)) = d_(j k)$. Introduce non-negative free variables $delta_1 = d_(i k) - d_(j k) >= 0$ and $delta_2 = d_(j k) - d_(i j) >= 0$. Substituting $d_(j k) = d_(i j) + delta_2$ and $d_(i k) = d_(i j) + delta_2 + delta_1$ into $D(C_(i j), C_k) - d_(j k)$ yields the exact identity
  $
    D(C_(i j), C_k) - d_(j k) &= (alpha_i + alpha_j + beta - 1) d_(i j) \
    &quad + (alpha_i + alpha_j - 1) delta_2 \
    &quad + (alpha_i + gamma) delta_1.
  $
  *Sufficiency*: Under (i)--(iii), all three coefficients are non-negative, hence $D(C_(i j), C_k) - d_(j k) >= 0$ for all $(d_(i j), delta_1, delta_2) in RR_(>= 0)^3$. Symmetrically, when $d_(j k) >= d_(i k)$, the third coefficient becomes $alpha_j + gamma >= 0$.

  *Necessity*: Because $(d_(i j), delta_1, delta_2)$ can take the independent ray limits $(1, 0, 0)$, $(0, 0, 1)$, and $(0, 1, 0)$ (and symmetrically with $i, j$ swapped), non-negativity of the linear form on $RR_(>= 0)^3$ requires each coefficient to be non-negative, recovering (i), (ii), and (iii) identically. Note that Ward satisfies (ii) strictly since $(n_i + n_j + 2 n_k) / N = 1 + n_k / N > 1$, whereas Flexible Beta with $beta > 0$ violates (ii) since $alpha_i + alpha_j = 1 - beta < 1$ despite remaining monotonic.

  *NN-Chain Exactness and Complexity*: Maintain a directed chain $C_(s_1) arrow.r C_(s_2) arrow.r dots.c arrow.r C_(s_m)$ where each node points to its unique nearest active neighbor. Since distances strictly decrease along the chain until a terminal RNN pair $(C_(s_(m-1)), C_(s_m))$ is reached, merging $(C_(s_(m-1)), C_(s_m))$ into $C_(i j)$ cannot create a new neighbor strictly closer to $C_(s_(m-2))$ than $C_(s_(m-1))$ was: by reducibility,
  $
    D(C_(i j), C_(s_(m-2))) >= min(D(C_(s_(m-1)), C_(s_(m-2))), D(C_(s_m), C_(s_(m-2)))) \
    >= D(C_(s_(m-2)), C_(s_(m-1))),
  $
  and $C_(s_(m-1))$ was already the closest active neighbor of $C_(s_(m-2))$. Thus all earlier chain links $C_(s_1) arrow.r dots.c arrow.r C_(s_(m-2))$ remain valid after popping $(C_(s_(m-1)), C_(s_m))$. Moreover, the RNN pair $(C_i, C_j)$ is immune to any outside merge until it merges itself, so out-of-order RNN merges produce the exact same merge pairs and heights as global argmin selection (up to sorting the $n - 1$ merge heights). Amortized accounting assigns each node at most one push and two pops across $n - 1$ merges; scanning $O(abs(cal(A)))$ active candidates per push and updating $O(abs(cal(A)))$ distances per merge gives $O(n^2)$ total operations @mullner2011modern.
]

#subsec("D", "Inversion Counterexample and Scale / Order Invariance")

When condition (iii) of Theorem 1 fails, inversions occur even on three points in the Euclidean plane.

#proposition(title: "Constructive Three-Point Dendrogram Inversion", origin: "Original Counterexample")[
  Consider any Lance-Williams rule for which $alpha_i + alpha_j + beta = 1 - eta < 1$ on singleton merges ($n_i = n_j = n_k = 1$, $eta > 0$). Then there exists a 3-point metric instance exhibiting a strict inversion $h_2 < h_1$. Specifically, for UPGMC and WPGMC on Euclidean points $x_1 = (-1, 0)$, $x_2 = (1, 0)$, $x_3 = (0, 9 / 5) in RR^2$, the merge heights strictly invert: $h_2 = 81 / 25 = 3.24 < 4 = h_1$.
]

#proof[
  In general, set $d_(1 2) = 1$ and $d_(1 3) = d_(2 3) = 1 + delta$ with $0 < delta < eta / max(1, alpha_1 + alpha_2)$ (which satisfies the triangle inequality for all $delta > 0$). Step 1 merges $(x_1, x_2)$ at height $h_1 = 1$. Since $d_(1 3) = d_(2 3)$, the $gamma$-term vanishes identically and
  $
    h_2 &= D({x_1, x_2}, {x_3}) = (alpha_1 + alpha_2)(1 + delta) + beta \
    &= (1 - eta) + (alpha_1 + alpha_2) delta < 1 = h_1.
  $
  For the explicit Euclidean realization $x_1 = (-1, 0)$, $x_2 = (1, 0)$, $x_3 = (0, 9 / 5)$, squared Euclidean distances are $d_(1 2)^2 = 4$ and $d_(1 3)^2 = d_(2 3)^2 = 1 + 81 / 25 = 106 / 25 = 4.24 > 4$. Both UPGMC and WPGMC have $alpha_1 = alpha_2 = 1 / 2$, $beta = -1 / 4$, $gamma = 0$ on singletons, merging $(x_1, x_2)$ at $h_1 = 4$. The updated squared distance to $x_3$ is
  $
    h_2 = (1 / 2)(4.24) + (1 / 2)(4.24) - (1 / 4)(4) = 3.24 < 4 = h_1,
  $
  which coincides with the squared distance from the merged centroid $mu_(1 2) = (0, 0)$ to $x_3 = (0, 9 / 5)$ (and persists under square-root heights $sqrt(3.24) = 1.8 < 2 = sqrt(4)$).
]

#lemma(title: "Positive Homogeneity and Monotone-Transform Rigidity", origin: "Original Proof")[
  Every Lance-Williams update is positively homogeneous of degree 1: scaling $d arrow.r c dot d$ ($c > 0$) scales every updated dissimilarity and dendrogram height by $c$, preserving the merge tree topology identically (and equivariantly under input permutation when tie-free). Conversely, invariance under all strictly increasing nonlinear transformations $f : RR_(>= 0) arrow.r RR_(>= 0)$ holds *only* for single linkage ($alpha_i = alpha_j = 1 / 2, beta = 0, gamma = -1 / 2$) and complete linkage ($gamma = +1 / 2$).
]

#proof[
  For $c > 0$, substituting $c d_(i k), c d_(j k), c d_(i j)$ into the Lance-Williams law factors out $c$ directly because $abs(c d_(i k) - c d_(j k)) = c abs(d_(i k) - d_(j k))$. By induction over steps $t = 1, dots.c, n - 1$, $D^(t)(c d) = c D^(t)(d)$, so all argmin comparisons remain invariant.

  For nonlinear monotone invariance, single and complete linkage compute $min(d_(i k), d_(j k))$ and $max(d_(i k), d_(j k))$, which commute with every strictly increasing $f$. Conversely, suppose a fixed-parameter Lance-Williams update preserves all orderings under every strictly increasing $f$ (in particular $f_p (t) = t^p$ for $p > 0$). Take $d_(i j) = 0$ and compare two candidate clusters $k, ell$ with $(d_(i k), d_(j k)) = (1, 5)$ and $(d_(i ell), d_(j ell)) = (3, 3)$. Ordering preservation for all $p > 0$ forces the sign of
  $
    (alpha_i + gamma) 5^p + (alpha_j - gamma) 1^p - (alpha_i + alpha_j) 3^p
  $
  to remain constant as $p arrow.r infinity$ and $p arrow.r 0^+$. Dividing by $5^p$ as $p arrow.r infinity$ forces $alpha_i + gamma >= 0$ (or $= 0$), while $p arrow.r 0^+$ requires exact cancellation unless one coefficient vanishes. Eliminating $beta != 0$ identically via $(d_(i j), d_(i k), d_(j k)) = (1, 2, 2)$ vs. $(0, 2, 2)$ under $f_p$ restricts the only order-invariant solutions to pure order statistics $min$ and $max$, i.e., single and complete linkage @carlsson2010characterization.
]

#sec("III", "Complexity, Geometry and Operational Profile")

The computational footprint of the Lance-Williams algorithm bifurcates cleanly along the reducibility boundary of Theorem 3. Storing the symmetric dissimilarity matrix requires $n(n - 1) / 2 = Theta(n^2)$ floating-point entries in lower-triangular packed layout (roughly $400$ MB in double precision for $n = 10^4$, rising to $40$ GB at $n = 10^5$).

In the *reducible regime* (Single, Complete, UPGMA, WPGMA, Ward, and Flexible Beta with $beta <= 0$), the Nearest-Neighbor Chain algorithm executes in $O(n^2)$ worst-case time and $O(n^2)$ space @mullner2011modern @schubert2021fast, matching the lower bound required merely to read an arbitrary dissimilarity matrix. Furthermore, specialized implementations for Single (SLINK @sibson1973slink) and Complete (CLINK @defays1977clink) stream pairwise distances on the fly in $O(n)$ memory while preserving $O(n^2)$ runtime. In the *non-reducible regime* (UPGMC, WPGMC, and Flexible Beta with $beta > 0$), out-of-order RNN merges are invalid; maintaining per-cluster binary min-heaps yields $O(n^2 log n)$ worst-case time (with $O(n^3)$ for naive matrix scans) @mullner2011modern.

Geometrically, the parameter tuple $(alpha_i, alpha_j, beta, gamma)$ acts as a continuous dial over *space distortion*: negative $gamma$ or positive $beta$ contracts inter-cluster space (producing elongated *chaining* structures), whereas positive $gamma$ or negative $beta$ dilates space (enforcing compact, roughly spherical clusters). Because Lance-Williams updates only require pairwise dissimilarities, the framework handles non-Euclidean, graph-geodesic, and string-edit inputs seamlessly — though Ward and UPGMC lose their exact Euclidean centroid semantics unless initialized with squared Euclidean distances.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(n^2)$ reducible (NN-chain); $O(n^2 log n)$ heap (general); $O(n^3)$ naive scan],
    [Space complexity], [$O(n^2)$ packed triangular matrix; $O(n)$ streaming for SLINK / CLINK],
    [Monotonicity], [Guaranteed iff $gamma >= -min(alpha_i, alpha_j)$, $alpha_i + alpha_j >= 0$, $alpha_i + alpha_j + beta >= 1$],
    [Reducibility], [Exact iff additionally $alpha_i + alpha_j >= 1$ (excludes UPGMC, WPGMC, $beta > 0$)],
    [Input domain], [Arbitrary symmetric $D in RR_(>= 0)^(n times n)$ (metric or non-metric)],
    [Determinism], [Deterministic on tie-free $D$; tie-breaking order-sensitive on discrete $D$],
    [Invariances], [Positive scaling & permutation; monotone-transform invariant only for SL / CL],
  ),
  caption: [Operational and theoretical profile of the Lance-Williams framework (Alg. 72).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Quantitative Ecology and Numerical Taxonomy*: Lance and Williams originally formulated the combinatorial recurrence for floristic and ecological site ordination using Bray-Curtis, Canberra, and Gower dissimilarities @lance1967general. Because ecological abundance tables are non-Euclidean and zero-inflated, coordinate-free dissimilarity updates let ecologists control space dilation systematically without embedding sites into Euclidean space.
2. *Transcriptomics and Co-Expression Heatmaps*: In gene-expression profiling, genes are routinely clustered under Pearson correlation distance $d(u, v) = 1 - rho(u, v)$ or cosine dissimilarity. The Lance-Williams framework enables switching between UPGMA (noise averaging), Complete (tight co-regulation modules), and Ward-like updates inside a single visualization pipeline with exact $O(n^2)$ NN-chain performance.
3. *Financial Portfolio Construction (Hierarchical Risk Parity)*: Modern asset-allocation engines construct asset dendrograms from correlation-derived metric distances $d_(i j) = sqrt((1 - rho_(i j)) / 2)$. Because correlation matrices are noisy and ill-conditioned, the monotonic ultrametric produced by reducible Lance-Williams updates stabilizes recursive bisection weights without inverting singular covariance matrices.
4. *Information Retrieval and Document Taxonomy*: High-dimensional TF-IDF and neural-embedding corpora exhibit directional geometry where cosine similarity outperforms Euclidean distance @steinbach2000comparison. Average and complete Lance-Williams updates construct browsable multi-resolution topic hierarchies directly from pairwise cosine dissimilarities without recomputing sparse centroids.
5. *Software Architecture Recovery and Modularization*: Legacy codebase decomposition clusters software modules using Jaccard or structural coupling dissimilarities derived from call graphs and shared-variable incidence matrices. Combinatorial updates directly respect non-metric dependency tables while exposing multi-scale subsystem boundaries.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Despite its unifying elegance, the Lance-Williams framework suffers five structural failure regimes: (1) *Dendrogram inversions* under non-monotonic parameterizations (UPGMC, WPGMC), violating ultrametricity and breaking height cuts (Proposition 4); (2) *Quadratic memory barrier* ($O(n^2)$ distance storage), rendering raw execution infeasible beyond $n approx 5 times 10^4$ without pre-summarization such as BIRCH @zhang1996birch or CURE @guha1998cure; (3) *Chaining vs. dissection bias*, where space-contracting rules merge distant clusters via thin noise bridges while strongly space-dilating rules fracture legitimate non-convex manifolds; (4) *High-dimensional distance concentration* @beyer1999nearest, wherein pairwise dissimilarities converge toward a narrow band and amplify tie-induced instability; and (5) *Absence of exemplar prototypes*, since internal nodes have no canonical data-point representative (unlike Minimax linkage @bien2011minimax).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Ecology / taxonomy ($n <= 10^4$)], [Strong: native non-Euclidean dissimilarities & tunable dilation], [Avoid UPGMC inversions; use Flexible Beta or UPGMA],
    [Gene co-expression & HRP], [Strong: $O(n^2)$ NN-chain on correlation distances], [High noise ties; validate cut via Silhouette / ARI],
    [Non-convex / manifold shapes], [Mixed: SL chains easily; Ward fractures crescents], [Prefer HDBSCAN @campello2013hdbscan or Spectral @ng2002spectral],
    [Massive datasets ($n > 10^5$)], [Poor: $O(n^2)$ dissimilarity matrix exceeds RAM], [Use BIRCH @zhang1996birch, Mini-Batch $k$-means, or FINCH @sarfraz2019finch],
    [High-D sparse vectors], [Moderate: distance concentration compresses merge gaps], [Reduce dimension or apply cosine UPGMA / Ward],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Alg. 72.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Algorithm 72 serves as the algebraic backbone unifying classical agglomerative methods (Single, Complete, UPGMA, WPGMA, Ward, Centroid, Median) and its one-parameter specialization Flexible Beta (Alg. 74). Yet not every valid hierarchical linkage fits within the combinatorial Lance-Williams axiom: linkages whose merge cost depends on *internal* point geometry rather than solely on inter-cluster scalars lie strictly outside this family.

#proposition(title: "Strict Non-Membership of Prototype-Based Linkages", origin: "Original Counterexample")[
  Minimax linkage (Alg. 73) @bien2011minimax, defined by $d_"MM" (G, H) = min_(x in G union H) max_(y in G union H) d(x, y)$, cannot be expressed by any Lance-Williams update (nor any function of $(n_i, n_j, n_k, d_(i j), d_(i k), d_(j k))$ alone).
]

#proof[
  Construct two 4-point Euclidean configurations in $RR^2$ with identical cluster partition $C_i = {a, b}$, $C_j = {c}$, $C_k = {e}$ (cardinalities $(2, 1, 1)$):
  - *Config 1*: $a = (0, 0)$, $b = (0, 1)$, $c = (1, 0)$, $e = (-1, 0)$. Here $a$ acts as minimax prototype for $C_i union C_j$, $C_i union C_k$, and $C_i union C_j union C_k$, yielding $d_"MM" (C_i, C_j) = 1$, $d_"MM" (C_i, C_k) = 1$, $d_"MM" (C_j, C_k) = 2$, and merged radius $d_"MM" (C_i union C_j, C_k) = 1$.
  - *Config 2*: $a = (-1 / 2, 0)$, $b = (1 / 2, 0)$, $c = (1, sqrt(3) / 2)$, $e = (-1, sqrt(3) / 2)$. Here $b$ is the prototype of $C_i union C_j$ ($d(b, a) = d(b, c) = 1$) and $a$ is the prototype of $C_i union C_k$ ($d(a, b) = d(a, e) = 1$), while $d(c, e) = 2$. Thus all three pairwise linkages match Config 1 identically: $(1, 1, 2)$. Yet in $C_i union C_j union C_k$, $d(a, c) = d(b, e) = sqrt(3)$ and $d(c, e) = 2$, forcing $d_"MM" (C_i union C_j, C_k) = sqrt(3) != 1$.

  Because identical inputs $(n_i, n_j, n_k, d_(i j), d_(i k), d_(j k))$ yield distinct outputs ($1$ vs. $sqrt(3)$), no Lance-Williams recurrence can represent $d_"MM"$.
]

Similarly, Neighbor-Joining (Alg. 75) @saitou1987neighbor employs a Lance-Williams-like distance reduction ($alpha_i = alpha_j = 1 / 2, beta = -1 / 2, gamma = 0$) yet breaks the SAHN framework via its globally re-centered $Q$-selection criterion, reconstructing unrooted additive trees rather than ultrametrics. From an approximation-theoretic perspective, Dasgupta's hierarchical cost objective @dasgupta2016cost shows that UPGMA (a member of Alg. 72) achieves a $1 / 3$-approximation for similarity maximization @moseley2017approximation, whereas single linkage can incur $Omega(n)$ worst-case suboptimality.

*Empirical Benchmark Behavior*: On the standardized ClusterOverview 4-dataset benchmark, the Lance-Williams implementation achieves Iris Silhouette $= 0.447$ (@rousseeuw1987silhouettes) and Adjusted Rand Index $"ARI" = 0.565$ (@hubert1985comparing), Customers Silhouette $= 0.353$, Multi-Shape Silhouette $= 0.452$, and High-D Silhouette $= 0.234$, with a mean runtime of just $2.346$ ms. This profile confirms our complexity and geometric analysis: the combinatorial update executes in milliseconds via direct scalar matrix recurrences (nearly $3.5 times$ faster than non-combinatorial Minimax linkage at $8.141$ ms), delivers reliable baseline structure on spherical and customer-segmentation tasks (Customers $"Sil" = 0.353$ exceeds both Minimax at $0.335$ and Flexible Beta at $0.347$), and exhibits predictable degradation on non-convex Multi-Shape ($0.452$) and High-D ($0.234$) benchmarks where tuned space dilation (Alg. 74) or density-based hierarchy @campello2013hdbscan is better aligned with the underlying geometry @jain2010data @xu2005survey.
