#import "../template/ieee.typ": *

#alg_entry(
  id: 157,
  name: "Cheng-Church Biclustering",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2000",
  paradigm: "Greedy Mean-Squared-Residue Submatrix Extraction",
  complexity: $O(K n p (n + p))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical partitioning algorithms cluster rows across *all* columns
simultaneously, assuming global similarity across the entire feature
space. In functional genomics, however, co-regulated genes often display
coherent expression patterns only over a localized subset of experimental
conditions, while a single gene may participate in multiple overlapping
biological pathways. Cheng and Church @cheng2000biclustering formalized
this joint row–column discovery task—termed _biclustering_—by seeking
maximum-volume submatrices of a real data matrix
$bold(A) = [a_(i j)] in bb(R)^(n times p)$ whose entries conform to an
additive two-way shifting model up to a prescribed variance tolerance
$delta >= 0$.

Let $cal(R) = {1, dots, n}$ index rows (objects / genes) and
$cal(C) = {1, dots, p}$ index columns (features / conditions). For any
nonempty submatrix $bold(A)_(I J)$ indexed by $I subset.eq cal(R)$ and
$J subset.eq cal(C)$, define its row means, column means, and grand mean:
$
  a_(i J) = 1 / |J| sum_(j in J) a_(i j),
  quad a_(I j) = 1 / |I| sum_(i in I) a_(i j),
  quad a_(I J) = 1 / (|I| |J|) sum_(i in I) sum_(j in J) a_(i j).
$
The _residue_ of entry $(i, j)$ in submatrix $(I, J)$ measures its
departure from a purely additive row-plus-column profile:
$
  r_(i j) (I, J) = a_(i j) - a_(i J) - a_(I j) + a_(I J).
$
Aggregating squared residues yields the _Mean Squared Residue_ (MSR)
score $H(I, J)$, together with per-row scores $d(i)$ and per-column
scores $e(j)$:
$
  H(I, J) = 1 / (|I| |J|) sum_(i in I) sum_(j in J) r_(i j) (I, J)^2, \
  d(i) = 1 / |J| sum_(j in J) r_(i j) (I, J)^2,
  quad e(j) = 1 / |I| sum_(i in I) r_(i j) (I, J)^2.
$

#definition(title: "Delta-Bicluster and Maximum-Volume Optimization")[
  Given $bold(A) in bb(R)^(n times p)$ and coherence threshold $delta >= 0$,
  a submatrix $(I, J) in 2^cal(R) times 2^cal(C)$ ($|I|, |J| >= 1$) is a
  _$delta$-bicluster_ if $H(I, J) <= delta$. The _Maximum-Volume
  $delta$-Bicluster Problem_ seeks $(I^ast, J^ast) = arg max_(H(I, J) <= delta) |I| |J|$
  (typically restricted to $|I|, |J| >= 2$ to exclude trivial zero-residue
  single rows or columns). To discover $K$ overlapping biclusters
  sequentially, Cheng-Church masks each discovered submatrix
  $I_k times J_k$ with i.i.d. uniform random values over
  $[min_(a,b) a_(a b), max_(a,b) a_(a b)]$ before extracting the next
  bicluster @cheng2000biclustering.
]

Because exact maximization of $|I||J|$ subject to $H(I, J) <= delta$ is
NP-hard (Theorem 3), Cheng and Church @cheng2000biclustering introduced
a four-phase greedy deletion–addition procedure: (1) _multiple-node
deletion_ removes batches of high-residue rows/columns exceeding
$alpha H(I, J)$ ($alpha > 1$) for computational speedup; (2) _single-node
deletion_ iteratively removes the worst row or column until
$H(I, J) <= delta$; (3) _node addition_ greedily reinserts rows (optionally
sign-inverted, $-a_(i j)$, to capture anti-correlated genes) and columns
whose mean squared residue against $(I, J)$ does not exceed $H(I, J)$;
and (4) _random masking_ overwrites $bold(A)_(I_k J_k)$ to prevent
rediscovery in subsequent extractions.

#algbox(
  title: "Cheng-Church Sequential Delta-Biclustering",
  inputs: [Matrix $bold(A) in bb(R)^(n times p)$, MSR threshold $delta >= 0$, batch factor $alpha > 1$, bicluster count $K$, minimum batch threshold $n_"min"$.],
  outputs: [Sequence of $K$ submatrices ${(I_k, J_k)}_(k=1)^K$ with $H(I_k, J_k) <= delta$.],
)[
  1. Set working matrix $tilde(bold(A)) <- bold(A)$ with global range $[L, U] = [min bold(A), max bold(A)]$.
  2. For $k = 1, dots, K$: initialize $I <- cal(R)$, $J <- cal(C)$.
  3. *Multiple-Node Deletion:* While $H(I, J) > delta$ and $|I| > n_"min"$: compute $d(i)$ on $I$ and remove all $i in I$ with $d(i) > alpha H(I, J)$; recompute $H(I, J)$ and, if $|J| > n_"min"$, remove all $j in J$ with $e(j) > alpha H(I, J)$. Exit loop if no node is removed.
  4. *Single-Node Deletion:* While $H(I, J) > delta$ and $|I|, |J| > 1$: compute $i^ast = arg max_(i in I) d(i)$ and $j^ast = arg max_(j in J) e(j)$; if $d(i^ast) >= e(j^ast)$, set $I <- I without {i^ast}$, else $J <- J without {j^ast}$.
  5. *Node Addition:* Recompute $a_(i J), a_(I j), a_(I J)$ on $(I, J)$: add every $j in cal(C) without J$ with $e(j) <= H(I, J)$; then add every $i in cal(R) without I$ with $d(i) <= H(I, J)$ (and optionally sign-inverted rows with $d^-(i) <= H(I, J)$).
  6. *Record and Mask:* Store $(I_k, J_k) <- (I, J)$; replace $tilde(a)_(i j) tilde.op "Unif"[L, U]$ independently for all $(i, j) in I_k times J_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five structural and algorithmic results for Cheng-Church
biclustering: the doubly-centered orthogonal projection representation,
exact algebraic deletion identities proving strict descent for both
single- and $alpha$-batch deletion, exact residue preservation under
node addition, deterministic polynomial-bit NP-hardness, and closed-form
failure bounds on multiplicative patterns and random-masking overlap.

#subsec("A", "Doubly-Centered Orthogonal Projection and ANOVA Kernel")

#lemma(title: "Doubly-Centered Projection Identity and Zero-Residue Kernel", origin: "Original Derivation")[
  Let $bold(P)_I = bold(I)_(|I|) - |I|^(-1) bold(1)_I bold(1)_I^top$ and
  $bold(P)_J = bold(I)_(|J|) - |J|^(-1) bold(1)_J bold(1)_J^top$ be
  orthogonal centering projectors. Then the residue matrix
  $bold(R)_(I J) = [r_(i j) (I, J)]$ equals $bold(P)_I bold(A)_(I J) bold(P)_J$,
  and
  $
    H(I, J) &= 1 / (|I| |J|) norm(bold(P)_I bold(A)_(I J) bold(P)_J)_F^2 \
    &= 1 / (|I| |J|) min_(mu, bold(alpha), bold(beta))
       sum_(i in I) sum_(j in J) (a_(i j) - mu - alpha_i - beta_j)^2.
  $
  Consequently, $H(I, J) = 0$ if and only if $bold(A)_(I J)$ admits an
  exact additive decomposition $a_(i j) = mu + alpha_i + beta_j$.
]

#proof[
  Evaluating entry $(i, j)$ of $bold(P)_I bold(A)_(I J) bold(P)_J$ gives
  $(bold(I) - |I|^(-1) bold(1)_I bold(1)_I^top) bold(A)_(I J) (bold(I) - |J|^(-1) bold(1)_J bold(1)_J^top)_(i j) = a_(i j) - a_(i J) - a_(I j) + a_(I J) = r_(i j) (I, J)$,
  so $H(I, J) = (|I| |J|)^(-1) norm(bold(P)_I bold(A)_(I J) bold(P)_J)_F^2$.
  Define the linear subspace of additive matrices
  $cal(U)_(I J) = {mu bold(1)_I bold(1)_J^top + bold(alpha) bold(1)_J^top + bold(1)_I bold(beta)^top : mu in bb(R), bold(alpha) in bb(R)^(|I|), bold(beta) in bb(R)^(|J|)}$.
  Because $bold(P)_I bold(1)_I = bold(0)$ and $bold(1)_J^top bold(P)_J = bold(0)^top$,
  every $bold(M) in cal(U)_(I J)$ satisfies $bold(P)_I bold(M) bold(P)_J = bold(0)$
  and is Frobenius-orthogonal to $bold(R)_(I J) = bold(P)_I bold(A)_(I J) bold(P)_J$:
  $
    chevron.l bold(R)_(I J), bold(alpha) bold(1)_J^top chevron.r_F
    = "Tr"(bold(P)_J bold(A)_(I J)^top bold(P)_I bold(alpha) bold(1)_J^top)
    = bold(1)_J^top bold(P)_J bold(A)_(I J)^top bold(P)_I bold(alpha) = 0,
  $
  and symmetrically for $bold(1)_I bold(beta)^top$ and $mu bold(1)_I bold(1)_J^top$.
  Furthermore, $bold(A)_(I J) - bold(R)_(I J) = (a_(i J) + a_(I j) - a_(I J))_(i,j) in cal(U)_(I J)$.
  By the orthogonal projection theorem in $(bb(R)^(|I| times |J|), norm(dot)_F)$,
  $bold(A)_(I J) - bold(R)_(I J)$ is the unique best Frobenius
  approximation to $bold(A)_(I J)$ in $cal(U)_(I J)$, so
  $min_(bold(M) in cal(U)_(I J)) norm(bold(A)_(I J) - bold(M))_F^2 = norm(bold(R)_(I J))_F^2 = |I| |J| H(I, J)$.
  Hence $H(I, J) = 0$ iff $bold(A)_(I J) in cal(U)_(I J)$.
]

#subsec("B", "Exact Deletion Identities and Strict Monotonicity")

#theorem(title: "Exact Single-Node Identity and Alpha-Batch Contraction Guarantee", origin: "Original Proof")[
  Assume $|I| >= 2$, $|J| >= 1$, and $H(I, J) > 0$.
  1. *(Exact Single-Node Identity)* Removing any single row $i^ast in I$
     transforms the Mean Squared Residue identically as
     $
       H(I without {i^ast}, J) = |I| / (|I| - 1)
       (H(I, J) - d(i^ast) / (|I| - 1)).
     $
     In particular, removing $i^ast = arg max_(i in I) d(i)$ guarantees
     strict contraction
     $H(I without {i^ast}, J) <= (1 - (|I| - 1)^(-2)) H(I, J) < H(I, J)$.
  2. *(Alpha-Batch Contraction)* For any threshold $alpha > 1$, let
     $S = {i in I : d(i) > alpha H(I, J)}$ with $s = |S| >= 1$. Then
     automatically $s < |I| / alpha < |I|$ (so $I without S != emptyset$)
     and simultaneous batch deletion strictly contracts $H$:
     $
       H(I without S, J) <= (1 - ((alpha - 1) s) / (|I| - s)) H(I, J) < H(I, J).
     $
  Symmetric identities hold identically for column deletions via $e(j)$.
]

#proof[
  For fixed $J$, define row-centered vectors
  $tilde(bold(a))_i = (a_(i j) - a_(i J))_(j in J) in bb(R)^(|J|)$ and
  their sample mean $bar(bold(c))(I) = |I|^(-1) sum_(i in I) tilde(bold(a))_i$.
  Since $(bar(bold(c))(I))_j = a_(I j) - a_(I J)$, we have
  $r_(i j) (I, J) = (tilde(bold(a))_i - bar(bold(c))(I))_j$, hence
  $|J| d(i) = norm(tilde(bold(a))_i - bar(bold(c))(I))_2^2$ and
  $|I| |J| H(I, J) = sum_(i in I) norm(tilde(bold(a))_i - bar(bold(c))(I))_2^2$.
  For any nonempty proper subset $S subset I$ ($|S| = s$, $I' = I without S$),
  the Huygens–Steiner parallel-axis identity on $I'$ around $bar(bold(c))(I')$
  gives
  $
    (|I| - s) |J| H(I', J) &= sum_(i in I') norm(tilde(bold(a))_i - bar(bold(c))(I))_2^2 \
    &- (|I| - s) norm(bar(bold(c))(I') - bar(bold(c))(I))_2^2.
  $
  The first term equals $|I| |J| H(I, J) - |J| sum_(i in S) d(i)$.

  For (1), set $S = {i^ast}$ ($s = 1$). Since
  $(|I| - 1) bar(bold(c))(I') = |I| bar(bold(c))(I) - tilde(bold(a))_(i^ast)$,
  we have $bar(bold(c))(I') - bar(bold(c))(I) = -(|I| - 1)^(-1) (tilde(bold(a))_(i^ast) - bar(bold(c))(I))$,
  so the second term equals $(|I| - 1)^(-1) |J| d(i^ast)$. Combining both
  terms yields
  $
    (|I| - 1) |J| H(I', J) = |I| |J| H(I, J)
    - |I| / (|I| - 1) |J| d(i^ast).
  $
  Dividing by $(|I| - 1) |J|$ establishes the exact identity. Because
  $|I|^(-1) sum_(i in I) d(i) = H(I, J)$, the maximizer $i^ast$ obeys
  $d(i^ast) >= H(I, J)$, giving
  $H(I', J) <= |I|(|I| - 2) / (|I| - 1)^2 H(I, J) = (1 - (|I| - 1)^(-2)) H(I, J)$.

  For (2), since $d(i) >= 0$ for all $i in I$,
  $|I| H(I, J) = sum_(i in I) d(i) >= sum_(i in S) d(i) > s alpha H(I, J)$,
  so $H(I, J) > 0$ forces $s alpha < |I|$, guaranteeing $I' != emptyset$.
  Dropping the non-positive shift term
  $-(|I| - s) norm(bar(bold(c))(I') - bar(bold(c))(I))_2^2 <= 0$ yields
  $(|I| - s) H(I', J) <= |I| H(I, J) - sum_(i in S) d(i) < (|I| - alpha s) H(I, J)$.
  Dividing by $|I| - s > 0$ yields
  $H(I', J) < (|I| - alpha s) / (|I| - s) H(I, J) = (1 - ((alpha - 1) s) / (|I| - s)) H(I, J)$.
  Column deletions follow identically by transposing $bold(A)$.
]

#subsec("C", "Node Addition Guarantee and NP-Hardness")

#theorem(title: "Residue Preservation under Row and Column Addition", origin: "Classical — Cited")[
  Let $(I, J)$ satisfy $H(I, J) <= delta$, and let $S subset.eq cal(R) without I$
  be any set of candidate rows whose residue scores against the current
  column means of $(I, J)$ satisfy $d(i) <= H(I, J)$ for all $i in S$
  (equivalently $d^-(i) <= H(I, J)$ for sign-inverted rows). Then
  $H(I union S, J) <= H(I, J) <= delta$ @cheng2000biclustering. For a
  single added row $i$, the exact update is
  $H(I union {i}, J) = |I| / (|I| + 1) H(I, J) + |I| / (|I| + 1)^2 d(i)$.
]

#proof[
  Let $I^+ = I union S$ with $s = |S| >= 1$. Using the row-centered
  representation of Theorem 1 (replacing $tilde(bold(a))_i$ by
  $-tilde(bold(a))_i$ for inverted rows), the parallel-axis identity on
  $I^+$ around $bar(bold(c))(I)$ yields
  $
    (|I| + s) |J| H(I^+, J) &= sum_(i in I^+) norm(tilde(bold(a))_i - bar(bold(c))(I))_2^2 \
    &- (|I| + s) norm(bar(bold(c))(I^+) - bar(bold(c))(I))_2^2.
  $
  The first sum equals $|I| |J| H(I, J) + |J| sum_(i in S) d(i) <= (|I| + s) |J| H(I, J)$.
  Because the second term is non-positive, dividing by $(|I| + s) |J|$
  immediately gives $H(I^+, J) <= H(I, J) <= delta$. For $S = {i}$,
  $bar(bold(c))(I^+) - bar(bold(c))(I) = (|I| + 1)^(-1) (tilde(bold(a))_i - bar(bold(c))(I))$,
  so the subtracted shift equals $(|I| + 1)^(-1) |J| d(i)$, yielding the
  exact single-row identity upon dividing by $(|I| + 1) |J|$. Symmetric
  reasoning holds for column additions $J union T$.
]

#theorem(title: "NP-Hardness of Maximum-Volume Delta-Biclustering", origin: "Classical — Cited")[
  Deciding whether $bold(A) in bb(R)^(n times p)$ contains a
  $delta$-bicluster $(I, J)$ with $|I|, |J| >= 2$ and volume
  $|I| |J| >= tau$ is NP-hard, even at $delta = 0$ @cheng2000biclustering
  @prelic2006bimax.
]

#proof[
  We reduce from the NP-hard *Maximum Edge Biclique* problem: given a
  bipartite graph $G = (U, W, E)$ ($|U| = n$, $|W| = p$) and integer
  $tau > max(n, p)$, decide whether $G$ has a biclique $(I, J)$
  ($I times J subset.eq E$) with $|I| |J| >= tau$. Construct
  $bold(A) in bb(R)^(n times p)$ deterministically by setting
  $a_(i j) = 0$ if $(i, j) in E$ and $a_(i j) = 2^((i - 1) p + j)$ if
  $(i, j) in.not E$. Each entry uses at most $n p + 1$ bits, so the
  total encoding size is $O((n p)^2)$, strictly polynomial in $(n, p)$.

  If $(I, J)$ is a biclique of $G$ with $|I| |J| >= tau$, then
  $bold(A)_(I J) = bold(0)$, so $H(I, J) = 0$ and $|I|, |J| >= 2$
  (as $|I| |J| >= tau > max(n, p)$). Conversely, suppose $(I, J)$ has
  $H(I, J) = 0$ and $|I| |J| >= tau > max(n, p)$, forcing $|I|, |J| >= 2$.
  By Lemma 1, $a_(i_1 j_1) - a_(i_1 j_2) - a_(i_2 j_1) + a_(i_2 j_2) = 0$
  for every $2 times 2$ submatrix of $bold(A)_(I J)$. Suppose for
  contradiction that $I times J$ contains at least one non-edge, and let
  $(i_1, j_1) in (I times J) without E$ maximize the exponent
  $e_"max" = (i_1 - 1) p + j_1$. Choose any $i_2 in I without {i_1}$ and
  $j_2 in J without {j_1}$. The remaining three entries are either $0$
  or distinct powers of $2$ strictly below $2^(e_"max")$, so
  $a_(i_1 j_2) + a_(i_2 j_1) - a_(i_2 j_2) <= 2^(e_"max" - 1) + 2^(e_"max" - 2) = (3/4) 2^(e_"max") < 2^(e_"max") = a_(i_1 j_1)$,
  contradicting zero $2 times 2$ contrast. Hence $I times J subset.eq E$,
  so $(I, J)$ is a biclique of size $|I| |J| >= tau$.
]

#subsec("D", "Failure on Multiplicative Co-Regulation and Masking Distortion")

#proposition(title: "Multiplicative Coherence Breakdown and Random-Masking Interference Bound", origin: "Original Proof")[
  1. *(Multiplicative Pattern Breakdown)* If $bold(A)_(I J)$ is a
     noiseless rank-1 multiplicative bicluster $a_(i j) = u_i v_j$ with
     empirical variances $sigma_u^2 = |I|^(-1) sum_(i in I) (u_i - bar(u))^2$
     and $sigma_v^2 = |J|^(-1) sum_(j in J) (v_j - bar(v))^2$, then its
     exact MSR is $H(I, J) = sigma_u^2 sigma_v^2$.
  2. *(Random-Masking Interference)* Suppose a noiseless additive
     bicluster $(I_2, J_2)$ overlaps a previously masked bicluster
     $(I_1, J_1)$ on fraction
     $rho = (|I_1 inter I_2| |J_1 inter J_2|) / (|I_2| |J_2|) in (0, 1)$,
     where masked entries are replaced by i.i.d. $Z_(i j) tilde.op "Unif"[L, U]$
     of variance $sigma_Z^2 = (U - L)^2 / 12$. Then the expected MSR of
     $(I_2, J_2)$ in the masked matrix $tilde(bold(A))$ obeys
     $
       bb(E)[H(tilde(bold(A))_(I_2 J_2))] >=
       rho (1 - 1 / |I_2|) (1 - 1 / |J_2|) (U - L)^2 / 12.
     $
]

#proof[
  For (1), direct substitution gives $a_(i J) = u_i bar(v)$,
  $a_(I j) = bar(u) v_j$, and $a_(I J) = bar(u) bar(v)$, so
  $r_(i j) = u_i v_j - u_i bar(v) - bar(u) v_j + bar(u) bar(v) = (u_i - bar(u))(v_j - bar(v))$.
  Squaring and averaging over $I times J$ factors identically into
  $H(I, J) = (|I|^(-1) sum_(i in I) (u_i - bar(u))^2)(|J|^(-1) sum_(j in J) (v_j - bar(v))^2) = sigma_u^2 sigma_v^2$.

  For (2), let $I_0 = I_1 inter I_2$, $J_0 = J_1 inter J_2$, $r = |I_2|$,
  and $c = |J_2|$. Write $tilde(bold(A))_(I_2 J_2) = bold(A)_(I_2 J_2) + bold(Delta)$
  where $Delta_(i j) = (Z_(i j) - a_(i j)) bb(I)((i, j) in I_0 times J_0)$.
  Because $bold(P)_(I_2) bold(A)_(I_2 J_2) bold(P)_(J_2) = bold(0)$,
  Lemma 1 gives $r c H(tilde(bold(A))_(I_2 J_2)) = norm(bold(P)_(I_2) bold(Delta) bold(P)_(J_2))_F^2$.
  Decompose $bold(Delta) = bb(E)[bold(Delta)] + bold(E)$, where
  $E_(i j) = (Z_(i j) - bb(E)[Z_(i j)]) bb(I)((i, j) in I_0 times J_0)$
  has zero mean and variance $sigma_Z^2$ independently on $I_0 times J_0$.
  By variance–bias orthogonality,
  $
    bb(E)[norm(bold(P)_(I_2) bold(Delta) bold(P)_(J_2))_F^2]
    >= bb(E)[norm(bold(P)_(I_2) bold(E) bold(P)_(J_2))_F^2] \
    = sum_(i in I_0) sum_(j in J_0) sigma_Z^2 ((bold(P)_(I_2)^2)_(i i))
      ((bold(P)_(J_2)^2)_(j j)).
  $
  Since $bold(P)_(I_2)$ is idempotent with diagonal entries $1 - 1 / r$
  and $bold(P)_(J_2)$ has diagonal entries $1 - 1 / c$, the double sum
  equals $sigma_Z^2 |I_0| |J_0| (1 - 1 / r)(1 - 1 / c)$. Dividing by
  $r c$ and substituting $rho = |I_0| |J_0| / (r c)$ and
  $sigma_Z^2 = (U - L)^2 / 12$ completes the proof.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In a single extraction, computing row means, column means, grand mean,
$H(I, J)$, and all scores $d(i), e(j)$ from scratch costs $O(|I| |J|)$
operations. Because deleting one row shifts every column mean $a_(I j)$,
re-evaluating all row scores $d(i)$ after a single deletion requires
$O(|I| |J|)$ time unless an $O(n^2)$ Gram matrix is cached. In the worst
case of $O(n + p)$ single-node deletions, one bicluster requires
$O(n p (n + p))$ time, and extracting $K$ biclusters costs
$O(K n p (n + p))$ time with $O(n p)$ working memory for the masked copy
$tilde(bold(A))$. When multiple-node deletion ($alpha > 1$, typically
$alpha = 1.2$) is active, each pass removes a geometric fraction of
excess-residue rows and columns, reducing typical practical runtime to
$O(K n p log(n + p))$ @cheng2000biclustering.

Geometrically, Cheng-Church does not partition $bb(R)^p$ into convex
Voronoi cells; instead it carves out axis-aligned local subspaces
$I times J$ lying within Frobenius radius $sqrt(|I| |J| delta)$ of the
$(|I| + |J| - 1)$-dimensional two-way additive affine subspace
$cal(U)_(I J)$. It is invariant to additive row shifts, additive column
shifts, and global matrix offsets, as well as row/column permutations,
but *not* invariant to non-uniform column rescaling or orthogonal feature
rotations.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(K n p (n + p))$ worst-case; $O(K n p log(n + p))$ typical with $alpha$-batch deletion],
    [Space Complexity], [$O(n p)$ (working matrix $tilde(bold(A))$ plus index bitmasks)],
    [Coherence Model], [Two-way additive shifting $a_(i j) approx mu + alpha_i + beta_j$ (zero interaction)],
    [Invariances], [Exact invariance to row shifts, column shifts, and matrix permutations],
    [Hyperparameters], [MSR tolerance $delta >= 0$, batch ratio $alpha > 1$, bicluster count $K$],
    [Overlap Handling], [Sequential uniform random masking (induces overlap distortion per Prop. 1)],
  ),
  caption: [Complexity, invariance, and operational profile of Cheng-Church biclustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *DNA Microarray and RNA-Seq Regulatory Module Discovery:* Under
   log-transformed expression levels, multiplicative fold-change
   co-regulation ($g_(i j) = u_i v_j$) transforms into an exact additive
   shifting pattern ($log g_(i j) = log u_i + log v_j$) with $H(I, J) = 0$
   (Lemma 1). Cheng-Church @cheng2000biclustering directly recovers
   condition-specific transcriptional modules while accommodating
   baseline gene-level and array-level offsets.
2. *Collaborative Filtering and Localized Preference Mining:* In
   user–item rating matrices, subgroups of users often evaluate a shared
   niche of items with consistent relative rankings up to per-user
   generosity offsets ($alpha_i$) and per-item popularity baselines
   ($beta_j$). Low-MSR submatrices correspond directly to coherent
   preference neighborhoods immune to user rating bias.
3. *Pharmacogenomics and Drug-Sensitivity Profiling:* Cancer cell line
   drug-response matrices (e.g., log-$"IC"_(50)$ across cell lines and
   compounds) contain localized pathway-driven sensitivity blocks where
   a subset of tumor lines responds consistently across a drug class.
   Additive MSR captures shared mechanism-of-action response profiles
   despite baseline tumor resistance shifts.
4. *Multi-Sensor Financial and Macroeconomic Regime Detection:* Asset
   return panels across sectors and time windows exhibit localized
   factor regimes in which a subset of equities co-move synchronously
   over specific trading epochs up to idiosyncratic drift and market-wide
   shocks—precisely matching $a_(i j) = mu + alpha_i + beta_j$.
5. *Text-Document Topic Submatrix Extraction:* In log-weighted
   document–term matrices, topical corpora share elevated relative
   co-occurrence patterns across specialized vocabularies despite
   differing document lengths (row offsets) and baseline term
   frequencies (column offsets).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure modes limit Cheng-Church in practice. First, *purely
multiplicative unlogged data* fails Proposition 1(1) because
$H = sigma_u^2 sigma_v^2 > 0$, shattering coherent scaling patterns.
Second, *random-masking interference* (Proposition 1(2)) injects uniform
noise variance $(U - L)^2 / 12$ into overlapping regions, corrupting
subsequently extracted overlapping biclusters @prelic2006bimax. Third,
*trivial constant / flat-gene bias*: a submatrix of near-zero variance
trivially achieves $H approx 0$ even when biologically uninformative,
requiring secondary row-variance filters. Fourth, *greedy myopic deletion*
can prematurely drop a core row or column during early high-residue
stages when competing biclusters inflate global means. Fifth, *low-column
tabular datasets* ($p <= 6$) lack enough columns for meaningful column
subsetting, collapsing biclustering into noisy row trimming.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Log-Gene Expression ($n, p >> 20$)], [Strong: log-fold changes are exact additive blocks], [Filter low-variance constant blocks; compare ISA @bergmann2003isa],
    [Heavily Overlapping Modules], [Weak: uniform masking corrupts shared entries], [Use additive layer model PLAID (Alg. 158) @lazzeroni2002plaid],
    [Multiplicative / Scaling Patterns], [Poor without log transform ($H = sigma_u^2 sigma_v^2$)], [Log-transform first, or use Spectral Biclustering @kluger2003spectral / FABIA @hochreiter2010fabia],
    [Binary Presence–Absence Matrices], [Moderate: additive MSR ignores 0/1 asymmetry], [Use Bimax @prelic2006bimax for exact maximal 1-submatrices],
    [Low-Dim Euclidean Clustering ($p <= 5$)], [Poor: shift-invariance contradicts metric proximity], [Use KFCM (Alg. 155), sFCM (Alg. 156), or subspace CLIQUE @agrawal1998clique],
  ),
  caption: [Domain suitability and failure-mode matrix for Cheng-Church biclustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Cheng-Church pioneered residual-based biclustering and anchors a rich
family of matrix co-clustering methods: PLAID (Alg. 158)
@lazzeroni2002plaid replaces destructive random masking with additive
layer subtraction; Spectral Biclustering @kluger2003spectral extracts
checkerboard multiplicative structures via SVD; ISA @bergmann2003isa
discovers high-signal signature modules by thresholded alternating
projections; FABIA @hochreiter2010fabia formulates sparse Bayesian
multiplicative factor analysis; and CLIQUE @agrawal1998clique performs
grid-based dense subspace clustering @kriegel2009subspacesurvey.

To explain why Cheng-Church exhibits low Euclidean cluster validity on
standard tabular benchmarks, we prove an exact geometric mismatch between
MSR coherence and Euclidean distance.

#proposition(title: "Unbounded Euclidean Divergence of Zero-MSR Biclusters", origin: "Original Derivation")[
  For any dimension $p >= 2$ and any constant $M > 0$, there exists a
  submatrix $bold(A)_(I J)$ with $J = {1, dots, p}$ achieving exact zero
  Mean Squared Residue ($H(I, J) = 0$) whose pairwise row Euclidean
  distances satisfy $norm(bold(a)_i - bold(a)_(i'))_2 >= M$ for all
  $i != i' in I$, while simultaneously admitting an external row
  $bold(a)_0 in.not I$ arbitrarily close to some $bold(a)_i$ in
  Euclidean distance yet violating $H(I union {0}, J) <= delta$.
]

#proof[
  Fix any template vector $bold(beta) in bb(R)^p$ and define rows
  $bold(a)_i = i (M / sqrt(p)) bold(1)_p + bold(beta)$ for $i = 1, dots, r$.
  Since $a_(i j) = alpha_i + beta_j$ with $alpha_i = i M / sqrt(p)$,
  Lemma 1 yields $H(I, J) = 0$ identically, yet for all $i != i'$,
  $norm(bold(a)_i - bold(a)_(i'))_2 = |i - i'| (M / sqrt(p)) norm(bold(1)_p)_2 = |i - i'| M >= M$.
  Conversely, pick $bold(a)_0 = bold(a)_1 + bold(eta)$ where
  $bold(eta) = (epsilon, -epsilon, 0, dots, 0)^top$ has small Euclidean
  norm $sqrt(2) epsilon$, while its zero-sum perturbation is orthogonal
  to $bold(1)_p$ and adds residue
  $d(0) = p^(-1) (2 epsilon^2) > delta$ whenever $epsilon > sqrt(p delta / 2)$.
  Thus zero-MSR coherence is geometrically orthogonal to Euclidean
  compactness.
]

*Empirical Benchmark Behavior (ClusterOverview Suite).* The four-dataset
ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and
Adjusted Rand Index @hubert1985comparing) directly illustrates
Proposition 2. Across all datasets, Cheng-Church achieves a very fast
mean runtime of $0.707 "ms"$—nearly $4 times$ faster than KFCM
($2.723 "ms"$)—due to rapid batch deletion on moderate matrix dimensions.
However, when projected to global hard row partitions for Euclidean
benchmarking, it scores $"Sil" = 0.175$, $"ARI" = 0.183$ on *Iris*,
$"Sil" = 0.116$ on *Customers*, $"Sil" = 0.404$ ($"ARI" = "N/A"$) on
*Multi-Shape*, and $"Sil" = 0.043$ ($"ARI" = "N/A"$) on *High-D*. These
scores are structurally expected rather than anomalous: (i) on *Iris*
($p = 4$) and *Customers*, additive shift-invariance groups rows with
parallel feature profiles separated by large Euclidean offsets (Prop. 2)
and deletes discriminative columns; (ii) on *Multi-Shape*, higher
silhouette ($0.404$) arises because separated coordinate blocks induce
localized constant-value submatrices; and (iii) on *High-D*, sequential
random masking (Prop. 1(2)) and partial column coverage leave many rows
in background residuals, underscoring that Cheng-Church must be evaluated
by submatrix coherence and biological enrichment rather than global
isotropic Silhouette @prelic2006bimax.
