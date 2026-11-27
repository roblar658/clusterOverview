#import "../template/ieee.typ": *

#alg_entry(
  id: 159,
  name: "Bimax (Binary Inclusion-Maximal Biclustering)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2006",
  paradigm: "Divide-and-Conquer Exact Binary Biclustering",
  complexity: $O(n m beta min(n, m))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Biclustering seeks simultaneous subsets of rows (objects, genes, or users) and columns (features, experimental conditions, or items) that exhibit coherent local structure across only a subset of dimensions, rather than partitioning objects globally across all features @cheng2000biclustering @prelic2006bimax. Whereas continuous formulations optimize heuristic variance or residue surrogates with no global optimality guarantees, *Bimax* (Binary inclusion-maximal biclustering) @prelic2006bimax introduces an exact combinatorial reference model on binarized matrices.

Let $X in bb(R)^(n times m)$ denote a real-valued data matrix with row index set $R = {1, dots, n}$ and column index set $C = {1, dots, m}$. A threshold binarization operator $phi_tau : bb(R)^(n times m) -> {0, 1}^(n times m)$ maps $X$ to a binary incidence matrix $E in {0, 1}^(n times m)$, typically via $E_(i j) = bold(1)[X_(i j) >= tau]$ (or two-sided thresholding $|X_(i j) - mu| >= tau$). For each row $i in R$ and column $j in C$, define the active support sets
$ C(i) = {j in C : E_(i j) = 1}, quad R(j) = {i in R : E_(i j) = 1}. $

#definition(title: "Inclusion-Maximal Binary Bicluster")[
  Given $E in {0, 1}^(n times m)$, a pair $(I, J) in 2^R times 2^C$ is a _valid all-ones bicluster_ if $I != emptyset$, $J != emptyset$, and $E_(i j) = 1$ for every $(i, j) in I times J$. A valid bicluster $(I, J)$ is _inclusion-maximal_ if no valid bicluster $(I', J') != (I, J)$ satisfies $I subset.eq I'$ and $J subset.eq J'$. We denote the complete set of inclusion-maximal all-ones biclusters of $E$ by $cal(B)^star (E)$ and its cardinality by $beta(E) = |cal(B)^star (E)|$.
]

To enumerate $cal(B)^star (E)$ exactly without generating duplicates or explicitly testing subset inclusion against previously discovered blocks, Bimax employs a recursive divide-and-conquer partitioning scheme (`Conquer`) parameterized by a current row set $U subset.eq R$, column set $V subset.eq C$, and a family of _mandatory column sets_ $cal(Z) = {Z_1, dots, Z_d}$ ($Z_ell subset.eq V$). Any bicluster $(I, J)$ emitted from a recursive state $(U, V, cal(Z))$ must satisfy $I subset.eq U$, $J subset.eq V$, and the mandatory intersection invariant $J inter Z_ell != emptyset$ for all $Z_ell in cal(Z)$.

At each non-uniform submatrix $E_(U, V)$, Bimax selects a _pivot row_ $r^star in U$ that is neither all-zero nor all-one on $V$, splits $V$ into the pivot's active columns $C_U = C(r^star) inter V$ and inactive columns $C_V = V without C_U$, and partitions the active rows of $U$ into three disjoint blocks: rows active solely in $C_U$ ($G_U$), rows active solely in $C_V$ ($G_V$), and bridging rows active in both $C_U$ and $C_V$ ($G_W$). Because $E_(i j) = 0$ on $G_U times C_V$ and $G_V times C_U$, any maximal bicluster of $U times V$ lies exclusively inside either $(G_U union G_W) times C_U$ or $(G_W union G_V) times V$ with at least one column in $C_V$.

#algbox(
  title: "Bimax Divide-and-Conquer Enumeration (Conquer)",
  inputs: [Binary matrix $E in {0, 1}^(n times m)$; optional minimum size bounds $(r_"min", c_"min")$; initial call $(U, V, cal(Z)) = (R, C, emptyset)$.],
  outputs: [Exact set $cal(B)^star (E)$ of all inclusion-maximal all-ones biclusters (filtered by $|I| >= r_"min"$, $|J| >= c_"min"$).],
)[
  1. *Feasibility Pruning:* Filter $U <- {i in U : C(i) inter V != emptyset "and" forall Z in cal(Z), C(i) inter V inter Z != emptyset}$ and $V <- {j in V : R(j) inter U != emptyset}$. If $U = emptyset$ or $V = emptyset$ (or $|U| < r_"min"$ or $|V| < c_"min"$), return $emptyset$.
  2. *Uniform Base Case:* If $E_(i j) = 1$ for all $(i, j) in U times V$, return ${(U, V)}$.
  3. *Pivot Selection:* Choose a row $r^star in U$ satisfying $0 < |C(r^star) inter V| < |V|$.
  4. *Column Partition:* Set $C_U <- C(r^star) inter V$ and $C_V <- V without C_U$.
  5. *Row Partition:* Partition $U$ into three disjoint subsets:
     $ G_U &= {i in U : C(i) inter V subset.eq C_U}, \
       G_V &= {i in U : C(i) inter V subset.eq C_V}, \
       G_W &= U without (G_U union G_V). $
  6. *Left Recursion (Subproblem A):* Set $cal(Z)_A <- {Z inter C_U : Z in cal(Z)}$. If every $Z' in cal(Z)_A$ is non-empty, compute $cal(M)_A <- "Conquer"(G_U union G_W, C_U, cal(Z)_A)$; else $cal(M)_A <- emptyset$.
  7. *Right Recursion (Subproblem B):* If $G_W = emptyset$, set $V_B <- C_V$ and $cal(Z)_B <- {Z inter C_V : Z in cal(Z)}$; else set $V_B <- V$ and $cal(Z)_B <- cal(Z) union {C_V}$. Compute $cal(M)_B <- "Conquer"(G_W union G_V, V_B, cal(Z)_B)$.
  8. *Aggregate:* Return $cal(M)_A union cal(M)_B$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Galois Connection and Concept Lattice Characterization")

We first establish an algebraic bijection between Bimax's combinatorial objective and the formal concept lattice of the binary relation $E$. Define the pair of derivation operators $(arrow.t, arrow.b)$ between $2^R$ and $2^C$ by
$ I^(arrow.t) = inter.big_(i in I) C(i), quad J^(arrow.b) = inter.big_(j in J) R(j), $
with the conventions $emptyset^(arrow.t) = C$ and $emptyset^(arrow.b) = R$.

#lemma(title: "Galois Closure Equivalence of Maximal Biclusters", origin: "Original Unified Proof")[
  The pair $(arrow.t, arrow.b)$ forms an antitone Galois connection between $(2^R, subset.eq)$ and $(2^C, subset.eq)$. A pair $(I, J)$ with $I != emptyset$ and $J != emptyset$ belongs to $cal(B)^star (E)$ if and only if $I = J^(arrow.b)$ and $J = I^(arrow.t)$. Equivalently, $I$ is a non-empty closed row set ($I = (I^(arrow.t))^(arrow.b)$ with $I^(arrow.t) != emptyset$) and $J = I^(arrow.t)$. Moreover, the partial order $(I_1, J_1) <= (I_2, J_2) <==> I_1 subset.eq I_2$ is anti-isomorphic to column inclusion ($J_2 subset.eq J_1$).
]

#proof[
  For any $I subset.eq R$ and $J subset.eq C$, we have $J subset.eq I^(arrow.t)$ iff $forall i in I, J subset.eq C(i)$, iff $forall (i, j) in I times J, E_(i j) = 1$, iff $forall j in J, I subset.eq R(j)$, iff $I subset.eq J^(arrow.b)$. Thus $J subset.eq I^(arrow.t) <==> I subset.eq J^(arrow.b)$, proving that $(arrow.t, arrow.b)$ is an antitone Galois connection. Standard consequences follow immediately: both maps are order-reversing ($I_1 subset.eq I_2 ==> I_2^(arrow.t) subset.eq I_1^(arrow.t)$), $I subset.eq (I^(arrow.t))^(arrow.b)$, $J subset.eq (J^(arrow.b))^(arrow.t)$, and $((I^(arrow.t))^(arrow.b))^(arrow.t) = I^(arrow.t)$.

  Suppose $(I, J) in cal(B)^star (E)$. Since $E_(i j) = 1$ on $I times J$, we have $I subset.eq J^(arrow.b)$ and $J subset.eq I^(arrow.t)$. Because $(I, I^(arrow.t))$ is a valid all-ones bicluster containing $(I, J)$, inclusion-maximality forces $J = I^(arrow.t)$; symmetrically $(J^(arrow.b), J)$ is valid and contains $(I, J)$, forcing $I = J^(arrow.b)$. Conversely, suppose $I = J^(arrow.b) != emptyset$ and $J = I^(arrow.t) != emptyset$. Then $J subset.eq I^(arrow.t)$ implies $E_(i j) = 1$ on $I times J$. If $(I', J')$ is any valid bicluster with $I subset.eq I'$ and $J subset.eq J'$, antitonicity yields $J' subset.eq (I')^(arrow.t) subset.eq I^(arrow.t) = J$ and $I' subset.eq (J')^(arrow.b) subset.eq J^(arrow.b) = I$, so $(I', J') = (I, J) in cal(B)^star (E)$. Finally, for any two closed pairs $(I_1, J_1), (I_2, J_2) in cal(B)^star (E)$, $I_1 subset.eq I_2 ==> J_2 = I_2^(arrow.t) subset.eq I_1^(arrow.t) = J_1$, and $J_2 subset.eq J_1 ==> I_1 = J_1^(arrow.b) subset.eq J_2^(arrow.b) = I_2$, establishing the order anti-isomorphism.
]

#subsec("B", "Soundness, Completeness, and Zero-Duplicate Guarantee")

The core algorithmic challenge in Bimax is proving that local base-case leaves $(U, V)$ are _globally_ inclusion-maximal in the full matrix $E$ without any post-hoc maximality check, and that no maximal bicluster is emitted twice across overlapping $G_W$ rows.

#theorem(title: "Exactness and Uniqueness of Bimax Enumeration", origin: "Classical — Cited")[
  For any binary matrix $E in {0, 1}^(n times m)$ with $r_"min" = c_"min" = 1$, `Conquer`$(R, C, emptyset)$ @prelic2006bimax terminates and returns every bicluster in $cal(B)^star (E)$ exactly once, emitting no non-maximal submatrices and no duplicates.
]

#proof[
  We first establish two structural loop invariants maintained at entry to every recursive call $(U, V, cal(Z))$ after Step 1:
  - *(I1 — Row Exclusion):* Every excluded row $i in R without U$ either satisfies $C(i) inter V = emptyset$ or there exists $Z in cal(Z)$ such that $C(i) inter V inter Z = emptyset$.
  - *(I2 — Column Exclusion):* Every excluded column $j in C without V$ either satisfies $R(j) inter U = emptyset$ or there exists an _anchor row_ $r in U$ that is all-ones on $V$ ($V subset.eq C(r)$) and satisfies $E_(r j) = 0$.

  At the root $(R, C, emptyset)$, $R without U$ contains only all-zero rows and $C without V$ contains only all-zero columns, so (I1) and (I2) hold. Assume both hold at a non-uniform call $(U, V, cal(Z))$ with pivot $r^star in U$. Because $C_U = C(r^star) inter V$ and $C_V = V without C_U$ are both non-empty, $r^star in G_U$ and $r^star$ is all-ones on $C_U$ with $E_(r^star j) = 0$ for all $j in C_V$. Crucially, observe the *anchor persistence property*: any existing anchor row $r in U$ satisfying $V subset.eq C(r)$ has $1$-entries in both $C_U$ and $C_V$, hence necessarily lies in $G_W$, and therefore survives into _both_ $G_U union G_W$ and $G_W union G_V$ (and remains all-ones on any subset of $V$).

  *Left child (Branch A):* Rows removed beyond $R without U$ are $G_V$ (whose rows satisfy $C(i) inter C_U = emptyset$) and rows failing Step 1 against $cal(Z)_A$, preserving (I1). Columns removed beyond $C without V$ are $C_V$ and all-zero columns on $U_A$. For every $j in C_V$, $E_(r^star j) = 0$; since $r^star$ is all-ones on $C_U$ and every $Z' in cal(Z)_A$ is non-empty, $C(r^star) inter C_U inter Z' = Z' != emptyset$, so $r^star$ survives Step 1 into $U_A$ as an anchor row covering all $j in C_V$. All prior anchors stay in $G_W inter U_A$ (also all-ones on $C_U$), preserving (I2).

  *Right child (Branch B):* Rows removed beyond $R without U$ are $G_U$ (and rows dropped in Step 1). Every $i in G_U$ satisfies $C(i) inter C_V = emptyset$. If $G_W != emptyset$, $C_V$ is added to $cal(Z)_B$, so $C(i) inter V_B inter C_V = emptyset$; if $G_W = emptyset$, $V_B = C_V$, so $C(i) inter V_B = emptyset$. Thus (I1) holds. Columns removed in Branch B are only all-zero columns on $U_B$ (including $C_U$ when $G_W = emptyset$), while every prior anchor row in $G_W$ satisfies $C(r) inter V = V$, hence meets $C_V$ and all $Z in cal(Z)$ and survives into $U_B$, preserving (I2).

  *Soundness:* Suppose Step 2 returns $(U, V)$ with $E_(U, V) = bold(1)$, $U != emptyset$, $V != emptyset$. By Step 1, $U != emptyset$ ensures $V inter Z != emptyset$ for all $Z in cal(Z)$. For any $i in R without U$, (I1) gives either $C(i) inter V = emptyset$ or $C(i) inter (V inter Z) = emptyset$ for some $Z in cal(Z)$; since $V$ and $V inter Z$ are non-empty, there exists $j in V$ with $E_(i j) = 0$, so $i$ cannot be added. For any $j in C without V$, (I2) gives either $R(j) inter U = emptyset$ or some $r in U$ with $E_(r j) = 0$, so $j$ cannot be added. Hence $(U, V) in cal(B)^star (E)$.

  *Completeness and Uniqueness:* Let $(I^star, J^star) in cal(B)^star (E)$. At the root, $I^star subset.eq U$, $J^star subset.eq V$, and $J^star inter Z != emptyset$ for all $Z in cal(Z)$ trivially. Inductively, at a split node $(U, V, cal(Z))$ satisfying these conditions, every $i in I^star$ has $J^star subset.eq C(i) inter V$, so $I^star$ survives Step 1. Exactly one of two mutually exclusive cases holds:
  1. $J^star inter C_V = emptyset$: Then $J^star subset.eq C_U$, so $(I^star, J^star)$ cannot satisfy $cal(Z)_B$ (and if $G_W = emptyset$, $J^star inter V_B = emptyset$). Meanwhile every $i in I^star$ hits $C_U$, so $I^star subset.eq G_U union G_W$, and $J^star inter (Z inter C_U) = J^star inter Z != emptyset$, routing $(I^star, J^star)$ uniquely to Branch A.
  2. $J^star inter C_V != emptyset$: Then $J^star subset.eq.not C_U$, excluding Branch A. Every $i in I^star$ hits $C_V$, so $I^star subset.eq G_W union G_V$, and $J^star$ satisfies $cal(Z)_B$ (if $G_W = emptyset$, $I^star subset.eq G_V$ implies $J^star subset.eq C_V = V_B$), routing $(I^star, J^star)$ uniquely to Branch B.

  Finally, termination holds because $r^star$ becomes all-ones on $C_U$ in Branch A (ineligible as a future pivot) and is deleted in Branch B, so the number of eligible pivot rows strictly decreases at every recursive step. At the unique leaf $(U, V)$ reached by $(I^star, J^star)$, $I^star subset.eq U$, $J^star subset.eq V$, and $E_(U, V) = bold(1)$ force $(I^star, J^star) = (U, V)$ by inclusion-maximality.
]

#subsec("C", "Tight Extremal Multiplicity and Output-Sensitive Complexity")

Because Bimax enumerates all maximal biclusters, its runtime depends critically on the extremal size of $cal(B)^star (E)$ as a function of $(n, m)$. We derive exact, attainable upper bounds on $beta(E)$ and then establish the output-sensitive time complexity.

#proposition(title: "Exact Extremal Bounds on Maximal Bicluster Multiplicity", origin: "Original Derivation")[
  For every binary matrix $E in {0, 1}^(n times m)$ with $n, m >= 1$, let $q = min(n, m)$. Then:
  1. $beta(E) <= 2^q - 1$, and this bound is attained with equality whenever $max(n, m) >= q + 1$.
  2. For square matrices ($n = m >= 2$), the sharp upper bound tightens to $beta(E) <= 2^n - 2$, attained by the complement of the identity matrix $E = bold(1)_n bold(1)_n^top - I_n$.
]

#proof[
  Without loss of generality assume $n <= m$, so $q = n$ (otherwise transpose $E$, which preserves $beta(E)$). By Lemma 1, the map $(I, J) |-> I$ from $cal(B)^star (E)$ into $2^R without {emptyset}$ is injective because $J = I^(arrow.t)$ is uniquely determined by $I$. Since $|2^R without {emptyset}| = 2^n - 1$, we immediately obtain $beta(E) <= 2^n - 1$.

  To prove tightness when $m >= n + 1$, construct $E in {0, 1}^(n times m)$ whose first column is $bold(1)_n$, columns $2, dots, n+1$ equal $bold(1)_n bold(1)_n^top - I_n$, and any remaining columns are all-zero. For any non-empty $I subset.eq {1, dots, n}$, column $1$ belongs to $I^(arrow.t)$, and for $k in {1, dots, n}$, column $k + 1$ belongs to $I^(arrow.t)$ iff $k in.not I$. Hence $I^(arrow.t) = {1} union {k + 1 : k in R without I} != emptyset$. Applying $(dot.c)^(arrow.b)$, row $i$ belongs to $(I^(arrow.t))^(arrow.b)$ iff $i in.not R without I$, i.e., $i in I$. Thus $(I^(arrow.t))^(arrow.b) = I$ for all $2^n - 1$ non-empty subsets $I$, giving $beta(E) = 2^n - 1$.

  Now suppose $n = m >= 2$ and assume for contradiction that $beta(E) = 2^n - 1$. Then every non-empty $I subset.eq R$ is closed with $I^(arrow.t) != emptyset$; in particular, $R^(arrow.t) != emptyset$, so at least one column $j_0 in C$ is identically $1$. Consequently, every $(I, J) in cal(B)^star (E)$ must satisfy $j_0 in I^(arrow.t) = J$. Because $(I, J) |-> J$ is also injective on $cal(B)^star (E)$, $beta(E)$ is bounded by the number of subsets of $C$ containing $j_0$, namely $2^(m - 1) = 2^(n - 1)$. However, $2^n - 1 > 2^(n - 1)$ for all $n >= 2$, a contradiction. Thus $beta(E) <= 2^n - 2$. For $E = bold(1)_n bold(1)_n^top - I_n$, every non-empty proper subset $I subset R$ satisfies $I^(arrow.t) = R without I != emptyset$ and $(R without I)^(arrow.b) = I$, whereas $R^(arrow.t) = emptyset$, yielding exactly $2^n - 2$ maximal biclusters.
]

#theorem(title: "Output-Sensitive Time Complexity of Bimax", origin: "Classical — Cited")[
  Assuming $E$ is oriented so that $n <= m$ (transposing in $O(n m)$ time if $n > m$), let $beta = beta(E)$ be the output cardinality. With Step 1 row-feasibility pruning, the recursion tree of `Conquer` has depth at most $n$ and at most $beta n + 1$ active calls, yielding total runtime $O(n m beta min(n, m))$ when mandatory-set tests run in $O(n m)$ per call (e.g., via $w$-bit word-parallel masks for $|cal(Z)| <= w$, or $O(n m beta min(n, m)^2)$ in the scalar RAM model) @prelic2006bimax. When all maximal biclusters are pairwise disjoint ($G_W = emptyset$ at every split), runtime reduces to $O(n m beta)$.
]

#proof[
  At every internal split node, the pivot row $r^star in U$ becomes all-ones on $C_U$ in Branch A (and remains all-ones on all column subsets of $C_U$ in every descendant of Branch A) and is excluded from $U_B = G_W union G_V$ in Branch B. Hence $r^star$ can never be chosen as a pivot again in either subtree. Since there are at most $n = min(n, m)$ rows initially, every root-to-leaf path contains at most $n$ splits, bounding recursion depth by $min(n, m)$.

  Next, we show that zero dead subtrees are traversed after Step 1 pruning. If $U != emptyset$ after Step 1, pick any surviving row $i_0 in U$ and let $J_0 = C(i_0) inter V$. By Step 1, $J_0 != emptyset$ and $J_0 inter Z != emptyset$ for all $Z in cal(Z)$. Let $I_0 = J_0^(arrow.b) inter U$ (which contains $i_0$). Then $(I_0, J_0)$ is a valid all-ones submatrix in $(U, V, cal(Z))$, so extending it to a maximal bicluster $(I^star, J^star) in cal(B)^star (E)$ satisfies $I^star subset.eq U$, $J^star subset.eq V$, and $J^star inter Z != emptyset$ for all $Z in cal(Z)$ (by invariants (I1)–(I2)). By Theorem 1, this subtree must emit $(I^star, J^star)$. Thus every active recursive call is an ancestor of at least one of the $beta$ output leaves. Since each leaf has at most $n$ ancestors, the total number of active calls is at most $beta n + 1 = O(beta min(n, m))$.

  At each call, partitioning $V$ into $(C_U, C_V)$ and $U$ into $(G_U, G_V, G_W)$ scans $U times V$ in $O(n m)$ time. In Branch B, $V_B = V$ does not shrink old mandatory sets, and the sole new constraint $C_V$ is already met by every row in $G_W union G_V$ in $O(1)$ extra work. In Branch A, updating $cal(Z)_A$ takes $O(m |cal(Z)|) <= O(n m)$ time, while testing $C(i) inter (Z inter C_U) != emptyset$ across $i in U_A$ and $Z in cal(Z)_A$ takes $O(n m ceil(|cal(Z)| / w))$ bitset operations ($O(n m)$ for $|cal(Z)| <= w$, or $O(n m min(n, m))$ scalar operations). Summing over $O(beta min(n, m))$ active calls gives $O(n m beta min(n, m))$ @prelic2006bimax. When $G_W = emptyset$ everywhere, $cal(Z) = emptyset$ identically, every internal split produces two non-empty children so the binary tree has at most $2 beta - 1$ nodes total, and $O(n m)$ work per node gives $O(n m beta)$.
]

#subsec("D", "Exponential Noise Sensitivity and Breakdown Analysis")

Although Bimax provides exact global enumeration on noise-free binary data, its hard all-ones requirement induces severe structural fragility under sparse bit-dropout noise.

#proposition(title: "Exponential Bicluster Shattering Under Sparse Bit Dropouts", origin: "Original Counterexample")[
  Let $E_0 = bold(1)_s bold(1)_s^top in {0, 1}^(s times s)$ be a single planted $s times s$ all-ones bicluster ($beta(E_0) = 1$). Suppose $k$ independent bit dropouts ($1 <= k < s$) flip the diagonal entries $(1, 1), dots, (k, k)$ from $1$ to $0$, yielding $E_k$. Then:
  1. $E_k$ shatters the single planted bicluster into _exactly_ $beta(E_k) = 2^k$ distinct inclusion-maximal biclusters, parameterized bijectively by subsets $A subset.eq {1, dots, k}$ as $I_A = A union T$, $J_A = (S without A) union T$, where $S = {1, dots, k}$ and $T = {k + 1, dots, s}$.
  2. Consequently, a bit-flip fraction as small as $k / s^2 = O(1 / s)$ reduces the maximum Jaccard area recovery of the ground-truth block by a factor of $(s - k) / s$ while multiplying output cardinality exponentially by $2^k$.
]

#proof[
  Let $S = {1, dots, k}$ and $T = {k + 1, dots, s}$ with $T != emptyset$ since $k < s$. For any $A subset.eq S$, set $I_A = A union T$ and $J_A = (S without A) union T$. For every $(i, j) in I_A times J_A$, if $i in S$ and $j in S$ then $i in A$ and $j in S without A$, so $i != j$ and $(E_k)_(i j) = 1$; if $i in T$ or $j in T$ then $(E_k)_(i j) = 1$ by construction. Thus $(I_A, J_A)$ is a valid all-ones bicluster of size $(s - k + |A|) times (s - |A|)$. No row $r in S without A$ can be added to $I_A$ because $r in J_A$ and $(E_k)_(r r) = 0$; symmetrically, no column $c in A$ can be added to $J_A$ because $c in I_A$ and $(E_k)_(c c) = 0$. Hence $(I_A, J_A) in cal(B)^star (E_k)$ for all $2^k$ subsets $A subset.eq S$.

  Conversely, let $(I, J) in cal(B)^star (E_k)$. Because every row $i in T$ is all-ones across $C$ and every column $j in T$ is all-ones across $R$, inclusion-maximality enforces $T subset.eq I$ and $T subset.eq J$. Let $A = I inter S$. For every $r in A$, $(E_k)_(r r) = 0$ forbids $r in J$, so $J inter S subset.eq S without A$. Since $(I, (S without A) union T)$ is valid, maximality forces $J inter S = S without A$, hence $(I, J) = (I_A, J_A)$. Thus $beta(E_k) = 2^k$. Each fragment has area $|I_A| |J_A| = (s - k + |A|)(s - |A|)$, maximized at $|A| = 0$ or $|A| = k$ with area $s(s - k)$, giving maximum Jaccard overlap $s(s - k) / s^2 = (s - k) / s$ relative to $E_0$, despite only $k$ corrupted entries out of $s^2$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Bimax's operational profile is dominated by the tension between exactness and output cardinality $beta(E)$. In dense matrices with scattered zeros, Proposition 2 shows that $beta(E)$ can grow as $2^(min(n, m))$, turning output-sensitive enumeration intractable unless constrained by minimum row/column support thresholds $(r_"min", c_"min")$ or an explicit cap on retained biclusters. Conversely, on sparse binary incidence matrices (edge density $p << 1$) or datasets with a modest number of clean co-active modules, the recursion depth collapses rapidly: Step 1 discards rows outside $C_U$ or $C_V$, shrinking submatrices geometrically at each split.

Memory complexity is remarkably lean: because `Conquer` is a depth-first recursion of depth at most $min(n, m)$, storing only the active row/column index masks and mandatory sets $cal(Z)$ along the current root-to-leaf stack requires $O(n m + min(n, m)(n + m))$ working space — strictly polynomial in $(n, m)$ even when $beta(E)$ is exponential (provided biclusters are streamed to disk or filtered on the fly). Bit-packing rows into $ceil(m / 64)$ 64-bit words converts set intersections $C(i) inter C_U$, $C(i) inter C_V$, and $C(i) inter Z$ into cache-friendly SIMD/AVX bitwise `AND` instructions, delivering up to $64 times$ practical throughput gains over byte arrays and fitting L1/L2 caches cleanly on moderate matrices. Geometrically, every Bimax bicluster corresponds to an axis-aligned hyper-rectangle (an exact Cartesian product $I times J$ of binary indicators) in the discretized feature space, affording zero invariance to arbitrary rotations or continuous feature rescaling after binarization, while being strictly invariant to arbitrary row permutations and column permutations.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Input Domain], [Binary $E in {0, 1}^(n times m)$ (real $X$ via threshold $tau$)],
    [Time Complexity], [$O(n m beta min(n, m))$ general; $O(n m beta)$ disjoint],
    [Working Space], [$O(n m + (n + m) min(n, m))$ stack (streamed)],
    [Optimality], [Exact global enumeration of $cal(B)^star (E)$; zero duplicates],
    [Overlap Structure], [Arbitrary row and column overlap permitted],
    [Invariances], [Row/column permutation invariant; monotone $tau$-equivariant],
    [Noise Breakdown], [$k$ diagonal dropouts split 1 block into $2^k$ fragments],
    [Key Hyperparameters], [Binarization threshold $tau$; minimum sizes $(r_"min", c_"min")$],
  ),
  caption: [Operational and complexity profile of Bimax @prelic2006bimax.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Differential Gene Expression and Functional Genomics:* In microarray and scRNA-seq perturbation screens, genes are discretized into differentially expressed ($1$) versus baseline ($0$) states across drug or knockout conditions @prelic2006bimax. Bimax's exact inclusion-maximal formulation guarantees that every complete co-activated transcription module is discovered without local-optimum misses, and permits genes belonging to multiple biological pathways to participate in overlapping biclusters simultaneously.
2. *Market Basket Analysis and Closed Frequent Itemset Mining:* In retail transaction matrices ($E_(i j) = 1$ iff customer $i$ purchased item $j$), Lemma 1 establishes that Bimax's maximal biclusters $(I, J)$ with $|I| >= r_"min"$ are identically the *closed frequent itemsets* together with their exact supporting transaction lists. Exact enumeration eliminates redundant sub-itemsets that share identical support.
3. *Role-Based Access Control (RBAC) and Permission Mining:* Given a user-permission binary matrix $E$, inclusion-maximal all-ones submatrices directly correspond to candidate security roles (maximal sets of permissions jointly held by a maximal user group). Bimax's completeness ensures no valid maximal role is overlooked during policy refactoring and least-privilege auditing.
4. *Bipartite Network and Ecological Interaction Analysis:* In plant-pollinator, host-pathogen, or user-tag bipartite graphs $G = (R union C, E)$, Bimax enumerates all maximal bicliques (complete bipartite subgraphs), exposing tightly coupled functional guilds and structural motifs with zero heuristic approximation error.
5. *Combinatorial Benchmark Reference for Biclustering Suites:* As established in @prelic2006bimax, Bimax serves as an exact deterministic baseline for calibrating stochastic or continuous biclustering algorithms (such as ISA @bergmann2003isa, Cheng-Church @cheng2000biclustering, and FABIA @hochreiter2010fabia) on synthetic constant-expression benchmarks where exact recall can be verified analytically.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Bimax encounters four principal failure regimes: (i) *dropout/noise fragmentation* (Proposition 2), where even $1"–"5%$ false-negative zeros shatter true modules into exponential families of partial submatrices; (ii) *binarization information loss*, where discretizing graded continuous patterns (e.g., linear scaling or coherent multiplicative slopes across conditions) destroys continuous rank structure; (iii) *combinatorial explosion* on dense, weakly structured binary matrices ($p approx 0.5"–"0.8$), where $beta(E)$ becomes infeasible unless aggressive $(r_"min", c_"min")$ pruning is enforced; and (iv) *lack of a global hard partition*, since overlapping biclusters can leave noisy objects unassigned and assign core objects to dozens of modules.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Binary Differential Omics], [Exact maximal co-activation; natural gene overlap], [Sensitive to dropouts; try FABIA @hochreiter2010fabia],
    [Closed Itemset / RBAC Mining], [Exact Galois concept bijection; zero redundancy], [Cap $(r_"min", c_"min")$ on dense tables],
    [Noisy Continuous Expression], [Weak; hard binarization discards amplitude ratios], [Use ISA @bergmann2003isa or Plaid @lazzeroni2002plaid],
    [Low-D Continuous Geometry], [Poor; axis-aligned binary cuts distort manifolds], [Use spectral or density clustering],
    [Bipartite Motif Enumeration], [Exact maximal biclique recovery in $O(n m beta min)$], [Exponential $beta$ when density $p > 0.5$],
  ),
  caption: [Domain suitability and failure-mode matrix for Bimax.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Bimax occupies a foundational combinatorial bridge linking formal concept analysis, bipartite graph theory, and axis-aligned subspace clustering @kriegel2009subspacesurvey.

#proposition(title: "Exact Reduction of Bimax to Maximal Bicliques and Binary CLIQUE", origin: "Original Proof")[
  Let $E in {0, 1}^(n times m)$ induce the bipartite graph $G_E = (R union.sq C, cal(E))$ with $(i, j) in cal(E) <==> E_(i j) = 1$, and let $cal(D)_E = {x_1, dots, x_n} subset {0, 1}^m$ be the rows of $E$ embedded in the $m$-dimensional binary hypercube under CLIQUE @agrawal1998clique with $xi = 2$ grid bins per coordinate and density threshold $tau_"den" = r_"min" / n$. Then:
  1. $(I, J) in cal(B)^star (E)$ iff $I union.sq J$ induces a non-trivial inclusion-maximal biclique (complete bipartite subgraph) in $G_E$.
  2. Every size-filtered Bimax bicluster $(I, J) in cal(B)^star (E)$ with $|I| >= r_"min"$ corresponds bijectively to a maximal dense $1$-anchored CLIQUE subspace unit $U_J = {x in {0, 1}^m : x_j = 1 "for all" j in J}$ whose point support $I = {i : x_i in U_J}$ cannot be preserved while adding any further coordinate $j' in.not J$.
]

#proof[
  For (1), a subset pair $(I, J)$ with $I subset.eq R$, $J subset.eq C$, $I, J != emptyset$ forms a biclique in $G_E$ iff every $i in I$ is adjacent to every $j in J$, i.e., $E_(i j) = 1$ on $I times J$. Maximality of the biclique under vertex inclusion on both sides coincides verbatim with Definition 1. For (2), under $2$-bin grid partitioning $\{0\}, \{1\}$ per feature, fixing bin $1$ along coordinates $J subset.eq C$ defines the subspace grid cell $U_J$ containing exactly $I = inter.big_(j in J) R(j) = J^(arrow.b)$ points. Its empirical selectivity satisfies $|I| / n >= r_"min" / n <==> |I| >= r_"min"$. By Lemma 1, $(I, J) in cal(B)^star (E)$ iff $I = J^(arrow.b)$ and $J = I^(arrow.t)$, which holds iff no $j' in.not J$ satisfies $I subset.eq R(j')$, i.e., iff $U_J$ is dimension-maximal among all grid cells sharing the same support $I$.
]

Relative to neighboring Cat. VII algorithms, Bimax contrasts sharply with continuous iterative thresholding in ISA (Alg. 160 @bergmann2003isa) and generative sparse factor analysis in FABIA (Alg. 161 @hochreiter2010fabia): Bimax sacrifices noise tolerance and multiplicative amplitude modeling in exchange for deterministic completeness and zero local-optima variance.

*Empirical Behavior on ClusterOverview Benchmarks:* On the standardized four-dataset ClusterOverview suite evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, Bimax attains a mean runtime of *0.215 ms* — the fastest among Algorithms 159–162 — reflecting the compact recursion depth $min(n, m)$ and word-fast binary partitioning on moderate matrices. On *Iris* (`Sil = 0.193`, `ARI = 0.420`), thresholding the four morphological dimensions cleanly isolates the well-separated _Iris setosa_ block (uniformly low petal length/width and high sepal width) into a crisp all-ones submatrix, achieving solid `ARI = 0.420` despite modest Euclidean Silhouette (`0.193`) caused by axis-aligned binarization cutting through the continuous _versicolor_–_virginica_ boundary. On *Customers* (`Sil = N/A`) and *Multi-Shape* (`Sil = N/A`, `ARI = N/A`), strict all-ones extraction fails to induce a valid multi-cluster hard partition across all rows (collapsing into a single dominant block or leaving unclustered / multiply-overlapping points), illustrating the structural mismatch between overlapping binary biclusters and global partition metrics. On *High-D* (`Sil = 0.095`, `ARI = N/A`), Bimax isolates genuine localized feature blocks spanning subsets of active dimensions; full-space Euclidean Silhouette remains low (`0.095`) because bicluster coherence holds strictly on the selected column subset $J subset C$ while distances over irrelevant coordinates $C without J$ remain unminimized.
