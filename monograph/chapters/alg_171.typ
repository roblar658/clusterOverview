#import "../template/ieee.typ": *

#alg_entry(
  id: 171,
  name: "CSPA (Cluster-based Similarity Partitioning)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2002",
  paradigm: "Co-Association Similarity Graph Partitioning",
  complexity: $O(n^2 (r + k))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $V = {1, dots, n}$ index a set of $n$ objects without requiring direct access to their original feature vectors. Instead, we are given a _cluster ensemble_ $Lambda = {lambda^((1)), dots, lambda^((r))}$ consisting of $r$ base partitions, where each labeling $lambda^((q)) : V -> {1, dots, k^((q))}$ assigns object $i$ to one of $k^((q))$ clusters (possibly omitting a subset of objects in distributed or partial ensembles). The _cluster ensemble problem_ formulated by Strehl and Ghosh @strehl2002cluster seeks a single consensus $k$-partition $lambda^* : V -> {1, dots, k}$ that shares maximal statistical information with the ensemble without accessing raw features or assuming label correspondence across base clusterers.

Formally, Strehl and Ghosh define the consensus quality criterion via the _Average Normalized Mutual Information_ ($"ANMI"$):
$ Phi^("ANMI") (Lambda, lambda) = 1 / r sum_(q=1)^r "NMI"(lambda, lambda^((q))), $
where the symmetric geometric-mean Normalized Mutual Information between two partitions $lambda^((a)), lambda^((b))$ with contingency counts $n_(h l) = abs({i : lambda^((a)) (i) = h and lambda^((b)) (i) = l})$ and marginals $n_h^((a)), n_l^((b))$ is
$ "NMI"(lambda^((a)), lambda^((b))) = (sum_(h, l) n_(h l) log((n dot n_(h l)) / (n_h^((a)) n_l^((b))))) / sqrt((sum_h n_h^((a)) log(n_h^((a)) / n))(sum_l n_l^((b)) log(n_l^((b)) / n))). $
Direct maximization of $Phi^("ANMI")$ over all $k$-partitions is combinatorially intractable. To obtain a principled, label-permutation-invariant surrogate, each base partition $lambda^((q))$ is encoded as a binary cluster-membership indicator matrix $H^((q)) in {0, 1}^(n times k^((q)))$ with $H^((q))_(i c) = bb(1)(lambda^((q)) (i) = c)$, and concatenated into the _hypergraph incidence matrix_
$ H = [ H^((1)), H^((2)), dots, H^((r)) ] in {0, 1}^(n times m), quad m = sum_(q=1)^r k^((q)). $

#definition(title: "CSPA Co-Association Similarity Matrix and Balanced Cut Objective")[
  The _Cluster-based Similarity Partitioning Algorithm (CSPA)_ @strehl2002cluster induces an $n times n$ pairwise co-association similarity matrix
  $ S = 1 / r H H^top in [0, 1]^(n times n), quad S_(i j) = 1 / r sum_(q=1)^r bb(1)(lambda^((q)) (i) = lambda^((q)) (j)). $
  Constructing the weighted undirected graph $G = (V, E, W)$ with edge weights $W_(i j) = S_(i j)$ for $i != j$ (and $W_(i i) = 0$), CSPA extracts a consensus $k$-partition $cal(C) = {C_1, dots, C_k}$ by minimizing the $k$-way graph cut
  $ "Cut"(cal(C); S) = 1 / 2 sum_(c=1)^k sum_(i in C_c) sum_(j in V without C_c) S_(i j), $
  subject to the cluster-size balance constraint $max_c abs(C_c) <= (1 + epsilon_"bal") n / k$ (solved in practice via multilevel $k$-way graph partitioning, METIS @karypis1998metis, or normalized spectral clustering @shi2000ncut).
]

Intuitively, every base clustering defines a coarse binary equivalence relation; averaging $r$ such relations constructs a refined similarity kernel $S$ whose entry $S_(i j)$ equals the empirical fraction of ensemble members that co-cluster objects $i$ and $j$. Re-clustering $G$ via graph partitioning then groups objects with high mutual co-association while discarding inconsistent boundary disagreements.

#algbox(
  title: "CSPA (Cluster-based Similarity Partitioning Algorithm)",
  inputs: [Cluster ensemble $Lambda = {lambda^((1)), dots, lambda^((r))}$ on $n$ objects, target cluster count $k$, balance tolerance $epsilon_"bal" >= 0$.],
  outputs: [Co-association matrix $S in [0, 1]^(n times n)$ and consensus partition $lambda^* : V -> {1, dots, k}$.],
)[
  + *Hypergraph encoding*: for each $q = 1, dots, r$, construct the binary membership matrix $H^((q)) in {0, 1}^(n times k^((q)))$ with $H^((q))_(i c) = bb(1)(lambda^((q)) (i) = c)$, and form $H = [H^((1)), dots, H^((r))] in {0, 1}^(n times m)$.
  + *Co-association construction*: compute the dense $n times n$ similarity matrix $S <- (1 / r) H H^top$ (equivalently accumulating $1 / r$ over all pairs within each cluster $C^((q))_c$).
  + *Similarity graph assembly*: form the weighted undirected graph $G = (V, E, W)$ with $W_(i j) = S_(i j)$ for $i != j$ (optionally zeroing self-loops and discarding zero-weight edges).
  + *Multi-way graph partitioning*: invoke multilevel $k$-way graph partitioning (METIS @karypis1998metis) on $G$ to minimize $"Cut"(cal(C); S)$ subject to $abs(C_c) <= (1 + epsilon_"bal") n / k$.
  + *Consensus output*: return the resulting $k$-way labeling $lambda^*$ and (when used inside a supra-consensus portfolio) evaluate $Phi^("ANMI") (Lambda, lambda^*)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal results characterizing the reproducing-kernel geometry of $S$, exact equivalence between CSPA cuts, Mirkin median partitioning, and kernel $K$-means, a constructive imbalance failure mode, a Davis–Kahan spectral recovery bound, and invariance properties.

#subsec("A", "Kernel Structure, Rank Bound, and Euclidean Isometric Embedding")

#theorem(title: "Gram-Kernel Geometry and Exact Rank Bound of CSPA Similarity", origin: "Original Proof")[
  Let $Lambda$ be a complete ensemble ($sum_c H^((q))_(i c) = 1$ for all $i, q$) with $m = sum_(q=1)^r k^((q))$ total base clusters, and let $S = (1 / r) H H^top$. Then:
  (a) $S$ is a positive semidefinite Mercer kernel on $V$ with unit diagonal $S_(i i) = 1$ and entries in $[0, 1]$.
  (b) The rank of $S$ obeys the tight linear bound
  $ "rank"(S) <= min(n, m - r + 1). $
  (c) The feature map $phi(i) = (1 / sqrt(r)) H_(i, :)^top in bb(R)^m$ embeds $V$ on the unit hypersphere $bb(S)^(m-1)$ such that $d_E (i, j) = sqrt(2 (1 - S_(i j)))$ is a Euclidean metric and $d_M (i, j) = 1 - S_(i j)$ is a pseudometric of negative type.
]
#proof[
  *Part (a):* For any $v in bb(R)^n$, $v^top S v = (1 / r) norm(H^top v)_2^2 >= 0$, proving $S succ.eq 0$. Since each complete base clustering assigns object $i$ to exactly one cluster, row $H_(i, :)$ contains exactly $r$ ones and $m - r$ zeros, giving $S_(i i) = (1 / r) norm(H_(i, :))_2^2 = 1$, and by Cauchy–Schwarz $0 <= S_(i j) <= sqrt(S_(i i) S_(j j)) = 1$.

  *Part (b):* Because each block $H^((q))$ partitions $V$, its columns sum to the all-ones vector: $H^((q)) bold(1)_(k^((q))) = bold(1)_n$ for every $q = 1, dots, r$. Define $r - 1$ vectors $u^((2)), dots, u^((r)) in bb(R)^m$ where $u^((q))$ places $+1$ on the indices of block $H^((q))$, $-1$ on the indices of block $H^((1))$, and $0$ elsewhere. Then $H u^((q)) = bold(1)_n - bold(1)_n = bold(0)$ for all $q = 2, dots, r$. Since the supports of the $+1$ blocks in $u^((2)), dots, u^((r))$ are pairwise disjoint, these $r - 1$ nullspace vectors are linearly independent. Hence $dim ker(H) >= r - 1$, and by rank–nullity,
  $ "rank"(S) = "rank"(H H^top) = "rank"(H) <= min(n, m - (r - 1)). $

  *Part (c):* By construction, $norm(phi(i))_2^2 = S_(i i) = 1$ and $chevron.l phi(i), phi(j) chevron.r = S_(i j)$. Expanding squared Euclidean distance in $bb(R)^m$ gives
  $ norm(phi(i) - phi(j))_2^2 = norm(phi(i))_2^2 + norm(phi(j))_2^2 - 2 chevron.l phi(i), phi(j) chevron.r = 2 (1 - S_(i j)). $
  Thus $sqrt(2 (1 - S_(i j)))$ inherits all Euclidean metric axioms. Moreover, for any vector $c in bb(R)^n$ with $sum_i c_i = 0$, $sum_(i, j) c_i c_j (1 - S_(i j)) = -c^top S c <= 0$, establishing that $1 - S_(i j)$ is a negative-type pseudometric (satisfying the triangle inequality as the average of equivalence pseudometrics $1 - M^((q))_(i j)$).
]

#subsec("B", "Equivalence to Mirkin Median Partitioning and Kernel K-Means")

Although CSPA is introduced heuristically as graph partitioning on $S$, we prove that its objective coincides algebraically with two classical optimization formulations: the Mirkin median partition and kernel $K$-means @dhillon2004kernel.

#theorem(title: "Exact Equivalence to Mirkin Median Partition and Kernel K-Means", origin: "Original Derivation")[
  Let $P in {0, 1}^(n times n)$ be the binary co-clustering matrix of a candidate $k$-partition $cal(C) = {C_1, dots, C_k}$, and let $M^((q)) = H^((q)) (H^((q)))^top$ be the $q$-th base co-clustering matrix. Define the ensemble-averaged Mirkin disagreement $J_"Mirkin" (cal(C); Lambda) = (1 / r) sum_(q=1)^r norm(P - M^((q)))_F^2$. Then:
  (a) For every partition $cal(C)$,
  $ J_"Mirkin" (cal(C); Lambda) = sum_(c=1)^k abs(C_c)^2 + 4 "Cut"(cal(C); S) - sum_(i, j=1)^n S_(i j). $
  In particular, over balanced $k$-partitions ($abs(C_c) = n / k$), minimizing the CSPA cut $"Cut"(cal(C); S)$ is *exact* Mirkin median partitioning.
  (b) Let $J_"KM" (cal(C); phi) = sum_(c=1)^k sum_(i in C_c) norm(phi(i) - mu_c)_2^2$ be the $K$-means objective in the CSPA feature embedding $phi(i) = (1 / sqrt(r)) H_(i, :)^top$. If $G$ is regular with uniform degree $d_i = sum_(j=1)^n S_(i j) = bar(d)$ (e.g., all base clusters have size $n / k$), then the Normalized Cut @shi2000ncut obeys the exact affine identity
  $ "Ncut"(cal(C); S) = 1 / bar(d) J_"KM" (cal(C); phi) + (k - n / bar(d)). $
]
#proof[
  *Part (a):* Because $P_(i j), M^((q))_(i j) in {0, 1}$, we have $P_(i j)^2 = P_(i j)$ and $(M^((q))_(i j))^2 = M^((q))_(i j)$. Expanding the squared Frobenius norm and averaging over $q = 1, dots, r$ yields
  $ J_"Mirkin" (cal(C); Lambda) & = 1 / r sum_(q=1)^r sum_(i, j=1)^n (P_(i j) + M^((q))_(i j) - 2 P_(i j) M^((q))_(i j)) \
    & = sum_(i, j=1)^n P_(i j) + sum_(i, j=1)^n S_(i j) - 2 sum_(i, j=1)^n P_(i j) S_(i j). $
  Observe that $sum_(i, j) P_(i j) = sum_(c=1)^k abs(C_c)^2$, and since $1 - P_(i j) = 1$ iff $i$ and $j$ lie in distinct clusters,
  $ 2 "Cut"(cal(C); S) = sum_(i, j=1)^n (1 - P_(i j)) S_(i j) = sum_(i, j=1)^n S_(i j) - sum_(i, j=1)^n P_(i j) S_(i j). $
  Solving for $sum_(i, j) P_(i j) S_(i j) = sum_(i, j) S_(i j) - 2 "Cut"(cal(C); S)$ and substituting into $J_"Mirkin"$ immediately produces
  $ J_"Mirkin" (cal(C); Lambda) = sum_(c=1)^k abs(C_c)^2 + 4 "Cut"(cal(C); S) - sum_(i, j=1)^n S_(i j). $
  Under strict balance $abs(C_c) = n / k$, the first term equals the constant $n^2 / k$ and the last term depends only on $Lambda$, proving exact objective equivalence.

  *Part (b):* For each cluster $C_c$ with centroid $mu_c = (1 / abs(C_c)) sum_(i in C_c) phi(i)$, expanding the within-cluster SSE and using $norm(phi(i))_2^2 = S_(i i) = 1$ gives
  $ J_"KM" (cal(C); phi) & = sum_(c=1)^k ( sum_(i in C_c) norm(phi(i))_2^2 - abs(C_c) norm(mu_c)_2^2 ) \
    & = n - sum_(c=1)^k 1 / abs(C_c) sum_(i, j in C_c) S_(i j). $
  When every vertex has degree $d_i = bar(d)$, cluster volume is $"vol"(C_c) = bar(d) abs(C_c)$, and $"cut"(C_c, V without C_c) = bar(d) abs(C_c) - sum_(i, j in C_c) S_(i j)$. Dividing by $"vol"(C_c)$ and summing over $c = 1, dots, k$ yields
  $ "Ncut"(cal(C); S) & = sum_(c=1)^k (1 - 1 / (bar(d) abs(C_c)) sum_(i, j in C_c) S_(i j)) \
    & = k - (n - J_"KM" (cal(C); phi)) / bar(d) = 1 / bar(d) J_"KM" (cal(C); phi) + (k - n / bar(d)). $
]

#subsec("C", "Constructive Imbalance Failure Under METIS Size Constraints")

Theorem 2(a) highlights a critical operational hazard: CSPA minimizes $"Cut"(cal(C); S)$ subject to balance constraints rather than directly minimizing $J_"Mirkin"$. When true clusters are strongly skewed in size, the balance constraint forces catastrophic splitting even when the ensemble is unanimous.

#proposition(title: "Unanimous-Ensemble Imbalance Failure of Balanced-Cut CSPA", origin: "Original Counterexample")[
  Consider an ensemble of $r >= 1$ identical (unanimous) bipartitions on $n$ objects with true clusters $C_1^*, C_2^*$ of sizes $(1 - alpha) n$ and $alpha n$, where $alpha in (0, 1/2)$. Under balanced $2$-way graph partitioning enforcing $max(abs(C_1), abs(C_2)) <= U = floor((1 + epsilon_"bal") n / 2)$, whenever $alpha < (1 - epsilon_"bal") / 2$, EVERY feasible CSPA solution suffers strictly positive cut
  $ "Cut"(cal(C); S) >= U ((1 - alpha) n - U) = Omega(n^2) $
  and misclassifies at least $(1 - alpha) n - U$ objects, even though the unanimous true partition attains zero cut.
]
#proof[
  Because all $r$ base clusterings agree with ${C_1^*, C_2^*}$, $S_(i j) = 1$ for every intra-cluster pair and $S_(i j) = 0$ across $C_1^*$ and $C_2^*$. In any feasible bipartition ${C_1, C_2}$, let $a = abs(C_1 inter C_1^*)$ and $b = abs(C_2 inter C_1^*) = (1 - alpha) n - a$. Feasibility requires $a <= abs(C_1) <= U$ and $b <= abs(C_2) <= U$. Since $alpha < (1 - epsilon_"bal") / 2$, we have $U < (1 - alpha) n$, forcing both $a > 0$ and $b >= (1 - alpha) n - U > 0$. Because every pair in $(C_1 inter C_1^*) times (C_2 inter C_1^*)$ carries weight $S_(i j) = 1$, the cut satisfies
  $ "Cut"(cal(C); S) >= a ((1 - alpha) n - a). $
  The concave quadratic $g(a) = a((1 - alpha) n - a)$ on the interval $[(1 - alpha) n - U, U]$ attains its minimum at the endpoints $a = U$ and $a = (1 - alpha) n - U$, yielding $"Cut"(cal(C); S) >= U ((1 - alpha) n - U)$. For instance, at $alpha = 0.1$ and standard METIS tolerance $epsilon_"bal" = 0.03$ ($U = 0.515 n$), at least $0.385 n$ objects ($38.5\%$ of the dataset) are forcibly detached from $C_1^*$ into the minority cluster despite $100\%$ ensemble agreement.
]

#subsec("D", "Spectral Subspace Recovery Rate for Noisy Ensembles")

When CSPA is coupled with spectral partitioning @ng2002spectral @vonluxburg2007tutorial on $S$, we derive an exact non-asymptotic recovery bound showing that the subspace error decays at the parametric rate $O(k / ((p - q) sqrt(r)))$, independent of $n$.

#theorem(title: "Dimension-Free Spectral Recovery of Planted Ensemble Partitions", origin: "Original Derivation")[
  Let $cal(C)^* = {C_1^*, dots, C_k^*}$ be a balanced planted $k$-partition ($abs(C_c^*) = n / k$). Suppose $Lambda$ consists of $r$ independent noisy base partitions with population co-association matrix
  $ S^* = bb(E)[S] = (p - q) P^* + q bold(1) bold(1)^top + (1 - p) I_n, $
  where $P^*$ is the block-diagonal indicator matrix of $cal(C)^*$ and $1 >= p > q >= 0$. Let $U^*, hat(U) in bb(R)^(n times k)$ have orthonormal columns spanning the top-$k$ eigenspaces of $S^*$ and $S$, respectively. For any failure probability $eta in (0, 1)$, if the ensemble size satisfies $r > k^2 / (eta (p - q)^2)$, then with probability at least $1 - eta$,
  $ norm(sin Theta(hat(U), U^*))_F <= k / ((p - q) sqrt(r eta)). $
]
#proof[
  *Step 1 (Population spectrum):* The indicator matrix $P^*$ has eigenvalue $n / k$ with multiplicity $k$ (eigenspace $cal(V)_k = "span"{bold(1)_(C_1^*), dots, bold(1)_(C_k^*)}$) and eigenvalue $0$ with multiplicity $n - k$. Since $bold(1) = sum_(c=1)^k bold(1)_(C_c^*) in cal(V)_k$, the matrix $q bold(1) bold(1)^top$ acts as $q n$ along $bold(1) / sqrt(n)$ and $0$ on $bold(1)^perp$. Therefore the exact eigenvalues of $S^*$ are:
  - $lambda_1 (S^*) = (p - q)(n / k) + q n + (1 - p)$ (multiplicity $1$),
  - $lambda_2 (S^*) = dots = lambda_k (S^*) = (p - q)(n / k) + (1 - p)$ (multiplicity $k - 1$),
  - $lambda_(k+1) (S^*) = dots = lambda_n (S^*) = 1 - p$ (multiplicity $n - k$).
  Thus the eigengap separating $U^*$ from its orthogonal complement is identically $gamma_k = lambda_k (S^*) - lambda_(k+1) (S^*) = (p - q) n / k$.

  *Step 2 (Perturbation energy):* Write $E = S - S^* = (1 / r) sum_(q'=1)^r (M^((q')) - S^*)$. Because $M^((q'))_(i i) = 1$ deterministically and the $r$ replicates are independent with $M^((q'))_(i j) in {0, 1}$ (variance at most $1 / 4$),
  $ bb(E)[norm(E)_F^2] = sum_(i != j) "Var"(S_(i j)) = 1 / r sum_(i != j) "Var"(M^((1))_(i j)) <= (n(n - 1)) / (4 r) < n^2 / (4 r). $
  By Markov's inequality on $norm(E)_F^2$, with probability at least $1 - eta$,
  $ norm(E)_2 <= norm(E)_F <= n / (2 sqrt(r eta)). $

  *Step 3 (Davis–Kahan $sin Theta$ bound @davis1970rotation):* Let $hat(U)_perp in bb(R)^(n times (n - k))$ span the bottom $n - k$ eigenvectors of $S$, with eigenvalues $hat(Lambda)_perp$, and let $S^* U^* = U^* Lambda_1^*$. Then
  $ hat(U)_perp^top E U^* = hat(U)_perp^top (S - S^*) U^* = hat(Lambda)_perp (hat(U)_perp^top U^*) - (hat(U)_perp^top U^*) Lambda_1^*. $
  By Weyl's inequality, $lambda_(k+1) (S) <= lambda_(k+1) (S^*) + norm(E)_2$. Under $r > k^2 / (eta (p - q)^2)$, on the high-probability event above we have $norm(E)_2 < (p - q) n / (2 k) = gamma_k / 2$, so every diagonal gap satisfies $lambda_b^* - hat(lambda)_(k+a) >= gamma_k - norm(E)_2 > gamma_k / 2$. Taking squared Frobenius norms yields
  $ norm(E)_F^2 >= norm(hat(U)_perp^top E U^*)_F^2 >= (gamma_k / 2)^2 norm(hat(U)_perp^top U^*)_F^2 = (gamma_k / 2)^2 norm(sin Theta(hat(U), U^*))_F^2. $
  Taking square roots and substituting $norm(E)_F <= n / (2 sqrt(r eta))$ and $gamma_k / 2 = (p - q) n / (2 k)$ cancels $n$ identically:
  $ norm(sin Theta(hat(U), U^*))_F <= (n / (2 sqrt(r eta))) / ((p - q) n / (2 k)) = k / ((p - q) sqrt(r eta)). $
]

#subsec("E", "Label-Permutation Invariance and Supra-Consensus Selection")

#proposition(title: "Label-Permutation Invariance and Non-Invariance to Correlated Clones", origin: "Classical — Cited")[
  As noted in Strehl and Ghosh's knowledge-reuse framework @strehl2002cluster, the CSPA similarity matrix $S$ is strictly invariant under arbitrary independent permutations $sigma_q in S_(k^((q)))$ of the cluster labels of each base partition, and naturally handles missing object labels via normalized denominators $S_(i j) = sum_(q in cal(Q)_(i j)) M^((q))_(i j) / abs(cal(Q)_(i j))$. However, CSPA is *not* invariant to duplicate or strongly correlated base clusterers: inserting $w$ copies of an identical noisy labeling $lambda_"bad"$ perturbs $S$ linearly by $(w / (r + w))(M_"bad" - S)$.
]
#proof[
  Relabeling base partition $q$ by a permutation $sigma_q$ multiplies $H^((q))$ on the right by a permutation matrix $Pi_q$, which is orthogonal ($Pi_q Pi_q^top = I_(k^((q)))$). Hence
  $ H^((q)) Pi_q (H^((q)) Pi_q)^top = H^((q)) Pi_q Pi_q^top (H^((q)))^top = H^((q)) (H^((q)))^top, $
  leaving $S = (1 / r) sum_q H^((q)) (H^((q)))^top$ unchanged. Conversely, appending $w$ copies of $M_"bad"$ to $Lambda$ yields $S' = (r S + w M_"bad") / (r + w) = S + (w / (r + w))(M_"bad" - S)$, shifting every entry linearly toward $M_"bad"$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Constructing $S = (1 / r) H H^top$ by iterating over pairs within each cluster $C^((q))_c$ requires $O(sum_(q=1)^r sum_(c=1)^(k^((q))) abs(C^((q))_c)^2) <= O(n^2 r)$ arithmetic operations (equivalently $O(n^2 r)$ via direct label comparisons, or $O(n^2 r / k)$ when base clusters are balanced). Partitioning the resulting dense graph with $abs(V) = n$ and up to $abs(E) = O(n^2)$ edges via multilevel $k$-way METIS @karypis1998metis requires $O(abs(E) + n log k) = O(n^2)$ time per coarsening-refinement pass, or $O(n^2 k)$ when executing $k$-way spectral / refinement passes, yielding overall complexity $O(n^2 (r + k))$ (bounded loosely by $O(n^2 k r)$ in Strehl and Ghosh's original analysis @strehl2002cluster).

The primary operational constraint is *space*: storing the dense co-association matrix $S$ requires $Theta(n^2)$ memory ($8 n^2$ bytes in 64-bit floating point). On small-to-moderate datasets ($n <= 5 times 10^3$), however, dense matrix multiplication $H H^top$ exploits BLAS Level-3 cache blocking and SIMD vectorization, making CSPA markedly faster in wall-clock time than pointer-heavy hypergraph coarsening (HGPA) or repeated resampling loops (Monti).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n^2 (r + k))$ (tight); $O(n^2 k r)$ worst-case bound @strehl2002cluster],
    [Space Complexity], [$Theta(n^2)$ for dense co-association matrix $S$],
    [Geometric Bias], [Pairwise co-occurrence smoothing; METIS favors balanced cluster sizes],
    [Hyperparameters], [Base ensemble $Lambda$, target $k$, balance tolerance $epsilon_"bal"$],
    [Determinism], [Exact $S$; stochasticity solely from METIS coarsening / spectral init],
    [Invariances], [Strictly label-permutation invariant; raw-feature agnostic (knowledge reuse)],
    [Statistical Rate], [Spectral subspace error $O(k / ((p - q) sqrt(r)))$, independent of $n$ (Thm. 4)],
  ),
  caption: [Operational and geometric profile of CSPA (Alg. 171).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Privacy-Preserving Distributed and Multi-View Clustering:* Because CSPA consumes only cluster-label vectors $Lambda$ and never accesses raw object features, independent organizations holding distinct views of the same entities (e.g., banking, telecom, and e-commerce views @cai2013multiview @kumar2011coregularized) can publish local cluster IDs and fuse them via $S = (1 / r) H H^top$ without disclosing sensitive feature tables.

*2. Multi-Algorithm Integration on Non-Convex Manifolds:* When individual base clusterers use complementary inductive biases (e.g., combining $K$-means over multiple $k^((q)) > k$, spectral clustering, and density-based runs), the pairwise co-association graph $S$ acts as a learned non-linear affinity kernel. Objects traversing the same curved manifold consistently co-cluster in fine-grained partitions, enabling graph cuts on $S$ to trace non-convex shapes seamlessly.

*3. Bioinformatic Multi-Platform Sample Stratification:* Integrating mRNA expression, DNA methylation, and proteomics clusterings on patient cohorts ($n approx 100$–$2000$) fits CSPA's sweet spot: $O(n^2)$ memory is negligible, heterogeneous cluster counts $k^((q))$ across assays require no cross-assay cluster alignment, and Theorem 4 guarantees noise cancellation across noisy modalities.

*4. Document and Image Segmentation Ensembles:* Combining multiple bag-of-words, embedding, or superpixel segmentations produces sparse-to-moderate similarity graphs $S$ where balanced $k$-way partitioning cleanly removes isolated mis-segmentations while enforcing coherent region granularity.

*5. Crowd-Sourced Human Taxonomy Aggregation:* When $r$ human annotators independently group $n$ items into idiosyncratic categories of varying cardinality $k^((q))$ without a shared codebook, CSPA's permutation invariance (Prop. 5) and equivalence to Mirkin median partitioning (Thm. 2) yield an objective consensus taxonomy directly from pairwise co-grouping frequencies.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CSPA faces three sharp structural limitations:
(i) _Quadratic $O(n^2)$ memory and time wall_: beyond $n approx 3 times 10^4$, explicit formation of $S$ exhausts RAM unless replaced by $k$-NN sparsification or switched to MCLA (Alg. 173).
(ii) _Balanced-cut distortion on skewed clusters_ (Prop. 3): METIS's default load-imbalance tolerance splits large dominant clusters even under $100\%$ ensemble agreement; replacing METIS with normalized spectral cuts or UPGMA mitigates this failure.
(iii) _Ensemble-correlation dominance_ (Prop. 5): including many near-identical weak clusterers can drown out a minority of high-accuracy base models unless weights $w_q$ are learned or base diversity is enforced.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Multi-View / Privacy Fusion], [Excellent: feature-free knowledge reuse via $S$], [Down-weight redundant views],
    [Non-Convex Manifolds ($n <= 10^4$)], [High: co-association kernel captures local connectivity], [Prefer spectral cut if sizes vary],
    [Multi-Omics Cohorts], [High: handles heterogeneous $k^((q))$ naturally], [Check cluster-size skew (Prop. 3)],
    [Imbalanced Anomaly / Fraud], [Poor under METIS balance: forcibly splits majority], [Use UPGMA (Alg. 170) or relaxed Ncut],
    [Large-Scale ($n > 5 times 10^4$)], [Infeasible: $Theta(n^2)$ dense matrix $S$], [Switch to MCLA (Alg. 173) or HGPA (Alg. 172)],
  ),
  caption: [Domain suitability and failure-mode matrix for CSPA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within Strehl and Ghosh's @strehl2002cluster ensemble triad (CSPA, HGPA, MCLA), CSPA is the sole _object-level_ pairwise similarity method: HGPA (Alg. 172) partitions the hypergraph $H$ directly without expanding hyperedges into cliques, and MCLA (Alg. 173) collapses the $m times m$ cluster-level meta-graph $H^top H$ instead of the $n times n$ object-level Gram matrix $H H^top$.

#proposition(title: "Dual Gram Matrix Symmetry Between CSPA and MCLA", origin: "Original Proof")[
  The nonzero eigenvalues of the CSPA object co-association matrix $S = (1 / r) H H^top in bb(R)^(n times n)$ coincide identically (including algebraic multiplicities) with the nonzero eigenvalues of the unnormalized MCLA cluster-overlap Gram matrix $G_M = (1 / r) H^top H in bb(R)^(m times m)$, and their top-$k$ singular subspaces are linked bijectively by $V_k = H^top U_k Sigma_k^(-1)$.
]
#proof[
  Let $H = U Sigma V^top$ be the thin singular value decomposition of $H in {0, 1}^(n times m)$ with $rho = "rank"(H) <= min(n, m - r + 1)$ positive singular values $sigma_1 >= dots >= sigma_rho > 0$. Then $S = (1 / r) U Sigma^2 U^top$ and $G_M = (1 / r) V Sigma^2 V^top$, so both share the exact nonzero spectrum ${sigma_1^2 / r, dots, sigma_rho^2 / r}$. For the leading $k <= rho$ singular triplets $(U_k, Sigma_k, V_k)$, multiplying $H^top = V Sigma U^top$ on the right by $U_k Sigma_k^(-1)$ yields $H^top U_k Sigma_k^(-1) = V_k$, establishing an exact spectral duality between CSPA's object-space partitioning and MCLA's cluster-space partitioning (with MCLA additionally applying Jaccard normalization to $G_M$).
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* Across the four standardized ClusterOverview datasets, CSPA delivers the strongest overall separation quality and fastest execution among the consensus algorithms (Algs. 170–173): *Iris* $"Sil" = 0.450$, $"ARI" = 0.563$; *Customers* $"Sil" = 0.367$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.303$ ($"ARI" = "N/A"$); with mean runtime $2.957$ ms. CSPA leads every peer in Silhouette across all four benchmarks—most dramatically on *Multi-Shape* ($0.712$ vs. $0.611$ for MCLA, $0.609$ for Monti, and $0.563$ for HGPA) and *High-D* ($0.303$ vs. $0.257$–$0.259$). This empirical dominance aligns tightly with Theorems 1, 2, and 4: dense pairwise co-association smoothing $S = (1 / r) H H^top$ preserves fine-grained non-convex manifold neighborhoods and maximizes the $O(n)$-separated cut signal across roughly balanced benchmark clusters, while BLAS-accelerated matrix multiplication keeps runtime under $3$ ms at benchmark sample sizes. On *Iris*, its slightly lower ARI ($0.563$) relative to Monti ($0.634$) reflects balanced-cut pressure (Prop. 3) across the unequal overlap of _versicolor_ and _virginica_.
