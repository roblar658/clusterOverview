#import "../template/ieee.typ": *

#alg_entry(
  id: 97,
  name: "Dirichlet Multinomial Mixture (DMM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2014",
  paradigm: "Collapsed Bayesian Mixture of Multinomials",
  complexity: $cal(O)(T K bar(L) n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Short-text corpora—such as microblog posts, search queries, clinical chief
complaints, and system log lines—typically contain only $N_i in [5, 30]$ word
tokens per document, rendering multi-topic mixed-membership models such as
Latent Dirichlet Allocation (LDA) @blei2003lda statistically under-determined
at the per-document level. The *Dirichlet Multinomial Mixture (DMM)* model,
popularized for short-text clustering through the Gibbs Sampling DMM (GSDMM)
formulation of Yin and Wang @yin2014dmm, enforces a parsimonious
*single-cluster-per-document* inductive bias while integrating out continuous
simplex parameters analytically via Dirichlet-Multinomial conjugacy.

Let $cal(D) = {bold(w)_1, ..., bold(w)_n}$ denote a corpus of $n$ count
documents over a finite vocabulary $cal(V) = {1, ..., V}$, where each
document $d_i$ is represented by its term-frequency vector
$bold(w)_i = (w_(i,1), ..., w_(i,V))^top in bb(Z)_(>= 0)^V$ with length
$N_i = sum_(v=1)^V w_(i,v) >= 1$ and $U_i = |{v : w_(i,v) > 0}|$ unique terms.
Given an upper bound $K in bb(Z)_(>= 1)$ on the cluster count and symmetric
Dirichlet concentration hyperparameters $alpha > 0$ and $beta > 0$, the DMM
generative hierarchy is:
$
  bold(pi) &~ "Dir"_K (alpha, ..., alpha), \
  bold(phi)_k &~ "Dir"_V (beta, ..., beta), quad k = 1, ..., K, \
  z_i | bold(pi) &~ "Cat"(bold(pi)), quad i = 1, ..., n, \
  bold(w)_i | z_i = k, bold(phi)_k &~ "Mult"(N_i, bold(phi)_k),
$
where $bold(pi) in Delta^(K-1)$ encodes cluster mixing proportions,
$bold(phi)_k in Delta^(V-1)$ denotes the unigram vocabulary distribution of
cluster $k$, and $z_i in {1, ..., K}$ is the latent assignment of document $i$.

#definition(title: "Collapsed Dirichlet-Multinomial Joint Model")[
  Let $bold(z) = (z_1, ..., z_n) in {1, ..., K}^n$ be a cluster assignment
  vector inducing sufficient statistics $m_k = sum_(i=1)^n bb(I)_(z_i = k)$
  (documents in cluster $k$), $n_(k,v) = sum_(i: z_i = k) w_(i,v)$ (occurrences
  of word $v$ in cluster $k$), and $n_k = sum_(v=1)^V n_(k,v)$ (total tokens in
  cluster $k$). Marginalizing over $(bold(pi), bold(phi)_1, ..., bold(phi)_K)$
  defines the *collapsed joint distribution*
  $
    p(cal(D), bold(z) | alpha, beta) =
    integral_(Delta^(K-1)) integral_((Delta^(V-1))^K)
    p(cal(D), bold(z), bold(pi), bold(Phi) | alpha, beta)
    dif bold(Phi) dif bold(pi),
  $
  whose discrete posterior $p(bold(z) | cal(D), alpha, beta)$ induces a
  compound Polya urn over both document-to-cluster counts and within-cluster
  word counts.
]

Because both Dirichlet priors are strictly conjugate to their corresponding
categorical and multinomial likelihoods, all $K V + K - 1$ continuous simplex
parameters can be integrated out in exact closed form, yielding a collapsed
Gibbs sampler (GSDMM @yin2014dmm) that operates purely on integer count tables
and automatically empties redundant clusters as it balances *completeness*
(preferring populous clusters) against *homogeneity* (preferring lexically
aligned clusters).

#algbox(
  title: "Collapsed Gibbs Sampling for DMM (GSDMM)",
  inputs: [Corpus $cal(D) = {bold(w)_i}_(i=1)^n$, upper bound $K$,
    hyperparameters $alpha > 0, beta > 0$, sweeps $T$],
  outputs: [MAP / posterior cluster assignments $bold(z) in {1, ..., K}^n$,
    effective cluster count $K_"eff"$],
)[
  + *Initialize:* Sample $z_i ~ "Uniform"({1, ..., K})$ for $i = 1, ..., n$;
    accumulate counts $m_k$, $n_k$, and $n_(k,v)$ for all $k, v$.
  + *For* sweep $t = 1, ..., T$ *do:*
    + *For* each document $i = 1, ..., n$ (in random order) *do:*
      + Remove $d_i$ from $c = z_i$: decrement $m_c <- m_c - 1$,
        $n_c <- n_c - N_i$, and $n_(c,v) <- n_(c,v) - w_(i,v)$ for $v in d_i$.
      + For each $k = 1, ..., K$, compute the unnormalized log-weight
        $
          log tilde(p)_k &= log(m_k + alpha)
          - sum_(r=0)^(N_i - 1) log(n_k + V beta + r) \
          &quad + sum_(v in d_i) sum_(j=0)^(w_(i,v) - 1)
          log(n_(k,v) + beta + j).
        $
      + Normalize $p_k = exp(log tilde(p)_k - M) / sum_(l=1)^K exp(log tilde(p)_l - M)$
        with $M = max_l log tilde(p)_l$; sample $z_i ~ "Cat"(p_1, ..., p_K)$.
      + Reinstate $d_i$ into $z_i$: increment $m_(z_i)$, $n_(z_i)$, and
        $n_(z_i, v)$ for all $v in d_i$.
  + *Return:* Final (or highest-joint-probability) $bold(z)$ and
    $K_"eff" = |{k : m_k > 0}|$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Collapsed Marginalization and Gibbs Conditional")

We first establish the exact closed-form collapsed joint density and full
conditional transition law without discarding repeated-token multiplicities
$w_(i,v) > 1$.

#theorem(
  title: "Exact Collapsed Joint Distribution and Gibbs Conditional",
  origin: "Original Derivation",
)[
  Let $C_i = N_i ! / product_(v=1)^V w_(i,v)!$ denote the multinomial
  combinatorial factor of document $i$. For any $alpha > 0$ and $beta > 0$,
  the exact collapsed joint probability of the DMM model is
  $
    p(cal(D), bold(z) | alpha, beta) &=
    (product_(i=1)^n C_i)
    (Gamma(K alpha) / Gamma(n + K alpha))
    (product_(k=1)^K Gamma(m_k + alpha) / Gamma(alpha)) \
    &quad times product_(k=1)^K [
      Gamma(V beta) / Gamma(n_k + V beta)
      product_(v=1)^V Gamma(n_(k,v) + beta) / Gamma(beta)
    ].
  $
  Consequently, the exact single-site Gibbs conditional distribution for $z_i$
  given $bold(z)_(-i)$ and $cal(D)$ is
  $
    p(z_i = k | bold(z)_(-i), cal(D)) prop
    (m_k^(-i) + alpha) / (n - 1 + K alpha) dot.op
    (product_(v in d_i) product_(j=0)^(w_(i,v) - 1) (n_(k,v)^(-i) + beta + j))
    / (product_(r=0)^(N_i - 1) (n_k^(-i) + V beta + r)),
  $
  where superscript $(-i)$ denotes sufficient statistics with $d_i$ removed.
]

#proof[
  Under the hierarchical specification, the uncollapsed joint density factors
  into an assignment integral over $bold(pi)$ and $K$ independent cluster
  integrals over $bold(phi)_1, ..., bold(phi)_K$:
  $
    p(cal(D), bold(z) | alpha, beta) =
    (product_(i=1)^n C_i) dot.op I_bold(pi) (bold(z)) dot.op
    product_(k=1)^K I_(bold(phi), k) (cal(D), bold(z)).
  $
  Using the multivariate Beta identity
  $integral_(Delta^(M-1)) product_(j=1)^M u_j^(a_j - 1) dif bold(u) = (product_(j=1)^M Gamma(a_j)) / Gamma(sum_(j=1)^M a_j)$
  for any $a_j > 0$, the mixing-proportion integral evaluates to
  $
    I_bold(pi) (bold(z)) &=
    Gamma(K alpha) / Gamma(alpha)^K
    integral_(Delta^(K-1)) product_(k=1)^K pi_k^(m_k + alpha - 1) dif bold(pi) \
    &= (Gamma(K alpha) / Gamma(n + K alpha))
    product_(k=1)^K Gamma(m_k + alpha) / Gamma(alpha),
  $
  since $sum_(k=1)^K m_k = n$. Identically, for each cluster $k in {1, ..., K}$,
  integrating $bold(phi)_k$ over $Delta^(V-1)$ yields
  $
    I_(bold(phi), k) (cal(D), bold(z)) &=
    Gamma(V beta) / Gamma(beta)^V
    integral_(Delta^(V-1)) product_(v=1)^V phi_(k,v)^(n_(k,v) + beta - 1)
    dif bold(phi)_k \
    &= (Gamma(V beta) / Gamma(n_k + V beta))
    product_(v=1)^V Gamma(n_(k,v) + beta) / Gamma(beta),
  $
  since $sum_(v=1)^V n_(k,v) = n_k$. Multiplying these factors establishes the
  first identity.

  To obtain the full conditional, note that
  $p(z_i = k | bold(z)_(-i), cal(D)) prop p(cal(D), z_i = k, bold(z)_(-i)) / p(cal(D)_(-i), bold(z)_(-i))$.
  Adding document $i$ to cluster $k$ modifies only the $k$-th cluster counts,
  replacing $(m_k^(-i), n_k^(-i), n_(k,v)^(-i))$ by
  $(m_k^(-i) + 1, n_k^(-i) + N_i, n_(k,v)^(-i) + w_(i,v))$, while incrementing
  the total document count from $n - 1$ to $n$. All factors for clusters
  $l != k$, as well as $C_i$ and the normalizing constants $Gamma(K alpha)$ and
  $Gamma(V beta)$, cancel identically. Applying the rising-factorial identity
  $Gamma(x + m) / Gamma(x) = product_(j=0)^(m-1) (x + j)$ for integer $m >= 1$
  to $m_k^(-i) + alpha$, $n_k^(-i) + V beta$, and $n_(k,v)^(-i) + beta$ yields
  $
    p(z_i = k | bold(z)_(-i), cal(D)) prop
    (m_k^(-i) + alpha) / (n - 1 + K alpha) dot.op
    (product_(v in d_i) product_(j=0)^(w_(i,v) - 1) (n_(k,v)^(-i) + beta + j))
    / (product_(r=0)^(N_i - 1) (n_k^(-i) + V beta + r)),
  $
  completing the derivation.
]

#subsec("B", "Exact Phase Transition for Automatic Cluster-Count Pruning")

A defining operational mechanism of GSDMM is its ability to converge to an
active subset of $K_"eff" <= K$ non-empty clusters without fixing the exact
number of topics beforehand @yin2014dmm. We derive the exact sharp boundary
governing when a lexically novel document spawns an empty cluster versus being
absorbed into an existing populous cluster.

#theorem(
  title: "Exact Empty-Cluster Re-Population Threshold",
  origin: "Original Proof",
)[
  Suppose cluster $k$ is currently empty ($m_k^(-i) = 0, n_k^(-i) = 0$) while
  cluster $c$ is occupied ($m_c^(-i) > 0, n_c^(-i) > 0$) and shares zero
  vocabulary overlap with document $d_i$ ($n_(c,v)^(-i) = 0$ for every
  $v in d_i$). Then the exact posterior odds of assigning $d_i$ to the empty
  cluster $k$ rather than the disjoint occupied cluster $c$ satisfy
  $
    (p(z_i = k | bold(z)_(-i), cal(D))) / (p(z_i = c | bold(z)_(-i), cal(D)))
    = alpha / (m_c^(-i) + alpha)
    product_(r=0)^(N_i - 1) (1 + n_c^(-i) / (V beta + r)).
  $
  In particular:
  + For a single-token document ($N_i = 1$), the empty cluster strictly
    dominates ($p_k > p_c$) if and only if the mean document length
    $bar(L)_c = n_c^(-i) / m_c^(-i)$ of cluster $c$ satisfies
    $bar(L)_c > (V beta) / alpha$.
  + For general length $N_i >= 1$, a sufficient condition for cluster spawning
    ($p_k > p_c$) is
    $(1 + n_c^(-i) / (V beta + N_i - 1))^(N_i) > 1 + m_c^(-i) / alpha$,
    whereas if $(1 + n_c^(-i) / (V beta))^(N_i) < 1 + m_c^(-i) / alpha$, then
    even a completely disjoint document is preferentially absorbed into the
    populous cluster $c$.
]

#proof[
  Because $n_(c,v)^(-i) = 0 = n_(k,v)^(-i)$ for every token $v in d_i$, the
  numerator product
  $product_(v in d_i) product_(j=0)^(w_(i,v) - 1) (beta + j)$
  is identical for clusters $k$ and $c$. Taking the ratio of the Gibbs
  conditionals from Theorem 1 cancels both the common denominator
  $n - 1 + K alpha$ and this numerator product identically, yielding
  $
    p_k / p_c &=
    alpha / (m_c^(-i) + alpha) dot.op
    (product_(r=0)^(N_i - 1) (n_c^(-i) + V beta + r))
    / (product_(r=0)^(N_i - 1) (V beta + r)) \
    &= alpha / (m_c^(-i) + alpha)
    product_(r=0)^(N_i - 1) (1 + n_c^(-i) / (V beta + r)).
  $
  Setting $N_i = 1$ gives $p_k / p_c > 1$ if and only if
  $alpha (V beta + n_c^(-i)) > (m_c^(-i) + alpha) V beta$, which simplifies
  after subtracting $alpha V beta$ from both sides to
  $alpha n_c^(-i) > m_c^(-i) V beta$, i.e., $bar(L)_c = n_c^(-i) / m_c^(-i) > (V beta) / alpha$.
  For arbitrary $N_i >= 1$, monotonicity of $r |-> 1 + n_c^(-i) / (V beta + r)$
  on $r in {0, ..., N_i - 1}$ bounds the product from below by
  $(1 + n_c^(-i) / (V beta + N_i - 1))^(N_i)$ and from above by
  $(1 + n_c^(-i) / (V beta))^(N_i)$. Comparing each bound to
  $(m_c^(-i) + alpha) / alpha = 1 + m_c^(-i) / alpha$ yields the stated
  sufficient and collapse regimes.
]

Theorem 2 quantitatively explains the empirically observed sensitivity of
$K_"eff"$ to the product $V beta$ in short-text clustering @yin2014dmm: when
$V beta$ is too large relative to $n_c^(-i)$ and $N_i$ is small, the
"completeness" factor $m_c^(-i) + alpha$ overwhelms lexical disjointness and
collapses distinct topics into a single mega-cluster.

#subsec("C", "Ergodicity and Geometric Convergence of the Gibbs Chain")

#theorem(
  title: "Irreducible Geometric Ergodicity of Collapsed Gibbs DMM",
  origin: "Classical — Cited",
)[
  For any finite corpus $cal(D)$ and hyperparameters $alpha > 0, beta > 0$, the
  single-site collapsed Gibbs Markov chain on the finite state space
  $cal(S) = {1, ..., K}^n$ has unique stationary distribution
  $pi(bold(z)) = p(bold(z) | cal(D), alpha, beta)$ and converges geometrically
  in total variation distance: there exist $C > 0$ and $rho in [0, 1)$ such
  that for every initial distribution $mu_0$ on $cal(S)$ and all $t >= 1$,
  $
    norm(mu_0 P_"GS"^t - pi)_"TV" <= C rho^t.
  $
  Moreover, $pi(bold(z))$ is strictly invariant under the symmetric label-permuting
  group $S_K$ acting on cluster indices @yin2014dmm @bishop2006pattern.
]

#proof[
  Each single-site update $P_i$ transitions $bold(z)$ to $bold(z)'$ (differing
  only at coordinate $i$) with probability $p(z_i' | bold(z)_(-i), cal(D))$.
  For any $bold(z), bold(z)'$ agreeing off coordinate $i$,
  $
    pi(bold(z)) P_i (bold(z), bold(z)')
    = p(bold(z)_(-i) | cal(D)) p(z_i | bold(z)_(-i), cal(D))
      p(z_i' | bold(z)_(-i), cal(D))
    = pi(bold(z)') P_i (bold(z)', bold(z)),
  $
  so each kernel $P_i$ satisfies detailed balance with respect to $pi$. Hence
  both the random-scan mixture $P_"RS" = n^(-1) sum_i P_i$ and the
  systematic-scan composition $P_"GS" = P_1 P_2 ... P_n$ leave $pi$ invariant.

  Because $alpha > 0$ and $beta > 0$, every numerator factor $m_k^(-i) + alpha$
  and $n_(k,v)^(-i) + beta + j$ in Theorem 1 is strictly positive, while every
  denominator factor is strictly finite. Thus
  $p_"min" := min_(i, bold(z), k) p(z_i = k | bold(z)_(-i), cal(D)) > 0$.
  Consequently, after one complete sweep $P_"GS"$, every state
  $bold(z)' in cal(S)$ is reachable from every state $bold(z) in cal(S)$ with
  strictly positive probability $P_"GS" (bold(z), bold(z)') >= p_"min"^n > 0$.
  This yields a uniform Doeblin minorization
  $P_"GS" (bold(z), dot.op) >= epsilon nu(dot.op)$ with $epsilon = |cal(S)| p_"min"^n > 0$
  and $nu$ uniform on $cal(S)$, ensuring irreducibility, aperiodicity, unique
  stationarity of $pi$, and geometric contraction
  $norm(mu_0 P_"GS"^t - pi)_"TV" <= (1 - epsilon)^t$ (setting $rho = 1 - epsilon < 1$).
  Finally, because the prior parameters $(alpha, ..., alpha)$ and
  $(beta, ..., beta)$ are identical across all $k$, any permutation
  $sigma in S_K$ leaves the product over $k = 1, ..., K$ in Theorem 1
  unchanged, giving $pi(sigma(bold(z))) = pi(bold(z))$.
]

#subsec("D", "Asymptotic Reduction to Regularized Multinomial KL-Means")

We connect DMM's probabilistic objective directly to Bregman divergence
clustering @banerjee2005clustering on the probability simplex.

#proposition(
  title: "Reduction of MAP-DMM to Length-Weighted KL-Bregman Clustering",
  origin: "Original Proof",
)[
  Let $hat(bold(p))_i = bold(w)_i / N_i in Delta^(V-1)$ denote the empirical
  unigram distribution of document $i$. Under the flat Dirichlet prior
  ($alpha = 1, beta = 1$), maximizing the uncollapsed MAP log-posterior
  $log p(cal(D), bold(z), bold(pi), bold(Phi))$ over $(bold(z), bold(pi), bold(Phi))$
  is strictly equivalent to minimizing the length-weighted Kullback-Leibler
  objective
  $
    cal(J)_"KL" (bold(z), bold(pi), bold(Phi)) =
    sum_(i=1)^n ( N_i D_"KL" (hat(bold(p))_i || bold(phi)_(z_i)) - log pi_(z_i) ).
  $
  For general $beta > 1$, the MAP estimate $hat(bold(phi))_k$ adds a
  cross-entropy Dirichlet barrier $(beta - 1) V D_"KL" (bold(u)_V || bold(phi)_k)$
  toward the uniform distribution $bold(u)_V = (1/V, ..., 1/V)^top$.
]

#proof[
  When $alpha = 1$ and $beta = 1$, the prior densities $p(bold(pi))$ and
  $p(bold(phi)_k)$ are constant on their simplices. Expanding the negative
  complete-data log-likelihood for document $i$ assigned to $z_i = k$:
  $
    -log p(bold(w)_i, z_i = k | bold(pi), bold(Phi)) &=
    -log C_i - log pi_k - sum_(v=1)^V w_(i,v) log phi_(k,v) \
    &= -log C_i - log pi_k - N_i sum_(v=1)^V hat(p)_(i,v) log phi_(k,v).
  $
  Adding and subtracting the empirical Shannon entropy
  $N_i H(hat(bold(p))_i) = -N_i sum_(v=1)^V hat(p)_(i,v) log hat(p)_(i,v)$
  (using $0 log 0 = 0$) yields
  $
    -N_i sum_(v=1)^V hat(p)_(i,v) log phi_(k,v)
    = N_i H(hat(bold(p))_i) + N_i D_"KL" (hat(bold(p))_i || bold(phi)_k).
  $
  Since $-log C_i + N_i H(hat(bold(p))_i)$ depends only on the fixed data
  $bold(w)_i$, summing over $i = 1, ..., n$ proves exact equivalence with
  minimizing $cal(J)_"KL"$. For $beta > 1$, the negative log-prior contributes
  $-(beta - 1) sum_(v=1)^V log phi_(k,v) = (beta - 1) V D_"KL" (bold(u)_V || bold(phi)_k) + "const"$,
  which corresponds to Laplace-smoothed Bregman centroids
  $hat(phi)_(k,v) = (n_(k,v) + beta - 1) / (n_k + V (beta - 1))$
  @banerjee2005clustering.
]

#subsec("E", "Constructive Breakdown on Multi-Topic Documents")

#proposition(
  title: "Irreducible Penalty of DMM Under Multi-Topic Mixtures",
  origin: "Original Counterexample",
)[
  Partition $cal(V)$ into $M >= 2$ disjoint sub-vocabularies
  $cal(V)_1, ..., cal(V)_M$ of equal size $V / M$, and generate a corpus
  comprising equal fractions of pure documents ($2 L$ tokens drawn uniformly
  from one $cal(V)_a$) and hybrid documents ($L$ tokens uniformly from $cal(V)_a$
  and $L$ uniformly from $cal(V)_b$, $a != b$). For $M = 2$ with $K = 2$
  clusters, every DMM parameter configuration incurs an asymptotic per-document
  excess KL risk of at least $L log(4 / 3) approx 0.288 L$ nats relative to LDA
  @blei2003lda. If $K$ is unrestricted, achieving zero excess risk forces DMM
  to allocate $K_"eff" = M + binom(M, 2) = Omega(M^2)$ clusters (and $2^M - 1$
  clusters for arbitrary non-empty topic subsets).
]

#proof[
  Consider $M = 2$ with equal $1 / 3$ population weights on pure-$cal(V)_1$,
  pure-$cal(V)_2$, and balanced hybrid documents of length $2 L$. Under LDA
  with $K = 2$, setting $bold(phi)_1 = "Unif"(cal(V)_1)$,
  $bold(phi)_2 = "Unif"(cal(V)_2)$, and per-document proportions
  $bold(theta)_i in {(1,0), (0,1), (1/2, 1/2)}$ achieves zero excess KL risk.
  In a $K = 2$ DMM, each document belongs to one cluster $bold(phi)_k$. By
  word symmetry within each half, an optimal cluster $k$ places total mass
  $q_k in (0, 1)$ uniformly on $cal(V)_1$ and $1 - q_k$ uniformly on $cal(V)_2$.
  The excess KL risks relative to the true unigram distributions are
  $-2 L log q_k$ on pure-$cal(V)_1$, $-2 L log(1 - q_k)$ on pure-$cal(V)_2$,
  and $2 L (-log 2 - (1/2) log(q_k (1 - q_k)))$ on hybrids.

  There are two non-degenerate partitions of the three document types into
  $K = 2$ clusters (up to symmetry):
  + *Hybrids paired with pure-$cal(V)_1$ in cluster 1, pure-$cal(V)_2$ in
    cluster 2:* Setting $q_2^* = 0$ incurs zero risk on pure-$cal(V)_2$, while
    cluster 1 minimizes the weighted excess risk
    $
      Delta R_1 (q_1) &= 1/3 (-2 L log q_1)
      + 1/3 (-2 L log 2 - L log q_1 - L log(1 - q_1)) \
      &= (L / 3) ( -3 log q_1 - log(1 - q_1) - 2 log 2 ).
    $
    Setting $dif Delta R_1 / dif q_1 = 0$ yields the unique minimizer
    $q_1^* = 3 / 4$, where
    $
      Delta R_1 (3 / 4) = (L / 3) ( 6 log 2 - 3 log 3 )
      = L log(4 / 3) approx 0.288 L.
    $
  + *Pure-$cal(V)_1$ and pure-$cal(V)_2$ merged in cluster 2, hybrids in
    cluster 1:* Setting $q_1^* = 1 / 2$ incurs zero risk on hybrids, while
    cluster 2 uniquely minimizes $-(2 L / 3)(log q_2 + log(1 - q_2))$ at
    $q_2^* = 1 / 2$, giving $Delta R_2 = (4 L / 3) log 2 approx 0.924 L > L log(4/3)$.

  Thus every $K = 2$ DMM configuration suffers at least $L log(4 / 3)$ nats of
  excess risk per document (matching the strictly positive $Omega(L)$ bias of
  collapsing all documents into one cluster). Eliminating this bias when $K$ is
  unrestricted requires assigning every distinct non-empty co-occurring subset
  of topics its own dedicated component, forcing $K_"eff" = M + binom(M, 2)$
  for pairwise hybrids and $2^M - 1$ for general multi-topic corpora.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In each collapsed Gibbs sweep, updating document $i$ requires evaluating
$log tilde(p)_k$ across all $K$ clusters (or $K_"eff" + 1$ when empty clusters
are collapsed into one representative state). Evaluating the token product for
document $i$ takes $cal(O)(N_i)$ arithmetic operations directly, or
$cal(O)(U_i)$ operations using cached $log Gamma$ evaluations. Across $n$
documents and $T$ sweeps, total runtime is $cal(O)(T K bar(L) n)$ where
$bar(L) = n^(-1) sum_i N_i$ is the mean document length. Memory complexity is
$cal(O)(K V + n bar(U))$ using dense count tables $n_(k,v)$, reducible to
$cal(O)("nnz"(bold(N)) + n bar(U))$ with sparse hash maps when $V$ exceeds
$10^5$. Geometrically, DMM partitions documents on the $(V-1)$-simplex via
smoothed discrete Polya / KL-divergence neighborhoods (Proposition 1) rather
than Euclidean balls, rendering it strictly invariant to vocabulary
permutations and asymptotically invariant to document ordering.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T K bar(L) n)$ collapsed Gibbs; $cal(O)(T K_"eff" bar(U) n)$ with empty-cluster pooling and $log Gamma$ cache],
    [Space Complexity], [$cal(O)(K V + n bar(U))$ dense; $cal(O)("nnz"(bold(N)) + n bar(U))$ sparse],
    [Cluster Geometry], [Simplex KL-divergence / Polya cones; single topic per document],
    [Model Selection], [Automatic pruning $K_"eff" <= K$ controlled by $(alpha, beta)$ (Theorem 2)],
    [Invariances], [Vocabulary-permutation invariant; asymptotically order-exchangeable],
    [Key Hyperparameters], [Upper bound $K$, cluster prior $alpha$, word smoothing $beta$],
  ),
  caption: [Operational and complexity profile of DMM / GSDMM.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Microblogging and Social-Media Event Clustering:* Tweets and status updates
  average $N_i approx 10$--$20$ words and overwhelmingly express a single
  event or topic. DMM's single-cluster-per-document inductive bias avoids the
  extreme sparsity failure of per-document topic vectors in LDA @yin2014dmm,
  while automatic pruning adapts $K_"eff"$ as breaking events emerge.
+ *Search Query Intent and Support-Ticket Triage:* Web queries ($N_i in [2, 6]$)
  and customer-service subject lines have near-zero term co-occurrence within a
  single record. Collapsed Dirichlet smoothing ($beta$) pools lexical counts at
  the cluster level, letting rare synonyms co-cluster through shared cluster
  statistics.
+ *Metagenomic 16S rRNA Microbiome Enterotyping:* Operational Taxonomic Unit
  (OTU) abundance tables consist of non-negative discrete read counts over
  thousands of microbial taxa with high overdispersion relative to a plain
  multinomial. Integrating out $bold(phi)_k ~ "Dir"_V (beta)$ captures
  biological count overdispersion naturally via the Dirichlet-Multinomial
  marginal likelihood.
+ *Biomedical Title and Clinical Chief-Complaint Stratification:* Emergency
  department chief complaints and PubMed article titles are brief, domain-specific
  token sequences whose counts obey Zipfian sparsity; DMM's KL geometry weights
  discriminative rare clinical terms far more faithfully than Euclidean $L_2$.
+ *AIOps Log-Template and Security-Alert Grouping:* Automated server logs and
  intrusion alerts consist of short tokenized messages with repeated template
  vocabularies; GSDMM's $cal(O)(T K_"eff" bar(L) n)$ updates support rapid
  clustering without pre-specifying exact fault counts.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DMM degrades in four identifiable regimes: (i) *long multi-topic documents*
(articles, books), where Proposition 2 proves an irreducible $Omega(L)$ excess
KL risk or quadratic/exponential cluster proliferation; (ii) *continuous
signed Euclidean features*, which violate discrete count assumptions and
require non-negative quantization that discards covariance orientation;
(iii) *mis-tuned $V beta$*, where Theorem 2 shows excessive $beta$ collapses
disjoint clusters while vanishing $beta$ fragments near-duplicates; and
(iv) *high-frequency stopword dominance*, since length weighting $N_i$ in
Proposition 1 lets uninformative frequent tokens dominate the KL sum unless
filtered.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Short texts / tweets], [Exact single-topic Polya match; automatic $K_"eff"$], [Tune $alpha, beta$ per Theorem 2; compare BTM @yan2013btm],
    [Microbiome OTU counts], [Dirichlet-Multinomial models count overdispersion], [Normalizing library size if $N_i$ varies by $10^3 times$],
    [Long multi-topic docs], [Poor; $Omega(L)$ bias or cluster explosion (Prop. 2)], [Use LDA @blei2003lda or CTM @blei2007ctm],
    [Continuous tabular data], [Weak; requires binning / non-negative rescaling], [Use GMM @fraley2002mclust or Skew-Normal @lin2009skewnormal],
    [Streaming huge $V$], [Good with sparse maps; exact online counts], [Monitor label switching across restarts],
  ),
  caption: [Domain suitability matrix for Dirichlet Multinomial Mixture.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Structurally, DMM occupies the exact *single-topic boundary* of the Bayesian
topic-modeling hierarchy: taking the per-document Dirichlet concentration
$alpha_"LDA" -> 0^+$ in LDA @blei2003lda forces every document's topic mixture
$bold(theta)_i$ onto a vertex of $Delta^(K-1)$, recovering DMM exactly. In the
opposite direction, letting $K -> infinity$ with prior $alpha_0 / K$ turns
DMM's finite Dirichlet prior on $bold(pi)$ into the Chinese Restaurant Process
used in the Infinite Gaussian Mixture Model (Alg. 99 @rasmussen2000igmm) and
Hierarchical Dirichlet Processes @teh2006hdp. Compared to Biterm Topic Model
(BTM @yan2013btm), which models unordered word pairs corpus-wide, DMM preserves
document-level assignment directly; compared to spherical $k$-means / von
Mises-Fisher mixtures @dhillon2001concept @banerjee2005vmf and AutoClass
@cheeseman1996autoclass, DMM replaces cosine or independent-feature distances
with exact Bayesian Polya-multinomial evidence.

On the ClusterOverview 4-dataset continuous benchmark (with non-negative
count-scaled features), DMM achieves Iris Silhouette $= 0.586$ and ARI $= 0.567$
@rousseeuw1987silhouettes @hubert1985comparing, Customers Silhouette $= 0.324$,
Multi-Shape Silhouette $= 0.712$ (ARI N/A), High-D Silhouette $= 0.287$
(ARI N/A), and a mean runtime of $16.789$ ms. This empirical profile reflects
our theoretical analysis directly: when clusters occupy well-separated positive
radial cones (Multi-Shape and the Setosa component of Iris), simplex KL
geometry yields strong geometric separation (Silhouette $0.712$ and $0.586$),
while the moderate Iris ARI ($0.567$, trailing Skew-Normal's $0.743$) reflects
DMM's inability to model correlated continuous covariance tilt between
overlapping Versicolor and Virginica classes. Notably, DMM attains the highest
High-D Silhouette ($0.287$) among Algs. 97--100 because Dirichlet smoothing
$beta$ regularizes sparse high-dimensional coordinates against distance
concentration @beyer1999nearest, while its $16.789$ ms runtime reflects the
per-sweep cost $cal(O)(T K bar(L) n)$ of iterative collapsed Gibbs resampling.
