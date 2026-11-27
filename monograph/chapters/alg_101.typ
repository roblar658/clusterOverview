#import "../template/ieee.typ": *

#alg_entry(
  id: 101,
  name: "Biterm Topic Model (BTM)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2013",
  paradigm: "Corpus-Level Biterm Topic Mixture",
  complexity: $O(T K |cal(B)|)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical probabilistic topic models such as Probabilistic Latent Semantic Analysis (pLSA) @hofmann1999plsa and Latent Dirichlet Allocation (LDA) @blei2003lda treat each document as a distinct mixture over $K$ latent topics and infer per-document topic proportions from document-level word co-occurrences. On short texts—such as microblog posts, search queries, titles, and clinical notes—each document contains only a handful of tokens, rendering document-level term-frequency histograms overwhelmingly sparse and causing severe posterior variance in document-specific topic mixtures. The Biterm Topic Model (BTM) @yan2013btm resolves this data-sparsity bottleneck by shifting the generative unit from individual tokens inside a document to unordered word pairs (_biterms_) aggregated across the entire corpus.

Let $cal(D) = \{d_1, dots, d_N\}$ denote a corpus of $N$ short documents over a finite vocabulary $cal(V) = \{1, dots, V\}$, where document $d$ has token sequence $(w_(d,1), dots, w_(d, L_d))$ of length $L_d >= 2$.

#definition(title: "Corpus Biterm Multiset and BTM Generative Process")[
  For each document $d in cal(D)$ and window width $W >= 2$ (with $W = L_d$ for full-document extraction), the _document biterm multiset_ is
  $ cal(B)_d := \{\{w_(d,i), w_(d,j)\} : 1 <= i < j <= L_d,\; j - i < W\}, $
  and the _corpus biterm multiset_ is the disjoint union $cal(B) := union.sq_(d=1)^N cal(B)_d$ of cardinality $M := |cal(B)|$. Given symmetric Dirichlet hyperparameters $alpha > 0$, $beta > 0$, and topic count $K >= 2$, the BTM generative process @yan2013btm specifies:
  + Draw global topic weights $bold(theta) tilde.op "Dir"(alpha bold(1)_K)$ on $Delta_(K-1)$.
  + For each topic $z in [K]$, draw word weights $bold(phi)_z tilde.op "Dir"(beta bold(1)_V)$ on $Delta_(V-1)$.
  + For each biterm $b_m = (w_(m,1), w_(m,2)) in cal(B)$, $m = 1, dots, M$:
    - Draw a single latent topic $z_m tilde.op "Cat"(bold(theta))$.
    - Draw both words conditionally independently from the same topic: $w_(m,1), w_(m,2) tilde.op "Cat"(bold(phi)_(z_m))$.
]

Integrating over the latent topic assignment $z_m in [K]$ yields the marginal probability of a biterm $b = (w_1, w_2)$ and the corpus log-likelihood:
$ P(b | bold(theta), bold(Phi)) &= sum_(z=1)^K theta_z phi_(z, w_1) phi_(z, w_2), \
  cal(L)(bold(theta), bold(Phi)) &= sum_(m=1)^M log ( sum_(z=1)^K theta_z phi_(z, w_(m,1)) phi_(z, w_(m,2)) ), $
where $bold(Phi) in [0,1]^(K times V)$ stacks row vectors $bold(phi)_z^top$. Unlike LDA, BTM introduces _zero_ per-document topic vectors during training, pooling co-occurrence statistics globally into $bold(theta)$ and $bold(Phi)$. After inferring posterior estimates $hat(bold(theta))$ and $hat(bold(Phi))$ via collapsed Gibbs sampling, document-level topic distributions are recovered via Bayes' rule over the empirical biterm distribution $P(b | d) = n_d (b) / |cal(B)_d|$:
$ P(z | d) &= sum_(b in cal(B)_d) P(z | b) P(b | d), \
  P(z | b) &= (hat(theta)_z hat(phi)_(z, w_1) hat(phi)_(z, w_2)) / (sum_(k=1)^K hat(theta)_k hat(phi)_(k, w_1) hat(phi)_(k, w_2)), $
and each document is assigned to its MAP cluster $c(d) = arg max_(z in [K]) P(z | d)$.

#algbox(
  title: "Biterm Topic Model (Collapsed Gibbs & Document Clustering)",
  inputs: [Corpus $cal(D)$, vocabulary size $V$, topics $K$, window $W$, priors $alpha, beta > 0$, sweeps $T$, burn-in $T_0 < T$.],
  outputs: [Topic proportions $hat(bold(theta))$, word distributions $hat(bold(Phi))$, document topic posteriors $\{P(z | d)\}$, hard clusters $\{c(d)\}$.],
)[
  + *Biterm Extraction:* Extract $cal(B)_d$ for every $d in cal(D)$ within window $W$ and form $cal(B) = union.sq_d cal(B)_d$ of size $M$.
  + *Initialization:* Draw $z_m^((0)) tilde.op "Unif"([K])$ for $m in [M]$; initialize topic counts $n_z := sum_m bb(I)[z_m = z]$ and topic-word counts $n_(z,w) := sum_m bb(I)[z_m = z](bb(I)[w_(m,1)=w] + bb(I)[w_(m,2)=w])$.
  + *Collapsed Gibbs Sweeps:* For $t = 1, dots, T$ and $m = 1, dots, M$ with $b_m = (w_1, w_2)$:
    - Decrement $n_(z_m)$, $n_(z_m, w_1)$, and $n_(z_m, w_2)$ by $1$.
    - For $z = 1, dots, K$, compute unnormalized weights $q_z$:
      $ q_z = (n_z + alpha) ((n_(z, w_1) + beta)(n_(z, w_2) + beta + bb(I)[w_1 = w_2])) / ((2 n_z + V beta)(2 n_z + V beta + 1)). $
    - Sample $z_m tilde.op "Cat"(q_1 / Z, dots, q_K / Z)$ with $Z = sum_k q_k$, and increment $n_(z_m)$, $n_(z_m, w_1)$, $n_(z_m, w_2)$ by $1$.
  + *Parameter Estimation:* Using post-burn-in counts ($t > T_0$), set $hat(theta)_z = (n_z + alpha)/(M + K alpha)$ and $hat(phi)_(z,w) = (n_(z,w) + beta)/(2 n_z + V beta)$.
  + *Document Clustering:* For each $d in cal(D)$, compute $P(z | d) = |cal(B)_d|^(-1) sum_(b in cal(B)_d) P(z | b)$ and return $c(d) = arg max_z P(z | d)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Collapsed Posterior and Self-Biterm Detailed Balance")

In raw short-text corpora, repeated terms within a document generate _self-biterms_ $(w, w)$. Dropping the duplicate-word increment in the Gibbs update breaks exact conjugacy. We derive the exact collapsed conditional for arbitrary biterm multisets.

#lemma(title: "Exact Collapsed Conditional and Self-Biterm Correction", origin: "Original Derivation")[
  Let $bold(z) = (z_1, dots, z_M) in [K]^M$, $n_z := sum_(m=1)^M bb(I)[z_m = z]$, and $n_(z,w) := sum_(m=1)^M bb(I)[z_m = z](bb(I)[w_(m,1)=w] + bb(I)[w_(m,2)=w])$. Integrating out $(bold(theta), bold(Phi))$ yields the exact collapsed joint distribution
  $ P(cal(B), bold(z) | alpha, beta) &= (Gamma(K alpha)) / (Gamma(M + K alpha)) product_(z=1)^K (Gamma(n_z + alpha)) / (Gamma(alpha)) \
    &times product_(z=1)^K ( (Gamma(V beta)) / (Gamma(2 n_z + V beta)) product_(w=1)^V (Gamma(n_(z,w) + beta)) / (Gamma(beta)) ). $
  Consequently, the exact full conditional for biterm $b_m = (w_1, w_2)$ given $bold(z)_(-m)$ is
  $ P(z_m = z | bold(z)_(-m), cal(B)) prop & (n_z^((-m)) + alpha) \
    &times ((n_(z, w_1)^((-m)) + beta)(n_(z, w_2)^((-m)) + beta + delta_(w_1, w_2))) / ((2 n_z^((-m)) + V beta)(2 n_z^((-m)) + V beta + 1)), $
  where $delta_(w_1, w_2) := bb(I)[w_1 = w_2]$. Omitting $delta_(w_1, w_2)$ when self-biterms occur violates detailed balance with respect to $P(bold(z) | cal(B), alpha, beta)$.
]

#proof[
  Conditioned on $(bold(theta), bold(Phi))$, the complete-data likelihood factorizes as
  $ P(cal(B), bold(z) | bold(theta), bold(Phi)) = ( product_(z=1)^K theta_z^(n_z) ) ( product_(z=1)^K product_(w=1)^V phi_(z,w)^(n_(z,w)) ), $
  because every biterm assigned to topic $z$ contributes exactly two word draws to $bold(phi)_z$, so $sum_(w=1)^V n_(z,w) = 2 n_z$. Multiplying by the independent Dirichlet priors $p(bold(theta)) = (Gamma(K alpha) / Gamma(alpha)^K) product_z theta_z^(alpha - 1)$ and $p(bold(Phi)) = product_z (Gamma(V beta) / Gamma(beta)^V) product_w phi_(z,w)^(beta - 1)$, and integrating over $Delta_(K-1)$ and $Delta_(V-1)^K$ using the Dirichlet normalizing integral $integral_(Delta_(r-1)) product_(j=1)^r u_j^(a_j - 1) dif bold(u) = (product_j Gamma(a_j)) / Gamma(sum_j a_j)$, yields the stated expression for $P(cal(B), bold(z) | alpha, beta)$.

  Now separate the terms depending on $z_m = z$ from leave-one-out counts $n_z^((-m))$ and $n_(z,w)^((-m))$. Since $P(z_m = z | bold(z)_(-m), cal(B)) prop P(cal(B), bold(z) | alpha, beta) / P(cal(B)_(-m), bold(z)_(-m) | alpha, beta)$, applying $Gamma(x + 1) = x Gamma(x)$ gives topic ratio $n_z^((-m)) + alpha$ and denominator word-count ratio $((2 n_z^((-m)) + V beta)(2 n_z^((-m)) + V beta + 1))^(-1)$ because $2 n_z = 2 n_z^((-m)) + 2$. In the numerator, if $w_1 != w_2$, both $n_(z, w_1)$ and $n_(z, w_2)$ increment by $1$, giving $(n_(z, w_1)^((-m)) + beta)(n_(z, w_2)^((-m)) + beta)$; if $w_1 = w_2$, then $n_(z, w_1) = n_(z, w_1)^((-m)) + 2$, giving $(n_(z, w_1)^((-m)) + beta)(n_(z, w_1)^((-m)) + beta + 1)$. Both cases unify via $+ delta_(w_1, w_2)$. For a single-site transition $z -> z'$ holding $bold(z)_(-m)$ fixed, detailed balance $pi(bold(z)) T(bold(z) -> bold(z)') = pi(bold(z)') T(bold(z)' -> bold(z))$ requires $T(z -> z') / T(z' -> z) = pi(bold(z)') / pi(bold(z))$. Omitting $delta_(w_1, w_2)$ when $w_1 = w_2$ alters this ratio by $(n_(z', w_1)^((-m)) + beta + 1)(n_(z, w_1)^((-m)) + beta) / [(n_(z, w_1)^((-m)) + beta + 1)(n_(z', w_1)^((-m)) + beta)] != 1$ whenever $n_(z, w_1)^((-m)) != n_(z', w_1)^((-m))$, violating detailed balance.
]

#subsec("B", "Exact Symmetric Nonnegative Factorization Equivalence")

Whereas pLSA and NMF factorize the document-term matrix @ding2008nmfplsa @lee1999nmf, BTM's corpus likelihood possesses an exact algebraic equivalence to symmetric low-rank factorization of the empirical word co-occurrence matrix under Kullback-Leibler divergence.

#theorem(title: "Equivalence to Symmetric KL Co-Occurrence Factorization", origin: "Original Proof")[
  Define the normalized symmetric empirical co-occurrence matrix $C in [0,1]^(V times V)$ with entries
  $ C_(u,v) := 1 / (2 M) sum_(m=1)^M (bb(I)[w_(m,1)=u, w_(m,2)=v] + bb(I)[w_(m,1)=v, w_(m,2)=u]), $
  and let $Q(bold(theta), bold(Phi)) := bold(Phi)^top "diag"(bold(theta)) bold(Phi) in [0,1]^(V times V)$. Then $C$ and $Q$ are valid probability distributions on $[V] times [V]$ ($sum_(u,v) C_(u,v) = sum_(u,v) Q_(u,v) = 1$), and for all $(bold(theta), bold(Phi)) in Delta_(K-1) times Delta_(V-1)^K$,
  $ 1 / M cal(L)(bold(theta), bold(Phi)) = -H(C) - D_"KL" (C || bold(Phi)^top "diag"(bold(theta)) bold(Phi)), $
  where $H(C) := -sum_(u,v) C_(u,v) log C_(u,v)$ is the empirical biterm entropy. Hence maximum-likelihood BTM is equivalent to minimizing $D_"KL" (C || Q)$ over rank-$<= K$ positive semidefinite stochastic factorizations $Q = bold(Phi)^top "diag"(bold(theta)) bold(Phi)$.
]

#proof[
  First, observe that every entry $Q_(u,v) = sum_(z=1)^K theta_z phi_(z,u) phi_(z,v) >= 0$ satisfies $Q_(u,v) = Q_(v,u)$ and
  $ sum_(u=1)^V sum_(v=1)^V Q_(u,v) = sum_(z=1)^K theta_z ( sum_(u=1)^V phi_(z,u) ) ( sum_(v=1)^V phi_(z,v) ) = sum_(z=1)^K theta_z = 1. $
  Likewise, $sum_(u,v) C_(u,v) = 1 / (2M) sum_(m=1)^M (1 + 1) = 1$, and $C_(u,v) = C_(v,u)$. Because $Q_(w_(m,1), w_(m,2)) = Q_(w_(m,2), w_(m,1))$, each summand in the normalized corpus log-likelihood rewrites symmetrically:
  $ 1 / M cal(L)(bold(theta), bold(Phi)) &= 1 / (2M) sum_(m=1)^M ( log Q_(w_(m,1), w_(m,2)) + log Q_(w_(m,2), w_(m,1)) ) \
    &= sum_(u=1)^V sum_(v=1)^V C_(u,v) log Q_(u,v). $
  Adding and subtracting $sum_(u,v) C_(u,v) log C_(u,v) = -H(C)$ (with $0 log 0 := 0$) gives
  $ 1 / M cal(L)(bold(theta), bold(Phi)) = -H(C) - sum_(u,v=1)^V C_(u,v) log (C_(u,v) / Q_(u,v)) = -H(C) - D_"KL" (C || Q). $
  Because $H(C)$ is independent of $(bold(theta), bold(Phi))$, maximizing $cal(L)$ is identically equivalent to minimizing $D_"KL" (C || bold(Phi)^top "diag"(bold(theta)) bold(Phi))$. Moreover, for any $bold(v) in bb(R)^V$, $bold(v)^top Q bold(v) = sum_(z=1)^K theta_z (bold(phi)_z^top bold(v))^2 >= 0$, confirming that $Q$ is symmetric positive semidefinite of rank at most $K$.
]

Theorem 1 reveals a sharp geometric characterization: BTM compresses the empirical word-pair graph $C$ onto the cone of symmetric completely positive stochastic matrices of cp-rank at most $K$ under Bregman/KL geometry @banerjee2005clustering.

#subsec("C", "Variance Reduction in the Short-Document Regime")

To formalize why BTM outperforms LDA @blei2003lda on short texts, we compare the exact Bayesian posterior uncertainty of topic-proportion estimation as the number of documents $N -> oo$ while document length $L$ remains fixed.

#theorem(title: "Exact Short-Text Posterior Variance Separation from LDA", origin: "Original Proof")[
  Consider a corpus of $N$ documents of equal length $L >= 2$ with full-document biterm extraction ($W = L$), so $M = N binom(L, 2)$.
  + In LDA with symmetric Dirichlet prior $alpha$ on each document's topic mixture $bold(theta)_d in Delta_(K-1)$, conditioned on token topic assignments $bold(z)_d$, the total posterior variance of $bold(theta)_d$ equals
    $ sum_(z=1)^K "Var"(theta_(d,z) | bold(z)_d) = (1 - ||tilde(bold(theta))_d||_2^2) / (L + K alpha + 1), $
    where $tilde(theta)_(d,z) = (n_(d,z) + alpha) / (L + K alpha)$. For fixed $L = O(1)$ and any non-degenerate document ($||tilde(bold(theta))_d||_2^2 <= 1 - eta$, $eta > 0$), this variance remains bounded below by $eta / (L + K alpha + 1) = Omega(1)$ as $N -> oo$.
  + In BTM conditioned on biterm assignments $bold(z)$, the total posterior variance of the corpus topic mixture $bold(theta)$ equals
    $ sum_(z=1)^K "Var"(theta_z | bold(z)) = (1 - ||tilde(bold(theta))||_2^2) / (N binom(L, 2) + K alpha + 1) <= 2 / (N L (L - 1)), $
    which vanishes at rate $O((N L^2)^(-1))$ as $N -> oo$ while eliminating all $N(K - 1)$ document-level nuisance parameters.
]

#proof[
  By Dirichlet-Multinomial conjugacy, if $bold(u) tilde.op "Dir"(a_1, dots, a_K)$ with $A := sum_(z=1)^K a_z$ and mean $tilde(u)_z = a_z / A$, each marginal follows $u_z tilde.op "Beta"(a_z, A - a_z)$, whose exact variance is
  $ "Var"(u_z) = (a_z (A - a_z)) / (A^2 (A + 1)) = (tilde(u)_z (1 - tilde(u)_z)) / (A + 1). $
  Summing over $z = 1, dots, K$ and using $sum_z tilde(u)_z = 1$ gives
  $ sum_(z=1)^K "Var"(u_z) = (1 - sum_(z=1)^K tilde(u)_z^2) / (A + 1) = (1 - ||tilde(bold(u))||_2^2) / (A + 1). $
  In LDA, conditioned on the $L$ token assignments $bold(z)_d$ of document $d$, the posterior $bold(theta)_d | bold(z)_d$ is $"Dir"(n_(d,1)+alpha, dots, n_(d,K)+alpha)$ with total concentration $A_"LDA" = L + K alpha$, independent of all other $N - 1$ documents given $bold(Phi)$. Substituting $A_"LDA"$ gives Claim 1; because $L$ is fixed, the lower bound $eta / (L + K alpha + 1)$ does not shrink with $N$.

  In BTM, conditioned on the $M = N binom(L, 2) = N L (L - 1) / 2$ biterm assignments $bold(z)$, the posterior $bold(theta) | bold(z)$ is $"Dir"(n_1 + alpha, dots, n_K + alpha)$ with total concentration $A_"BTM" = M + K alpha$. Substituting $A_"BTM"$ and noting $1 - ||tilde(bold(theta))||_2^2 < 1$ and $K alpha + 1 > 0$ yields
  $ sum_(z=1)^K "Var"(theta_z | bold(z)) = (1 - ||tilde(bold(theta))||_2^2) / (N binom(L, 2) + K alpha + 1) <= 1 / (N binom(L, 2)) = 2 / (N L (L - 1)), $
  establishing Claim 2.
]

#subsec("D", "Geometric Ergodicity of Collapsed Gibbs Sampling")

#theorem(title: "Uniform Doeblin Minorization and Geometric Ergodicity", origin: "Classical — Cited")[
  For any finite biterm multiset $cal(B)$ of size $M >= 1$, topic count $K >= 2$, vocabulary size $V >= 2$, and strictly positive hyperparameters $alpha > 0, beta > 0$, the systematic-scan collapsed Gibbs Markov chain $\{bold(z)^((t))\} _ (t >= 0)$ on $[K]^M$ is irreducible, aperiodic, and geometrically ergodic @yan2013btm @bishop2006pattern. Specifically, letting $pi(bold(z)) := P(bold(z) | cal(B), alpha, beta)$ denote the exact collapsed posterior and $P^t$ the $t$-sweep transition kernel,
  $ || P^t (bold(z)^((0)), dot.op) - pi(dot.op) ||_"TV" <= (1 - delta)^t, $
  for all initial states $bold(z)^((0)) in [K]^M$, with explicit Doeblin constant
  $ delta := ( (alpha beta^2) / ((M + K alpha)(2 M + V beta)(2 M + V beta + 1)) )^M > 0. $
]

#proof[
  At any single-biterm update $m in [M]$ from any configuration $bold(z)_(-m)$, the unnormalized weight $q_z$ of Lemma 1 satisfies $n_z^((-m)) >= 0$, $n_(z, w_r)^((-m)) >= 0$, $2 n_z^((-m)) <= 2 M$, and
  $ q_z >= (alpha beta^2) / ((2 M + V beta)(2 M + V beta + 1)) =: c_"min" > 0. $
  Meanwhile, for each $k in [K]$, every word-count ratio is strictly bounded above by $1$ because $(n_(k, w_1)^((-m)) + beta) <= 2 n_k^((-m)) + V beta$ and $(n_(k, w_2)^((-m)) + beta + delta_(w_1, w_2)) <= 2 n_k^((-m)) + V beta + 1$. Thus the normalization constant satisfies
  $ Z = sum_(k=1)^K q_k <= sum_(k=1)^K (n_k^((-m)) + alpha) = M - 1 + K alpha < M + K alpha. $
  Consequently, every single-site conditional probability is uniformly bounded below by
  $ P(z_m = z | bold(z)_(-m), cal(B)) = q_z / Z >= (alpha beta^2) / ((M + K alpha)(2 M + V beta)(2 M + V beta + 1)) =: epsilon_"min". $
  In one complete systematic sweep over $m = 1, dots, M$, the transition probability to any target $bold(z)' in [K]^M$ from any starting state $bold(z)$ satisfies $P(bold(z), bold(z)') >= epsilon_"min"^M = delta$. Because $pi(bold(z)') <= 1$ for every $bold(z)' in [K]^M$, we obtain the uniform Doeblin minorization $P(bold(z), A) >= delta pi(A)$ for all $bold(z) in [K]^M$ and $A subset.eq [K]^M$ (with $delta <= K^(-M) < 1$). Since each single-site Gibbs update preserves $pi$ by Lemma 1, their composition $P$ has unique invariant distribution $pi$, and the residual kernel $Q := (1 - delta)^(-1)(P - delta pi)$ satisfies $P^t (bold(z)^((0)), dot.op) - pi = (1 - delta)^t (Q^t (bold(z)^((0)), dot.op) - pi)$, yielding $||P^t (bold(z)^((0)), dot.op) - pi||_"TV" <= (1 - delta)^t$.
]

#subsec("E", "Document-Length Skew and Breakdown Counterexample")

Because unwindowed extraction generates $binom(L_d, 2) = Theta(L_d^2)$ biterms per document, BTM is highly sensitive to heterogeneous document lengths—a structural vulnerability formalized below.

#proposition(title: "Quadratic Length-Skew Dominance and Zero Asymptotic Breakdown", origin: "Original Counterexample")[
  Let $cal(D)$ consist of $N_s >= 1$ short documents of length $L_s >= 2$ supported on vocabulary subset $cal(V)_1 subset cal(V)$, together with a single long outlier document $d_0$ of length $L_0$ supported on a disjoint subset $cal(V)_2 subset cal(V) without cal(V)_1$. Under full-document extraction ($W = oo$), let $rho := |cal(B)_(d_0)| / |cal(B)|$ be the fraction of corpus biterms generated by $d_0$. Whenever
  $ L_0 >= sqrt(N_s) L_s + 1, $
  we have $rho > 1 / 2$, and for $K = 2$ the unique global MLE of BTM assigns mixing weight $hat(theta)_(z_2) = rho > 1 / 2$ to the single outlier document $d_0$. As $L_0 / (sqrt(N_s) L_s) -> oo$, $hat(theta)_(z_2) -> 1$, proving that unwindowed BTM has finite-sample document-length breakdown point $1 / (N_s + 1) -> 0$.
]

#proof[
  The short documents contribute $M_s = N_s L_s (L_s - 1) / 2$ biterms, while $d_0$ contributes $M_0 = L_0 (L_0 - 1) / 2$ biterms. If $L_0 >= sqrt(N_s) L_s + 1$, then
  $ L_0 (L_0 - 1) >= (sqrt(N_s) L_s + 1)(sqrt(N_s) L_s) > N_s L_s^2 > N_s L_s (L_s - 1), $
  so $M_0 > M_s$ and $rho = M_0 / (M_0 + M_s) > 1 / 2$. Because $cal(V)_1 inter cal(V)_2 = emptyset$, the empirical co-occurrence matrix $C$ of Theorem 1 is block-diagonal, $C = "diag"((1 - rho) C_1, rho C_2)$, with $sum_(u,v in cal(V)_r) (C_r)_(u,v) = 1$. Any $K = 2$ factorization $Q = theta_1 bold(phi)_1 bold(phi)_1^top + theta_2 bold(phi)_2 bold(phi)_2^top$ minimizing $D_"KL" (C || Q)$ puts zero mass on cross-blocks $cal(V)_1 times cal(V)_2$ by allocating disjoint supports $"supp"(bold(phi)_1) subset.eq cal(V)_1$ and $"supp"(bold(phi)_2) subset.eq cal(V)_2$. Differentiating $sum_(u,v in cal(V)_1) C_(u,v) log(theta_1 phi_(1,u) phi_(1,v)) + sum_(u,v in cal(V)_2) C_(u,v) log(theta_2 phi_(2,u) phi_(2,v))$ subject to $theta_1 + theta_2 = 1$ gives $hat(theta)_1 = 1 - rho$ and $hat(theta)_2 = rho$. Hence a single contaminated document of length $L_0 >> sqrt(N_s) L_s$ captures $hat(theta)_2 -> 1$ of the corpus topic mass, yielding breakdown point $1 / N -> 0$.
]

Bounding the sliding window $W = O(1)$ restricts each document to at most $O(W L_d)$ biterms, reducing length influence from quadratic to linear; alternatively, weighting each biterm in $cal(B)_d$ by $|cal(B)_d|^(-1)$ equalizes per-document contribution.

#sec("III", "Complexity, Geometry and Operational Profile")

In each collapsed Gibbs sweep, visiting all $M = |cal(B)|$ biterms and computing $K$ unnormalized topic weights takes $O(K)$ scalar arithmetic operations per biterm, giving total training complexity $O(T K M)$ over $T$ sweeps. With bounded sliding window $W$, $M <= sum_(d=1)^N W L_d = O(N W overline(L))$; without windowing, $M = O(N overline(L^2))$, which explodes quadratically on long documents. Post-training document inference scans each document's biterm multiset $cal(B)_d$ once across $K$ topics in $O(K M)$ total time.

Space complexity is $O(M + K V)$: $O(M)$ integers store the corpus biterm index pairs $(w_(m,1), w_(m,2))$ and current topic labels $z_m$, while $O(K V)$ integers store the dense topic-word count table $n_(z,w)$ and $O(K)$ stores $n_z$. Because BTM updates only _two_ columns ($w_1, w_2$) of $n_(z,w)$ per biterm, storing $n_(z,w)$ in word-major layout (`V x K`) ensures that all $K$ topic counts for $w_1$ and $w_2$ reside in two contiguous cache lines, eliminating cache thrashing during inner loops. When $M$ exceeds RAM, BTM admits online and incremental mini-batch Gibbs variants (oBTM / iBTM @yan2013btm) whose memory scales only with $O(K V)$ plus the current batch size.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T K |cal(B)|)$ training; $O(K |cal(B)|)$ inference ($|cal(B)| = O(N W overline(L))$ windowed, $O(N overline(L^2))$ full).],
    [Space Complexity], [$O(|cal(B)| + K V)$ in-core; $O(K V)$ streaming mini-batch (no $O(N K)$ doc-topic matrix).],
    [Cluster Geometry], [Simplex topics in $Delta_(V-1)$ inducing symmetric rank-$K$ KL factorization of co-occurrence $C$.],
    [Invariances], [Invariant to document ordering and word order within window $W$; sensitive to length skew.],
    [Hyperparameters], [Topic count $K$, Dirichlet priors $(alpha, beta)$, window size $W$, Gibbs sweeps $T$.],
    [Convergence], [Geometrically ergodic Markov chain (Theorem 3); non-convex multimodal posterior.],
  ),
  caption: [Operational and geometric profile of the Biterm Topic Model (BTM).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

BTM excels whenever observations are ultra-short discrete bags of co-occurring symbols where single-sample sparsity defeats document-level latent mixtures:

+ *Microblogging and Social-Media Event Clustering:* Posts on platforms such as X/Twitter or Weibo average 10--20 content tokens. Pooling biterms across millions of posts stabilizes topic-word co-occurrence estimates (Theorem 2), allowing coherent clustering of breaking news, hashtags, and crisis reports without relying on heuristic tweet aggregation @yan2013btm.
+ *Search-Query Intent Taxonomy:* Web search queries typically contain only 2--5 terms, yielding $binom(L, 2) in [1, 10]$ biterms per query. Because LDA's per-query Dirichlet posterior suffers $Omega(1)$ variance at $L <= 5$, BTM's corpus-wide co-occurrence factorization (Theorem 1) directly recovers latent query intents for sponsored-search matching and query rewriting.
+ *Short Clinical Notes and ICD Code Co-Occurrence:* Emergency-room chief complaints, triage notes, and billing-code baskets per patient encounter are short discrete multisets. Modeling pairwise symptom-diagnosis biterms uncovers clinical phenotypes and comorbidity topics while naturally assigning soft multi-topic memberships $P(z | d)$ to multi-morbid patients.
+ *E-Commerce Product Titles and Market-Basket Profiling:* Catalog titles and small shopping carts contain sparse co-occurring attribute tokens or SKU IDs. By factoring the symmetric pairwise co-occurrence graph $C$ under KL divergence, BTM clusters products into coherent semantic departments even when individual descriptions share few exact tokens.
+ *Bioinformatics Short-Read $k$-mer and Motif Clustering:* Metagenomic amplicon reads and short regulatory sequences can be tokenized into overlapping $k$-mers within a local window $W$. BTM captures conserved local $k$-mer co-occurrence patterns across millions of short reads without requiring a high-dimensional per-read mixture parameter.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

BTM exhibits four distinct failure regimes: (i) *Long or heterogeneous documents* suffer quadratic biterm blow-up ($O(L_d^2)$) and outlier dominance (Proposition 1) unless windowed or reweighted; (ii) *Continuous dense tabular data* lack sparse discrete co-occurrence structure—naive discretization or non-negative feature scaling creates near-complete biterm graphs that collapse posterior topics; (iii) *Very large vocabularies without pruning* inflate $O(K V)$ memory and dilute $V beta$ smoothing; and (iv) *Topic-correlation structure* is not modeled because $bold(theta)$ draws independent per-biterm topics rather than logistic-normal correlations @blei2007ctm.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Tweets & Queries ($L <= 25$)], [Optimal: eliminates $N(K-1)$ sparse nuisance parameters.], [Prune stop-words; compare GSDMM @yin2014dmm.],
    [Long Articles ($L >= 500$)], [Poor unwindowed: $O(L^2)$ biterm blow-up (Prop. 1).], [Enforce window $W <= 15$ or use LDA @blei2003lda.],
    [E-Commerce & Baskets], [Strong: KL factorization of pairwise co-occurrence.], [Down-weight mega-baskets by $|cal(B)_d|^(-1)$.],
    [Dense Continuous Tabular], [Mismatch: dense co-occurrence induces topic collapse.], [Use GMM / MCLUST @fraley2002mclust or SNOB @wallace1968snob.],
    [Unknown Topic Count $K$], [Fixed-$K$ parametric model; requires external selection.], [Use HDP @teh2006hdp or held-out perplexity.],
  ),
  caption: [Domain suitability and failure-mode matrix for BTM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

BTM occupies a principled bridge between unigram document-mixture models (pLSA @hofmann1999plsa, LDA @blei2003lda), single-topic document mixtures (GSDMM @yin2014dmm), and symmetric nonnegative matrix factorization @lee1999nmf @ding2008nmfplsa.

#proposition(title: "Marginal Unigram Consistency and Reduction to GSDMM", origin: "Original Derivation")[
  Let $Q = bold(Phi)^top "diag"(bold(theta)) bold(Phi)$ be BTM's symmetric biterm distribution (Theorem 1).
  + *Marginal Unigram Equivalence:* Summing $Q_(u,v)$ over $v in [V]$ recovers the corpus-level unigram topic mixture $q_u := sum_(v=1)^V Q_(u,v) = sum_(z=1)^K theta_z phi_(z,u)$, matching the marginal corpus word distribution of pLSA and LDA while preserving second-order off-diagonal coupling $Q_(u,v) - q_u q_v = "Cov"_(z tilde.op bold(theta)) (phi_(z,u), phi_(z,v))$.
  + *Tied-Biterm Reduction to GSDMM:* Restricting all biterms within document $d$ to share a single latent topic $z_m = c_d$ for all $m in cal(B)_d$ transforms BTM's conditional document likelihood into
    $ P(cal(B)_d | c_d = z, bold(Phi)) = product_(i=1)^(L_d) phi_(z, w_(d,i))^(L_d - 1), $
    which equals the GSDMM @yin2014dmm document likelihood raised to the length-dependent power $L_d - 1$.
]

#proof[
  For Claim 1, using $sum_(v=1)^V phi_(z,v) = 1$,
  $ sum_(v=1)^V Q_(u,v) = sum_(z=1)^K theta_z phi_(z,u) sum_(v=1)^V phi_(z,v) = sum_(z=1)^K theta_z phi_(z,u) = bb(E)_(z)[phi_(z,u)]. $
  Thus $Q_(u,v) - q_u q_v = bb(E)_z [phi_(z,u) phi_(z,v)] - bb(E)_z [phi_(z,u)] bb(E)_z [phi_(z,v)] = "Cov"_(z tilde.op bold(theta)) (phi_(z,u), phi_(z,v))$. For Claim 2, under full-document extraction ($W = L_d$), each token occurrence $w_(d,i)$ participates in exactly $L_d - 1$ pairs $\{w_(d,i), w_(d,j)\}$ ($j != i$). Tying $z_m = c_d = z$ across all $binom(L_d, 2)$ biterms in $cal(B)_d$ yields
  $ P(cal(B)_d | c_d = z, bold(Phi)) = product_(1 <= i < j <= L_d) phi_(z, w_(d,i)) phi_(z, w_(d,j)) = product_(i=1)^(L_d) phi_(z, w_(d,i))^(L_d - 1), $
  which is the $(L_d - 1)$-th power of the unigram multinomial likelihood $product_(i=1)^(L_d) phi_(z, w_(d,i))$ used in GSDMM @yin2014dmm.
]

*Empirical Behavior on Continuous Tabular Benchmarks.* The ClusterOverview 4-dataset benchmark provides an instructive stress test of BTM outside its native sparse-text regime (evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing):
- *Iris:* $"Sil" = "N/A"$, $"ARI" = 0.000$.
- *Customers:* $"Sil" = 0.150$.
- *Multi-Shape:* $"Sil" = 0.344$, $"ARI" = "N/A"$.
- *High-D:* $"Sil" = 0.162$, $"ARI" = "N/A"$.
- *Mean Runtime:* $6.008$ ms.

These results directly mirror our theoretical analysis. Because BTM requires sparse pairwise co-occurrences over a discrete vocabulary, adapting it to dense, low-dimensional continuous tabular data like Iris ($d = 4$) yields nearly uniform pairwise co-occurrence graphs in which every feature pair co-occurs in every sample; consequently, the symmetric rank-$K$ KL factorization (Theorem 1) collapses onto a single dominant topic, leaving Silhouette undefined (`N/A`) and `"ARI" = 0.000`. On datasets exhibiting richer positive multi-feature contrasts (Customers, Multi-Shape, High-D), BTM extracts nontrivial cluster partitions (`Sil` up to $0.344$), though Euclidean separation remains modest relative to metric or Gaussian mixture models @jain2010data @wallace1968snob. Conversely, its $6.008$ ms mean runtime confirms the computational efficiency of corpus-level biterm pooling without per-document iterative inference.
