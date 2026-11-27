#import "../template/ieee.typ": *

#alg_entry(
  id: 84,
  name: "Probabilistic Latent Semantic Analysis (pLSA)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "1999",
  paradigm: "Aspect Model / Multinomial Matrix Factorization via EM",
  complexity: $cal(O)(T K N_"nnz")$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Probabilistic Latent Semantic Analysis (pLSA)—originally introduced as Probabilistic Latent Semantic Indexing (pLSI) @hofmann1999plsa—provides a statistical foundation for Latent Semantic Analysis (LSA) by replacing algebraic truncated SVD with a discrete latent-variable generative model of co-occurrence data. Given a non-negative count matrix $bold(N) = (n(d, w)) in bb(R)_(>= 0)^(M times V)$ recording co-occurrences between $M$ documents $d in cal(D) = \{1, dots.h, M\}$ and $V$ terms $w in cal(V) = \{1, dots.h, V\}$, let $n(d) = sum_(w=1)^V n(d, w) > 0$, total count $N_"tot" = sum_(d=1)^M n(d)$, and sparse support $Omega = \{(d, w) : n(d, w) > 0\}$ with $|Omega| = N_"nnz"$.

pLSA posits the *aspect model*: a discrete latent variable (aspect or topic) $z in \{z_1, dots.h, z_K\}$ renders documents and words conditionally independent. Two canonical parameterizations exist @hofmann1999plsa:
1. *Asymmetric form:* $P(d, w) = P(d) sum_(k=1)^K P(w | z_k) P(z_k | d)$, with parameters $bold(Theta) = (P(z_k | d)) in bb(R)^(M times K)$ (rows in $Delta^(K-1)$) and $bold(B) = (P(w | z_k)) in bb(R)^(K times V)$ (rows in $Delta^(V-1)$).
2. *Symmetric form:* $P(d, w) = sum_(k=1)^K P(z_k) P(d | z_k) P(w | z_k)$, with $bold(p)_Z in Delta^(K-1)$, $bold(U) = (P(d | z_k)) in bb(R)^(M times K)$ (columns in $Delta^(M-1)$), and $bold(V) = (P(w | z_k)) in bb(R)^(V times K)$ (columns in $Delta^(V-1)$).

Maximizing the multinomial log-likelihood of the observed co-occurrence counts reduces (up to the empirical document-marginal term $sum_d n(d) log P(d)$) to maximizing
$
  ell(bold(Theta), bold(B)) = sum_(d=1)^M sum_(w=1)^V n(d, w) log (sum_(k=1)^K P(z_k | d) P(w | z_k)) .
$
Because direct MLE over $Omega(M K)$ parameters readily overfits finite corpora, Hofmann @hofmann1999plsa introduced Tempered Expectation-Maximization (TEM), a deterministic annealing generalization of EM controlled by an inverse temperature $beta in (0, 1]$.

#definition(title: "Tempered EM Variational Free Energy for pLSA")[
  Let $Q = \{Q(z_k | d, w)\}_((d, w) in Omega)$ denote auxiliary posterior distributions over $z_k$ for each observed pair $(d, w) in Omega$. At inverse temperature $beta in (0, 1]$, the tempered free energy is
  $
    cal(F)_beta (Q, bold(Theta), bold(B)) &= -sum_((d, w) in Omega) n(d, w) sum_(k=1)^K Q(z_k | d, w) log [P(z_k | d) P(w | z_k)] \
    &quad + 1 / beta sum_((d, w) in Omega) n(d, w) sum_(k=1)^K Q(z_k | d, w) log Q(z_k | d, w) .
  $
  At $beta = 1$, $-cal(F)_1 (Q, bold(Theta), bold(B))$ equals the standard EM variational lower bound on $ell(bold(Theta), bold(B))$; for $beta < 1$, the amplified entropy term smooths latent posteriors to curb overfitting.
]

#algbox(
  title: "EM and Tempered EM (TEM) Algorithm for pLSA",
  inputs: [Sparse count matrix $bold(N) in bb(R)_(>= 0)^(M times V)$ with support $Omega$, topics $K$, annealing factor $eta in (0, 1]$, tolerance $tau > 0$.],
  outputs: [Topic-word matrix $hat(bold(B)) = (P(w | z_k))$, document-topic weights $hat(bold(Theta)) = (P(z_k | d))$.],
)[
  1. *Initialize* $bold(Theta)^((0)) > bold(0)$ and $bold(B)^((0)) > bold(0)$ with strictly positive simplex rows; set $beta <- 1$.
  2. *For* $t = 0, 1, dots.h, T - 1$:
  3. #h(1em) *Tempered E-step:* For every nonzero pair $(d, w) in Omega$ and topic $k in \{1, dots.h, K\}$, compute
     $
       Q(z_k | d, w) = ([P^((t))(z_k | d) P^((t))(w | z_k)]^beta) / (sum_(j=1)^K [P^((t))(z_j | d) P^((t))(w | z_j)]^beta) .
     $
  4. #h(1em) *M-step:* Accumulate expected topic counts over $(d, w) in Omega$ and normalize:
     $
       P^((t+1))(w | z_k) &= (sum_(d : (d, w) in Omega) n(d, w) Q(z_k | d, w)) / (sum_((d', w') in Omega) n(d', w') Q(z_k | d', w')) , \
       P^((t+1))(z_k | d) &= (sum_(w : (d, w) in Omega) n(d, w) Q(z_k | d, w)) / n(d) .
     $
  5. #h(1em) *Convergence / Tempering:* If $beta = 1$ and $|ell^((t+1)) - ell^((t))| < tau$, *break*; in TEM mode, if held-out perplexity stops decreasing, anneal $beta <- eta beta$ and stop when no further validation gain occurs.
  6. *Return* $hat(bold(B))$ and $hat(bold(Theta))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove the algebraic equivalence of symmetric and asymmetric aspect models, the exact KL-divergence equivalence between pLSA and Non-Negative Matrix Factorization (NMF), EM monotonicity alongside simplex-face absorption, a constructive finite-sample overfitting lower bound, and global concavity of out-of-sample fold-in projection.

#subsec("A", "Symmetric-Asymmetric Parameterization Equivalence")

#theorem(title: "Exact Model Equivalence of Symmetric and Asymmetric pLSA", origin: "Original Proof")[
  Let $cal(M)_K^"asym"$ and $cal(M)_K^"sym"$ denote the families of joint distributions $bold(P) in Delta^(M V - 1)$ expressible in the asymmetric and symmetric aspect forms, respectively, with strictly positive marginals $P(d) > 0$ and $P(z_k) > 0$. Then $cal(M)_K^"asym" = cal(M)_K^"sym"$, both models have identical parameter dimension $p_"pLSA" = M K + K V - K - 1$, and at every MLE of the joint log-likelihood the document marginal satisfies $hat(P)(d) = n(d) / N_"tot"$ identically.
]

#proof[
  Given any asymmetric model $(bold(p)_D, bold(Theta), bold(B)) in cal(M)_K^"asym"$ with $p_(D, d) = P(d) > 0$, define the symmetric parameters via Bayes' rule:
  $
    P(z_k) &:= sum_(d=1)^M P(d) P(z_k | d) , \
    P(d | z_k) &:= (P(d) P(z_k | d)) / P(z_k) ,
  $
  keeping $P(w | z_k)$ unchanged. Clearly $sum_k P(z_k) = sum_d P(d) (sum_k P(z_k | d)) = 1$ and $sum_d P(d | z_k) = 1$, while
  $
    sum_(k=1)^K P(z_k) P(d | z_k) P(w | z_k) = P(d) sum_(k=1)^K P(z_k | d) P(w | z_k) ,
  $
  showing $cal(M)_K^"asym" subset.eq cal(M)_K^"sym"$. Conversely, given any symmetric model $(bold(p)_Z, bold(U), bold(V)) in cal(M)_K^"sym"$, define $P(d) := sum_k P(z_k) P(d | z_k)$ and $P(z_k | d) := (P(z_k) P(d | z_k)) / P(d)$, which obey $sum_d P(d) = 1$ and $sum_k P(z_k | d) = 1$ and reproduce the identical joint probability $P(d, w)$, proving $cal(M)_K^"sym" subset.eq cal(M)_K^"asym"$.

  Counting independent degrees of freedom gives $(M - 1) + M (K - 1) + K (V - 1) = M K + K V - K - 1$ for the asymmetric form, and $(K - 1) + K (M - 1) + K (V - 1) = M K + K V - K - 1$ for the symmetric form. Finally, the joint log-likelihood decouples as
  $
    sum_(d, w) n(d, w) log P(d, w) = sum_(d=1)^M n(d) log P(d) + ell(bold(Theta), bold(B)) .
  $
  Because $bold(p)_D in Delta^(M-1)$ enters only the first strictly concave term, Lagrange multipliers yield the unique global maximizer $hat(P)(d) = n(d) / N_"tot"$ independently of $(bold(Theta), bold(B))$.
]

#subsec("B", "Equivalence to Kullback-Leibler NMF")

#theorem(title: "Objective and Stationary-Point Equivalence of pLSA and KL-NMF", origin: "Classical — Cited")[
  Let $hat(bold(P)) = bold(N) / N_"tot" in Delta^(M V - 1)$ be the empirical normalized co-occurrence matrix. Maximizing the joint pLSA log-likelihood over $bold(P) = bold(U) bold(Lambda) bold(V)^top in cal(M)_K^"sym"$ is strictly equivalent to minimizing the generalized KL (I-divergence) Non-Negative Matrix Factorization objective @lee2001algorithms @ding2008nmfplsa:
  $
    min_(bold(W) in bb(R)_(>= 0)^(M times K), bold(H) in bb(R)_(>= 0)^(K times V)) D_I (hat(bold(P)) || bold(W) bold(H)) ,
  $
  where $D_I (bold(A) || bold(C)) = sum_(d, w) [A_(d w) log(A_(d w) / C_(d w)) - A_(d w) + C_(d w)]$. In particular, every stationary point $(bold(W)^*, bold(H)^*)$ of KL-NMF satisfies $sum_(d, w) (bold(W)^* bold(H)^*)_(d w) = 1$ and corresponds bijectively (up to diagonal rescaling) to a stationary point of pLSA-EM.
]

#proof[
  Consider any nonnegative factorization $bold(C) = bold(W) bold(H)$ with $bold(W) >= bold(0)$, $bold(H) >= bold(0)$ and nonzero columns/rows. Define column sums $s_k = sum_(d=1)^M W_(d k) > 0$ and row sums $r_k = sum_(w=1)^V H_(k w) > 0$, and normalize $U_(d k) = W_(d k) / s_k$, $V_(w k) = H_(k w) / r_k$, and $lambda_k = s_k r_k > 0$. Setting total mass $Z = sum_(k=1)^K lambda_k = sum_(d, w) (bold(W) bold(H))_(d w)$ and $P(z_k) = lambda_k / Z$ decomposes $bold(W) bold(H) = Z dot bold(P)$ with $bold(P) in cal(M)_K^"sym"$ and $sum_(d, w) P_(d w) = 1$. Substituting $bold(C) = Z bold(P)$ into $D_I (hat(bold(P)) || Z bold(P))$ gives
  $
    D_I (hat(bold(P)) || Z bold(P)) &= sum_(d, w) hat(P)_(d w) log hat(P)_(d w) - sum_(d, w) hat(P)_(d w) log P_(d w) \
    &quad - log Z - 1 + Z .
  $
  Minimizing the convex scalar map $Z |-> Z - log Z - 1$ over $Z > 0$ gives the unique global minimizer $Z^* = 1$ with minimum value $0$! At $Z = 1$,
  $
    D_I (hat(bold(P)) || bold(P)) = -H(hat(bold(P))) - 1 / N_"tot" sum_(d, w) n(d, w) log P_(d w) ,
  $
  where $-H(hat(bold(P)))$ is constant. Thus minimizing $D_I (hat(bold(P)) || bold(W) bold(H))$ over unconstrained nonnegative factors is identical to maximizing the pLSA log-likelihood over $bold(P) in cal(M)_K^"sym"$, and every KKT stationary point of KL-NMF automatically satisfies $Z^* = sum_(d, w) (bold(W)^* bold(H)^*)_(d w) = 1$ @ding2008nmfplsa. Matching the KKT conditions $W_(d k) (sum_w (hat(P)_(d w) / (bold(W) bold(H))_(d w)) H_(k w) - sum_w H_(k w)) = 0$ with the pLSA-EM fixed-point equations establishes identical stationary sets.
]

#subsec("C", "EM Monotonicity and Simplex-Boundary Absorption")

#lemma(title: "Monotonic Ascent and Boundary Absorption of pLSA-EM", origin: "Original Proof")[
  For standard pLSA-EM ($beta = 1$) initialized with $bold(Theta)^((0)) > bold(0)$ and $bold(B)^((0)) > bold(0)$:
  1. *(Strict ascent):* At every iteration $t >= 0$,
     $
       &ell(bold(Theta)^((t+1)), bold(B)^((t+1))) - ell(bold(Theta)^((t)), bold(B)^((t))) \
       &>= sum_((d, w) in Omega) n(d, w) D_"KL" (Q^((t))_(d w) || Q^((t+1))_(d w)) >= 0 .
     $
  2. *(Boundary absorption):* Every coordinate face $P(w | z_k) = 0$ or $P(z_k | d) = 0$ is strictly invariant under pLSA-EM: if $P^((t))(w | z_k) = 0$ (or $P^((t))(z_k | d) = 0$) at any $t$, then the same coordinate remains identically $0$ for all $s >= t$.
]

#proof[
  *Part 1.* At $beta = 1$, Definition 1 gives $ell(bold(Theta), bold(B)) = -cal(F)_1 (Q, bold(Theta), bold(B)) + sum_((d, w) in Omega) n(d, w) D_"KL" (Q_(d w) || P(bold(z) | d, w; bold(Theta), bold(B)))$. Setting $Q^((t))_(d w) (z_k) = P(z_k | d, w; bold(Theta)^((t)), bold(B)^((t)))$ zeros the KL gap so $-cal(F)_1 (Q^((t)), bold(Theta)^((t)), bold(B)^((t))) = ell(bold(Theta)^((t)), bold(B)^((t)))$. Fixing $Q^((t))$, maximizing $-cal(F)_1$ decouples into $K$ concave multinomial problems over $bold(beta)_k in Delta^(V-1)$ and $M$ concave multinomial problems over $bold(theta)_d in Delta^(K-1)$:
  $
    max_(bold(theta)_d in Delta^(K-1)) sum_(k=1)^K (sum_w n(d, w) Q^((t)) (z_k | d, w)) log theta_(d, k) ,
  $
  and analogously for $bold(beta)_k$. Solving the Lagrange conditions uniquely yields the M-step formulas, guaranteeing $-cal(F)_1 (Q^((t)), bold(Theta)^((t+1)), bold(B)^((t+1))) >= -cal(F)_1 (Q^((t)), bold(Theta)^((t)), bold(B)^((t)))$. Adding the non-negative KL term at $(bold(Theta)^((t+1)), bold(B)^((t+1)))$ yields claim 1, and since $ell <= 0$ is bounded above, convergence to a stationary value follows @dempster1977em @wu1983convergence.

  *Part 2.* Suppose $P^((t))(w_0 | z_(k_0)) = 0$. In the E-step at iteration $t$, the numerator for $(d, w_0, k_0)$ is $P^((t))(z_(k_0) | d) P^((t))(w_0 | z_(k_0)) = 0$ for every document $d$, hence $Q^((t))(z_(k_0) | d, w_0) = 0$ for all $d$. Substituting into the M-step gives
  $
    P^((t+1))(w_0 | z_(k_0)) = (sum_d n(d, w_0) dot 0) / (sum_(d', w') n(d', w') Q^((t))(z_(k_0) | d', w')) = 0 .
  $
  An identical argument holds whenever $P^((t))(z_(k_0) | d_0) = 0$: then $Q^((t))(z_(k_0) | d_0, w) = 0$ for all $w$, forcing $P^((t+1))(z_(k_0) | d_0) = 0$. Induction over $s >= t$ establishes exact absorption on all boundary faces.
]

#subsec("D", "Overfitting Divergence and Parameter-to-Sample Scaling")

#theorem(title: "Non-Vanishing Parameter Ratio and Exact Interpolation Witness", origin: "Original Proof")[
  Let $bar(N) = N_"tot" / M$ denote the mean document length of a corpus of $M$ documents.
  1. *(Persistent parameter ratio):* Under unregularized pLSA, the ratio of free parameters to total observed tokens satisfies
     $
       p_"pLSA" / N_"tot" = (M K + K V - K - 1) / (M bar(N)) > K / bar(N) - 1 / (M bar(N)) ,
     $
     which remains bounded below by $K / bar(N) > 0$ even as corpus size $M -> infinity$ for fixed $bar(N)$.
  2. *(Interpolation degeneracy):* Whenever $K >= M$, the global MLE achieves exact empirical memorization $hat(P)(w | d) = n(d, w) / n(d)$ with zero training I-divergence, while assigning exact probability $0$ (infinite held-out perplexity) to any unseen word $w$ in document $d$.
]

#proof[
  *Part 1.* Substituting $p_"pLSA" = M K + K V - K - 1$ from Theorem 1 and $N_"tot" = M bar(N)$ gives
  $
    p_"pLSA" / N_"tot" = K / bar(N) + (K (V - 1) - 1) / (M bar(N)) .
  $
  For $K >= 1, V >= 2$, $K (V - 1) >= 1$, so $p_"pLSA" / N_"tot" >= K / bar(N)$ for all $M >= 1$. Consequently, unlike parametric models with fixed dimension (such as LDA, where $p_"LDA" / N_"tot" = K V / (M bar(N)) -> 0$ as $M -> infinity$), pLSA never enters the classical fixed-dimension asymptotic regime as more documents arrive @blei2003lda.

  *Part 2.* Suppose $K >= M$. Set $P(z_k | d) = delta_(k d)$ for $d in \{1, dots.h, M\}$ and $P(w | z_d) = n(d, w) / n(d)$ for all $w in cal(V)$ (with arbitrary simplex vectors for $k > M$). Then for every document $d$ and word $w$,
  $
    sum_(k=1)^K P(z_k | d) P(w | z_k) = P(w | z_d) = n(d, w) / n(d) .
  $
  Because $D_I (hat(bold(P)) || bold(P)) >= 0$ with equality iff $bold(P) = hat(bold(P))$, this memorization solution attains the global maximum of $ell(bold(Theta), bold(B))$ (Theorem 2). However, for any vocabulary term $w^*$ not occurring in document $d$ ($n(d, w^*) = 0$), the fitted model predicts $hat(P)(w^* | d) = 0$, causing the log-likelihood of any held-out continuation of document $d$ containing $w^*$ to diverge to $-infinity$ (perplexity $+infinity$).
]

#subsec("E", "Global Concavity of Out-of-Sample Fold-In Projection")

#proposition(title: "Global Concavity of pLSA Document Fold-In", origin: "Original Derivation")[
  Given a fixed topic-word matrix $bold(B) in bb(R)_(> 0)^(K times V)$ learned from a training corpus, folding in an unseen document with count vector $bold(n)_"new" in bb(R)_(>= 0)^V$ by maximizing
  $
    f(bold(theta)) = sum_(w=1)^V n_"new" (w) log (sum_(k=1)^K theta_k B_(k w)) quad "subject to" quad bold(theta) in Delta^(K-1)
  $
  is a globally concave optimization problem. Consequently, holding $bold(B)$ fixed in Algorithm 84 yields fold-in EM iterates that converge to a global maximizer of $f$ from any strictly positive initialization $bold(theta)^((0)) > bold(0)$.
]

#proof[
  Let $p_w (bold(theta)) = sum_(k=1)^K theta_k B_(k w) > 0$ for all $bold(theta) in Delta^(K-1)$ since $bold(B) > bold(0)$. Differentiating $f$ twice with respect to $bold(theta)$ yields gradient $(nabla f(bold(theta)))_j = sum_w n_"new" (w) B_(j w) / p_w (bold(theta))$ and Hessian entries
  $
    (nabla^2 f(bold(theta)))_(j l) = -sum_(w=1)^V n_"new" (w) (B_(j w) B_(l w)) / (p_w (bold(theta)))^2 .
  $
  For any vector $bold(v) in bb(R)^K$, the quadratic form equals
  $
    bold(v)^top nabla^2 f(bold(theta)) bold(v) = -sum_(w=1)^V (n_"new" (w)) / (p_w (bold(theta)))^2 (sum_(k=1)^K v_k B_(k w))^2 <= 0 ,
  $
  establishing that $nabla^2 f(bold(theta)) prec.eq bold(0)$ everywhere on $Delta^(K-1)$. Because $Delta^(K-1)$ is compact and convex and $f$ is concave, every KKT stationary point of $f$ on $Delta^(K-1)$ is a global maximizer. By Lemma 1, fold-in EM is monotonically non-decreasing; because $bold(theta)^((0)) > bold(0)$ stays strictly positive until converging to the boundary only where KKT slackness holds, every limit point satisfies the concave KKT conditions and attains $max_(bold(theta) in Delta^(K-1)) f(bold(theta))$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational Complexity and Speed Advantage.* Unlike LDA (Alg. 83), pLSA has *no inner variational loop* and requires *no transcendental digamma or exponential evaluations* when $beta = 1$. Each EM iteration visits only the $N_"nnz"$ nonzero entries of the sparse matrix $bold(N)$: the E-step computes $K$ scalar products $P(z_k | d) P(w | z_k)$ and normalizes in $cal(O)(K N_"nnz")$ flops, and the M-step accumulates weighted sums into $K V$ word slots and $M K$ document slots in $cal(O)(K N_"nnz")$ flops. Total time over $T$ iterations is $cal(O)(T K N_"nnz")$ and memory is $cal(O)(K (M + V) + N_"nnz")$—streaming posteriors $Q(z_k | d, w)$ on the fly avoids storing an $N_"nnz" times K$ tensor altogether.

*Simplex-Conical Geometry vs. SVD.* Truncated SVD / LSA approximates $bold(N)$ in Frobenius norm using orthogonal singular vectors with mixed positive and negative signs, hindering probabilistic interpretation. By contrast, pLSA minimizes the KL / I-divergence (Theorem 2) over a *non-negative conical factorization* $hat(bold(P)) approx "diag"(bold(p)_D) bold(Theta) bold(B)$, representing every normalized document vector inside the convex polytope $"conv"(\{bold(beta)_1, dots.h, bold(beta)_K\}) subset.eq Delta^(V-1)$. Because there is no Dirichlet prior pulling $bold(theta)_d$ toward the barycenter, pLSA can represent sharp sparse vertex allocations cleanly when tempered or early-stopped.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Objective], [Maximize $ell = sum_(d, w) n(d, w) log sum_k P(z_k | d) P(w | z_k)$],
    [Equivalent Loss], [Generalized KL / I-divergence NMF $D_I (hat(bold(P)) || bold(W) bold(H))$],
    [Free Parameters], [$M K + K V - K - 1$ (grows linearly with $M$)],
    [Time Complexity], [$cal(O)(T K N_"nnz")$ (single-loop multiply-add; no digamma)],
    [Space Complexity], [$cal(O)(K (M + V) + N_"nnz")$ with on-the-fly E-step],
    [Cluster Geometry], [Non-negative low-rank simplex cone inside $Delta^(V-1)$],
    [Boundary Behavior], [Zero coordinates are absorbing (Lemma 1); init $> 0$],
    [Regularization], [Tempered EM ($beta < 1$) or held-out early stopping],
  ),
  caption: [Operational and theoretical profile of pLSA (Alg. 84).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Information Retrieval and Ad-Hoc Query Disambiguation (pLSI):* Classical keyword retrieval suffers from synonymy (distinct words expressing one concept) and polysemy (one word spanning multiple senses). By folding queries into the low-dimensional aspect space via concave projection (Proposition 1) and smoothing $P(w | d)$ over $K$ topics, pLSA matches semantically related documents without exact lexical overlap @hofmann1999plsa.
2. *Implicit-Feedback Collaborative Filtering:* User-item interaction tables (clicks, plays, purchases) are vast, sparse, non-negative count matrices. The symmetric aspect form $P(u, i) = sum_k P(z_k) P(u | z_k) P(i | z_k)$ models latent taste communities directly in $cal(O)(T K N_"nnz")$ time without fabricating negative ratings.
3. *Visual Bag-of-Words Scene Categorization:* Quantizing local image descriptors (SIFT or patch tokens) into visual words transforms image collections into sparse image-word histograms. pLSA separates foreground object parts from recurring background textures through soft aspect responsibilities $Q(z_k | d, w)$ per spatial patch.
4. *Web Search-Log and Clickstream Intent Clustering:* Query-URL click graphs form bipartite co-occurrence matrices with millions of sparse entries. Streaming pLSA factorizes user search intents in linear pass time over $N_"nnz"$, uncovering latent informational vs. navigational session clusters.
5. *Compositional Microbiome and Gene-Condition Profiling:* Metagenomic species-abundance tables and gene-expression count matrices are inherently non-negative and compositional (proportional across samples). Minimizing I-divergence (Theorem 2) aligns closely with Poisson/multinomial sequencing noise, decomposing samples into interpretable microbial enterotypes or transcriptional programs.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four operational pitfalls must be managed: (i) *overfitting without tempering*—as proved in Theorem 3, unregularized EM with large $K$ memorizes training counts and drives unseen-word likelihood to zero (apply TEM $beta < 1$, validation early stopping, or LDA Alg. 83 @blei2003lda); (ii) *lack of a generative document prior*—pLSA assigns no prior probability to unseen documents without an auxiliary fold-in step; (iii) *non-identifiability and local optima*—the simplex polytope $"conv"(\{bold(beta)_k\})$ is generally non-unique without anchor-word separability, and Lemma 1's boundary absorption locks premature zeros permanently; and (iv) *negative or strongly correlated continuous data*—Euclidean metric geometry with negative signs violates the multinomial count assumption (use GMM Alg. 81).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Large sparse co-occurrence matrices], [Fast $cal(O)(T K N_"nnz")$ KL-NMF factorization], [Use TEM ($beta < 1$) or validation early stop],
    [Strict generative new-doc scoring], [No prior on $bold(theta)_"new"$; fold-in only], [LDA (Alg. 83) @blei2003lda or HDP @teh2006hdp],
    [Ultra-short texts / microblogs], [Severe memorization risk ($K / bar(N)$ large)], [BTM @yan2013btm or DMM @yin2014dmm],
    [Continuous correlated Gaussian clouds], [Loses covariance rotation and radial scale], [GMM (Alg. 81) or DPGMM (Alg. 82)],
    [Dense signed spectral embeddings], [Requires non-negative matrix entries], [Spectral clustering @ng2002spectral or LSA/SVD],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for pLSA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Unified Hub: Latent Class Analysis, pLSA, KL-NMF, and LDA", origin: "Original Unified Proof")[
  Let $hat(bold(P)) = bold(N) / N_"tot"$ be a normalized contingency table. Across Category IV and matrix-factorization clustering:
  1. pLSA is the two-way dyadic specialization of Lazarsfeld's Latent Class Analysis (LCA) @lazarsfeld1968latent on $(d, w)$ pairs;
  2. pLSA and KL-NMF @lee2001algorithms @ding2008nmfplsa minimize the identical divergence $D_I (hat(bold(P)) || bold(W) bold(H))$ (Theorem 2), with pLSA preserving probability normalization at every step;
  3. Placing $bold(theta)_d tilde.op "Dir"(bold(alpha))$ and marginalizing $bold(theta)_d$ converts pLSA into LDA (Alg. 83 @blei2003lda), reducing parameter growth from $Omega(M K)$ to $cal(O)(K V)$ while strictly bounding empirical overfitting.
]

#proof[
  In Lazarsfeld's LCA @lazarsfeld1968latent, $J$ observed categorical attributes $(X_1, dots.h, X_J)$ are conditionally independent given a discrete latent class $Z in \{1, dots.h, K\}$: $P(X_1, dots.h, X_J) = sum_(k=1)^K P(Z = k) product_(j=1)^J P(X_j | Z = k)$. Setting $J = 2$ with $X_1 = d in cal(D)$ and $X_2 = w in cal(V)$ reproduces the exact symmetric pLSA model $P(d, w) = sum_k P(z_k) P(d | z_k) P(w | z_k)$. Claim 2 follows directly from Theorem 2, where the KL-NMF multiplicative updates $W_(d k) <- W_(d k) (sum_w (hat(P)_(d w) / (bold(W) bold(H))_(d w)) H_(k w)) / (sum_w H_(k w))$ coincide with pLSA's E- and M-step composition whenever rows of $bold(H)$ are renormalized to sum to $1$ after each update. Claim 3 follows from comparing Theorem 3 ($p_"pLSA" = M K + K V - K - 1$) with LDA's hierarchical marginalization over $bold(theta)_d tilde.op "Dir"(bold(alpha))$ (Alg. 83), whose fixed global parameter count $K V$ eliminates the $M K$ degrees of freedom responsible for Theorem 3's interpolation degeneracy.
]

*Empirical Benchmark Behavior (ClusterOverview 4-Dataset Suite).* On the standardized four-dataset evaluation, pLSA records: *Iris* Silhouette $= 0.298$, ARI $= 0.507$; *Customers* Silhouette $= 0.311$; *Multi-Shape* Silhouette $= 0.706$ (ARI N/A); *High-D* Silhouette $= 0.287$ (ARI N/A); and mean runtime $2.291$ ms. Three sharp empirical findings validate the theory above. First, pLSA achieves the *highest High-D Silhouette ($0.287$) among all four Category IV models* (surpassing LDA at $0.232$, DPGMM at $0.208$, and GMM at $0.181$): in high-dimensional positive-scaled spaces, pLSA's unsmoothed low-rank KL-NMF decomposition (Theorem 2) extracts crisp sparse additive components without either the $cal(O)(K d^2)$ covariance noise of Gaussian mixtures or the Dirichlet barycentric shrinkage of LDA. Second, on *Iris* and *Customers*, pLSA consistently outperforms LDA (Iris ARI $0.507$ vs. $0.429$; Customers Sil $0.311$ vs. $0.245$) while matching LDA identically on *Multi-Shape* ($0.706$ vs. $0.706$), showing that unregularized document-topic weights adapt more flexibly on low-dimensional tabular profiles—though still trailing full-covariance GMM on *Iris* (ARI $0.950$) where ellipsoidal correlation dominates. Third, its $2.291$ ms mean runtime is $10.7 times$ faster than LDA ($24.476$ ms) and $2.25 times$ faster than DPGMM ($5.162$ ms), confirming the exact single-loop $cal(O)(T K N_"nnz")$ efficiency advantage analyzed in Section III.
