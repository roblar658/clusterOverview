#import "../template/ieee.typ": *

#alg_entry(
  id: 83,
  name: "Latent Dirichlet Allocation (LDA)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2003",
  paradigm: "Hierarchical Bayesian Generative Topic Model",
  complexity: $cal(O)(T_"out" (T_"in" K N_"nnz" + K V))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Latent Dirichlet Allocation (LDA) @blei2003lda is a three-level hierarchical Bayesian generative model for mixed-membership (soft) clustering of discrete count data—most prominently text corpora, where documents simultaneously express multiple latent thematic topics. Unlike single-assignment mixture models (Algs. 81--82) in which each observation belongs to one latent component, LDA permits every document to draw a continuous topic-proportion vector on the simplex $Delta^(K-1)$ and allows individual tokens within the same document to originate from different topics.

Consider a corpus $cal(D) = \{bold(w)_1, dots.h, bold(w)_M\}$ of $M$ documents over a discrete vocabulary $cal(V) = \{1, dots.h, V\}$. Document $d$ consists of $N_d$ word tokens $bold(w)_d = (w_(d, 1), dots.h, w_(d, N_d))$ with $w_(d, n) in cal(V)$ and total corpus length $N_"tot" = sum_(d=1)^M N_d$; equivalently, sufficient bag-of-words counts are $c_(d, v) = sum_(n=1)^(N_d) bb(1)[w_(d, n) = v]$ with $N_"nnz" = sum_(d=1)^M |\{v : c_(d, v) > 0\}|$ nonzero entries. Given topic count $K$, Dirichlet hyperparameter $bold(alpha) in bb(R)_(> 0)^K$, and topic-word distributions $bold(B) = (beta_(k, v)) in bb(R)_(> 0)^(K times V)$ with rows $bold(beta)_k in Delta^(V-1)$ (optionally with symmetric prior $bold(beta)_k tilde.op "Dir"(bold(eta))$), LDA generates each document $d in \{1, dots.h, M\}$ as:
1. Draw document topic mixture $bold(theta)_d tilde.op "Dir"(bold(alpha))$ on $Delta^(K-1)$.
2. For each token $n in \{1, dots.h, N_d\}$, draw latent topic $z_(d, n) tilde.op "Cat"(bold(theta)_d)$ and observed word $w_(d, n) tilde.op "Cat"(bold(beta)_(z_(d, n)))$.

Marginalizing $(bold(theta)_d, bold(z)_d)$ gives the exact document marginal likelihood:
$
  p(bold(w)_d | bold(alpha), bold(B)) &= integral_(Delta^(K-1)) p(bold(theta)_d | bold(alpha)) \
  &quad dot (product_(n=1)^(N_d) sum_(k=1)^K theta_(d, k) beta_(k, w_(d, n))) dif bold(theta)_d ,
$
where the Dirichlet prior density is
$
  p(bold(theta) | bold(alpha)) = (Gamma(sum_(k=1)^K alpha_k)) / (product_(k=1)^K Gamma(alpha_k)) product_(k=1)^K theta_k^(alpha_k - 1) .
$
Because this integral couples $bold(theta)_d$ and $bold(B)$ inside a product of sums, empirical Bayes estimation maximizes a tractable variational Evidence Lower Bound (ELBO) @blei2003lda @bishop2006pattern.

#definition(title: "Document-Level Mean-Field Variational ELBO for LDA")[
  For each document $d$, introduce free variational Dirichlet parameters $bold(gamma)_d in bb(R)_(> 0)^K$ and token-level multinomial distributions $bold(phi)_(d, n) in Delta^(K-1)$ via the factorized surrogate
  $
    q(bold(theta)_d, bold(z)_d | bold(gamma)_d, bold(phi)_d) = q(bold(theta)_d | bold(gamma)_d) product_(n=1)^(N_d) q(z_(d, n) | bold(phi)_(d, n)) .
  $
  The document ELBO $cal(L)_d (bold(gamma)_d, bold(phi)_d; bold(alpha), bold(B)) := bb(E)_q [log p(bold(theta)_d, bold(z)_d, bold(w)_d | bold(alpha), bold(B))] - bb(E)_q [log q(bold(theta)_d, bold(z)_d)]$ expands explicitly as
  $
    cal(L)_d &= log Gamma(sum_(k=1)^K alpha_k) - sum_(k=1)^K log Gamma(alpha_k) \
    &quad + sum_(k=1)^K (alpha_k - gamma_(d, k)) (psi(gamma_(d, k)) - psi(Gamma_d)) \
    &quad - log Gamma(Gamma_d) + sum_(k=1)^K log Gamma(gamma_(d, k)) \
    &quad + sum_(n=1)^(N_d) sum_(k=1)^K phi_(d, n, k) [psi(gamma_(d, k)) - psi(Gamma_d)] \
    &quad + sum_(n=1)^(N_d) sum_(k=1)^K phi_(d, n, k) log (beta_(k, w_(d, n)) / phi_(d, n, k)) ,
  $
  where $Gamma_d = sum_(j=1)^K gamma_(d, j)$ and $psi(x) = (dif / dif x) log Gamma(x)$ is the digamma function.
]

#algbox(
  title: "Variational EM Algorithm for Latent Dirichlet Allocation",
  inputs: [Corpus counts $\{c_(d, v)\}$, topics $K$, smoothing $eta >= 0$, tolerances $tau_"in", tau_"out" > 0$.],
  outputs: [Topic-word matrix $hat(bold(B)) in bb(R)^(K times V)$, prior $hat(bold(alpha))$, document mixtures $hat(bold(gamma))_d$.],
)[
  1. *Initialize* $bold(alpha)^((0)) in bb(R)_(> 0)^K$ and $bold(B)^((0))$ with strictly positive rows $bold(beta)_k^((0)) in Delta^(V-1)$.
  2. *Repeat* (outer Variational EM iteration $t = 0, 1, dots.h, T_"out" - 1$):
  3. #h(1em) *Variational E-step:* For each document $d = 1, dots.h, M$, initialize $gamma_(d, k) <- alpha_k + N_d / K$ and iterate up to $T_"in"$ times:
     - For each unique active word $v$ with $c_(d, v) > 0$ and topic $k$:
       $
         phi_(d, v, k) prop beta_(k, v) exp(psi(gamma_(d, k)) - psi(sum_(j=1)^K gamma_(d, j))) .
       $
     - Update $gamma_(d, k) <- alpha_k + sum_(v : c_(d, v) > 0) c_(d, v) phi_(d, v, k)$ until $||bold(gamma)_d^("new") - bold(gamma)_d||_1 < tau_"in"$.
  4. #h(1em) *Variational M-step:* Update topic-word distributions (with optional Dirichlet smoothing $eta >= 0$):
     $
       beta_(k, v) <- (eta + sum_(d=1)^M c_(d, v) phi_(d, v, k)) / (V eta + sum_(v'=1)^V sum_(d=1)^M c_(d, v') phi_(d, v', k)) .
     $
     Update $bold(alpha)$ by linear-time Newton-Raphson (Proposition 2) on $sum_(d=1)^M cal(L)_d$.
  5. *Until* relative corpus ELBO improvement falls below $tau_"out"$.
  6. *Return* $hat(bold(B))$, $hat(bold(alpha))$, and normalized mixtures $hat(theta)_(d, k) = hat(gamma)_(d, k) / sum_j hat(gamma)_(d, j)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove the combinatorial intractability of exact posterior inference, establish a self-contained strict positive-definiteness theorem for the Dirichlet Fisher matrix that proves unique variational stationarity, derive the exact $cal(O)(K)$ Newton step for $bold(alpha)$, formalize de Finetti exchangeability and simplex geometry, and prove LDA's parameter-scaling generalization advantage over pLSA.

#subsec("A", "Combinatorial Intractability of Exact Inference")

#theorem(title: "Exponential Monomial Explosion of the Exact LDA Marginal", origin: "Original Proof")[
  For any document $bold(w)_d$ of length $N_d$ with strictly positive $bold(B) > bold(0)$, the exact marginal likelihood admits the finite expansion
  $
    p(bold(w)_d | bold(alpha), bold(B)) &= (Gamma(A_0)) / (Gamma(A_0 + N_d)) sum_(bold(m) in cal(M)(N_d, K)) (product_(k=1)^K (Gamma(alpha_k + m_k)) / (Gamma(alpha_k))) \
    &quad dot C_(bold(w)_d) (bold(m); bold(B)) ,
  $
  where $A_0 = sum_k alpha_k$, $cal(M)(N_d, K) = \{bold(m) in bb(Z)_(>= 0)^K : sum_k m_k = N_d\}$, and $C_(bold(w)_d) (bold(m); bold(B)) = sum_(bold(z) : n(bold(z)) = bold(m)) product_(n=1)^(N_d) beta_(z_n, w_(d, n))$. The number of distinct coupled Dirichlet modes equals $binom(N_d + K - 1, K - 1)$ and the direct sum contains $K^(N_d)$ assignment paths, rendering exact marginalization exponential in $min(K, N_d)$.
]

#proof[
  Expand the product inside the Dirichlet integral using the distributive law over all $K^(N_d)$ topic assignment sequences $bold(z) = (z_1, dots.h, z_(N_d)) in \{1, dots.h, K\}^(N_d)$:
  $
    product_(n=1)^(N_d) (sum_(k=1)^K theta_k beta_(k, w_(d, n))) &= sum_(bold(z) in \{1, dots.h, K\}^(N_d)) (product_(k=1)^K theta_k^(m_k (bold(z)))) \
    &quad dot (product_(n=1)^(N_d) beta_(z_n, w_(d, n))) ,
  $
  where $m_k (bold(z)) = sum_(n=1)^(N_d) bb(1)[z_n = k]$ counts tokens assigned to topic $k$. Grouping sequences with identical count vector $bold(m) = (m_1, dots.h, m_K) in cal(M)(N_d, K)$ produces the multinomial generating coefficient $C_(bold(w)_d) (bold(m); bold(B))$. Integrating each monomial $product_k theta_k^(alpha_k + m_k - 1)$ over the simplex $Delta^(K-1)$ via the exact Dirichlet normalization identity
  $
    integral_(Delta^(K-1)) product_(k=1)^K theta_k^(alpha_k + m_k - 1) dif bold(theta) &= (product_(k=1)^K Gamma(alpha_k + m_k)) / (Gamma(A_0 + sum_(k=1)^K m_k))
  $
  immediately yields the stated formula. By stars-and-bars combinatorics, $|cal(M)(N_d, K)| = binom(N_d + K - 1, K - 1)$, which grows as $(N_d + K - 1)^(K-1) / (K - 1)!$ for fixed $K$ and exponentially as $Omega(2^(min(K, N_d)))$ when both $K$ and $N_d$ scale. Moreover, even for a single count vector $bold(m)$, evaluating $C_(bold(w)_d) (bold(m); bold(B))$ requires summing products of entries from distinct rows of $bold(B)$ subject to row-sum multiplicities $bold(m)$—a permanent-like multilinear sum that does not factor across $n$—demonstrating why exact posterior inference $p(bold(theta)_d | bold(w)_d, bold(alpha), bold(B))$ is a mixture of exponentially many coupled Dirichlet densities.
]

#subsec("B", "Strict Positivity of the Dirichlet Fisher Matrix and Unique Variational Updates")

#theorem(title: "Positive Definiteness of the Dirichlet Fisher Matrix and Unique E-Step Fixed Point", origin: "Original Proof")[
  For any $K >= 2$ and any $bold(gamma) in bb(R)_(> 0)^K$ with sum $S = sum_(k=1)^K gamma_k$, define the trigamma matrix
  $
    bold(H)(bold(gamma)) := "diag"(psi'(gamma_1), dots.h, psi'(gamma_K)) - psi'(S) bold(1) bold(1)^top ,
  $
  where $psi'(x) = (dif / dif x) psi(x) = sum_(m=0)^infinity (x + m)^(-2)$. Then:
  1. $bold(H)(bold(gamma)) in bb(S)_(succ 0)^K$ (strictly positive definite everywhere on $bb(R)_(> 0)^K$);
  2. Holding $bold(phi)_d$ fixed, the stationarity condition $nabla_(bold(gamma)_d) cal(L)_d = bold(0)$ has the unique closed-form solution $gamma_(d, k)^* = alpha_k + sum_(n=1)^(N_d) phi_(d, n, k)$, which strictly maximizes $cal(L)_d$ over $bold(gamma)_d in bb(R)_(> 0)^K$.
]

#proof[
  *Part 1.* We first prove that $h(x) := x psi'(x)$ is strictly decreasing on $(0, infinity)$. Using the integral representation $psi'(x) = integral_0^infinity (t e^(-x t)) / (1 - e^(-t)) dif t$ and substituting $u = x t$ ($dif t = dif u / x$), we obtain
  $
    h(x) = x psi'(x) = integral_0^infinity (u e^(-u)) / (x (1 - e^(-u / x))) dif u .
  $
  For any fixed $u > 0$, set $y = u / x > 0$; then $x (1 - e^(-u / x)) = u (1 - e^(-y)) / y = u integral_0^1 e^(-s y) dif s$. As $x$ increases, $y$ strictly decreases, so $e^(-s y)$ strictly increases for every $s in (0, 1)$, making the denominator $x (1 - e^(-u / x))$ strictly increasing in $x$. Consequently, the integrand is strictly decreasing in $x$ for all $u > 0$, proving that $h(x)$ is strictly decreasing on $(0, infinity)$.

  Therefore $g(x) := 1 / psi'(x)$ satisfies that $g(x) / x = 1 / h(x)$ is strictly increasing on $(0, infinity)$. Because $K >= 2$ and every $gamma_k > 0$, we have $gamma_k < S = sum_(j=1)^K gamma_j$, whence $g(gamma_k) / gamma_k < g(S) / S$ for every $k$. Multiplying by $gamma_k$ and summing over $k = 1, dots.h, K$ yields strict superadditivity:
  $
    sum_(k=1)^K 1 / psi'(gamma_k) &= sum_(k=1)^K g(gamma_k) \
    &< (g(S) / S) sum_(k=1)^K gamma_k = 1 / psi'(S) ,
  $
  so $psi'(S) sum_(k=1)^K psi'(gamma_k)^(-1) < 1$. Now let $bold(v) in bb(R)^K backslash \{bold(0)\}$ and $D_(k k) = psi'(gamma_k) > 0$. By Cauchy-Schwarz,
  $
    (bold(1)^top bold(v))^2 &= (sum_(k=1)^K D_(k k)^(-1 / 2) (D_(k k)^(1 / 2) v_k))^2 \
    &<= (sum_(k=1)^K 1 / D_(k k)) (sum_(k=1)^K D_(k k) v_k^2) < (1 / psi'(S)) bold(v)^top bold(D) bold(v) .
  $
  Rearranging yields $bold(v)^top bold(H)(bold(gamma)) bold(v) = bold(v)^top bold(D) bold(v) - psi'(S) (bold(1)^top bold(v))^2 > 0$, proving $bold(H)(bold(gamma)) in bb(S)_(succ 0)^K$.

  *Part 2.* Differentiating $cal(L)_d$ from Definition 1 with respect to $gamma_(d, k)$ (noting that $(partial / partial gamma_(d, k)) (log Gamma(gamma_(d, k)) - log Gamma(Gamma_d)) = psi(gamma_(d, k)) - psi(Gamma_d)$ cancels the derivative of the linear coefficient) gives the exact vector gradient
  $
    nabla_(bold(gamma)_d) cal(L)_d = bold(H)(bold(gamma)_d) (bold(alpha) + sum_(n=1)^(N_d) bold(phi)_(d, n) - bold(gamma)_d) .
  $
  Because $bold(H)(bold(gamma)_d)$ is strictly positive definite everywhere on $bb(R)_(> 0)^K$, it is nonsingular, so $nabla_(bold(gamma)_d) cal(L)_d = bold(0)$ holds if and only if $bold(gamma)_d = bold(alpha) + sum_(n=1)^(N_d) bold(phi)_(d, n) := bold(gamma)_d^*$. Evaluating the Hessian at $bold(gamma)_d^*$ cancels all third-order tetragamma terms (which multiply the vanishing factor $bold(gamma)_d^* - bold(gamma)_d^* = bold(0)$), leaving $nabla^2_(bold(gamma)_d) cal(L)_d (bold(gamma)_d^*) = -bold(H)(bold(gamma)_d^*) prec bold(0)$. Because $bold(gamma)_d^*$ is the unique critical point in $bb(R)_(> 0)^K$ and $cal(L)_d -> -infinity$ on the boundary of $bb(R)_(> 0)^K$, it is the unique global maximizer for fixed $bold(phi)_d$. Optimizing $cal(L)_d$ over $bold(phi)_(d, n) in Delta^(K-1)$ is strictly concave by negative entropy, uniquely yielding $phi_(d, n, k) prop beta_(k, w_(d, n)) exp(psi(gamma_(d, k)) - psi(Gamma_d))$.
]

#subsec("C", "Linear-Time Exact Newton Step for Hyperparameters")

#proposition(title: "Exact O(K) Newton-Raphson Step for Dirichlet Concentration", origin: "Original Derivation")[
  Let $cal(L)(bold(alpha)) = sum_(d=1)^M cal(L)_d$ denote the corpus ELBO as a function of $bold(alpha) in bb(R)_(> 0)^K$, with gradient components
  $
    g_k &= M [psi(sum_(j=1)^K alpha_j) - psi(alpha_k)] \
    &quad + sum_(d=1)^M [psi(gamma_(d, k)) - psi(Gamma_d)] .
  $
  Then $cal(L)(bold(alpha))$ is strictly concave on $bb(R)_(> 0)^K$, and the exact Newton ascent step $bold(Delta) = -(nabla^2 cal(L)(bold(alpha)))^(-1) bold(g)$ is computed in $cal(O)(K)$ operations without matrix inversion via
  $
    Delta_k = (g_k + chi) / D_(k k) , quad chi = (sum_(j=1)^K g_j / D_(j j)) / (c_0^(-1) - sum_(j=1)^K D_(j j)^(-1)) ,
  $
  where $D_(k k) = M psi'(alpha_k) > 0$, $c_0 = M psi'(sum_j alpha_j) > 0$, and the denominator is strictly positive everywhere.
]

#proof[
  Differentiating $g_k$ gives the Hessian $nabla^2 cal(L)(bold(alpha)) = -bold(D) + c_0 bold(1) bold(1)^top = -M bold(H)(bold(alpha))$, where $bold(H)(bold(alpha))$ is the Dirichlet Fisher matrix of Theorem 2. By Theorem 2, $bold(H)(bold(alpha)) in bb(S)_(succ 0)^K$, so $nabla^2 cal(L)(bold(alpha)) prec bold(0)$ everywhere on $bb(R)_(> 0)^K$, establishing strict concavity and certifying $c_0^(-1) - sum_j D_(j j)^(-1) = M^(-1) [psi'(sum_j alpha_j)^(-1) - sum_j psi'(alpha_j)^(-1)] > 0$. Applying the Sherman-Morrison inversion identity to $(bold(D) - c_0 bold(1) bold(1)^top)^(-1)$ yields
  $
    bold(Delta) &= (bold(D) - c_0 bold(1) bold(1)^top)^(-1) bold(g) \
    &= bold(D)^(-1) bold(g) + (bold(1)^top bold(D)^(-1) bold(g)) / (c_0^(-1) - bold(1)^top bold(D)^(-1) bold(1)) bold(D)^(-1) bold(1) .
  $
  Reading off the $k$-th coordinate gives $Delta_k = (g_k + chi) / D_(k k)$, which requires only two scalar sums over $j = 1, dots.h, K$, achieving exact $cal(O)(K)$ complexity instead of $cal(O)(K^3)$ @blei2003lda.
]

#subsec("D", "De Finetti Exchangeability and Topic-Simplex Geometry")

#proposition(title: "Infinite Token Exchangeability and Convex-Hull Geometry", origin: "Classical — Cited")[
  Under LDA @blei2003lda, every document's token sequence $(w_(d, 1), dots.h, w_(d, N_d))$ is infinitely exchangeable (invariant under all finite permutations $sigma in S_(N_d)$), and its conditional unigram distribution $bold(p)_d = bold(B)^top bold(theta)_d$ lies in the $(K - 1)$-dimensional topic simplex $"conv"(\{bold(beta)_1, dots.h, bold(beta)_K\}) subset.eq Delta^(V-1)$. Moreover, the model is invariant under simultaneous permutation of topic indices $k in \{1, dots.h, K\}$.
]

#proof[
  For any permutation $sigma in S_(N_d)$, the product inside the marginal integral satisfies $product_(n=1)^(N_d) sum_k theta_(d, k) beta_(k, w_(d, sigma(n))) = product_(n=1)^(N_d) sum_k theta_(d, k) beta_(k, w_(d, n))$ because multiplication in $bb(R)$ is commutative; hence $p(bold(w)_d | bold(alpha), bold(B))$ depends solely on the word-count histogram $\{c_(d, v)\}$, matching de Finetti's representation of infinitely exchangeable sequences as mixtures of i.i.d. categorical draws with mixing measure $bold(theta)_d tilde.op "Dir"(bold(alpha))$ @blei2003lda. Conditioned on $bold(theta)_d in Delta^(K-1)$, the marginal probability of drawing word $v$ is $"Pr"(w_(d, n) = v | bold(theta)_d, bold(B)) = sum_(k=1)^K theta_(d, k) beta_(k, v) = (bold(B)^top bold(theta)_d)_v$, which sweeps out the convex hull of the $K$ topic vectors $\{bold(beta)_k\}$ inside $Delta^(V-1)$. Permuting columns of $bold(theta)_d$, coordinates of $bold(alpha)$, and rows of $bold(B)$ leaves the sum over $k$ identical.
]

#subsec("E", "Parameter-Scaling Generalization Advantage")

#proposition(title: "Corpus-Independent Parameter Complexity and Out-of-Sample Support", origin: "Original Proof")[
  The global parameter count of empirical Bayes LDA is $p_"LDA" = K + K (V - 1) = K V$, independent of corpus size $M$, in contrast to $p_"pLSA" = M K + K V - K - 1 = Omega(M K)$ for pLSA (Alg. 84). Moreover, whenever $bold(B) > bold(0)$ (guaranteed by any smoothing $eta > 0$), LDA defines a strictly positive normalized probability distribution over every finite unseen document $bold(w)_"new" in cal(V)^N$, satisfying $sum_(bold(w) in cal(V)^N) p(bold(w) | bold(alpha), bold(B)) = 1$ and finite held-out perplexity.
]

#proof[
  In LDA, the document proportions $bold(theta)_d$ are integrated out as random variables drawn from $"Dir"(bold(alpha))$ rather than stored as free parameters; only $bold(alpha) in bb(R)_(> 0)^K$ ($K$ parameters) and $K$ simplex rows $bold(beta)_k in Delta^(V-1)$ ($K (V - 1)$ parameters) are estimated globally, giving $p_"LDA" = K V$ regardless of $M$. For any fixed length $N$ and unseen document $bold(w)_"new" in cal(V)^N$, set $beta_min := min_(k, v) beta_(k, v) > 0$. Then $sum_(k=1)^K theta_k beta_(k, w_n) >= beta_min sum_k theta_k = beta_min$ for all $bold(theta) in Delta^(K-1)$, so $p(bold(w)_"new" | bold(alpha), bold(B)) >= beta_min^N > 0$, ensuring finite perplexity $exp(-N^(-1) log p(bold(w)_"new")) <= beta_min^(-1) < infinity$. Summing over $bold(w) in cal(V)^N$ and exchanging finite sum and Dirichlet integral gives
  $
    &sum_(bold(w) in cal(V)^N) p(bold(w) | bold(alpha), bold(B)) \
    &= integral_(Delta^(K-1)) p(bold(theta) | bold(alpha)) (sum_(v=1)^V sum_(k=1)^K theta_k beta_(k, v))^N dif bold(theta) = 1 ,
  $
  since $sum_v beta_(k, v) = 1$ and $sum_k theta_k = 1$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational Complexity and Sparse Aggregation.* Because identical word tokens within document $d$ share identical variational distributions $bold(phi)_(d, v)$, the variational E-step iterates only over the $N_"nnz"$ unique nonzero document-term pairs rather than all $N_"tot"$ raw tokens. Each outer iteration runs up to $T_"in"$ inner document updates costing $cal(O)(K N_"nnz")$ arithmetic plus digamma/exp evaluations, followed by an $cal(O)(K N_"nnz" + K V)$ M-step and $cal(O)(K)$ Newton hyperparameter step (Proposition 2), yielding total complexity $cal(O)(T_"out" (T_"in" K N_"nnz" + K V))$ and memory $cal(O)(K V + M K + K max_d |"supp"(bold(w)_d)|)$. While asymptotically linear in $N_"nnz"$ and $K$, the nested inner loop and transcendental $psi(x)$ and $exp(x)$ calls give LDA a substantially higher constant factor than single-loop pLSA (Alg. 84) or closed-form GMM (Alg. 81). Collapsed Gibbs Sampling (CGS) avoids variational approximation bias by sampling $z_(d, n)$ directly in $cal(O)(T K N_"tot")$ time, at the cost of Monte Carlo variance and harder convergence diagnostics.

*Simplex Geometry and Sparsity Control.* The Dirichlet concentration $bold(alpha)$ acts as a geometric sparsity dial on $Delta^(K-1)$: when $alpha_k < 1$, prior density diverges at the vertices and edges of $Delta^(K-1)$, concentrating each document onto a few active topics (crisp soft clustering); when $alpha_k > 1$, prior mass concentrates near the simplex barycenter, producing diffuse mixtures. Likewise, $eta < 1$ encourages sharp, non-overlapping topic vocabularies $bold(beta)_k$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Objective], [Maximize variational ELBO $sum_d cal(L)_d (bold(gamma)_d, bold(phi)_d; bold(alpha), bold(B))$],
    [Global Parameters], [$K V$ (constant in corpus size $M$, Proposition 4)],
    [Time Complexity], [$cal(O)(T_"out" (T_"in" K N_"nnz" + K V))$ (variational EM)],
    [Space Complexity], [$cal(O)(K V + M K)$ global + $cal(O)(K |"supp"(bold(w)_d)|)$ per doc],
    [Cluster Geometry], [Mixed-membership simplex $"conv"(\{bold(beta)_k\}) subset.eq Delta^(V-1)$],
    [Invariance], [Token-order exchangeable; topic-permutation invariant],
    [Input Domain], [Non-negative count / frequency vectors (bag-of-words)],
    [Key Controls], [Topic count $K$, concentration $bold(alpha)$, word prior $eta$],
  ),
  caption: [Operational and theoretical profile of Latent Dirichlet Allocation (Alg. 83).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Scientific Literature Indexing and Patent Landscaping:* Research articles and patents naturally span multiple intersecting methodologies (e.g., optimization + genomics + hardware). LDA's mixed-membership simplex representation models multi-topic articles natively while corpus-independent parameterization (Proposition 4) enables principled perplexity scoring and fold-in of newly published papers.
2. *Population Genetics Admixture Inference:* Multilocus allele-count genotypes of admixed individuals are mathematically isomorphic to documents composed of allele tokens drawn from ancestral population proportions $bold(theta)_d tilde.op "Dir"(bold(alpha))$. LDA's hierarchical Dirichlet-multinomial model directly parallels the classical genetic admixture model, recovering latent ancestral allele frequencies $bold(beta)_k$.
3. *Single-Cell Epigenomics (scATAC-seq) and Transcriptomic Modules:* Sparse chromatin-accessibility peak counts across thousands of cells exhibit extreme zero inflation and overlapping regulatory gene programs. Representing each cell as a mixture of cis-regulatory topics $bold(theta)_d$ denoises dropout noise and maps co-accessible enhancers to cell-state pathways.
4. *Electronic Health Record (EHR) Clinical Phenotyping:* Patient visit histories encoded as sparse counts of ICD diagnosis codes, medications, and laboratory orders reflect concurrent chronic comorbidities. LDA decomposes multimorbid patients into interpretable disease-burden weights $bold(theta)_d$ without forcing mutually exclusive clinical labels.
5. *E-Commerce Market-Basket and Digital Humanities Analysis:* Customer purchase histories and historical archives consist of high-dimensional sparse item or lexical histograms ($V approx 10^4$--$10^5$). Additive linear-in-$V$ parameterization avoids quadratic covariance explosion while yielding human-auditable top-word topic summaries.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

LDA exhibits five documented failure modes: (i) *short-text sparsity*—on microblogs or titles ($N_d < 15$), per-document word co-occurrence is too scarce to separate $bold(theta)_d$ from the prior $bold(alpha)$, degrading coherence (use Biterm Topic Model @yan2013btm or Dirichlet Multinomial Mixture @yin2014dmm); (ii) *inter-topic correlation*—the Dirichlet prior enforces weak negative covariance $"Cov"(theta_j, theta_k) < 0$ and cannot capture positively co-occurring topic pairs (use Correlated Topic Models @blei2007ctm); (iii) *continuous real-valued / signed geometry*—applying LDA to dense Euclidean features requires ad-hoc non-negative shifting/discretization and discards radial scale (use GMM Alg. 81 or DPGMM Alg. 82); (iv) *fixed topic number $K$*—mis-specifying $K$ fragments or merges themes (use Hierarchical Dirichlet Process HDP-LDA @teh2006hdp); and (v) *syntax blindness*—bag-of-words exchangeability ignores word order and negation.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Long sparse text / count corpora ($N_d >= 50$)], [Exact hierarchical Dirichlet-multinomial generative fit], [Tune $K, bold(alpha), eta$; use HDP @teh2006hdp if $K$ unknown],
    [Short texts / tweets ($N_d <= 15$)], [Severe per-document sparsity; noisy $bold(gamma)_d$], [BTM @yan2013btm or DMM @yin2014dmm],
    [Strongly correlated themes], [Dirichlet forbids positive topic correlation], [Correlated Topic Model (CTM) @blei2007ctm],
    [Dense continuous Euclidean data], [Poor; discards radial norm and covariance tilt], [GMM (Alg. 81) or DPGMM (Alg. 82)],
    [Latency-critical large count matrices], [Inner loop + digamma adds $10 times$ overhead], [pLSA (Alg. 84) or KL-NMF @lee2001algorithms],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for LDA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Uninformative-Prior Reduction of LDA to pLSA", origin: "Original Derivation")[
  Consider the MAP estimation problem for $(bold(theta)_1, dots.h, bold(theta)_M, bold(B))$ under the LDA joint density $p(bold(W), bold(Theta) | bold(alpha), bold(B)) = product_(d=1)^M p(bold(theta)_d | bold(alpha)) product_(v=1)^V (sum_k theta_(d, k) beta_(k, v))^(c_(d, v))$. Setting the uniform Dirichlet hyperparameter $bold(alpha) = bold(1)$ (flat prior on $Delta^(K-1)$) and replacing posterior marginalization over $bold(theta)_d$ with point maximization recovers the exact pLSA log-likelihood objective (Alg. 84) @hofmann1999plsa @blei2003lda.
]

#proof[
  When $alpha_k = 1$ for all $k in \{1, dots.h, K\}$, the Dirichlet prior reduces to the constant density $p(bold(theta)_d | bold(1)) = (K - 1)!$ everywhere on $Delta^(K-1)$. Taking the logarithm of the joint density with respect to $(bold(Theta), bold(B))$ yields
  $
    log p(bold(W), bold(Theta) | bold(1), bold(B)) &= M log((K - 1)!) \
    &quad + sum_(d=1)^M sum_(v=1)^V c_(d, v) log (sum_(k=1)^K theta_(d, k) beta_(k, v)) .
  $
  Dropping the constant $M log((K - 1)!)$ leaves identically the asymmetric pLSA log-likelihood $ell_"pLSA" (bold(Theta), bold(B))$ up to the fixed empirical document-length constant $sum_d N_d log(N_d / N_"tot")$. Thus LDA strictly generalizes pLSA by placing a learnable non-uniform prior $"Dir"(bold(alpha))$ on $bold(theta)_d$ and integrating over $bold(theta)_d$ rather than point-estimating $M (K - 1)$ document parameters.
]

Within Category IV, LDA sits at the center of probabilistic topic modeling: removing the Dirichlet prior and marginalization recovers pLSA (Alg. 84 @hofmann1999plsa), which in turn matches KL-NMF @ding2008nmfplsa; replacing the Dirichlet with a logistic-normal prior yields Correlated Topic Models @blei2007ctm; placing a Dirichlet Process over the topic inventory yields HDP-LDA @teh2006hdp; and restricting each document to a single topic ($bold(theta)_d$ one-hot) yields Dirichlet Multinomial Mixture (DMM) @yin2014dmm.

*Empirical Benchmark Behavior (ClusterOverview 4-Dataset Suite).* On the four standardized datasets (evaluated after non-negative feature adaptation for count-based topic models), LDA achieves: *Iris* Silhouette $= 0.182$, ARI $= 0.429$; *Customers* Silhouette $= 0.245$; *Multi-Shape* Silhouette $= 0.706$ (ARI N/A); *High-D* Silhouette $= 0.232$ (ARI N/A); and mean runtime $24.476$ ms. These results illuminate LDA's inductive bias with surgical precision. On low-dimensional dense Euclidean benchmarks (*Iris* and *Customers*), LDA underperforms Gaussian mixtures (GMM ARI $= 0.950$, Customers Sil $= 0.332$) because its bag-of-counts simplex likelihood (Proposition 3) responds only to relative feature proportions and is further blurred by Dirichlet prior smoothing—even lagging unsmoothed pLSA (Iris ARI $= 0.507$, Customers Sil $= 0.311$). By sharp contrast, on *High-D*, LDA achieves Silhouette $= 0.232$, *outperforming both GMM ($0.181$) and DPGMM ($0.208$)* because its linear $cal(O)(K V)$ simplex parameterization avoids the quadratic $cal(O)(K d^2)$ covariance overparameterization of Gaussian models! Finally, its $24.476$ ms mean runtime is $10.7 times$ slower than pLSA ($2.291$ ms) and $17.5 times$ slower than GMM ($1.400$ ms), directly confirming the nested inner-loop and digamma cost quantified in Section III.
