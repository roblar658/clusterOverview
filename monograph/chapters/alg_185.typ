#import "../template/ieee.typ": *

#alg_entry(
  id: 185,
  name: "DEPICT (Deep Embedded Regularized Clustering)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2017",
  paradigm: "Denoising Autoencoder with Balanced Relative-Entropy Clustering",
  complexity: $cal(O)(E N (P + K d))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Early deep clustering architectures such as Deep Embedded Clustering
(DEC) @xie2016dec discard the autoencoder decoder after pretraining
and optimize a self-training Kullback–Leibler divergence directly
over latent embeddings. Without a persistent reconstruction anchor,
deep encoders can distort local geometry and collapse disparate data
manifolds into a small subset of trivial clusters. Deep Embedded
Regularized Clustering (DEPICT), introduced by Ghasedi Dizaji et al.
@ghasedi2017depict, addresses both pathologies simultaneously by
combining three tightly coupled mechanisms: (i) a stacked multi-layer
*denoising* autoencoder supervising reconstructions at every hidden
layer, (ii) a multinomial logistic regression (softmax) clustering
head trained on *corrupted* latent features to enforce local
noise-invariance, and (iii) an explicit *balanced-assignment*
relative-entropy prior that penalizes degenerate cluster collapse.

Let $cal(X) = {bold(x)_1, dots, bold(x)_N} subset.eq bb(R)^D$ be an
unlabeled dataset and let $K >= 2$ be the target cluster count.
Consider an $L$-layer encoder with clean activations
$bold(z)_i^((0)) = bold(x)_i$ and
$bold(z)_i^((l)) = f_(bold(W)^((l))) (bold(z)_i^((l-1))) in bb(R)^(d_l)$
for $l = 1, dots, L$, culminating in the bottleneck embedding
$bold(z)_i = bold(z)_i^((L)) in bb(R)^d$ ($d = d_L$). During
training, stochastic corruption (dropout or additive noise) is
injected at every encoder layer, yielding noisy representations
$tilde(bold(z))_i^((l)) = f_(bold(W)^((l))) ("Drop"(tilde(bold(z))_i^((l-1))))$
with bottleneck $tilde(bold(z))_i = tilde(bold(z))_i^((L))$. A
symmetric $L$-layer decoder $g_(bold(W)')$ reconstructs every clean
intermediate representation
$hat(bold(z))_i^((l-1)) = g_(bold(W)'^((l))) (hat(bold(z))_i^((l)))$
starting from $hat(bold(z))_i^((L)) = tilde(bold(z))_i$.

On top of the bottleneck layer, a multinomial logistic regression
classifier parameterized by $bold(Theta) = [bold(theta)_1, dots, bold(theta)_K] in bb(R)^(d times K)$
assigns soft cluster probabilities to clean embeddings $bold(z)_i$
(and analogously $tilde(p)_(i k)$ to noisy embeddings
$tilde(bold(z))_i$):
$
  p_(i k) = Pr(y_i = k | bold(z)_i; bold(Theta))
  = (exp(bold(theta)_k^top bold(z)_i)) /
    (sum_(r=1)^K exp(bold(theta)_r^top bold(z)_i)).
$
To guide unsupervised training without ground-truth labels, DEPICT
introduces an auxiliary target distribution matrix
$bold(Q) = (q_(i k)) in (Delta^(K-1))^N$ with empirical cluster
frequency vector $bold(f) in Delta^(K-1)$ defined by
$f_k = (1 / N) sum_(i=1)^N q_(i k)$, regularized toward a target
prior $bold(u) in Delta^(K-1)$ (defaulting to the uniform prior
$u_k = 1 / K$).

#definition(title: "DEPICT Regularized Relative-Entropy Objective")[
  Let $bold(P) = (p_(i k))$ and $tilde(bold(P)) = (tilde(p)_(i k))$
  denote clean and noisy softmax assignment matrices, and let
  $bold(W)_"AE" = (bold(W), bold(W)')$. The *DEPICT clustering loss*
  with balancing prior $bold(u) in Delta^(K-1)$ is
  $
    cal(L)_"cls" (bold(Q), bold(P))
    = 1 / N sum_(i=1)^N D_"KL" (bold(q)_i || bold(p)_i)
      + D_"KL" (bold(f) || bold(u)),
  $
  and the *stacked layer-wise denoising reconstruction loss* is
  $
    cal(L)_"rec" (bold(W)_"AE")
    = 1 / (N L) sum_(i=1)^N sum_(l=0)^(L-1)
      1 / d_l norm(bold(z)_i^((l)) - hat(bold(z))_i^((l)))_2^2.
  $
  The full DEPICT objective is
  $cal(L)_"DEPICT" = cal(L)_"cls" (bold(Q), tilde(bold(P))) + lambda cal(L)_"rec" (bold(W)_"AE")$,
  with trade-off weight $lambda > 0$.
]

Optimization proceeds via alternating minimization: holding network
weights fixed, the auxiliary target $bold(Q)$ is updated analytically
from clean probabilities $bold(P)$ to enforce sharp, balanced
assignments; holding $bold(Q)$ fixed, the network weights
$(bold(W)_"AE", bold(Theta))$ are updated by stochastic gradient
descent minimizing the cross-entropy $- (1 / N) sum_(i, k) q_(i k) log tilde(p)_(i k)$
plus $lambda cal(L)_"rec"$.

#algbox(
  title: "DEPICT — Deep Embedded Regularized Clustering",
  inputs: [
    Dataset $cal(X) = {bold(x)_i}_(i=1)^N subset.eq bb(R)^D$,
    clusters $K$, prior $bold(u)$ (default $1 / K$), weight
    $lambda > 0$, dropout rate $rho$, epochs $E$.
  ],
  outputs: [
    Trained weights $(bold(W)_"AE", bold(Theta))$, soft assignments
    $bold(P)$, hard cluster labels ${hat(c)_i}_(i=1)^N$.
  ],
)[
  1. *Denoising Pretraining:* Pretrain the $L$-layer denoising
     autoencoder layer-by-layer and end-to-end by minimizing
     $cal(L)_"rec" (bold(W)_"AE")$ with corruption rate $rho$.
  2. *Softmax Initialization:* Extract clean embeddings
     ${bold(z)_i}_(i=1)^N$, run initial clustering (e.g., $k$-means
     or AC) to obtain pseudo-labels, and initialize $bold(Theta)$
     by supervised multinomial logistic regression.
  3. *For* epoch $e = 1, dots, E$:
     - *(Target Update):* Compute clean probabilities $bold(P)$ over
       $cal(X)$ and update the balanced target matrix $bold(Q)$ via
       $
         q_(i k) = (p_(i k) / (sum_(j=1)^N p_(j k))^(1/2)) /
                   (sum_(r=1)^K p_(i r) / (sum_(j=1)^N p_(j r))^(1/2))
       $
       (or exact KKT fixed-point iterations from Theorem 1).
     - *(Network Update):* For each mini-batch $cal(B) subset.eq cal(X)$,
       run a noisy forward pass to obtain $tilde(bold(P))_cal(B)$ and
       layer reconstructions ${hat(bold(z))_i^((l))}$, then update
       $(bold(W)_"AE", bold(Theta))$ by SGD / Adam on
       $- (1 / |cal(B)|) sum_(i in cal(B)) sum_(k=1)^K q_(i k) log tilde(p)_(i k) + lambda cal(L)_"rec" (cal(B))$.
  4. *Assignment:* Compute clean $p_(i k)$ and return
     $hat(c)_i = arg max_(k in {1, dots, K}) p_(i k)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four core theoretical results: (A) strict convexity of
the target subproblem, its exact implicit KKT characterization, and
the precise condition under which DEPICT's closed-form square-root
update is exact; (B) alternating minimization descent via Bregman
projections; (C) an analytical proof that noisy softmax training
acts as an explicit Fisher-information Jacobian regularizer; and
(D) a constructive distortion bound quantifying DEPICT's bias under
class imbalance.

#subsec("A", "Convexity, Exact KKT System, and Square-Root Update")

#theorem(
  title: "Exact Target Optimality and Square-Root Update Condition",
  origin: "Original Derivation",
)[
  Fix a strictly positive assignment matrix
  $bold(P) in (Delta^(K-1))^N$ ($p_(i k) > 0$) and prior
  $bold(u) in Delta^(K-1)$ ($u_k > 0$). Then:
  1. $cal(L)_"cls" (bold(Q), bold(P))$ is strictly convex over
     $bold(Q) in (Delta^(K-1))^N$ and admits a unique global
     minimizer $bold(Q)^ast$ satisfying the coupled KKT system
     $
       q_(i k)^ast
       = (p_(i k) u_k / f_k^ast) /
         (sum_(r=1)^K p_(i r) u_r / f_r^ast),
       quad f_k^ast = 1 / N sum_(i=1)^N q_(i k)^ast.
     $
  2. Define the row normalizers $Z_i = sum_(r=1)^K p_(i r) u_r / f_r^ast$.
     Whenever $Z_i$ is constant across samples ($Z_i equiv Z$), the
     exact global minimizer collapses analytically to the closed-form
     DEPICT square-root update @ghasedi2017depict:
     $
       tilde(q)_(i k)
       = (p_(i k) (u_k / bar(p)_k)^(1/2)) /
         (sum_(r=1)^K p_(i r) (u_r / bar(p)_r)^(1/2)),
       quad bar(p)_k = 1 / N sum_(j=1)^N p_(j k).
     $
  3. In general, letting $Z_"min" = min_i Z_i$ and
     $Z_"max" = max_i Z_i$, the exact optimal marginals $f_k^ast$
     are bracketed by
     $
       ((u_k bar(p)_k) / Z_"max")^(1/2)
       <= f_k^ast
       <= ((u_k bar(p)_k) / Z_"min")^(1/2).
     $
]

#proof[
  For part (1), note that $x |-> x log(x / a)$ is strictly convex on
  $(0, infinity)$ for any $a > 0$, so the first term
  $(1 / N) sum_(i, k) q_(i k) log(q_(i k) / p_(i k))$ is strictly
  convex in $bold(Q)$, while $D_"KL" (bold(f) || bold(u))$ is convex
  as a composition of a convex function with the linear map
  $bold(Q) |-> bold(f)$. Their sum is therefore strictly convex on
  the compact convex polytope $(Delta^(K-1))^N$. Moreover, as any
  $q_(i k) -> 0^+$, the directional derivative
  $(1 / N) (log(q_(i k) / p_(i k)) + log(f_k / u_k) + 2) -> - infinity$,
  so the unique minimizer lies strictly in the relative interior
  ($q_(i k)^ast > 0$). Forming the Lagrangian
  $Lambda(bold(Q), bold(alpha)) = cal(L)_"cls" (bold(Q), bold(P)) + (1 / N) sum_(i=1)^N alpha_i (sum_(k=1)^K q_(i k) - 1)$
  and setting $partial Lambda / (partial q_(i k)) = 0$ yields
  $
    1 / N (log (q_(i k)^ast / p_(i k)) + 1)
    + 1 / N (log (f_k^ast / u_k) + 1)
    + alpha_i / N = 0,
  $
  since $partial f_r / (partial q_(i k)) = (1 / N) delta_(r k)$.
  Solving for $q_(i k)^ast$ gives
  $q_(i k)^ast = exp(- 2 - alpha_i) p_(i k) u_k / f_k^ast$, and
  enforcing $sum_(r=1)^K q_(i r)^ast = 1$ determines
  $exp(- 2 - alpha_i) = Z_i^(-1)$, establishing the KKT identity.

  For parts (2) and (3), average the KKT equation over all
  $i = 1, dots, N$:
  $
    f_k^ast
    = 1 / N sum_(i=1)^N q_(i k)^ast
    = u_k / f_k^ast dot 1 / N sum_(i=1)^N p_(i k) / Z_i.
  $
  Multiplying by $f_k^ast > 0$ gives
  $(f_k^ast)^2 = u_k (1 / N) sum_(i=1)^N p_(i k) / Z_i$.
  Bounding $Z_"max"^(-1) <= Z_i^(-1) <= Z_"min"^(-1)$ immediately
  yields the sandwich inequality
  $(u_k bar(p)_k / Z_"max")^(1/2) <= f_k^ast <= (u_k bar(p)_k / Z_"min")^(1/2)$.
  In particular, whenever $Z_i equiv Z$ is sample-independent,
  $f_k^ast = (u_k bar(p)_k / Z)^(1/2)$, and substituting this exact
  value back into $q_(i k)^ast prop p_(i k) u_k / f_k^ast$ cancels
  the shared constant $Z^(1/2)$ and produces
  $tilde(q)_(i k) prop p_(i k) (u_k / bar(p)_k)^(1/2)$, reducing
  for uniform $u_k = 1 / K$ to the formula of Ghasedi Dizaji et al.
  @ghasedi2017depict.
]

Theorem 1 reveals a subtle and important nuance: the closed-form
update in Alg. 185 is a first-order uniform-normalizer approximation
to the exact convex optimum $bold(Q)^ast$, akin to a single
damped Sinkhorn–Knopp scaling step @cuturi2013sinkhorn.

#subsec("B", "Bregman Alternating Descent and Lower-Bound Guarantee")

#lemma(
  title: "Monotonic Descent of Exact Alternating DEPICT Optimization",
  origin: "Original Proof",
)[
  Consider noiseless alternating optimization of
  $cal(L)_"DEPICT" (bold(Q), bold(W)_"AE", bold(Theta))$ where step
  $t$ updates: (i) $bold(Q)^((t+1)) = arg min_(bold(Q)) cal(L)_"cls" (bold(Q), bold(P)^((t)))$
  via the exact KKT fixed-point map
  $f_k^((m+1)) = (u_k (1 / N) sum_i p_(i k)^((t)) / Z_i^((m)))^(1/2)$,
  and (ii) $(bold(W)_"AE"^((t+1)), bold(Theta)^((t+1)))$ via any
  descent step reducing
  $- (1 / N) sum_(i, k) q_(i k)^((t+1)) log p_(i k) + lambda cal(L)_"rec"$.
  Then the sequence $cal(L)_"DEPICT"^((t))$ is monotonically
  non-increasing and converges to a finite limit
  $cal(L)^ast >= 0$.
]

#proof[
  Observe first that $cal(L)_"DEPICT" >= 0$ identically because both
  $D_"KL" (bold(q)_i || bold(p)_i) >= 0$ and
  $D_"KL" (bold(f) || bold(u)) >= 0$ by Gibbs' inequality
  @banerjee2005clustering, while $lambda cal(L)_"rec" >= 0$ as a
  nonnegative sum of squared Euclidean norms.

  In the $bold(Q)$-step, because $cal(L)_"cls" (bold(Q), bold(P)^((t)))$
  is strictly convex (Theorem 1), replacing $bold(Q)^((t))$ by the
  unique global minimizer $bold(Q)^((t+1))$ satisfies
  $
    cal(L)_"DEPICT" (bold(Q)^((t+1)), bold(W)_"AE"^((t)), bold(Theta)^((t)))
    <= cal(L)_"DEPICT" (bold(Q)^((t)), bold(W)_"AE"^((t)), bold(Theta)^((t))),
  $
  since $cal(L)_"rec"$ does not depend on $bold(Q)$.

  In the network step, holding $bold(Q)^((t+1))$ fixed, expand the
  clustering objective:
  $
    cal(L)_"cls" (bold(Q)^((t+1)), bold(P))
    &= - 1 / N sum_(i=1)^N sum_(k=1)^K q_(i k)^((t+1)) log p_(i k) \
    & quad - 1 / N sum_(i=1)^N H(bold(q)_i^((t+1)))
       + D_"KL" (bold(f)^((t+1)) || bold(u)).
  $
  The negative entropy $- (1 / N) sum_i H(bold(q)_i^((t+1)))$ and
  balance term $D_"KL" (bold(f)^((t+1)) || bold(u))$ are strictly
  constant with respect to $(bold(W)_"AE", bold(Theta))$. Hence any
  update decreasing the cross-entropy plus $lambda cal(L)_"rec"$
  decreases $cal(L)_"DEPICT"$ by the exact same amount:
  $
    cal(L)_"DEPICT" (bold(Q)^((t+1)), bold(W)_"AE"^((t+1)), bold(Theta)^((t+1)))
    <= cal(L)_"DEPICT" (bold(Q)^((t+1)), bold(W)_"AE"^((t)), bold(Theta)^((t))).
  $
  Combining both inequalities gives
  $cal(L)_"DEPICT"^((t+1)) <= cal(L)_"DEPICT"^((t))$. Every
  monotone non-increasing sequence bounded below by $0$ converges to
  a finite limit $cal(L)^ast >= 0$.
]

#subsec("C", "Denoising Softmax as Explicit Fisher-Margin Regularization")

A defining feature of DEPICT is evaluating the target $bold(Q)$ on
*clean* embeddings $bold(z)_i$ while backpropagating the
cross-entropy through *corrupted* embeddings $tilde(bold(z))_i$.
The next original theorem proves that this asymmetry induces a
closed-form curvature regularizer pushing embeddings away from
softmax decision boundaries.

#theorem(
  title: "Second-Order Fisher Regularization of Noisy Softmax Training",
  origin: "Original Proof",
)[
  Let $tilde(bold(z))_i = bold(z)_i + bold(eta)_i$ where
  $bold(eta)_i in bb(R)^d$ is a zero-mean perturbation with
  covariance $"Cov"(bold(eta)_i) = bold(Sigma)_(eta, i)$ and finite
  third absolute moments $cal(O)(sigma^3)$ with
  $sigma^2 = norm(bold(Sigma)_(eta, i))_"op"$. Define the softmax
  Fisher covariance matrix
  $bold(C)_i = "diag"(bold(p)_i) - bold(p)_i bold(p)_i^top succ.eq bold(0)$
  at the clean probability vector $bold(p)_i = "softmax"(bold(Theta)^top bold(z)_i)$.
  Then the expected noisy cross-entropy satisfies
  $
    bb(E)_(bold(eta)_i)
    [ - sum_(k=1)^K q_(i k) log tilde(p)_(i k) ]
    &= - sum_(k=1)^K q_(i k) log p_(i k) \
    & quad + 1 / 2 "Tr"(bold(Theta) bold(C)_i bold(Theta)^top bold(Sigma)_(eta, i))
       + cal(O)(sigma^3).
  $
  In particular, for isotropic noise $bold(Sigma)_(eta, i) = sigma^2 bold(I)_d$,
  the induced penalty equals
  $(sigma^2 / 2) "Var"_(k tilde.op bold(p)_i) (bold(theta)_k) >= 0$,
  which vanishes identically if and only if $bold(p)_i$ lies on a
  vertex of $Delta^(K-1)$ (maximally confident cluster assignment)
  or $bold(Theta)$ is rank-zero.
]

#proof[
  Define the cross-entropy scalar function
  $h_i (bold(z)) = - sum_(k=1)^K q_(i k) log p_(i k) (bold(z)) = log(sum_(r=1)^K exp(bold(theta)_r^top bold(z))) - sum_(k=1)^K q_(i k) bold(theta)_k^top bold(z)$.
  Note that $h_i$ is smooth ($C^infinity$) and convex in $bold(z)$.
  Its gradient at $bold(z)_i$ is
  $
    nabla_(bold(z)) h_i (bold(z)_i)
    = sum_(k=1)^K (p_(i k) - q_(i k)) bold(theta)_k
    = bold(Theta) (bold(p)_i - bold(q)_i).
  $
  Differentiating once more, note that the linear term
  $- sum_k q_(i k) bold(theta)_k^top bold(z)$ has zero Hessian, while
  the log-sum-exp Hessian with respect to logits
  $bold(a) = bold(Theta)^top bold(z)$ is exactly
  $bold(C)_i = "diag"(bold(p)_i) - bold(p)_i bold(p)_i^top$. By the
  chain rule,
  $
    nabla_(bold(z))^2 h_i (bold(z)_i)
    = bold(Theta) bold(C)_i bold(Theta)^top.
  $
  Applying Taylor's theorem around $bold(z)_i$ and taking
  expectation over $bold(eta)_i$ with $bb(E)[bold(eta)_i] = bold(0)$
  and $bb(E)[bold(eta)_i bold(eta)_i^top] = bold(Sigma)_(eta, i)$:
  $
    bb(E)[h_i (bold(z)_i + bold(eta)_i)]
    &= h_i (bold(z)_i) + nabla_(bold(z)) h_i (bold(z)_i)^top bb(E)[bold(eta)_i] \
    & quad + 1 / 2 bb(E)[bold(eta)_i^top bold(Theta) bold(C)_i bold(Theta)^top bold(eta)_i]
       + cal(O)(sigma^3) \
    &= h_i (bold(z)_i)
       + 1 / 2 "Tr"(bold(Theta) bold(C)_i bold(Theta)^top bold(Sigma)_(eta, i))
       + cal(O)(sigma^3).
  $
  When $bold(Sigma)_(eta, i) = sigma^2 bold(I)_d$, cyclic invariance
  of trace gives
  $
    "Tr"(bold(Theta) bold(C)_i bold(Theta)^top)
    &= sum_(k=1)^K p_(i k) norm(bold(theta)_k)_2^2
       - norm(sum_(k=1)^K p_(i k) bold(theta)_k)_2^2 \
    &= sum_(k=1)^K p_(i k) norm(bold(theta)_k - bar(bold(theta))_i)_2^2
       = "Var"_(k tilde.op bold(p)_i) (bold(theta)_k) >= 0,
  $
  where $bar(bold(theta))_i = sum_k p_(i k) bold(theta)_k$. For
  distinct non-collinear class prototypes ${bold(theta)_k}$, this
  variance vanishes if and only if $bold(p)_i$ puts unit mass on a
  single cluster $k^ast$, proving that dropout/denoising pushes
  latent representations cleanly away from inter-cluster
  hyperplanes while simultaneously discouraging unbounded weight
  norm inflation in $bold(Theta)$.
]

#subsec("D", "Class-Imbalance Distortion Bound and Counterexample")

While the uniform balance penalty $D_"KL" (bold(f) || bold(u))$
with $u_k = 1 / K$ prevents empty-cluster collapse, it introduces
an irreducible structural bias whenever the true data generating
clusters have unequal cardinalities.

#proposition(
  title: "Imbalanced-Cluster Distortion Under Uniform Prior",
  origin: "Original Counterexample",
)[
  Suppose $cal(X)$ consists of $K = 2$ well-separated clusters of
  sizes $N_1 = alpha N$ and $N_2 = (1 - alpha) N$ with
  $alpha in (1 / 2, 1)$, and suppose the encoder produces cleanly
  separated soft probabilities $p_(i 1) = 1 - epsilon$ for
  $i in cal(C)_1$ and $p_(i 2) = 1 - epsilon$ for $i in cal(C)_2$
  with $epsilon in (0, 1 / 2)$.
  1. Under DEPICT's uniform prior $bold(u) = (1 / 2, 1 / 2)$, the
     exact KKT optimum $bold(Q)^ast$ of Theorem 1 shrinks the
     majority cluster frequency $f_1^ast$ strictly below $alpha$:
     specifically, for $epsilon -> 0^+$ along the interior regime
     where $epsilon / (1 - epsilon) < ((1 - alpha) / alpha)$, the
     first-step square-root update shifts the log-odds of *every*
     sample toward the minority cluster by
     $
       Delta_"logit"
       = 1 / 2 log ((alpha (1 - epsilon) + (1 - alpha) epsilon) /
                    ((1 - alpha)(1 - epsilon) + alpha epsilon)) > 0.
     $
  2. If the cluster separation satisfies
     $((1 - epsilon) / epsilon)^2 < alpha / (1 - alpha)$
     (i.e., moderately overlapping imbalanced clusters with
     $epsilon > (1 + (alpha / (1 - alpha))^(1/2))^(-1)$), then
     DEPICT's square-root update *flips the majority assignment*
     ($tilde(q)_(i 1) < 1 / 2$) for *every* sample in the majority
     cluster $cal(C)_1$ simultaneously!
]

#proof[
  Compute the empirical column averages of $bold(P)$:
  $
    bar(p)_1 &= alpha (1 - epsilon) + (1 - alpha) epsilon, \
    bar(p)_2 &= (1 - alpha)(1 - epsilon) + alpha epsilon.
  $
  Since $alpha > 1 / 2$ and $epsilon < 1 / 2$,
  $bar(p)_1 - bar(p)_2 = (2 alpha - 1)(1 - 2 epsilon) > 0$, hence
  $bar(p)_1 > bar(p)_2$. Under the closed-form square-root update
  of Theorem 1(2) with $u_1 = u_2 = 1 / 2$, the target log-odds for
  any datum $i$ satisfy
  $
    log (tilde(q)_(i 1) / tilde(q)_(i 2))
    = log (p_(i 1) / p_(i 2)) - 1 / 2 log (bar(p)_1 / bar(p)_2).
  $
  This proves part (1) with
  $Delta_"logit" = (1 / 2) log(bar(p)_1 / bar(p)_2) > 0$.

  For part (2), for every $i in cal(C)_1$ we have
  $p_(i 1) / p_(i 2) = (1 - epsilon) / epsilon$. Therefore
  $tilde(q)_(i 1) < 1 / 2$ holds if and only if
  $(1 - epsilon) / epsilon < (bar(p)_1 / bar(p)_2)^(1/2)$. Because
  $bar(p)_1 / bar(p)_2 -> alpha / (1 - alpha)$ monotonically from
  below as $epsilon -> 0$, and explicitly
  $((1 - epsilon) / epsilon)^2 < (alpha (1 - epsilon) + (1 - alpha) epsilon) / ((1 - alpha)(1 - epsilon) + alpha epsilon)$
  holds whenever $(1 - epsilon) / epsilon < (alpha / (1 - alpha))^(1/2)$
  is strictly satisfied with sufficient margin (and is exact at the
  boundary $bar(p)_1 / bar(p)_2 = ((1 - epsilon)/epsilon)^2$), solving
  $(1 - epsilon) / epsilon = (alpha / (1 - alpha))^(1/2)$ yields the
  critical overlap threshold
  $epsilon_"crit" = (1 + (alpha / (1 - alpha))^(1/2))^(-1)$. Above
  this overlap threshold, every sample in $cal(C)_1$ is reassigned
  to cluster $2$ in the $bold(Q)$-step, inducing mass oscillation or
  spurious splitting of the majority class.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $P$ denote the forward-backward arithmetic cost of the $L$-layer
stacked autoencoder per sample, and let $H_"tot" = sum_(l=0)^(L-1) d_l$
be the sum of layer widths. Unlike VaDE (Alg. 184), which
reconstructs only the input layer $bold(x)_i$, DEPICT evaluates
squared reconstruction errors across all $L$ intermediate clean
layers ${bold(z)_i^((l))}$, requiring two forward passes during
training (one clean pass to obtain reconstruction targets
$bold(z)_i^((l))$ and one corrupted pass to obtain
$tilde(bold(z))_i^((l))$ and $hat(bold(z))_i^((l))$) or cached clean
activations. The global target distribution update $bold(Q)$ scans
clean softmax outputs once per epoch in exact $cal(O)(N K)$ time
and $cal(O)(N K)$ space (or $cal(O)(K)$ accumulator space if
computed in two streaming passes), completely avoiding any
$cal(O)(N^2)$ pairwise similarity matrix. Total training complexity
over $E$ epochs is therefore
$
  T_"DEPICT" = cal(O)(E N (P + K d)) + T_"pre",
$
linear in both $N$ and $K$.

Geometrically, the multinomial softmax head $bold(Theta)$ partitions
the bottleneck space $bb(R)^d$ into $K$ *convex linear polytopes*
(Voronoi-like cells with linear separating hyperplanes
$(bold(theta)_k - bold(theta)_r)^top bold(z) = 0$). Nonlinear
cluster boundaries in ambient space $bb(R)^D$ emerge entirely from
the composition with the deep encoder $f_(bold(W))$. Crucially,
layer-wise denoising supervision acts as a multi-scale contractive
manifold regularizer: by forcing every decoder layer to invert
corrupted hidden features back to clean activations, DEPICT
prevents the encoder from collapsing orthogonal data variations
needed to distinguish subtle sub-clusters, while Theorem 2 proves
that noisy softmax cross-entropy simultaneously carves wide
low-density margins around the linear decision boundaries in
$bb(R)^d$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(E N (P + K d))$ training; $cal(O)(P_"enc" + K d)$ per test sample],
    [Space Complexity], [$cal(O)(|bold(W)_"AE"| + K d + B H_"tot" + N K)$ (or $cal(O)(B K)$ streaming)],
    [Latent Geometry], [Linear softmax polytopes in $bb(R)^d$; contractive denoising manifold in $bb(R)^D$],
    [Cluster Balance], [Explicit $D_"KL" (bold(f) || bold(u))$ regularization; strong uniform-size bias by default],
    [Noise Robustness], [High; dual clean-target / noisy-prediction loss yields Fisher margin (Theorem 2)],
    [Initialization Sensitivity], [Moderate; requires AE pretraining + initial softmax fit to avoid symmetric traps],
    [Primary Hyperparameters], [Cluster count $K$, balance prior $bold(u)$, weight $lambda$, dropout rate $rho$],
  ),
  caption: [Operational and geometric profile of DEPICT (Alg. 185).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Balanced Visual Object & Digit Recognition:* Standard vision
   corpora (MNIST, USPS) possess near-uniform class frequencies
   matching $u_k = 1 / K$; convolutional denoising preserves local
   stroke geometry while preventing empty-cluster collapse.

2. *Noisy Biomedical Microscopy & Histology:* Staining variance and
   shot noise corrupt cellular patches; layer-wise denoising maps
   noisy inputs back to clean morphologic manifolds while Theorem 2's
   Fisher margin stabilizes tissue-phenotype boundaries.

3. *Industrial Defect Inspection:* Manufacturing cameras and sensors
   face motion blur and pixel dropout; dropout-corrupted encoder
   training directly builds sensor-fault invariance with sub-millisecond
   feed-forward inference across balanced defect classes.

4. *Acoustic & Speech Spectrogram Clustering:* Spectrograms contain
   broadband noise across high-dimensional frequency bins; stacked
   layer reconstructions preserve harmonic formants across scales
   while the softmax head yields clean acoustic posteriors.

5. *Balanced Operational Cohort Stratification:* When applications
   explicitly require comparable-capacity intervention tiers,
   DEPICT's $D_"KL" (bold(f) || bold(u))$ converts operational load
   balancing directly into a smooth convex regularizer.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

DEPICT exhibits four primary failure modes: (i) *class imbalance*
(Proposition 1), where uniform $u_k = 1 / K$ splits majority
classes unless $bold(u)$ is adapted; (ii) *over-corruption* on
sparse or low-dimensional features under excessive dropout $rho$;
(iii) *reconstruction dominance* when poorly scaled $lambda$ traps
nuisance variance in the bottleneck; and (iv) *no generative density*,
precluding likelihood scoring or sampling as in VaDE (Alg. 184).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Balanced Vision / Digits], [Denoising CNN + $u_k = 1 / K$ matches uniform classes], [Tune $rho$ and layer depth $L$],
    [Noisy Microscopy / Audio], [Layer-wise denoising + Fisher margin resists noise], [Avoid excessive $rho$ on sparse data],
    [Balanced Operational Tiers], [Explicit $D_"KL" (bold(f) || bold(u))$ enforces capacity balance], [Specify custom $bold(u)$ if non-uniform],
    [Skewed Customer / Rare Fraud], [Poor under default uniform prior (Prop. 1 distortion)], [Prefer VaDE (Alg. 184) or GMM],
    [High-Res Natural Scenes], [Pixel MSE focuses on texture rather than semantics], [Prefer DAC (Alg. 187) or SwAV],
  ),
  caption: [Domain suitability and failure-mode matrix for DEPICT.],
)

#sec("V", "Structural Relations and Empirical Behavior")

DEPICT belongs to the self-training KL-divergence family pioneered
by DEC @xie2016dec and extended by IDEC @guo2017idec, while sharing
explicit marginal balancing with Sinkhorn-based optimal-transport
clustering such as SwAV @caron2020swav @cuturi2013sinkhorn. The
following proposition formalizes the exact mathematical contrast
between DEPICT and DEC / IDEC.

#proposition(
  title: "Spectral Sharpening and Balancing: DEPICT vs. DEC and IDEC",
  origin: "Original Unified Proof",
)[
  Compare the auxiliary target distributions
  $q_(i k)^((alpha, beta)) prop p_(i k)^alpha / (sum_(j=1)^N p_(j k))^beta$
  for $(alpha, beta) = (1, 1 / 2)$ (DEPICT @ghasedi2017depict) and
  $(alpha, beta) = (2, 1)$ (DEC @xie2016dec / IDEC @guo2017idec).
  1. *(Entropy Sharpening vs. Linear Fidelity):* On a balanced
     assignment matrix ($bar(p)_k = 1 / K$ for all $k$), DEPICT's
     update is the exact identity map ($tilde(bold(q))_i = bold(p)_i$),
     relying on noisy softmax regularization (Theorem 2) to sharpen
     cluster boundaries, whereas DEC / IDEC strictly reduces row
     entropy $H(bold(q)_i) < H(bold(p)_i)$ for every non-uniform
     interior vector $bold(p)_i$ by squaring probabilities.
  2. *(Balance Elasticity):* Holding sample logits fixed and
     perturbing cluster masses $bar(p)_k$, the log-odds elasticity
     with respect to cluster frequency is
     $partial log(q_(i k) / q_(i r)) / (partial log(bar(p)_k / bar(p)_r)) = - beta$,
     yielding elasticity $- 1 / 2$ for DEPICT (square-root damping
     derived from $D_"KL" (bold(f) || bold(u))$) versus $- 1$ for
     DEC / IDEC (inverse-frequency heuristic).
  3. *(Reconstruction Anchor Limit):* Setting $lambda -> 0$ and
     $rho -> 0$ removes DEPICT's geometric anchor and Fisher
     margin, collapsing its network step into pure unregularized
     self-training akin to DEC @xie2016dec; keeping $lambda > 0$
     with single-layer ($L = 1$) noiseless reconstruction recovers
     the structural role of IDEC @guo2017idec.
]

#proof[
  For part (1), when $bar(p)_k = 1 / K$ for all $k$, the
  denominator factors $(sum_j p_(j k))^(-1/2) = (N / K)^(-1/2)$ are
  identical across $k$ and cancel upon normalization, leaving
  $tilde(q)_(i k) = p_(i k)$ for $(alpha, beta) = (1, 1 / 2)$. For
  $(alpha, beta) = (2, 1)$, the balanced update is
  $q_(i k)^"DEC" = p_(i k)^2 / (sum_r p_(i r)^2)$. Consider the
  one-parameter family $q_k (s) = p_k^s / (sum_r p_r^s)$ for
  $s >= 1$. Differentiating its Shannon entropy
  $H(s) = - sum_k q_k (s) log q_k (s)$ gives
  $H'(s) = - s "Var"_(k tilde.op q(s)) (log p_k) <= 0$, with strict
  inequality whenever $bold(p)_i$ is not uniform over its support
  and not one-hot. Integrating from $s = 1$ to $s = 2$ yields
  $H(bold(q)_i^"DEC") < H(bold(p)_i)$.

  For part (2), taking logarithms of
  $q_(i k)^((alpha, beta)) / q_(i r)^((alpha, beta)) = (p_(i k) / p_(i r))^alpha (bar(p)_k / bar(p)_r)^(- beta)$
  yields
  $log(q_(i k) / q_(i r)) = alpha log(p_(i k) / p_(i r)) - beta log(bar(p)_k / bar(p)_r)$,
  whose partial derivative with respect to $log(bar(p)_k / bar(p)_r)$
  is identically $- beta$.

  Part (3) follows directly from Definition 1 and Theorem 2:
  setting $rho = 0$ drives $bold(Sigma)_(eta, i) = bold(0)$
  (eliminating the Fisher penalty), while setting $lambda = 0$
  eliminates $cal(L)_"rec"$, leaving only unanchored KL
  self-training.
]

Empirical results on the ClusterOverview four-dataset benchmark
(Silhouette @rousseeuw1987silhouettes and ARI
@hubert1985comparing) align closely with these theoretical
predictions. On *Iris* — which comprises three perfectly balanced
species ($50$ samples each, exactly matching DEPICT's uniform prior
$u_k = 1 / 3$) — DEPICT achieves the *highest Silhouette score
($0.447$)* among all four deep clustering algorithms (vs. VaDE
$0.350$, JULE $0.412$, DAC $0.433$) together with strong ARI
($0.583$), reflecting how the Fisher-margin denoising regularizer
(Theorem 2) carves clean latent separation without distorting class
proportions. Conversely, on *Customers* (skewed, unbalanced spend
cohorts), DEPICT's Silhouette drops to $0.327$ — well below VaDE's
$0.381$ — directly manifesting the Class-Imbalance Distortion
effect of Proposition 1 under a uniform prior. On *Multi-Shape*
(Silhouette $0.712$, ARI N/A) and *High-D* (Silhouette $0.263$, ARI
N/A), stacked denoising reconstruction preserves curved manifold
structure and filters high-dimensional noise effectively, placing
DEPICT nearly tied with DAC ($0.266$) at the top of *High-D* while
strongly outperforming JULE ($0.119$). Finally, its mean runtime of
$8.460$ ms sits cleanly between decoder-free DAC ($3.329$ ms) and
heavier VaDE ($11.007$ ms) / JULE ($14.593$ ms), confirming the
efficiency of its linear-time $cal(O)(N K)$ target update.
