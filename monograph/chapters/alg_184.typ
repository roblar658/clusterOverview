#import "../template/ieee.typ": *

#alg_entry(
  id: 184,
  name: "VaDE (Variational Deep Embedding)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2017",
  paradigm: "Deep Generative Variational Mixture Model",
  complexity: $cal(O)(E N (P + K d))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical Gaussian Mixture Models (GMMs) @mclachlan2000finite
provide a principled probabilistic framework for density estimation
and soft clustering, yet their expressive capacity degrades
severely when observations $bold(x)_i in bb(R)^D$ lie on curved,
low-dimensional nonlinear manifolds embedded in high-dimensional
ambient spaces. Conversely, deep Variational Autoencoders (VAEs)
@kingma2014vae learn smooth nonlinear latent representations via
amortized variational inference, but their standard unimodal
isotropic prior $p(bold(z)) = cal(N)(bold(0), bold(I)_d)$
explicitly discourages multimodal cluster structure in the latent
space. Variational Deep Embedding (VaDE), introduced by Jiang et
al. @jiang2017vade, unifies both paradigms inside a single directed
probabilistic graphical model by replacing the unimodal VAE prior
with a $K$-component Gaussian mixture over the latent manifold.

Let $cal(X) = {bold(x)_1, dots, bold(x)_N} subset.eq bb(R)^D$ be an
unlabeled dataset and let $K >= 2$ denote the target cluster count.
VaDE posits a three-stage generative process parameterized by
decoder weights $theta$ and latent mixture parameters
$Psi = {pi_c, bold(mu)_c, bold(sigma)_c^2}_(c=1)^K$, where
$bold(pi) in Delta^(K-1)$ lies in the $(K-1)$-simplex
($pi_c > 0$, $sum_(c=1)^K pi_c = 1$), $bold(mu)_c in bb(R)^d$, and
$bold(sigma)_c^2 in bb(R)_(>0)^d$ encodes a diagonal component
covariance $"diag"(bold(sigma)_c^2)$:
1. Draw a discrete cluster assignment $c tilde.op "Cat"(bold(pi))$
   with prior probability $p_Psi (c) = pi_c$.
2. Draw a latent representation from the selected Gaussian mode:
   $
     bold(z) | c tilde.op cal(N)(bold(mu)_c, "diag"(bold(sigma)_c^2)).
   $
3. Generate the observation $bold(x)$ through a nonlinear decoder
   network $g_theta : bb(R)^d -> bb(R)^D$ defining the likelihood
   $p_theta (bold(x) | bold(z))$ — either a factorized Bernoulli
   model $"Bern"(bold(mu)_x (bold(z)))$ for binary or
   $[0, 1]$-normalized data, or a diagonal Gaussian
   $cal(N)(bold(mu)_x (bold(z)), "diag"(bold(sigma)_x^2 (bold(z))))$
   for continuous features.

The resulting joint density factorizes cleanly as
$
  p_(theta, Psi) (bold(x), bold(z), c)
  = p_theta (bold(x) | bold(z)) p_Psi (bold(z) | c) p_Psi (c).
$
Because integrating out $(bold(z), c)$ to obtain the marginal
log-likelihood $log p_(theta, Psi) (bold(x))$ is intractable for
nonlinear decoders $g_theta$, VaDE introduces an amortized
variational posterior $q_phi (bold(z), c | bold(x))$ parameterized
by an inference (encoder) network $f_phi : bb(R)^D -> bb(R)^(2 d)$
that outputs mean and log-variance vectors
$(tilde(bold(mu))_phi (bold(x)), log tilde(bold(sigma))_phi^2 (bold(x)))$.

#definition(title: "VaDE Variational Objective (ELBO)")[
  Under the mean-field variational factorization
  $q_phi (bold(z), c | bold(x)) = q_phi (bold(z) | bold(x)) q(c | bold(x))$
  with Gaussian recognition density
  $q_phi (bold(z) | bold(x)) = cal(N)(tilde(bold(mu))_phi (bold(x)), "diag"(tilde(bold(sigma))_phi^2 (bold(x))))$
  and categorical responsibility $gamma_c (bold(x)) = q(c | bold(x))$,
  the *VaDE Evidence Lower Bound* (ELBO) for a single datum
  $bold(x)$ is
  $
    cal(L)_"ELBO" (theta, phi, Psi; bold(x))
    &= bb(E)_(q_phi (bold(z), c | bold(x)))
       [ log (p_(theta, Psi) (bold(x), bold(z), c)) / (q_phi (bold(z), c | bold(x))) ] \
    &= bb(E)_(q_phi (bold(z) | bold(x))) [ log p_theta (bold(x) | bold(z)) ] \
    & quad - D_"KL" (q_phi (bold(z), c | bold(x)) || p_Psi (bold(z), c)).
  $
  The empirical objective over $cal(X)$ is
  $hat(cal(L))(theta, phi, Psi) = (1 / N) sum_(i=1)^N cal(L)_"ELBO" (theta, phi, Psi; bold(x)_i)$.
]

Applying the reparameterization trick @kingma2014vae, latent
samples are drawn as
$bold(z)^((l)) = tilde(bold(mu))_phi (bold(x)_i) + tilde(bold(sigma))_phi (bold(x)_i) dot.o bold(epsilon)^((l))$
with $bold(epsilon)^((l)) tilde.op cal(N)(bold(0), bold(I)_d)$ and
Hadamard product $dot.o$, enabling low-variance Stochastic Gradient
Variational Bayes (SGVB) backpropagation through the Monte Carlo
reconstruction term while evaluating the GMM-prior KL divergence in
closed form. Because deep generative models suffer from early
posterior collapse when trained from random weights, VaDE employs
stacked autoencoder pretraining followed by classical GMM EM
@dempster1977em in the pretrained latent space to initialize $Psi$.

#algbox(
  title: "VaDE — Variational Deep Embedding",
  inputs: [
    Dataset $cal(X) = {bold(x)_i}_(i=1)^N subset.eq bb(R)^D$,
    clusters $K$, latent dimension $d$, MC samples $L$,
    epochs $E$, learning rate $eta$.
  ],
  outputs: [
    Parameters $(theta, phi, Psi)$, soft assignments
    ${gamma_(i c)}$, hard cluster labels ${hat(c)_i}_(i=1)^N$.
  ],
)[
  1. *Pretraining:* Train a deterministic autoencoder $(f_phi, g_theta)$
     on $cal(X)$ (MSE / BCE) to initialize $tilde(bold(mu))_phi$ and $g_theta$.
  2. *Prior Initialization:* Fit a $K$-component diagonal GMM via EM
     @dempster1977em on ${tilde(bold(mu))_phi (bold(x)_i)}$ to initialize $Psi$.
  3. *For* epoch $e = 1, dots, E$ *and mini-batch* $cal(B) subset.eq cal(X)$:
     - Compute $(tilde(bold(mu))_i, tilde(bold(sigma))_i^2) = f_phi (bold(x)_i)$
       and sample $bold(z)_i^((l)) = tilde(bold(mu))_i + tilde(bold(sigma))_i dot.o bold(epsilon)_i^((l))$,
       $bold(epsilon)_i^((l)) tilde.op cal(N)(bold(0), bold(I)_d)$ ($l = 1, dots, L$).
     - Evaluate optimal responsibilities $gamma_(i c)$ (Thm. 1 / Prop. 1)
       and update $(theta, phi, Psi)$ via Adam on the SGVB ELBO,
       enforcing $bold(pi) in Delta^(K-1)$ and $bold(sigma)_c^2 > bold(0)$.
  4. *Assignment:* Return $hat(c)_i = arg max_(c in {1, dots, K}) gamma_(i c)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing VaDE: (A) its
exact analytical ELBO decomposition and the variational gap between
mean-field and structured posteriors, (B) exact closed-form M-step
monotonicity and automatic covariance regularization, (C) a
constructive counterexample characterizing posterior collapse and
its instability threshold, and (D) exact latent gauge symmetry.

#subsec("A", "Analytical ELBO and Variational Factorization Gap")

#theorem(
  title: "Exact Closed-Form VaDE ELBO and Optimal Posterior",
  origin: "Classical — Cited",
)[
  Let $q_phi (bold(z), c | bold(x)) = q_phi (bold(z) | bold(x)) q(c | bold(x))$
  with $q_phi (bold(z) | bold(x)) = cal(N)(tilde(bold(mu)), "diag"(tilde(bold(sigma))^2))$
  and $gamma_c = q(c | bold(x))$ @jiang2017vade @kingma2014vae. Then:
  1. For any valid $(theta, phi, Psi)$ and any $bold(gamma) in Delta^(K-1)$,
     $log p_(theta, Psi) (bold(x)) >= cal(L)_"ELBO" (theta, phi, Psi; bold(x))$,
     with exact gap $D_"KL" (q_phi (bold(z), c | bold(x)) || p_(theta, Psi) (bold(z), c | bold(x)))$.
  2. The analytical ELBO takes the closed form
     $
       cal(L)_"ELBO"
       &= bb(E)_(q_phi) [ log p_theta (bold(x) | bold(z)) ]
          + sum_(c=1)^K gamma_c log (pi_c / gamma_c) \
       & quad - 1 / 2 sum_(c=1)^K gamma_c sum_(j=1)^d
          (log sigma_(c, j)^2
           + (tilde(sigma)_j^2 + (tilde(mu)_j - mu_(c, j))^2) / sigma_(c, j)^2) \
       & quad + 1 / 2 sum_(j=1)^d (1 + log tilde(sigma)_j^2).
     $
  3. Under the enriched structured posterior
     $q_"struct" (bold(z), c | bold(x)) = q_phi (bold(z) | bold(x)) q(c | bold(z))$,
     the unique functional maximizer over $q(c | bold(z))$ is the
     exact latent mixture posterior
     $
       q^ast (c | bold(z)) = p_Psi (c | bold(z))
       = (pi_c p_Psi (bold(z) | c)) / (sum_(k=1)^K pi_k p_Psi (bold(z) | k)).
     $
]

#proof[
  For part (1), apply Jensen's inequality to the marginal
  log-likelihood:
  $
    log p_(theta, Psi) (bold(x))
    &= log sum_(c=1)^K integral p_(theta, Psi) (bold(x), bold(z), c) dif bold(z) \
    &= cal(L)_"ELBO" + D_"KL" (q_phi (bold(z), c | bold(x)) || p_(theta, Psi) (bold(z), c | bold(x))) \
    &>= cal(L)_"ELBO",
  $
  since Kullback–Leibler divergence is strictly nonnegative.

  For part (2), expand the joint prior KL divergence using the
  factorization $q_phi (bold(z), c | bold(x)) = q_phi (bold(z) | bold(x)) gamma_c$:
  $
    - D_"KL" (q_phi (bold(z), c | bold(x)) || p_Psi (bold(z), c))
    &= bb(E)_(q_phi) [ sum_(c=1)^K gamma_c log p_Psi (bold(z) | c) ] \
    & quad + sum_(c=1)^K gamma_c log (pi_c / gamma_c)
       - bb(E)_(q_phi) [ log q_phi (bold(z) | bold(x)) ].
  $
  For each coordinate $j in {1, dots, d}$, since
  $z_j tilde.op cal(N)(tilde(mu)_j, tilde(sigma)_j^2)$ under $q_phi$,
  we have
  $bb(E)_(q_phi) [(z_j - mu_(c, j))^2] = tilde(sigma)_j^2 + (tilde(mu)_j - mu_(c, j))^2$.
  Substituting the diagonal Gaussian log-density
  $log p_Psi (bold(z) | c) = - (d / 2) log(2 pi) - (1 / 2) sum_j (log sigma_(c, j)^2 + (z_j - mu_(c, j))^2 / sigma_(c, j)^2)$
  and the differential entropy of $q_phi$,
  $- bb(E)_(q_phi) [log q_phi (bold(z) | bold(x))] = (d / 2) log(2 pi) + (1 / 2) sum_j (1 + log tilde(sigma)_j^2)$,
  the $(d / 2) log(2 pi)$ terms cancel identically because
  $sum_(c=1)^K gamma_c = 1$, yielding the stated closed-form ELBO.

  For part (3), write the structured KL divergence via the chain
  rule for KL divergence @bishop2006pattern:
  $
    D_"KL" (q_phi (bold(z) | bold(x)) q(c | bold(z)) || p_(theta, Psi) (bold(z), c | bold(x)))
    &= D_"KL" (q_phi (bold(z) | bold(x)) || p_(theta, Psi) (bold(z) | bold(x))) \
    & quad + bb(E)_(q_phi) [ D_"KL" (q(c | bold(z)) || p_Psi (c | bold(z))) ].
  $
  Note that $p_(theta, Psi) (c | bold(z), bold(x)) = p_Psi (c | bold(z))$
  because $bold(x) perp c | bold(z)$ by $d$-separation in the
  generative graph $c -> bold(z) -> bold(x)$. The first term is
  independent of $q(c | bold(z))$, while the second term is uniquely
  minimized to zero $q_phi$-almost everywhere if and only if
  $q(c | bold(z)) = p_Psi (c | bold(z))$.
]

In practice, implementations of VaDE @jiang2017vade often substitute
Monte Carlo evaluations of $p_Psi (c | bold(z))$ for $gamma_c$ inside
the mean-field ELBO. The next original proposition establishes the
exact relationship and strict ordering between the mean-field
optimal responsibility and the structured posterior.

#proposition(
  title: "Exact Variational Gap Between Mean-Field and Structured VaDE",
  origin: "Original Derivation",
)[
  Fix $q_phi (bold(z) | bold(x)) = cal(N)(tilde(bold(mu)), "diag"(tilde(bold(sigma))^2))$
  and mixture parameters $Psi$.
  1. The unique global maximizer of the mean-field ELBO over
     $bold(gamma) in Delta^(K-1)$ is the variance-penalized
     geometric-mean posterior
     $
       gamma_c^ast
       = (pi_c exp(bb(E)_(q_phi) [ log p_Psi (bold(z) | c) ])) /
         (sum_(k=1)^K pi_k exp(bb(E)_(q_phi) [ log p_Psi (bold(z) | k) ])).
     $
  2. Denote the resulting maximized mean-field and structured ELBOs
     by $cal(L)_"MF" (bold(x))$ and $cal(L)_"struct" (bold(x))$.
     Then
     $
       0 <= cal(L)_"struct" (bold(x)) - cal(L)_"MF" (bold(x))
       = bb(E)_(q_phi (bold(z) | bold(x)))
         [ D_"KL" (bold(gamma)^ast || p_Psi (dot | bold(z))) ],
     $
     with equality if and only if $p_Psi (c | bold(z))$ is
     $q_phi$-almost surely constant in $bold(z)$ (e.g., when
     $tilde(bold(sigma)) -> bold(0)$ or all component densities
     $p_Psi (bold(z) | c)$ coincide).
]

#proof[
  For part (1), isolate the terms in $cal(L)_"ELBO"$ depending on
  $bold(gamma)$:
  $
    J(bold(gamma))
    = sum_(c=1)^K gamma_c (log pi_c + bb(E)_(q_phi) [ log p_Psi (bold(z) | c) ] - log gamma_c).
  $
  Define unnormalized weights
  $omega_c = pi_c exp(bb(E)_(q_phi) [ log p_Psi (bold(z) | c) ]) > 0$
  and $Z_"MF" = sum_(k=1)^K omega_k$. Rewriting $J(bold(gamma))$
  gives
  $
    J(bold(gamma))
    = sum_(c=1)^K gamma_c log (omega_c / gamma_c)
    = log Z_"MF" - D_"KL" (bold(gamma) || bold(gamma)^ast),
  $
  where $gamma_c^ast = omega_c / Z_"MF"$. By strict convexity of
  $D_"KL"$, $J(bold(gamma))$ is uniquely maximized at
  $bold(gamma) = bold(gamma)^ast$ with maximum value $log Z_"MF"$.

  For part (2), note that the structured ELBO evaluated at
  $q^ast (c | bold(z)) = p_Psi (c | bold(z))$ replaces
  $bb(E)_(q_phi) [ sum_c gamma_c^ast log (p_Psi (bold(z), c) / gamma_c^ast) ]$
  with $bb(E)_(q_phi) [ log p_Psi (bold(z)) ]$, where
  $p_Psi (bold(z)) = sum_(c=1)^K p_Psi (bold(z), c)$. Subtracting
  the two bounds cancels the reconstruction expectation and
  differential entropy identically:
  $
    cal(L)_"struct" - cal(L)_"MF"
    &= bb(E)_(q_phi) [ log p_Psi (bold(z))
       - sum_(c=1)^K gamma_c^ast log (p_Psi (bold(z), c) / gamma_c^ast) ] \
    &= bb(E)_(q_phi) [ sum_(c=1)^K gamma_c^ast
       log (gamma_c^ast / (p_Psi (bold(z), c) / p_Psi (bold(z)))) ] \
    &= bb(E)_(q_phi) [ D_"KL" (bold(gamma)^ast || p_Psi (dot | bold(z))) ] >= 0.
  $
  Because $q_phi$ has full support on $bb(R)^d$ for
  $tilde(bold(sigma)) > bold(0)$, the expectation vanishes if and
  only if $p_Psi (c | bold(z)) = gamma_c^ast$ for Lebesgue-almost
  every $bold(z) in bb(R)^d$, which holds if and only if
  $p_Psi (bold(z) | c)$ is independent of $c$ (or in the
  deterministic encoder limit $tilde(bold(sigma)) -> bold(0)$).
]

#subsec("B", "Closed-Form Prior M-Step Monotonicity and Variance Floor")

While standard VaDE updates $Psi$ via stochastic gradient ascent,
analyzing the exact block-coordinate update of $Psi$ reveals a
crucial regularization mechanism absent in classical GMM EM
@mclachlan2000finite @wu1983convergence.

#lemma(
  title: "Exact Prior M-Step Monotonicity and Singularity Prevention",
  origin: "Original Proof",
)[
  Hold encoder and decoder parameters $(theta, phi)$ and
  responsibilities ${gamma_(i c)}$ fixed with effective counts
  $N_c = sum_(i=1)^N gamma_(i c) > 0$. The empirical VaDE ELBO
  $hat(cal(L))(theta, phi, Psi)$ is uniquely maximized over $Psi$ by
  $
    pi_c^+ &= N_c / N,
    quad bold(mu)_c^+ = 1 / N_c sum_(i=1)^N gamma_(i c) tilde(bold(mu))_i, \
    (sigma_(c, j)^2)^+
    &= 1 / N_c sum_(i=1)^N gamma_(i c)
       ((tilde(mu)_(i, j) - mu_(c, j)^+)^2 + tilde(sigma)_(i, j)^2).
  $
  Furthermore, if the encoder variance satisfies
  $inf_(i, j) tilde(sigma)_(i, j)^2 >= tau^2 > 0$, then
  $(sigma_(c, j)^2)^+ >= tau^2$ unconditionally for every cluster
  $c$ and coordinate $j$, preventing covariance-singularity blowup
  even when a single datum dominates component $c$.
]

#proof[
  Using Theorem 1(2), the terms of $N hat(cal(L))$ depending on
  $bold(pi)$ are $sum_(c=1)^K N_c log pi_c$. Maximizing over
  $bold(pi) in Delta^(K-1)$ via a Lagrange multiplier yields
  $pi_c^+ = N_c / N$, which is the unique global maximizer by
  strict concavity of $log$.

  For each component $c$ and coordinate $j$, let
  $v_(c, j) = sigma_(c, j)^2 > 0$ and isolate the terms depending on
  $(mu_(c, j), v_(c, j))$:
  $
    Phi_(c, j) (mu, v)
    = - 1 / 2 (N_c log v
      + 1 / v sum_(i=1)^N gamma_(i c)
        ((tilde(mu)_(i, j) - mu)^2 + tilde(sigma)_(i, j)^2)).
  $
  Setting $partial Phi_(c, j) / (partial mu) = 0$ yields the unique
  stationary point $mu_(c, j)^+ = (1 / N_c) sum_i gamma_(i c) tilde(mu)_(i, j)$,
  with $partial^2 Phi_(c, j) / (partial mu^2) = - N_c / v < 0$
  everywhere. Substituting $mu_(c, j)^+$ into $Phi_(c, j)$ and
  setting $S_(c, j) = sum_i gamma_(i c) ((tilde(mu)_(i, j) - mu_(c, j)^+)^2 + tilde(sigma)_(i, j)^2) > 0$,
  the function $h(v) = - (1 / 2)(N_c log v + S_(c, j) / v)$ has
  derivative $h'(v) = (S_(c, j) - N_c v) / (2 v^2)$, which is
  strictly positive on $(0, S_(c, j) / N_c)$ and strictly negative
  on $(S_(c, j) / N_c, infinity)$. Hence $v^+ = S_(c, j) / N_c$ is
  the unique global maximizer, guaranteeing
  $hat(cal(L))(theta, phi, Psi^+) >= hat(cal(L))(theta, phi, Psi)$
  with strict inequality unless $Psi = Psi^+$.

  Finally, because $(tilde(mu)_(i, j) - mu_(c, j)^+)^2 >= 0$ and
  $tilde(sigma)_(i, j)^2 >= tau^2$ for all $(i, j)$,
  $
    (sigma_(c, j)^2)^+
    >= 1 / N_c sum_(i=1)^N gamma_(i c) tilde(sigma)_(i, j)^2
    >= tau^2.
  $
  Thus even if $gamma_(i c) -> 0$ for all $i != i_0$ (a classical
  singleton-collapse regime where standard GMM EM @dempster1977em
  drives $sigma_(c, j)^2 -> 0$ and log-likelihood to $+infinity$),
  the VaDE ELBO remains strictly bounded above and nonsingular.
]

#subsec("C", "Posterior Collapse Counterexample and Destabilization")

Despite Lemma 1's protection against variance singularities, VaDE
exhibits a distinct neural failure mode inherited from VAEs:
*posterior collapse*, in which an expressive decoder ignores
$bold(z)$ and collapses all cluster responsibilities to the prior.

#theorem(
  title: "Posterior Collapse Counterexample and SNR Threshold",
  origin: "Original Counterexample",
)[
  Consider continuous data $cal(X) subset.eq bb(R)^D$ with zero
  sample mean and covariance $bold(Sigma)_X$, modeled by a linear
  isotropic Gaussian decoder
  $p_theta (bold(x) | bold(z)) = cal(N)(bold(W) bold(z) + bold(b), sigma_x^2 bold(I)_D)$
  with $bold(W) in bb(R)^(D times d)$ and fixed observation variance
  $sigma_x^2 > 0$.
  1. For any $K >= 2$ and any $bold(pi) in Delta^(K-1)$, the
     degenerate configuration $bold(W) = bold(0)$, $bold(b) = bold(0)$,
     $bold(mu)_c = bold(0)$, $bold(sigma)_c^2 = bold(1)$ for all $c$,
     and $tilde(bold(mu))_phi (bold(x)) equiv bold(0)$,
     $tilde(bold(sigma))_phi^2 (bold(x)) equiv bold(1)$ is an exact
     stationary point of $hat(cal(L))$ with constant responsibilities
     $gamma_(i c) = pi_c$ and zero mutual information
     $I_q (bold(x); c) = 0$.
  2. When $sigma_x^2 > lambda_"max" (bold(Sigma)_X)$, this collapsed
     point is a *strict local maximizer* of the empirical marginal
     log-likelihood over $(bold(W), bold(b))$ at identical mixture
     modes, locking VaDE into trivial clustering ($"ARI" = 0$).
     Conversely, whenever $sigma_x^2 < lambda_"max" (bold(Sigma)_X)$,
     the collapsed point is a strict saddle point and gradient
     ascent escapes along the leading eigenvector of $bold(Sigma)_X$.
]

#proof[
  For part (1), at the degenerate configuration we have
  $q_phi (bold(z) | bold(x)_i) = p_Psi (bold(z) | c) = cal(N)(bold(0), bold(I)_d)$
  for every $i, c$. By Proposition 1(1),
  $gamma_(i c)^ast = pi_c$ independently of $bold(x)_i$, so
  $I_q (bold(x); c) = (1 / N) sum_i D_"KL" (bold(gamma)_i || bold(pi)) = 0$.
  Since $D_"KL" (q_phi || p_Psi) = 0$ attains its global minimum
  $0$, all first-order derivatives of the KL penalty with respect to
  $(tilde(bold(mu))_i, tilde(bold(sigma))_i^2, Psi)$ vanish. For the
  reconstruction term
  $R = - (1 / (2 N sigma_x^2)) sum_i (norm(bold(x)_i - bold(b) - bold(W) tilde(bold(mu))_i)_2^2 + "Tr"(bold(W)^top bold(W) "diag"(tilde(bold(sigma))_i^2)))$,
  differentiating at $bold(W) = bold(0), bold(b) = bold(0), tilde(bold(mu))_i = bold(0)$
  gives
  $nabla_(bold(W)) R = (1 / (N sigma_x^2)) sum_i bold(x)_i tilde(bold(mu))_i^top - (1 / (N sigma_x^2)) bold(W) sum_i (tilde(bold(mu))_i tilde(bold(mu))_i^top + "diag"(tilde(bold(sigma))_i^2)) = bold(0)$,
  $nabla_(bold(b)) R = (1 / (N sigma_x^2)) sum_i bold(x)_i = bold(0)$,
  and $nabla_(tilde(bold(mu))_i) R = (1 / sigma_x^2) bold(W)^top bold(x)_i = bold(0)$.
  Hence every first-order KKT condition holds identically.

  For part (2), with identical mixture modes
  $p_Psi (bold(z)) = cal(N)(bold(0), bold(I)_d)$, the exact
  marginal distribution of $bold(x)$ under the linear Gaussian
  decoder is $cal(N)(bold(b), bold(C))$ with
  $bold(C) = bold(W) bold(W)^top + sigma_x^2 bold(I)_D$. Let
  $bold(W) = epsilon bold(u) bold(v)^top$ be a rank-1 perturbation
  with unit vectors $bold(u) in bb(R)^D, bold(v) in bb(R)^d$. By the
  matrix determinant lemma and Sherman–Morrison identity,
  $
    log det bold(C)
    &= D log sigma_x^2 + log(1 + epsilon^2 / sigma_x^2), \
    "Tr"(bold(C)^(-1) bold(Sigma)_X)
    &= 1 / sigma_x^2 "Tr"(bold(Sigma)_X)
       - (epsilon^2 / (sigma_x^4 + sigma_x^2 epsilon^2)) bold(u)^top bold(Sigma)_X bold(u).
  $
  Expanding the average marginal log-likelihood
  $ell(epsilon) = - (1 / 2)(D log(2 pi) + log det bold(C) + "Tr"(bold(C)^(-1) bold(Sigma)_X))$
  to order $epsilon^2$ yields
  $
    ell(epsilon) - ell(0)
    = epsilon^2 / (2 sigma_x^2)
      ((bold(u)^top bold(Sigma)_X bold(u)) / sigma_x^2 - 1)
      + cal(O)(epsilon^4).
  $
  The leading quadratic coefficient is strictly negative for all
  unit directions $bold(u)$ whenever $sigma_x^2 > lambda_"max" (bold(Sigma)_X)$,
  making $bold(W) = bold(0)$ a strict local maximizer of the exact
  marginal likelihood — and hence also of the ELBO, which is bounded
  above by $ell$ and agrees with $ell$ at $bold(W) = bold(0)$.
  Conversely, when $sigma_x^2 < lambda_"max" (bold(Sigma)_X)$, taking
  $bold(u)$ as the principal eigenvector of $bold(Sigma)_X$ yields
  strictly positive curvature, destabilizing collapse.
]

Theorem 2 rigorously explains the empirical necessity of
autoencoder pretraining in Alg. 184: pretraining initializes
$bold(W)$ and $tilde(bold(mu))_phi$ far outside the basin of
attraction of the collapsed stationary point, while ensuring the
effective observation noise variance stays well below the signal
spectrum.

#subsec("D", "Latent Gauge Symmetry and Hyperoctahedral Identifiability")

#proposition(
  title: "Gauge Group Reduction of Diagonal-Covariance VaDE",
  origin: "Original Proof",
)[
  Let $cal(T)_"diag"$ be the group of latent transformations
  $T(bold(z)) = bold(P) (bold(D) bold(z) + bold(b))$, where
  $bold(P)$ is a $d times d$ permutation matrix,
  $bold(D) = "diag"(a_1, dots, a_d)$ with $a_j != 0$, and
  $bold(b) in bb(R)^d$. For every $T in cal(T)_"diag"$ and cluster
  permutation $tau in S_K$, the VaDE marginal likelihood and ELBO
  are strictly invariant under the coupled reparameterization
  $
    (tilde(bold(mu))_phi, tilde(bold(sigma))_phi)
    &|-> (bold(P)(bold(D) tilde(bold(mu))_phi + bold(b)),
          bold(P) |bold(D)| tilde(bold(sigma))_phi), \
    (pi_c, bold(mu)_c, bold(sigma)_c)
    &|-> (pi_(tau(c)), bold(P)(bold(D) bold(mu)_(tau(c)) + bold(b)),
          bold(P) |bold(D)| bold(sigma)_(tau(c))), \
    g_theta &|-> g_theta compose T^(-1).
  $
  Unlike an isotropic Gaussian VAE @kingma2014vae — whose prior is
  invariant under the full orthogonal group $O(d)$ — diagonal VaDE
  breaks continuous rotational symmetry down to the discrete
  hyperoctahedral group $B_d$ of signed coordinate permutations
  whenever component means or variances are non-isotropic.
]

#proof[
  Under the change of variables $bold(z)' = T(bold(z))$, the
  Jacobian determinant is $|det(bold(P) bold(D))| = product_(j=1)^d |a_j|$.
  Both the transformed encoder density $q'_phi (bold(z)' | bold(x))$
  and the transformed mixture component density
  $p'_Psi (bold(z)' | c')$ acquire the identical reciprocal Jacobian
  factor $(product_j |a_j|)^(-1)$ because $bold(P) bold(D)$ maps
  diagonal Gaussians strictly to diagonal Gaussians with variances
  permuted and scaled by $a_j^2$. Therefore the density ratio
  $p'_Psi (bold(z)' | c') / q'_phi (bold(z)' | bold(x))$ is
  pointwise identical to $p_Psi (bold(z) | tau(c')) / q_phi (bold(z) | bold(x))$.
  Likewise, $g'_theta (bold(z)') = g_theta (T^(-1) (T(bold(z)))) = g_theta (bold(z))$
  preserves $p_theta (bold(x) | bold(z))$ identically, and relabeling
  summation index $c |-> tau(c)$ leaves $cal(L)_"ELBO"$ and
  $log p_(theta, Psi) (bold(x))$ invariant.

  By contrast, for a non-monomial orthogonal rotation
  $bold(Q) in O(d) backslash B_d$, the rotated covariance
  $bold(Q) "diag"(bold(sigma)_c^2) bold(Q)^top$ acquires nonzero
  off-diagonal entries whenever $bold(sigma)_c^2$ has distinct
  coordinates; projecting onto the diagonal family strictly
  increases the KL divergence by
  $(1 / 2) sum_j log(("diag"(bold(Q) "diag"(bold(sigma)_c^2) bold(Q)^top))_j / sigma_(c, j)^2) > 0$
  by Hadamard's determinant inequality. Thus only $B_d$ preserves
  the diagonal variational family without loss.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $P$ denote the floating-point multiply-accumulate count of a
single forward-backward pass through both encoder $f_phi$ and
decoder $g_theta$ (for an $M$-layer MLP with hidden width $H$,
$P = cal(O)(D H + M H^2 + H d)$). Evaluating the analytical
GMM-prior KL term and exact responsibilities $gamma_(i c)$ requires
computing $K$ diagonal Mahalanobis terms in $bb(R)^d$ per sample,
costing $cal(O)(K d)$ time and $cal(O)(K d)$ parameter storage.
Across $E$ epochs over $N$ samples with $L$ Monte Carlo draws per
sample (typically $L = 1$ in SGVB @kingma2014vae), total training
time is
$
  T_"VaDE" = cal(O)(E N (L P + K d)) + T_"pre",
$
where $T_"pre" = cal(O)(E_"pre" N P + I_"EM" N K d)$ accounts for
deterministic autoencoder pretraining and initial latent GMM EM.
Crucially, unlike graph-based deep clustering algorithms (e.g.,
JULE, Alg. 186), VaDE has *strictly linear complexity in $N$*, zero
$cal(O)(N^2)$ pairwise affinity storage, and peak memory
$cal(O)(B (D + H + K + d) + |theta| + |phi| + K d)$ for mini-batch
size $B$.

Geometrically, VaDE partitions the latent space $bb(R)^d$ via
*quadratic Mahalanobis decision surfaces*: the log-odds between
components $c_1$ and $c_2$ under Proposition 1(1) are quadratic in
$tilde(bold(mu))_phi (bold(x))$ and linear in
$tilde(bold(sigma))_phi^2 (bold(x))$. Pullback through the nonlinear
encoder $f_phi$ bends these quadratic latent regions into curved,
non-convex ambient clusters in $bb(R)^D$, while the generative
decoder $g_theta$ endows VaDE with true out-of-sample density
scoring and conditional class generation $bold(x) tilde.op p_theta (bold(x) | bold(z))$
with $bold(z) tilde.op cal(N)(bold(mu)_c, "diag"(bold(sigma)_c^2))$.
Because latent covariances are diagonal, VaDE relies on $f_phi$ to
locally disentangle correlated input features along coordinate axes
(Proposition 2).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(E N (L P + K d))$ training; $cal(O)(P + K d)$ per test sample],
    [Space Complexity], [$cal(O)(|theta| + |phi| + K d + B(D + H + K))$ mini-batch RAM],
    [Latent Geometry], [Axis-aligned quadratic ellipsoids in $bb(R)^d$; nonlinear manifold pullback in $bb(R)^D$],
    [Cluster Balance], [Unconstrained simplex prior $bold(pi) in Delta^(K-1)$ adapts naturally to skewed class sizes],
    [Generative Capability], [Full conditional sampling $c -> bold(z) -> bold(x)$ and ELBO-based density scoring],
    [Initialization Sensitivity], [High without pretraining (Theorem 2 collapse); stabilized by AE + GMM warm start],
    [Primary Hyperparameters], [Cluster count $K$, latent dimension $d$, architecture $(f_phi, g_theta)$, pretraining epochs],
  ),
  caption: [Operational and geometric profile of VaDE (Alg. 184).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Single-Cell RNA Sequencing (scRNA-seq):* Zero-inflated dropout
   noise, curved trajectories, and rare subpopulations suit VaDE's
   nonlinear manifold and unconstrained simplex weights $bold(pi)$,
   which preserve rare cell states under NB / ZINB decoders without
   forcing uniform cluster sizes.

2. *Medical Imaging Subtype Discovery:* Radiomic and histology
   cohorts benefit from calibrated posteriors $gamma_(i c)$ (flagging
   ambiguous cases via high entropy) and synthesizable anatomical
   prototypes $g_theta (bold(mu)_c)$ guaranteed by Theorem 1.

3. *Industrial Telemetry & Predictive Maintenance:* Machinery
   regimes exhibit heterogeneous variances; per-cluster
   $bold(sigma)_c^2$ accommodates tight nominal modes alongside broad
   fault transients, while the ELBO scores out-of-distribution faults.

4. *Customer & Portfolio Segmentation:* Skewed transaction cohorts
   violate uniform-balance assumptions. Unlike DEPICT (Alg. 185),
   VaDE's adaptive $pi_c$ and component scales $bold(sigma)_c^2$ fit
   heavy-tailed commercial segments without cardinality distortion.

5. *Astrophysical Spectral Taxonomy:* Million-scale sky surveys
   demand linear-time $cal(O)(N)$ clustering with Gaussian noise
   likelihoods and physically inspectable synthetic cluster templates
   $g_theta (bold(mu)_c)$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

VaDE fails in four regimes: (i) *posterior collapse* (Theorem 2)
when $sigma_x^2$ is too large or autoregressive decoders ignore
$bold(z)$; (ii) *reconstruction distraction* on textured scenes
where pixel MSE prioritizes nuisance background over object class;
(iii) *small-sample overfitting* when $N << |theta| + |phi|$; and
(iv) *fixed $K$*, requiring nonparametric priors @blei2006variational
when cluster cardinality is unknown.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [scRNA-seq & Omics], [Nonlinear manifold + flexible $bold(pi)$ preserves rare states], [Use NB/ZINB decoder likelihood],
    [Tabular / Telemetry], [Heteroscedastic $bold(sigma)_c^2$ + ELBO anomaly score], [Overfits if $N < 500$; use GMM],
    [Medical Phenotyping], [Calibrated $gamma_(i c)$ + inspectable prototypes $g_theta (bold(mu)_c)$], [Monitor collapse via $I_q (bold(x); c)$],
    [Complex Natural Images], [Moderate; pixel MSE emphasizes background over semantics], [Prefer DAC (Alg. 187) or SCAN],
    [Streaming / Dynamic $K$], [Online SGD possible, but fixed $K$ limits new modes], [Use DenStream or DP-GMM],
  ),
  caption: [Domain suitability and failure-mode matrix for VaDE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

VaDE occupies a central unifying position connecting classical
probabilistic mixture models, generative VAEs, and deterministic
autoencoder-based deep clustering methods (DEC @xie2016dec, IDEC
@guo2017idec, DCN @yang2017dcn). The following unified reduction
formalizes these exact parametric limits.

#proposition(
  title: "Exact Parametric Reductions of VaDE",
  origin: "Original Unified Proof",
)[
  Let $hat(cal(L))_"VaDE" (theta, phi, Psi)$ denote the empirical
  VaDE ELBO with isotropic Gaussian observation variance $sigma_x^2$.
  1. *(Standard VAE Limit):* Setting $K = 1$, $pi_1 = 1$,
     $bold(mu)_1 = bold(0)$, and $bold(sigma)_1^2 = bold(1)$
     recovers the exact vanilla Gaussian VAE ELBO of Kingma and
     Welling @kingma2014vae.
  2. *(Classical Diagonal GMM Limit):* Restricting $f_phi$ and
     $g_theta$ to the identity map ($d = D$,
     $tilde(bold(mu))_phi (bold(x)) = bold(x)$) and taking the
     vanishing-noise limit $tilde(bold(sigma))_phi -> bold(0)$
     reduces the exact coordinate ascent of Proposition 1(1) and
     Lemma 1 identically to classical diagonal-covariance GMM EM
     @dempster1977em up to an additive constant independent of $Psi$.
  3. *(Deep $k$-Means / DCN Limit):* Constraining $pi_c = 1 / K$,
     $bold(sigma)_c^2 = beta^(-1) bold(1)$, and deterministic
     encoder variance $tilde(bold(sigma))_phi^2 = epsilon^2 bold(1)$
     with $beta -> infinity$, $epsilon -> 0$, and $beta epsilon^2 -> 0$
     reduces $- 2 hat(cal(L))_"VaDE"$ (up to constants) to the
     hard-assignment Deep Clustering Network (DCN) objective
     @yang2017dcn:
     $
       1 / N sum_(i=1)^N
       (1 / sigma_x^2 norm(bold(x)_i - g_theta (tilde(bold(mu))_i))_2^2
        + beta min_(c in {1, dots, K}) norm(tilde(bold(mu))_i - bold(mu)_c)_2^2).
     $
]

#proof[
  Part (1) follows immediately from Theorem 1(2) upon substituting
  $K = 1$, $gamma_1 = pi_1 = 1$, $mu_(1, j) = 0$, and
  $sigma_(1, j)^2 = 1$, yielding
  $bb(E)_(q_phi) [log p_theta (bold(x) | bold(z))] - (1 / 2) sum_j (tilde(sigma)_j^2 + tilde(mu)_j^2 - 1 - log tilde(sigma)_j^2)$,
  which is the exact Gaussian VAE ELBO @kingma2014vae.

  For part (2), when $tilde(bold(mu))_i = bold(x)_i$ and
  $tilde(bold(sigma))_i -> bold(0)$, the reconstruction term and
  entropy term $(1 / 2) sum_j (1 + log tilde(sigma)_(i, j)^2)$ are
  independent of $Psi$, while
  $bb(E)_(q_phi) [log p_Psi (bold(z) | c)] -> log p_Psi (bold(x)_i | c)$.
  Consequently, Proposition 1(1) yields exact GMM posterior
  responsibilities
  $gamma_(i c) = (pi_c p_Psi (bold(x)_i | c)) / (sum_k pi_k p_Psi (bold(x)_i | k))$,
  and Lemma 1 reduces to the exact diagonal GMM M-step updates
  @dempster1977em @mclachlan2000finite.

  For part (3), substituting $pi_c = 1 / K$,
  $sigma_(c, j)^2 = beta^(-1)$, and $tilde(sigma)_(i, j)^2 = epsilon^2$
  into Proposition 1(1) gives
  $gamma_(i c)^ast = exp(- (beta / 2) norm(tilde(bold(mu))_i - bold(mu)_c)_2^2) / (sum_k exp(- (beta / 2) norm(tilde(bold(mu))_i - bold(mu)_k)_2^2))$.
  As $beta -> infinity$, this softmax concentrates on the nearest
  centroid $c^ast (i) = arg min_c norm(tilde(bold(mu))_i - bold(mu)_c)_2^2$,
  driving the responsibility entropy $- sum_c gamma_(i c)^ast log gamma_(i c)^ast -> 0$.
  Simultaneously, $beta epsilon^2 -> 0$ eliminates the variance-floor
  term in Theorem 1(2) and makes the reconstruction expectation
  converge to the deterministic squared error
  $- (1 / (2 sigma_x^2)) norm(bold(x)_i - g_theta (tilde(bold(mu))_i))_2^2$.
  Multiplying $- hat(cal(L))_"VaDE"$ by $2$ and dropping terms
  independent of $(theta, phi, {bold(mu)_c})$ yields the exact DCN
  objective @yang2017dcn.
]

Empirical evaluation on the ClusterOverview four-dataset benchmark
(quantified via Silhouette width @rousseeuw1987silhouettes and
Adjusted Rand Index @hubert1985comparing) illuminates every aspect
of this theory. On *Multi-Shape* (non-convex geometric manifolds),
VaDE achieves Silhouette $0.718$ (ARI N/A), tying JULE (Alg. 186)
for the highest manifold separation among deep architectures and
demonstrating how nonlinear decoding pulls curved ambient shapes
into clean Gaussian latent clusters. On *Customers* (skewed,
unbalanced tabular spend features), VaDE achieves Silhouette
$0.381$ — outperforming DEPICT ($0.327$), JULE ($0.335$), and DAC
($0.197$) — directly confirming the advantage of its unconstrained
simplex weights $bold(pi)$ and heteroscedastic component variances
$bold(sigma)_c^2$. On *Iris*, VaDE records Silhouette $0.350$ and
ARI $0.579$: the KL contraction toward overlapping Versicolor and
Virginica prior modes slightly reduces Euclidean margin compared
with discriminative sharpeners (DEPICT $0.447$, DAC $0.433$) while
maintaining comparable label recovery. On *High-D*, VaDE attains
Silhouette $0.245$ (ARI N/A), avoiding the severe $K_s$-NN
degradation exhibited by JULE ($0.119$) because latent compression
filters ambient noise. Finally, its mean runtime of $11.007$ ms
reflects the modest linear-time overhead of joint encoder–decoder
SGVB passes relative to decoder-free DAC ($3.329$ ms), while
remaining substantially faster than graph-agglomerative JULE
($14.593$ ms).
