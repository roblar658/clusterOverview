#import "../template/ieee.typ": *

#alg_entry(
  id: 190,
  name: "SwAV (Swapping Assignments between Views)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2020",
  paradigm: "Online Optimal-Transport Swapped Prototype Clustering",
  complexity: $O(E N (C_"net" + M (r K + L_"sk" K)))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Instance contrastive methods such as InfoNCE @oord2018infonce require massive negative-pair banks, whereas offline DeepCluster @caron2018deepcluster requires full-dataset passes and global K-Means sweeps. *Swapping Assignments between Views (SwAV)* @caron2020swav unifies contrastive view invariance with online deep clustering by comparing views indirectly via $K$ learnable unit prototypes while enforcing equipartitioned minibatch assignments through entropic optimal transport @cuturi2013sinkhorn.

Let $cal(X) = {bold(x)_1, ..., bold(x)_N} subset.eq bb(R)^d$ and view transforms $cal(T)$. An encoder $f_theta : bb(R)^d -> bb(S)^(r-1)$ maps each augmented view $bold(x)_i^((v)) = g_v (bold(x)_i)$ ($v = 1, ..., M$) to $ell_2$-normalized features $bold(z)_i^((v)) = f_theta (bold(x)_i^((v)))$. Let $bold(C) = [bold(c)_1, ..., bold(c)_K] in bb(R)^(r times K)$ hold unit prototypes ($||bold(c)_k||_2 = 1$). At temperature $tau > 0$, view $v$ yields softmax predictions
$
  p_k (bold(z)_i^((v))) = (exp(chevron.l bold(z)_i^((v)), bold(c)_k chevron.r / tau)) / (sum_(k'=1)^K exp(chevron.l bold(z)_i^((v)), bold(c)_(k') chevron.r / tau)).
$
To compute target codes *online* without collapse, SwAV stacks minibatch features $bold(Z)^((v)) in bb(R)^(r times B)$ (optionally padded with a FIFO feature queue to length $B_Q >= B$) and projects onto the uniform transportation polytope.

#definition(title: "Entropic Equipartition Polytope and SwAV Objective")[
  Let $cal(U)(B, K) = {bold(Q) in bb(R)_(>= 0)^(K times B) : bold(Q) bold(1)_B = (1 / K) bold(1)_K, bold(Q)^top bold(1)_K = (1 / B) bold(1)_B}$ and $H(bold(Q)) = - sum_(k, i) Q_(k i) ln Q_(k i)$. For smoothing $epsilon > 0$ and score matrix $bold(S)^((v)) = bold(C)^top bold(Z)^((v)) in [-1, 1]^(K times B)$, the *online Sinkhorn code matrix* is
  $
    bold(Q)^((v)*) = arg max_(bold(Q) in cal(U)(B, K)) ( "Tr"(bold(Q)^top bold(S)^((v))) + epsilon H(bold(Q)) ),
  $
  with column-normalized codes $tilde(bold(q))_i^((v)) = B bold(q)_i^((v)*) in Delta^(K-1)$. For global views $(1, 2)$, the *SwAV swapped loss* is
  $
    cal(L)_"SwAV" (theta, bold(C)) &= 1 / B sum_(i=1)^B ( ell(bold(z)_i^((1)), tilde(bold(q))_i^((2))) \
      &quad + ell(bold(z)_i^((2)), tilde(bold(q))_i^((1))) ),
  $
  where $ell(bold(z), tilde(bold(q))) = - sum_(k=1)^K tilde(q)_k ln p_k (bold(z))$, with stop-gradient on $tilde(bold(q))_i^((v))$.
]

In *multi-crop* SwAV @caron2020swav, only $2$ full-resolution global views compute target codes while $M - 2$ cheap local crops predict them, multiplying invariance constraints at low compute cost.

#algbox(
  title: "SwAV: Online Swapped Prototype Clustering via Sinkhorn-Knopp",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$; prototypes $K$; temperature $tau$; smoothing $epsilon > 0$; steps $L_"sk"$; views $M >= 2$; batch $B$; queue $B_Q$; step $eta$; epochs $E$.],
  outputs: [Encoder $f_theta$, unit prototypes $bold(C)$, hard labels $hat(y)_i = arg max_k chevron.l f_theta (bold(x)_i), bold(c)_k chevron.r$.],
)[
  + Initialize $f_theta$ and unit prototypes $bold(C) in bb(R)^(r times K)$ ($||bold(c)_k||_2 = 1$).
  + *For* epoch $e = 1, ..., E$ *do*:
    + *For* each minibatch $cal(B)$ of size $B$ *do*:
      + Compute $ell_2$-normalized embeddings $bold(z)_i^((v)) = f_theta (bold(x)_i^((v)))$ for $v = 1, ..., M$.
      + *For* global views $v in {1, 2}$ *do* (stop-gradient):
        + Form $bold(S)^((v)) = bold(C)^top [bold(Z)^((v)), bold(Z)_"queue"]$, $bold(K) = exp(bold(S)^((v)) / epsilon)$.
        + Set $bold(u) <- bold(1)_K$, $bold(v) <- bold(1)_(B_Q)$; *for* $l = 1, ..., L_"sk"$: $bold(u) <- (1 / K) (bold(K) bold(v))^(-1)$, $bold(v) <- (1 / B_Q) (bold(K)^top bold(u))^(-1)$.
        + Extract codes $tilde(bold(q))_i^((v)) <- B_Q u_k K_(k i) v_i$ for $i in cal(B)$.
      + Descend $(1 / (2 B)) sum_(i in cal(B)) sum_(t in {1, 2}) sum_(v != t) ell(bold(z)_i^((v)), tilde(bold(q))_i^((t)))$ on $(theta, bold(C))$ and re-normalize $bold(c)_k <- bold(c)_k / ||bold(c)_k||_2$.
  + *Return* $(f_theta, bold(C))$ and $hat(y)_i = arg max_k chevron.l f_theta (bold(x)_i), bold(c)_k chevron.r$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Unique Sinkhorn Optimum and Birkhoff Linear Convergence")

SwAV's inner assignment is entropic optimal transport @cuturi2013sinkhorn, guaranteeing a unique interior scaling optimum and geometric Hilbert-metric contraction.

#theorem(title: "Unique Diagonal-Scaling Optimum and Hilbert Contraction Rate", origin: "Classical — Cited")[
  For any $bold(S) = bold(C)^top bold(Z) in [-1, 1]^(K times B)$ and $epsilon > 0$ @cuturi2013sinkhorn @caron2020swav:
  + $J(bold(Q)) = "Tr"(bold(Q)^top bold(S)) + epsilon H(bold(Q))$ is strictly concave on $cal(U)(B, K)$ with unique strictly positive maximizer
    $
      bold(Q)^* = "diag"(bold(u)^*) bold(K) "diag"(bold(v)^*),
      quad bold(K) = exp(bold(S) / epsilon),
    $
    for $bold(u)^* in bb(R)_(>0)^K, bold(v)^* in bb(R)_(>0)^B$ unique up to $(c bold(u)^*, c^(-1) bold(v)^*)$.
  + Alternating Sinkhorn--Knopp iterations contract the Hilbert projective metric $d_H (bold(u), bold(u)') = ln (max_k (u_k / u'_k) dot max_k (u'_k / u_k))$ linearly:
    $
      d_H (bold(u)^((l)), bold(u)^*) &<= kappa(bold(K))^(2 l) d_H (bold(u)^((0)), bold(u)^*), \
      kappa(bold(K)) &<= tanh(1 / epsilon) < 1.
    $
]

#proof[
  *Part 1.* Because $cal(U)(B, K)$ is nonempty, compact, and convex and $-H(bold(Q))$ is strictly convex with boundary gradient diverging to $-infinity$, the unique maximizer lies strictly in the interior ($Q^*_(k i) > 0$). Forming the Lagrangian with multipliers $bold(alpha) in bb(R)^K, bold(beta) in bb(R)^B$,
  $
    Lambda(bold(Q), bold(alpha), bold(beta)) &= sum_(k, i) (Q_(k i) S_(k i) - epsilon Q_(k i) ln Q_(k i)) \
      &quad + bold(alpha)^top (bold(Q) bold(1)_B - 1 / K bold(1)_K) \
      &quad + bold(beta)^top (bold(Q)^top bold(1)_K - 1 / B bold(1)_B),
  $
  and setting $partial Lambda / (partial Q_(k i)) = 0$ gives $S_(k i) - epsilon (ln Q_(k i) + 1) + alpha_k + beta_i = 0$, hence
  $
    Q^*_(k i) &= exp(alpha_k / epsilon - 1 / 2) exp(S_(k i) / epsilon) exp(beta_i / epsilon - 1 / 2) \
      &= u^*_k K_(k i) v^*_i.
  $

  *Part 2.* Elementwise inversion is an exact isometry of $d_H$, since $d_H (bold(w)^(-1), (bold(w)')^(-1)) = d_H (bold(w), bold(w)')$. By the Birkhoff--Hopf theorem, the positive map $bold(v) |-> bold(K) bold(v)$ contracts $d_H$ by $kappa(bold(K)) = tanh(Delta(bold(K)) / 4)$, with projective diameter
  $
    Delta(bold(K)) &= max_(k, k', i, i') ln (K_(k i) K_(k' i')) / (K_(k i') K_(k' i)) \
      &= 1 / epsilon max_(k, k', i, i') (S_(k i) + S_(k' i') - S_(k i') - S_(k' i)).
  $
  Because $S_(k i) in [-1, 1]$, we have $Delta(bold(K)) <= 4 / epsilon$ and $kappa(bold(K)) <= tanh(1 / epsilon) < 1$. Each full iteration applies both $bold(K)$ and $bold(K)^top$ (sharing diameter $Delta(bold(K))$), contracting $d_H$ by $kappa(bold(K))^2$.
]

At practical settings ($epsilon = 0.05, L_"sk" = 3$ @caron2020swav), $L_"sk" = 3$ acts as an approximate projection whose final column scaling guarantees $tilde(bold(q))_i in Delta^(K-1)$ while softly balancing prototypes.

#subsec("B", "Exact Variational Bregman Decomposition and Swapped Discrepancy")

We prove an exact variational KL identity showing that SwAV's stop-gradient update corresponds to alternating Bregman / EM descent augmented by a nonnegative cross-view invariance discrepancy.

#theorem(title: "Exact Variational KL Decomposition and Nonnegative Swapped Discrepancy", origin: "Original Proof")[
  Suppose $epsilon = tau > 0$, and define the view-$v$ self-prediction energy
  $
    cal(F)_"self" (bold(Z)^((v)), bold(Q)^((v))) := tau / B sum_(i=1)^B D_"KL" (B bold(q)_i^((v)) || bold(p)(bold(z)_i^((v)))).
  $
  + Minimizing $cal(F)_"self"$ over $bold(Q)^((v)) in cal(U)(B, K)$ is identical to SwAV's Sinkhorn code step, while minimizing over $(theta, bold(C))$ at fixed $bold(Q)^((v))$ (stop-gradient) is identical to cross-entropy minimization.
  + Define $cal(F)_"swap" := (tau / B) sum_(i=1)^B [D_"KL" (tilde(bold(q))_i^((2)) || bold(p)_i^((1))) + D_"KL" (tilde(bold(q))_i^((1)) || bold(p)_i^((2)))]$. At Sinkhorn-optimal codes $bold(Q)^((1)*), bold(Q)^((2)*)$,
    $
      cal(F)_"swap" &= cal(F)_"self"^((1)) + cal(F)_"self"^((2)) \
        &quad + "Tr"( (bold(Q)^((1)*) - bold(Q)^((2)*))^top (bold(S)^((1)) - bold(S)^((2))) ),
    $
    where the cross-view discrepancy is strictly nonnegative and vanishes iff $bold(Q)^((1)*) = bold(Q)^((2)*)$.
]

#proof[
  *Part 1.* Setting $tilde(q)_(k i) = B Q_(k i)$ and $Z_i = sum_(k'=1)^K exp(S_(k' i) / tau)$, expand
  $
    cal(F)_"self" &= tau / B sum_(i=1)^B sum_(k=1)^K tilde(q)_(k i) (ln tilde(q)_(k i) - S_(k i) / tau + ln Z_i) \
      &= - "Tr"(bold(Q)^top bold(S)) - tau H(bold(Q)) + tau ln B + tau / B sum_(i=1)^B ln Z_i,
  $
  using $(1 / B) sum_(i, k) B Q_(k i) ln(B Q_(k i)) = -H(bold(Q)) + ln B$. Because $tau ln B + (tau / B) sum_i ln Z_i$ is independent of $bold(Q)$, minimizing $cal(F)_"self"$ over $cal(U)(B, K)$ maximizes $"Tr"(bold(Q)^top bold(S)) + tau H(bold(Q))$. Holding $bold(Q)$ fixed eliminates $-tau H(bold(Q))$ from parameter gradients, leaving $tau$ times cross-entropy.

  *Part 2.* Expanding $cal(F)_"swap"$ cancels all entropy terms $-tau H(bold(Q)^((v)*))$ and log-partition terms $(tau / B) sum_i ln Z_i^((v))$ against $cal(F)_"self"^((1)) + cal(F)_"self"^((2))$, leaving
  $
    cal(F)_"swap" - cal(F)_"self"^((1)) - cal(F)_"self"^((2))
    &= "Tr"((bold(Q)^((1)*))^top bold(S)^((1))) + "Tr"((bold(Q)^((2)*))^top bold(S)^((2))) \
    &quad - "Tr"((bold(Q)^((2)*))^top bold(S)^((1))) - "Tr"((bold(Q)^((1)*))^top bold(S)^((2))) \
    &= "Tr"((bold(Q)^((1)*) - bold(Q)^((2)*))^top (bold(S)^((1)) - bold(S)^((2)))).
  $
  By optimality of $bold(Q)^((1)*)$ over $cal(U)(B, K)$ against feasible $bold(Q)^((2)*)$,
  $
    "Tr"((bold(Q)^((1)*))^top bold(S)^((1))) + tau H(bold(Q)^((1)*)) \
    >= "Tr"((bold(Q)^((2)*))^top bold(S)^((1))) + tau H(bold(Q)^((2)*)),
  $
  strictly whenever $bold(Q)^((1)*) != bold(Q)^((2)*)$. Summing with the symmetric optimality bound for $bold(Q)^((2)*)$ cancels $tau H$ and proves $"Tr"((bold(Q)^((1)*) - bold(Q)^((2)*))^top (bold(S)^((1)) - bold(S)^((2)))) >= 0$, with equality iff $bold(Q)^((1)*) = bold(Q)^((2)*)$.
]

#subsec("C", "Rank-Collapse Immunity and Equiangular Simplex Prototypes")

The transportation polytope strictly blocks constant-representation collapse regardless of prototype geometry.

#proposition(title: "Exact Rank-1 Collapse Immunity and Simplex Prototype Margin", origin: "Original Proof")[
  + *(Collapse Barrier).* If $bold(z)_i^((v)) equiv bold(z)_0 in bb(S)^(r-1)$, then for every $bold(C)$ and $epsilon > 0$, the unique Sinkhorn optimum is $bold(Q)^* = (1 / (K B)) bold(1)_K bold(1)_B^top$, forcing $cal(L)_"SwAV" >= 2 ln K$.
  + *(Equiangular Simplex Margin).* For $K <= r + 1$, placing prototypes at regular $(K - 1)$-simplex vertices ($chevron.l bold(c)_k, bold(c)_(k') chevron.r = -1 / (K - 1)$) with aligned invariant classes achieves
    $
      cal(L)_"SwAV" ("simplex") = 2 ln ( 1 + (K - 1) exp(- K / ((K - 1) tau)) ) < 2 ln K,
    $
    vanishing exponentially as $tau -> 0^+$.
]

#proof[
  *Part 1.* Under $bold(z)_i^((v)) equiv bold(z)_0$, $bold(S) = bold(a) bold(1)_B^top$ ($bold(a) = bold(C)^top bold(z)_0$). For every $bold(Q) in cal(U)(B, K)$, row equipartition yields
  $
    "Tr"(bold(Q)^top bold(S)) = bold(a)^top bold(Q) bold(1)_B = 1 / K bold(a)^top bold(1)_K,
  $
  constant over $cal(U)(B, K)$. Maximizing $"Tr"(bold(Q)^top bold(S)) + epsilon H(bold(Q))$ reduces to maximizing $H(bold(Q))$, uniquely attained at $(1 / (K B)) bold(1)_K bold(1)_B^top$ ($tilde(bold(q))_i = (1 / K) bold(1)_K$). Thus per view,
  $
    ell(bold(z)_0, 1 / K bold(1)_K) = ln K + D_"KL" (1 / K bold(1)_K || bold(p)(bold(z)_0)) >= ln K,
  $
  giving $cal(L)_"SwAV" >= 2 ln K$.

  *Part 2.* For a regular unit simplex ($sum_k bold(c)_k = bold(0)$), expanding $||sum_k bold(c)_k||_2^2 = K + K(K - 1) rho = 0$ gives $rho = -1 / (K - 1)$. Aligned features yield
  $
    p_c (bold(c)_c) &= (1 + (K - 1) exp(- (1 - rho) / tau))^(-1) \
      &= (1 + (K - 1) exp(- K / ((K - 1) tau)))^(-1).
  $
  Taking $-2 ln p_c (bold(c)_c)$ gives the stated bound, strictly below $2 ln K$ for every $tau > 0$.
]

#subsec("D", "Equipartition Distortion and Queue Necessity under Small Batches")

Enforcing exact uniform row marginals on every minibatch introduces a tight Total-Variation error floor whenever empirical batch composition deviates from uniformity.

#theorem(title: "Total-Variation Error Floor and Finite-Batch Queue Bound", origin: "Original Counterexample")[
  Let minibatch $cal(B)$ of size $B$ have ground-truth class proportions $hat(bold(pi))_B$ and $bold(u)_K = (1 / K) bold(1)_K$.
  + *(Imbalance Floor).* For every feasible $tilde(bold(Q)) = B bold(Q)$ ($bold(Q) in cal(U)(B, K)$) and every permutation $sigma in S_K$,
    $
      min_(sigma in S_K) 1 / B sum_(i=1)^B (1 - tilde(q)_(sigma(y_i), i))
      >= d_"TV" (hat(bold(pi))_B, bold(u)_K) := 1 / 2 ||hat(bold(pi))_B - bold(u)_K||_1.
    $
  + *(Balanced Finite-Batch Floor).* Even when population classes are perfectly balanced ($bold(pi) = bold(u)_K$), an i.i.d. batch of size $B$ without a queue obeys
    $
      bb(E)[d_"TV" (hat(bold(pi))_B, bold(u)_K)] >= max( 1 - B / K, (1 - 1 / K)^B ).
    $
    At $B = K$, at least $(1 - 1 / K)^K approx 36.8\%$ of codes are forced erroneous; a FIFO queue of size $B_Q = m K$ shrinks this empty-class floor exponentially as $e^(-m)$.
]

#proof[
  *Part 1.* For any $sigma in S_K$, let $M_k = (1 / B) sum_(i : sigma(y_i) = k) tilde(q)_(k, i)$. Both $tilde(q)_(k, i) <= 1$ and $sum_i tilde(q)_(k, i) = B / K$ give $M_k <= min(hat(pi)_(B, sigma^(-1)(k)), 1 / K)$. Summing over $k$ via $min(a, b) = (a + b - |a - b|) / 2$ yields
  $
    1 / B sum_(i=1)^B tilde(q)_(sigma(y_i), i) &= sum_(k=1)^K M_k <= sum_(k=1)^K min(hat(pi)_(B, k), 1 / K) \
      &= 1 - d_"TV" (hat(bold(pi))_B, bold(u)_K).
  $
  Subtracting from $1$ proves the tight Total-Variation error floor.

  *Part 2.* Using $d_"TV" (hat(bold(pi))_B, bold(u)_K) = sum_(k=1)^K max(1 / K - hat(pi)_(B, k), 0)$, whenever $B < K$ at least $K - B$ classes are empty, giving $d_"TV" >= 1 - B / K$. For any $B >= 1$ under $B hat(bold(pi))_B tilde.op "Multinomial"(B, bold(u)_K)$, $max(1 / K - hat(pi)_(B, k), 0) >= (1 / K) bb(1)_(hat(pi)_(B, k) = 0)$, so taking expectations yields
  $
    bb(E)[d_"TV" (hat(bold(pi))_B, bold(u)_K)] >= 1 / K sum_(k=1)^K "Pr"(hat(pi)_(B, k) = 0) = (1 - 1 / K)^B.
  $
  Setting $B_Q = m K$ gives $(1 - 1 / K)^(m K) <= e^(-m)$, explaining mathematically why SwAV requires a feature queue whenever $B$ is not much larger than $K$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Per minibatch of size $B$ with $M$ views, network passes cost $O(B M C_"net")$, prototype dot products $bold(C)^top bold(Z)^((v))$ cost $O(B M r K)$, and $L_"sk"$ Sinkhorn iterations on $K times B_Q$ over $2$ global views cost $O(L_"sk" B_Q K)$. Total training time is $O(E N (C_"net" + M (r K + L_"sk" K)))$ (for $B_Q = O(B)$), and inference costs $O(N (C_"net" + r K))$. Space scales as $O(|theta| + r K + B_Q (r + K))$, decoupling memory from dataset size $N$ and scaling linearly in $K$. Geometrically, embeddings and prototypes lie on $bb(S)^(r-1)$, making SwAV an online, equipartitioned deep relative of spherical K-Means @dhillon2001concept and von Mises--Fisher mixtures @banerjee2005vmf.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(E N (C_"net" + M (r K + L_"sk" K)))$ ($L_"sk" approx 3$); $O(N (C_"net" + r K))$ inference.],
    [Space Complexity], [$O(|theta| + r K + B_Q (r + K))$; no full-dataset cache or $k$-NN index.],
    [Inductive Bias], [Hyperspherical conical prototypes on $bb(S)^(r-1)$ with equipartitioned codes.],
    [Collapse Resistance], [Strict online barrier $cal(L)_"SwAV" >= 2 ln K$ (Proposition II.C).],
    [Batch / Queue Regime], [Needs $B_Q >> K$ (Theorem II.D) to suppress empty-prototype floor $(1 - 1 / K)^(B_Q)$.],
    [Multi-Crop Efficiency], [Amortizes global-view Sinkhorn codes across $M - 2$ cheap local crops.],
    [Hyperparameters], [Prototypes $K$, temperature $tau$, Sinkhorn $epsilon$, steps $L_"sk"$, queue $B_Q$, views $M$.],
    [Streaming / Online], [Fully online minibatch algorithm; scales to billion-scale web streams.],
  ),
  caption: [Operational and geometric profile of SwAV (Alg. 190).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Web-Scale Self-Supervised Vision Pretraining:* On ImageNet-1K and billion-image corpora (IG-1B), SwAV scales to $K = 3000$ prototypes online without huge negative banks or offline K-Means sweeps, matching supervised transfer accuracy via multi-crop views @caron2020swav.
+ *Continuous Satellite and Streaming Video Indexing:* Online prototype updates and FIFO feature queues maintain balanced semantic hash codes in $O(B_Q K)$ memory when dataset size $N$ is unbounded.
+ *Bioacoustic and Speech Representation Learning:* Pairing temporal crops and spectral masks discovers discrete phonetic or species call prototypes on $bb(S)^(r-1)$ while blocking collapse onto silence frames.
+ *Cryo-EM and Molecular Conformation Discovery:* Noisy 2D macromolecular projections or subgraph perturbations align angularly onto equiangular structural prototypes without $O(B^2)$ pairwise contrastive matrices.
+ *Distributed and Federated Representation Learning:* Workers synchronize compact $r times K$ prototype gradients and batch statistics rather than shipping raw negative samples across nodes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Three failure regimes arise: (1) *strong class imbalance* incurs the tight $d_"TV" (bold(pi), bold(u)_K)$ error floor (Theorem II.D) unless prototypes are over-parameterized ($K >> K_"true"$); (2) *small batches without a queue* ($B <= K$) suffer severe empty-prototype sampling noise; and (3) *very small $epsilon$ with $L_"sk" = 3$* risks exponential overflow and incomplete marginal balancing.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Large Balanced Vision / Streams], [Excellent: online, linear in $K$, multi-crop boosts invariance.], [Keep queue $B_Q >= 5 K$ for small batches.],
    [High-D Angular Data], [Strong: $bb(S)^(r-1)$ normalization resists norm drift.], [Tune $tau$ and $epsilon$ jointly (Theorem II.B).],
    [Skewed Direct $K$-Clustering], [Weak: exact equipartition forces $d_"TV"$ floor (Theorem II.D).], [Overcluster ($K >> K_"true"$) or use SCAN.],
    [Small Tabular Sets ($N < 500$)], [Moderate: small $B$ amplifies batch-level equipartition noise.], [Use full-batch Sinkhorn or SCAN (Alg. 189).],
    [Static Single-View Data], [Limited: swapping needs a meaningful view generator $cal(T)$.], [Use DeepCluster (Alg. 191) or DEC.],
  ),
  caption: [Domain suitability and failure-mode matrix for SwAV (Alg. 190).],
)

#sec("V", "Structural Relations and Empirical Behavior")

SwAV bridges offline DeepCluster @caron2018deepcluster, spherical K-Means @dhillon2001concept @banerjee2005vmf, and InfoNCE @oord2018infonce via exact zero-entropy reduction.

#proposition(title: "Zero-Entropy Reduction to Balanced Spherical K-Means and DeepCluster", origin: "Original Derivation")[
  For full-batch ($B = N$), single-view ($bold(Z)^((1)) = bold(Z)^((2)) = bold(Z)$) SwAV with fixed $f_theta$ and unit prototypes $bold(C)$:
  + As $epsilon -> 0^+$, Sinkhorn optimal transport converges to cardinality-balanced spherical K-Means @dhillon2001concept @banerjee2005vmf:
    $
      lim_(epsilon -> 0^+) bold(Q)^*(epsilon) in arg min_(bold(Q) in cal(U)(N, K)) 1 / 2 sum_(i=1)^N sum_(k=1)^K Q_(k i) ||bold(z)_i - bold(c)_k||_2^2.
    $
  + Dropping row equipartition $bold(Q) bold(1)_N = (1 / K) bold(1)_K$ as $epsilon -> 0^+$ recovers DeepCluster's offline K-Means assignment step @caron2018deepcluster (Alg. 191), while persistent prototypes $bold(C)$ remove epoch-wise head resets.
]

#proof[
  Since $||bold(z)_i||_2 = ||bold(c)_k||_2 = 1$, $(1 / 2)||bold(z)_i - bold(c)_k||_2^2 = 1 - S_(k i)$, so
  $
    1 / 2 sum_(i=1)^N sum_(k=1)^K Q_(k i) ||bold(z)_i - bold(c)_k||_2^2 = 1 - "Tr"(bold(Q)^top bold(S)).
  $
  Maximizing $"Tr"(bold(Q)^top bold(S)) + epsilon H(bold(Q))$ over $cal(U)(N, K)$ minimizes spherical K-Means distortion minus $epsilon H(bold(Q))$; since $0 <= H(bold(Q)) <= ln(N K)$, limit points as $epsilon -> 0^+$ solve exact equal-size spherical K-Means @cuturi2013sinkhorn. Removing the row constraint decouples samples into unconstrained Voronoi assignments $arg max_k chevron.l bold(z)_i, bold(c)_k chevron.r$, recovering DeepCluster's assignment step @caron2018deepcluster.
]

*Empirical Behavior on the ClusterOverview Benchmark.* On the 4-dataset ClusterOverview benchmark @rousseeuw1987silhouettes @hubert1985comparing, SwAV attains *Iris* $"Sil" = 0.334$, $"ARI" = 0.372$; *Customers* $"Sil" = 0.284$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.292$ ($"ARI" =$ N/A), with mean runtime *6.479 ms*. SwAV achieves strong *High-D* ($0.292$, matching IIC's $0.298$ and nearly tripling DeepCluster's $0.109$) and *Multi-Shape* ($0.712$) Silhouette because $bb(S)^(r-1)$ normalization and entropic optimal transport stabilize angular clusters against Euclidean concentration @beyer1999nearest. Its *Customers* Silhouette ($0.284$) exceeds both IIC ($0.257$) and DeepCluster ($0.259$), whereas lower *Iris* ARI ($0.372$) and Silhouette ($0.334$) reflect Theorem II.D's finite-batch equipartition variance on small datasets. Its $6.479$ ms runtime reflects the lightweight $L_"sk"$-step Sinkhorn loop.
