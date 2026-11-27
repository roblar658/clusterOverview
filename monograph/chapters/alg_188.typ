#import "../template/ieee.typ": *

#alg_entry(
  id: 188,
  name: "IIC (Invariant Information Clustering)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2019",
  paradigm: "Paired-View Mutual Information Maximization",
  complexity: $O(E N (C_"net" + K^2))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical deep clustering pipelines such as DeepCluster @caron2018deepcluster and DEC @xie2016dec alternate external pseudo-labeling with discriminative training, exposing representations to empty-cluster collapse unless auxiliary heuristics are imposed. *Invariant Information Clustering (IIC)* @ji2019iic removes external clustering stages by directly training a neural softmax head to maximize discrete mutual information between soft cluster assignments of paired semantically equivalent views.

Let $cal(X) = {bold(x)_1, ..., bold(x)_N} subset.eq bb(R)^d$ be drawn i.i.d. from $P_X$, and let $cal(T)$ be a label-preserving stochastic transformation family with law $P_cal(T)$ (geometric warps, photometric jitter, or spatial offsets). Given $bold(x) tilde.op P_X$, IIC draws paired views $bold(x)^((1)) = g_1 (bold(x))$ and $bold(x)^((2)) = g_2 (bold(x))$ with $g_1, g_2 tilde.op P_cal(T)$ (or pairs $bold(x)$ with $bold(x)' = g(bold(x))$). A neural network $Phi_theta : bb(R)^d -> Delta^(K-1)$ maps each view onto the probability simplex
$
  Delta^(K-1) = {bold(p) in bb(R)^K_(>= 0) : bold(1)_K^top bold(p) = 1},
$
yielding $bold(p)_i = Phi_theta (bold(x)_i^((1)))$ and $bold(p)'_i = Phi_theta (bold(x)_i^((2)))$ as conditional cluster distributions $P(Z = c | bold(x)_i^((1)))$ and $P(Z' = c' | bold(x)_i^((2)))$.

#definition(title: "Symmetrized Joint Matrix and Generalized IIC Objective")[
  For a minibatch $cal(B)$ of size $B$, the raw empirical joint matrix is $bold(P)_"raw" (theta) = (1 / B) sum_(i in cal(B)) bold(p)_i (bold(p)'_i)^top in Delta^(K times K)$, symmetrized as
  $
    bold(P)(theta) = 1 / 2 (bold(P)_"raw" (theta) + bold(P)_"raw" (theta)^top),
  $
  with identical marginals $P_c = sum_(c'=1)^K P_(c c')$. For marginal-entropy weight $lambda >= 1$, the *generalized IIC objective* is
  $
    I_lambda (bold(P)) &= sum_(c=1)^K sum_(c'=1)^K P_(c c') ln (P_(c c') / (P_c^lambda P_(c')^lambda)) \
                       &= I(Z; Z') + 2(lambda - 1) H(Z),
  $
  where $H(Z) = - sum_(c=1)^K P_c ln P_c$ and $I(Z; Z') = 2 H(Z) - H(bold(P))$.
]

Setting $lambda = 1$ recovers exact mutual information $I(Z; Z')$; setting $lambda > 1$ boosts marginal entropy @ji2019iic. In practice, IIC attaches an auxiliary *overclustering head* $tilde(Phi)_phi : bb(R)^d -> Delta^(K'-1)$ ($K' > K$) sharing backbone $f_omega$ so that fine-grained nuisance modes enrich shared features.

#algbox(
  title: "Invariant Information Clustering (IIC) with Overclustering",
  inputs: [Dataset $cal(X) subset.eq bb(R)^d$; counts $K < K'$; law $P_cal(T)$; batch $B$; weight $lambda >= 1$; step $eta$; epochs $E$; floor $epsilon > 0$.],
  outputs: [Encoder $f_omega$, head $Phi_theta$, hard labels $hat(y)_i = arg max_c [Phi_theta (bold(x)_i)]_c$.],
)[
  + Initialize backbone $f_omega$, primary head $Phi_theta$ ($K$), and auxiliary head $tilde(Phi)_phi$ ($K'$).
  + *For* epoch $e = 1, ..., E$ *do*:
    + *For* each minibatch $cal(B)$ of size $B >> K^2$ *do*:
      + Sample $bold(x)_i^((1)) = g_(i, 1) (bold(x)_i)$, $bold(x)_i^((2)) = g_(i, 2) (bold(x)_i)$ ($g_(i, 1), g_(i, 2) tilde.op P_cal(T)$).
      + Compute $bold(p)_i = Phi_theta (bold(x)_i^((1)))$ and $bold(p)'_i = Phi_theta (bold(x)_i^((2)))$.
      + Form $bold(P)_"raw" = (1 / B) sum_(i in cal(B)) bold(p)_i (bold(p)'_i)^top$, set $bold(P) <- (bold(P)_"raw" + bold(P)_"raw"^top) / 2$, and clip entries below $epsilon$.
      + Compute $P_c <- sum_(c'=1)^K P_(c c')$ and descend $-I_lambda (bold(P))$ on $(omega, theta)$.
      + Periodically form $tilde(bold(P)) in Delta^(K' times K')$ via $tilde(Phi)_phi$ and ascend $I_lambda (tilde(bold(P)))$ on $(omega, phi)$.
  + *Return* $(f_omega, Phi_theta)$ and $hat(y)_i = arg max_c [Phi_theta (bold(x)_i)]_c$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Positive Semidefiniteness, Collapse Immunity, and Global Optimality")

Under i.i.d. paired augmentations, the population joint matrix $bold(P)^*$ is positive semidefinite—a spectral constraint that rules out off-diagonal permutation solutions (such as involutions swapping labels across views) and enforces a unique diagonal global maximizer.

#theorem(title: "Automatic Collapse Immunity and Unique Diagonal Global Maximizer", origin: "Original Proof")[
  Let $bold(x) tilde.op P_X$ and $g_1, g_2 tilde.op P_cal(T)$ independently given $bold(x)$. Set $bar(bold(p))(bold(x)) = bb(E)_g [Phi_theta (g(bold(x)))]$ and $bold(P)^* = bb(E)_(bold(x), g_1, g_2) [Phi_theta (g_1 (bold(x))) Phi_theta (g_2 (bold(x)))^top]$. For any $lambda >= 1$:
  + $bold(P)^*$ is symmetric and positive semidefinite ($bold(P)^* succ.eq 0$).
  + If $bold(P)^*$ has $m <= K$ active clusters ($P^*_c > 0$), then $I_lambda (bold(P)^*) <= (2 lambda - 1) ln m$, with strict gap $(2 lambda - 1) ln(K / m)$ below $(2 lambda - 1) ln K$; full collapse ($m = 1$) gives $I_lambda (bold(P)^*) = 0$.
  + $I_lambda (bold(P)^*) = (2 lambda - 1) ln K$ iff $bold(P)^* = (1 / K) bold(I)_K$, which holds iff $Phi_theta (g(bold(x))) = Phi_theta (bold(x)) in {bold(e)_1, ..., bold(e)_K}$ a.s. with $P^*_c = 1 / K$.
]

#proof[
  *Step 1 (PSD structure).* Conditional on $bold(x)$, $Phi_theta (g_1 (bold(x)))$ and $Phi_theta (g_2 (bold(x)))$ are i.i.d. with mean $bar(bold(p))(bold(x))$. By the tower property,
  $
    bold(P)^* &= bb(E)_(bold(x)) [ bb(E)_(g_1) [Phi_theta (g_1 (bold(x)))] \, bb(E)_(g_2) [Phi_theta (g_2 (bold(x)))]^top ] \
              &= bb(E)_(bold(x)) [ bar(bold(p))(bold(x)) bar(bold(p))(bold(x))^top ].
  $
  Thus $bold(v)^top bold(P)^* bold(v) = bb(E)_(bold(x)) [(bold(v)^top bar(bold(p))(bold(x)))^2] >= 0$ for all $bold(v) in bb(R)^K$, so $bold(P)^* succ.eq 0$.

  *Step 2 (Entropy bound).* By symmetry $H(bold(P)^*) = H(Z) + H(Z' | Z)$, giving
  $
    I_lambda (bold(P)^*) &= (2 lambda - 1) H(Z) - H(Z' | Z) \
                         &<= (2 lambda - 1) ln m <= (2 lambda - 1) ln K,
  $
  since $H(Z' | Z) >= 0$ and $H(Z) <= ln m$ on $m$ active clusters. At $m = 1$, $H(Z) = H(Z' | Z) = 0$, so $I_lambda (bold(P)^*) = 0$.

  *Step 3 (Uniqueness of diagonal optimum).* Achieving $(2 lambda - 1) ln K$ forces $m = K$, $P^*_c = 1 / K$, and $H(Z' | Z) = 0$, so every row of $bold(P)^*$ has at most one positive entry. If $P^*_(c c') > 0$ for some $c != c'$, then $bold(P)^* succ.eq 0$ implies the $2 times 2$ principal minor satisfies $P^*_(c c) P^*_(c' c') >= (P^*_(c c'))^2 > 0$, forcing $P^*_(c c) > 0$. Row $c$ then has at least two positive entries ($P^*_(c c), P^*_(c c')$), giving $H(Z' | Z = c) > 0$, a contradiction. Thus $bold(P)^* = (1 / K) bold(I)_K$.

  *Step 4 (Exact invariance).* Since $bar(p)_c (bold(x)) in [0, 1]$, $P^*_(c c) = bb(E)[bar(p)_c^2] <= bb(E)[bar(p)_c] = P^*_c$, with equality iff $bar(p)_c (bold(x)) in {0, 1}$ a.s., hence $bar(bold(p))(bold(x)) in {bold(e)_1, ..., bold(e)_K}$ a.s. Because $Phi_theta (g(bold(x))) in Delta^(K-1)$, its conditional expectation is an extreme point $bold(e)_c$ iff $Phi_theta (g(bold(x))) = bold(e)_c$ a.s. over $g tilde.op P_cal(T)$.
]

#subsec("B", "Exact Gradient Dynamics and Finite-Batch Plug-In Bias")

Because $I_lambda (bold(P))$ is nonlinear in the batch-averaged matrix $bold(P)$, finite minibatches incur a systematic upward estimation bias. We derive both the exact gradient and the exact $O(K^2 / B)$ plug-in bias.

#lemma(title: "Exact Log-Ratio Gradient and Finite-Batch Plug-In MI Bias", origin: "Original Derivation")[
  Let $bold(P) = (1 / (2 B)) sum_(i=1)^B (bold(p)_i (bold(p)'_i)^top + bold(p)'_i bold(p)_i^top) > 0$.
  + The gradient with respect to $bold(p)_i$ is
    $
      nabla_(bold(p)_i) I_lambda (bold(P)) &= 1 / B bold(G) bold(p)'_i, \
      G_(c c') &= ln (P_(c c') / (P_c^lambda P_(c')^lambda)) + 1 - 2 lambda,
    $
    where the constant $(1 - 2 lambda)$ lies in the softmax Jacobian nullspace.
  + For $B$ i.i.d. hard cluster pairs from a symmetric interior law $bold(P)^* > 0$, the symmetrized plug-in MI estimator ($lambda = 1$) has exact bias
    $
      bb(E)[I(hat(bold(P))_B)] - I(bold(P)^*) &= 1 / (2 B) ( (K(K - 1)) / 2 + 1 \
        &quad - sum_(c=1)^K P^*_(c c) / P^*_c ) + O(B^(-2)).
    $
    At uniform initialization ($P^*_(c c') = 1 / K^2$), this equals $K(K - 1) / (4 B) + O(B^(-2))$.
]

#proof[
  *Part 1.* Treating $P_(a b)$ as unconstrained with row sums $P_a$ and column sums $Q_b$, differentiating $I_lambda (bold(P))$ yields
  $
    (partial I_lambda) / (partial P_(a b)) &= (ln P_(a b) + 1) - lambda (ln P_a + 1) \
      &quad - lambda (ln Q_b + 1) = G_(a b),
  $
  using $Q_b = P_b$. By $(partial P_(a b)) / (partial p_(i, c)) = (delta_(a c) p'_(i, b) + delta_(b c) p'_(i, a)) / (2 B)$ and symmetry $G_(a c) = G_(c a)$,
  $
    (partial I_lambda) / (partial p_(i, c)) &= 1 / (2 B) sum_(b=1)^K G_(c b) p'_(i, b) + 1 / (2 B) sum_(a=1)^K G_(a c) p'_(i, a) \
      &= 1 / B (bold(G) bold(p)'_i)_c.
  $
  The softmax Jacobian $J_"sm" (bold(p)_i) = "diag"(bold(p)_i) - bold(p)_i bold(p)_i^top$ annihilates $bold(1)_K$, so $(1 - 2 lambda)$ vanishes in logit gradients.

  *Part 2.* Let $bold(Delta) = hat(bold(P))_"raw" - bold(P)^*$, $bold(Delta)_"sym" = (bold(Delta) + bold(Delta)^top) / 2$, and $Delta_"sym", c = (Delta_c + Delta'_c) / 2$. Expanding $I(bold(P)^* + bold(Delta)_"sym")$ to second order gives
  $
    bb(E)[I(hat(bold(P))_B)] - I(bold(P)^*) &= 1 / 2 sum_(c, c') ("Var"(Delta_"sym", c c')) / P^*_(c c') \
      &quad - sum_(c=1)^K ("Var"(Delta_"sym", c)) / P^*_c + O(B^(-2)).
  $
  For multinomial counts, $"Var"(Delta_(c c')) = P^*_(c c')(1 - P^*_(c c')) / B$ and $"Cov"(Delta_(c c'), Delta_(c' c)) = -(P^*_(c c'))^2 / B$ ($c != c'$). Thus diagonal terms satisfy $"Var"(Delta_"sym", c c) / P^*_(c c) = (1 - P^*_(c c)) / B$, while off-diagonals satisfy
  $
    "Var"(Delta_"sym", c c') / P^*_(c c') &= (2 P^*_(c c')(1 - P^*_(c c')) - 2 (P^*_(c c'))^2) / (4 B P^*_(c c')) \
      &= (1 - 2 P^*_(c c')) / (2 B).
  $
  Summing over $K$ diagonals and $K(K - 1)$ off-diagonals using $sum_(c, c') P^*_(c c') = 1$ yields $(K(K + 1) / 2 - 1) / B$. For marginals, $"Cov"(Delta_c, Delta'_c) = (P^*_(c c) - (P^*_c)^2) / B$, giving
  $
    "Var"(Delta_"sym", c) / P^*_c &= (2 P^*_c (1 - P^*_c) + 2(P^*_(c c) - (P^*_c)^2)) / (4 B P^*_c) \
      &= (1 + P^*_(c c) / P^*_c - 2 P^*_c) / (2 B).
  $
  Summing over $c$ yields $(K - 2 + sum_c P^*_(c c) / P^*_c) / (2 B)$. Combining both terms in the second-order expansion gives
  $
    1 / (2 B) ( (K(K + 1)) / 2 - 1 - K + 2 - sum_(c=1)^K P^*_(c c) / P^*_c ) \
    = 1 / (2 B) ( (K(K - 1)) / 2 + 1 - sum_(c=1)^K P^*_(c c) / P^*_c ).
  $
  At $P^*_(c c') = 1 / K^2$, $sum_c P^*_(c c) / P^*_c = 1$, leaving $K(K - 1) / (4 B)$.
]

Lemma 1 explains why IIC requires $B >> K^2$ (or repeated view pairing): otherwise, the $O(K^2 / B)$ plug-in noise overwhelms early semantic gradients @ji2019iic.

#subsec("C", "Overclustering Information Decomposition and Contrastive Bounds")

Auxiliary overclustering heads ($K' > K$) empirically improve representation quality @ji2019iic by preserving intra-class sub-structure inside the shared backbone.

#proposition(title: "Exact Overclustering Information Decomposition and InfoNCE Link", origin: "Classical — Cited")[
  Let $(tilde(Z), tilde(Z)') in {1, ..., K'}^2$ be auxiliary cluster variables and $Z = pi(tilde(Z))$, $Z' = pi(tilde(Z)')$ be coarsened $K$-class variables under a surjection $pi : {1, ..., K'} -> {1, ..., K}$ @ji2019iic @strehl2002cluster @oord2018infonce.
  + The mutual information decomposes identically as
    $
      I(tilde(Z); tilde(Z)') &= I(Z; Z') + I(Z; tilde(Z)' | Z') \
        &quad + I(tilde(Z); Z' | Z) + I(tilde(Z); tilde(Z)' | Z, Z') >= I(Z; Z').
    $
  + For shared features $bold(z) = f_omega (bold(x)^((1)))$ and $bold(z)' = f_omega (bold(x)^((2)))$, the Data Processing Inequality yields $I(Z; Z') <= I(tilde(Z); tilde(Z)') <= I(bold(z); bold(z)')$, paralleling the InfoNCE variational lower bound @oord2018infonce.
]

#proof[
  Since $Z = pi(tilde(Z))$ and $Z' = pi(tilde(Z)')$ are deterministic functions, $I(tilde(Z); tilde(Z)') = I(Z, tilde(Z); Z', tilde(Z)')$. Expanding by the chain rule first on the left,
  $
    I(Z, tilde(Z); Z', tilde(Z)') &= I(Z; Z', tilde(Z)') \
      &quad + I(tilde(Z); Z', tilde(Z)' | Z),
  $
  and then on the right gives
  $
    I(Z; Z', tilde(Z)') &= I(Z; Z') + I(Z; tilde(Z)' | Z'), \
    I(tilde(Z); Z', tilde(Z)' | Z) &= I(tilde(Z); Z' | Z) \
      &quad + I(tilde(Z); tilde(Z)' | Z, Z').
  $
  Summing yields the exact four-term identity; nonnegativity of conditional mutual information implies $I(tilde(Z); tilde(Z)') >= I(Z; Z')$. Applying the Data Processing Inequality twice along the Markov chain $tilde(Z) <- bold(z) <-> bold(z)' -> tilde(Z)'$ yields $I(tilde(Z); tilde(Z)') <= I(bold(z); bold(z)')$, showing that overclustering strictly tightens the discrete MI surrogate whenever intra-class nuisance modes carry $I(tilde(Z); tilde(Z)' | Z, Z') > 0$.
]

#subsec("D", "Imbalanced-Prior Distortion and Constructive Counterexample")

Because Theorem 1 rewards maximal marginal entropy $H(Z) = ln K$, standard IIC exhibits a provable distortion under skewed class priors.

#theorem(title: "Majority-Splitting Counterexample under Class Imbalance", origin: "Original Counterexample")[
  Let $K = 2$ true semantic classes $cal(C)_1, cal(C)_2$ have proportions $pi_1 = 1 - epsilon$ and $pi_2 = epsilon$ with $epsilon in (0, 1 / 3)$. Suppose $cal(C)_1 = A_1 union A_2$ splits into equal-mass sub-modes (each $(1 - epsilon) / 2$) such that augmentations preserve semantic class exactly while leaking symmetrically between $A_1$ and $A_2$ with probability $delta in [0, 1 / 2)$. For any $lambda >= 1$:
  + The true partition achieves $I_lambda (bold(P)_"true") = (2 lambda - 1) h_2 (epsilon)$, with $h_2 (u) = -u ln u - (1 - u) ln(1 - u)$.
  + The majority-splitting partition $(A_1, A_2 union cal(C)_2)$ strictly beats the true partition ($I_lambda (bold(P)_"split") > I_lambda (bold(P)_"true")$) whenever
    $
      h_2 (delta) < (2 lambda - 1) ( h_2 ((1 - epsilon) / 2) - h_2 (epsilon) ).
    $
  At $lambda = 1$ and $epsilon = 0.05$, any leak $delta <= 0.19$ forces IIC to merge the minority class and split the majority class.
]

#proof[
  *Part 1.* Under the true partition, $bold(P)_"true" = "diag"(1 - epsilon, epsilon)$, giving $H(Z' | Z) = 0$ and $I_lambda (bold(P)_"true") = (2 lambda - 1) h_2 (epsilon)$.

  *Part 2.* Under $C'_1 = A_1$ and $C'_2 = A_2 union cal(C)_2$, the symmetric joint matrix has entries
  $
    P'_(1 1) &= (1 - epsilon) / 2 (1 - delta), \
    P'_(1 2) = P'_(2 1) &= (1 - epsilon) / 2 delta, \
    P'_(2 2) &= epsilon + (1 - epsilon) / 2 (1 - delta).
  $
  Its marginals are $(1 - epsilon) / 2$ and $(1 + epsilon) / 2$, so $H(Z') = h_2 ((1 - epsilon) / 2) > h_2 (epsilon)$ since $epsilon < (1 - epsilon) / 2 < 1 / 2$. For conditional entropy, $P(Z'_2 = 2 | Z'_1 = 1) = delta$ gives $H(Z'_2 | Z'_1 = 1) = h_2 (delta)$, while $P(Z'_2 = 1 | Z'_1 = 2) = ((1 - epsilon) / (1 + epsilon)) delta < delta < 1 / 2$ gives $H(Z'_2 | Z'_1 = 2) < h_2 (delta)$, hence $H(Z'_2 | Z'_1) < h_2 (delta)$. Thus
  $
    I_lambda (bold(P)_"split") &= (2 lambda - 1) h_2 ((1 - epsilon) / 2) - H(Z'_2 | Z'_1) \
      &> (2 lambda - 1) h_2 ((1 - epsilon) / 2) - h_2 (delta),
  $
  exceeding $I_lambda (bold(P)_"true")$ whenever the stated inequality holds. At $lambda = 1$ and $epsilon = 0.05$, $h_2 (0.475) - h_2 (0.05) approx 0.6919 - 0.1985 = 0.4934 > h_2 (0.19) approx 0.4862$. Increasing $lambda > 1$ scales the gap by $(2 lambda - 1)$, widening the failure regime.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Let $C_"net"$ denote the per-sample forward-and-backward pass cost of backbone $f_omega$ and head $Phi_theta$. Per minibatch of size $B$, computing paired probabilities costs $O(B C_"net")$, forming $bold(P)_"raw" = (1 / B) bold(P)_1 bold(P)_2^top$ costs $O(B K^2)$, evaluating the log-ratio gradient $bold(G)$ costs $O(K^2)$, and backpropagating $(1 / B) bold(G) bold(p)'_i$ costs $O(B K^2)$. Over $E$ epochs, total training time is $O(E N (C_"net" + K^2))$ (plus $O(E N (K')^2)$ for overclustering). Memory is $O(|theta| + B (r + K) + K^2)$, requiring no memory bank @oord2018infonce or $k$-NN index @vangansbeke2020scan. However, Lemma 1 requires $B = Omega(K^2)$ to control the $K(K - 1) / (4 B)$ plug-in bias. Geometrically, IIC learns smooth non-convex neural boundaries aligned with $cal(T)$-orbits and achieves exact cluster-permutation invariance (within the axiom trade-offs of @kleinberg2002impossibility).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(E N (C_"net" + K^2))$ ($+ O(E N (K')^2)$ overclustering); $O(N C_"net")$ inference.],
    [Space Complexity], [$O(|theta| + B (r + K) + K^2)$; no $k$-NN index or feature queue.],
    [Inductive Bias], [Equipartitioned, $cal(T)$-orbit-invariant discrete neural partitions.],
    [Collapse Resistance], [Provably collapse-immune ($I_lambda = 0$ vs. $(2 lambda - 1) ln K$, Theorem II.A).],
    [Batch-Size Coupling], [Requires $B >> K^2$ so plug-in MI bias $Theta(K^2 / B)$ stays negligible.],
    [Class-Balance Prior], [Implicit uniform marginal pressure; degrades on skewed classes (Theorem II.D).],
    [Hyperparameters], [Counts $K, K'$, weight $lambda >= 1$, batch $B$, transformation family $cal(T)$.],
    [Streaming / Online], [Native minibatch SGD; directly supports online arrival streams.],
  ),
  caption: [Operational and geometric profile of IIC (Alg. 188).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Unsupervised Natural Image Categorization:* On balanced visual benchmarks (CIFAR-10, STL-10, MNIST), geometric and Sobel-filtered photometric augmentations preserve object identity while erasing pixel shortcuts; maximizing $I(Z; Z')$ isolates semantic invariance end-to-end @ji2019iic.
+ *Unsupervised Dense Semantic Segmentation:* Pairing spatially adjacent pixel patches $(bold(x)_(u, v), bold(x)_(u + delta_u, v + delta_v))$ in satellite (Potsdam) and street scenes (COCO-Stuff) maps IIC directly to pixel-level clustering without manual augmentations @ji2019iic.
+ *Multi-Sensor Remote Sensing Land-Cover Mapping:* Cross-date or cross-spectral satellite tiles form natural invariant pairs; IIC's PSD joint objective aligns cross-modal predictions onto discrete land-cover types without collapsing onto dominant terrain.
+ *Digital Histopathology Tile Phenotyping:* Pairing tissue patches under $D_4$ rotations and stain perturbations separates morphological phenotypes while auxiliary overclustering ($K' > K$) absorbs slide-level staining artifacts.
+ *Multi-View Tabular and Acoustic Clustering:* With calibrated feature dropout, additive noise, or spectrogram masking, IIC provides a lightweight single-stage neural clustering primitive returning direct posterior probabilities.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Three mathematical regimes govern failure: (1) *skewed class priors* trigger majority-splitting collapse (Theorem II.D); (2) *large $K >= 100$* inflates the $O(K^2)$ matrix and $K(K - 1) / (4 B)$ plug-in bias (Lemma 1); and (3) *misspecified augmentations* that leave color shortcuts unmasked (e.g., raw RGB without Sobel filtering @ji2019iic) cause $I(Z; Z')$ to lock onto low-level texture invariants.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Balanced Visual Objects], [Strong: geometric/Sobel orbits match semantics; overclustering helps.], [Avoid raw RGB shortcuts; compare SCAN (Alg. 189).],
    [Dense Pixel Segmentation], [Strong: spatial adjacency gives natural parameter-free pairs.], [Watch tiny foreground classes; tune $lambda$.],
    [Long-Tailed / Skewed Data], [Poor: marginal entropy splits majority modes (Theorem II.D).], [Prefer SCAN (Alg. 189) or GMM / VaDE.],
    [High-$K$ Taxonomies ($K > 200$)], [Weak: $O(K^2 / B)$ plug-in bias limits small-batch training.], [Prefer SwAV (Alg. 190) with feature queue.],
    [Tabular / Sensor Data], [Moderate: fast $O(K^2)$ step with calibrated noise views.], [Validate perturbations; compare DEC / K-Means++.],
  ),
  caption: [Domain suitability and failure-mode matrix for IIC (Alg. 188).],
)

#sec("V", "Structural Relations and Empirical Behavior")

IIC occupies the pure information-theoretic vertex connecting paired-view neural clustering to Contrastive Clustering @li2021contrastive and SwAV @caron2020swav (Alg. 190).

#proposition(title: "Structural Equivalence Between IIC, SwAV, and Contrastive Clustering", origin: "Original Derivation")[
  Let $bold(P)_"raw" = (1 / B) sum_(i=1)^B bold(p)_i (bold(p)'_i)^top$ on minibatch $cal(B)$.
  + Constraining soft assignments onto SwAV's equipartition polytope $(1 / B) sum_i bold(p)_i = (1 / B) sum_i bold(p)'_i = (1 / K) bold(1)_K$ @caron2020swav fixes $H(Z) = ln K$, reducing maximizing $I_lambda (bold(P))$ identically to minimizing joint entropy $H(bold(P))$.
  + Writing cluster column vectors $bold(u)_c = (p_(1, c), ..., p_(B, c))^top in bb(R)^B_(>= 0)$, entry $P_("raw", c c') = (1 / B) chevron.l bold(u)_c, bold(u)'_(c') chevron.r$ is an unnormalized cross-view cluster inner product; replacing $I(bold(P))$ by cosine InfoNCE over ${(bold(u)_c, bold(u)'_c)}_(c=1)^K$ plus marginal entropy recovers Contrastive Clustering @li2021contrastive.
]

#proof[
  Under exact equipartition $P_c = 1 / K$, $H(Z) = ln K$ is constant, so $I_lambda (bold(P)) = 2 lambda ln K - H(bold(P))$ is maximized iff $H(bold(P))$ is minimized (uniquely at $(1 / K) bold(I)_K$ by Theorem 1). Entrywise expansion gives $P_("raw", c c') = (1 / B) sum_i p_(i, c) p'_(i, c') = (1 / B) chevron.l bold(u)_c, bold(u)'_(c') chevron.r$, normalizing overlaps by $ell_1$ column sums $B P_c$ rather than the $ell_2$ norms of cosine InfoNCE in @li2021contrastive.
]

*Empirical Behavior on the ClusterOverview Benchmark.* On the standardized 4-dataset ClusterOverview benchmark @rousseeuw1987silhouettes @hubert1985comparing, IIC attains *Iris* $"Sil" = 0.360$, $"ARI" = 0.531$; *Customers* $"Sil" = 0.257$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.298$ ($"ARI" =$ N/A), with mean runtime *4.143 ms*—the fastest among Algs. 188--191, reflecting its single-stage $O(K^2)$ contingency update without $k$-NN indexing, Sinkhorn loops, or offline K-Means sweeps. Crucially, IIC achieves the highest *High-D* Silhouette ($0.298$) of the quartet (far exceeding DeepCluster's $0.109$ and SCAN's $0.211$), since optimizing a $K times K$ probability matrix avoids Euclidean distance concentration @beyer1999nearest. Strong *Multi-Shape* separation ($0.712$) confirms flexible non-convex boundaries, while lower *Customers* Silhouette ($0.257$ vs. SCAN's $0.319$) aligns with Theorem II.D's uniform-marginal bias on skewed tabular segments.
