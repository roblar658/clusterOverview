#import "../template/ieee.typ": *

#alg_entry(
  id: 126,
  name: "SLPA (Speaker-Listener Label Propagation)",
  category: "Cat. V: Graph, Spectral & Network",
  year: "2011",
  paradigm: "Memory-Based Stochastic Label Propagation",
  complexity: $cal(O)(T m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical Label Propagation (LPA) @raghavan2007lpa detects disjoint communities in near-linear time by iteratively replacing each vertex's label with the majority label of its neighbors. Because LPA discards all prior label history at every step, it cannot represent overlapping memberships and can oscillate indefinitely on bipartite subgraphs. Xie, Szymanski, and Liu @xie2011slpa introduced *SLPA* (Speaker-Listener Label Propagation Algorithm), extending LPA into a stochastic *memory-based* information-diffusion process that preserves the full empirical history of labels observed at each vertex and extracts arbitrary overlapping covers via post-hoc frequency thresholding.

Let $G = (V, E)$ be an undirected (or directed/weighted) graph with $n = abs(V)$ vertices, $m = abs(E)$ edges, and neighbor sets $N(i) = {j in V : (i, j) in E}$. Every vertex $i in V$ maintains a dynamic multiset (memory) $M_i (t)$ of size $t + 1$ at discrete iteration $t = 0, 1, ..., T$. Initially, $M_i (0) = {i}$ holds vertex $i$'s own unique label. During each asynchronous sweep $t$, vertices are visited in a uniform random permutation $pi_t$: the active vertex $i$ plays the *listener*, while each neighbor $j in N(i)$ plays a *speaker* that broadcasts one label sampled from its empirical memory distribution.

#definition(title: "Speaker-Listener Memory Dynamics and Thresholded Cover")[
  Let $N_i^((t)) (c)$ denote the multiplicity of label $c in {1, ..., n}$ in vertex $i$'s memory $M_i (t)$ after $t$ iterations, defining the empirical label distribution $hat(p)_i^((t)) (c) = N_i^((t)) (c) / (t + 1)$. When vertex $i$ is updated as listener at step $t$:
  1. *Speaking rule*: Each neighbor $j in N(i)$ independently draws one label $L_j ~ hat(p)_j^((t - 1))$ (using the most recently updated memory of $j$ in the asynchronous sweep).
  2. *Listening rule*: Listener $i$ computes the received vote counts $V_i (c) = sum_(j in N(i)) bb(I)[L_j = c]$, selects a maximizing label $ell_i^(*) in arg max_c V_i (c)$ (breaking ties uniformly at random), and appends $ell_i^(*)$ to its memory: $M_i (t) = M_i (t - 1) union.sq {ell_i^(*)}$.
  3. *Post-processing cover*: Given probability threshold $r in (0, 1]$, vertex $i$'s overlapping community assignment is
     $
       cal(C)_i (r) = { c in {1, ..., n} : hat(p)_i^((T)) (c) >= r } ,
     $
     (assigning $arg max_c hat(p)_i^((T)) (c)$ if $cal(C)_i (r) = emptyset$ so every vertex belongs to at least one community).
]

Equivalently, SLPA can be interpreted as a coupled network of generalized Pólya-like reinforced urn processes interacting through nonlinear majority-vote filters: vertices deep inside a cohesive community quickly accumulate a single dominant label, whereas overlapping boundary vertices receive competing speaker emissions from multiple adjacent communities and preserve a multi-modal memory histogram.

#algbox(
  title: "SLPA (Speaker-Listener Label Propagation)",
  inputs: [Graph $G = (V, E)$, iteration horizon $T$ (default $T = 20$), post-processing threshold $r in (0, 1]$.],
  outputs: [Overlapping community assignment ${cal(C)_i (r)}_(i in V)$ and community sets ${C_c}$.],
)[
  1. *Initialize memory*: For each $i in V$, initialize memory array $M_i <- [i]$ of length $1$.
  2. *Asynchronous speaker-listener sweeps*: For $t = 1, ..., T$:
     - Draw a uniform random permutation $pi_t$ of $V$.
     - For each listener $i$ in order $pi_t$:
       - Initialize empty vote map $V_i$.
       - For each speaker $j in N(i)$, sample index $r ~ "Unif"{0, ..., abs(M_j) - 1}$, emit $ell_j = M_j [r]$, and increment $V_i [ell_j] <- V_i [ell_j] + 1$.
       - Select $ell_i^(*) in arg max_c V_i [c]$ (uniform random tie-breaking) and append $ell_i^(*)$ to $M_i$.
  3. *Threshold empirical distributions*: For each $i in V$, compute $hat(p)_i (c) = N_i (c) / (T + 1)$ and set $cal(C)_i (r) = {c : hat(p)_i (c) >= r}$ (fallback to mode if empty).
  4. *Assemble communities*: Form $C_c = {i in V : c in cal(C)_i (r)}$, discarding empty sets and optionally merging nested duplicate subsets.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four mathematical results for SLPA: an exact deterministic cardinality bound and disjointness threshold, a mean-field majority-amplification theorem characterizing when overlapping bridges survive versus vanish, a finite-sample Azuma-Hoeffding iteration bound for reliable overlap recovery, and a formal proof contrasting Cesàro oscillation damping against monochromatic collapse.

#subsec("A", "Deterministic Overlap Cap and Disjointness Threshold")

The post-processing parameter $r in (0, 1]$ acts as a sharp algebraic control knob governing maximum multi-community multiplicity per vertex.

#proposition(title: "Exact Overlap Cardinality Bound and Disjointness Threshold", origin: "Original Proof")[
  For every graph $G$, every horizon $T >= 1$, and every threshold $r in (0, 1]$, the SLPA membership set $cal(C)_i (r)$ of any vertex $i in V$ satisfies:
  1. $abs(cal(C)_i (r)) <= floor(1 / r)$ (without fallback, and $max(1, floor(1 / r))$ with mode fallback).
  2. Whenever $r > 1 / 2$, every vertex belongs to *at most one* community ($abs(cal(C)_i (r)) <= 1$), guaranteeing a strictly disjoint partition.
  3. At $r = 1 / 2$, a vertex belongs to two communities if and only if $T$ is odd and its memory splits exactly equally ($N_i (c_1) = N_i (c_2) = (T + 1) / 2$).
]

#proof[
  Because $hat(p)_i^((T))$ is a probability distribution over ${1, ..., n}$, summing over the retained label set $cal(C)_i (r)$ yields
  $
    1 = sum_(c = 1)^n hat(p)_i^((T)) (c) >= sum_(c in cal(C)_i (r)) hat(p)_i^((T)) (c) >= r abs(cal(C)_i (r)) .
  $
  Dividing by $r > 0$ gives $abs(cal(C)_i (r)) <= 1 / r$, and integrality of $abs(cal(C)_i (r))$ tightens this to $floor(1 / r)$. When $r > 1 / 2$, we have $1 / r < 2$, so $floor(1 / r) = 1$, forcing disjoint communities (and with fallback to the unique or tie-broken mode when $max_c hat(p)_i^((T)) (c) < r$, every vertex receives exactly one label). Finally, when $r = 1 / 2$, retaining two distinct labels $c_1 != c_2$ requires $hat(p)_i^((T)) (c_1) >= 1 / 2$ and $hat(p)_i^((T)) (c_2) >= 1 / 2$, which forces equality $N_i (c_1) = N_i (c_2) = (T + 1) / 2$; because $N_i (c)$ is an integer, $T + 1$ must be even, i.e., $T$ must be odd.
]

#subsec("B", "Majority-Vote Amplification and Bridge-Overlap Survival Law")

Understanding why SLPA preserves *symmetric* overlaps while potentially erasing *asymmetric* bridges requires analyzing the nonlinear majority-vote operator acting on speaker distributions.

#theorem(title: "Majority Amplification Derivative and Asymmetric Bridge Erasure", origin: "Original Derivation")[
  Let $Phi_d (p) = "Pr"("Binom"(d, p) > d / 2) + (1 / 2) "Pr"("Binom"(d, p) = d / 2)$ denote the probability that a listener of degree $d$ selects label $a$ when $d$ speakers emit label $a$ independently with probability $p$ and label $b$ with probability $1 - p$. Then:
  1. *Sharp S-curve instability at $p = 1 / 2$*: For every odd $d >= 3$, $Phi_d$ has fixed points ${0, 1 / 2, 1}$ and exact derivative
     $
       Phi'_d (p) = d binom(d - 1, (d - 1) / 2) (p (1 - p))^((d - 1) / 2) ,
     $
     giving $Phi'_d (1 / 2) = d binom(d - 1, (d - 1) / 2) 2^(-(d - 1)) = sqrt(2 d / pi) (1 + cal(O)(1 / d)) > 1$. Thus $p = 1 / 2$ is an *unstable* fixed point in homogeneous cores, driving rapid consensus toward $0$ or $1$.
  2. *Exact survival of symmetric bridges*: Suppose a bridge vertex $w$ has $d_A$ neighbors in a stabilized pure-$a$ core and $d_B$ neighbors in a stabilized pure-$b$ core (with no other neighbors). If $d_A = d_B$, then $hat(p)_w^((T)) (a) -> 1 / 2$ almost surely as $T -> infinity$, with $sqrt(T) (hat(p)_w^((T)) (a) - 1 / 2) -> cal(N)(0, 1 / 4)$ in distribution, recovering both communities for any $r < 1 / 2$.
  3. *Deterministic erasure of asymmetric pure-core bridges*: Under the same pure-core setting, if $d_A > d_B$, then listener $w$ appends label $a$ deterministically at every step, yielding $hat(p)_w^((T)) (a) = T / (T + 1) -> 1$ and $hat(p)_w^((T)) (b) = 0$; hence majority listening completely erases strictly asymmetric overlap between pure cores.
]

#proof[
  For part (1), when $d = 2 k + 1$ is odd ($k >= 1$), there are no ties and $Phi_d (p)$ equals the regularized incomplete beta function $I_p (k + 1, k + 1)$:
  $
    Phi_d (p) = sum_(j = k + 1)^(2 k + 1) binom(2 k + 1, j) p^j (1 - p)^(2 k + 1 - j) .
  $
  Differentiating term-by-term with respect to $p$ yields
  $
    Phi'_d (p) &= sum_(j = k + 1)^d j binom(d, j) p^(j - 1) (1 - p)^(d - j) - sum_(j = k + 1)^(d - 1) (d - j) binom(d, j) p^j (1 - p)^(d - j - 1) .
  $
  Using the identity $(d - j) binom(d, j) = (j + 1) binom(d, j + 1)$ and shifting the index $j' = j + 1$ in the second sum cancels every term except $j = k + 1$, leaving
  $
    Phi'_d (p) = (k + 1) binom(2 k + 1, k + 1) p^k (1 - p)^k = d binom(d - 1, k) (p (1 - p))^k .
  $
  Evaluating at $p = 1 / 2$ gives $Phi'_d (1 / 2) = d binom(2 k, k) 2^(-2 k)$. Since $binom(2 k, k) >= 4^k / sqrt(pi (k + 1 / 2))$ (and by Stirling's formula $binom(2 k, k) 2^(-2 k) = (pi k)^(-1/2) (1 + cal(O)(1 / k))$), we obtain $Phi'_d (1 / 2) = sqrt(2 d / pi) (1 + cal(O)(1 / d))$, which strictly exceeds $1$ for all $d >= 3$ (at $d = 3$, $Phi'_3 (1 / 2) = 3 / 2 > 1$). Because $Phi_d (p) > p$ on $(1 / 2, 1)$ and $Phi_d (p) < p$ on $(0, 1 / 2)$, the fixed point $1 / 2$ is repelling while $0$ and $1$ are locally attracting.

  For part (2), when $d_A = d_B$ and all $A$-neighbors emit $a$ while all $B$-neighbors emit $b$, listener $w$ faces an exact $d_A$-to-$d_A$ tie at every iteration $t >= 1$. Uniform random tie-breaking appends i.i.d. indicators $Y_t ~ "Bern"(1 / 2)$ for label $a$, so $hat(p)_w^((T)) (a) = (1 / (T + 1)) sum_(t = 1)^T Y_t$ (assuming $w$'s initial self-label differs from $a$ and $b$). By the Strong Law of Large Numbers and Central Limit Theorem, $hat(p)_w^((T)) (a) -> 1 / 2$ almost surely with asymptotic variance $1 / (4 T)$.

  For part (3), when $d_A >= d_B + 1$, the received vote count satisfies $V_w (a) = d_A > d_B = V_w (b)$ deterministically at every step, so $w$ appends $a$ every round and never records $b$ in $M_w (T)$, giving $hat(p)_w^((T)) (b) = 0$.
]

Theorem 1 clarifies an important practical nuance: in real networks, SLPA still recovers asymmetric bridges ($d_A > d_B$) *only* when transient burn-in fluctuations or mixed speaker memories in the surrounding overlap zone keep the winning probability of label $b$ above $r$.

#subsec("C", "Finite-Sample Concentration and Horizon Complexity")

Because SLPA terminates after a finite horizon $T$ (commonly $T = 20$ @xie2011slpa), choosing $T$ too small causes high Monte Carlo variance in $hat(p)_i^((T))$ and misses genuine overlaps. The following original concentration bound quantifies the exact minimum horizon $T$.

#lemma(title: "Azuma-Hoeffding Horizon Bound for Stable Overlap Recovery", origin: "Original Proof")[
  Suppose after a burn-in of $T_0$ iterations, the filtration $cal(F)_(t - 1)$ satisfies $"Pr"(ell_w^((t)) = c | cal(F)_(t - 1)) >= q_c$ for all $t > T_0$ at an overlapping vertex $w$ and true community $c$, with margin $Delta = q_c - r > 0$. Then for every horizon $T > T_0 (1 + r / Delta) + r / Delta$, the probability that SLPA fails to recover community $c$ at vertex $w$ obeys
  $
    "Pr"(c in.not cal(C)_w (r)) <= exp(- 2 ((T - T_0) Delta - r (T_0 + 1))^2 / (T - T_0)) .
  $
  In particular, recovering $M_"ov"$ overlapping memberships simultaneously with probability at least $1 - delta$ requires
  $
    T >= T_0 + (2 r (T_0 + 1)) / Delta + (2 / Delta^2) log(M_"ov" / delta) .
  $
]

#proof[
  Let $Y_t = bb(I)[ell_w^((t)) = c]$ for $t = T_0 + 1, ..., T$ and let $mu_t = bb(E)[Y_t | cal(F)_(t - 1)] >= q_c$. Then $Z_t = Y_t - mu_t$ is a martingale difference sequence taking values in $[-mu_t, 1 - mu_t]$, an interval of conditional width $1$. Because $N_w^((T)) (c) >= sum_(t = T_0 + 1)^T Y_t$, the miss event ${N_w^((T)) (c) < r (T + 1)}$ implies
  $
    sum_(t = T_0 + 1)^T Z_t < r (T + 1) - sum_(t = T_0 + 1)^T mu_t <= r (T + 1) - (T - T_0) q_c = -s ,
  $
  where $s = (T - T_0) Delta - r (T_0 + 1) > 0$ under the stated condition on $T$. By the Hoeffding-Azuma lemma for bounded martingale differences on intervals of length $1$, $bb(E)[e^(-lambda Z_t) | cal(F)_(t - 1)] <= e^(lambda^2 / 8)$ for all $lambda > 0$. Applying Markov's inequality and optimizing $lambda = 4 s / (T - T_0)$ gives
  $
    "Pr"(sum_(t = T_0 + 1)^T Z_t < -s) <= exp(- 2 s^2 / (T - T_0)) .
  $
  Setting $T - T_0 >= 2 r (T_0 + 1) / Delta$ ensures $s >= (T - T_0) Delta / 2$, so the exponent is at most $-(T - T_0) Delta^2 / 2$. Solving $M_"ov" exp(-(T - T_0) Delta^2 / 2) <= delta$ via the union bound yields the stated horizon bound.
]

#subsec("D", "Oscillation Damping vs Monochromatic Collapse")

SLPA resolves the classic period-2 oscillation of synchronous LPA @raghavan2007lpa via time-averaged memories, yet inherits LPA's absorbing monochromatic states on weakly modular networks.

#theorem(title: "Cesàro Oscillation Damping and Absorbing Monochromatic Flood", origin: "Classical — Cited")[
  Let $hat(p)_i^((T))$ denote the empirical memory distribution of vertex $i$ in SLPA @xie2011slpa.
  1. *Deterministic Cesàro damping*: For every graph $G$ (including complete bipartite graphs $K_(a, b)$ where synchronous LPA @raghavan2007lpa cycles perpetually), every vertex $i in V$, and every step $T >= 0$,
     $
       norm(hat(p)_i^((T + 1)) - hat(p)_i^((T)))_1 <= 2 / (T + 2) -> 0 .
     $
  2. *Monochromatic absorption*: On any connected graph $G$, every homogeneous configuration in which $M_i (t) = {c, ..., c}$ for all $i in V$ is an absorbing state of the SLPA Markov chain. Moreover, on weakly modular graphs with external mixing fraction $mu > 1 / 2$, the majority instability (Theorem 1) drives giant-community flooding unless $T$ is kept finite (@raghavan2007lpa, @xie2011slpa).
]

#proof[
  For part (1), let $e_ell$ be the standard basis vector of the label $ell = ell_i^((T + 1))$ appended at step $T + 1$. By definition of the empirical frequency vector,
  $
    hat(p)_i^((T + 1)) = ((T + 1) hat(p)_i^((T)) + e_ell) / (T + 2) ,
  $
  and therefore
  $
    hat(p)_i^((T + 1)) - hat(p)_i^((T)) = (e_ell - hat(p)_i^((T))) / (T + 2) .
  $
  Taking the $ell_1$ norm directly yields
  $
    norm(hat(p)_i^((T + 1)) - hat(p)_i^((T)))_1 &= (1 / (T + 2)) ((1 - hat(p)_i^((T)) (ell)) + sum_(c != ell) hat(p)_i^((T)) (c)) \
                                               &= (2 (1 - hat(p)_i^((T)) (ell))) / (T + 2) <= 2 / (T + 2) .
  $
  Thus step-to-step changes vanish at rate $cal(O)(1 / T)$, smoothing bipartite alternations into stable fractional memberships rather than perpetual label flips.

  For part (2), if every vertex memory contains only label $c$, then every speaker emits $c$ with probability $1$ and every listener receives $V_i (c) = "deg"(i)$ votes for $c$ and $0$ for all other labels, appending $c$ with probability $1$; hence monochromatic states are absorbing. When external mixing $mu > 1 / 2$, vertices inside a planted cluster receive more external than internal votes, so by Theorem 1(1) ($Phi'_d (1 / 2) > 1$) a globally advantaged label amplifies super-linearly across community boundaries and floods the graph into a single giant cluster as $T -> infinity$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time and Space Complexity.* Storing each vertex's memory $M_j$ as a contiguous array of length $t$ allows every speaker emission `M_j[rand(0..t-1)]` to execute in exact $cal(O)(1)$ time without binary searching a cumulative distribution. One full listener-speaker sweep visits every edge twice (for undirected graphs), requiring $cal(O)(sum_(i in V) "deg"(i)) = cal(O)(m)$ arithmetic and hash-map vote operations. Across $T$ iterations (empirically $T = 20$ suffices on million-node graphs @xie2011slpa), total propagation time is $cal(O)(T m)$, plus $cal(O)(T n)$ for histogram thresholding—making SLPA one of the fastest overlapping community detectors available. Memory consumption is $cal(O)(T n + m)$ for array-backed memories and adjacency lists, and every speaker draw has excellent sequential cache locality when memory arrays are stored contiguously.

*Operational Characteristics.* SLPA requires neither the community count $K$ nor global graph statistics (such as total edge count $m$ or matrix factorizations), operating purely through local neighbor exchanges. Its sole structural hyperparameters are the iteration horizon $T$ (trading Monte Carlo variance against giant-label flood risk) and the overlap threshold $r in (0, 0.5]$ (controlling maximum overlap multiplicity via Proposition 1).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$cal(O)(T m + T n) = cal(O)(T m)$ on connected graphs (near-linear in edges).],
    [Space complexity], [$cal(O)(T n + m)$ for array-backed vertex memories and adjacency lists.],
    [Cluster geometry], [Arbitrary sparse graph modules and overlapping covers; no convexity assumption.],
    [Overlap / noise], [Native overlapping covers controlled by $r$; no explicit outlier/homeless class.],
    [Scalability], [Purely local neighbor message passing; scales easily to $10^7$+ edges.],
    [Hyperparameters], [Horizon $T$ (default $20$) and probability threshold $r in [0.01, 0.5]$.],
    [Determinism], [Stochastic due to random schedules, speaker draws, and tie-breaking.],
  ),
  caption: [Operational and theoretical profile of SLPA.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Massive Online Social and Communication Graphs.* Social platforms contain tens of millions of accounts with multi-group affiliation where $cal(O)(n^2)$ statistical tests or iterative matrix factorizations are expensive. SLPA's $cal(O)(T m)$ local propagation discovers overlapping interest groups in seconds without knowing $K$ in advance @xie2011slpa.
2. *Protein-Protein Interaction and Metabolic Multi-Module Discovery.* Biological networks feature sparse topologies with multi-functional proteins bridging pathways. SLPA's memory histograms capture multi-pathway participation directly through empirical label frequencies without requiring dense cliques as in CPM @palla2005cpm or BigCLAM @yang2013bigclam.
3. *Streaming and Dynamic Social Network Monitoring.* Because SLPA updates vertex memories incrementally via local speaker broadcasts, replacing the unbounded memory with a sliding window adapts naturally to evolving edge streams without full recomputation.
4. *Financial Fraud-Ring and Shared-Mule Detection.* Transaction networks frequently exhibit overlapping illicit rings sharing common mule accounts or shell entities. SLPA's ability to surface vertices with split empirical memories ($0.2 <= hat(p)_i (c) <= 0.5$) flags multi-ring intermediaries at linear transaction-graph scale.
5. *Decentralized Wireless and IoT Self-Organization.* In ad-hoc sensor and peer-to-peer networks lacking a central coordinator, SLPA requires zero global state: every node communicates exclusively with one-hop neighbors to self-organize into overlapping routing zones and gateway nodes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SLPA fails or degrades in four identifiable regimes: (i) *high-mixing / weak-community networks* ($mu > 0.5$), where absorbing monochromatic collapse (Theorem 2) produces a single trivial giant cluster; (ii) *strongly asymmetric overlapping bridges* ($d_A >> d_B$), where pure majority listening erases minority memberships (Theorem 1(3)); (iii) *unreplicated single-run deployments*, because random asynchronous permutations induce noticeable run-to-run variance (ARI/NMI fluctuations) requiring consensus aggregation; and (iv) *networks with severe background noise*, where SLPA assigns every noise vertex to some community rather than isolating homeless vertices as OSLOM does @lancichinetti2011oslom.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Massive social networks], [Excellent: $cal(O)(T m)$ linear scaling + overlap.], [Keep $T approx 20$ to avoid giant-cluster flood.],
    [Sparse $k$-NN manifolds], [Strong: local majority voting respects neighborhoods.], [Run multi-seed consensus to reduce variance.],
    [Decentralized / P2P], [Excellent: zero global state; 1-hop messages.], [Use bounded sliding memory window.],
    [High-noise / weak SBM], [Weak: floods into giant cluster when $mu > 0.5$.], [Prefer OSLOM @lancichinetti2011oslom or spectral/SBM methods.],
    [Asymmetric overlap], [Moderate: majority filter erases minority bridge.], [Prefer BigCLAM @yang2013bigclam or weighted speaking rules.],
  ),
  caption: [Domain suitability matrix for SLPA.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SLPA sits directly at the core of the label-propagation family alongside classical LPA @raghavan2007lpa, Chinese Whispers @biemann2006chinese, and Fluid Communities @pares2018fluidc, while serving as the linear-time algorithmic counterpart to BigCLAM @yang2013bigclam and OSLOM @lancichinetti2011oslom.

#proposition(title: "Exact Memory-Window Reduction to Classical LPA", origin: "Original Proof")[
  Define *windowed SLPA* with buffer length $W in bb(N) union {infinity}$ by retaining only the $W$ most recent labels in $M_i (t)$ (discarding the oldest entry once $abs(M_i) > W$). Then:
  1. For $W = 1$, every speaker distribution collapses to a point mass $hat(p)_j^((t)) = e_(ell_j^((t)))$, speaker emission becomes deterministic, and windowed SLPA reduces *identically* to asynchronous Raghavan-Albert-Kumara LPA @raghavan2007lpa (with $r = 1$ recovering the exact disjoint label partition).
  2. For $W = infinity$, windowed SLPA equals full SLPA @xie2011slpa; intermediate finite $W$ interpolates smoothly between memoryless greedy propagation and full Cesàro-damped overlapping history, bounding space at $cal(O)(W n + m)$.
]

#proof[
  When $W = 1$, initialization sets $M_i (0) = [i]$. At every update of vertex $j$, appending the new majority label $ell_j^(*)$ and retaining the latest $W = 1$ entry yields $M_j = [ell_j^(*)]$, so $hat(p)_j (ell_j^(*)) = 1$ and $hat(p)_j (c) = 0$ for $c != ell_j^(*)$. Consequently, when neighbor $i$ listens to $j$, speaker $j$ emits $ell_j^(*)$ with probability $1$, so $V_i (c) = sum_(j in N(i)) bb(I)[ell_j^(*) = c]$ equals the exact neighborhood label count of asynchronous LPA @raghavan2007lpa. Selecting $ell_i^(*) in arg max_c V_i (c)$ with uniform random tie-breaking matches the LPA transition kernel identically, and post-processing any $r in (0, 1]$ returns the singleton ${ell_i^((T))}$. Conversely, setting $W = infinity$ retains every appended label, recovering Definition 1 identically.
]

Unlike Fluid Communities @pares2018fluidc—which enforces a fixed community count $K$ and density conservation to prevent giant-cluster collapse—and unlike Walktrap @pons2006walktrap or MCL @vandongen2000mcl—which rely on Markov random-walk similarity matrices—SLPA achieves linear edge complexity via discrete majority filtering. Compared to OSLOM @lancichinetti2011oslom and BigCLAM @yang2013bigclam, SLPA trades formal probabilistic/statistical guarantees for superior raw speed and zero requirement to specify $K$.

*Empirical Behavior on ClusterOverview Benchmarks.* On the 4-dataset ClusterOverview benchmark suite (evaluated on $k$-NN graphs via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), SLPA records *Iris* Sil $= 0.341$, ARI $= 0.319$; *Customers* Sil $= 0.382$; *Multi-Shape* Sil $= 0.510$ (ARI N/A); and *High-D* Sil $= 0.303$ (ARI N/A), with the fastest mean runtime of all four algorithms in this section ($2.628$ ms). This performance profile aligns tightly with our theoretical findings. The $2.628$ ms mean runtime directly reflects SLPA's $cal(O)(T m)$ complexity and $cal(O)(1)$ array speaker sampling. Remarkably, SLPA achieves the *highest Silhouette* among Algs. 123–126 on both *Iris* ($0.341$) and *High-D* ($0.303$): in high-dimensional $k$-NN graphs where distance concentration flattens global configuration-null contrast (reducing OSLOM to $0.175$) and sparse neighborhoods starve clique factorization (reducing BigCLAM to $0.065$), SLPA's purely local majority-voting memory preserves compact neighborhood consensus. Conversely, on *Iris* its ARI ($0.319$) trails OSLOM ($0.451$) because majority amplification across the overlapping Versicolor-Virginica boundary can merge adjacent classes (Theorem 1(3)), and on *Multi-Shape* ($0.510$ vs OSLOM's $0.780$) random propagation schedules occasionally split elongated thin manifolds before labels traverse their full geodesic diameter within finite $T$.
