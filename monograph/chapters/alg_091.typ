#import "../template/ieee.typ": *

#alg_entry(
  id: 91,
  name: "Hierarchical Dirichlet Process (HDP)",
  category: "Cat. IV: Distribution, Statistical & Topic",
  year: "2006",
  paradigm: "Nonparametric Bayesian Hierarchical Mixture",
  complexity: $cal(O)(T N K)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Parametric mixture models and Latent Dirichlet Allocation (LDA)
@blei2003lda require specifying the cluster or topic count $K$ in advance,
forcing costly cross-validation sweeps across candidate cardinalities.
Conversely, assigning an independent Dirichlet Process (DP) mixture
@rasmussen2000igmm to each group of a grouped dataset (such as multiple
documents, patients, or experimental batches) prevents distinct groups
from sharing mixture components whenever the base measure $H$ is
non-atomic. The *Hierarchical Dirichlet Process (HDP)* introduced by Teh
et al. @teh2006hdp resolves both obstacles simultaneously by coupling
group-level Dirichlet Processes through a shared random discrete base
measure $G_0$ drawn from a top-level Dirichlet Process.

#definition(title: "Hierarchical Dirichlet Process and CRF Representation")[
  Let $(Theta, cal(B))$ be a measurable parameter space equipped with a
  probability base measure $H$, concentration parameters $gamma > 0$ and
  $alpha_0 > 0$, and $J$ data groups
  $bold(x)_j = (x_(j 1), dots, x_(j n_j))$ totaling $N = sum_(j=1)^J n_j$
  observations. The *Hierarchical Dirichlet Process* mixture model is
  defined hierarchically by
  $ G_0 mid(|) gamma, H &tilde.op "DP"(gamma, H), \
    G_j mid(|) alpha_0, G_0 &tilde.op "DP"(alpha_0, G_0) quad (j = 1, dots, J), \
    theta_(j i) mid(|) G_j &tilde.op G_j, quad x_(j i) mid(|) theta_(j i) tilde.op F(theta_(j i)). $
  Under Sethuraman's stick-breaking construction,
  $G_0 = sum_(k=1)^infinity beta_k delta_(phi_k)$ with
  $phi_k tilde.op H$ i.i.d. and $bold(beta) tilde.op "GEM"(gamma)$, while
  each group measure shares the exact atom set ${phi_k}_(k=1)^infinity$
  via $G_j = sum_(k=1)^infinity pi_(j k) delta_(phi_k)$ with
  $bold(pi)_j mid(|) alpha_0, bold(beta) tilde.op "DP"(alpha_0, bold(beta))$.
  Marginalizing $G_0$ and ${G_j}$ induces the *Chinese Restaurant
  Franchise (CRF)*: each customer $(j, i)$ selects table $t_(j i)$ in
  restaurant $j$ with weight $n_(j t)$ (existing table) or $alpha_0$ (new
  table), and each table $(j, t)$ selects global dish $k_(j t)$ with
  weight $m_(dot k)$ (tables serving dish $k$ across the franchise) or
  $gamma$ (new dish $phi_(K+1) tilde.op H$).
]

Integrating out $G_j$ while explicitly maintaining the global weights
$bold(beta)$ yields the *direct-assignment* collapsed Gibbs sampler
@teh2006hdp, which tracks latent atom assignments $z_(j i) in {1, dots, K}$
(for $K$ active global atoms) and auxiliary table counts $m_(j k)$
representing the number of CRF tables in group $j$ assigned to atom $k$.
For conjugate likelihood-prior pairs $(F, H)$ (such as
Multinomial-Dirichlet for documents or Normal-Inverse-Wishart for
continuous vectors), the marginal likelihood
$f_k^(-j i) (x_(j i)) = integral F(x_(j i) mid(|) phi) p(phi mid(|) cal(X)_k^(-j i), H) dif phi$
admits a closed-form expression, allowing rapid collapsed MCMC or
truncated mean-field variational inference @blei2006variational.

#algbox(
  title: "Direct-Assignment Collapsed Gibbs Sampler for HDP",
  inputs: [Grouped data ${x_(j i)}$ ($j = 1, dots, J$, $i = 1, dots, n_j$); conjugate pair $(F, H)$; concentrations $gamma, alpha_0 > 0$; sweeps $T$.],
  outputs: [Active atom count $K$; posterior assignments ${z_(j i)}$; global weights $bold(beta)$ and group proportions ${bold(pi)_j}$.],
)[
  + *Initialize:* Assign all observations to $K = 1$ atom ($z_(j i) = 1$), set table counts $m_(j 1) = 1$, and draw $bold(beta) = (beta_1, beta_u) tilde.op "Dir"(J, gamma)$.
  + *For* $t = 1, dots, T$ *do:*
  + #h(1em) *Sample assignments $z_(j i)$:* For each $(j, i)$, decrement $n_(j z_(j i))$ and compute conditional probabilities for $k = 1, dots, K$ and new atom $k_"new" = K + 1$:
    $ "Pr"(z_(j i) = k mid(|) dots) &prop (n_(j k)^(-j i) + alpha_0 beta_k) f_k^(-j i) (x_(j i)), \
      "Pr"(z_(j i) = k_"new" mid(|) dots) &prop alpha_0 beta_u f_0 (x_(j i)). $
    Draw $z_(j i)$; if $z_(j i) = k_"new"$, split $beta_u$ via $b tilde.op "Beta"(1, gamma)$, set $beta_(K+1) = b beta_u$, $beta_u <- (1 - b) beta_u$, and increment $K$. Remove any empty atom and absorb its weight into $beta_u$.
  + #h(1em) *Sample auxiliary table counts $m_(j k)$:* For each $(j, k)$ with $n_(j k) > 0$, set $m_(j k) = 0$ and for $r = 1, dots, n_(j k)$ increment $m_(j k)$ with probability $alpha_0 beta_k / (alpha_0 beta_k + r - 1)$.
  + #h(1em) *Resample global weights $bold(beta)$:* Let $m_(dot k) = sum_(j=1)^J m_(j k)$; draw
    $ (beta_1, dots, beta_K, beta_u) tilde.op "Dir"(m_(dot 1), dots, m_(dot K), gamma). $
  + *Return* posterior MAP/consensus partition ${z_(j i)}$, $K$, and $bold(beta)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of the HDP: exact
cross-group correlation and almost-sure atom sharing, the coupled
stick-breaking construction, an original double-logarithmic bound on
expected global cluster growth, an original table-sharpened total
variation truncation bound for variational HDP, and exact parametric and
flat-DP limits.

#subsec("A", "Almost-Sure Atom Sharing and Cross-Group Correlation")

#proposition(title: "Atom Sharing and Exact Cross-Group Measure Correlation", origin: "Original Derivation")[
  Let $(G_0, G_1, dots, G_J)$ follow an HDP with concentrations
  $gamma, alpha_0 > 0$ and base probability measure $H$. Then:
  (i) $"Pr"("supp"(G_j) subset.eq "supp"(G_0)) = 1$ for every $j$, and if
  $H$ is non-atomic, distinct groups share at least one drawn parameter
  with strictly positive probability
  $"Pr"(theta_(j 1) = theta_(j' 1)) = 1 / (gamma + 1) > 0$ ($j != j'$),
  whereas independent DPs satisfy $"Pr"(theta_(j 1) = theta_(j' 1)) = 0$;
  (ii) for every measurable set $A in cal(B)$ with $0 < H(A) < 1$ and
  every $j != j'$,
  $ bb(E)[G_j (A)] &= H(A), \
    "Var"(G_j (A)) &= H(A)(1 - H(A)) frac(alpha_0 + gamma + 1, (alpha_0 + 1)(gamma + 1)), \
    "Cov"(G_j (A), G_(j') (A)) &= frac(H(A)(1 - H(A)), gamma + 1), $
  yielding the set-independent inter-group correlation
  $ "Corr"(G_j (A), G_(j') (A)) = frac(alpha_0 + 1, alpha_0 + gamma + 1). $
]

#proof[
  Since $G_0 tilde.op "DP"(gamma, H)$ is almost surely purely atomic
  ($G_0 = sum_(k=1)^infinity beta_k delta_(phi_k)$), any draw
  $G_j mid(|) G_0 tilde.op "DP"(alpha_0, G_0)$ places probability $1$ on
  $"supp"(G_0)$, proving the first claim of (i). Conditioned on $G_0$,
  $theta_(j 1)$ and $theta_(j' 1)$ ($j != j'$) are independent draws from
  $G_0$ (since $bb(E)[G_j (dot) mid(|) G_0] = G_0 (dot)$), so
  $ "Pr"(theta_(j 1) = theta_(j' 1)) = bb(E)[sum_(k=1)^infinity beta_k^2] = frac(1, gamma + 1), $
  using the standard GEM second-moment identity
  $bb(E)[sum_k beta_k^2] = bb(E)[(beta'_1)^2] / (1 - bb(E)[(1 - beta'_1)^2]) = (2 / ((gamma + 1)(gamma + 2))) / (1 - gamma / (gamma + 2)) = 1 / (gamma + 1)$.

  For (ii), set $p = H(A)$. By DP marginal properties,
  $G_0 (A) tilde.op "Beta"(gamma p, gamma (1 - p))$ has mean $p$ and
  variance $p(1 - p) / (gamma + 1)$, while given $G_0$,
  $G_j (A) mid(|) G_0 tilde.op "Beta"(alpha_0 G_0 (A), alpha_0 (1 - G_0 (A)))$
  independently across $j$. Hence
  $bb(E)[G_j (A)] = bb(E)[G_0 (A)] = p$, and
  $ bb(E)[G_0 (A)(1 - G_0 (A))] &= p - p^2 - "Var"(G_0 (A)) \
    &= p(1 - p) frac(gamma, gamma + 1). $
  Applying the law of total variance gives
  $ "Var"(G_j (A)) &= bb(E)[frac(G_0 (A)(1 - G_0 (A)), alpha_0 + 1)] + "Var"(G_0 (A)) \
    &= p(1 - p) (frac(gamma, (alpha_0 + 1)(gamma + 1)) + frac(1, gamma + 1)) \
    &= p(1 - p) frac(alpha_0 + gamma + 1, (alpha_0 + 1)(gamma + 1)). $
  For $j != j'$, conditional independence given $G_0$ kills
  $bb(E)["Cov"(G_j (A), G_(j') (A) mid(|) G_0)]$, leaving
  $"Cov"(G_j (A), G_(j') (A)) = "Var"(G_0 (A)) = p(1 - p) / (gamma + 1)$.
  Dividing covariance by variance cancels $p(1 - p) / (gamma + 1)$ and
  yields $(alpha_0 + 1) / (alpha_0 + gamma + 1)$.
]

#subsec("B", "Coupled Stick-Breaking Construction of Group Measures")

#theorem(title: "Coupled Stick-Breaking Representation of Group Weights", origin: "Classical — Cited")[
  Let $bold(beta) tilde.op "GEM"(gamma)$ and define independent stick
  fractions
  $ tilde(pi)_(j k) mid(|) bold(beta) tilde.op "Beta"(alpha_0 beta_k, alpha_0 (1 - sum_(l=1)^k beta_l)) quad (k >= 1), $
  with $pi_(j k) := tilde(pi)_(j k) product_(l=1)^(k-1) (1 - tilde(pi)_(j l))$
  @teh2006hdp. Then $sum_(k=1)^infinity pi_(j k) = 1$ almost surely and
  $G_j = sum_(k=1)^infinity pi_(j k) delta_(phi_k)$ satisfies
  $G_j mid(|) alpha_0, G_0 tilde.op "DP"(alpha_0, G_0)$.
]

#proof[
  Fix $G_0 = sum_(k=1)^infinity beta_k delta_(phi_k)$ and any integer
  $K >= 1$. Define the tail masses $B_K := 1 - sum_(l=1)^K beta_l$ and
  $R_(j K) := 1 - sum_(l=1)^K pi_(j l) = product_(l=1)^K (1 - tilde(pi)_(j l))$.
  We verify by induction on $K$ that
  $(pi_(j 1), dots, pi_(j K), R_(j K)) tilde.op "Dir"(alpha_0 beta_1, dots, alpha_0 beta_K, alpha_0 B_K)$.
  For $K = 1$, $pi_(j 1) = tilde(pi)_(j 1) tilde.op "Beta"(alpha_0 beta_1, alpha_0 B_1)$,
  which is $"Dir"(alpha_0 beta_1, alpha_0 B_1)$. Assuming the statement for
  $K - 1$, we have
  $(pi_(j K), R_(j K)) = R_(j, K-1) (tilde(pi)_(j K), 1 - tilde(pi)_(j K))$
  where $tilde(pi)_(j K) tilde.op "Beta"(alpha_0 beta_K, alpha_0 B_K)$ is
  independent of $(pi_(j 1), dots, pi_(j, K-1), R_(j, K-1))$ and
  $alpha_0 beta_K + alpha_0 B_K = alpha_0 B_(K-1)$. By the neutrality
  (splitting) property of the Dirichlet distribution, splitting the last
  component $R_(j, K-1)$ by an independent
  $"Beta"(alpha_0 beta_K, alpha_0 B_K)$ yields
  $"Dir"(alpha_0 beta_1, dots, alpha_0 beta_K, alpha_0 B_K)$. Moreover,
  $bb(E)[R_(j K) mid(|) bold(beta)] = B_K -> 0$ as $K -> infinity$, so
  $R_(j K) -> 0$ almost surely by monotone convergence. By the Dirichlet
  aggregation property, for every finite measurable partition
  $(A_1, dots, A_r)$ of $Theta$, $(G_j (A_1), dots, G_j (A_r)) tilde.op "Dir"(alpha_0 G_0 (A_1), dots, alpha_0 G_0 (A_r))$,
  which is the defining property of $"DP"(alpha_0, G_0)$.
]

#subsec("C", "Double-Logarithmic Growth of Global Cluster Count")

Whereas a flat Dirichlet Process generates $cal(O)(gamma log N)$ clusters
from $N$ samples, coupling groups through a two-level hierarchy sharply
compresses the growth rate of distinct global clusters.

#theorem(title: "Double-Logarithmic Upper Bound on Expected Global HDP Clusters", origin: "Original Proof")[
  Let $K_N$ be the number of distinct global atoms occupied by
  $N = sum_(j=1)^J n_j$ observations ($n_j >= 1$) drawn from an HDP with
  parameters $(gamma, alpha_0)$ and non-atomic $H$. Then:
  (i) the expected table count $M_j$ in group $j$ satisfies
  $ alpha_0 log(1 + frac(n_j, alpha_0)) <= bb(E)[M_j] <= 1 + alpha_0 log(1 + frac(n_j - 1, alpha_0)); $
  (ii) defining $bar(M) := J + alpha_0 sum_(j=1)^J log(1 + (n_j - 1) / alpha_0)$,
  the expected number of distinct global clusters obeys the exact
  non-asymptotic bound
  $ bb(E)[K_N] <= 1 + gamma log(1 + frac(bar(M) - 1, gamma)). $
  In particular, for fixed $J$ and balanced group sizes $n_j = N / J -> infinity$,
  $bb(E)[K_N] = cal(O)(gamma log(J alpha_0 log(N / J)))$ grows
  double-logarithmically in $N$.
]

#proof[
  In the CRF representation of group $j$, customer $i in {1, dots, n_j}$
  opens a new table with probability $alpha_0 / (alpha_0 + i - 1)$
  independently of earlier table counts, so
  $ bb(E)[M_j] = sum_(r=0)^(n_j - 1) frac(alpha_0, alpha_0 + r). $
  Because $x |-> alpha_0 / (alpha_0 + x)$ is strictly decreasing on
  $[0, infinity)$, integral comparison yields
  $ integral_0^(n_j) frac(alpha_0, alpha_0 + x) dif x &<= sum_(r=0)^(n_j - 1) frac(alpha_0, alpha_0 + r) \
    &<= 1 + integral_0^(n_j - 1) frac(alpha_0, alpha_0 + x) dif x, $
  evaluating directly to (i).

  Let $M = sum_(j=1)^J M_j >= J >= 1$ be the total number of tables across
  all groups. Conditioned on $M$, the $M$ tables draw their global dishes
  from a top-level CRP with concentration $gamma$, so applying the same
  upper integral bound conditionally gives
  $ bb(E)[K_N mid(|) M] = sum_(r=0)^(M-1) frac(gamma, gamma + r) <= 1 + gamma log(1 + frac(M - 1, gamma)) =: g(M). $
  Observe that $g'(m) = gamma / (gamma + m - 1) > 0$ and
  $g''(m) = -gamma / (gamma + m - 1)^2 < 0$ on $[1, infinity)$, so $g$ is
  strictly increasing and strictly concave. Applying the tower property,
  Jensen's inequality, and part (i) yields
  $ bb(E)[K_N] = bb(E)[bb(E)[K_N mid(|) M]] <= bb(E)[g(M)] <= g(bb(E)[M]), $
  and substituting $bb(E)[M] = sum_(j=1)^J bb(E)[M_j] <= J + alpha_0 sum_(j=1)^J log(1 + (n_j - 1) / alpha_0)$
  produces (ii).
]

#subsec("D", "Table-Sharpened Truncation Bound for Variational HDP")

Practical variational inference @blei2006variational truncates $G_0$ to
$T$ global atoms by merging the residual stick mass
$R_T := 1 - sum_(k=1)^T beta_k = product_(k=1)^T (1 - beta'_k)$ into atom
$T$. While naive sample-level coupling incurs a factor of $N$, the CRF
table bottleneck yields an exponentially sharper bound in $N / J$.

#proposition(title: "CRF-Sharpened Total Variation Truncation Bound", origin: "Original Derivation")[
  Let $p_infinity$ denote the marginal joint density of all $N$ grouped
  observations under the infinite HDP, and let $p_T$ denote the density
  obtained by truncating $G_0$ at $T$ sticks ($G_0^(T) = sum_(k=1)^(T-1) beta_k delta_(phi_k) + (beta_T + R_T) delta_(phi_T)$).
  Then the total variation distance satisfies
  $ frac(1, 2) norm(p_infinity - p_T)_1 &<= (J + alpha_0 sum_(j=1)^J log(1 + frac(n_j - 1, alpha_0))) \
    &quad times (frac(gamma, gamma + 1))^T. $
  Consequently, achieving total variation error at most $epsilon in (0, 1)$
  requires only
  $T = cal(O)(gamma log((J alpha_0 log(N / J)) / epsilon))$ global sticks.
]

#proof[
  Couple the infinite and truncated processes by sharing the same CRF
  table partition (total table count $M = sum_(j=1)^J M_j$) and the same
  top-level stick lengths $(beta_1, beta_2, dots)$ and atoms
  $(phi_1, phi_2, dots)$. Whenever every one of the $M$ tables draws a
  dish index in ${1, dots, T}$, the realized parameter vector
  ${theta_(j i)}$ and observation likelihood under $p_infinity$ and $p_T$
  are identical. By the coupling inequality for total variation,
  $ frac(1, 2) norm(p_infinity - p_T)_1 <= "Pr"(exists t in {1, dots, M} : k_t > T). $
  Conditioned on $(M, G_0)$, the $M$ table dishes are i.i.d. draws from
  $G_0$, each landing beyond index $T$ with probability $R_T$, so by the
  union bound
  $"Pr"(exists t <= M : k_t > T mid(|) M, G_0) <= M R_T$.
  Crucially, in the CRF, table creation depends solely on $(alpha_0, n_1, dots, n_J)$
  and is statistically independent of $G_0$. Therefore
  $ bb(E)[M R_T] = bb(E)[M] dot bb(E)[R_T]. $
  Since $1 - beta'_k tilde.op "Beta"(gamma, 1)$ are i.i.d. with mean
  $gamma / (gamma + 1)$, we have $bb(E)[R_T] = (gamma / (gamma + 1))^T$.
  Combining this with Theorem 2(i) yields the stated bound; setting the
  right-hand side $<= epsilon$ and solving for $T$ gives
  $T >= log(bb(E)[M] / epsilon) / log(1 + 1 / gamma) = cal(O)(gamma log(bb(E)[M] / epsilon))$.
]

#subsec("E", "Exact Limits to LDA and Flat Dirichlet Processes")

#proposition(title: "Parametric LDA and Degenerate-Hierarchy Limits of HDP", origin: "Classical — Cited")[
  The HDP strictly encompasses three canonical models as exact limits
  @teh2006hdp @blei2003lda @rasmussen2000igmm:
  (i) *Finite symmetric limit:* replacing $G_0$ by $bold(beta) tilde.op "Dir"(gamma / L, dots, gamma / L)$
  converges weakly in partition law to HDP as $L -> infinity$, whereas
  fixing $bold(beta) = (1 / L, dots, 1 / L)$ recovers symmetric $L$-topic
  LDA with Dirichlet prior $alpha = alpha_0 / L$;
  (ii) *Pooled flat DP ($alpha_0 -> infinity$):* $G_j -> G_0$ in
  $L^2$-measure for every measurable $A$, collapsing all $J$ groups into
  a single DP mixture;
  (iii) *Decoupled independent DPs ($gamma -> infinity$):* $G_0 -> H$
  weakly and $"Corr"(G_j (A), G_(j') (A)) -> 0$, eliminating cross-group
  sharing.
]

#proof[
  For (i), under $bold(beta) = (1 / L, dots, 1 / L)$, the group mixing
  distribution $bold(pi)_j mid(|) alpha_0, bold(beta) tilde.op "Dir"(alpha_0 / L, dots, alpha_0 / L)$
  coincides identically with the symmetric document-topic prior of LDA
  @blei2003lda, while the weak convergence of symmetric $L$-component
  Dirichlet-Dirichlet hierarchies to HDP follows from Kingman's
  Poisson-Dirichlet limit @teh2006hdp. For (ii), Proposition 1 gives
  $bb(E)[(G_j (A) - G_0 (A))^2] = bb(E)["Var"(G_j (A) mid(|) G_0)] = H(A)(1 - H(A)) gamma / ((alpha_0 + 1)(gamma + 1)) -> 0$
  as $alpha_0 -> infinity$. For (iii), $"Var"(G_0 (A)) = H(A)(1 - H(A)) / (gamma + 1) -> 0$
  as $gamma -> infinity$, so $G_0 (A) -> H(A)$ in probability and
  $"Corr"(G_j (A), G_(j') (A)) = (alpha_0 + 1) / (alpha_0 + gamma + 1) -> 0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In direct-assignment collapsed Gibbs sampling, each sweep over $N$
observations evaluates $K + 1$ candidate atoms in $cal(O)(N K)$ scalar
likelihood calls (multiplied by $d$ for continuous $d$-dimensional
conjugate Gaussian atoms, or $cal(O)(1)$ per token for multinomial text
counts). Sampling auxiliary table counts $m_(j k)$ via Antoniak Bernoulli
trials visits $N$ total increments across all $(j, k)$ pairs, costing
$cal(O)(N)$ per sweep; resampling $bold(beta)$ costs $cal(O)(K)$. Thus
one full Gibbs sweep runs in $cal(O)(N K)$ time for discrete/low-$d$
conjugate observations (and $cal(O)(N K d)$ for diagonal-covariance
continuous vectors, or $cal(O)(N K d^2)$ for full Normal-Inverse-Wishart
atoms). Truncated stochastic variational inference @blei2006variational
replaces Markov-chain autocorrelation with deterministic coordinate ascent
in $cal(O)(T_"var" N T_"trunc")$, where Proposition 2 guarantees small
$T_"trunc"$ even for massive $N$.

Because $K$ is random rather than fixed, HDP adapts partition complexity
automatically to the data geometry: multimodal or non-convex shapes can
be approximated by unions of locally coherent atoms without prescribing
$K$. However, posterior inference depends critically on the concentration
hyperparameters $(gamma, alpha_0)$, which are typically endowed with
vague Gamma hyperpriors and resampled via auxiliary Beta variables
@teh2006hdp.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time complexity], [$cal(O)(T N K)$ discrete / $cal(O)(T N K d)$ diagonal continuous],
    [Space complexity], [$cal(O)(N + J K + K d)$ for sufficient statistics & tables],
    [Cluster count $K$], [Learned automatically; $bb(E)[K_N] = cal(O)(gamma log(J alpha_0 log(N/J)))$],
    [Statistical sharing], [Exact cross-group correlation $(alpha_0 + 1) / (alpha_0 + gamma + 1)$],
    [Conjugacy], [Closed-form collapsed updates for exponential families],
    [Inference options], [CRF / direct Gibbs MCMC; truncated variational Bayes],
    [Permutation / scale], [Exchangeable within/across groups; scale follows base $F$],
    [Hyperparameters], [Top $gamma$, group $alpha_0$, base prior $H$ (Gamma hyperpriors)],
  ),
  caption: [Operational and statistical profile of the Hierarchical Dirichlet Process.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Multi-Corpus and Streaming Topic Modeling:* Large document collections
  spanning evolving news archives or scientific repositories continuously
  introduce new latent themes while sharing core vocabulary across
  documents @teh2006hdp @blei2003lda. HDP eliminates fixed-$K$ topic
  selection, borrowing statistical strength across short documents via
  the shared global measure $G_0$ while spawning new topics on demand.
+ *Multi-Center Biomedical Cohort Clustering:* Patient records collected
  across distinct hospitals exhibit site-specific prevalence shifts
  ($G_j$) over shared underlying biological disease subtypes ($G_0$).
  HDP explicitly separates shared cluster parameters ${phi_k}$ from
  hospital-specific mixture weights $bold(pi)_j$, avoiding the confounding
  of batch prevalence with cluster identity.
+ *Speaker Diarization and Acoustic Segmentation (HDP-HMM):* Extending
  HDP so that each hidden state $j$ owns a transition distribution
  $G_j tilde.op "DP"(alpha_0, G_0)$ yields the infinite Hidden Markov
  Model @teh2006hdp, automatically inferring an unknown number of speakers
  in audio streams while coupling transition targets across states.
+ *Information Retrieval and Multi-View Image Scene Parsing:* Images
  represented as bags of local visual patches share recurring object
  parts across scenes ($G_0$) in scene-dependent proportions ($G_j$).
  Nonparametric atom sharing lets rare visual categories benefit from
  global parameter estimation without distorting dominant scene mixtures.
+ *Population Genetics and Haplotype Sharing:* Individuals sampled across
  geographically structured subpopulations share ancestral haplotype
  blocks in varying proportions. The CRF genealogy naturally models
  recurrent migration and founder effects with double-logarithmic global
  allele discovery (Theorem 2).

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four operational pitfalls govern HDP in practice:
(1) *Hyperparameter sensitivity and over-splitting:* On small, noisy
datasets, misspecified $(gamma, alpha_0)$ or narrow base-measure scale
$H$ can fracture unimodal overlapping classes into multiple redundant
atoms, depressing ARI despite strong local cohesion.
(2) *MCMC mixing and label switching:* Collapsed Gibbs samplers can mix
slowly between split and merged modes on high-dimensional continuous data,
requiring split-merge Metropolis-Hastings moves or variational Bayes
@blei2006variational.
(3) *Non-conjugate likelihood overhead:* Outside exponential families
with conjugate $H$, evaluating $f_k^(-j i)$ requires auxiliary-variable
MCMC or Laplace approximations, inflating runtime substantially.
(4) *Single-group degeneracy:* When data exhibit no genuine grouping
structure ($J = 1$), the second hierarchy level is superfluous and adds
variance relative to a flat DP mixture @rasmussen2000igmm.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Multi-document text corpora], [Adaptive topic count; exact cross-doc sharing], [For very short texts, consider BTM / DMM],
    [Multi-site clinical cohorts], [Decouples site prevalence from shared subtypes], [Place Gamma hyperpriors on $gamma, alpha_0$],
    [Sequential audio / speech], [HDP-HMM learns unknown state cardinality], [Add sticky self-transition bias for stability],
    [Small single-group tabular], [Extra hierarchy level can over-split classes], [Use flat DP-GMM, BIC-GMM, or LCA (Alg. 90)],
    [Ultra-low-latency streaming], [Collapsed MCMC sweeps too slow per batch], [Use online variational HDP or spherical $k$-means],
  ),
  caption: [Domain suitability and failure-mode matrix for HDP.],
)

#sec("V", "Structural Relations and Empirical Behavior")

HDP occupies the nonparametric apex of the probabilistic topic and
mixture family in Cat. IV. By Proposition 3, it strictly subsumes
symmetric LDA @blei2003lda and flat infinite mixture models
@rasmussen2000igmm, while standing in structural contrast to Correlated
Topic Models @blei2007ctm (which model topic covariance via logistic-normal
priors at fixed $K$), pLSA @hofmann1999plsa (which lacks document-level
generative priors), and short-text single-cluster models such as BTM
@yan2013btm and DMM @yin2014dmm. Unlike NMF (Alg. 92), which performs
deterministic low-rank factorization without uncertainty quantification,
HDP delivers full posterior distributions over both cluster cardinality
and component parameters.

In the ClusterOverview 4-dataset benchmark (evaluated via Silhouette
@rousseeuw1987silhouettes and ARI @hubert1985comparing), HDP incurs the
highest computational cost among Algs. 89--92 (mean runtime *5.778 ms*,
reflecting iterative stochastic sampling and dynamic atom management),
yet delivers the *strongest Silhouette scores across three of the four
benchmarks*:
- *Multi-Shape:* Silhouette $= 0.782$ ($"ARI" =$ N/A) --- the highest
  score in this group (exceeding movMF $0.709$, LCA $0.709$, and NMF
  $0.710$). Allowing the DP hierarchy to instantiate an adaptive number of
  local components lets HDP conform closely to multi-modal geometric
  substructures rather than forcing rigid fixed-$K$ partitions.
- *Customers:* Silhouette $= 0.380$ --- again leading all four models
  (vs. movMF $0.351$, NMF $0.266$, LCA $0.182$), as nonparametric atom
  creation cleanly isolates skewed heavy-tail wholesale spenders without
  corrupting core cluster centroids.
- *High-D:* Silhouette $= 0.320$ ($"ARI" =$ N/A) --- leading the group
  by regularizing component estimates through the conjugate base prior $H$
  in high ambient dimensions.
- *Iris:* Silhouette $= 0.267$, $"ARI" = 0.466$. Here the classic
  nonparametric trade-off surfaces clearly: on a small ($N = 150$)
  three-class dataset with overlapping _Versicolor_-_Virginica_ margins,
  the infinite prior over-splits or merges borderline modes relative to
  oracle-$K$ parametric mixtures, lowering both Iris Silhouette and ARI.

These empirical diagnostics corroborate the theory directly: HDP excels
whenever cluster multiplicity is unknown, heavy-tailed, or hierarchically
structured across groups @teh2006hdp @blei2006variational, while fixed-$K$
parametric models remain sharper on small datasets with known cardinality.
