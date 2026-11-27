#import "../template/ieee.typ": *

#alg_entry(
  id: 200,
  name: "ACO-Clustering (Ant Colony Optimization)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2004",
  paradigm: "Pheromone-Guided Combinatorial Partitioning",
  complexity: $O(T_"max" n (R (k + d) + L I_"ls" k d))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Given $cal(X) = {x_1, dots, x_n} subset bb(R)^d$ and $k in {2, dots, n - 1}$ (with at least $k$ distinct points), hard partitional clustering seeks a surjective assignment string $s = (s_1, dots, s_n) in Omega_k subset.eq {1, dots, k}^n$, inducing non-empty clusters $C_j (s) = {i : s_i = j}$ and centroids $mu_j (s) = abs(C_j (s))^(-1) sum_(i in C_j (s)) x_i$, that minimizes the within-cluster sum of squared errors (SSE):
$ J(s) = sum_(j=1)^k sum_(i in C_j (s)) norm(x_i - mu_j (s))_2^2. $
Minimizing $J(s)$ over $Omega_k$ is NP-hard even for $k = 2$ or $d = 2$ @aloise2009np, and deterministic refinement such as Lloyd's algorithm @lloyd1982least can take exponentially many steps @vattani2011kmeans while stalling in suboptimal Voronoi basins. Ant Colony Optimization for clustering (ACO-Clustering) @shelokar2004aco performs stochastic global search via stigmergic positive feedback: a colony of $R$ artificial ants samples candidate strings from a shared pheromone trail matrix $T^((t)) in bb(R)_+^(n times k)$, refines elite strings with local search, and reinforces trails inversely with SSE @jain2010data @xu2005survey.

Two ant-inspired clustering paradigms exist. In *spatial ant-sorting*, mobile agents traverse a 2D toroidal grid and stochastically pick up or drop objects based on local neighborhood similarity $f(i) in [0, 1]$, via $p_"pick" (i) = (k_1 / (k_1 + f(i)))^2$ and $p_"drop" (i) = (f(i) / (k_2 + f(i)))^2$, forming spatial heaps without fixing $k$ upfront; however, grid sorting lacks an explicit objective and mixes slowly. The *pheromone-matrix formulation* @shelokar2004aco optimizes $J(s)$ directly over $Omega_k$, couples cleanly with local refinement, and admits the exact convergence and stability analysis below.

#definition(title: "Bounded Elitist ACO-Clustering System")[
  An *ACO-Clustering system* on $cal(X)$ is the stochastic process ${(T^((t)), s^(ast (t)))}_(t >= 0)$ with colony size $R >= 1$, elite count $L in {1, dots, R}$, evaporation $rho in (0, 1)$, exploitation threshold $q_0 in [0, 1)$, exponents $alpha > 0, beta >= 0$, deposit scale $Q > 0$, stabilizer $epsilon > 0$, and bounds $0 <= tau_"min" < tau_"max" <= infinity$:
  1. *Pheromone Matrix*: $T^((t)) = [tau_(i j)^((t))] in [tau_"min", tau_"max"]^(n times k)$, initialized to $tau_(i j)^((0)) = tau_0 in [tau_"min", tau_"max"]$.
  2. *Transition Kernel*: At step $t$, each ant $r <= R$ independently assigns every $i <= n$ by drawing $u ~ "Uniform"(0, 1)$. Given visibility $eta_(i j)^((t)) > 0$ ($eta_(i j)^((t)) = 1$ when $beta = 0$ @shelokar2004aco, or $(1 + norm(x_i - hat(mu)_j^((t)))_2^2)^(-1)$ from incumbent centroids), set $s_i^((r)) = arg max_j (tau_(i j)^((t)))^alpha (eta_(i j)^((t)))^beta$ if $u <= q_0$; else sample $s_i^((r)) = j$ with
     $ p_(i j)^((t)) = ((tau_(i j)^((t)))^alpha (eta_(i j)^((t)))^beta) / (sum_(ell=1)^k (tau_(i ell)^((t)))^alpha (eta_(i ell)^((t)))^beta). $
     Any empty cluster is repaired deterministically by relocating the point farthest from its centroid so $s^((r)) in Omega_k$.
  3. *Local Search*: The $L$ lowest-SSE ants undergo monotone refinement $tilde(s)^((ell)) = cal(L)(s^((ell))) in Omega_k$ with $J(tilde(s)^((ell))) <= J(s^((ell)))$.
  4. *Evaporation and Deposit*: With $Pi_([a, b]) (z) = min(b, max(a, z))$,
     $ tau_(i j)^((t+1)) &= Pi_([tau_"min", tau_"max"]) ((1 - rho) tau_(i j)^((t)) + sum_(ell=1)^L Delta tau_(i j)^((ell, t))), \
       Delta tau_(i j)^((ell, t)) &= Q / (J(tilde(s)^((ell))) + epsilon) bb(I)[tilde(s)_i^((ell)) = j]. $
]

Setting $tau_"min" = 0, tau_"max" = infinity, alpha = 1, beta = 0$ recovers the unclipped update of @shelokar2004aco; enforcing $0 < tau_"min" < tau_"max" < infinity$ acts as a max-min safeguard ensuring persistent exploration.

#algbox(
  title: "ACO-Clustering (Bounded Elitist Formulation)",
  inputs: [$cal(X) = {x_1, dots, x_n} subset bb(R)^d$, $k$, $R$, $L$, $(rho, q_0, alpha, beta, Q, epsilon, tau_0, tau_"min", tau_"max", T_"max")$.],
  outputs: [Incumbent partition $s^ast in Omega_k$ and centroids ${mu_j (s^ast)}_(j=1)^k$.],
)[
  1. Set $tau_(i j)^((0)) <- tau_0$, $J^ast <- infinity$, $s^(ast (0)) <- emptyset$; if $beta > 0$, seed $hat(mu)_j^((0))$ via $k$-means++ @arthur2007kmeanspp.
  2. *For* $t = 0, 1, dots, T_"max" - 1$:
  3. #h(1em) *For* $r = 1, dots, R$ (parallel): sample $s^((r))$ via $(T^((t)), eta^((t)))$, repair empty clusters, and compute $mu_j (s^((r)))$ and $J(s^((r)))$.
  4. #h(1em) Retain the $L$ lowest-SSE strings and apply local search $tilde(s)^((ell)) <- cal(L)(s^((ell)))$ (up to $I_"ls"$ Lloyd / Hartigan-Wong sweeps).
  5. #h(1em) Let $ell^ast = arg min_(ell <= L) J(tilde(s)^((ell)))$. If $J(tilde(s)^((ell^ast))) < J^ast$, update $J^ast$, $s^(ast (t+1)) <- tilde(s)^((ell^ast))$, and $hat(mu)_j^((t+1))$; else keep $s^(ast (t))$.
  6. #h(1em) Update $T^((t+1))$ via bounded evaporation and elitist deposit; stop early if $J^ast$ stalls for $T_"stall"$ steps.
  7. *Return* $s^ast = s^(ast (T))$ and ${mu_j (s^ast)}_(j=1)^k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Throughout, let $J_"opt" = min_(s in Omega_k) J(s)$ and let $J_"tot" = sum_(i=1)^n norm(x_i - bar(x))_2^2$ denote total variance around $bar(x) = n^(-1) sum_(i=1)^n x_i$.

#subsec("A", "Pheromone Envelope and Geometric Relaxation")

#lemma(title: "Invariant Pheromone Envelope and Geometric Relaxation", origin: "Original Derivation")[
  Under unclipped dynamics ($tau_"min" = 0, tau_"max" = infinity$) with $rho in (0, 1)$, $L >= 1$, $Q > 0$, $epsilon > 0$, let $M := (L Q) / (rho (J_"opt" + epsilon))$ and $m := (L Q) / (rho (J_"tot" + epsilon))$. For all $i <= n$, $j <= k$, and $t >= 0$:
  1. $0 <= tau_(i j)^((t)) <= (1 - rho)^t tau_0 + (1 - (1 - rho)^t) M <= max(tau_0, M)$.
  2. The row sum $S_i^((t)) := sum_(j=1)^k tau_(i j)^((t))$ satisfies
     $ (1 - rho)^t k tau_0 + (1 - (1 - rho)^t) m &<= S_i^((t)) \
       &<= (1 - rho)^t k tau_0 + (1 - (1 - rho)^t) M. $
  3. If $tilde(s)^((ell)) = s^ast$ for all $ell <= L$ and $t >= t_0$, then $tau_(i j)^((t))$ converges at exact rate $(1 - rho)$ to $tau_(i j)^infinity = (L Q) / (rho (J(s^ast) + epsilon)) bb(I)[s_i^ast = j]$.
]

#proof[
  For any $s in Omega_k$, expanding $x_i - bar(x) = (x_i - mu_j (s)) + (mu_j (s) - bar(x))$ on each $C_j (s)$ and using $sum_(i in C_j (s)) (x_i - mu_j (s)) = 0$ gives the variance identity $J_"tot" = J(s) + sum_(j=1)^k abs(C_j (s)) norm(mu_j (s) - bar(x))_2^2 >= J(s)$, so $J_"opt" <= J(tilde(s)^((ell))) <= J_"tot"$.

  *Claim 1*: Since $bb(I)[tilde(s)_i^((ell)) = j] in {0, 1}$, we have $0 <= sum_(ell=1)^L Delta tau_(i j)^((ell, t)) <= (L Q) / (J_"opt" + epsilon) = rho M$. By induction on $t$ (base $t = 0$: $tau_(i j)^((0)) = tau_0$),
  $ tau_(i j)^((t+1)) &<= (1 - rho) [(1 - rho)^t tau_0 + (1 - (1 - rho)^t) M] + rho M \
    &= (1 - rho)^(t+1) tau_0 + (1 - (1 - rho)^(t+1)) M <= max(tau_0, M). $
  Non-negativity holds since $tau_0 >= 0$, $1 - rho > 0$, and $Delta tau >= 0$.

  *Claim 2*: Since each elite assigns datum $i$ to one unique cluster, $sum_(j=1)^k bb(I)[tilde(s)_i^((ell)) = j] = 1$, giving $S_i^((t+1)) = (1 - rho) S_i^((t)) + sum_(ell=1)^L Q / (J(tilde(s)^((ell))) + epsilon)$. Because the deposit sum lies in $[rho m, rho M]$ and $S_i^((0)) = k tau_0$, applying the same linear induction to both endpoints yields the bilateral band.

  *Claim 3*: For $t >= t_0$, $tau_(i j)^((t+1)) = (1 - rho) tau_(i j)^((t)) + c_(i j)$ with $c_(i j) = (L Q) / (J(s^ast) + epsilon) bb(I)[s_i^ast = j]$. Subtracting $tau_(i j)^infinity = c_(i j) / rho$ yields $tau_(i j)^((t)) - tau_(i j)^infinity = (1 - rho)^(t - t_0) (tau_(i j)^((t_0)) - tau_(i j)^infinity)$, decaying geometrically at rate $(1 - rho)$.
]

#subsec("B", "Global Convergence and Hitting-Time Guarantee")

#theorem(title: "Almost-Sure Global SSE Convergence with Hitting-Time Bound", origin: "Original Unified Proof")[
  Suppose trails satisfy $0 < tau_"min" <= tau_(i j)^((t)) <= tau_"max" < infinity$, $q_0 in [0, 1)$, $0 < eta_"min" <= eta_(i j)^((t)) <= eta_"max" < infinity$ ($eta_"min" = eta_"max" = 1$ if $beta = 0$), and incumbent retention is monotone. Define
  $ p_"min" := (1 - q_0) (tau_"min"^alpha eta_"min"^beta) / (tau_"min"^alpha eta_"min"^beta + (k - 1) tau_"max"^alpha eta_"max"^beta) > 0. $
  Then: (1) for every $T >= 1$, $"Pr"(J(s^(ast (T))) > J_"opt") <= (1 - k! p_"min"^n)^(R T) <= exp(-R T k! p_"min"^n)$; (2) $"Pr"(lim_(T -> infinity) J(s^(ast (T))) = J_"opt") = 1$; and (3) the hitting time $tau_"hit" := inf {t >= 1 : J(s^(ast (t))) = J_"opt"}$ satisfies
  $ bb(E)[tau_"hit"] <= (1 - (1 - k! p_"min"^n)^R)^(-1) <= (k! p_"min"^n)^(-1). $
]

#proof[
  Fix $s^ast in cal(S)^ast := arg min_(s in Omega_k) J(s)$. For any permutation $pi in S_k$, $pi compose s^ast$ induces the same partition with $J(pi compose s^ast) = J_"opt"$. Since $s^ast$ is surjective, $pi_1 compose s^ast = pi_2 compose s^ast$ implies $pi_1 = pi_2$, so $abs(cal(S)^ast) >= k!$.

  Let $cal(F)_t = sigma(T^((0)), dots, T^((t)), s^(ast (0)), dots, s^(ast (t)))$ and $w_(i ell) = (tau_(i ell)^((t)))^alpha (eta_(i ell)^((t)))^beta$. Because $z |-> z / (z + C)$ increases strictly for $z, C > 0$, the exploration branch gives for every $(i, j)$:
  $ "Pr"(s_i^((r)) = j mid(|) cal(F)_t) >= (1 - q_0) w_(i j) / (w_(i j) + sum_(ell != j) w_(i ell)) >= p_"min". $
  Note $p_"min" < 1 / k$, so $k! p_"min"^n < k! / k^n <= 1$ for $n >= k$. Every $s in cal(S)^ast$ is surjective, so empty-cluster repair leaves $s$ unchanged. By coordinate independence given $cal(F)_t$ and disjointness over $cal(S)^ast$,
  $ "Pr"(s^((r)) in cal(S)^ast mid(|) cal(F)_t) = sum_(s in cal(S)^ast) product_(i=1)^n "Pr"(s_i^((r)) = s_i mid(|) cal(F)_t) >= k! p_"min"^n. $
  Any ant hitting $cal(S)^ast$ attains the global minimum $J_"opt"$, automatically entering the top $L$; monotonicity $J(cal(L)(s)) <= J(s)$ forces $J(s^(ast (t+1))) = J_"opt"$. By independence across the $R$ ants, $"Pr"(J(s^(ast (t+1))) > J_"opt" mid(|) cal(F)_t) <= (1 - k! p_"min"^n)^R$. Since $J(s^(ast (t)))$ is non-increasing, iterating the tower property over $t = T - 1, dots, 0$ and applying $1 - u <= e^(-u)$ proves Claim 1.

  For Claim 2, $sum_(T=1)^infinity exp(-R T k! p_"min"^n) < infinity$, so the first Borel-Cantelli lemma implies $J(s^(ast (T))) > J_"opt"$ occurs only finitely often almost surely; monotonicity yields $J(s^(ast (T))) -> J_"opt"$ almost surely. For Claim 3, setting $q_R := 1 - (1 - k! p_"min"^n)^R >= k! p_"min"^n$ gives $bb(E)[tau_"hit"] = sum_(T=0)^infinity "Pr"(tau_"hit" > T) <= sum_(T=0)^infinity (1 - q_R)^T = q_R^(-1) <= (k! p_"min"^n)^(-1)$.
]

While Theorem 2 guarantees global optimality asymptotically, the $p_"min"^(-n)$ scaling reflects independent coordinate sampling; practical speed therefore depends on enlarging the local-search attraction basin $cal(B)(cal(S)^ast) = {s in Omega_k : cal(L)(s) in cal(S)^ast}$ via $cal(L)$.

#subsec("C", "Label-Symmetry Interference and Premature Lock-in")

#theorem(title: "Label-Permutation Interference and Unclipped Premature Lock-in", origin: "Original Counterexample")[
  Consider pure-pheromone ACO-Clustering ($alpha = 1, beta = 0, rho in (0, 1), q_0 in [0, 1)$).
  1. *(Symmetry Washout)* If the $L = m dot.c k!$ elites at step $t$ consist of $m >= 1$ copies of all permutations ${pi compose s^ast : pi in S_k}$, every entry $(i, j)$ receives identical deposit $(m (k - 1)! Q) / (J(s^ast) + epsilon)$. Thus uniform $T^((t))$ stays strictly uniform ($p_(i j)^((t+1)) = 1 / k$) even when every elite achieves $J_"opt"$.
  2. *(Unclipped Stagnation)* Let $tau_"min" = 0, tau_"max" = infinity, tau_0 > 0, Q > 0, epsilon >= 0$, with $cal(L)$ being Lloyd @lloyd1982least or single-point relocation @hartigan1979algorithm. On the $4$-point rectangle $x_1 = (0, 0), x_2 = (0, delta), x_3 = (1, 0), x_4 = (1, delta)$ with $delta in (1 / sqrt(2), 1)$ and $k = 2$, $s_"sub" = (1, 2, 1, 2)$ has $J(s_"sub") = 1 > delta^2 = J_"opt"$ yet traps the colony permanently with strictly positive probability:
     $ "Pr"(forall t >= 0 : J(s^(ast (t))) = 1 > J_"opt") > 0. $
]

#proof[
  *Part 1*: For fixed $(i, j)$, exactly $(k - 1)!$ permutations $pi in S_k$ satisfy $pi(s_i^ast) = j$, each achieving identical SSE $J(s^ast)$. Summing over $m$ full copies of $S_k$ yields $(m (k - 1)! Q) / (J(s^ast) + epsilon)$, independent of $(i, j)$, leaving any uniform $T^((t))$ strictly uniform.

  *Part 2*: The vertical split $s^ast = (1, 1, 2, 2)$ has centroids $(0, delta / 2), (1, delta / 2)$ and $J(s^ast) = delta^2 < 1$, whereas $s_"sub" = (1, 2, 1, 2)$ has centroids $(1/2, 0), (1/2, delta)$ and $J(s_"sub") = 1$. Under Lloyd's step, each vertex lies at squared distance $1/4$ from its centroid versus $1/4 + delta^2$ from the other, so $s_"sub"$ is strictly stationary. Under single-point relocation, moving $x_1$ to $C_2$ gives $C_1' = {x_3}$ (SSE $0$) and $C_2' = {x_1, x_2, x_4}$ with centroid $(1/3, 2 delta / 3)$ and SSE
  $ J' = 2 / 3 (1 + delta^2) > 2 / 3 (1 + 1 / 2) = 1 = J(s_"sub"), $
  as $delta > 1 / sqrt(2)$; by vertex symmetry, $cal(L)(s_"sub") = s_"sub"$.

  Let $E_t := {forall r <= R, s^((r, t)) = s_"sub"}$. At $t = 0$, every ant exploring ($u > q_0$) and drawing $s_"sub"$ yields $"Pr"(E_0) >= ((1 - q_0) / 2)^(4 R) > 0$. Conditioned on $E_0 inter dots inter E_(t-1)$ ($t >= 1$), every prior elite equals $s_"sub"$. Setting $c := (L Q) / (rho (1 + epsilon)) > 0$, Lemma 1 yields $tau_(i, s_("sub", i))^((t)) = (1 - rho)^t tau_0 + (1 - (1 - rho)^t) c > (1 - rho)^t tau_0 = tau_(i, 3 - s_("sub", i))^((t))$. Thus the $q_0$-branch selects $s_("sub", i)$ deterministically, and the per-coordinate deviation probability is
  $ q_t = (1 - q_0) ((1 - rho)^t tau_0) / (2 (1 - rho)^t tau_0 + (1 - (1 - rho)^t) c). $
  Since $1 - (1 - rho)^t >= rho$ for $t >= 1$, we have $0 < q_t < 1 / 2$ and $q_t <= ((1 - q_0) tau_0) / (rho c) (1 - rho)^t$. Using $1 - u >= e^(-2 u)$ on $[0, 1/2]$ (as $dif / (dif u)[ln(1 - u) + 2 u] = (1 - 2 u) / (1 - u) >= 0$), independence across ants and coordinates yields
  $ "Pr"(inter.big_(t=1)^infinity E_t mid(|) E_0) &= product_(t=1)^infinity (1 - q_t)^(4 R) >= exp(-8 R sum_(t=1)^infinity q_t) \
    &>= exp(-(8 R (1 - q_0) tau_0 (1 - rho)) / (rho^2 c)) > 0. $
  Multiplying by $"Pr"(E_0) > 0$ proves permanent suboptimal lock-in with strictly positive probability, establishing the necessity of $tau_"min" > 0$.
]

#subsec("D", "Local-Search Descent and Geometric Invariances")

#proposition(title: "Exact Relocation Descent and Finite Monotone Refinement", origin: "Classical — Cited")[
  For $cal(L)$ using Hartigan-Wong relocation @hartigan1979algorithm or Lloyd's step @lloyd1982least inside ACO-Clustering @shelokar2004aco, moving $x_i$ from $C_a$ ($n_a = abs(C_a) >= 2$) to $C_b$ ($n_b = abs(C_b) >= 1, b != a$) changes $J(s)$ by
  $ Delta J(i : a -> b) = n_b / (n_b + 1) norm(x_i - mu_b)_2^2 - n_a / (n_a - 1) norm(x_i - mu_a)_2^2. $
  Accepting moves with $Delta J < 0$ (or Lloyd steps) terminates finitely with $J(cal(L)(s)) <= J(s)$, though global SSE minimization is NP-hard @aloise2009np.
]

#proof[
  Removing $x_i$ shifts $mu_a$ to $mu_a' = mu_a - (n_a - 1)^(-1) (x_i - mu_a)$. By the parallel-axis theorem on $C_a without {i}$, $S_a - norm(x_i - mu_a)_2^2 = S_a' + (n_a - 1) norm(mu_a' - mu_a)_2^2 = S_a' + (n_a - 1)^(-1) norm(x_i - mu_a)_2^2$, so $S_a' - S_a = -(n_a / (n_a - 1)) norm(x_i - mu_a)_2^2$. Adding $x_i$ to $C_b$ shifts $mu_b$ to $mu_b' = mu_b + (n_b + 1)^(-1) (x_i - mu_b)$, yielding $S_b' - S_b = n_b norm(mu_b - mu_b')_2^2 + norm(x_i - mu_b')_2^2 = (n_b / (n_b + 1)) norm(x_i - mu_b)_2^2$. Summing gives $Delta J$. Strict SSE decrease on every move and finiteness of $Omega_k$ ($abs(Omega_k) <= k^n$) guarantee finite termination with $J(cal(L)(s)) <= J(s)$.
]

#proposition(title: "Isometric Invariance, Scale Sensitivity, and Breakdown Point", origin: "Original Proof")[
  For Euclidean ACO-Clustering with fixed hyperparameters: (1) any rigid motion $x_i' = U x_i + b$ ($U in O(d), b in bb(R)^d$) leaves the law of ${(T^((t)), s^(ast (t)))}_(t >= 0)$ strictly invariant; (2) uniform scaling $x_i |-> c x_i$ ($c > 0$) scales deposits as $O(c^(-2))$ as $c -> infinity$ and saturates them at $Q / epsilon$ as $c -> 0^+$, whereas setting $Q = hat(Q) J_"tot"$, $epsilon = hat(epsilon) J_"tot"$, and $eta_(i j) = (1 + n norm(x_i - hat(mu)_j)_2^2 / J_"tot")^(-1)$ restores exact similarity invariance; (3) the centroid breakdown point is $epsilon_n^ast = 1 / n$.
]

#proof[
  (1) Centroid linearity yields $mu_j (s; cal(X)') = U mu_j (s; cal(X)) + b$, so orthogonality preserves all distances $norm(x_i' - mu_j')_2 = norm(x_i - mu_j)_2$, keeping $J$, $eta$, $cal(L)$, and $Delta tau$ identical. (2) Under $x_i |-> c x_i$, both $J(s)$ and $J_"tot"$ scale by $c^2$; substituting $Q = hat(Q) J_"tot"$ and $epsilon = hat(epsilon) J_"tot"$ gives $Delta tau = hat(Q) / (J(tilde(s)^((ell))) / J_"tot" + hat(epsilon))$ and cancels $c^2$ identically in $eta_(i j)$. (3) Displacing $x_1$ to $norm(x_1)_2 = M -> infinity$ with $R_cal(X) = max_(i >= 2) norm(x_i)_2$ gives, for $j_1 = s_1$ and $m = abs(C_(j_1)) <= n$,
  $ norm(mu_(j_1) (s))_2 >= M / m - (m - 1) / m R_cal(X) >= M / n - R_cal(X) -> infinity, $
  so corrupting a single datum drives the centroids unbounded: $epsilon_n^ast = 1 / n$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Per iteration $t < T_"max"$, ACO-Clustering executes four stages: (i) *ant construction* costs $O(R n (k + d))$ for direct row-CDF sampling plus centroid / SSE evaluation (or $O(n k + R n d)$ via alias sampling when $beta = 0$, plus $O(n k d)$ if $beta > 0$); (ii) *elite selection* takes $O(R log R)$ (or $O(R)$ via quickselect); (iii) *local search* applies up to $I_"ls"$ sweeps to $L$ elites in $O(L I_"ls" n k d)$; (iv) *pheromone update* takes $O(n (k + L))$. Total time is $O(T_"max" n (R (k + d) + L I_"ls" k d))$. By comparison, spatial ant-sorting runs in $O(T_"walk" R_"ant" s^2 d)$ on grid window $s^2$, avoiding $k$-dependence but needing $T_"walk" >> n^2$ steps to coalesce heaps.

Space complexity is $O(n (d + k + L))$ when streaming ants through an elite heap (or $O(n (d + k + R))$ batched). The $n times k$ matrix $T^((t))$ ($8 n k$ bytes in FP64) stays wholly L3-cache-resident for $n <= 10^4, k <= 32$ ($<= 2.56$ MB), and all $R$ ant constructions are *embarrassingly parallel* over read-only $T^((t))$. Geometrically, Euclidean SSE yields convex Voronoi cells biased toward isotropic, equal-variance clusters.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T_"max" n (R (k + d) + L I_"ls" k d))$],
    [Space Complexity], [$O(n (d + k + L))$ streaming; $O(n (d + k + R))$ batched],
    [Decision Geometry], [Convex Euclidean Voronoi cells (spherical bias)],
    [Parallelism], [Embarrassingly parallel across $R$ ants and $L$ elites],
    [Hyperparameters], [$k, R, L, rho, q_0, tau_"min", tau_"max", I_"ls"$ (high tuning load)],
    [Invariances], [Rigid-motion invariant; scale-invariant if $Q, epsilon prop J_"tot"$],
    [Outlier Robustness], [Breakdown point $epsilon_n^ast = 1 / n$ (non-robust quadratic loss)],
    [Convergence], [Almost-sure global SSE optimum iff $tau_"min" > 0, q_0 < 1$],
  ),
  caption: [Operational, computational, and geometric profile of ACO-Clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Chemometrics and Analytical Spectroscopy*: Introduced in analytical chemistry @shelokar2004aco, ACO-Clustering excels on IR, NMR, and chromatographic spectra ($n ~ 10^2$–$10^3$) where multi-modal SSE landscapes trap single-start $k$-means; stigmergic exploration plus local search consistently recovers compact spectral classes.
2. *Wireless Sensor Network (WSN) Energy-Aware Zoning*: Partitioning $n$ sensor nodes into $k$ cluster-head zones often adds non-smooth battery-balance and hop penalties to Euclidean dispersion. Because ACO reinforces via black-box fitness $Q / (J(s) + epsilon)$, it optimizes composite discrete costs without closed-form gradients.
3. *Bioinformatics Gene Co-Expression Modules*: On expression profiles of $n ~ 10^2$–$10^3$ genes across conditions, pheromone consensus smooths noisy local SSE plateaus where deterministic initializations fluctuate markedly across replicates.
4. *Operations Research Zoning and Customer Segmentation*: In logistics catchment design and CRM segmentation, Euclidean SSE models service-radius transport cost directly while cache-resident $O(n k)$ trails deliver sub-$10$ ms multi-basin optimization.
5. *Exploratory Document and Visual Collection Sorting*: When $k$ is unknown and only cosine or Jaccard similarity is available, spatial ant-sorting organizes collections into emergent 2D density heaps analogous to Self-Organizing Maps @kohonen1982som.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure modes follow directly from Section II: (1) *Symmetry Washout (Theorem 3.1)* cancels trail signal unless elite labels are canonicalized (e.g., lexicographic centroid sorting or Hungarian matching to $s^(ast (t))$); (2) *Unclipped Collapse (Theorem 3.2)* traps the colony permanently when $tau_"min" = 0$; (3) *Non-Convex Manifolds* are sliced erroneously by convex Voronoi boundaries (use DBSCAN @ester1996dbscan or Spectral Clustering @ng2002spectral); (4) *High-Dimensional Concentration* @beyer1999nearest flattens $J(tilde(s)^((ell))) / J_"tot" -> 1$, eroding pheromone contrast; and (5) *Massive Scale / Outliers* hit the $O(n k)$ memory wall and $1 / n$ breakdown point (prefer mini-batch $k$-means @sculley2010webscale, coresets @harpeled2004coresets, trimmed $k$-means @cuesta1997trimmed, or PAM @kaufman1990finding).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Chemometrics & spectra], [Multi-modal SSE; moderate $n$; spherical modes], [Canonicalize labels; set $tau_"min" > 0$],
    [WSN / constrained zoning], [Black-box fitness $Q / (J + epsilon)$ handles penalties], [Scale-normalize $Q, epsilon$ with $J_"tot"$],
    [Gene co-expression], [Consensus trail memory stabilizes noisy basins], [Apply PCA before high-$d$ Euclidean],
    [Non-convex manifolds], [Poor: SSE forces convex Voronoi cuts], [DBSCAN @ester1996dbscan or Spectral @ng2002spectral],
    [Web-scale / heavy outliers], [Poor: $O(n k)$ memory; breakdown $1 / n$], [Mini-batch @sculley2010webscale; PAM @kaufman1990finding],
  ),
  caption: [Domain suitability matrix and failure-regime alternatives for ACO-Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

To position ACO-Clustering alongside Genetic $K$-Means (GKA) @krishna1999gka, PSO-Clustering @vandermerwe2003pso, and cluster ensembles @strehl2002cluster, we prove that its sampling kernel is an online discounted weighted consensus ensemble.

#proposition(title: "Discounted Consensus Ensemble Representation and Reductions", origin: "Original Proof")[
  Under unclipped pure-pheromone dynamics ($alpha = 1, beta = 0, tau_"min" = 0, tau_"max" = infinity$) with $w_(ell, u) := (J(tilde(s)^((ell, u))) + epsilon)^(-1)$, define $W_(i j)^((u)) := sum_(ell=1)^L w_(ell, u) bb(I)[tilde(s)_i^((ell, u)) = j]$ and $W^((u)) := sum_(ell=1)^L w_(ell, u)$. Then for all $t >= 1$,
  $ p_(i j)^((t)) = ((1 - rho)^t tau_0 + Q sum_(u=0)^(t-1) (1 - rho)^(t-1-u) W_(i j)^((u))) / (k (1 - rho)^t tau_0 + Q sum_(u=0)^(t-1) (1 - rho)^(t-1-u) W^((u))). $
  Hence: (i) as $rho -> 1^-$, $p_(i j)^((t)) = W_(i j)^((t-1)) / W^((t-1))$ equals the inverse-SSE weighted consensus frequency of the prior $L$ elites; (ii) for $rho -> 1^-$ and $L = 1$, $p_(i j)^((t)) = bb(I)[tilde(s)_i^((1, t-1)) = j]$, reducing to deterministic single-elite local search @lloyd1982least.
]

#proof[
  Unrolling $tau_(i j)^((t)) = (1 - rho) tau_(i j)^((t-1)) + Q W_(i j)^((t-1))$ from $u = 0$ to $t - 1$ gives $tau_(i j)^((t)) = (1 - rho)^t tau_0 + Q sum_(u=0)^(t-1) (1 - rho)^(t-1-u) W_(i j)^((u))$. Because $sum_(j=1)^k W_(i j)^((u)) = W^((u))$, summing over $j <= k$ yields the exact row sum $S_i^((t)) = k (1 - rho)^t tau_0 + Q sum_(u=0)^(t-1) (1 - rho)^(t-1-u) W^((u))$, and dividing gives $p_(i j)^((t))$. Letting $rho -> 1^-$ annihilates $(1 - rho)^t$ and all terms $u < t - 1$, leaving $W_(i j)^((t-1)) / W^((t-1))$, which collapses to $bb(I)[tilde(s)_i^((1, t-1)) = j]$ when $L = 1$.
]

Thus, whereas GKA @krishna1999gka evolves discrete strings via biased mutation and Fuzzy $c$-Means @bezdek1981fcm optimizes smooth memberships analytically, ACO-Clustering compresses elite history into a factorized consensus table @strehl2002cluster. Conversely, PSO-Clustering @vandermerwe2003pso searches continuous centroid space $bb(R)^(k times d)$, avoiding ACO's $O(n k)$ trail table and $k!$-label interference (Theorem 3.1) while forfeiting ACO's direct extensibility to non-centroid discrete losses.

*Empirical Behavior (ClusterOverview 4-Dataset Benchmark)*: Evaluated via Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing, ACO-Clustering attains *Iris* $"Sil" = 0.447, "ARI" = 0.583$; *Customers* $"Sil" = 0.365$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.303$ ($"ARI" = "N/A"$); with *mean runtime* $6.771$ ms. These results mirror our analysis: on *Multi-Shape* ($"Sil" = 0.712$), hybrid colony search plus local refinement isolates well-separated low-SSE modes cleanly; on *Iris* ($"Sil" = 0.447, "ARI" = 0.583$), overlapping _versicolor_ / _virginica_ boundaries create shallow SSE gradients where finite-budget label-symmetry interference (Theorem 3.1) trails tuned $k$-means++ @arthur2007kmeanspp; on *Customers* ($"Sil" = 0.365$), heavy-tailed spending exposes Proposition 5's $1 / n$ breakdown sensitivity; on *High-D* ($"Sil" = 0.303$), distance concentration @beyer1999nearest flattens inverse-SSE weights $w_(ell, u)$ (Proposition 6); and the $6.771$ ms runtime confirms Section III's cache-resident operational profile.
