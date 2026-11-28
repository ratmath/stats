# Changelog

## Unreleased

- Docs: the browser bundle is rebuilt against Core 0.12 (it still rounded 65
  to 60 at the tens and could not mix `floor()` results with Rationals), is
  built without identifier minification (Bun renames classes even with
  `keepNames`), and the docs build now runs every code box against the
  emitted bundle after a mixed Integer/Rational check.
- Distinct outcomes are no longer merged when their text runs together: the
  internal outcome key now length-prefixes tuple and set members, so
  `finite([["a,sb", "c"], ["a", "b,sc"]])` has two outcomes of probability
  1/2 each instead of one of probability 1. Values inside plain-object
  outcomes are keyed by type too (`{ x: "1/2" }` and `{ x: 1/2 }` differ;
  `{ x: 5 }` and `{ x: 5n }` agree). Found by the external QA review.
- Requires Core 0.12: the `@ratmath/core` peer range is now `>=0.12.0` (Core 0.12 returns Integers from `floor`/`ceil`/`trunc`/`round` and rounds half up by default).
- `fromDensity` no longer presents an enclosure of the mass that contains 1 as
  a proof: `checks.massProven` (true only for an exact mass) and
  `checks.massEnclosure` are added, and `{ requireExactMass: true }` refuses
  such densities. (1 + 10⁻⁵⁰)·cos(x)/sin(1) on [0, 1] was accepted as "proven".
- `sampleSize` returns `n` as an exact bigint above 2^53 − 1 instead of a
  rounded Number, and with no effect (mu = mu0, p = p0) returns n = 1 for a
  target power up to alpha instead of throwing.
- Products are implicit: `product`, `power` and `dice` keep their independent
  parts instead of listing every tuple (seven dice used to take about 190 MB
  and nine were refused). One outcome's probability multiplies the parts,
  `marginal` returns the parts, `prob(event)` and `map(f)` stream the outcomes
  one at a time with integer weights and divide once, `prob(event, given)`
  counts both in one pass, and `sample` reads each part off one integer.
  Answers, outcome order and seeded samples are unchanged; listing (`entries`,
  `outcomes`, `table`, `given`) is still limited to two million outcomes.
  New: `count` (bigint), `factors`, `each()`, and `sum()` (the total of tuple
  coordinates, by convolution for a product). Docs: a "Many dice" page.
  Found by Astra's dice storage report.
- Review of batches A3 and A10 (Course Test B1): `fromDensity` on an unbounded
  support now says its nonnegativity check stops at the window it proved
  ("not checked beyond that"), instead of claiming a leading-term argument it
  does not make.
- `fromDensity(f, support, { normalize })`: a continuous distribution from a
  density formula or a list of pieces, proven nonnegative with mass 1 (exactly
  for polynomial pieces, which keep the CDF and moments exact).
- `transform(dist, g)` / `dist.transform(g)`: the distribution of g(X) for g
  proven monotone on the support.
- `orderStatistic(dist, k, n)` / `dist.orderStatistic(k, n)`: the k-th smallest
  of n draws, exact for uniform and polynomial densities.
- Docs: a "Your own densities" page. Course Test batch A10 (ratmath/shelf#10).
- `meanInterval` and `varianceInterval` take a summary in place of data:
  `{ n, mean, sd }`, `{ n, mean, variance }`, `{ n, variance }`, `{ n, sd }`.
- `power({ test: "z" | "proportion", ... })`: closed-form power as an
  enclosure (the proportion test by the normal approximation).
- `sampleSize(...)`: the smallest whole n proven to meet a margin of error
  (proportion or mean) or a target power.
- Docs: a "Power and sample size" page.
- Quartiles the school way: `quartiles`, `iqr`, `fiveNumber` and `quantile`
  take `{ method: "halves" }` (medians of the lower and upper halves, the
  median left out when n is odd, as in most textbooks and the TI-84) and
  `{ method: "tukey" }` (the median kept in both halves). Found by the Course
  Test (ratmath/shelf#7).
- Finite distributions accept Core's `Integer` (what ``R`1` `` gives) wherever
  they need a number: `mean(f)` and friends no longer refuse it, and an
  `Integer` outcome is the same outcome as the equal `Rational`.

- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

## 0.3.0

### Added

- Rank tests: `ranks`, `signTest`, `wilcoxonTest` (with the Hodges–Lehmann estimate), `mannWhitneyTest`, `kruskalWallisTest`, `spearmanTest`, `kendallTest`. Null distributions are counted exactly with ties (midranks), with a certified normal approximation for large samples; `method` says which.
- Kolmogorov–Smirnov: `ksTest` for one sample (against a certified Continuous or an exact cdf) and two samples (exact by counting lattice paths); `kolmogorovCdf` (Marsaglia–Tsang–Wang, exact), `smirnovTail`, `kolmogorovTail`.
- `bootstrap`: percentile and basic intervals, exact over every resample for small samples, otherwise seeded.
- Markov chains: `markovChain` with n-step distributions, classes and periods, stationary distributions, absorption probabilities and steps (the fundamental matrix), hitting times, mean return times and seeded paths, all exact.
- Conjugate updating: `betaPrior`, `gammaPrior` (Poisson or exponential), `normalPrior` (known variance), `dirichletPrior`, with exact posterior parameters and predictive probabilities and certified credible intervals.
- Tests compare with published values (R's examples, Marsaglia–Tsang–Wang) and with brute-force counts over every sign pattern, split and path.
- Docs: pages on rank tests, Kolmogorov–Smirnov, the bootstrap, Markov chains and Bayesian updating.

## 0.2.0

### Added

- Inference: `zTest`, `tTest` (one-sample, paired, Welch, pooled), `meanInterval`, `varianceInterval`, `varianceTest`, `fTest`, `proportionTest` (one and two samples, Wilson intervals), `chiSquareTest` (goodness of fit and independence, Yates' correction), `correlationTest`. Statistics are exact where possible (t², X², F); p-values are certified; interval ends are rounded outward.
- `fit`: maximum likelihood or moments for ten families.
- Regression: `linearModel` with exact coefficients, coefficient tests and intervals, prediction intervals, leverage and Cook's distance; `polyFit`, `dummies`, `compareModels`, `anova`.
- Tests that check identities between procedures exactly, and seeded experiments that check coverage and size.
- Docs: pages on tests for means, variances and proportions, chi-square tests, fitting distributions, testing the tests, linear regression and analysis of variance; a two-column panel on the explore page.

## 0.1.0

### Added

- Finite probability: `finite`, `die`, `dice`, `coin`, `deck`, `draw`, `bayes` and the `Finite` class, exact in fractions.
- Discrete distributions: `bernoulli`, `binomial`, `geometric`, `negativeBinomial`, `hypergeometric`, `discreteUniform` (exact) and `poisson` (certified).
- Descriptive statistics: means, variances, moments, quantiles by nine named rules, histograms, correlation, `leastSquares`, `summary`; enclosures for interval data; floating-point versions under `fast`.
- Exact tests: `binomialTest`, `fisherTest`, `permutationTest` (exact, or simulated and labelled when the deals cannot be counted), `proportionInterval`.
- Continuous distributions with certified values and quantiles: `normal`, `lognormal`, `gamma`, `chiSquare`, `exponential`, `weibull`, `beta`, `studentT`, `fDist`, `uniform`, `cauchy`.
- `show` and `compare` for enclosures; comparisons can answer "too close to call".
- Seeded sampling from every distribution, and exact rational draws.
- `model`: data from a known truth with noise, unequal spread, contamination, error in x, rounding, censoring and missing values.
- Experiments: `howOften`, `coverage`, `rejectionRate`, `powerCurve`, `estimatorStudy`, each carrying its seed and sampling error.
- Docs site with a page per topic.
