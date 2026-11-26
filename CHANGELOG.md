# Changelog

## Unreleased

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
