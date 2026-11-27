// @ratmath/stats: probability and statistics for RatMath.
//
// Three kinds of answer, never mixed up:
//   exact      rationals: finite probability, descriptive statistics, exact tests
//   certified  enclosures: values of continuous distributions, to any precision
//   simulated  results of seeded random experiments, with their sampling error

export { StatsError, toQ, toValue, show, compare, sqrt, toNumber } from "./src/num.js";
export { Finite, finite, die, dice, coin, deck, draw, bayes } from "./src/finite.js";
export { Discrete, bernoulli, binomial, geometric, negativeBinomial, hypergeometric, discreteUniform, poisson } from "./src/discrete-dists.js";
export {
  count, sum, mean, weightedMean, variance, sd, moment, skewness, kurtosis, min, max, range,
  quantile, median, quartiles, iqr, fiveNumber, mad, trimmedMean, frequencies, modes, histogram,
  covariance, correlation, leastSquares, summary, fast,
} from "./src/describe.js";
export { Continuous, normal, lognormal, gamma, chiSquare, exponential, weibull, beta, studentT, fDist, uniform, cauchy } from "./src/continuous.js";
export { TestResult, binomialTest, fisherTest, permutationTest, proportionInterval } from "./src/tests.js";
export { Model, Dataset, model } from "./src/model.js";
export {
  zTest, tTest, meanInterval, varianceInterval, varianceTest, fTest, proportionTest, chiSquareTest, correlationTest, fit,
} from "./src/inference.js";
export { power, sampleSize } from "./src/power.js";
export { fromDensity, transform, orderStatistic } from "./src/density.js";
export { LinearModel, linearModel, polyFit, dummies, compareModels, anova } from "./src/regression.js";
export { SimulatedRate, SimulatedEstimates, howOften, coverage, rejectionRate, powerCurve, estimatorStudy } from "./src/experiment.js";
export {
  ranks, signTest, wilcoxonTest, mannWhitneyTest, kruskalWallisTest, spearmanTest, kendallTest,
  ksTest, kolmogorovCdf, kolmogorovTail, smirnovTail, bootstrap, Bootstrap,
} from "./src/rank.js";
export { MarkovChain, markovChain } from "./src/markov.js";
export { betaPrior, gammaPrior, normalPrior, dirichletPrior } from "./src/conjugate.js";
export { rng, Random } from "@ratmath/discrete";
export { plan } from "./src/plan.js";
