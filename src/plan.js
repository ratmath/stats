// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Nothing here is built yet; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 1,
    "contents": "Finite probability, discrete distributions, descriptive statistics, exact tests",
    "estimate": "1,100"
  },
  {
    "phase": 2,
    "contents": "Continuous distributions, certified",
    "estimate": "1,200, plus 500 in analytic"
  },
  {
    "phase": 3,
    "contents": "Seeded generation, models with noise and contamination, experiments",
    "estimate": "1,300"
  },
  {
    "phase": 4,
    "contents": "Estimation, confidence intervals, the standard tests",
    "estimate": "1,000"
  },
  {
    "phase": 5,
    "contents": "Regression and analysis of variance",
    "estimate": "1,100"
  },
  {
    "phase": 6,
    "contents": "Rank tests, permutation tests, the bootstrap",
    "estimate": "700"
  },
  {
    "phase": 7,
    "contents": "Markov chains; conjugate Bayesian updating",
    "estimate": "700"
  },
  {
    "phase": 8,
    "contents": "Graphics hooks",
    "estimate": "800 in graphics"
  }
].map(Object.freeze));
