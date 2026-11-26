// What this package is still planned to contain, phase by phase, from the
// RatMath library roadmap. Phases 1 to 3 are built; a phase is removed from
// this list when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
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
