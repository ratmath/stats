// What this package is still planned to contain, phase by phase, from the
// RatMath library roadmap. Phases 1 to 5 are built; a phase is removed from
// this list when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
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
