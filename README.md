# @ratmath/stats

Probability and statistics for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: set up, nothing built yet.** This repository holds the package
skeleton and the docs site. The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#6-a-probability-and-statistics-package-ratmathstats).

Probability and statistics in one library. Every answer is labelled exact, certified or simulated, and data can be generated from a model you choose, disturbed with noise from another distribution, and handed to the tests and regressions to see whether they recover what was put in.

## Planned

| Phase | Contents |
|---|---|
| 1 | Finite probability, discrete distributions, descriptive statistics, exact tests |
| 2 | Continuous distributions, certified |
| 3 | Seeded generation, models with noise and contamination, experiments |
| 4 | Estimation, confidence intervals, the standard tests |
| 5 | Regression and analysis of variance |
| 6 | Rank tests, permutation tests, the bootstrap |
| 7 | Markov chains; conjugate Bayesian updating |
| 8 | Graphics hooks |

It will depend on: Core, analytic (the distribution functions), linalg (regression) and discrete (counting and the seeded generator).

## Layout

- `index.js`, `src/`: the library. For now it exports only `plan`, the list
  above, so that the tests and the docs build have something to load.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, then one
  page per topic as topics are built). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/stats.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
