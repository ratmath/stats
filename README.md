# @ratmath/stats

Probability and statistics for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: phases 1 to 5 of 8 are built.** Docs: <https://ratmath.com/stats/>.

Probability and statistics in one library, with three kinds of answer that are
never mixed up:

- **exact**: fractions. Finite probability, discrete distributions,
  descriptive statistics, exact tests.
- **certified**: enclosures, as narrow as you ask. Continuous distributions
  and their quantiles. Comparing a certified p-value with a level can answer
  "too close to call".
- **simulated**: seeded random experiments that carry their seed and their
  own sampling error.

Standard inference sits on top: z, t, chi-square and F tests, proportions,
fitted distributions, linear regression and analysis of variance. Statistics
such as t², X² and F stay exact fractions, p-values are certified, and seeded
experiments in the test suite check that intervals cover and tests keep their
stated size.

Data can be generated from a model you choose, disturbed with noise,
contamination, error in x, rounding, censoring and missing values, and handed
to a procedure to see whether it recovers what was put in.

```js
import { dice, binomialTest, normal, show, model, estimatorStudy, fast } from "@ratmath/stats";

dice(2).prob(([a, b]) => a + b === 7);          // 1/6
binomialTest(14, 20, "1/2").pValue;             // 15115/131072
show(normal().quantile("0.975"), 18);           // "1.95996398454005424…"

const line = model({ truth: [2, 0.5], x: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], noise: normal(0, 1) });
String(estimatorStudy({ model: line, n: 50, truth: 0.5, runs: 2000, seed: 1,
  estimator: (data) => fast.leastSquares(data.x, data.y).slope }));
// "mean 0.49993, sd 0.049765, bias -0.000066071 ± 0.0021811, rmse 0.049753 (simulated: 2000 runs, seed 1)"
```

## Phases

| Phase | Contents | |
|---|---|---|
| 1 | Finite probability, discrete distributions, descriptive statistics, exact tests | built |
| 2 | Continuous distributions, certified | built |
| 3 | Seeded generation, models with noise and contamination, experiments | built |
| 4 | Estimation, confidence intervals, the standard tests | built |
| 5 | Regression and analysis of variance | built |
| 6 | Rank tests, permutation tests, the bootstrap | planned |
| 7 | Markov chains; conjugate Bayesian updating | planned |
| 8 | Graphics hooks | planned |

The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#6-a-probability-and-statistics-package-ratmathstats).

It depends on Core, analytic (the certified special functions) and discrete
(counting and the seeded generator, exported again here as `rng`). Until those
are published, analytic and discrete are vendored as git submodules under
`vendor/` and used as workspaces: clone with `--recurse-submodules`.

## Layout

- `index.js`, `src/`: the library.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, and one
  page per topic). `bun run docs:build` writes `docs/`, typesets the
  mathematics and runs the code box at the end of every page.
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
