// Estimation, confidence intervals and the standard tests.
//
// Each test statistic is computed exactly from rational data. A z or t
// statistic is a square root, so its square is kept exact and the p-value is
// taken from the square: P(|T| ≥ |t|) = P(F(1, ν) ≥ t²) and
// P(|Z| ≥ |z|) = P(χ²₁ ≥ z²). The p-values and the ends of the intervals
// are certified enclosures, so a decision at a level is either settled or
// reported as too close to call.
import { Rational, RationalInterval } from "@ratmath/core";
import { enclose } from "@ratmath/analytic";
import { refuse, toQ, Q, ZERO, ONE, TWO, HALF, add, sub, mul, div, neg, sqrt, lowOf, highOf, interval, isQ, whole, probability, min as minQ, max as maxQ } from "./num.js";
import { values } from "./describe.js";
import { TestResult, proportionInterval } from "./tests.js";
import { normal, studentT, chiSquare, fDist } from "./continuous.js";
import { crossing } from "./solve.js";

const tidy = (v) => (v instanceof RationalInterval && v.low.equals(v.high) ? v.low : v);
const ALTERNATIVES = new Set(["two-sided", "less", "greater"]);
export const alternativeOf = (a) => {
  if (!ALTERNATIVES.has(a)) refuse("bad-argument", `the alternative must be "two-sided", "less" or "greater", got "${a}"`);
  return a;
};
export const levelOf = (level) => {
  const c = probability(level, "the level");
  if (c.numerator === 0n || c.equals(ONE)) refuse("bad-argument", "the level must be strictly between 0 and 1");
  return c;
};
/** Exact data: rationals only. */
export function exactValues(data, what = "the data") {
  const xs = values(data, what);
  if (xs.some((v) => !isQ(v))) refuse("needs-exact", `${what} must be exact numbers here, not intervals`);
  return xs;
}
const meanOf = (xs) => xs.reduce((s, v) => s.add(v), ZERO).divide(Q(xs.length));
/** Σ(x − x̄)². */
const sumSquares = (xs, m = meanOf(xs)) => xs.reduce((s, v) => {
  const d = v.subtract(m);
  return s.add(d.multiply(d));
}, ZERO);
const sign = (q) => (q.numerator > 0n ? 1 : q.numerator < 0n ? -1 : 0);
const clamp01 = (v) => {
  if (isQ(v)) return v.numerator < 0n ? ZERO : v.compareTo(ONE) > 0 ? ONE : v;
  const lo = v.low.numerator < 0n ? ZERO : v.low;
  const hi = v.high.compareTo(ONE) > 0 ? ONE : v.high;
  return tidy(interval(lo, hi));
};
/** The signed square root of an exact square: the statistic z or t from z² or t². */
export function signedRoot(square, s, bits) {
  const r = sqrt(square, bits);
  return s < 0 ? neg(r) : r;
}
/**
 * The p-value from the tail of a squared statistic. `two` is P(S ≥ s²) for
 * the square, which is the two-sided p-value; one side is half of it or one
 * minus half of it, depending on the direction the data point.
 */
export function sidedP(two, s, alternative) {
  if (alternative === "two-sided") return two;
  if (s === 0) return HALF;
  const half = div(two, TWO);
  const towards = alternative === "greater" ? s > 0 : s < 0;
  return clamp01(towards ? half : sub(ONE, half));
}
/** Twice the smaller tail, at most 1: the two-sided p-value of an asymmetric statistic. */
function twoTails(lower, upper) {
  const lo = minQ(lowOf(lower), lowOf(upper)).multiply(TWO);
  const hi = minQ(highOf(lower), highOf(upper)).multiply(TWO);
  return clamp01(tidy(interval(lo, hi)));
}
const tailP = (dist, x, alternative, bits) => {
  if (alternative === "greater") return dist.sf(x, { bits });
  if (alternative === "less") return dist.cdf(x, { bits });
  return twoTails(dist.cdf(x, { bits }), dist.sf(x, { bits }));
};
const kindOf = (p) => (isQ(p) ? "exact" : "certified");

/**
 * estimate ± q·√se², with q the quantile of `dist` that leaves (1 − level)/2
 * in each tail; one-sided when `alternative` is "less" or "greater" (the open
 * end is null). The ends are rationals, rounded outward from their enclosures.
 */
export function centredInterval(estimate, se2, dist, { level = "0.95", alternative = "two-sided", bits = 64, method } = {}) {
  const conf = levelOf(level);
  const alt = alternativeOf(alternative);
  const tail = alt === "two-sided" ? ONE.subtract(conf).divide(TWO) : ONE.subtract(conf);
  const q = dist.quantile(ONE.subtract(tail), { bits });
  const half = highOf(mul(q, sqrt(se2, bits)));
  return {
    low: alt === "less" ? null : lowOf(sub(estimate, half)),
    high: alt === "greater" ? null : highOf(add(estimate, half)),
    level: conf, estimate, method,
  };
}

// ---- z and t -------------------------------------------------------------------------------

/** The t test from its pieces: estimate, null value, se² and degrees of freedom. */
function tResult(name, estimate, nullValue, se2, df, alternative, bits, extra = {}) {
  const diff = estimate.subtract(nullValue);
  if (se2.numerator === 0n) refuse("undefined", "the standard error is 0 (every value is the same), so the t statistic is not defined");
  const t2 = diff.multiply(diff).divide(se2);
  const s = sign(diff);
  const pValue = sidedP(fDist(1, df).sf(t2, { bits }), s, alternative);
  return new TestResult({
    test: name, kind: kindOf(pValue), statistic: signedRoot(t2, s, bits), statisticSquared: t2, df, estimate, null: nullValue,
    standardError: sqrt(se2, bits), alternative, pValue,
    interval: (options = {}) => centredInterval(estimate, se2, studentT(df), { bits, method: "t", ...options }),
    ...extra,
  });
}

/**
 * The z test for a mean when the standard deviation σ is known.
 *   zTest(data, { sigma, mu = 0, alternative })
 */
export function zTest(data, { sigma, mu = 0, alternative = "two-sided", bits = 64 } = {}) {
  const xs = exactValues(data);
  if (sigma === undefined) refuse("bad-argument", "the z test needs the known standard deviation σ; without it use tTest");
  const s = toQ(sigma, "σ");
  if (s.numerator <= 0n) refuse("bad-argument", "σ must be positive");
  const m0 = toQ(mu, "μ0");
  const alt = alternativeOf(alternative);
  const n = xs.length;
  const xbar = meanOf(xs);
  const se2 = s.multiply(s).divide(Q(n));
  const d = xbar.subtract(m0);
  const z2 = d.multiply(d).divide(se2);
  const pValue = sidedP(chiSquare(1).sf(z2, { bits }), sign(d), alt);
  return new TestResult({
    test: "z test", kind: kindOf(pValue), statistic: signedRoot(z2, sign(d), bits), statisticSquared: z2, n, estimate: xbar, null: m0,
    standardError: sqrt(se2, bits), alternative: alt, pValue,
    interval: (options = {}) => centredInterval(xbar, se2, normal(0, 1), { bits, method: "z", ...options }),
  });
}

const isOptions = (v) => v !== null && typeof v === "object" && !Array.isArray(v) && typeof v[Symbol.iterator] !== "function";

/**
 * Student's t tests.
 *   tTest(x, { mu = 0 })                       one sample: is the mean μ0?
 *   tTest(x, y, { paired: true })              the mean of the differences x − y
 *   tTest(x, y)                                two samples, Welch's test (unequal variances)
 *   tTest(x, y, { equalVariance: true })       two samples, pooled variance
 * The statistic is t; statisticSquared is t² exactly. Welch's degrees of
 * freedom are the exact rational of the Welch–Satterthwaite formula.
 */
export function tTest(x, y, options) {
  if (isOptions(y) || y === undefined) [y, options] = [undefined, y];
  const { mu = 0, paired = false, equalVariance = false, alternative = "two-sided", bits = 64 } = options ?? {};
  const alt = alternativeOf(alternative);
  const m0 = toQ(mu, "μ0");
  const xs = exactValues(x, "x");
  if (y === undefined) {
    if (paired) refuse("bad-argument", "a paired test needs two lists");
    return oneSample("one-sample t test", xs, m0, alt, bits);
  }
  const ys = exactValues(y, "y");
  if (paired) {
    if (xs.length !== ys.length) refuse("bad-argument", `paired data need lists of the same length (${xs.length} and ${ys.length})`);
    return oneSample("paired t test", xs.map((v, i) => v.subtract(ys[i])), m0, alt, bits, { pairs: xs.length });
  }
  const [n1, n2] = [xs.length, ys.length];
  if (n1 < 2 || n2 < 2) refuse("no-data", "each sample needs at least two values");
  const [m1, m2] = [meanOf(xs), meanOf(ys)];
  const [v1, v2] = [sumSquares(xs, m1).divide(Q(n1 - 1)), sumSquares(ys, m2).divide(Q(n2 - 1))];
  const estimate = m1.subtract(m2);
  if (equalVariance) {
    const sp2 = v1.multiply(Q(n1 - 1)).add(v2.multiply(Q(n2 - 1))).divide(Q(n1 + n2 - 2));
    const se2 = sp2.multiply(Q(1, n1).add(Q(1, n2)));
    return tResult("two-sample t test (pooled)", estimate, m0, se2, Q(n1 + n2 - 2), alt, bits, { n: [n1, n2], means: [m1, m2], variances: [v1, v2], pooledVariance: sp2 });
  }
  const [a, b] = [v1.divide(Q(n1)), v2.divide(Q(n2))];
  const se2 = a.add(b);
  if (se2.numerator === 0n) refuse("undefined", "both samples are constant, so the standard error is 0");
  const df = se2.multiply(se2).divide(a.multiply(a).divide(Q(n1 - 1)).add(b.multiply(b).divide(Q(n2 - 1))));
  return tResult("Welch two-sample t test", estimate, m0, se2, df, alt, bits, { n: [n1, n2], means: [m1, m2], variances: [v1, v2] });
}
function oneSample(name, xs, m0, alt, bits, extra = {}) {
  const n = xs.length;
  if (n < 2) refuse("no-data", "a t test needs at least two values");
  const m = meanOf(xs);
  const s2 = sumSquares(xs, m).divide(Q(n - 1));
  return tResult(name, m, m0, s2.divide(Q(n)), Q(n - 1), alt, bits, { n, variance: s2, ...extra });
}

/**
 * A confidence interval for a mean: t when σ is unknown, z when it is given.
 * { low, high, level, estimate, method }; one-sided with { alternative }.
 */
export function meanInterval(data, { level = "0.95", sigma, alternative = "two-sided", bits = 64 } = {}) {
  const xs = exactValues(data);
  const n = xs.length;
  const m = meanOf(xs);
  if (sigma !== undefined) {
    const s = toQ(sigma, "σ");
    if (s.numerator <= 0n) refuse("bad-argument", "σ must be positive");
    return centredInterval(m, s.multiply(s).divide(Q(n)), normal(0, 1), { level, alternative, bits, method: "z" });
  }
  if (n < 2) refuse("no-data", "a t interval needs at least two values");
  return centredInterval(m, sumSquares(xs, m).divide(Q(n - 1)).divide(Q(n)), studentT(n - 1), { level, alternative, bits, method: "t" });
}

// ---- variances -----------------------------------------------------------------------------

/**
 * The chi-square interval for a variance from normal data:
 * ((n − 1)s²/χ²(1 − α/2), (n − 1)s²/χ²(α/2)). Also `sd`, the interval for σ.
 */
export function varianceInterval(data, { level = "0.95", bits = 64 } = {}) {
  const xs = exactValues(data);
  const n = xs.length;
  if (n < 2) refuse("no-data", "the variance needs at least two values");
  const conf = levelOf(level);
  const tail = ONE.subtract(conf).divide(TWO);
  const ss = sumSquares(xs);
  const d = chiSquare(n - 1);
  const upper = d.quantile(ONE.subtract(tail), { bits });
  const lower = d.quantile(tail, { bits });
  const low = ss.divide(highOf(upper));
  const high = ss.divide(lowOf(lower));
  return { low, high, level: conf, estimate: ss.divide(Q(n - 1)), method: "chi-square", sd: { low: lowOf(sqrt(low, bits)), high: highOf(sqrt(high, bits)) } };
}
/** The chi-square test of H₀: σ² = σ0² from normal data. The statistic (n − 1)s²/σ0² is exact. */
export function varianceTest(data, variance0, { alternative = "two-sided", bits = 64 } = {}) {
  const xs = exactValues(data);
  const n = xs.length;
  if (n < 2) refuse("no-data", "the variance needs at least two values");
  const v0 = toQ(variance0, "σ0²");
  if (v0.numerator <= 0n) refuse("bad-argument", "σ0² must be positive");
  const alt = alternativeOf(alternative);
  const ss = sumSquares(xs);
  const X = ss.divide(v0);
  const pValue = tailP(chiSquare(n - 1), X, alt, bits);
  return new TestResult({
    test: "chi-square test for a variance", kind: kindOf(pValue), statistic: X, df: n - 1, n, estimate: ss.divide(Q(n - 1)), null: v0, alternative: alt, pValue,
    interval: (options = {}) => varianceInterval(xs, { bits, ...options }),
  });
}
/**
 * The F test of H₀: σ₁²/σ₂² = ratio for two normal samples. The statistic
 * (s₁²/s₂²)/ratio is exact; the interval is for σ₁²/σ₂².
 */
export function fTest(x, y, { ratio = 1, alternative = "two-sided", bits = 64 } = {}) {
  const xs = exactValues(x, "x");
  const ys = exactValues(y, "y");
  const [n1, n2] = [xs.length, ys.length];
  if (n1 < 2 || n2 < 2) refuse("no-data", "each sample needs at least two values");
  const r0 = toQ(ratio, "the ratio");
  if (r0.numerator <= 0n) refuse("bad-argument", "the ratio must be positive");
  const alt = alternativeOf(alternative);
  const v1 = sumSquares(xs).divide(Q(n1 - 1));
  const v2 = sumSquares(ys).divide(Q(n2 - 1));
  if (v2.numerator === 0n) refuse("undefined", "the second sample is constant, so the ratio of variances is not defined");
  const est = v1.divide(v2);
  const F = est.divide(r0);
  const dist = fDist(n1 - 1, n2 - 1);
  const pValue = tailP(dist, F, alt, bits);
  return new TestResult({
    test: "F test for two variances", kind: kindOf(pValue), statistic: F, df: [n1 - 1, n2 - 1], estimate: est, null: r0, alternative: alt, pValue,
    interval({ level = "0.95" } = {}) {
      const conf = levelOf(level);
      const tail = ONE.subtract(conf).divide(TWO);
      return { low: est.divide(highOf(dist.quantile(ONE.subtract(tail), { bits }))), high: est.divide(lowOf(dist.quantile(tail, { bits }))), level: conf, estimate: est, method: "F" };
    },
  });
}

// ---- proportions ---------------------------------------------------------------------------

/**
 * The z test for proportions.
 *   proportionTest(k, n, p0)                 one sample: z² = (k − n·p0)²/(n·p0(1 − p0))
 *   proportionTest([k1, k2], [n1, n2])       two samples, with the pooled proportion
 * { correct: true } applies Yates's continuity correction (one sample).
 * For small samples binomialTest and fisherTest are exact.
 */
export function proportionTest(successes, trials, p0 = "1/2", { alternative = "two-sided", correct = false, bits = 64 } = {}) {
  if (isOptions(p0)) [p0, { alternative = "two-sided", correct = false, bits = 64 } = {}] = ["1/2", p0];
  const alt = alternativeOf(alternative);
  if (Array.isArray(successes)) {
    if (!Array.isArray(trials) || successes.length !== 2 || trials.length !== 2) refuse("bad-argument", "two samples need [k1, k2] and [n1, n2]");
    const [n1, n2] = trials.map((t) => whole(t, "the number of trials", { min: 1 }));
    const [k1, k2] = successes.map((k) => whole(k, "the number of successes"));
    if (k1 > n1 || k2 > n2) refuse("bad-argument", "more successes than trials");
    const [p1, p2] = [Q(k1, n1), Q(k2, n2)];
    const pooled = Q(k1 + k2, n1 + n2);
    if (pooled.numerator === 0n || pooled.equals(ONE)) refuse("undefined", "every trial had the same outcome, so the pooled variance is 0");
    const d = p1.subtract(p2);
    const se2 = pooled.multiply(ONE.subtract(pooled)).multiply(Q(1, n1).add(Q(1, n2)));
    const z2 = d.multiply(d).divide(se2);
    const pValue = sidedP(chiSquare(1).sf(z2, { bits }), sign(d), alt);
    const wald = p1.multiply(ONE.subtract(p1)).divide(Q(n1)).add(p2.multiply(ONE.subtract(p2)).divide(Q(n2)));
    return new TestResult({
      test: "two-sample z test for proportions", kind: kindOf(pValue), statistic: signedRoot(z2, sign(d), bits), statisticSquared: z2,
      estimate: d, proportions: [p1, p2], pooled, null: ZERO, alternative: alt, pValue,
      interval: (options = {}) => centredInterval(d, wald, normal(0, 1), { bits, method: "wald", ...options }),
    });
  }
  const n = whole(trials, "the number of trials", { min: 1 });
  const k = whole(successes, "the number of successes");
  if (k > n) refuse("bad-argument", "more successes than trials");
  const p = probability(p0, "p0");
  if (p.numerator === 0n || p.equals(ONE)) refuse("bad-argument", "p0 must be strictly between 0 and 1");
  let d = Q(k).subtract(Q(n).multiply(p));
  if (correct) {
    const a = d.abs().subtract(HALF);
    d = a.numerator <= 0n ? ZERO : sign(d) < 0 ? a.negate() : a;
  }
  const z2 = d.multiply(d).divide(Q(n).multiply(p).multiply(ONE.subtract(p)));
  const pValue = sidedP(chiSquare(1).sf(z2, { bits }), sign(d), alt);
  return new TestResult({
    test: `z test for a proportion${correct ? " (continuity corrected)" : ""}`, kind: kindOf(pValue), statistic: signedRoot(z2, sign(d), bits), statisticSquared: z2,
    n, estimate: Q(k, n), null: p, alternative: alt, pValue,
    interval: (options = {}) => proportionInterval(k, n, { method: "wilson", bits, ...options }),
  });
}

// ---- chi-square ----------------------------------------------------------------------------

const isTable = (o) => Array.isArray(o) && o.length > 0 && o.every((r) => Array.isArray(r));
/**
 * Pearson's chi-square tests.
 *   chiSquareTest(counts, { p })             goodness of fit to probabilities p (default: all equal)
 *   chiSquareTest(counts, { expected })      … or to expected counts with the same total
 *   chiSquareTest([[…], […]])                independence in a table of counts
 * Options: estimated (parameters estimated from the data, subtracted from the
 * degrees of freedom), correct (Yates, 2×2 tables only). X² is exact; the
 * result lists the expected counts and each cell's contribution, and
 * `smallExpected` is true when an expected count is below 5, where the
 * chi-square approximation is poor.
 */
export function chiSquareTest(observed, { p, expected, estimated = 0, correct = false, bits = 64 } = {}) {
  const est = whole(estimated, "the number of estimated parameters");
  if (isTable(observed)) {
    const O = observed.map((r) => r.map((v) => Q(whole(v, "a count"))));
    const c = O[0].length;
    if (O.some((r) => r.length !== c)) refuse("bad-argument", "every row of the table needs the same number of cells");
    if (O.length < 2 || c < 2) refuse("bad-argument", "a test of independence needs at least two rows and two columns");
    const rows = O.map((r) => r.reduce((s, v) => s.add(v), ZERO));
    const cols = O[0].map((_, j) => O.reduce((s, r) => s.add(r[j]), ZERO));
    const N = rows.reduce((s, v) => s.add(v), ZERO);
    if (rows.some((v) => v.numerator === 0n) || cols.some((v) => v.numerator === 0n)) refuse("undefined", "a row or column of the table is empty, so its expected counts are 0");
    const E = rows.map((r) => cols.map((col) => r.multiply(col).divide(N)));
    const yates = correct && O.length === 2 && c === 2;
    if (correct && !yates) refuse("bad-argument", "Yates's correction is for 2×2 tables");
    const contributions = O.map((r, i) => r.map((o, j) => cell(o, E[i][j], yates)));
    const X2 = contributions.flat().reduce((s, v) => s.add(v), ZERO);
    const df = (O.length - 1) * (c - 1) - est;
    return chiResult(`chi-square test of independence${yates ? " (Yates)" : ""}`, X2, df, E, contributions, bits);
  }
  const O = values(observed, "the counts").map((v) => Q(whole(v, "a count")));
  const k = O.length;
  if (k < 2) refuse("bad-argument", "a goodness-of-fit test needs at least two categories");
  const N = O.reduce((s, v) => s.add(v), ZERO);
  if (N.numerator === 0n) refuse("no-data", "every count is 0");
  let E;
  if (expected !== undefined) {
    E = values(expected, "the expected counts").map((v) => toQ(v, "an expected count"));
    const total = E.reduce((s, v) => s.add(v), ZERO);
    if (!total.equals(N)) refuse("bad-argument", `the expected counts add up to ${total}, not to the ${N} observed`);
  } else {
    const ps = p === undefined ? O.map(() => Q(1, k)) : values(p, "the probabilities").map((v) => probability(v, "a probability"));
    const total = ps.reduce((s, v) => s.add(v), ZERO);
    if (!total.equals(ONE)) refuse("bad-argument", `the probabilities add up to ${total}, not 1`);
    E = ps.map((q) => q.multiply(N));
  }
  if (E.length !== k) refuse("bad-argument", "one expected value for each count");
  if (E.some((e) => e.numerator <= 0n)) refuse("bad-argument", "every expected count must be positive");
  const contributions = O.map((o, i) => cell(o, E[i], false));
  const X2 = contributions.reduce((s, v) => s.add(v), ZERO);
  return chiResult("chi-square goodness-of-fit test", X2, k - 1 - est, E, contributions, bits);
}
function cell(o, e, yates) {
  let d = o.subtract(e).abs();
  if (yates) d = d.compareTo(HALF) > 0 ? d.subtract(HALF) : ZERO;
  return d.multiply(d).divide(e);
}
function chiResult(name, X2, df, expected, contributions, bits) {
  if (df < 1) refuse("bad-argument", "no degrees of freedom are left");
  const pValue = chiSquare(df).sf(X2, { bits });
  const small = expected.flat().some((e) => e.compareTo(Q(5)) < 0);
  return new TestResult({ test: name, kind: kindOf(pValue), statistic: X2, df, expected, contributions, smallExpected: small, alternative: "greater", pValue });
}

// ---- correlation ---------------------------------------------------------------------------

/**
 * The test of zero correlation for paired normal data: t² = r²(n − 2)/(1 − r²)
 * on n − 2 degrees of freedom, exact since r² is rational. The interval is
 * Fisher's: tanh(atanh r ± z/√(n − 3)).
 */
export function correlationTest(x, y, { alternative = "two-sided", bits = 64 } = {}) {
  const xs = exactValues(x, "x");
  const ys = exactValues(y, "y");
  if (xs.length !== ys.length) refuse("bad-argument", `x has ${xs.length} values and y has ${ys.length}`);
  const n = xs.length;
  if (n < 3) refuse("no-data", "a correlation test needs at least three pairs");
  const alt = alternativeOf(alternative);
  const [mx, my] = [meanOf(xs), meanOf(ys)];
  const sxy = xs.reduce((s, v, i) => s.add(v.subtract(mx).multiply(ys[i].subtract(my))), ZERO);
  const [sxx, syy] = [sumSquares(xs, mx), sumSquares(ys, my)];
  if (sxx.numerator === 0n || syy.numerator === 0n) refuse("undefined", "the correlation is not defined when one variable is constant");
  const r2 = sxy.multiply(sxy).divide(sxx.multiply(syy));
  const s = sign(sxy);
  const r = signedRoot(r2, s, bits);
  let pValue;
  let t2 = null;
  if (r2.equals(ONE)) pValue = alt === "two-sided" || (alt === "greater") === s > 0 ? ZERO : ONE;
  else {
    t2 = r2.multiply(Q(n - 2)).divide(ONE.subtract(r2));
    pValue = sidedP(fDist(1, n - 2).sf(t2, { bits }), s, alt);
  }
  return new TestResult({
    test: "test of zero correlation", kind: kindOf(pValue), statistic: t2 === null ? null : signedRoot(t2, s, bits), statisticSquared: t2, df: n - 2, n,
    estimate: r, rSquared: r2, null: ZERO, alternative: alt, pValue,
    interval: (options = {}) => fisherInterval(r, n, { bits, ...options }),
  });
}
function fisherInterval(r, n, { level = "0.95", bits = 64 } = {}) {
  if (n < 4) refuse("no-data", "Fisher's interval needs at least four pairs");
  const conf = levelOf(level);
  const [rl, rh] = [lowOf(r), highOf(r)];
  if (rh.compareTo(ONE) >= 0 || rl.compareTo(ONE.negate()) <= 0) return { low: rl, high: rh, level: conf, estimate: r, method: "Fisher z" };
  const atanh = (v) => div(enclose.log(ONE.add(v).divide(ONE.subtract(v)), bits + 8), TWO);
  const q = normal(0, 1).quantile(ONE.subtract(ONE.subtract(conf).divide(TWO)), { bits });
  const half = highOf(div(q, sqrt(Q(n - 3), bits + 8)));
  const zl = lowOf(atanh(rl)).subtract(half);
  const zh = highOf(atanh(rh)).add(half);
  return { low: lowOf(enclose.tanh(zl, bits)), high: highOf(enclose.tanh(zh, bits)), level: conf, estimate: r, method: "Fisher z" };
}

// ---- estimation ----------------------------------------------------------------------------

const FAMILIES = ["normal", "exponential", "poisson", "bernoulli", "binomial", "geometric", "uniform", "gamma", "beta", "lognormal"];
const positiveData = (xs, family) => {
  if (xs.some((v) => v.numerator <= 0n)) refuse("bad-argument", `the ${family} family needs positive data`);
};
const wholeData = (xs, family) => {
  if (xs.some((v) => v.denominator !== 1n || v.numerator < 0n)) refuse("bad-argument", `the ${family} family needs whole, non-negative data`);
};
/**
 * Fit a family of distributions to data, by maximum likelihood ("mle", the
 * default) or the method of moments ("moments").
 *   fit(data, "normal")  →  { family, method, n, params, distribution }
 * Parameters are exact when they can be (all the closed forms); the gamma
 * maximum-likelihood shape and the lognormal parameters involve logarithms
 * and are certified enclosures, and then `distribution` is null.
 * Families: normal, exponential, poisson, bernoulli, binomial ({ trials }),
 * geometric (counting trials, or failures with { failures: true }), uniform,
 * gamma, beta (moments only), lognormal (maximum likelihood only).
 */
export function fit(data, family, { method = "mle", trials, failures = false, bits = 64 } = {}) {
  if (!FAMILIES.includes(family)) refuse("bad-argument", `unknown family "${family}"; use one of ${FAMILIES.join(", ")}`);
  if (method !== "mle" && method !== "moments") refuse("bad-argument", `the method is "mle" or "moments", got "${method}"`);
  const xs = exactValues(data);
  const n = xs.length;
  const m = meanOf(xs);
  const v = sumSquares(xs, m).divide(Q(n)); // the divisor n: both methods agree on it
  const out = (params, distribution = null) => ({ family, method, n, params, distribution });
  switch (family) {
    case "normal": {
      const sigma = sqrt(v, bits);
      return out({ mu: m, variance: v, sigma }, isQ(sigma) && sigma.numerator > 0n ? normal(m, sigma) : null);
    }
    case "exponential":
      positiveData(xs, family);
      return out({ rate: m.reciprocal() });
    case "poisson":
      wholeData(xs, family);
      return out({ lambda: m });
    case "bernoulli":
      if (xs.some((x) => !x.equals(ZERO) && !x.equals(ONE))) refuse("bad-argument", "Bernoulli data are 0s and 1s");
      return out({ p: m });
    case "binomial": {
      const t = whole(trials, "the number of trials", { min: 1 });
      wholeData(xs, family);
      if (xs.some((x) => x.compareTo(Q(t)) > 0)) refuse("bad-argument", `a count is above the ${t} trials`);
      return out({ trials: t, p: m.divide(Q(t)) });
    }
    case "geometric":
      wholeData(xs, family);
      if (!failures && xs.some((x) => x.numerator === 0n)) refuse("bad-argument", "counting trials, every value is at least 1 (use { failures: true } to count failures)");
      return out({ p: failures ? ONE.divide(ONE.add(m)) : m.reciprocal(), failures });
    case "uniform": {
      if (method === "mle") {
        const sorted = [...xs].sort((a, b) => a.compareTo(b));
        return out({ a: sorted[0], b: sorted[n - 1] });
      }
      const half = sqrt(v.multiply(Q(3)), bits);
      return out({ a: tidy(sub(m, half)), b: tidy(add(m, half)) });
    }
    case "gamma": {
      positiveData(xs, family);
      if (v.numerator === 0n) refuse("undefined", "the data are constant");
      if (method === "moments") return out({ shape: m.multiply(m).divide(v), rate: m.divide(v) });
      return out(gammaMle(xs, m, bits));
    }
    case "beta": {
      if (method === "mle") refuse("not-available", "the beta maximum-likelihood equations have no closed form and are not solved here; use { method: \"moments\" }");
      if (xs.some((x) => x.numerator <= 0n || x.compareTo(ONE) >= 0)) refuse("bad-argument", "beta data lie strictly between 0 and 1");
      if (v.numerator === 0n) refuse("undefined", "the data are constant");
      const common = m.multiply(ONE.subtract(m)).divide(v).subtract(ONE);
      if (common.numerator <= 0n) refuse("undefined", "the sample variance is too large for any beta distribution");
      return out({ a: m.multiply(common), b: ONE.subtract(m).multiply(common) });
    }
    case "lognormal": {
      if (method === "moments") refuse("not-available", "fit the lognormal by maximum likelihood (the moments of the logs)");
      positiveData(xs, family);
      const logs = xs.map((x) => tidy(enclose.log(x, bits + 8)));
      const mu = div(logs.reduce(add, ZERO), Q(n));
      const s2 = div(logs.reduce((s, l) => add(s, intervalSquare(sub(l, mu))), ZERO), Q(n));
      return out({ mu: tidy(mu), variance: tidy(s2), sigma: tidy(sqrt(s2, bits)) });
    }
  }
  /* c8 ignore next */
  return refuse("bad-argument", family);
}
const intervalSquare = (a) => {
  if (isQ(a)) return a.multiply(a);
  const [l, h] = [a.low.multiply(a.low), a.high.multiply(a.high)];
  if (a.low.numerator <= 0n && a.high.numerator >= 0n) return interval(ZERO, maxQ(l, h));
  return interval(minQ(l, h), maxQ(l, h));
};
/**
 * The gamma maximum-likelihood shape k solves log k − ψ(k) = s with
 * s = log x̄ − mean(log x) > 0; the left side decreases from ∞ to 0, so the
 * root is enclosed by bisection on certified values at both ends of s.
 */
function gammaMle(xs, m, bits) {
  const n = xs.length;
  const logs = xs.map((x) => tidy(enclose.log(x, bits + 16)));
  const s = tidy(sub(tidy(enclose.log(m, bits + 16)), div(logs.reduce(add, ZERO), Q(n))));
  if (highOf(s).numerator <= 0n) refuse("undefined", "the data are too close to constant for a gamma fit");
  // G(k) = ψ(k) − log k increases from −∞ to 0; solve G(k) = −s
  const G = (k, b) => sub(tidy(enclose.digamma(k, b)), tidy(enclose.log(k, b)));
  const sf = (lowOf(s).toNumber() + highOf(s).toNumber()) / 2;
  const guess = (3 - sf + Math.sqrt((sf - 3) ** 2 + 24 * sf)) / (12 * sf); // Minka's starting value
  if (lowOf(s).numerator <= 0n) refuse("too-close", "log x̄ − mean(log x) cannot be told apart from 0 at this precision; ask for more bits");
  // the larger s, the smaller k
  const kLow = lowOf(crossing(G, highOf(s).negate(), { guess, lo: ZERO, bits }));
  const kHigh = highOf(crossing(G, lowOf(s).negate(), { guess, lo: ZERO, bits }));
  const shape = tidy(interval(minQ(kLow, kHigh), maxQ(kLow, kHigh)));
  return { shape, rate: tidy(div(shape, m)) };
}
