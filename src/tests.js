// Exact tests: the p-value is a rational number computed by counting, with
// no approximation by a continuous distribution.
import { Rational } from "@ratmath/core";
import { combinations, binomial as choose, rng as makeRng } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, TWO, Q, whole, probability, compare, show, lowOf, highOf, interval, isQ, sqrt as sqrtAny } from "./num.js";
import { binomial, hypergeometric } from "./discrete-dists.js";
import { beta, normal } from "./continuous.js";
import { values } from "./describe.js";

const ALTERNATIVES = new Set(["two-sided", "less", "greater"]);
const alternativeOf = (a) => {
  if (!ALTERNATIVES.has(a)) refuse("bad-argument", `the alternative must be "two-sided", "less" or "greater", got "${a}"`);
  return a;
};

/**
 * The outcome of a test.
 *   test, statistic, pValue, alternative, kind ("exact", "certified" or "simulated")
 *   decision(alpha)   "reject", "keep", or "too close to call" when the
 *                     p-value is an enclosure that contains alpha
 *   reject(alpha)     true or false; refuses when it is too close to call
 */
export class TestResult {
  constructor(fields) {
    Object.assign(this, fields);
    Object.freeze(this);
  }
  decision(alpha = "0.05") {
    const where = compare(this.pValue, probability(alpha, "the level α"));
    return where === "above" ? "keep" : where === "too close to call" ? where : "reject";
  }
  reject(alpha = "0.05") {
    const d = this.decision(alpha);
    if (d === "too close to call") refuse("too-close", `the p-value ${show(this.pValue, 12)} cannot be told apart from ${alpha}; ask for more precision`);
    return d === "reject";
  }
  toString() {
    return `${this.test}: p = ${show(this.pValue, 6)} (${this.kind}, ${this.alternative})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/**
 * An interval for a binomial proportion from k successes in n trials.
 *   method "exact"    Clopper and Pearson's interval: its coverage is never
 *                     below the level. The ends are quantiles of beta
 *                     distributions, certified and rounded outward.
 *   method "wilson"   Wilson's score interval (certified)
 *   method "wald"     the textbook p̂ ± z·√(p̂(1−p̂)/n) (certified); its
 *                     coverage can be far below the level
 * Returns { low, high, level, method, estimate }.
 */
export function proportionInterval(successes, trials, { level = "0.95", method = "exact", bits = 64 } = {}) {
  const n = whole(trials, "the number of trials", { min: 1 });
  const k = whole(successes, "the number of successes");
  if (k > n) refuse("bad-argument", "more successes than trials");
  const conf = probability(level, "the level");
  if (conf.numerator === 0n || conf.equals(ONE)) refuse("bad-argument", "the level must be strictly between 0 and 1");
  const tail = ONE.subtract(conf).divide(TWO);
  const estimate = Q(k, n);
  if (method === "exact") {
    const key = `${k}|${n}|${conf}|${bits}`;
    let ends = KEPT.get(key);
    if (!ends) {
      const low = k === 0 ? ZERO : lowOf(beta(k, n - k + 1).quantile(tail, { bits }));
      const high = k === n ? ONE : highOf(beta(k + 1, n - k).quantile(ONE.subtract(tail), { bits }));
      ends = [low, high];
      if (KEPT.size >= 4096) KEPT.clear();
      KEPT.set(key, ends);
    }
    return { low: ends[0], high: ends[1], level: conf, method, estimate };
  }
  const z = zFor(ONE.subtract(tail), bits);
  if (method === "wald") {
    const half = mulIv(z, sqrtIv(estimate.multiply(ONE.subtract(estimate)).divide(Q(n)), bits));
    return { low: estimate.subtract(highOf(half)), high: estimate.add(highOf(half)), level: conf, method, estimate };
  }
  if (method === "wilson") {
    // (p̂ + z²/2n ± z·√(p̂(1−p̂)/n + z²/4n²)) / (1 + z²/n), evaluated at both ends of the enclosure of z
    const at = (zq, sign) => {
      const z2 = zq.multiply(zq);
      const centre = estimate.add(z2.divide(Q(2 * n)));
      const root = sqrtIv(estimate.multiply(ONE.subtract(estimate)).divide(Q(n)).add(z2.divide(Q(4 * n * n))), bits);
      const denom = ONE.add(z2.divide(Q(n)));
      return sign < 0 ? centre.subtract(zq.multiply(highOf(root))).divide(denom) : centre.add(zq.multiply(highOf(root))).divide(denom);
    };
    const lo = [at(lowOf(z), -1), at(highOf(z), -1)].reduce((a, b) => (a.compareTo(b) < 0 ? a : b));
    const hi = [at(lowOf(z), 1), at(highOf(z), 1)].reduce((a, b) => (a.compareTo(b) > 0 ? a : b));
    return { low: lo.numerator < 0n ? ZERO : lo, high: hi.compareTo(ONE) > 0 ? ONE : hi, level: conf, method, estimate };
  }
  return refuse("bad-argument", `unknown method "${method}"; use "exact", "wilson" or "wald"`);
}
const STANDARD = normal(0, 1);
const zFor = (p, bits) => STANDARD.quantile(p, { bits });
const KEPT = new Map();
const sqrtIv = (q, bits) => sqrtAny(q, bits);
const mulIv = (a, b) => {
  // both non-negative
  const lo = lowOf(a).multiply(lowOf(b));
  const hi = highOf(a).multiply(highOf(b));
  return lo.equals(hi) ? lo : interval(lo, hi);
};

/**
 * The exact binomial test of H₀: p = p0, from k successes in n trials.
 * The two-sided p-value adds up every outcome that is no more likely than
 * the one observed.
 */
export function binomialTest(successes, trials, p0 = "1/2", { alternative = "two-sided" } = {}) {
  const n = whole(trials, "the number of trials", { min: 1 });
  const k = whole(successes, "the number of successes");
  if (k > n) refuse("bad-argument", "more successes than trials");
  const p = probability(p0, "p0");
  const alt = alternativeOf(alternative);
  const d = binomial(n, p);
  let pValue;
  if (alt === "less") pValue = d.cdf(k);
  else if (alt === "greater") pValue = d.prob(k, n);
  else {
    const observed = d.pmf(k);
    pValue = ZERO;
    for (let i = 0; i <= n; i++) {
      const q = d.pmf(i);
      if (q.compareTo(observed) <= 0) pValue = pValue.add(q);
    }
  }
  return new TestResult({
    test: "exact binomial test", kind: "exact", statistic: k, n, estimate: Q(k, n), null: p, alternative: alt, pValue,
    interval: (options) => proportionInterval(k, n, options),
  });
}

/**
 * Fisher's exact test on a 2×2 table [[a, b], [c, d]]: with the margins
 * fixed, a follows a hypergeometric distribution under independence.
 * "greater" means the odds ratio ad/bc is above 1.
 */
export function fisherTest(table, { alternative = "two-sided" } = {}) {
  if (!Array.isArray(table) || table.length !== 2 || table.some((r) => !Array.isArray(r) || r.length !== 2)) refuse("bad-argument", "fisherTest needs a 2×2 table [[a, b], [c, d]]");
  const [[a, b], [c, d]] = table.map((r) => r.map((v) => whole(v, "a count")));
  const alt = alternativeOf(alternative);
  const n = a + b + c + d;
  if (n === 0) refuse("no-data", "the table is empty");
  const dist = hypergeometric(n, a + c, a + b); // a: marked items among the first row's draws
  let pValue;
  if (alt === "less") pValue = dist.cdf(a);
  else if (alt === "greater") pValue = dist.prob(a, dist.max);
  else {
    const observed = dist.pmf(a);
    pValue = ZERO;
    for (let i = dist.min; i <= dist.max; i++) {
      const q = dist.pmf(i);
      if (q.compareTo(observed) <= 0) pValue = pValue.add(q);
    }
  }
  return new TestResult({
    test: "Fisher's exact test", kind: "exact", statistic: a, table: [[a, b], [c, d]], alternative: alt, pValue,
    oddsRatio: b * c === 0 ? (a * d === 0 ? null : Infinity) : Q(a * d, b * c),
    null: dist,
  });
}

const STATISTICS = {
  mean: (xs) => xs.reduce((s, v) => s.add(v), ZERO).divide(Q(xs.length)),
  median(xs) {
    const s = [...xs].sort((p, q) => p.compareTo(q));
    const m = s.length >> 1;
    return s.length % 2 ? s[m] : s[m - 1].add(s[m]).divide(TWO);
  },
};
/**
 * A permutation test for two groups: if the groups do not differ, every way
 * of dealing the pooled values into groups of these sizes is equally likely.
 * The statistic is statistic(x) − statistic(y), with "mean" (the default),
 * "median", or a function of a list of Rationals.
 *
 * When the number of deals is at most `limit` they are all counted and the
 * p-value is exact. Otherwise `resamples` random deals are drawn from `seed`
 * and the answer is simulated: it carries the seed and an exact interval for
 * the p-value it estimates.
 *
 * With { paired: true } the two lists are pairs, and the deals are the 2ⁿ
 * ways of swapping within pairs.
 */
export function permutationTest(x, y, { statistic = "mean", alternative = "two-sided", paired = false, limit = 200000, resamples = 9999, seed = 1 } = {}) {
  const xs = values(x, "x");
  const ys = values(y, "y");
  if (xs.some((v) => !isQ(v)) || ys.some((v) => !isQ(v))) refuse("needs-exact", "a permutation test needs exact data, not intervals");
  const alt = alternativeOf(alternative);
  const stat = typeof statistic === "function" ? (list) => toQ(statistic(list), "the statistic") : STATISTICS[statistic];
  if (!stat) refuse("bad-argument", `unknown statistic "${statistic}"; use "mean", "median" or a function`);
  if (paired && xs.length !== ys.length) refuse("bad-argument", "paired data need lists of the same length");
  const diff = (p, q) => stat(p).subtract(stat(q));
  const observed = diff(xs, ys);
  const extreme = (d) => (alt === "greater" ? d.compareTo(observed) >= 0 : alt === "less" ? d.compareTo(observed) <= 0 : d.abs().compareTo(observed.abs()) >= 0);
  const pooled = [...xs, ...ys];
  const n = xs.length;
  const total = paired ? 1n << BigInt(n) : choose(pooled.length, n);
  const base = { test: paired ? "paired permutation test" : "permutation test", statistic: observed, alternative: alt, arrangements: total };

  const deal = (chosen) => {
    // chosen: the positions (0-based, increasing) that go to the first group
    const first = [];
    const second = [];
    let c = 0;
    for (let i = 0; i < pooled.length; i++) {
      if (c < chosen.length && chosen[c] === i) {
        first.push(pooled[i]);
        c++;
      } else second.push(pooled[i]);
    }
    return diff(first, second);
  };
  const swap = (mask) => {
    const first = [];
    const second = [];
    for (let i = 0; i < n; i++) {
      const flip = (mask >> BigInt(i)) & 1n;
      first.push(flip ? ys[i] : xs[i]);
      second.push(flip ? xs[i] : ys[i]);
    }
    return diff(first, second);
  };

  if (total <= BigInt(limit)) {
    let count = 0n;
    if (paired) {
      for (let mask = 0n; mask < total; mask++) if (extreme(swap(mask))) count++;
    } else {
      const positions = Array.from({ length: pooled.length }, (_, i) => i);
      for (const chosen of combinations(positions, n)) if (extreme(deal(chosen))) count++;
    }
    return new TestResult({ ...base, kind: "exact", pValue: new Rational(count, total), count });
  }
  const B = whole(resamples, "resamples", { min: 1 });
  const g = makeRng(seed).fork("permutation test");
  let hits = 0;
  const positions = Array.from({ length: pooled.length }, (_, i) => i);
  for (let r = 0; r < B; r++) {
    const d = paired ? swap(g.bits(n)) : deal(g.sample(positions, n).sort((p, q) => p - q));
    if (extreme(d)) hits++;
  }
  const ci = proportionInterval(hits, B, { level: "0.99" });
  return new TestResult({
    ...base, kind: "simulated", seed, resamples: B, hits,
    // counting the observed arrangement itself keeps the test valid
    pValue: Q(hits + 1, B + 1),
    pValueInterval: { low: ci.low, high: ci.high, level: ci.level },
  });
}
