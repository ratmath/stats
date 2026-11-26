// Conjugate Bayesian updating. For these pairs of prior and likelihood the
// posterior is in the prior's family, with parameters that are simple sums of
// the data, so updating is exact arithmetic. Posterior means, variances and
// modes are exact; credible intervals are certified quantiles; predictive
// probabilities for counts are exact.
//
//   beta prior          × binomial (or Bernoulli, geometric) data   → beta
//   gamma prior (rate)  × Poisson counts or exponential times       → gamma
//   normal prior        × normal data with known variance           → normal
//   Dirichlet prior     × multinomial counts                        → Dirichlet
import { Rational } from "@ratmath/core";
import { binomial as choose } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, TWO, Q, whole, probability, show, isQ } from "./num.js";
import { values } from "./describe.js";
import { beta, gamma, normal } from "./continuous.js";

const positive = (v, what) => {
  const q = toQ(v, what);
  if (q.numerator <= 0n) refuse("bad-argument", `${what} must be positive`);
  return q;
};
const exactList = (data, what) => {
  const xs = values(data, what);
  if (xs.some((v) => !isQ(v))) refuse("needs-exact", `${what} must be exact values`);
  return xs;
};
/** Γ(a + k)/Γ(a) for a whole k: the rising factorial (a)ₖ. */
const rising = (a, k) => {
  let s = ONE;
  for (let i = 0; i < k; i++) s = s.multiply(a.add(Q(i)));
  return s;
};

class Belief {
  constructor(fields) {
    // keep getters as getters: the distribution is built only when asked for
    Object.defineProperties(this, Object.getOwnPropertyDescriptors(fields));
    Object.freeze(this);
  }
  /** An equal-tailed credible interval, certified: the level/2 and 1 − level/2 quantiles of the posterior. */
  interval(level = "0.95", { bits = 64 } = {}) {
    const L = probability(level, "the level");
    const a = ONE.subtract(L).divide(TWO);
    const d = this.distribution;
    return Object.freeze({ low: d.quantile(a, { bits }), high: d.quantile(ONE.subtract(a), { bits }), level: L });
  }
  /** P(θ ≤ x) under this belief, certified. */
  cdf(x, options) {
    return this.distribution.cdf(x, options);
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/**
 * A beta prior for a probability p: Beta(a, b), as if a − 1 successes and
 * b − 1 failures had been seen. Beta(1, 1) is uniform.
 *   .update({ successes, failures })   or { successes, trials }, or a list of 0/1 outcomes
 *   .mean, .variance, .mode (exact), .distribution, .interval(level)
 *   .predictive(m)   the beta-binomial law of the successes in m more trials, exact
 */
export function betaPrior(a = 1, b = 1, history = []) {
  const A = positive(a, "a");
  const B = positive(b, "b");
  const s = A.add(B);
  const mode = A.compareTo(ONE) > 0 && B.compareTo(ONE) > 0 ? A.subtract(ONE).divide(s.subtract(TWO)) : null;
  return new Belief({
    family: "beta", params: { a: A, b: B }, history,
    mean: A.divide(s),
    variance: A.multiply(B).divide(s.multiply(s).multiply(s.add(ONE))),
    mode,
    get distribution() {
      return beta(A, B);
    },
    update(data) {
      let k;
      let f;
      if (Array.isArray(data)) {
        const xs = exactList(data, "the outcomes");
        if (xs.some((x) => !x.equals(ZERO) && !x.equals(ONE))) refuse("bad-argument", "outcomes must be 0 or 1");
        k = xs.filter((x) => x.equals(ONE)).length;
        f = xs.length - k;
      } else {
        k = whole(data?.successes ?? 0, "successes");
        f = data?.failures !== undefined ? whole(data.failures, "failures") : data?.trials !== undefined ? whole(data.trials, "trials") - k : refuse("bad-argument", "give failures or trials");
        if (f < 0) refuse("bad-argument", "more successes than trials");
      }
      return betaPrior(A.add(Q(k)), B.add(Q(f)), [...history, { successes: k, failures: f }]);
    },
    predictive(m) {
      const M = whole(m, "m");
      const pmf = (k) => new Rational(choose(M, k)).multiply(rising(A, k)).multiply(rising(B, M - k)).divide(rising(s, M));
      return Object.freeze({ name: `beta-binomial(${M}, ${A}, ${B})`, min: 0, max: M, pmf: (k) => (k < 0 || k > M ? ZERO : pmf(k)), table: Array.from({ length: M + 1 }, (_, k) => pmf(k)) });
    },
    toString() {
      return `Beta(${A}, ${B}): mean ${show(A.divide(s), 6)}`;
    },
  });
}

/**
 * A gamma prior for a rate λ: Gamma(shape α, rate β).
 *   likelihood "poisson"       .update(counts): α + Σx, β + n
 *   likelihood "exponential"   .update(times):  α + n,  β + Σx
 *   .predictive()  for Poisson: the negative binomial law of the next count, as a pmf, exact
 */
export function gammaPrior(shape, rate, { likelihood = "poisson" } = {}, history = []) {
  const a = positive(shape, "the shape α");
  const b = positive(rate, "the rate β");
  if (likelihood !== "poisson" && likelihood !== "exponential") refuse("bad-argument", `the likelihood must be "poisson" or "exponential"`);
  return new Belief({
    family: "gamma", likelihood, params: { shape: a, rate: b }, history,
    mean: a.divide(b),
    variance: a.divide(b.multiply(b)),
    mode: a.compareTo(ONE) >= 0 ? a.subtract(ONE).divide(b) : ZERO,
    get distribution() {
      return gamma(a, b);
    },
    update(data) {
      const xs = exactList(data, likelihood === "poisson" ? "the counts" : "the times");
      const total = xs.reduce((s, v) => s.add(v), ZERO);
      if (likelihood === "poisson") {
        if (xs.some((x) => x.denominator !== 1n || x.numerator < 0n)) refuse("bad-argument", "Poisson counts are whole numbers ≥ 0");
        return gammaPrior(a.add(total), b.add(Q(xs.length)), { likelihood }, [...history, { n: xs.length, sum: total }]);
      }
      if (xs.some((x) => x.numerator < 0n)) refuse("bad-argument", "times are not negative");
      return gammaPrior(a.add(Q(xs.length)), b.add(total), { likelihood }, [...history, { n: xs.length, sum: total }]);
    },
    predictive() {
      if (likelihood !== "poisson") refuse("bad-argument", "the exact predictive here is for Poisson counts");
      // P(X = k) = Γ(α + k)/(Γ(α) k!) · (β/(β+1))^α · (1/(β+1))^k; exact when α is whole
      if (a.denominator !== 1n) refuse("needs-exact", "the predictive probabilities are exact for a whole shape α");
      const p = b.divide(b.add(ONE));
      const q = ONE.divide(b.add(ONE));
      const pa = p.pow(a.numerator);
      let fact = ONE;
      const pmf = (k) => {
        if (k < 0) return ZERO;
        fact = ONE;
        for (let i = 2; i <= k; i++) fact = fact.multiply(Q(i));
        return rising(a, k).divide(fact).multiply(pa).multiply(q.pow(BigInt(k)));
      };
      return Object.freeze({ name: `negative binomial(${a}, ${p})`, min: 0, max: null, pmf, mean: a.divide(b) });
    },
    toString() {
      return `Gamma(shape ${a}, rate ${b}): mean ${show(a.divide(b), 6)}`;
    },
  });
}

/**
 * A normal prior N(μ₀, τ₀²) for the mean of normal data with known variance
 * σ². Precisions add: 1/τ₁² = 1/τ₀² + n/σ², and the posterior mean is the
 * precision-weighted average of μ₀ and the sample mean. Give the prior by its
 * variance (`variance`) or standard deviation (`sd`), and σ² as `dataVariance`.
 * Exact while the variances are rational.
 */
export function normalPrior(mean, { variance, sd, dataVariance, dataSd } = {}, history = []) {
  const m = toQ(mean, "the prior mean");
  const v = variance !== undefined ? positive(variance, "the prior variance") : sd !== undefined ? positive(sd, "the prior sd").pow(2n) : refuse("bad-argument", "give the prior's variance or sd");
  const s2 = dataVariance !== undefined ? positive(dataVariance, "the data variance σ²") : dataSd !== undefined ? positive(dataSd, "the data sd σ").pow(2n) : refuse("bad-argument", "give the known data variance σ² (dataVariance) or sd (dataSd)");
  return new Belief({
    family: "normal", params: { mean: m, variance: v, dataVariance: s2 }, history,
    mean: m, variance: v, mode: m,
    get distribution() {
      const r = exactRoot(v);
      if (!r) refuse("needs-exact", "the normal distribution here needs a rational standard deviation; the mean and variance are exact all the same");
      return normal(m, r);
    },
    update(data) {
      const xs = exactList(data, "the data");
      const n = Q(xs.length);
      const total = xs.reduce((s, x) => s.add(x), ZERO);
      const precision = ONE.divide(v).add(n.divide(s2));
      const v1 = ONE.divide(precision);
      const m1 = v1.multiply(m.divide(v).add(total.divide(s2)));
      return normalPrior(m1, { variance: v1, dataVariance: s2 }, [...history, { n: xs.length, sum: total }]);
    },
    toString() {
      return `Normal(mean ${show(m, 6)}, variance ${show(v, 6)})`;
    },
  });
}
const exactRoot = (q) => {
  const r = (x) => {
    if (x < 0n) return null;
    let y = BigInt(Math.floor(Math.sqrt(Number(x))));
    while (y * y > x) y--;
    while ((y + 1n) * (y + 1n) <= x) y++;
    return y * y === x ? y : null;
  };
  const a = r(q.numerator);
  const b = r(q.denominator);
  return a !== null && b !== null ? new Rational(a, b) : null;
};

/**
 * A Dirichlet prior for the probabilities of k categories, with weights α
 * (a list, or an object keyed by category). Updating adds the counts.
 *   .mean (exact, by category), .marginal(category) a beta belief,
 *   .predictive()  the probability that the next observation is each category
 */
export function dirichletPrior(alpha, history = []) {
  const keys = Array.isArray(alpha) ? alpha.map((_, i) => i) : Object.keys(alpha ?? {});
  if (keys.length < 2) refuse("bad-argument", "a Dirichlet prior needs at least two categories");
  const A = keys.map((k) => positive(Array.isArray(alpha) ? alpha[k] : alpha[k], `α for ${k}`));
  const total = A.reduce((s, a) => s.add(a), ZERO);
  const byKey = (f) => Object.fromEntries(keys.map((k, i) => [k, f(A[i], i)]));
  return new Belief({
    family: "dirichlet", params: byKey((a) => a), history,
    mean: byKey((a) => a.divide(total)),
    variance: byKey((a) => a.multiply(total.subtract(a)).divide(total.multiply(total).multiply(total.add(ONE)))),
    mode: A.every((a) => a.compareTo(ONE) > 0) ? byKey((a) => a.subtract(ONE).divide(total.subtract(Q(keys.length)))) : null,
    get distribution() {
      refuse("bad-argument", "a Dirichlet belief has no single distribution here; use marginal(category)");
    },
    marginal(key) {
      const i = keys.findIndex((k) => String(k) === String(key));
      if (i < 0) refuse("bad-argument", `${key} is not a category`);
      return betaPrior(A[i], total.subtract(A[i]));
    },
    update(counts) {
      const c = keys.map((k) => whole(Array.isArray(counts) ? counts[k] ?? 0 : counts?.[k] ?? 0, `the count for ${k}`));
      const next = keys.map((k, i) => A[i].add(Q(c[i])));
      return dirichletPrior(Array.isArray(alpha) ? next : Object.fromEntries(keys.map((k, i) => [k, next[i]])), [...history, Object.fromEntries(keys.map((k, i) => [k, c[i]]))]);
    },
    predictive() {
      return byKey((a) => a.divide(total));
    },
    interval() {
      refuse("bad-argument", "use marginal(category).interval(level) for a Dirichlet belief");
    },
    toString() {
      return `Dirichlet(${keys.map((k, i) => `${k}: ${A[i]}`).join(", ")})`;
    },
  });
}
