// Closed-form power and sample sizes for the z test of a mean and the
// one-proportion z test (normal approximation), the formulas textbooks use.
// Quantiles and normal probabilities are certified enclosures, so the power is
// an interval certain to contain the true value, and a sample size is the
// smallest whole n whose power (or margin) is proven to meet the target.
import { Rational, RationalInterval } from "@ratmath/core";
import { refuse, toQ, Q, ZERO, ONE, TWO, add, sub, mul, div, neg, sqrt, lowOf, highOf, interval, isQ, probability, whole } from "./num.js";
import { normal } from "./continuous.js";
import { alternativeOf, levelOf } from "./inference.js";

const Z = normal(0, 1);
const tidy = (v) => (v instanceof RationalInterval && v.low.equals(v.high) ? v.low : v);
const clamp01 = (v) => {
  const lo = lowOf(v).numerator < 0n ? ZERO : lowOf(v);
  const hi = highOf(v).compareTo(ONE) > 0 ? ONE : highOf(v);
  return tidy(interval(lo, hi));
};
const alphaOf = (alpha) => {
  const a = probability(alpha, "alpha");
  if (a.numerator === 0n || a.equals(ONE)) refuse("bad-argument", "alpha must be strictly between 0 and 1");
  return a;
};
const positive = (v, what) => {
  const q = toQ(v, what);
  if (q.numerator <= 0n) refuse("bad-argument", `${what} must be positive`);
  return q;
};

/**
 * Power of the z test at n: the probability of rejecting H₀ when the truth is
 * the alternative value. Returns an enclosure.
 */
function zPower({ mu0, mu, sigma, n, alpha, alternative, bits }) {
  // δ = (μ − μ0)√n/σ: the shift of the statistic's mean
  const delta = div(mul(sub(mu, mu0), sqrt(Q(n), bits)), sigma);
  return powerFromShift(delta, alpha, alternative, bits, (x) => x);
}

/** P(reject) for a statistic that is N(shift·…) — generic: rejection cut(s) c and the statistic's law under the truth. */
function powerFromShift(delta, alpha, alternative, bits, scale) {
  if (alternative === "greater") {
    const c = Z.quantile(ONE.subtract(alpha), { bits });
    return clamp01(Z.sf(sub(scale(c), delta), { bits }));
  }
  if (alternative === "less") {
    const c = Z.quantile(ONE.subtract(alpha), { bits });
    return clamp01(Z.cdf(sub(neg(scale(c)), delta), { bits }));
  }
  const c = Z.quantile(ONE.subtract(alpha.divide(TWO)), { bits });
  return clamp01(add(Z.sf(sub(scale(c), delta), { bits }), Z.cdf(sub(neg(scale(c)), delta), { bits })));
}

/** Power of the one-proportion z test (normal approximation): cut-offs from p0, spread under p. */
function proportionPower({ p0, p, n, alpha, alternative, bits }) {
  const se0 = sqrt(p0.multiply(ONE.subtract(p0)).divide(Q(n)), bits);
  const se1 = sqrt(p.multiply(ONE.subtract(p)).divide(Q(n)), bits);
  // reject when p̂ > p0 + c·se0 (and/or below p0 − c·se0); p̂ ~ N(p, se1²)
  const delta = div(sub(p, p0), se1);
  const ratio = div(se0, se1);
  return powerFromShift(delta, alpha, alternative, bits, (c) => mul(c, ratio));
}

/**
 * Closed-form power.
 *   power({ test: "z", mu0, mu, sigma, n, alpha = 0.05, alternative = "two-sided" })
 *   power({ test: "proportion", p0, p, n, alpha, alternative })   normal approximation
 * Returns { power, test, n, alpha, alternative, kind }: `power` is a Rational
 * or an interval certain to contain it.
 */
export function power(spec = {}) {
  const { test = "z", alpha = "0.05", alternative = "two-sided", bits = 64 } = spec;
  const a = alphaOf(alpha);
  const alt = alternativeOf(alternative);
  const n = whole(spec.n, "n", { min: 1 });
  if (test === "z") {
    const value = zPower({ mu0: toQ(spec.mu0 ?? 0, "mu0"), mu: toQ(spec.mu, "mu"), sigma: positive(spec.sigma, "sigma"), n, alpha: a, alternative: alt, bits });
    return { power: value, test: "z", n, alpha: a, alternative: alt, kind: isQ(value) ? "exact" : "certified" };
  }
  if (test === "proportion") {
    const p0 = openProbability(spec.p0, "p0");
    const p = openProbability(spec.p, "p");
    const value = proportionPower({ p0, p, n, alpha: a, alternative: alt, bits });
    return { power: value, test: "proportion", n, alpha: a, alternative: alt, kind: "normal approximation", note: "the normal approximation to the binomial, as textbooks use it" };
  }
  refuse("unsupported", `closed-form power is for test "z" or "proportion"; for others use powerCurve, which simulates`);
}

function openProbability(v, what) {
  const p = probability(v, what);
  if (p.numerator === 0n || p.equals(ONE)) refuse("bad-argument", `${what} must be strictly between 0 and 1`);
  return p;
}

/** The smallest whole number ≥ v, settled from enclosures at rising precision. */
function ceilOf(make) {
  for (let bits = 64; bits <= 4096; bits *= 2) {
    const v = make(bits);
    const lo = BigInt(lowOf(v).ceil());
    const hi = BigInt(highOf(v).ceil());
    if (lo === hi) return { n: lo, value: v };
  }
  refuse("precision", "the sample size is too close to a whole number to settle");
}

/** A sample size: a Number while it is exact as one, else the exact bigint. */
const count = (n) => (n <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(n) : n);

/**
 * Sample sizes.
 *   sampleSize({ margin, p = 1/2, level = 0.95 })            a proportion: n ≥ z²p(1 − p)/m²
 *   sampleSize({ margin, sigma, level })                     a mean with known σ: n ≥ (zσ/m)²
 *   sampleSize({ power, mu0, mu, sigma, alpha, alternative }) the z test reaching that power
 *   sampleSize({ power, test: "proportion", p0, p, alpha, alternative })
 * Returns { n, ... }: the smallest whole n that is proven to meet the target,
 * a Number, or a bigint when it is above 2^53 − 1 (a Number would round it).
 * With no effect (mu = mu0, or p = p0) the power is alpha for every n, so
 * n = 1 when the target power is at most alpha and a refusal otherwise.
 */
export function sampleSize(spec = {}) {
  if (spec.power !== undefined) return sampleSizeForPower(spec);
  const margin = positive(spec.margin, "margin");
  const conf = levelOf(spec.level ?? "0.95");
  const zq = (bits) => Z.quantile(ONE.subtract(ONE.subtract(conf).divide(TWO)), { bits });
  if (spec.sigma !== undefined) {
    const sigma = positive(spec.sigma, "sigma");
    const { n, value } = ceilOf((bits) => {
      const t = div(mul(zq(bits), sigma), margin);
      return mul(t, t);
    });
    return { n: count(n), bound: value, kind: "mean", formula: "n ≥ (z·σ/m)²", level: conf, margin, sigma };
  }
  const p = spec.p === undefined ? Q(1, 2) : probability(spec.p, "p");
  const { n, value } = ceilOf((bits) => {
    const t = div(zq(bits), margin);
    return mul(mul(t, t), p.multiply(ONE.subtract(p)));
  });
  return { n: count(n), bound: value, kind: "proportion", formula: "n ≥ z²·p(1 − p)/m²", level: conf, margin, p };
}

/** With no effect the power is the size alpha for every n: n = 1 reaches a target up to alpha, nothing reaches more. */
function noEffect(target, alpha, at, info, why) {
  if (target.compareTo(alpha) > 0) refuse("bad-argument", `${why}: the power is alpha = ${alpha} for every sample size, so it never reaches ${target}`);
  return { n: 1, power: at(1), target, ...info, kind: "smallest n with proven power", note: `${why}: the power is alpha for every n` };
}

function sampleSizeForPower(spec) {
  const target = openProbability(spec.power, "power");
  const { test = "z", alpha = "0.05", alternative = "two-sided", bits = 96 } = spec;
  const a = alphaOf(alpha);
  const alt = alternativeOf(alternative);
  let at;
  if (test === "z") {
    const mu0 = toQ(spec.mu0 ?? 0, "mu0"), mu = toQ(spec.mu, "mu"), sigma = positive(spec.sigma, "sigma");
    at = (n) => zPower({ mu0, mu, sigma, n, alpha: a, alternative: alt, bits });
    if (mu.equals(mu0)) return noEffect(target, a, at, { test, alpha: a, alternative: alt }, "mu equals mu0");
  } else if (test === "proportion") {
    const p0 = openProbability(spec.p0, "p0"), p = openProbability(spec.p, "p");
    at = (n) => proportionPower({ p0, p, n, alpha: a, alternative: alt, bits });
    if (p.equals(p0)) return noEffect(target, a, at, { test, alpha: a, alternative: alt }, "p equals p0");
  } else refuse("unsupported", `sample sizes for power are for test "z" or "proportion"`);
  // power rises with n: double until it is proven to reach the target, then bisect
  const reaches = (n) => {
    const v = at(n);
    if (lowOf(v).compareTo(target) >= 0) return true;
    if (highOf(v).compareTo(target) < 0) return false;
    refuse("precision", `the power at n = ${n} is too close to the target to settle`);
  };
  let hi = 1;
  while (!reaches(hi)) {
    hi *= 2;
    if (hi > 2 ** 40) refuse("budget", "the power does not reach the target for any reasonable n (wrong direction for the alternative?)");
  }
  let lo = Math.floor(hi / 2);
  while (hi - lo > 1) {
    const mid = Math.floor((lo + hi) / 2);
    if (reaches(mid)) hi = mid;
    else lo = mid;
  }
  return { n: hi, power: at(hi), target, test, alpha: a, alternative: alt, kind: "smallest n with proven power" };
}
