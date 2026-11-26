// The continuous distributions. Density, cumulative probability and quantile
// are certified enclosures (rational intervals that contain the true value),
// built on the special functions of @ratmath/analytic. Each distribution
// also has a floating-point side, `fast`, for previews and simulation.
import { Rational, RationalInterval } from "@ratmath/core";
import { enclose, constants } from "@ratmath/analytic";
import { refuse, toQ, toValue, ZERO, ONE, TWO, HALF, Q, add, sub, mul, div, neg, square, sqrt, lowOf, highOf, interval, isQ, probability } from "./num.js";
import * as F from "./float.js";
import { crossing } from "./solve.js";

const tidy = (v) => (v instanceof RationalInterval && v.low.equals(v.high) ? v.low : v);
const iexp = (v, bits) => tidy(enclose.exp(v, bits));
const ilog = (v, bits) => tidy(enclose.log(v, bits));
const ierf = (v, bits) => tidy(enclose.erf(v, bits));
const ipi = (bits) => constants.pi(bits);
const clamp01 = (v) => {
  if (isQ(v)) return v;
  const lo = v.low.numerator < 0n ? ZERO : v.low;
  const hi = v.high.compareTo(ONE) > 0 ? ONE : v.high;
  return tidy(new RationalInterval(lo, hi));
};
const positive = (v, what) => {
  const q = toQ(v, what);
  if (q.numerator <= 0n) refuse("bad-argument", `${what} must be positive, got ${q}`);
  return q;
};
/** a uniform double strictly between 0 and 1, from 53 random bits */
export const u01 = (rng) => ((rng.u32() >>> 5) * 67108864 + (rng.u32() >>> 6) + 0.5) / 9007199254740992;
const stdNormal = (rng) => F.normalQuantile(u01(rng));
/** Marsaglia and Tsang's sampler for the gamma distribution with rate 1. */
function gammaDraw(rng, shape) {
  if (shape < 1) return gammaDraw(rng, shape + 1) * u01(rng) ** (1 / shape);
  const d = shape - 1 / 3;
  const c = 1 / Math.sqrt(9 * d);
  for (;;) {
    const z = stdNormal(rng);
    const v = (1 + c * z) ** 3;
    if (v <= 0) continue;
    const u = u01(rng);
    if (Math.log(u) < 0.5 * z * z + d - d * v + d * Math.log(v)) return d * v;
  }
}

/**
 * A continuous distribution.
 *   pdf(x)  the density            cdf(x)  P(X ≤ x)        sf(x)  P(X > x)
 *   prob(a, b)  P(a ≤ X ≤ b)       quantile(p)  the x with P(X ≤ x) = p
 *   mean, variance, sd()           (null where they do not exist)
 * Each takes { bits } (default 96 for values, 64 for quantiles) and returns
 * an exact rational where there is one and a certified enclosure otherwise.
 *
 *   fast.pdf, fast.cdf, fast.quantile    the same in floating point
 *   sample(rng), samples(n, rng)         floating-point draws from a seeded generator
 *   sampleExact(rng, { bits })           an exact uniform draw pushed through the certified quantile
 */
export class Continuous {
  constructor(spec) {
    Object.assign(this, spec);
    this.kind = "certified";
    Object.defineProperty(this, "_kept", { value: new Map() });
    this.fast = Object.freeze({
      pdf: spec._fpdf,
      cdf: spec._fcdf,
      quantile: spec._fquantile ?? ((p) => F.quantileBy(spec._fcdf, p, spec.low === null ? -Infinity : spec.low.toNumber(), spec.high === null ? Infinity : spec.high.toNumber(), spec._start ?? 0)),
    });
    Object.freeze(this);
  }
  #where(x) {
    const q = toQ(x, "x");
    if (this.low !== null && q.compareTo(this.low) <= 0) return [q, -1];
    if (this.high !== null && q.compareTo(this.high) >= 0) return [q, 1];
    return [q, 0];
  }
  /** The density at x. */
  pdf(x, { bits = 96 } = {}) {
    const [q, where] = this.#where(x);
    if (where !== 0) {
      const edge = where < 0 ? this.low : this.high;
      if (!q.equals(edge)) return ZERO;
      return this._edge ? this._edge(where, bits) : ZERO;
    }
    return this._pdf(q, bits);
  }
  /** P(X ≤ x). An interval x gives the enclosure over the interval. */
  cdf(x, { bits = 96 } = {}) {
    const v = toValue(x, "x");
    if (!isQ(v)) return interval(lowOf(this.cdf(v.low, { bits })), highOf(this.cdf(v.high, { bits })));
    const [q, where] = this.#where(v);
    if (where < 0) return ZERO;
    if (where > 0) return ONE;
    return clamp01(this._cdf(q, bits));
  }
  /** P(X > x). */
  sf(x, { bits = 96 } = {}) {
    const v = toValue(x, "x");
    if (!isQ(v)) return interval(lowOf(this.sf(v.high, { bits })), highOf(this.sf(v.low, { bits })));
    const [q, where] = this.#where(v);
    if (where < 0) return ONE;
    if (where > 0) return ZERO;
    return clamp01(this._sf ? this._sf(q, bits) : sub(ONE, this._cdf(q, bits)));
  }
  /** P(a ≤ X ≤ b). */
  prob(a, b, { bits = 96 } = {}) {
    const lo = toQ(a, "a");
    const hi = toQ(b, "b");
    if (hi.compareTo(lo) < 0) refuse("bad-argument", "prob(a, b) needs a ≤ b");
    return clamp01(sub(this.cdf(hi, { bits: bits + 2 }), this.cdf(lo, { bits: bits + 2 })));
  }
  /** The x with P(X ≤ x) = p, enclosed to about 2^(−bits) of its size. */
  quantile(p, { bits = 64 } = {}) {
    const level = probability(p, "the level");
    if (level.numerator === 0n) return this.low ?? -Infinity;
    if (level.equals(ONE)) return this.high ?? Infinity;
    // quantiles are asked for again and again in simulations, so the last few hundred are kept
    const key = `${level}|${bits}`;
    const kept = this._kept.get(key);
    if (kept !== undefined) return kept;
    const q = this._quantile
      ? tidy(this._quantile(level, bits))
      : crossing((x, b) => this._cdf(x, b), level, { guess: this.fast.quantile(level.toNumber()), lo: this.low, hi: this.high, bits });
    if (this._kept.size >= 512) this._kept.clear();
    this._kept.set(key, q);
    return q;
  }
  median(options) {
    return this.quantile(HALF, options);
  }
  /** The standard deviation, where the variance exists. */
  sd({ bits = 96 } = {}) {
    return this.variance === null ? null : sqrt(this.variance, bits);
  }
  /** One floating-point draw. */
  sample(rng) {
    return this._draw(rng);
  }
  samples(n, rng) {
    const out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = this._draw(rng);
    return out;
  }
  /**
   * An exact rational draw: a uniform fraction with `bits` random bits is
   * pushed through the certified quantile, and the midpoint of the enclosure
   * (which is within 2^(−bits) of the true quantile, relative to its size) is
   * returned. The same seed gives the same fractions on every machine.
   */
  sampleExact(rng, { bits = 32 } = {}) {
    const u = new Rational((rng.bits(bits) << 1n) | 1n, 1n << BigInt(bits + 1));
    const q = this.quantile(u, { bits: bits + 4 });
    if (isQ(q)) return q;
    const m = q.low.add(q.high).divide(TWO);
    const size = Math.max(Math.abs(m.toNumber()), 2 ** -40);
    const grid = 1n << BigInt(Math.max(0, bits + 4 - Math.ceil(Math.log2(size))));
    return new Rational(m.multiply(new Rational(grid)).round(), grid);
  }
  samplesExact(n, rng, options) {
    return Array.from({ length: n }, () => this.sampleExact(rng, options));
  }
  toString() {
    return `${this.name}(${Object.values(this.params).join(", ")})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

// ---- the normal family ---------------------------------------------------------------------

/** Φ(z) for a rational or interval z. */
function phiCertified(z, bits) {
  const root = enclose.sqrt(HALF, bits + 8);
  return div(add(ONE, ierf(mul(z, root), bits + 4)), TWO);
}
function phiDensity(z, bits) {
  // e^(−z²/2)/√(2π)
  const root = sqrt(mul(ipi(bits + 8), TWO), bits + 8);
  return div(iexp(neg(div(square(z), TWO)), bits + 4), root);
}
/** The certified standard normal quantile. */
function zQuantile(p, bits) {
  return crossing((x, b) => phiCertified(x, b), p, { guess: F.normalQuantile(p.toNumber()), bits });
}

/** The normal distribution with mean μ and standard deviation σ (not the variance). */
export function normal(mu = 0, sigma = 1) {
  const m = toQ(mu, "the mean μ");
  const s = positive(sigma, "the standard deviation σ");
  const [mf, sf] = [m.toNumber(), s.toNumber()];
  const z = (x) => x.subtract(m).divide(s);
  return new Continuous({
    name: "normal", params: { mu: m, sigma: s }, low: null, high: null, mean: m, variance: s.multiply(s),
    _pdf: (x, bits) => div(phiDensity(z(x), bits), s),
    _cdf: (x, bits) => phiCertified(z(x), bits),
    _sf: (x, bits) => phiCertified(z(x).negate(), bits),
    _quantile(p, bits) {
      if (p.equals(HALF)) return m;
      return add(m, mul(zQuantile(p, bits), s));
    },
    _fpdf: (x) => Math.exp(-0.5 * ((x - mf) / sf) ** 2) / (sf * Math.sqrt(2 * Math.PI)),
    _fcdf: (x) => F.phi((x - mf) / sf),
    _fquantile: (p) => mf + sf * F.normalQuantile(p),
    _draw: (rng) => mf + sf * stdNormal(rng),
  });
}

/** The lognormal distribution: log X is normal with mean μ and standard deviation σ. */
export function lognormal(mu = 0, sigma = 1) {
  const m = toQ(mu, "μ");
  const s = positive(sigma, "σ");
  const [mf, sf] = [m.toNumber(), s.toNumber()];
  const z = (x, bits) => div(sub(ilog(x, bits + 8), m), s);
  const s2 = s.multiply(s);
  const mean = iexp(m.add(s2.divide(TWO)), 96);
  return new Continuous({
    name: "lognormal", params: { mu: m, sigma: s }, low: ZERO, high: null, _start: Math.exp(mf),
    mean, variance: mul(sub(iexp(s2, 96), ONE), iexp(m.multiply(TWO).add(s2), 96)),
    _pdf: (x, bits) => div(phiDensity(z(x, bits), bits), s.multiply(x)),
    _cdf: (x, bits) => phiCertified(z(x, bits), bits),
    _quantile: (p, bits) => iexp(add(m, mul(p.equals(HALF) ? ZERO : zQuantile(p, bits + 4), s)), bits + 4),
    _fpdf: (x) => (x <= 0 ? 0 : Math.exp(-0.5 * ((Math.log(x) - mf) / sf) ** 2) / (x * sf * Math.sqrt(2 * Math.PI))),
    _fcdf: (x) => (x <= 0 ? 0 : F.phi((Math.log(x) - mf) / sf)),
    _fquantile: (p) => Math.exp(mf + sf * F.normalQuantile(p)),
    _draw: (rng) => Math.exp(mf + sf * stdNormal(rng)),
  });
}

// ---- the gamma family ----------------------------------------------------------------------

function gammaLike(name, params, shape, rate) {
  const [af, rf] = [shape.toNumber(), rate.toNumber()];
  const lg = F.lgamma(af);
  return new Continuous({
    name, params, low: ZERO, high: null, _start: af / rf,
    mean: shape.divide(rate), variance: shape.divide(rate.multiply(rate)),
    _pdf(x, bits) {
      const y = rate.multiply(x);
      const lead = sub(sub(mul(ilog(y, bits + 8), shape.subtract(ONE)), y), tidy(enclose.lgamma(shape, bits + 8)));
      return mul(iexp(lead, bits + 4), rate);
    },
    _edge(where) {
      if (where > 0) return ZERO;
      const c = shape.compareTo(ONE);
      if (c < 0) refuse("unbounded", `the ${name} density is unbounded at 0`);
      return c === 0 ? rate : ZERO;
    },
    _cdf: (x, bits) => tidy(enclose.gammaP(shape, rate.multiply(x), bits)),
    _fpdf: (x) => (x <= 0 ? 0 : rf * Math.exp((af - 1) * Math.log(rf * x) - rf * x - lg)),
    _fcdf: (x) => F.gammaP(af, rf * x),
    _draw: (rng) => gammaDraw(rng, af) / rf,
  });
}
/** The gamma distribution with the given shape and rate (mean shape/rate). */
export const gamma = (shape, rate = 1) => {
  const a = positive(shape, "the shape");
  const r = positive(rate, "the rate");
  return gammaLike("gamma", { shape: a, rate: r }, a, r);
};
/** The chi-square distribution with k degrees of freedom. */
export const chiSquare = (k) => {
  const d = positive(k, "the degrees of freedom");
  return gammaLike("chiSquare", { df: d }, d.divide(TWO), HALF);
};
/** The exponential distribution with the given rate (mean 1/rate). */
export function exponential(rate = 1) {
  const r = positive(rate, "the rate");
  const rf = r.toNumber();
  return new Continuous({
    name: "exponential", params: { rate: r }, low: ZERO, high: null, mean: r.reciprocal(), variance: r.multiply(r).reciprocal(),
    _pdf: (x, bits) => mul(iexp(r.multiply(x).negate(), bits + 4), r),
    _edge: (where) => (where < 0 ? r : ZERO),
    _cdf: (x, bits) => sub(ONE, iexp(r.multiply(x).negate(), bits)),
    _sf: (x, bits) => iexp(r.multiply(x).negate(), bits),
    _quantile: (p, bits) => div(neg(ilog(ONE.subtract(p), bits + 4)), r),
    _fpdf: (x) => (x < 0 ? 0 : rf * Math.exp(-rf * x)),
    _fcdf: (x) => (x <= 0 ? 0 : -Math.expm1(-rf * x)),
    _fquantile: (p) => -Math.log1p(-p) / rf,
    _draw: (rng) => -Math.log(u01(rng)) / rf,
  });
}
/** The Weibull distribution with the given shape k and scale λ. */
export function weibull(shape, scale = 1) {
  const k = positive(shape, "the shape");
  const lam = positive(scale, "the scale");
  const [kf, lf] = [k.toNumber(), lam.toNumber()];
  // (x/λ)^k, exactly when k is a whole number
  const power = (x, bits) => (k.denominator === 1n && k.numerator < 64n ? x.divide(lam).pow(Number(k.numerator)) : iexp(mul(ilog(x.divide(lam), bits + 8), k), bits + 4));
  const g = (j) => tidy(enclose.gamma(ONE.add(Q(j).divide(k)), 104));
  const g1 = g(1);
  return new Continuous({
    name: "weibull", params: { shape: k, scale: lam }, low: ZERO, high: null, _start: lf,
    mean: mul(g1, lam), variance: mul(sub(g(2), square(g1)), lam.multiply(lam)),
    _pdf(x, bits) {
      const t = power(x, bits);
      return mul(div(mul(t, iexp(neg(t), bits + 4)), x), k);
    },
    _edge(where) {
      if (where > 0) return ZERO;
      const c = k.compareTo(ONE);
      if (c < 0) refuse("unbounded", "the weibull density is unbounded at 0");
      return c === 0 ? lam.reciprocal() : ZERO;
    },
    _cdf: (x, bits) => sub(ONE, iexp(neg(power(x, bits)), bits)),
    _sf: (x, bits) => iexp(neg(power(x, bits)), bits),
    _quantile(p, bits) {
      const t = neg(ilog(ONE.subtract(p), bits + 8)); // (x/λ)^k
      return mul(iexp(div(ilog(t, bits + 8), k), bits + 4), lam);
    },
    _fpdf: (x) => (x <= 0 ? 0 : (kf / lf) * (x / lf) ** (kf - 1) * Math.exp(-((x / lf) ** kf))),
    _fcdf: (x) => (x <= 0 ? 0 : -Math.expm1(-((x / lf) ** kf))),
    _fquantile: (p) => lf * (-Math.log1p(-p)) ** (1 / kf),
    _draw: (rng) => lf * (-Math.log(u01(rng))) ** (1 / kf),
  });
}

// ---- the beta family -----------------------------------------------------------------------

const lnBeta = (a, b, bits) => sub(add(enclose.lgamma(a, bits), enclose.lgamma(b, bits)), enclose.lgamma(a.add(b), bits));
const lnBetaF = (a, b) => F.lgamma(a) + F.lgamma(b) - F.lgamma(a + b);

/** The beta distribution on [0, 1] with parameters a and b. */
export function beta(a, b) {
  const A = positive(a, "a");
  const B = positive(b, "b");
  const [af, bf] = [A.toNumber(), B.toNumber()];
  const s = A.add(B);
  const lb = lnBetaF(af, bf);
  return new Continuous({
    name: "beta", params: { a: A, b: B }, low: ZERO, high: ONE, _start: af / (af + bf),
    mean: A.divide(s), variance: A.multiply(B).divide(s.multiply(s).multiply(s.add(ONE))),
    _pdf(x, bits) {
      const lead = sub(add(mul(ilog(x, bits + 8), A.subtract(ONE)), mul(ilog(ONE.subtract(x), bits + 8), B.subtract(ONE))), lnBeta(A, B, bits + 8));
      return iexp(lead, bits + 4);
    },
    _edge(where) {
      const [own, other] = where < 0 ? [A, B] : [B, A];
      const c = own.compareTo(ONE);
      if (c < 0) refuse("unbounded", `the beta density is unbounded at ${where < 0 ? 0 : 1}`);
      return c === 0 ? other : ZERO;
    },
    _cdf: (x, bits) => tidy(enclose.betaI(x, A, B, bits)),
    _fpdf: (x) => (x <= 0 || x >= 1 ? 0 : Math.exp((af - 1) * Math.log(x) + (bf - 1) * Math.log1p(-x) - lb)),
    _fcdf: (x) => F.betaI(x, af, bf),
    _draw(rng) {
      const g1 = gammaDraw(rng, af);
      return g1 / (g1 + gammaDraw(rng, bf));
    },
  });
}

/** Student's t distribution with ν degrees of freedom. */
export function studentT(df) {
  const v = positive(df, "the degrees of freedom");
  const vf = v.toNumber();
  const half = v.divide(TWO);
  const lead = F.lgamma((vf + 1) / 2) - F.lgamma(vf / 2) - 0.5 * Math.log(vf * Math.PI);
  const tail = (x, bits) => div(tidy(enclose.betaI(v.divide(v.add(x.multiply(x))), half, HALF, bits)), TWO); // P(T > |x|)
  return new Continuous({
    name: "studentT", params: { df: v }, low: null, high: null,
    mean: v.compareTo(ONE) > 0 ? ZERO : null,
    variance: v.compareTo(TWO) > 0 ? v.divide(v.subtract(TWO)) : null,
    _pdf(x, bits) {
      const lg = sub(enclose.lgamma(v.add(ONE).divide(TWO), bits + 8), enclose.lgamma(half, bits + 8));
      const body = mul(ilog(ONE.add(x.multiply(x).divide(v)), bits + 8), v.add(ONE).divide(TWO));
      return div(iexp(sub(lg, body), bits + 4), sqrt(mul(ipi(bits + 8), v), bits + 8));
    },
    _cdf: (x, bits) => (x.numerator === 0n ? HALF : x.numerator > 0n ? sub(ONE, tail(x, bits)) : tail(x, bits)),
    _sf: (x, bits) => (x.numerator === 0n ? HALF : x.numerator > 0n ? tail(x, bits) : sub(ONE, tail(x, bits))),
    _fpdf: (x) => Math.exp(lead - ((vf + 1) / 2) * Math.log1p((x * x) / vf)),
    _fcdf(x) {
      const t = 0.5 * F.betaI(vf / (vf + x * x), vf / 2, 0.5);
      return x > 0 ? 1 - t : t;
    },
    _draw: (rng) => stdNormal(rng) / Math.sqrt((2 * gammaDraw(rng, vf / 2)) / vf),
  });
}

/** The F distribution with d1 and d2 degrees of freedom. */
export function fDist(d1, d2) {
  const a = positive(d1, "the numerator degrees of freedom");
  const b = positive(d2, "the denominator degrees of freedom");
  const [af, bf] = [a.toNumber(), b.toNumber()];
  const [ha, hb] = [a.divide(TWO), b.divide(TWO)];
  const lb = lnBetaF(af / 2, bf / 2);
  const four = Q(4);
  return new Continuous({
    name: "fDist", params: { d1: a, d2: b }, low: ZERO, high: null, _start: 1,
    mean: b.compareTo(TWO) > 0 ? b.divide(b.subtract(TWO)) : null,
    variance: b.compareTo(four) > 0 ? TWO.multiply(b).multiply(b).multiply(a.add(b).subtract(TWO)).divide(a.multiply(b.subtract(TWO).pow(2)).multiply(b.subtract(four))) : null,
    _pdf(x, bits) {
      const ax = a.multiply(x);
      const lead = sub(
        sub(add(mul(ilog(ax, bits + 8), ha), mul(ilog(b, bits + 8), hb)), mul(ilog(ax.add(b), bits + 8), ha.add(hb))),
        lnBeta(ha, hb, bits + 8),
      );
      return div(iexp(lead, bits + 4), x);
    },
    _edge(where) {
      if (where > 0) return ZERO;
      const c = a.compareTo(TWO);
      if (c < 0) refuse("unbounded", "the F density is unbounded at 0");
      return c === 0 ? ONE : ZERO;
    },
    _cdf: (x, bits) => tidy(enclose.betaI(a.multiply(x).divide(a.multiply(x).add(b)), ha, hb, bits)),
    _fpdf: (x) => (x <= 0 ? 0 : Math.exp((af / 2) * Math.log(af * x) + (bf / 2) * Math.log(bf) - ((af + bf) / 2) * Math.log(af * x + bf) - lb) / x),
    _fcdf: (x) => (x <= 0 ? 0 : F.betaI((af * x) / (af * x + bf), af / 2, bf / 2)),
    _draw: (rng) => gammaDraw(rng, af / 2) / af / (gammaDraw(rng, bf / 2) / bf),
  });
}

// ---- the rest ------------------------------------------------------------------------------

/** The uniform distribution on [a, b]. Everything about it is exact. */
export function uniform(a = 0, b = 1) {
  const lo = toQ(a, "a");
  const hi = toQ(b, "b");
  if (hi.compareTo(lo) <= 0) refuse("bad-argument", "uniform(a, b) needs a < b");
  const w = hi.subtract(lo);
  const [lf, wf] = [lo.toNumber(), w.toNumber()];
  return new Continuous({
    name: "uniform", params: { a: lo, b: hi }, low: lo, high: hi, mean: lo.add(hi).divide(TWO), variance: w.multiply(w).divide(Q(12)),
    _pdf: () => w.reciprocal(),
    _edge: () => w.reciprocal(),
    _cdf: (x) => x.subtract(lo).divide(w),
    _quantile: (p) => lo.add(p.multiply(w)),
    _fpdf: (x) => (x < lf || x > lf + wf ? 0 : 1 / wf),
    _fcdf: (x) => Math.min(1, Math.max(0, (x - lf) / wf)),
    _fquantile: (p) => lf + p * wf,
    _draw: (rng) => lf + wf * u01(rng),
  });
}

/** The Cauchy distribution with the given location and scale. It has no mean and no variance. */
export function cauchy(location = 0, scale = 1) {
  const x0 = toQ(location, "the location");
  const g = positive(scale, "the scale");
  const [xf, gf] = [x0.toNumber(), g.toNumber()];
  const z = (x) => x.subtract(x0).divide(g);
  const atanOverPi = (t, bits) => div(tidy(enclose.atan(t, bits + 8)), ipi(bits + 8));
  return new Continuous({
    name: "cauchy", params: { location: x0, scale: g }, low: null, high: null, mean: null, variance: null,
    _pdf: (x, bits) => div(ONE, mul(ipi(bits + 8), g.multiply(ONE.add(z(x).pow(2))))),
    _cdf: (x, bits) => add(HALF, atanOverPi(z(x), bits)),
    _sf: (x, bits) => sub(HALF, atanOverPi(z(x), bits)),
    _fpdf: (x) => 1 / (Math.PI * gf * (1 + ((x - xf) / gf) ** 2)),
    _fcdf: (x) => 0.5 + Math.atan((x - xf) / gf) / Math.PI,
    _fquantile: (p) => xf + gf * Math.tan(Math.PI * (p - 0.5)),
    _draw: (rng) => xf + gf * Math.tan(Math.PI * (u01(rng) - 0.5)),
  });
}
