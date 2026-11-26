// The named discrete distributions. Probabilities are exact rationals for
// rational parameters; the Poisson distribution involves e^(−λ), so its
// probabilities are certified enclosures.
import { Rational, RationalInterval } from "@ratmath/core";
import { enclose } from "@ratmath/analytic";
import { binomial as choose } from "@ratmath/discrete";
import { refuse, toQ, ZERO, ONE, Q, sqrt, whole, probability, add, sub, lowOf, highOf, interval, toNumber } from "./num.js";
import { Finite } from "./finite.js";

const B = (n) => BigInt(n);
const wholeK = (k) => {
  const q = typeof k === "number" ? (Number.isInteger(k) ? k : NaN) : typeof k === "bigint" ? Number(k) : (() => {
    const r = toQ(k, "k");
    return r.denominator === 1n ? Number(r.numerator) : NaN;
  })();
  return q;
};

/**
 * A distribution on whole numbers.
 *   pmf(k)  P(X = k)        cdf(k)  P(X ≤ k)        sf(k)  P(X > k)
 *   prob(a, b)  P(a ≤ X ≤ b)        quantile(p)  the least k with P(X ≤ k) ≥ p
 *   mean, variance, sd(), mode      sample(rng), samples(n, rng)
 *   toFinite()  the same distribution as a finite table (bounded support only)
 * `kind` is "exact" when the probabilities are rationals and "certified"
 * when they are enclosures.
 */
export class Discrete {
  constructor(spec) {
    Object.assign(this, spec);
    Object.freeze(this);
  }
  /** P(X = k); 0 off the support. */
  pmf(k, options) {
    const n = wholeK(k);
    if (!Number.isInteger(n) || n < this.min || (this.max !== null && n > this.max)) return ZERO;
    return this._pmf(n, options);
  }
  /** P(X ≤ k). */
  cdf(k, options) {
    const q = typeof k === "number" ? Math.floor(k) : Number(toQ(k, "k").floor());
    if (q < this.min) return ZERO;
    if (this.max !== null && q >= this.max) return ONE;
    if (this._cdf) return this._cdf(q, options);
    let s = ZERO;
    for (let i = this.min; i <= q; i++) s = add(s, this._pmf(i, options));
    return s;
  }
  /** P(X > k). */
  sf(k, options) {
    return sub(ONE, this.cdf(k, options));
  }
  /** P(a ≤ X ≤ b). */
  prob(a, b, options) {
    const lo = Math.max(this.min, Math.ceil(typeof a === "number" ? a : toNumber(toQ(a, "a"))));
    const hi0 = Math.floor(typeof b === "number" ? b : toNumber(toQ(b, "b")));
    const hi = this.max === null ? hi0 : Math.min(this.max, hi0);
    let s = ZERO;
    for (let i = lo; i <= hi; i++) s = add(s, this._pmf(i, options));
    return s;
  }
  /** The least k with P(X ≤ k) ≥ p. */
  quantile(p, options) {
    const level = probability(p, "the level");
    if (level.equals(ONE)) return this.max === null ? Infinity : this.max;
    let s = ZERO;
    for (let k = this.min; ; k++) {
      if (this.max !== null && k >= this.max) return this.max;
      s = add(s, this._pmf(k, options));
      if (lowOf(s).compareTo(level) >= 0) return k;
      if (highOf(s).compareTo(level) >= 0) refuse("too-close", `P(X ≤ ${k}) cannot be told apart from ${level} at this precision; pass more bits`);
      if (k - this.min > 10000000) refuse("budget", "the quantile lies too far out");
    }
  }
  median(options) {
    return this.quantile(Q(1, 2), options);
  }
  /** The standard deviation: exact when the variance is a perfect square, otherwise certified. */
  sd({ bits = 96 } = {}) {
    return sqrt(this.variance, bits);
  }
  /**
   * One draw from a seeded generator. Where the distribution is built from
   * coin tosses or draws the sampler does exactly that; otherwise it inverts
   * the cumulative probabilities at a uniform point with 64 random bits.
   */
  sample(rng) {
    if (this._sample) return this._sample(rng);
    const u = new Rational((rng.bits(64) << 1n) | 1n, 1n << 65n);
    let s = ZERO;
    for (let k = this.min; ; k++) {
      if (this.max !== null && k >= this.max) return this.max;
      s = add(s, this._pmf(k, { bits: 80 }));
      if (lowOf(s).compareTo(u) >= 0) return k;
      if (highOf(s).compareTo(u) >= 0) {
        // undecided at 80 bits: settle it with more
        let t = ZERO;
        for (let i = this.min; i <= k; i++) t = add(t, this._pmf(i, { bits: 160 }));
        if (lowOf(t).compareTo(u) >= 0) return k;
      }
    }
  }
  samples(n, rng) {
    return Array.from({ length: n }, () => this.sample(rng));
  }
  /** The distribution as a finite table, for a bounded support with exact probabilities. */
  toFinite() {
    if (this.max === null) refuse("unbounded", `the ${this.name} distribution has infinitely many values; it has no finite table`);
    if (this.kind !== "exact") refuse("not-exact", `the ${this.name} distribution has certified, not exact, probabilities`);
    const rows = [];
    for (let k = this.min; k <= this.max; k++) rows.push([k, this._pmf(k)]);
    return new Finite(rows);
  }
  toString() {
    return `${this.name}(${Object.values(this.params).join(", ")})`;
  }
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }
}

/** One trial: 1 with probability p, otherwise 0. */
export function bernoulli(p) {
  const q = probability(p, "p");
  return new Discrete({
    name: "bernoulli", params: { p: q }, kind: "exact", min: 0, max: 1,
    mean: q, variance: q.multiply(ONE.subtract(q)), mode: q.compareTo(Q(1, 2)) > 0 ? 1 : 0,
    _pmf: (k) => (k === 1 ? q : ONE.subtract(q)),
    _sample: (rng) => (rng.bool(q) ? 1 : 0),
  });
}

/** The number of successes in n independent trials, each succeeding with probability p. */
export function binomial(n, p) {
  const N = whole(n, "the number of trials n");
  const q = probability(p, "p");
  const c = ONE.subtract(q);
  const mode = Number(q.multiply(Q(N + 1)).floor());
  return new Discrete({
    name: "binomial", params: { n: N, p: q }, kind: "exact", min: 0, max: N,
    mean: q.multiply(Q(N)), variance: q.multiply(c).multiply(Q(N)), mode: Math.min(N, mode),
    _pmf: (k) => new Rational(choose(N, k)).multiply(q.pow(k)).multiply(c.pow(N - k)),
    _sample(rng) {
      let s = 0;
      for (let i = 0; i < N; i++) if (rng.bool(q)) s++;
      return s;
    },
  });
}

/**
 * The number of the trial on which the first success comes (1, 2, 3, …).
 * With { failures: true }, the number of failures before it (0, 1, 2, …).
 */
export function geometric(p, { failures = false } = {}) {
  const q = probability(p, "p");
  if (q.numerator === 0n) refuse("bad-argument", "a geometric distribution needs p > 0");
  const c = ONE.subtract(q);
  const shift = failures ? 0 : 1;
  return new Discrete({
    name: "geometric", params: failures ? { p: q, failures: true } : { p: q }, kind: "exact", min: shift, max: null,
    mean: failures ? c.divide(q) : q.reciprocal(), variance: c.divide(q.multiply(q)), mode: shift,
    _pmf: (k) => c.pow(k - shift).multiply(q),
    _cdf: (k) => ONE.subtract(c.pow(k - shift + 1)),
    _sample(rng) {
      let k = shift;
      while (!rng.bool(q)) k++;
      return k;
    },
  });
}

/** The number of failures before the r-th success. */
export function negativeBinomial(r, p) {
  const R = whole(r, "the number of successes r", { min: 1 });
  const q = probability(p, "p");
  if (q.numerator === 0n) refuse("bad-argument", "a negative binomial distribution needs p > 0");
  const c = ONE.subtract(q);
  return new Discrete({
    name: "negativeBinomial", params: { r: R, p: q }, kind: "exact", min: 0, max: null,
    mean: c.multiply(Q(R)).divide(q), variance: c.multiply(Q(R)).divide(q.multiply(q)),
    mode: R > 1 ? Number(c.multiply(Q(R - 1)).divide(q).floor()) : 0,
    _pmf: (k) => new Rational(choose(k + R - 1, k)).multiply(q.pow(R)).multiply(c.pow(k)),
    _sample(rng) {
      let fails = 0;
      for (let wins = 0; wins < R; ) {
        if (rng.bool(q)) wins++;
        else fails++;
      }
      return fails;
    },
  });
}

/** The number of marked items among n drawn without replacement from N items of which K are marked. */
export function hypergeometric(N, K, n) {
  const pop = whole(N, "the population size N");
  const marked = whole(K, "the number marked K");
  const drawn = whole(n, "the number drawn n");
  if (marked > pop || drawn > pop) refuse("bad-argument", "K and n cannot exceed N");
  const total = choose(pop, drawn);
  const lo = Math.max(0, drawn - (pop - marked));
  const hi = Math.min(drawn, marked);
  const mean = pop === 0 ? ZERO : Q(drawn * marked, pop);
  const variance = pop <= 1 ? ZERO : mean.multiply(Q(pop - marked, pop)).multiply(Q(pop - drawn, pop - 1));
  return new Discrete({
    name: "hypergeometric", params: { N: pop, K: marked, n: drawn }, kind: "exact", min: lo, max: hi,
    mean, variance, mode: Math.min(hi, Math.max(lo, Math.floor(((drawn + 1) * (marked + 1)) / (pop + 2)))),
    _pmf: (k) => new Rational(choose(marked, k) * choose(pop - marked, drawn - k), total),
    _sample(rng) {
      let left = marked;
      let got = 0;
      for (let i = 0; i < drawn; i++) {
        if (rng.int(pop - i) < left) {
          got++;
          left--;
        }
      }
      return got;
    },
  });
}

/** A whole number from a to b, each equally likely. */
export function discreteUniform(a, b) {
  const lo = whole(a, "a", { min: -Infinity });
  const hi = whole(b, "b", { min: -Infinity });
  if (hi < lo) refuse("bad-argument", "discreteUniform(a, b) needs a ≤ b");
  const n = hi - lo + 1;
  return new Discrete({
    name: "discreteUniform", params: { a: lo, b: hi }, kind: "exact", min: lo, max: hi,
    mean: Q(lo + hi, 2), variance: Q(B(n) * B(n) - 1n, 12), mode: lo,
    _pmf: () => Q(1, n),
    _cdf: (k) => Q(k - lo + 1, n),
    _sample: (rng) => rng.between(lo, hi),
  });
}

/**
 * The Poisson distribution with mean λ. Its probabilities are e^(−λ)·λᵏ/k!,
 * which are not rational, so they are returned as certified enclosures;
 * pass { bits } to pmf, cdf and the rest for a different precision.
 */
export function poisson(lambda) {
  const L = toQ(lambda, "the mean λ");
  if (L.numerator <= 0n) refuse("bad-argument", `a Poisson distribution needs λ > 0, got ${L}`);
  const factorials = [1n];
  const fact = (k) => {
    while (factorials.length <= k) factorials.push(factorials[factorials.length - 1] * B(factorials.length));
    return factorials[k];
  };
  const eNeg = (bits) => enclose.exp(L.negate(), bits + 8);
  const scaled = (coefficient, bits) => {
    const e = eNeg(bits);
    return interval(e.low.multiply(coefficient), e.high.multiply(coefficient));
  };
  return new Discrete({
    name: "poisson", params: { lambda: L }, kind: "certified", min: 0, max: null,
    mean: L, variance: L, mode: Number(L.floor()),
    _pmf: (k, { bits = 96 } = {}) => scaled(L.pow(k).divide(new Rational(fact(k))), bits),
    _cdf(k, { bits = 96 } = {}) {
      let s = ZERO;
      let t = ONE;
      for (let i = 0; i <= k; i++) {
        s = s.add(t);
        t = t.multiply(L).divide(Q(i + 1));
      }
      const v = scaled(s, bits);
      return v.high.compareTo(ONE) > 0 ? new RationalInterval(v.low, ONE) : v;
    },
  });
}
